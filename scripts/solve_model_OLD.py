#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import pypsa
from __helpers import replace_su, average_every_nhours
import numpy as np
import os
from pypsa.descriptors import get_switchable_as_dense as as_dense
import math

# In[2]:


def solve_rolling_horizon(m, o, time_aggregation, horizon):
    
    o.storage_units.cyclic_state_of_charge=False
    o.stores.e_cyclic=False

    o.stores_t.e_min_pu = (
        (m.stores_t.e/m.stores.e_nom)
        .shift(1)
        .resample(time_aggregation).asfreq()
        .shift(-1)
        .fillna(0)
    ).clip(upper=.999).reindex(o.snapshots).fillna(0)

    o.stores.e_initial = m.stores_t.e.iloc[-1]

    o.madd("Generator", 
           o.buses.query("carrier == 'AC'").index + " load", 
           bus = o.buses.query("carrier == 'AC'").index, 
           p_nom = 1e6, 
           marginal_cost = 10e3
          )

    o.generators_t["r"] = pd.DataFrame(
        index=o.snapshots,
        columns=o.generators.index,
        data=np.nan,
    )
    
    o.links_t["r"] = pd.DataFrame(
        index=o.snapshots,
        columns=o.links.index,
        data=np.nan,
    )
    
    o.buses_t["reserve_price"] = pd.DataFrame(
        index=o.snapshots,
        columns=o.buses.index,
        data=np.nan,
    )
    
    for i in range(len(o.snapshots)//horizon):

        start= i* horizon
        end = min((i+1)*horizon, len(o.snapshots))
        
        print("optimizing time period between " + str(o.snapshots[start]) + " and " + str(o.snapshots[end-1]))
        
        snapshots = o.snapshots[start:end]
              
        o.optimize(solver_name=solver_name,
                   snapshots = snapshots,
                   assign_all_duals=True,
                   extra_functionality=balancing_market,
                   **solver_options
              )
        
        reserve_generators = o.model.solution["Generator-r"].to_dataframe()["Generator-r"].unstack()
        o.generators_t.r.loc[reserve_generators.index, reserve_generators.columns] = reserve_generators
        
        reserve_storage = o.model.solution["Link-r"].to_dataframe()["Link-r"].unstack()
        o.links_t.r.loc[reserve_storage.index, reserve_storage.columns] = reserve_storage
        
        reserve_price = o.model.dual["Bus-reserve_balance"].to_dataframe()["Bus-reserve_balance"].unstack(0)
        o.buses_t.reserve_price.loc[reserve_price.index, reserve_price.columns] =  reserve_price
        
        
        o.stores.e_initial = o.stores_t.e.loc[snapshots].iloc[-1]

def get_ordc_integral(excess_reserve, mean, std):
    
    return (
        -0.5*excess_reserve 
        + 0.5*(
            (excess_reserve - mean)*math.erf((excess_reserve-mean)/(std*math.sqrt(2)))
            + math.sqrt(2)*std*math.exp(-(excess_reserve-mean)**2/(2*std**2))/math.sqrt(math.pi)
        )
    )
    


def determine_cutoff(zone):
    cost_function = pd.Series([get_ordc_integral(i, **ordc_parameters.loc[zone]) for i in range(0,40000,10)], 
                              index = range(0,40000,10))

    return cost_function.diff()[cost_function.diff()/pd.Series(cost_function.index).diff().median()<=-0.00001].index[-1]
    

def get_break_points(n_ints, zone, cut_off):
    
    break_points = pd.DataFrame(index=range(n_ints+1))
    break_points["cap_break_points"] = list(np.arange(0, cut_off[zone], cut_off[zone]/n_ints)) + [2*cut_off[zone]]
    break_points["welfare_break_points"] = [get_ordc_integral(i, **ordc_parameters.loc[zone]) for i in break_points.cap_break_points]
    break_points.set_index(pd.Series(zone,index=break_points.index), inplace=True,append=True)
    break_points = break_points.reorder_levels([1,0])
    
    return break_points


def get_linear_function(p, q):
    m = (p[1] - q[1])/(p[0] - q[0])
    b = (p[0]*q[1] - p[1]*q[0])/(p[0] - q[0])
    return m, b
    
def get_line_parameters(ordc_parameters, n_ints):
    cut_off = pd.Series(index=ordc_parameters.index, data=[determine_cutoff(i) for i in ordc_parameters.index])

    break_points = pd.DataFrame()
    for zone in ordc_parameters.index:

        break_points = pd.concat([break_points, get_break_points(n_ints, zone, cut_off)])

    line_parameters = pd.DataFrame()
    for zone in ordc_parameters.index:
        j = 0    
        for i in break_points.loc[zone].index[:-1]:
            m, b = get_linear_function(break_points.loc[zone, i], break_points.loc[zone, i+1])
            line_parameters = pd.concat([line_parameters, pd.Series(index=["a", "b"], data=[m, b], name=(zone, j))],axis=1)
            j+=1

    return line_parameters.T.set_index((i for i in line_parameters.columns))


def balancing_market(n, snapshots):
    
    m = n.model

    reserve_gens = n.generators.loc[n.generators.bus.map(n.buses.carrier) == "electricity"]
    reserve_gens = reserve_gens.query("carrier in @reserve_participation_generators")

    coord_reserve_gens = reserve_gens.index
    coord_reserve_gens.name = "Generator-fix"

    gen_r = m.add_variables(lower=0, name="Generator-r", coords = [snapshots, coord_reserve_gens])

    lhs = m.constraints["Generator-fix-p-upper"].lhs
    rhs = m.constraints["Generator-fix-p-upper"].rhs
    m.remove_constraints("Generator-fix-p-upper")
    m.add_constraints(lhs + gen_r <= rhs, name= "Generator-fix-p-upper")

    committable = reserve_gens.query("committable == True").copy()

    committable_grouper = pd.Series(committable.index, committable.index)
    committable_grouper.name="Generator-com"

    lhs = m.constraints["Generator-com-p-upper"].lhs
    rhs = m.constraints["Generator-com-p-upper"].rhs

    lhs += gen_r.loc[:, committable.index].groupby(committable_grouper.to_xarray()).sum()
    m.remove_constraints("Generator-com-p-upper")
    m.add_constraints( lhs <= rhs, name="Generator-com-p-upper")

    gen_status = m["Generator-status"]

    reserve_gen_status = gen_status.loc[:, committable.index]

    m.add_constraints(
        gen_r.loc[:, committable.index].groupby(committable_grouper.to_xarray()).sum()
        <= (
            reserve_gen_status*
            (
                committable.groupby(committable_grouper.to_xarray()).ramp_limit_up.sum()
                .multiply(delivery_time_reserves)
                .multiply(committable.groupby(committable_grouper.to_xarray()).p_nom.sum())     
            )
        ),
        name = "Generator-com-ramp_upper"
    )

    res_shedding = m.add_variables(
        lower=0, 
        name="Bus-reserve_shedding", 
        coords=[snapshots, reserve_requirements.index]
    )

    storage = n.links.loc[n.links.type == "discharging"]

    storage_link_r = m.add_variables(lower = 0, name="Link-r", coords = m["Link-p"].loc[:, storage.index].coords)

    link_p = m["Link-p"]

    m.add_constraints(
        link_p.loc[:, storage.index] + storage_link_r <= n.links.p_nom.to_xarray(),
        name = "Link-capacity_upper"
    );

    store_e = m["Store-e"]

    store_link_map = n.stores.set_index("bus",append=True).reset_index("Store").reindex(storage.bus0)["Store"]
    store_link_map.index = storage.index
    store_link_map.index.name = "Link"
    store_link_map.name="Store"

    m.add_constraints(
        storage_link_r.groupby(store_link_map.to_xarray()).sum() <= store_e,
        name="Link-reserve_reservoir_constraint"
    )

    gen_grouper = reserve_gens.bus
    gen_grouper.name = "Bus"

    storage_grouper = storage.bus1
    storage_grouper.name = "Bus"

    if balancing_market_design == "ORDC":

        excess_r = m.add_variables(
            lower = 0,
            upper= 40,
            name = "Bus-excess_r",
            coords = [snapshots, reserve_requirements.index]
        )

        m.add_constraints(
            gen_r.groupby(gen_grouper.to_xarray()).sum() 
            + storage_link_r.groupby(storage_grouper.to_xarray()).sum()
            + res_shedding
            - excess_r
            >= 0,
            name="Bus-reserve_balance"
        )

        
        a = linear_ordc_approximation["a"]
        b = linear_ordc_approximation["b"]

        a.index.names = ["Bus", "break_point"]
        b.index.names = ["Bus", "break_point"]


        ordc_cost_term = m.add_variables(coords=excess_r.coords, name="Bus-ORDC_cost")
        
        m.add_constraints(
            ordc_cost_term >= (a.to_xarray()*excess_r + b.to_xarray())*VOLL,
            name="Bus-ORDC"
        )

        obj = m.objective
        m.add_objective(obj + res_shedding.sum()*1e5 + ordc_cost_term.sum(), overwrite=True) 

    else:

        m.add_constraints(
            gen_r.groupby(gen_grouper.to_xarray()).sum() 
            + storage_link_r.groupby(storage_grouper.to_xarray()).sum()
            + res_shedding
            >= reserve_requirements.reindex(n.buses.query("carrier == 'electricity'").index).dropna().to_xarray(),
            name="Bus-reserve_balance"
        )

        obj = m.objective

        m.add_objective(obj + res_shedding.sum()*1e5, overwrite=True)


target_year = int(snakemake.params.ty)

solver_name = snakemake.config["solving"]["solver"]["name"]
options = snakemake.config["solving"]["solver"]["options"]
solver_options = snakemake.config["solving"]["solver_options"][options]

delivery_time_reserves = snakemake.config["delivery_time_reserves"] # in minutes
delivery_time_reserves = delivery_time_reserves/60

ordc_parameters = pd.read_hdf(snakemake.input.ordc_parameters, "ordc_parameters")

years = snakemake.params.years

VOLL = snakemake.config["VOLL"]

solved_network = snakemake.output.solved_network

n = pypsa.Network(snakemake.input.network)

reserve_requirements = pd.read_excel("data/pemmdb.xlsx", "Reserve Requirements",index_col=[0,1])[["FCR (MW)", "FRR (MW)"]].sum(axis=1)#.loc[:, target_year]
reserve_requirements = reserve_requirements.unstack(1).reindex(years, axis=1).interpolate(axis=1)[target_year]
reserve_requirements.index.name="Bus"
reserve_requirements = reserve_requirements[reserve_requirements>0]

reserve_participation_generators = snakemake.config["balancing"]["participation"]

for su in n.storage_units.index:
    replace_su(n, su);

storage_preopt_aggregation = snakemake.config["storage_preopt_aggregation"]


n.madd("Generator", 
       n.buses.query("carrier == 'electricity'").index + " load-shedding", 
       bus = n.buses.query("carrier == 'electricity'").index, 
       p_nom = 1e6, 
       marginal_cost = VOLL
      )
      
m = average_every_nhours(n, storage_preopt_aggregation)
m.stores.e_cyclic = True
m.generators["p_min_pu"] = 0
m.optimize(
    solver_name="highs",
)

dispatchable = ['CCGT', 'OCGT', 'oil', 'biomass', 'other', 'lignite', 'nuclear', 'coal']

n.generators.loc[n.generators.query("carrier in @dispatchable").index, "committable"] = True
n.generators.shut_down_cost = n.generators.start_up_cost


solve_rolling_horizon(m, n, "h", int(storage_preopt_aggregation[:-1]))

dirname = os.path.dirname(solved_network)

if not os.path.isdir(dirname):
    os.makedirs(dirname)

n.export_to_netcdf(solved_network)

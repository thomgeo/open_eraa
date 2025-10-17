#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import pypsa
from __helpers import replace_su, average_every_nhours
import numpy as np
import os
import math
from pypsa.descriptors import get_switchable_as_dense as as_dense

def solve_rolling_horizon(m, o, time_aggregation, horizon):

    horizon = int(storage_preopt_aggregation.split("h")[0])
    
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
    
    for i in range(len(o.snapshots)//horizon):

        start= i* horizon
        end = min((i+1)*horizon, len(o.snapshots))
        
        print("optimizing time period between " + str(o.snapshots[start]) + " and " + str(o.snapshots[end-1]))
        
        snapshots = o.snapshots[start:end]
              
        o.optimize(solver_name=solver_name,
                   snapshots = snapshots,
                   assign_all_duals=True,
                   extra_functionality=extra_functionality,
                   linearized_unit_commitment = True,
                   **solver_options
              )
        
        
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


# In[6]:


def get_break_points(n_ints, zone, cut_off):
    
    break_points = pd.DataFrame(index=range(n_ints+1))
    break_points["cap_break_points"] = list(np.arange(0, cut_off[zone], cut_off[zone]/n_ints)) + [2*cut_off[zone]]
    break_points["welfare_break_points"] = [get_ordc_integral(i, **ordc_parameters.loc[zone]) for i in break_points.cap_break_points]
    break_points.set_index(pd.Series(zone,index=break_points.index), inplace=True,append=True)
    break_points = break_points.reorder_levels([1,0])
    
    return break_points


# In[7]:


def get_linear_function(p, q):
    m = (p[1] - q[1])/(p[0] - q[0])
    b = (p[0]*q[1] - p[1]*q[0])/(p[0] - q[0])
    return m, b


# In[8]:


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


# In[44]:


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

    lhs += gen_r.loc[:, committable.index].groupby(committable_grouper).sum()
    m.remove_constraints("Generator-com-p-upper")
    m.add_constraints( lhs <= rhs, name="Generator-com-p-upper")

    gen_status = m["Generator-status"]

    reserve_gen_status = gen_status.loc[:, committable.index]

    m.add_constraints(
        gen_r.loc[:, committable.index].groupby(committable_grouper).sum()
        <= (
            reserve_gen_status*
            (
                committable.groupby(committable_grouper).ramp_limit_up.sum()
                .multiply(delivery_time_reserves)
                .multiply(committable.groupby(committable_grouper).p_nom.sum())     
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

    if ORDC == True:

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
        m.add_objective(obj + res_shedding.sum()*reserve_scarcity_price + ordc_cost_term.sum(), overwrite=True) 

    else:

        m.add_constraints(
            gen_r.groupby(gen_grouper.to_xarray()).sum() 
            + storage_link_r.groupby(storage_grouper.to_xarray()).sum()
            + res_shedding
            >= reserve_requirements.reindex(n.buses.query("carrier == 'electricity'").index).dropna().to_xarray(),
            name="Bus-reserve_balance"
        )

        obj = m.objective

        m.add_objective(obj + res_shedding.sum()*reserve_scarcity_price, overwrite=True)


# In[10]:


def flow_based_market_coupling(n, snapshots)   :
    
    m = n.model
    
    ac_core = n.links.query("(carrier == 'AC') & (bus0 in @ptdf_core.columns) & (bus1 in @ptdf_core.columns) ")

    ac_bus0 = ac_core.bus0.copy()
    ac_bus0.name = "bus"

    ac_bus1 = ac_core.bus1.copy()
    ac_bus1.name = "bus"

    net_position_core = m.add_variables(coords=[snapshots, ptdf_core.columns], name="Bus-net_position")
    
    m.add_constraints(
        net_position_core  == (
            m["Link-p"].loc[:, ac_core.index].groupby(ac_bus0.to_xarray()).sum()
            - m["Link-p"].loc[:, ac_core.index].groupby(ac_bus1.to_xarray()).sum()
        ),
        name="Bus-net_position_zones"
    );

    net_position_alegro = m.add_variables(
        coords=[
            snapshots, 
            pd.Index(["Allegro"], name="Link")
        ],
        name="Net_position_alegro"
    )

    m.add_constraints(
        net_position_alegro == m["Link-p"].loc[:, "BE00-DE00-DC"] - m["Link-p"].loc[:, "DE00-BE00-DC"],
        name="Bus-alegro"
    );

    ahc_flows = n.links[
        (
            (n.links.bus0.isin(ptdf_core.columns) & ~n.links.bus1.isin(ptdf_core.columns))
            | (n.links.bus1.isin(ptdf_core.columns) & ~n.links.bus0.isin(ptdf_core.columns))
        ) & n.links.carrier.isin(['AC', 'DC'])

    ]

    ahc0 = pd.Series(ahc_flows.index.str[:-3], ahc_flows.index)
    ahc0 = ahc0[ahc0.isin(ptdf_ahc.columns)]
    ahc1 = pd.Series([i.split("-")[1] + "-" + i.split("-")[0] for i in ahc_flows.index], ahc_flows.index)
    ahc1 = ahc1[ahc1.isin(ptdf_ahc.columns)]
    ahc0.name = "Flow"
    ahc1.name = "Flow"

    external_flow = m.add_variables(
        coords = [snapshots, ptdf_ahc.columns],
        name="Link-external_flow"
    );

    m.add_constraints(
        external_flow == m["Link-p"].loc[:,ahc0.index].groupby(ahc0.to_xarray()).sum() - m["Link-p"].loc[:,ahc1.index].groupby(ahc1.to_xarray()).sum(),
        name="Link-external_flows"
    );

    m.add_constraints(
        (
            (ptdf_ahc.loc[snapshots].stack().to_xarray()*external_flow).sum("Flow")  
            + (ptdf_core.loc[snapshots].stack().to_xarray()*net_position_core).sum("bus")
            + (net_position_alegro*ptdf_dc.loc[snapshots].ALEGRO.to_xarray()).sum("Link")
        ) <= ram.loc[snapshots].to_xarray(),
        name="Link-RAM_constraint"
    )


# In[11]:


def build_flow_based_data(excel_file):

    ptdf_core = pd.read_excel(excel_file, "PTDF", index_col=[0,1,2], header=[0,1])["PTDF_SZ"]
    ptdf_ahc = pd.read_excel(excel_file, "PTDF", index_col=[0,1,2], header=[0,1])["PTDF*_AHC,SZ"]
    ptdf_dc = pd.read_excel(excel_file, "PTDF", index_col=[0,1,2], header=[0,1])["PTDF_EvFB"]
    
    domain_assign = pd.read_excel(excel_file, "FB Domain Assignment", index_col=[0]).iloc[:, 4:]
    domain_assign.columns = [int(i.split("_")[1]) for i in domain_assign.columns]
    domain_assign = domain_assign[climate_year]

    ptdf_core = ptdf_core.reset_index(2,drop=True).unstack(1).reindex(domain_assign.values)
    ptdf_core.index = n.snapshots
    ptdf_core = ptdf_core.stack(1)

    ptdf_ahc = ptdf_ahc.reset_index(2,drop=True).unstack(1).reindex(domain_assign.values)
    ptdf_ahc.index = n.snapshots
    ptdf_ahc = ptdf_ahc.stack(1)

    ptdf_dc = ptdf_dc.reset_index(2, drop=True).unstack(1).reindex(domain_assign.values)
    ptdf_dc.index = n.snapshots
    ptdf_dc = ptdf_dc.stack(1)


    ram_sheets = [i for i in excel_file.sheet_names if "RAM" in i]

    base_years = pd.Series([int(i.split("_")[1]) for i in ram_sheets], ram_sheets) - target_year

    ram_sheet = base_years[base_years <=0].idxmax()
    ram = pd.read_excel(excel_file, ram_sheet, index_col=[0], skiprows=3)
    ram.columns = ram.columns.astype(int)
    ram = ram.T.reindex(domain_assign.values)
    ram.index = n.snapshots
    ram = ram.stack()

    ptdf_core.index.names = ["snapshot", "CNEC"]
    ptdf_core.columns.name = "bus"
    ptdf_ahc.index.names = ["snapshot", "CNEC"]
    ptdf_ahc.columns.name = "Flow"
    ptdf_dc.index.names = ["snapshot", "CNEC"]
    ram.index.names = ["snapshot", "CNEC"]
    
    ptdf_core = ptdf_core.reindex(
        pd.MultiIndex.from_product(
            [ptdf_core.index.levels[0], ptdf_core.index.levels[1]]
        ),
        fill_value=0
    )
    
    ptdf_ahc = ptdf_ahc.reindex(
        pd.MultiIndex.from_product(
            [ptdf_ahc.index.levels[0], ptdf_ahc.index.levels[1]]
        ),
        fill_value=0
    )
    ptdf_dc = ptdf_dc.reindex(
        pd.MultiIndex.from_product(
            [ptdf_dc.index.levels[0], ptdf_dc.index.levels[1]]
        ), 
        fill_value=0
    )
    
    ptdf_ahc = ptdf_ahc[ptdf_ahc.sum()[ptdf_ahc.sum()!=0].index]
    ptdf_ahc = ptdf_ahc[[i for i in ptdf_ahc.columns if i.split("-")[0][:4] in n.buses.index and i.split("-")[1][:4] in n.buses.index]]
    ptdf_ahc = ptdf_ahc.groupby(ptdf_ahc.columns.str[:9],axis=1).sum()
    
    return ptdf_core, ptdf_ahc, ptdf_dc, ram


# In[12]:


def extra_functionality(n, snapshots):
    balancing_market(n, snapshots)
    if FBMC == True:
        flow_based_market_coupling(n, snapshots)

def market_revenue(n):
    
    spot_prices = n.buses_t.marginal_price.reindex(n.generators.bus.values, axis=1)
    spot_prices.columns = n.generators.index

    spot_payments = spot_prices.subtract(
        n.generators.marginal_cost
    ).multiply(n.generators_t.p).sum()

    reserve_prices = n.buses_t.mu_reserve_balance.reindex(n.generators.bus.values, axis=1)
    reserve_prices.columns = n.generators.index

    reserve_payments = reserve_prices.multiply(n.generators_t.r).sum()

    return (spot_payments + reserve_payments)


def apply_lm_and_cs(n):

    o = n.copy()

    o_snapshots = n.snapshots[n.generators_t.p.filter(like="load-shedding").sum(axis=1)>0]
    o.set_snapshots(o_snapshots)

    gen_techs = ['CCGT', 'OCGT', 'biomass', 'coal', 'lignite', 'nuclear', 'oil', 'other', 'onwind', 'solar', 'ROR', 'offwind', 'CSP', 'DSR'] 
    generators = o.generators.query("carrier in @gen_techs")

    o.generators_t.p_max_pu.loc[:, generators.index] = (
        o.generators_t.p
        .div(o.generators.p_nom.add(1e-6))
        [generators.index]
        .clip(lower=0)
    )

    discharger = n.links[n.links.type == "discharging"]
    charger = n.links[n.links.type == "charging"]

    discharge_gens = pd.DataFrame()
    discharge_gens["p_nom"] = n.links_t.p1[discharger.index].abs().max()
    discharge_gens["bus"] = o.links.loc[discharge_gens.index,"bus1"]

    discharge_p_max_pu = o.links_t.p1[discharger.index].abs().div(
        o.links_t.p1[discharger.index].abs().max().add(1e-6)
    ).clip(lower=0)

    o.madd(
        "Generator",
        discharge_gens.index,
        **discharge_gens,
        p_max_pu = discharge_p_max_pu
    )


    charging_load_t = o.links_t.p0[charger.index].groupby(o.links.bus0, axis=1).sum()
    charging_load_t.columns = charging_load_t.columns + " charging"
    charging_load = pd.DataFrame(index = charging_load_t.columns)
    charging_load["bus"] = charging_load.index.str[:4]
    charging_load = charging_load[charging_load_t.sum()>0]
    charging_load_t = charging_load_t[charging_load.index]

    o.madd("Load", charging_load.index, **charging_load, p_set = charging_load_t)
    o = o[o.buses.query("carrier == 'electricity'").index]

    o.generators_t.p.filter(like="load-shed").sum(axis=1)
    o.generators.p_min_pu = 0
    o.generators.committable = False

    gen_bus = o.generators.p_nom.multiply(o.generators_t.p_max_pu)[generators.index].groupby(generators.bus, axis=1).sum()
    discharge_bus = o.generators.p_nom.multiply(o.generators_t.p_max_pu)[discharge_gens.index].groupby(discharge_gens.bus, axis=1).sum()

    ls_active = gen_bus.add(discharge_bus, fill_value=0).subtract(
        o.loads_t.p_set.groupby(o.loads.bus, axis=1).sum()
    )<0


    load_shedding = o.generators.filter(like="load-shedding", axis=0)
    ls_active.columns = ls_active.columns.map(pd.Series(load_shedding.index, load_shedding.bus))

    o.generators_t.p_max_pu.loc[:,ls_active.columns] = ls_active.astype(float)

    o.optimize.optimize_with_rolling_horizon(solver_name=solver_name, extra_functionality=flow_based_market_coupling, **solver_options)

    n.generators_t.p.loc[o.snapshots, load_shedding.index] = o.generators_t.p[load_shedding.index]

    load_shedding_bus = o.generators_t.p[load_shedding.index].groupby(load_shedding.bus, axis=1).sum()
    correct_prices = n.buses_t.marginal_price.loc[o.snapshots][load_shedding_bus>0].stack()
    correct_prices.loc[correct_prices.index] = VOLL
    correct_prices = correct_prices.unstack(1)

    n.buses_t.marginal_price.loc[correct_prices.index, correct_prices.columns] = correct_prices


target_year = int(snakemake.params.ty)

solver_name = snakemake.config["solving"]["solver"]["name"]
options = snakemake.config["solving"]["solver"]["options"]
solver_options = snakemake.config["solving"]["solver_options"][options]

climate_year=int(snakemake.params.cy)

revenues = snakemake.output.revenues
lole = snakemake.output.lole
solved_network = snakemake.output.solved_network
n = pypsa.Network(snakemake.input.network)
years = snakemake.params.years

reserve_requirements = pd.read_excel("data/pemmdb.xlsx", "Reserve Requirements",index_col=[0,1])[["FCR (MW)", "FRR (MW)"]].sum(axis=1)#.loc[:, target_year]
reserve_requirements = reserve_requirements.unstack(1).reindex(years, axis=1).interpolate(axis=1)[target_year]
reserve_requirements.index.name = "Bus"
reserve_requirements = reserve_requirements[reserve_requirements>0]

flow_based = pd.ExcelFile(snakemake.input.fb_domains)
ptdf_core, ptdf_ahc, ptdf_dc, ram = build_flow_based_data(flow_based)


VOLL = float(snakemake.config["VOLL"])
reserve_scarcity_price = float(snakemake.config["balancing"]["price_cap"])
FBMC = snakemake.config["FBMC"]

reserve_participation_generators = snakemake.config["balancing"]["participation"]

ORDC = snakemake.config["balancing"]["ORDC"]["active"]

for su in n.storage_units.index:
    replace_su(n, su);
    
print(n.links.query("type=='charging'").efficiency)

storage_preopt_aggregation = snakemake.config["storage_preopt_aggregation"]

delivery_time_reserves = snakemake.config["delivery_time_reserves"] # in minutes
delivery_time_reserves = delivery_time_reserves/60

n.madd("Generator", 
       n.buses.query("carrier == 'electricity'").index + " load-shedding", 
       bus = n.buses.query("carrier == 'electricity'").index, 
       p_nom = 1e6, 
       marginal_cost = VOLL
      )
      
n.madd("Generator", 
       n.buses.query("carrier == 'electricity'").index + " gen-shedding", 
       bus = n.buses.query("carrier == 'electricity'").index, 
       p_nom = 1e6, 
       p_min_pu= -1,
       p_max_pu = 0,
       marginal_cost = -500
      )

ordc_parameters = pd.read_hdf(snakemake.input.ordc_parameters, "ordc_parameters")

m = average_every_nhours(n, storage_preopt_aggregation)

m.generators["p_min_pu"] = 0
m.stores.e_cyclic = True
m.generators.committable = False

m.optimize(
    solver_name=solver_name,
    **solver_options
)


# In[36]:


dispatchable = ['CCGT', 'OCGT', 'oil', 'biomass', 'other', 'lignite', 'nuclear', 'coal']
n.generators.loc[n.generators.query("carrier in @dispatchable").index, "committable"] = True
n.generators.shut_down_cost = n.generators.start_up_cost


# In[37]:


linear_ordc_approximation = get_line_parameters(ordc_parameters, n_ints = int(snakemake.config["balancing"]["ORDC"]["n_ints"]))

if FBMC == True:
    n.links.loc[n.links.bus0.isin(ptdf_core.columns) & n.links.bus1.isin(ptdf_core.columns), "p_nom"]

solve_rolling_horizon(m, n, "h", storage_preopt_aggregation)

apply_lm_and_cs(n)

dirname = os.path.dirname(solved_network)

if not os.path.isdir(dirname):
    os.makedirs(dirname)

n.export_to_netcdf(solved_network)

if not os.path.isdir(os.path.dirname(revenues)):
    os.makedirs(os.path.dirname(revenues))

market_revenue(n).to_hdf(revenues, "revenues")

if not os.path.isdir(os.path.dirname(lole)):
    os.makedirs(os.path.dirname(lole))

load_shed = n.generators.filter(like = "load-shedding", axis=0)
(n.generators_t.p[load_shed.index]>0).sum().to_frame().to_csv(lole)
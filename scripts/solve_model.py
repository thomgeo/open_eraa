#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pypsa
from scripts.__helpers import replace_su, average_every_nhours
import numpy as np
import os


def get_flow_based_groupers(PTDF_zone, PTDF_DC, PTDF_AHC):

    zones = PTDF_zone.columns

    AC = n.links.query("(carrier == 'AC') & (bus0 in @zones) & (bus1 in @zones)")
    AC_bus0 = AC.bus0.copy()
    AC_bus0.name = "name"
    AC_bus1 = AC.bus1.copy()
    AC_bus1.name = "name"
    
    DC = pd.Series(PTDF_DC.columns + "-DC", PTDF_DC.columns)
    DC = DC.to_frame("direct")
    DC["opposite"] = [i.split("-")[1] + "-" + i.split("-")[0] + "-DC" for i in DC.index]
    grouper_DC_direct = pd.Series(DC.index, DC.direct, name="name")
    grouper_DC_opposite = pd.Series(DC.index, DC.opposite, name="name")
    grouper_DC_direct.index.name="name"
    grouper_DC_opposite.index.name="name"
    
    
    AHC = PTDF_AHC.columns
    
    AHC_matching = pd.Series(index  = AHC)
    
    for i in AHC:
        matched_links = [j for j in n.links.index if "-".join(j.split("-")[:2]) ==i]
        if len(matched_links) >0:
            AHC_matching.loc[i] = matched_links[0]
    
    AHC_matching.dropna(inplace=True)
    AHC_matching = AHC_matching.to_frame("direct")
    
    for i in AHC_matching.index:
        direct_link = AHC_matching.loc[i, "direct"]
    
        bus0, bus1, carrier = direct_link.split("-")
        
        opposite = "{bus1}-{bus0}-{carrier}".format(bus0=bus0, bus1=bus1, carrier=carrier)
        if opposite in n.links.index:
            AHC_matching.loc[i, "opposite"] = opposite
    
    if "UK00-FR00_1" in AHC_matching.index:
    
        AHC_matching.loc["UK00-FR00_1", "opposite"] = "UK00-FR00_1-DC"
        AHC_matching.loc["UK00-FR00_2", "opposite"] = "UK00-FR00_2-DC"
    
    grouper_AHC_direct = pd.Series(AHC_matching.index, AHC_matching.direct.values)
    grouper_AHC_direct.index.names = ["name"]
    
    grouper_AHC_opposite = pd.Series(AHC_matching.index, AHC_matching.opposite.values)
    grouper_AHC_opposite.index.names = ["name"]

    return AC_bus0, AC_bus1, AC, grouper_DC_direct, grouper_DC_opposite, DC, grouper_AHC_direct, grouper_AHC_opposite, AHC_matching, 


def PTDF_timeseries(PTDF, domain_assignment, nordic=False):
    
    if nordic==False:
        PTDF = PTDF.unstack(1).reindex(domain_assignment.values)
    else:
        PTDF = PTDF.unstack(1).reindex([1 for i in range(len(n.snapshots))])

    
    PTDF.index = n.snapshots
    PTDF = PTDF.stack(1)

    return PTDF


def prepare_PTDFs(PTDF, domain_assignment, nordic=False):
    
    PTDF_zone = PTDF_timeseries(PTDF["PTDF_SZ"], domain_assignment, nordic=nordic)
    PTDF_DC = PTDF_timeseries( PTDF["PTDF_EvFB"], domain_assignment, nordic=nordic)
    PTDF_AHC = PTDF_timeseries(PTDF["PTDF*_AHC,SZ"], domain_assignment, nordic=nordic)
    
    PTDF_zone.columns.name = "name"
    PTDF_DC.columns.name = "name"
    PTDF_AHC.columns.name = "name"

    PTDF_zone.index.names=["snapshot", "CNEC"]
    PTDF_DC.index.names=["snapshot", "CNEC"]
    PTDF_AHC.index.names=["snapshot", "CNEC"]
    
    return PTDF_zone, PTDF_DC, PTDF_AHC


def prepare_for_uc():
    
    n.optimize.add_load_shedding(sign=1, marginal_cost=VOLL, buses=n.loads.bus.unique())

    """
    n.add(
        "Generator", 
        n.buses.query("carrier == 'electricity'").index + " gen-shedding", 
        bus = n.buses.query("carrier == 'electricity'").index, 
        p_nom = 1e6, 
        p_min_pu= -1,
        p_max_pu = 0,
        marginal_cost = -500,
        carrier="gen-shedding",
        ramp_limit_up = 1., 
        ramp_limit_down = 1.
    )
    """
    
    #n.generators.ramp_limit_up = n.generators.ramp_limit_up.clip(upper=1.).fillna(1.)
    #n.generators.ramp_limit_down = n.generators.ramp_limit_down.clip(upper=1.).fillna(1.)
    
    n.generators.shut_down_cost = n.generators.start_up_cost
    
    n.generators.up_time_before = n.generators.min_up_time
    n.generators.down_time_before = 0


# In[6]:


def flow_based_market_coupling(n, snapshots):
    
    
    core = PTDF_zone_core.columns
    nordic = PTDF_zone_nordic.columns
    
    core_groupers = get_flow_based_groupers(PTDF_zone_core, PTDF_DC_core, PTDF_AHC_core)
    
    AC_core_bus0 = core_groupers[0]
    AC_core_bus1 = core_groupers[1]
    AC_core = core_groupers[2]
    grouper_DC_core_direct = core_groupers[3]
    grouper_DC_core_opposite = core_groupers[4] 
    DC_core = core_groupers[5]
    grouper_AHC_core_direct  = core_groupers[6]
    grouper_AHC_core_opposite = core_groupers[7]
    AHC_core_matching = core_groupers[8]

    nordic_groupers = get_flow_based_groupers(PTDF_zone_nordic, PTDF_DC_nordic, PTDF_AHC_nordic)
    
    AC_nordic_bus0 = nordic_groupers[0]
    AC_nordic_bus1 = nordic_groupers[1]
    AC_nordic = nordic_groupers[2]
    grouper_DC_nordic_direct= nordic_groupers[3]
    grouper_DC_nordic_opposite = nordic_groupers[4]
    DC_nordic = nordic_groupers[5]
    grouper_AHC_nordic_direct = nordic_groupers[6]
    grouper_AHC_nordic_opposite = nordic_groupers[7]
    AHC_nordic_matching = nordic_groupers[8]
    
    m = n.model
    
    net_position_zone_core = m.add_variables(coords=[snapshots, core], name="Bus-zonal_NP_core")
    net_position_AHC_core = m.add_variables(coords=[snapshots, AHC_core_matching.index], name="Link-AHC_NP_core")
    net_position_DC_core = m.add_variables(coords = [snapshots, DC_core.index], name="Link-DC_NP_core")
    
    net_position_zone_nordic = m.add_variables(coords=[snapshots, nordic], name="Bus-zonal_NP_nordic")
    net_position_AHC_nordic = m.add_variables(coords=[snapshots, AHC_nordic_matching.index], name="Link-AHC_NP_nordic")
    net_position_DC_nordic = m.add_variables(coords = [snapshots, DC_nordic.index], name="Link-DC_NP_nordic")
    
    m.add_constraints(
        net_position_zone_core  == (
            m["Link-p"].loc[:, AC_core.index].groupby(AC_core_bus0.to_xarray()).sum()
            - m["Link-p"].loc[:, AC_core.index].groupby(AC_core_bus1.to_xarray()).sum()
        ),
        name="Bus-zonal_net_position_core"
    )
    
    m.add_constraints(
        net_position_DC_core == (
            m["Link-p"].loc[:, DC_core.direct.values].groupby(grouper_DC_core_direct).sum() 
            - m["Link-p"].loc[:, DC_core.opposite.values].groupby(grouper_DC_core_opposite).sum() 
        ),
        name="Link-DC_net_position_core"
    )
    
    
    m.add_constraints(
        net_position_AHC_core == (
            m["Link-p"].loc[:, AHC_core_matching.direct.values].groupby(grouper_AHC_core_direct).sum()
            - m["Link-p"].loc[:, AHC_core_matching.opposite.values].groupby(grouper_AHC_core_opposite).sum()
        ), name="Link-AHC_net_position_core"
    )
    
    m.add_constraints(
        net_position_zone_nordic  == (
            m["Link-p"].loc[:, AC_nordic.index].groupby(AC_nordic_bus0.to_xarray()).sum()
            - m["Link-p"].loc[:, AC_nordic.index].groupby(AC_nordic_bus1.to_xarray()).sum()
        ),
        name="Bus-zonal_net_position_nordic"
    )
    
    m.add_constraints(
        net_position_DC_nordic == (
            m["Link-p"].loc[:, DC_nordic.direct.values].groupby(grouper_DC_nordic_direct).sum() 
            - m["Link-p"].loc[:, DC_nordic.opposite.values].groupby(grouper_DC_nordic_opposite).sum() 
        ),
        name="Link-DC_net_position_nordic"
    )
    
    
    m.add_constraints(
        net_position_AHC_nordic == (
            m["Link-p"].loc[:, AHC_nordic_matching.direct.values].groupby(grouper_AHC_nordic_direct).sum()
            - m["Link-p"].loc[:, AHC_nordic_matching.opposite.values].groupby(grouper_AHC_nordic_opposite).sum()
        ), name="Link-AHC_net_position_nordic"
    )
    
    m.add_constraints(
        (net_position_DC_core*PTDF_DC_core.loc[snapshots].stack().to_xarray()).sum("name")
        + (net_position_zone_core*PTDF_zone_core.loc[snapshots].stack().to_xarray()).sum("name")
        + (net_position_AHC_core*PTDF_AHC_core.loc[snapshots].stack().to_xarray()).sum("name")
        <= RAM_core.loc[snapshots].stack().sort_index().to_xarray(),
        name="CNEC_capacity_constraint_core"
    )
    
    m.add_constraints(
        (net_position_DC_nordic*PTDF_DC_nordic.loc[snapshots].stack().to_xarray()).sum("name")
        + (net_position_zone_nordic*PTDF_zone_nordic.loc[snapshots].stack().to_xarray()).sum("name")
        + (net_position_AHC_nordic*PTDF_AHC_nordic.loc[snapshots].stack().to_xarray()).sum("name")
        <= RAM_nordic.loc[snapshots].stack().sort_index().to_xarray(),
        name="CNEC_capacity_constraint_nordic"
    )


def storage_targets(n, snapshots):

    m = n.model
    
    last_snapshot = snapshots[-1]
    soc = m["StorageUnit-state_of_charge"]
    
    soc_deviation = m.add_variables(lower=0, coords=[n.storage_units.index], name="StorageUnit-deviation_target")
    
    m.add_constraints(
        soc_deviation >= soc_target.loc[last_snapshot] - soc.loc[last_snapshot]  
    )

    obj = m.objective.expression 
    
    m.add_objective(
        obj + storage_target_deviation_penalty*(soc_deviation).sum("name"),
        overwrite=True
    )


def balancing_market(n, snapshots):
    
        
    m = n.model
    
    FRR_buses = reserve_requirements.xs("FRR", level=1).index  
    FCR_buses = reserve_requirements.xs("FCR", level=1).index

    FCR_gens = n.generators.loc[n.generators.bus.isin(FCR_buses)].query("carrier in @reserve_participation_generators")
    FRR_gens = n.generators.loc[n.generators.bus.isin(FRR_buses)].query("carrier in @reserve_participation_generators")

 
    gen_FRR = m.add_variables(lower=0, name="Generator-FRR", coords = [snapshots, FRR_gens.index])
    gen_FCR = m.add_variables(lower=0, name="Generator-FCR", coords = [snapshots, FCR_gens.index])
    non_committable = m.constraints["Generator-fix-p-upper"].coords["name"].to_index()
    non_committable_FRR = non_committable.intersection(FRR_gens.index)
    non_committable_FCR = non_committable.intersection(FCR_gens.index)

    if len(non_committable_FRR.union(non_committable_FCR)) >0:

        lhs = m.constraints["Generator-fix-p-upper"].lhs 
        rhs = m.constraints["Generator-fix-p-upper"].rhs
        m.remove_constraints("Generator-fix-p-upper")
        
        m.add_constraints(
            (lhs + gen_FRR.loc[:, non_committable_FRR] + gen_FCR.loc[:, non_committable_FCR])<= rhs,
            name="Generator-fix-p-upper"
        ) 

    lhs = m.constraints["Generator-com-p-upper"].lhs
    rhs = m.constraints["Generator-com-p-upper"].rhs
    m.remove_constraints("Generator-com-p-upper")
    

    committable_FRR = FRR_gens.query("committable == True").copy()
    committable_FCR = FCR_gens.query("committable == True").copy()

    m.add_constraints(
        lhs + gen_FRR.loc[:, committable_FRR.index] + gen_FCR.loc[:, committable_FCR.index] <= rhs,
        name="Generator-com-p-upper"
    )

    FRR_status = m["Generator-status"].loc[:, committable_FRR.index]
    FCR_status = m["Generator-status"].loc[:, committable_FCR.index]
    
    m.add_constraints(
        gen_FRR.loc[:, committable_FRR.index]<=FRR_status*committable_FRR.ramp_limit_up*delivery_time_FRR*committable_FRR.p_nom,
        name="Generator-FRR_ramping_limit"
    )
    
    m.add_constraints(
        gen_FCR.loc[:, committable_FCR.index]<=FCR_status*committable_FCR.ramp_limit_up*delivery_time_FCR*committable_FCR.p_nom,
            name="Generator-FCR_ramping_limit"
    )

    FRR_storages = n.storage_units.query("bus in @FRR_buses").copy()
    FCR_storages = n.storage_units.query("bus in @FCR_buses").copy()

    storage_FRR = m.add_variables(lower=0, coords=[snapshots, FRR_storages.index], name="StorageUnit-FRR")
    storage_FCR = m.add_variables(lower=0, coords=[snapshots, FCR_storages.index], name="StorageUnit-FCR")
    
    lhs = m.constraints["StorageUnit-fix-p_dispatch-upper"].lhs
    rhs = m.constraints["StorageUnit-fix-p_dispatch-upper"].rhs
    
    m.remove_constraints("StorageUnit-fix-p_dispatch-upper")
    
    m.add_constraints(lhs + storage_FRR + storage_FCR <= rhs, name="StorageUnit-fix-p_dispatch-upper")
    
    m.add_constraints(
        storage_FCR + storage_FRR <= m["StorageUnit-state_of_charge"],
        name="StorageUnit-reservoir_reserve_constraint"
    )

    FRR_buses.name = "name"
    FCR_buses.name = "name"

    FRR_shedding = m.add_variables(
        lower=0, 
        name="Bus-FRR_shedding", 
        coords=[snapshots, FRR_buses]
    )

    FCR_shedding = m.add_variables(
        lower=0, 
        name="Bus-FCR_shedding", 
        coords=[snapshots, FCR_buses]
    )

    FRR_gens.bus.name="name"
    FCR_gens.bus.name="name"

    FRR_storages.bus.name="name"
    FCR_storages.bus.name="name"

    m.add_constraints(
        gen_FRR.groupby(FRR_gens.bus.to_xarray()).sum()
        + storage_FRR.groupby(FRR_storages.bus).sum()
        + FRR_shedding == reserve_requirements.xs("FRR", level=1).to_xarray(),
        name="Bus-FRR_balance"
    )

    m.add_constraints(
        gen_FCR.groupby(FCR_gens.bus.to_xarray()).sum()
        + storage_FCR.groupby(FCR_storages.bus).sum()
        + FCR_shedding == reserve_requirements.xs("FCR", level=1).to_xarray(),
        name="Bus-FCR_balance"
    )

    obj = m.objective.expression
    
    m.add_objective(obj + (FRR_shedding.sum() + FCR_shedding.sum())*reserve_scarcity_price, overwrite=True)


def extra_functionality(n, snapshots):
    balancing_market(n, snapshots)
    storage_targets(n, snapshots)
    if FBMC == True:
        flow_based_market_coupling(n, snapshots)

def market_revenue(n):
    
    spot_prices_gens = n.buses_t.marginal_price.reindex(n.generators.bus.values, axis=1)
    spot_prices_gens.columns = n.generators.index
    
    spot_payments_gens = spot_prices_gens.subtract(
        n.generators.marginal_cost
    ).multiply(n.generators_t.p).sum()
    
    FRR_price_gens = n.buses_t.mu_FRR_balance.reindex(n.generators.bus, axis=1)
    FCR_price_gens = n.buses_t.mu_FCR_balance.reindex(n.generators.bus, axis=1)
    FRR_price_gens.columns = n.generators.index
    FCR_price_gens.columns = n.generators.index
    FRR_payments_gens = FRR_price_gens.multiply(n.generators_t.FRR).sum()
    FCR_payments_gens = FCR_price_gens.multiply(n.generators_t.FCR).sum()
    
    generator_revenues = pd.concat([spot_payments_gens, FRR_payments_gens, FCR_payments_gens], keys=["Spot", "FRR", "FCR"], axis=1)
    
    spot_prices_storages = n.buses_t.marginal_price.reindex(n.storage_units.bus, axis=1)
    spot_prices_storages.columns = n.storage_units.index
    spot_payments_storages = spot_prices_storages.multiply(n.storage_units_t.p).sum()
    FRR_price_storages = n.buses_t.mu_FRR_balance.reindex(n.storage_units.bus, axis=1)
    FCR_price_storages = n.buses_t.mu_FCR_balance.reindex(n.storage_units.bus, axis=1)
    FRR_price_storages.columns = n.storage_units.index
    FCR_price_storages.columns = n.storage_units.index
    FRR_payments_storages = FRR_price_storages.multiply(n.storage_units_t.FRR).sum()
    FCR_payments_storages = FCR_price_storages.multiply(n.storage_units_t.FCR).sum()
    
    storage_revenues = pd.concat([spot_payments_storages, FRR_payments_storages, FCR_payments_storages], keys=["Spot", "FRR", "FCR"], axis=1)
    
    return generator_revenues, storage_revenues

target_year = int(snakemake.params.ty)
climate_year=int(snakemake.params.cy)
weather_scenario = f"WS{climate_year:02}"

solver_name = snakemake.config["solving"]["solver"]["name"]
options = snakemake.config["solving"]["solver"]["options"]
solver_options = snakemake.config["solving"]["solver_options"][options]

revenues = snakemake.output.revenues
lole = snakemake.output.lole
solved_network = snakemake.output.solved_network

n = pypsa.Network(snakemake.input.network)

reserve_requirements = pd.read_csv(snakemake.input.reserve_requirements, index_col=[0,1,2, 3]).Value
reserve_requirements = reserve_requirements.loc["ERAA 2025 post-CfE", target_year]
reserve_requirements.index.names = ["name", "Category"]

RAM_core = pd.read_hdf(snakemake.input.core_domain, "RAM")
PTDF_core = pd.read_hdf(snakemake.input.core_domain, "PTDF")
domain_assignment_core = pd.read_hdf(snakemake.input.core_domain, "domain_assignment")

RAM_nordic = pd.read_hdf(snakemake.input.nordic_domain, "RAM")
PTDF_nordic = pd.read_hdf(snakemake.input.nordic_domain, "PTDF")
domain_assignment_nordic = pd.read_hdf(snakemake.input.nordic_domain, "domain_assignment")

domain_assignment_nordic = domain_assignment_nordic.loc[target_year][weather_scenario].sort_index()

if 29 in domain_assignment_nordic.loc[2].index.remove_unused_levels().levels[0]:
    domain_assignment_nordic = domain_assignment_nordic.drop(domain_assignment_nordic.loc[[2], 29, :].index)

domain_assignment_nordic.index = n.snapshots

domain_assignment_core = domain_assignment_core.loc[target_year][weather_scenario]
domain_assignment_core.sort_index(inplace=True)

if 29 in domain_assignment_core.loc[2].index.remove_unused_levels().levels[0]:
    domain_assignment_core = domain_assignment_core.drop(domain_assignment_core.loc[[2], 29, :].index)
    
domain_assignment_core.index = n.snapshots

RAM_core = RAM_core.T.reindex(domain_assignment_core.values)
RAM_core.index = n.snapshots
RAM_nordic = RAM_nordic.T.reindex(domain_assignment_nordic.values)
RAM_nordic.index = n.snapshots

PTDF_zone_core, PTDF_DC_core, PTDF_AHC_core = prepare_PTDFs(PTDF_core, domain_assignment_core)

PTDF_zone_nordic, PTDF_DC_nordic, PTDF_AHC_nordic = prepare_PTDFs(PTDF_nordic, domain_assignment_nordic, nordic=True)
PTDF_DC_nordic.drop("FI00-NON1", axis=1, inplace=True) # wrong assignment as DC

VOLL = float(snakemake.config["VOLL"])
reserve_scarcity_price = float(snakemake.config["balancing"]["price_cap"])
storage_target_deviation_penalty = float(snakemake.config["storage"]["target_deviation_penalty"])
FBMC = snakemake.config["FBMC"]
reserve_participation_generators = snakemake.config["balancing"]["participation"]
storage_preopt_aggregation = snakemake.config["storage"]["preopt_aggregation"]
rolling_horizon_window = int(snakemake.config["rolling_horizon"]["window"])
rolling_horizon_overlap = int(snakemake.config["rolling_horizon"]["overlap"])

delivery_time_FRR = float(snakemake.config["balancing"]["delivery_time"]["FRR"])/60 
delivery_time_FCR = float(snakemake.config["balancing"]["delivery_time"]["FCR"])/60 

m = average_every_nhours(n, "{}h".format(storage_preopt_aggregation))
m.generators["p_min_pu"] = 0
m.storage_units.cyclic_state_of_charge = True
m.generators.committable = False

FBMC_buses = PTDF_zone_core.columns.to_list() + PTDF_zone_nordic.columns.to_list() 

if FBMC == True:
    n.links.loc[n.links.query("bus1 in @FBMC_buses and bus0 in @FBMC_buses and (carrier == 'AC')").index, "p_nom"] *=5

m.optimize.add_load_shedding(sign=1, marginal_cost=VOLL)

status, condition = m.optimize(
    solver_name=solver_name,
    **solver_options
)

if "infeasible" in condition:
    raise RuntimeError("Solving status 'infeasible'")

soc_target = m.storage_units_t.state_of_charge.copy()
soc_target = soc_target.reindex(n.snapshots).shift(storage_preopt_aggregation-1).dropna()
n.storage_units.state_of_charge_initial = m.storage_units_t.state_of_charge.iloc[-1]
prepare_for_uc()

n.optimize.optimize_with_rolling_horizon(
    solver_name=solver_name,
    **solver_options,
    horizon=rolling_horizon_window,
    overlap=rolling_horizon_overlap,
    linearized_unit_commitment=True,
    extra_functionality=extra_functionality,
    assign_all_duals=True,
)

loss_of_load = (n.generators_t.p.filter(like="load shedding", axis=1)>1).sum().groupby(n.generators.filter(like="load shedding", axis=0).bus).sum()
dirname = os.path.dirname(solved_network)

if not os.path.isdir(dirname):
    os.makedirs(dirname)

n.export_to_netcdf(solved_network)

generator_revenues, storage_revenues = market_revenue(n)

if not os.path.exists(os.path.dirname(revenues)):
    os.makedirs(os.path.dirname(revenues))

if not os.path.exists(os.path.dirname(lole)):
    os.makedirs(os.path.dirname(lole))

generator_revenues.to_hdf(revenues, "generators")
storage_revenues.to_hdf(revenues, "storages")

loss_of_load.to_hdf(lole, "LOLE")


#!/usr/bin/env python
# coding: utf-8

# In[43]:


import pandas as pd
import powerplantmatching as pm
import pypsa
import pycountry
import os 


# In[2]:


def prepare_capacity_table(capacity):
    
    capacity_matching = pd.Series(
        ["nuclear", "lignite", "coal", "gas", "oil", "other", "ROR", "ROR", "reservoir", "reservoir", "PHS",
         "reservoir (pumping)", "PHS (pumping)", "onwind", "offwind", "CSP", "solar", "biomass", "biomass", "battery"],
        capacity.index
    )

    capacity = capacity.groupby(capacity_matching, ).sum()
    
    pm_plants = pm.powerplants()
    gas_share = pd.DataFrame()
    gas_share["CCGT"] = pm_plants.groupby(["Fueltype", "Technology", "Country"]).sum(numeric_only=True).loc["Natural Gas", "CCGT",:].Capacity.div(
        pm_plants.groupby(["Fueltype", "Technology", "Country"]).sum(numeric_only=True).loc["Natural Gas", :].Capacity.groupby(level=1).sum(),
        fill_value=0
    )

    gas_share.index = [pycountry.countries.lookup(i).alpha_2 for i in gas_share.index]

    gas_share["OCGT"] = 1- gas_share["CCGT"]

    gas_share = gas_share.reindex(capacity.columns.str[:2]).fillna(gas_share.mean())
    gas_share.index = capacity.columns

    capacity = pd.concat([
        capacity.drop("gas"),
        gas_share.multiply(capacity.loc["gas"],axis=0).T
    ])
    
    return capacity


# In[3]:

def build_capacity_table():

    capacity_years = pd.DataFrame()

    for sheet in plant_sheets.index:

        capacity_raw = pd.read_excel(excel_file, sheet, header=1, index_col=1).dropna(how="all", axis=1).iloc[:20]

        capacity_years[sheet[2:]] = prepare_capacity_table(capacity_raw).stack().astype(float)

    capacity_years.columns = capacity_years.columns.astype(int)
    capacity_years = capacity_years.reindex(years,axis=1)
    capacity_years.loc[distributed_resources, :] = capacity_years.loc[distributed_resources, :].interpolate(axis=1)

    return capacity_years.fillna(method="ffill", axis=1)[year].unstack(1)


def build_thermal_properties(properties):
    
    for col in properties.columns:
        try: 
            properties[col] = properties[col].astype(float)
        except:
            None

    properties.set_index(properties.index.remove_unused_levels(), inplace=True)

    properties_matching = pd.Series(
        ["OCGT", "coal", "oil", "oil", "lignite", "nuclear", "oil"],
        properties.index.levels[0]
    )

    properties_matching = properties_matching.reindex([i[0] for i in properties.index])

    properties_matching.index = properties.index

    properties_matching.loc[[i for i in properties_matching.index if "CCGT" in i[1]]] = "CCGT"

    properties_transformed = properties.iloc[:, 1:].groupby(
        properties_matching.reindex(properties.index).values
    ).mean()

    missing_properties = properties.reindex([("Gas", "conventional old 2"), ("Gas", "conventional old 2")]).copy()
    missing_properties.index = ["biomass", "other"]
    missing_properties.loc["biomass","CO2 emission factor"] = 0

    properties_transformed = pd.concat([properties_transformed, missing_properties[properties_transformed.columns]])

    properties_transformed.set_index(pd.Series("existing", properties_transformed.index),append=True, inplace=True)

    properties_new = properties.loc["Gas", ["CCGT new", "OCGT new"], :].set_index(pd.MultiIndex.from_product([["CCGT", "OCGT"], ["new"]]))

    properties_transformed = pd.concat([properties_transformed, properties_new[properties_transformed.columns]])
    properties_transformed.index.names=["carrier", "invest_status"]

    return properties_transformed


def build_inflows(inflows):
    base_year_hydro = pd.Series(
        inflows.index.levels[4].astype(int), 
        inflows.index.levels[4]).subtract(year).abs().idxmin()

    inflows = inflows.loc[:, str(climate_year), :, :, base_year_hydro]

    inflow_grouper = pd.Series(["ROR", "PHS", "reservoir", "reservoir", "ROR"], inflows.index.levels[1])

    inflows = inflows.unstack(1).groupby(inflow_grouper, axis=1).sum().multiply(1e3) # conversion to MWh
    
    inflows = inflows.unstack(0).stack(0)
    inflows.index = [" ".join(i) for i in inflows.index]
    
    inflows = inflows.T
    
    inflows.index = snapshots
    
    return inflows


# In[6]:


def add_existing_storage():  
 
    storage = capacity.loc[["PHS", "reservoir", "battery"]].unstack().reset_index().copy()
    storage.columns = ["bus", "carrier", "p_nom"]
    storage.index = storage.bus + " " + storage.carrier

    p_min_pu = capacity.loc[["PHS (pumping)", "reservoir (pumping)",]].div(capacity.loc[["PHS", "reservoir"]].add(1e-5).values)
    p_min_pu.index = ["PHS", "reservoir"]

    storage["p_min_pu"] = (
        p_min_pu.unstack()
        .reindex(pd.MultiIndex.from_arrays([storage.bus, storage.carrier]))
    ).values

    storage.loc[storage.carrier=="battery", "p_min_pu"] = -1

    storage = storage[storage.p_nom >0]

    storage_inflows = inflows[storage.loc[storage.carrier=="reservoir"].index]

    storage_capacity_raw = pd.read_excel(excel_file, sheet, index_col=1, header=28).dropna(how="all", axis=1)

    storage_capacity= storage_capacity_raw.groupby(
        ["ROR", "reservoir", "reservoir", "PHS", "battery"]
    ).sum()

    storage["max_hours"] = storage_capacity.unstack().reindex(
        pd.MultiIndex.from_arrays(
            [storage.bus, storage.carrier]  
        )
    ).div(
        storage.p_nom.values
    ).values

    map_storage = pd.Series( ["battery inverter", "PHS", "hydro"], ["battery", "PHS", "reservoir"])

    storage["efficiency"] = technology_data.loc[:, "efficiency", :].reindex(storage.carrier.map(map_storage)).value.values

    n.madd(
        "StorageUnit",
        storage.index,
        **storage,
        inflow=storage_inflows,
        invest_status = "existing"
    )


def prepare_commodity_prices(commodity_prices):

    commodity_prices = commodity_prices_raw.iloc[:, 2:].groupby(commodity_prices.index).mean()
    commodity_prices = commodity_prices.groupby(["co2", "CCGT", "coal", "oil", "hydrogen", "oil", "lignite", "nuclear", "oil"]).mean()
       
    converter = pd.DataFrame(3.6 , commodity_prices.index, commodity_prices.columns)
    converter["CO2"] = 3.6/1000
    converter.loc["co2"] = 1
    
    commodity_prices = converter*commodity_prices

    to_add = commodity_prices.reindex(["CCGT", "CCGT", "CCGT"])
    to_add.index = ["OCGT", "biomass", "other"]
    to_add.loc[["other"]] = 0
    to_add.loc[["biomass"]] = biomass_price
    to_add.loc[["biomass"], "CO2"] = 0
    commodity_prices = pd.concat([commodity_prices, to_add])
    
    return commodity_prices

def add_existing_dispatchables(legacy):
    
    legacy["efficiency"] = properties["Standard efficiency in NCV terms"].loc[:, "existing"].reindex(legacy.carrier).values

    legacy["start_up_cost"] = properties["Start-up fix cost (e.g. wear) warm start"].loc[:, "existing"].reindex(legacy.carrier, level=1).fillna(0).values
    legacy["ramp_limit_up"] = properties["Ramp up rate % of max output power / min"].loc[:, "existing"].reindex(legacy.carrier).values*60
    legacy["ramp_limit_down"] = properties["Ramp down rate % of max output power / min"].loc[:, "existing"].reindex(legacy.carrier).values*60
    legacy["p_min_pu"] = properties["Minimum stable generation (% of max power)"].loc[:, "existing"].reindex(legacy.carrier).values
    legacy["min_up_time"] = properties["Min Time on"].loc[:, "existing"].reindex(legacy.carrier).values
    legacy["min_down_time"] = properties["Min Time off"].loc[:, "existing"].reindex(legacy.carrier).values

    legacy["marginal_cost"] = commodity_prices.reindex(legacy.carrier)[base_year].div(
        legacy.efficiency.values,
    ).add(
        (
            commodity_prices.reindex(
                legacy.carrier)["CO2"]
            .multiply(
                commodity_prices.loc["co2", base_year]
            )
        )
    ).add(
        properties.loc[:, "existing",:].reindex(legacy.carrier)["Variable O&M cost"]
    ).values

    legacy = legacy.loc[(legacy.entry <= year) & (legacy.exit > year)]
    
    n.madd(
        "Generator",
        legacy.index,
        **legacy,
        invest_status = "existing"
    )


def add_renewables():
    
    res = capacity.loc[["onwind","offwind", "solar", "CSP", "ROR"]].unstack().copy()
    res = res.reset_index()
    res.columns = ["bus", "carrier", "p_nom"]
    res.index = res.bus + " " + res.carrier

    res = res[res.p_nom>0]

    vre = ["onwind", "offwind", "solar", "CSP"]

    base_year_res = int(plant_sheets[plant_sheets >= year].subtract(year).idxmin()[2:])
    
    p_max_pu = pd.concat(
        [pd.read_hdf("resources/res_profile.h5", tech).loc[:, str(climate_year), str(base_year_res), :] for tech in vre],
        axis=1
    )


    p_max_pu.columns = vre
    p_max_pu = p_max_pu.unstack(0).stack(0)
    p_max_pu.index = [" ".join(i) for i in p_max_pu.index]
    p_max_pu.columns = snapshots
    p_max_pu = p_max_pu.T

    p_max_pu = pd.concat(
        [p_max_pu, inflows[res.filter(like="ROR", axis=0).index].div(res.filter(like="ROR", axis=0).p_nom)],
        axis=1
    )

    n.madd(
        "Generator",
        res.index,
        **res,
        p_max_pu = p_max_pu[res.index],
        invest_status="policy"
    )

def group_luxembourg(demand, links):
    
    demand_grouper = pd.Series(demand.columns, demand.columns)
    demand_grouper.loc[demand_grouper.index.str[:2] == "LU"] = "LUG1"
    demand = demand.groupby(demand_grouper, axis=1).sum()

    links.loc[links.bus0.str[:2] == "LU", "bus0"] = "LUG1"
    links.loc[links.bus1.str[:2] == "LU", "bus1"] = "LUG1"
    
    return demand, links

def add_dsr():
    
    p_nom_dsr.columns = ["band " + i.split(" ")[2] for i in p_nom_dsr.columns]

    marginal_cost_dsr.columns = ["band " + i.split(" ")[8] for i in marginal_cost_dsr.columns]

    dsr = pd.concat([p_nom_dsr.stack(), marginal_cost_dsr.stack()],axis=1)

    dsr = dsr.loc[:, base_year, :]

    dsr = dsr.set_index((i[0] + " dsr " + i[1] for i in dsr.index), append=True).reset_index([0,1]).drop("level_1",axis=1)

    dsr.columns = ["bus", "p_nom", "marginal_cost"]
    
    dsr["carrier"] = "DSR"

    n.madd(
        "Generator",
        dsr.index,
        **dsr,
        invest_status = "existing"
    )
    
    
def add_dispatchable_investment_options():

    zones_for_investment = capacity.sum()[capacity.loc[dispatchable].sum() >0 ].index # no investment in offshores zones etc.


    new_dispatchables = pd.DataFrame(
        0.01, 
        index= pd.MultiIndex.from_product(
            [zones_for_investment,technologies_for_investment, range(years[0], year +1)], 
            names=["bus", "carrier", "invest_year"]), 
        columns=["p_nom"]
    )


    new_dispatchables.reset_index(inplace=True)

    new_dispatchables.index = new_dispatchables.bus + " " + new_dispatchables.carrier + " new " + new_dispatchables.invest_year.astype(str) 

    new_dispatchables["efficiency"] = properties.loc[:, "new", :]["Standard efficiency in NCV terms"].reindex(new_dispatchables.carrier).values

    new_dispatchables["start_up_cost"] = properties["Start-up fix cost (e.g. wear) warm start"].loc[:, "new"].reindex(new_dispatchables.carrier, level=1).fillna(0).values

    new_dispatchables["ramp_limit_up"] = properties["Ramp up rate % of max output power / min"].loc[:, "new"].reindex(new_dispatchables.carrier).values*60
    new_dispatchables["ramp_limit_down"] = properties["Ramp down rate % of max output power / min"].loc[:, "new"].reindex(new_dispatchables.carrier).values*60
    new_dispatchables["p_min_pu"] = properties["Minimum stable generation (% of max power)"].loc[:, "new"].reindex(new_dispatchables.carrier).values
    new_dispatchables["min_up_time"] = properties["Min Time on"].loc[:, "new"].reindex(new_dispatchables.carrier).values
    new_dispatchables["min_down_time"] = properties["Min Time off"].loc[:, "new"].reindex(new_dispatchables.carrier).values

    new_dispatchables["marginal_cost"] = commodity_prices.reindex(new_dispatchables.carrier)[base_year].div(
            new_dispatchables.efficiency.values,
        ).add(
            (
                commodity_prices.reindex(
                    new_dispatchables.carrier)["CO2"]
                .multiply(
                    commodity_prices.loc["co2", base_year]
                )
            )
        ).add(
            properties.loc[:, "existing",:].reindex(new_dispatchables.carrier)["Variable O&M cost"]
        ).values

    n.madd(
        "Generator",
        new_dispatchables.index,
        **new_dispatchables,
        invest_status = "new"
    )

def set_investment_bounds():
    
    n.generators.loc[n.generators.invest_status == "existing", "p_nom_max"] = n.generators.loc[n.generators.invest_status == "existing", "p_nom"].clip(0.01)
    n.generators.loc[n.generators.invest_status == "existing", "p_nom"] = n.generators.loc[n.generators.invest_status == "existing", "p_nom"].clip(0.01) 
    n.generators.p_nom_min = 0.01
    n.generators.loc[n.generators.invest_status == "policy", "p_nom_max"] = n.generators.loc[n.generators.invest_status == "policy", "p_nom"]
    n.generators.loc[n.generators.invest_status == "policy", "p_nom_min"] = n.generators.loc[n.generators.invest_status == "policy", "p_nom"]
    
commodity_prices_raw = pd.read_excel(snakemake.input.commodity_prices, index_col = 0)

year = int(snakemake.params.ty)
climate_year = int(snakemake.params.cy)

years = snakemake.params.years

save_path = snakemake.output.network

biomass_price = snakemake.config["biomass_price"]

dispatchable = ["CCGT", "OCGT", 'biomass', 'coal', 'lignite', 'nuclear','oil',  'other', ]
distributed_resources = ["onwind", "offwind", "CSP","solar", "battery"]

technologies_for_investment = snakemake.config["investment"]["technologies_new_investment"]

inflows_raw = pd.read_hdf(snakemake.input.inflow, "inflow")
excel_file = pd.ExcelFile(snakemake.input.pemmdb)

snapshots = pd.date_range(start="2010-01-01", freq="h", periods=8760)#.strftime('%m-%d %H:%M:%S')
#snapshots = snapshots = pd.period_range(start="01-01 00:00:00", freq="h", periods=8760)

plant_sheets = [i for i in excel_file.sheet_names if "TY" in i]
plant_sheets = pd.Series([int(i[2:]) for i in plant_sheets], plant_sheets )
sheet = plant_sheets[plant_sheets <= year].subtract(year).idxmax()

base_year = int(sheet[2:])

technology_data = pd.read_csv(snakemake.input.technology_data, index_col=[0,1])

demand = pd.read_hdf(snakemake.input.demand)
demand = demand.loc[:, climate_year, :, str(base_year), :].unstack(1)
demand.index = snapshots
demand.drop("TR00",axis=1, inplace=True)

links = pd.read_hdf(snakemake.input.ntc, "p_nom")
links_p_max_pu = pd.read_hdf(snakemake.input.ntc, "p_max_pu")
links = links[base_year].unstack(1)
links_p_max_pu = links_p_max_pu[base_year].unstack(2)
links_p_max_pu.index = snapshots
links.dropna(inplace=True)
links_p_max_pu.dropna(axis=1, inplace=True)
links["carrier"] = [i[-2:] for i in links.index]

demand, links = group_luxembourg(demand, links)

commodity_prices = prepare_commodity_prices(commodity_prices_raw)

properties_raw = pd.read_excel(excel_file, "Thermal Properties", index_col=[2,3], header=3).dropna(how="all").iloc[1:, 2:].dropna(how="all", axis=1).iloc[:24, :12]
properties_raw2 = pd.read_excel(excel_file, "Thermal Properties", index_col=[2,3], skiprows=35, header=[0,3]).iloc[:, 2:].dropna(how="all", axis=1).dropna(how="all")
properties_raw2.columns = [" ".join(i) for i in properties_raw2.columns]
properties_raw = pd.concat([properties_raw, properties_raw2],axis=1)

properties = build_thermal_properties(properties_raw)

legacy = pd.read_hdf(snakemake.input.legacy_capacities)


n = pypsa.Network()
n.set_snapshots(snapshots)

capacity = build_capacity_table()

p_nom_dsr = pd.read_excel(excel_file, "Explicit DSR", index_col = [0,1]).iloc[:, :8]
marginal_cost_dsr = pd.read_excel(excel_file, "Explicit DSR", index_col = [0,1]).iloc[:, 8:16]

buses = (
    capacity.sum()[capacity.abs().sum() >0].index
    .union(demand.columns)
    .union(links.bus0.unique())
    .union(links.bus1.unique())
)

n.madd(
    "Bus", 
    buses, 
    carrier = "electricity", 
    country = buses.str[:2]
)

n.madd(
    "Load", 
    demand.columns,
    bus=demand.columns,
    p_set = demand
)

inflows = build_inflows(inflows_raw)
add_existing_storage()

add_existing_dispatchables(legacy)

add_renewables()

add_dsr()

add_dispatchable_investment_options()

set_investment_bounds()

n.madd(
    "Link",
    links.index,
    bus0 = links.bus0,
    bus1 = links.bus1,
    p_nom = links.p_nom,
    p_max_pu = links_p_max_pu,
    carrier = links.carrier,
)

n.generators.loc[(n.generators.bus == "CY00") & (n.generators.carrier == "CCGT"), "p_min_pu"] = 0 # remove minimum load of CCGT in Cyprus as this can exceed actual load.

dirname = os.path.dirname(save_path)  

if not os.path.exists(dirname):
    os.makedirs(dirname)

n.export_to_netcdf(save_path)



#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import powerplantmatching as pm
import pycountry
import os 
import zipfile
import numpy as np


# In[3]:


def build_capacity_table(yf):

    capacity_years = pd.DataFrame()
    capacity_disp = pd.DataFrame()
 

    for year, df in yf.items():
   
        capacity_years[year], capacity_disp[year] = [
        df.stack().astype(float) for df in prepare_capacity_table(yf[year])
]

    capacity_disp.columns = capacity_disp.columns.astype(int)
    capacity_disp = capacity_disp.reindex(years,axis=1)
    capacity_years.columns = capacity_years.columns.astype(int)
    capacity_years = capacity_years.reindex(years,axis=1)
    #capacity_years.loc[distributed_resources, :] = capacity_years.loc[distributed_resources, :].interpolate(axis=1)

    return [capacity_disp.fillna(method="ffill", axis=1), capacity_years.fillna(method="ffill", axis=1)]


def build_legacy_caps():

    existing = pd.DataFrame()

    for i in capacity.loc[dispatchable].index:

        capacity_evolution = capacity.loc[i]

        capacity_changes = capacity_evolution.diff()

        commissioning = capacity_changes[capacity_changes>0]
        decommissioning = capacity_changes[capacity_changes<0]

        commissioning.index.names=["entry"]
        commissioning.name = "p_nom"

        legacy = pd.Series(name="p_nom")
        for iyear in decommissioning.index:
            legacy[iyear] = - decommissioning.loc[iyear]

        legacy[years[-1] +1] = capacity_evolution.loc[years[0]]- legacy.sum()

        legacy.index.name = "exit"

        legacy = legacy.to_frame().reset_index()
        legacy["entry"] = years[0]
        legacy.index = i[1] + " " + i[0] + " exit " + legacy.exit.astype(str)
        legacy["carrier"] = i[0]
        legacy["bus"] = i[1]

        commissioning = commissioning.to_frame().reset_index()
        commissioning.index = i[1] + " " + i[0] + " entry " + commissioning.entry.astype(str)
        commissioning["exit"] = np.inf #years[-1] + 1
        commissioning["carrier"] = i[0]
        commissioning["bus"] = i[1]

        existing = pd.concat([existing, legacy, commissioning])
        existing = existing[["bus", "p_nom", "carrier", "entry", "exit"]]
        existing = existing[existing.p_nom >0]
    

    legacy_properties = technological_parameters.drop(["new", "CCS"], level=1).groupby(level=0).mean()
    
    legacy_properties.min_down_time = legacy_properties.min_down_time.round().astype(int)
    legacy_properties.min_up_time = legacy_properties.min_up_time.round().astype(int)

    existing["efficiency"] = legacy_properties["efficiency"].reindex(existing.carrier).values
    existing["start_up_cost"] = legacy_properties.start_up_fix_cost.reindex(existing.carrier).values
    existing["ramp_limit_up"] = legacy_properties.ramp_limit_up.reindex(existing.carrier).values

    existing["ramp_limit_down"] = legacy_properties.ramp_limit_down.reindex(existing.carrier).values

    existing["p_min_pu"] = legacy_properties.p_min_pu.reindex(existing.carrier).values

    existing["min_up_time"] = legacy_properties.min_up_time.reindex(existing.carrier).values
    existing["min_down_time"] = legacy_properties.min_down_time.reindex(existing.carrier).values

    existing["var_om"] = legacy_properties.var_OM.reindex(existing.carrier).fillna(0).values

    existing["p_nom_max"] = existing["p_nom"]
    existing["p_nom_min"] = 0.01
    existing["invest_status"] = "existing"

    return existing
    
def build_new_investments():

    zones_for_investment = (
        capacity.groupby(level=1).sum().sum(axis=1)
        [capacity.loc[dispatchable].groupby(level=1).sum().sum(axis=1)>0]
        .index
    ) # no investment in offshores zones etc.


    new = pd.DataFrame(
        0.01, 
        index= pd.MultiIndex.from_product(
            [zones_for_investment,technologies_for_investment, years], 
            names=["bus", "carrier", "entry"]), 
        columns=["p_nom"]
    )

    new.reset_index(inplace=True)

    new.index = new.bus + " " + new.carrier + " new " + new.entry.astype(str)

    new["efficiency"] = technological_parameters.efficiency.loc[:, "new", :].reindex(new.carrier).values
    new["start_up_cost"] = technological_parameters.start_up_fix_cost.loc[:, "new"].reindex(new.carrier).values
    new["ramp_limit_up"] = technological_parameters.ramp_limit_up.loc[:, "new"].reindex(new.carrier).values    
    new["ramp_limit_down"] = technological_parameters.ramp_limit_down.loc[:, "new"].reindex(new.carrier).values    
    new["p_min_pu"] = technological_parameters.p_min_pu.loc[:, "new"].reindex(new.carrier).values    
    new["min_up_time"] = technological_parameters.min_up_time.loc[:, "new"].reindex(new.carrier).values    
    new["min_down_time"] = technological_parameters.min_down_time.loc[:, "new"].reindex(new.carrier).values    
    new["var_om"] = technological_parameters.var_OM.loc[:, "new"].reindex(new.carrier).values
    new["exit"] = np.inf
    new["p_nom_max"] = np.inf
    new["p_nom_min"] = 0.01
    new["invest_status"] = "new"
    
    return new

"""
def build_thermal_properties(properties):
    
    for col in properties.columns:
        try: 
            properties[col] = properties[col].astype(float)
        except:
            None

    properties.set_index(properties.index.remove_unused_levels(), inplace=True)

    properties_matching = pd.Series(
        ["OCGT", "coal", "oil", "hydrogen", "oil", "lignite", "nuclear", "oil"],
        properties.index.levels[0]
    )

    properties_matching = properties_matching.reindex([i[0] for i in properties.index])

    properties_matching.index = properties.index

    mask = (
    properties_matching.index.get_level_values(1).str.contains("CCGT") &
    (properties_matching.index.get_level_values(0) != "Hydrogen"))

    # Update only the entries matching the mask
    properties_matching.loc[mask] = "CCGT"
    
    coal_index = properties_matching[properties_matching == "hydrogen"].index
    
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

"""
    
def prepare_capacity_table(capacity):
      
    pm_plants = pm.powerplants()
    gas_share = pd.DataFrame()
    gas_share["CCGT"] = pm_plants.groupby(["Fueltype", "Technology", "Country"]).sum(numeric_only=True).loc["Natural Gas", "CCGT",:].Capacity.div(
        pm_plants.groupby(["Fueltype", "Technology", "Country"]).sum(numeric_only=True).loc["Natural Gas", :].Capacity.groupby(level=1).sum(),
        fill_value=0
    )


    index = ["biomass", "coal", 'hydrogen', 'lignite', 'nuclear','oil','CCGT','OCGT']

    gas_share.index = [pycountry.countries.lookup(i).alpha_2 for i in gas_share.index]

    gas_share["OCGT"] = 1- gas_share["CCGT"]

    gas_share = gas_share.reindex(capacity.columns.str[:2]).fillna(gas_share.mean())
    gas_share.index = capacity.columns

    capacity = pd.concat([
        capacity.drop("gas"),
        gas_share.multiply(capacity.loc["gas"],axis=0).T
    ])
    
    capacity_disp = capacity.reindex(index=index).dropna(how='all').fillna(0)
    
    return [capacity,capacity_disp]


technological_parameters = pd.read_hdf(snakemake.input.technology_parameters)

#Different methodology to get data 
with zipfile.ZipFile(snakemake.input.pemmdb) as zip_f:
    zip_f.extractall(snakemake.params.data_folder)
        
#with zipfile.ZipFile(snakemake.input.thermal) as zip_f:
#    zip_f.extractall(snakemake.params.data_folder)
    
with zipfile.ZipFile(snakemake.input.dsr) as zip_f:
    zip_f.extractall(snakemake.params.data_folder)
    
#file2 = snakemake.params.data_folder + "Common data/Common Data.xlsx"

all_columns = pd.read_csv(
    snakemake.params.data_folder + "Other data/Explicit DSR detailed.csv",
    nrows=0
).columns.tolist()

columns_to_use = all_columns[:-1]

properties_raw2 = pd.read_csv(
    snakemake.params.data_folder + "Other data/Explicit DSR detailed.csv",
    usecols=columns_to_use,
    header=0
)

years = properties_raw2['TARGET_YEAR'].unique()

DSR = {
  year: group.rename(columns={
      "MARKET_NODE": "bus",
      "DA activation price (EUR/MWh)": "marginal_cost",
      "Max hours dispatched per day [h]": "hours",
      "MAX CAP (MW)": "p_nom"
  })
  for year, group in properties_raw2.groupby("TARGET_YEAR")
}

DSR = pd.concat(DSR, names=["TARGET_YEAR"]).fillna(0)

DSR.to_hdf(snakemake.output.dsr, "dsr")

#excel_file2 = pd.ExcelFile(file2)

generation = pd.read_csv(snakemake.params.data_folder + "Dashboard_raw_data/GenerationCapacities.csv", header=0)

filtered_rows = generation[
    generation.apply(
        lambda row: row.astype(str).str.contains("ERAA 2025 pre-CfE").any(), axis=1
    )
]

tech_names = sorted(filtered_rows["Technology"].unique())

tech_mapping = pd.Series(
    [
        "battery", "battery", "biomass", "PHS", "DSR", "DSR", "electrolyser", "geothermal",
        "DSR", "coal", "oil", "hydrogen", "oil", "lignite", "marine", "gas", "nuclear",
        "PHS Open", "pondage", "P2H", "hydro", "ROR", "oil", "biomass", "solar-ind", "solar-rsd",
        "solar-fix", "solar-track", "CSP-stor", "CSP", "biomass", "offwind", "offwind", "onwind"
    ],
    index=tech_names
).to_dict()

filtered_rows["Technology"] = filtered_rows["Technology"].map(tech_mapping)

filtered_rows .loc[
    (filtered_rows["Technology"] == "ROR") & (filtered_rows["Market_Node"] == "UKNI"),
    "Value"
] = 0

all_caps = {
    year: (
        group.groupby(["Technology", "Market_Node"])["Value"]
             .sum()
    )
    for year, group in filtered_rows.groupby("Target year")
}

all_caps = pd.concat(all_caps, names=["Target year"]).fillna(0)

all_caps.to_hdf(snakemake.output.all_capacities, "capacities")

allCap = {
    year: (
        group.groupby(["Technology", "Market_Node"])["Value"]
             .sum()
             .unstack(fill_value=0)  
    )
    for year, group in filtered_rows.groupby("Target year")
}

#distributed_resources = ["onwind", "offwind", "CSP","solar", "battery"]
dispatchable = ["CCGT", "OCGT", 'biomass', 'coal','lignite', 'nuclear','oil','hydrogen']
technologies_for_investment = snakemake.config["power_plants"]["technologies_new_investment"]

[capacity,investcap] = build_capacity_table(allCap)

investcap.to_hdf(snakemake.output.investcap, "investcap")

existing = build_legacy_caps()
new = build_new_investments()

dispatchable_plants = pd.concat([existing, new])

#dispatchable_plants.loc[(dispatchable_plants.bus == "CY00") & (dispatchable_plants.carrier == "CCGT"), "p_min_pu"] = 0 # lift p_min_pu in cyprus as it can lead to infeasibilities

dispatchable_plants.to_hdf(snakemake.output.dispatchable_capacities, "dispatchable")


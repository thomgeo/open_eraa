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


# In[43]:

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
    
    existing["efficiency"] = properties["Standard efficiency in NCV terms"].loc[:, "existing"].reindex(existing.carrier).values

    existing["start_up_cost"] = properties["Start-up fix cost (e.g. wear) warm start"].loc[:, "existing"].reindex(existing.carrier, level=1).fillna(0).values
    existing["ramp_limit_up_real"] = properties["Ramp up rate % of max output power / min"].loc[:, "existing"].reindex(existing.carrier).values*60
    existing["ramp_limit_up"] = existing["ramp_limit_up_real"].clip(upper=1.)
    existing["ramp_limit_down_real"] = properties["Ramp down rate % of max output power / min"].loc[:, "existing"].reindex(existing.carrier).values*60
    existing["ramp_limit_down"] = existing["ramp_limit_down_real"].clip(upper=1.0)
    existing["p_min_pu"] = properties["Minimum stable generation (% of max power)"].loc[:, "existing"].reindex(existing.carrier).values
    existing["min_up_time"] = properties["Min Time on"].loc[:, "existing"].reindex(existing.carrier).values
    existing["min_down_time"] = properties["Min Time off"].loc[:, "existing"].reindex(existing.carrier).values
    existing["var_om"] = properties.loc[:, "existing",:].reindex(existing.carrier)["Variable O&M cost"].values
    existing["p_nom_max"] = existing["p_nom"]
    existing["p_nom_min"] = 0.01
    existing["invest_status"] = "existing"
        
    return existing

# In[4]:

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

    new["efficiency"] = properties.loc[:, "new", :]["Standard efficiency in NCV terms"].reindex(new.carrier).values

    new["start_up_cost"] = properties["Start-up fix cost (e.g. wear) warm start"].loc[:, "new"].reindex(new.carrier, level=1).fillna(0).values

    new["ramp_limit_up_real"] = properties["Ramp up rate % of max output power / min"].loc[:, "new"].reindex(new.carrier).values*60
    existing["ramp_limit_up"] = existing["ramp_limit_up_real"].clip(upper=1.)
    new["ramp_limit_down_real"] = properties["Ramp down rate % of max output power / min"].loc[:, "new"].reindex(new.carrier).values*60
    existing["ramp_limit_down"] = existing["ramp_limit_down_real"].clip(upper=1.)
    new["p_min_pu"] = properties["Minimum stable generation (% of max power)"].loc[:, "new"].reindex(new.carrier).values
    new["min_up_time"] = properties["Min Time on"].loc[:, "new"].reindex(new.carrier).values
    new["min_down_time"] = properties["Min Time off"].loc[:, "new"].reindex(new.carrier).values
    new["exit"] = np.inf
    new["p_nom_max"] = np.inf
    new["p_nom_min"] = 0.01
    new["var_om"] = properties.loc[:, "new",:].reindex(new.carrier)["Variable O&M cost"].values
    new["invest_status"] = "new"
    
    return new
    
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

years = range(2028, 2036)

#Different methodology to get data 
with zipfile.ZipFile(snakemake.input.pemmdb) as zip_f:
    zip_f.extractall(snakemake.params.data_folder)
        
with zipfile.ZipFile(snakemake.input.thermal) as zip_f:
    zip_f.extractall(snakemake.params.data_folder)
    
with zipfile.ZipFile(snakemake.input.dsr) as zip_f:
    zip_f.extractall(snakemake.params.data_folder)
    
file2 = snakemake.params.data_folder + "Common data/Common Data.xlsx"

all_columns = pd.read_csv(
    snakemake.params.data_folder + "Other data/Explicit DSR detailed.csv",
    nrows=0
).columns.tolist()

# Drop the last column
columns_to_use = all_columns[:-1]

# Read CSV with only the selected columns
properties_raw2 = pd.read_csv(
    snakemake.params.data_folder + "Other data/Explicit DSR detailed.csv",
    usecols=columns_to_use,
    header=0
)

allDSR = {
  year: group.rename(columns={
      "MARKET_NODE": "bus",
      "DA activation price (EUR/MWh)": "marginal_cost",
      "Max hours dispatched per day [h]": "hours",
      "MAX CAP (MW)": "p_nom"
  })
  for year, group in properties_raw2.groupby("TARGET_YEAR")
}

allDSR = pd.concat(allDSR, names=["TARGET_YEAR"]).fillna(0)

allDSR.to_hdf(snakemake.output.dsr, "dsr")

excel_file2 = pd.ExcelFile(file2)

properties_raw = pd.read_csv(snakemake.params.data_folder + "Dashboard_raw_data/GenerationCapacities.csv", header=0)

filtered_rows = properties_raw[
    properties_raw.apply(
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

allCap = {
    year: (
        group.groupby(["Technology", "Market_Node"])["Value"]
             .sum()
    )
    for year, group in filtered_rows.groupby("Target year")
}

# Combine into one DataFrame with hierarchical indices
allCap_df = pd.concat(allCap, names=["Target year"]).fillna(0)

allCap_df.to_hdf(snakemake.output.all_capacities, "capacities")
allCap_df.to_csv('allCapdf.csv', index=False)

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

properties_raw = pd.read_excel(excel_file2, "Common Data", index_col=[2,3], skiprows=10, header=0).dropna(how="all").iloc[2:, 1:].dropna(how="all", axis=1).iloc[:27, 1:17]

properties_raw2 = pd.read_excel(excel_file2, "Common Data", index_col=[2,3], skiprows=44, header=[0,3]).iloc[:, 1:].dropna(how="all", axis=1).dropna(how="all").iloc[:27, 1:17]

properties_raw2.columns = [" ".join(i) for i in properties_raw2.columns]
properties_raw = pd.concat([properties_raw, properties_raw2],axis=1)

mask = properties_raw.index.get_level_values(0).str.contains("fuel cell", case=False) | \
       properties_raw.index.get_level_values(1).str.contains("fuel cell", case=False)
properties_raw = properties_raw[~mask]
properties_raw = properties_raw.fillna(0)

#properties_raw.to_csv('propertiesprev.csv', index=False)

properties = build_thermal_properties(properties_raw)

#properties.to_csv('propertiesfinal.csv', index=False)

[capacity,investcap] = build_capacity_table(allCap)

investcap.to_csv('InvestCap.csv', index=False)
investcap.to_hdf(snakemake.output.investcap, "investcap")

capacity.to_csv('ManipulatedCap.csv', index=False)

existing = build_legacy_caps()
new = build_new_investments()

dispatchable_plants = pd.concat([existing, new])

#dispatchable_plants.loc[(dispatchable_plants.bus == "CY00") & (dispatchable_plants.carrier == "CCGT"), "p_min_pu"] = 0 # lift p_min_pu in cyprus as it can lead to infeasibilities

dispatchable_plants.to_hdf(snakemake.output.dispatchable_capacities, "dispatchable")

#dispatchable_plants.to_csv('output_pandas.csv', index=False)

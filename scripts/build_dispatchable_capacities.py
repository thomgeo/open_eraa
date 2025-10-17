#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import powerplantmatching as pm
import pycountry
import os 
import numpy as np


# In[3]:


def build_capacity_table():

    capacity_years = pd.DataFrame()

    for sheet in plant_sheets.index:

        capacity_raw = pd.read_excel(excel_file, sheet, header=1, index_col=1).dropna(how="all", axis=1).iloc[:20]

        capacity_years[sheet[2:]] = prepare_capacity_table(capacity_raw).stack().astype(float)

    capacity_years.columns = capacity_years.columns.astype(int)
    capacity_years = capacity_years.reindex(years,axis=1)
    capacity_years.loc[distributed_resources, :] = capacity_years.loc[distributed_resources, :].interpolate(axis=1)
    

    return capacity_years.fillna(method="ffill", axis=1)


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
    existing["ramp_limit_up"] = properties["Ramp up rate % of max output power / min"].loc[:, "existing"].reindex(existing.carrier).values*60
    existing["ramp_limit_down"] = properties["Ramp down rate % of max output power / min"].loc[:, "existing"].reindex(existing.carrier).values*60
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

    new["ramp_limit_up"] = properties["Ramp up rate % of max output power / min"].loc[:, "new"].reindex(new.carrier).values*60
    new["ramp_limit_down"] = properties["Ramp down rate % of max output power / min"].loc[:, "new"].reindex(new.carrier).values*60
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

years = range(2025, 2034)


excel_file = pd.ExcelFile(snakemake.input.pemmdb)


plant_sheets = [i for i in excel_file.sheet_names if "TY" in i]
plant_sheets = pd.Series([int(i[2:]) for i in plant_sheets], plant_sheets )


distributed_resources = ["onwind", "offwind", "CSP","solar", "battery"]
dispatchable = ["CCGT", "OCGT", 'biomass', 'coal', 'lignite', 'nuclear','oil',  'other', ]
technologies_for_investment = snakemake.config["investment"]["technologies_new_investment"]

properties_raw = pd.read_excel(excel_file, "Thermal Properties", index_col=[2,3], header=3).dropna(how="all").iloc[1:, 2:].dropna(how="all", axis=1).iloc[:24, :12]
properties_raw2 = pd.read_excel(excel_file, "Thermal Properties", index_col=[2,3], skiprows=35, header=[0,3]).iloc[:, 2:].dropna(how="all", axis=1).dropna(how="all")
properties_raw2.columns = [" ".join(i) for i in properties_raw2.columns]
properties_raw = pd.concat([properties_raw, properties_raw2],axis=1)

properties = build_thermal_properties(properties_raw)

capacity = build_capacity_table()

existing = build_legacy_caps()
new = build_new_investments()

dispatchable_plants = pd.concat([existing, new])

dispatchable_plants.loc[(dispatchable_plants.bus == "CY00") & (dispatchable_plants.carrier == "CCGT"), "p_min_pu"] = 0 # lift p_min_pu in cyprus as it can lead to infeasibilities

dispatchable_plants.to_hdf(snakemake.output.dispatchable_capacities, "dispatchable")

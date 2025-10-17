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
        commissioning["exit"] = years[-1] + 1
        commissioning["carrier"] = i[0]
        commissioning["bus"] = i[1]

        existing = pd.concat([existing, legacy, commissioning])
        existing = existing[["bus", "p_nom", "carrier", "entry", "exit"]]
        existing = existing[existing.p_nom >0]

    existing.exit.replace(years[-1] +1, np.nan)
    
    return existing


# In[4]:


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
technologies_for_investment = ["OCGT", "CCGT"]


capacity = build_capacity_table()

existing = build_legacy_caps()

existing.to_hdf(snakemake.output.legacy_capacities, "legacy")

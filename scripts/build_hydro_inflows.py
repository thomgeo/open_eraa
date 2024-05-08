#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import zipfile
import glob
import os 

unzip_folder = snakemake.params.unzip_folder
save_hdf = snakemake.output.inflow



with zipfile.ZipFile(snakemake.input.pecd22) as zip_f:
    zip_f.extract("Climate Data/Hydro Inflows.zip", unzip_folder)

with zipfile.ZipFile(unzip_folder + "Climate Data/Hydro Inflows.zip") as zip_f:
    zip_f.extractall(unzip_folder)


files = glob.glob(unzip_folder + "/Hydro Inflows/*")

inflow = pd.DataFrame()
for file in files:

    print(file)
    
    excel_file = pd.ExcelFile(file)
    inflow_zone = pd.DataFrame()

    for sheet in ['Run of River', 'Pondage', 'Reservoir', 'Pump storage - Open Loop', 'Pump Storage - Closed Loop']:

        inflow_zone_tech = pd.read_excel(file, sheet, header=12).iloc[:,17:17+35]
        inflow_zone_tech.columns = range(1982, 2017)
        inflow_zone_tech.columns = inflow_zone_tech.columns.astype(str)
        inflow_zone_tech.index = inflow_zone_tech.index*24*7

        scaler = pd.Series(inflow_zone_tech.index.to_list() + [8760]).diff().shift(-1).dropna()
        scaler.index *=(24*7)

        inflow_zone_hourly = inflow_zone_tech.reindex(range(8760)).fillna(method="ffill")
        inflow_zone_hourly = inflow_zone_hourly.div(scaler.reindex(inflow_zone_hourly.index).fillna(method="ffill"), axis=0)    
        inflow_zone[sheet] = inflow_zone_hourly.stack()

    inflow_zone = inflow_zone.stack()
    inflow_zone.index.names = ["hour", "climate year", "technology"]
    inflow_zone.name = (file.split("_")[1], file.split("_")[-1][:4])

    inflow = pd.concat([inflow, inflow_zone], axis=1)


# In[19]:


inflow = inflow.T.set_index((i for i in inflow.columns)).T.stack(0)

dirname = os.path.dirname(save_hdf)

if not os.path.exists(dirname):
    os.mkdir(dirname)
    
    
inflow.stack().to_hdf(save_hdf, "inflow")


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import glob
import zipfile
import os


save_hdf=snakemake.output.demand


with zipfile.ZipFile(snakemake.input.demand_files) as zip_f:
    zip_f.extractall(snakemake.params.demand_folder)

files = glob.glob(snakemake.params.demand_folder + "Demand data/Demand timeseries/*")


#demand = [i for i in files if "Demand" in i]

demand = pd.DataFrame()

# structural change to collect 2024 data
for file in files:
    
    print("\n" + file + "\n")

    #demand_excel = pd.ExcelFile(file)

    demand_target_year = pd.DataFrame()

    zone = os.path.splitext(os.path.basename(file))[0].split("_")[0]
        
    demand_target_year[zone] = pd.read_csv(file, index_col=[0,3], header=0).reset_index(drop=True).stack()

    demand_target_year = demand_target_year.stack()

    target_year = os.path.splitext(os.path.basename(file))[0].split("_")[3]

    demand_target_year.name=target_year
    
    demand = pd.concat([demand, demand_target_year], axis=1)

demand = demand.stack()

demand.index.names = ["hour", "climate year", "zone", "target year"]

dirname = os.path.dirname(save_hdf)

if not os.path.exists(dirname):
    os.mkdir(dirname)

demand.to_hdf(save_hdf, "demand")

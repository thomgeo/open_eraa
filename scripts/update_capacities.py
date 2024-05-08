#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import pypsa
import os


save_path = snakemake.output.new_network

n = pypsa.Network(snakemake.input.old_network)

target_year = int(snakemake.params.ty)
climate_year = int(snakemake.params.cy)

capacity_table = pd.read_csv(snakemake.input.capacity_table, index_col=[0,1])
capacity_table = capacity_table.loc[:, target_year,:]
update = capacity_table.index.intersection(n.generators.index)
n.generators.loc[update, "p_nom"] = capacity_table.loc[update, "p_nom"]

os.makedirs(os.path.dirname(save_path), exist_ok=True)

print(save_path)

n.export_to_netcdf(save_path)

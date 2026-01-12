#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import zipfile
import glob
import os

def get_ntc_data(df):
    
    p_nom = df.max().reset_index([1,2])
    p_nom.columns = ["bus0", "bus1", "p_nom"]
    p_max_pu = df.div(df.max().add(1e-6))
    p_max_pu = p_max_pu.T.reset_index([1,2], drop=True).T
    
    return p_nom, p_max_pu

def get_ntc_tables(excel_file):
    
    hvac = pd.read_excel(excel_file, "HVAC", header=[7, 8, 9], index_col=[0, 1]).iloc[5:]
    hvdc = pd.read_excel(excel_file, "HVDC", header=[7, 8, 9], index_col=[0, 1]).iloc[5:]

    p_nom_ac, p_max_pu_ac = get_ntc_data(hvac)
    p_nom_ac.index = p_nom_ac.index + "-AC"
    p_max_pu_ac.columns = p_max_pu_ac.columns + "-AC"
    
    p_nom_ac = p_nom_ac[~p_nom_ac.index.duplicated(keep="first")]
    p_max_pu_ac = p_max_pu_ac.loc[:, ~p_max_pu_ac.columns.duplicated(keep="first")]

    p_nom_dc, p_max_pu_dc = get_ntc_data(hvdc)
    p_nom_dc.index = p_nom_dc.index + "-DC"
    p_max_pu_dc.columns = p_max_pu_dc.columns + "-DC"

    p_nom = pd.concat([p_nom_ac, p_nom_dc])
    p_max_pu = pd.concat([p_max_pu_ac, p_max_pu_dc], axis=1)

    p_nom = p_nom[p_nom.p_nom >0]
    p_max_pu = p_max_pu[p_nom.index]
    
    return p_nom, p_max_pu

def split_uk_france_interconnection(p_nom, p_max_pu):
    
    p_nom_uk_fr = p_nom.loc["UK00-FR00-DC"].copy()
    p_nom_uk_fr.loc["p_nom"] = 2000
    p_nom_uk_fr = pd.concat([p_nom_uk_fr, p_nom_uk_fr], keys=["UK00-FR00_1-DC", "UK00-FR00_2-DC"])
    p_nom.drop("UK00-FR00-DC", inplace=True)
    p_nom = pd.concat([p_nom, p_nom_uk_fr]).sort_index()
    
    p_max_pu_uk_fr = p_max_pu.loc[:, :, ["UK00-FR00-DC"], :].copy()
    
    p_max_pu_uk_fr.index = p_max_pu_uk_fr.index.remove_unused_levels()
    
    p_max_pu_uk_fr_1 = p_max_pu_uk_fr.copy()
    p_max_pu_uk_fr_2 = p_max_pu_uk_fr.copy()
    
    p_max_pu_uk_fr_1.index = p_max_pu_uk_fr_1.index.set_levels(["UK00-FR00_1-DC"], level=2)
    p_max_pu_uk_fr_2.index = p_max_pu_uk_fr_2.index.set_levels(["UK00-FR00_2-DC"], level=2)
    
    p_max_pu = pd.concat([
        p_max_pu.drop("UK00-FR00-DC", level=2),
        pd.concat([p_max_pu_uk_fr_1, p_max_pu_uk_fr_2])
    ]).sort_index()

    p_nom_fr_uk = p_nom.loc["FR00-UK00-DC"].copy()
    p_nom_fr_uk.loc["p_nom"] = 2000
    p_nom_fr_uk = pd.concat([p_nom_fr_uk, p_nom_fr_uk], keys=["FR00-UK00_1-DC", "FR00-UK00_2-DC"])
    p_nom.drop("FR00-UK00-DC", inplace=True)
    p_nom = pd.concat([p_nom, p_nom_fr_uk]).sort_index()
    
    p_max_pu_fr_uk = p_max_pu.loc[:, :, ["FR00-UK00-DC"], :].copy()
    
    p_max_pu_fr_uk.index = p_max_pu_fr_uk.index.remove_unused_levels()
    
    p_max_pu_fr_uk_1 = p_max_pu_fr_uk.copy()
    p_max_pu_fr_uk_2 = p_max_pu_fr_uk.copy()
    
    p_max_pu_fr_uk_1.index = p_max_pu_fr_uk_1.index.set_levels(["FR00-UK00_1-DC"], level=2)
    p_max_pu_fr_uk_2.index = p_max_pu_fr_uk_2.index.set_levels(["FR00-UK00_2-DC"], level=2)
    
    p_max_pu = pd.concat([
        p_max_pu.drop("FR00-UK00-DC", level=2),
        pd.concat([p_max_pu_fr_uk_1, p_max_pu_fr_uk_2])
    ]).sort_index()
    

    return p_nom, p_max_pu    
folder = snakemake.params.folder
ntc = snakemake.input.ntc

save_hdf = snakemake.output.save_hdf

with zipfile.ZipFile(ntc) as zip_f:
    zip_f.extractall(folder)

p_nom = pd.DataFrame()
p_max_pu = pd.DataFrame()

for file in sorted([i for i in glob.glob("data/NTCs/*") if "NTCs" in i]):

    excel_file = pd.ExcelFile(file)

    p_nom_year, p_max_pu_year = get_ntc_tables(excel_file)

    p_nom_year, p_max_pu_year = get_ntc_tables(excel_file)
    
    p_nom_year = p_nom_year.stack()
    p_nom_year.name = int(os.path.splitext(file.split("TY")[-1])[0])
    p_nom = pd.concat([p_nom, p_nom_year], axis=1)
    
    p_max_pu_year = p_max_pu_year.stack()
    p_max_pu_year.name= int(os.path.splitext(file.split("TY")[-1])[0])
    p_max_pu = pd.concat([p_max_pu, p_max_pu_year],axis=1)
    
for col in p_max_pu.columns:
    p_max_pu[col] = p_max_pu[col].astype(float)

p_nom = p_nom.set_index(i for i in p_nom.index)
p_max_pu = p_max_pu.set_index((i for i in p_max_pu.index))

p_nom, p_max_pu = split_uk_france_interconnection(p_nom, p_max_pu)

p_nom.to_hdf(save_hdf, "p_nom")
p_max_pu.to_hdf(save_hdf, "p_max_pu")





#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import zipfile
import os

def extract_FB_data(excel_file):

    RAM_sheets = [i for i in excel_file.sheet_names if "RAM" in i]
    PTDF_sheets = [i for i in excel_file.sheet_names if "PTDF" in i]
    
    RAM = []
    for sheet in RAM_sheets:
        RAM.append(pd.read_excel(excel_file, sheet, index_col=0))
    
    RAM = pd.concat(RAM, keys=RAM_sheets)
    RAM.index = RAM.index.set_levels(RAM.index.levels[0].str.replace("RAM ", "").astype(int), level=0)
    RAM.index.names = ["year", "CNEC"]
    
    PTDF = []
    for sheet in PTDF_sheets:
        PTDF.append(pd.read_excel(excel_file, sheet, index_col=[0,1], header=[0,1]))
    
    PTDF = pd.concat(PTDF, keys = PTDF_sheets).iloc[:, 1:]
    PTDF.index = PTDF.index.set_levels(PTDF.index.levels[0].str.replace("PTDF ", "").astype(int), level=0)
    
    domain_assignment = pd.read_excel(excel_file, "Domain Assignment", index_col=[0,1,2,3])

    return RAM, PTDF, domain_assignment


fb_zip = snakemake.input.zip_file
fb_folder=  snakemake.params.folder

year = 2030
core_data = snakemake.output.core_domain
nordic_data = snakemake.output.nordic_domain

with zipfile.ZipFile(fb_zip) as zip_f:
    zip_f.extractall(fb_folder)

core = pd.ExcelFile(fb_folder + 'FB domains/FB-Domain-CORE_Merged.xlsx')
nordic = pd.ExcelFile(fb_folder + 'FB domains/FB-Domain-NORDIC_Merged.xlsx')

RAM_nordic, PTDF_nordic, domain_assignment_nordic = extract_FB_data(nordic)
RAM_core, PTDF_core, domain_assignment_core = extract_FB_data(core)
RAM_nordic.drop("CNEC_Name", axis=1, inplace=True)

os.makedirs(os.path.dirname(core_data), exist_ok=True)


RAM_core.loc[year].to_hdf(core_data, key="RAM")
PTDF_core.loc[year].to_hdf(core_data, key="PTDF")
domain_assignment_core.to_hdf(core_data, key="domain_assignment")

RAM_nordic.loc[year].to_hdf(nordic_data, key="RAM")
PTDF_nordic.loc[year].to_hdf(nordic_data, key="PTDF")
domain_assignment_nordic.to_hdf(nordic_data, key="domain_assignment")

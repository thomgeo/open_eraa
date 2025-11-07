#!/usr/bin/env python
# coding: utf-8

# In[10]:


import pandas as pd
import numpy as np
import glob
import zipfile
import os

with zipfile.ZipFile(snakemake.input.pecd) as zip_f:
    zip_f.extractall(snakemake.params.data_folder)
    
save_hdf = snakemake.output.res_profile

#line changed to 2025
all_files = glob.glob(snakemake.params.data_folder + "PECD - RES/Capacity Factors_250716/*")

onshore = [i for i in all_files if "Onshore" in i]
offshore = [i for i in all_files if "Offshore" in i]
PVtrack = [i for i in all_files if "tracking" in i]
PVfixed = [i for i in all_files if "fixed" in i]
PVRroof = [i for i in all_files if "residential" in i]
PVIroof = [i for i in all_files if "industrial" in i]
csp = [i for i in all_files if "CSP_no" in i]
cspS = [i for i in all_files if "7h_dis" in i]

#missing_csp = ["GR03","PT00"]
#csp = csp + [i for i in all_files if "CSP_withStorage_7h_preDispatch" in i and any(zone in i for zone in missing_csp)]
# previous instruction changed
#csp = csp + [i for i in all_files if "CSP_withStorage_7h_dispatched" in i and any(zone in i for zone in missing_csp)]

def get_availabilities(files):
    
    availability_factors = pd.DataFrame()

    for file in files:

        info = pd.read_csv(file, index_col=0, header=1).iloc[:3,0]

        availability_factors_zone = pd.read_csv(file, index_col=[0,1], header=10).reset_index(drop=True).stack()

        availability_factors_zone.name = (info["Target Year"], info["PECD Zone"])

        availability_factors = pd.concat([availability_factors, availability_factors_zone], axis=1)

    availability_factors = availability_factors.T.set_index((i for i in availability_factors.columns)).T.stack([0,1])

    availability_factors.index.names = ["hour", "climate year", "target year", "zone"]

    return availability_factors

onshore_availabilities = get_availabilities(onshore)
PVtrack_availabilities = get_availabilities(PVtrack)
PVfixed_availabilities = get_availabilities(PVfixed)
PVRroof_availabilities = get_availabilities(PVRroof)
PVIroof_availabilities = get_availabilities(PVIroof)
offshore_availabilities = get_availabilities(offshore)
csp_availabilities = get_availabilities(csp)
csps_availabilities=get_availabilities(cspS)

dirname = os.path.dirname(save_hdf)

if not os.path.exists(dirname):
    os.mkdir(dirname)

onshore_availabilities.to_hdf(save_hdf, "onwind")
offshore_availabilities.to_hdf(save_hdf, "offwind")
PVtrack_availabilities.to_hdf(save_hdf, "PVtrack")
PVfixed_availabilities.to_hdf(save_hdf, "PVfixed")
PVRroof_availabilities.to_hdf(save_hdf, "PVRroof")
PVIroof_availabilities.to_hdf(save_hdf, "PVIroof")
csp_availabilities.to_hdf(save_hdf, "CSP")
csps_availabilities.to_hdf(save_hdf, "CSPS")



#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import zipfile
import glob
import os 
import time 

save_hdf = snakemake.output.inflow

all_files = glob.glob(snakemake.params.data_folder + "PECD - RES/Hydro Inflows_250704/*")

HRI = [i for i in all_files if "HRI" in i]
HRR = [i for i in all_files if "HRR" in i]
HOL = [i for i in all_files if "HOL" in i]
HPI = [i for i in all_files if "HPI" in i]

def get_inflows(files,tech):

   inflow_zone = pd.DataFrame()

   for file in files:
 
  
      inflow_zone_tech = pd.read_csv(file, header=0).iloc[:,1:]    
      inflow_zone_tech.columns = range(1, 37)
      inflow_zone_tech.columns = inflow_zone_tech.columns.astype(str)
      scale=8760/len(inflow_zone_tech.index)
      
      inflow_zone_hourly = pd.DataFrame(index=np.arange(8760), columns=inflow_zone_tech.columns, dtype=float)

      for i in range(len(inflow_zone_tech.index)):
          start = int(i * scale)
          end = int((i + 1) * scale)
          inflow_zone_hourly.iloc[start:end] = inflow_zone_tech.iloc[i] / (end - start)
      
      inflow_zone_hourly.fillna(0, inplace=True)   

      inflow_zone_hourly=inflow_zone_hourly.stack();
      inflow_zone_hourly.name = (file.split("/")[-1].split("_")[0], os.path.splitext(file.split("_")[-1])[0]) 
   
      inflow_zone = pd.concat([inflow_zone, inflow_zone_hourly], axis=1)
      
   
   inflow_zone = inflow_zone.T.set_index((i for i in inflow_zone.columns)).T.stack([0,1])
   inflow_zone.index.names = ["hour", "climate year", "target year", "zone"]
   inflow_zone.name = tech 
      
   return inflow_zone
   
inflow = pd.DataFrame()


HRI_inflow = get_inflows(HRI, "hydro")
inflow = pd.concat([inflow, HRI_inflow], axis=1)
HRR_inflow = get_inflows(HRR, "ROR")
inflow = pd.concat([inflow, HRR_inflow], axis=1)
HOL_inflow = get_inflows(HOL, "PHS open")
inflow = pd.concat([inflow, HOL_inflow], axis=1)
HPI_inflow = get_inflows(HPI, "pondage")
inflow = pd.concat([inflow, HPI_inflow], axis=1)

inflow = inflow.T.set_index((i for i in inflow.columns)).T.stack(0)

dirname = os.path.dirname(save_hdf)

if not os.path.exists(dirname):
    os.mkdir(dirname)

inflow.to_hdf(save_hdf2, "inflow")  
HRR_inflow[:,"1","AT00","2030"].to_csv('scaledinflow.csv', index=False)     



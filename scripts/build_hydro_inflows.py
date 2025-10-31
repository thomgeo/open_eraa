#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import zipfile
import glob
import os 

save_hdf = snakemake.output.inflow

all_files = glob.glob(snakemake.params.data_folder + "PECD - RES/Hydro Inflows_250704/*")

HRI = [i for i in all_files if "HRI" in i]
HRR = [i for i in all_files if "HRR" in i]
HOL = [i for i in all_files if "HOL" in i]
HPI = [i for i in all_files if "HPI" in i]

def get_inflows(files):

   inflow_zone = pd.DataFrame()

   for file in files:

    
  
      inflow_zone_tech = pd.read_csv(file, header=0).iloc[:,1:]
      
        
      
      inflow_zone_tech.columns = range(1, 37)
      inflow_zone_tech.columns = inflow_zone_tech.columns.astype(str)
      #print(inflow_zone_tech.index)
      inflow_zone_tech.index = inflow_zone_tech.index*24*7
      
      scaler = pd.Series(inflow_zone_tech.index.to_list() + [8760]).diff().shift(-1).dropna()
      scaler.index *=(24*7)
   
      inflow_zone_hourly = inflow_zone_tech.reindex(range(8760)).fillna(method="ffill")
      #print(inflow_zone_hourly.index)
      inflow_zone_hourly = inflow_zone_hourly.div(scaler.reindex(inflow_zone_hourly.index).fillna(method="ffill"), axis=0)    

      inflow_zone_hourly=inflow_zone_hourly.stack();
      inflow_zone_hourly.name = (file.split("/")[-1].split("_")[0], os.path.splitext(file.split("_")[-1])[0]) 
   
      inflow_zone = pd.concat([inflow_zone, inflow_zone_hourly], axis=1)
      
   
   inflow_zone = inflow_zone.T.set_index((i for i in inflow_zone.columns)).T.stack([0,1])
   inflow_zone.index.names = ["hour", "climate year", "target year", "zone"]
   #print(inflow_zone.index.nlevels)   
      
   return inflow_zone
   
inflow = pd.DataFrame()


HRI_inflow = get_inflows(HRI)
#inflow = pd.concat([inflow, HRI_inflow], axis=1)
HRR_inflow = get_inflows(HRR)
#inflow = pd.concat([inflow, HRR_inflow], axis=1)
HOL_inflow = get_inflows(HOL)
#inflow = pd.concat([inflow, HOL_inflow], axis=1)
HPI_inflow = get_inflows(HPI)
#inflow = pd.concat([inflow, HPI_inflow], axis=1)

inflow = inflow.T.set_index((i for i in inflow.columns)).T.stack(0)

dirname = os.path.dirname(save_hdf)

if not os.path.exists(dirname):
    os.mkdir(dirname)
      
HRI_inflow.to_hdf(save_hdf, "hydro")
HRR_inflow.to_hdf(save_hdf, "ROR")
HOL_inflow.to_hdf(save_hdf, "PHS open")
HPI_inflow.to_hdf(save_hdf, "pondage")



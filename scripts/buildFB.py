# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 11:00:58 2025

@author: lourehu
"""

import pandas as pd
import numpy as np
import glob
import zipfile
import pypsa
import os
import xlsxwriter
import time


def build_FB_table(excel_file,excel_file2):

    ptdf_corey = pd.DataFrame()
    ptdf_ahcy = pd.DataFrame()
    ptdf_dcy = pd.DataFrame()
    ramy = pd.DataFrame()
    
    core2= pd.DataFrame()
    ptdf2= pd.DataFrame()
    ptdf2= pd.DataFrame()
    
    ptdf_sheets = [i for i in excel_file.sheet_names if "PTDF" in i]
    #print(ptdf_sheets)
    ptdf_sheets = pd.Series([int(i[5:]) for i in ptdf_sheets], ptdf_sheets) 
    years= range(int(ptdf_sheets[0]),int(ptdf_sheets[len(ptdf_sheets)-1]+1))


    for idx,sheet in enumerate(ptdf_sheets.index):
        #print(idx)
        #time.sleep(2)
                   
        ptdf_corey[sheet], ptdf_ahcy[sheet], ptdf_dcy[sheet], domain_assign = [
            df.stack().astype(str) for df in build_flow_based_data(excel_file, sheet,years[idx])
        ]
        print("passou assign1")
        time.sleep(1)
        
        
        core2, ptdf2, ptdf2, domain_assign2 = [
            df.stack().astype(float) for df in build_flow_based_data(excel_file2, sheet,years[idx])
        ]
        
        ptdf_corey[sheet] = pd.concat([float(ptdf_corey[sheet]), core2])
        ptdf_ahcy[sheet]  = pd.concat([float(ptdf_ahcy[sheet]), ahc2])
        ptdf_dcy[sheet]   = pd.concat([float(ptdf_dcy[sheet]), dc2])
        print("passou assign2")
        time.sleep(1)
    
    print("vai para o RAM")    
    ram_sheets = [i for i in excel_file.sheet_names if "RAM" in i]
    ram_sheets = pd.Series([int(i[4:]) for i in ram_sheets], ram_sheets) 
    
    for sheet in ram_sheets.index:    
        ram = pd.read_excel(excel_file, sheet, index_col=[0], skiprows=1)
        ram.columns = ram.columns.astype(int)
        ram = ram.T.reindex(domain_assign.values)
        ram.index = n.snapshots
        ram = ram.stack()
        ram.index.names = ["snapshot", "CNEC"]
        ram2 = pd.read_excel(excel_file2, sheet, index_col=[0], skiprows=1)
        ram2.columns = ram.columns.astype(int)
        ram2 = ram2.T.reindex(domain_assign2.values)
        ram2.index2 = n.snapshots
        ram2 = ram2.stack()
        ram2.index.names = ["snapshot", "CNEC"]
                
        ramy[sheet] = pd.concat([ram, ram2])
    
    print("passou RAM")
    ptdf_corey.columns = ptdf_corey.columns.astype(int)
    ptdf_ahcy.columns = ptdf_ahcy.columns.astype(int)
    ptdf_dcy.columns = ptdf_dcy.columns.astype(int)
    ramy.columns = ramy.columns.astype(int)
    ptdf_corey.columns = ptdf_corey.reindex(years,axis=1)
    ptdf_ahcy.columns = ptdf_ahcy.reindex(years,axis=1)
    ptdf_dcy.columns = ptdf_dcy.reindex(years,axis=1)
    ramy.columns = ramy.columns.reindex(years,axis=1)
    #capacity_years.columns = capacity_years.columns.astype(int)
    #capacity_years = capacity_years.reindex(years,axis=1)
    #capacity_years.loc[distributed_resources, :] = capacity_years.loc[distributed_resources, :].interpolate(axis=1)

    return [ptdf_core.fillna(method="ffill", axis=1)[year].unstack(1),ptdf_ahc.fillna(method="ffill", axis=1)[year].unstack(1),ptdf_dc.fillna(method="ffill", axis=1)[year].unstack(1),ram.fillna(method="ffill", axis=1)[year].unstack(1)]



def build_flow_based_data(excel_file, ptdf_sheet,year):
 
    
    #print(ptdf_sheets)
    
    #base_years = pd.Series([int(i.split(" ")[1]) for i in ptdf_sheets], ptdf_sheets) - snakemake.params.target_year
    
    #print(base_years)
    
    #ptdf_sheet = base_years[base_years <=0].idxmax()

    ptdf_core = pd.read_excel(excel_file, ptdf_sheet, index_col=[0,1,2], header=[0,1])["PTDF_SZ"]
    ptdf_ahc = pd.read_excel(excel_file, ptdf_sheet, index_col=[0,1,2], header=[0,1])["PTDF*_AHC,SZ"]
    ptdf_dc = pd.read_excel(excel_file, ptdf_sheet, index_col=[0,1,2], header=[0,1])["PTDF_EvFB"]
  
    
    domain_assign = pd.read_excel(excel_file, "Domain Assignment", index_col=None)
    
    domain_assigny = (
        domain_assign[domain_assign["Year"] == year]
        .sort_values(by=["Month", "Day", "Hour"])
        .head(8761)
        .reset_index(drop=True)
    )
    
    new_header = domain_assigny.iloc[0]  
    domain_assigny = domain_assigny[1:]  
    domain_assigny = domain_assigny.reset_index(drop=True)
    domain_assigny.columns = new_header  
    
    domain_assigny = domain_assigny.iloc[:, 4:]
    
    domain_assigny.to_csv("domain_assigny.values.csv")
    
    
    domain_assigny.columns = [
        int(col.replace("WS", "")) if isinstance(col, str) and col.startswith("WS") else col
        for col in domain_assigny.columns
    ]

    domain_assigny.to_csv("domain_assigny.values1.csv", index=False)
    
    print(len(ptdf_core))
    
    ptdf_core.index.to_frame().to_csv("ptdfindex.csv", index=False)
    ptdf_core = ptdf_core.reset_index(2,drop=True).unstack(1).reindex(domain_assigny.values)
    
    time.sleep(2)
    ptdf_core.index = n.snapshots
    ptdf_core = ptdf_core.stack(1)
    
    domain_assigny.to_csv("domain_assigny.values2.csv", index=False)
    

    ptdf_ahc = ptdf_ahc.reset_index(2,drop=True).unstack(1).reindex(domain_assigny.values)
    ptdf_ahc.index = n.snapshots
    ptdf_ahc = ptdf_ahc.stack(1)

    ptdf_dc = ptdf_dc.reset_index(2, drop=True).unstack(1).reindex(domain_assigny.values)
    ptdf_dc.index = n.snapshots
    ptdf_dc = ptdf_dc.stack(1)


    ptdf_core.index.names = ["snapshot", "CNEC"]
    ptdf_core.columns.name = "bus"
    ptdf_ahc.index.names = ["snapshot", "CNEC"]
    ptdf_ahc.columns.name = "Flow"
    ptdf_dc.index.names = ["snapshot", "CNEC"]
    
    
    ptdf_core = ptdf_core.reindex(
        pd.MultiIndex.from_product(
            [ptdf_core.index.levels[0], ptdf_core.index.levels[1]]
        ),
        fill_value=0
    )
    print("passou reindex2")
    ptdf_core.to_csv("ptdf.csv", index=False)
    
    ptdf_ahc = ptdf_ahc.reindex(
        pd.MultiIndex.from_product(
            [ptdf_ahc.index.levels[0], ptdf_ahc.index.levels[1]]
        ),
        fill_value=0
    )
    ptdf_dc = ptdf_dc.reindex(
        pd.MultiIndex.from_product(
            [ptdf_dc.index.levels[0], ptdf_dc.index.levels[1]]
        ), 
        fill_value=0
    )
    
    ptdf_ahc = ptdf_ahc[ptdf_ahc.sum()[ptdf_ahc.sum()!=0].index]
    ptdf_ahc = ptdf_ahc[[i for i in ptdf_ahc.columns if i.split("-")[0][:4] in n.buses.index and i.split("-")[1][:4] in n.buses.index]]
    ptdf_ahc = ptdf_ahc.groupby(ptdf_ahc.columns.str[:9],axis=1).sum()
    
    #print("passou ahc")
    #time.sleep(2)
    
    return ptdf_core, ptdf_ahc, ptdf_dc, domain_assigny

n = pypsa.Network(str(snakemake.input.network))
#snapshots = pd.date_range(start="2010-01-01", freq="h", periods=8760)
#n.set_snapshots(snapshots)

with zipfile.ZipFile(snakemake.input.FB_files) as zip_f:
    zip_f.extractall(snakemake.params.data_folder)
                            
save_hdf = snakemake.output.FB   
    
#fb_domainsC = glob.glob(snakemake.params.data_folder + "FB domains/FB-Domain-CORE_Merged.xlsx")

flow_based = pd.ExcelFile(snakemake.params.data_folder + "FB domains/FB-Domain-CORE_Merged.xlsx")
flow_basedN = pd.ExcelFile(snakemake.params.data_folder + "FB domains/FB-Domain-NORDIC_Merged.xlsx")

ptdf_core, ptdf_ahc, ptdf_dc, ram = build_FB_table(flow_based,flow_basedN)

# Create a new Excel file and add a worksheet
workbook = xlsxwriter.Workbook('FB.xlsx')
worksheet = workbook.add_worksheet('Core')
worksheet.write(ptdf_core)
worksheet = workbook.add_worksheet('AHC')
worksheet.write(ptdf_ahc)
worksheet = workbook.add_worksheet('DC')
worksheet.write(ptdf_dc)
worksheet = workbook.add_worksheet('RAM')
worksheet.write(ram)
workbook.close()

dirname = os.path.dirname(save_hdf)

if not os.path.exists(dirname):
    os.mkdir(dirname)

ptdf_core.to_hdf(save_hdf, "Core")
ptdf_ahc.to_hdf(save_hdf, "AHC")
ptdf_dc.to_hdf(save_hdf, "DC")
ram.to_hdf(save_hdf, "RAM")


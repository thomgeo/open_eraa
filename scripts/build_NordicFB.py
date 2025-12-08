import pandas as pd
import zipfile
import glob
import os
import pypsa
import time

def build_flow_based_data(excel_file):

    ptdf_sheets = [i for i in excel_file.sheet_names if "PTDF" in i]

    ptdf_core = pd.read_excel(excel_file, ptdf_sheets[0], index_col=[0,1,2], header=[0,1])["PTDF_SZ"]
    ptdf_ahc = pd.read_excel(excel_file, ptdf_sheets[0], index_col=[0,1,2], header=[0,1])["PTDF*_AHC,SZ"]
    ptdf_dc = pd.read_excel(excel_file, ptdf_sheets[0], index_col=[0,1,2], header=[0,1])["PTDF_EvFB"]
       
    ptdf_core.index.to_frame().to_csv("ptdfindex.csv", index=False)
    
    ptdf_core = ptdf_core.reset_index(2,drop=True).unstack(1).reindex(pd.Series(1, index=range(8760)))
    ptdf_core.index = n.snapshots
    ptdf_core = ptdf_core.stack(1)

    ptdf_ahc = ptdf_ahc.reset_index(2,drop=True).unstack(1).reindex(pd.Series(1, index=range(8760)))
    ptdf_ahc.index = n.snapshots
    ptdf_ahc = ptdf_ahc.stack(1)

    ptdf_dc = ptdf_dc.reset_index(2, drop=True).unstack(1).reindex(pd.Series(1, index=range(8760)))
    ptdf_dc.index = n.snapshots
    ptdf_dc = ptdf_dc.stack(1)

    ptdf_core.index.names = ["snapshot", "CNEC"]
    ptdf_core.columns.name = "bus"
    ptdf_ahc.index.names = ["snapshot", "CNEC"]
    ptdf_ahc.columns.name = "Flow"
    ptdf_dc.index.names = ["snapshot", "CNEC"]
    ptdf_dc.columns.name = "Link"
    
    ptdf_core = ptdf_core.reindex(
        pd.MultiIndex.from_product(
            [ptdf_core.index.levels[0], ptdf_core.index.levels[1]]
        ),
        fill_value=0
    )
    ptdf_core = ptdf_core.rename(lambda x: int(x) * 100000, level="CNEC")
    
    ptdf_ahc = ptdf_ahc.reindex(
        pd.MultiIndex.from_product(
            [ptdf_ahc.index.levels[0], ptdf_ahc.index.levels[1]]
        ),
        fill_value=0
    )
    ptdf_ahc = ptdf_ahc.rename(lambda x: int(x) * 100000, level="CNEC")
    
    
    ptdf_dc = ptdf_dc.reindex(
        pd.MultiIndex.from_product(
            [ptdf_dc.index.levels[0], ptdf_dc.index.levels[1]]
        ), 
        fill_value=0
    )
    ptdf_dc = ptdf_dc.rename(lambda x: int(x) * 100000, level="CNEC")
    #ptdf_dc.columns.to_frame().to_csv("ptdfdc.csv", index=False)
    
    ptdf_ahc = ptdf_ahc[ptdf_ahc.sum()[ptdf_ahc.sum()!=0].index]
    ptdf_ahc = ptdf_ahc[[i for i in ptdf_ahc.columns if i.split("-")[0][:4] in n.buses.index and i.split("-")[1][:4] in n.buses.index]]
    ptdf_ahc = ptdf_ahc.groupby(ptdf_ahc.columns.str[:9],axis=1).sum()
    
    ptdf_dc = ptdf_dc[ptdf_dc.sum()[ptdf_dc.sum()!=0].index]
    ptdf_dc = ptdf_dc[[i for i in ptdf_dc.columns if i.split("-")[0][:4] in n.buses.index and i.split("-")[1][:4] in n.buses.index]]
    ptdf_dc = ptdf_dc.groupby(ptdf_dc.columns.str[:9],axis=1).sum()
    
    return ptdf_core, ptdf_ahc, ptdf_dc
    
folder = snakemake.params.folder
NFB = snakemake.output.NFB
n = pypsa.Network(str(snakemake.input.network))

with zipfile.ZipFile(snakemake.input.fb_domains) as zip_f:
    zip_f.extractall(snakemake.params.folder)
    
fb_domains = pd.ExcelFile(glob.glob(snakemake.params.folder + "FB domains/FB-Domain-NORDIC_Merged.xlsx")[0])

ptdf_core, ptdf_ahc, ptdf_dc = build_flow_based_data(fb_domains)

ptdf_dc.columns.to_frame().to_csv("ptdfdc2.csv", index=False)

dirname = os.path.dirname(NFB)

if not os.path.exists(dirname):
    os.mkdir(dirname)
    
ptdf_core.to_hdf(NFB, "ptdf_core")
ptdf_ahc.to_hdf(NFB, "ptdf_ahc")
ptdf_dc.to_hdf(NFB, "ptdf_dc")


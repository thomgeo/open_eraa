#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import geopandas as gpd
import entsoe
from urllib.error import HTTPError
import zipfile


def load_bidding_zones():
    
    gdfs: list[gpd.GeoDataFrame] = []
    
    for area in entsoe.Area:
        name = area.name
        try:
            url = f"https://raw.githubusercontent.com/EnergieID/entsoe-py/refs/heads/master/entsoe/geo/geojson/{name}.geojson"
            gdfs.append(gpd.read_file(url))
        except HTTPError:
            continue
    
    shapes = pd.concat(gdfs, ignore_index=True)  # type: ignore

    shapes.set_index("zoneName",inplace=True)
    
    return shapes


def align_with_ERAA_classification():
    
    unique_zones = pd.unique([i for i in bidding_zones.index if len(i) == 2])
    
    for zone in unique_zones:
        bidding_zones.loc[zone, "zone"] = zone + "00"
    
    bidding_zones.loc["DK_1", "zone"] = "DKW1"
    bidding_zones.loc["DK_2", "zone"] = "DKE1"
    bidding_zones.loc["IT_CALA", "zone"] = "ITCA"
    bidding_zones.loc["IT_CNOR", "zone"] = "ITCN"
    bidding_zones.loc["IT_CSUD", "zone"] = "ITCS"
    bidding_zones.loc["IT_NORD", "zone"] = "ITN1"
    bidding_zones.loc["IT_SARD", "zone"] = "ITSA"
    bidding_zones.loc["IT_SICI", "zone"] = "ITSI"
    bidding_zones.loc["IT_SUD", "zone"] = "ITS1"
    bidding_zones.loc["NO_4", "zone"] = "NON1"
    bidding_zones.loc["NO_3", "zone"] = "NOM1"
    bidding_zones.loc["NO_2", "zone"] = "NOS2"
    bidding_zones.loc["NO_5", "zone"] = "NOS3"
    bidding_zones.loc["NO_1", "zone"] = "NOS1"
    bidding_zones.loc["SE_1", "zone"] = "SE01"
    bidding_zones.loc["SE_2", "zone"] = "SE02"
    bidding_zones.loc["SE_3", "zone"] = "SE03"
    bidding_zones.loc["SE_4", "zone"] = "SE04"


def extract_shape_files():
    with zipfile.ZipFile(snakemake.input.nuts24) as zip_f:
        zip_f.extract("NUTS_RG_01M_2024_3035_LEVL_0.geojson",snakemake.params.extract_to)
    
    with zipfile.ZipFile(snakemake.input.nuts13) as zip_f:
        zip_f.extract("NUTS_RG_01M_2013_3035_LEVL_1.geojson",snakemake.params.extract_to)

    with zipfile.ZipFile(snakemake.input.moldova) as zip_f:
        zip_f.extract("mda_admbnda_unhcr_20220510_SHP/mda_admbnda_adm0_unhcr_20220510.shp",snakemake.params.extract_to)

    with zipfile.ZipFile("data/moldova_regions.zip") as zip_f:
        zip_f.extract("mda_admbnda_unhcr_20220510_SHP/mda_admbnda_adm0_unhcr_20220510.shp",snakemake.params.extract_to)
        zip_f.extract("mda_admbnda_unhcr_20220510_SHP/mda_admbnda_adm0_unhcr_20220510.sbn",snakemake.params.extract_to)
        zip_f.extract("mda_admbnda_unhcr_20220510_SHP/mda_admbnda_adm0_unhcr_20220510.sbx",snakemake.params.extract_to)
        zip_f.extract("mda_admbnda_unhcr_20220510_SHP/mda_admbnda_adm0_unhcr_20220510.prj",snakemake.params.extract_to)
        zip_f.extract("mda_admbnda_unhcr_20220510_SHP/mda_admbnda_adm0_unhcr_20220510.dbf",snakemake.params.extract_to)
        zip_f.extract("mda_admbnda_unhcr_20220510_SHP/mda_admbnda_adm0_unhcr_20220510.cpg",snakemake.params.extract_to)
        zip_f.extract("mda_admbnda_unhcr_20220510_SHP/mda_admbnda_adm0_unhcr_20220510.shx",snakemake.params.extract_to)


# In[5]:


def add_missing_zones():
    
    moldova = gpd.read_file("{}mda_admbnda_unhcr_20220510_SHP/mda_admbnda_adm0_unhcr_20220510.shp".format(snakemake.params.extract_to))
    
    nuts0 = gpd.read_file("{}NUTS_RG_01M_2024_3035_LEVL_0.geojson".format(snakemake.params.extract_to))
    nuts0 = nuts0.to_crs("EPSG:4326")
    nuts0.set_index("NUTS_ID", inplace=True)
    
    nuts1 = gpd.read_file("{}NUTS_RG_01M_2013_3035_LEVL_1.geojson".format(snakemake.params.extract_to))
    nuts1 = nuts1.to_crs("EPSG:4326")
    nuts1.set_index("NUTS_ID", inplace=True)
    
    bidding_zones.drop("DE_LU", inplace=True)
    bidding_zones.set_index("zone", inplace=True)
    
    bidding_zones.loc["DE00", "geometry"] = nuts0.loc["DE", "geometry"]
    bidding_zones.loc["LUG1", "geometry"] = nuts0.loc["LU", "geometry"]
    bidding_zones.loc["IE00", "geometry"] = nuts0.loc["IE", "geometry"]
    #bidding_zones.loc["TR00", "geometry"] = nuts0.loc["TR", "geometry"]
    bidding_zones.loc["UA00", "geometry"] = nuts0.loc["UA", "geometry"]
    bidding_zones.loc["UK00", "geometry"] = nuts1.query("CNTR_CODE == 'UK'").drop("UKN").unary_union
    bidding_zones.loc["UKNI", "geometry"] = nuts1.loc["UKN"].geometry
    bidding_zones.loc["UK00", "geometry"] = nuts1.filter(like="UK", axis=0).drop("UKN").unary_union
    bidding_zones.loc["MD00", "geometry"] = moldova.loc[0, "geometry"]


bidding_zones = load_bidding_zones()
align_with_ERAA_classification()
extract_shape_files()
add_missing_zones()
bidding_zones.to_file(snakemake.output.bidding_zones)


# In[ ]:





# In[ ]:





# In[ ]:





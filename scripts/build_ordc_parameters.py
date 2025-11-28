#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from entsoe import EntsoePandasClient
from entsoe.exceptions import InvalidBusinessParameterError, NoMatchingDataError
from requests import HTTPError
import time

def get_balancing(zone_list):
    
    unavailable_zones = []
    
    balancing = pd.DataFrame()

    for zone in zone_list:
        
        print("Extracting balancing of zone " + zone)
        
        for business_type in business_types:
            
            bal = pd.DataFrame()

            try: 
                act = client.query_activated_balancing_energy(zone, start=start, end=end, business_type=business_type)            
                act = act.resample("1h").sum().unstack()            
                act.name = business_type
                bal = pd.concat([bal, act],axis=1)

            except: # (HTTPError, NoMatchingDataError, InvalidBusinessParameterError):
                unavailable_zones.append(zone)
        
        if len(bal)>0:
            bal = bal.sum(axis=1)
            bal.name = zone 
            balancing = pd.concat([balancing, bal], axis=1)
            
    if len(unavailable_zones) >0:
        print("Unable to retrieve " + ", ".join(unavailable_zones) + " from ENTSO-E transparency" )
    
    
    balancing.reset_index(0, drop=True, inplace=True)

    balancing.loc[
        ("Down", balancing.index.levels[1]), balancing.columns
    ] = -balancing.loc[("Down", balancing.index.levels[1]), balancing.columns]

    return balancing
            


api_key = snakemake.config["entso-e"]["api_key"]
client = EntsoePandasClient(api_key=api_key)


start = pd.Timestamp(snakemake.config["balancing"]["start_extraction"], tz="Europe/Brussels")
end = pd.Timestamp(snakemake.config["balancing"]["end_extraction"], tz="Europe/Brussels")

ordc_hdf = snakemake.output.ordc_hdf

rr_hdf = snakemake.output.reserve_requirements

business_types = [ 'A96',  'A97']

reserve_requirements = pd.read_csv(snakemake.input.pemmdb,index_col=[2,1])

reserve_requirements = reserve_requirements[reserve_requirements['data_version'] == "ERAA 2024"]

#reserve_requirements.to_csv("reserve_requirements0.csv")

reserve_requirements = reserve_requirements.groupby(["MARKET_NODE","YEAR"])["Value"].sum()

#reserve_requirements.to_csv("reserve_requirements1.csv")

reserve_requirements = reserve_requirements[reserve_requirements >0]

if isinstance(reserve_requirements.index, pd.MultiIndex):
    reserve_requirements.index = reserve_requirements.index.remove_unused_levels()

reserve_requirements.to_hdf(rr_hdf, "reserve_requirements")
reserve_requirements.to_csv("reserve_requirements.csv")

bz_map = pd.Series(
    reserve_requirements.index.levels[0].str[:2],
    index = reserve_requirements.index.levels[0]
)

#ERAA now considers 5 zones in Norway
#bz_map.drop(["UKNI", "NOM1", "NOS0", "NON1"], inplace=True) # norway to be built separately because ERAA contains only 3 NO zones 
                                                            # while ENTSO-E transparency features 5 zones.

bz_map.loc[
    ["DKW1", "DKE1", "SE01", "SE02", "SE03", "SE04", "DE00", 'ITCA', 'ITCN', 'ITCS', 'ITN1', 'ITS1', 'ITSA', 'ITSI',"NOS1", "NOS2", "NOS3","NOM1", "NON1"]
] = ["DK_1","DK_2", 'SE_1', 'SE_2', 'SE_3', 'SE_4', 'DE_LU','IT_CALA', 
     'IT_CNOR', 'IT_CSUD',  'IT_NORD',  'IT_SUD' ,  'IT_SARD', 'IT_SICI','NO_1', 'NO_2', 'NO_3', 'NO_4', 'NO_5']

#norway = pd.Series(
#    ["NOS0", "NOS0", "NOM1", "NON1", "NOS0"],
#    ['NO_1', 'NO_2', 'NO_3', 'NO_4', 'NO_5'] 
#)


balancing = get_balancing(bz_map.values)

balancing_mean = balancing.stack().groupby(
    pd.Series(bz_map.index, bz_map.values)
    .reindex([i[2] for i in balancing.stack().index]).values
).mean()


balancing_std = balancing.stack().groupby(
    pd.Series(bz_map.index, bz_map.values)
    .reindex([i[2] for i in balancing.stack().index]).values
).std()


#balancing_norway = get_balancing(norway.index)

#balancing_mean = pd.concat(
#    [
#        balancing_mean, 
#        balancing_norway.stack().groupby(
#            norway.reindex([i[2] for i in balancing_norway.stack().index]).values
#        ).mean()
#    ]
#)


# In[19]:


#balancing_std = pd.concat(
#    [
#        balancing_std, 
#        balancing_norway.stack().groupby(
#            norway.reindex([i[2] for i in balancing_norway.stack().index]).values
#        ).std()
#    ]
#)



balancing_mean = balancing_mean.reindex(reserve_requirements.index.levels[0], fill_value=0)
balancing_std = balancing_std.reindex(reserve_requirements.index.levels[0])

missing_zones = balancing_std.loc[balancing_std.isna()].index

scale_missing = reserve_requirements.loc[missing_zones].groupby(level=0).mean().div(
    reserve_requirements.groupby(level=0).mean().mean()
)

balancing_std.loc[missing_zones] = scale_missing.multiply(balancing_std.mean())

ordc_parameters = pd.concat([balancing_mean, balancing_std], axis=1)
ordc_parameters.columns = ["mean", "std"]
ordc_parameters.to_hdf(ordc_hdf, "ordc_parameters")

# In[ ]:





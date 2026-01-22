#!/usr/bin/env python
# coding: utf-8

# In[43]:


import pandas as pd
import powerplantmatching as pm
import numpy as np
import pypsa
import pycountry
import os 

import time


# In[2]:


def build_inflows(inflows):

    base_year_hydro = pd.Series(
    inflows.index.levels[3].astype(int), 
    inflows.index.levels[3]).subtract(year).abs().idxmin()

    inflows = inflows.loc[:, str(climate_year),:, base_year_hydro,:]
    
    inflow_grouper = pd.Series(["hydro", "ROR", "PHS Open", "pondage"],inflows.index.levels[2])

    inflows = inflows.unstack(2).groupby(inflow_grouper, axis=1).sum().multiply(1) 
    inflows = inflows.unstack(0).stack(0)
    
    inflows.index = [" ".join(i) for i in inflows.index]
    
    inflows = inflows.T
    
    inflows.index = snapshots  
   
    #inflows.to_csv('inflows.csv', index=False)
   
    return inflows


# In[6]:


def add_existing_storage(all_c):  
 
 
    hydro = pd.read_csv(snakemake.input.hydrodata,header=0)
    
    all_c = all_c.reset_index(name="Value")
    
    closest_year = all_c.loc[(all_c['Target year'] - base_year).abs().idxmin(), 'Target year']

    all_cy = all_c[all_c['Target year'] == closest_year]

    all_cdata = (all_cy
      .pivot(index="Technology", columns="Market_Node", values=all_cy.columns[-1]) 
      )

    storage = all_cdata.loc[["PHS","PHS Open","hydro","pondage","battery"]].unstack().reset_index().copy()
    storage.columns = ["bus", "carrier", "p_nom"]
    storage.index = storage.bus + " " + storage.carrier

    hydro = hydro[hydro.apply(lambda row: row.astype(str).str.contains("ERAA 2025 pre-CfE").any(), axis=1)]

    target_year = hydro.loc[hydro["TARGET_YEAR"].sub(base_year).abs().idxmin(), "TARGET_YEAR"]
    hydro_y = hydro[hydro["TARGET_YEAR"] == target_year]
    
    hydro_pump = -(
    hydro_y.loc[
        hydro_y["PEMMDB_PLANT_TYPE"].isin(["Closed loop pumping", "Open loop pumping"])
    ]
    .pivot(index="PEMMDB_PLANT_TYPE", columns="MARKET_NODE", values="MAX PUMPING CAP (MW)")
    .fillna(0)
    .rename(index={
        "Closed loop pumping": "PHS",
        "Open loop pumping": "PHS Open"
    })
)
    hydro_pump=hydro_pump.drop(columns=["UA00"], errors='ignore').loc[:, (hydro_pump != 0).any(axis=0)]  

    hydro_cap=all_cdata.loc[all_cdata.index.isin(["PHS", "PHS Open"])]

    p_min_pu = hydro_pump.div(hydro_cap.loc[:, hydro_cap.notna().any(axis=0)].add(1e-5).values).fillna(0)
    p_min_pu.index = ["PHS", "PHS Open"]
    
    valid_countries = p_min_pu.columns[(p_min_pu != 0).any(axis=0)]
    p_min_pu_filtered = p_min_pu[valid_countries]
    p_min_pu_stacked = p_min_pu_filtered.stack().fillna(0) 
    storage_index = pd.MultiIndex.from_arrays([storage.carrier, storage.bus])
    storage["p_min_pu"] = p_min_pu_stacked.reindex(storage_index).values 

    storage.loc[storage.carrier=="battery", "p_min_pu"] = -1
    storage.loc[storage.carrier=="hydro", "p_min_pu"] = 0
    storage.loc[storage.carrier=="pondage", "p_min_pu"] = 0

    storage = storage[storage.p_nom >0]
    
    storage_inflows = inflows[storage.loc[storage.carrier=="PHS Open"].index] 
    storage_inflows = pd.concat([storage_inflows,inflows[storage.loc[storage.carrier == "hydro"].index]],axis=1)
    storage_inflows = pd.concat([storage_inflows,inflows[storage.loc[storage.carrier == "pondage"].index]],axis=1)
        
    hydro_stored = 1000000*(hydro_y.loc[hydro_y["PEMMDB_PLANT_TYPE"].isin(["Closed loop pumping", "Open loop pumping", "Reservoir", "Pondage"])]
      .pivot(index="PEMMDB_PLANT_TYPE", columns="MARKET_NODE", values="Storage Capacity [TWh]")
      .fillna(0)
      .rename(index={
          "Closed loop pumping": "PHS",
          "Open loop pumping": "PHS Open",
          "Reservoir": "hydro",
          "Pondage": "pondage"                    
      })
)       

    
    #storage_inflows.index.to_series().to_csv('storage_inflows_index.csv', index=True)
        
    battery = pd.read_csv(snakemake.input.battery,header=0)
    battery = battery[battery.apply(lambda row: row.astype(str).str.contains("ERAA 2025 pre-CfE").any(), axis=1)]
    #battery.to_csv('battery.csv', index=False)
    
    battery_y = battery[battery["TARGET_YEAR"] == base_year]
    
    #battery_y.to_csv('battery_y.csv', index=False)
    
    battery_by_node = battery_y.groupby("MARKET_NODE")["STORAGE CAPACITY (MWh)"].sum().reset_index()
    
    #battery_by_node.to_csv('battery_y.csv', index=False)
    battery_by_node["Technology"] = "battery"
    battery_pivot = battery_by_node.pivot(index="Technology", columns="MARKET_NODE", values="STORAGE CAPACITY (MWh)").fillna(0)
    
    #battery_pivot.to_csv('battery_pivot.csv', index=False)
    
    storage_capacity=pd.concat([hydro_stored, battery_pivot])
    
    #storage_capacity.to_csv('storage_capacity.csv', index=False)
      
    valid_countries = storage_capacity.columns[(storage_capacity != 0).any(axis=0)]
    stored_filtered = storage_capacity[valid_countries]
        
    stored_stacked = stored_filtered.stack() 
    
    storage_index = pd.MultiIndex.from_arrays([storage.carrier, storage.bus])

    storage["max_hours"] = stored_stacked.reindex(storage_index).clip(lower=5000).div(storage.p_nom.values).values 
       
    map_storage = pd.Series( ["battery inverter", "PHS", "hydro","hydro","hydro"], ["battery", "PHS", "PHS Open","hydro","pondage"])

    storage["efficiency_store"] = technology_data.loc[:, "efficiency", :].reindex(storage.carrier.map(map_storage)).value.values
    
    #storage.to_csv('storage.csv', index=False)

    n.add(
        "StorageUnit",
        storage.index,
        **storage,
        inflow=storage_inflows.reindex(storage.index, axis=1, fill_value=0.),
        invest_status = "existing"
    )


def prepare_commodity_prices(commodity_prices):
    
    fuels = ["co2", "gas", "coal", "oil","lignite", "hydrogen", "nuclear"]
    C02= [57,0, 94, 0,101, 0, 85.3] #gas, co2, coal, hydrogen, lignite, nuclear, oil
        
    def assign_fuel_group(fuel_name):
      fuel_name_lower = str(fuel_name).lower()
      for f in fuels:
          if f.lower() in fuel_name_lower:
              return "CCGT" if f.lower() == "gas" else f
      return None
    
    commodity_prices["Fuel"] = commodity_prices["Fuel"].apply(assign_fuel_group)  
    
    commodity_prices_filtered = commodity_prices[commodity_prices["Fuel"].notna()]
    
    commodity_prices["Fuel"] = commodity_prices["Fuel"].replace("gas", "CCGT")
    
    commodity_prices = (
    commodity_prices_filtered
    .groupby(["Year", "Fuel"], as_index=False)["Value"]
    .mean()
    )

     
    commodity_prices = commodity_prices.pivot(index="Fuel", columns="Year", values="Value").reset_index()
     
    commodity_prices["CO2"] = C02
    
    commodity_prices = commodity_prices.set_index("Fuel")
    
    commodity_prices = commodity_prices.apply(pd.to_numeric, errors="coerce")
      
    converter = pd.DataFrame(3.6 , commodity_prices.index, commodity_prices.columns)
    converter["CO2"] = 3.6/1000 
    converter.loc["co2"] = 1

    commodity_prices = converter*commodity_prices


    to_add = commodity_prices.reindex(["CCGT", "CCGT"])
    to_add.index = ["OCGT", "biomass"]
    to_add.loc[["biomass"]] = biomass_price
    to_add.loc[["biomass"], "CO2"] = 0
    commodity_prices = pd.concat([commodity_prices, to_add])
    #commodity_prices.to_csv('commodity_prices.csv', index=False)
    
    return commodity_prices

def add_renewables():
       
    res = capacity.loc[["CSP","CSP-stor","onwind","offwind","solar-fix","solar-ind", "solar-rsd", "solar-track","ROR"]].unstack().copy()
    res = res.reset_index()
    res.columns = ["bus", "carrier", "p_nom"]
    res.index = res.bus + " " + res.carrier

    res = res[res.p_nom>0]
    
    vre=["CSP","CSP-stor","onwind","offwind","solar-fix","solar-ind", "solar-rsd", "solar-track"]
   
    p_max_pu = pd.concat(
        [pd.read_hdf(snakemake.input.res_profile, tech).loc[:, "WS{:02}".format(climate_year), str(base_year), :] for tech in vre],
        axis=1
    )
    
    p_max_pu.columns = vre
    p_max_pu = p_max_pu.unstack(0).stack(0)
    p_max_pu.index = [" ".join(i) for i in p_max_pu.index]
    p_max_pu.columns = snapshots
    p_max_pu = p_max_pu.T
    p_max_pu = pd.concat(
        [p_max_pu, inflows[res.filter(like="ROR", axis=0).index].div(res.filter(like="ROR", axis=0).p_nom).clip(upper=1)],
        axis=1
    )

    n.add(
        "Generator",
        res.index,
        **res,
        p_max_pu=p_max_pu[res.index], 
        invest_status="policy"
    )

def add_dispatchables():
         
    plants = dispatchable_plants.query("(entry <= @year) & exit > @year")

    plants["marginal_cost"] = (
        commodity_prices
        .reindex(plants.carrier)[base_year].div(
            plants.efficiency.values,
        ).add(
            (
                commodity_prices.reindex(
                    plants.carrier)["CO2"]
                .multiply(
                    commodity_prices.loc["co2", base_year]
                )
            )
        ).add(plants.var_om.values).values
    )
    
    plants.loc[:, "committable"]= True
    
    """
    p_max_puA = pd.read_hdf(snakemake.input.maintenance, key=f"maintenance{base_year}") 
    
      
    p_max_puA.index=snapshots
    p_max_pu_columns_aligned = pd.DataFrame(0, index=snapshots, columns=plants.index)
    
    p_min_pu=pd.DataFrame(np.tile(plants["p_min_pu"].values, (len(snapshots), 1)), index=snapshots, columns=plants.index)
    
    common_cols = plants.index.intersection(p_max_puA.columns)
    
    p_max_pu_columns_aligned.loc[:, common_cols] = p_max_puA[common_cols].astype(int)

    p_max_pu =1-p_max_pu_columns_aligned
    
    p_min_pu=np.minimum(p_min_pu, p_max_pu)
    
    #p_min_pu.to_csv('p_min_pu1.csv', index=True)

    plants.drop(columns="p_min_pu", inplace=True)
    """
    
    n.add(
        "Generator",
        plants.index,
        **plants,
        #p_max_pu=p_max_pu[plants.index], 
        #p_min_pu=p_min_pu[plants.index], 
    )

def group_luxembourg(demand, links):
    
    demand_grouper = pd.Series(demand.columns, demand.columns)
    demand_grouper.loc[demand_grouper.index.str[:2] == "LU"] = "LUG1"
    demand = demand.groupby(demand_grouper, axis=1).sum()

    links.loc[links.bus0.str[:2] == "LU", "bus0"] = "LUG1"
    links.loc[links.bus1.str[:2] == "LU", "bus1"] = "LUG1"
    
    return demand, links

def add_dsr():
    
    
    dsr=pd.read_hdf(snakemake.input.dsr)
    
    dsr = dsr.loc[base_year, :]
    
    dsr.unstack(1).columns = ["bus", "p_nom", "marginal_cost"]
    
    dsr=dsr.drop(columns=["hours","TARGET_YEAR"])
    
    dsr["carrier"] = "DSR"
    
    dsr.index = [f"{bus} {carrier} {i}" for i, (bus, carrier) in enumerate(zip(dsr.bus, dsr.carrier))]
    
    print(dsr.index)

    n.add(
        "Generator",
        dsr.index,
        **dsr,
        invest_status = "existing"
    )
    

def set_investment_bounds():
    
    n.generators.loc[n.generators.invest_status == "policy", "p_nom_max"] = n.generators.loc[n.generators.invest_status == "policy", "p_nom"]
    n.generators.loc[n.generators.invest_status == "policy", "p_nom_min"] = n.generators.loc[n.generators.invest_status == "policy", "p_nom"]
    
def links(PTDF_core, PTDF_nordic):
    
    DC_core=PTDF_core["PTDF_EvFB"].columns
    AHC_core=PTDF_core["PTDF*_AHC,SZ"].columns
    DC_nordic=PTDF_nordic["PTDF_EvFB"].columns
    AHC_nordic=PTDF_nordic["PTDF*_AHC,SZ"].columns
    
    
    links = pd.read_hdf("resources/ntcs.h5", "p_nom")
    links = links[base_year].unstack(1)
    links = links.iloc[:len(snapshots)]
    links.dropna(inplace=True)
    links["carrier"] = [i[-2:] for i in links.index]
    #links.to_csv('Linkp_nom0.csv', index=True)
    links_df=links.copy()
    
    links_df["sorted_link"] = links[["bus0", "bus1"]].apply(lambda x: "-".join(sorted([x[0], x[1]])), axis=1)
    links_df = links_df.groupby(["sorted_link", "carrier"]).agg(
        forward_pnom=("p_nom", "first"),
        reverse_pnom=("p_nom", "last")
    ).reset_index()
    
    #links_df.to_csv('LinkResult.csv', index=True)
    
    links_p_max_pu = pd.read_hdf("resources/ntcs.h5", "p_max_pu")  
    links_p_max_pu = links_p_max_pu[base_year].unstack(2)
    links_p_max_pu = links_p_max_pu.dropna(how='all', axis=0)
    links_p_max_pu = links_p_max_pu.iloc[:len(snapshots)]  
    links_p_max_pu.index = snapshots
    links_p_max_pu = links_p_max_pu.dropna(how='all',axis=1)

    links_p_min_pu=links_p_max_pu.copy()

    
    for index, row in links_df.iterrows():
        bus0, bus1 =row["sorted_link"].split("-")
        carrier = row["carrier"]
        reverse_link = f"{bus1}-{bus0}-{carrier}"
        direct_link = f"{bus0}-{bus1}-{carrier}"
        aux2=f"{bus1}-{bus0}"

        if [j for j in DC_core if aux2 ==j]or[j for j in DC_nordic if aux2 ==j]or[j for j in AHC_core if aux2 ==j]or[j for j in AHC_nordic if aux2 ==j]:
            links_p_min_pu[reverse_link]=-links_p_min_pu[direct_link]*row["forward_pnom"]/row["reverse_pnom"]
            links_p_min_pu = links_p_min_pu.drop(columns=direct_link)
            links_p_max_pu = links_p_max_pu.drop(columns=direct_link)
            links=links.drop(index=direct_link)
        else:
            if direct_link in links_p_min_pu.columns and reverse_link in links_p_min_pu.columns:
                links_p_min_pu[direct_link]=-links_p_min_pu[reverse_link]*row["reverse_pnom"]/row["forward_pnom"]
                links_p_min_pu = links_p_min_pu.drop(columns=reverse_link)
                links_p_max_pu = links_p_max_pu.drop(columns=reverse_link)
                links=links.drop(index=reverse_link)
            elif direct_link in links_p_min_pu.columns:
                links_p_min_pu[direct_link]=0
            elif reverse_link in links_p_min_pu.columns:
                links_p_min_pu[reverse_link]=0
    
    links_p_max_pu = links_p_max_pu.drop(columns="UK00-FR00_1-DC")
    links_p_max_pu = links_p_max_pu.drop(columns="UK00-FR00_2-DC")
    links_p_min_pu = links_p_min_pu.drop(columns="UK00-FR00_1-DC")
    links_p_min_pu = links_p_min_pu.drop(columns="UK00-FR00_2-DC")
    links_p_min_pu["FR00-UK00_1-DC"] = -links_p_min_pu["FR00-UK00_1-DC"]
    links_p_min_pu["FR00-UK00_2-DC"] = -links_p_min_pu["FR00-UK00_2-DC"]
    links=links.drop(index="UK00-FR00_1-DC")
    links=links.drop(index="UK00-FR00_2-DC")
    
    #links.to_csv('Linkp_nom.csv', index=True)
#    links_p_max_pu.to_csv('links_p_max_pu1.csv', index=True)
#    links_p_min_pu.to_csv('links_p_min_pu.csv', index=True)
    
    return links_p_max_pu, links_p_min_pu, links    
commodity_prices_raw = pd.read_csv(snakemake.input.commodity_prices,header=0)

all_cap=pd.read_hdf(snakemake.input.all_capacities)

year = int(snakemake.params.ty)
climate_year = int(snakemake.params.cy)

years = snakemake.params.years

save_path = snakemake.output.network

biomass_price = snakemake.config["biomass_price"]

dispatchable = ["CCGT", "OCGT", 'biomass', 'coal', 'lignite', 'nuclear','oil','hydrogen']
distributed_resources = ["onwind", "offwind", "CSP","f", "battery"]

technologies_for_investment = snakemake.config["power_plants"]["technologies_new_investment"]

inflows_raw = pd.read_hdf(snakemake.input.inflow, "inflow")

snapshots = pd.date_range(start="2010-01-01", freq="h", periods=8760)

base_year = pd.Series(
all_cap.index.levels[0].astype(int), 
all_cap.index.levels[0]).subtract(year).abs().idxmin()

technology_data = pd.read_csv(snakemake.input.technology_data, index_col=[0,1])
demand = pd.read_hdf(snakemake.input.demand)

demand = demand.loc[:, "WS{:02}".format(climate_year), :, str(base_year), :].unstack(1)
demand.index = snapshots

PTDF_core = pd.read_hdf(snakemake.input.core_domain, "PTDF")
PTDF_nordic = pd.read_hdf(snakemake.input.nordic_domain, "PTDF")

links_p_max_pu,links_p_min_pu, links=links(PTDF_core, PTDF_nordic)

demand, links = group_luxembourg(demand, links)

commodity_prices = prepare_commodity_prices(commodity_prices_raw[commodity_prices_raw.apply(lambda row: row.astype(str).str.contains("ERAA 2025 post-CfE").any(), axis=1)])

if snakemake.config["power_plants"]["aggregation_level"] == "none":
  dispatchable_plants = pd.read_hdf(snakemake.input.individual_plants, key='detailed')
elif snakemake.config["power_plants"]["aggregation_level"] == "small":
  dispatchable_plants = pd.read_hdf(snakemake.input.individual_plants, key='small_aggregated')
else: 
  dispatchable_plants = pd.read_hdf(snakemake.input.dispatchable_plants)

n = pypsa.Network()
n.set_snapshots(snapshots)

capacity=pd.read_hdf(snakemake.input.investcap)[base_year].unstack(1)

buses = (
    capacity.sum()[capacity.abs().sum() >0].index
    .union(demand.columns)
    .union(links.bus0.unique())
    .union(links.bus1.unique())
)

n.add(
    "Bus", 
    buses, 
    carrier = "electricity", 
    country = buses.str[:2]
)

n.add(
    "Load", 
    demand.columns,
    bus=demand.columns,
    p_set = demand
)

inflows = build_inflows(inflows_raw)

add_existing_storage(all_cap)

add_dispatchables()

add_renewables()

add_dsr()

set_investment_bounds()


links_p_max_pu = links_p_max_pu.reindex(
      index=n.snapshots,   
      columns=links.index, 
      fill_value=1
)

links_p_min_pu = links_p_min_pu.reindex(
      index=n.snapshots,   
      columns=links.index, 
      fill_value=1, 
)

n.add(
    "Link",
    links.index,
    bus0 = links.bus0,
    bus1 = links.bus1,
    p_nom = links.p_nom,
    p_max_pu = links_p_max_pu.reindex(links.index, fill_value=1, axis=1),
    p_min_pu = links_p_min_pu.reindex(links.index, fill_value=1, axis=1),
    carrier = links.carrier,
)

#n.generators.loc[(n.generators.bus == "CY00") & (n.generators.carrier == "CCGT"), "p_min_pu"] = 0 # remove minimum load of CCGT in Cyprus as this can exceed actual load.

dirname = os.path.dirname(save_path)  

if not os.path.exists(dirname):
    os.makedirs(dirname)


n.export_to_netcdf(save_path)

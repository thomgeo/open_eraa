#!/usr/bin/env python
# coding: utf-8

# In[1]:


import powerplantmatching as pm
import linopy
import numpy as np
import pandas as pd
import geopandas as gpd
import os


# In[2]:


def prepare_power_plant_table():

    res = ['Wind', 'Solar', "Hydro"]
    
    match_carriers = {
        'Nuclear': "nuclear",
        'Hard Coal': "coal",
        'Lignite': "lignite",
        'Natural Gas': np.nan,
        'Oil': "oil",
        'Solid Biomass': "biomass",
        'Other': "other",
        'Biogas': "biomass",
        'Waste': "biomass",
        'Geothermal': "geothermal"
    }
    
    power_plants = pm.powerplants()
    power_plants = power_plants.query("(Fueltype not in @res)")
    power_plants = power_plants[~(power_plants.DateOut < 2005)]
    
    power_plant_locations = gpd.GeoDataFrame(index=power_plants.index, geometry=gpd.points_from_xy(power_plants.lon, power_plants.lat), crs="EPSG:4326")
    power_plant_locations = power_plant_locations.to_crs("EPSG:3035")
    
    plant_to_zone = gpd.sjoin_nearest(
        power_plant_locations,
        bidding_zones,
        how="left"
    )
    
    plant_to_zone = plant_to_zone.groupby(plant_to_zone.index).last()
    power_plants["bus"] = power_plants.index.map(plant_to_zone.zone)

    power_plants["carrier"] = power_plants.Fueltype.map(match_carriers)

    for plant in power_plants.index[power_plants.carrier.isna()]:
        if power_plants.loc[plant, "Technology"] == "CCGT":
            power_plants.loc[plant, "carrier"] = "CCGT"
        else:
            power_plants.loc[plant, "carrier"] = "OCGT"
    
    power_plants.loc[(power_plants.bus=="ES00") & (power_plants.Set == "CHP"), "carrier"] = "other"
    
    return power_plants



def get_age_matching():
    
    age_matching = pd.DataFrame(index=technology_parameters.index, dtype=str)

    age_matching["from"] = np.nan
    age_matching["to"] = np.nan
    
    idx = pd.IndexSlice
    
    age_matching.loc[idx[:, "old 1"], "from"] = 1900
    age_matching.loc[idx[:, "old 1"], "to"] = 1980
    age_matching.loc[idx[:, "old 2"], "from"] = 1980
    age_matching.loc[idx[["coal", "lignite", "oil", "OCGT"], "old 2"], "to"] = 2010
    age_matching.loc[idx[["CCGT"], "old 2"], "to"] = 2000
    age_matching.loc[idx[["coal", "lignite", "oil", "OCGT"], "new"], "from"] = 2010
    age_matching.loc[idx[:, "new"], "to"] = 2040
    age_matching.loc[idx[:, "present 1"], "from"] = 2000
    age_matching.loc[idx[:, "present 1"], "to"] = 2010
    age_matching.loc[idx[:, "present 2"], "from"] = 2010
    age_matching.loc[idx[:, "present 2"], "to"] = 2020
    age_matching.loc[idx[["CCGT"], "new"], "from"] = 2020
    age_matching.drop(["old"], level=1, inplace=True)
    age_matching.drop("CCS", level=1, inplace=True)
    age_matching.drop(('oil','-'), inplace=True)
    age_matching.loc[["nuclear", "hydrogen", "biomass", "other"], :] = 1900, 2040

    return age_matching
    
def add_technical_properties_existing(eraa_plants_detailed, detailed_plants, properties, age_matching):

    idx = pd.IndexSlice
    
    eraa_plants_detailed["age_class"] = np.nan
    
    for plant in detailed_plants.index:
    
        datein = power_plants.loc[plant, "DateIn"]
        carrier = power_plants.loc[plant, "carrier"]
        
        if (carrier in age_matching.index.levels[0]) and not np.isnan(datein):
            eraa_plants_detailed.loc[plant, "age_class"] = age_matching[(age_matching["from"] <= datein)&(age_matching["to"] > datein)].loc[carrier].index[0]
        elif (carrier in age_matching.index.levels[0]) and ("old 2" in age_matching.loc[carrier].index):
            eraa_plants_detailed.loc[plant, "age_class"] = "old 2"    
    
    eraa_plants_detailed.age_class = eraa_plants_detailed.age_class.fillna("new")
    all_plants.loc[all_plants.carrier.isin(["biomass", "nuclear", "other"]), "age_class"] = "-"
    eraa_plants_detailed.loc[eraa_plants_detailed.carrier =="hydrogen", "age_class"] = "CCGT"
    
    eraa_plants_detailed['invest_status'] = "existing"
    eraa_plants_detailed['p_nom_max'] = eraa_plants_detailed.p_nom
    eraa_plants_detailed['p_nom_min'] = 0.01
    
    eraa_plants_detailed[['start_up_cost', 'ramp_limit_up', 'ramp_limit_down', 'p_min_pu', 'min_up_time', 'min_down_time', 'var_om' ]] = np.nan
    
    start_up_cost = []
    ramp_limit_up = []
    ramp_limit_down = []
    p_min_pu = []
    min_up_time = []
    min_down_time = [] 
    var_om = [] 
    efficiency = []
    
    for plant in eraa_plants_detailed.index:
    
        carrier = eraa_plants_detailed.loc[plant, "carrier"]
        age_class = eraa_plants_detailed.loc[plant, "age_class"]
        
        start_up_cost.append(properties.loc[idx[carrier, age_class], "start_up_fix_cost"])
        ramp_limit_up.append(properties.loc[idx[carrier, age_class], "ramp_limit_up"])
        ramp_limit_down.append(properties.loc[idx[carrier, age_class], "ramp_limit_down"])
        p_min_pu.append(properties.loc[idx[carrier, age_class], "p_min_pu"])
        min_up_time.append(properties.loc[idx[carrier, age_class], "min_up_time"])
        min_down_time.append(properties.loc[idx[carrier, age_class], "min_down_time"])
        var_om.append(properties.loc[idx[carrier, age_class], "var_OM"])
        efficiency.append(properties.loc[idx[carrier, age_class], "efficiency"])
    
    eraa_plants_detailed.start_up_cost = start_up_cost
    
    """
    eraa_plants_detailed.ramp_limit_down_real = ramp_limit_down
    eraa_plants_detailed.ramp_limit_down = eraa_plants_detailed.ramp_limit_down_real.clip(upper=1.) 
    eraa_plants_detailed.ramp_limit_up_real = ramp_limit_up
    eraa_plants_detailed.ramp_limit_up = eraa_plants_detailed.ramp_limit_up_real.clip(upper=1.)
    """
    
    eraa_plants_detailed.ramp_limit_down = ramp_limit_down
    eraa_plants_detailed.ramp_limit_up = ramp_limit_up
    eraa_plants_detailed.p_min_pu = p_min_pu
    eraa_plants_detailed.min_up_time = min_up_time
    eraa_plants_detailed.min_down_time = min_down_time
    eraa_plants_detailed.var_om = var_om
    eraa_plants_detailed.efficiency = efficiency

    eraa_plants_detailed.p_nom = eraa_plants_detailed.p_nom.astype(float)
    eraa_plants_detailed.p_nom_max = eraa_plants_detailed.p_nom_max.astype(float)

    return eraa_plants_detailed
    
def build_new_investments(df, properties):

    technologies_for_investment = ["OCGT", "CCGT"]
    years = df.entry.unique()
    zones_for_investment = df.bus.unique()

    new = pd.DataFrame(
        0.01, 
        index= pd.MultiIndex.from_product(
            [zones_for_investment,technologies_for_investment, years], 
            names=["bus", "carrier", "entry"]), 
        columns=["p_nom"]
    )

    new.reset_index(inplace=True)
    new.index = new.bus + " " + new.carrier + " new " + new.entry.astype(str)
    new["efficiency"] = properties.loc[:, "new", :]["efficiency"].reindex(new.carrier).values
    new["start_up_cost"] = properties["start_up_fix_cost"].loc[:, "new"].reindex(new.carrier, level=1).fillna(0).values
    new["ramp_limit_up"] = properties["ramp_limit_up"].loc[:, "new"].reindex(new.carrier).values*60
    new["ramp_limit_down"] = properties["ramp_limit_down"].loc[:, "new"].reindex(new.carrier).values*60
    new["p_min_pu"] = properties["p_min_pu"].loc[:, "new"].reindex(new.carrier).values
    new["min_up_time"] = properties["min_up_time"].loc[:, "new"].reindex(new.carrier).values
    new["min_down_time"] = properties["min_down_time"].loc[:, "new"].reindex(new.carrier).values
    new["exit"] = target_years[-1] + 1
    new["p_nom_max"] = np.inf
    new["p_nom_min"] = 0.01
    new["var_om"] = properties.loc[:, "new",:].reindex(new.carrier)["var_OM"].values
    new["invest_status"] = "new"
    new["age_class"] = "new"
    
    return pd.concat([df, new])
    
def aggregate_small(eraa_plants_detailed):

    small_plants = eraa_plants_detailed[eraa_plants_detailed.p_nom.round()<=de_minimis]
    
    agg_func = pd.Series(
        "mean",
        index = ['p_nom', 'efficiency',  'p_nom_max', 'p_nom_min', 'start_up_cost', 'ramp_limit_up',
                 'ramp_limit_down', 'p_min_pu', 'min_up_time', 'min_down_time', 'var_om', "age_class"]
    )
    
    agg_func.loc[["p_nom", "p_nom_max"]] = "sum"
    agg_func.loc["age_class"] = pd.Series.mode
    
    aggregated_small = small_plants.groupby(["bus", "carrier", "entry", "exit", 'invest_status', "cluster"]).agg(agg_func).reset_index()
    small_index = pd.Index(aggregated_small.bus + " small " + aggregated_small.carrier + " ")
    aggregated_small.index = small_index + pd.Series(small_index).to_frame().groupby(0).cumcount().astype(str)
    
    eraa_plants = eraa_plants_detailed.drop(small_plants.index)
    return pd.concat([eraa_plants, aggregated_small])

def match_to_cluster(bus, carrier):
    
    capacity = power_plants.query("(carrier == @carrier) and (bus == @bus)").Capacity.copy()
    clusters = aggregated_plants.query("bus == @bus and carrier == @carrier").copy()
    clusters = clusters.filter(like="exit", axis=0)
    
    if (len(capacity)>0) and (len(clusters)>0):
        
        m = linopy.Model()    
        clusters.index.name = "cluster"      
        cluster_matching = m.add_variables(binary=True, coords=[clusters.index, capacity.index], name="cluster_matching")       
        matched_capacity = m.add_variables(lower=0, coords=[clusters.index], name="matched_capacity")
        
        m.add_constraints(
            matched_capacity <= clusters.p_nom,
            name="upper_limit_is_total"
        )
                
        m.add_constraints(
            (cluster_matching*capacity).sum("id") == matched_capacity,
            name="match_plants_to_cluster"
        )
        
        m.add_constraints(
            cluster_matching.sum("cluster") <= 1,
            name="unique_assignment"
        )
        
        m.add_objective((matched_capacity**2).sum() - (2*clusters.p_nom*matched_capacity).sum())       
        m.solve(solver_name="cplex", **solver_options)   
        
        return m.solution["cluster_matching"].to_series(), m.termination_condition
    
    else:
        return pd.Series([]), "optimal"


def match_existing_plants():
    
    matched_plants=[]
    status = pd.Series(index= pd.unique(pd.MultiIndex.from_frame(aggregated_plants[["bus", "carrier"]])), dtype=str)
    for i in status.index:
        bus, carrier = i
        print(bus, carrier)
        matched, status[i] = match_to_cluster(bus, carrier)
        matched_plants.append(matched)
    
    if (status != "optimal").sum() == 0:
        print("Large decomposition finished successfully")
    else:
        raise ValueError('Not all decompositions converged')
    
    matched_plants = pd.concat(matched_plants)
    matched_plants.index = pd.MultiIndex.from_tuples(matched_plants.index)
    matched_plants = matched_plants.round(0).astype(int)
    
    cluster_association = matched_plants[matched_plants == 1].reset_index(0).level_0
    cluster_association.name="cluster"
    
    detailed_plants = power_plants[["bus", "carrier", "Capacity", "Efficiency"]].copy()
    detailed_plants.columns = detailed_plants.columns.str.replace("Efficiency", "efficiency").str.replace("Capacity", "p_nom")
    detailed_plants = detailed_plants.reindex(cluster_association.index)
    detailed_plants["cluster"] = cluster_association
    detailed_plants.sort_index(inplace=True)
    detailed_plants["invest_status"] = "existing"
    detailed_plants["entry"] = target_years[0]
    detailed_plants["exit"] = [int(i[1]) for i in cluster_association.str.split("exit ")]
    
    return detailed_plants, matched_plants


def build_generic_plants():
    
    unmatched_capacity = aggregated_plants.query("invest_status == 'existing'").p_nom.subtract(
        matched_plants.unstack(0).multiply(power_plants.loc[matched_plants.index.levels[1], "Capacity"], axis=0).sum(),
        fill_value=0
    )
    
    typical_size_unmatched = typical_size.reindex(aggregated_plants.reindex(unmatched_capacity.index).carrier.values)
    typical_size_unmatched.index = unmatched_capacity.index
    n_generic = (unmatched_capacity//typical_size_unmatched).astype(int)
    residual = unmatched_capacity.subtract(n_generic.multiply(typical_size_unmatched))
    
    generic_plants = []
    
    for cluster in n_generic[(n_generic>0)|(residual>0)].index:    
        generic_unit = aggregated_plants.loc[cluster][["bus", "carrier", "efficiency","p_nom", "entry", "exit", "invest_status"]].copy()
        generic_unit.p_nom = float(typical_size[generic_unit.carrier])
        generic_unit["cluster"] = cluster
        generic_unit["efficiency"] = np.nan
    
        add_on = " ".join(cluster.split(" ")[-2:])
        for number in range(n_generic[cluster]):
            unit_to_add = generic_unit.copy()
            
            unit_to_add.name = "{bus} {carrier} {add_on} {number}".format(
                bus=generic_unit["bus"], carrier=generic_unit["carrier"], number=number, add_on=add_on
            )
            
            generic_plants.append(unit_to_add)
    
        residual_unit = generic_unit.copy()
        
        residual_unit.name = "{bus} {carrier} {add_on} {number}".format(
            bus=generic_unit["bus"], carrier=generic_unit["carrier"], number=n_generic[cluster], add_on=add_on
        )
        
        residual_unit.p_nom = residual[cluster]
        generic_plants.append(residual_unit)
    
    generic_plants = pd.concat(generic_plants, axis=1).T

    return generic_plants



technology_parameters = pd.read_hdf(snakemake.input.technology_parameters)

typical_size = pd.Series(snakemake.params.typical_size)
aggregated_plants = pd.read_hdf(snakemake.input.dispatchable_capacities)
bidding_zones = gpd.read_file(snakemake.input.bidding_zones)
bidding_zones = bidding_zones.to_crs("EPSG:3035")
de_minimis = float(snakemake.params.de_minimis)

solver = snakemake.config["solving"]["solver_decomposition"]["name"]
options = snakemake.config["solving"]["solver_decomposition"]["options"]
solver_options = snakemake.config["solving"]["solver_options"][options]

power_plants = prepare_power_plant_table()
target_years = sorted(aggregated_plants.entry.unique())

detailed_plants, matched_plants = match_existing_plants()
generic_plants = build_generic_plants()
all_plants = pd.concat([detailed_plants, generic_plants])

age_matching = get_age_matching()
all_plants = add_technical_properties_existing(all_plants, detailed_plants, technology_parameters, age_matching)
small_aggregated = aggregate_small(all_plants)
all_plants = build_new_investments(all_plants, technology_parameters)
small_aggregated = build_new_investments(small_aggregated, technology_parameters)

hdf_file = snakemake.output.individual_power_plants
os.makedirs(os.path.dirname(hdf_file), exist_ok=True)
all_plants.to_hdf(hdf_file, key="detailed")
small_aggregated.to_hdf(hdf_file, key="small_aggregated")


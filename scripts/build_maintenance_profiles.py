#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import linopy

def prepare_input_data(bus, year):
    
    availability_factors = pd.concat([pd.read_hdf("resources/res_profile.h5", carrier) for carrier in res_carriers], keys=res_carriers)
    res_caps = all_caps.loc[year, res_carriers,bus, :].reset_index([0,2], drop=True)
    
    country_plants = plants[(plants.entry<=year)&(plants.exit>year)&(plants.bus==bus)]
    total_res_generation = availability_factors.loc[:, :, :, str(year), bus].unstack(0).multiply(res_caps).sum(axis=1)
    
    country_plants.index.name = "plant"
        
    max_demand = demand.loc[:, :, bus, str(year)].groupby(level=0).max()
    min_res = total_res_generation.groupby(level=0).min()
    max_residual_demand_per_week = (max_demand - min_res).groupby(max_demand.index//(24*7)).max()
    max_residual_demand_per_week.index = max_residual_demand_per_week.index*7
    max_residual_demand_per_week = max_residual_demand_per_week.reindex(day).ffill()
    
    outage_days_country_plants = number_of_days.reindex(pd.MultiIndex.from_arrays([country_plants.carrier, country_plants.age_class]))
    outage_days_country_plants.index = country_plants.index
    share_winter_country_plants = share_winter.reindex(pd.MultiIndex.from_arrays([country_plants.carrier, country_plants.age_class]))
    share_winter_country_plants.index = country_plants.index

    return country_plants, max_residual_demand_per_week, outage_days_country_plants, share_winter_country_plants, 


def optimize_maintenance_scheduling(max_residual_demand_per_week):

    country_plants_p_nom = country_plants.p_nom.div(1e3) 
    max_residual_demand_per_week_normed = max_residual_demand_per_week.div(1e3)
    
    m = linopy.Model()

    plant = country_plants.index
    
    maintenance = m.add_variables(binary=True, coords=[day, plant], name="maintenance")
    
    daily_margin = m.add_variables(coords=[day], name="daily_margin")
    
    days_winter = share_winter_country_plants.multiply(outage_days_country_plants).sort_values().round()
    
    m.add_constraints(
        maintenance.loc[winter_time_days].sum("day")>=days_winter,
        name="winter_maintenance_days"
    )
    
    m.add_constraints(
        maintenance.sum("day") == outage_days_country_plants,
        name="maintenance_target"
    )

    m.add_constraints(
        maintenance.sum("plant").loc[364] == 0,
        name="no_maintenance_silvester"
    )
    
    m.add_constraints(
        daily_margin == ((1 - maintenance)*country_plants_p_nom).sum("plant") - max_residual_demand_per_week_normed + constant,
        name="daily_margin"
    )
    
    m.add_objective(daily_margin**2)

    m.solve(solver_name=solver_name, **solver_options)

    return m

technology_parameters = pd.read_hdf(snakemake.input.technology_parameters)

plants = pd.read_hdf(snakemake.input.power_plants, "detailed")
demand = pd.read_hdf(snakemake.input.demand)

solver_name = snakemake.config["solving"]["solver_maintenance"]["name"]
options = snakemake.config["solving"]["solver_maintenance"]["options"]
solver_options = snakemake.config["solving"]["solver_options"][options]

res_carriers = [ "onwind",  "offwind", "solar-track", "solar-fix", "solar-rsd",  "solar-ind",  "CSP", "CSP-stor"]

save_hdf = snakemake.output.maintenance_profiles
all_caps = pd.read_hdf(snakemake.input.all_caps)

winter_time_days = list(range(59)) + list(range(365-31, 365))
constant = 1000
day = pd.Index(range(365), name="day")
target_years = all_caps.index.levels[0]

number_of_days = technology_parameters.maintenance_n_days.copy()
share_winter = technology_parameters.maintenance_share_winter.copy()

maintenance_profiles = []
status = []
key=[]
for bus in plants.bus.unique():
    for year in target_years:
        print(bus, year)
        country_plants, max_residual_demand_per_week, outage_days_country_plants, share_winter_country_plants = prepare_input_data(bus, year)
        m =  optimize_maintenance_scheduling(max_residual_demand_per_week)
        maintenance_profiles.append(m.solution["maintenance"].to_series())
        status.append(m.termination_condition)
        key.append((bus, year))


index = pd.MultiIndex.from_product([plants.bus.unique(), target_years])
maintenance_profiles = pd.concat(maintenance_profiles, keys=index)
maintenance_profiles = maintenance_profiles.reset_index(0, drop=True)



for year in target_years:
        
    maintenance_profiles_year = maintenance_profiles.loc[year].unstack(1)
    maintenance_profiles_year.index = maintenance_profiles_year.index * 24
    maintenance_profiles_year = maintenance_profiles_year.reindex(range(8760)).ffill()

    maintenance_profiles_year = 1 - maintenance_profiles_year.astype(int)
    
    maintenance_profiles_year.to_hdf(save_hdf, key="maintenance{}".format(year))


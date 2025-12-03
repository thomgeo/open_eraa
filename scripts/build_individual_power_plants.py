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


# In[3]:


def run_disaggregation(bus, carrier):

    weight_float = 4
    
    existing_large = power_plants.query("(carrier == @carrier) and (bus == @bus)").Capacity.copy()
    weights_large = power_plants.loc[existing_large.index, ("DateIn", "DateRetrofit")].max(axis=1).subtract(2025).div(100).add(1.5).fillna(1)
    
    m = linopy.Model()
    
    years = caps_per_year.columns
    years.name ="year"
    
    status_existing_large = m.add_variables(binary=True, coords=[years, existing_large.index], name="status_existing_large")

    added_float = m.add_variables(lower=0, coords=[years], name="added_float")
    retired_float = m.add_variables(lower=0, coords=[years], name="retired_float")

    m.add_constraints(
        (
            (status_existing_large*existing_large).sum("id")
            + added_float - retired_float
            == caps_per_year.loc[bus, carrier]
            
        ),
        name="capacity_balance"
    )
    
    m.add_objective(
        (
            (- status_existing_large)*existing_large*weights_large 
            + retired_float*weight_float 
            + added_float*weight_float 
        ).sum()    
    )
        
    m.add_constraints(
        added_float >= added_float.shift({"year":1}),
        name="new_small_plants_constraint"
    )
    
    m.add_constraints(
        added_float >= retired_float,
        name="less_retirement_than_addition"
    )
    
    
    m.add_constraints(
        retired_float >= retired_float.shift({"year":1}),
        name="retired_small_plants_constraint"
    )
    
    m.add_constraints(
        status_existing_large.shift({"year":-1}) <=status_existing_large,
        name="decommissioning_large_plants_constraint"
    )

    m.solve(
        solver_name="cplex", 
        **solver_options
    )

    return m


# In[4]:


def decompose_residual_cap(i, residual_cap, blocks, block_sizes):
    
    m = linopy.Model()
    years = residual_cap.columns
    
    vin = m.add_variables(binary=True, coords = [years, blocks], name="entry")
    vout = m.add_variables(binary=True, coords=[years, blocks], name="exit")
    total_cap_year = m.add_variables(lower=0, name="total_cap_year", coords=[years])
    
    m.add_constraints(vin >= vin.shift({"year":1}))
    m.add_constraints(vout >= vout.shift({"year":1}))
    m.add_constraints(vin >= vout)  
    m.add_constraints(((vin-vout)*block_sizes).sum("block") == residual_cap.loc[i])
    
    m.add_objective(total_cap_year*1)
    m.solve(solver_name="cplex", **solver_options)

    return m


# In[5]:


def decompose_power_plant_table():

    status = pd.Series(index=caps_per_year.index)
    
    status_existing_large = []
    added_float = []
    retired_float = []
    
    for i in caps_per_year.index:
    
        print(i)
        
        m = run_disaggregation(i[0], i[1])
        status.loc[i] = m.termination_condition
        status_existing_large.append(m.solution["status_existing_large"].to_series())
        added_float.append(m.solution["added_float"].to_series())
        retired_float.append(m.solution["retired_float"].to_series())
    
    if (status != "optimal").sum() == 0:
        print("Large decomposition finished successfully")
    else:
        raise ValueError('Not all decompositions converged')
    
    eraa_plants_detailed = power_plants[["bus", "carrier", "Capacity", "Efficiency"]].copy()
    eraa_plants_detailed.columns = eraa_plants_detailed.columns.str.replace("Efficiency", "efficiency").str.replace("Capacity", "p_nom")
    
    status_existing_large = pd.concat(status_existing_large, axis=0).unstack(0).round().astype(int)
    
    fully_retired = status_existing_large[~status_existing_large.sum(axis=1).astype("bool")].index
    eraa_plants_detailed.drop(fully_retired, inplace=True)
    status_existing_large.drop(fully_retired,errors="ignore", inplace=True)
    
    eraa_plants_detailed["entry"] = target_years[0]
    eraa_plants_detailed["exit"] = np.nan
    
    eraa_plants_detailed = eraa_plants_detailed.loc[status_existing_large.index]
    eraa_plants_detailed.loc[status_existing_large[target_years[-1]] == 1, "exit"] = target_years[-1] + 1
    
    for i in eraa_plants_detailed[eraa_plants_detailed.exit.isna()].index:
        eraa_plants_detailed.loc[i, "exit"] = status_existing_large.loc[i].iloc[1:][status_existing_large.loc[i].diff().iloc[1:].round(1) == -1].index[0]
    
    added_float = pd.concat(added_float, keys=caps_per_year.index).unstack(2)
    retired_float = pd.concat(retired_float, keys=caps_per_year.index).unstack(2)
    capacity_balance = added_float.subtract(retired_float)
    
    n_plants = capacity_balance.unstack(0).div(typical_size,axis=0).stack(1).reorder_levels([1,0]).apply(np.floor).sort_index()
    
    change_in_plants = n_plants.diff(axis=1)
    change_in_plants[target_years[0]] = n_plants[target_years[0]]
    
    
    residual_cap = capacity_balance.subtract(
        n_plants.unstack(0)
        .multiply(typical_size, axis=0)
        .stack(1).reorder_levels([1,0]).sort_index()
    ).round()
    
    residual_cap = residual_cap.round(-1)
    
    change_resid_cap = residual_cap.diff(axis=1)
    change_resid_cap[target_years[0]] = residual_cap[target_years[0]]
    change_resid_cap[target_years[-1] + 1] = -change_resid_cap.sum(axis=1)
    
    status = pd.Series(index = residual_cap.index)
    
    entry = []
    exit = []
    block = []
    
    for i in residual_cap.index:
        
        block_sizes = change_resid_cap.loc[i].abs()[change_resid_cap.loc[i].abs() > 0 ].values
        
        block_sizes = np.append(block_sizes,
            change_resid_cap.loc[i].abs().sort_values().diff().abs()[change_resid_cap.loc[i].abs().sort_values().diff().abs()>0].values
        )
        
        blocks = pd.Index(range(len(block_sizes)), name="block")
        
        block_sizes = pd.Series(block_sizes, blocks)
        
        m = decompose_residual_cap(i, residual_cap, blocks, block_sizes)
    
        status.loc[i] = m.termination_condition
        entry.append(m.solution["entry"].to_series())
        exit.append(m.solution["exit"].to_series())
        block.append(block_sizes)
    
    if (status != "optimal").sum() == 0:
        print("Residual decomposition finished successfully")
    else:
        raise ValueError('Not all decompositions converged')
    
    entry = pd.concat(entry, keys=residual_cap.index).unstack("year")
    exit = pd.concat(exit, keys=residual_cap.index).unstack("year")
    entry_years = pd.Series(np.nan, index=entry.index)
    exit_years = pd.Series(np.nan, index=exit.index)
    
    for i in entry.index:
        if entry.loc[i].max()>0:
            entry_years.loc[i] = entry.loc[i][entry.loc[i]>0].index[0]
        if exit.loc[i].max()>0:
            exit_years.loc[i] = exit.loc[i][exit.loc[i]>0].index[0]
    
    
    residual_plants = pd.concat(block, keys=residual_cap.index).to_frame("p_nom")
    residual_plants["entry"] = entry_years
    residual_plants["exit"] = exit_years
    residual_plants = residual_plants[~residual_plants.entry.isna()]
    residual_plants.exit = residual_plants.exit.fillna(target_years[-1] + 1)
    entry_plants = change_in_plants.clip(lower=0).astype(int)
    exit_plants = change_in_plants.clip(upper=0).astype(int)
    exit_plants[target_years[-1] + 1] = - entry_plants.sum(axis=1).add(exit_plants.sum(axis=1))
    
    for i in change_in_plants[change_in_plants.abs().sum(axis=1)>0].index:
        
    
        bus = i[0]
        carrier = i[1]
        p_nom = typical_size[carrier]
    
        entry_years = entry_plants.loc[i][entry_plants.loc[i]>0].cumsum()
        exit_years = -exit_plants.loc[i][exit_plants.loc[i]!=0].cumsum()
        
        entry_years = pd.Series(
            entry_years.index, entry_years.values
        ).reindex(
            range(1, entry_years.astype(int).max()+1)
        ).bfill().astype(int)
        
        entry_years.index -=1
        
        exit_years = pd.Series(
            exit_years.index, exit_years.values
        ).reindex(
            range(1, exit_years.astype(int).max()+1)
        ).bfill().astype(int)
        
        exit_years.index -=1
    
        for number in entry_years.index:
    
            entry = entry_years.loc[number]
            exit = exit_years.loc[number]
            name = "{bus} {carrier} {number}".format(bus=bus, carrier=carrier, number=number)
    
            add_plant = pd.Series(
                [bus, carrier, p_nom, np.nan, entry, exit],
                index = eraa_plants_detailed.columns,
                name=name    
            )
    
            eraa_plants_detailed = pd.concat([eraa_plants_detailed, add_plant.to_frame().T],axis=0)
    
    
        
    
    for i in residual_plants.index:
    
        bus = i[0]
        carrier = i[1]
        number = i[2]
        p_nom = residual_plants.loc[i, "p_nom"]
        entry = residual_plants.loc[i, "entry"]
        exit = residual_plants.loc[i, "exit"]
    
        name = "{bus} small {carrier} {number}".format(bus=bus, carrier=carrier, number=number)
    
        add_plant = pd.Series(
                [bus, carrier, p_nom, np.nan, entry, exit],
                index = eraa_plants_detailed.columns,
                name=name    
            )
    
        eraa_plants_detailed = pd.concat([eraa_plants_detailed, add_plant.to_frame().T],axis=0)

    eraa_plants_detailed.exit = eraa_plants_detailed.exit.astype(int)
    eraa_plants_detailed.entry = eraa_plants_detailed.entry.astype(int)

    eraa_plants_detailed = eraa_plants_detailed[eraa_plants_detailed.entry != eraa_plants_detailed.exit]

    eraa_plants_detailed.columns.name = ""

    return eraa_plants_detailed, status_existing_large


# In[6]:


def build_thermal_properties():
    
    properties_raw = pd.read_excel(excel_file, "Common Data", index_col=[2,3], skiprows=10, header=0).dropna(how="all").iloc[2:, 1:].dropna(how="all", axis=1).iloc[:27, 1:17]
    
    properties_raw2 = pd.read_excel(excel_file, "Common Data", index_col=[2,3], skiprows=44, header=[0,3]).iloc[:, 1:].dropna(how="all", axis=1).dropna(how="all").iloc[:27, 1:17]
    
    properties_raw2.columns = [" ".join(i) for i in properties_raw2.columns]
    properties_raw = pd.concat([properties_raw, properties_raw2],axis=1)
    
    mask = properties_raw.index.get_level_values(0).str.contains("fuel cell", case=False) | \
           properties_raw.index.get_level_values(1).str.contains("fuel cell", case=False)
    properties_raw = properties_raw[~mask]
    properties_raw = properties_raw.fillna(0)
    
    properties = properties_raw.copy()
    
    properties.index = pd.MultiIndex.from_tuples([(i[0], i[1].replace("New", "new")) for i in properties.index])
    
    
    for col in properties.columns:
        try: 
            properties[col] = properties[col].astype(float)
        except:
            None
    
    properties.set_index(properties.index.remove_unused_levels(), inplace=True)
    
    properties_matching = pd.Series(
        ["OCGT", "coal", "oil", "hydrogen", "oil", "lignite", "nuclear", "oil"],
        properties.index.levels[0]
    )
    
    properties_matching = properties_matching.reindex([i[0] for i in properties.index])
    
    properties_matching.index = properties.index
    
    mask = (
    properties_matching.index.get_level_values(1).str.contains("CCGT") &
    (properties_matching.index.get_level_values(0) != "Hydrogen"))
    
    properties_matching.loc[mask] = "CCGT"
    
    age_categorization = pd.Series([i[1] for i in properties_matching.index], properties_matching.index, name="age")
    age_categorization = age_categorization.str.replace("CCGT ", "").str.replace("OCGT ", "").str.replace("conventional ", "")
    
    properties_matching = properties_matching.to_frame("carrier")
    properties_matching["age"] = age_categorization
    properties.index = pd.MultiIndex.from_arrays([properties_matching.carrier, properties_matching.age])
    
    properties.drop(("OCGT", "old"), inplace=True)
    
    age_matching = pd.DataFrame(index=properties.index, dtype=str)
    
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
    #age_matching.drop("hydrogen", level=0, inplace=True)
    age_matching.drop(('oil','-'), inplace=True)
    age_matching.loc["nuclear", "-"] = 1900, 2040

    return properties, age_matching
    
def add_technical_properties_existing(properties, age_matching):

    idx = pd.IndexSlice
    
    eraa_plants_detailed["age_class"] = np.nan
    
    for plant in status_existing_large.index:
    
        datein = power_plants.loc[plant, "DateIn"]
        carrier = power_plants.loc[plant, "carrier"]
        
        if (carrier in age_matching.index.levels[0]) and not np.isnan(datein):
            eraa_plants_detailed.loc[plant, "age_class"] = age_matching[(age_matching["from"] <= datein)&(age_matching["to"] > datein)].loc[carrier].index[0]
        elif (carrier in age_matching.index.levels[0]) and ("old 2" in age_matching.loc[carrier].index):
            eraa_plants_detailed.loc[plant, "age_class"] = "old 2"    
    
    eraa_plants_detailed.age_class = eraa_plants_detailed.age_class.fillna("new")
    eraa_plants_detailed.loc[eraa_plants_detailed.carrier =="nuclear", "age_class"] = "-"
    eraa_plants_detailed.loc[eraa_plants_detailed.carrier =="hydrogen", "age_class"] = "CCGT"
    
    missing_properties = properties.reindex([("OCGT", "old 2"), ("OCGT", "old 2")]).copy()
    missing_properties.index = [("biomass", "new"), ("other", "new")]
    missing_properties.loc[:,"CO2 emission factor"] = 0
    properties = pd.concat([properties, missing_properties])
    
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
        
        start_up_cost.append(properties.loc[idx[carrier, age_class], "Start-up fix cost (e.g. wear) warm start"])
        ramp_limit_up.append(properties.loc[idx[carrier, age_class], "Ramp up rate % of max output power / min"]*60)
        ramp_limit_down.append(properties.loc[idx[carrier, age_class], "Ramp down rate % of max output power / min"]*60)
        p_min_pu.append(properties.loc[idx[carrier, age_class], "Minimum stable generation (% of max power)"])
        min_up_time.append(properties.loc[idx[carrier, age_class], "Min Time on"])
        min_down_time.append(properties.loc[idx[carrier, age_class], "Min Time off"])
        var_om.append(properties.loc[idx[carrier, age_class], "Variable O&M cost"])
        efficiency.append(properties.loc[idx[carrier, age_class], "Standard efficiency in NCV terms"])
    
    eraa_plants_detailed.start_up_cost = start_up_cost
    
    eraa_plants_detailed.ramp_limit_down_real = ramp_limit_down
    eraa_plants_detailed.ramp_limit_down = pd.Series(ramp_limit_down).clip(upper=1.).values 
    eraa_plants_detailed.ramp_limit_up_real = ramp_limit_up
    eraa_plants_detailed.ramp_limit_up = pd.Series(ramp_limit_up).clip(upper=1.).values
    eraa_plants_detailed.p_min_pu = p_min_pu
    eraa_plants_detailed.min_up_time = min_up_time
    eraa_plants_detailed.min_down_time = min_down_time
    eraa_plants_detailed.var_om = var_om
    eraa_plants_detailed.efficiency = efficiency
    
    eraa_plants_detailed.p_nom = eraa_plants_detailed.p_nom.astype(float)
    eraa_plants_detailed.p_nom_max = eraa_plants_detailed.p_nom_max.astype(float)
    
def build_new_investments(df):

    technologies_for_investment = ["OCGT", "CCGT"]
    years = eraa_plants_detailed.entry.unique()
    zones_for_investment = eraa_plants.bus.unique()

    new = pd.DataFrame(
        0.01, 
        index= pd.MultiIndex.from_product(
            [zones_for_investment,technologies_for_investment, years], 
            names=["bus", "carrier", "entry"]), 
        columns=["p_nom"]
    )

    new.reset_index(inplace=True)
    new.index = new.bus + " " + new.carrier + " new " + new.entry.astype(str)
    new["efficiency"] = properties.loc[:, "new", :]["Standard efficiency in NCV terms"].reindex(new.carrier).values
    new["start_up_cost"] = properties["Start-up fix cost (e.g. wear) warm start"].loc[:, "new"].reindex(new.carrier, level=1).fillna(0).values
    new["ramp_limit_up_real"] = properties["Ramp up rate % of max output power / min"].loc[:, "new"].reindex(new.carrier).values*60
    new["ramp_limit_up"] = new["ramp_limit_up_real"].clip(upper=1.)
    new["ramp_limit_down_real"] = properties["Ramp down rate % of max output power / min"].loc[:, "new"].reindex(new.carrier).values*60
    new["ramp_limit_down"] = new["ramp_limit_down_real"].clip(upper=1.)
    new["p_min_pu"] = properties["Minimum stable generation (% of max power)"].loc[:, "new"].reindex(new.carrier).values
    new["min_up_time"] = properties["Min Time on"].loc[:, "new"].reindex(new.carrier).values
    new["min_down_time"] = properties["Min Time off"].loc[:, "new"].reindex(new.carrier).values
    new["exit"] = target_years[-1] + 1
    new["p_nom_max"] = np.inf
    new["p_nom_min"] = 0.01
    new["var_om"] = properties.loc[:, "new",:].reindex(new.carrier)["Variable O&M cost"].values
    new["invest_status"] = "new"
    new["age_class"] = "new"
    
    return pd.concat([df, new])
    
def aggregate_small():

    small_plants = eraa_plants_detailed[eraa_plants_detailed.p_nom.round()<=de_minimis]
    
    agg_func = pd.Series(
        "mean",
        index = ['p_nom', 'efficiency',  'p_nom_max', 'p_nom_min', 'start_up_cost', 'ramp_limit_up',
                 'ramp_limit_down', 'p_min_pu', 'min_up_time', 'min_down_time', 'var_om', "age_class"]
    )
    
    agg_func.loc[["p_nom", "p_nom_max"]] = "sum"
    agg_func.loc["age_class"] = pd.Series.mode
    
    aggregated_small = small_plants.groupby(["bus", "carrier", "entry", "exit", 'invest_status']).agg(agg_func).reset_index()
    small_index = pd.Index(aggregated_small.bus + " small " + aggregated_small.carrier + " ")
    aggregated_small.index = small_index + pd.Series(small_index).to_frame().groupby(0).cumcount().astype(str)
    
    eraa_plants = eraa_plants_detailed.drop(small_plants.index)
    return pd.concat([eraa_plants, aggregated_small])


excel_file = pd.ExcelFile(snakemake.input.technology_data)

typical_size = pd.Series(snakemake.params.typical_size)
solver_options = {"mip.tolerances.mipgap" : 0.01}
aggregated_plants = pd.read_hdf(snakemake.input.dispatchable_capacities)
bidding_zones = gpd.read_file(snakemake.input.bidding_zones)
bidding_zones = bidding_zones.to_crs("EPSG:3035")
de_minimis = snakemake.params.de_minimis

solver = snakemake.config["solving"]["solver_decomposition"]["name"]
options = snakemake.config["solving"]["solver_decomposition"]["options"]
solver_options = snakemake.config["solving"]["solver_options"][options]

power_plants = prepare_power_plant_table()
target_years = sorted(aggregated_plants.entry.unique())

caps_per_year = pd.DataFrame()

for year in target_years:
    caps_per_year[year] = aggregated_plants.query("(entry <= @year) and (exit > @year)").groupby(["bus", "carrier"]).p_nom.sum()

caps_per_year.fillna(0, inplace=True)

eraa_plants_detailed, status_existing_large = decompose_power_plant_table()
properties, age_matching = build_thermal_properties()
add_technical_properties_existing(properties, age_matching)
eraa_plants = aggregate_small()
eraa_plants_detailed = build_new_investments(eraa_plants_detailed)
eraa_plants = build_new_investments(eraa_plants)

hdf_file = snakemake.output.individual_power_plants

os.makedirs(os.path.dirname(hdf_file), exist_ok=True)

eraa_plants_detailed.to_hdf(hdf_file, key="detailed")
eraa_plants.to_hdf(hdf_file, key="small_aggregated")


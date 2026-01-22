
import pandas as pd
import zipfile


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

    
    properties = properties[['Standard efficiency in NCV terms', 'CO2 emission factor', 'Variable O&M cost', 'Min Time on', 
            'Min Time off', 'Start-up fuel consumption - warm start', 'Start-up fix cost (e.g. wear) warm start',
            'Unavailability %', 'Unavailability Days','Unavailability number of days','Unavailability % of annual number of days',
            'Minimum stable generation (% of max power)',"Ramp up rate % of max output power / min","Ramp down rate % of max output power / min"]]

    properties.columns = ["efficiency", "emission_factor", "var_OM", "min_up_time", "min_down_time",  
                          "start_up_fuel_consumption", "start_up_fix_cost", "forced_outage_share", "forced_outage_n_days" ,"maintenance_n_days", "maintenance_share_winter", "p_min_pu", "ramp_limit_up", "ramp_limit_down"]
    
    properties.loc[:, "ramp_limit_up"] *= 60
    properties.loc[:, "ramp_limit_down"] *= 60

    missing_properties = properties.reindex([("OCGT", "old 2"), ("OCGT", "old 2")]).copy()
    missing_properties.index = pd.MultiIndex.from_product([["biomass", "other"], ["-"]])
    missing_properties.emission_factor = 0

    properties = pd.concat([properties, missing_properties])


    return properties


# In[4]:


def build_default_data():

    data = []
    
    for sheet in invest_excel_default.sheet_names[1:]:
        data.append(pd.read_excel(invest_excel_default, sheet , index_col=[0,1,2]))
    
    data = pd.concat(data, keys=invest_excel_default.sheet_names[1:])
    
    data = data.reset_index("Node", drop=True).drop("Additional Comments", axis=1)
    
    age_class = ['new', 'new', 'new', 'old 1', 'old 2', 'present 1', 'present 2', 'new', 'old', 'old 1', 'old 2', 'CCGT', 'OCGT', 'new', 'old 1', 'old 2', 'old 1', 'old 2', '-', 
     'new', 'old','old 1', 'old 2', '-', 'new', 'old']
    
    matching = pd.Series(
        [
            "battery", "DSR", "CCGT", "CCGT","CCGT","CCGT","CCGT", "OCGT", "OCGT","OCGT","OCGT", "hydrogen", "hydrogen", "coal", "coal","coal",
            "oil","oil", "oil", "lignite","lignite","lignite", "lignite", "nuclear", "oil", "oil"
        ],
        index = data.index.levels[1]
    ).to_frame("carrier")
    
    matching["age_class"] = age_class
    
    data = data.unstack([0,2])
    data.index = pd.MultiIndex.from_frame(matching)
    data = data.stack([1,2])
    data.index.names = ['carrier', 'age_class', "indicator", 'reference_technology']

    return data


# In[5]:


def build_country_data():
    
    data = []
    
    for sheet in invest_excel_country.sheet_names[2:]:
        data.append(pd.read_excel(invest_excel_country, sheet , index_col=[0,1,2]).dropna(how="any", axis=1))
    
    data = pd.concat(data, keys=invest_excel_country.sheet_names[2:])
    data.drop("Source", axis=1, inplace=True)
    data = data.unstack([0,1,3])
    age_class = ['new', 'new', 'new', 'new', 'CCGT']
    
    matching = pd.Series(
            [
                "battery", "DSR", "CCGT",  "OCGT", "hydrogen"
            ],
            index = data.index
        ).to_frame("carrier")
    
    matching["age_class"] = age_class
    data.index = pd.MultiIndex.from_frame(matching)
    data = data.stack([1,2,3])
    
    return data


# In[6]:


data_folder = snakemake.params.data_folder

with zipfile.ZipFile(snakemake.input.thermal) as zip_f:
    zip_f.extractall(data_folder)


# In[7]:


with zipfile.ZipFile(snakemake.input.invest) as zip_f:
    zip_f.extractall(data_folder)

excel_file = pd.ExcelFile(data_folder + "Common data/Common Data.xlsx")
invest_excel_default = pd.ExcelFile(data_folder + "Economic and technical investment parameters/Economic and technical investment parameters_Default.xlsx")
invest_excel_country = pd.ExcelFile(data_folder + "Economic and technical investment parameters/Economic and technical investment parameters_Country Specific.xlsx")

properties = build_thermal_properties()
default = build_default_data()
country = build_country_data()

default.to_hdf(snakemake.output.economic_data, key="default")
country.to_hdf(snakemake.output.economic_data, key="country")
properties.to_hdf(snakemake.output.technology_parameters, key="technology_parameters")

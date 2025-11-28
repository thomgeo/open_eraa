from shutil import copyfile, move, rmtree, unpack_archive

configfile:	"config/config.yaml"
configfile:	"technology-data/config.yaml"

module technology_data:
    snakefile:  "technology-data/Snakefile"
    config:     config
    prefix:     "technology-data"


target_years = range(config["time_horizon"]["start"], config["time_horizon"]["end"])
climate_years = range(1, 2)

ruleorder: base_model > update_capacities

rule retrieve_all_data:
	input:
		#"data/pecd.zip",
		#"data/fb_domains.zip",
		#"data/pemmdb.zip",
		#"data/thermal.zip",
		#"resources/res_profile.h5",
		#"resources/demand.h5",
		#"resources/inflow.h5",	
		#"data/ntc.zip",
		#"resources/dispatchable_capacities.h5",
		#"resources/capacity_tables/individual_plants.h5",
		#"resources/maintenance_profiles.h5",
		#"resources/ntcs.h5",
		#"resources/dsr.h5",
		#"resources/all_capacities.h5",
		#"resources/investcap.h5",
		#"resources/networks/0/cy{cy}_ty{ty}.nc",
		#expand("resources/networks/0/cy{cy}_ty{ty}.nc", cy = climate_years, ty=target_years),
		#"data/ordc_parameters.h5",
		#"data/reserve_requirements.h5",
		#"data/NFB.h5",
		expand("results/lole/0/cy{cy}_ty{ty}.csv",cy = climate_years, ty=target_years)


rule retrieve_ntc:
	input:
		storage("https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/NTCs.zip"),
	output:
		protected("data/ntc.zip"),
	resources:	
		mem_mb=5000
	retries:
		2
	run:	
		move(input[0], output[0])

rule retrieve_pecd:
	input:
		storage("https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/PECD-RES.zip")
	output:	
		protected("data/pecd.zip"),
	resources:	
		mem_mb=5000
	retries:
		2
	run:	
		move(input[0], output[0])

rule retrieve_fb_domains:
	input:	storage("https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/FB_domains.zip"),
	output:	protected("data/fb_domains.zip")
	run:	move(input[0], output[0])

rule retrieve_pemmdb:
	input:	storage("https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/Dashboard_raw_data.zip"),
	output:	protected("data/pemmdb.zip"),
	resources:
		mem_mb=5000,
	retries:
		2
	run:	
		move(input[0], output[0])

rule retrieve_thermaldata:
	input:
		storage("https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/Common_data.zip"),
	output:
		protected("data/thermal.zip"),
	resources:
		mem_mb=5000,
	retries:
		2
	run:
		move(input[0], output[0])

rule retrieve_demand:
	input:
		storage("https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/Demand_data.zip")
	output:
		protected("data/demand.zip"),
	resources:
		mem_mb=5000,
	retries:
		2
	log:
		"logs/retrieve_demand.log"
	run:
		move(input[0], output[0])


rule retrieve_nuts24_shapes:
	input:
		storage("https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/ref-nuts-2024-01m.geojson.zip")
	output:
		protected("data/nuts24_regions.zip")
	retries:
		2
	run:
		move(input[0], output[0])
			
rule retrieve_nuts13_shapes:
	input:
		storage("https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/ref-nuts-2013-01m.geojson.zip")
	output:
		protected("data/nuts13_regions.zip")
	retries:
		2
	run:
		move(input[0], output[0])


rule retrieve_moldova_shapes:
	input:
		storage("https://data.humdata.org/dataset/3cd53554-3ad7-4aae-b862-9bbbc6fa3bfc/resource/e93ce536-41e6-41ed-a3b4-71268c1d677e/download/mda_admbnda_unhcr_20220510_shp.zip")
	output:
		protected("data/moldova_regions.zip")
	retries:
		2
	run:
		move(input[0], output[0])


rule retrieve_DSR:
	input:	storage("https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/Other_data.zip"),
	output:
		protected("data/dsr.zip"),
	resources:
		mem_mb=5000,
	retries:
		2
	log:
		"logs/retrieve_demand.log"
	run:
		move(input[0], output[0])

rule build_availability_factors:
	input:
		pecd = "data/pecd.zip",
	params:	
		data_folder="data/",
	output:
		res_profile = "resources/res_profile.h5",
	script:
		"scripts/build_availability_factors.py"

rule build_hydro_inflows:
	params:	
		data_folder = "data/"
	output:
		inflow = "resources/inflow.h5"
	script:
		"scripts/build_hydro_inflows.py"


rule buildFB:
	input:	FB_files = "data/fb_domains.zip",
		network = expand("resources/networks/0/cy{cy}_ty{ty}.nc", cy=climate_years, ty=target_years)
	params:	data_folder = "data/",
		climate_year=climate_years,
		target_year=target_years
	output:	FB = "resources/FB.h5",
	script:	"scripts/buildFB.py"

rule build_NordicFB:
	input:	fb_domains = "data/fb_domains.zip",
		network = expand("resources/networks/0/cy{cy}_ty{ty}.nc", cy=climate_years, ty=target_years)
	params:	folder = "data/",
	output:	NFB = "data/NFB.h5",
	script:	"scripts/build_NordicFB.py"

rule build_ntc:
	input:	ntc="data/ntc.zip",
	params:	folder = "data/"
	output:	save_hdf= "resources/ntcs.h5"
	script:	"scripts/build_ntc.py"

rule build_dispatchable_capacities:
	input:	
		dsr = "data/dsr.zip",
		pemmdb = "data/pemmdb.zip",
		thermal = "data/thermal.zip",
	params:	data_folder = "data/"
	output:	dispatchable_capacities = "resources/dispatchable_capacities.h5",
		all_capacities = "resources/all_capacities.h5",
		dsr = "resources/dsr.h5",
		investcap = "resources/investcap.h5",
	script:	"scripts/build_dispatchable_capacities.py"

rule build_demand:
	input:	
		demand_files = "data/demand.zip",
	params:
		demand_folder = "data/demand/",
	output:
		demand = "resources/demand.h5",
	script:	"scripts/build_demand.py"

rule build_bidding_zones:
	input:
		nuts24 = "data/nuts24_regions.zip",
		nuts13 = "data/nuts13_regions.zip",
		moldova = "data/moldova_regions.zip",
	params:
		extract_to = "data/nuts_regions/",
	output:
		bidding_zones = "resources/bidding_zones.geojson"
	script:	"scripts/build_bidding_zones.py"


rule build_individual_power_plants:
	input:	technology_data="data/Common data/Common Data.xlsx",
			bidding_zones = "resources/bidding_zones.geojson",
			dispatchable_capacities = "resources/dispatchable_capacities.h5",
	params:
			typical_size = config["power_plants"]["typical_size"],
			de_minimis = config["power_plants"]["de_minimis"],
	output:	
			individual_power_plants = "resources/capacity_tables/individual_plants.h5",
	script:	"scripts/build_individual_power_plants.py"

rule base_model:
	input:	
		commodity_prices = "data/Dashboard_raw_data/Commodity Prices.csv",
		battery = "data/Dashboard_raw_data/Batteries additional information.csv",
		dsr = "resources/dsr.h5",
		hydrodata = "data/Dashboard_raw_data/Hydro additional information.csv",
		inflow = "resources/inflow.h5",
		#pemmdb = "data/pemmdb.xlsx",
		demand = "resources/demand.h5",
		ntc = "resources/ntcs.h5",
		res_profile = "resources/res_profile.h5",
		technology_data = "technology-data/outputs/costs_2025.csv",
		dispatchable_plants = "resources/dispatchable_capacities.h5",
		individual_plants = "resources/capacity_tables/individual_plants.h5",
		all_capacities = "resources/all_capacities.h5",
		investcap = "resources/investcap.h5",
		maintenance="resources/maintenance_profiles.h5",
		#FB = "resources/FB.h5",
	params:
		ty = "{ty}",
		cy = "{cy}",
		years = target_years
	output:	
		network= "resources/networks/0/cy{cy}_ty{ty}.nc"
	script:	"scripts/base_model.py"	

rule solve_model:
	input:	network = "resources/networks/0/cy{cy}_ty{ty}.nc",
		reserve_requirements="data/reserve_requirements.h5",
		ordc_parameters = "data/ordc_parameters.h5",
		fb_domains = "data/fb_domains.zip",
		NFB = "data/NFB.h5",
	params:
		ty = "{ty}",
		cy = "{cy}",
		years = target_years,
		data_folder = "data/",
	output:
		solved_network= "results/networks/0/cy{cy}_ty{ty}.nc",
		revenues = "results/revenues/0/cy{cy}_ty{ty}.h5",
		lole = "results/lole/0/cy{cy}_ty{ty}.csv",
	script:	"scripts/solve_model.py"	


rule build_ordc_parameters:
	input:	pemmdb = "data/Dashboard_raw_data/Reserve requirements.csv",
	output:	ordc_hdf="data/ordc_parameters.h5",
		reserve_requirements="data/reserve_requirements.h5",
	script:	"scripts/build_ordc_parameters.py"

rule build_initial_capacity_table:
	input:	networks = expand("resources/networks/0/cy{cy}_ty{ty}.nc", cy=climate_years[0], ty=target_years)
	output:	initial_capacity_table = "resources/capacity_tables/0.csv"
	script: "scripts/build_initial_capacity_table.py"

rule build_maintenance_profiles:
	input:	power_plants="resources/capacity_tables/individual_plants.h5",
		demand = "resources/demand.h5",
		common_data = "data/Common data/Common Data.xlsx",
		all_caps = "resources/all_capacities.h5",
	output:	maintenance_profiles = "resources/maintenance_profiles.h5"
	script:	"scripts/build_maintenance_profiles.py"

rule update_capacities:
	input:	old_network = lambda w: ("resources/networks/{iter}/cy{cy}_ty{ty}.nc").format(iter = int(w.iter)-1, cy = w.cy, ty = w.ty),
		capacity_table = "resources/capacity_tables/{iter}.csv"
	params:	ty = "{ty}",
		cy = "{cy}",
	output:	new_network= "resources/networks/{iter}/cy{cy}_ty{ty}.nc"	
	script: "scripts/update_capacities.py"

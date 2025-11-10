from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from shutil import copyfile, move, rmtree, unpack_archive

HTTP = HTTPRemoteProvider()

configfile: "config/config.yaml"

target_years = range(config["time_horizon"]["start"], config["time_horizon"]["end"])
climate_years = range(1, 2)

ruleorder: base_model > update_capacities
#ruleorder: build_initial_capacity_table > capacity_adjustment

rule retrieve_all_data:
	input:	
		"data/pecd.zip",
		"data/fb_domains.zip",
		"data/pemmdb.zip",
		"data/thermal.zip",
		"resources/res_profile.h5",
		"resources/demand.h5",
		"resources/inflow.h5",	
		"data/ntc.zip",
		"resources/ntcs.h5",
		#"resources/FB.h5",
<<<<<<< HEAD
		"resources/dispatchable_capacities.h5",
		"resources/bidding_zones.geojson"
||||||| 3d004cc
		"resources/dispatchable_capacities.h5"	
=======
		"resources/dispatchable_capacities.h5",
		"resources/dsr.h5",
		"resources/all_capacities.h5",
		"resources/investcap.h5",
		"resources/networks/0/cy{cy}_ty{ty}.nc",
>>>>>>> main

subworkflow technology_data:
	workdir: "../technology-data"
	snakefile: "../technology-data/Snakefile"
	configfile: "../technology-data/config.yaml"


rule all:
	input:	expand("resources/networks/0/cy{cy}_ty{ty}.nc", cy = climate_years, ty=target_years)
		#"resources/capacity_tables/50.csv",
		#expand("resources/networks/0/cy{cy}_ty{ty}.nc", cy = climate_years, ty=target_years),


#2024 all links changed. PECD-RES

rule retrieve_pecd:
        input:
            HTTP.remote(
		"https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/PECD-RES.zip",
                keep_local=True,
                static=True,
            ),
        output:
            protected("data/pecd.zip"),
        resources:
            mem_mb=5000,
        retries: 2
	log:	"logs/retrieve_pecd.log"
        run:
            move(input[0], output[0])

rule retrieve_fb_domains:
        input:
            HTTP.remote(
		"https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/FB_domains.zip",
                keep_local=True,
                static=True,
            ),
        output:
            protected("data/fb_domains.zip"),
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])

rule retrieve_pemmdb:
        input:
            HTTP.remote(
		"https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/Dashboard_raw_data.zip",
                keep_local=True,
                static=True,
            ),
        output:
            protected("data/pemmdb.zip"),
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])

rule retrieve_thermaldata:
        input:
            HTTP.remote(
		"https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/Common_data.zip",
                keep_local=True,
                static=True,
            ),
        output:
            protected("data/thermal.zip"),
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])

rule retrieve_demand:
        input:
            HTTP.remote(
		"https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/Demand_data.zip",
                keep_local=True,
                static=True,
            ),
        output:
            protected("data/demand.zip"),
        resources:
            mem_mb=5000,
        retries: 2
	log:	"logs/retrieve_demand.log"
        run:
            move(input[0], output[0])

<<<<<<< HEAD
rule retrieve_nuts24_shapes:
		input:
			HTTP.remote(
		"https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/ref-nuts-2024-01m.geojson.zip",
			keep_local=True,
			static=True,
			)
		output:
			protected("data/nuts24_regions.zip")
		retries: 2
		run:
			move(input[0], output[0])
			
rule retrieve_nuts13_shapes:
		input:
			HTTP.remote(
		"https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/ref-nuts-2013-01m.geojson.zip",
			keep_local=True,
			static=True,
			)
		output:
			protected("data/nuts13_regions.zip")
		retries: 2
		run:
			move(input[0], output[0])




rule retrieve_moldova_shapes:
		input:
			HTTP.remote(
		"https://data.humdata.org/dataset/3cd53554-3ad7-4aae-b862-9bbbc6fa3bfc/resource/e93ce536-41e6-41ed-a3b4-71268c1d677e/download/mda_admbnda_unhcr_20220510_shp.zip",
			keep_local=True,
			static=True,
			)
		output:
			protected("data/moldova_regions.zip")
		retries: 2
		run:
			move(input[0], output[0])



||||||| 3d004cc
=======

rule retrieve_DSR:
        input:
            HTTP.remote(
		"https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/Other_data.zip",
                keep_local=True,
                static=True,
            ),
        output:
            protected("data/dsr.zip"),
        resources:
            mem_mb=5000,
        retries: 2
	log:	"logs/retrieve_demand.log"
        run:
            move(input[0], output[0])

>>>>>>> main
#rule retrieve_pecd22:
#        input:
#            HTTP.remote(
#		"https://eepublicdownloads.azureedge.net/clean-documents/sdc-documents/ERAA/2022/data-for-publication/Climate%20Data.zip",
#                keep_local=True,
#                static=True,
#            ),
#        output:
#            protected("data/pecd22.zip"),
#        resources:
#            mem_mb=5000,
#        retries: 2
#	log:	"logs/retrieve_pecd22.log"
#        run:
#            move(input[0], output[0])

rule retrieve_ntc:
        input:
            HTTP.remote(
		"https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/ERAA_2025/Preliminary_input_data/NTCs.zip",
                keep_local=True,
                static=True,
            ),
        output:
            protected("data/ntc.zip"),
        resources:
            mem_mb=5000,
        retries: 2
	log:	"logs/retrieve_ntc.log"
        run:
            move(input[0], output[0])

rule build_availability_factors:
	input:	pecd = "data/pecd.zip",
	params:	data_folder="data/",
	output:	res_profile = "resources/res_profile.h5",
	script:	"scripts/build_availability_factors.py"

#rule checkh5:
#	input:	inflow = "resources/inflow.h5",
#	params:	data_folder="data/",
#	output:	dispatchable_capacities = "resources/dispatchable_capacities.h5",
#	script:	"scripts/checkh5.py"

rule build_demand:
	input:	demand_files = "data/demand.zip",
	params:	demand_folder = "data/demand/",
	output:	demand = "resources/demand.h5",
	script:	"scripts/build_demand.py"

rule buildFB:
	input:	FB_files = "data/fb_domains.zip",
	params:	data_folder = "data/",
		climate_year=climate_years,
		target_year=target_years
	output:	FB = "resources/FB.h5",
	script:	"scripts/buildFB.py"

rule build_hydro_inflows:
	params:	data_folder = "data/"
	output: inflow = "resources/inflow.h5"
	script:	"scripts/build_hydro_inflows.py"

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

rule base_model:
	input:	
		commodity_prices = "data/Dashboard_raw_data/Commodity Prices.csv",
		battery = "data/Dashboard_raw_data/Batteries additional information.csv",
		dsr = "resources/dsr.h5",
		hydrodata = "data/Dashboard_raw_data/Hydro additional information.csv",
		inflow = "resources/inflow.h5",
		pemmdb = "data/ERAA2023 PEMMDB Generation.xlsx",
		demand = "resources/demand.h5",
		ntc = "resources/ntcs.h5",
		res_profile = "resources/res_profile.h5",
		technology_data = "outputs/costs_2025.csv",
		dispatchable_plants = "resources/dispatchable_capacities.h5",
		all_capacities = "resources/all_capacities.h5",
		investcap = "resources/investcap.h5",
		#FB = "resources/FB.h5",
	params:
		ty = "{ty}",
		cy = "{cy}",
		years = target_years
	output:	
		network= "resources/networks/0/cy{cy}_ty{ty}.nc"
	script:	"scripts/base_model.py"	


rule build_initial_capacity_table:
    input:  networks = expand("resources/networks/0/cy{cy}_ty{ty}.nc", cy=climate_years[0], ty=target_years)
    output: initial_capacity_table = "resources/capacity_tables/0.csv"
    script: "scripts/build_initial_capacity_table.py"

rule solve_model:
	input:	network = "resources/networks/{iter}/cy{cy}_ty{ty}.nc",
		pemmdb = "data/pemmdb.xlsx",
		ordc_parameters = "data/ordc_parameters.h5",
		fb_domains = "data/fb_domains.zip"
	params:
		ty = "{ty}",
		cy = "{cy}",
		years = target_years,
	output:
		solved_network= "results/networks/{iter}/cy{cy}_ty{ty}.nc",
		revenues = "results/revenues/{iter}/cy{cy}_ty{ty}.h5",
		lole = "results/lole/{iter}/cy{cy}_ty{ty}.csv",
	script:	"scripts/solve_model.py"	


rule build_ordc_parameters:
	input:	pemmdb = "data/Dashboard_raw_data/Reserve requirements.csv",
	output:	ordc_hdf="data/ordc_parameters.h5",
	script:	"scripts/build_ordc_parameters.py"

#rule capacity_adjustment:
#	input:	demand = "resources/demand.h5",
#		revenue_files = lambda w: expand("results/revenues/{iter}/cy{cy}_ty{ty}.h5", cy = climate_years, ty = target_years, iter = int(w.iter)-1),
#		previous_capacity_table = lambda w: ("resources/capacity_tables/{iter}.csv").format(iter = int(w.iter)-1),
#		costs = technology_data("outputs/costs_2025.csv"),
#		capacity_adjustment_size = "data/capacity_adjustement_size.csv"
#	params:	iteration = "{iter}"
#	output:	next_capacity_table = "resources/capacity_tables/{iter}.csv",
#		revenue_ratios_existing = "results/revenue_ratios/{iter}_existing.csv",
#		revenue_rarios_new = "results/revenue_ratios/{iter}_new.csv",
#	script:	"scripts/capacity_adjustment.py"
		
rule update_capacities:
	input:	old_network = lambda w: ("resources/networks/{iter}/cy{cy}_ty{ty}.nc").format(iter = int(w.iter)-1, cy = w.cy, ty = w.ty),
		capacity_table = "resources/capacity_tables/{iter}.csv"
	params:	ty = "{ty}",
		cy = "{cy}",
	output:	new_network= "resources/networks/{iter}/cy{cy}_ty{ty}.nc"	
	script: "scripts/update_capacities.py"		

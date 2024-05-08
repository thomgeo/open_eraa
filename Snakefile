from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from shutil import copyfile, move, rmtree, unpack_archive

HTTP = HTTPRemoteProvider()

configfile: "config/config.yaml"

target_years = range(config["time_horizon"]["start"], config["time_horizon"]["end"])
climate_years = range(1982, 1987)

ruleorder: base_model > update_capacities

#rule retrieve_all_data:
#	input:	"data/pecd22.zip",
#		"data/fb_domains.xlsx",
#		"data/pemmdb.xlsx",
#		"resources/res_profile.h5",
#		"resources/demand.h5",
#		"resources/inflow.h5",	
#		"data/ntc.zip",
#		"resources/ntcs.h5"	

subworkflow technology_data:
	workdir: "../technology-data"
	snakefile: "../technology-data/Snakefile"
	configfile: "../technology-data/config.yaml"


rule all:
	input:	"resources/capacity_tables/10.csv"
#		expand("results/networks/0/cy{cy}_ty{ty}.nc", cy = climate_years, ty=target_years),


rule retrieve_pecd:
        input:
            HTTP.remote(
		"https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/2023/PECD.zip",
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
		"https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/2023/FB-Domain-CORE_Merged.xlsx",
                keep_local=True,
                static=True,
            ),
        output:
            protected("data/fb_domains.xlsx"),
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])

rule retrieve_pemmdb:
        input:
            HTTP.remote(
		"https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/2023/ERAA2023%20PEMMDB%20Generation.xlsx",
                keep_local=True,
                static=True,
            ),
        output:
            protected("data/pemmdb.xlsx"),
        resources:
            mem_mb=5000,
        retries: 2
        run:
            move(input[0], output[0])


rule retrieve_demand:
        input:
            HTTP.remote(
		"https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/2023/Demand%20Dataset.zip",
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

rule retrieve_pecd22:
        input:
            HTTP.remote(
		"https://eepublicdownloads.azureedge.net/clean-documents/sdc-documents/ERAA/2022/data-for-publication/Climate%20Data.zip",
                keep_local=True,
                static=True,
            ),
        output:
            protected("data/pecd22.zip"),
        resources:
            mem_mb=5000,
        retries: 2
	log:	"logs/retrieve_pecd22.log"
        run:
            move(input[0], output[0])

rule retrieve_ntc:
        input:
            HTTP.remote(
		"https://eepublicdownloads.blob.core.windows.net/public-cdn-container/clean-documents/sdc-documents/ERAA/2023/Net%20Transfer%20Capacities.zip",
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

rule build_demand:
	input:	demand_files = "data/demand.zip",
	params:	demand_folder = "data/demand/",
	output:	demand = "resources/demand.h5",
	script:	"scripts/build_demand.py"

rule build_hydro_inflows:
	input:	pecd22 = "data/pecd22.zip",
	params:	unzip_folder = "data/hydro/"
	output: inflow = "resources/inflow.h5"
	script:	"scripts/build_hydro_inflows.py"

rule build_ntc:
	input:	ntc="data/ntc.zip",
	params:	folder = "data/"
	output:	save_hdf= "resources/ntcs.h5"
	script:	"scripts/build_ntc.py"

rule base_model:
	input:	
		commodity_prices = "data/Fuel and CO2 prices_ERAA2023.xlsx",
		inflow = "resources/inflow.h5",
		pemmdb = "data/pemmdb.xlsx",
		demand = "resources/demand.h5",
		ntc = "resources/ntcs.h5",
		res_profile = "resources/res_profile.h5",
	params:
		ty = "{ty}",
		cy = "{cy}",
		years = target_years
	output:	
		network= "resources/networks/0/cy{cy}_ty{ty}.nc"
	script:	"scripts/base_model.py"	

rule solve_model:
	input:	network = "resources/networks/{iter}/cy{cy}_ty{ty}.nc",
		pemmdb = "data/pemmdb.xlsx",
		ordc_parameters = "data/ordc_parameters.h5",
		fb_domains = "data/fb_domains.xlsx"
	params:
		ty = "{ty}",
		cy = "{cy}",
		years = target_years,
	output:
		solved_network= "results/networks/{iter}/cy{cy}_ty{ty}.nc",
		revenues = "results/revenues/{iter}/cy{cy}_ty{ty}.nc",
	script:	"scripts/solve_model.py"	


rule build_ordc_parameters:
	input:	pemmdb = "data/pemmdb.xlsx",
	output:	ordc_hdf="data/ordc_parameters.h5",
	script:	"scripts/build_ordc_parameters.py"

rule capacity_adjustment:
	input:	demand = "resources/demand.h5",
		revenue_files = lambda w: expand("results/revenues/{iter}/cy{cy}_ty{ty}.nc", cy = climate_years, ty = target_years, iter = int(w.iter)-1),
		previous_capacity_table = lambda w: ("resources/capacity_tables/{iter}.csv").format(iter = int(w.iter)-1),
		costs = technology_data("outputs/costs_2025.csv"),
	output:	next_capacity_table = "resources/capacity_tables/{iter}.csv"
	script:	"scripts/capacity_adjustment.py"
		
rule update_capacities:
	input:	old_network = lambda w: ("resources/networks/{iter}/cy{cy}_ty{ty}.nc").format(iter = int(w.iter)-1, cy = w.cy, ty = w.ty),
		capacity_table = "resources/capacity_tables/{iter}.csv"
	params:	ty = "{ty}",
		cy = "{cy}",
	output:	new_network= "resources/networks/{iter}/cy{cy}_ty{ty}.nc"	
	script: "scripts/update_capacities.py"		
# %%
import pandas as pd
import numpy as np
import pypsa
import glob
import os


iteration = 0
save_path = snakemake.output.initial_capacity_table
network_list = snakemake.input.networks

sorted(network_list)

networks = pd.Series(
    index = [i.split("ty")[1][:4] for i in network_list]
)

networks.index = networks.index.astype(int)
networks.index.name = "target_year"

for i, path in enumerate(network_list):
    networks.iloc[i] = pypsa.Network(path)

networks.sort_index(inplace=True)

generation_capacities = []
storage_capacities = []
for year in networks.index:
    generation_capacities.append(networks.loc[year].generators)
    storage_capacities.append(networks.loc[year].storage_units)

generation_capacities = generation_capacities.reorder_levels(["name", "target_year"]).sort_index()
storage_capacities = storage_capacities.reorder_levels(["name", "target_year"]).sort_index()   

os.makedirs(os.path.dirname(save_path), exist_ok=True)

generation_capacities.to_hdf(save_path, key="generators")
storage_capacities.to_hdf(save_path, key="storages")
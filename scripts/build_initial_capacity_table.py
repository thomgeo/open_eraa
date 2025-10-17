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

for i, path in enumerate(network_list):
    networks.iloc[i] = pypsa.Network(path)

networks.sort_index(inplace=True)
capacities = pd.DataFrame()

for year in networks.index:
    capacities = pd.concat(
        [
            capacities,
            (
                networks[year].generators
                .set_index(
                    pd.Series(year, networks[year].generators.index),
                    append=True
                )
            )
        ]
    )

os.makedirs(os.path.dirname(save_path), exist_ok=True)

capacities.to_csv(save_path)

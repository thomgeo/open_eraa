
import pandas as pd
import glob
import numpy as np
import os

def profits_new_plants(average_expected_revenue):
    
    revenue_lifetime_new = pd.Series(index = average_expected_revenue.filter(like="new").columns)

    lifetime_new = costs.loc[:, "lifetime", :].reindex([i.split(" ")[1] for i in revenue_lifetime_new.index]).value
    lifetime_new.index = revenue_lifetime_new.index

    for plant in revenue_lifetime_new.index:

        tech = initial_capacity_table.loc[plant, "carrier"].iloc[0]
        entry = initial_capacity_table.loc[plant, "entry"].iloc[0]

        initial_capacity_table.loc[plant]

        lifetime = lifetime_new.loc[plant]

        revenue_horizon = average_expected_revenue[plant].loc[entry: min(end_year, entry+lifetime)].sum()
        revenue_beyond = average_expected_revenue.loc[end_year,plant]*max(lifetime - (end_year - entry), 0)
        revenue_lifetime_new.loc[plant] = revenue_horizon + revenue_beyond

    invest_cost_new = costs.loc[:, "investment", :].reindex([i.split(" ")[1] for i in revenue_lifetime_new.index]).value.multiply(1e3) #converstion to EUR/MW
    invest_cost_new.index = revenue_lifetime_new.index




    fixed_om_new = (
        costs.loc[:, "FOM", :]
        .reindex([i.split(" ")[1] for i in revenue_lifetime_new.index])
        .value
        .multiply(
            costs.loc[:, "investment", :]
            .reindex([i.split(" ")[1] for i in revenue_lifetime_new.index])
            .value
            .multiply(1e3)
        )
    )

    fixed_om_new.index = revenue_lifetime_new.index

    total_fixed_new = invest_cost_new + fixed_om_new.multiply(lifetime_new)

    return revenue_lifetime_new.subtract(total_fixed_new)


# In[ ]:


def profits_existing_plants(average_expected_revenue):
    

    revenue_lifetime_existing = pd.DataFrame(
        index = average_expected_revenue.index, 
        columns = [i for i in average_expected_revenue if "new" not in i]
    )

    relative_fixed_om = costs.loc[:, "FOM", :].reindex([i.split(" ")[1] for i in revenue_lifetime_existing]).value

    relative_fixed_om.index = revenue_lifetime_existing.columns

    invest_cost = (
        costs.loc[:, "investment", :]
        .reindex([i.split(" ")[1] for i in revenue_lifetime_existing])
        .value
        .multiply(1e3) #converstion to EUR/MW
    ) 

    invest_cost.index = revenue_lifetime_existing.columns
    total_fixed_om = relative_fixed_om.multiply(invest_cost).divide(100) # O&M given in %
    

    total_fixed = pd.DataFrame(
            index = average_expected_revenue.index, 
            columns = [i for i in average_expected_revenue if "new" not in i]
        )

    for plant in revenue_lifetime_existing.columns:
        
        tech = initial_capacity_table.loc[plant, "carrier"].unique()[0]

        entry = int(initial_capacity_table.loc[plant, "entry"].unique()[0])
        exit = int(initial_capacity_table.loc[plant, "exit"].unique()[0])

        for year in range(max(start_year, entry), min(exit, end_year+1)):
            revenue_lifetime_existing.loc[year, plant] = average_expected_revenue.loc[start_year:year, plant].fillna(0).sum()
            total_fixed.loc[year, plant] = (year + 1 - entry)*total_fixed_om.loc[plant]

    revenue_lifetime_existing = revenue_lifetime_existing.astype(float)
    total_fixed = total_fixed.astype(float)

    return revenue_lifetime_existing.subtract(total_fixed)

        
        
iteration = int(snakemake.params.iteration)

adjustment_factor = snakemake.config["EVA"]["adjustment_factor_initial"]

offset = snakemake.config["EVA"]["offset"]
eva_technologies = snakemake.config["EVA"]["eva_technologies"]

demand = pd.read_hdf(snakemake.input.demand)
revenue_files = snakemake.input.revenue_files
initial_capacity_table = pd.read_csv(snakemake.input.previous_capacity_table, index_col=[0,1])

save_path_rr_exist = snakemake.output.revenue_ratios_existing
save_path_rr_new = snakemake.output.revenue_rarios_new

capacity_change = pd.read_csv(snakemake.input.capacity_adjustment_size, index_col=0).loc[iteration, "capacity_change"]


revenues = []
climate_years = []
target_years = []

for file in revenue_files:
    revenues.append(pd.read_hdf(file))

    climate_years.append(int(file.split("cy")[1][:4]))
    target_years.append(int(file.split("ty")[1][:4]))

revenues = pd.concat(revenues,axis=1).T

revenues.index = pd.MultiIndex.from_arrays([climate_years, target_years])

costs = pd.read_csv(snakemake.input.costs, index_col=[0,1])

costs = costs.reindex(eva_technologies, level=0)

average_expected_revenue = revenues.groupby(level=1).mean()

average_expected_revenue = average_expected_revenue[[i for i in average_expected_revenue.columns if i.split(" ")[1] in eva_technologies]]

average_expected_revenue = average_expected_revenue.div(initial_capacity_table.p_nom.unstack(0)[average_expected_revenue.columns])

start_year = revenues.index.levels[1][0]
end_year = revenues.index.levels[1][-1]

profits_new = profits_new_plants(average_expected_revenue)
profits_existing = profits_existing_plants(average_expected_revenue)

capacity_adjustment_new = profits_new*adjustment_factor
capacity_adjustment_existing = profits_existing*adjustment_factor

p_nom_initial = initial_capacity_table.p_nom.unstack(0)
p_nom_min = initial_capacity_table.p_nom_min.unstack(0)
p_nom_max = initial_capacity_table.p_nom_max.unstack(0)

#initial_capacity_changes = initial_capacity_table.p_nom.unstack(0).subtract(
#    initial_capacity_table.loc[:, start_year, :].p_nom
#)


next_capacities_new = pd.concat(
    [
        p_nom_initial[capacity_adjustment_new.index].add(capacity_adjustment_new).unstack(),
        p_nom_min[capacity_adjustment_new.index].unstack()
    ],
    axis=1
).max(axis=1).dropna()

capacity_adjustment_existing = (
    capacity_adjustment_existing
    .reindex(capacity_adjustment_existing.index[::-1])
    .cummax()
    .sort_index()
)

next_capacities_existing = pd.concat(
    [
        p_nom_initial[capacity_adjustment_existing.columns].add(capacity_adjustment_existing).unstack(),
        p_nom_max[capacity_adjustment_existing.columns].unstack()
    ],
    axis=1
).min(axis=1)

next_capacities_existing = pd.concat(
    [
        next_capacities_existing,
        p_nom_min[capacity_adjustment_existing.columns].unstack(),
    ],
    axis=1
).max(axis=1)

next_capacity_table = initial_capacity_table.copy()
next_capacities_existing.dropna(inplace=True)
next_capacities_new.dropna(inplace=True)
next_capacity_table.loc[next_capacities_existing.index, "p_nom"] = next_capacities_existing

next_capacity_table.loc[next_capacities_new.index, "p_nom"] = next_capacities_new

next_capacity_table.to_csv(snakemake.output.next_capacity_table)

os.makedirs(os.path.dirname(save_path_rr_exist), exist_ok=True)

profits_existing.to_csv(save_path_rr_exist)
profits_new.to_frame().to_csv(save_path_rr_new)

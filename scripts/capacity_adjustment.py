
import pandas as pd
import glob
import numpy as np

def scaling(input_data, offset, adjustment_factor):
    return (
        np.log(
            (input_data - offset)
            .clip(1e-6))
        .multiply(adjustment_factor)
    )


def capacity_adjustment_new_plants(average_expected_revenue):
    
    revenue_lifetime_new = pd.DataFrame(index = average_expected_revenue.index, columns = average_expected_revenue.filter(like="new").columns)

    for plant in average_expected_revenue.filter(like="new").columns:

        tech = plant.split(" ")[1]
        lifetime = costs.loc[tech, "lifetime"].value


        for year in revenues.index.levels[1]:

            years_beyond_horizon = lifetime - (end_year - year)
            revenue_lifetime_new.loc[year, plant] = (
                average_expected_revenue[plant].loc[year:].fillna(0).sum() 
                + average_expected_revenue[plant].fillna(0).loc[end_year]*years_beyond_horizon
            )

    revenue_lifetime_new = revenue_lifetime_new.astype(float)

    invest_cost_new = costs.loc[:, "investment", :].reindex([i.split(" ")[1] for i in revenue_lifetime_new]).value.multiply(1e3) #converstion to EUR/MW
    invest_cost_new.index = revenue_lifetime_new.columns

    lifetime_new = costs.loc[:, "lifetime", :].reindex([i.split(" ")[1] for i in revenue_lifetime_new]).value
    lifetime_new.index = revenue_lifetime_new.columns

    fixed_om_new = (
        costs.loc[:, "FOM", :]
        .reindex([i.split(" ")[1] for i in revenue_lifetime_new])
        .value
        .multiply(
            costs.loc[:, "investment", :]
            .reindex([i.split(" ")[1] for i in revenue_lifetime_new])
            .value
            .multiply(1e3)
        )
    )

    fixed_om_new.index = revenue_lifetime_new.columns

    total_fixed_new = invest_cost_new + fixed_om_new.multiply(lifetime_new)

    demand_scaling_factor = (
        demand.groupby(level=2).mean()
        .reindex(revenue_lifetime_new.columns.str[:4])
    )

    demand_scaling_factor.index = revenue_lifetime_new.columns

    return scaling(
        revenue_lifetime_new.div(total_fixed_new),
        offset,
        adjustment_factor    
    ).multiply(
        demand_scaling_factor
    )


def capacity_adjustment_existing_plants(average_expected_revenue):

    revenue_lifetime_existing = pd.DataFrame(
        index = average_expected_revenue.index, 
        columns = [i for i in average_expected_revenue if "new" not in i]
    )

    for plant in [i for i in average_expected_revenue if "new" not in i]:

        tech = plant.split(" ")[1]

        for year in revenues.index.levels[1]:
            revenue_lifetime_existing.loc[year, plant] = average_expected_revenue.loc[start_year:year, plant].fillna(0).sum()

    revenue_lifetime_existing = revenue_lifetime_existing.astype(float)

    fixed_om = costs.loc[:, "FOM", :].reindex([i.split(" ")[1] for i in revenue_lifetime_existing]).value

    fixed_om.index = revenue_lifetime_existing.columns

    invest_cost = (
        costs.loc[:, "investment", :]
        .reindex([i.split(" ")[1] for i in revenue_lifetime_existing])
        .value
        .multiply(1e3) #converstion to EUR/MW
    ) 

    invest_cost.index = revenue_lifetime_existing.columns

    demand_scaling_factor = (
            demand.groupby(level=2).mean()
            .reindex(revenue_lifetime_existing.columns.str[:4])
        )

    demand_scaling_factor.index = revenue_lifetime_existing.columns

    total_fixed_existing = fixed_om.multiply(invest_cost)

    return scaling(
        revenue_lifetime_existing.divide(total_fixed_existing),
        offset,
        adjustment_factor    
    ).multiply(
        demand_scaling_factor
    )
        

adjustment_factor = snakemake.config["EVA"]["adjustment_factor"]
offset = snakemake.config["EVA"]["offset"]
eva_technologies = snakemake.config["EVA"]["eva_technologies"]

demand = pd.read_hdf(snakemake.input.demand)
revenue_files = snakemake.input.revenue_files
initial_capacity_table = pd.read_csv(snakemake.input.previous_capacity_table, index_col=[0,1])

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

start_year = revenues.index.levels[1][0]
end_year = revenues.index.levels[1][-1]

capacity_adjustment_new = capacity_adjustment_new_plants(average_expected_revenue)

capacity_adjustment_existing = capacity_adjustment_existing_plants(average_expected_revenue)

p_nom_initial = initial_capacity_table.p_nom.unstack(0)

p_nom_min = initial_capacity_table.p_nom_min.unstack(0)
p_nom_max = initial_capacity_table.p_nom_max.unstack(0)

initial_capacity_changes = initial_capacity_table.p_nom.unstack(0).subtract(
    initial_capacity_table.loc[:, start_year, :].p_nom
)


next_capacities_new = pd.concat(
    [
        p_nom_initial[capacity_adjustment_new.columns].add(capacity_adjustment_new).unstack(),
        p_nom_min[capacity_adjustment_new.columns].unstack()
    ],
    axis=1
).max(axis=1).unstack(1).cummax().stack()

capacity_adjustment_existing = (
    capacity_adjustment_existing.subtract(
        initial_capacity_changes[capacity_adjustment_existing.columns]
    ).reindex(capacity_adjustment_existing.index[::-1])
    .cummax()
    .sort_index()
    .add(initial_capacity_changes[capacity_adjustment_existing.columns])
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



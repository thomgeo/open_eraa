#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd


# ## Load data

year = int(snakemake.params.ty)

scenario = int(snakemake.params.op)

technology_parameters = pd.read_hdf(snakemake.input.technology_parameters)

plants = pd.read_hdf(snakemake.input.power_plants, "detailed")

maintenance_profile = pd.read_hdf(snakemake.input.maintenance_profiles, key="maintenance{}".format(year))

save_hdf = snakemake.output.outage_patterns



# # Prepare and quality check input data

outage_data = plants.join(technology_parameters[['forced_outage_share', 'forced_outage_n_days', 'maintenance_n_days']], on=['carrier', 'age_class'])
outage_data.forced_outage_n_days = outage_data.forced_outage_n_days.astype(int)

# Check if all plants have forced outage parameters
if outage_data[outage_data.forced_outage_share.isnull() | outage_data.forced_outage_n_days.isnull()].shape[0] != 0:
    raise ValueError(
                "There are plants with missing forced outage parameters!"
            )



# Check if number of planned maintenance days match target vales in technology parameters

if outage_data[abs(outage_data['maintenance_n_days'].subtract(
    -maintenance_profile.sum().subtract(
        maintenance_profile.index.shape[0]).div(24))) > 5].shape[0] != 0:
    raise ValueError(
                "The number of maintenance days deviates by more than 5 days from the target number!"
            )



#Check if all generators have a cluster information attached

idx = [i for i in outage_data.index if 'new' not in str(i)]

if outage_data.loc[idx][outage_data.loc[idx].cluster.isnull()].shape[0] > 0:
    raise ValueError(
                "The are generators without a cluster!"
            )


# ## Simulate forced outages

# Outages for generators detailed
key=[]

np.random.seed(year*scenario)

hours_per_year = len(maintenance_profile.index)
ls_generators = maintenance_profile.columns
np_forcedoutages = np.ones((len(ls_generators), hours_per_year), dtype=int)

for i, gen in enumerate(ls_generators):

    np_forcedoutages[i] = maintenance_profile[gen]

    duration = outage_data.loc[gen].forced_outage_n_days * 24 # convert to hours
    num_maintenance_hours = outage_data.loc[gen].maintenance_n_days * 24
    num_outages_year = outage_data.loc[gen].forced_outage_share * (hours_per_year - num_maintenance_hours) / duration
    num_drawn_outages = int(num_outages_year * 2)
    mean_interval = (hours_per_year - num_maintenance_hours)/num_outages_year

    # Generate inter-arrival times (in hours) using exponential distribution
    times = np.cumsum(np.random.exponential(mean_interval, num_drawn_outages))

    # Only keep outage start times within the year
    outage_starts = times[times < hours_per_year].astype(int)

    # No forced outages during planned maintenance
    outage_starts = outage_starts[maintenance_profile[gen][outage_starts] != 0]

    for start in outage_starts:
        end = int(min(start + duration, hours_per_year))
        np_forcedoutages[i, start:end] = 0

    key.append((gen))

pd_unavailability_profiles = pd.DataFrame(np_forcedoutages, index=key).T
pd_unavailability_profiles.index.name = 'hour'



# Outage profile for aggregates as average outages weighted by generator capacity

pd_unavailability_profiles_cluster = pd_unavailability_profiles.multiply(
    outage_data.p_nom[pd_unavailability_profiles.columns]).T.groupby(outage_data.cluster).sum().T.div(
        outage_data.p_nom[pd_unavailability_profiles.columns].groupby(outage_data.cluster).sum())



# # Quality check and readout

# Check if mean outage days do not deviate more than 10% from target given in technology parameters 

mean_offdays_carrier_target = (outage_data.groupby(outage_data.carrier).maintenance_n_days.mean() +
    outage_data.forced_outage_share.groupby(outage_data.carrier).mean().multiply(
        365 - outage_data.groupby(outage_data.carrier).maintenance_n_days.mean()
    ))

mean_offdays_carrier = (1 - pd_unavailability_profiles).sum().groupby(outage_data.carrier).mean() / 24

if (1 - mean_offdays_carrier.div(mean_offdays_carrier_target).abs() > 0.1).any():
    raise ValueError(
                "The outage rate deviates more than 10% from target!"
            )


# Check if mean outage days for clusters do not deviate more than 10% from target given in technology parameters 

mean_offdays_cluster_carrier = (1 - pd_unavailability_profiles_cluster).sum().groupby(
    outage_data.groupby('cluster').carrier.max()).mean() / 24

if ((1 - mean_offdays_cluster_carrier.div(mean_offdays_carrier)).abs() > 0.1).any():
    raise ValueError(
                "The clustered outage rate deviates more than 10% from target!"
            )


pd_unavailability_profiles.to_hdf(save_hdf, key="detailed")
pd_unavailability_profiles_cluster.to_hdf(save_hdf, key="aggregated")

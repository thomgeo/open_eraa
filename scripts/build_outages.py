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
                "The number of maintenance days deviates by more than 5 from the target number!"
            )


# Check if all generators have a cluster information attached

# if outage_data[outage_data.cluster.isnull()].shape[0] > 0:
#     raise ValueError(
#                 "The are generators without a cluster!"
#             )


# ## Simulate forced outages

# Outages for single generator
key=[]

np.random.seed(year*scenario)

hours_per_year = len(maintenance_profile.index)
ls_generators = maintenance_profile.columns #[86, 95] #
np_forcedoutages = np.ones((len(ls_generators), hours_per_year), dtype=int)

for i, gen in enumerate(ls_generators):

    duration = outage_data.loc[gen].forced_outage_n_days * 24 # convert to hours
    num_outages_year = outage_data.loc[gen].forced_outage_share * hours_per_year / duration
    num_drawn_outages = int(num_outages_year * 2)
    mean_interval = hours_per_year/num_outages_year

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


# Outage profile for clusters as average outages weighted by generator capacity
pd_unavailability_profiles_cluster = pd_unavailability_profiles.multiply(
    outage_data.p_nom).T.groupby(outage_data.cluster).sum().T.div(
        outage_data.p_nom.groupby(outage_data.cluster).sum())



# # Quality check and readout

# Check if outage rate does not deviate more than 5pp from target given in technology parameters 
if ((1 - pd_unavailability_profiles).sum().div(len(pd_unavailability_profiles)).groupby(
    outage_data.carrier).mean().subtract(
        outage_data.forced_outage_share.groupby(outage_data.carrier).mean()).abs() > 0.05).any():
    raise ValueError(
                "The outage rate deviates more than 5pp from target in technology parameter"
            )


pd_unavailability_profiles.to_hdf(save_hdf, key="single")
pd_unavailability_profiles_cluster.to_hdf(save_hdf, key="cluster")


import pandas as pd
import numpy
import pypsa
from pypsa.descriptors import get_switchable_as_dense as as_dense

def average_every_nhours(n, offset):
    m = n.copy(with_time=False)

    snapshot_weightings = n.snapshot_weightings.resample(offset).sum()
    m.set_snapshots(snapshot_weightings.index)
    m.snapshot_weightings = snapshot_weightings

    for c in n.iterate_components():
        pnl = getattr(m, c.list_name + "_t")
        for k, df in c.pnl.items():
            if not df.empty:
                pnl[k] = df.resample(offset).mean()

    return m


def replace_su(network, su_to_replace):

    su = network.storage_units.loc[su_to_replace]

    bus_name = "{} {}".format(su["bus"], su["carrier"])
    link_1_name = "{} discharger".format(su_to_replace, su["carrier"])
    link_2_name = "{} charger".format(su_to_replace, su["carrier"])
    store_name = "{} store".format(su_to_replace, su["carrier"])
    gen_name = "{} inflow".format(su_to_replace)

    if bus_name not in network.buses.index:
        network.add("Bus", bus_name, carrier=su["carrier"])

    # dispatch link
    network.add(
        "Link",
        link_1_name,
        bus0=bus_name,
        bus1=su["bus"],
        capital_cost=su["capital_cost"] * su["efficiency_dispatch"],
        p_nom=su["p_nom"] / su["efficiency_dispatch"],
        p_nom_extendable=su["p_nom_extendable"],
        p_nom_max=su["p_nom_max"] / su["efficiency_dispatch"],
        p_nom_min=su["p_nom_min"] / su["efficiency_dispatch"],
        p_max_pu=su["p_max_pu"],
        marginal_cost=su["marginal_cost"] * su["efficiency_dispatch"],
        efficiency=su["efficiency_dispatch"],
        carrier = su["carrier"],
        type="discharging",
    )

    # store link
    network.add(
        "Link",
        link_2_name,
        bus1=bus_name,
        bus0=su["bus"],
        p_nom=su["p_nom"],
        p_nom_extendable=su["p_nom_extendable"],
        p_nom_max=su["p_nom_max"],
        p_nom_min=su["p_nom_min"],
        p_max_pu=-su["p_min_pu"],
        efficiency=su["efficiency_store"],
        carrier = su["carrier"],
        type="charging"
    )

    if (
        su_to_replace in network.storage_units_t.state_of_charge_set.columns
        and (
            ~pd.isnull(network.storage_units_t.state_of_charge_set[su_to_replace])
        ).any()
    ):
        e_max_pu = pd.Series(data=1.0, index=network.snapshots)
        e_min_pu = pd.Series(data=0.0, index=network.snapshots)
        non_null = ~pd.isnull(
            network.storage_units_t.state_of_charge_set[su_to_replace]
        )
        e_max_pu[non_null] = network.storage_units_t.state_of_charge_set[su_to_replace][
            non_null
        ]
        e_min_pu[non_null] = network.storage_units_t.state_of_charge_set[su_to_replace][
            non_null
        ]
    else:
        e_max_pu = 1.0
        e_min_pu = 0.0

    network.add(
        "Store",
        store_name,
        bus=bus_name,
        e_nom=su["p_nom"] * su["max_hours"],
        e_nom_min=su["p_nom_min"] / su["efficiency_dispatch"] * su["max_hours"],
        e_nom_max=su["p_nom_max"] / su["efficiency_dispatch"] * su["max_hours"],
        e_nom_extendable=su["p_nom_extendable"],
        e_max_pu=e_max_pu,
        e_min_pu=e_min_pu,
        standing_loss=su["standing_loss"],
        e_cyclic=su["cyclic_state_of_charge"],
        e_initial=su["state_of_charge_initial"],
        carrier = su["carrier"]
    )


    # inflow from a variable generator, which can be curtailed (i.e. spilled)
    inflow_max = as_dense(network, "StorageUnit", "inflow").max()[su_to_replace]
    
    if inflow_max == 0.0:
        inflow_pu = 0.0
    else:
        inflow_pu = network.storage_units_t.inflow[su_to_replace] / inflow_max

    if inflow_max >0:
        network.add(
            "Generator",
            gen_name,
            bus=bus_name,
            carrier="rain",
            p_nom=inflow_max,
            p_max_pu=inflow_pu,
        )
        
    network.add(
        "Generator",
        gen_name + " loss-of-charge",
        bus=bus_name,
        carrier= su["carrier"],
        p_nom=su["p_nom"] * su["max_hours"],
        marginal_cost = 100000,
        type="inflow",
        )

    network.remove("StorageUnit", su_to_replace)

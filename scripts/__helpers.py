import pandas as pd
import numpy as np
import pypsa
from pypsa.descriptors import get_switchable_as_dense as as_dense

def average_every_nhours(n, offset):
    m = n.copy(snapshots=[])

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

    bus_name = list(map(' '.join, zip(su["bus"], su["carrier"])))
    link_1_name = [s1 + " converter " + s2 + " to AC" for s1, s2 in zip(su_to_replace, su["carrier"])]
    link_2_name = [s1 + " converter AC to " + s2 for s1, s2 in zip(su_to_replace, su["carrier"])]
    store_name = [s1 + " store " + s2 for s1, s2 in zip(su_to_replace, su["carrier"])]
    gen_name = [s1 + " inflow" for s1 in list(su_to_replace)]
 

    buses = pd.DataFrame(index = bus_name, data=su["carrier"].values, columns=["carrier"])
    buses = buses[~buses.index.duplicated()]

    network.add("Bus", buses.index, **buses)

    dispatch_links = pd.DataFrame(index=link_1_name)
    dispatch_links["bus0"] = bus_name
    dispatch_links["bus1"]=list(su["bus"])
    dispatch_links["capital_cost"]=list(su["capital_cost"] * su["efficiency_dispatch"])
    dispatch_links["p_nom"]=list(su["p_nom"] / su["efficiency_dispatch"])
    dispatch_links["p_nom_extendable"]=list(su["p_nom_extendable"])
    dispatch_links["p_nom_max"]=list(su["p_nom_max"] / su["efficiency_dispatch"])
    dispatch_links["p_nom_min"]=list(su["p_nom_min"] / su["efficiency_dispatch"])
    dispatch_links["p_max_pu"]=list(su["p_max_pu"])
    dispatch_links["carrier"]=list(su["carrier"])
    dispatch_links["marginal_cost"]=list(su["marginal_cost"] * su["efficiency_dispatch"])
    dispatch_links["efficiency"]=list(su["efficiency_dispatch"])

 
    # dispatch link
    network.add( "Link", dispatch_links.index, **dispatch_links)
 

    store_links = pd.DataFrame(index=link_2_name)

    store_links["bus1"]=bus_name
    store_links["bus0"]=list(su["bus"])
    store_links["p_nom"]=list(su["p_nom"])
    store_links["p_nom_extendable"]=list(su["p_nom_extendable"])
    store_links["p_nom_max"]=list(su["p_nom_max"])
    store_links["p_nom_min"]=list(su["p_nom_min"])
    store_links["p_max_pu"]=list(-su["p_min_pu"])
    store_links["carrier"]=list(su["carrier"])
    store_links["efficiency"]=list(su["efficiency_store"])

    # store link
    network.add("Link", store_links.index, **store_links)

    e_max_pu = pd.DataFrame(data=np.ones((len(network.snapshots), len(bus_name))), index=network.snapshots, columns= store_name)
    e_min_pu = pd.DataFrame(data=np.zeros((len(network.snapshots), len(bus_name))), index=network.snapshots, columns= store_name)
 
    idx_stateofcharge = su_to_replace[su_to_replace.isin(network.storage_units_t.state_of_charge_set.columns)]
    idx_stateofcharge = idx_stateofcharge[(~pd.isnull(network.storage_units_t.state_of_charge_set[idx_stateofcharge])).any() == True]
 
    e_max_pu[idx_stateofcharge] = network.storage_units_t.state_of_charge_set[idx_stateofcharge]
    e_min_pu[idx_stateofcharge] = network.storage_units_t.state_of_charge_set[idx_stateofcharge]
 
    e_max_pu[e_max_pu.isnull()] = 1
    e_min_pu[e_max_pu.isnull()] = 0
 

    stores = pd.DataFrame(index=store_name)

    stores["bus"]=bus_name
    stores["e_nom"]=list(su["p_nom"] * su["max_hours"])
    stores["e_nom_min"]=list(su["p_nom_min"] / su["efficiency_dispatch"] * su["max_hours"])
    stores["e_nom_max"]=list(su["p_nom_max"] / su["efficiency_dispatch"] * su["max_hours"])
    stores["e_nom_extendable"]=list(su["p_nom_extendable"])
    stores["carrier"]=list(su["carrier"])
    stores["standing_loss"]=list(su["standing_loss"])
    stores["e_cyclic"]=list(su["cyclic_state_of_charge"])
    stores["e_initial"]=list(su["state_of_charge_initial"])

    network.add( "Store", stores.index, **stores, e_max_pu=e_max_pu, e_min_pu=e_min_pu)
 

    # inflow from a variable generator, which can be curtailed (i.e. spilled)
    s_inflowmax = as_dense(network, "StorageUnit", "inflow").max()[su_to_replace]
    s_inflowmax_idx = s_inflowmax[s_inflowmax > 0].index
    su_inflowmax = network.storage_units.loc[s_inflowmax_idx]
    inflow_max = list(s_inflowmax[s_inflowmax_idx])
    inflow_pu = network.storage_units_t.inflow[s_inflowmax_idx] / s_inflowmax[s_inflowmax_idx]

    bus_name_inflowmax = list(map(' '.join, zip(su_inflowmax["bus"], su_inflowmax["carrier"])))
    gen_name_inflowmax = [s1 + " inflow" for s1 in list(s_inflowmax_idx)]

    inflow_pu.columns = gen_name_inflowmax
 

    inflow = pd.DataFrame(index = gen_name_inflowmax)

    inflow["bus"]=bus_name_inflowmax
    inflow["carrier"]=["inflow"] * len(bus_name_inflowmax)
    inflow["p_nom"]=inflow_max


    network.add("Generator",inflow.index, **inflow, p_max_pu=inflow_pu)
 

    loss_of_charge = pd.DataFrame(index = [s1 + " loss-of-charge" for s1 in list(gen_name)])

    loss_of_charge["bus"]=bus_name
    loss_of_charge["carrier"]=["storage-content"] * len(bus_name)
    loss_of_charge["p_nom"]=list(su["p_nom"] * su["max_hours"])
    loss_of_charge["marginal_cost"] = 20000

 
    network.add("Generator", loss_of_charge.index, **loss_of_charge)
 
    network.remove("StorageUnit", su_to_replace)

"""

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
"""
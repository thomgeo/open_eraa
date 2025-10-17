import pypsa
import pandas as pd


def solve_partial_model(n, n_horizon, horizon, solver_name, solver_options, path, e_initial, e_final_path):

    e_final = pd.read_csv(e_final_path, index_col=0, parse_dates=True)
    
    def storage_targets(n, snapshots):
        m = n.model
    
        rhs = m.constraints["Store-energy_balance"].rhs
        rhs = rhs.to_dataframe()
        rhs.loc[snapshots[-1], "rhs"] = e_final.loc[snapshots[-1]].values    
        m.constraints["Store-energy_balance"].rhs = rhs.rhs.to_xarray()
    
    p = n.copy()
    
    p.set_snapshots(n.snapshots[n_horizon*horizon:horizon*(n_horizon+1)])

    p.stores.e_initial = e_initial.loc[p.snapshots[0]]
    
    p.optimize(solver_name="cplex", 
               assign_all_duals=True,
               extra_functionality=storage_targets
              );

    p.export_to_netcdf(path)

    return p.model.termination_condition
# Weekly BYM model with covariates
# Snow-transition models fitted to the full data
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import coo_matrix, csr_matrix, diags, bmat, block_diag
from sksparse.cholmod import cholesky
from polyagamma import random_polyagamma
import pyreadr
import pickle
from joblib import Parallel, delayed
from pathlib import Path
import multiprocessing

# Set this path to the local data and results directory.
BASE_DIR = Path("path/to/snow/data-and-results")

burn = 3000
thin = 2
tot_save = 1000
total_iters = burn + thin * tot_save
n_chains = 1
period = 52
# Data
def load_data():
    snow = pyreadr.read_r(BASE_DIR / "snow_cleaned_full.Rda")
    snow = list(snow.values())[0].reset_index(drop=True)

    coords = snow.iloc[:, :2].to_numpy()
    y = snow.iloc[:, 2:].to_numpy()

    return coords, y
# BUILD DATASET
def build_dataset(event):

    coords, y = load_data()

    gdf = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(coords[:,0], coords[:,1]),
        crs="EPSG:4326"
    ).to_crs("+proj=aeqd +lat_0=90 +lon_0=-100")

    xy = np.vstack([gdf.geometry.x, gdf.geometry.y]).T / 1e6

    W = (squareform(pdist(xy)) <= 0.22).astype(int)
    np.fill_diagonal(W, 0)
    W = csr_matrix(W)

    S, TT = y.shape
    print(f"Using FULL S = {S}")

    deg = np.array(W.sum(axis=1)).flatten()
    Q_car = diags(deg) - W
    I_S = diags(np.ones(S))

    t_full = np.arange(1, TT+1)

    # >>> CHANGED: scale t
    t_scaled_full = (t_full - t_full.mean()) / t_full.std()

    # >>> CHANGED: square first, then scale
    t2_full = t_full**2
    t2_scaled_full = (t2_full - t2_full.mean()) / t2_full.std()

    if event == "p01":
        loc_mask = (y[:, :-1] == 0)
        kappa_builder = lambda ny: ny - 0.5
    else:
        loc_mask = (y[:, :-1] == 1)
        kappa_builder = lambda ny: (1 - ny) - 0.5

    row_idx, time_idx = np.where(loc_mask)
    N = len(row_idx)

    next_y = y[row_idx, time_idx+1]
    kappa = kappa_builder(next_y)

    t_raw = time_idx + 1

    # >>> CHANGED
    t_scaled = t_scaled_full[time_idx]
    t2_scaled = t2_scaled_full[time_idx]

    week_idx = (t_raw - 1) % 52

    # >>> CHANGED: cov4 -> cov5
    cov5 = np.column_stack([
        np.ones(N),
        np.cos(2*np.pi*t_raw/period),
        np.sin(2*np.pi*t_raw/period),
        t_scaled,
        t2_scaled
    ])

    # >>> CHANGED
    K_base = 5
    K_total = 2 * K_base
    # covariates (UNCHANGED)
    lat = (coords[:,1] - coords[:,1].mean()) / coords[:,1].std()

    no_nbs = np.array([
        57,170,236,269,343,685,946,947,989,
        1037,1084,1090,1109,1118,1127,1176,1203
    ]) - 1

    elev_raw = pd.read_csv(BASE_DIR/"curr_elev.csv").iloc[:,3].to_numpy()
    nnbs_elev = pd.read_csv(BASE_DIR/"nnbs_elev.csv", sep="\t").iloc[:,2].to_numpy()

    mask = np.ones(S, dtype=bool)
    mask[no_nbs] = False

    elev_all = np.zeros(S)
    elev_all[mask] = elev_raw
    elev_all[no_nbs] = nnbs_elev
    elev = (elev_all - elev_all.mean()) / elev_all.std()

    snow_temp = pyreadr.read_r(BASE_DIR/"snow_temp_full.Rda")
    temp = list(snow_temp.values())[0].iloc[:,2:].to_numpy()
    temp_scaled = (temp - temp.mean()) / temp.std()

    # interactions (UNCHANGED)
    t_lat = t_scaled * lat[row_idx]
    t_elev = t_scaled * elev[row_idx]
    t_temp = t_scaled * temp_scaled[row_idx, time_idx]
    # X_eta
    eta_dim = K_total*S + 3

    rows_e, cols_e, vals_e = [], [], []

    for i in range(N):
        s = row_idx[i]

        # >>> CHANGED
        for k in range(K_base):
            xval = cov5[i,k]
            rows_e += [i,i]
            cols_e += [(2*k)*S+s, (2*k+1)*S+s]
            vals_e += [xval,xval]

        rows_e += [i,i,i]
        cols_e += [K_total*S, K_total*S+1, K_total*S+2]
        vals_e += [t_lat[i],t_elev[i],t_temp[i]]

    X_eta_base = coo_matrix((vals_e,(rows_e,cols_e)),
                            shape=(N,eta_dim)).tocsr()

    # Block indices
    rows_all, cols_all = X_eta_base.nonzero()
    sp_mask = cols_all < K_total*S

    rows_sp = rows_all[sp_mask]
    cols_sp = cols_all[sp_mask]

    j = cols_sp // S
    w = week_idx[rows_sp]
    tau_indices = j*52 + w
    # X_tau
    tau_dim = K_total * 52

    rows_t, cols_t, vals_t = [], [], []

    for i in range(N):
        w = week_idx[i]

        # >>> CHANGED
        for k in range(K_base):
            xval = cov5[i,k]
            rows_t += [i,i]
            cols_t += [(2*k)*52+w, (2*k+1)*52+w]
            vals_t += [xval,xval]

    X_tau_base = coo_matrix((vals_t,(rows_t,cols_t)),
                            shape=(N,tau_dim)).tocsr()
    # PRIOR
    Q_blocks = []
    for k in range(K_base):
        Q_blocks.append(Q_car)
        Q_blocks.append(I_S)

    Q_spatial = bmat(
        [[Q_blocks[i] if i==j else None for j in range(K_total)]
         for i in range(K_total)],
        format="csr"
    )

    Q_factor = diags(np.ones(3)/9)
    Q_eta = block_diag((Q_spatial, Q_factor), format="csr")

    tau_prior_prec = diags(np.ones(tau_dim)/9)

    return (
        X_eta_base, X_tau_base, Q_eta, tau_prior_prec,
        row_idx, week_idx, kappa,
        eta_dim, tau_dim, K_total, S,
        sp_mask, tau_indices
    )
# MCMC (UNCHANGED)
def run_chain(chain_id, data, event):

    (
        X_eta_base, X_tau_base, Q_eta, tau_prior_prec,
        row_idx, week_idx, kappa,
        eta_dim, tau_dim, K_total, S,
        sp_mask, tau_indices
    ) = data

    np.random.seed(2000 + chain_id)

    curr_eta = np.zeros(eta_dim)
    curr_tau = np.ones(tau_dim)

    all_eta = np.zeros((eta_dim, tot_save), dtype=np.float32)
    all_tau = np.zeros((tau_dim, tot_save), dtype=np.float32)

    save_idx = 0

    for it in tqdm(range(total_iters),
                   desc=f"{event}-Chain {chain_id}",
                   position=chain_id):

        X_tilde_eta = X_eta_base.copy()
        X_tilde_eta.data[sp_mask] *= curr_tau[tau_indices]

        psi = X_tilde_eta @ curr_eta
        omega = random_polyagamma(1, psi)

        XtOmega = X_tilde_eta.T.multiply(omega)
        post_prec_eta = (XtOmega @ X_tilde_eta + Q_eta).tocsc()

        factor = cholesky(post_prec_eta, mode="simplicial")
        mu = factor.solve_A(X_tilde_eta.T @ kappa)

        z = np.random.randn(eta_dim)
        z = z / np.sqrt(factor.D())
        z = factor.solve_Lt(z)
        z = factor.apply_Pt(z)

        curr_eta = mu + z

        X_tilde_tau = X_tau_base.copy()

        rows_nz, cols_nz = X_tilde_tau.nonzero()
        j = cols_nz // 52
        s = row_idx[rows_nz]
        eta_indices = j*S + s

        X_tilde_tau.data *= curr_eta[eta_indices]

        gamma = curr_eta[K_total*S:K_total*S+3]
        X_factor = X_eta_base[:, K_total*S:K_total*S+3]

        c = X_factor @ gamma
        residual = kappa - omega*c

        XtOmega = X_tilde_tau.T.multiply(omega)
        post_prec_tau = (XtOmega @ X_tilde_tau + tau_prior_prec).tocsc()

        factor = cholesky(post_prec_tau, mode="simplicial")
        mu = factor.solve_A(X_tilde_tau.T @ residual)

        z = np.random.randn(tau_dim)
        z = z / np.sqrt(factor.D())
        z = factor.solve_Lt(z)
        z = factor.apply_Pt(z)

        curr_tau = mu + z

        if it >= burn and (it-burn)%thin==0:
            all_eta[:,save_idx] = curr_eta
            all_tau[:,save_idx] = curr_tau
            save_idx += 1
            if save_idx == tot_save:
                break

    with open(BASE_DIR / f"{event}_weekly_covt2_chain{chain_id}.pkl","wb") as f:
        pickle.dump({"eta":all_eta,"tau":all_tau},f)

    return chain_id
# MAIN
def main():

    for event in ["p01", "p10"]:
        print(f"\nBuilding dataset for {event} ...")
        data = build_dataset(event)

        print(f"Running {event} ...")

        Parallel(n_jobs=n_chains, backend="loky", batch_size=1)(
            delayed(run_chain)(i, data, event) for i in range(n_chains)
        )

    print("All done.")


if __name__ == "__main__":
    multiprocessing.set_start_method("spawn", force=True)
    main()

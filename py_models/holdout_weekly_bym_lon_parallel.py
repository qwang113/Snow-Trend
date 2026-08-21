# Weekly BYM model with covariates
# + NA vs Non-NA longitude (separate scaling)
# Spatial model with aligned tau indices
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

# Set these paths for the local system.
BASE_DIR = Path("path/to/snow/data")
OUTPUT_DIR = Path("path/to/snow/results") / "holdout_38_14"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

burn = 10000
thin = 2
tot_save = 5000
total_iters = burn + thin * tot_save
n_chains = 10
period = 52
SEED_BASE = 2000
N_TRAIN_YEARS = 38
N_TEST_YEARS = 14
TRAIN_WEEKS = N_TRAIN_YEARS * period
TEST_WEEKS = N_TEST_YEARS * period
PRED_THIN = 15

CHAIN_START = 0
CHAIN_END   = 10
# Data
def load_data():
    snow = pyreadr.read_r(BASE_DIR / "snow_cleaned_full.Rda")
    snow = list(snow.values())[0].reset_index(drop=True)

    coords = snow.iloc[:, :2].to_numpy()
    y_full = snow.iloc[:, 2:].to_numpy()
    assert y_full.shape[1] == (N_TRAIN_YEARS + N_TEST_YEARS) * period
    y = y_full[:, :TRAIN_WEEKS]

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
    t_scaled_full = (t_full - t_full.mean()) / t_full.std()

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
    t_scaled = t_scaled_full[time_idx]

    week_idx = (t_raw - 1) % 52

    cov4 = np.column_stack([
        np.ones(N),
        np.cos(2*np.pi*t_raw/period),
        np.sin(2*np.pi*t_raw/period),
        t_scaled
    ])

    K_base = 4
    K_total = 8
    # covariates
    lat = (coords[:,1] - coords[:,1].mean()) / coords[:,1].std()

    lon_raw = coords[:,0]
    region = np.zeros(S, dtype=int)
    region[lon_raw < -30] = 0
    region[lon_raw >= -30] = 1

    lon_na = np.zeros(S)
    lon_euas = np.zeros(S)

    mask_na = (region == 0)
    mask_euas = (region == 1)

    lon_na[mask_na] = (
        (lon_raw[mask_na] - lon_raw[mask_na].mean()) /
        max(lon_raw[mask_na].std(), 1e-6)
    )

    lon_euas[mask_euas] = (
        (lon_raw[mask_euas] - lon_raw[mask_euas].mean()) /
        max(lon_raw[mask_euas].std(), 1e-6)
    )

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
    temp_full = list(snow_temp.values())[0].iloc[:,2:].to_numpy()
    temp = temp_full[:, :TRAIN_WEEKS]
    temp_scaled = (temp - temp.mean()) / temp.std()

    t_lon_NA   = t_scaled * lon_na[row_idx]
    t_lon_EUAS = t_scaled * lon_euas[row_idx]
    t_lat  = t_scaled * lat[row_idx]
    t_elev = t_scaled * elev[row_idx]
    t_temp = t_scaled * temp_scaled[row_idx, time_idx]
    # X_eta
    eta_dim = K_total*S + 5

    rows_e, cols_e, vals_e = [], [], []

    for i in range(N):
        s = row_idx[i]

        for k in range(K_base):
            xval = cov4[i,k]
            rows_e += [i,i]
            cols_e += [(2*k)*S+s, (2*k+1)*S+s]
            vals_e += [xval,xval]

        rows_e += [i,i,i,i,i]
        cols_e += [K_total*S, K_total*S+1, K_total*S+2, K_total*S+3, K_total*S+4]
        vals_e += [t_lon_NA[i], t_lon_EUAS[i], t_lat[i], t_elev[i], t_temp[i]]

    coo = coo_matrix((vals_e,(rows_e,cols_e)), shape=(N,eta_dim))
    print("\n===== COO DEBUG CHECK =====")

    idx_check = np.random.choice(len(coo.data), size=20, replace=False)

    for k in idx_check:
        r = coo.row[k]
        c = coo.col[k]
        v = coo.data[k]

        is_spatial = (c < K_total*S)

        print(f"\n[k={k}] row={r}, col={c}, val={v:.4f}, spatial={is_spatial}")

        if is_spatial:
            j = c // S
            w = week_idx[r]
            tau = j*52 + w

            print(f"   → j={j}, week={w}, tau_idx={tau}")

            base_k = j // 2
            print(f"   → base_k={base_k}, cov4={cov4[r, base_k]:.4f}")

        else:
            gamma_idx = c - K_total*S
            print(f"   → gamma_idx={gamma_idx}")

            expected = [
                t_lon_NA[r],
                t_lon_EUAS[r],
                t_lat[r],
                t_elev[r],
                t_temp[r]
            ][gamma_idx]

            print(f"   → expected={expected:.4f}")
    sp_mask = coo.col < K_total*S
    j = coo.col[sp_mask] // S
    w = week_idx[coo.row[sp_mask]]
    tau_indices = j*52 + w

    X_eta_base = coo.tocsr()
    # X_tau
    tau_dim = K_total * 52

    rows_t, cols_t, vals_t = [], [], []

    for i in range(N):
        w = week_idx[i]
        for k in range(K_base):
            xval = cov4[i,k]
            rows_t += [i,i]
            cols_t += [(2*k)*52+w, (2*k+1)*52+w]
            vals_t += [xval,xval]

    X_tau_base = coo_matrix((vals_t,(rows_t,cols_t)), shape=(N,tau_dim)).tocsr()

    Q_blocks = []
    for k in range(K_base):
        Q_blocks.append(Q_car)
        Q_blocks.append(I_S)

    Q_spatial = bmat(
        [[Q_blocks[i] if i==j else None for j in range(K_total)]
         for i in range(K_total)],
        format="csr"
    )

    Q_factor = diags(np.ones(5)/9)
    Q_eta = block_diag((Q_spatial, Q_factor), format="csr")
    tau_prior_prec = diags(np.ones(tau_dim)/9)

    return (
        X_eta_base, X_tau_base, Q_eta, tau_prior_prec,
        row_idx, week_idx, kappa,
        eta_dim, tau_dim, K_total, S,
        sp_mask, tau_indices
    )
# MCMC
def run_chain(chain_id, data, event):

    (
        X_eta_base, X_tau_base, Q_eta, tau_prior_prec,
        row_idx, week_idx, kappa,
        eta_dim, tau_dim, K_total, S,
        sp_mask, tau_indices
    ) = data

    chain_seed = SEED_BASE + chain_id
    np.random.seed(chain_seed)
    pg_rng = np.random.default_rng(chain_seed)

    curr_eta = np.zeros(eta_dim)
    curr_tau = np.ones(tau_dim)

    all_eta = np.zeros((eta_dim, tot_save), dtype=np.float32)
    all_tau = np.zeros((tau_dim, tot_save), dtype=np.float32)

    save_idx = 0

    for it in tqdm(
        range(total_iters),
        desc=f"{event}-Chain {chain_id}",
        position=chain_id,
        leave=True,
        dynamic_ncols=False
    ):

        X_tilde_eta = X_eta_base.copy()
        X_tilde_eta.data[sp_mask] *= curr_tau[tau_indices]

        psi = X_tilde_eta @ curr_eta
        omega = random_polyagamma(1, psi, random_state=pg_rng)

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

        gamma = curr_eta[K_total*S:K_total*S+5]
        X_factor = X_eta_base[:, K_total*S:K_total*S+5]

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

    with open(OUTPUT_DIR / f"{event}_weekly_cov+lon_train38_chain{chain_id}.pkl","wb") as f:
        pickle.dump({"eta":all_eta,"tau":all_tau},f)

    return chain_id


def predict_chain(chain_id):
    snow = pyreadr.read_r(BASE_DIR / "snow_cleaned_full.Rda")
    snow = list(snow.values())[0].reset_index(drop=True)
    coords = snow.iloc[:, :2].to_numpy()
    y_full = snow.iloc[:, 2:].to_numpy()
    S = y_full.shape[0]

    with open(OUTPUT_DIR / f"p01_weekly_cov+lon_train38_chain{chain_id}.pkl", "rb") as f:
        d01 = pickle.load(f)
    with open(OUTPUT_DIR / f"p10_weekly_cov+lon_train38_chain{chain_id}.pkl", "rb") as f:
        d10 = pickle.load(f)

    eta01_all = d01["eta"][:, ::PRED_THIN]
    eta10_all = d10["eta"][:, ::PRED_THIN]
    assert eta01_all.shape[1] == eta10_all.shape[1] == 334
    gamma01 = eta01_all[8*S:, :]
    gamma10 = eta10_all[8*S:, :]
    eta01 = eta01_all[:8*S, :].reshape(8, S, 334)
    eta10 = eta10_all[:8*S, :].reshape(8, S, 334)
    tau01 = d01["tau"][:, ::PRED_THIN].reshape(8, period, 334)
    tau10 = d10["tau"][:, ::PRED_THIN].reshape(8, period, 334)

    lon_raw = coords[:, 0]
    mask_na = lon_raw < -30
    mask_euas = ~mask_na
    lon_na = np.zeros(S)
    lon_euas = np.zeros(S)
    lon_na[mask_na] = (lon_raw[mask_na] - lon_raw[mask_na].mean()) / max(lon_raw[mask_na].std(), 1e-6)
    lon_euas[mask_euas] = (lon_raw[mask_euas] - lon_raw[mask_euas].mean()) / max(lon_raw[mask_euas].std(), 1e-6)
    lat = (coords[:, 1] - coords[:, 1].mean()) / coords[:, 1].std()

    no_nbs = np.array([
        57,170,236,269,343,685,946,947,989,
        1037,1084,1090,1109,1118,1127,1176,1203
    ]) - 1
    elev_raw = pd.read_csv(BASE_DIR / "curr_elev.csv").iloc[:, 3].to_numpy()
    nnbs_elev = pd.read_csv(BASE_DIR / "nnbs_elev.csv", sep="\t").iloc[:, 2].to_numpy()
    keep = np.ones(S, dtype=bool)
    keep[no_nbs] = False
    elev_all = np.zeros(S)
    elev_all[keep] = elev_raw
    elev_all[no_nbs] = nnbs_elev
    elev = (elev_all - elev_all.mean()) / elev_all.std()

    snow_temp = pyreadr.read_r(BASE_DIR / "snow_temp_full.Rda")
    temp_full = list(snow_temp.values())[0].iloc[:, 2:].to_numpy()
    temp_train = temp_full[:, :TRAIN_WEEKS]
    temp_scaled = (temp_full - temp_train.mean()) / temp_train.std()

    y_prev_test = y_full[:, TRAIN_WEEKS - 1:-1]
    assert y_prev_test.shape == (S, TEST_WEEKS)
    train_t = np.arange(1, TRAIN_WEEKS + 1)
    transition_t = np.arange(TRAIN_WEEKS, TRAIN_WEEKS + TEST_WEEKS)
    transition_scaled = (transition_t - train_t.mean()) / train_t.std(ddof=0)

    output = OUTPUT_DIR / f"pred_weekly_bym_lon_train38_test14_chain{chain_id}.npy"
    p_snow = np.lib.format.open_memmap(
        output, mode="w+", dtype=np.float32,
        shape=(334, S, TEST_WEEKS)
    )

    for j, (t_raw, t_scaled) in enumerate(zip(transition_t, transition_scaled)):
        week = (t_raw - 1) % period
        cov = np.array([
            1.0, 1.0,
            np.cos(2 * np.pi * t_raw / period),
            np.cos(2 * np.pi * t_raw / period),
            np.sin(2 * np.pi * t_raw / period),
            np.sin(2 * np.pi * t_raw / period),
            t_scaled, t_scaled
        ])
        phi01 = np.sum(cov[:, None, None] * eta01 * tau01[:, week, None, :], axis=0)
        phi10 = np.sum(cov[:, None, None] * eta10 * tau10[:, week, None, :], axis=0)
        z = (lon_na, lon_euas, lat, elev, temp_scaled[:, TRAIN_WEEKS + j])
        phi01 += t_scaled * sum(z[k][:, None] * gamma01[k] for k in range(5))
        phi10 += t_scaled * sum(z[k][:, None] * gamma10[k] for k in range(5))
        p01 = 1.0 / (1.0 + np.exp(-phi01))
        p10 = 1.0 / (1.0 + np.exp(-phi10))
        p_snow[:, :, j] = np.where(
            y_prev_test[:, j, None] == 0, p01, 1.0 - p10
        ).T.astype(np.float32)

    p_snow.flush()
    return chain_id
# MAIN
def main():
    for event in ["p01", "p10"]:
        print(f"\nBuilding dataset for {event} ...")
        data = build_dataset(event)

        print(f"Running {event} ...")

        Parallel(
            n_jobs=CHAIN_END-CHAIN_START,
            backend="loky",
            batch_size=1
        )(
            delayed(run_chain)(i, data, event)
            for i in range(CHAIN_START, CHAIN_END)
        )

    print("Running rolling one-step holdout predictions ...")
    Parallel(
        n_jobs=CHAIN_END-CHAIN_START,
        backend="loky",
        batch_size=1
    )(
        delayed(predict_chain)(i)
        for i in range(CHAIN_START, CHAIN_END)
    )

    print("All training and prediction done.")

if __name__ == "__main__":
    multiprocessing.set_start_method("spawn", force=True)
    main()

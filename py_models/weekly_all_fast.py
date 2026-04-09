# ================================================================
# WEEKLY BYM (NO COVARIATES, Double PolyGamma)
# p01 + p10 version (FULL DATA, TOP1 optimized)
# ================================================================

import numpy as np
import geopandas as gpd
from tqdm import tqdm
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import coo_matrix, csr_matrix, diags, bmat
from sksparse.cholmod import cholesky
from polyagamma import random_polyagamma
import pyreadr
import pickle
from joblib import Parallel, delayed
from pathlib import Path
import multiprocessing

BASE_DIR = Path(r"D:\77\Research\temp\snow")

burn = 10000
thin = 2
tot_save = 5000
total_iters = burn + thin * tot_save
n_chains = 10
period = 52

# ================================================================
# LOAD DATA
# ================================================================
def load_data():
    snow = pyreadr.read_r(BASE_DIR / "snow_cleaned_full.Rda")
    snow = list(snow.values())[0].reset_index(drop=True)

    coords = snow.iloc[:, :2].to_numpy()
    y = snow.iloc[:, 2:].to_numpy()

    return coords, y

# ================================================================
# BUILD DATASET (FULL)
# ================================================================
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

    # ============================================================
    # X_eta (NO COVARIATES)
    # ============================================================
    eta_dim = K_total*S

    rows_e, cols_e, vals_e = [], [], []

    for i in range(N):
        s = row_idx[i]

        for k in range(K_base):
            xval = cov4[i,k]
            rows_e += [i,i]
            cols_e += [(2*k)*S+s, (2*k+1)*S+s]
            vals_e += [xval,xval]

    X_eta_base = coo_matrix((vals_e,(rows_e,cols_e)),
                            shape=(N,eta_dim)).tocsr()

    # ===== TOP1 OPTIMIZATION =====
    rows_all, cols_all = X_eta_base.nonzero()

    j = cols_all // S
    w = week_idx[rows_all]
    tau_indices = j*52 + w

    # ============================================================
    # X_tau
    # ============================================================
    tau_dim = K_total * 52

    rows_t, cols_t, vals_t = [], [], []

    for i in range(N):
        w = week_idx[i]

        for k in range(K_base):
            xval = cov4[i,k]
            rows_t += [i,i]
            cols_t += [(2*k)*52+w, (2*k+1)*52+w]
            vals_t += [xval,xval]

    X_tau_base = coo_matrix((vals_t,(rows_t,cols_t)),
                            shape=(N,tau_dim)).tocsr()

    # ============================================================
    # PRIOR
    # ============================================================
    Q_blocks = []
    for k in range(K_base):
        Q_blocks.append(Q_car)
        Q_blocks.append(I_S)

    Q_eta = bmat(
        [[Q_blocks[i] if i==j else None for j in range(K_total)]
         for i in range(K_total)],
        format="csr"
    )

    tau_prior_prec = diags(np.ones(tau_dim)/9)

    return (
        X_eta_base, X_tau_base, Q_eta, tau_prior_prec,
        row_idx, week_idx, kappa,
        eta_dim, tau_dim, K_total, S,
        tau_indices
    )

# ================================================================
# MCMC
# ================================================================
def run_chain(chain_id, data, event):

    (
        X_eta_base, X_tau_base, Q_eta, tau_prior_prec,
        row_idx, week_idx, kappa,
        eta_dim, tau_dim, K_total, S,
        tau_indices
    ) = data

    np.random.seed(2000 + chain_id)

    curr_eta = np.zeros(eta_dim)
    curr_tau = np.ones(tau_dim)

    all_eta = np.zeros((eta_dim, tot_save), dtype=np.float32)
    all_tau = np.zeros((tau_dim, tot_save), dtype=np.float32)

    save_idx = 0

    base_rows, base_cols = X_eta_base.nonzero()
    j_idx = base_cols // S

    # ===== INIT CHOLESKY (FAST) =====
    # eta
    X_tilde_eta = X_eta_base.copy()
    tau_idx = j_idx*52 + week_idx[base_rows]
    X_tilde_eta.data *= curr_tau[tau_idx]

    psi = X_tilde_eta @ curr_eta
    omega = random_polyagamma(1, psi)

    XtOmega = X_tilde_eta.T.multiply(omega)
    post_prec_eta = (XtOmega @ X_tilde_eta + Q_eta).tocsc()

    factor_eta = cholesky(post_prec_eta, mode="simplicial")

    # tau
    X_tilde_tau = X_tau_base.copy()
    XtOmega = X_tilde_tau.T.multiply(omega)
    post_prec_tau = (XtOmega @ X_tilde_tau + tau_prior_prec).tocsc()

    factor_tau = cholesky(post_prec_tau, mode="simplicial")
    for it in tqdm(range(total_iters),
                   desc=f"{event}-Chain {chain_id}",
                   position=chain_id):

        # ====================================================
        # psi
        # ====================================================
        X_tilde_eta = X_eta_base.copy()
        tau_idx = j_idx*52 + week_idx[base_rows]
        X_tilde_eta.data *= curr_tau[tau_idx]

        psi = X_tilde_eta @ curr_eta

        # ====================================================
        # PG
        # ====================================================
        omega = random_polyagamma(1, psi)

        # ====================================================
        # ETA update
        # ====================================================
        XtOmega = X_tilde_eta.T.multiply(omega)
        post_prec_eta = (XtOmega @ X_tilde_eta + Q_eta).tocsc()

        factor_eta.cholesky_inplace(post_prec_eta)
        mu = factor_eta.solve_A(X_tilde_eta.T @ kappa)

        z = np.random.randn(eta_dim)
        z = z / np.sqrt(factor_eta.D())
        z = factor_eta.solve_Lt(z)
        z = factor_eta.apply_Pt(z)

        curr_eta = mu + z

        # ====================================================
        # TAU update
        # ====================================================
        X_tilde_tau = X_tau_base.copy()

        rows_nz, cols_nz = X_tilde_tau.nonzero()
        j = cols_nz // 52
        s = row_idx[rows_nz]
        eta_indices = j*S + s

        X_tilde_tau.data *= curr_eta[eta_indices]

        residual = kappa

        XtOmega = X_tilde_tau.T.multiply(omega)
        post_prec_tau = (XtOmega @ X_tilde_tau + tau_prior_prec).tocsc()

        factor_tau.cholesky_inplace(post_prec_tau)
        mu = factor_tau.solve_A(X_tilde_tau.T @ residual)

        z = np.random.randn(tau_dim)
        z = z / np.sqrt(factor_tau.D())
        z = factor_tau.solve_Lt(z)
        z = factor_tau.apply_Pt(z)

        curr_tau = mu + z

        # ====================================================
        # SAVE
        # ====================================================
        if it >= burn and (it-burn)%thin==0:
            all_eta[:,save_idx] = curr_eta
            all_tau[:,save_idx] = curr_tau
            save_idx += 1
            if save_idx == tot_save:
                break

    with open(BASE_DIR / f"{event}_weekly_chain{chain_id}.pkl","wb") as f:
        pickle.dump({"eta":all_eta,"tau":all_tau},f)

    return chain_id

# ================================================================
# MAIN
# ================================================================
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
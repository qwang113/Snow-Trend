import numpy as np
import geopandas as gpd
import pyreadr
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix, coo_matrix, diags, bmat
from scipy.sparse.csgraph import connected_components
from tqdm import tqdm
from pathlib import Path
from multiprocessing import Pool, current_process
from sksparse.cholmod import cholesky
from polyagamma import random_polyagamma
import os

# ============================================================
# CONFIG
# ============================================================
DIST_TH = 0.22
period = 52

burn = 10000
thin = 2
tot_save = 5000
total_iters = burn + tot_save * thin

prior_prec = 1.0 / 25.0
N_CHAINS = 10

BASE_DIR = Path(r"D:\77\Research\temp\snow")

# ============================================================
# 1️⃣ Load data
# ============================================================
snow = pyreadr.read_r("snow_cleaned_full.Rda")
snow = list(snow.values())[0].reset_index(drop=True)

coords_all = snow.iloc[:, :2].to_numpy()
y_all = snow.iloc[:, 2:].to_numpy()

S_full, TT = y_all.shape

# ============================================================
# 2️⃣ USE FULL DATA (no filtering)
# ============================================================

coords = coords_all
y = y_all
S = S_full

print("Using FULL data")
print("S =", S)

# ============================================================
# 3️⃣ global trend
# ============================================================
t_full = np.arange(1, TT + 1)
t_trend_full = (t_full - t_full.mean()) / t_full.std(ddof=0)

def run_chain(args):
    chain_id, event_name = args

    np.random.seed(1234 + chain_id)

    # --------------------------------------------------------
    # select events
    # --------------------------------------------------------
    if event_name == "p01":
        loc_mask = (y[:, :-1] == 0)
        kappa_builder = lambda ny: ny - 0.5
    else:
        loc_mask = (y[:, :-1] == 1)
        kappa_builder = lambda ny: (1 - ny) - 0.5

    loc = np.where(loc_mask)
    pairs = np.column_stack(loc)
    pairs = pairs[np.lexsort((pairs[:,0], pairs[:,1]))]
    pairs[:,1] += 1

    row_idx  = pairs[:,0]
    time_idx = pairs[:,1] - 1
    N = len(row_idx)

    next_y = y[pairs[:,0], pairs[:,1]]
    kappa  = kappa_builder(next_y)

    t_raw   = time_idx + 1
    t_trend = t_trend_full[time_idx]

    covariates = np.column_stack([
        np.ones(N),
        np.cos(2*np.pi*t_raw / period),
        np.sin(2*np.pi*t_raw / period),
        t_trend
    ])

    K = covariates.shape[1]
    theta_dim = K * S

    # --------------------------------------------------------
    # build X
    # --------------------------------------------------------
    rows, cols, vals = [], [], []

    for i in range(N):
        s = row_idx[i]
        for k in range(K):
            rows.append(i)
            cols.append(s + k*S)
            vals.append(covariates[i, k])

    X = coo_matrix((vals,(rows,cols)), shape=(N,theta_dim)).tocsr()

    # --------------------------------------------------------
    # prior
    # --------------------------------------------------------
    blocks = [
        [prior_prec * diags(np.ones(S)) if i == j else None
         for j in range(K)]
        for i in range(K)
    ]
    curr_prec = bmat(blocks, format="csr")

    # --------------------------------------------------------
    # MCMC
    # --------------------------------------------------------
    curr_theta = np.zeros(theta_dim)
    all_theta = np.zeros((theta_dim, tot_save))
    save_idx = 0

    pbar = tqdm(
        range(total_iters),
        desc=f"{event_name} chain {chain_id}",
        position=chain_id,
        leave=True
    )

    for it in pbar:

        phi = X @ curr_theta
        omega = random_polyagamma(1, phi, size=N)

        XtOmega = X.T.multiply(omega)
        post_prec = (XtOmega @ X + curr_prec).tocsc()

        factor = cholesky(post_prec, mode="simplicial")

        rhs = X.T @ kappa
        mu = factor.solve_A(rhs)

        z = np.random.randn(theta_dim)
        z = z / np.sqrt(factor.D())
        z = factor.solve_Lt(z)
        z = factor.apply_Pt(z)

        curr_theta = mu + z

        if it >= burn and (it - burn) % thin == 0:
            all_theta[:, save_idx] = curr_theta
            save_idx += 1
            if save_idx == tot_save:
                break

    save_path = BASE_DIR / f"{event_name}_ind_all_chain{chain_id}.npz"
    np.savez_compressed(save_path, all_theta=all_theta)

    return f"{event_name} chain {chain_id} done"

if __name__ == "__main__":

    for event in ["p01", "p10"]:

        print(f"\nRunning {event} ...")

        args = [(i, event) for i in range(N_CHAINS)]

        with Pool(processes=N_CHAINS) as pool:
            results = pool.map(run_chain, args)

        print(results)

    print("\nAll done.")
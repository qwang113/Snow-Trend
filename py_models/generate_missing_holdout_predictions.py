"""Generate only the missing holdout prediction arrays from saved MCMC chains.

Inputs and outputs use a user-configured local holdout directory. The script
never overwrites an existing valid prediction. A temporary .partial.npy file is atomically renamed after the
expected shape and dtype have been verified.
"""

from __future__ import annotations

import argparse
import gc
import os
import pickle
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr


os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

# Set these paths for the local system.
BASE_DIR = Path("path/to/snow/data")
CHAIN_DIR = Path("path/to/snow/results") / "holdout_38_14"
DEFAULT_OUTPUT_DIR = CHAIN_DIR

PERIOD = 52
N_TRAIN_YEARS = 38
N_TEST_YEARS = 14
TRAIN_WEEKS = N_TRAIN_YEARS * PERIOD
TEST_WEEKS = N_TEST_YEARS * PERIOD
PRED_THIN = 15
N_DRAWS = 334
N_CELLS = 1618
SEED_BASE = 2000
EXPECTED_SHAPE = (N_DRAWS, N_CELLS, TEST_WEEKS)


def load_snow_data():
    snow = pyreadr.read_r(BASE_DIR / "snow_cleaned_full.Rda")
    snow = list(snow.values())[0].reset_index(drop=True)
    coords = snow.iloc[:, :2].to_numpy()
    y_full = snow.iloc[:, 2:].to_numpy()
    expected_weeks = (N_TRAIN_YEARS + N_TEST_YEARS) * PERIOD
    if y_full.shape != (N_CELLS, expected_weeks):
        raise ValueError(f"Unexpected snow shape: {y_full.shape}")
    return coords, y_full


def valid_output(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        arr = np.load(path, mmap_mode="r", allow_pickle=False)
        ok = arr.shape == EXPECTED_SHAPE and arr.dtype == np.float32
        del arr
        return ok
    except Exception:
        return False


def common_time_inputs(y_full):
    y_prev_test = y_full[:, TRAIN_WEEKS - 1 : -1]
    if y_prev_test.shape != (N_CELLS, TEST_WEEKS):
        raise ValueError(f"Unexpected previous-state shape: {y_prev_test.shape}")
    train_t = np.arange(1, TRAIN_WEEKS + 1)
    transition_t = np.arange(TRAIN_WEEKS, TRAIN_WEEKS + TEST_WEEKS)
    transition_scaled = (transition_t - train_t.mean()) / train_t.std(ddof=0)
    return y_prev_test, transition_t, transition_scaled


def load_bym_chains(chain_id: int):
    with open(CHAIN_DIR / f"p01_weekly_train38_chain{chain_id}.pkl", "rb") as f:
        d01 = pickle.load(f)
    with open(CHAIN_DIR / f"p10_weekly_train38_chain{chain_id}.pkl", "rb") as f:
        d10 = pickle.load(f)

    eta01 = d01["eta"][:, ::PRED_THIN].reshape(8, N_CELLS, N_DRAWS)
    tau01 = d01["tau"][:, ::PRED_THIN].reshape(8, PERIOD, N_DRAWS)
    eta10 = d10["eta"][:, ::PRED_THIN].reshape(8, N_CELLS, N_DRAWS)
    tau10 = d10["tau"][:, ::PRED_THIN].reshape(8, PERIOD, N_DRAWS)
    return eta01, tau01, eta10, tau10


def load_bym_plus_chains(chain_id: int):
    with open(
        CHAIN_DIR / f"p01_weekly_cov+lon_train38_chain{chain_id}.pkl", "rb"
    ) as f:
        d01 = pickle.load(f)
    with open(
        CHAIN_DIR / f"p10_weekly_cov+lon_train38_chain{chain_id}.pkl", "rb"
    ) as f:
        d10 = pickle.load(f)

    eta01_all = d01["eta"][:, ::PRED_THIN]
    eta10_all = d10["eta"][:, ::PRED_THIN]
    if eta01_all.shape[1] != N_DRAWS or eta10_all.shape[1] != N_DRAWS:
        raise ValueError("Unexpected number of posterior draws")

    gamma01 = eta01_all[8 * N_CELLS :, :]
    gamma10 = eta10_all[8 * N_CELLS :, :]
    if gamma01.shape != (5, N_DRAWS) or gamma10.shape != (5, N_DRAWS):
        raise ValueError(f"Unexpected gamma shapes: {gamma01.shape}, {gamma10.shape}")

    eta01 = eta01_all[: 8 * N_CELLS, :].reshape(8, N_CELLS, N_DRAWS)
    eta10 = eta10_all[: 8 * N_CELLS, :].reshape(8, N_CELLS, N_DRAWS)
    tau01 = d01["tau"][:, ::PRED_THIN].reshape(8, PERIOD, N_DRAWS)
    tau10 = d10["tau"][:, ::PRED_THIN].reshape(8, PERIOD, N_DRAWS)
    return eta01, tau01, gamma01, eta10, tau10, gamma10


def prediction_paths(model: str, chain_id: int, output_dir: Path):
    if model == "bym":
        name = f"pred_weekly_bym_train38_test14_chain{chain_id}.npy"
    elif model == "bym_plus":
        name = f"pred_weekly_bym_lon_train38_test14_chain{chain_id}.npy"
    else:
        raise ValueError(f"Unknown model: {model}")
    final = output_dir / name
    partial = output_dir / f"{name}.partial.npy"
    return final, partial


def generate_one(model: str, chain_id: int, output_dir_text: str):
    # Matches the original chain-specific seed convention. Prediction is
    # deterministic given the saved chains, but setting this explicitly keeps
    # the execution configuration aligned with the fitting scripts.
    np.random.seed(SEED_BASE + chain_id)
    output_dir = Path(output_dir_text)
    output_dir.mkdir(parents=True, exist_ok=True)
    final, partial = prediction_paths(model, chain_id, output_dir)

    if valid_output(final):
        return model, chain_id, "skipped-valid", str(final)
    if final.exists():
        raise ValueError(f"Existing output is invalid; refusing to overwrite: {final}")
    if partial.exists():
        partial.unlink()

    coords, y_full = load_snow_data()
    y_prev_test, transition_t, transition_scaled = common_time_inputs(y_full)

    if model == "bym":
        eta01, tau01, eta10, tau10 = load_bym_chains(chain_id)
        gamma01 = gamma10 = covariates = None
    else:
        eta01, tau01, gamma01, eta10, tau10, gamma10 = load_bym_plus_chains(
            chain_id
        )

        lon_raw = coords[:, 0]
        mask_na = lon_raw < -30
        mask_euas = ~mask_na
        lon_na = np.zeros(N_CELLS)
        lon_euas = np.zeros(N_CELLS)
        lon_na[mask_na] = (lon_raw[mask_na] - lon_raw[mask_na].mean()) / max(
            lon_raw[mask_na].std(), 1e-6
        )
        lon_euas[mask_euas] = (
            lon_raw[mask_euas] - lon_raw[mask_euas].mean()
        ) / max(lon_raw[mask_euas].std(), 1e-6)
        lat = (coords[:, 1] - coords[:, 1].mean()) / coords[:, 1].std()

        no_nbs = np.array(
            [
                57,
                170,
                236,
                269,
                343,
                685,
                946,
                947,
                989,
                1037,
                1084,
                1090,
                1109,
                1118,
                1127,
                1176,
                1203,
            ]
        ) - 1
        elev_raw = pd.read_csv(BASE_DIR / "curr_elev.csv").iloc[:, 3].to_numpy()
        nnbs_elev = (
            pd.read_csv(BASE_DIR / "nnbs_elev.csv", sep="\t").iloc[:, 2].to_numpy()
        )
        keep = np.ones(N_CELLS, dtype=bool)
        keep[no_nbs] = False
        elev_all = np.zeros(N_CELLS)
        elev_all[keep] = elev_raw
        elev_all[no_nbs] = nnbs_elev
        elev = (elev_all - elev_all.mean()) / elev_all.std()

        snow_temp = pyreadr.read_r(BASE_DIR / "snow_temp_full.Rda")
        temp_full = list(snow_temp.values())[0].iloc[:, 2:].to_numpy()
        temp_train = temp_full[:, :TRAIN_WEEKS]
        temp_scaled = (temp_full - temp_train.mean()) / temp_train.std()
        covariates = (lon_na, lon_euas, lat, elev, temp_scaled)

    p_snow = np.lib.format.open_memmap(
        partial, mode="w+", dtype=np.float32, shape=EXPECTED_SHAPE
    )

    for j, (t_raw, t_scaled) in enumerate(zip(transition_t, transition_scaled)):
        week = (t_raw - 1) % PERIOD
        cov = np.array(
            [
                1.0,
                1.0,
                np.cos(2 * np.pi * t_raw / PERIOD),
                np.cos(2 * np.pi * t_raw / PERIOD),
                np.sin(2 * np.pi * t_raw / PERIOD),
                np.sin(2 * np.pi * t_raw / PERIOD),
                t_scaled,
                t_scaled,
            ]
        )
        phi01 = np.sum(
            cov[:, None, None] * eta01 * tau01[:, week, None, :], axis=0
        )
        phi10 = np.sum(
            cov[:, None, None] * eta10 * tau10[:, week, None, :], axis=0
        )

        if model == "bym_plus":
            lon_na, lon_euas, lat, elev, temp_scaled = covariates
            z = (lon_na, lon_euas, lat, elev, temp_scaled[:, TRAIN_WEEKS + j])
            phi01 += t_scaled * sum(z[k][:, None] * gamma01[k] for k in range(5))
            phi10 += t_scaled * sum(z[k][:, None] * gamma10[k] for k in range(5))

        p01 = 1.0 / (1.0 + np.exp(-phi01))
        p10 = 1.0 / (1.0 + np.exp(-phi10))
        p_snow[:, :, j] = np.where(
            y_prev_test[:, j, None] == 0, p01, 1.0 - p10
        ).T.astype(np.float32)

    p_snow.flush()
    del p_snow
    gc.collect()

    if not valid_output(partial):
        raise ValueError(f"Generated output failed validation: {partial}")
    if valid_output(final):
        partial.unlink()
        return model, chain_id, "skipped-arrived-during-run", str(final)
    if final.exists():
        raise ValueError(f"A conflicting invalid output appeared during generation: {final}")
    os.replace(partial, final)
    return model, chain_id, "generated", str(final)


def required_jobs(output_dir: Path):
    jobs = [("bym", chain) for chain in range(10)]
    jobs.extend(("bym_plus", chain) for chain in range(7, 10))
    return [
        job
        for job in jobs
        if not valid_output(prediction_paths(job[0], job[1], output_dir)[0])
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--list-only", action="store_true")
    args = parser.parse_args()

    jobs = required_jobs(args.output_dir)
    print(f"Output directory: {args.output_dir}")
    print(f"Jobs to generate ({len(jobs)}): {jobs}")
    if args.list_only or not jobs:
        return

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(generate_one, model, chain, str(args.output_dir)): (model, chain)
            for model, chain in jobs
        }
        for future in as_completed(futures):
            model, chain = futures[future]
            try:
                print("DONE", future.result(), flush=True)
            except Exception as exc:
                print(f"FAILED {model} chain {chain}: {exc}", flush=True)
                raise


if __name__ == "__main__":
    main()

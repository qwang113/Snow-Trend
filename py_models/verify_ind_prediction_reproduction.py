"""Regenerate one IND holdout prediction and compare it with the downloaded file."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pyreadr


# Set these paths for the local system.
BASE_DIR = Path("path/to/snow/data")
RESULT_DIR = Path("path/to/snow/results") / "holdout_38_14"
CHAIN_ID = 0
SEED_BASE = 1234
PERIOD = 52
TRAIN_WEEKS = 38 * PERIOD
TEST_WEEKS = 14 * PERIOD
PRED_THIN = 15
EXPECTED_SHAPE = (334, 1618, TEST_WEEKS)

REFERENCE = RESULT_DIR / f"pred_ind_train38_test14_chain{CHAIN_ID}.npy"
REGENERATED = RESULT_DIR / f"pred_ind_train38_test14_chain{CHAIN_ID}.verify.npy"


def sha256(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def regenerate():
    # Matches the original IND chain seed convention. Prediction itself is
    # deterministic given the saved posterior chain arrays.
    np.random.seed(SEED_BASE + CHAIN_ID)

    snow = pyreadr.read_r(BASE_DIR / "snow_cleaned_full.Rda")
    snow = list(snow.values())[0].reset_index(drop=True)
    y_full = snow.iloc[:, 2:].to_numpy()
    s, total_weeks = y_full.shape
    if (s, total_weeks) != (1618, 52 * PERIOD):
        raise ValueError(f"Unexpected snow data shape: {y_full.shape}")

    with np.load(
        RESULT_DIR / f"p01_ind_train38_chain{CHAIN_ID}.npz", allow_pickle=False
    ) as d01:
        theta01 = d01["all_theta"][:, ::PRED_THIN]
    with np.load(
        RESULT_DIR / f"p10_ind_train38_chain{CHAIN_ID}.npz", allow_pickle=False
    ) as d10:
        theta10 = d10["all_theta"][:, ::PRED_THIN]

    if theta01.shape[1] != 334 or theta10.shape[1] != 334:
        raise ValueError(f"Unexpected retained draws: {theta01.shape}, {theta10.shape}")
    theta01 = theta01.reshape(4, s, 334)
    theta10 = theta10.reshape(4, s, 334)

    y_prev_test = y_full[:, TRAIN_WEEKS - 1 : -1]
    if y_prev_test.shape != (s, TEST_WEEKS):
        raise ValueError(f"Unexpected previous-state shape: {y_prev_test.shape}")

    t_full = np.arange(1, TRAIN_WEEKS + 1)
    transition_t = np.arange(TRAIN_WEEKS, TRAIN_WEEKS + TEST_WEEKS)
    transition_scaled = (transition_t - t_full.mean()) / t_full.std(ddof=0)

    output = np.lib.format.open_memmap(
        REGENERATED, mode="w+", dtype=np.float32, shape=EXPECTED_SHAPE
    )
    for j, (t_raw, t_scaled) in enumerate(zip(transition_t, transition_scaled)):
        x = np.array(
            [
                1.0,
                np.cos(2 * np.pi * t_raw / PERIOD),
                np.sin(2 * np.pi * t_raw / PERIOD),
                t_scaled,
            ]
        )
        phi01 = np.sum(x[:, None, None] * theta01, axis=0)
        phi10 = np.sum(x[:, None, None] * theta10, axis=0)
        p01 = 1.0 / (1.0 + np.exp(-phi01))
        p10 = 1.0 / (1.0 + np.exp(-phi10))
        output[:, :, j] = np.where(
            y_prev_test[:, j, None] == 0, p01, 1.0 - p10
        ).T.astype(np.float32)
    output.flush()
    del output


def compare():
    reference = np.load(REFERENCE, mmap_mode="r", allow_pickle=False)
    regenerated = np.load(REGENERATED, mmap_mode="r", allow_pickle=False)
    if reference.shape != regenerated.shape or reference.dtype != regenerated.dtype:
        raise ValueError(
            f"Metadata mismatch: reference={reference.shape}/{reference.dtype}, "
            f"regenerated={regenerated.shape}/{regenerated.dtype}"
        )

    max_abs = 0.0
    abs_sum = 0.0
    exact_count = 0
    close_count = 0
    total = reference.size
    for start in range(0, reference.shape[2], 26):
        stop = min(start + 26, reference.shape[2])
        left = np.asarray(reference[:, :, start:stop])
        right = np.asarray(regenerated[:, :, start:stop])
        difference = np.abs(left - right)
        max_abs = max(max_abs, float(difference.max()))
        abs_sum += float(difference.sum(dtype=np.float64))
        exact_count += int(np.count_nonzero(left == right))
        close_count += int(np.count_nonzero(np.isclose(left, right, rtol=1e-6, atol=1e-7)))

    ref_hash = sha256(REFERENCE)
    regen_hash = sha256(REGENERATED)
    print(f"reference:   {REFERENCE}")
    print(f"regenerated: {REGENERATED}")
    print(f"shape: {reference.shape}")
    print(f"dtype: {reference.dtype}")
    print(f"max_abs_error: {max_abs:.12g}")
    print(f"mean_abs_error: {abs_sum / total:.12g}")
    print(f"exact_equal_fraction: {exact_count / total:.12%}")
    print(f"allclose_fraction: {close_count / total:.12%}")
    print(f"reference_sha256:   {ref_hash}")
    print(f"regenerated_sha256: {regen_hash}")
    print(f"byte_identical: {ref_hash == regen_hash}")


def main():
    if REGENERATED.exists():
        REGENERATED.unlink()
    try:
        regenerate()
        compare()
    finally:
        if REGENERATED.exists():
            REGENERATED.unlink()


if __name__ == "__main__":
    main()

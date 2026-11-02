"""Axisymmetric cells and bounded-memory initial position generation."""
import math
import numpy as np

CELL_DTYPE = np.dtype([
    ("cell_id", "<u8"), ("r_min_mm", "<f8"), ("r_max_mm", "<f8"),
    ("z_min_mm", "<f8"), ("z_max_mm", "<f8"),
    ("r_center_mm", "<f8"), ("z_center_mm", "<f8"),
])
INITIAL_DTYPE = np.dtype([
    ("event", "<u4"), ("cell_id", "<u8"),
    ("x0", "<f8"), ("y0", "<f8"), ("z0", "<f8"),
])

def edges(low, high, size):
    ratio = (high - low) / size
    nearest = round(ratio)
    n = max(1, nearest if math.isclose(ratio, nearest, rel_tol=1e-12, abs_tol=1e-12) else math.ceil(ratio))
    if n > 10_000_000:
        raise ValueError("La rejilla tiene demasiadas divisiones.")
    values = low + np.arange(n + 1) * size
    values[-1] = high
    return values

def cells(config):
    r = edges(*config["bounds"]["r"], config["grid"]["size_mm"])
    z = edges(*config["bounds"]["z"], config["grid"]["size_mm"])
    total = (len(r)-1) * (len(z)-1)
    if total * config["grid"]["electrons_per_cell"] > 2**32:
        raise ValueError("El formato event de SENTINEL admite como máximo 2**32 electrones.")
    table = np.empty(total, dtype=CELL_DTYPE)
    i = 0
    for r0, r1 in zip(r[:-1], r[1:]):
        for z0, z1 in zip(z[:-1], z[1:]):
            table[i] = (i, r0, r1, z0, z1, (r0+r1)/2, (z0+z1)/2)
            i += 1
    return table

def initial_batches(table, per_cell, rng=None, batch_size=8192):
    rng = rng if rng is not None else np.random.default_rng()
    total = len(table) * per_cell
    for start in range(0, total, batch_size):
        stop = min(total, start + batch_size)
        event = np.arange(start, stop, dtype=np.uint64)
        cell = table[event // per_cell]
        positions = np.zeros(stop-start, dtype=INITIAL_DTYPE)
        positions["event"] = event
        positions["cell_id"] = cell["cell_id"]
        positions["x0"] = np.sqrt(cell["r_min_mm"]**2 + rng.random(stop-start)
                                   * (cell["r_max_mm"]**2 - cell["r_min_mm"]**2))
        positions["z0"] = cell["z_min_mm"] + rng.random(stop-start) * (cell["z_max_mm"]-cell["z_min_mm"])
        yield positions

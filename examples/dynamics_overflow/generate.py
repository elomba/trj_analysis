#!/usr/bin/env python3
"""Create a small, reproducible LAMMPS-style NetCDF dynamics fixture."""

import pathlib
import numpy as np
from netCDF4 import Dataset

out = pathlib.Path(__file__).resolve().parent / "run.nc"
rng = np.random.default_rng(20261108)
nframe, natom = 700, 256
side = np.float32(12.86)

with Dataset(out, "w", format="NETCDF4_CLASSIC") as nc:
    nc.createDimension("frame", None)
    nc.createDimension("atom", natom)
    nc.createDimension("spatial", 3)
    nc.createDimension("cell_spatial", 3)
    nc.createDimension("cell_angular", 3)
    nc.createDimension("label", 10)
    nc.Conventions = "AMBER"
    nc.ConventionVersion = "1.0"
    nc.program = "LAMMPS"

    for name, dimensions, values in (
        ("spatial", ("spatial",), b"xyz"),
        ("cell_spatial", ("cell_spatial",), b"xyz"),
    ):
        var = nc.createVariable(name, "S1", dimensions)
        var[:] = np.frombuffer(values, dtype="S1")
    angles = nc.createVariable("cell_angular", "S1", ("cell_angular", "label"))
    angles[:] = np.array([list(s.ljust(10)) for s in ("alpha", "beta", "gamma")], dtype="S1")

    time = nc.createVariable("time", "f4", ("frame",))
    time.units = "lj"
    time.scale_factor = np.float32(0.01)
    time.set_auto_maskandscale(False)
    origin = nc.createVariable("cell_origin", "f4", ("frame", "cell_spatial"))
    lengths = nc.createVariable("cell_lengths", "f4", ("frame", "cell_spatial"))
    cell_angles = nc.createVariable("cell_angles", "f4", ("frame", "cell_angular"))
    ids = nc.createVariable("id", "i4", ("frame", "atom"))
    types = nc.createVariable("type", "i4", ("frame", "atom"))
    coords = nc.createVariable("coordinates", "f4", ("frame", "atom", "spatial"))
    coords.units = "lj"
    velocities = nc.createVariable("velocities", "f4", ("frame", "atom", "spatial"))
    velocities.units = "lj"

    particle_ids = np.arange(1, natom + 1, dtype=np.int32)
    particle_types = np.where(particle_ids <= natom // 2, 1, 2).astype(np.int32)
    for frame in range(nframe):
        time[frame] = frame + 1
        origin[frame] = 0
        lengths[frame] = side
        cell_angles[frame] = 90
        ids[frame] = particle_ids
        types[frame] = particle_types
        coords[frame] = rng.uniform(0, side, (natom, 3)).astype("f4")
        v = rng.normal(0, np.sqrt(1.5), (natom, 3)).astype("f4")
        v -= v.mean(axis=0)
        v *= np.float32(np.sqrt(1.5 / np.mean(v * v)))
        velocities[frame] = v

print(out)

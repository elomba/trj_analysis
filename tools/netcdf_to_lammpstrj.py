#!/usr/bin/env python3
import sys
import numpy as np
from netCDF4 import Dataset


def write_lammpstrj(ncfile, outfile, start=0, end=None, stride=1):

    nc = Dataset(ncfile, "r")

    coords = nc.variables["coordinates"]
    times  = nc.variables["time"]
    scale = times.getncattr("scale_factor")
    units = times.getncattr("units")

    print("scale_factor =", scale)
    print("units =", units)
    if "cell_lengths" in nc.variables:
        cell_lengths = nc.variables["cell_lengths"]
    else:
        cell_lengths = None

    natoms = coords.shape[1]
    nframes = coords.shape[0]

    if end is None or end > nframes:
        end = nframes

    # atom ids
    if "atom_id" in nc.variables:
        atom_ids = nc.variables["atom_id"][:]
    else:
        atom_ids = np.arange(1, natoms + 1)

    # atom types
    if "atom_type" in nc.variables:
        atom_types = nc.variables["atom_type"][:]
    else:
        atom_types = np.ones(natoms, dtype=int)

    with open(outfile, "w") as f:

        for t in range(start, end, stride):

            xyz = coords[t]
            time = int(times[t]/scale)


            f.write("ITEM: TIMESTEP\n")
            f.write(f"{time}\n")

            f.write("ITEM: NUMBER OF ATOMS\n")
            f.write(f"{natoms}\n")

            if cell_lengths is not None:
                Lx, Ly, Lz = cell_lengths[t]
                xlo, xhi = 0.0, Lx
                ylo, yhi = 0.0, Ly
                zlo, zhi = 0.0, Lz
            else:
                xlo = ylo = zlo = 0.0
                xhi = yhi = zhi = 1.0

            f.write("ITEM: BOX BOUNDS pp pp pp\n")
            f.write(f"{xlo} {xhi}\n")
            f.write(f"{ylo} {yhi}\n")
            f.write(f"{zlo} {zhi}\n")

            f.write("ITEM: ATOMS id type x y z\n")

            for i in range(natoms):
                x, y, z = xyz[i]
                f.write(f"{atom_ids[i]} {atom_types[i]} {x} {y} {z}\n")

    nc.close()


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Usage:")
        print("python netcdf_to_lammpstrj.py input.nc output.lammpstrj [start end stride]")
        sys.exit()

    ncfile = sys.argv[1]
    outfile = sys.argv[2]

    start = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    end   = int(sys.argv[4]) if len(sys.argv) > 4 else None
    stride= int(sys.argv[5]) if len(sys.argv) > 5 else 1

    write_lammpstrj(ncfile, outfile, start, end, stride)


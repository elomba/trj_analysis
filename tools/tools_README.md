# Tools

This directory contains **utility scripts for trajectory manipulation and conversion**.

Typical utilities include:

- conversion of LAMMPS dump files
- preprocessing trajectories
- preparing trajectories for analysis

## Example: LAMMPS dump → NetCDF

```
python lammpstrj_to_lammps_netcdf.py trajectory.lammpstrj trajectory.nc
```

Dependencies:

- Python
- NumPy
- netCDF4

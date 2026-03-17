# Tools

This directory contains **utility scripts for trajectory manipulation,conversion, calculation of structural properties**.

Typical utilities include:

- conversion of LAMMPS dump files (script names self explanatory)
    - LAMMPS trj to pdb
    - LAMMPS trj to netcdf
    - LAMMPS netcdf to trj 
- preprocessing trajectories
    - Extract frames from a netcdf
- Structural properties 
  - **fast_grfreud_netcdf.py** g(r) using Freud python package (https://freud.readthedocs.io/en/latest/)
  - **fast_sq_freud_netcdf.py** S(Q) (using Q-space sampling) with Freud package
  - **rerun_grSALRlj.lmp** LAMMPS rerun script to calculate g(r) from a lammpstrj dump
  - **rerun_grSALRlj_ord.lmp** LAMMPS rerun script to calculate $Q_l$ bond orientational order parameters

## Example: LAMMPS dump → NetCDF

```
python lammpstrj_to_lammps_netcdf.py trajectory.lammpstrj trajectory.nc
```

Dependencies:

- Python
- NumPy
- netCDF4
- freud

## Name: 
Trajectory analysis

## Description
 This program performs an analysis of a netcdf trajectory file from LAMMPS.

 At this stage it computes:

- Pair distribution functions
- Particle&cluster fluctuation analysis (for hyperuniformity tests)
- Static & dynamic structure factors
- Cluster analysis (average cluster profiles, cluster-cluster rdf's and sq's) size and
  radii distributions, and a trajectory file with the evolution of cluster  com's 
- Dynamics (position and velocity time correlation functions)
- Kinetic energy (if velocities present in trajectory file) global and per cluster
- Potential energy (if compute pe/peratom pressent in the trajectory file)
- Stress tensor (pressure)  (if compute stress/atom pressent in the trajectory file)

-   **Important notice:** in LAMMPS script the following computes must be included in the dump
                     in order to compute potential energies and pressures
             
      compute stress all stress/atom NULL

      compute ener all pe/atom
      
      dump trj1 all netcdf ${Ndump} run.nc  id type x y z vx vy vz c_stress[*] c_ener

-    In the first INPUT namelist optional character variables "ener_name" and "press_name" 
     refer to the  names of the computes

*Restrictions*

    - The trajectory file **MUST** us NETDCF format 
    - Only orthogonal simulation cells
    - The number of particles MUST be constant (NpT simulations are allowed, but structure
      factors will be affected by minor errors due to changes in \Delta Q = 2\pi/L)
    - The mimimum dump information to process is 
      dump trj1 all netcdf ${Ndump} run.nc  id type x y z      

## Input 
The input is provided as a set of namelist directives (see attached example) 

## Installation
A Makefile is included (-mno-avx512f can be removed from compilation options if AVX512 instruction set present in the CPU)

Requires NVIDIA CUDA SDK >= 11.6, netcdf v 4.9, FFTW3.

## Authors and acknowledgment
 A. Diaz-Pozuelo & E. Lomba, CSIC-Madrid/USC Santiago de Compostela, October 2024 

## License
Not licensed yet.
Creative Commons Non commercial CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/)
## Project status
Version 0.2.3 ready. Awaiting for publication


## Name: 
Trajectory analysis

## Description
 This program performs an analysis of a netcdf trajectory file from LAMMPS.
  At this stage it computes:

- Pair distribution functions
- Particle&cluster fluctuation analysis (for hyperuniformity tests)
- Static & dynamic structure factors
- Cluster analysis (average cluster profiles, cluster-cluster rdf's and sq's) size and radii distributions, and a trajectory file with the evolution of cluster  com's 
- Dynamics (position and velocity time correlation functions)
- Kinetic energy (if velocities present in trajectory file) global and per cluster
- Potential energy. Only when potential table files LAMMPS style (using RSQ tabulation) are available as 
  ulm.dat (l,m denote the species)

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
Version 0.2 ready. Awaiting for publication


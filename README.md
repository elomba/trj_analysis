## Name: 
Trajectoryanalysis

## Description
 This program performs an analysis of a netcdf trajectory file from LAMMPS.
  At this stage it computes:

- Pair distribution functions
- Static structure factors
- Cluster analysis (average cluster profiles, cluster-cluster rdf's and sq's) size and radii distributions, and a trajectory file with the evolution of cluster  com's 
- Dynamics (position and velocity correlation functions)

## Input 
The input is provided as a set of namelist directives (see attached example) 

## Installation

Requires NVIDIA CUDA SDK >= 11.0, netcdf v 4.9, FFTW3, source code distribution includes a Makefile

## Authors and acknowledgment
 A. Diaz-Pozuelo & E. Lomba, Madrid February 2024 

## License
Not licensed yet.

## Project status
Under current developement


## Name: 
Trajectory analysis

## Description
 This program performs an analysis of a netcdf trajectory file from LAMMPS.

 At this stage it computes:

- Pair distribution functions
- Particle&cluster fluctuation analysis (for hyperuniformity tests)
- Static & dynamic structure factors
- Orientational order parameters, averaged, per particles, and profiles across 
  clusters
- Cluster analysis (average cluster profiles, cluster-cluster rdf's and sq's) size and
  radii distributions, and a trajectory file with the evolution of cluster  com's 
- Dynamics (position and velocity time correlation functions, viscosity if stress tensor
- is computed)
- Kinetic energy (if velocities present in trajectory file) global and per cluster
- Potential energy (if compute pe/peratom pressent in the trajectory file)
- Stress tensor (pressure)  (if compute stress/atom pressent in the trajectory file)
- Particle density and charge density profiles for systems confined along a space direction

-   **Important notice:** in LAMMPS script the following computes must be 
    included in the dump in order to compute potential energies and pressures
      
```
      compute stress all stress/atom NULL
      compute ener all pe/atom
      dump trj1 all netcdf ${Ndump} run.nc  id type x y z vx vy vz q c_stress[*] c_ener
```


-    In the first INPUT namelist optional character variables
     "ener_name" and "press_name" refer to the  names of the computes

*Restrictions*

    - The trajectory file **MUST** us NETDCF format 
    - Only orthogonal simulation cells
    - The number of particles MUST be constant (NpT simulations are allowed, 
      but structure factors will be affected by minor errors due to changes in \Delta Q = 2\pi/L)
    - Molecules are analyzed in terms of their atoms, so far no internal
      degrees of freedom 
      taken into account
    - LAMMPS units must be "real" or "lj" (unit conversion soon to be implemented)
    - The mimimum dump information to process is 
      `dump trj1 all netcdf ${Ndump} run.nc  id type x y z`   

## Usage   
        trj_analysis.exe input.nml (input file with sequence of namelists)
## Input 
!
    ! namelist /INPUT/ log_output_file, trj_input_file, ndim, nsp, norder,nthread, &
    !       ncfs_from_to, rcl, rdf_sq_cl_dyn_sqw_conf_ord, nqw, ener_name, press_name, 
    !       potnbins, potengmargin, periodic(ndim)
    !       Name of log file, name of netcdf trajectory file, no. of dimensions (2,3), no. of species,
    !       no. of CUDA threads (default 128), no. of configurations-start-end, modules to run 
    !       NN distance for cluster analysis and order parameter computation
    !       (logical vars: RDF, structure factor, cluster analysis, dynamics, dynamic S(q,w),
    !       no. of Q's for dynamic analysis, name of compute for potential energy/atom in LAMMPS,
    !       name of stress/atom, no. of bins for energy histograms (def. 100), extra margins 
    !       in energy histograms (def. 0))
    !       periodic(ndim): defaults to true, logical var to indicate non-periodic dimension
    ! namelist /INPUT_SP/ sp_types_selected, sp_labels, mat
    !          IDs of selected species (if nsp<ntypes in trajectory), character labels, atomic mass
    ! namelist /INPUT_RDF/ deltar, rcrdf, nrandom
    !          grid in RDF calculation, cut-off for rdf (default half box size), no. os random origins 
    !          for calculation of local number fluctuations
    ! namelist /INPUT_SQ/ qmax, qmin, bsc
    !          Max value of Q for S(Q), max value for full calculations (all Qs for 0<Q<=qmin), 
    !          scattering lengths 
    ! namelist /INPUT_CL/ dcl, jmin, minclsize, ndrclus
    !          Geometric clustering distance, grid for cluster distribution, minimum cluster size for analysis, no. of bins for cluster profiles
    ! namelist /INPUT_ORD/ orderp, print_orderp, rclcl, nnbond
    !     Orientational order parameters  to be computed
    !    (def. .true.), print average per particle order parameters 
    !    (def. .false.), distance between clusters, for neighbor search
    !    in order parameter calculation  (def. rcl), if nnbond=0 (default)
    !    all NN are used up to the cutoff to compute, otherwise the number of 
    !    neighbors must be the same as the order of the parameter 
    !
    !          
    ! namelist /INPUT_CONF/ idir, pwall, pwallp
    !          Direction of confinement (1,2,3->x,y,z), position of left wall, position of right wall
    ! namelist /INPUT_DYN/ nbuffer, tmax, tmaxp, tlimit, jump
    !           Number of buffers (time origins) for dynamic correlation analysis
    !           buffers are separated by jump configurations (def. 1), maximum time for 
    !           correlation functions (at tmax a window function is applied for FFTs, if omitted tmax=tlimit), if omitted all
    !           t values are used. For viscosity a window function is used at tmaxp (if omitted =tmax)
    !           Averages stored over tlimit only (in ps), and then origin for calculation is shifted if tlimit is
    !           omitted is calculated as a function of the trajectory length, no. buffers, and jump.
    ! namelist /INPUT_SQW/ qw, tmqw
    !          values of Q to compute F(Q,t), Fs(Q,t) and S(Q,w),Ss(Q,w) maximum times for F(Q,t),
    !          if omitted tmax is used
    !
    !
    !    OUTPUT FILES:
    !      * Thermodynamics
    !      - thermo_run.dat (Instantaneous values of thermod. quantities)
    !      * Dynamics
    !      - dyn.dat (msd, <v(t)v(0)>)
    !      - dynw.dat w, Z(w)
    !      - fkt.dat (F(Q_i,t) for nqw Qs)
    !      - fskt.dat (F_self(Q_i,t) for nqw Qs)
    !      - sqw.dat   (S(q,w), S_self(Q,w))
    !      - viscor.dat (<p_xy(t)p_xy(0)  eta(t) (shear viscosity integral))
    !
    !      * Structure
    !      - gmixsim.dat (g_ab(r), g_clcl(r) in cluster analysis on)
    !      - sq.dat  S_NN(Q), (S_cc(Q) , S_11, S_12, S_22 in binary systems)
    !      - sqcl.dat Cluster-cluster S(Q)
    !      - sqmix.dat S_ii (i<=nsp)
    !      - sqw.dat  (S(Q_i,w), S_self(Q_i,w) for nwq Qs)
    !      - order.dat Average Steinhardt order parameters <Q_l> (l=1,norder)
    !      - order_per_mol.dat Average per particle Steinhardt order parameters <q_l>
    !      - last_clconf.dat Last configuration with cluster particle IDs (defined by size)
    !
    !      * Cluster analysis
    !      - rhoprof.dat Average cluster density profile (only for finite clusters) 
    !      - radii.dat Distribution of cluster gyration radii
    !      - clustdistr.dat Distribution of cluster particle size
    !      - distUcl_N.dat Distribution of cluster internal energies per particle
    !      - distUcltot.dat Distribution of cluster internal energies (total)
    !      - clusevol.dat , conf no., no. of clusters, % of particles in clusters
    !      - fshape.dat, cluster shape (deviation from sphere and cylnder)
    !                    (0,0) perfect sphere, (>0,0) perfect cylindero
   !      - ordprof_clust.dat Steinhardt order parameter density profiles (within clusters): averaged over all clusters across slices
    !      - ordprof_clcum.dat Steinhardt order parameter density profiles (cummulative)
    !      - order_per_cl.dat Average Steinhardt order parameters per cluster <Q_l> (l=1,norder)
    !      - centers.lammpstrj trajectory of clusters centers of mass (to be visualized with Ovito)
    !                          Particle no. not constant along the trajectory !!
    ! 
    
## Installation
A Makefile is included (-mno-avx512f can be removed from compilation options if AVX512 instruction set present in the CPU)

Requires NVIDIA CUDA SDK >= 11.6, netcdf v 4.9, FFTW3.

## Authors and acknowledgment
 A. Diaz-Pozuelo & E. Lomba, CSIC-Madrid/USC Santiago de Compostela, August 2025 

## License
Creative Commons Non commercial CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/)
## Project status
Version 0.6 ready. Awaiting for publication


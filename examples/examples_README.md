
# Examples

This directory contains **reference examples illustrating the use of the `trj_analysis` package** for the analysis of molecular simulation trajectories.

The examples cover several types of systems commonly encountered in soft-matter and molecular simulations:

- **SALR systems** (short-range attraction / long-range repulsion) exhibiting cluster formation
- **Binary Lennard-Jones mixtures**
- **Confined systems** (e.g. ionic-liquid-like systems near walls)

Each example includes:

- simulation input files (typically for **LAMMPS**)
- analysis configuration files (`input.nml`)
- example output files produced by the analysis routines

These examples can be used to:

- verify the correct installation of the code
- understand the workflow of the analysis tools
- reproduce typical structural and dynamical observables.

---

# Directory structure

```
examples
│
├── 2D
│   └── salr
│
└── 3D
    ├── salr
    ├── LJmix
    └── Conf
```

The directory is divided according to **system dimensionality and model type**.

---

# Typical workflow

```
LAMMPS simulation
        │
        │ trajectory file
        ▼
   trj_analysis
        │
        ├── structural analysis
        │       g(r), S(q)
        │
        ├── cluster analysis
        │       cluster size distribution
        │       cluster radii
        │       cluster shapes
        │
        ├── order parameters
        │
        └── dynamic observables
                F(k,t)
                dynamical correlations
```

Each example directory contains the necessary files to perform this analysis chain.

---

# 2D systems

## 2D/salr

Example of a **two-dimensional SALR system** exhibiting cluster formation.

### Simulation input

| File | Description |
|-----|-------------|
| `in.pysalr` | LAMMPS input script for the SALR system |
| `py_pot.py` | Python implementation of the SALR potential |
| `input.nml` | Analysis control file |

### Structural observables

| File | Description |
|-----|-------------|
| `grlmp.dat` | Radial distribution function |
| `sq.dat` | Static structure factor |
| `sqmix.dat` | Partial structure factor |
| `sqcl.dat` | Structure factor of clusters |
| `gmixsim.dat` | Cross pair correlation |

### Cluster analysis

| File | Description |
|-----|-------------|
| `clustdistr.dat` | Cluster size distribution |
| `clusevol.dat` | Cluster time evolution |
| `radii.dat` | Cluster radii |
| `fshape.dat` | Cluster shape factor |

---

# 3D systems

## 3D/salr

Three-dimensional SALR system showing cluster formation and dynamic heterogeneity.

### Simulation input

| File | Description |
|-----|-------------|
| `in.pysalr` | LAMMPS simulation input |
| `py_pot.py` | SALR potential implementation |
| `system.data` | Initial configuration |

### Structural observables

| File | Description |
|-----|-------------|
| `sq.dat` | Static structure factor |
| `sqmix.dat` | Partial structure factor |
| `sqcl.dat` | Cluster structure factor |

### Dynamic correlation functions

| File | Description |
|-----|-------------|
| `fkt.dat` | Intermediate scattering function |
| `fskt.dat` | Self intermediate scattering function |

---

## 3D/salr/dyn

Three-dimensional SALR system evaluation of dynamic properties

### Simulation input

| File | Description |
|-----|-------------|
| `in.pysalr` | LAMMPS simulation input |
| `py_pot.py` | SALR potential implementation |
| `system.data` | Initial configuration |

### Structural observables

| File | Description |
|-----|-------------|
| `sq.dat` | Static structure factor |
| `sqmix.dat` | Partial structure factor |
| `sqcl.dat` | Cluster structure factor |

### Dynamic correlation functions

| File | Description |
|-----|-------------|
| `fkt.dat` | Intermediate scattering function |
| `fskt.dat` | Self intermediate scattering function |
| `dyn.dat`| Position and velocity self correlation functions |
| `dynw.dat`| Position and velocity self correlation functions (fequency domain) |

## 3D/LJmix

Example of a **binary Lennard-Jones mixture**.

### Simulation input

| File | Description |
|-----|-------------|
| `in.2LJ` | LAMMPS simulation script |
| `forcefield.lj` | Lennard-Jones force field |
| `input.nml` | Analysis configuration |

### Outputs

| File | Description |
|-----|-------------|
| `gmixsim.dat` | Pair correlations between species |
| `sq.dat` | Static structure factor |
| `sqmix.dat` | Partial structure factor |
| `order.dat` | Structural order parameter |

---

## 3D/Conf

Example of a **confined system**.

### Simulation input

| File | Description |
|-----|-------------|
| `MWwall_cV1.lmp` | LAMMPS simulation input |
| `data.atoms` | Atomic configuration |
| `forcefield.electrode` | Interaction parameters |
| `input.nml` | Analysis configuration |

### Output observables

| File | Description |
|-----|-------------|
| `densprof.dat` | Density profile along confinement direction |
| `qdens.dat` | Charge density profile |
| `sq.dat` | Static structure factor |
| `sqmix.dat` | Partial structure factor |

---

# Running the examples

A typical analysis run is performed with:

```bash
trj_analysis input.nml
```

The parameters controlling the analysis are specified in the **`input.nml`** file included in each example directory.

---

# Purpose of the example outputs

The output files included in this directory correspond to **reference results** obtained with the analysis routines. They can be used to:

- validate a new installation
- verify numerical consistency
- illustrate the expected format of the output files.

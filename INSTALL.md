# Installation Guide for Trajectory Analysis for LAMMPS

This document provides complete instructions for installing prerequisites, compiling the **NetCDF-Fortran** library with the NVIDIA HPC SDK (`nvfortran`), and building the `trj_analysis` toolkit for GPU-accelerated molecular dynamics trajectory analysis.

---

## Table of Contents

1. [Hardware & Software Prerequisites](#1-hardware--software-prerequisites)
2. [The NetCDF-Fortran & nvfortran Requirement](#2-the-netcdf-fortran--nvfortran-requirement)
3. [Automated Installation: `install_netcdff_nvfortran.sh`](#3-automated-installation-install_netcdff_nvfortransh)
4. [Manual NetCDF Installation](#4-manual-netcdf-installation)
   - [4.1 Install Base System Dependencies](#41-install-base-system-dependencies)
   - [4.2 Compile NetCDF-C with nvc](#42-compile-netcdf-c-with-nvc)
   - [4.3 Compile NetCDF-Fortran with nvfortran](#43-compile-netcdf-fortran-with-nvfortran)
5. [Configuring and Compiling `trj_analysis`](#5-configuring-and-compiling-trj_analysis)
   - [5.1 Environment Variables](#51-environment-variables)
   - [5.2 Building with Make](#52-building-with-make)
   - [5.3 GPU Target Architectures](#53-gpu-target-architectures)
6. [HPC Cluster Environment (e.g., CSIC Ladon / Slurm)](#6-hpc-cluster-environment-eg-csic-ladon--slurm)
7. [Testing and Verification](#7-testing-and-verification)
8. [Troubleshooting & FAQ](#8-troubleshooting--faq)

---

## 1. Hardware & Software Prerequisites

### 1.1 Hardware
- **NVIDIA GPU** with Compute Capability $\ge 7.5$ (Turing, Ampere, Ada Lovelace, Hopper, Blackwell):
  - Turing: `cc75` (e.g., RTX 20xx, T4)
  - Ampere: `cc80`, `cc86` (e.g., A100, RTX 30xx)
  - Ada Lovelace: `cc89` (e.g., RTX 40xx, L40)
  - Hopper: `cc90` (e.g., H100)
  - Blackwell: `cc120` (e.g., RTX Pro 4500, B200)

### 1.2 Compilers & Toolchain
- **NVIDIA HPC SDK** (Version 22.x through 26.x or newer):
  - `nvfortran` (CUDA Fortran compiler)
  - `nvc` (C compiler)
  - `nvcc` (CUDA C++ compiler with C++17 support)
  *Download from [NVIDIA Developer](https://developer.nvidia.com/hpc-sdk) or use cluster modules.*
- **GNU Make** ($\ge 4.0$)

### 1.3 External Libraries
- **NetCDF-C** ($\ge 4.9.0$) and **NetCDF-Fortran** ($\ge 4.6.0$) compiled with `nvfortran` (see [Section 2](#2-the-netcdf-fortran--nvfortran-requirement))
- **HDF5** ($\ge 1.10.x$, required by NetCDF-4)
- **FFTW3** (`libfftw3.so` and `fftw3.h`, single- and double-precision)
- **BLAS & LAPACK** (provided by system, OpenBLAS, or NVIDIA HPC SDK math libraries)

---

## 2. The NetCDF-Fortran & nvfortran Requirement

> [!IMPORTANT]
> **Why distribution packages (apt/dnf) cannot be used directly:**
> Fortran module interface files (`.mod`) are **compiler-specific and ABI-incompatible**. Pre-built NetCDF-Fortran packages from Ubuntu (`libnetcdff-dev`) or RHEL/Fedora (`netcdf-fortran-devel`) were compiled with `gfortran`. When `nvfortran` attempts to compile `use netcdf`, it will fail with:
> ```text
> PGF90-F-0004-Corrupt or unrecognized .mod file .../netcdf.mod
> ```
> Therefore, **NetCDF-Fortran must be compiled with `nvfortran`** so that a compatible `netcdf.mod` and shared library (`libnetcdff.so`) are produced.

---

## 3. Automated Installation: `install_netcdff_nvfortran.sh`

The repository provides a script, [`install_netcdff_nvfortran.sh`](file:///Users/elomba/trj_analysis/install_netcdff_nvfortran.sh), that automates downloading, configuring, compiling, and installing NetCDF-C and NetCDF-Fortran with `nvfortran` and `nvc`.

### Features:
- **No root / sudo required**: Installs into user space (`$HOME/software/netcdf-nvfortran` by default) or any custom directory via `--prefix`.
- **Auto-detection**: Detects `nvfortran`, `nvc`, existing NetCDF-C, and HDF5 installations.
- **Official Unidata Releases**: Downloads clean source tarballs with pre-built configure scripts.

### Quick Start:

```bash
# 1. Ensure your compiler is available (or load NVHPC module on a cluster)
nvfortran --version

# 2. Run the installer (installs to $HOME/software/netcdf-nvfortran)
./install_netcdff_nvfortran.sh

# Or specify a custom prefix and parallel jobs:
./install_netcdff_nvfortran.sh --prefix $HOME/local/netcdf -j 8
```

### If NetCDF-C is already installed on your system:
If your system or cluster already has NetCDF-C (built with `-fPIC`), you can compile only NetCDF-Fortran against it:

```bash
./install_netcdff_nvfortran.sh \
  --prefix $HOME/local/netcdf-nvfortran \
  --netcdf-c-dir /path/to/existing/netcdf-c \
  -j 8
```

At the end of the build, the script outputs the exact environment variables needed to build `trj_analysis`.

---

## 4. Manual NetCDF Installation

If you prefer to compile NetCDF manually or need fine-grained control over flags, follow these steps.

### 4.1 Install Base System Dependencies

#### On Ubuntu / Debian:
```bash
sudo apt update
sudo apt install -y build-essential m4 wget curl libhdf5-dev zlib1g-dev libfftw3-dev liblapack-dev libblas-dev
```

#### On Rocky Linux / AlmaLinux / RHEL / Fedora:
```bash
sudo dnf install -y gcc gcc-c++ make m4 wget curl hdf5-devel zlib-devel fftw-devel lapack-devel blas-devel
```

---

### 4.2 Compile NetCDF-C with nvc

Set an installation directory (e.g., `$HOME/software/netcdf-nvfortran`):

```bash
export PREFIX="$HOME/software/netcdf-nvfortran"
mkdir -p /tmp/build_netcdf && cd /tmp/build_netcdf

# Download and extract NetCDF-C
wget https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
tar -xzf netcdf-c-4.9.2.tar.gz
cd netcdf-c-4.9.2

# Configure and install
CC=nvc \
CFLAGS="-O3 -fPIC" \
./configure \
  --prefix="${PREFIX}" \
  --enable-netcdf-4 \
  --enable-shared \
  --disable-dap \
  --disable-doxygen

make -j$(nproc)
make install
```

---

### 4.3 Compile NetCDF-Fortran with nvfortran

Now compile the Fortran library pointing to the NetCDF-C installation:

```bash
cd /tmp/build_netcdf

# Download and extract NetCDF-Fortran
wget https://downloads.unidata.ucar.edu/netcdf-fortran/4.6.1/netcdf-fortran-4.6.1.tar.gz
tar -xzf netcdf-fortran-4.6.1.tar.gz
cd netcdf-fortran-4.6.1

# Point compiler to NetCDF-C
export CC=nvc
export FC=nvfortran
export F90=nvfortran
export F77=nvfortran
export CFLAGS="-O3 -fPIC"
export FFLAGS="-O3 -fPIC"
export FCFLAGS="-O3 -fPIC"
export CPPFLAGS="-I${PREFIX}/include"
export LDFLAGS="-L${PREFIX}/lib -L${PREFIX}/lib64"
export LD_LIBRARY_PATH="${PREFIX}/lib:${PREFIX}/lib64:${LD_LIBRARY_PATH}"

# Configure, compile, and install
./configure \
  --prefix="${PREFIX}" \
  --enable-shared \
  --disable-doxygen

make -j$(nproc)
make install
```

Verify that `netcdf.mod` exists:
```bash
ls -l ${PREFIX}/include/netcdf.mod
```

---

## 5. Configuring and Compiling `trj_analysis`

### 5.1 Environment Variables

The `Makefile` requires several environment variables pointing to your compiler, NetCDF, FFTW, and math libraries:

| Variable | Description | Example Path |
| :--- | :--- | :--- |
| `NVBIN` | Path to NVHPC compiler binaries (`nvfortran`, `nvcc`) | `/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/compilers/bin` |
| `NVINCLUDE` | Path to NVHPC header files and Fortran modules | `/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/compilers/include` |
| `NVLIBS` | Path to NVHPC dynamic runtime libraries | `/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/compilers/lib` |
| `NETCDFINC` | Path to NetCDF include dir (containing `netcdf.mod`) | `$HOME/software/netcdf-nvfortran/include` |
| `NETCDFLIB` | Path to NetCDF libraries (containing `libnetcdff.so`) | `$HOME/software/netcdf-nvfortran/lib` |
| `FFTWINC` | *(Optional)* Path to `fftw3.h` if not in system include path | `/usr/include` or module include path |
| `FFTWLIB` | *(Optional)* Path to `libfftw3.so` if not in system library path | `/usr/lib64` or module lib path |

#### Setting Environment Variables (bash / zsh):

```bash
# Example for a standard custom NetCDF installation:
export NETCDFINC="$HOME/software/netcdf-nvfortran/include"
export NETCDFLIB="$HOME/software/netcdf-nvfortran/lib"

# If NVHPC is in /opt/nvidia/hpc_sdk:
NV_ROOT="/opt/nvidia/hpc_sdk/Linux_x86_64/25.3/compilers"
export NVBIN="$NV_ROOT/bin"
export NVINCLUDE="$NV_ROOT/include"
export NVLIBS="$NV_ROOT/lib"
export PATH="$NVBIN:$PATH"

# Ensure runtime dynamic linker can find NetCDF-Fortran and NVHPC:
export LD_LIBRARY_PATH="$NETCDFLIB:$NVLIBS:$LD_LIBRARY_PATH"
```

---

### 5.2 Building with Make

From the root directory of the repository:

```bash
# Clean previous build artifacts
make clean

# Compile trj_analysis
make
```

The executable will be placed in `bin/trj_analysis`.

You can also override paths directly on the command line:

```bash
make NETCDFINC=/path/to/include NETCDFLIB=/path/to/lib
```

---

### 5.3 GPU Target Architectures

The `Makefile` compiles GPU kernels targeting multiple architectures using the `-gpu` flag:

```makefile
FCOPTS = -O3 -gpu=cc75,cc80,cc86,cc89,cc90,cc120,maxregcount:96 -cudalib=curand
LKOPTS = -cuda -gpu=cc75,cc80,cc86,cc89,cc90,cc120 -c++libs -lnetcdff -lfftw3 -llapack -lblas
```

If your compiler does not support newer architectures (e.g. `cc120` on older NVHPC versions) or you want to optimize strictly for your specific GPU, adjust `FCOPTS`, `LKOPTS`, and the `nvcc` flags in `Makefile`.

---

## 6. HPC Cluster Environment (e.g., CSIC Ladon / Slurm)

On modern scientific clusters running Environment Modules (Lmod) and Slurm (such as the CSIC Ladon cluster):

### 6.1 Loading Required Modules

Check available modules:
```bash
module avail NVHPC
module avail netCDF
module avail FFTW
```

Load the compatible compiler and libraries:
```bash
# Load NVHPC SDK (providing nvfortran, nvc, nvcc)
module load NVHPC/25.3-CUDA-12.8.0

# Load NetCDF modules compiled with NVHPC
module load netCDF/4.9.2-NVHPC-25.3-CUDA-12.8.0
module load netCDF-Fortran/4.6.1-NVHPC-25.3-CUDA-12.8.0
```

On clusters where these modules are available, the module system automatically sets `NETCDFINC`, `NETCDFLIB`, `NVBIN`, `NVINCLUDE`, and `NVLIBS`.

### 6.2 Running on GPU Compute Nodes via Slurm

To launch an interactive GPU session:
```bash
srun -p gpu_v4_short --gres=gpu:1 -n1 --pty bash
```

Or execute directly with `srun`:
```bash
srun -p gpu_v4_short --gres=gpu:1 -n1 ./bin/trj_analysis input.nml 0
```

---

## 7. Testing and Verification

### 7.1 Binary Check

Run the compiled executable without arguments:
```bash
./bin/trj_analysis
```

Expected output:
```text
FORTRAN STOP
!!! Error: You must specify at least ONE argument: the name of the 
    input file e.g.: trj_analisys.exe input_file.nml
    Second argument (optional) is the GPU device number to use, (default 0)
```

Verify that all dynamic dependencies are resolved:
```bash
ldd bin/trj_analysis
```
Make sure `libnetcdff.so`, `libfftw3.so`, `libcudart.so`, and `libcudafor.so` point to valid library paths and are not listed as `not found`.

### 7.2 Example Verification

Verify the installation using one of the reference examples provided in [`examples/`](file:///Users/elomba/trj_analysis/examples):

```bash
cd examples/3D/RTIL+H20_rigid
../../bin/trj_analysis input.nml 0
```

The run will read the trajectory file (`run.nc`), compute structural/dynamic quantities, and produce output data files (`sq.dat`, `gmixsim*.dat`, `trj_analysis.log`). Compare the results with the reference `.dat` files in the example directory.

---

## 8. Troubleshooting & FAQ

### Q1: `Cannot read module file 'netcdf.mod' because it was created by a different compiler`
**Cause**: `netcdf.mod` was created with `gfortran` or another compiler instead of `nvfortran`.  
**Fix**: Use [`install_netcdff_nvfortran.sh`](file:///Users/elomba/trj_analysis/install_netcdff_nvfortran.sh) to compile NetCDF-Fortran with `nvfortran`, and ensure `NETCDFINC` points to that directory.

### Q2: `error while loading shared libraries: libnetcdff.so.7: cannot open shared object file`
**Cause**: The runtime dynamic linker cannot locate the NetCDF-Fortran shared library.  
**Fix**: Add the directory containing `libnetcdff.so` to `LD_LIBRARY_PATH`:
```bash
export LD_LIBRARY_PATH="/path/to/netcdf-nvfortran/lib:$LD_LIBRARY_PATH"
```

### Q3: `nvfortran-Error-Unknown switch: -gpu=cc120`
**Cause**: Your version of NVIDIA HPC SDK predates support for architecture `cc120` (Blackwell).  
**Fix**: In [`Makefile`](file:///Users/elomba/trj_analysis/Makefile), remove `cc120` from `FCOPTS` and `LKOPTS` (leaving e.g. `cc75,cc80,cc86,cc89,cc90`).

### Q4: `undefined reference to curandCreateGenerator` or `curand...`
**Cause**: Missing CUDA random library link flags.  
**Fix**: Ensure `-cudalib=curand` is present in `FCOPTS` and `-cuda` is in `LKOPTS`.

### Q5: `CUDA driver version is insufficient for CUDA runtime version`
**Cause**: The host system's NVIDIA display driver is older than the CUDA toolkit version used by the HPC SDK.  
**Fix**: Update the host NVIDIA GPU driver, or select an older NVHPC SDK / CUDA version compatible with the installed driver.

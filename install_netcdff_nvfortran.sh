#!/usr/bin/env bash
# ==============================================================================
# Script: install_netcdff_nvfortran.sh
# Purpose: Download, build, and install NetCDF-C and NetCDF-Fortran using the
#          NVIDIA HPC SDK (nvc / nvfortran).
#
# Can be run by regular (non-root) users on workstations or HPC clusters.
# ==============================================================================

set -e

# Default software versions
DEFAULT_NC_C_VERSION="4.9.2"
DEFAULT_NC_F_VERSION="4.6.1"

# Default paths
if [ "$EUID" -eq 0 ]; then
    DEFAULT_INSTALL_PREFIX="/usr/local/netcdf-nvfortran"
else
    DEFAULT_INSTALL_PREFIX="$HOME/software/netcdf-nvfortran"
fi
DEFAULT_CACHE_DIR="$HOME/.cache/netcdf_nvfortran_build"

# Configuration variables
INSTALL_PREFIX="$DEFAULT_INSTALL_PREFIX"
CACHE_DIR="$DEFAULT_CACHE_DIR"
NC_C_VERSION="$DEFAULT_NC_C_VERSION"
NC_F_VERSION="$DEFAULT_NC_F_VERSION"
BUILD_JOBS=""
CUSTOM_NC_C_DIR=""
CUSTOM_HDF5_DIR=""
CLEAN_BUILD=0

# Usage / Help function
usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Downloads, compiles, and installs NetCDF-C and NetCDF-Fortran using the
NVIDIA HPC SDK compilers (nvc, nvfortran). Does NOT require root/sudo.

Options:
  --prefix <DIR>          Installation prefix directory
                          [Default: $DEFAULT_INSTALL_PREFIX]
  --netcdf-c-dir <DIR>    Use existing NetCDF-C installation instead of building it
  --hdf5-dir <DIR>        Path to existing HDF5 installation (for NetCDF-C build)
  --cache-dir <DIR>       Directory to store downloaded archives and build trees
                          [Default: $DEFAULT_CACHE_DIR]
  -j, --jobs <N>          Number of parallel make jobs [Default: auto-detect]
  --clean                 Remove build cache directory before building
  -h, --help              Show this help message and exit

Examples:
  # Standard user installation (installs to ~/software/netcdf-nvfortran):
  $(basename "$0")

  # Custom prefix with 8 parallel jobs:
  $(basename "$0") --prefix /opt/netcdf-nvfortran -j 8

  # Build NetCDF-Fortran only, using existing NetCDF-C:
  $(basename "$0") --prefix \$HOME/local --netcdf-c-dir /usr/local

EOF
}

# Parse command-line arguments
while [ $# -gt 0 ]; do
    case "$1" in
        --prefix)
            INSTALL_PREFIX="$2"
            shift 2
            ;;
        --netcdf-c-dir)
            CUSTOM_NC_C_DIR="$2"
            shift 2
            ;;
        --hdf5-dir)
            CUSTOM_HDF5_DIR="$2"
            shift 2
            ;;
        --cache-dir)
            CACHE_DIR="$2"
            shift 2
            ;;
        -j|--jobs)
            BUILD_JOBS="$2"
            shift 2
            ;;
        --clean)
            CLEAN_BUILD=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: Unknown argument '$1'" >&2
            echo "Run '$(basename "$0") --help' for usage." >&2
            exit 1
            ;;
    esac
done

# Determine parallel jobs
if [ -z "$BUILD_JOBS" ]; then
    if command -v nproc >/dev/null 2>&1; then
        BUILD_JOBS=$(nproc)
    elif command -v sysctl >/dev/null 2>&1; then
        BUILD_JOBS=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
    else
        BUILD_JOBS=4
    fi
fi

echo "======================================================================"
echo " NetCDF (C & Fortran) Installer with NVIDIA HPC SDK (nvc / nvfortran)"
echo "======================================================================"
echo " Installation Prefix : $INSTALL_PREFIX"
echo " Download/Build Cache: $CACHE_DIR"
echo " Parallel Jobs       : $BUILD_JOBS"
echo " NetCDF-C Version    : $NC_C_VERSION"
echo " NetCDF-Fortran Vers.: $NC_F_VERSION"
echo "======================================================================"

# Check compiler availability
if ! command -v nvfortran >/dev/null 2>&1; then
    # Look for common NVHPC install locations
    POSSIBLE_NVHPC=$(ls -d /opt/nvidia/hpc_sdk/Linux_x86_64/*/compilers/bin 2>/dev/null | sort -V | tail -n 1 || true)
    if [ -n "$POSSIBLE_NVHPC" ] && [ -x "$POSSIBLE_NVHPC/nvfortran" ]; then
        echo "--> Detected nvfortran at: $POSSIBLE_NVHPC"
        export PATH="$POSSIBLE_NVHPC:$PATH"
        NV_LIB_DIR="$(dirname "$POSSIBLE_NVHPC")/lib"
        [ -d "$NV_LIB_DIR" ] && export LD_LIBRARY_PATH="$NV_LIB_DIR:$LD_LIBRARY_PATH"
    else
        echo "Error: 'nvfortran' compiler not found in PATH." >&2
        echo "Please load your NVHPC environment module (e.g. 'module load NVHPC')" >&2
        echo "or add the compiler binary directory to your PATH." >&2
        exit 1
    fi
fi

if ! command -v nvc >/dev/null 2>&1; then
    echo "Warning: 'nvc' compiler not found in PATH. Checking for gcc fallback for C compilation..."
    if command -v gcc >/dev/null 2>&1; then
        CC_EXEC="gcc"
    else
        echo "Error: Neither 'nvc' nor 'gcc' was found." >&2
        exit 1
    fi
else
    CC_EXEC="nvc"
fi

FC_EXEC="nvfortran"

echo "--> Using C compiler      : $(command -v $CC_EXEC) ($($CC_EXEC --version 2>&1 | head -n 1))"
echo "--> Using Fortran compiler: $(command -v $FC_EXEC) ($($FC_EXEC --version 2>&1 | head -n 1))"

# Check required utilities
for tool in curl wget tar make; do
    if command -v "$tool" >/dev/null 2>&1; then
        FETCH_TOOL="$tool"
        break
    fi
done

download_file() {
    local url="$1"
    local dest="$2"
    if [ -f "$dest" ]; then
        echo "--> Found cached archive: $(basename "$dest")"
        return 0
    fi
    echo "--> Downloading: $url"
    if command -v curl >/dev/null 2>&1; then
        curl -fSL "$url" -o "$dest"
    elif command -v wget >/dev/null 2>&1; then
        wget -q "$url" -O "$dest"
    else
        echo "Error: Neither curl nor wget is available for downloading." >&2
        exit 1
    fi
}

if [ "$CLEAN_BUILD" -eq 1 ] && [ -d "$CACHE_DIR" ]; then
    echo "--> Cleaning build cache directory: $CACHE_DIR"
    rm -rf "$CACHE_DIR"
fi

mkdir -p "$CACHE_DIR"
mkdir -p "$INSTALL_PREFIX"

# Base compiler flags for PIC and Optimization
COMMON_CFLAGS="-O3 -fPIC"
COMMON_FFLAGS="-O3 -fPIC"

# -----------------------------------------------------------------------------
# 1. NetCDF-C Build (if not using existing installation)
# -----------------------------------------------------------------------------
if [ -n "$CUSTOM_NC_C_DIR" ]; then
    echo ""
    echo "======================================================================"
    echo " Step 1: Using Existing NetCDF-C Installation"
    echo "======================================================================"
    NC_C_PREFIX="$CUSTOM_NC_C_DIR"
    echo "--> NetCDF-C prefix: $NC_C_PREFIX"
    if [ ! -f "$NC_C_PREFIX/include/netcdf.h" ]; then
        echo "Error: Could not find $NC_C_PREFIX/include/netcdf.h" >&2
        exit 1
    fi
else
    echo ""
    echo "======================================================================"
    echo " Step 1: Downloading & Compiling NetCDF-C v$NC_C_VERSION"
    echo "======================================================================"

    NC_C_TAR="$CACHE_DIR/netcdf-c-${NC_C_VERSION}.tar.gz"
    NC_C_URL="https://downloads.unidata.ucar.edu/netcdf-c/${NC_C_VERSION}/netcdf-c-${NC_C_VERSION}.tar.gz"
    download_file "$NC_C_URL" "$NC_C_TAR"

    NC_C_SRC="$CACHE_DIR/netcdf-c-${NC_C_VERSION}"
    if [ ! -d "$NC_C_SRC" ]; then
        echo "--> Extracting $(basename "$NC_C_TAR")..."
        tar -xzf "$NC_C_TAR" -C "$CACHE_DIR"
    fi

    # Detect HDF5 if not specified
    HDF5_CONFIG_ARG=""
    NC_C_CPPFLAGS=""
    NC_C_LDFLAGS=""

    if [ -n "$CUSTOM_HDF5_DIR" ]; then
        HDF5_CONFIG_ARG="--with-hdf5=$CUSTOM_HDF5_DIR"
        NC_C_CPPFLAGS="-I${CUSTOM_HDF5_DIR}/include"
        if [ -d "${CUSTOM_HDF5_DIR}/lib64" ]; then
            NC_C_LDFLAGS="-L${CUSTOM_HDF5_DIR}/lib64"
        else
            NC_C_LDFLAGS="-L${CUSTOM_HDF5_DIR}/lib"
        fi
    elif [ -n "$EBROOTHDF5" ]; then
        # EasyBuild / Lmod HDF5 environment variable
        HDF5_CONFIG_ARG="--with-hdf5=$EBROOTHDF5"
        NC_C_CPPFLAGS="-I${EBROOTHDF5}/include"
        NC_C_LDFLAGS="-L${EBROOTHDF5}/lib -L${EBROOTHDF5}/lib64"
    elif [ -f "/usr/include/hdf5.h" ] || [ -f "/usr/include/hdf5/serial/hdf5.h" ]; then
        if [ -d "/usr/include/hdf5/serial" ]; then
            NC_C_CPPFLAGS="-I/usr/include/hdf5/serial"
            NC_C_LDFLAGS="-L/usr/lib/x86_64-linux-gnu/hdf5/serial"
        fi
    fi

    echo "--> Configuring NetCDF-C..."
    cd "$NC_C_SRC"
    
    CC="$CC_EXEC" \
    CFLAGS="$COMMON_CFLAGS" \
    CPPFLAGS="$NC_C_CPPFLAGS $CPPFLAGS" \
    LDFLAGS="$NC_C_LDFLAGS $LDFLAGS" \
    ./configure \
        --prefix="$INSTALL_PREFIX" \
        --disable-dap \
        --enable-netcdf-4 \
        --enable-shared \
        --disable-doxygen \
        $HDF5_CONFIG_ARG

    echo "--> Compiling NetCDF-C (make -j$BUILD_JOBS)..."
    make -j"$BUILD_JOBS"

    echo "--> Installing NetCDF-C to $INSTALL_PREFIX..."
    make install

    NC_C_PREFIX="$INSTALL_PREFIX"
fi

# Determine NetCDF-C library path (lib or lib64)
NC_C_LIBDIR="$NC_C_PREFIX/lib"
if [ -d "$NC_C_PREFIX/lib64" ] && [ ! -d "$NC_C_PREFIX/lib" ]; then
    NC_C_LIBDIR="$NC_C_PREFIX/lib64"
elif [ -d "$NC_C_PREFIX/lib64" ] && [ -f "$NC_C_PREFIX/lib64/libnetcdf.so" ]; then
    NC_C_LIBDIR="$NC_C_PREFIX/lib64"
fi

# -----------------------------------------------------------------------------
# 2. NetCDF-Fortran Build
# -----------------------------------------------------------------------------
echo ""
echo "======================================================================"
echo " Step 2: Downloading & Compiling NetCDF-Fortran v$NC_F_VERSION"
echo "======================================================================"

NC_F_TAR="$CACHE_DIR/netcdf-fortran-${NC_F_VERSION}.tar.gz"
NC_F_URL="https://downloads.unidata.ucar.edu/netcdf-fortran/${NC_F_VERSION}/netcdf-fortran-${NC_F_VERSION}.tar.gz"
download_file "$NC_F_URL" "$NC_F_TAR"

NC_F_SRC="$CACHE_DIR/netcdf-fortran-${NC_F_VERSION}"
if [ ! -d "$NC_F_SRC" ]; then
    echo "--> Extracting $(basename "$NC_F_TAR")..."
    tar -xzf "$NC_F_TAR" -C "$CACHE_DIR"
fi

echo "--> Configuring NetCDF-Fortran with nvfortran..."
cd "$NC_F_SRC"

export CC="$CC_EXEC"
export FC="$FC_EXEC"
export F90="$FC_EXEC"
export F77="$FC_EXEC"
export CFLAGS="$COMMON_CFLAGS"
export FFLAGS="$COMMON_FFLAGS"
export FCFLAGS="$COMMON_FFLAGS"

# Provide paths to NetCDF-C headers and libraries
export CPPFLAGS="-I${NC_C_PREFIX}/include $CPPFLAGS"
export LDFLAGS="-L${NC_C_LIBDIR} $LDFLAGS"
export LD_LIBRARY_PATH="${NC_C_LIBDIR}:${LD_LIBRARY_PATH}"

./configure \
    --prefix="$INSTALL_PREFIX" \
    --enable-shared \
    --disable-doxygen

echo "--> Compiling NetCDF-Fortran (make -j$BUILD_JOBS)..."
make -j"$BUILD_JOBS"

echo "--> Installing NetCDF-Fortran to $INSTALL_PREFIX..."
make install

# Locate Fortran library directory
NC_F_LIBDIR="$INSTALL_PREFIX/lib"
if [ -d "$INSTALL_PREFIX/lib64" ] && [ ! -d "$INSTALL_PREFIX/lib" ]; then
    NC_F_LIBDIR="$INSTALL_PREFIX/lib64"
elif [ -d "$INSTALL_PREFIX/lib64" ] && [ -f "$INSTALL_PREFIX/lib64/libnetcdff.so" ]; then
    NC_F_LIBDIR="$INSTALL_PREFIX/lib64"
fi

# -----------------------------------------------------------------------------
# 3. Verification & Setup Instructions
# -----------------------------------------------------------------------------
echo ""
echo "======================================================================"
echo " Installation Summary & Verification"
echo "======================================================================"

if [ -f "$INSTALL_PREFIX/include/netcdf.mod" ]; then
    echo "[OK] Fortran module found: $INSTALL_PREFIX/include/netcdf.mod"
else
    echo "[WARNING] netcdf.mod not found in $INSTALL_PREFIX/include!" >&2
fi

if [ -f "$NC_F_LIBDIR/libnetcdff.so" ] || [ -f "$NC_F_LIBDIR/libnetcdff.dylib" ] || [ -f "$NC_F_LIBDIR/libnetcdff.a" ]; then
    echo "[OK] Fortran library found: $NC_F_LIBDIR"
else
    echo "[WARNING] libnetcdff was not found in $NC_F_LIBDIR!" >&2
fi

echo ""
echo "======================================================================"
echo " Build Environment Variables for trj_analysis"
echo "======================================================================"
echo ""
echo "To compile trj_analysis with this installation, export the following:"
echo ""
echo "  export NETCDFINC=\"$INSTALL_PREFIX/include\""
echo "  export NETCDFLIB=\"$NC_F_LIBDIR\""
echo "  export PATH=\"$INSTALL_PREFIX/bin:\$PATH\""
echo "  export LD_LIBRARY_PATH=\"$NC_F_LIBDIR:\$LD_LIBRARY_PATH\""
echo ""
echo "You can also add those lines to your ~/.bashrc or job scripts."
echo "======================================================================"

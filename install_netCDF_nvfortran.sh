#!/bin/bash

# ==============================================================================
# INSTALADOR UNIVERSAL: NVHPC 26.1 + NETCDF (C & FORTRAN) - CACHE & PIC FIXED
# ==============================================================================

if [ "$EUID" -ne 0 ]; then
  echo "Error: Este script debe ejecutarse con sudo o como root."
  exit 1
fi

# --- Configuración de Versiones ---
NC_C_VERSION="4.9.2"
NC_F_VERSION="4.6.1"
NV_VER_LONG="2026_261"
NV_VER_SHORT="26.1"
DEFAULT_INSTALL_DIR="/usr/local/netcdf-nvfortran"
NV_BASE_PATH="/opt/nvidia/hpc_sdk/Linux_x86_64/$NV_VER_SHORT/compilers"

# Carpeta de descargas persistente para evitar re-descargas
DL_DIR="$HOME/netcdf_downloads"
mkdir -p "$DL_DIR"

echo "=========================================================="
echo " CONFIGURACIÓN DE INSTALACIÓN"
echo "=========================================================="
read -p "Ruta de instalación para NetCDF [$DEFAULT_INSTALL_DIR]: " USER_PATH
INSTALL_DIR=${USER_PATH:-$DEFAULT_INSTALL_DIR}

# --- 1. Gestión de Dependencias ---
if command -v apt-get &> /dev/null; then
    apt-get update && apt-get install -y wget make m4 g++ libhdf5-dev curl libcurl4-openssl-dev zlib1g-dev libxml2-dev
    HDF5_DIR="/usr"
elif command -v dnf &> /dev/null; then
    dnf install -y wget make m4 gcc-c++ hdf5-devel curl libcurl-devel zlib-devel libxml2-devel
    HDF5_DIR="/usr"
fi

# --- 2. Verificación de NVHPC ---
if [ -f "$NV_BASE_PATH/bin/nvfortran" ]; then
    echo "-> NVHPC detectado en $NV_BASE_PATH. Saltando instalación."
    export PATH="$NV_BASE_PATH/bin:$PATH"
    export LD_LIBRARY_PATH="$NV_BASE_PATH/lib:$LD_LIBRARY_PATH"
else
    echo "--- NVHPC no detectado. Procesando instalación ---"
    NV_TAR="$DL_DIR/nvhpc_${NV_VER_LONG}.tar.gz"
    if [ ! -f "$NV_TAR" ]; then
        echo "Descargando NVHPC (~15GB)..."
        wget "https://developer.download.nvidia.com/hpc-sdk/26.1/nvhpc_${NV_VER_LONG}_Linux_x86_64_cuda_multi.tar.gz" -O "$NV_TAR"
    else
        echo "-> Usando instalador NVHPC ya descargado: $NV_TAR"
    fi

    NV_TMP="/tmp/nvhpc_unpack"
    mkdir -p "$NV_TMP" && tar -xpzf "$NV_TAR" -C "$NV_TMP"
    export NVHPC_INSTALL_DIR=/opt/nvidia/hpc_sdk
    export NVHPC_SILENT=true
    "$NV_TMP/nvhpc_${NV_VER_LONG}_Linux_x86_64_cuda_multi/install"

    export PATH="$NV_BASE_PATH/bin:$PATH"
    export LD_LIBRARY_PATH="$NV_BASE_PATH/lib:$LD_LIBRARY_PATH"
    rm -rf "$NV_TMP"
fi

# --- 3. Preparación del Entorno (Fuerza Bruta para PIC) ---
export CC="nvc -fPIC"
export FC="nvfortran -fPIC"
export F90="nvfortran -fPIC"
export F77="nvfortran -fPIC"
export CPP="nvc -E"
export CFLAGS="-O3"
export FFLAGS="-O3"
export CPPFLAGS="-I/usr/include/hdf5/serial"
export LDFLAGS="-L/usr/lib/x86_64-linux-gnu/hdf5/serial"

# --- 4. Compilación de NetCDF-C ---
echo "--- Paso 1/2: NetCDF-C v$NC_C_VERSION ---"
NC_C_TAR="$DL_DIR/netcdf-c-${NC_C_VERSION}.tar.gz"
if [ ! -f "$NC_C_TAR" ]; then
    wget "https://github.com/Unidata/netcdf-c/archive/v${NC_C_VERSION}.tar.gz" -O "$NC_C_TAR"
else
    echo "-> Usando NetCDF-C ya descargado."
fi

WORK_DIR="/tmp/netcdf_compile_$(date +%s)"
mkdir -p "$WORK_DIR" && cd "$WORK_DIR"
tar -xf "$NC_C_TAR" && cd "netcdf-c-${NC_C_VERSION}"

./configure --prefix="${INSTALL_DIR}" --disable-dap --enable-netcdf-4 --with-hdf5="${HDF5_DIR}"
make -j$(nproc)
make install

# --- 5. Compilación de NetCDF-Fortran ---
echo "--- Paso 2/2: NetCDF-Fortran v$NC_F_VERSION ---"
NC_F_TAR="$DL_DIR/netcdf-fortran-${NC_F_VERSION}.tar.gz"
if [ ! -f "$NC_F_TAR" ]; then
    wget "https://github.com/Unidata/netcdf-fortran/archive/v${NC_F_VERSION}.tar.gz" -O "$NC_F_TAR"
else
    echo "-> Usando NetCDF-Fortran ya descargado."
fi

cd "$WORK_DIR"
tar -xf "$NC_F_TAR" && cd "netcdf-fortran-${NC_F_VERSION}"

# Inyectar rutas del NetCDF-C recién instalado
export LDFLAGS="-L${INSTALL_DIR}/lib $LDFLAGS"
export CPPFLAGS="-I${INSTALL_DIR}/include $CPPFLAGS"
export LD_LIBRARY_PATH="${INSTALL_DIR}/lib:$LD_LIBRARY_PATH"

./configure --prefix="${INSTALL_DIR}" --enable-shared
make -j$(nproc)
make install

# --- 6. Limpieza ---
rm -rf "$WORK_DIR"

echo "=========================================================="
echo " ¡INSTALACIÓN COMPLETADA!"
echo "=========================================================="
echo "Descargas guardadas en: $DL_DIR"
echo "NetCDF instalado en:    $INSTALL_DIR"
echo "----------------------------------------------------------"
echo "Añade esto a tu ~/.bashrc:"
echo "export PATH=$NV_BASE_PATH/bin:$INSTALL_DIR/bin:\$PATH"
echo "export LD_LIBRARY_PATH=$NV_BASE_PATH/lib:$INSTALL_DIR/lib:\$LD_LIBRARY_PATH"
echo "=========================================================="

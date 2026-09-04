# Makefile to build trj_analysis.
#
# The following environment variables must be defined :
#
# NVBIN path of nvfortran binaries
# NVINCLUDE path to nvfortran header files and modules
# NVLIBS path to nvfortran dynamic libraries
# NETCDFINC path to NETCDF and NETCDFF includes
# NETCDFLIB path to NETCDF and NETCDFF libraries
# FFTWINC path to FFTW3 include files (might be in system's path)
# FFTWLIB path to FFTW2 dynamic libs (might be in system's path)
#
#
PATH := $(NVBIN):$(PATH)

# ==========================================
# Directories & Structure
# ==========================================
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# VPATH tells Make to look inside these subdirectories for source files automatically
VPATH = $(SRC_DIR)/core:$(SRC_DIR)/io:$(SRC_DIR)/modules:$(SRC_DIR)/main

# ==========================================
# Compilers and Flags
# ==========================================
FC = nvfortran
# NOTE: Added -module $(OBJ_DIR) to redirect .mod generation and -I$(OBJ_DIR) to read them
FCOPTS = -O3 -gpu=cc75,cc80,cc86,cc89,cc90,cc120,maxregcount:96 -cudalib=curand
FCINC = -I$(NETCDFINC) -I$(NVINCLUDE) $(if $(FFTWINC),-I$(FFTWINC)) -I$(OBJ_DIR) -module $(OBJ_DIR)

F90 = $(FC)
F90OPTS = -O3  
F90INC = $(FCINC)
F90LIBS = $(FCLIBS)

LKOPTS = -cuda -gpu=cc75,cc80,cc86,cc89,cc90,cc120 -c++libs -lnetcdff -lfftw3 -llapack -lblas 
LKLIBS = -L$(NETCDFLIB) -L$(NVLIBS) $(if $(FFTWLIB),-L$(FFTWLIB))

CC = nvcc

# ==========================================
# Targets and Objects
# ==========================================
TARGET = $(BIN_DIR)/trj_analysis

# Base flat list of objects
OBJ_FILES = precision.o thrust.o common.o sorts.o input.o moltools.o netcdf.o cells.o \
	sq.o rdf.o densprof.o thermo.o clusters.o order.o log.o fftwlib.o dynamics.o \
	util.o trj_analysis.o ex-scan.o

# Translate raw object files to be placed in the obj/ directory
OBJS = $(patsubst %.o, $(OBJ_DIR)/%.o, $(OBJ_FILES))

# ==========================================
# Build Rules
# ==========================================
.PHONY: all directories clean

all: directories $(TARGET)

# Create target directories
directories:
	@mkdir -p $(OBJ_DIR) $(BIN_DIR)

# Link everything to create the executable
$(TARGET): $(OBJS)
	$(FC) $(LKOPTS) $(LKINC) -o $@ $(OBJS) $(LKLIBS)

# Compile standalone C++ / CUDA file
$(OBJ_DIR)/ex-scan.o: ex-scan.cu
	$(CC) -O3 --std c++17 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_89,code=sm_89 -gencode arch=compute_90,code=sm_90 -gencode arch=compute_90,code=compute_90 -c $< -o $@

# ------------------------------------------
# Fortran Module Dependencies (CRITICAL)
# ------------------------------------------
# 1. precision.f90 is the root dependency. It builds first.
$(OBJ_DIR)/precision.o: precision.f90
	$(F90) -c $(F90OPTS) $(F90INC) $(F90LIBS) $< -o $@

# 2. sorts.f90 needs precision.o
$(OBJ_DIR)/sorts.o: sorts.f90 $(OBJ_DIR)/precision.o
	$(F90) -c $(F90OPTS) $(F90INC) $(F90LIBS) $< -o $@

# 3. common.cuf needs precision.o
$(OBJ_DIR)/common.o: common.cuf $(OBJ_DIR)/precision.o
	$(FC) -c $(FCOPTS) $(FCINC) $(FCLIBS) $< -o $@

# 4. input.f90 needs common.o and sorts.o
$(OBJ_DIR)/input.o: input.f90 $(OBJ_DIR)/common.o $(OBJ_DIR)/sorts.o
	$(F90) -c $(F90OPTS) $(F90INC) $(F90LIBS) $< -o $@

# 5. Generic rules for the remaining .cuf files
$(OBJ_DIR)/%.o: %.cuf $(OBJ_DIR)/common.o $(OBJ_DIR)/input.o $(OBJ_DIR)/precision.o
	$(FC) -c $(FCOPTS) $(FCINC) $(FCLIBS) $< -o $@

# 6. Generic rules for the remaining .f90 files
$(OBJ_DIR)/%.o: %.f90 $(OBJ_DIR)/common.o $(OBJ_DIR)/input.o $(OBJ_DIR)/precision.o
	$(F90) -c $(F90OPTS) $(F90INC) $(F90LIBS) $< -o $@

# Clean up build artifacts
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

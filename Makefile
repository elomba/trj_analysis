# Makefile to build trj_analysis.
# A word of warning: when making changes in commo
# Specify main local directory
#
PATH := /usr/local/nvidia/hpc_sdk/Linux_x86_64/2024/compilers/bin:$(PATH)
PROG_DIR = /usr/local
# Specify directory to Netcdf libs. The must be compiled with NVIDIA's nvfortran !!
NETCDF = nv_netcdf
NETCDF_INC = $(PROG_DIR)/$(NETCDF)/include
NETCDF_LIB = $(PROG_DIR)/$(NETCDF)/lib64

FC = nvfortran
FCOPTS = -O2 -mavx2 -mno-avx512f -gpu=cc60,cc70,cc80 -cudalib=curand
FCINC = -I$(NETCDF_INC)
FCLIBS = -L$(NETCDF_LIB)

F90 = $(FC)
F90OPTS = -O2 -C -mavx2 -mno-avx512f
F90INC = $(FCINC)
F90LIBS = $(FCLIBS)

LKOPTS =  -cuda -c++libs -gpu=cc60,cc70,cc80 -lnetcdff -lfftw3 -llapack -lblas 
LKINC = -I$(NETCDF_INC)
LKLIBS = -L$(NETCDF_LIB)

CC = nvcc

EXE = -o trj_analysis_n.exe

OBJ = precision.o thrust.o common.o input.o netcdf.o cells.o \
	sq.o rdf.o densprof.o thermo.o clusters.o log.o fftwlib.o dynamics.o\
	util.o trj_analysis.o ex-scan.o

%.o : %.mod
.SUFFIXES : .cuf .f90

%.o: %.cuf common.cuf precision.f90 
	$(FC) -c $(FCOPTS) $(FCINC) $(FCLIBS) $< -o $@

%.o: %.f90 common.cuf precision.f90 
	$(F90) -c $(F90OPTS) $(F90INC) -c $(F90LIBS) $< -o $@

all: $(OBJ) 
	$(FC) $(LKOPTS) $(LKINC) $(EXE) $(OBJ) $(LKLIBS)

ex-scan.o:
	$(CC) -O4 --std c++14 -c ex-scan.cu
common.o: precision.o common.cuf 
	$(FC) -c $(FCOPTS) $(FCINC) $(FCLIBS) common.cuf
precision.o: precision.f90
	$(F90) -c $(F90OPTS) $(F90INC) -c $(F90LIBS) precision.f90

clean:
	rm -f *.o *.mod 

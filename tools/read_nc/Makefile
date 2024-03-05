FC = nvfortran
F90 = nvfortran
OBJ_CORE = general.o netcdflib_v0.o read_nc.o 

FCOPTS = $(INC) -O3 
F90OPTS =  $(INC) -O3
LKOPTS = -L/usr/local/nv_netcdf/lib64 
LKLIBS = -lnetcdff
INC = -I /usr/local/nv_netcdf/include


%.o : %.mod

.SUFFIXES : .f90

	$(FC)  -c $(FCOPTS) $< -o $@
%.o: %.f90
	$(F90)  $(INC) -c $(F9OPTS) $< -o $@

all: $(OBJ_CORE) 
	$(FC) $(LKOPTS) -o read_nc.exe $(OBJ_CORE)  $(LKLIBS)

clean:
	rm -f *.o *.mod *.exe


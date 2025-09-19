TEST = centers
MODS = common.f90 distances.f90 io.f90 trajectories.f90
FC = gfortran
EXE = exe
FCFLAGS = -Ofast

all: build 

build: $(TEST).f90
	$(FC) $(FCFLAGS) $(MODS) -o $(TEST).$(EXE) $<


clean:
	@echo 'Cleaning up...'
	@rm -rf *.$(EXE) *.dwf *.pdb *.mod *.o prof

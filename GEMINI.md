# Gemini Configuration for Trajectory Analysis for LAMMPS

## Project Context
- **Language Stack:** CUDA Fortran (`.cuf`), Fortran 90 (`.f90`), Python (tools), Bash (scripts).
- **Core Domain:** Molecular dynamics trajectory analysis, structural/thermodynamic/dynamic calculations on LAMMPS NetCDF files.
- **Hardware Target:** NVIDIA GPUs (CUDA).

## AI Assistant Directives

### 1. Code Style & Conventions
- **Fortran:** Adhere to modern Fortran standards (F90 and later). Use `implicit none` in all program units, modules, and subroutines. Ensure strict type safety.
- **CUDA Fortran:** Maintain efficient host-device memory transfers. Emphasize performance and parallelism when modifying or creating `.cuf` kernels. Keep array indexing and memory allocations consistent with existing patterns (e.g., using `thrust` or direct CUDA bindings if applicable).
- **Python:** For scripts in `tools/`, use PEP 8 style conventions. Rely on standard data science and MD analysis libraries (like `numpy`, `netCDF4`, `MDAnalysis`, `freud`, etc.) as established by the existing scripts.

### 2. File Organization
- `.cuf` files are for CUDA Fortran modules and kernels.
- `.f90` files are for standard CPU-bound Fortran routines, utilities, and I/O.
- Python tools should remain in the `tools/` directory.

### 3. Build & Dependencies
- Respect the existing `Makefile`. When adding new files, always update the `Makefile` accordingly.
- The project relies on NetCDF and FFTW libraries. Do not change the build configuration without verifying paths to these libraries.

### 4. General Best Practices
- **Performance:** Always consider the memory layout and coalesced memory access when modifying CUDA kernels.
- **Verification:** Before making significant code changes, check the `examples/` directory to understand how the toolkit is currently utilized.
- **Documentation:** Briefly document the intent of any new subroutine or module. Maintain existing comments and CITATION information.

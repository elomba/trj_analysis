# Species order regression

`generate.py` copies the first 20 frames of `reruns/fine/d01_bp20/run.nc` into two
AMBER/LAMMPS NetCDF files. Case `a` retains the original atom order (first type 2);
case `b` puts all 2041 type-1 atoms before the 875 type-2 atoms. Each atom's `id`,
`type`, coordinates, velocities, energy, and stress move together. The script copies
dimensions, global and variable attributes, including the time `scale_factor`.

Generate the inputs on the node with access to the source trajectory and Python
packages `netCDF4` and `numpy`:

```sh
cd /sif/users/e.lomba/Trabajo/IPHS-Mie/trj_analysis_species_order
python tools/test_species_order/generate.py
```

Run each binary on a GPU node, from separate output directories so files from one
run cannot mask another. The original binary is reached through `~/bin/trj_analysis`;
the fixed binary is `bin/trj_analysis_spfix` in this copy:

```sh
cd /sif/users/e.lomba/Trabajo/IPHS-Mie/trj_analysis_species_order
for case in a b; do
  (cd tools/test_species_order/$case && ~/bin/trj_analysis input.nml > old.log 2>&1 &&
   mkdir -p ../old/$case && cp gmixsim.dat sq.dat ../old/$case/)
done
python tools/test_species_order/compare.py tools/test_species_order/old --expect-swapped

for case in a b; do
  (cd tools/test_species_order/$case && ../../../bin/trj_analysis_spfix input.nml > new.log 2>&1)
done
python tools/test_species_order/compare.py
```

Expected: the old binary's `g_11`/`g_22` and `S_11`/`S_22` swap between `a` and
`b`; the fixed binary's partials agree within the comparison tolerance. For the
fixed binary, `g_11` has the larger contact peak because species 1 is the majority
trajectory type 1. Both fixed logs should include `species 1 -> trajectory type 1
(2041 atoms)` and `species 2 -> trajectory type 2 (875 atoms)`.

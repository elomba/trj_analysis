#!/usr/bin/env python3
"""Make two 20-frame LAMMPS/AMBER NetCDF inputs differing only in atom order."""
from pathlib import Path
import argparse
import re

import netCDF4
import numpy as np


HERE = Path(__file__).resolve().parent
SOURCE = Path('/sif/users/e.lomba/Trabajo/IPHS-Mie/reruns/fine/d01_bp20/run.nc')
INPUT = Path('/sif/users/e.lomba/Trabajo/IPHS-Mie/reruns/fine/d01_bp20/input.nml')
ATOM_VARIABLES = {'id', 'type', 'coordinates', 'velocities', 'c_ener', 'c_stress'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, default=SOURCE)
    parser.add_argument('--input', type=Path, default=INPUT)
    parser.add_argument('--output', type=Path, default=HERE)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with netCDF4.Dataset(args.source) as src:
        assert len(src.dimensions['frame']) >= 20
        types = np.asarray(src.variables['type'][0, :])
        assert types[0] == 2, 'case a must begin with trajectory type 2'
        assert np.count_nonzero(types == 1) == 2041
        assert np.count_nonzero(types == 2) == 875
        permutations = {'a': np.arange(len(types)), 'b': np.argsort(types, kind='stable')}
        for case, order in permutations.items():
            folder = args.output / case
            folder.mkdir(exist_ok=True)
            with netCDF4.Dataset(folder / 'run.nc', 'w', format=src.file_format) as dst:
                dst.setncatts({name: src.getncattr(name) for name in src.ncattrs()})
                for name, dim in src.dimensions.items():
                    dst.createDimension(name, None if dim.isunlimited() else len(dim))
                for name, var in src.variables.items():
                    out = dst.createVariable(name, var.datatype, var.dimensions)
                    out.setncatts({attr: var.getncattr(attr) for attr in var.ncattrs() if attr != '_FillValue'})
                    if 'frame' in var.dimensions:
                        data = var[:20, ...]
                        if name in ATOM_VARIABLES:
                            data = np.take(data, order, axis=var.dimensions.index('atom'))
                        out[:20, ...] = data
                    else:
                        out[...] = var[...]

            input_text = args.input.read_text()
            input_text = re.sub(r'(?im)^\s*ncfs_from_to\s*=.*$', '  ncfs_from_to=19 2 20', input_text)
            input_text = re.sub(r'(?im)^\s*nbuffer\s*=.*$', 'nbuffer=4', input_text)
            input_text = re.sub(r'(?im)^\s*jump\s*=.*$', 'jump=1', input_text)
            input_text = re.sub(r'(?im)^\s*tmax\s*=.*$', 'tmax=0.005', input_text)
            input_text = re.sub(r'(?im)^\s*tlimit\s*=.*$', 'tlimit=0.005', input_text)
            input_text = re.sub(r'(?im)^\s*rdf_sq_cl_dyn_sqw_conf_ord\s*=.*$',
                                '  rdf_sq_cl_dyn_sqw_conf_ord=T T F F F F F', input_text)
            (folder / 'input.nml').write_text(input_text)

    with netCDF4.Dataset(args.output / 'a/run.nc') as a, netCDF4.Dataset(args.output / 'b/run.nc') as b:
        assert a.variables['type'][0, 0] == 2 and b.variables['type'][0, 0] == 1
        for name in ATOM_VARIABLES:
            assert np.array_equal(a.variables[name][:20, order], b.variables[name][:20, :])
        assert np.array_equal(a.variables['time'][:], b.variables['time'][:])
    print(f'Wrote and verified {args.output / "a/run.nc"} and {args.output / "b/run.nc"}')


if __name__ == '__main__':
    main()

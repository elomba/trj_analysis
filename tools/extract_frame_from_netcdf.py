import netCDF4 as nc
import numpy as np
import os


def write_lammps_slice(input_file: str, output_file: str, frame_indices, compression: bool = False):
    """
    Core engine: extract specific frames from a LAMMPS/AMBER NetCDF trajectory,
    preserving all dimensions, variables, and attributes exactly.
    """

    if os.path.exists(output_file):
        os.remove(output_file)

    with nc.Dataset(input_file) as src, nc.Dataset(output_file, "w", format="NETCDF4") as dst:

        # ---------------------------------------------------------------
        # 1️⃣  Copy all dimensions
        # ---------------------------------------------------------------
        for name, dim in src.dimensions.items():
            length = len(dim) if not dim.isunlimited() else None
            dst.createDimension(name, length)

        # ---------------------------------------------------------------
        # 2️⃣  Copy all variables (with data + attributes)
        # ---------------------------------------------------------------
        for name, var in src.variables.items():

            # create variable (reuse fill_value if present)
            fill_value = var._FillValue if "_FillValue" in var.ncattrs() else None
            kwargs = {}
            if compression and var.dtype.kind in {"f", "i"}:
                kwargs.update(dict(zlib=True, complevel=4, shuffle=True))
            new = dst.createVariable(
                name, var.dtype, var.dimensions,
                **({} if fill_value is None else {"fill_value": fill_value}),
                **kwargs
            )

            # copy attributes
            for a in var.ncattrs():
                if a != "_FillValue":
                    new.setncattr(a, getattr(var, a))

            # copy data, subsetting 'frame' if present
            if "frame" in var.dimensions:
                frame_axis = var.dimensions.index("frame")
                slicer = [slice(None)] * len(var.dimensions)
                data = []
                for idx in frame_indices:
                    slicer[frame_axis] = idx
                    data.append(var[tuple(slicer)])
                new[:] = np.stack(data, axis=frame_axis)
            else:
                new[:] = var[:]

        # ---------------------------------------------------------------
        # 3️⃣  Copy global attributes (program, version, etc.)
        # ---------------------------------------------------------------
        for a in src.ncattrs():
            dst.setncattr(a, getattr(src, a))

        dst.Conventions = "AMBER"
        dst.ConventionVersion = "1.0"

    print(f"✅  Wrote {len(frame_indices)} frame(s) → {output_file}")


def lammps_subset(input_file: str, output_file: str, every: int = 1, start: int = 0, compression: bool = False):
    """
    Extract every Nth frame from a LAMMPS NetCDF trajectory.
    Keeps exact attributes and structure.

    Parameters
    ----------
    input_file : str
        Original NetCDF trajectory (e.g. 'run.nc')
    output_file : str
        Output subset file (e.g. 'run_subset.nc')
    every : int
        Keep every Nth frame (default 1 keeps all).
    start : int
        Starting frame index.
    compression : bool
        Compress float/int variables using zlib.
    """
    with nc.Dataset(input_file) as f:
        nframes = len(f.dimensions["frame"])
    frame_indices = list(range(start, nframes, every))
    print(f"Extracting {len(frame_indices)} of {nframes} frames "
          f"({frame_indices[:5]}{'...' if len(frame_indices) > 5 else ''})")
    write_lammps_slice(input_file, output_file, frame_indices, compression=compression)
    
write_lammps_slice("run.nc", "run_subset.nc", frame_indices=[0])
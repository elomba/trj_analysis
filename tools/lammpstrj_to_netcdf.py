import numpy as np
from netCDF4 import Dataset
import sys

if len(sys.argv) < 3:
    print("Usage: python lammpstrj_to_lammps_netcdf.py dump.lammpstrj output.nc")
    sys.exit()

infile = sys.argv[1]
outfile = sys.argv[2]

# ----------------------------------------------------
# detect number of atoms
# ----------------------------------------------------

with open(infile) as f:
    for line in f:
        if "ITEM: NUMBER OF ATOMS" in line:
            natoms = int(next(f))
            break

print("Atoms:",natoms)

# ----------------------------------------------------
# create NetCDF file
# ----------------------------------------------------

nc = Dataset(outfile,"w",format="NETCDF4")

# dimensions
nc.createDimension("frame",None)
nc.createDimension("atom",natoms)
nc.createDimension("spatial",3)
nc.createDimension("cell_spatial",3)
nc.createDimension("cell_angular",3)
nc.createDimension("label",10)

# spatial labels
spatial = nc.createVariable("spatial","S1",("spatial",))
spatial[:] = np.array(list("xyz"),dtype="S1")

cell_spatial = nc.createVariable("cell_spatial","S1",("cell_spatial",))
cell_spatial[:] = np.array(list("xyz"),dtype="S1")

cell_ang = nc.createVariable("cell_angular","S1",("cell_angular","label"))
names = ["alpha","beta","gamma"]

for i,n in enumerate(names):
    arr = np.zeros(10,dtype="S1")
    arr[:len(n)] = np.array(list(n),dtype="S1")
    cell_ang[i,:] = arr

# variables
timev = nc.createVariable("time","f4",("frame",))
timev.setncattr("units","femtosecond")
timev.setncattr("scale_factor",np.float32(2.0))

cell_origin = nc.createVariable("cell_origin","f4",("frame","cell_spatial"))
cell_origin.setncattr("units","angstrom")

cell_lengths = nc.createVariable("cell_lengths","f4",("frame","cell_spatial"))
cell_lengths.setncattr("units","angstrom")

cell_angles = nc.createVariable("cell_angles","f4",("frame","cell_angular"))
cell_angles.setncattr("units","degree")

idv = nc.createVariable("id","i4",("frame","atom"))
typev = nc.createVariable("type","i4",("frame","atom"))

coord = nc.createVariable("coordinates","f4",("frame","atom","spatial"))
coord.setncattr("units","angstrom")

velv = nc.createVariable("velocities","f4",("frame","atom","spatial"))
velv.setncattr("units","angstrom/femtosecond")

# global attributes
nc.setncattr("Conventions","AMBER")
nc.setncattr("ConventionVersion","1.0")
nc.setncattr("program","LAMMPS")
nc.setncattr("programVersion","converted_dump")

# ----------------------------------------------------
# streaming conversion
# ----------------------------------------------------

frame = 0

with open(infile) as f:

    while True:

        line = f.readline()
        if not line:
            break

        if "ITEM: TIMESTEP" in line:

            timestep = float(f.readline())

            f.readline()                 # ITEM NUMBER OF ATOMS
            nat = int(f.readline())

            f.readline()                 # BOX BOUNDS

            xlo,xhi = map(float,f.readline().split())
            ylo,yhi = map(float,f.readline().split())
            zlo,zhi = map(float,f.readline().split())

            box = [xhi-xlo,yhi-ylo,zhi-zlo]

            header = f.readline()
            cols = header.split()[2:]

            raw = np.fromfile(f,count=nat*len(cols),sep=' ')
            raw = raw.reshape(nat,len(cols))

            idcol = cols.index("id")
            order = np.argsort(raw[:,idcol])
            raw = raw[order]

            ids = raw[:,cols.index("id")].astype(np.int32)
            types = raw[:,cols.index("type")].astype(np.int32)

            pos = raw[:,[cols.index("x"),cols.index("y"),cols.index("z")]]

            idv[frame,:] = ids
            typev[frame,:] = types
            coord[frame,:,:] = pos

            if "vx" in cols:
                vel = raw[:,[cols.index("vx"),cols.index("vy"),cols.index("vz")]]
                velv[frame,:,:] = vel

            timev[frame] = timestep / 2.0
            cell_lengths[frame,:] = box
            cell_origin[frame,:] = 0.0
            cell_angles[frame,:] = [90.0,90.0,90.0]

            frame += 1

            if frame % 10 == 0:
                print("Processed frame",frame)

nc.close()

print("Finished writing:",outfile)

import numpy as np
import freud
import sys
from netCDF4 import Dataset

trajfile = sys.argv[1]

bins = 200

# open NetCDF trajectory
nc = Dataset(trajfile)

coords = nc.variables["coordinates"]      # (frame, atom, xyz)
cell_lengths = nc.variables["cell_lengths"]

nframes = coords.shape[0]
natoms = coords.shape[1]

rdf = None

for frame in range(nframes):
    print(" Processing frame ", frame)
    positions = coords[frame,:,:].astype(np.float32)

    Lx, Ly, Lz = cell_lengths[frame]

    box = freud.box.Box(Lx=Lx, Ly=Ly, Lz=Lz)

    if rdf is None:
#        rmax = 0.49 * min(Lx, Ly, Lz)
        rmax = 150.0
        rdf = freud.density.RDF(
            bins=bins,
            r_max=rmax
        )

    # fast neighbor search
    aq = freud.locality.AABBQuery(box, positions)

    nlist = aq.query(
        positions,
        dict(mode='ball', r_max=rdf.r_max, exclude_ii=True)
    ).toNeighborList()

    rdf.compute(system=(box, positions),
                neighbors=nlist,
                reset=False)

    if frame % 10 == 0:
        print("processed frame", frame)

r = rdf.bin_centers
g = rdf.rdf

np.savetxt("gr.dat", np.column_stack((r,g)), header="r g(r)")

print("Finished. Frames processed:", nframes)

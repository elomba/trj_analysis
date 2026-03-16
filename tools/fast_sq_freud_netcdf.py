import numpy as np
import freud
import sys
from netCDF4 import Dataset

trajfile = sys.argv[1]

# parameters
bins = 60
kmax = 0.5

# open trajectory
nc = Dataset(trajfile)

coords = nc.variables["coordinates"]      # (frame, atom, xyz)
cell_lengths = nc.variables["cell_lengths"]

nframes = coords.shape[0]
natoms = coords.shape[1]

sq = None

for frame in range(nframes):

    positions = coords[frame,:,:].astype(np.float32)

    Lx, Ly, Lz = cell_lengths[frame]

    box = freud.box.Box(Lx=Lx, Ly=Ly, Lz=Lz)

    kmin = 2*np.pi/max(Lx,Ly,Lz)

    if sq is None:

        sq = freud.diffraction.StaticStructureFactorDirect(
            bins=bins,
            k_min=kmin,
            k_max=kmax
        )

    sq.compute((box, positions), reset=False)

    if frame % 10 == 0:
        print("processed frame", frame)

k = sq.bin_centers
S = sq.S_k

np.savetxt("Sq.dat", np.column_stack((k, S)), header="k S(k)")

print("Finished. Frames processed:", nframes)

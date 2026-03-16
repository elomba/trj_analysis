import MDAnalysis as mda

u = mda.Universe("system.data", format="DATA")
u.atoms.write("topology.pdb")

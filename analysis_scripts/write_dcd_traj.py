import MDAnalysis

import sys
trajectory=sys.argv[1]

u = MDAnalysis.Universe(trajectory)

from MDAnalysis.coordinates.DCD import DCDWriter
d=DCDWriter(filename="WAT-E-pos-1.dcd",n_atoms=2100,convert_units=True,step=1)
for ts in u.trajectory:
    d.write(u)
d.close()

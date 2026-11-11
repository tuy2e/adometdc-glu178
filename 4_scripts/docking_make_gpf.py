import sys, numpy as np
rec=sys.argv[1]; name=sys.argv[2]
types=sorted({l[77:79].strip() for l in open(rec) if l.startswith(('ATOM','HETATM'))})
c=np.load('put_xyz.npy').mean(0)
lig=['C','N','HD']
g=f"""npts 60 60 60
gridfld {name}.maps.fld
spacing 0.375
receptor_types {' '.join(types)}
ligand_types {' '.join(lig)}
receptor {rec}
gridcenter {c[0]:.3f} {c[1]:.3f} {c[2]:.3f}
smooth 0.5
"""
for t in lig: g+=f"map {name}.{t}.map\n"
g+=f"elecmap {name}.e.map\ndsolvmap {name}.d.map\ndielectric -0.1465\n"
open(f'{name}.gpf','w').write(g)

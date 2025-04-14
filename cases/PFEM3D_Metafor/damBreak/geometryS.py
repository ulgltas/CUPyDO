import os, gmsh
from gmsh import model as sh
gmsh.initialize()

# Mesh Parameters

L = 0.146
W = 0.012
H = 0.08

N = 32
M = 5

# Points list

p = list()

p.append(sh.occ.addPoint(2*L, 0, 0))
p.append(sh.occ.addPoint(2*L+W, 0, 0))
p.append(sh.occ.addPoint(2*L+W, H, 0))
p.append(sh.occ.addPoint(2*L, H, 0))

# Lines list

l = list()

l.append(sh.occ.addLine(p[0], p[1]))
l.append(sh.occ.addLine(p[1], p[2]))
l.append(sh.occ.addLine(p[2], p[3]))
l.append(sh.occ.addLine(p[3], p[0]))

# Physical surface

k = sh.occ.addCurveLoop(l)
s = sh.occ.addPlaneSurface([k])
sh.occ.synchronize()

sh.mesh.setTransfiniteCurve(l[0], M)
sh.mesh.setTransfiniteCurve(l[1], N)
sh.mesh.setTransfiniteCurve(l[2], M)
sh.mesh.setTransfiniteCurve(l[3], N)

sh.mesh.setTransfiniteSurface(s)
sh.mesh.setRecombine(2, s)

# Physical boundary

sh.addPhysicalGroup(2, [s], name='Solid')
sh.addPhysicalGroup(1, l[1:], name='FSI')
sh.addPhysicalGroup(1, l[0:1], name='Base')

# Write the mesh file

sh.mesh.generate(2)
gmsh.write(f'{os.path.dirname(__file__)}/geometryS.msh')
gmsh.fltk.run()
gmsh.finalize()
import os, gmsh
from gmsh import model as sh
gmsh.initialize()

# Mesh Parameters

L = 0.146
W = 0.012
H = 0.08
d = 0.003

# Points list

p = list()

p.append(sh.occ.addPoint(0, 0, 0, d))
p.append(sh.occ.addPoint(4*L, 0, 0, d))
p.append(sh.occ.addPoint(4*L, 4*L, 0, d))
p.append(sh.occ.addPoint(0, 3*L, 0, d))
p.append(sh.occ.addPoint(L, 0, 0, d))
p.append(sh.occ.addPoint(L, 2*L, 0, d))
p.append(sh.occ.addPoint(0, 2*L, 0, d))
p.append(sh.occ.addPoint(2*L, 0, 0, d))
p.append(sh.occ.addPoint(2*L, H, 0, d))
p.append(sh.occ.addPoint(2*L+W, H, 0, d))
p.append(sh.occ.addPoint(2*L+W, 0, 0, d))

# Lines list

l = list()
h = list()

l.append(sh.occ.addLine(p[3], p[6]))
l.append(sh.occ.addLine(p[6], p[0]))
l.append(sh.occ.addLine(p[0], p[4]))
l.append(sh.occ.addLine(p[4], p[5]))
l.append(sh.occ.addLine(p[5], p[6]))
l.append(sh.occ.addLine(p[4], p[7]))
l.append(sh.occ.addLine(p[10], p[1]))
l.append(sh.occ.addLine(p[1], p[2]))

h.append(sh.occ.addLine(p[7], p[8]))
h.append(sh.occ.addLine(p[8], p[9]))
h.append(sh.occ.addLine(p[9], p[10]))

# Physical surface

k = sh.occ.addCurveLoop(l[1:5])
s = sh.occ.addPlaneSurface([k])
sh.occ.synchronize()

A = 5
B = 32

sh.mesh.setTransfiniteCurve(h[0], B)
sh.mesh.setTransfiniteCurve(h[1], A)
sh.mesh.setTransfiniteCurve(h[2], B)

# Physical boundary

sh.addPhysicalGroup(2, [s], name='Fluid')
sh.addPhysicalGroup(1, l[0:3]+l[5:], name='Reservoir')
sh.addPhysicalGroup(1, h, name='FSI')

# Write the mesh

sh.mesh.generate(2)
gmsh.write(f'{os.path.dirname(__file__)}/geometryF.msh')
gmsh.fltk.run()
gmsh.finalize()
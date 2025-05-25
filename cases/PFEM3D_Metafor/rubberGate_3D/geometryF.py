import os, gmsh
from gmsh import model as sh
gmsh.initialize()

# Mesh Parameters

H = 0.079
S = 0.005
D = 0.14
W = 0.05
L1 = 0.1
L2 = 0.1

d = 2.5e-3
E = 1e-5

# Points list

p = list()

p.append(sh.occ.addPoint(0, -W, 0, d))
p.append(sh.occ.addPoint(L1, -W, 0, d))
p.append(sh.occ.addPoint(L1+S+L2, -W, 0, d))
p.append(sh.occ.addPoint(0, -W, D, d))
p.append(sh.occ.addPoint(L1, -W, D, d))
p.append(sh.occ.addPoint(L1, -W, H, d))
p.append(sh.occ.addPoint(L1, -W, E, d))
p.append(sh.occ.addPoint(L1+S+L2, -W, H, d))

p.append(sh.occ.addPoint(0, W, 0, d))
p.append(sh.occ.addPoint(L1, W, 0, d))
p.append(sh.occ.addPoint(L1+S+L2, W, 0, d))
p.append(sh.occ.addPoint(0, W, D, d))
p.append(sh.occ.addPoint(L1, W, D, d))
p.append(sh.occ.addPoint(L1, W, H, d))
p.append(sh.occ.addPoint(L1, W, E, d))
p.append(sh.occ.addPoint(L1+S+L2, W, H, d))

# Lines list

l = list()

l.append(sh.occ.addLine(p[6], p[5]))
l.append(sh.occ.addLine(p[5], p[4]))
l.append(sh.occ.addLine(p[4], p[3]))
l.append(sh.occ.addLine(p[3], p[0]))
l.append(sh.occ.addLine(p[0], p[1]))
l.append(sh.occ.addLine(p[1], p[6]))

l.append(sh.occ.addLine(p[14], p[13]))
l.append(sh.occ.addLine(p[13], p[12]))
l.append(sh.occ.addLine(p[12], p[11]))
l.append(sh.occ.addLine(p[11], p[8]))
l.append(sh.occ.addLine(p[8], p[9]))
l.append(sh.occ.addLine(p[9], p[14]))

l.append(sh.occ.addLine(p[0], p[8]))
l.append(sh.occ.addLine(p[1], p[9]))
l.append(sh.occ.addLine(p[2], p[10]))
l.append(sh.occ.addLine(p[3], p[11]))
l.append(sh.occ.addLine(p[4], p[12]))
l.append(sh.occ.addLine(p[5], p[13]))
l.append(sh.occ.addLine(p[6], p[14]))

h = list()

h.append(sh.occ.addLine(p[1], p[2]))
h.append(sh.occ.addLine(p[2], p[7]))
h.append(sh.occ.addLine(p[7], p[5]))

h.append(sh.occ.addLine(p[9], p[10]))
h.append(sh.occ.addLine(p[10], p[15]))
h.append(sh.occ.addLine(p[15], p[13]))

k = list()
q = list()

k.append(sh.occ.addCurveLoop([l[0], l[1], l[2], l[3], l[4], l[5]]))
k.append(sh.occ.addCurveLoop([l[6], l[7], l[8], l[9], l[10], l[11]]))
k.append(sh.occ.addCurveLoop([l[10], l[12], l[4], l[13]]))
k.append(sh.occ.addCurveLoop([l[2], l[16], l[8], l[15]]))
k.append(sh.occ.addCurveLoop([l[3], l[15], l[9], l[12]]))
k.append(sh.occ.addCurveLoop([l[13], l[5], l[18], l[11]]))
k.append(sh.occ.addCurveLoop([l[0], l[18], l[6], l[17]]))
k.append(sh.occ.addCurveLoop([l[1], l[17], l[7], l[16]]))

q.append(sh.occ.addCurveLoop([h[0], l[14], h[3], l[13]]))
q.append(sh.occ.addCurveLoop([h[0], h[1], h[2], l[0], l[5]]))
q.append(sh.occ.addCurveLoop([h[3], h[4], h[5], l[6], l[11]]))

# Physical volume

s = list()
u = list()

for a in k: s.append(sh.occ.addPlaneSurface([a]))
for a in q: u.append(sh.occ.addPlaneSurface([a]))
sh.occ.synchronize()

h = sh.occ.addSurfaceLoop(s)
v = sh.occ.addVolume([h])
sh.occ.synchronize()

# Physical boundary

sh.addPhysicalGroup(3, [v], name='Fluid')
sh.addPhysicalGroup(2, s[6:7], name='FSI')
sh.addPhysicalGroup(2, s[0:2]+s[4:5]+s[7:]+u, name='Reservoir')
sh.addPhysicalGroup(2, s[2:3], name='Bottom')

# Mesh characteristic size

sh.mesh.field.add('Distance', 1)
sh.mesh.field.setNumber(1, 'Sampling', 1e3)
sh.mesh.field.setNumbers(1, 'SurfacesList', s[0:2]+s[3:]+u)

sh.mesh.field.add('MathEval', 2)
sh.mesh.field.setString(2, 'F', f'Min({d}+4*F1^2, 1e-2)')

sh.mesh.field.setAsBackgroundMesh(2)
gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)

# Write the mesh

sh.mesh.generate(3)
gmsh.write(f'{os.path.dirname(__file__)}/geometryF.msh')
gmsh.fltk.run()
gmsh.finalize()
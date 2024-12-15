import os, gmsh
from gmsh import model as sh
gmsh.initialize()

# Mesh Parameters

H = 0.4
B = 0.25
L1 = 0.15
L2 = 0.2

BS = 0.02
HS = 0.1
W = 0.005

d = 2e-3
N = 40
M = 9
P = 4

# Points list

p = list()

p.append(sh.occ.addPoint(0, 0, 0, d))
p.append(sh.occ.addPoint(L1+L2, 0, 0, d))
p.append(sh.occ.addPoint(L1+L2, B, 0, d))
p.append(sh.occ.addPoint(0, B, 0, d))
p.append(sh.occ.addPoint(0, 0, H, d))
p.append(sh.occ.addPoint(L1+L2, 0, H, d))
p.append(sh.occ.addPoint(L1+L2, B, H, d))
p.append(sh.occ.addPoint(0, B, H, d))

p.append(sh.occ.addPoint(L1, (B-BS)/2, 0, d))
p.append(sh.occ.addPoint(L1+W, (B-BS)/2, 0, d))
p.append(sh.occ.addPoint(L1+W, (B+BS)/2, 0, d))
p.append(sh.occ.addPoint(L1, (B+BS)/2, 0, d))
p.append(sh.occ.addPoint(L1, (B-BS)/2, HS, d))
p.append(sh.occ.addPoint(L1+W, (B-BS)/2, HS, d))
p.append(sh.occ.addPoint(L1+W, (B+BS)/2, HS, d))
p.append(sh.occ.addPoint(L1, (B+BS)/2, HS, d))

p.append(sh.occ.addPoint(0, (B-BS)/2, 0, d))
p.append(sh.occ.addPoint(0, (B+BS)/2, 0, d))
p.append(sh.occ.addPoint(L1+L2, (B-BS)/2, 0, d))
p.append(sh.occ.addPoint(L1+L2, (B+BS)/2, 0, d))

# Lines list

l = list()

l.append(sh.occ.addLine(p[0], p[1]))
l.append(sh.occ.addLine(p[2], p[3]))
l.append(sh.occ.addLine(p[4], p[5]))
l.append(sh.occ.addLine(p[5], p[6]))
l.append(sh.occ.addLine(p[6], p[7]))
l.append(sh.occ.addLine(p[7], p[4]))

l.append(sh.occ.addLine(p[8], p[9]))
l.append(sh.occ.addLine(p[9], p[10]))
l.append(sh.occ.addLine(p[10], p[11]))
l.append(sh.occ.addLine(p[8], p[11]))
l.append(sh.occ.addLine(p[12], p[13]))
l.append(sh.occ.addLine(p[13], p[14]))
l.append(sh.occ.addLine(p[14], p[15]))
l.append(sh.occ.addLine(p[15], p[12]))

l.append(sh.occ.addLine(p[0], p[4]))
l.append(sh.occ.addLine(p[1], p[5]))
l.append(sh.occ.addLine(p[2], p[6]))
l.append(sh.occ.addLine(p[3], p[7]))
l.append(sh.occ.addLine(p[8], p[12]))
l.append(sh.occ.addLine(p[9], p[13]))
l.append(sh.occ.addLine(p[10], p[14]))
l.append(sh.occ.addLine(p[11], p[15]))

l.append(sh.occ.addLine(p[16], p[8]))
l.append(sh.occ.addLine(p[17], p[11]))
l.append(sh.occ.addLine(p[18], p[9]))
l.append(sh.occ.addLine(p[19], p[10]))

l.append(sh.occ.addLine(p[16], p[0]))
l.append(sh.occ.addLine(p[16], p[17]))
l.append(sh.occ.addLine(p[17], p[3]))

l.append(sh.occ.addLine(p[1], p[18]))
l.append(sh.occ.addLine(p[18], p[19]))
l.append(sh.occ.addLine(p[19], p[2]))

# Surfaces list

k = list()
s = list()
sh.occ.synchronize()

k.append(sh.occ.addCurveLoop([l[26], l[27], l[28], l[17], l[5], l[14]]))
k.append(sh.occ.addCurveLoop([l[29], l[30], l[31], l[16], l[3], l[15]]))
k.append(sh.occ.addCurveLoop([l[0], l[29], l[24], l[6], l[22], l[26]]))
k.append(sh.occ.addCurveLoop([l[23], l[8], l[25], l[31], l[1], l[28]]))

k.append(sh.occ.addCurveLoop([l[22], l[9], l[23], l[27]]))
k.append(sh.occ.addCurveLoop([l[24], l[30], l[25], l[7]]))
k.append(sh.occ.addCurveLoop([l[9], l[18], l[13], l[21]]))
k.append(sh.occ.addCurveLoop([l[7], l[19], l[11], l[20]]))
k.append(sh.occ.addCurveLoop([l[6], l[19], l[10], l[18]]))
k.append(sh.occ.addCurveLoop([l[8], l[20], l[12], l[21]]))

k.append(sh.occ.addCurveLoop([l[10], l[11], l[12], l[13]]))
k.append(sh.occ.addCurveLoop([l[1], l[16], l[4], l[17]]))
k.append(sh.occ.addCurveLoop([l[0], l[15], l[2], l[14]]))
k.append(sh.occ.addCurveLoop([l[2], l[3], l[4], l[5]]))

for a in k: s.append(sh.occ.addPlaneSurface([a]))
sh.occ.synchronize()

sh.mesh.setTransfiniteCurve(l[18], N)
sh.mesh.setTransfiniteCurve(l[19], N)
sh.mesh.setTransfiniteCurve(l[20], N)
sh.mesh.setTransfiniteCurve(l[21], N)

sh.mesh.setTransfiniteCurve(l[7], M)
sh.mesh.setTransfiniteCurve(l[9], M)
sh.mesh.setTransfiniteCurve(l[11], M)
sh.mesh.setTransfiniteCurve(l[13], M)

sh.mesh.setTransfiniteCurve(l[6], P)
sh.mesh.setTransfiniteCurve(l[8], P)
sh.mesh.setTransfiniteCurve(l[10], P)
sh.mesh.setTransfiniteCurve(l[12], P)

for a in s[6:11]: sh.mesh.setTransfiniteSurface(a)

# Volumes list

h = sh.occ.addSurfaceLoop(s)
v = sh.occ.addVolume([h])
sh.occ.synchronize()

# Physical surface

sh.addPhysicalGroup(3, [v], name='Fluid')
sh.addPhysicalGroup(2, s[6:11], name='FSI')
sh.addPhysicalGroup(2, s[0:1]+s[11:14], name='Inlet')
sh.addPhysicalGroup(2, s[2:6], name='Bottom')
sh.addPhysicalGroup(2, s[1:2], name='Outlet')

# Mesh characteristic size

sh.mesh.field.add('Distance', 1)
sh.mesh.field.setNumber(1, 'Sampling', 1e3)
sh.mesh.field.setNumbers(1, 'SurfacesList', s[6:11])

sh.mesh.field.add('MathEval', 2)
sh.mesh.field.setString(2, 'F', f'Min({d}+0.1*F1, 3e-2)')

sh.mesh.field.setAsBackgroundMesh(2)
gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)

# Write the mesh

sh.mesh.generate(3)
gmsh.write(f'{os.path.dirname(__file__)}/geometryF.msh')
gmsh.fltk.run()
gmsh.finalize()
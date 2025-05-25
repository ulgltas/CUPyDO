import os, gmsh
from gmsh import model as sh
gmsh.initialize()

# Mesh Parameters

H = 0.079
S = 0.005
W = 0.05
L1 = 0.1

E = 1e-5
N = 60
R = 75
M = 5

# Points list

p = list()

p.append(sh.occ.addPoint(L1, -W, H))
p.append(sh.occ.addPoint(L1+S, -W, H))
p.append(sh.occ.addPoint(L1, -W, E))
p.append(sh.occ.addPoint(L1+S, -W, E))

p.append(sh.occ.addPoint(L1, W, H))
p.append(sh.occ.addPoint(L1+S, W, H))
p.append(sh.occ.addPoint(L1, W, E))
p.append(sh.occ.addPoint(L1+S, W, E))

# Lines list

l = list()

l.append(sh.occ.addLine(p[3], p[1]))
l.append(sh.occ.addLine(p[1], p[0]))
l.append(sh.occ.addLine(p[0], p[2]))
l.append(sh.occ.addLine(p[2], p[3]))

l.append(sh.occ.addLine(p[7], p[5]))
l.append(sh.occ.addLine(p[5], p[4]))
l.append(sh.occ.addLine(p[4], p[6]))
l.append(sh.occ.addLine(p[6], p[7]))

l.append(sh.occ.addLine(p[0], p[4]))
l.append(sh.occ.addLine(p[1], p[5]))
l.append(sh.occ.addLine(p[2], p[6]))
l.append(sh.occ.addLine(p[3], p[7]))

k = list()
s = list()

k.append(sh.occ.addCurveLoop([l[0], l[1], l[2], l[3]]))
k.append(sh.occ.addCurveLoop([l[4], l[5], l[6], l[7]]))
k.append(sh.occ.addCurveLoop([l[0], l[9], l[4], l[11]]))
k.append(sh.occ.addCurveLoop([l[2], l[10], l[6], l[8]]))
k.append(sh.occ.addCurveLoop([l[1], l[8], l[5], l[9]]))
k.append(sh.occ.addCurveLoop([l[3], l[10], l[7], l[11]]))

for a in k: s.append(sh.occ.addPlaneSurface([a]))
sh.occ.synchronize()

sh.mesh.setTransfiniteCurve(l[0], N)
sh.mesh.setTransfiniteCurve(l[2], N)
sh.mesh.setTransfiniteCurve(l[4], N)
sh.mesh.setTransfiniteCurve(l[6], N)

sh.mesh.setTransfiniteCurve(l[1], M)
sh.mesh.setTransfiniteCurve(l[3], M)
sh.mesh.setTransfiniteCurve(l[5], M)
sh.mesh.setTransfiniteCurve(l[7], M)

sh.mesh.setTransfiniteCurve(l[8], R)
sh.mesh.setTransfiniteCurve(l[9], R)
sh.mesh.setTransfiniteCurve(l[10], R)
sh.mesh.setTransfiniteCurve(l[11], R)

# Physical surface

for a in s: sh.mesh.setTransfiniteSurface(a)
for a in s: sh.mesh.setRecombine(2, a)

h = sh.occ.addSurfaceLoop(s)
v = sh.occ.addVolume([h])

sh.occ.synchronize()
sh.mesh.setTransfiniteVolume(v)
sh.mesh.setRecombine(3, v)

# Physical boundary

sh.addPhysicalGroup(3, [v], name='Solid')
sh.addPhysicalGroup(2, s[3:4], name='FSI')
sh.addPhysicalGroup(2, s[4:5], name='Base')

# Write the mesh file

sh.mesh.generate(3)
gmsh.write(f'{os.path.dirname(__file__)}/geometryS.msh')
gmsh.fltk.run()
gmsh.finalize()
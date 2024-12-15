import os, gmsh
from gmsh import model as sh
gmsh.initialize()

# Mesh Parameters

H = 0.079
S = 0.005
D = 0.14
L1 = 0.1
L2 = 0.1

d = 1e-3
E = 1e-5
N = 120

# Points list

p = list()

p.append(sh.occ.addPoint(0, 0, 0, d))
p.append(sh.occ.addPoint(L1, 0, 0, d))
p.append(sh.occ.addPoint(L1+S+L2, 0, 0, d))
p.append(sh.occ.addPoint(0, D, 0, d))
p.append(sh.occ.addPoint(L1, D, 0, d))
p.append(sh.occ.addPoint(L1, H, 0, d))
p.append(sh.occ.addPoint(L1, E, 0, d))

# Lines list

l = list()

l.append(sh.occ.addLine(p[0], p[1]))
l.append(sh.occ.addLine(p[1], p[6]))
l.append(sh.occ.addLine(p[5], p[6]))
l.append(sh.occ.addLine(p[4], p[5]))
l.append(sh.occ.addLine(p[4], p[3]))
l.append(sh.occ.addLine(p[3], p[0]))

h = list()
h.append(sh.occ.addLine(p[1], p[2]))

# Physical surface

k = sh.occ.addCurveLoop(l)
s = sh.occ.addPlaneSurface([k])
sh.occ.synchronize()

sh.mesh.setTransfiniteCurve(l[2], N)

# Physical boundary

sh.addPhysicalGroup(2, [s], name='Fluid')
sh.addPhysicalGroup(1, l[2:3], name='FSI')
sh.addPhysicalGroup(1, l[0:1]+l[5:6], name='Reservoir')
sh.addPhysicalGroup(1, h+l[3:4], name='Refine')

# Mesh characteristic size

sh.mesh.field.add('Distance', 1)
sh.mesh.field.setNumber(1, 'Sampling', 1e4)
sh.mesh.field.setNumbers(1, 'CurvesList', l[2:5]+h)

sh.mesh.field.add('MathEval', 2)
sh.mesh.field.setString(2, 'F', f'Min({d}+0.1*F1, 5e-3)')

sh.mesh.field.setAsBackgroundMesh(2)
gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)

# Write the mesh

sh.mesh.generate(2)
gmsh.write(f'{os.path.dirname(__file__)}/geometryF.msh')
gmsh.fltk.run()
gmsh.finalize()
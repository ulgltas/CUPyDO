import os,gmsh
from gmsh import model as sh
gmsh.initialize()

# |---------------------------|
# |   Mesh Size Parameters    |
# |---------------------------|

L = 1
HF = 0.2
d = 0.01
N = 101

# |----------------------------------|
# |   Points and Lines Definition    |
# |----------------------------------|

p = list()

p.append(sh.occ.addPoint(L,0,0,d))
p.append(sh.occ.addPoint(0,0,0,d))
p.append(sh.occ.addPoint(L,HF,0,d))
p.append(sh.occ.addPoint(0,HF,0,d))

# Lines List

l = list()

l.append(sh.occ.addLine(p[0],p[1]))
l.append(sh.occ.addLine(p[0],p[2]))
l.append(sh.occ.addLine(p[2],p[3]))
l.append(sh.occ.addLine(p[3],p[1]))

# |------------------------------------|
# |   Physical Surface and Boundary    |
# |------------------------------------|

k = sh.occ.addCurveLoop(l)
s = sh.occ.addPlaneSurface([k])

sh.occ.synchronize()
sh.mesh.setTransfiniteCurve(l[2],N)

# Boundaries

sh.addPhysicalGroup(2,[s],name='Fluid')
sh.addPhysicalGroup(1,[l[2]],name='FSInterface')
sh.addPhysicalGroup(1,[l[0],l[1],l[3]],name='Wall')

# |--------------------------|
# |   Write the Mesh File    |
# |--------------------------|

sh.mesh.generate(2)
gmsh.model.mesh.reverse()
gmsh.write(os.path.dirname(__file__)+'/geometryF.msh')
gmsh.fltk.run()
gmsh.finalize()
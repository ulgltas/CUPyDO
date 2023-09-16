import os,gmsh
from gmsh import model as sh
gmsh.initialize()

# |---------------------------|
# |   Mesh Size Parameters    |
# |---------------------------|

L = 1
HF = 0.2
HS = 0.1
N = 101
M = 10

# |----------------------------------|
# |   Points and Lines Definition    |
# |----------------------------------|

p = list()

p.append(sh.occ.addPoint(0,HF,0))
p.append(sh.occ.addPoint(L,HF,0))
p.append(sh.occ.addPoint(L,HF+HS,0))
p.append(sh.occ.addPoint(0,HF+HS,0))

# Lines List

l = list()

l.append(sh.occ.addLine(p[0],p[1]))
l.append(sh.occ.addLine(p[1],p[2]))
l.append(sh.occ.addLine(p[2],p[3]))
l.append(sh.occ.addLine(p[3],p[0]))

# |------------------------------------|
# |   Physical Surface and Boundary    |
# |------------------------------------|

k = sh.occ.addCurveLoop(l)
s = sh.occ.addPlaneSurface([k])
sh.occ.synchronize()

sh.mesh.setTransfiniteCurve(l[0],N)
sh.mesh.setTransfiniteCurve(l[1],M)
sh.mesh.setTransfiniteCurve(l[2],N)
sh.mesh.setTransfiniteCurve(l[3],M)

sh.mesh.setTransfiniteSurface(s)
sh.mesh.setRecombine(2,s)

# Physical Boundary

sh.addPhysicalGroup(2,[s],name='Solid')
sh.addPhysicalGroup(1,[l[0]],name='FSInterface')
sh.addPhysicalGroup(1,[l[1],l[3]],name='Clamped')
sh.addPhysicalGroup(1,[l[2]],name='Top')

# |--------------------------|
# |   Write the Mesh File    |
# |--------------------------|

sh.mesh.generate(2)
gmsh.write(os.path.dirname(__file__)+'/geometryS.msh')
gmsh.fltk.run()
gmsh.finalize()
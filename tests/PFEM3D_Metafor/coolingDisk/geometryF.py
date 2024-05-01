import os,gmsh
from gmsh import model as sh
gmsh.initialize()

# |---------------------------|
# |   Mesh Size Parameters    |
# |---------------------------|

L = 0.3
HF = 0.25
HS = 0.03
R = 0.025

d = HF/33
N = 13

# |----------------------------------|
# |   Points and Lines Definition    |
# |----------------------------------|

p = list()

p.append(sh.occ.addPoint(L,0,0,d))
p.append(sh.occ.addPoint(0,0,0,d))
p.append(sh.occ.addPoint(L,HF,0,d))
p.append(sh.occ.addPoint(0,HF,0,d))
p.append(sh.occ.addPoint(0,HF+HS+R,0,d))
p.append(sh.occ.addPoint(L,HF+HS+R,0,d))

p.append(sh.occ.addPoint(L/2,HF+HS,0,d))
p.append(sh.occ.addPoint(L/2,HF+HS+R,0,d))
p.append(sh.occ.addPoint(L/2,HF+HS-R,0,d))

# Lines List

l = list()
h = list()
r = list()
u = list()

l.append(sh.occ.addLine(p[0],p[1]))
l.append(sh.occ.addLine(p[0],p[2]))
l.append(sh.occ.addLine(p[2],p[3]))
l.append(sh.occ.addLine(p[3],p[1]))
l.append(sh.occ.addLine(p[3],p[4]))
l.append(sh.occ.addLine(p[2],p[5]))

h.append(sh.occ.addCircleArc(p[7],p[6],p[8]))
h.append(sh.occ.addCircleArc(p[8],p[6],p[7]))

# |------------------------------------|
# |   Physical Surface and Boundary    |
# |------------------------------------|

k = sh.occ.addCurveLoop(l[:4])
s = sh.occ.addPlaneSurface([k])
sh.occ.synchronize()

sh.mesh.setTransfiniteCurve(h[0],N)
sh.mesh.setTransfiniteCurve(h[1],N)

# Physical Boundary

sh.addPhysicalGroup(2,[s],name='Fluid')
sh.addPhysicalGroup(1,h+r+u,name='FSInterface')
sh.addPhysicalGroup(1,[l[2]],name='FreeSurface')
sh.addPhysicalGroup(1,l[:2]+l[3:],name='Wall')

# |--------------------------|
# |   Write the Mesh File    |
# |--------------------------|

sh.mesh.generate(2)
gmsh.model.mesh.reverse()
gmsh.write(os.path.dirname(__file__)+'/geometryF.msh')
gmsh.fltk.run()
gmsh.finalize()
#!/usr/bin/env python

import gmsh
import numpy as np

factory = gmsh.model.occ

gmsh.initialize()

gmsh.model.add("Cylinder")

gmsh.logger.start()

# # cyl = factory.addCylinder(0, 0, 0, 0, 0, 5, 1, tag=1)
cyl = factory.addBox(0, 0, 0, 1, 1, 2, tag=1)
dt_cyl = [(3,cyl)]

sph = factory.addSphere(1, 0, 1, 0.5)
dt_sph = [(3,sph)]

# disk = factory.addDisk(0,0,0, 1,1)
# dt_disk= [(2,disk)]


factory.synchronize()

gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 0.2)

b_cyl = gmsh.model.getBoundary(dt_cyl,True,False,False)
b_sph = gmsh.model.getBoundary(dt_sph,True,False,False)

print("b_cyl: ", b_cyl)
# print("b_cyl_fused: ", b_cyl_fused)
print("b_sph: ", b_sph)

# dt_int, _ = factory.fragment([(2,b_cyl_fused)], b_sph)
dt_int, _ = factory.fragment(b_cyl, b_sph)
# dt_int, _ = factory.cut(b_sph, b_cyl)
# dt_int, _ = factory.fragment(dt_disk, b_sph)
# dt_int, _ = factory.fragment(b_sph, dt_disk)

print("dt_int:", dt_int)
# dt_int = [ x for x in dt_int if x not in b_cyl ]
# dt_int = [ x for x in dt_int if x not in b_sph ]
# print("dt_int:", dt_int)

eps = 1e-3
dt_ent = factory.getEntitiesInBoundingBox(0-eps, 0-eps, 0-eps, 1+eps, 1+eps, 5+eps, dim=2)

print("dt_ent:", dt_ent)

## fragment and inside box
dt_new = list(set(dt_int).intersection(set(dt_ent)))
print("dt_new:", dt_new)

## fragment and inside box and not part of original surfaces
dt_newer = list(set(dt_new).difference(set(b_cyl)))
print("dt_newer:",dt_newer)


t_int = [t for _,t in dt_int]
t_ent = [t for _,t in dt_ent]
t_new = [t for _,t in dt_new]
t_newer = [t for _,t in dt_newer]



# print(t_int)

# import sys; sys.exit()

# factory.remove(dt_disk)
factory.remove(dt_cyl)
factory.remove(dt_sph)

factory.synchronize()

b_cyl_bbox = []

for d,t in b_cyl:
    b_cyl_bbox.append(gmsh.model.getBoundingBox(d,t))

print("")
print(b_cyl_bbox)
print("")

for dim, tag in dt_newer:
    rdt = gmsh.model.getBoundary([(dim,tag)],False,False,True)
    bbox = gmsh.model.getBoundingBox(dim,tag)
    bbox_cmp = []
    for _bbox in b_cyl_bbox:
        bbox_cmp.append(np.allclose(bbox, _bbox))
    print("tag:",tag, "cmp:", bbox_cmp)
    if any(bbox_cmp):
        print("skipping")
        continue
    for d, t in rdt:
        coord = []
        coord = gmsh.model.getValue(d, t, coord)
        pcoord = gmsh.model.getParametrization(dim,tag,coord)
        curv = gmsh.model.getCurvature(dim, tag, pcoord)
        if (curv == 0):
            print(dim, tag, curv, bbox)



# import sys; sys.exit()
# boundary_tags = [ y for _,y in boundaries]

# gmsh.model.addPhysicalGroup(2,[t_int[1]], 1)
# gmsh.model.addPhysicalGroup(2,t_newer, 1)
gmsh.model.addPhysicalGroup(2,[9,11], 1)
gmsh.model.setPhysicalName(2,1,"abc")

gmsh.option.setNumber('Print.GeoOnlyPhysicals', 1)

gmsh.model.mesh.generate(3)

gmsh.write("int.vtk")

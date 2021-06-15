#!/usr/bin/env python

import gmsh

gmsh.initialize()

gmsh.model.add("Cylinder")

gmsh.logger.start()

gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, 5, 1, tag=1)
dt_cyl = [(3,1)]

gmsh.model.occ.synchronize()

gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 0.2)

boundaries = gmsh.model.getBoundary(dt_cyl,False,False,False)

boundary_tags = [ y for _,y in boundaries]

gmsh.model.addPhysicalGroup(2, [boundary_tags[0]], 3)
gmsh.model.setPhysicalName(2,3,"Wall")

gmsh.model.addPhysicalGroup(2, [boundary_tags[1]], 2)
gmsh.model.setPhysicalName(2,2,"Outlet")

gmsh.model.addPhysicalGroup(2, [boundary_tags[2]], 1)
gmsh.model.setPhysicalName(2,1,"Inlet")

gmsh.model.addPhysicalGroup(3, [1], 4)
gmsh.model.setPhysicalName(3,4,"Volume")

gmsh.option.setNumber('Print.GeoOnlyPhysicals', 1)

gmsh.model.mesh.generate(3)

gmsh.write("cyl.msh2")

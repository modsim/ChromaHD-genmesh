#!/usr/bin/env python3

"""
@pre: export PYTHONPATH=$PYTHONPATH:/usr/local/lib64
"""

import gmsh

gmsh.initialize()

gmsh.option.setNumber("General.Terminal", 1)

gmsh.merge("../output/OB-0.03.msh2")
gmsh.plugin.setNumber("MeshVolume", "Physical", -1)
gmsh.plugin.setNumber("MeshVolume", "Dimension", 3)
gmsh.plugin.run("MeshVolume")

gmsh.plugin.setNumber("MeshVolume", "Physical", 5)
gmsh.plugin.setNumber("MeshVolume", "Dimension", 3)
gmsh.plugin.run("MeshVolume")

gmsh.plugin.setNumber("MeshVolume", "Physical", 6)
gmsh.plugin.setNumber("MeshVolume", "Dimension", 3)
gmsh.plugin.run("MeshVolume")


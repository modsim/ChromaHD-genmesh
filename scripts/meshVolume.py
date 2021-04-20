#!/usr/bin/env python3

"""
@pre: export PYTHONPATH=$PYTHONPATH:/usr/local/lib64
or
@pre: export PYTHONPATH=$PYTHONPATH:/usr/local/lib
"""

import gmsh
import sys

gmsh.initialize()

gmsh.option.setNumber("General.Terminal", 1)

gmsh.merge(sys.argv[1])

print("[Column]")
gmsh.plugin.setNumber("MeshVolume", "PhysicalGroup", -1)
gmsh.plugin.setNumber("MeshVolume", "Dimension", 3)
gmsh.plugin.run("MeshVolume")

print("[Interstitial]")
gmsh.plugin.setNumber("MeshVolume", "PhysicalGroup", 5)
gmsh.plugin.setNumber("MeshVolume", "Dimension", 3)
gmsh.plugin.run("MeshVolume")

print("[Packed Bed]")
gmsh.plugin.setNumber("MeshVolume", "PhysicalGroup", 6)
gmsh.plugin.setNumber("MeshVolume", "Dimension", 3)
gmsh.plugin.run("MeshVolume")



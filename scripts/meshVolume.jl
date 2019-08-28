push!(LOAD_PATH, "/usr/local/lib64")
import gmsh

gmsh.initialize(ARGS)

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


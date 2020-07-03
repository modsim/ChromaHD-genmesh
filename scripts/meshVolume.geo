Merge "../output/OB-0.03.msh2" ;

Plugin(MeshVolume).Dimension = 3 ;
Plugin(MeshVolume).PhysicalGroup = -1 ;
Plugin(MeshVolume).Run ;

Plugin(MeshVolume).Dimension = 3 ;
Plugin(MeshVolume).PhysicalGroup = 5 ;
Plugin(MeshVolume).Run ;

Plugin(MeshVolume).PhysicalGroup = 6 ;
Plugin(MeshVolume).Dimension = 3 ;
Plugin(MeshVolume).Run ;

Merge "../output/OB-0.03.msh2" ;

Plugin(MeshVolume).Dimension = 3 ;
Plugin(MeshVolume).Physical = 5 ;
Plugin(MeshVolume).Run ;

Plugin(MeshVolume).Physical = 6 ;
Plugin(MeshVolume).Dimension = 3 ;
Plugin(MeshVolume).Run ;

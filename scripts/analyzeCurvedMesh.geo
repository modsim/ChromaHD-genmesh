Merge "../output/VC-mono-0.04-bri/VC-mono-0.04-bri.msh2" ;                                                                                                                   

Plugin(AnalyseCurvedMesh).JacobianDeterminant = 1 ;
Plugin(AnalyseCurvedMesh).IGEMeasure = 1 ;
Plugin(AnalyseCurvedMesh).ICNMeasure = 1 ;
Plugin(AnalyseCurvedMesh).DimensionOfElements = 3 ;

Plugin(AnalyseCurvedMesh).Run ;

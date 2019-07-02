all: 
	g++ Parameters.cpp PackedBed.cpp Files.cpp Bead.cpp main.cpp -L/usr/local/lib64 -lgmsh -o mesher  -O3 

# ibt:
# 	g++ Parameters.cpp PackedBed.cpp Files.cpp Bead.cpp main.cpp -L../../../tools/gmsh/lib -lgmsh -o mesher -O3 -I../../../tools/gmsh/include -std=c++11 -L../../../tools/occt/lib -lgmp -lTKXDESTEP -lTKXDEIGES -lTKXCAF -lTKLCAF -lTKVCAF -lTKCAF -lTKV3d -lTKService -lTKCDF -lfreetype -lTKSTEP -lTKSTEP209 -lTKSTEPAttr -lTKSTEPBase -lTKIGES -lTKXSBase -lTKOffset -lTKFeat -lTKFillet -lTKBool -lTKMesh -lTKHLR -lTKBO -lTKPrim -lTKShHealing -lTKTopAlgo -lTKGeomAlgo -lTKBRep -lTKGeomBase -lTKG3d -lTKG2d -lTKMath -lTKernel -ldl -lblas -lpthread -lrt -lgfortran

ibt:
	g++ Parameters.cpp PackedBed.cpp Files.cpp Bead.cpp main.cpp -L../../../tools/gmsh/lib -lgmsh -o mesher -O3 -I../../../tools/gmsh/include -std=c++11 -L../../../tools/occt/lib

all: git-check
	g++ Parameters.cpp PackedBed.cpp Model.cpp Files.cpp Bead.cpp main.cpp -L/usr/local/lib64 -lgmsh -o mesher  -O3 

ibt: git-check
	g++ Parameters.cpp PackedBed.cpp Model.cpp Files.cpp Bead.cpp main.cpp -L../../../tools/gmsh/lib -lgmsh -o mesher -O3 -I../../../tools/gmsh/include -std=c++11 -L../../../tools/occt/lib

git-check:
	@echo "#ifndef VERSION_H" > version.h
	@echo "#define VERSION_H" >>version.h
	@echo  >> version.h
	@echo "#define GITCOMMIT \"$$(git rev-parse --verify HEAD)\"" >>version.h
	@echo "#define GITSTATE \"$$(git diff --quiet || echo '(dirty!)')\""  >>version.h
	@echo  >> version.h
	@echo "#endif /* VERSION_H */" >>version.h

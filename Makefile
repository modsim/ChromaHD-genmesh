HOST=$(shell hostname)
FLAGS=-O3
CXX := $(if $(CXX),$(CXX),g++)

PREFIX=$(HOME)/local

SRC=Parameters.cpp PackedBed.cpp Geometry.cpp Model.cpp Files.cpp Bead.cpp main.cpp 
OBJ=$(SRC:.cpp=.o)
EXE=genmesh

## Set GMSH_ROOT and OCCT_ROOT before running make
INC=-I${GMSH_ROOT}/include
LIB=-L${GMSH_ROOT}/lib -lgmsh -L${OCCT_ROOT}/lib

.PHONY: all local ibt clean install

all: build 

clean: 
	@rm -rf *.o $(EXE)

build: git-check $(OBJ)
	$(CXX)  $(OBJ) $(FLAGS) $(INC) $(LIB)  -o $(EXE)
	@echo "Compiled $(EXE) !" 

%.o: %.cpp
	$(CXX) -c -o $@ $< $(FLAGS) $(INC)

git-check:
	@echo "#ifndef VERSION_H" > version.h
	@echo "#define VERSION_H" >>version.h
	@echo  >>version.h
	@echo "#define GITCOMMIT \"$$(git rev-parse --verify HEAD)\"" >>version.h
	@echo "#define GITSTATE \"$$(git diff --quiet || echo '(dirty!)')\""  >>version.h
	@echo "#define GMSHVERSION \"$$($$GMSH_ROOT/bin/gmsh --version 2>&1 > /dev/null)\"" >>version.h
	@echo  >>version.h
	@echo "#endif /* VERSION_H */" >>version.h

install: build
	cp genmesh $(PREFIX)/bin	

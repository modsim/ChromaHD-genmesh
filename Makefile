HOST=$(shell hostname)
CC=g++
FLAGS=-O3

SRC=Parameters.cpp PackedBed.cpp Model.cpp Files.cpp Bead.cpp main.cpp
OBJ=$(SRC:.cpp=.o)
EXE=genmesh

ifeq (IBT918, $(findstring IBT918, $(HOST)))
	INC=
	LIB=-lgmsh -L/usr/local/lib64
else
	INC=-I../../../tools/gmsh/include
	LIB=-L../../../tools/gmsh/lib -lgmsh -L../../../tools/occt/lib
endif

.PHONY: all local ibt clean

all: build 

clean: 
	@rm -rf *.o $(EXE)

build: git-check $(OBJ)
	$(CC)  $(OBJ) $(FLAGS) $(INC) $(LIB)  -o $(EXE)
	@echo "Compiled $(EXE) !" 

%.o: %.cpp
	$(CC) -c -o $@ $< $(FLAGS) $(INC)

git-check:
	@echo "#ifndef VERSION_H" > version.h
	@echo "#define VERSION_H" >>version.h
	@echo  >> version.h
	@echo "#define GITCOMMIT \"$$(git rev-parse --verify HEAD)\"" >>version.h
	@echo "#define GITSTATE \"$$(git diff --quiet || echo '(dirty!)')\""  >>version.h
	@echo  >> version.h
	@echo "#endif /* VERSION_H */" >>version.h

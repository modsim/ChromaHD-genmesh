#!/bin/bash

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
## This is an imperfect script and must be used mostly as a GUIDE to installing the dependencies for genmesh. 
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

set -e

NTHREADS=48
BASE_DIR=${PWD}

TCL_LINK="https://prdownloads.sourceforge.net/tcl/tcl8.6.9-src.tar.gz"
TK_LINK="https://prdownloads.sourceforge.net/tcl/tk8.6.9.1-src.tar.gz"
OCCT_LINK="https://www.opencascade.com/sites/default/files/private/occt/OCC_7.3.0_release/opencascade-7.3.0.tgz"

install_tcl()
{
    wget TCL_LINK -O tcl-src.tar.gz
    tar xf tcl-src.tar.gz
    cd tcl8.6.9/unix
    ./configure --enable-gcc --enable-shared --enable-thrds --prefix=$BASE_DIR/tcl --enable-64bit
    make -j $NTHREADS
    make install
    cd -
}

install_tk()
{
    wget TK_LINK -O tk-src.tar.gz
    tar xf tk-src.tar.gz
    cd tk8.6.9/unix
    ./configure --enable-gcc --enable-shared --enable-thrds --prefix=$BASE_DIR/tk --enable-64bit --with-tcl=$BASE_DIR/tcl/lib                           
    make -j $NTHREADS
    make install
    cd -
}

install_occt()
{
    wget OCCT_LINK -O occt-src.tar.gz
    tar xf occt-src.tar.gz
    mkdir -p opencascade-7.3.0/build
    cd opencascade-7.3.0/build
    cmake -DCMAKE_INSTALL_PREFIX=$BASE_DIR/occt -D3RDPARTY_TCL_LIBRARY_DIR=$BASE_DIR/tcl/lib/ -D3RDPARTY_TK_LIBRARY_DIR=$BASE_DIR/tk/lib -D3RDPARTY_TK_INCLUDE_DIR=$BASE_DIR/tk/include/ -D3RDPARTY_TCL_INCLUDE_DIR=$BASE_DIR/tcl/include/ ..
    make -j $NTHREADS
    make install
    cd -

}

install_gmsh()
{
    git clone https://gitlab.onelab.info/gmsh/gmsh.git gmsh-git
    mkdir gmsh-git/build
    cd gmsh-git/build
    cmake -DCMAKE_PREFIX_PATH=$BASE_DIR/occt -DCMAKE_INSTALL_PREFIX=$BASE_DIR/gmsh-git-installed ..
    make -j $NTHREADS
    make install
    cd -
    ln -s $BASE_DIR/gmsh-git-installed $BASE_DIR/gmsh
}

check_pre()
{
    #TODO: check headers also
    #TODO: exit if not found
    [[ $(ldconfig -p | grep libXmu.so) -gt 1 ]] && echo "Found libXmu"
    [[ $(ldconfig -p | grep libXi.so) -gt 1 ]] && echo "Found libXi"
    [[ $(ldconfig -p | grep libfreetype.so) -gt 1 ]] && echo "Found libfreetype"

}

check_pre
install_tcl
install_tk
install_occt
install_gmsh

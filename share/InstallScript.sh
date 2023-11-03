#!/bin/bash

#------------------------------------------------------------

# set spack and moris installation directory
export WORKSPACE=$HOME/codes_new

# define developper mode: 0 for users; 1 for developpers
export DEVELOPPER_MODE=0

# define blas implementation INTEL_MKL: 0 for no; 1 for yes
export INTEL_MKL=0

# set compiler and version
export COMPILER='gcc@11.4.0'

# set directory for temporary files during built
export TMPDIR=/tmp

#------------------------------------------------------------

if [ -d $WORKSPACE ];then
    echo "WORKSPACE $WORKSPACE already exits; remove or rename it first before installing moris"
    exit
fi

mkdir $WORKSPACE
cd $WORKSPACE

#------------------------------------------------------------

# get spack

git clone https://github.com/spack/spack.git

#------------------------------------------------------------

. $WORKSPACE/spack/share/spack/setup-env.sh

#------------------------------------------------------------

if [ $DEVELOPPER_MODE = "1" ];then
    git clone git@github.com:kkmaute/moris
else
    git clone https://github.com/kkmaute/moris
fi

#------------------------------------------------------------

spack create --name moris --skip-editor 
spack create --name gcmma --skip-editor 
spack create --name lbfgs --skip-editor 
spack create --name snopt --skip-editor 

#------------------------------------------------------------

cp $WORKSPACE/moris/share/spack/trilinos_package.py  $WORKSPACE/spack/var/spack/repos/builtin/packages/trilinos/package.py
cp $WORKSPACE/moris/share/spack/moris_package.py     $WORKSPACE/spack/var/spack/repos/builtin/packages/moris/package.py
cp $WORKSPACE/moris/share/spack/gcmma_package.py     $WORKSPACE/spack/var/spack/repos/builtin/packages/gcmma/package.py
cp $WORKSPACE/moris/share/spack/lbfgs_package.py     $WORKSPACE/spack/var/spack/repos/builtin/packages/lbfgs/package.py
cp $WORKSPACE/moris/share/spack/snopt_package.py     $WORKSPACE/spack/var/spack/repos/builtin/packages/snopt/package.py

#------------------------------------------------------------

spack env create -d .

spack env activate .

#------------------------------------------------------------

spack compiler find

#------------------------------------------------------------

if [ $DEVELOPPER_MODE = "1" ];then
    spack add moris+pardiso+mumps
else
    spack add moris+openblas~petsc~slepc~pardiso~mumps~gcmma~lbfgs~snopt
fi

spack develop --path $WORKSPACE/moris moris@main

if [ $INTEL_MKL = "0" ];then
    spack add openblas
fi

if [ $DEVELOPPER_MODE = "1" ];then
    spack add doxygen
    spack add llvm~gold
fi

spack add openmpi fabrics=auto 

spack add python 

#------------------------------------------------------------

sed -i -e 's/unify: true/unify: when_possible/g' ./spack.yaml

#------------------------------------------------------------

spack concretize -f -U

#------------------------------------------------------------

spack install python %"$COMPILER"

spack install openmpi %"$COMPILER"

if [ $INTEL_MKL = "0" ];then
    spack install openblas %"$COMPILER"
fi

#------------------------------------------------------------

if [ $DEVELOPPER_MODE = "1" ];then
    spack install --only dependencies moris %"$COMPILER"
    
    spack install doxygen %"$COMPILER"

    spack install llvm %"$COMPILER"
else
    spack install moris %"$COMPILER"
fi

#------------------------------------------------------------

spack clean

#------------------------------------------------------------

if [ $DEVELOPPER_MODE = "1" ];then
    tcsh $WORKSPACE/moris/share/spack/make_moris_cshrc.sh
else
    tcsh $WORKSPACE/moris/share/spack/make_moris_cshrc.sh view
fi

sed -rn 's/^\s*setenv\s+(\S+)\s+/export \1=/p' ~/.cshrc_moris > ~/.bashrc_moris


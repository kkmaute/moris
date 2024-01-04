#!/bin/bash
#------------------------------------------------------------

# Installation script for MORIS 

# tested on: ubuntu; OpenSUSE 15.5 

# For special cases, such as building openmpi with scheduler, 
# see moris/share/isntall/Install_Step_By_Step.txt

#------------------------------------------------------------
# edit the following lines
#------------------------------------------------------------

# set spack and moris installation directory
export WORKSPACE=$HOME/codes

# define developper mode: 0 for users; 1 for developpers
export DEVELOPPER_MODE=0

# define blas and lapack provider:
# options are:   amd, intel-mkl, intel-oneapi-mkl, openblas, netlib
export BLASLAPACK=openblas

# define whether petsc and slepc should be installed: 0 for no; 1 for yes
export PETSC=0

# define whether pardiso solver should be installed: 0 for no; 1 for yes
export PARDISO=0

# define whether mumps solver should be installed: 0 for no; 1 for yes
export MUMPS=0

# define wether optimization should be enabled: 0 for no; 1 for yes
# note: optimization requires access to github optimization packages
#       under gitbub.com/kkmaute: gcmma, lbfgs,and snopt
export OPTALG=0

# set compiler and version (see out of gcc --version)
export COMPILER='gcc@11.3.0'

# set number of processors used for installation (0: automatically detected)
# use lower number than available processors if build process runs out of memory
export NUMPROC=0

# set directory for temporary files during built
export TMPDIR=/tmp

# verbose mode for spack builds: 0 for no; 1 for yes
export VERBOSE=0

#------------------------------------------------------------
# no need to edit script below
#------------------------------------------------------------

echo ""                                      > moris_config.log
echo "MORIS installation parameters:"       >> moris_config.log
echo ""                                     >> moris_config.log
echo `date`                                 >> moris_config.log
echo ""                                     >> moris_config.log
echo "WORKSPACE          $WORKSPACE"        >> moris_config.log
echo "DEVELOPPER_MODE    $DEVELOPPER_MODE"  >> moris_config.log
echo "BLASLAPACK         $BLASLAPACK"       >> moris_config.log
echo "PETSC              $PETSC"            >> moris_config.log
echo "PARDISO            $PARDISO"          >> moris_config.log
echo "MUMPS              $MUMPS"            >> moris_config.log
echo "OPTALG             $OPTALG"           >> moris_config.log
echo "COMPILER           $COMPILER"         >> moris_config.log
echo "NUMPROC            $NUMPROC"          >> moris_config.log
echo "TMPDIR             $TMPDIR"           >> moris_config.log
echo "VERBOSE            $VERBOSE"          >> moris_config.log
echo ""

if [ "$VERBOSE" = "1" ];then
    cat moris_config.log
    echo "Press return to continue; to abort press ctrl+c"
    echo ""
    read ans
    echo ""
fi

ORGDIR=`pwd`

#------------------------------------------------------------

export blaspro=na

if [ $PETSC = "1" ];then
    export petopt='+petsc+slepc'
else
    export petopt='~petsc~slepc'
fi

if [ $PARDISO = "1" ];then
    export mklpro=intel-oneapi-mkl
    export paropt='+pardiso'
else
    export paropt='~pardiso'
fi

if [ $MUMPS = "1" ];then
    export mumopt='+mumps'
else
    export mumopt='~mumps'
fi

if [ $OPTALG = "1" ];then
    export optopt='+gcmma+lbfgs+snopt'
else
    export optopt='~gcmma~lbfgs~snopt'
fi

#------------------------------------------------------------

if [ $BLASLAPACK = "amd" ];then
    export blaspro=amdblis
    export lapackpro=amdlibflame
    export sclpackpro=amdscalapack
fi

if [ $BLASLAPACK = "intel-mkl" ];then
    export blaspro=intel-mkl
    export lapackpro=intel-mkl
    export sclpackpro=intel-mkl
    export mklpro=intel-mkl
fi

if [ $BLASLAPACK = "intel-oneapi-mkl" ];then
    export blaspro=intel-oneapi-mkl
    export lapackpro=intel-oneapi-mkl
    export sclpackpro=intel-oneapi-mkl
    export mklpro=intel-oneapi-mkl
fi

if [ $BLASLAPACK = "openblas" ];then
    export blaspro=openblas
    export lapackpro=openblas
    export sclpackpro=netlib-scalapack
fi

if [ $BLASLAPACK = "netlib" ];then
    export blaspro=netlib-lapack
    export lapackpro=netlib-lapack
    export sclpackpro=netlib-scalapack
fi

if [ $blaspro = "na" ];then
    echo "Error - incorrect blas and lapack provider"
    exit
fi

#------------------------------------------------------------

if [ -d $WORKSPACE ];then
    echo "WORKSPACE $WORKSPACE already exits; remove or rename it first before installing moris"
    exit
fi

mkdir $WORKSPACE
cd $WORKSPACE

mv $ORGDIR/moris_config.log .

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

spack compiler list | grep gcc | grep '@' | awk -v compiler=$COMPILER '{ if ( match($0,compiler) == 0) {cmd="spack compiler rm "$0; system(cmd)}}'

spack compiler list 

ret=`spack compiler list | grep gcc | grep '@' | awk -v compiler=$COMPILER 'BEGIN{n=0}{ if ( match($0,compiler) > 0) {n=n+1}}END{print n}'`

if [ ! $ret = "1" ];then
    echo "Error - $COMPILER not available"
    exit
fi
    
#------------------------------------------------------------

echo "  packages:"                      >> spack.yaml
echo "    all:"                         >> spack.yaml
echo "      providers:"                 >> spack.yaml
echo "        blas: [$blaspro]"         >> spack.yaml
echo "        lapack: [$lapackpro]"     >> spack.yaml
echo "        scalapack: [$sclpackpro]" >> spack.yaml
echo "        mkl: [$mklpro]"           >> spack.yaml

#------------------------------------------------------------

spack add moris$petopt$paropt$mumopt$optopt %"$COMPILER"

spack develop --path $WORKSPACE/moris moris@main

if [ $DEVELOPPER_MODE = "1" ];then
    spack add doxygen %"$COMPILER"
    spack add llvm~gold %"$COMPILER"
fi

spack add openmpi %"$COMPILER" fabrics=auto 

spack add python %"$COMPILER"

#------------------------------------------------------------

sed -i -e 's/unify: true/unify: when_possible/g' ./spack.yaml

#------------------------------------------------------------

spack concretize -f -U >> moris_config.log

#------------------------------------------------------------

SOPTION=""
if [ $VERBOSE = "1" ];then
   SOPTION="-v"
fi

if [ ! $NUMPROC = "0" ];then
  SOPTION="$SOPTION -j $NUMPROC"
fi

#------------------------------------------------------------

spack install $SOPTION python %"$COMPILER"

#------------------------------------------------------------

isOpenSUSE=`grep NAME /etc/os-release | head -1 | awk -F '=' '{ if ( match($2,"openSUSE") > 1 ) {print 1}else{print 0}}'`

if [ $isOpenSUSE = "1" ];then

  echo "fixing python installation for OpenSUSE"

  export PYIDIR=`find spack/opt/spack/ -type d -name "python-*"`

  export PYLVERS=`ls $PYIDIR/include`
  
  cp -R $PYIDIR/lib64/$PYLVERS/lib-dynload $PYIDIR/lib/$PYLVERS/.
fi

#------------------------------------------------------------

spack install $SOPTION openmpi %"$COMPILER"

#------------------------------------------------------------

if [ $DEVELOPPER_MODE = "1" ];then
    spack install $SOPTION --only dependencies moris %"$COMPILER"
    
    spack install $SOPTION doxygen %"$COMPILER"

    spack install $SOPTION llvm %"$COMPILER"
else
    spack install $SOPTION moris %"$COMPILER"
fi

#------------------------------------------------------------

spack clean

#------------------------------------------------------------

if [ $DEVELOPPER_MODE = "1" ];then
    bash $WORKSPACE/moris/share/spack/make_moris_resource.sh
else
    bash $WORKSPACE/moris/share/spack/make_moris_resource.sh view
fi

#------------------------------------------------------------

if [ $DEVELOPPER_MODE = "1" ];then
    echo ""
    echo ""
    echo "Developers: follow build instructions under moris/share/install/Build_Instructions.txt"
    echo "            to build moris in either debug or release (opt) mode"
    echo ""
    echo ""
fi

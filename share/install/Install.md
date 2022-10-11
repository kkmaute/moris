# Installation

*------------------------------------------------------------
* The following document describes the installation of MORIS. This installation
* is based on spack.
*
* For installation of MORIS on a standard workstation under Linux, follow the
* instructions below. These instructions have been tested for OPENSuse 15.x.
*
* For installation of MORIS on a cluster system, see Install_Cluter.md
*
* Note: The installation process uses tcsh
*------------------------------------------------------------

*------------------------------------------------------------
* Prerequisites on OS installation:
*
* The following packages should be installed; 
* if install them with space (see option below)
*
* patch 
* pkgconf (or pkg-config)
* makeinfo (part of texinfo package)
*
* To check whether the above packages are installed, excute 
* the command:  

which patch
which pkgconf
which pkg-config
which makeinfo
*------------------------------------------------------------

*------------------------------------------------------------
* For users with an existing MORIS installation, it is 
* strongly recommended to remove all environment variables
* related to MORIS from their environment source files (i.e.,
* .cshrc) before performing the MORIS installation described 
* below.
*
*------------------------------------------------------------

*------------------------------------------------------------
* basic setup
*
* MORIS and spack along with third party libraries will be installed 
* in in a workspace, e.g. $HOME/codes
*------------------------------------------------------------

* set workspace; change directory name as needed

setenv WORKSPACE $HOME/codes

*------------------------------------------------------------

* create and enter workspace

mkdir $WORKSPACE
cd $WORKSPACE

*------------------------------------------------------------

* get spack

git clone https://github.com/spack/spack.git

*------------------------------------------------------------

* set spack root directory and add to PATH variable

setenv SPACK_ROOT $WORKSPACE/spack
setenv PATH $PATH/:$SPACK_ROOT/bin

*------------------------------------------------------------

* source spack environment

source $SPACK_ROOT/share/spack/setup-env.csh

*------------------------------------------------------------

* let spack find installed compilers; see $HOME/.spack/linux/compilers.yaml

spack compiler find

*------------------------------------------------------------

* check what compilers have been found and remove unwanted;
* recommendation remove all but the one compiler you want to build moris with

spack compiler list

spack compiler rm <unwanted compiler>

*------------------------------------------------------------

* let spack find installed external package
* after commend has been exectued, see $HOME/.spack/packages.yaml

spack external find

*------------------------------------------------------------

* if texinfo is installed through your OS but makeinfo is not found
* remove texinfo from $HOME/.spack/packages.yaml

*------------------------------------------------------------

* get moris (main branch)

git clone git@github.com:kkmaute/moris 

*------------------------------------------------------------

* add moris specific packages to spack; ignore warnings

spack create --name gcmma --skip-editor
spack create --name lbfgs --skip-editor
spack create --name snopt --skip-editor
spack create --name moris --skip-editor

*------------------------------------------------------------

* add package files; overwrite trilinos package file to enable mkl and pardiso

cp $WORKSPACE/moris/share/spack/gcmma_package.py     $WORKSPACE/spack/var/spack/repos/builtin/packages/gcmma/package.py
cp $WORKSPACE/moris/share/spack/lbfgs_package.py     $WORKSPACE/spack/var/spack/repos/builtin/packages/lbfgs/package.py
cp $WORKSPACE/moris/share/spack/snopt_package.py     $WORKSPACE/spack/var/spack/repos/builtin/packages/snopt/package.py
cp $WORKSPACE/moris/share/spack/trilinos_package.py  $WORKSPACE/spack/var/spack/repos/builtin/packages/trilinos/package.py
cp $WORKSPACE/moris/share/spack/moris_package.py     $WORKSPACE/spack/var/spack/repos/builtin/packages/moris/package.py

*------------------------------------------------------------

* create spack environment and get spack ready to install moris

spack env create -d .

spack env activate .

spack add moris

spack develop --path $WORKSPACE/moris moris@main

spack add openblas

*------------------------------------------------------------

* if patch, pkgconf, and/or texinfo are not installed on your OS

spack add patch
spack add pkgconf
spack add texinfo

*------------------------------------------------------------

* finalize installation configuration

spack concretize -f -U

*------------------------------------------------------------

* if patch, pkgconf, and/or texinfo are not installed on your OS

spack install patch
spack install pkgconf
spack install texinfo

*------------------------------------------------------------

* install moris dependencies - by default in debug version

spack install --only dependencies moris

spack install openblas

*------------------------------------------------------------

* remove left over build directories

spack clean

*------------------------------------------------------------

* create the following resource file and source it as part of your .cshrc

tcsh $WORKSPACE/moris/share/spack/make_moris_cshrc.sh
        
*------------------------------------------------------------
*
* Important: .cshrc_moris needs to be sourced when working with moris
*
* Recommendation: include it in your .cshrc file
*
*------------------------------------------------------------

* build moris optimized version and run tests

source $HOME/.cshrc_moris

cd $MORISROOT

mkdir build_opt

cd build_opt

cmake -DBUILD_ALL=ON -DMORIS_USE_EXAMPLES=ON ..
cmake -DBUILD_ALL=ON -DMORIS_USE_EXAMPLES=ON ..

make -j 16

ctest
                                                                           
*------------------------------------------------------------

* build opt version and run test

source $HOME/.cshrc_moris

cd $MORISROOT

mkdir build_dbg

cd build_dbg

cmake -DBUILD_ALL=ON -DMORIS_HAVE_DEBUG=ON -DMORIS_HAVE_SYMBOLIC=ON -DMORIS_HAVE_SYMBOLIC_STRONG=ON -DMORIS_USE_EXAMPLES=ON ..
cmake -DBUILD_ALL=ON -DMORIS_HAVE_DEBUG=ON -DMORIS_HAVE_SYMBOLIC=ON -DMORIS_HAVE_SYMBOLIC_STRONG=ON -DMORIS_USE_EXAMPLES=ON ..

make -j 16

ctest


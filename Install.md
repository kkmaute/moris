# Installation

*------------------------------------------------------------
* The following document describes the installation of MORIS. This installation
* is based on spack.
*
* For installation of MORIS on a standard workstation under Linux, follow the
* instructions below.
*
* For installation of MORIS on a cluster system, see Install_Cluter.md
*
* Note: The installation process uses tcsh
*------------------------------------------------------------

*------------------------------------------------------------
* Prerequisites on OS installation:
*
* The following packages should be installed; if not spack will do it
*
* Modules 
* patch 
* makeinfo
* pkgconf
*
* sudo: zypper in Modules patch makeinfo
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
* if you do not have root access, enter <ctrl>+d when asked for root password
* after commend has been exectued, see $HOME/.spack/packages.yaml

spack external find --all

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

spack concretize -f

*------------------------------------------------------------

* install moris dependencies - by default in debug version

spack install --only dependencies moris

spack install openblas
*
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

source $HOME/.cshrc_moris

*------------------------------------------------------------

* build moris optimized version and run tests
* you may have to set the compilers manually (CC,CXX,F77,FC)

cd $WORKSPACE/moris

mkdir build_opt

cd build_opt

cmake -DBUILD_ALL=ON -DMORIS_USE_EXAMPLES=ON ..

make -j 16

ctest
                                                                           
*------------------------------------------------------------

* build opt version and run test

source $HOME/.cshrc_moris

cd $WORKSPACE/moris

mkdir build_dbg

cd build_dbg

cmake -DBUILD_ALL=ON -DMORIS_HAVE_DEBUG=ON -DMORIS_HAVE_SYMBOLIC=ON -DMORIS_HAVE_SYMBOLIC_STRONG=ON -DMORIS_USE_EXAMPLES=ON ..

make -j 16

ctest


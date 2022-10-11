# Installation on clusters using gcc compilers

*------------------------------------------------------------
* The following document describes the installation of MORIS. This installation
* is based on spack.
*
* This doucment describes the installation of MORIS on a HPC cluster system;
* using gcc; for installing MORIS on a cluster using an intel complier, 
* see Install_Cluster_intel.md ; for installing MORIS on a workstation with 
* Linux OS, see Install.md.
*
* Note: The installation process uses tcsh
*------------------------------------------------------------

*------------------------------------------------------------
* Prerequisites on OS installation:
*
* Load the modules for compiler, mpi implementation etc. 
* that you want to use for your installation, e.g.
* 
* gcc/10.2.0
* openmpi/4.0.5
* perl/5.24.0
* zlib/1.2.11 
* szip/2.1.1
* lmod/6.3.7 
* autotools/2.71 
* cmake/3.20.2   
*------------------------------------------------------------

* Check that the following software is installed; if it is not
* installed use the steps defined below
*
* patch 
* pkgconf
* makeinfo (part of texinfo package)
*
* To check whether the above packages are installed, excute 
* the command:  which patch; which pkgconf; which makeinfo
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

* optional: set tmp directory

setenv TMP /home/maute/work/tmp
rm -r -f $TMP
mkdir $TMP

*------------------------------------------------------------

module load gcc/10.2.0
module load openmpi/4.0.5
module load perl/5.24.0
module load zlib/1.2.11 
module load szip/2.1.1
module load lmod/6.3.7 
module load autotools/2.71 
module load cmake/3.20.2   

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

spack external find --all

*------------------------------------------------------------

* if texinfo is installed through your OS but makeinfo is not found
* remove texinfo from $HOME/.spack/packages.yaml

*------------------------------------------------------------

* you may also have to remove the following packages from 
* ~/.spack/packages.yaml and spack.yaml  as the pre-installed 
* version may cause issues

* hwloc                   as external version causes issues in trilinos
* libtool (older version) as it causes issues in mpfr
* autoconf (older version)

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

*add the lines below to 

  spack/var/spack/repos/builtin/packages/gmp/package.py
  spack/var/spack/repos/builtin/packages/mpfr/package.py

    def autoreconf(self, spec, prefix):
        sh = which('sh')
        sh('-c', 'autoreconf -fi')

*------------------------------------------------------------

* create spack environment and get spack ready to install moris

spack env create -d .

spack env activate .

spack add moris

spack develop --path $WORKSPACE/moris moris@main

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

*------------------------------------------------------------

* remove left over build directories

spack clean

*------------------------------------------------------------

* create the following resource file and source it as part of your .cshrc

tcsh $WORKSPACE/moris/share/spack/make_moris_cluster_gcc_cshrc.sh
        
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

cmake -DMORIS_USE_MKL=ON -DMORIS_USE_OPENBLAS=OFF  -DBUILD_ALL=ON -DMORIS_USE_EXAMPLES=ON ..
cmake -DMORIS_USE_MKL=ON -DMORIS_USE_OPENBLAS=OFF  -DBUILD_ALL=ON -DMORIS_USE_EXAMPLES=ON ..

make -j 16

ctest
                                                                           
*------------------------------------------------------------

* build opt version and run test

source $HOME/.cshrc_moris

cd $WORKSPACE/moris

mkdir build_dbg

cd build_dbg

cmake -DMORIS_USE_MKL=ON -DMORIS_USE_OPENBLAS=OFF -DBUILD_ALL=ON -DMORIS_HAVE_DEBUG=ON -DMORIS_HAVE_SYMBOLIC=ON -DMORIS_HAVE_SYMBOLIC_STRONG=ON -DMORIS_USE_EXAMPLES=ON ..
cmake -DMORIS_USE_MKL=ON -DMORIS_USE_OPENBLAS=OFF -DBUILD_ALL=ON -DMORIS_HAVE_DEBUG=ON -DMORIS_HAVE_SYMBOLIC=ON -DMORIS_HAVE_SYMBOLIC_STRONG=ON -DMORIS_USE_EXAMPLES=ON ..

make -j 16

ctest


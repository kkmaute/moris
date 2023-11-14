# Installation on clusters using intel compilers

*------------------------------------------------------------
* The following document describes the installation of MORIS. This installation
* is based on spack.
*
* This doucment describes the installation of MORIS on a HPC cluster system using 
* an intel compiler; for installing MORIS on a cluster using gcc, 
* see Install_Cluster_gcc.md ; for installing MORIS on a workstation with 
* Linux OS, see Install.md.;
*
* Note: The installation process uses tcsh
*------------------------------------------------------------

*------------------------------------------------------------
* Prerequisites on OS installation:
*
* Load the modules for compiler, mpi implementation etc. 
* that you want to use for your installation, e.g.

module load intel/20.2 
module load impi/19.8
module load mkl/20.2
module load perl/5.24.0
module load zlib/1.2.11 
module load szip/2.1.1
module load lmod/6.3.7 
module load autotools/2.71 
module load cmake/3.20.2   

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

* as needed to avoid compilation issues: edit file   and add 
* entries under modules and LD_LIBRARY_PATH
 
    modules: [intel/20.2, impi/19.8]
    
    environment: 
        append_path:
            LD_LIBRARY_PATH: /curc/sw/intel/20.2/compilers_and_libraries_2020.2.254/linux/compiler/lib/intel64/:/usr/lib64/psm2-compat/

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

*to force the use of machine stalled mpi and mkl, 
* add the following lines to ~/.spack/packages 

  mpi:
    buildable: False
    require:
    - one_of: [
      intel-mpi
    ]
  intel-mpi:
    externals:
    - spec: intel-mpi@19.8
      prefix: /curc/sw/intel/20.2/compilers_and_libraries_2020.2.254/linux/mpi/intel64/
    buildable: False
  blas:
    buildable: False
    require:
    - one_of: [
      intel-mkl
    ]
  lapack:
    buildable: False
    require:
    - one_of: [
      intel-mkl
    ]
  mkl:
    buildable: False
    require:
    - one_of: [
      intel-mkl
    ]    
  intel-mkl:
    externals:
    - spec: intel-mkl@20.2
      prefix: /curc/sw/intel/20.2/compilers_and_libraries_2020.2.254/linux/mkl
    buildable: False

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

tcsh $WORKSPACE/moris/share/spack/make_cluster_moris_cshrc.sh
        
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

cmake -DMORIS_USE_INTEL=ON -DMORIS_USE_MKL=ON -DMORIS_USE_OPENBLAS=OFF -DBUILD_ALL=ON -DMORIS_HAVE_DEBUG=ON -DMORIS_HAVE_SYMBOLIC=ON -DMORIS_HAVE_SYMBOLIC_STRONG=ON -DMORIS_USE_EXAMPLES=ON ..

make -j 16

ctest


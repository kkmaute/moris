#!/bin/bash

#-----------------------------------------------------------------
# script to generate resource files for bash and tcsh shells
#-----------------------------------------------------------------

if [ ! $WORKSPACE ];then
    echo ""
    echo "Environment variable WORKSPACE needs to be defined"
    echo ""
    exit
fi

if [ -f $HOME/.cshrc_moris ];then
    echo ""
    echo "saving $HOME/.cshrc_moris to $HOME/.cshrc_moris.org"
    mv $HOME/.cshrc_moris  $HOME/.cshrc_moris.org
fi
if [ -f $HOME/.bashrc_moris ];then
    echo ""
    echo "saving $HOME/.bashrc_moris to $HOME/.bashrc_moris.org"
    mv $HOME/.bashrc_moris $HOME/.bashrc_moris.org
fi

echo ""
echo "creating $HOME/.cshrc_moris"

cd $WORKSPACE
export SPACK_ROOT=$WORKSPACE/spack
export PATH=$PATH/:$SPACK_ROOT/bin
. $SPACK_ROOT/share/spack/setup-env.sh
spack env activate .

export SPACKCOMP=`spack compiler list | tail -1`

export  CC=`spack compiler info $SPACKCOMP | grep 'cc ='  | awk -F = '{print $2}' | xargs ls`
export CXX=`spack compiler info $SPACKCOMP | grep 'cxx =' | awk -F = '{print $2}' | xargs ls`
export  FC=`spack compiler info $SPACKCOMP | grep 'fc ='  | awk -F = '{print $2}' | xargs ls`
export F77=`spack compiler info $SPACKCOMP | grep 'f77 =' | awk -F = '{print $2}' | xargs ls`

export GCCLIB=`spack compiler info $SPACKCOMP | grep 'cc ='  | awk -F = '{split($2,a,"/bin/");print a[1]}'`

export GCMMA_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "gcmma" )           {n=1}}END{print n}'`
export SNOPT_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "snopt" )           {n=1}}END{print n}'`
export LBFGS_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "lbfgs" )           {n=1}}END{print n}'`
export PETSC_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "petsc" )           {n=1}}END{print n}'`
export ABLIS_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "amdblis" )         {n=1}}END{print n}'`
export AFLAM_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "amdlibflame" )     {n=1}}END{print n}'`
export OAMKL_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "intel-oneapi-mkl" ){n=1}}END{print n}'`
export  IMKL_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "intel-mkl" )       {n=1}}END{print n}'`
export OBLAS_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "openblas" )        {n=1}}END{print n}'`
export SLEPC_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "slepc" )           {n=1}}END{print n}'`
export  DOXY_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "doxygen" )         {n=1}}END{print n}'`
export CLANG_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "llvm" )            {n=1}}END{print n}'`
export NINJA_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "ninja" )           {n=1}}END{print n}'`
export  ARBX_INSTALLED=`spack find | awk -F '@' 'BEGIN{n=0}{ if ( $1 == "arborx" )          {n=1}}END{print n}'`

export Trilinos_DIR=`spack location --install-dir trilinos`

if [ "$1" = "view" ];then
export MORISROOT=`spack location --install-dir moris`
echo 'setenv PATH $PATH/:'"$MORISROOT/bin/"                                    >> $HOME/.cshrc_moris
else
echo "setenv MORISROOT      $WORKSPACE/moris"                                  >> $HOME/.cshrc_moris
echo 'setenv MORISBUILDDBG  build_dbg'                                         >> $HOME/.cshrc_moris
echo 'setenv MORISBUILDOPT  build_opt'                                         >> $HOME/.cshrc_moris
echo 'setenv MORISOUTPUT    $MORISROOT/$MORISBUILDDBG/'                        >> $HOME/.cshrc_moris
echo ""                                                                        >> $HOME/.cshrc_moris
echo 'setenv MRD $MORISROOT/build_dbg/projects/mains/moris'                    >> $HOME/.cshrc_moris
echo 'setenv MRO $MORISROOT/build_opt/projects/mains/moris'                    >> $HOME/.cshrc_moris
echo 'setenv PATH $PATH/:$MORISROOT/share/scripts/'                            >> $HOME/.cshrc_moris
fi

echo ""                                                                        >> $HOME/.cshrc_moris
echo "setenv MPI_HOME"         `spack location --install-dir openmpi`          >> $HOME/.cshrc_moris
echo ""                                                                        >> $HOME/.cshrc_moris
echo "setenv Armadillo_DIR"    `spack location --install-dir armadillo`        >> $HOME/.cshrc_moris
echo "setenv Eigen3_DIR"       `spack location --install-dir eigen`            >> $HOME/.cshrc_moris
echo "setenv Boost_DIR"        `spack location --install-dir boost`            >> $HOME/.cshrc_moris
echo "setenv Boost_ROOT"       `spack location --install-dir boost`            >> $HOME/.cshrc_moris
echo "setenv ARPACK_DIR"       `spack location --install-dir arpack-ng`        >> $HOME/.cshrc_moris
echo "setenv SUPERLU_DIR"      `spack location --install-dir superlu`          >> $HOME/.cshrc_moris
echo "setenv SuperLU_DIST_DIR" `spack location --install-dir superlu-dist`     >> $HOME/.cshrc_moris
echo "setenv HDF5_DIR"         `spack location --install-dir hdf5`             >> $HOME/.cshrc_moris
echo "setenv NETCDF_DIR "      `spack location --install-dir netcdf-c`         >> $HOME/.cshrc_moris
echo "setenv ZLIB_DIR "        `spack location --install-dir zlib-ng`          >> $HOME/.cshrc_moris
echo "setenv SSL_DIR  "        `spack location --install-dir openssl`          >> $HOME/.cshrc_moris
echo "setenv CMAKE_DIR  "      `spack location --install-dir cmake`            >> $HOME/.cshrc_moris
echo "setenv Trilinos_DIR       $Trilinos_DIR"                                 >> $HOME/.cshrc_moris
echo ""                                                                        >> $HOME/.cshrc_moris

if [ $GCMMA_INSTALLED == "1" ];then
echo "setenv GCMMA_DIR"        `spack location --install-dir gcmma`            >> $HOME/.cshrc_moris
fi
if [ $SNOPT_INSTALLED == "1" ];then
echo "setenv SNOPT_DIR"        `spack location --install-dir snopt`            >> $HOME/.cshrc_moris
fi
if [ $LBFGS_INSTALLED == "1" ];then
echo "setenv LBFGSB_DIR"       `spack location --install-dir lbfgs`            >> $HOME/.cshrc_moris
fi
if [ $PETSC_INSTALLED == "1" ];then
export PETSC_DIR=`spack location --install-dir petsc`
echo "setenv PETSC_DIR"        $PETSC_DIR                                      >> $HOME/.cshrc_moris
fi
if [ $SLEPC_INSTALLED == "1" ];then
echo "setenv SLEPC_DIR"        `spack location --install-dir slepc`            >> $HOME/.cshrc_moris
fi
if [ $OAMKL_INSTALLED == "1" ];then
echo "setenv MKL_DIR" \
                 `spack location --install-dir intel-oneapi-mkl`"/mkl/latest"  >> $HOME/.cshrc_moris
fi
if [ $IMKL_INSTALLED == "1" ];then
echo "setenv MKL_DIR"          `spack location --install-dir intel-mkl`"/mkl"  >> $HOME/.cshrc_moris
fi
if [ $ABLIS_INSTALLED == "1" ];then
echo "setenv AMDBLIS_DIR"      `spack location --install-dir amdblis`          >> $HOME/.cshrc_moris
fi
if [ $AFLAM_INSTALLED == "1" ];then
echo "setenv AMDLIBFLAME_DIR"  `spack location --install-dir amdlibflame`      >> $HOME/.cshrc_moris
fi
if [ $OBLAS_INSTALLED == "1" ];then
echo "setenv OPENBLAS_DIR"     `spack location --install-dir openblas`         >> $HOME/.cshrc_moris
fi
if [ $DOXY_INSTALLED == "1" ];then
export DOXYGEN_DIR=`spack location --install-dir doxygen`
echo "setenv DOXYGEN_DIR"     $DOXYGEN_DIR                                     >> $HOME/.cshrc_moris
fi
if [ $CLANG_INSTALLED == "1" ];then
export CLANG_DIR=`spack location --install-dir llvm`
echo "setenv CLANG_DIR"       $CLANG_DIR                                       >> $HOME/.cshrc_moris
fi

if [ $NINJA_INSTALLED == "1" ];then
export NINJA_DIR=`spack location --install-dir ninja`
echo "setenv NINJA_DIR"       $NINJA_DIR                                       >> $HOME/.cshrc_moris
fi

if [ $ARBX_INSTALLED == "1" ];then
export ARBX_DIR=`spack location --install-dir arborx`
echo "setenv ARBX_DIR"       $ARBX_DIR                                         >> $HOME/.cshrc_moris
fi

echo ""                                                                        >> $HOME/.cshrc_moris
echo 'setenv PATH $MPI_HOME/bin/:$PATH'                                        >> $HOME/.cshrc_moris 
echo 'setenv PATH $NETCDF_DIR/bin/:$PATH'                                      >> $HOME/.cshrc_moris 
echo 'setenv PATH $Trilinos_DIR/bin/:$PATH'                                    >> $HOME/.cshrc_moris 
echo 'setenv PATH $CMAKE_DIR/bin/:$PATH'                                       >> $HOME/.cshrc_moris 

if [ $DOXY_INSTALLED == "1" ];then
echo 'setenv PATH $DOXYGEN_DIR/bin/:$PATH'                                     >> $HOME/.cshrc_moris 
fi

if [ $CLANG_INSTALLED == "1" ];then
echo 'setenv PATH $CLANG_DIR/bin/:$PATH'                                       >> $HOME/.cshrc_moris 
fi

if [ $NINJA_INSTALLED == "1" ];then
echo 'setenv PATH $NINJA_DIR/bin/:$PATH'                                       >> $HOME/.cshrc_moris 
fi

echo ""                                                                        >> $HOME/.cshrc_moris
echo "setenv LD_LIBRARY_PATH $GCCLIB/lib64"                                    >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$ZLIB_DIR/lib'                  >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$SSL_DIR/lib64'                 >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$MPI_HOME/lib64'                >> $HOME/.cshrc_moris 
echo ""                                                                        >> $HOME/.cshrc_moris
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$Armadillo_DIR/lib64'           >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$Boost_DIR/lib'                 >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$ARPACK_DIR/lib64'              >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$SUPERLU_DIR/lib'               >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$SuperLU_DIST_DIR/lib'          >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$Trilinos_DIR/lib'              >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$HDF5_DIR/lib'                  >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$NETCDF_DIR/lib64'              >> $HOME/.cshrc_moris 
echo ""                                                                        >> $HOME/.cshrc_moris

if [ $GCMMA_INSTALLED == "1" ];then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$GCMMA_DIR/lib'                 >> $HOME/.cshrc_moris 
fi
if [ $SNOPT_INSTALLED == "1" ];then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$SNOPT_DIR/lib'                 >> $HOME/.cshrc_moris 
fi
if [ $LBFGS_INSTALLED == "1" ];then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$LBFGSB_DIR/lib'                >> $HOME/.cshrc_moris 
fi
if [ $PETSC_INSTALLED == "1" ];then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$PETSC_DIR/lib'                 >> $HOME/.cshrc_moris 
fi
if [ $SLEPC_INSTALLED == "1" ];then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$SLEPC_DIR/lib'                 >> $HOME/.cshrc_moris 
fi
if [ $OBLAS_INSTALLED == "1" ];then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$OPENBLAS_DIR/lib'              >> $HOME/.cshrc_moris 
fi
if [ $ABLIS_INSTALLED == "1" ];then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$AMDBLIS_DIR/lib'               >> $HOME/.cshrc_moris
fi
if [ $AFLAM_INSTALLED == "1" ];then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$AMDLIBFLAME_DIR/lib'           >> $HOME/.cshrc_moris
fi
if [ $IMKL_INSTALLED == "1" ];then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$MKL_DIR/lib'                   >> $HOME/.cshrc_moris 
fi
if [ $OAMKL_INSTALLED == "1" ];then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$MKL_DIR/lib'                   >> $HOME/.cshrc_moris 
fi

echo ""                                                                        >> $HOME/.cshrc_moris
echo "setenv OMPI_MCA_rmaps_base_oversubscribe 1"                              >> $HOME/.cshrc_moris
echo "setenv OMP_NUM_THREADS 1"                                                >> $HOME/.cshrc_moris
echo "setenv OMPI_MCA_btl ^tcp"                                                >> $HOME/.cshrc_moris

echo ""                                                                        >> $HOME/.cshrc_moris
echo "setenv CC  $CC"                                                          >> $HOME/.cshrc_moris
echo "setenv CXX $CXX"                                                         >> $HOME/.cshrc_moris
echo "setenv FC  $FC"                                                          >> $HOME/.cshrc_moris
echo "setenv F77 $F77"                                                         >> $HOME/.cshrc_moris
echo ""                                                                        >> $HOME/.cshrc_moris

if [ $PETSC_INSTALLED == "1" ];then
export      GFORTLIB=`ldd $PETSC_DIR/lib/libpetsc.so | grep gfortran | awk '{print $1}'`
export GFORTLIB_PATH=`ldd $PETSC_DIR/lib/libpetsc.so | grep gfortran | awk '{print $3}' | xargs dirname`
else
export      GFORTLIB=`ldd $Trilinos_DIR/lib/libexodus.so | grep gfortran | awk '{print $1}'`
export GFORTLIB_PATH=`ldd $Trilinos_DIR/lib/libexodus.so | grep gfortran | awk '{print $3}' | xargs dirname`
fi

echo "setenv GFORTLIB $GFORTLIB"                                               >> $HOME/.cshrc_moris
echo "setenv GFORTLIB_PATH $GFORTLIB_PATH"                                     >> $HOME/.cshrc_moris

echo ""
echo "creating $HOME/.bashrc_moris"

sed -rn 's/^\s*setenv\s+(\S+)\s+/export \1=/p' $HOME/.cshrc_moris > $HOME/.bashrc_moris
echo ""

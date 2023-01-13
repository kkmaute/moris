#!/bin/tcsh

rm -f $HOME/.bashrc_moris

cd $WORKSPACE
setenv SPACK_ROOT $WORKSPACE/spack
setenv PATH $PATH/:$SPACK_ROOT/bin
source $SPACK_ROOT/share/spack/setup-env.csh
spack env activate .

setenv SPACKCOMP `spack compiler list | tail -1`
setenv CC        `spack compiler info $SPACKCOMP | grep 'cc ='  | awk -F = '{print $2}' | xargs ls`
setenv CXX       `spack compiler info $SPACKCOMP | grep 'cxx =' | awk -F = '{print $2}' | xargs ls`
setenv FC        `spack compiler info $SPACKCOMP | grep 'fc ='  | awk -F = '{print $2}' | xargs ls`
setenv F77       `spack compiler info $SPACKCOMP | grep 'f77 =' | awk -F = '{print $2}' | xargs ls`

setenv GCCLIB    `spack compiler info $SPACKCOMP | grep 'cc ='  | awk -F = '{split($2,a,"/bin/");print a[1]}'`
setenv PETSC_DIR `spack location --install-dir petsc`

echo "export MORISROOT=$WORKSPACE/moris"                                       >> $HOME/.bashrc_moris
echo 'export MORISBUILDDBG=build_dbg'                                          >> $HOME/.bashrc_moris
echo 'export MORISBUILDOPT=build_opt'                                          >> $HOME/.bashrc_moris
echo 'export MORISOUTPUT=$MORISROOT/$MORISBUILDDBG/'                           >> $HOME/.bashrc_moris
echo ""                                                                        >> $HOME/.bashrc_moris
echo 'setenv MRD $MORISROOT/build_dbg/projects/mains/moris'                    >> $HOME/.bashrc_moris
echo 'setenv MRO $MORISROOT/build_opt/projects/mains/moris'                    >> $HOME/.bashrc_moris
echo ""                                                                        >> $HOME/.bashrc_moris
echo "export MPI_HOME=`spack location --install-dir openmpi`"                  >> $HOME/.bashrc_moris
echo ""                                                                        >> $HOME/.bashrc_moris
echo "export PETSC_DIR=$PETSC_DIR"                                             >> $HOME/.bashrc_moris
echo ""                                                                        >> $HOME/.bashrc_moris
echo "export Armadillo_DIR=`spack location --install-dir armadillo`"           >> $HOME/.bashrc_moris
echo "export Eigen3_DIR=`spack location --install-dir eigen`"                  >> $HOME/.bashrc_moris
echo "export BOOST_DIR=`spack location --install-dir boost`"                   >> $HOME/.bashrc_moris
echo "export BOOST_ROOT=`spack location --install-dir boost`"                  >> $HOME/.bashrc_moris
echo "export GCMMA_DIR=`spack location --install-dir gcmma`"                   >> $HOME/.bashrc_moris
echo "export SNOPT_DIR=`spack location --install-dir snopt`"                   >> $HOME/.bashrc_moris
echo "export LBFGSB_DIR=`spack location --install-dir lbfgs`"                  >> $HOME/.bashrc_moris
echo "export ARPACK_DIR=`spack location --install-dir arpack-ng`"              >> $HOME/.bashrc_moris
echo "export SUPERLU_DIR=`spack location --install-dir superlu`"               >> $HOME/.bashrc_moris
echo "export SuperLU_DIST_DIR=`spack location --install-dir superlu-dist`"     >> $HOME/.bashrc_moris
echo "export Trilinos_DIR=`spack location --install-dir trilinos`"             >> $HOME/.bashrc_moris
echo "export HDF5_DIR=`spack location --install-dir hdf5`"                     >> $HOME/.bashrc_moris
echo "export MKL_DIR=`spack location --install-dir intel-mkl`/mkl"             >> $HOME/.bashrc_moris
echo "export NETCDF_DIR=`spack location --install-dir netcdf-c`"               >> $HOME/.bashrc_moris
echo "export ZLIB_DIR=`spack location --install-dir zlib`"                     >> $HOME/.bashrc_moris
echo "export SSL_DIR=`spack location --install-dir openssl`"                   >> $HOME/.bashrc_moris
echo "export CMAKE_DIR=`spack location --install-dir cmake`"                   >> $HOME/.bashrc_moris
echo ""                                                                        >> $HOME/.bashrc_moris
echo 'export ZLIB_LIBRARY_DIR=$ZLIB_DIR/lib'                                   >> $HOME/.bashrc_moris 
echo 'export SSL_LIBRARY_DIR=$SSL_DIR/lib'                                     >> $HOME/.bashrc_moris 
echo ""                                                                        >> $HOME/.bashrc_moris
echo 'export PATH=$PATH/:$MORISROOT/share/scripts/'                            >> $HOME/.bashrc_moris 
echo 'export PATH=$MPI_HOME/bin/:$PATH'                                        >> $HOME/.bashrc_moris 
echo 'export PATH=$NETCDF_DIR/bin/:$PATH'                                      >> $HOME/.bashrc_moris 
echo 'export PATH=$Trilinos_DIR/bin/:$PATH'                                    >> $HOME/.bashrc_moris 
echo 'export PATH=$CMAKE_DIR/bin/:$PATH'                                       >> $HOME/.bashrc_moris 
echo ""                                                                        >> $HOME/.bashrc_moris
echo 'export LD_LIBRARY_PATH=$GCCLIB/lib64'                                    >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$MPI_HOME/lib64'                >> $HOME/.bashrc_moris 
echo ""                                                                        >> $HOME/.bashrc_moris
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$OPENBLAS_DIR/lib'              >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$Armadillo_DIR/lib64'           >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$BOOST_DIR/lib'                 >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$GCMMA_DIR/lib'                 >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$SNOPT_DIR/lib'                 >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$LBFGSB_DIR/lib'                >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$ARPACK_DIR/lib64'              >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$SUPERLU_DIR/lib'               >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$SuperLU_DIST_DIR/lib'          >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$Trilinos_DIR/lib'              >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$PETSC_DIR/lib'                 >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$HDF5_DIR/lib'                  >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$MKL_DIR/lib'                   >> $HOME/.bashrc_moris 
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:$NETCDF_DIR/lib64'              >> $HOME/.bashrc_moris
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH/:/usr/lib64/psm2-compat'         >> $HOME/.bashrc_moris
echo ""                                                                        >> $HOME/.bashrc_moris
echo "export CC=$CC"                                                           >> $HOME/.bashrc_moris
echo "export CXX=$CXX"                                                         >> $HOME/.bashrc_moris
echo "export FC=$FC"                                                           >> $HOME/.bashrc_moris
echo "export F77=$F77"                                                         >> $HOME/.bashrc_moris

setenv GFORTLIB      `ldd $PETSC_DIR/lib/libpetsc.so | grep libgfortran | awk '{print $3}' | awk -F "/" '{print $NF}'`
setenv GFORTLIB_PATH `ldd $PETSC_DIR/lib/libpetsc.so | grep libgfortran | awk '{print $3}' | awk -F libgfortran '{print $1}'`

echo ""
echo "export GFORTLIB=$GFORTLIB"                                               >> $HOME/.bashrc_moris
echo "export GFORTLIB_PATH=$GFORTLIB_PATH"                                     >> $HOME/.bashrc_moris
echo ""
echo "export OMPI_MCA_rmaps_base_oversubscribe=1"                              >> $HOME/.bashrc_moris
echo "export OMP_NUM_THREADS=1"                                                >> $HOME/.bashrc_moris
echo "export OMPI_MCA_mtl=psm2"                                                >> $HOME/.bashrc_moris


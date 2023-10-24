#!/bin/tcsh

rm -f $HOME/.cshrc_moris

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

setenv GCMMA_INSTALLED `spack find gcmma     | awk 'BEGIN{n=0}{ n=n+1; if ($2 == "No" && n == 1) {print 0} else { if ( n == 1 ) {print 1}}}'`
setenv SNOPT_INSTALLED `spack find snopt     | awk 'BEGIN{n=0}{ n=n+1; if ($2 == "No" && n == 1) {print 0} else { if ( n == 1 ) {print 1}}}'`
setenv LBFGS_INSTALLED `spack find lbfgs     | awk 'BEGIN{n=0}{ n=n+1; if ($2 == "No" && n == 1) {print 0} else { if ( n == 1 ) {print 1}}}'`
setenv PETSC_INSTALLED `spack find petsc     | awk 'BEGIN{n=0}{ n=n+1; if ($2 == "No" && n == 1) {print 0} else { if ( n == 1 ) {print 1}}}'`
setenv   MKL_INSTALLED `spack find intel-mkl | awk 'BEGIN{n=0}{ n=n+1; if ($2 == "No" && n == 1) {print 0} else { if ( n == 1 ) {print 1}}}'`
setenv OBLAS_INSTALLED `spack find openblas  | awk 'BEGIN{n=0}{ n=n+1; if ($2 == "No" && n == 1) {print 0} else { if ( n == 1 ) {print 1}}}'`
setenv SLEPC_INSTALLED `spack find slepc     | awk 'BEGIN{n=0}{ n=n+1; if ($2 == "No" && n == 1) {print 0} else { if ( n == 1 ) {print 1}}}'`
setenv  DOXY_INSTALLED `spack find doxygen   | awk 'BEGIN{n=0}{ n=n+1; if ($2 == "No" && n == 1) {print 0} else { if ( n == 1 ) {print 1}}}'`
setenv CLANG_INSTALLED `spack find llvm      | awk 'BEGIN{n=0}{ n=n+1; if ($2 == "No" && n == 1) {print 0} else { if ( n == 1 ) {print 1}}}'`

setenv Trilinos_DIR `spack location --install-dir trilinos`

echo "setenv MORISROOT      $WORKSPACE/moris"                                  >> $HOME/.cshrc_moris
echo 'setenv MORISBUILDDBG  build_dbg'                                         >> $HOME/.cshrc_moris
echo 'setenv MORISBUILDOPT  build_opt'                                         >> $HOME/.cshrc_moris
echo 'setenv MORISOUTPUT    $MORISROOT/$MORISBUILDDBG/'                        >> $HOME/.cshrc_moris
echo ""                                                                        >> $HOME/.cshrc_moris
echo 'setenv MRD $MORISROOT/build_dbg/projects/mains/moris'                    >> $HOME/.cshrc_moris
echo 'setenv MRO $MORISROOT/build_opt/projects/mains/moris'                    >> $HOME/.cshrc_moris
echo ""                                                                        >> $HOME/.cshrc_moris
echo "setenv MPI_HOME"         `spack location --install-dir openmpi`          >> $HOME/.cshrc_moris
echo ""                                                                        >> $HOME/.cshrc_moris
echo "setenv Armadillo_DIR"    `spack location --install-dir armadillo`        >> $HOME/.cshrc_moris
echo "setenv Eigen3_DIR"       `spack location --install-dir eigen`            >> $HOME/.cshrc_moris
echo "setenv BOOST_DIR"        `spack location --install-dir boost`            >> $HOME/.cshrc_moris
echo "setenv BOOST_ROOT"       `spack location --install-dir boost`            >> $HOME/.cshrc_moris
echo "setenv ARPACK_DIR"       `spack location --install-dir arpack-ng`        >> $HOME/.cshrc_moris
echo "setenv SUPERLU_DIR"      `spack location --install-dir superlu`          >> $HOME/.cshrc_moris
echo "setenv SuperLU_DIST_DIR" `spack location --install-dir superlu-dist`     >> $HOME/.cshrc_moris
echo "setenv HDF5_DIR"         `spack location --install-dir hdf5`             >> $HOME/.cshrc_moris
echo "setenv NETCDF_DIR "      `spack location --install-dir netcdf-c`         >> $HOME/.cshrc_moris
echo "setenv ZLIB_DIR "        `spack location --install-dir zlib`             >> $HOME/.cshrc_moris
echo "setenv SSL_DIR  "        `spack location --install-dir openssl`          >> $HOME/.cshrc_moris
echo "setenv CMAKE_DIR  "      `spack location --install-dir cmake`            >> $HOME/.cshrc_moris
echo "setenv Trilinos_DIR       $Trilinos_DIR"                                 >> $HOME/.cshrc_moris
echo ""                                                                        >> $HOME/.cshrc_moris
echo 'setenv ZLIB_LIBRARY_DIR  $ZLIB_DIR/lib'                                  >> $HOME/.cshrc_moris 
echo 'setenv SSL_LIBRARY_DIR   $SSL_DIR/lib'                                   >> $HOME/.cshrc_moris 
echo ""                                                                        >> $HOME/.cshrc_moris

if ( $GCMMA_INSTALLED == "1" ) then
echo "setenv GCMMA_DIR"        `spack location --install-dir gcmma`            >> $HOME/.cshrc_moris
endif
if ( $SNOPT_INSTALLED == "1" ) then
echo "setenv SNOPT_DIR"        `spack location --install-dir snopt`            >> $HOME/.cshrc_moris
endif
if ( $LBFGS_INSTALLED == "1" ) then
echo "setenv LBFGSB_DIR"       `spack location --install-dir lbfgs`            >> $HOME/.cshrc_moris
endif
if ( $PETSC_INSTALLED == "1" ) then
setenv PETSC_DIR    `spack location --install-dir petsc`
echo "setenv PETSC_DIR"        $PETSC_DIR                                      >> $HOME/.cshrc_moris
endif
if ( $SLEPC_INSTALLED == "1" ) then
setenv SLEPC_DIR    `spack location --install-dir slepc`
echo "setenv SLEPC_DIR"        $SLEPC_DIR                                      >> $HOME/.cshrc_moris
endif
if ( $MKL_INSTALLED == "1" ) then
echo "setenv MKL_DIR"          `spack location --install-dir intel-mkl`"/mkl"  >> $HOME/.cshrc_moris
endif
if ( $OBLAS_INSTALLED == "1" ) then
setenv OPENBLAS_DIR `spack location --install-dir openblas`"/lib"
echo "setenv OPENBLAS_DIR"     $OPENBLAS_DIR                                   >> $HOME/.cshrc_moris
endif
if ( $DOXY_INSTALLED == "1" ) then
setenv DOXYGEN_DIR `spack location --install-dir doxygen`
echo "setenv DOXYGEN_DIR"     $DOXYGEN_DIR                                     >> $HOME/.cshrc_moris
endif
if ( $CLANG_INSTALLED == "1" ) then
setenv CLANG_DIR `spack location --install-dir llvm`
echo "setenv CLANG_DIR"       $CLANG_DIR                                       >> $HOME/.cshrc_moris
endif

echo ""                                                                        >> $HOME/.cshrc_moris
echo 'setenv PATH $MPI_HOME/bin/:$PATH'                                        >> $HOME/.cshrc_moris 
echo 'setenv PATH $NETCDF_DIR/bin/:$PATH'                                      >> $HOME/.cshrc_moris 
echo 'setenv PATH $Trilinos_DIR/bin/:$PATH'                                    >> $HOME/.cshrc_moris 
echo 'setenv PATH $CMAKE_DIR/bin/:$PATH'                                       >> $HOME/.cshrc_moris 
if ( $DOXY_INSTALLED == "1" ) then
echo 'setenv PATH $DOXYGEN_DIR/bin/:$PATH'                                     >> $HOME/.cshrc_moris 
endif
if ( $CLANG_INSTALLED == "1" ) then
echo 'setenv PATH $CLANG_DIR/bin/:$PATH'                                       >> $HOME/.cshrc_moris 
endif

echo ""                                                                        >> $HOME/.cshrc_moris
echo "setenv LD_LIBRARY_PATH $GCCLIB/lib64"                                    >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$MPI_HOME/lib64'                >> $HOME/.cshrc_moris 
echo ""                                                                        >> $HOME/.cshrc_moris
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$Armadillo_DIR/lib64'           >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$BOOST_DIR/lib'                 >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$ARPACK_DIR/lib64'              >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$SUPERLU_DIR/lib'               >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$SuperLU_DIST_DIR/lib'          >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$Trilinos_DIR/lib'              >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$HDF5_DIR/lib'                  >> $HOME/.cshrc_moris 
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$NETCDF_DIR/lib64'              >> $HOME/.cshrc_moris 
echo ""                                                                        >> $HOME/.cshrc_moris

if ( $GCMMA_INSTALLED == "1" ) then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$GCMMA_DIR/lib'                 >> $HOME/.cshrc_moris 
endif
if ( $SNOPT_INSTALLED == "1" ) then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$SNOPT_DIR/lib'                 >> $HOME/.cshrc_moris 
endif
if ( $LBFGS_INSTALLED == "1" ) then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$LBFGSB_DIR/lib'                >> $HOME/.cshrc_moris 
endif
if ( $PETSC_INSTALLED == "1" ) then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$PETSC_DIR/lib'                 >> $HOME/.cshrc_moris 
endif
if ( $SLEPC_INSTALLED == "1" ) then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$SLEPC_DIR/lib'                 >> $HOME/.cshrc_moris 
endif
if ( $MKL_INSTALLED == "1" ) then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$MKL_DIR/lib'                   >> $HOME/.cshrc_moris 
endif
if ( $OBLAS_INSTALLED == "1" ) then
echo 'setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$OPENBLAS_DIR/'                 >> $HOME/.cshrc_moris 
endif

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

if ( $PETSC_INSTALLED == "1" ) then
setenv GFORTLIB      `ldd $PETSC_DIR/lib/libpetsc.so | grep gfortran | awk '{print $1}'`
setenv GFORTLIB_PATH `ldd $PETSC_DIR/lib/libpetsc.so | grep gfortran | awk '{print $3}' | xargs dirname`
else
setenv GFORTLIB      `ldd $Trilinos_DIR/lib/libexodus.so | grep gfortran | awk '{print $1}'`
setenv GFORTLIB_PATH `ldd $Trilinos_DIR/lib/libexodus.so | grep gfortran | awk '{print $3}' | xargs dirname`
endif

echo "setenv GFORTLIB $GFORTLIB"                                               >> $HOME/.cshrc_moris
echo "setenv GFORTLIB_PATH $GFORTLIB_PATH"                                     >> $HOME/.cshrc_moris

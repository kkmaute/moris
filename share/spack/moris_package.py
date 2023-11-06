# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install moris
#
# You can edit this file again by typing:
#
#     spack edit moris
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *

import os
import re
import subprocess

class Moris(CMakePackage):
    """MORIS"""

    git      = "ssh://git@github.com/kkmaute/moris"

    maintainers = ['kkmaute']

    version('main', branch='main', submodules=True, preferred=True)

    variant("default",  default=True,  description='Compile with default setting')
    variant("petsc",    default=True,  description="Compile with support for petsc")
    variant("slepc",    default=True,  description="Compile with support for slepc")
    variant("pardiso",  default=True,  description="Compile with support for pardiso solver")
    variant("gcmma",    default=True,  description="Compile with support for gcmma algorithm")
    variant("snopt",    default=True,  description="Compile with support for snopt algorithm")
    variant("lbfgs",    default=True,  description="Compile with support for lbfgs algorithm")
    variant("mumps",    default=False, description="Compile with support for mumps solver")
    variant("openblas", default=False, description="Compile with support for openblas")
    variant("mkl",      default=False, description="Compile with support for intel-oneapi-mkl")
    variant("lapack",   default=False, description="Compile with support generic blas and lapack")
    variant("tests",    default=False, description="Compile with unit tests")
    variant("examples", default=False, description="Compile with examples")
    
    depends_on('armadillo@9.800.3')
    depends_on('arpack-ng@3.8.0')
    depends_on('boost+filesystem+log+serialization+system+thread+timer')
    depends_on('eigen')
    depends_on('hdf5')

    depends_on('gcmma',             when="+gcmma")
    depends_on('snopt',             when="+snopt")
    depends_on('netlib-lapack',     when="+lapack")
    depends_on('openblas',          when="+openblas")
    depends_on('intel-oneapi-mkl',  when="+mkl")
    depends_on('intel-oneapi-mkl',  when="+pardiso")
    depends_on('lbfgs',             when="+lbfgs")
 
    depends_on('mpi')
 
    depends_on('superlu')
    depends_on('superlu-dist')
   
    depends_on('trilinos@13.4')
    depends_on('trilinos+boost+hdf5+mpi+suite-sparse+superlu-dist+amesos+anasazi+aztec+belos+chaco+epetra+exodus+ifpack+ifpack2+ml+rol+stk+zoltan2')
    depends_on('trilinos+pardiso', when="+pardiso")
    depends_on('trilinos+mumps',   when="+mumps")

    depends_on('petsc@3.17.4',                       when="+petsc")
    depends_on('petsc+mpi+metis+hypre+suite-sparse', when="+petsc")
    #depends_on('petsc+mkl-pardiso',                  when="+petsc +pardiso")
    depends_on('petsc+mumps',                        when="+petsc +mumps")

    depends_on('slepc',                              when="+slepc")

    conflicts('openblas',   when='+pardiso')
    conflicts('openblas',   when='+mkl')
    conflicts('openblas',   when='+lapack')
    
    conflicts('mkl',        when='+openblas')
    conflicts('mkl',        when='+lapack')
    
    conflicts('+lapack',     when='+pardiso')
    conflicts('+lapack',     when='+mkl')
    conflicts('+lapack',     when='+openblas')

    def cmake_args(self):
        spec = self.spec
        options = []

        options.extend(['-DBUILD_ALL=ON'])
        options.extend(['-DMORIS_HAVE_DEBUG=ON'])
        options.extend(['-DMORIS_HAVE_SYMBOLIC=ON'])
        options.extend(['-DMORIS_HAVE_SYMBOLIC_STRONG=ON'])
        options.extend(['-DMORIS_HAVE_TRILINOS_NEW_CMAKE=OFF'])

        if '+mkl' in spec:
            options.extend(['-DMORIS_USE_MKL=ON'])
            options.extend(['-DMORIS_USE_OPENBLAS=OFF'])
            options.extend(['-DMORIS_USE_LAPACK=OFF'])
            
        if '+openblas' in spec:
            options.extend(['-DMORIS_USE_MKL=OFF'])
            options.extend(['-DMORIS_USE_OPENBLAS=ON'])
            options.extend(['-DMORIS_USE_LAPACK=OFF'])

        if '+lapack' in spec:
            options.extend(['-DMORIS_USE_MKL=OFF'])
            options.extend(['-DMORIS_USE_OPENBLAS=OFF'])
            options.extend(['-DMORIS_USE_LAPACK=ON'])

        if '-petsc' in spec:
            options.extend([ '-DMORIS_HAVE_PETSC=OFF' ])

        if '-slepc' in spec:
            options.extend([ '-DMORIS_HAVE_SLEPC=OFF' ])

        if '-gcmma' in spec:
            options.extend([ '-DMORIS_HAVE_GCMMA=OFF' ])

        if '-lbfgs' in spec:
            options.extend([ '-DMORIS_HAVE_LBFGS=OFF' ])

        if '-snopt' in spec:
            options.extend([ '-DMORIS_HAVE_SNOPT=OFF' ])

        if '-mumps' in spec:
            options.extend([ '-DMORIS_USE_MUMPS=OFF' ])

        if '-pardiso' in spec:
            options.extend([ '-DMORIS_USE_PARDISO=OFF' ])

        if '+pardiso' in spec:
            options.extend([ '-DMORIS_USE_MKL=ON' ])
            options.extend([ '-DMORIS_USE_OPENBLAS=OFF' ])

        if '+tests' in spec:
            options.extend([ '-DMORIS_USE_TESTS=ON' ])

        if '-tests' in spec:
            options.extend([ '-DMORIS_USE_TESTS=OFF' ])

        if '+examples' in spec:
            options.extend([ '-DMORIS_USE_EXAMPLES=ON' ])

        if '-examples' in spec:
            options.extend([ '-DMORIS_USE_EXAMPLES=OFF' ])

        return options
    
    def setup_build_environment(self, env):
        if '-petsc' in self.spec:
            # Get gfortran library of space fortran compiler
            fc = Executable(self.compiler.fc)
            gfortlib_path = fc("--print-file-name", "libgfortran.so", output=str
 	                ).rsplit("/",1)[0]
            gfortlib_name = fc("--print-file-name", "libgfortran.so", output=str
 	                ).rsplit("/",1)[1]
            gfortlib_name = gfortlib_name.replace('\n','')
            env.append_path('LD_LIBRARY_PATH', gfortlib_path)
            env.set('GFORTLIB_PATH', gfortlib_path)
            env.set('GFORTLIB', gfortlib_name)
        else:
            # Get gfortran library that petsc links with
            petsc_libs = self.spec['petsc'].libs
            ldd_command = ['ldd', petsc_libs[0]]
            proc_result = subprocess.run(ldd_command, stdout=subprocess.PIPE, check=True)
            gfortlib_path_regex = r'(/[\S]*libgfortran.so(\.\d+)*)'
            gfortlib_result = re.search(gfortlib_path_regex, proc_result.stdout.decode('utf-8'))
            if gfortlib_result:
                gfortlib_path, gfortlib_name = os.path.split(gfortlib_result.group(0)) 
                env.append_path('LD_LIBRARY_PATH', gfortlib_path)
                env.set('GFORTLIB_PATH', gfortlib_path)
                env.set('GFORTLIB', gfortlib_name)
            else:
                print('WARNING: Could not determine which gfortran library to use. moris may not find the fortran runtime library correctly.')

        if '+openblas' in self.spec:
            env.set('OPENBLAS_DIR', self.spec['blas'].prefix)

        if '+lapack' in self.spec:
            env.set('LAPACK_DIR', self.spec['blas'].prefix)

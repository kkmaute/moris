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

    git = "ssh://git@github.com/kkmaute/moris"

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
    variant("tests",    default=False, description="Compile with unit tests")
    variant("examples", default=False, description="Compile with examples")
    variant("debug",    default=False, description="Compile with debug support")

    depends_on('arborx@1.6')
    depends_on('armadillo')
    depends_on('arpack-ng')
    depends_on('boost+filesystem+log+serialization+system+thread+timer')
    depends_on('eigen')
    depends_on('hdf5')

    depends_on('gcmma', when="+gcmma")
    depends_on('snopt', when="+snopt")
    depends_on('lbfgs', when="+lbfgs")
 
    depends_on('mpi')
 
    depends_on('superlu')
    depends_on('superlu-dist@8.1')
    depends_on('suite-sparse@5.13.0')
   
    depends_on('trilinos@15.1.1')
    depends_on('trilinos+boost+hdf5+mpi+suite-sparse+superlu-dist+amesos+anasazi+aztec+belos+chaco+epetra+exodus+ifpack+ifpack2+ml+rol+stk+zoltan2')
    depends_on('trilinos+pardiso', when="+pardiso")
    depends_on('trilinos~pardiso', when="-pardiso")
    depends_on('trilinos+mumps',   when="+mumps")
    depends_on('trilinos~mumps',   when="-mumps")

    depends_on('petsc',                              when="+petsc")
    depends_on('petsc+mpi+metis+hypre+suite-sparse', when="+petsc")
    #depends_on('petsc+mkl-pardiso',                 when="+petsc +pardiso")
    depends_on('petsc+mumps',                        when="+petsc +mumps")
    depends_on('petsc~mumps',                        when="+petsc -mumps")

    depends_on('slepc', when="+slepc")

    def cmake_args(self):
        spec = self.spec
        options = []

        options.extend(['-DBUILD_ALL=ON'])
        options.extend(['-DMORIS_HAVE_DEBUG=OFF'])
        options.extend(['-DMORIS_HAVE_SYMBOLIC=OFF'])
        options.extend(['-DMORIS_HAVE_SYMBOLIC_STRONG=OFF'])

        if '+debug' in spec:
            options.extend(['-DMORIS_HAVE_DEBUG=ON'])
            options.extend(['-DMORIS_HAVE_SYMBOLIC=ON'])
            options.extend(['-DMORIS_HAVE_SYMBOLIC_STRONG=ON'])

        if self.spec["blas"].name in ["intel-mkl", "intel-oneapi-mkL"]:
            options.extend(['-DMORIS_USE_MKL=ON'])
            options.extend(['-DMORIS_USE_OPENBLAS=OFF'])
            options.extend(['-DMORIS_USE_LAPACK=OFF'])
            options.extend(['-DMORIS_USE_ACML=OFF'])
           
        if self.spec["blas"].name in ["openblas"]: 
            options.extend(['-DMORIS_USE_MKL=OFF'])
            options.extend(['-DMORIS_USE_OPENBLAS=ON'])
            options.extend(['-DMORIS_USE_LAPACK=OFF'])
            options.extend(['-DMORIS_USE_ACML=OFF'])

        if self.spec["blas"].name in ["amdblis"]:
            options.extend(['-DMORIS_USE_MKL=OFF'])
            options.extend(['-DMORIS_USE_OPENBLAS=OFF'])
            options.extend(['-DMORIS_USE_LAPACK=ON'])
            options.extend(['-DMORIS_USE_ACML=ON'])

        if self.spec["blas"].name in ["netlib-lapack"]:
            options.extend(['-DMORIS_USE_MKL=OFF'])
            options.extend(['-DMORIS_USE_OPENBLAS=OFF'])
            options.extend(['-DMORIS_USE_LAPACK=ON'])
            options.extend(['-DMORIS_USE_ACML=OFF'])
        
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
        openssl_lib_path = self.spec['openssl'].prefix.lib
        if not os.path.exists(openssl_lib_path):
            openssl_lib_path = self.spec['openssl'].prefix.lib64
        env.append_path('LD_LIBRARY_PATH', openssl_lib_path)

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

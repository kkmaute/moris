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


class Moris(CMakePackage):
    """MORIS"""

    git = "ssh://git@github.com/kkmaute/moris"
    
    maintainers = ['kkmaute']

    version('main', branch='main', submodules=True, preferred=True)

    variant("default",  default=True, description="Compile with default setting")
    variant("petsc",    default=True, description="Compile with support for petsc")
    variant("pardiso",  default=True, description="Compile with support for pardiso solver")
    variant("mumps",    default=True, description="Compile with support for mumps solver")
    variant("gcmma",    default=True, description="Compile with support for gcmma algorithm")
    variant("snopt",    default=True, description="Compile with support for snopt algorithm")
    variant("lbfgs",    default=True, description="Compile with support for lbfgs algorithm")
    variant("openblas", default=False,description="Compile with support for openblas")
    
    depends_on('armadillo@9.800.3')

    depends_on('arpack-ng@3.8.0')

    depends_on('boost+filesystem+log+serialization+system+thread+timer')

    depends_on('eigen')

    depends_on('gcmma', when="+gcmma")
    
    depends_on('hdf5')

    depends_on('intel-mkl', when="+pardiso")

    depends_on('lbfgs', when="+lbfgs")
 
    depends_on('mpi')
 
    depends_on('superlu')
    
    depends_on('snopt', when="+snopt")

    depends_on('trilinos@13.4')
    depends_on('trilinos+boost+hdf5+mpi+suite-sparse+superlu-dist+amesos+anasazi+aztec+belos+chaco+epetra+exodus+ifpack+ifpack2+ml+rol+stk+zoltan2')
    depends_on('trilinos+pardiso', when="+pardiso")
    depends_on('trilinos+mumps',   when="+mumps")
    
    depends_on('petsc@3.17.4',                       when="+petsc")
    depends_on('petsc+mpi+metis+hypre+suite-sparse', when="+petsc")
    depends_on('petsc+mkl-pardiso',                  when="+petsc +pardiso")
    depends_on('petsc+mumps',                        when="+petsc +mumps")
 
    conflicts('+openblas', when='+pardiso')
 
    def cmake_args(self):
    
        spec = self.spec
        options = []

        options.extend([ '-DBUILD_ALL=ON -DMORIS_HAVE_DEBUG=ON -DMORIS_HAVE_SYMBOLIC=ON -DMORIS_HAVE_SYMBOLIC_STRONG=ON -DMORIS_USE_EXAMPLES=ON' ])
        
        if '-petsc' in spec:
            options.extend([ '-DMORIS_HAVE_PETSC=OFF' ])
        
        if '-gcmma' in spec:
            options.extend([ '-DMORIS_HAVE_GCMMA=OFF' ])

        if '-lbfgs' in spec:
            options.extend([ '-DMORIS_HAVE_LBFGS=OFF' ])

        if '-snopt' in spec:
            options.extend([ '-DMORIS_HAVE_SNOPT=OFF' ])

        if '-mumps' in spec:
            options.extend([ '-DMORIS_USE_MUMPS=OFF' ])

        if '+pardiso' in spec:
            options.extend([ '-DMORIS_USE_OPENBLAS=OFF' ])
            options.extend([ '-DMORIS_USE_MKL=ON' ])

        if '-pardiso' in spec:
            options.extend([ '-DMORIS_USE_PARDISO=OFF' ])

       if '-openblas' in spec:
            options.extend([ '-DMORIS_USE_OPENBLAS=OFF' ])
            options.extend([ '-DMORIS_USE_MKL=ON' ])

        return options

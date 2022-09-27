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

    #git      = "ssh://git@github.com/kkmaute/moris"
    git      = "ssh://pegasus:/home/maute/codes/moris"

    maintainers = ['kmaute']

    version('main', branch='main', submodules=True, preferred=True)

    variant( 'default',     default=True,    description='Compile with default setting'          )

    depends_on('armadillo@9.800.3')
    depends_on('armadillo+hdf5')

    depends_on('arpack-ng@3.8.0')

    depends_on('boost+filesystem+log+serialization+system+thread+timer')

    depends_on('eigen')

    depends_on('gcmma')

    depends_on('hdf5')

    depends_on('lbfgs')
 
    depends_on('mpi')
 
    depends_on('superlu')
    
    depends_on('snopt')

    depends_on('trilinos@13.4')
    depends_on('trilinos+boost+hdf5+mpi+mumps+suite-sparse+superlu-dist+pardiso+amesos+anasazi+aztec+belos+chaco+epetra+exodus+ifpack+ifpack2+ml+rol+stk+zoltan2')

    depends_on('petsc@3.17.4')
    depends_on('petsc+mpi+metis+hypre+mumps+mkl-pardiso+suite-sparse')
 
    def cmake_args(self):
    
        spec = self.spec
        options = []

        options.extend([ '-DBUILD_ALL=ON -DMORIS_HAVE_DEBUG=ON -DMORIS_HAVE_SYMBOLIC=ON -DMORIS_HAVE_SYMBOLIC_STRONG=ON -DMORIS_USE_EXAMPLES=ON' ])

        return options

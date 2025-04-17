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
#     spack install lbfgs
#
# You can edit this file again by typing:
#
#     spack edit lbfgs
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack.package import *


class Lbfgs(Package):
    """FIXME: Put a proper description of your package here."""

    # FIXME: Add a proper url for your package's homepage here.
    git = "ssh://git@github.com/kkmaute/lbfgs"

    version('main', branch='main', submodules=True, preferred=True)

    depends_on("gmake", type="build")
    depends_on("fortran", type="build")
    
    # FIXME: Add a list of GitHub accounts to
    # notify when the package is updated.
    # maintainers = ["github_user1", "github_user2"]

    # FIXME: Add proper versions here.
    # version("1.2.4")

    # FIXME: Add dependencies if required.
    # depends_on("foo")

    def install(self, spec, prefix):
        # FIXME: Unknown build system
                
        make('lib', parallel=False)
        make('prefix={0}'.format(self.prefix), parallel=False)


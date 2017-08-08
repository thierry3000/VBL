##############################################################################
# Copyright (c) 2013-2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# This file is part of Spack.
# Created by Todd Gamblin, tgamblin@llnl.gov, All rights reserved.
# LLNL-CODE-647188
#
# For details, see https://github.com/llnl/spack
# Please also see the NOTICE and LICENSE files for our notice and the LGPL.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License (as
# published by the Free Software Foundation) version 2.1, February 1999.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
# conditions of the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
##############################################################################

# T.F.
# on my maschine this file resides at: 
# /localdisk/thierry/local/spack/var/spack/repos/builtin/packages/vbl-master/package.py
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install vbl-master
#
# You can edit this file again by typing:
#
#     spack edit vbl-master
#
# See the Spack documentation for more information on packaging.
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
from spack import *


class Vbl(CMakePackage):
    """FIXME: Put a proper description of your package here."""
    homepage = "http://http://vbl.ts.infn.it/SiteVBL/index.html"
    url      = "https://github.com/edymil/VBL/archive/master.zip"

    # version('2014-10-08', git='https://github.com/example-project/example.git',commit='9d38cd4e2c94c3cea97d0e2924814acc')
    version('develop', git='https://github.com/edymil/VBL.git',branch='master')
    depends_on('cgal')
    depends_on('boost')

    def cmake_args(self):
        args = ['-DCMAKE_BUILD_TYPE=Release',]
        return args
    #def install(self, spec, prefix):
        # FIXME: Unknown build system
        # make()
    #    make('install')

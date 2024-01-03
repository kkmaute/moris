#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# Trilinos libraries and includes -----------------------------------------
# -------------------------------------------------------------------------
#
# CMake example that uses FIND_PACKAGE(Trilinos ...) to build your C++
# application with Trilinos.  You should know a little bit about CMake
# before reading this example; in particular, you should know how to
# add C++ source files and header files to your project.
#

# Your "do-configure" script that invokes CMake should set
# TRILINOS_PATH to the path to your Trilinos install.
# You do _not_ need to edit this line.

if(NOT TARGET ${MORIS}::trilinos)
	set(MORIS_TRILINOS_TPLS
		"arpack"
		)
	
	foreach(TPL ${MORIS_TRILINOS_TPLS})
		include(${TPL}_new)
		list(APPEND MORIS_TRILINOS_TPLS_TARGETS ${MORIS}::${TPL})
	endforeach()

	add_library(${MORIS}::trilinos INTERFACE IMPORTED GLOBAL)
	target_link_libraries(${MORIS}::trilinos INTERFACE ${MORIS_TRILINOS_TPLS_TARGETS})
endif()


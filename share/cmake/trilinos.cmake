# -----------------------------------------------------------------------------
# Trilinos libraries and includes ------------------------------------------------
# -----------------------------------------------------------------------------
#
# CMake example that uses FIND_PACKAGE(Trilinos ...) to build your C++
# application with Trilinos.  You should know a little bit about CMake
# before reading this example; in particular, you should know how to
# add C++ source files and header files to your project.
#

# Your "do-configure" script that invokes CMake should set
# TRILINOS_PATH to the path to your Trilinos install.
# You do _not_ need to edit this line.
#> what? ^

if(NOT DEFINED ENV{Trilinos_DIR})
    message(FATAL_ERROR 
        "\nPlease set the Trilinos_DIR environment variable. It should be the absolute path to the Trilinos library (ex: /lib/trilinos/gcc-openmpi).\n" )
endif()

if ( NOT MORIS_HAVE_DEBUG )
    set(TRILINOS_PATH $ENV{Trilinos_DIR})
else()
    if(DEFINED ENV{Trilinos_DEBUG_DIR})
        set(TRILINOS_PATH "$ENV{Trilinos_DEBUG_DIR}")
    else()
        message(WARNING 
            "\nMORIS will use the release version of Trilinos unless the Trilinos_DEBUG_DIR environment variable is set.\n" )
        set(TRILINOS_PATH "$ENV{Trilinos_DIR}")
    endif()
endif()

MESSAGE("\nLooking for ${TRILINOS_PATH}\n\n")


if(NOT "${TRILINOS_PATH}" MATCHES "${Trilinos_DIR}")
    set(Trilinos_DIR Trilinos_DIR-NOTFOUND CACHE PATH 
        "Force CMake to find_package(Trilinos ...)"
        FORCE )
endif()

FIND_PACKAGE(Trilinos PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos ${TRILINOS_PATH} NO_CMAKE_ENVIRONMENT_PATH NO_SYSTEM_ENVIRONMENT_PATH)
    

# If FIND_PACKAGE successfully found your Trilinos install, it will
# set the Boolean flag Trilinos_FOUND.  The following IF statement
# fails with a FATAL_ERROR if Trilinos was not found.  If it _was_
# found, it prints out the values of some Trilinos configuration
# details.  You may find them useful for building your application
# that uses Trilinos.
IF(Trilinos_FOUND)
    message("\nFound Trilinos! Details can be found in a config file somewhere...")

#     MESSAGE("\nFound Trilinos!  Here are the details: ")
#         #message("GET LOST, TRILINOS!\n  - <3 MORIS")
#     MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
#     MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
#     MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
#     MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
#     MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
#     MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
#     MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
#     MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
#     MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
#     MESSAGE("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
#     MESSAGE("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
#     MESSAGE("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
#     MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
#     MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
#     MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS = ${Trilinos_Fortran_COMPILER_FLAGS}")
#     MESSAGE("   Trilinos_LINKER = ${Trilinos_LINKER}")
#     MESSAGE("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
#     MESSAGE("   Trilinos_AR = ${Trilinos_AR}")
#     MESSAGE("End of Trilinos details\n")

    list(APPEND MORIS_INCDIRS ${Trilinos_INCLUDE_DIRS})

    # need to add trilinos libraries explicitly otherwise acml.so 
    # is automatically added
    
    foreach( lib ${Trilinos_LIBRARIES})
        list(APPEND MORIS_LDLIBS "${TRILINOS_PATH}/lib/lib${lib}.a")
    endforeach(lib)

    list(APPEND MORIS_INCDIRS ${Trilinos_TPL_INCLUDE_DIRS})
    list(APPEND MORIS_LDLIBS  ${Trilinos_TPL_LIBRARIES})
    
ELSE()
  MESSAGE(FATAL_ERROR "Could not find Trilinos!")
ENDIF()

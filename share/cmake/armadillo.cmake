# -------------------------------------------------------------------------
# Armadillo libraries -----------------------------------------------------
# -------------------------------------------------------------------------

if(NOT ARMADILLO_FOUND_ONCE)
    set(ARMADILLO_ENV_VARS
        $ENV{ARMADILLODIR}
        $ENV{ARMADILLO_DIR}
        $ENV{Armadillo_DIR}
        $ENV{ARMADILLO_ROOT}
        $ENV{Armadillo_ROOT}
        $ENV{ARMADILLO_PATH}
        $ENV{Armadillo_PATH} )

    find_package(Armadillo REQUIRED HINTS ${ARMADILLO_ENV_VARS})
    
    set(ARMADILLO_INCLUDE_DIRS ${ARMADILLO_INCLUDE_DIRS}
        CACHE PATH "Armadillo include directories." )
    set(ARMADILLO_LIBRARY_DIRS ${ARMADILLO_LIBRARY_DIRS}
        CACHE PATH "Armadillo library directories." )
    set(ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARIES}
        CACHE FILEPATH "Armadillo libraries." )
    
    mark_as_advanced(Armadillo_DIR
        ARMADILLO_INCLUDE_DIRS
        ARMADILLO_LIBRARY_DIRS
        ARMADILLO_LIBRARIES )
    
    if(ARMADILLO_FOUND)
        set(ARMADILLO_FOUND_ONCE TRUE CACHE INTERNAL "Armadillo was found.")
    endif()
endif()

# list(APPEND MORIS_DEFINITIONS "-DMORIS_USE_ARMA")
# list(APPEND MORIS_INCDIRS "${ARMADILLO_INCLUDE_DIRS}")
# list(APPEND MORIS_LDFLAGS "${ARMADILLO_LIBRARY_DIRS}")
# list(APPEND MORIS_LDLIBS "-larmadillo")

add_definitions("-DMORIS_USE_ARMA")
include_directories(${ARMADILLO_INCLUDE_DIRS})
link_directories(${ARMADILLO_LIBRARY_DIRS})
list(APPEND MORIS_ARMADILLO_EIGEN_LIBS "-larmadillo")

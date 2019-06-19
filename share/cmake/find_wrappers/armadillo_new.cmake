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
    
    set(MORIS_ARMADILLO_INCLUDE_DIRS ${ARMADILLO_INCLUDE_DIRS}
        CACHE PATH "Armadillo include directories." )
    #set(MORIS_ARMADILLO_LIBRARY_DIRS ${ARMADILLO_LIBRARY_DIRS}
    #    CACHE PATH "Armadillo library directories." )
    #set(MORIS_ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARIES}
    #    CACHE FILEPATH "Armadillo libraries." )
    set(MORIS_ARMADILLO_LIBRARIES "${ARMADILLO_LIBRARY_DIRS}/libarmadillo.a"
        CACHE FILEPATH "Armadillo libraries." )
    set(MORIS_ARMADILLO_DEFINITIONS "-DMORIS_USE_ARMA"
    	CACHE INTERNAL "Moris preprocessor definitions for Armadillo.")
    
    mark_as_advanced(Armadillo_DIR
        MORIS_ARMADILLO_INCLUDE_DIRS
        MORIS_ARMADILLO_LIBRARIES
        MORIS_ARMADILLO_DEFINITIONS )
    
    if(Armadillo_FOUND)
        set(ARMADILLO_FOUND_ONCE TRUE CACHE INTERNAL "Armadillo was found.")
    endif()
    
    message(STATUS "ARMADILLO_INCLUDE_DIRS: ${ARMADILLO_INCLUDE_DIRS}")
	message(STATUS "ARMADILLO_LIBRARY_DIRS: ${ARMADILLO_LIBRARY_DIRS}")
	message(STATUS "ARMADILLO_LIBRARIES: ${ARMADILLO_LIBRARIES}")
endif()

if(NOT TARGET armadillo)
	set(MORIS_ARMADILLO_TPLS
		${ACML_LAPACK_MKL}
		"arpack"
		)
	
	foreach(TPL ${MORIS_ARMADILLO_TPLS})
		include(${MORIS_TPL_DIR}/${TPL}_new.cmake)
	endforeach()

	add_library(armadillo STATIC IMPORTED GLOBAL)
	set_target_properties(armadillo PROPERTIES
		IMPORTED_LOCATION ${MORIS_ARMADILLO_LIBRARIES})
	target_link_libraries(armadillo INTERFACE ${MORIS_ARMADILLO_TPLS})
endif()

#add_definitions("-DMORIS_USE_ARMA")
#include_directories(${MORIS_ARMADILLO_INCLUDE_DIRS})
#link_directories(${MORIS_ARMADILLO_LIBRARY_DIRS})
#list(APPEND MORIS_ARMADILLO_EIGEN_LIBS "-larmadillo")

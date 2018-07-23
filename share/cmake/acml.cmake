# -------------------------------------------------------------------------
# ACML libraries ----------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT ACML_FOUND_ONCE)
    find_package(ACML)
    if(ACML_FOUND)
        set(ACML_FOUND_ONCE TRUE CACHE INTERNAL "ACML was found." FORCE)
    endif()
endif()

# find_package(ACML)

message(STATUS "ACML_LIBRARIES: ${ACML_LIBRARIES}")

# list(APPEND MORIS_DEFINITIONS "-DMORIS_HAVE_ACML")
# list(APPEND MORIS_LDLIBS ${ACML_LIBRARIES})

add_definitions("-DMORIS_HAVE_ACML")
set(MORIS_ACML_LIBS ${ACML_LIBRARIES})
set(MORIS_MATH_LIBS ${ACML_LIBRARIES})

# -------------------------------------------------------------------------
# MKL libraries -----------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT MKL_FOUND_ONCE)
    find_package(MKL)
    if(MKL_FOUND)
        set(MKL_FOUND_ONCE TRUE CACHE INTERNAL "MKL was found.")
    endif()
endif()

message(STATUS "MKL_LIBRARIES: ${MKL_LIBRARIES}")

# list(APPEND MORIS_DEFINITIONS "-DMORIS_HAVE_MKL")
# list(APPEND MORIS_INCDIRS ${MKL_INCLUDE_DIRS})
# list(APPEND MORIS_LDLIBS ${MKL_LIBRARIES})

if(MORIS_USE_MKL)
    add_definitions("-DMORIS_HAVE_MKL")
endif()
include_directories(${MKL_INCLUDE_DIRS})
set(MORIS_MKL_LIBS ${MKL_LIBRARIES})
set(MORIS_MATH_LIBS ${MKL_LIBRARIES})

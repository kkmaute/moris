# -------------------------------------------------------------------------
# SNOPT libraries ---------------------------------------------------------
# -------------------------------------------------------------------------

find_package(SNOPT)

message(STATUS "SNOPT_LIBRARIES: ${SNOPT_LIBRARIES}")

list(APPEND MORIS_LDFLAGS "${SNOPT_LIBRARY_DIRS}")
list(APPEND MORIS_LDLIBS "${SNOPT_LIBRARIES}")

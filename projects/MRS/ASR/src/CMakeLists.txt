# Assert Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    assert.hpp
    fn_assert.hpp
    fn_log_assert.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert Target
add_library(${ASR}-target INTERFACE)
target_include_directories(${ASR}-target INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert Library

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert Config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert Installation

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${MRS}/${ASR} )


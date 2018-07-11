# MORIS mains -------------------------------------------------------------
# -------------------------------------------------------------------------

# Note: These execiutables are to be replaced by project-specific mains
#       ...eventaually.

# Set executable output location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# List mains
set(MAINS
    main.cpp )

# List include directories
set(MAIN_INCLUDES
    ${COM}
    ${IOS}
    ${LNA})

foreach(MAIN_INCLUDE ${MAIN_INCLUDES})
    include_directories(../${MAIN_INCLUDE}/src)
endforeach()

# List executable dependencies
set(MAIN_DEPENDENCIES
    ${COM}-lib)

# Create and link main executables
foreach(MAIN ${MAINS})
    string(REPLACE ".cpp" "-${EXECUTABLE}" EXE_NAME ${MAIN})
    add_executable(${EXE_NAME} ${MAIN})
    target_link_libraries(${EXE_NAME} ${MAIN_DEPENDENCIES})
    list(APPEND MAIN_EXECUTABLES ${EXE_NAME})
endforeach()

# Install executables
install(TARGETS ${MAIN_EXECUTABLES}
    DESTINATION ${EXE_INSTALL_DIR} )


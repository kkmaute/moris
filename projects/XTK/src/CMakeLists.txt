# XTK Source --------------------------------------------------------------                                                                                                                                                                                                    
# -------------------------------------------------------------------------                                                                                                                                                                                                    

# Set version                                                                                                                                                                                                                                                                  
set(${XTK}_VERSION ${MORIS_VERSION})                                                                                                                                                                                                                                           

# List header files                                                                                                                                                                                                                                                            
set(HEADERS)                                                                                                                                                                                                                                                                   

# List header dependencies                                                                                                                                                                                                                                                     
set(XTK_INCLUDES                                                                                                                                                                                                                                                               
    "algorithms"                                                                                                                                                                                                                                                               
    "assert"                                                                                                                                                                                                                                                                   
    "containers"                                                                                                                                                                                                                                                               
    "core"                                                                                                                                                                                                                                                                     
    "geomeng"                                                                                                                                                                                                                                                                  
    "geometry"                                                                                                                                                                                                                                                                 
    "ios"                                                                                                                                                                                                                                                                      
    "linalg"                                                                                                                                                                                                                                                                   
    "mesh"                                                                                                                                                                                                                                                                     
    "tools"                                                                                                                                                                                                                                                                    
    "topology"                                                                                                                                                                                                                                                                 
    "xtk"                                                                                                                                                                                                                                                                      
    )                                                                                                                                                                                                                                                                          

# Include depedency directories                                                                                                                                                                                                                                                
foreach(XTK_INCLUDE ${XTK_INCLUDES})                                                                                                                                                                                                                                           
    include_directories(${MORIS_PACKAGE_DIR}/XTK/src/${XTK_INCLUDE})                                                                                                                                                                                                           
endforeach()                                                                                                                                                                                                                                                                   

include_directories(${MORIS_PACKAGE_DIR}/XTK/include)                                                                                                                                                                                                                          

include_directories(.) 


# List header dependencies
set(XTK_INCLUDES
    ${LINALG}
    )


# Include depedency directories
foreach(XTK_INCLUDE ${XTK_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${XTK_INCLUDE}/src)
endforeach()
                                                                                                                                                                                                                                                        

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                                                                                                                                    
# XTK Library                                                                                                                                                                                                                                                                  

# List source files                                                                                                                                                                                                                                                            
# set(LIB_SOURCES
#       linalg/cl_XTK_Linalg_Enums.${XTK_SOURCE_EXT)

# List library dependencies
# set(LIB_DEPENDENCIES "")

# Set the output path for library
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
#add_library(${XTK}-lib STATIC ${LIB_SOURCES})
#target_link_libraries(${XTK}-lib ${LIB_DEPENDENCIES})
#set_target_properties(${XTK}-lib PROPERTIES OUTPUT_NAME ${XTK})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# XTK Executable
# List source files
set(EXE_SOURCES
        main.cpp)

# List library dependencies
set(EXE_DEPENDENCIES 
    ${MORIS_LDLIBS} 
    ${MORIS_ARMADILLO_EIGEN_LIBS} 	
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_SUPERLU_LIBS})
 
# Set the output path for library
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${bin})

# Create library
add_executable(${XTK} ${EXE_SOURCES})
target_link_libraries(${XTK} ${EXE_DEPENDENCIES})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# XTK Config Files

# # Build *Config and *ConfigVersion files
# include(CMakePackageConfigHelpers)
# 
# configure_package_config_file(
#     ${XTK}Config.cmake.in
#     ${${XTK}_CONFIG_BUILD_DIR}/${XTK}Config.cmake
#     INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${XTK}
#     PATH_VARS LIB_INSTALL_DIR )
# 
# write_basic_package_version_file(
#     ${${XTK}_CONFIG_BUILD_DIR}/${XTK}ConfigVersion.cmake
#     COMPATIBILITY AnyNewerVersion )
# 
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# # XTK Installation
# 
# # Install *Config and *ConfigVersion files
# install(FILES ${${XTK}_CONFIG_BUILD_DIR}/${XTK}Config.cmake 
#     ${${XTK}_CONFIG_BUILD_DIR}/${XTK}ConfigVersion.cmake
#     DESTINATION ${CONFIG_INSTALL_DIR}/${XTK} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${XTK} )

# Install library
install(TARGETS ${XTK}
    DESTINATION ${EXE_INSTALL_DIR}/${XTK}
    EXPORT ${XTK}Targets )

# Export library in to *Targets file
install(EXPORT ${XTK}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${XTK} )
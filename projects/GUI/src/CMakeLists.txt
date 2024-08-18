cmake_minimum_required(VERSION 3.5)
# Set executable output location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

list(REMOVE_DUPLICATES MAIN_TPL_DEPENDENCIES)
foreach(MAIN_TPL ${MAIN_TPL_DEPENDENCIES})
    include(${MORIS_TPL_DIR}/${MAIN_TPL}_new.cmake)
endforeach()

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)

message("QT version: ${QT_VERSION_MAJOR}")
  
set(PROJECT_SOURCES
       main_gui.cpp
       main.cpp 
       cl_line_edit.cpp
       cl_combo_box.cpp
       cl_bool_combo_box.cpp
       TestWindow.cpp
       cl_double_spin_box.cpp
       cl_int_spin_box.cpp
       cl_pair_box.cpp
       cl_tree_widget_item.cpp
       fn_read_file_dialog.cpp
       fn_write_file_dialog.cpp
       cl_group_box.cpp
       )

set(HEADERS
       main_gui.hpp
       cl_line_edit.hpp
       cl_combo_box.hpp
       cl_bool_combo_box.hpp
       TestWindow.hpp
       cl_double_spin_box.hpp
       cl_int_spin_box.hpp
       cl_pair_box.hpp
       cl_tree_widget_item.hpp
       fn_read_file_dialog.hpp
       fn_write_file_dialog.hpp
       cl_group_box.hpp

)
 
 # List include directories
set(MAIN_INCLUDES
    ${COM}
    ${IOS}
    ${WRK}
    ${LINALG}
        )

foreach(MAIN_INCLUDE ${MAIN_INCLUDES})
    include_directories(../${MAIN_INCLUDE}/src)
endforeach()
 
 # List executable dependencies
set(MAIN_DEPENDENCIES
    ${HMR}-lib
    ${COM}-lib
    ${IOS}-lib
    ${GEN}-lib
    ${WRK}-lib
    ${MORIS_PETSC_LIBS}
    ${MORIS_BOOST_LIBS}
    ${MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_LDLIBS}
    ${MORIS_BASE_LIBS}
    )

string(REPLACE ".cpp" "-${EXECUTABLE}" EXE_NAME ${PROJECT_SOURCES})

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(moris_gui
        MANUAL_FINALIZATION
        ${PROJECT_SOURCES}
    )
else()
    add_executable(moris_gui
        ${PROJECT_SOURCES}
    )
endif()

target_link_libraries(moris_gui PRIVATE Qt${QT_VERSION_MAJOR}::Widgets)
target_link_libraries(moris_gui PRIVATE ${MAIN_DEPENDENCIES})

# following lines commented out as cause error on linux
#if(QT_VERSION_MAJOR EQUAL 6)
#    qt_finalize_executable(moris_gui)
#endif()

list(APPEND MAIN_EXECUTABLES moris_gui)

# Install executables
install(TARGETS ${MAIN_EXECUTABLES}
    DESTINATION ${EXE_INSTALL_DIR} )

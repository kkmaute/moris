
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${BIN}")

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 REQUIRED COMPONENTS Test Core Widgets)

set(TEST_SOURCES
    test_main.cpp
    ut_gui_startup.cpp
    ut_cl_line_edit.cpp
    ut_cl_int_spin_box.cpp
    ut_cl_double_spin_box.cpp
    ut_cl_combo_box.cpp
    ut_cl_bool_combo_box.cpp
)

add_executable(${GUI}-test  ${TEST_SOURCES})
# target_compile_options(${GUI}-test PRIVATE -fno-qt-genetrate-moc)

# set_property(SOURCE ut_gui_startup.cpp
#     PROPERTY SKIP_AUTOMOC ON)

target_link_libraries(${GUI}-test
    PRIVATE
        gui_core
        Qt6::Test # could update this such that it is able to detect the major version and chose that
        Qt6::Core
        Qt6::Gui
        Qt6::Widgets
        test-libs
)


enable_testing()

add_test(NAME ${GUI}-test
         COMMAND $<TARGET_FILE:${GUI}-test>)
# List header files
set(HEADERS

    cl_FEM_Set.hpp
    cl_FEM_NodeProxy.hpp
    cl_FEM_Node_Base.hpp
    cl_FEM_Node.hpp
    cl_FEM_Field.hpp        
    fn_FEM_Side_Coordinate_Map.hpp
    cl_FEM_Property.hpp
    cl_FEM_Set_User_Info.hpp
    cl_FEM_Phase_User_Info.hpp
    cl_FEM_Model.hpp
    fn_FEM_Check.hpp
    fn_FEM_FD_Scheme.hpp
    )

# List library source files
#set(LIB_SOURCES
#cl_FEM_Set.cpp
#cl_FEM_Property.cpp
#cl_FEM_Model.cpp
#cl_FEM_Field.cpp 
#    )
set(LIB_SOURCES
#cl_FEM_Set.cpp
#cl_FEM_Property.cpp
cl_FEM_Model.cpp
#cl_FEM_Field.cpp 
    )

# Create library

#include(${CMAKE_SOURCE_DIR}/projects/tgt_includes.cmake)
list(TRANSFORM LIB_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
list(TRANSFORM HEADERS PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

#target_sources(morisCombinedLib PUBLIC ${LIB_SOURCES} ${HEADERS})
target_sources(morisCombinedLib PUBLIC ${HEADERS})

target_include_directories(morisCombinedLib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include> )


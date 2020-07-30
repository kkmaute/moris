#ifndef MORIS_FN_CREATE_GEOMETRY_HPP
#define MORIS_FN_CREATE_GEOMETRY_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_Param_List.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace ge
    {
        /**
         * Higher-level call for creating a cell of geometries, which ensures that the ADVs are resized properly
         *
         * @param aGeometryParameterLists Parameter lists for creating geometry classes
         * @param aADVs Reference to the initial ADV vector
         * @param aLibrary Pointer to library for loading user-defined functions
         * @return Pointers to created Geometry classes
         */
        Cell<std::shared_ptr<Geometry>> create_geometries(
                Cell<ParameterList>         aGeometryParameterLists,
                Matrix<DDRMat>&             aADVs,
                std::shared_ptr<Library_IO> aLibrary = nullptr);

        /**
         * Creates an instance of the specified Geometry class and returns a shared pointer to it.
         *
         * @param aGeometryParameterList Parameter list for creating a geometry class
         * @param aADVs Reference to the initial ADV vector
         * @param aLibrary Pointer to library for loading user-defined functions
         * @return Pointer to specific Geometry class
         */
        std::shared_ptr<Geometry> create_geometry(
                ParameterList               aGeometryParameterList,
                Matrix<DDRMat>&             aADVs,
                std::shared_ptr<Library_IO> aLibrary = nullptr);

        /**
         * Sets the geometry variables which depend on ADVs. Used by create_geometry().
         *
         * @param aGeometryParameterList Parameter list for creating a geometry class
         * @param aNumADVs The number of total ADVs
         * @param aGeometryVariableIndices Indices of geometry variables to fill with ADVs
         * @param aADVIndices Indices of ADVs for filling the geometry variables
         */
        void set_geometry_variable_inputs(
                ParameterList  aGeometryParameterList,
                uint           aNumADVs,
                Matrix<DDUMat>& aGeometryVariableIndices,
                Matrix<DDUMat>& aADVIndices);

    }
}

#endif //MORIS_FN_CREATE_GEOMETRY_HPP

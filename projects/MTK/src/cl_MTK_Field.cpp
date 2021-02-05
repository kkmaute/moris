#include <iostream>
#include <cstdio>

// HD5 c-interface
#include "hdf5.h"

#include "fn_save_matrix_to_binary_file.hpp"

#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Field.hpp"
#include "linalg_typedefs.hpp"

#include "HDF5_Tools.hpp"

#include "cl_MTK_Mesh.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace mtk
    {

        //------------------------------------------------------------------------------

        Field::Field(sint        aDiscretizationMeshIndex,
                     std::string aName)
                : mDiscretizationMeshIndex(aDiscretizationMeshIndex)
                , mName(aName)
        {
        }

        //------------------------------------------------------------------------------

        Field::~Field()
        {
        }

        //------------------------------------------------------------------------------

        const Matrix<DDRMat>& Field::get_nodal_values(mtk::Mesh* aMesh)
        {
            mNodalValues.set_size(aMesh->get_num_nodes(), 1);
            for (uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++)
            {
                mNodalValues(tNodeIndex) = this->get_field_value(tNodeIndex, aMesh->get_node_coordinate(tNodeIndex));
            }
            return mNodalValues;
        }

        //------------------------------------------------------------------------------

        uint Field::get_number_of_dimensions()
        {
            return mNumberOfDimensions;
        }

        //------------------------------------------------------------------------------

        uint Field::get_discretization_mesh_index()
        {
            return mDiscretizationMeshIndex;
        }

        //------------------------------------------------------------------------------

        std::string Field::get_name()
        {
            return mName;
        }

        //------------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */

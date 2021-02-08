#include "cl_MTK_Discrete_Field.hpp"
#include "fn_dot.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris
{
    namespace mtk
    {

        //--------------------------------------------------------------------------------------------------------------

        Discrete_Field::Discrete_Field(
                Mesh_Pair aMeshPair,
                uint      aDiscretizationMeshIndex)
                : Field(aDiscretizationMeshIndex)
                , mMeshPair(aMeshPair)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Discrete_Field::~Discrete_Field()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Discrete_Field::get_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            return mNodalValues(aNodeIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Discrete_Field::set_nodal_values(const Matrix<DDRMat>& aNodalValues)
        {
            MORIS_ERROR(aNodalValues.length() == mMeshPair.mInterpolationMesh->get_num_nodes(),
                        "B-spline coefficients set to a field must match the number of coefficients on the mesh.");
            mNodalValues = aNodalValues;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Discrete_Field::get_nodal_values()
        {
            return mNodalValues;
        }

        //--------------------------------------------------------------------------------------------------------------

        Mesh_Pair Discrete_Field::get_mesh_pair()
        {
            return mMeshPair;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}


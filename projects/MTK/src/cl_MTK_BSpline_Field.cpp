#include "cl_MTK_BSpline_Field.hpp"
#include "fn_dot.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris
{
    namespace mtk
    {

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Field::BSpline_Field(
                Mesh_Pair aMeshPair,
                uint      aDiscretizationMeshIndex)
                : Field(aDiscretizationMeshIndex)
                , mMeshPair(aMeshPair)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Field::~BSpline_Field()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void BSpline_Field::set_coefficients(const Matrix<DDRMat>& aCoefficients)
        {
            MORIS_ERROR(aCoefficients.length() == mMeshPair.mInterpolationMesh->get_num_coeffs(get_discretization_mesh_index()),
                        "B-spline coefficients set to a field must match the number of coefficients on the mesh.");
            mCoefficients = aCoefficients;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& BSpline_Field::get_nodal_values()
        {
            return Field::get_nodal_values(mMeshPair.mInterpolationMesh);
        }

        //--------------------------------------------------------------------------------------------------------------

        real BSpline_Field::get_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            // Get mesh
            Mesh* tMesh = mMeshPair.mInterpolationMesh;

            // Get B-spline information
//            sint tNodeID = tMesh->get_glb_entity_id_from_entity_loc_index(
//                    aNodeIndex,
//                    EntityRank::NODE,
//                    this->get_discretization_mesh_index());
            Matrix<IndexMat> tBSplineIndices = tMesh->get_bspline_inds_of_node_loc_ind(aNodeIndex, this->get_discretization_mesh_index());
            Matrix<DDRMat> tMatrix = tMesh->get_t_matrix_of_node_loc_ind(aNodeIndex, this->get_discretization_mesh_index());

            // Multiply T-matrix
            real tFieldValue = 0.0;
            for (uint tBSpline = 0; tBSpline < tBSplineIndices.length(); tBSpline++)
            {
                tFieldValue += tMatrix(tBSpline) * mCoefficients(tBSplineIndices(tBSpline));
            }

            return tFieldValue;
        }

        //--------------------------------------------------------------------------------------------------------------

        Mesh_Pair BSpline_Field::get_mesh_pair()
        {
            return mMeshPair;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}


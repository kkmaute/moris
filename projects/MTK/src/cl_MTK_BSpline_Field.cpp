#include "cl_MTK_BSpline_Field.hpp"
#include "fn_dot.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace mtk
    {

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Field::BSpline_Field(
                std::shared_ptr<mtk::Mesh_Manager> aMeshManager,
                uint                               aMeshIndex,
                uint                               aDiscretizationMeshIndex )
                : Field(aDiscretizationMeshIndex)
                , mMeshManager(aMeshManager)
                , mMeshIndex(aMeshIndex)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        BSpline_Field::~BSpline_Field()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void BSpline_Field::set_coefficients(const Matrix<DDRMat>& aCoefficients)
        {
            MORIS_ERROR(aCoefficients.length() == mMeshManager->get_interpolation_mesh(mMeshIndex)->get_num_coeffs(get_discretization_mesh_index()),
                        "B-spline coefficients set to a field must match the number of coefficients on the mesh.");
            mCoefficients = aCoefficients;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& BSpline_Field::get_nodal_values()
        {
            return Field::get_nodal_values(mMeshManager->get_interpolation_mesh(mMeshIndex));
        }

        //--------------------------------------------------------------------------------------------------------------

        real BSpline_Field::get_field_value(
                uint                  aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            // Get mesh
            Mesh* tMesh = mMeshManager->get_interpolation_mesh(mMeshIndex);

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

        std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > BSpline_Field::get_mesh_pair()
        {
            return std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> >( mMeshIndex, mMeshManager );
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}


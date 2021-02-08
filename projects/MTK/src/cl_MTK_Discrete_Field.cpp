#include "cl_MTK_Discrete_Field.hpp"
#include "fn_dot.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace mtk
    {

        //--------------------------------------------------------------------------------------------------------------

        Discrete_Field::Discrete_Field(
                std::shared_ptr<mtk::Mesh_Manager> aMeshManager,
                uint                               aMeshIndex,
                uint                               aDiscretizationMeshIndex )
                : Field(aDiscretizationMeshIndex)
                , mMeshManager(aMeshManager)
                , mMeshIndex(aMeshIndex)
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
            MORIS_ERROR(aNodalValues.length() == mMeshManager->get_interpolation_mesh(mMeshIndex)->get_num_nodes(),
                        "B-spline coefficients set to a field must match the number of coefficients on the mesh.");
            mNodalValues = aNodalValues;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Discrete_Field::get_nodal_values()
        {
            return mNodalValues;
        }

        //--------------------------------------------------------------------------------------------------------------

        std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> > Discrete_Field::get_mesh_pair()
        {
            return std::pair< moris_index, std::shared_ptr<mtk::Mesh_Manager> >( mMeshIndex, mMeshManager );
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}


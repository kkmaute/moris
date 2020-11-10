#include "cl_GEN_Field_Discrete_Interpolation.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Field_Discrete_Interpolation::Field_Discrete_Interpolation(mtk::Mesh* aMesh)
                : mMesh(aMesh)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Field_Discrete_Interpolation::get_field_value(
                uint aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            return this->get_field_value(mMesh->get_base_node_index(aNodeIndex));
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Field_Discrete_Interpolation::get_field_sensitivities(
                uint aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            return this->get_field_sensitivities(mMesh->get_base_node_index(aNodeIndex));
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Field_Discrete_Interpolation::get_determining_adv_ids(
                uint aNodeIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            return this->get_determining_adv_ids(mMesh->get_base_node_index(aNodeIndex));
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Field_Discrete_Interpolation::get_determining_adv_ids(uint aNodeIndex)
        {
            return Field::get_determining_adv_ids(aNodeIndex, {{}});
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

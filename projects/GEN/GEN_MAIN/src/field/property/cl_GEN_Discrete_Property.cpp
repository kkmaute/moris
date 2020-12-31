#include "cl_GEN_Discrete_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        real Discrete_Property::get_field_value(uint aNodeIndex)
        {
            return *mFieldVariables(aNodeIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Discrete_Property::get_field_sensitivities(uint aNodeIndex)
        {
            return mSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Discrete_Property::get_determining_adv_ids(uint aNodeIndex)
        {
            return {{(sint)aNodeIndex}};
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

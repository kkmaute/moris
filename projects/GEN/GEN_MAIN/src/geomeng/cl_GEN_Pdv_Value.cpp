#include "cl_GEN_Pdv_Value.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Pdv_Value::Pdv_Value(real aValue)
        : mValue(aValue)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Pdv_Value::get_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            return mValue;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Value::get_sensitivity(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

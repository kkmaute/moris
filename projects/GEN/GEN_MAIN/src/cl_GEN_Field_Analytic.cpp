#include "cl_GEN_Field_Analytic.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Field_Analytic::Field_Analytic()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Field_Analytic::get_field_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            return this->get_field_value(aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Field_Analytic::get_field_sensitivities(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            return this->get_field_sensitivities(aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

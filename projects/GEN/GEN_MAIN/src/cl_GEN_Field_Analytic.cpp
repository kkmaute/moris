//
// Created by christopherson on 5/19/20.
//

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

        real Field_Analytic::get_field_value(
                uint                  aIndex,
                const Matrix<DDRMat>& aCoordinates)
        {
            return this->get_field_value(aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field_Analytic::evaluate_sensitivities(
                uint                  aIndex,
                const Matrix<DDRMat>& aCoordinates,
                Matrix<DDRMat>&       aSensitivities)
        {
            this->evaluate_sensitivities(aCoordinates, aSensitivities);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

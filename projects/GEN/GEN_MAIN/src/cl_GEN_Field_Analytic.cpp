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

        real Field_Analytic::evaluate_field_value(      uint            aIndex,
                                                  const Matrix<DDRMat>& aCoordinates)
        {
            return this->evaluate_field_value(aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field_Analytic::evaluate_all_sensitivities(      uint            aIndex,
                                                        const Matrix<DDRMat>& aCoordinates,
                                                              Matrix<DDRMat>& aSensitivities)
        {
            this->evaluate_all_sensitivities(aCoordinates, aSensitivities);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
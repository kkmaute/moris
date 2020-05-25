//
// Created by christopherson on 5/19/20.
//

#include "cl_GEN_Field_Discrete.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Field_Discrete::Field_Discrete()
        {

        }

        //--------------------------------------------------------------------------------------------------------------

        real Field_Discrete::evaluate_field_value(      uint            aIndex,
                                                  const Matrix<DDRMat>& aCoordinates)
        {
            return this->evaluate_field_value(aIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Field_Discrete::evaluate_all_sensitivities(      uint            aIndex,
                                                        const Matrix<DDRMat>& aCoordinates,
                                                              Matrix<DDRMat>& aSensitivities)
        {
            this->evaluate_all_sensitivities(aIndex, aSensitivities);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
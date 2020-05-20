//
// Created by christopherson on 5/19/20.
//

#include "cl_GEN_Pdv_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Pdv_Property::Pdv_Property(std::shared_ptr<Property> aPropertyPointer)
        : Pdv(0.0),
          mProperty(aPropertyPointer)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real Pdv_Property::get_value(uint aNodeIndex, Matrix<DDRMat> aCoordinates)
        {
            return mProperty->evaluate_field_value(aNodeIndex, aCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

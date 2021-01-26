#include "cl_GEN_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Property::Property(Property_Parameters aParameters)
                : mParameters(aParameters)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Property::Property(std::shared_ptr<Property> aProperty)
                : mParameters(aProperty->mParameters)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

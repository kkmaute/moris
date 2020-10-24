#include "cl_GEN_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Property::Property(Cell<std::shared_ptr<Property>> aPropertyDependencies)
                           : mPropertyDependencies(aPropertyDependencies)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

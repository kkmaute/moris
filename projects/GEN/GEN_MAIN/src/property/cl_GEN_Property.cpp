#include "cl_GEN_Property.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Property::Property(Cell<std::shared_ptr<Field>> aFieldDependencies)
                : mFieldDependencies(aFieldDependencies)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

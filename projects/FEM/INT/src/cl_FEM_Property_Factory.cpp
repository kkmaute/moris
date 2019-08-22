#include "assert.hpp"
#include "cl_FEM_Property_Factory.hpp"                               //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Property * Property_Factory::create_property( Property_Type aPropertyType )
        {
            Property * tProperty = nullptr;

            switch( aPropertyType )
            {
                default:
                    MORIS_ERROR( false, " Property_Factory::create_Property - No property type specified. " );
                    break;
            }
            return tProperty;
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

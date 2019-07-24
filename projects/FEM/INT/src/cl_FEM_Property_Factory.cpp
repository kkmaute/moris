#include "assert.hpp"
#include "cl_FEM_Property_Factory.hpp"                               //FEM/INT/src
#include "cl_FEM_Property_Temp_Dirichlet.hpp"                                    //FEM/INT/src

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
                case ( Property_Type::TEMP_DIRICHLET ):
                    tProperty = new Property_Temp_Dirichlet();
                    break;

                default:
                    MORIS_ERROR( false, " Property_Factory::create_Property - No property type specified. " );
                    break;
            }
            return tProperty;
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

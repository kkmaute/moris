#include "assert.hpp"
#include "cl_FEM_Element_Factory.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Element_Factory::Element_Factory(){}

//------------------------------------------------------------------------------

        Element_Factory::~Element_Factory(){}

//------------------------------------------------------------------------------

        MSI::Equation_Object * Element_Factory::create_element( Element_Type         aElementType,
                                                                mtk::Cell          * aCell,
                                                                Cell< IWG* >       & aIWGs,
                                                                Cell< Node_Base* > & aNodes )
        {
            MSI::Equation_Object * tElement;

            switch( aElementType )
            {
                case ( Element_Type::UNDEFINED ):
                    tElement = new Element( aCell, aIWGs, aNodes );
                    break;

                default:
                    MORIS_ERROR( false, "Element_Factory::create_element - No element type specified" );
                    break;
            }
            return tElement;
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

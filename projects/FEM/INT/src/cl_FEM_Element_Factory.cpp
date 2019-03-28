#include "assert.hpp"
#include "cl_FEM_Element_Factory.hpp"      //FEM/INT/src
#include "cl_FEM_Element_Bulk.hpp"         //FEM/INT/src
#include "cl_FEM_Element_Sideset.hpp"      //FEM/INT/src
#include "cl_FEM_Element_Time_Sideset.hpp" //FEM/INT/src

#include "cl_MSI_Equation_Object.hpp" //FEM/MSI/src

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
                                                                mtk::Cell           * aCell,
                                                                Cell< IWG* >        & aIWGs,
                                                                Cell< Node_Base* >  & aNodes )
        {
            MSI::Equation_Object * tElement = nullptr;

            switch( aElementType )
            {
                case ( fem::Element_Type::BULK ):
                    tElement = new Element_Bulk( aCell, aIWGs, aNodes );
                    break;

                case ( fem::Element_Type::SIDESET ):
                    tElement = new Element_Sideset( aCell, aIWGs, aNodes );
                    break;

                case ( fem::Element_Type::TIME_SIDESET ):
                    tElement = new Element_Time_Sideset( aCell, aIWGs, aNodes );
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

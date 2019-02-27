#include "assert.hpp"
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Factory::IWG_Factory(){}

//------------------------------------------------------------------------------

        IWG_Factory::~IWG_Factory(){}

//------------------------------------------------------------------------------

        Cell< IWG* > IWG_Factory::create_IWGs( Element_Type aElementType )
        {
            Cell< IWG* > tIWGs;

            switch( aElementType )
            {
                case ( Element_Type::L2 ):
                    tIWGs.resize( 1, nullptr );
                    tIWGs( 0 ) = new IWG_L2();
                    break;

                case ( Element_Type::HJ ):
                    tIWGs.resize( 1, nullptr );
                    tIWGs( 0 ) = new IWG_Helmholtz_Bulk();
                    break;

                default:
                    MORIS_ERROR( false, " IWG_Factory::create_IWGs - No element type specified. " );
                    break;
            }
            return tIWGs;
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

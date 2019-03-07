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

        IWG* IWG_Factory::create_IWGs( IWG_Type aIWGType )
        {
            IWG* tIWG;

            switch( aIWGType )
            {
                case ( IWG_Type::L2 ):
                    tIWG = new IWG_L2();
                    break;

                case ( IWG_Type::HELMHOLTZ ):
                    tIWG = new IWG_Helmholtz_Bulk2();
                    break;

                case ( IWG_Type::HJ ):
                    tIWG = new IWG_Hamilton_Jacobi_Bulk2();
                    break;

                case ( IWG_Type::LSNORMAL ):
                    tIWG = new IWG_LSNormal_Bulk();
                    break;

                case ( IWG_Type::OLSSON ):
                    tIWG = new IWG_Olsson_CLS_Bulk();
                    break;

                case ( IWG_Type::SPATIALDIFF ):
                    tIWG = new IWG_Spatial_Diffusion_Bulk();
                    break;

                default:
                    MORIS_ERROR( false, " IWG_Factory::create_IWGs - No IWG type specified. " );
                    break;
            }
            return tIWG;
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

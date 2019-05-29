#include "assert.hpp"
#include "cl_FEM_IWG_Factory.hpp"                               //FEM/INT/src
#include "cl_FEM_IWG_L2.hpp"                                    //FEM/INT/src
#include "cl_FEM_IWG_Helmholtz_Bulk.hpp"                        //FEM/INT/src
#include "cl_FEM_IWG_Helmholtz_Bulk2.hpp"                       //FEM/INT/src
#include "cl_FEM_IWG_Helmholtz_Interface.hpp"                   //FEM/INT/src
#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk.hpp"                  //FEM/INT/src
#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk2.hpp"                 //FEM/INT/src
#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk_Test.hpp"             //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Bulk.hpp"      //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Dirichlet.hpp" //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Neumann.hpp"   //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Ghost.hpp"     //FEM/INT/src
#include "cl_FEM_IWG_LSNormal_Bulk.hpp"                         //FEM/INT/src
#include "cl_FEM_IWG_Olsson_CLS_Bulk.hpp"                       //FEM/INT/src
#include "cl_FEM_IWG_Olsson_CLS_Interface.hpp"                  //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG* IWG_Factory::create_IWGs( IWG_Type aIWGType )
        {
            IWG* tIWG = nullptr;

            switch( aIWGType )
            {
                case ( IWG_Type::L2 ):
                    tIWG = new IWG_L2();
                    break;

                case ( IWG_Type::HELMHOLTZ ):
                    tIWG = new IWG_Helmholtz_Bulk();
                    break;

                case ( IWG_Type::HJ ):
                    tIWG = new IWG_Hamilton_Jacobi_Bulk2();
                    break;

                case ( IWG_Type::HJTEST ):
                    tIWG = new IWG_Hamilton_Jacobi_Bulk_Test();
                    break;

                case ( IWG_Type::LSNORMAL ):
                    tIWG = new IWG_LSNormal_Bulk();
                    break;

                case ( IWG_Type::OLSSON ):
                    tIWG = new IWG_Olsson_CLS_Bulk();
                    break;

                case ( IWG_Type::SPATIALDIFF_BULK ):
                    tIWG = new IWG_Isotropic_Spatial_Diffusion_Bulk();
                    break;

                case ( IWG_Type::SPATIALDIFF_DIRICHLET ):
                    tIWG = new IWG_Isotropic_Spatial_Diffusion_Dirichlet();
                    break;

                case ( IWG_Type::SPATIALDIFF_NEUMANN ):
                    tIWG = new IWG_Isotropic_Spatial_Diffusion_Neumann();
                    break;

                case ( IWG_Type::SPATIALDIFF_GHOST ):
                    tIWG = new IWG_Isotropic_Spatial_Diffusion_Ghost();
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

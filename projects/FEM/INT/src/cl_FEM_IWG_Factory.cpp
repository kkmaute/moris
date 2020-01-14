#include <memory>
#include "assert.hpp"
#include "cl_FEM_IWG_Factory.hpp"                               //FEM/INT/src
#include "cl_FEM_IWG_L2.hpp"                                    //FEM/INT/src
//#include "cl_FEM_IWG_Helmholtz_Bulk.hpp"                        //FEM/INT/src
//#include "cl_FEM_IWG_Helmholtz_Bulk2.hpp"                       //FEM/INT/src
//#include "cl_FEM_IWG_Helmholtz_Interface.hpp"                   //FEM/INT/src
//#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk.hpp"                  //FEM/INT/src
//#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk2.hpp"                 //FEM/INT/src
#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk_Test.hpp"             //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Bulk.hpp"      //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Dirichlet.hpp" //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Neumann.hpp"   //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Ghost.hpp"     //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost.hpp"     //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Interface.hpp" //FEM/INT/src
//#include "cl_FEM_IWG_LSNormal_Bulk.hpp"                         //FEM/INT/src
//#include "cl_FEM_IWG_Olsson_CLS_Bulk.hpp"                       //FEM/INT/src
//#include "cl_FEM_IWG_Olsson_CLS_Interface.hpp"                  //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.hpp"           //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Bulk.hpp"  //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.hpp"      //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Dirichlet.hpp"      //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Interface.hpp"        //FEM/INT/src
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Neumann.hpp"        //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        std::shared_ptr< IWG > IWG_Factory::create_IWG( IWG_Type aIWGType )
        {
            switch( aIWGType )
            {
                case ( IWG_Type::L2 ):
                    return std::make_shared< IWG_L2 >();

//                case ( IWG_Type::HELMHOLTZ ):
//                    return std::make_shared< IWG_Helmholtz_Bulk >();
//
//                case ( IWG_Type::HJ ):
//                    return std::make_shared< IWG_Hamilton_Jacobi_Bulk2 >();
//
                case ( IWG_Type::HJTEST ):
                    return std::make_shared< IWG_Hamilton_Jacobi_Bulk_Test >();

//                case ( IWG_Type::LSNORMAL ):
//                    return std::make_shared< IWG_LSNormal_Bulk >();
//
//                case ( IWG_Type::OLSSON ):
//                    return std::make_shared< IWG_Olsson_CLS_Bulk >();

                case ( IWG_Type::SPATIALDIFF_BULK ):
                    return std::make_shared< IWG_Isotropic_Spatial_Diffusion_Bulk >();

                case ( IWG_Type::SPATIALDIFF_DIRICHLET ):
                    return std::make_shared< IWG_Isotropic_Spatial_Diffusion_Dirichlet >();

                case ( IWG_Type::SPATIALDIFF_NEUMANN ):
                    return std::make_shared< IWG_Isotropic_Spatial_Diffusion_Neumann >();

                case ( IWG_Type::SPATIALDIFF_GHOST ):
                    return std::make_shared< IWG_Isotropic_Spatial_Diffusion_Ghost >();

                case ( IWG_Type::SPATIALDIFF_VW_GHOST ):
                    return std::make_shared< IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost >();

                case ( IWG_Type::SPATIALDIFF_INTERFACE ):
                    return std::make_shared< IWG_Isotropic_Spatial_Diffusion_Interface >();

                case ( IWG_Type::STRUC_LINEAR_BULK ):
                    return std::make_shared< IWG_Isotropic_Struc_Linear_Bulk >();

                case ( IWG_Type::STRUC_LINEAR_PRESSURE_BULK ):
                    return std::make_shared< IWG_Isotropic_Struc_Linear_Pressure_Bulk >();

                case ( IWG_Type::STRUC_LINEAR_DIRICHLET ):
                    return std::make_shared< IWG_Isotropic_Struc_Linear_Dirichlet >();

                case ( IWG_Type::STRUC_LINEAR_PRESSURE_DIRICHLET ):
                    return std::make_shared< IWG_Isotropic_Struc_Linear_Pressure_Dirichlet >();

                case ( IWG_Type::STRUC_LINEAR_INTERFACE ):
                    return std::make_shared< IWG_Isotropic_Struc_Linear_Interface >();

                case ( IWG_Type::STRUC_LINEAR_NEUMANN ):
                    return std::make_shared< IWG_Isotropic_Struc_Linear_Neumann >();

                default:
                    MORIS_ERROR( false, " IWG_Factory::create_IWGs - No IWG type specified. " );
                    return nullptr;
            }
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

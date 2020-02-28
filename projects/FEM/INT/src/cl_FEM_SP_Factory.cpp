#include "assert.hpp"
#include "cl_FEM_SP_Factory.hpp"            //FEM/INT/src
#include "cl_FEM_SP_Dirichlet_Nitsche.hpp"  //FEM/INT/src
#include "cl_FEM_SP_Ghost_Displacement.hpp" //FEM/INT/src
#include "cl_FEM_SP_Ghost_Virtual_Work.hpp" //FEM/INT/src
#include "cl_FEM_SP_Nitsche_Interface.hpp" //FEM/INT/src
#include "cl_FEM_SP_Master_Weight_Interface.hpp" //FEM/INT/src
#include "cl_FEM_SP_Slave_Weight_Interface.hpp" //FEM/INT/src
#include "cl_FEM_SP_Reciprocal_Total_Volume.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        std::shared_ptr< Stabilization_Parameter > SP_Factory::create_SP( fem::Stabilization_Type aStabilizationType )
        {
            switch( aStabilizationType )
            {
                case ( fem::Stabilization_Type::DIRICHLET_NITSCHE ):
                    return std::make_shared< SP_Dirichlet_Nitsche >();

                case ( fem::Stabilization_Type::GHOST_DISPL ):
                    return std::make_shared< SP_Ghost_Displacement >();

                case ( fem::Stabilization_Type::GHOST_VW ):
                    return std::make_shared< SP_Ghost_Virtual_Work >();

                case ( fem::Stabilization_Type::NITSCHE_INTERFACE ):
                    return std::make_shared< SP_Nitsche_Interface >();

                case ( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE ):
                    return std::make_shared< SP_Master_Weight_Interface >();

                case ( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE ):
                    return std::make_shared< SP_Slave_Weight_Interface >();

                case ( fem::Stabilization_Type::RECIPROCAL_TOTAL_VOLUME ):
                    return std::make_shared< SP_Reciprocal_Total_Volume >();

                default:
                    MORIS_ERROR( false, " SP_Factory::create_SP - No stabilization type specified. " );
                    return nullptr;
            }
        }
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

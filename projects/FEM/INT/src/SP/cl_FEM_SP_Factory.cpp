#include "assert.hpp"
//FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_SP_Dirichlet_Nitsche.hpp"
#include "cl_FEM_SP_Ghost_Displacement.hpp"
#include "cl_FEM_SP_Ghost_Virtual_Work.hpp"
#include "cl_FEM_SP_Nitsche_Interface.hpp"
#include "cl_FEM_SP_Reciprocal_Total_Volume.hpp"
#include "cl_FEM_SP_Incompressible_Flow.hpp"
#include "cl_FEM_SP_Viscous_Ghost.hpp"
#include "cl_FEM_SP_Convective_Ghost.hpp"
#include "cl_FEM_SP_Pressure_Ghost.hpp"
#include "cl_FEM_SP_Time_Velocity_Ghost.hpp"
#include "cl_FEM_SP_Velocity_Dirichlet_Nitsche.hpp"
#include "cl_FEM_SP_Velocity_SlipBoundary_Nitsche.hpp"
#include "cl_FEM_SP_Compressible_Velocity_Dirichlet_Nitsche.hpp"
#include "cl_FEM_SP_Compressible_Dirichlet_Nitsche.hpp"
#include "cl_FEM_SP_SUPG_Advection.hpp"
#include "cl_FEM_SP_YZBeta_Advection.hpp"
#include "cl_FEM_SP_GGLS_Diffusion.hpp"
#include "cl_FEM_SP_SUPG_Spalart_Allmaras_Turbulence.hpp"
#include "cl_FEM_SP_Turbulence_Dirichlet_Nitsche.hpp"
#include "cl_FEM_SP_Stab_Penalty_Contact.hpp"
#include "cl_FEM_SP_Penalty_Contact.hpp"
#include "cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.hpp"
#include "cl_FEM_SP_Measure.hpp"


namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        std::shared_ptr< Stabilization_Parameter > SP_Factory::create_SP(
                fem::Stabilization_Type aStabilizationType )
        {
            switch( aStabilizationType )
            {
                case fem::Stabilization_Type::DIRICHLET_NITSCHE :
                    return std::make_shared< SP_Dirichlet_Nitsche >();

                case fem::Stabilization_Type::GGLS_DIFFUSION :
                    return std::make_shared< SP_GGLS_Diffusion >();

                case fem::Stabilization_Type::GHOST_DISPL :
                    return std::make_shared< SP_Ghost_Displacement >();

                case fem::Stabilization_Type::GHOST_VW :
                    return std::make_shared< SP_Ghost_Virtual_Work >();

                case fem::Stabilization_Type::NITSCHE_INTERFACE :
                    return std::make_shared< SP_Nitsche_Interface >();

                case fem::Stabilization_Type::RECIPROCAL_TOTAL_VOLUME :
                    return std::make_shared< SP_Reciprocal_Total_Volume >();

                case fem::Stabilization_Type::INCOMPRESSIBLE_FLOW :
                    return std::make_shared< SP_Incompressible_Flow >();

                case fem::Stabilization_Type::VISCOUS_GHOST :
                    return std::make_shared< SP_Viscous_Ghost >();

                case fem::Stabilization_Type::CONVECTIVE_GHOST :
                    return std::make_shared< SP_Convective_Ghost >();

                case fem::Stabilization_Type::PRESSURE_GHOST :
                    return std::make_shared< SP_Pressure_Ghost >();

                case fem::Stabilization_Type::TIME_VELOCITY_GHOST :
                    return std::make_shared< SP_Time_Velocity_Ghost >();

                case fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE :
                    return std::make_shared< SP_Velocity_Dirichlet_Nitsche >();

                case fem::Stabilization_Type::VELOCITY_SLIPBOUNDARY_NITSCHE :
                    return std::make_shared< SP_Velocity_SlipBoundary_Nitsche >();

                case fem::Stabilization_Type::COMPRESSIBLE_VELOCITY_DIRICHLET_NITSCHE :
                    return std::make_shared< SP_Compressible_Velocity_Dirichlet_Nitsche >();

                case fem::Stabilization_Type::COMPRESSIBLE_DIRICHLET_NITSCHE :
                    return std::make_shared< SP_Compressible_Dirichlet_Nitsche >();

                case fem::Stabilization_Type::SUPG_ADVECTION :
                    return std::make_shared< SP_SUPG_Advection >();

                case fem::Stabilization_Type::YZBETA_ADVECTION :
                    return std::make_shared< SP_YZBeta_Advection >();

                case fem::Stabilization_Type::SUPG_SPALART_ALLMARAS_TURBULENCE :
                    return std::make_shared< SP_SUPG_Spalart_Allmaras_Turbulence >();

                case fem::Stabilization_Type::TURBULENCE_DIRICHLET_NITSCHE :
                    return std::make_shared< SP_Turbulence_Dirichlet_Nitsche >();

                case fem::Stabilization_Type::SPALART_ALLMARAS_NITSCHE_INTERFACE :
                    return std::make_shared< SP_Spalart_Allmaras_Nitsche_Interface >();

                case fem::Stabilization_Type::PENALTY_CONTACT :
                    return std::make_shared< SP_Penalty_Contact >();

                case fem::Stabilization_Type::STAB_PENALTY_CONTACT :
                    return std::make_shared< SP_Stab_Penalty_Contact >();

                case fem::Stabilization_Type::MEASURE :
                    return std::make_shared< SP_Measure >();

                default:
                    MORIS_ERROR( false, " SP_Factory::create_SP - No stabilization type specified. " );
                    return nullptr;
            }
        }
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

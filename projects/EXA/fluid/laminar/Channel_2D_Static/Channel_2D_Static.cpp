#include <string>
#include <iostream>
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

extern bool gInletVelocityBCFlag;
extern bool gInletPressureBCFlag;

#ifdef  __cplusplus
extern "C"
{
#endif
    //------------------------------------------------------------------------------
    namespace moris
    {

        bool sUseGhost          = true;

        std::string sFluid      = "HMR_dummy_c_p10,HMR_dummy_n_p10";
        std::string sFluidGhost = "ghost_p10";
        std::string sInlet      = "iside_b0_10_b1_8";
        std::string sWalls      = "iside_b0_10_b1_2,iside_b0_10_b1_14";

        moris::real tPlaneBottom = -0.5;                  /* y bottom plane (m) */
        moris::real tPlaneTop    =  0.5;                  /* y top plane    (m) */
        moris::real tPlaneLeft   = -2.5;                  /* x left plane   (m) */
        moris::real tPlaneRight  =  2.5;                  /* x right plane  (m) */
        moris::real tChannelRadius = ( tPlaneTop - tPlaneBottom ) / 2.0; /* channel radius  (m) */

        //Material Parameters
        moris::real tFluidDensity   = 1.0; /* Fluid density   () */
        moris::real tFluidViscosity = 1.0; /* Fluid viscosity () */

        // Boundary Conditions
        moris::real tInletPressure  = 20.0;  /* Inlet pressure  () */
        moris::real tGammaNitsche   = 50.0;  /* Penalty for Dirichlet BC */
        moris::real tGammaGPmu      = 0.05;  /* Penalty for ghost viscosity */
        moris::real tGammaGPu       = 0.05;  /* Penalty for ghost velocity */
        moris::real tGammaGPp       = 0.005; /* Penalty for ghost pressure */

        // Constant function for properties
        void Func_Const(
                moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                moris::fem::Field_Interpolator_Manager         * aFIManager )
        {
            aPropMatrix = aParameters( 0 );
        }

        // Inlet velocity function
        void Func_Inlet_U(
                moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                moris::fem::Field_Interpolator_Manager         * aFIManager )
        {
            // unpack parameters
            real tRadiusChannel = aParameters( 0 )( 0 );
            real tYChannel      = aParameters( 1 )( 0 );

            // get position in space
            real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

            // set size for aPropMatrix
            aPropMatrix.set_size( 2, 1, 0.0 );

            // velocity along x direction
            aPropMatrix( 0 ) = - ( tY - ( tYChannel + tRadiusChannel ) ) *
                    ( tY - ( tYChannel - tRadiusChannel ) ) /
                    ( 2.0 * std::pow( tRadiusChannel, 2.0 ) );
        }

        // Output criterion function
        bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
        {
            return true;
        }

        moris::real Func_Plane(
                const moris::Matrix< DDRMat >     & aCoordinates,
                const moris::Cell< moris::real* > & aGeometryParameters )
        {
            moris::real tXNormal = *( aGeometryParameters( 0 ) );
            moris::real tYNormal = *( aGeometryParameters( 1 ) );

            moris::real tXCenter = *( aGeometryParameters( 2 ) );
            moris::real tYCenter = *( aGeometryParameters( 3 ) );

            moris::real aReturnValue =
                    tXNormal * ( aCoordinates( 0 ) - tXCenter ) + tYNormal * ( aCoordinates( 1 ) - tYCenter );
            return aReturnValue;
        }

        void OPTParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 1 );
            tParameterlist( 0 ).resize( 1 );

            tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

            tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false);
        }

        void HMRParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 1 );
            tParameterlist( 0 ).resize( 1 );

            tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

            tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", "200,200");
            tParameterlist( 0 )( 0 ).set( "processor_decomposition_method",   1 );
            tParameterlist( 0 )( 0 ).set( "processor_dimensions",             "2,1");
            tParameterlist( 0 )( 0 ).set( "domain_dimensions",                "10.0,10.0");
            tParameterlist( 0 )( 0 ).set( "domain_offset",                    "-4.63,-4.74");
            tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  "1,2,3,4");
            tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           "0");

            tParameterlist( 0 )( 0 ).set( "lagrange_orders",  "1" );
            tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
            tParameterlist( 0 )( 0 ).set( "bspline_orders",   "1" );
            tParameterlist( 0 )( 0 ).set( "bspline_pattern",  "0" );

            tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0") ;

            tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
            tParameterlist( 0 )( 0 ).set( "refinement_buffer",  3 );
            tParameterlist( 0 )( 0 ).set( "staircase_buffer",   3 );
            tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );
            tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

            tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

            tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
            tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

            tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 0 );
        }

        void XTKParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 1 );
            tParameterlist( 0 ).resize( 1 );

            tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
            tParameterlist( 0 )( 0 ).set( "decompose",                 true );
            tParameterlist( 0 )( 0 ).set( "decomposition_type",        "conformal") ;
            tParameterlist( 0 )( 0 ).set( "enrich",                    true );
            tParameterlist( 0 )( 0 ).set( "basis_rank",                "bspline") ;
            tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices",       "0") ;
            tParameterlist( 0 )( 0 ).set( "ghost_stab",                true );
            tParameterlist( 0 )( 0 ).set( "multigrid",                 false );
            tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh",    false );
            tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", false );
        }

        void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 3 );
            tParameterlist( 0 ).resize( 1 );

            // Main GEN parameter list
            tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();

            // init geometry counter
            uint tGeoCounter = 0;

            // Bottom plane
            tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane");
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0.0,1.0,0.0," + std::to_string(tPlaneBottom)) ;
            tGeoCounter++;

            // Top plane
            tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane");
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0.0,1.0,0.0,"+ std::to_string(tPlaneTop));
            tGeoCounter++;

            // Left plane
            tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane");
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "1.0,0.0," + std::to_string(tPlaneLeft) + ",0.0");
            tGeoCounter++;

            // Right plane
            tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane");
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "1.0,0.0," + std::to_string(tPlaneRight) + ",0.0");
        }

        void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
        {
            // create a cell of cell of parameter list for fem
            tParameterList.resize( 8 );

            //------------------------------------------------------------------------------
            // fill the property part of the parameter list

            // init property counter
            uint tPropCounter = 0;

            // create viscosity property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropViscosity") ;
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::to_string(tFluidViscosity) );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
            tPropCounter++;

            // create density property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropDensity") ;
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::to_string(tFluidDensity) );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
            tPropCounter++;

            // create zero velocity property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropZeroU") ;
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "0.0;0.0") ;
            tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
            tPropCounter++;

            // create init velocity property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropInletU") ;
            tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Inlet_U") ;
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::to_string(tChannelRadius) + "/0.0" );
            tPropCounter++;

            // create pressure property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            "InletPressure") ;
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::to_string(tInletPressure) );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
            tPropCounter++;

            // create total pressure property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            "PropInletUpwind") ;
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      "1.0") ;
            tParameterList( 0 )( tPropCounter ).set( "value_function",           "Func_Const") ;
            tPropCounter++;

            //------------------------------------------------------------------------------
            // fill the constitutive model part of the parameter list

            // init CM counter
            uint tCMCounter = 0;

            // create fluid constitutive model
            tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
            tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMFluid") ;
            tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
            tParameterList( 1 )( tCMCounter ).set( "properties",        "PropViscosity,Viscosity;PropDensity,Density") ;
            tCMCounter++;

            //------------------------------------------------------------------------------
            // fill the stabilization parameter part of the parameter list

            // init SP counter
            uint tSPCounter = 0;

            // create SUPG stabilization parameter for Navier-Stokes
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPSUPGNS") ;
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters",     "36.0") ;
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       "PropViscosity,Viscosity;PropDensity,Density") ;
            tParameterList( 2 )( tSPCounter ).set( "master_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
            tSPCounter++;

            // create Nitsche stabilization parameter for velocity
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPNitscheU") ;
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::to_string(tGammaNitsche) +"/1.0" );
            tParameterList( 2 )( tSPCounter ).set( "master_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       "PropViscosity,Viscosity;PropDensity,Density") ;
            tSPCounter++;

            // create Ghost stabilization parameter for viscous velocity
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  "SPGPViscous") ;
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::VISCOUS_GHOST ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::to_string(tGammaGPmu) );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",   "PropViscosity,Viscosity") ;
            tSPCounter++;

            // create Ghost stabilization parameter for velocity
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGPVelocity") ;
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::CONVECTIVE_GHOST ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::to_string(tGammaGPu) );
            tParameterList( 2 )( tSPCounter ).set( "master_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       "PropDensity,Density") ;
            tSPCounter++;

            // create Ghost stabilization parameter for pressure
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      "SPGPPressure") ;
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::PRESSURE_GHOST ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::to_string(tGammaGPp) +"/1.0" );
            tParameterList( 2 )( tSPCounter ).set( "master_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       "PropViscosity,Viscosity;PropDensity,Density") ;
            tSPCounter++;

            //------------------------------------------------------------------------------
            // fill the IWG part of the parameter list

            // init IWG counter
            uint tIWGCounter = 0;

            // create incompressible NS velocity bulk IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGVelocityBulk") ;
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VX,VY") ;
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "VX,VY;P") ;
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,IncompressibleFluid") ;
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPSUPGNS,IncompressibleFlow") ;
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             sFluid );
            tIWGCounter++;

            // create incompressible NS pressure bulk IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGPressureBulk") ;
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "P") ;
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "VX,VY;P") ;
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,IncompressibleFluid") ;
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPSUPGNS,IncompressibleFlow") ;
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             sFluid );
            tIWGCounter++;

            if (gInletPressureBCFlag)
            {
                // create inlet total pressure BC
                tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
                tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInletImposedPressure") ;
                tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_IMPOSED_PRESSURE ) );
                tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VX,VY") ;
                tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "VX,VY;P") ;
                tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "InletPressure,Pressure") ;
                tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             sInlet );
                tIWGCounter++;
            }

            if (gInletVelocityBCFlag)
            {
                // create incompressible NS velocity Dirichlet IWG for inlet
                tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
                tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInletVelocity") ;
                tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) );
                tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VX,VY") ;
                tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "VX,VY;P") ;
                tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropInletU,Dirichlet") ;
                tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,IncompressibleFluid") ;
                tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPNitscheU,DirichletNitsche") ;
                tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             sInlet );
                tIWGCounter++;

                // create incompressible NS pressure Dirichlet IWG for inlet
                tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
                tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGInletPressure") ;
                tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) );
                tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "P") ;
                tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "VX,VY;P") ;
                tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropInletU,Dirichlet") ;
                tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,IncompressibleFluid") ;
                tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             sInlet );
                tIWGCounter++;
            }

            // create incompressible NS velocity Dirichlet IWG for walls
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGZeroVelocity") ;
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VX,VY") ;
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "VX,VY;P") ;
            tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropZeroU,Dirichlet") ;
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,IncompressibleFluid") ;
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPNitscheU,DirichletNitsche") ;
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             sWalls );
            tIWGCounter++;

            // create incompressible NS pressure Dirichlet IWG for walls
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGZeroPressure") ;
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "P") ;
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "VX,VY;P") ;
            tParameterList( 3 )( tIWGCounter ).set( "master_properties",          "PropZeroU,Dirichlet") ;
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", "CMFluid,IncompressibleFluid") ;
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             sWalls );
            tIWGCounter++;

            if( sUseGhost )
            {
                // create Ghost stabilization viscous IWG
                tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
                tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGGPViscous") ;
                tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
                tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VX,VY") ;
                tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "VX,VY;P") ;
                tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     "VX,VY;P") ;
                tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGPViscous,GhostSP") ;
                tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             sFluidGhost );
                tIWGCounter++;

                // create Ghost stabilization convective IWG
                tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
                tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGGPConvective") ;
                tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
                tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "VX,VY") ;
                tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "VX,VY;P") ;
                tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     "VX,VY;P") ;
                tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGPVelocity,GhostSP") ;
                tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             sFluidGhost );
                tIWGCounter++;

                // create Ghost stabilization pressure IWG
                tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
                tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   "IWGGPPressure") ;
                tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
                tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               "P") ;
                tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    "VX,VY;P") ;
                tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     "VX,VY;P") ;
                tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   "SPGPPressure,GhostSP") ;
                tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             sFluidGhost );
                tIWGCounter++;
            }

            //------------------------------------------------------------------------------
            // fill the IQI part of the parameter list

            // init IQI counter
            uint tIQICounter = 0;

            // create parameter list for IQI 1
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkVX") ;
            tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
            tParameterList( 4 )( tIQICounter ).set( "dof_quantity",               "VX,VY");
            tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    "VX,VY") ;
            tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
            tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             sFluid );
            tIQICounter++;

            // create parameter list for IQI 2
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkVY") ;
            tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
            tParameterList( 4 )( tIQICounter ).set( "dof_quantity",               "VX,VY");
            tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    "VX,VY") ;
            tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      1 );
            tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             sFluid );
            tIQICounter++;

            // create parameter list for IQI 3
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   "IQIBulkP") ;
            tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
            tParameterList( 4 )( tIQICounter ).set( "dof_quantity",               "P");
            tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    "P") ;
            tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
            tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             sFluid );
            tIQICounter++;

            //------------------------------------------------------------------------------
            // fill the computation part of the parameter list
            tParameterList( 5 ).resize( 1 );
            tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
        }

        void SOLParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 7 );
            for( uint Ik = 0; Ik < 7; Ik ++ )
            {
                tParameterlist( Ik ).resize( 1 );
            }

            tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

            tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

            tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
            tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop",    1e-06 );
            tParameterlist( 2 )( 0 ).set("NLA_relaxation_parameter", 1.0 );
            tParameterlist( 2 )( 0 ).set("NLA_max_iter",             50 );

            tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 0 ).set("NLA_DofTypes", "VX,VY;P") ;

            tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();

            tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
            tParameterlist( 5 )( 0 ).set("TSA_DofTypes",            "VX,VY;P") ;
            tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec",  "VX,0.5;VY,0.0;P,0.0") ;
            tParameterlist( 5 )( 0 ).set("TSA_Output_Indices",      "0") ;
            tParameterlist( 5 )( 0 ).set("TSA_Output_Crteria",      "Output_Criterion") ;

            tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
        }

        void MSIParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 1 );
            tParameterlist( 0 ).resize( 1 );

            tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        }

        void VISParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 1 );
            tParameterlist( 0 ).resize( 1 );

            tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
            tParameterlist( 0 )( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "./", "Channel_2D_Static.exo" ) );
            tParameterlist( 0 )( 0 ).set( "Mesh_Type"  , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
            tParameterlist( 0 )( 0 ).set( "Set_Names"  , sFluid );
            tParameterlist( 0 )( 0 ).set( "Field_Names",  "VX,VY,P" ) ;
            tParameterlist( 0 )( 0 ).set( "Field_Type" ,  "NODAL,NODAL,NODAL" ) ;
            tParameterlist( 0 )( 0 ).set( "IQI_Names"  ,  "IQIBulkVX,IQIBulkVY,IQIBulkP" ) ;
            tParameterlist( 0 )( 0 ).set( "Save_Frequency",1);
        }

        void MORISGENERALParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {

        }

        //------------------------------------------------------------------------------
    }

    //------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

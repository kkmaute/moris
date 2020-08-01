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
#include "cl_PRM_FEM_Parameters.hpp"
#include "cl_PRM_MSI_Parameters.hpp"
#include "cl_PRM_SOL_Parameters.hpp"
#include "cl_PRM_VIS_Parameters.hpp"
#include "cl_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "cl_PRM_XTK_Parameters.hpp"
#include "cl_PRM_OPT_Parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"


#ifdef  __cplusplus
extern "C"
{
#endif
    //------------------------------------------------------------------------------
    namespace moris
    {

        // Geometry Parameters
        moris::real tBottomPlane         = 0.0;  /* Bottom plane y (m) */
        moris::real tTopPlane            = 0.41; /* Top plane y (m) */
        moris::real tLeftPlane           = 0.0; /* Top plane y (m) */
        moris::real tRightPlane          = 1.845; /* Top plane y (m) */
        moris::real tMid                 = 0.205; /* Mid y (m) */
        moris::real tWallThickness       = 0.025; /* Wall thickness (m) */
        moris::real tWallBottomPlane     = tMid - tWallThickness; /* Mid y (m) */
        moris::real tWallTopPlane        = tMid + tWallThickness; /* Mid y (m) */
        moris::real tRadiusTopChannel    = ( tTopPlane - tWallTopPlane ) / 2;
        moris::real tRadiusBottomChannel = ( tWallBottomPlane - tBottomPlane ) / 2;
        moris::real tYcBottomChannel     = tBottomPlane + tRadiusBottomChannel;
        moris::real tYcTopChannel        = tTopPlane - tRadiusTopChannel;

        // Property constant function
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
            real tDir           = aParameters( 2 )( 0 );

            // get position in space
            real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

            // set size for aPropMatrix
            aPropMatrix.set_size( 2, 1, 0.0 );

            // velocity along x direction
            aPropMatrix( 0, 0 ) = - tDir * ( tY - ( tYChannel + tRadiusChannel ) ) * ( tY - ( tYChannel - tRadiusChannel ) ) / ( 2.0 * std::pow( tRadiusChannel, 2.0 ) );
        }

        // Output criterion function
        bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
        {
            return true;
        }

        // Plane geometry function
        moris::real Func_Plane(
                const moris::Matrix< DDRMat >     & aCoordinates,
                const moris::Cell< moris::real* > & aGeometryParameters )
        {
            moris::real tXNormal = *( aGeometryParameters( 0 ) );
            moris::real tXCenter = *( aGeometryParameters( 1 ) );
            moris::real tYNormal = *( aGeometryParameters( 2 ) );
            moris::real tYCenter = *( aGeometryParameters( 3 ) );

            moris::real aReturnValue =
                    tXNormal * ( aCoordinates( 0 ) - tXCenter ) +
                    tYNormal * ( aCoordinates( 1 ) - tYCenter );
            return aReturnValue;
        }

        // Dummy sensitivity function
        moris::Matrix<DDRMat> Func_Dummy(
                const moris::Matrix< DDRMat >     & aCoordinates,
                const moris::Cell< moris::real* > & aGeometryParameters )
        {
            moris::Matrix< DDRMat > dummy;
            return dummy;
        }

        // OPT parameter list
        void OPTParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 1 );
            tParameterlist( 0 ).resize( 1 );

            tParameterlist( 0 )( 0 ) = prm::create_opt_problem_parameter_list();

            tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false);
        }

        // HMR parameter list
        void HMRParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 1 );
            tParameterlist( 0 ).resize( 1 );

            tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

            tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", std::string("62,31"));
            tParameterlist( 0 )( 0 ).set( "domain_dimensions",                std::string("3,1"));
            tParameterlist( 0 )( 0 ).set( "domain_offset",                    std::string("-0.835,-0.186"));
            tParameterlist( 0 )( 0 ).set( "domain_sidesets",                  std::string("1,2,3,4"));
            tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes",           std::string("0"));

            tParameterlist( 0 )( 0 ).set( "lagrange_orders",  std::string("1" ));
            tParameterlist( 0 )( 0 ).set( "lagrange_pattern", std::string("0" ));
            tParameterlist( 0 )( 0 ).set( "bspline_orders",   std::string("1" ));
            tParameterlist( 0 )( 0 ).set( "bspline_pattern",  std::string("0" ));

            tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", std::string("0") );

            tParameterlist( 0 )( 0 ).set( "truncate_bsplines",  1 );
            tParameterlist( 0 )( 0 ).set( "refinement_buffer",  3 );
            tParameterlist( 0 )( 0 ).set( "staircase_buffer",   3 );
            tParameterlist( 0 )( 0 ).set( "initial_refinement", 0 );

            tParameterlist( 0 )( 0 ).set( "use_multigrid",  0 );
            tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

            tParameterlist( 0 )( 0 ).set( "adaptive_refinement_level", 0 );
        }

        // XTK parameter list
        void XTKParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 1 );
            tParameterlist( 0 ).resize( 1 );

            tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
            tParameterlist( 0 )( 0 ).set( "decompose",                 true );
            tParameterlist( 0 )( 0 ).set( "decomposition_type",        std::string("conformal") );
            tParameterlist( 0 )( 0 ).set( "enrich",                    true );
            tParameterlist( 0 )( 0 ).set( "basis_rank",                std::string("bspline") );
            tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices",       std::string("0") );
            tParameterlist( 0 )( 0 ).set( "ghost_stab",                true );
            tParameterlist( 0 )( 0 ).set( "multigrid",                 false );
            tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
            tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh",    true );
            tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", false );
        }

        // GEN parameter list
        void GENParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 3 );
            tParameterlist( 0 ).resize( 1 );

            tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();

            // init geometry counter
            uint tGeoCounter = 0;

            // Geometry parameter lists
            tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane");
            tParameterlist( 1 )( tGeoCounter ).set( "sensitivity_function_name", "Func_Dummy");
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0,0,1,0");
            tGeoCounter++;

            tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane");
            tParameterlist( 1 )( tGeoCounter ).set( "sensitivity_function_name", "Func_Dummy");
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0,0,1,0.41");
            tGeoCounter++;

            tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane");
            tParameterlist( 1 )( tGeoCounter ).set( "sensitivity_function_name", "Func_Dummy");
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0,0,1,0.18");
            tGeoCounter++;

            tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane");
            tParameterlist( 1 )( tGeoCounter ).set( "sensitivity_function_name", "Func_Dummy");
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "0,0,1,0.23");
            tGeoCounter++;

            tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane");
            tParameterlist( 1 )( tGeoCounter ).set( "sensitivity_function_name", "Func_Dummy");
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "1,0,0,0");
            tGeoCounter++;

            tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
            tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Func_Plane");
            tParameterlist( 1 )( tGeoCounter ).set( "sensitivity_function_name", "Func_Dummy");
            tParameterlist( 1 )( tGeoCounter ).set( "constant_parameters", "1,1.845,0,0");
            tGeoCounter++;
        }

        // FEM parameter list
        void FEMParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterList )
        {
            // create a cell of cell of parameter list for fem
            tParameterList.resize( 5 );

            //------------------------------------------------------------------------------
            // fill the property part of the parameter list

            // init property counter
            uint tPropCounter = 0;

            // fluid viscosity property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropFluidViscosity") );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("1.0") );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
            tPropCounter++;

            // fluid density property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropFluidDensity") );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("1.0") );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
            tPropCounter++;

            // fluid conductivity property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropFluidConductivity") );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("500.0") );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
            tPropCounter++;

            // fluid heat capacity property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropFluidHeatCapacity") );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("0.001") );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
            tPropCounter++;

            // solid density property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropSolidDensity") );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("2.0") );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
            tPropCounter++;

            // solid conductivity property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropSolidConductivity") );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("1.0") );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
            tPropCounter++;

            // solid heat capacity property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropSolidHeatCapacity") );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("0.001") );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
            tPropCounter++;

            // inlet velocity top channel property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropInletVelocityTop") );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string(
                    std::to_string(tRadiusTopChannel) + "/" +
                    std::to_string(tYcTopChannel) + "/-1.0" ) );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Inlet_U") );
            tPropCounter++;

            // inlet velocity bottom channel property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropInletVelocityBottom") );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string(
                    std::to_string(tRadiusBottomChannel) + "/" +
                    std::to_string(tYcBottomChannel) + "/1.0" ) );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Inlet_U") );
            tPropCounter++;

            // zero velocity property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropZeroVelocity") );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("0.0;0.0") );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
            tPropCounter++;

            // inlet temperature top channel property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropInletTempTop") );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("1.0") );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
            tPropCounter++;

            // inlet temperature bottom channel property
            tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
            tParameterList( 0 )( tPropCounter ).set( "property_name",            std::string("PropInletTempBottom") );
            tParameterList( 0 )( tPropCounter ).set( "function_parameters",      std::string("0.0") );
            tParameterList( 0 )( tPropCounter ).set( "value_function",           std::string("Func_Const") );
            tPropCounter++;

            //------------------------------------------------------------------------------
            // fill the constitutive model part of the parameter list

            // init CM counter
            moris::uint tCMCounter = 0;

            // fluid constitutive model
            tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
            tParameterList( 1 )( tCMCounter ).set( "constitutive_name", std::string("CMFluid") );
            tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
            tParameterList( 1 )( tCMCounter ).set( "properties",        std::string("PropFluidViscosity,Viscosity;PropFluidDensity,Density") );
            tCMCounter++;

            // fluid diffusion constitutive model
            tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
            tParameterList( 1 )( tCMCounter ).set( "constitutive_name", std::string("CMFluidDiffusion") );
            tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
            tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            tParameterList( 1 )( tCMCounter ).set( "properties",        std::string("PropFluidConductivity,Conductivity;PropFluidDensity,Density;PropFluidHeatCapacity,HeatCapacity") );
            tCMCounter++;

            // solid diffusion constitutive model
            tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
            tParameterList( 1 )( tCMCounter ).set( "constitutive_name", std::string("CMSolidDiffusion") );
            tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
            tParameterList( 1 )( tCMCounter ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            tParameterList( 1 )( tCMCounter ).set( "properties",        std::string("PropSolidConductivity,Conductivity;PropSolidDensity,Density;PropSolidHeatCapacity,HeatCapacity") );
            tCMCounter++;

            //------------------------------------------------------------------------------
            // fill the stabilization parameter part of the parameter list

            // init SP counter
            moris::uint tSPCounter = 0;

            // NS SUPG stabilization parameter
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPSUPGVelocity") );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("36.0") );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropFluidViscosity,Viscosity;PropFluidDensity,Density") );
            tParameterList( 2 )( tSPCounter ).set( "master_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
            tSPCounter++;

            // advection SUPG stabilization parameter
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPSUPGTemp") );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::SUPG_ADVECTION ) );
            tParameterList( 2 )( tSPCounter ).set( "master_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropFluidConductivity,Conductivity") );
            tSPCounter++;

            // nitsche for fluid velocity stabilization parameter
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPNitscheVelocity") );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("100.0/1.0") );
            tParameterList( 2 )( tSPCounter ).set( "master_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropFluidViscosity,Viscosity;PropFluidDensity,Density") );
            tSPCounter++;

            // nitsche for fluid temperature stabilization parameter 4
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPNitscheFluidTemp") );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("100.0") );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropFluidConductivity,Material") );
            tSPCounter++;

            // nitsche thermal interface stabilization parameter
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  std::string("SPInterfaceNitsche") );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string("100.0") );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",   std::string("PropFluidConductivity,Material") );
            tParameterList( 2 )( tSPCounter ).set( "slave_properties",    std::string("PropSolidConductivity,Material") );
            tSPCounter++;

            // master thermal interface stabilization parameter
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  std::string("SPInterfaceMasterWeight") );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE ) );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",   std::string("PropFluidConductivity,Material") );
            tParameterList( 2 )( tSPCounter ).set( "slave_properties",    std::string("PropSolidConductivity,Material") );
            tSPCounter++;

            // slave thermal interface stabilization parameter
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  std::string("SPInterfaceSlaveWeight") );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE ) );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",   std::string("PropFluidConductivity,Material") );
            tParameterList( 2 )( tSPCounter ).set( "slave_properties",    std::string("PropSolidConductivity,Material") );
            tSPCounter++;

            // ghost viscous stabilization parameter
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",  std::string("SPGPViscous") );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::VISCOUS_GHOST ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string("0.05") );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",   std::string("PropFluidViscosity,Viscosity") );
            tSPCounter++;

            // ghost convective stabilization parameter
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPGPVelocity") );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::CONVECTIVE_GHOST ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("0.05") );
            tParameterList( 2 )( tSPCounter ).set( "master_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropFluidDensity,Density") );
            tSPCounter++;

            // ghost fluid pressure stabilization parameter
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPGPPressure") );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::PRESSURE_GHOST ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("0.05/1.0") );
            tParameterList( 2 )( tSPCounter ).set( "master_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropFluidViscosity,Viscosity;PropFluidDensity,Density") );
            tSPCounter++;

            // ghost fluid temperature stabilization parameter
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPGPFluidTemp") );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("0.05") );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropFluidConductivity,Material") );
            tSPCounter++;

            // ghost solid temperature stabilization parameter
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name",      std::string("SPGPSolidTemp") );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type",      static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters",     std::string("0.05") );
            tParameterList( 2 )( tSPCounter ).set( "master_properties",       std::string("PropSolidConductivity,Material") );
            tSPCounter++;

            //------------------------------------------------------------------------------
            // fill the IWG part of the parameter list

            // init IWG counter
            moris::uint tIWGCounter = 0;

            // NS incompressible velocity bulk IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGVelocityBulk") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("VX,VY") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMFluid,IncompressibleFluid") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPSUPGVelocity,IncompressibleFlow") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46") );
            tIWGCounter++;

            // NS incompressible pressure bulk IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGPressureBulk") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("P") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMFluid,IncompressibleFluid") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPSUPGVelocity,IncompressibleFlow") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46") );
            tIWGCounter++;

            // fluid diffusion bulk IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGFluidDiffusionBulk") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMFluidDiffusion,Diffusion") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46") );
            tIWGCounter++;

            // solid diffusion bulk IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGSolidDiffusionBulk") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMSolidDiffusion,Diffusion") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p42,HMR_dummy_c_p42") );
            tIWGCounter++;

            // fluid advection bulk IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGAdvectionBulk") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::ADVECTION_BULK ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMFluidDiffusion,Diffusion") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPSUPGTemp,SUPG") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",           std::string("HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46") );
            tIWGCounter++;

            // top inlet fluid velocity IWG (velocity part)
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGInVelocityTop") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("VX,VY") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropInletVelocityTop,Dirichlet") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMFluid,IncompressibleFluid") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPNitscheVelocity,DirichletNitsche") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("iside_b0_46_b1_47") );
            tIWGCounter++;

            // top inlet fluid velocity IWG (pressure part)
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGInPressureTop") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("P") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropInletVelocityTop,Dirichlet") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMFluid,IncompressibleFluid") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("iside_b0_46_b1_47") );
            tIWGCounter++;

            // bottom inlet fluid velocity IWG (velocity part)
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGInVelocityBottom") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("VX,VY") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropInletVelocityBottom,Dirichlet") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMFluid,IncompressibleFluid") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPNitscheVelocity,DirichletNitsche") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("iside_b0_34_b1_32") );
            tIWGCounter++;

            // bottom inlet fluid velocity IWG (pressure part)
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGInPressureBottom") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("P") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropInletVelocityBottom,Dirichlet") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMFluid,IncompressibleFluid") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("iside_b0_34_b1_32") );
            tIWGCounter++;

            // wall fluid velocity IWG (velocity part)
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGWallVelocity") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("VX,VY") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropZeroVelocity,Dirichlet") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMFluid,IncompressibleFluid") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPNitscheVelocity,DirichletNitsche") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("iside_b0_34_b1_2,iside_b0_34_b1_42,iside_b0_46_b1_42,iside_b0_46_b1_62") );
            tIWGCounter++;

            // wall fluid velocity IWG (pressure part)
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGWallPressure") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("P") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropZeroVelocity,Dirichlet") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMFluid,IncompressibleFluid") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("iside_b0_34_b1_2,iside_b0_34_b1_42,iside_b0_46_b1_42,iside_b0_46_b1_62") );
            tIWGCounter++;

            // top inlet fluid temperature IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGInletTempTop") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropInletTempTop,Dirichlet") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMFluidDiffusion,Diffusion") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPNitscheFluidTemp,DirichletNitsche") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("iside_b0_46_b1_47") );
            tIWGCounter++;

            // bottom inlet fluid temperature IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGInletTempBottom") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_properties",          std::string("PropInletTempBottom,Dirichlet") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models", std::string("CMFluidDiffusion,Diffusion") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPNitscheFluidTemp,DirichletNitsche") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("iside_b0_34_b1_32") );
            tIWGCounter++;

            // fluid/solid thermal interface IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                        std::string("IWGInterfaceFluidSolid") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                        static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",                    std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",         std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",          std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_constitutive_models",      std::string("CMFluidDiffusion,Diffusion") );
            tParameterList( 3 )( tIWGCounter ).set( "slave_constitutive_models",       std::string("CMSolidDiffusion,Diffusion") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",        std::string("SPInterfaceNitsche,NitscheInterface;SPInterfaceMasterWeight,MasterWeightInterface;SPInterfaceSlaveWeight,SlaveWeightInterface") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",                  std::string("dbl_iside_p0_34_p1_42,dbl_iside_p0_46_p1_42") );
            tIWGCounter++;

            // ghost viscous IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGGPViscous") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("VX,VY") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPGPViscous,GhostSP") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("ghost_p34,ghost_p46") );
            tIWGCounter++;

            // ghost convective IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGGPConvective") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("VX,VY") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPGPVelocity,GhostSP") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("ghost_p34,ghost_p46") );
            tIWGCounter++;

            // ghost pressure IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGGPPressure") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("P") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPGPPressure,GhostSP") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("ghost_p34,ghost_p46") );
            tIWGCounter++;

            // ghost fluid temperature IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGGPFluidTemp") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPGPFluidTemp,GhostSP") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("ghost_p34,ghost_p46") );
            tIWGCounter++;

            // ghost solid temperature IWG
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name",                   std::string("IWGGPSolidTemp") );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual",               std::string("TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "master_dof_dependencies",    std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "slave_dof_dependencies",     std::string("VX,VY;P;TEMP") );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters",   std::string("SPGPSolidTemp,GhostSP") );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names",             std::string("ghost_p42") );
            tIWGCounter++;

            //------------------------------------------------------------------------------
            // fill the IQI part of the parameter list

            // init IQI counter
            moris::uint tIQICounter = 0;

            // VX IQI
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIBulkVX") );
            tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
            tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::VX ) );
            tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("VX,VY") );
            tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
            tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42") );
            tIQICounter++;

            // VY IQI
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIBulkVY") );
            tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
            tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::VY ) );
            tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("VX,VY") );
            tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      1 );
            tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42") );
            tIQICounter++;

            // P IQI
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIBulkP") );
            tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
            tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::P ) );
            tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("P") );
            tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
            tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42") );
            tIQICounter++;

            // TEMP IQI
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( tIQICounter ).set( "IQI_name",                   std::string("IQIBulkTEMP") );
            tParameterList( 4 )( tIQICounter ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::DOF ) );
            tParameterList( 4 )( tIQICounter ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::TEMP ) );
            tParameterList( 4 )( tIQICounter ).set( "master_dof_dependencies",    std::string("TEMP") );
            tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index",      0 );
            tParameterList( 4 )( tIQICounter ).set( "mesh_set_names",             std::string("HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42") );
            tIQICounter++;
        }

        // SOL parameter list
        void SOLParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 7 );
            for( uint Ik = 0; Ik < 7; Ik ++)
            {
                tParameterlist( Ik ).resize( 1 );
            }

            tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

            tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

            tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
            tParameterlist( 2 )( 0 ).set("NLA_rel_res_norm_drop", 1e-07 );

            tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 0 ).set("NLA_DofTypes"      , std::string("VX,VY;P;TEMP") );

            tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
            tParameterlist( 4 )( 0 ).set("TSA_Num_Time_Steps",20 );
            tParameterlist( 4 )( 0 ).set("TSA_Time_Frame",    0.4 );

            tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
            tParameterlist( 5 )( 0 ).set("TSA_DofTypes"       , std::string("VX,VY;P;TEMP") );
            tParameterlist( 5 )( 0 ).set("TSA_Initialize_Sol_Vec" , std::string("VX,1E-4;VY,1E-4;P,0.0;TEMP,0.0") );
            tParameterlist( 5 )( 0 ).set("TSA_Output_Indices" , std::string("0") );
            tParameterlist( 5 )( 0 ).set("TSA_Output_Crteria" , std::string("Output_Criterion") );
            tParameterlist( 5 )( 0 ).set("TSA_time_level_per_type", std::string("VX,2;VY,2;P,2;TEMP,2") );

            tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
        }

        // MSI parametr list
        void MSIParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 1 );
            tParameterlist( 0 ).resize( 1 );

            tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        }

        // VIS parameter list
        void VISParameterList( moris::Cell< moris::Cell< ParameterList > > & tParameterlist )
        {
            tParameterlist.resize( 1 );
            tParameterlist( 0 ).resize( 1 );

            tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
            tParameterlist( 0 )( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "./", "2_Channels_Temp.exo" ) );
            tParameterlist( 0 )( 0 ).set( "Mesh_Type"  , static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
            tParameterlist( 0 )( 0 ).set( "Set_Names"  , std::string( "HMR_dummy_n_p34,HMR_dummy_c_p34,HMR_dummy_n_p46,HMR_dummy_c_p46,HMR_dummy_n_p42,HMR_dummy_c_p42" ) );
            tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "VX,VY,P,TEMP" ) );
            tParameterlist( 0 )( 0 ).set( "Field_Type" , std::string( "NODAL,NODAL,NODAL,NODAL" ) );
            tParameterlist( 0 )( 0 ).set( "Output_Type", std::string( "VX,VY,P,TEMP" ) );
        }

        //------------------------------------------------------------------------------
    }

    //------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

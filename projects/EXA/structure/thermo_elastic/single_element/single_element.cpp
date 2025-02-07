/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * single_element.cpp
 *
 */

#include <string>
#include <iostream>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "parameters.hpp"
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"
#include "fn_stringify_matrix.hpp"

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    /* ------------------------------------------------------------------------ */
    // General

    // file name
    std::string tName = "single_element";

    // activate thermal/structural portion
    bool tHaveThermo    = true;
    bool tHaveStruct    = true;
    bool tHaveStaggered = true;

    bool tHaveMixedSolverConfig = true;

    bool tUseBelos = false;

    bool tUseGGLS = true;

    // Pressure load
    std::string tPressure = "100000.0";

    // CTE
    std::string tThermalExpansion = "0.1";

    /* ------------------------------------------------------------------------ */
    // Solver config

    int         tNLA_max_iter             = 3;
    moris::real tNLA_rel_res_norm_drop    = 1.0e-10;
    moris::real tNLA_relaxation_parameter = 1.0;

    int         tTSA_Num_Time_Steps = 3;
    moris::real tTSA_Time_Frame     = 0.3;

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    // Interpolation order
    std::string tOrder = "1";

    Vector< uint > tNumElemsPerDim = { 1, 1 };
    Vector< real > tDomainDims     = { 0.01, 0.05 };

    /* ------------------------------------------------------------------------ */
    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    // Bulk Phases
    std::string tBulk = "HMR_dummy_n_p0";

    // boundaries
    std::string tBottom = "SideSet_1_n_p0";
    std::string tRight  = "SideSet_2_n_p0";
    std::string tTop    = "SideSet_3_n_p0";
    std::string tLeft   = "SideSet_4_n_p0";

    /* ------------------------------------------------------------------------ */
    // material parameters, kg is scaled with a factor 1e-6

    // Copper material properties
    std::string tDensity       = "8960.0";
    std::string tHeatCapacity  = "1460.0";
    std::string tConductivity  = "370.0";
    std::string tYoungsModulus = "120e9";
    std::string tPoissonRatio  = "0.34";

    /* ------------------------------------------------------------------------ */
    // boundary conditions

    // bedding to suppress RBM
    std::string tBedding = std::to_string( 1.0 * 1.0e-5 );

    // Initial Temperature
    moris::real tInitialTemp   = 1.0;    // in deg C
    moris::real tInitialStruct = 0.0;
    std::string tReferenceTemp = "0.0";

    // Heat flux
    std::string tHeatLoad = "1.0e4";

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName  = tName + ".exo";
    std::string tLibraryName     = tName + ".so";
    std::string tHDF5Path        = tName + ".hdf5";
    std::string tGENOutputFile   = "GEN_" + tName + ".exo";
    bool        tOutputCriterion = true;

    //------------------------------------------------------------------------------
    //-------------------------------- FUNCTIONS -----------------------------------
    //------------------------------------------------------------------------------

    /* ------------------------------------------------------------------------ */
    // GEOMETRY (LEVEL-SET) FUNCTIONS
    /* ------------------------------------------------------------------------ */

    //-----------------------------------------------------------------------------

    // Back Wall
    moris::real Back_Wall(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        real tVal = aCoordinates( 0 ) - 100.0;

        // clean return value to return non-zero value
        return std::abs( tVal ) < 1.0e-8 ? 1.0e-8 : tVal;
    }

    /* ------------------------------------------------------------------------ */
    // PROPERTY FUNCTIONS (incl. INITIAL & BOUNDARY CONDITIONS)
    /* ------------------------------------------------------------------------ */

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    // initial temperature
    void
    Func_Initial_Condition(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = { { tInitialTemp } };
    }

    // initial structure
    void
    Func_Initial_Condition_Struct(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = { { tInitialStruct }, { tInitialStruct } };
    }

    void
    Func_Neumann_U(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // aPropMatrix = { { aParameters( 0 )( 0 ) },{ 0.0 } };
        aPropMatrix = { { 0.0 }, { aParameters( 0 )( 0 ) } };
    }

    /* ------------------------------------------------------------------------ */
    // DUMMY FUNCTIONS
    /* ------------------------------------------------------------------------ */

    // Output criterion for VIS mesh
    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return tOutputCriterion;
    }

    // Dummy function for unused sensitivities if needed
    moris::Matrix< DDRMat > Func_Dummy_Sensitivity(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&     aGeometryParameters )
    {
        moris::Matrix< DDRMat > aReturnValue = { { 0.0 } };
        return aReturnValue;
    }

    /* ------------------------------------------------------------------------ */
    // PARAMETER LISTS
    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists.set( "domain_dimensions", tDomainDims );
        aParameterLists.set( "lagrange_output_meshes", std::string( "0" ) );

        aParameterLists.set( "lagrange_orders", tOrder );
        aParameterLists.set( "lagrange_pattern", std::string( "0" ) );
        aParameterLists.set( "bspline_orders", tOrder );
        aParameterLists.set( "bspline_pattern", std::string( "0" ) );


        aParameterLists.set( "use_number_aura", 1 );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", std::string( "conformal" ) );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", std::string( "bspline" ) );
        aParameterLists.set( "enrich_mesh_indices", std::string( "0" ) );
        aParameterLists.set( "ghost_stab", false );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "output_mesh_file", tGENOutputFile );

        // Dummy Geometry
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Back_Wall" );
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - STRUCTURE (ni-w-alloy?)
        //------------------------------------------------------------------------------

        // Density Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", tDensity );
        aParameterLists.set( "value_function", "Func_Const" );

        // Heat Capacity Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropHeatCapacity" );
        aParameterLists.set( "function_parameters", tHeatCapacity );
        aParameterLists.set( "value_function", "Func_Const" );

        // Conductivity Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", tConductivity );
        aParameterLists.set( "value_function", "Func_Const" );

        // Youngs Modulus Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungsModulus" );
        aParameterLists.set( "function_parameters", tYoungsModulus );
        aParameterLists.set( "value_function", "Func_Const" );

        // Poisson Ratio Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropPoissonRatio" ) );
        aParameterLists.set( "function_parameters", tPoissonRatio );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        // CTE for Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropThermalExpansion" ) );
        aParameterLists.set( "function_parameters", tThermalExpansion );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        //------------------------------------------------------------------------------
        // OTHER MATERIAL PARAMETERS
        //------------------------------------------------------------------------------

        // properties for bedding (supression for RBMs, both Shell and PCM)
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropBedding" );
        aParameterLists.set( "function_parameters", tBedding );
        aParameterLists.set( "value_function", "Func_Const" );

        // Phase Change Temperature of PCM
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropPCTemp_Dummy" ) );
        aParameterLists.set( "function_parameters", "1.0e12" );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        // Phase Change Temperature Range of PCM
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropPCconst_Dummy" ) );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        // Cubic Phase State Function
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropPhaseStateFnct" ) );
        aParameterLists.set( "function_parameters", std::string( "2.0" ) );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        // Dummy latent heat for non-pc material
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropLatentHeat_Dummy" ) );
        aParameterLists.set( "function_parameters", std::string( "0.0" ) );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropPCTemp_Dummy" ) );
        aParameterLists.set( "function_parameters", std::string( "10000.0" ) );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        //------------------------------------------------------------------------------
        // BOUNDARY CONDITIONS
        //------------------------------------------------------------------------------

        // heat flux from outside
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropImposedFlux" ) );
        aParameterLists.set( "function_parameters", tHeatLoad );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        // reference temperature for thermal expansion
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropReferenceTemp" ) );
        aParameterLists.set( "function_parameters", tReferenceTemp );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        // pressure load
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropTraction" ) );
        aParameterLists.set( "function_parameters", tPressure );
        aParameterLists.set( "value_function", std::string( "Func_Neumann_U" ) );

        // pressure load
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropPressure" ) );
        aParameterLists.set( "function_parameters", tPressure );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        // Dirichlet structure
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropDirichletStruct" ) );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        // time continuity weights
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropWeightCurrent" ) );
        aParameterLists.set( "function_parameters", std::string( "100.0" ) );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropWeightPrevious" ) );
        aParameterLists.set( "function_parameters", std::string( "100.0" ) );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        // Initial Temperature
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropInitialCondition" ) );
        aParameterLists.set( "value_function", std::string( "Func_Initial_Condition" ) );

        // Initial Structure
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropInitialConditionStruct" ) );
        aParameterLists.set( "value_function", std::string( "Func_Initial_Condition_Struct" ) );

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // DIFFUSION
        //------------------------------------------------------------------------------

        if ( tHaveThermo )
        {
            // create parameter list for constitutive model - shell - 1
            aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
            aParameterLists.set( "constitutive_name", std::string( "CMDiffusion" ) );
            aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
            aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            aParameterLists.set( "properties", std::string( "PropConductivity , Conductivity;" ) + std::string( "PropDensity      , Density;" ) + std::string( "PropHeatCapacity , HeatCapacity" ) );
            }

        //------------------------------------------------------------------------------
        // LINEAR ELASTICITY
        //------------------------------------------------------------------------------

        if ( tHaveStruct )
        {
            // linear elasticity - shell - 1
            aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
            aParameterLists.set( "constitutive_name", "CMStrucLinIso" );
            aParameterLists.set( "model_type",  fem::Model_Type::PLANE_STRESS ) ;
            aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;

            if ( tHaveThermo )
            {
                aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
                aParameterLists.set( "properties", std::string( "PropYoungsModulus,    YoungsModulus;" ) + std::string( "PropPoissonRatio,     PoissonRatio;" ) + std::string( "PropThermalExpansion, CTE;" ) + std::string( "PropReferenceTemp,    ReferenceTemperature" ) );
            }
            else
            {
                aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
                aParameterLists.set( "properties", std::string( "PropYoungsModulus,  YoungsModulus;" ) + std::string( "PropPoissonRatio,   PoissonRatio" ) );
            }
            }

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // GGLS
        //------------------------------------------------------------------------------

        if ( tHaveThermo )
        {
            // create parameter list for GGLS stabilization parameter for Skin
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", std::string( "SPGGLSDiffusion" ) );
            aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GGLS_DIFFUSION ) ;
            aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            aParameterLists.set( "leader_properties", std::string( "PropConductivity      , Conductivity;" ) + std::string( "PropDensity           , Density;" ) + std::string( "PropHeatCapacity      , HeatCapacity;" ) + std::string( "PropLatentHeat_Dummy  , LatentHeat;" ) + std::string( "PropPCTemp_Dummy      , PCTemp;" ) + std::string( "PropPhaseStateFnct    , PhaseStateFunction;" ) + std::string( "PropPCconst_Dummy     , PhaseChangeConst" ) );
            }

        //------------------------------------------------------------------------------
        // NITSCHE DIRICHLET
        //------------------------------------------------------------------------------

        if ( tHaveStruct )
        {
            // Displacements - Shell - back wall
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPNitscheStruc" );
            aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
            aParameterLists.set( "function_parameters", "100.0" );
            aParameterLists.set( "leader_properties", "PropYoungsModulus,Material" );
            }

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        //------------------------------------------------------------------------------
        // BULK IWGs
        //------------------------------------------------------------------------------

        // diffusion - Shell
        if ( tHaveThermo )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGDiffusionBulk" ) );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
            aParameterLists.set( "dof_residual", std::string( "TEMP" ) );
            aParameterLists.set( "leader_dof_dependencies", std::string( "TEMP" ) );
            aParameterLists.set( "leader_constitutive_models", std::string( "CMDiffusion,Diffusion" ) );
            aParameterLists.set( "mesh_set_names", tBulk );

            if ( tUseGGLS )
            {
                aParameterLists.set( "stabilization_parameters", std::string( "SPGGLSDiffusion,GGLSParam" ) );
            }
        }

        // linear elasticity - Shell
        if ( tHaveStruct )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGStructShell" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
            aParameterLists.set( "dof_residual", "UX,UY" );

            if ( tHaveThermo )
            {
                aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
            }
            else
            {
                aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            }

            aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
            // aParameterLists.set( "leader_properties",          "PropBedding,Bedding" );
            aParameterLists.set( "mesh_set_names", tBulk );
            }

        //------------------------------------------------------------------------------
        // NEUMANN BCs - IWGs
        //------------------------------------------------------------------------------

        // heat flux on outside of Shell
        if ( tHaveThermo )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGInletFlux" ) );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
            aParameterLists.set( "dof_residual", std::string( "TEMP" ) );

            if ( tHaveStruct )
            {
                aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
            }
            else
            {
                aParameterLists.set( "leader_dof_dependencies", "TEMP" );
            }

            aParameterLists.set( "leader_properties", std::string( "PropImposedFlux,Neumann" ) );
            aParameterLists.set( "mesh_set_names", tTop );
            }

        // pressure pulling on outside of Shell
        if ( tHaveStruct )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGNeumannPressure" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
            aParameterLists.set( "dof_residual", "UX,UY" );

            if ( tHaveThermo )
            {
                aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
            }
            else
            {
                aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            }

            // aParameterLists.set( "leader_properties",          "PropTraction,Traction");
            aParameterLists.set( "leader_properties", "PropPressure,Pressure" );
            aParameterLists.set( "mesh_set_names", tLeft );
            }

        //------------------------------------------------------------------------------
        // DIRICHLET BCS - IWGs
        //------------------------------------------------------------------------------

        // displacements - shell - back wall
        if ( tHaveStruct )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGDirichletStruct" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
            aParameterLists.set( "dof_residual", "UX,UY" );

            if ( tHaveThermo )
            {
                aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
            }
            else
            {
                aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            }

            aParameterLists.set( "leader_properties", "PropDirichletStruct,Dirichlet" );
            aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
            aParameterLists.set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
            aParameterLists.set( "mesh_set_names", tBottom );
            }

        //------------------------------------------------------------------------------
        // IWGs - TIME CONTINUITY
        //------------------------------------------------------------------------------

        // Time continuity
        if ( tHaveThermo )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGTimeContinuityTemp" ) );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
            aParameterLists.set( "dof_residual", std::string( "TEMP" ) );

            if ( tHaveStruct )
            {
                aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
            }
            else
            {
                aParameterLists.set( "leader_dof_dependencies", "TEMP" );
            }

            aParameterLists.set( "leader_properties", std::string( "PropWeightCurrent,WeightCurrent;" ) + std::string( "PropWeightPrevious,WeightPrevious;" ) + std::string( "PropInitialCondition,InitialCondition" ) );
            aParameterLists.set( "mesh_set_names", tBulk );
            aParameterLists.set( "time_continuity", true );
            }

        // Time continuity Structure
        //         if ( tHaveStruct )
        //         {
        //             aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        //             aParameterLists.set( "IWG_name",                   std::string("IWGTimeContinuityStruct") );
        //             aParameterLists.set( "IWG_type",                    fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        //             aParameterLists.set( "dof_residual",               std::string("UX,UY") );
        //
        //             if ( tHaveThermo ) {
        //                 aParameterLists.set( "leader_dof_dependencies",    "UX,UY;TEMP");
        //             } else {
        //                 aParameterLists.set( "leader_dof_dependencies",    "UX,UY");
        //             }
        //
        //             aParameterLists.set( "leader_properties", std::string("PropWeightCurrent,WeightCurrent;")      +
        //                                                                          std::string("PropWeightPrevious,WeightPrevious;")    +
        //                                                                          std::string("PropInitialConditionStruct,InitialCondition") );
        //             aParameterLists.set( "mesh_set_names",             tBulk );
        //             aParameterLists.set( "time_continuity",            true );
        //             //         }

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        if ( tHaveThermo )
        {
            // Nodal Temperature IQI
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
            aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
            aParameterLists.set( "dof_quantity", "TEMP" );
            aParameterLists.set( "leader_dof_dependencies", std::string( "TEMP" ) );
            aParameterLists.set( "vectorial_field_index", 0 );
            aParameterLists.set( "mesh_set_names", tBulk );
            }

        // Volume IQI - TotalDomain - use once to find total volume to compute max dof
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQITotalVolume" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::VOLUME ) ;

        if ( tHaveThermo )
        {
            aParameterLists.set( "leader_dof_dependencies", "TEMP" );
        }
        else if ( tHaveStruct )
        {
            aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        }
        else if ( tHaveThermo && tHaveStruct )
        {
            aParameterLists.set( "leader_dof_dependencies", "UX,UY;TEMP" );
        }

        aParameterLists.set( "mesh_set_names", tBulk );

        if ( tHaveStruct )
        {
            // X-displacement
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "IQIBulkDISPX" );
            aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
            aParameterLists.set( "dof_quantity", "UX,UY" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists.set( "vectorial_field_index", 0 );
            aParameterLists.set( "mesh_set_names", tBulk );

            // Y-displacement
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "IQIBulkDISPY" );
            aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
            aParameterLists.set( "dof_quantity", "UX,UY" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists.set( "vectorial_field_index", 1 );
            aParameterLists.set( "mesh_set_names", tBulk );
            }

        // create computation parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        // Thermo Only solver config
        if ( !tHaveStruct && tHaveThermo )
        {
            std::cout << "Thermo ONLY !!! \n"
                      << std::flush;

            if ( tUseBelos )
            {
                aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::BELOS_IMPL );
            }
            else
            {
                aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
            }

            aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

            aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
            aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
            aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
            aParameterLists.set( "NLA_max_iter", tNLA_max_iter );
            aParameterLists.set( "NLA_combined_res_jac_assembly", true );

            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_DofTypes", "TEMP" );

            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_DofTypes", "TEMP" );
            aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );

            aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
            aParameterLists.set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
            aParameterLists.set( "TSA_Time_Frame", tTSA_Time_Frame );

            aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
            aParameterLists.set( "TSA_DofTypes", "TEMP" );
            aParameterLists.set( "TSA_Initialize_Sol_Vec", "TEMP,1.0" );
            aParameterLists.set( "TSA_Output_Indices", "0" );
            aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
            aParameterLists.set( "TSA_time_level_per_type", "TEMP,2" );
        }    // end: thermo only

        // Struct Only solver config
        else if ( tHaveStruct && !tHaveThermo )
        {
            std::cout << "Struct ONLY !!! \n"
                      << std::flush;

            if ( tUseBelos )
                aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::BELOS_IMPL );
            else
                aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

            aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

            aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
            aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
            aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
            aParameterLists.set( "NLA_max_iter", tNLA_max_iter );
            aParameterLists.set( "NLA_combined_res_jac_assembly", true );

            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_DofTypes", "UX,UY" );

            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_DofTypes", "UX,UY" );
            aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );

            aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
            aParameterLists.set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
            aParameterLists.set( "TSA_Time_Frame", tTSA_Time_Frame );

            aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
            aParameterLists.set( "TSA_DofTypes", "UX,UY" );
            aParameterLists.set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0" );
            aParameterLists.set( "TSA_Output_Indices", "0" );
            aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
            aParameterLists.set( "TSA_time_level_per_type", "UX,2;UY,2" );
        }    // end: struct only

        else if ( tHaveStruct && tHaveThermo && !tHaveStaggered )
        {
            std::cout << "NON-Staggered SOLVE !!! \n"
                      << std::flush;

            if ( tUseBelos )
            {
                aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::BELOS_IMPL );
            }
            else
            {
                aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
            }

            aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

            aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
            aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
            aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
            aParameterLists.set( "NLA_max_iter", tNLA_max_iter );
            aParameterLists.set( "NLA_combined_res_jac_assembly", true );

            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_DofTypes", "UX,UY;TEMP" );

            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_DofTypes", "UX,UY;TEMP" );
            aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );

            aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
            aParameterLists.set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
            aParameterLists.set( "TSA_Time_Frame", tTSA_Time_Frame );

            aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
            aParameterLists.set( "TSA_DofTypes", "UX,UY;TEMP" );
            aParameterLists.set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0;TEMP,1.0" );
            aParameterLists.set( "TSA_Output_Indices", "0" );
            aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
            if ( tHaveMixedSolverConfig )
                aParameterLists.set( "TSA_time_level_per_type", "UX,1;UY,1;TEMP,2" );
            else
                aParameterLists.set( "TSA_time_level_per_type", "UX,2;UY,2;TEMP,2" );
        }    // end: coupled & NOT staggered

        // Coupled Thermo Structural with Staggered solver config
        else if ( tHaveStruct && tHaveThermo && tHaveStaggered )
        {
            std::cout << "Staggered SOLVE !!! \n"
                      << std::flush;

            if ( tUseBelos )
            {
                aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::BELOS_IMPL );
            }
            else
            {
                aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
            }

            aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

            // ----------------------------------------------------------

            aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
            aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
            aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
            aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
            aParameterLists.set( "NLA_max_iter", tNLA_max_iter );
            aParameterLists.set( "NLA_combined_res_jac_assembly", true );

            aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
            aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;

            // ----------------------------------------------------------

            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
            aParameterLists.set( "NLA_DofTypes", "UX,UY" );

            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
            aParameterLists.set( "NLA_DofTypes", "TEMP" );

            aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
            aParameterLists.set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;
            aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "1,0" );
            aParameterLists.set( "NLA_DofTypes", "UX,UY;TEMP" );
            aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );

            // ----------------------------------------------------------

            aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
            aParameterLists.set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
            aParameterLists.set( "TSA_Time_Frame", tTSA_Time_Frame );
            aParameterLists.set( "TSA_Nonlinear_Solver", 2 );

            // ----------------------------------------------------------

            aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
            aParameterLists.set( "TSA_DofTypes", "UX,UY;TEMP" );
            aParameterLists.set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0;TEMP,1.0" );
            aParameterLists.set( "TSA_Output_Indices", "0" );
            aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
            if ( tHaveMixedSolverConfig )
                aParameterLists.set( "TSA_time_level_per_type", "UX,1;UY,1;TEMP,2" );
            else
                aParameterLists.set( "TSA_time_level_per_type", "UX,2;UY,2;TEMP,2" );
        }    // end: coupled & staggered

        else
        {
            std::cout << "Input file: incomplete physics configuration. Check what physics are used. \n"
                      << std::flush;
        }

        aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_save_operator_to_matlab", "jacobian.dat" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", tBulk );

        if ( tHaveStruct and tHaveThermo )
        {
            aParameterLists.set( "Field_Names", "TEMP,UX,UY,VOLUME" );
            aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,GLOBAL" );
            aParameterLists.set( "IQI_Names", "IQIBulkTEMP,IQIBulkDISPX,IQIBulkDISPY,IQITotalVolume" );
        }
        else if ( tHaveStruct )
        {
            aParameterLists.set( "Field_Names", "UX,UY,VOLUME" );
            aParameterLists.set( "Field_Type", "NODAL,NODAL,GLOBAL" );
            aParameterLists.set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQITotalVolume" );
        }
        else if ( tHaveThermo )
        {
            aParameterLists.set( "Field_Names", "TEMP,VOLUME" );
            aParameterLists.set( "Field_Type", "NODAL,GLOBAL" );
            aParameterLists.set( "IQI_Names", "IQIBulkTEMP,IQITotalVolume" );
        }

        aParameterLists.set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif

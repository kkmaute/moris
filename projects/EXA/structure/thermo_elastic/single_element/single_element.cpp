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
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
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

    std::string tNumElemsPerDim = "1, 1";
    std::string tDomainDims     = "0.01, 0.05";
    std::string tDomainOffset   = "0.0,0.0";
    std::string tDomainSidesets = "1,2,3,4";

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
    OPTParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        tParameterlist( 0 ).set( "is_optimization_problem", false );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        tParameterlist( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 ).set( "domain_dimensions", tDomainDims );
        tParameterlist( 0 ).set( "domain_offset", tDomainOffset );
        tParameterlist( 0 ).set( "domain_sidesets", tDomainSidesets );
        tParameterlist( 0 ).set( "lagrange_output_meshes", std::string( "0" ) );

        tParameterlist( 0 ).set( "lagrange_orders", tOrder );
        tParameterlist( 0 ).set( "lagrange_pattern", std::string( "0" ) );
        tParameterlist( 0 ).set( "bspline_orders", tOrder );
        tParameterlist( 0 ).set( "bspline_pattern", std::string( "0" ) );

        tParameterlist( 0 ).set( "truncate_bsplines", 1 );

        tParameterlist( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 ).set( "severity_level", 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        tParameterlist( 0 ).set( "decompose", true );
        tParameterlist( 0 ).set( "decomposition_type", std::string( "conformal" ) );
        tParameterlist( 0 ).set( "enrich", true );
        tParameterlist( 0 ).set( "basis_rank", std::string( "bspline" ) );
        tParameterlist( 0 ).set( "enrich_mesh_indices", std::string( "0" ) );
        tParameterlist( 0 ).set( "ghost_stab", false );
        tParameterlist( 0 ).set( "multigrid", false );
        tParameterlist( 0 ).set( "verbose", true );
        tParameterlist( 0 ).set( "print_enriched_ig_mesh", true );
        tParameterlist( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
        tParameterlist( 0 ).set( "output_mesh_file", tGENOutputFile );


        // Dummy Geometry
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Back_Wall" );
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Vector< Submodule_Parameter_Lists >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - STRUCTURE (ni-w-alloy?)
        //------------------------------------------------------------------------------

        // Density Shell
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropDensity" );
        tParameterList( 0 ).set( "function_parameters", tDensity );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // Heat Capacity Shell
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropHeatCapacity" );
        tParameterList( 0 ).set( "function_parameters", tHeatCapacity );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // Conductivity Shell
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropConductivity" );
        tParameterList( 0 ).set( "function_parameters", tConductivity );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // Youngs Modulus Shell
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropYoungsModulus" );
        tParameterList( 0 ).set( "function_parameters", tYoungsModulus );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // Poisson Ratio Shell
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropPoissonRatio" ) );
        tParameterList( 0 ).set( "function_parameters", tPoissonRatio );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Const" ) );

        // CTE for Shell
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropThermalExpansion" ) );
        tParameterList( 0 ).set( "function_parameters", tThermalExpansion );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Const" ) );

        //------------------------------------------------------------------------------
        // OTHER MATERIAL PARAMETERS
        //------------------------------------------------------------------------------

        // properties for bedding (supression for RBMs, both Shell and PCM)
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropBedding" );
        tParameterList( 0 ).set( "function_parameters", tBedding );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // Phase Change Temperature of PCM
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropPCTemp_Dummy" ) );
        tParameterList( 0 ).set( "function_parameters", "1.0e12" );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Const" ) );

        // Phase Change Temperature Range of PCM
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropPCconst_Dummy" ) );
        tParameterList( 0 ).set( "function_parameters", "100.0" );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Const" ) );

        // Cubic Phase State Function
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropPhaseStateFnct" ) );
        tParameterList( 0 ).set( "function_parameters", std::string( "2.0" ) );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Const" ) );

        // Dummy latent heat for non-pc material
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropLatentHeat_Dummy" ) );
        tParameterList( 0 ).set( "function_parameters", std::string( "0.0" ) );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Const" ) );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropPCTemp_Dummy" ) );
        tParameterList( 0 ).set( "function_parameters", std::string( "10000.0" ) );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Const" ) );

        //------------------------------------------------------------------------------
        // BOUNDARY CONDITIONS
        //------------------------------------------------------------------------------

        // heat flux from outside
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropImposedFlux" ) );
        tParameterList( 0 ).set( "function_parameters", tHeatLoad );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Const" ) );

        // reference temperature for thermal expansion
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropReferenceTemp" ) );
        tParameterList( 0 ).set( "function_parameters", tReferenceTemp );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Const" ) );

        // pressure load
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropTraction" ) );
        tParameterList( 0 ).set( "function_parameters", tPressure );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Neumann_U" ) );

        // pressure load
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropPressure" ) );
        tParameterList( 0 ).set( "function_parameters", tPressure );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Const" ) );

        // Dirichlet structure
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropDirichletStruct" ) );
        tParameterList( 0 ).set( "function_parameters", "0.0;0.0" );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Const" ) );

        // time continuity weights
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropWeightCurrent" ) );
        tParameterList( 0 ).set( "function_parameters", std::string( "100.0" ) );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Const" ) );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropWeightPrevious" ) );
        tParameterList( 0 ).set( "function_parameters", std::string( "100.0" ) );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Const" ) );

        // Initial Temperature
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropInitialCondition" ) );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Initial_Condition" ) );

        // Initial Structure
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropInitialConditionStruct" ) );
        tParameterList( 0 ).set( "value_function", std::string( "Func_Initial_Condition_Struct" ) );

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // DIFFUSION
        //------------------------------------------------------------------------------

        if ( tHaveThermo )
        {
            // create parameter list for constitutive model - shell - 1
            tParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
            tParameterList( 1 ).set( "constitutive_name", std::string( "CMDiffusion" ) );
            tParameterList( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
            tParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            tParameterList( 1 ).set( "properties", std::string( "PropConductivity , Conductivity;" ) + std::string( "PropDensity      , Density;" ) + std::string( "PropHeatCapacity , HeatCapacity" ) );
            }

        //------------------------------------------------------------------------------
        // LINEAR ELASTICITY
        //------------------------------------------------------------------------------

        if ( tHaveStruct )
        {
            // linear elasticity - shell - 1
            tParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
            tParameterList( 1 ).set( "constitutive_name", "CMStrucLinIso" );
            tParameterList( 1 ).set( "model_type",  fem::Model_Type::PLANE_STRESS ) ;
            tParameterList( 1 ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;

            if ( tHaveThermo )
            {
                tParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY;TEMP", "Displacement,Temperature" ) );
                tParameterList( 1 ).set( "properties", std::string( "PropYoungsModulus,    YoungsModulus;" ) + std::string( "PropPoissonRatio,     PoissonRatio;" ) + std::string( "PropThermalExpansion, CTE;" ) + std::string( "PropReferenceTemp,    ReferenceTemperature" ) );
            }
            else
            {
                tParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
                tParameterList( 1 ).set( "properties", std::string( "PropYoungsModulus,  YoungsModulus;" ) + std::string( "PropPoissonRatio,   PoissonRatio" ) );
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
            tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 ).set( "stabilization_name", std::string( "SPGGLSDiffusion" ) );
            tParameterList( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GGLS_DIFFUSION ) ;
            tParameterList( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
            tParameterList( 2 ).set( "leader_properties", std::string( "PropConductivity      , Conductivity;" ) + std::string( "PropDensity           , Density;" ) + std::string( "PropHeatCapacity      , HeatCapacity;" ) + std::string( "PropLatentHeat_Dummy  , LatentHeat;" ) + std::string( "PropPCTemp_Dummy      , PCTemp;" ) + std::string( "PropPhaseStateFnct    , PhaseStateFunction;" ) + std::string( "PropPCconst_Dummy     , PhaseChangeConst" ) );
            }

        //------------------------------------------------------------------------------
        // NITSCHE DIRICHLET
        //------------------------------------------------------------------------------

        if ( tHaveStruct )
        {
            // Displacements - Shell - back wall
            tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 ).set( "stabilization_name", "SPNitscheStruc" );
            tParameterList( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
            tParameterList( 2 ).set( "function_parameters", "100.0" );
            tParameterList( 2 ).set( "leader_properties", "PropYoungsModulus,Material" );
            }

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        //------------------------------------------------------------------------------
        // BULK IWGs
        //------------------------------------------------------------------------------

        // diffusion - Shell
        if ( tHaveThermo )
        {
            tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( 3 ).set( "IWG_name", std::string( "IWGDiffusionBulk" ) );
            tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
            tParameterList( 3 ).set( "dof_residual", std::string( "TEMP" ) );
            tParameterList( 3 ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
            tParameterList( 3 ).set( "leader_constitutive_models", std::string( "CMDiffusion,Diffusion" ) );
            tParameterList( 3 ).set( "mesh_set_names", tBulk );

            if ( tUseGGLS )
            {
                tParameterList( 3 ).set( "stabilization_parameters", std::string( "SPGGLSDiffusion,GGLSParam" ) );
            }
        }

        // linear elasticity - Shell
        if ( tHaveStruct )
        {
            tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( 3 ).set( "IWG_name", "IWGStructShell" );
            tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
            tParameterList( 3 ).set( "dof_residual", "UX,UY" );

            if ( tHaveThermo )
            {
                tParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
            }
            else
            {
                tParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            }

            tParameterList( 3 ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
            // tParameterList( 3 ).set( "leader_properties",          "PropBedding,Bedding" );
            tParameterList( 3 ).set( "mesh_set_names", tBulk );
            }

        //------------------------------------------------------------------------------
        // NEUMANN BCs - IWGs
        //------------------------------------------------------------------------------

        // heat flux on outside of Shell
        if ( tHaveThermo )
        {
            tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( 3 ).set( "IWG_name", std::string( "IWGInletFlux" ) );
            tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
            tParameterList( 3 ).set( "dof_residual", std::string( "TEMP" ) );

            if ( tHaveStruct )
            {
                tParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
            }
            else
            {
                tParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
            }

            tParameterList( 3 ).set( "leader_properties", std::string( "PropImposedFlux,Neumann" ) );
            tParameterList( 3 ).set( "mesh_set_names", tTop );
            }

        // pressure pulling on outside of Shell
        if ( tHaveStruct )
        {
            tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( 3 ).set( "IWG_name", "IWGNeumannPressure" );
            tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
            tParameterList( 3 ).set( "dof_residual", "UX,UY" );

            if ( tHaveThermo )
            {
                tParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
            }
            else
            {
                tParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            }

            // tParameterList( 3 ).set( "leader_properties",          "PropTraction,Traction");
            tParameterList( 3 ).set( "leader_properties", "PropPressure,Pressure" );
            tParameterList( 3 ).set( "mesh_set_names", tLeft );
            }

        //------------------------------------------------------------------------------
        // DIRICHLET BCS - IWGs
        //------------------------------------------------------------------------------

        // displacements - shell - back wall
        if ( tHaveStruct )
        {
            tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( 3 ).set( "IWG_name", "IWGDirichletStruct" );
            tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
            tParameterList( 3 ).set( "dof_residual", "UX,UY" );

            if ( tHaveThermo )
            {
                tParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
            }
            else
            {
                tParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            }

            tParameterList( 3 ).set( "leader_properties", "PropDirichletStruct,Dirichlet" );
            tParameterList( 3 ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
            tParameterList( 3 ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
            tParameterList( 3 ).set( "mesh_set_names", tBottom );
            }

        //------------------------------------------------------------------------------
        // IWGs - TIME CONTINUITY
        //------------------------------------------------------------------------------

        // Time continuity
        if ( tHaveThermo )
        {
            tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( 3 ).set( "IWG_name", std::string( "IWGTimeContinuityTemp" ) );
            tParameterList( 3 ).set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
            tParameterList( 3 ).set( "dof_residual", std::string( "TEMP" ) );

            if ( tHaveStruct )
            {
                tParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
            }
            else
            {
                tParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
            }

            tParameterList( 3 ).set( "leader_properties", std::string( "PropWeightCurrent,WeightCurrent;" ) + std::string( "PropWeightPrevious,WeightPrevious;" ) + std::string( "PropInitialCondition,InitialCondition" ) );
            tParameterList( 3 ).set( "mesh_set_names", tBulk );
            tParameterList( 3 ).set( "time_continuity", true );
            }

        // Time continuity Structure
        //         if ( tHaveStruct )
        //         {
        //             tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        //             tParameterList( 3 ).set( "IWG_name",                   std::string("IWGTimeContinuityStruct") );
        //             tParameterList( 3 ).set( "IWG_type",                    fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        //             tParameterList( 3 ).set( "dof_residual",               std::string("UX,UY") );
        //
        //             if ( tHaveThermo ) {
        //                 tParameterList( 3 ).set( "leader_dof_dependencies",    "UX,UY;TEMP");
        //             } else {
        //                 tParameterList( 3 ).set( "leader_dof_dependencies",    "UX,UY");
        //             }
        //
        //             tParameterList( 3 ).set( "leader_properties", std::string("PropWeightCurrent,WeightCurrent;")      +
        //                                                                          std::string("PropWeightPrevious,WeightPrevious;")    +
        //                                                                          std::string("PropInitialConditionStruct,InitialCondition") );
        //             tParameterList( 3 ).set( "mesh_set_names",             tBulk );
        //             tParameterList( 3 ).set( "time_continuity",            true );
        //             //         }

        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        if ( tHaveThermo )
        {
            // Nodal Temperature IQI
            tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
            tParameterList( 4 ).set( "IQI_name", "IQIBulkTEMP" );
            tParameterList( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
            tParameterList( 4 ).set( "dof_quantity", "TEMP" );
            tParameterList( 4 ).set( "leader_dof_dependencies", std::string( "TEMP" ) );
            tParameterList( 4 ).set( "vectorial_field_index", 0 );
            tParameterList( 4 ).set( "mesh_set_names", tBulk );
            }

        // Volume IQI - TotalDomain - use once to find total volume to compute max dof
        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQITotalVolume" );
        tParameterList( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;

        if ( tHaveThermo )
        {
            tParameterList( 4 ).set( "leader_dof_dependencies", "TEMP" );
        }
        else if ( tHaveStruct )
        {
            tParameterList( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        }
        else if ( tHaveThermo && tHaveStruct )
        {
            tParameterList( 4 ).set( "leader_dof_dependencies", "UX,UY;TEMP" );
        }

        tParameterList( 4 ).set( "mesh_set_names", tBulk );

        if ( tHaveStruct )
        {
            // X-displacement
            tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
            tParameterList( 4 ).set( "IQI_name", "IQIBulkDISPX" );
            tParameterList( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
            tParameterList( 4 ).set( "dof_quantity", "UX,UY" );
            tParameterList( 4 ).set( "leader_dof_dependencies", "UX,UY" );
            tParameterList( 4 ).set( "vectorial_field_index", 0 );
            tParameterList( 4 ).set( "mesh_set_names", tBulk );

            // Y-displacement
            tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
            tParameterList( 4 ).set( "IQI_name", "IQIBulkDISPY" );
            tParameterList( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
            tParameterList( 4 ).set( "dof_quantity", "UX,UY" );
            tParameterList( 4 ).set( "leader_dof_dependencies", "UX,UY" );
            tParameterList( 4 ).set( "vectorial_field_index", 1 );
            tParameterList( 4 ).set( "mesh_set_names", tBulk );
            }

        // create computation parameter list
        tParameterList( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 8 );


        // Thermo Only solver config
        if ( !tHaveStruct && tHaveThermo )
        {
            std::cout << "Thermo ONLY !!! \n"
                      << std::flush;


            if ( tUseBelos )
            {
                tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
            }
            else
            {
                tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
            }

            tParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

            tParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
            tParameterlist( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
            tParameterlist( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
            tParameterlist( 2 ).set( "NLA_max_iter", tNLA_max_iter );
            tParameterlist( 2 ).set( "NLA_combined_res_jac_assembly", true );

            tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            tParameterlist( 3 ).set( "NLA_DofTypes", "TEMP" );

            tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            tParameterlist( 3 ).set( "NLA_DofTypes", "TEMP" );
            tParameterlist( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0" );

            tParameterlist( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
            tParameterlist( 4 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
            tParameterlist( 4 ).set( "TSA_Time_Frame", tTSA_Time_Frame );

            tParameterlist( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
            tParameterlist( 5 ).set( "TSA_DofTypes", "TEMP" );
            tParameterlist( 5 ).set( "TSA_Initialize_Sol_Vec", "TEMP,1.0" );
            tParameterlist( 5 ).set( "TSA_Output_Indices", "0" );
            tParameterlist( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );
            tParameterlist( 5 ).set( "TSA_time_level_per_type", "TEMP,2" );
        }    // end: thermo only

        // Struct Only solver config
        else if ( tHaveStruct && !tHaveThermo )
        {
            std::cout << "Struct ONLY !!! \n"
                      << std::flush;


            if ( tUseBelos )
                tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
            else
                tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

            tParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

            tParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
            tParameterlist( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
            tParameterlist( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
            tParameterlist( 2 ).set( "NLA_max_iter", tNLA_max_iter );
            tParameterlist( 2 ).set( "NLA_combined_res_jac_assembly", true );

            tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            tParameterlist( 3 ).set( "NLA_DofTypes", "UX,UY" );

            tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            tParameterlist( 3 ).set( "NLA_DofTypes", "UX,UY" );
            tParameterlist( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0" );

            tParameterlist( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
            tParameterlist( 4 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
            tParameterlist( 4 ).set( "TSA_Time_Frame", tTSA_Time_Frame );

            tParameterlist( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
            tParameterlist( 5 ).set( "TSA_DofTypes", "UX,UY" );
            tParameterlist( 5 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0" );
            tParameterlist( 5 ).set( "TSA_Output_Indices", "0" );
            tParameterlist( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );
            tParameterlist( 5 ).set( "TSA_time_level_per_type", "UX,2;UY,2" );
        }    // end: struct only

        else if ( tHaveStruct && tHaveThermo && !tHaveStaggered )
        {
            std::cout << "NON-Staggered SOLVE !!! \n"
                      << std::flush;


            if ( tUseBelos )
            {
                tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
            }
            else
            {
                tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
            }

            tParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

            tParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
            tParameterlist( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
            tParameterlist( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
            tParameterlist( 2 ).set( "NLA_max_iter", tNLA_max_iter );
            tParameterlist( 2 ).set( "NLA_combined_res_jac_assembly", true );

            tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            tParameterlist( 3 ).set( "NLA_DofTypes", "UX,UY;TEMP" );

            tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            tParameterlist( 3 ).set( "NLA_DofTypes", "UX,UY;TEMP" );
            tParameterlist( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0" );

            tParameterlist( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
            tParameterlist( 4 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
            tParameterlist( 4 ).set( "TSA_Time_Frame", tTSA_Time_Frame );

            tParameterlist( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
            tParameterlist( 5 ).set( "TSA_DofTypes", "UX,UY;TEMP" );
            tParameterlist( 5 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0;TEMP,1.0" );
            tParameterlist( 5 ).set( "TSA_Output_Indices", "0" );
            tParameterlist( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );
            if ( tHaveMixedSolverConfig )
                tParameterlist( 5 ).set( "TSA_time_level_per_type", "UX,1;UY,1;TEMP,2" );
            else
                tParameterlist( 5 ).set( "TSA_time_level_per_type", "UX,2;UY,2;TEMP,2" );
        }    // end: coupled & NOT staggered

        // Coupled Thermo Structural with Staggered solver config
        else if ( tHaveStruct && tHaveThermo && tHaveStaggered )
        {
            std::cout << "Staggered SOLVE !!! \n"
                      << std::flush;

            if ( tUseBelos )
            {
                tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
            }
            else
            {
                tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
            }

            tParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

            // ----------------------------------------------------------

            tParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
            tParameterlist( 2 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
            tParameterlist( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
            tParameterlist( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
            tParameterlist( 2 ).set( "NLA_max_iter", tNLA_max_iter );
            tParameterlist( 2 ).set( "NLA_combined_res_jac_assembly", true );

            tParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
            tParameterlist( 2 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;

            // ----------------------------------------------------------

            tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            tParameterlist( 3 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
            tParameterlist( 3 ).set( "NLA_DofTypes", "UX,UY" );

            tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            tParameterlist( 3 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
            tParameterlist( 3 ).set( "NLA_DofTypes", "TEMP" );

            tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
            tParameterlist( 3 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;
            tParameterlist( 3 ).set( "NLA_Sub_Nonlinear_Solver", "1,0" );
            tParameterlist( 3 ).set( "NLA_DofTypes", "UX,UY;TEMP" );
            tParameterlist( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );

            // ----------------------------------------------------------

            tParameterlist( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
            tParameterlist( 4 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
            tParameterlist( 4 ).set( "TSA_Time_Frame", tTSA_Time_Frame );
            tParameterlist( 4 ).set( "TSA_Nonlinear_Solver", 2 );

            // ----------------------------------------------------------

            tParameterlist( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
            tParameterlist( 5 ).set( "TSA_DofTypes", "UX,UY;TEMP" );
            tParameterlist( 5 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0;TEMP,1.0" );
            tParameterlist( 5 ).set( "TSA_Output_Indices", "0" );
            tParameterlist( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );
            if ( tHaveMixedSolverConfig )
                tParameterlist( 5 ).set( "TSA_time_level_per_type", "UX,1;UY,1;TEMP,2" );
            else
                tParameterlist( 5 ).set( "TSA_time_level_per_type", "UX,2;UY,2;TEMP,2" );
        }    // end: coupled & staggered

        else
        {
            std::cout << "Input file: incomplete physics configuration. Check what physics are used. \n"
                      << std::flush;
        }

        tParameterlist( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        tParameterlist( 6 ).set( "SOL_save_operator_to_matlab", "jacobian.dat" );

        tParameterlist( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    void
    VISParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        tParameterlist( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        tParameterlist( 0 ).set( "Set_Names", tBulk );

        if ( tHaveStruct and tHaveThermo )
        {
            tParameterlist( 0 ).set( "Field_Names", "TEMP,UX,UY,VOLUME" );
            tParameterlist( 0 ).set( "Field_Type", "NODAL,NODAL,NODAL,GLOBAL" );
            tParameterlist( 0 ).set( "IQI_Names", "IQIBulkTEMP,IQIBulkDISPX,IQIBulkDISPY,IQITotalVolume" );
        }
        else if ( tHaveStruct )
        {
            tParameterlist( 0 ).set( "Field_Names", "UX,UY,VOLUME" );
            tParameterlist( 0 ).set( "Field_Type", "NODAL,NODAL,GLOBAL" );
            tParameterlist( 0 ).set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY,IQITotalVolume" );
        }
        else if ( tHaveThermo )
        {
            tParameterlist( 0 ).set( "Field_Names", "TEMP,VOLUME" );
            tParameterlist( 0 ).set( "Field_Type", "NODAL,GLOBAL" );
            tParameterlist( 0 ).set( "IQI_Names", "IQIBulkTEMP,IQITotalVolume" );
        }

        tParameterlist( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif

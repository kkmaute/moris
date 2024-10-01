/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Contact_Sensitivity_Test.cpp
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
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

// global variable for interpolation order
extern moris::uint gLevelSetInterpolationOrder;
extern moris::uint gFEMInterpolationOrder;
extern moris::uint gLagrMeshInterpolationOrder;

// problem dimension: 2D or 3D
extern moris::uint gDim;

// use Bspline parameterization or primitives
extern bool gUseBspline;

// test case index
extern moris::uint gTestCaseIndex;

/* To visualize results in paraview and check continuity of normal displacements use
 * calculator in paraview with:
 *
 *                    U_X*sin(75*3.142/180)+U_Y*cos(75*3.142/180)
 */

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------------

namespace moris
{
    //------------------------------------------------------------------------------

    // Main problem parameters

    std::string tName = "Contact_Sensitivity_Test";

    bool tIs3D     = gDim == 3 ? true : false;
    bool tIsOpt    = true;
    bool tUseGhost = true;

    // x-position of interface
    real tInterfaceX = 0.638;

    // hole radius
    real tHoleRadius = 0.166;

    // Nitsche interface penalty
    real tNitschePenalty = 10.0;

    // Prescribed load or displacement
    bool tPrescribedTraction = false;

    // use Contact
    bool tUseContact = true;

    // subset of ADVs for FD check
    std::string tAdvIndicesForFD = gUseBspline ? ( gDim == 3 ? "7,15,21" : "3,4,5" ) : "1";

    // FD in adjoint
    real tFEMFdEpsilon = 1.0e-5;

    // FD step size in sweep
    std::string tFDsweep = "1.0e-5";

    // Number of constraints
    uint tNumConstraints = 4;

    // background mesh parameters
    std::string tNumElementsPerDir = tIs3D ? "1,1,1" : "1,1";
    std::string tDimensions        = tIs3D ? "1,1,1" : "1,1";
    std::string tOffSet            = tIs3D ? "0.0,0.0,0.0" : "0.0,0.0";

    bool tUseBsplineForLevelset = gUseBspline;    // use simple plane to define intersection
    real tPlaneTilde            = 75;             // tilde of plane in degree; 90 being straight up
    bool tUseMultiLinear        = false;          // use multi linear intersection computation

    int tLevelsetOrder = gLevelSetInterpolationOrder;
    int tDispOrder     = gFEMInterpolationOrder;
    int tLagMeshOrder  = gLagrMeshInterpolationOrder;

    int tLevelsetInitialRef = 1;
    int tDispInitialRef     = 3;

    int tRefineBuffer = 0;

    // note: pattern 0 - Levelset field  pattern 1 - displacement field
    std::string tLagrangeOrder   = std::to_string( tLagMeshOrder );
    std::string tBsplineOrder    = std::to_string( tLevelsetOrder ) + "," + std::to_string( tDispOrder );
    std::string tInitialRef      = std::to_string( tLevelsetInitialRef ) + "," + std::to_string( tDispInitialRef );
    std::string tLagrangePattern = tLevelsetInitialRef > tDispInitialRef ? "0" : "1";

    //------------------------------------------------------------------------------
    // Derived problem parameters

    std::string tProblemConfig = "_" + std::to_string( gTestCaseIndex );

    std::string tOutputFileName = tName + tProblemConfig + ".exo";
    std::string tLibraryName    = tName + ".so";
    std::string tGENOutputFile  = tLagMeshOrder == 3 ? "" : tName + tProblemConfig + "_GEN" + ".exo";
    std::string tHDF5FileName   = tName + tProblemConfig + "_SEN" + ".hdf5";

    std::string tMaterial1Sets = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tMaterial2Sets = "HMR_dummy_n_p1,HMR_dummy_c_p1";

    std::string tTotalDomainSets = tMaterial1Sets + "," + tMaterial2Sets;

    std::string tMaterial1Ghost = "ghost_p2";
    std::string tMaterial2Ghost = "ghost_p1";

    std::string tLoadSSsets   = "SideSet_2_n_p1,SideSet_2_c_p1";
    std::string tSupportSSets = "SideSet_4_n_p2,SideSet_4_c_p2";

    std::string tMaterial12SSets = "iside_b0_2_b1_1";

    std::string tMaterial12DSets = "dbl_iside_p0_2_p1_1";    // watch:from Material 1 (p2) to Material 2 (p1)

    std::string tInterfaces = tLoadSSsets + "," + tSupportSSets + "," + tMaterial12SSets + "," + tMaterial12DSets;

    std::string tDofStrg = tIs3D ? "UX,UY,UZ" : "UX,UY";

    std::string tDirichletStr = tIs3D ? "0.0;0.0;0.0" : "0.0;0.0";

    std::string tDirichletLoadStr = tIs3D ? "-0.1;0.0;0.0" : "-0.1;0.0";

    std::string tNeumannStr = tIs3D ? "1.0;0.0;0.0" : "1.0;0.0";

    //------------------------------------------------------------------------------

    // Hole pattern
    real
    Plane_2D3D(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&    aGeometryParameters )
    {
        real tLSval = aCoordinates( 0 ) - tInterfaceX;

        // clean return value to return non-zero value
        return -tLSval;
    }

    //------------------------------------------------------------------------------

    void
    tLevelSetFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );
    }

    void
    tDerLevelSetFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->N();
    }

    //------------------------------------------------------------------------------

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //------------------------------------------------------------------------------

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( tNumConstraints, 1, 1 );

        return tConstraintTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );

        real obj0 = aCriteria( 0 );
        real obj1 = aCriteria( 1 );
        real obj2 = aCriteria( 2 );
        real obj3 = aCriteria( 3 );

        tObjectives( 0, 0 ) = obj0 + obj1 + obj2 + obj3;

        std::cout << "% --------------------------------- % \n";
        std::cout << "Objective                 = " << tObjectives( 0, 0 ) << " \n";
        std::cout << "Strain Energy (Material1) = " << aCriteria( 0 ) << " ( " << obj0 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Strain Energy (Material2) = " << aCriteria( 1 ) << " ( " << obj1 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Volume                    = " << aCriteria( 2 ) << " ( " << obj2 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Perimeter of Material 1-2 = " << aCriteria( 3 ) << " ( " << obj3 / tObjectives( 0, 0 ) << " )\n";
        std::cout << " \n";

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, tNumConstraints );

        tConstraints( 0 ) = aCriteria( 0 );
        tConstraints( 1 ) = aCriteria( 1 );
        tConstraints( 2 ) = aCriteria( 2 );
        tConstraints( 3 ) = aCriteria( 3 );

        return tConstraints;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );

        return tDObjectiveDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.size(), 0.0 );

        tDObjectiveDCriteria( 0 ) = 1.0;
        tDObjectiveDCriteria( 1 ) = 1.0;
        tDObjectiveDCriteria( 2 ) = 1.0;
        tDObjectiveDCriteria( 3 ) = 1.0;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( tNumConstraints, aADVs.size(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( tNumConstraints, aCriteria.size(), 0.0 );

        tDConstraintDCriteria( 0, 0 ) = 1.0;

        tDConstraintDCriteria( 1, 1 ) = 1.0;

        tDConstraintDCriteria( 2, 2 ) = 1.0;

        tDConstraintDCriteria( 3, 3 ) = 1.0;

        return tDConstraintDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 0 );
        tParameterlist( 2 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = moris::prm::create_opt_problem_parameter_list();
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", tIsOpt );
        tParameterlist( 0 )( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 )( 0 ).set( "library", tLibraryName );
        tParameterlist( 0 )( 0 ).set( "restart_file", "" );

        tParameterlist( 2 )( 0 ) = moris::prm::create_sweep_parameter_list();
        tParameterlist( 2 )( 0 ).set( "hdf5_path", tHDF5FileName );
        tParameterlist( 2 )( 0 ).set( "num_evaluations_per_adv", "1" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_type", "all" );
        tParameterlist( 2 )( 0 ).set( "finite_difference_epsilons", tFDsweep );
        tParameterlist( 2 )( 0 ).set( "finite_difference_adv_indices", tAdvIndicesForFD );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElementsPerDir );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", tDimensions );
        tParameterlist( 0 )( 0 ).set( "domain_offset", tOffSet );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", tLagrangeOrder );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", tLagrangePattern );

        tParameterlist( 0 )( 0 ).set( "bspline_orders", tBsplineOrder );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0,1" );

        tParameterlist( 0 )( 0 ).set( "initial_refinement", tInitialRef );
        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0,1" );
        tParameterlist( 0 )( 0 ).set( "use_advanced_T_matrix_scheme", 1 );

        tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0,1" );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 )( 0 ).set( "refinement_buffer", tRefineBuffer );
        tParameterlist( 0 )( 0 ).set( "staircase_buffer", tRefineBuffer );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose", true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 )( 0 ).set( "enrich", true );
        tParameterlist( 0 )( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", "0,1" );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );

        tParameterlist( 0 )( 0 ).set( "ghost_stab", tUseGhost );
        tParameterlist( 0 )( 0 ).set( "visualize_ghost", tUseGhost );

        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 3 );

        tParameterlist( 0 ).push_back( moris::prm::create_gen_parameter_list() );
        tParameterlist( 0 )( 0 ).set( "IQI_types",
                "IQIBulkStrainEnergy_Material1", "IQIBulkStrainEnergy_Material2", "IQIBulkVolume_Material1", "IQIPerimeter_InterfaceMaterial12" );
        tParameterlist( 0 )( 0 ).set( "output_mesh_file", tGENOutputFile );
        tParameterlist( 0 )( 0 ).set( "time_offset", 10.0 );

        Matrix< DDUMat > tPhaseMap( 4, 1, 0 );
        tPhaseMap( 0 ) = 0;
        tPhaseMap( 1 ) = 0;
        tPhaseMap( 2 ) = 1;
        tPhaseMap( 3 ) = 2;
        tParameterlist( 0 )( 0 ).set( "phase_table", moris::ios::stringify( tPhaseMap ) );

        tParameterlist( 0 )( 0 ).set( "print_phase_table", true );

        // init geometry counter
        uint tGeoCounter   = 0;
        uint tParamCounter = 0;

        const real pi = std::acos( -1 );

        if ( tIs3D )
        {
            tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::SPHERE ) );
            tParameterlist( 1 )( tGeoCounter ).set( "center_x", 0.5 );
            tParameterlist( 1 )( tGeoCounter ).set( "center_y", 0.5 );
            tParameterlist( 1 )( tGeoCounter ).set( "center_z", 0.5 );
            tParameterlist( 1 )( tGeoCounter ).set( "radius", tHoleRadius );
        }
        else
        {
            tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
            tParameterlist( 1 )( tGeoCounter ).set( "center_x", 0.5 );
            tParameterlist( 1 )( tGeoCounter ).set( "center_y", 0.5 );
            tParameterlist( 1 )( tGeoCounter ).set( "radius", tHoleRadius );
        }
        tGeoCounter++;

        // initialize geometry
        if ( tUseBsplineForLevelset )
        {
            tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            tParameterlist( 1 )( tGeoCounter ).set( "name", "ADVfield" );
            tParameterlist( 1 )( tGeoCounter ).set( "center_x", tInterfaceX );
            tParameterlist( 1 )( tGeoCounter ).set( "center_y", 0.0 );
            tParameterlist( 1 )( tGeoCounter ).set( "normal_x", -1.0 * std::sin( tPlaneTilde / 180.0 * pi ) );
            tParameterlist( 1 )( tGeoCounter ).set( "normal_y", -1.0 * std::cos( tPlaneTilde / 180.0 * pi ) );
            tParameterlist( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tUseMultiLinear );

            if ( tIsOpt )
            {
                tParameterlist( 1 )( tGeoCounter ).set( "discretization_mesh_index", 0 );
                tParameterlist( 1 )( tGeoCounter ).set( "discretization_lower_bound", -2.0 );
                tParameterlist( 1 )( tGeoCounter ).set( "discretization_upper_bound", 2.0 );
            }
            tGeoCounter++;

            // Levelset property
            tParameterlist( 2 ).push_back( moris::prm::create_gen_property_parameter_list( gen::Field_Type::SCALED_FIELD ) );

            tParameterlist( 2 )( tParamCounter ).set( "name", "LevelsetField" );
            tParameterlist( 2 )( tParamCounter ).set( "dependencies", "ADVfield" );
            tParameterlist( 2 )( tParamCounter ).set( "scaling_factor", 1.0 );
            tParameterlist( 2 )( tParamCounter ).set( "pdv_type", "LS1" );
            tParameterlist( 2 )( tParamCounter ).set( "pdv_mesh_set_names", tTotalDomainSets );
            tParamCounter++;
        }
        else
        {
            tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            tParameterlist( 1 )( tGeoCounter ).set( "center_x", tInterfaceX * 0.5, tInterfaceX, tInterfaceX / 0.5 );
            tParameterlist( 1 )( tGeoCounter ).set( "center_y", 0.0 );
            tParameterlist( 1 )( tGeoCounter ).set( "normal_x", -1.0 * std::sin( tPlaneTilde / 180.0 * pi ) );
            tParameterlist( 1 )( tGeoCounter ).set( "normal_y", -1.0 * std::cos( tPlaneTilde / 180.0 * pi ) );
            tGeoCounter++;

            // Levelset property
            tParameterlist( 2 ).push_back( moris::prm::create_gen_property_parameter_list( gen::Field_Type::CONSTANT ) );

            tParameterlist( 2 )( tParamCounter ).set( "name", "LevelsetField" );
            tParameterlist( 2 )( tParamCounter ).set( "constant", 1.0 );
            tParameterlist( 2 )( tParamCounter ).set( "pdv_type", "LS1" );
            tParameterlist( 2 )( tParamCounter ).set( "pdv_mesh_set_names", tTotalDomainSets );

            if ( tIsOpt )
            {
                tParameterlist( 2 )( tParamCounter ).set( "discretization_mesh_index", 0 );
                tParameterlist( 2 )( tParamCounter ).set( "discretization_lower_bound", -2.0 );
                tParameterlist( 2 )( tParamCounter ).set( "discretization_upper_bound", 2.0 );
            }

            tParamCounter++;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( Vector< Vector< Parameter_List > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------
        // init property counter
        uint tPropCounter = 0;

        // create parameter list for property 1
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs1" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "2.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungs2" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 2
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropBedding" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0e-6" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichletU" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDirichletStr );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 4
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDirichletULoad" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDirichletLoadStr );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 10
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropTraction" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tNeumannStr );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // create parameter list for property 7
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisson1" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropPoisson2" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.3" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ) = prm::create_property_parameter_list();
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropLevelSet" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "1.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "tLevelSetFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_derivative_functions", "tDerLevelSetFunc" );
        tParameterList( 0 )( tPropCounter ).set( "dv_dependencies", "LS1" );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // init CM counter
        uint tCMCounter = 0;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso_Material1" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungs1,YoungsModulus;PropPoisson1,PoissonRatio" );
        tCMCounter++;

        // create parameter list for constitutive model 1
        tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso_Material2" );
        tParameterList( 1 )( tCMCounter ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        tParameterList( 1 )( tCMCounter ).set( "properties", "PropYoungs2,YoungsModulus;PropPoisson2,PoissonRatio" );
        tCMCounter++;

        //------------------------------------------------------------------------------
        // init SP counter
        uint tSPCounter = 0;

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheDirichletBC" );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::to_string( tNitschePenalty ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungs1,Material" );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPNitscheMaterial12Interface" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::to_string( tNitschePenalty ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", std::string( "PropYoungs1,Material" ) );
        tParameterList( 2 )( tSPCounter ).set( "follower_properties", std::string( "PropYoungs2,Material" ) );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPGhost_Material1" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "0.005" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", std::string( "PropYoungs1,Material" ) );
        tSPCounter++;

        tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_name", std::string( "SPGhost_Material2" ) );
        tParameterList( 2 )( tSPCounter ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        tParameterList( 2 )( tSPCounter ).set( "function_parameters", std::string( "0.005" ) );
        tParameterList( 2 )( tSPCounter ).set( "leader_properties", std::string( "PropYoungs1,Material" ) );
        tSPCounter++;

        //------------------------------------------------------------------------------
        // init IWG counter
        uint tIWGCounter = 0;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkU_Material1" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso_Material1,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropBedding,Bedding" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tMaterial1Sets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGBulkU_Material2" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso_Material2,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropBedding,Bedding" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tMaterial2Sets );
        tIWGCounter++;

        tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletU" );
        tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichletU,Dirichlet" );
        tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso_Material1,ElastLinIso" );
        tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tSupportSSets );
        tIWGCounter++;

        if ( tPrescribedTraction )
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGTraction" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropTraction,Traction" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tLoadSSsets );
            tIWGCounter++;
        }
        else
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletULoad" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichletULoad,Dirichlet" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso_Material1,ElastLinIso" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tLoadSSsets );
            tIWGCounter++;
        }

        if ( tUseContact )
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGMaterial12Interface" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_CONTACT_UNSYMMETRIC_NITSCHE ) ;
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", std::string( "CMStrucLinIso_Material1,ElastLinIso" ) );
            tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models", std::string( "CMStrucLinIso_Material2,ElastLinIso" ) );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPNitscheMaterial12Interface,NitscheInterface" ) );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tMaterial12DSets );
            tIWGCounter++;
        }
        else
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGMaterial12Interface" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", std::string( "CMStrucLinIso_Material1,ElastLinIso" ) );
            tParameterList( 3 )( tIWGCounter ).set( "follower_constitutive_models", std::string( "CMStrucLinIso_Material2,ElastLinIso" ) );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPNitscheMaterial12Interface,NitscheInterface" ) );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tMaterial12DSets );
            tIWGCounter++;
        }

        if ( tUseGhost )
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGGhostMaterial1" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGhost_Material1,GhostSP" ) );
            tParameterList( 3 )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tMaterial1Ghost );
            tIWGCounter++;

            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGGhostMaterial2" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "follower_dof_dependencies", tDofStrg );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", std::string( "SPGhost_Material2,GhostSP" ) );
            tParameterList( 3 )( tIWGCounter ).set( "ghost_order", (uint)tDispOrder );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tMaterial2Ghost );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // init IQI counter
        uint tIQICounter = 0;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUX" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomainSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUY" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "dof_quantity", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomainSets );
        tIQICounter++;

        if ( tIs3D )
        {
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkUZ" );
            tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
            tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", tDofStrg );
            tParameterList( 4 )( tIQICounter ).set( "dof_quantity", tDofStrg );
            tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 2 );
            tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomainSets );
            tIQICounter++;
        }

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQILevelSet" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::PROPERTY ) ;
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropLevelSet,Property" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tTotalDomainSets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy_Material1" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso_Material1,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMaterial1Sets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkStrainEnergy_Material2" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "leader_constitutive_models", "CMStrucLinIso_Material2,Elast" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMaterial2Sets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume_Material1" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropLevelSet,Density" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMaterial1Sets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkVolume_Material2" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        tParameterList( 4 )( tIQICounter ).set( "leader_properties", "PropLevelSet,Density" );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMaterial2Sets );
        tIQICounter++;

        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIPerimeter_InterfaceMaterial12" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", tDofStrg );
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tMaterial12SSets );
        tIQICounter++;

        // create computation  parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
        tParameterList( 5 )( 0 ).set( "print_physics_model", false );

        tParameterList( 5 )( 0 ).set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        tParameterList( 5 )( 0 ).set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        ///*
        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

#ifdef MORIS_USE_MUMPS
        tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Mumps" );
#else
        tParameterlist( 0 )( 0 ).set( "Solver_Type", "Amesos_Superludist" );
#endif

        //*/

        /*
        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL );

        // Solver type: GMRES, Flexible GMRES, Block CG , PseudoBlockCG, Stochastic CG, Recycling GMRES, Recycling CG, MINRES, LSQR, TFQMR
        //              Pseudoblock TFQMR, Seed GMRES, Seed CG
        tParameterlist( 0 )( 0 ).set( "Solver Type" ,  "GMRES" );

        // Diagnostics: Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails
        sint tVerbosity = Belos::Errors; // + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails;
        tParameterlist( 0 )( 0 ).set( "Verbosity" , tVerbosity );

        // Maximum number of blocks in Krylov factorization
        tParameterlist( 0 )( 0 ).set( "Num Blocks", 500   );

        // Block size to be used by iterative solver
        tParameterlist( 0 )( 0 ).set( "Block Size", 1   );

        // Allowable Belos solver iterations
        tParameterlist( 0 )( 0 ).set( "Maximum Iterations" , 500 );

        // Allowable Belos solver iterations
        //tParameterlist( 0 )( 0 ).set( "Maximum Restarts" ,  );

        // Convergence criteria
        tParameterlist( 0 )( 0 ).set( "Convergence Tolerance" ,  1e-12 );

        // Preconditioner
        //tParameterlist( 0 )( 0 ).set( "ifpack_prec_type",  "ILU");
        //tParameterlist( 0 )( 0 ).set( "fact: level-of-fill",  3 );

        //tParameterlist( 0 )( 0 ).set( "ifpack_prec_type",  "ILUT");
        //tParameterlist( 0 )( 0 ).set( "fact: ilut level-of-fill", 1.0 );
        //tParameterlist( 0 )( 0 ).set( "fact: drop tolerance", 1e-1 );

        //tParameterlist( 0 )( 0 ).set( "ifpack_prec_type",  "Amesos");
        //tParameterlist( 0 )( 0 ).set( "amesos: solver type", "Amesos_Umfpack");

        tParameterlist( 0 )( 0 ).set( "ml_prec_type",  "NSSA");
         */

        /*
        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AZTEC_IMPL );

        //options are: AZ_gmres, AZ_gmres_condnum, AZ_cg, AZ_cg_condnum, AZ_cgs, AZ_tfqmr, AZ_bicgstab
        tParameterlist( 0 )( 0 ).set( "AZ_solver" ,  AZ_gmres );

            // Allowable Aztec solver iterations
        tParameterlist( 0 )( 0 ).set( "AZ_max_iter", 500   );

            // Allowable Aztec iterative residual
        tParameterlist( 0 )( 0 ).set( "rel_residual" , 1e-08 );

        // set Az_conv -convergence criteria
        // options are AZ_r0, AZ_rhs, AZ_Anorm, AZ_noscaled, AZ_sol
        tParameterlist( 0 )( 0 ).set( "AZ_conv" ,  AZ_r0 );

        // set Az_diagnostic parameters
        // Set whether or not diagnostics for every linear iteration are printed or not. options are AZ_all, AZ_none
        tParameterlist( 0 )( 0 ).set( "AZ_diagnostics" ,  AZ_all );

        // set AZ_output options
        // options are AZ_all, AZ_none, AZ_warnings, AZ_last, AZ_summary
        tParameterlist( 0 )( 0 ).set( "AZ_output" ,  AZ_all );

        // Determines the submatrices factored with the domain decomposition algorithms
        // Option to specify with how many rows from other processors each processor\u2019s local submatrix is augmented.
        tParameterlist( 0 )( 0 ).set( "AZ_overlap" , 1 );

        // Determines how overlapping subdomain results are combined when different processors have computed different values for the same unknown.
        // Options are AZ_standard, AZ_symmetric
        tParameterlist( 0 )( 0 ).set( "AZ_type_overlap" , AZ_standard );

        // Determines whether RCM reordering will be done in conjunction with domain decomposition incomplete factorizations.
        // Option to enable (=1) or disable (=0) the Reverse Cuthill\u2013McKee (RCM) algorithm to reorder system equations for smaller bandwidth
        tParameterlist( 0 )( 0 ).set( "AZ_reorder" , 1 );

        // Use preconditioner from a previous Iterate() call
        // Option are AZ_calc, AZ_recalc, AZ_reuse
        tParameterlist( 0 )( 0 ).set( "AZ_pre_calc" , AZ_calc );

        // Determines  whether  matrix  factorization  information will be kept after this solve
        // for example for preconditioner_recalculation
        tParameterlist( 0 )( 0 ).set( "AZ_keep_info" , 0 );

        //--------------------------GMRES specific solver parameters--------------------------------------------------------------------------
        // Set AZ_kspace
        // Krylov subspace size for restarted GMRES
        // Setting mKrylovSpace larger improves the robustness, decreases iteration count, but increases memory consumption.
        // For very difficult problems, set it equal to the maximum number of iterations.
        tParameterlist( 0 )( 0 ).set( "AZ_kspace" ,500 );

        // Set AZ_orthog
        //AZ_classic or AZ_modified
        tParameterlist( 0 )( 0 ).set( "AZ_orthog" , AZ_classic );

        // Set AZ_rthresh
        // Parameter used to modify the relative magnitude of the diagonal entries of the matrix that is used to compute
        // any of the incomplete factorization preconditioners
        tParameterlist( 0 )( 0 ).set( "AZ_rthresh" ,  0.0 );

        // Set AZ_athresh
        // Parameter used to modify the absolute magnitude of the diagonal entries of the matrix that is used to compute
        // any of the incomplete factorization preconditioners
        tParameterlist( 0 )( 0 ).set( "AZ_athresh" ,  0.0 );

        //--------------------------Preconsitioner specific parameters--------------------------------------------------------------------------
        // Determine which preconditioner is used
        // Options are AZ_none, AZ_Jacobi, AZ_sym_GS, AZ_Neumann, AZ_ls, AZ_dom_decomp,
        tParameterlist( 0 )( 0 ).set( "AZ_precond" ,  AZ_dom_decomp );

        // Set preconditioner subdomain solve - direct solve or incomplete
        // Options are AZ_lu, AZ_ilut, , AZ_rilu, AZ_bilu, AZ_icc
        tParameterlist( 0 )( 0 ).set( "AZ_subdomain_solve" ,  AZ_ilut );

        // Set preconditioner polynomial order - polynomial preconditioning, Gauss-Seidel, Jacobi
        tParameterlist( 0 )( 0 ).set( "AZ_poly_ord" ,  3 );

        // Set drop tolerance - for LU, ILUT
        tParameterlist( 0 )( 0 ).set(  "AZ_drop" ,  1.0e-12 );

        // Set level of graph fill in - for ilu(k), icc(k), bilu(k)
        tParameterlist( 0 )( 0 ).set( "AZ_graph_fill" ,  3 );

        // Set ilut fill
        tParameterlist( 0 )( 0 ).set( "AZ_ilut_fill" ,  5.0 );

        // Set Damping or relaxation parameter used for RILU
        tParameterlist( 0 )( 0 ).set( "AZ_omega" ,  1.0 );

        // Set external preconditioner
        tParameterlist( 0 )( 0 ).set( "ifpack_prec_type",  "ILU");
        tParameterlist( 0 )( 0 ).set( "fact: level-of-fill",  3 );

        tParameterlist( 0 )( 0 ).set( "prec_reuse" ,     false );
         */

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set( "NLA_combined_res_jac_assembly", true );
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", 1e-9 );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", 1.00 );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", 20 );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", tDofStrg );

        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", tDofStrg );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "UX", 1 );
        tParameterlist( 0 )( 0 ).set( "UY", 1 );
        if ( tIs3D )
        {
            tParameterlist( 0 )( 0 ).set( "UZ", 1 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        tParameterlist( 0 )( 0 ).set( "Set_Names", tTotalDomainSets + "," + tInterfaces );

        if ( tIs3D )
        {
            tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "UX,UY,UZ,StrainEnergyMaterial1,StrainEnergyMaterial2,VolumeMaterial2,PerimeterMaterial12" ) );
            tParameterlist( 0 )( 0 ).set( "Field_Type", std::string( "NODAL,NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL,GLOBAL" ) );
            tParameterlist( 0 )( 0 ).set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkUZ,IQIBulkStrainEnergy_Material1,IQIBulkStrainEnergy_Material2,"
                                                                    "IQIBulkVolume_Material1,IQIPerimeter_InterfaceMaterial12" ) );
        }
        else
        {
            tParameterlist( 0 )( 0 ).set( "Field_Names", std::string( "UX,UY,,Levelset,StrainEnergyMaterial1,StrainEnergyMaterial2,VolumeMaterial2,PerimeterMaterial12" ) );
            tParameterlist( 0 )( 0 ).set( "Field_Type", std::string( "NODAL,NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL,GLOBAL" ) );
            tParameterlist( 0 )( 0 ).set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQILevelSet,IQIBulkStrainEnergy_Material1,IQIBulkStrainEnergy_Material2,"
                                                                    "IQIBulkVolume_Material1,IQIPerimeter_InterfaceMaterial12" ) );
        }

        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
        tParameterlist( 0 )( 0 ).set( "Time_Offset", 10.0 );
    }

    void
    MORISGENERALParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
    {
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif

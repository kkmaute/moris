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
#include "parameters.hpp"
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
    Vector< uint > tNumElementsPerDir( 2 + tIs3D, 1 );
    Vector< real > tDimensions( 2 + tIs3D, 1.0 );

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
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", tIsOpt );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tLibraryName );
        aParameterLists.set( "restart_file", "" );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::SWEEP );
        aParameterLists.set( "hdf5_path", tHDF5FileName );
        aParameterLists.set( "num_evaluations_per_adv", "1" );
        aParameterLists.set( "finite_difference_type", "all" );
        aParameterLists.set( "finite_difference_epsilons", tFDsweep );
        aParameterLists.set( "finite_difference_adv_indices", tAdvIndicesForFD );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", tNumElementsPerDir );
        aParameterLists.set( "domain_dimensions", tDimensions );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", tLagrangeOrder );
        aParameterLists.set( "lagrange_pattern", tLagrangePattern );

        aParameterLists.set( "bspline_orders", tBsplineOrder );
        aParameterLists.set( "bspline_pattern", "0,1" );

        aParameterLists.set( "initial_refinement", tInitialRef );
        aParameterLists.set( "initial_refinement_pattern", "0,1" );
        aParameterLists.set( "use_advanced_T_matrix_scheme", 1 );

        aParameterLists.set( "lagrange_to_bspline", "0,1" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", tRefineBuffer );
        aParameterLists.set( "staircase_buffer", tRefineBuffer );

        aParameterLists.set( "use_number_aura", 1 );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0,1" );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );

        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "visualize_ghost", tUseGhost );

        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists.set( "IQI_types",
                "IQIBulkStrainEnergy_Material1", "IQIBulkStrainEnergy_Material2", "IQIBulkVolume_Material1", "IQIPerimeter_InterfaceMaterial12" );
        aParameterLists.set( "output_mesh_file", tGENOutputFile );
        aParameterLists.set( "time_offset", 10.0 );

        Matrix< DDUMat > tPhaseMap( 4, 1, 0 );
        tPhaseMap( 0 ) = 0;
        tPhaseMap( 1 ) = 0;
        tPhaseMap( 2 ) = 1;
        tPhaseMap( 3 ) = 2;
        aParameterLists.set( "phase_table", moris::ios::stringify( tPhaseMap ) );

        aParameterLists.set( "print_phase_table", true );

        const real pi = std::acos( -1 );

        if ( tIs3D )
        {
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::SPHERE ) );
            aParameterLists.set( "center_x", 0.5 );
            aParameterLists.set( "center_y", 0.5 );
            aParameterLists.set( "center_z", 0.5 );
            aParameterLists.set( "radius", tHoleRadius );
        }
        else
        {
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
            aParameterLists.set( "center_x", 0.5 );
            aParameterLists.set( "center_y", 0.5 );
            aParameterLists.set( "radius", tHoleRadius );
        }

        // initialize geometry
        if ( tUseBsplineForLevelset )
        {
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            aParameterLists.set( "name", "ADVfield" );
            aParameterLists.set( "center_x", tInterfaceX );
            aParameterLists.set( "center_y", 0.0 );
            aParameterLists.set( "normal_x", -1.0 * std::sin( tPlaneTilde / 180.0 * pi ) );
            aParameterLists.set( "normal_y", -1.0 * std::cos( tPlaneTilde / 180.0 * pi ) );
            aParameterLists.set( "use_multilinear_interpolation", tUseMultiLinear );

            if ( tIsOpt )
            {
                aParameterLists.set( "discretization_mesh_index", 0 );
                aParameterLists.set( "discretization_lower_bound", -2.0 );
                aParameterLists.set( "discretization_upper_bound", 2.0 );
            }

            // Levelset property
            aParameterLists( GEN::PROPERTIES ).add_parameter_list( gen::Field_Type::SCALED_FIELD );

            aParameterLists.set( "name", "LevelsetField" );
            aParameterLists.set( "dependencies", "ADVfield" );
            aParameterLists.set( "scaling_factor", 1.0 );
            aParameterLists.set( "pdv_type", "LS1" );
            aParameterLists.set( "pdv_mesh_set_names", tTotalDomainSets );
        }
        else
        {
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            aParameterLists.set( "center_x", tInterfaceX * 0.5, tInterfaceX, tInterfaceX / 0.5 );
            aParameterLists.set( "center_y", 0.0 );
            aParameterLists.set( "normal_x", -1.0 * std::sin( tPlaneTilde / 180.0 * pi ) );
            aParameterLists.set( "normal_y", -1.0 * std::cos( tPlaneTilde / 180.0 * pi ) );

            // Levelset property
            aParameterLists( GEN::PROPERTIES ).add_parameter_list( gen::Field_Type::CONSTANT );

            aParameterLists.set( "name", "LevelsetField" );
            aParameterLists.set( "constant", 1.0 );
            aParameterLists.set( "pdv_type", "LS1" );
            aParameterLists.set( "pdv_mesh_set_names", tTotalDomainSets );

            if ( tIsOpt )
            {
                aParameterLists.set( "discretization_mesh_index", 0 );
                aParameterLists.set( "discretization_lower_bound", -2.0 );
                aParameterLists.set( "discretization_upper_bound", 2.0 );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // create parameter list for property 1
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungs1" );
        aParameterLists.set( "function_parameters", "2.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungs2" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropBedding" );
        aParameterLists.set( "function_parameters", "1.0e-6" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletU" );
        aParameterLists.set( "function_parameters", tDirichletStr );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletULoad" );
        aParameterLists.set( "function_parameters", tDirichletLoadStr );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 10
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropTraction" );
        aParameterLists.set( "function_parameters", tNeumannStr );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 7
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson1" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson2" );
        aParameterLists.set( "function_parameters", "0.3" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSet" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "tLevelSetFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerLevelSetFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso_Material1" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs1,YoungsModulus;PropPoisson1,PoissonRatio" );

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso_Material2" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs2,YoungsModulus;PropPoisson2,PoissonRatio" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheDirichletBC" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", std::to_string( tNitschePenalty ) );
        aParameterLists.set( "leader_properties", "PropYoungs1,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPNitscheMaterial12Interface" ) );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", std::to_string( tNitschePenalty ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungs1,Material" ) );
        aParameterLists.set( "follower_properties", std::string( "PropYoungs2,Material" ) );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPGhost_Material1" ) );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", std::string( "0.005" ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungs1,Material" ) );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPGhost_Material2" ) );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", std::string( "0.005" ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungs1,Material" ) );

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkU_Material1" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Material1,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tMaterial1Sets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkU_Material2" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Material2,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tMaterial2Sets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichletU" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_properties", "PropDirichletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Material1,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tSupportSSets );

        if ( tPrescribedTraction )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGTraction" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
            aParameterLists.set( "dof_residual", tDofStrg );
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "leader_properties", "PropTraction,Traction" );
            aParameterLists.set( "mesh_set_names", tLoadSSsets );
            }
        else
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGDirichletULoad" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
            aParameterLists.set( "dof_residual", tDofStrg );
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "leader_properties", "PropDirichletULoad,Dirichlet" );
            aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Material1,ElastLinIso" );
            aParameterLists.set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
            aParameterLists.set( "mesh_set_names", tLoadSSsets );
            }

        if ( tUseContact )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGMaterial12Interface" ) );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_CONTACT_UNSYMMETRIC_NITSCHE ) ;
            aParameterLists.set( "dof_residual", tDofStrg );
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "follower_dof_dependencies", tDofStrg );
            aParameterLists.set( "leader_constitutive_models", std::string( "CMStrucLinIso_Material1,ElastLinIso" ) );
            aParameterLists.set( "follower_constitutive_models", std::string( "CMStrucLinIso_Material2,ElastLinIso" ) );
            aParameterLists.set( "stabilization_parameters", std::string( "SPNitscheMaterial12Interface,NitscheInterface" ) );
            aParameterLists.set( "mesh_set_names", tMaterial12DSets );
            }
        else
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGMaterial12Interface" ) );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
            aParameterLists.set( "dof_residual", tDofStrg );
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "follower_dof_dependencies", tDofStrg );
            aParameterLists.set( "leader_constitutive_models", std::string( "CMStrucLinIso_Material1,ElastLinIso" ) );
            aParameterLists.set( "follower_constitutive_models", std::string( "CMStrucLinIso_Material2,ElastLinIso" ) );
            aParameterLists.set( "stabilization_parameters", std::string( "SPNitscheMaterial12Interface,NitscheInterface" ) );
            aParameterLists.set( "mesh_set_names", tMaterial12DSets );
            }

        if ( tUseGhost )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGGhostMaterial1" ) );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual", tDofStrg );
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "follower_dof_dependencies", tDofStrg );
            aParameterLists.set( "stabilization_parameters", std::string( "SPGhost_Material1,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tMaterial1Ghost );

            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGGhostMaterial2" ) );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual", tDofStrg );
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "follower_dof_dependencies", tDofStrg );
            aParameterLists.set( "stabilization_parameters", std::string( "SPGhost_Material2,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tMaterial2Ghost );
            }

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkUX" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomainSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkUY" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "dof_quantity", tDofStrg );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tTotalDomainSets );

        if ( tIs3D )
        {
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "IQIBulkUZ" );
            aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "dof_quantity", tDofStrg );
            aParameterLists.set( "vectorial_field_index", 2 );
            aParameterLists.set( "mesh_set_names", tTotalDomainSets );
            }

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQILevelSet" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::PROPERTY ) ;
        aParameterLists.set( "leader_properties", "PropLevelSet,Property" );
        aParameterLists.set( "mesh_set_names", tTotalDomainSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy_Material1" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Material1,Elast" );
        aParameterLists.set( "mesh_set_names", tMaterial1Sets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy_Material2" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Material2,Elast" );
        aParameterLists.set( "mesh_set_names", tMaterial2Sets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume_Material1" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists.set( "leader_properties", "PropLevelSet,Density" );
        aParameterLists.set( "mesh_set_names", tMaterial1Sets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume_Material2" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists.set( "leader_properties", "PropLevelSet,Density" );
        aParameterLists.set( "mesh_set_names", tMaterial2Sets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIPerimeter_InterfaceMaterial12" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "mesh_set_names", tMaterial12SSets );

        // create computation  parameter list
        aParameterLists( FEM::COMPUTATION );
        aParameterLists.set( "print_physics_model", false );

        aParameterLists.set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        aParameterLists.set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        ///*
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

#ifdef MORIS_USE_MUMPS
        aParameterLists.set( "Solver_Type", "Amesos_Mumps" );
#else
        aParameterLists.set( "Solver_Type", "Amesos_Superludist" );
#endif

        //*/

        /*
        aParameterLists( 0 ).push_back( add_parameter_list( sol::SolverType::BELOS_IMPL ) );

        // Solver type: GMRES, Flexible GMRES, Block CG , PseudoBlockCG, Stochastic CG, Recycling GMRES, Recycling CG, MINRES, LSQR, TFQMR
        //              Pseudoblock TFQMR, Seed GMRES, Seed CG
        aParameterLists.set( "Solver Type" ,  "GMRES" );

        // Diagnostics: Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails
        sint tVerbosity = Belos::Errors; // + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails;
        aParameterLists.set( "Verbosity" , tVerbosity );

        // Maximum number of blocks in Krylov factorization
        aParameterLists.set( "Num Blocks", 500   );

        // Block size to be used by iterative solver
        aParameterLists.set( "Block Size", 1   );

        // Allowable Belos solver iterations
        aParameterLists.set( "Maximum Iterations" , 500 );

        // Allowable Belos solver iterations
        //aParameterLists.set( "Maximum Restarts" ,  );

        // Convergence criteria
        aParameterLists.set( "Convergence Tolerance" ,  1e-12 );

        // Preconditioner
        //aParameterLists.set( "ifpack_prec_type",  "ILU");
        //aParameterLists.set( "fact: level-of-fill",  3 );

        //aParameterLists.set( "ifpack_prec_type",  "ILUT");
        //aParameterLists.set( "fact: ilut level-of-fill", 1.0 );
        //aParameterLists.set( "fact: drop tolerance", 1e-1 );

        //aParameterLists.set( "ifpack_prec_type",  "Amesos");
        //aParameterLists.set( "amesos: solver type", "Amesos_Umfpack");

        aParameterLists.set( "ml_prec_type",  "NSSA");
         */

        /*
        aParameterLists( 0 ).push_back( add_parameter_list( sol::SolverType::AZTEC_IMPL ) );

        //options are: AZ_gmres, AZ_gmres_condnum, AZ_cg, AZ_cg_condnum, AZ_cgs, AZ_tfqmr, AZ_bicgstab
        aParameterLists.set( "AZ_solver" ,  AZ_gmres );

            // Allowable Aztec solver iterations
        aParameterLists.set( "AZ_max_iter", 500   );

            // Allowable Aztec iterative residual
        aParameterLists.set( "rel_residual" , 1e-08 );

        // set Az_conv -convergence criteria
        // options are AZ_r0, AZ_rhs, AZ_Anorm, AZ_noscaled, AZ_sol
        aParameterLists.set( "AZ_conv" ,  AZ_r0 );

        // set Az_diagnostic parameters
        // Set whether or not diagnostics for every linear iteration are printed or not. options are AZ_all, AZ_none
        aParameterLists.set( "AZ_diagnostics" ,  AZ_all );

        // set AZ_output options
        // options are AZ_all, AZ_none, AZ_warnings, AZ_last, AZ_summary
        aParameterLists.set( "AZ_output" ,  AZ_all );

        // Determines the submatrices factored with the domain decomposition algorithms
        // Option to specify with how many rows from other processors each processor\u2019s local submatrix is augmented.
        aParameterLists.set( "AZ_overlap" , 1 );

        // Determines how overlapping subdomain results are combined when different processors have computed different values for the same unknown.
        // Options are AZ_standard, AZ_symmetric
        aParameterLists.set( "AZ_type_overlap" , AZ_standard );

        // Determines whether RCM reordering will be done in conjunction with domain decomposition incomplete factorizations.
        // Option to enable (=1) or disable (=0) the Reverse Cuthill\u2013McKee (RCM) algorithm to reorder system equations for smaller bandwidth
        aParameterLists.set( "AZ_reorder" , 1 );

        // Use preconditioner from a previous Iterate() call
        // Option are AZ_calc, AZ_recalc, AZ_reuse
        aParameterLists.set( "AZ_pre_calc" , AZ_calc );

        // Determines  whether  matrix  factorization  information will be kept after this solve
        // for example for preconditioner_recalculation
        aParameterLists.set( "AZ_keep_info" , 0 );

        //--------------------------GMRES specific solver parameters--------------------------------------------------------------------------
        // Set AZ_kspace
        // Krylov subspace size for restarted GMRES
        // Setting mKrylovSpace larger improves the robustness, decreases iteration count, but increases memory consumption.
        // For very difficult problems, set it equal to the maximum number of iterations.
        aParameterLists.set( "AZ_kspace" ,500 );

        // Set AZ_orthog
        //AZ_classic or AZ_modified
        aParameterLists.set( "AZ_orthog" , AZ_classic );

        // Set AZ_rthresh
        // Parameter used to modify the relative magnitude of the diagonal entries of the matrix that is used to compute
        // any of the incomplete factorization preconditioners
        aParameterLists.set( "AZ_rthresh" ,  0.0 );

        // Set AZ_athresh
        // Parameter used to modify the absolute magnitude of the diagonal entries of the matrix that is used to compute
        // any of the incomplete factorization preconditioners
        aParameterLists.set( "AZ_athresh" ,  0.0 );

        //--------------------------Preconsitioner specific parameters--------------------------------------------------------------------------
        // Determine which preconditioner is used
        // Options are AZ_none, AZ_Jacobi, AZ_sym_GS, AZ_Neumann, AZ_ls, AZ_dom_decomp,
        aParameterLists.set( "AZ_precond" ,  AZ_dom_decomp );

        // Set preconditioner subdomain solve - direct solve or incomplete
        // Options are AZ_lu, AZ_ilut, , AZ_rilu, AZ_bilu, AZ_icc
        aParameterLists.set( "AZ_subdomain_solve" ,  AZ_ilut );

        // Set preconditioner polynomial order - polynomial preconditioning, Gauss-Seidel, Jacobi
        aParameterLists.set( "AZ_poly_ord" ,  3 );

        // Set drop tolerance - for LU, ILUT
        aParameterLists.set(  "AZ_drop" ,  1.0e-12 );

        // Set level of graph fill in - for ilu(k), icc(k), bilu(k)
        aParameterLists.set( "AZ_graph_fill" ,  3 );

        // Set ilut fill
        aParameterLists.set( "AZ_ilut_fill" ,  5.0 );

        // Set Damping or relaxation parameter used for RILU
        aParameterLists.set( "AZ_omega" ,  1.0 );

        // Set external preconditioner
        aParameterLists.set( "ifpack_prec_type",  "ILU");
        aParameterLists.set( "fact: level-of-fill",  3 );

        aParameterLists.set( "prec_reuse" ,     false );
         */

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_combined_res_jac_assembly", true );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1e-9 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.00 );
        aParameterLists.set( "NLA_max_iter", 20 );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", tDofStrg );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", tDofStrg );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "UX", 1 );
        aParameterLists.set( "UY", 1 );
        if ( tIs3D )
        {
            aParameterLists.set( "UZ", 1 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", tTotalDomainSets + "," + tInterfaces );

        if ( tIs3D )
        {
            aParameterLists.set( "Field_Names", std::string( "UX,UY,UZ,StrainEnergyMaterial1,StrainEnergyMaterial2,VolumeMaterial2,PerimeterMaterial12" ) );
            aParameterLists.set( "Field_Type", std::string( "NODAL,NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL,GLOBAL" ) );
            aParameterLists.set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkUZ,IQIBulkStrainEnergy_Material1,IQIBulkStrainEnergy_Material2,"
                                                                    "IQIBulkVolume_Material1,IQIPerimeter_InterfaceMaterial12" ) );
        }
        else
        {
            aParameterLists.set( "Field_Names", std::string( "UX,UY,,Levelset,StrainEnergyMaterial1,StrainEnergyMaterial2,VolumeMaterial2,PerimeterMaterial12" ) );
            aParameterLists.set( "Field_Type", std::string( "NODAL,NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL,GLOBAL" ) );
            aParameterLists.set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQILevelSet,IQIBulkStrainEnergy_Material1,IQIBulkStrainEnergy_Material2,"
                                                                    "IQIBulkVolume_Material1,IQIPerimeter_InterfaceMaterial12" ) );
        }

        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 10.0 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif

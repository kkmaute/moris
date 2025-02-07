/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Level_Set_Beam_SIMP_Hole_Seeding.cpp
 *
 */

#include <string>
#include <iostream>
#include <sstream>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "parameters.hpp"
#include "fn_PRM_MORIS_GENERAL_Parameters.hpp"
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"
#include "fn_stringify_matrix.hpp"

#include "AztecOO.h"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------------

namespace moris
{
    std::string
    moris_to_string( real tValue )
    {
        std::ostringstream streamObj;

        // Set precision
        streamObj << std::scientific;
        streamObj << std::setprecision( 15 );

        // Add value to stream
        streamObj << tValue;

        // Get string from output string stream
        return streamObj.str();
    }

    std::string tName = "Level_Set_Beam_SIMP_Hole_Seeding";

    std::string tLibraryName   = tName + ".so";
    std::string tGENOutputFile = "GEN_" + tName + ".exo";

    std::string tOutputFileName = tName + ".exo";

    bool tIs3D             = false;
    bool tIsOpt            = true;
    bool tUseGhost         = true;
    bool tUseAbsoluteValue = true;

    bool tUseDensityShift = true;

    static real tDensityShift          = 0.0;
    uint        tDensityShiftStart     = 20;
    uint        tDensityShiftIntervall = 100;

    static uint tItarationCounter = 0;

    //-------------------------------
    // Opt constant_parameters

    real tMMAPenalty  = 5.0;
    real tMMAStepSize = 0.021;
    int  tMMAMaxIter  = 20;

    real tInitialStrainEnergy = 1437.01;
    real tMaxMass             = 400;

    real tPerimeterPenalty = 0.02;
    real tInitialPerimeter = 60.0;

    real tRegularizationPenalty = 0.02;    // 0.02
    real tInitialRegularization = 60;

    //-------------------------------

    real phi_sh = 0.5;
    real phi_rt = 3.0;

    int tDispOrder = 1;

    real tBsplineLimitTop    = 1.0;
    real tBsplineLimitBottom = 0.0;

    real tElementEdgeLength = 1.0;

    real tInitialDensity = 0.0001;

    moris::real tBSplineLimit = 1.0;
    moris::real tPhiBandwidth = 3.0 * tElementEdgeLength;
    moris::real tPhiGradient  = tBSplineLimit * std::log( 199.0 ) / ( 2.0 * tPhiBandwidth );
    moris::real tPhiGamma     = 2.0 * std::log( 10.0 ) / std::pow( 2.0 / ( std::exp( -2.0 * tPhiBandwidth * tPhiGradient / tBSplineLimit ) + 1.0 ) - 1.0, 2.0 );

    /* ------------------------------------------------------------------------ */
    // material parameters

    // conductivity
    std::string tConductivity = "100.0";

    // capacity
    std::string tCapacityTheta = "0.001";
    std::string tCapacityPhi   = "0.0";

    // density
    std::string tDensityTheta = "1.0";
    std::string tDensityPhi   = "0.0";

    // prescribed theta on interface
    std::string tPrescTheta = "1.0";

    // prescribed phi on interface
    std::string tPrescPhi = "0.0";

    //------------------------------------------------------------------------------

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }
    //-------------------------------
    std::string tBulkSets = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tVoidSets = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    std::string tDirichletSets = "SideSet_4_n_p1,SideSet_4_n_p1";
    std::string tLoadSets      = "SideSet_2_n_p1,SideSet_2_n_p1";

    std::string tInterfaceVoidSets = "iside_b0_1_b1_0";

    std::string tVoidInterfaceSets = "iside_b0_0_b1_1";

    std::string tInteriorGhost = "ghost_p1";

    std::string tTotalDomain = tBulkSets + "," + tVoidSets;

    std::string tTotalDomainAGhost = tUseGhost ? tTotalDomain + "," + tInteriorGhost : tTotalDomain;

    void
    Func_Neumann_U(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        if ( aFIManager->get_IG_geometry_interpolator()->valx()( 1 ) < 0.2 )
        {
            if ( tIs3D )
            {
                aPropMatrix = { { 0.0 }, { aParameters( 0 )( 0 ) }, { 0.0 } };
            }
            else
            {
                aPropMatrix = { { 0.0 }, { aParameters( 0 )( 0 ) } };
            }
        }
        else
        {
            if ( tIs3D )
            {
                aPropMatrix = { { 0.0 }, { 0.0 }, { 0.0 } };
            }
            else
            {
                aPropMatrix = { { 0.0 }, { 0.0 } };
            }
        }
    }

    real
    Const_Geometry(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aGeometryParameters )
    {
        return 0.8;
    }

    //------------------------------------------------------------------------------
    void
    compute_density_shift()
    {
        if ( tUseDensityShift )
        {
            if ( tItarationCounter >= tDensityShiftStart && tItarationCounter < ( tDensityShiftStart + tDensityShiftIntervall ) )
            {
                real tPow     = std::pow( ( ( (real)tItarationCounter - (real)tDensityShiftStart ) / ( (real)tDensityShiftIntervall ) ), 2.0 );
                tDensityShift = tInitialDensity + ( 1 - tInitialDensity ) * tPow;
            }
            else if ( tItarationCounter < tDensityShiftStart )
            {
                tDensityShift = 0.0;
            }
            else
            {
                tDensityShift = 1.0;
            }
        }
    }

    //------------------------------------------------------------------------------
    void
    tYoungsFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        const moris::Matrix< moris::DDRMat >& tHCT  = aParameters( 0 );
        moris::real                           tBeta = aParameters( 0 )( 1 );

        real tLevelSet = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        real tDensity = ( tLevelSet - phi_sh ) / ( 1 - phi_sh );

        if ( tDensity < 0 ) { tDensity = 0.0001; }
        if ( tDensity > 1 ) { tDensity = 1; }

        tDensity = tDensityShift + ( 1 - tDensityShift ) * tDensity;

        aPropMatrix = tHCT * std::pow( tDensity, tBeta );
    }

    //------------------------------------------------------------------------------
    void
    tDerYoungsFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        MORIS_ERROR( false, "Do not need this one" );
        const moris::Matrix< moris::DDRMat >& tHCT  = aParameters( 0 );
        moris::real                           tBeta = aParameters( 0 )( 1 );

        real tLevelSet = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        real tDensity = ( tLevelSet - phi_sh ) / ( 1 - phi_sh );

        if ( tDensity < 0 ) { tDensity = 0.0001; }
        if ( tDensity > 1 ) { tDensity = 1; }

        // FIXME density shift missing

        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->N() * tBeta * tHCT( 0 ) * std::pow( tDensity, tBeta - 1 ) / ( 1 - phi_sh );
    }

    //------------------------------------------------------------------------------

    void
    tDensityFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        real tLevelSet = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        real tDensity = ( tLevelSet - phi_sh ) / ( 1 - phi_sh );

        tDensity = tDensityShift + ( 1 - tDensityShift ) * tDensity;

        if ( tDensity < 0 ) { tDensity = 0.0001; }
        if ( tDensity > 1 ) { tDensity = 1; }

        aPropMatrix.set_size( 1, 1, tDensity );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    tDerDensityFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        MORIS_ERROR( false, "Do not need this one" );

        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->N() / ( 1.0 - phi_sh );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    tLevelSetFuncReal(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // return absolute value of level set function
        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val();
    }

    /* ------------------------------------------------------------------------ */

    // Level set function defining property in FEM
    void
    tLevelSetFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // Matrix< DDRMat > tCoords = aFIManager->get_IG_geometry_interpolator()->valx();

        // bool tBool = is_in_vicinity_of_load( tCoords );

        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        value = value - phi_sh;

        // return PDV derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        // return absolute value of level set function
        aPropMatrix = factor * value;
    }

    /* ------------------------------------------------------------------------ */

    // Derivative of level set function with respect to PDV
    void
    tDerLevelSetFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        value = value - phi_sh;

        // return PDV derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->N();
    }

    /* ------------------------------------------------------------------------ */

    // Spatial derivative of level set function defining property in FEM
    void
    tLevelSetGradxFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        value = value - phi_sh;

        // return spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->gradx( 1 );
    }

    /* ------------------------------------------------------------------------ */

    // Derivative of spatial derivative of level set function with respect to PDV
    void
    tDerLevelSetGradxFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        value = value - phi_sh;

        // return PDV derivative of spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsoluteValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->dnNdxn( 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );

        return tConstraintTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );

        real obj1 = aCriteria( 0 ) / tInitialStrainEnergy;
        // real obj2 = aCriteria( 1 ) / tInitialStrainEnergy;
        real obj2 = tPerimeterPenalty * aCriteria( 2 ) / tInitialPerimeter;
        real obj3 = tRegularizationPenalty * aCriteria( 3 ) / tInitialRegularization;
        // real obj5 = tRegularizationPenalty* aCriteria( 5 ) ;

        tObjectives( 0, 0 ) = obj1 + obj2 + obj3;    // + obj2 + obj3 + obj4 + obj5;

        std::cout << "% --------------------------------- % \n";
        std::cout << "Objective                = " << tObjectives( 0, 0 ) << " \n";
        std::cout << "Strain Energy            = " << aCriteria( 0 ) << " ( " << obj1 / tObjectives( 0, 0 ) << " )\n";
        // std::cout << "Strain Energy (Interior) = " << aCriteria( 1 )      << " ( " <<  obj2 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Perimeter                = " << aCriteria( 2 ) << " ( " << obj2 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Regularization            = " << aCriteria( 3 ) << " ( " << obj3 / tObjectives( 0, 0 ) << " )\n";
        // std::cout << "H1Error PDV              = " << aCriteria( 5 )      << " ( " <<  aCriteria( 1 ) / tObjectives( 0, 0 ) << " )\n";
        std::cout << " \n";

        std::cout << "min ADV                  = " << aADVs.min() << " \n";
        std::cout << "max ADV                  = " << aADVs.max() << " \n"
                  << std::flush;

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, 1 );
        tConstraints( 0 ) = aCriteria( 1 ) / tMaxMass - 1.0;

        std::cout << "Volume     = " << aCriteria( 1 ) << " \n";
        std::cout << "Constraint = " << tConstraints( 0 ) << " \n";
        std::cout << "% --------------------------------- % \n"
                  << std::flush;

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

        tDObjectiveDCriteria( 0 ) = 1.0 / tInitialStrainEnergy;
        // tDObjectiveDCriteria( 1 ) = 1.0 / tInitialStrainEnergy;
        tDObjectiveDCriteria( 2 ) = tPerimeterPenalty / tInitialPerimeter;
        tDObjectiveDCriteria( 3 ) = tRegularizationPenalty / tInitialRegularization;
        // tDObjectiveDCriteria( 5 ) = tRegularizationPenalty* 1.0;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria(
            const Vector< real >& aADVs,
            const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( 1, aCriteria.size(), 0.0 );

        tDConstraintDCriteria( 1 ) = 1.0 / tMaxMass;

        return tDConstraintDCriteria;
    }

    //------------------------------------------------------------------------------

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        compute_density_shift();

        tItarationCounter++;

        return true;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", tIsOpt );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tLibraryName );
        aParameterLists.set( "restart_file", "" );
        aParameterLists.set( "reinitialize_interface_iter", 10000 );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::GCMMA );
        aParameterLists.set( "step_size", tMMAStepSize );
        aParameterLists.set( "penalty", tMMAPenalty );
        aParameterLists.set( "max_its", tMMAMaxIter );    // Maximum number of iterations
        aParameterLists.set( "restart_index", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", 60, 20 );
        aParameterLists.set( "domain_dimensions", 60.0, 20.0 );
        aParameterLists.set( "domain_offset", "0.0,0.0" );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "2,1" );
        aParameterLists.set( "lagrange_pattern", "0,1" );

        aParameterLists.set( "bspline_orders", "1,2" );
        aParameterLists.set( "bspline_pattern", "0,1" );

        aParameterLists.set( "initial_refinement", "1,0" );
        aParameterLists.set( "initial_refinement_pattern", "0,1" );

        aParameterLists.set( "lagrange_to_bspline", "0,1;-1" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", 1 );
        aParameterLists.set( "staircase_buffer", 1 );

        aParameterLists.set( "use_number_aura", 1 );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );

        // aParameterLists.set( "lagrange_mesh_output_file_name", "HMRLagrangeMesh.vtk" );

        aParameterLists.set( "use_refine_low_level_elements", false );
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
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists.set( "IQI_types", "IQIBulkStrainEnergy", "IQIBulkVolume", "IQIPerimeter_InterfaceVoid", "IQIHeatMethodPenalty" );
        // aParameterLists.set("output_mesh_file", tGENOutputFile );
        aParameterLists.set( "time_offset", 10.0 );
        // FIXME     this has to change
        // aParameterLists.set("PDV_types"         , "LS1");

        // initialize fins as swiss cheese geometry
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "isocontour_tolerance", 10e-14 );
        aParameterLists.set( "isocontour_threshold", 0.5 );
        aParameterLists.set( "field_function_name", "Const_Geometry" );
        aParameterLists.set( "name", "Level_Set_Field" );
        aParameterLists.set( "number_of_refinements", 0, 0 );
        aParameterLists.set( "refinement_mesh_index", 0, 1 );
        aParameterLists.set( "use_multilinear_interpolation", false );

        if ( tIsOpt )
        {
            aParameterLists.set( "discretization_mesh_index", 1 );
            aParameterLists.set( "discretization_lower_bound", tBsplineLimitBottom );
            aParameterLists.set( "discretization_upper_bound", tBsplineLimitTop );
        }

        uint tParamCounter = 0;
        aParameterLists( GEN::PROPERTIES ).add_parameter_list( gen::Field_Type::SCALED_FIELD );
        aParameterLists.set( "name", "LvL_Set_Field" );
        aParameterLists.set( "dependencies", "Level_Set_Field" );
        aParameterLists.set( "scaling_factor", 1.0 );
        aParameterLists.set( "pdv_type", "LS1" );
        // aParameterLists.set("discretization_mesh_index",   -1);
        // aParameterLists.set("discretization_lower_bound", 0.001);
        // aParameterLists.set("discretization_upper_bound", 1.0);
        aParameterLists.set( "pdv_mesh_set_names", tTotalDomainAGhost );
        tParamCounter++;

        aParameterLists( GEN::PROPERTIES ).add_parameter_list( gen::Field_Type::SCALED_FIELD );
        aParameterLists.set( "name", "Density_Field" );
        aParameterLists.set( "dependencies", "Level_Set_Field" );
        aParameterLists.set( "scaling_factor", 1.0 );
        aParameterLists.set( "pdv_type", "DENSITY" );
        aParameterLists.set( "pdv_mesh_set_names", tTotalDomainAGhost );
        tParamCounter++;
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
        aParameterLists.set( "value_function", "tDensityFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerDensityFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", "1.0,3.0" );
        aParameterLists.set( "value_function", "tYoungsFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerYoungsFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropBedding" );
        aParameterLists.set( "function_parameters", "1.0e-8" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletU" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 10
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropTraction" );
        aParameterLists.set( "function_parameters", "10.0" );
        aParameterLists.set( "value_function", "Func_Neumann_U" );

        // create parameter list for property 7
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // common properties for theta and phi problems

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", tConductivity );
        aParameterLists.set( "value_function", "Func_Const" );

        // properties for Theta
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensityTheta" );
        aParameterLists.set( "function_parameters", tDensityTheta );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacityTheta" );
        aParameterLists.set( "function_parameters", tCapacityTheta );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPrescTheta" );
        aParameterLists.set( "function_parameters", tPrescTheta );
        aParameterLists.set( "value_function", "Func_Const" );

        // properties for phi problem

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensityPhi" );
        aParameterLists.set( "function_parameters", tDensityPhi );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacityPhi" );
        aParameterLists.set( "function_parameters", tCapacityPhi );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPrescPhi" );
        aParameterLists.set( "function_parameters", tPrescPhi );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropEigenStrainPhi" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSetConst" );
        aParameterLists.set( "value_function", "Func_Const" );
        aParameterLists.set( "function_parameters", "1.0" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSetGradxConst" );
        aParameterLists.set( "value_function", "Func_Const" );
        aParameterLists.set( "function_parameters", "1.0;1.0" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSet" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "tLevelSetFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerLevelSetFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSetReal" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "tLevelSetFuncReal" );
        aParameterLists.set( "dv_derivative_functions", "tDerLevelSetFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLevelSetGradx" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "tLevelSetGradxFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerLevelSetGradxFunc" );
        aParameterLists.set( "dv_dependencies", "LS1" );

        // time continuity weights
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightCurrent" );
        aParameterLists.set( "function_parameters", "10.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWeightPrevious" );
        aParameterLists.set( "function_parameters", "10.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // initial condition
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropInitialCondition" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        // create parameter list for constitutive model - Theta problem
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusionTheta" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );

        // create parameter list for constitutive model - Phi problem
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusionPhi" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        aParameterLists.set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );

        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheDirichletBC" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPGhost" ) );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", std::string( "0.005" ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungs,Material" ) );
        aParameterLists.set( "follower_properties", std::string( "PropYoungs,Material" ) );

        // create parameter list for ghost stabilization parameter for theta and phi problems
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPTemp" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        // create parameter list for DBC on interface for theta problem
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // create IWG  - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionThetaBulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceInnerTheta" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tVoidInterfaceSets );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceOuterTheta" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tInterfaceVoidSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTimeContinuityTheta" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );
        aParameterLists.set( "time_continuity", true );

        // create IWG - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionOuterBulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceInnerPhi" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tVoidInterfaceSets );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceOuterPhi" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tInterfaceVoidSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichletU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropDirichletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tLoadSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTraction" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropTraction,Traction" );
        aParameterLists.set( "mesh_set_names", tDirichletSets );

        if ( tUseGhost )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGGhost" ) );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "UX,UY" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY;THETA;PHID" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY;THETA;PHID" );
            // aParameterLists.set( "leader_dv_dependencies",     "LS1") ;
            // aParameterLists.set( "follower_dv_dependencies",      "LS1") ;
            aParameterLists.set( "stabilization_parameters", std::string( "SPGhost,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tInteriorGhost );
        }

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkUX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkUY" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,Elast" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIPerimeter_InterfaceVoid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "mesh_set_names", tInterfaceVoidSets );

        // smooth fuction
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQILevelSet" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropLevelSetReal,Property" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQILevelSetHeatMethod" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropLevelSet,Property" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTHETA" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", "THETA" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // Nodal PHID IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkPHID" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "PHID" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // H1 Error if reference is constant
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIH1ErrorConst" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::H1_ERROR );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "leader_properties", "PropLevelSetConst,L2_Reference;PropLevelSetGradxConst,H1S_Reference" );
        aParameterLists.set( "function_parameters", "1.0 / 1.0 / 1.0" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // H1 Error if reference is design dependent
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIH1Error" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::H1_ERROR );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "leader_properties", "PropLevelSet,L2_Reference;PropLevelSetGradx,H1S_Reference" );
        aParameterLists.set( "function_parameters", "1.0 / 1.0 / 1.0" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // H1 Error if reference is design dependent
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIHeatMethodPenalty" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::HEAT_METHOD_PENALTY );
        aParameterLists.set( "dof_quantity", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "leader_properties", "PropLevelSet,L2_Reference;PropLevelSetGradx,H1S_Reference" );
        aParameterLists.set( "function_parameters", moris_to_string( tBSplineLimit ) + " / " + moris_to_string( tPhiGradient ) + " / " + moris_to_string( tPhiGamma ) + " / 0.1 / 1.0 / 1.0 / 1.0 " );
        aParameterLists.set( "mesh_set_names", tVoidSets );

        // create computation  parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

#ifdef MORIS_USE_MUMPS
        aParameterLists.set( "Solver_Type", "Amesos_Mumps" );
#else
        aParameterLists.set( "Solver_Type", "Amesos_Superludist" );
#endif

        /*
        // Solver type: GMRES, Flexible GMRES, Block CG , PseudoBlockCG, Stochastic CG, Recycling GMRES, Recycling CG, MINRES, LSQR, TFQMR
        //              Pseudoblock TFQMR, Seed GMRES, Seed CG
        aParameterLists.set( "Solver Type" ,  "GMRES" );

        // Diagnostics: Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails
        sint tVerbosity = Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails;
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
        aParameterLists.set( "ifpack_prec_type",  "ILU");
        aParameterLists.set( "fact: level-of-fill", 10  );

        //aParameterLists.set( "ifpack_prec_type",  "ILUT");
        //aParameterLists.set( "fact: ilut level-of-fill", 15.0 );
        //aParameterLists.set( "fact: drop tolerance", 1e-12 );
        */

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 5 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 5 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 1
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists.set( "NLA_max_iter", 5 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-7 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 5 );

        //------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();    // nonlinear solver index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );    // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_DofTypes", "THETA" );
        aParameterLists.set( "NLA_Secondary_DofTypes", "" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();    // nonlinear solver index 1
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "1" );    // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_DofTypes", "PHID" );
        aParameterLists.set( "NLA_Secondary_DofTypes", "" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "3" );    // set nonlinear algorithm with index 0
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );
        aParameterLists.set( "NLA_Secondary_DofTypes", "" );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();    // nonlinear solver index 2
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "2" );    // set nonlinear algorithm with index 1.
        aParameterLists.set( "NLA_Sub_Nonlinear_Solver", "0,1,2" );       // set sub nonlinear solvers with index 0 and 1
        aParameterLists.set( "NLA_DofTypes", "THETA;PHID;UX,UY" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Nonlinear_Solver", 3 );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "THETA;PHID;UX,UY" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_save_operator_to_matlab", "Mat.dat" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "UX", 0 );
        aParameterLists.set( "UY", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tTotalDomain );

        aParameterLists.set( "Field_Names", std::string( "UX,UY,VOL,PHID,THETA,LVLSET,HEATMETHOD,LEVELSETHEAT" ) );
        aParameterLists.set( "Field_Type", std::string( "NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL,NODAL" ) );
        aParameterLists.set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkVolume,IQIBulkPHID,IQIBulkTHETA,IQILevelSet,IQIHeatMethodPenalty,IQILevelSetHeatMethod" ) );

        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 10.0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( Parameter_List( "" ) );
        prm::create_remeshing_parameterlist( aParameterLists( 0 )( 0 ) );
        aParameterLists.set( "mode", "ab_initio" );
        aParameterLists.set( "remeshing_field_names", "Box,Level_Set_Field" );
        aParameterLists.set( "remeshing_levels_of_refinement", "0,0;0,0" );
        aParameterLists.set( "remeshing_refinement_pattern", "0,1;0,1" );
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif

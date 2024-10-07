/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Levelset_Boxbeam_Restart.cpp
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
#include "fn_PRM_MORIS_GENERAL_Parameters.hpp"
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
    //------------------------------------------------------------------------------

    // Main problem parameters

    std::string tName = "Levelset_Boxbeam_Restart";

    bool tIs3D     = false;
    bool tIsOpt    = true;
    bool tUseGhost = true;

    bool tUseAbsolutValue = true;

    // wall thickness
    real tWallThickness = 0.05;

    // background mesh parameters
    std::string tNumElementsPerDir = tIs3D ? "45,15,15" : "30,10";
    std::string tDimensions        = tIs3D ? "3,1,1" : "3,1";
    std::string tOffSet            = tIs3D ? "0.0,0.0,0.0" : "0.0,0.0";
    std::string tSideSets          = tIs3D ? "1,2,3,4,5,6" : "1,2,3,4";

    int tDispOrder = 1;

    // Hole Seeding parameters
    sint tNumHolesY = 3;
    sint tNumHolesZ = 1 * tNumHolesY;
    sint tNumHolesX = 2 * tNumHolesY + 1;

    real tHoleRadius   = tIs3D ? 0.2771 / tNumHolesY : 0.2771 / tNumHolesY;
    real tHoleExponent = 6.0;
    real tHoleScaling  = 1.0;

    real tHoleXdim = 3.0;
    real tHoleYdim = 1.0;
    real tHoleZdim = 1.0;

    real tHoleXOrg = 0.0;
    real tHoleYOrg = 0.0;
    real tHoleZOrg = 0.0;

    bool tHoleOffsetRow = false;

    // Nitsche penalty
    std::string tNitschePenalty = "100.0";

    // optimization parameters
    real tInitialStrainEnergy = tIs3D ? 1.49721 + 1.47892 : 2.04222e+00 + 2.36298e+00;
    real tInitialPerimeter    = tIs3D ? 27.3111 : 14.0771;

    real tInitialRegularization = 0.243945;
    real tPerimeterPenalty      = 0.05;
    real tRegularizationPenalty = 0.01;

    // prescribed gradient when using non-pdv level set field
    std::string tLevelSetGradxConstant = tIs3D ? "1.0;1.0;1.0" : "1.0;1.0";

    real tMaxMass = tIs3D ? 1.5 : 1.0;

    real tMMAPenalty  = 5.0;
    real tMMAStepSize = 0.05;
    int  tMMAMaxIter  = 5;

    real tBsplineLimit = tHoleRadius;

    // other mesh depedendent parameters
    real tElementEdgeLength = 1.0 / tNumHolesX / pow( 2, 3 );
    real tLoadLimitY        = std::floor( 0.2 / tElementEdgeLength ) * tElementEdgeLength;

    //------------------------------------------------------------------------------
    // Derived problem paramters

    std::string tOutputFileName = tName + ".exo";
    std::string tLibraryName    = "Levelset_Boxbeam_Restart.so";
    std::string tGENOutputFile  = "GEN_" + tName + ".exo";

    std::string tFrameSets    = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tInteriorSets = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tVoidSets     = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    std::string tTotalDomainSets = tFrameSets + "," + tInteriorSets;

    std::string tFrameGhost    = "ghost_p2";
    std::string tInteriorGhost = "ghost_p1";

    std::string tFrameLoadSSets    = "SideSet_2_n_p2,SideSet_2_c_p2";
    std::string tFrameSupportSSets = "SideSet_4_n_p2,SideSet_4_c_p2";
    std::string tFrameFreeSSets =
            std::string( "SideSet_1_n_p2,SideSet_1_c_p2" ) + std::string( "SideSet_3_n_p2,SideSet_3_c_p2" );

    std::string tInterfaceVoidSSets  = "iside_b0_2_b1_0,iside_b0_1_b1_0";
    std::string tInterfaceVoidSSets2 = "iside_b0_0_b1_1,iside_b0_0_b1_2";

    std::string tInterfaceVoidSSets3 = "iside_b0_1_b1_0";
    std::string tInterfaceVoidSSets4 = "iside_b0_0_b1_1";

    std::string tFrameInteriorDSets = "dbl_iside_p0_1_p1_2";

    std::string tDofStrg = tIs3D ? "UX,UY,UZ" : "UX,UY";

    std::string tHeatDofs = "THETA;PHID";

    std::string tDofStrgAll = tIs3D ? "UX,UY,UZ;THETA;PHID" : "UX,UY;THETA;PHID";

    std::string tDirichletStr = tIs3D ? "0.0;0.0;0.0" : "0.0;0.0";

    std::string tTotalDomain = tInteriorSets + "," + tVoidSets;

    std::string tTotalDomain1 = tFrameSets + "," + tInteriorSets + "," + tVoidSets;

    /* ------------------------------------------------------------------------ */
    // material parameters

    // conductivity
    std::string tConductivity = "1.0";

    // capacity
    std::string tCapacityTheta = "0.1";
    std::string tCapacityPhi   = "0.0";

    // density
    std::string tDensityTheta = "1.0";
    std::string tDensityPhi   = "0.0";

    // prescribed theta on interface
    std::string tPrescTheta = "1.0";

    // prescribed phi on interface
    std::string tPrescPhi = "0.0";

    /* ------------------------------------------------------------------------ */

    /* ------------------------------------------------------------------------ */

    void
    tLevelSetFuncReal(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // return absolute value of level set function
        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val();
    }

    // Level set function defining property in FEM
    void
    tLevelSetFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        // return PDV derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsolutValue )
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
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        // return PDV derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsolutValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->N();
    }

    /* ------------------------------------------------------------------------ */

    // Spatial derivative of level set function defining property in FEM
    void
    tLevelSetGradxFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        // return spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsolutValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->gradx( 1 );
    }

    /* ------------------------------------------------------------------------ */

    // Derivative of spatial derivative of level set function with respect to PDV
    void
    tDerLevelSetGradxFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // get value of design level set function
        real value = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val()( 0 );

        // return PDV derivative of spatial derivative of absolute value of level set function
        real factor = 1.0;
        if ( tUseAbsolutValue )
        {
            factor = value > 0.0 ? 1.0 : -1.0;
        }

        aPropMatrix = factor * aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->dnNdxn( 1 );
    }

    //------------------------------------------------------------------------------

    // Hole pattern
    real
    Box_2D3D(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&    aGeometryParameters )
    {
        real tBoxExponent = 24.0;

        Matrix< DDRMat > tCenter   = { { 1.5, 0.5, 0.5 } };
        Matrix< DDRMat > tDimOUter = { { 1.5 - tWallThickness, 0.5 - tWallThickness, 0.5 - tWallThickness } };

        real XContrib = std::pow( ( tCenter( 0 ) - aCoordinates( 0 ) ) / tDimOUter( 0 ), tBoxExponent );
        real YContrib = std::pow( ( tCenter( 1 ) - aCoordinates( 1 ) ) / tDimOUter( 1 ), tBoxExponent );

        real tLSval;

        if ( tIs3D )
        {
            real tBoxScaling = pow( tDimOUter( 0 ) * tDimOUter( 1 ) * tDimOUter( 2 ), 1.0 / 3.0 );

            real ZContrib = std::pow( ( tCenter( 2 ) - aCoordinates( 2 ) ) / tDimOUter( 2 ), tHoleExponent );

            tLSval = tBoxScaling * ( 1.0 - std::pow( XContrib + YContrib + ZContrib, 1.0 / tBoxExponent ) );
        }
        else
        {
            real tBoxScaling = pow( tDimOUter( 0 ) * tDimOUter( 1 ), 1.0 / 2.0 );

            tLSval = tBoxScaling * ( 1.0 - std::pow( XContrib + YContrib, 1.0 / tBoxExponent ) );
        }

        // clean return value to return non-zero value
        return -tLSval;
    }

    //------------------------------------------------------------------------------

    // Hole pattern
    real
    Hole_Pattern_2D3D(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&    aGeometryParameters )
    {
        Matrix< DDRMat > tDelta  = { { tHoleXdim / tNumHolesX }, { tHoleYdim / tNumHolesY }, { tHoleZdim / tNumHolesZ } };
        Matrix< DDRMat > tOrigin = { { tHoleXOrg + tDelta( 0 ) / 2.0 }, { tHoleYOrg + tDelta( 1, 0 ) / 2.0 }, { tHoleZOrg + tDelta( 2, 0 ) / 2.0 } };
        Matrix< DDRMat > tOffset = { { tDelta( 0, 0 ) / 2.0 }, { 0.0 }, { tDelta( 2, 0 ) / 2.0 } };

        bool tOddRow = true;

        real tLSval = -1e20;
        real dist   = 0.0;

        for ( int iz = -2; iz < tNumHolesZ + 2; ++iz )
        {
            for ( int iy = -2; iy < tNumHolesY + 2; ++iy )
            {
                for ( int ix = -2; ix < tNumHolesX + 2; ++ix )
                {

                    real XCenter = tOrigin( 0 ) + (real)ix * tDelta( 0 ) + tOddRow * tOffset( 0 );
                    real YCenter = tOrigin( 1 ) + (real)iy * tDelta( 1 ) + tOddRow * tOffset( 1 );

                    real XContrib = std::pow( XCenter - aCoordinates( 0 ), tHoleExponent );
                    real YContrib = std::pow( YCenter - aCoordinates( 1 ), tHoleExponent );

                    if ( tIs3D )
                    {
                        real ZCenter  = tOrigin( 2 ) + (real)iz * tDelta( 2 ) + tOddRow * tOffset( 2 );
                        real ZContrib = std::pow( ZCenter - aCoordinates( 2 ), tHoleExponent );

                        dist = tHoleScaling * ( tHoleRadius - std::pow( XContrib + YContrib + ZContrib, 1.0 / tHoleExponent ) );
                    }
                    else
                    {
                        dist = tHoleScaling * ( tHoleRadius - std::pow( XContrib + YContrib, 1.0 / tHoleExponent ) );
                    }
                    tLSval = std::max( tLSval, dist );
                }

                tOddRow = ( !tOddRow && tHoleOffsetRow );
            }
        }

        // clean return value to return non-zero value
        return -tLSval;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Func_Neumann_U(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
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
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );

        return tConstraintTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );

        real obj1 = aCriteria( 0 ) / tInitialStrainEnergy;
        real obj2 = aCriteria( 1 ) / tInitialStrainEnergy;
        real obj3 = tPerimeterPenalty * aCriteria( 3 ) / tInitialPerimeter;
        real obj4 = tRegularizationPenalty * aCriteria( 4 ) / tInitialRegularization;
        // real obj5 = tRegularizationPenalty* aCriteria( 5 ) ;

        tObjectives( 0, 0 ) = obj1 + obj2 + obj3 + obj4;    // + obj5;

        std::cout << "% --------------------------------- % \n";
        std::cout << "Objective                = " << tObjectives( 0, 0 ) << " \n";
        std::cout << "Strain Energy (Frame)    = " << aCriteria( 0 ) << " ( " << obj1 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Strain Energy (Interior) = " << aCriteria( 1 ) << " ( " << obj2 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Perimeter                = " << aCriteria( 3 ) << " ( " << obj3 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "H1Error Const            = " << aCriteria( 4 ) << " ( " << aCriteria( 0 ) / tObjectives( 0, 0 ) << " )\n";
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
        tConstraints( 0 ) = aCriteria( 2 ) / tMaxMass - 1.0;

        std::cout << "Volume     = " << aCriteria( 2 ) << " \n";
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
        tDObjectiveDCriteria( 1 ) = 1.0 / tInitialStrainEnergy;
        tDObjectiveDCriteria( 3 ) = tPerimeterPenalty / tInitialPerimeter;
        tDObjectiveDCriteria( 4 ) = tRegularizationPenalty / tInitialRegularization;
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

        tDConstraintDCriteria( 2 ) = 1.0 / tMaxMass;

        return tDConstraintDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", tIsOpt );
        aParameterLists( 0 ).set( "problem", "user_defined" );
        aParameterLists( 0 ).set( "library", tLibraryName );
        aParameterLists( 0 ).set( "restart_file", "ADV_Alg_0_Iter_11.hdf5" );
        aParameterLists( 0 ).set( "reinitialize_interface_iter", 50 );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_gcmma_parameter_list() );
        aParameterLists( 2 ).set( "step_size", tMMAStepSize );
        aParameterLists( 2 ).set( "penalty", tMMAPenalty );
        aParameterLists( 2 ).set( "max_its", tMMAMaxIter );    // Maximum number of iterations
        aParameterLists( 2 ).set( "restart_index", 10 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", tNumElementsPerDir );
        aParameterLists( 0 ).set( "domain_dimensions", tDimensions );
        aParameterLists( 0 ).set( "domain_offset", tOffSet );
        aParameterLists( 0 ).set( "domain_sidesets", tSideSets );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        aParameterLists( 0 ).set( "lagrange_orders", "1,1" );
        aParameterLists( 0 ).set( "lagrange_pattern", "0,1" );

        aParameterLists( 0 ).set( "bspline_orders", "1,1" );
        aParameterLists( 0 ).set( "bspline_pattern", "0,1" );

        aParameterLists( 0 ).set( "initial_refinement", "0,0" );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0,1" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0,1;-1" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", 1 );
        aParameterLists( 0 ).set( "staircase_buffer", 1 );

        aParameterLists( 0 ).set( "use_number_aura", 1 );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );

        aParameterLists( 0 ).set( "write_lagrange_output_mesh_to_exodus", "HMRLagrangeMesh.exo" );

        aParameterLists( 0 ).set( "use_refine_low_level_elements", false );

        aParameterLists( 0 ).set( "restart_refinement_pattern_file", "HMR_Background_Refinement_Iter_11.hdf5" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0,1" );
        aParameterLists( 0 ).set( "ghost_stab", tUseGhost );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", false );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "IQI_types", "IQIBulkStrainEnergy_Frame", "IQIBulkStrainEnergy_Interior", "IQIBulkVolume_Interior", "IQIPerimeter_InterfaceVoid", "IQIHeatMethodPenalty" );
        aParameterLists( 0 ).set( "output_mesh_file", tGENOutputFile );
        aParameterLists( 0 ).set( "time_offset", 10.0 );

        Matrix< DDUMat > tPhaseMap( 4, 1, 0 );
        tPhaseMap( 1 ) = 1;
        tPhaseMap( 2 ) = 2;
        tPhaseMap( 3 ) = 2;
        aParameterLists( 0 ).set( "phase_table", moris::ios::stringify( tPhaseMap ) );

        aParameterLists( 0 ).set( "print_phase_table", true );

        // outer frame
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Box_2D3D" );
        aParameterLists( 1 ).set( "number_of_refinements", 1, 0 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0, 1 );
        aParameterLists( 1 ).set( "isocontour_tolerance", 10e-14 );
        aParameterLists( 1 ).set( "name", "Box" );

        // initialize fins as swiss cheese geometry
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Hole_Pattern_2D3D" );
        aParameterLists( 1 ).set( "name", "Level_Set_Field" );
        aParameterLists( 1 ).set( "number_of_refinements", 1, 0 );
        aParameterLists( 1 ).set( "refinement_mesh_index", 0, 1 );

        if ( tIsOpt )
        {
            aParameterLists( 1 ).set( "discretization_mesh_index", 1 );
            aParameterLists( 1 ).set( "discretization_lower_bound", -tBsplineLimit );
            aParameterLists( 1 ).set( "discretization_upper_bound", tBsplineLimit );
        }

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_gen_property_parameter_list( gen::Field_Type::SCALED_FIELD ) );
        uint tParamCounter = 0;

        aParameterLists( 2 ).set( "name", "LevelsetField" );
        aParameterLists( 2 ).set( "dependencies", "Level_Set_Field" );
        aParameterLists( 2 ).set( "scaling_factor", 1.0 );
        aParameterLists( 2 ).set( "pdv_type", "LS1" );
        aParameterLists( 2 ).set( "pdv_mesh_set_names", tTotalDomain1 );
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
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensity" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropYoungs" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropBedding" );
        aParameterLists( 0 ).set( "function_parameters", "1.0e-6" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 5
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropFlux" );
        aParameterLists( 0 ).set( "function_parameters", "10.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDirichletU" );
        aParameterLists( 0 ).set( "function_parameters", tDirichletStr );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 10
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropTraction" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Neumann_U" );

        // create parameter list for property 7
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPoisson" );
        aParameterLists( 0 ).set( "function_parameters", "0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------
        // common properties for theta and phi problems

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropConductivity" );
        aParameterLists( 0 ).set( "function_parameters", tConductivity );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // properties for Theta

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensityTheta" );
        aParameterLists( 0 ).set( "function_parameters", tDensityTheta );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropCapacityTheta" );
        aParameterLists( 0 ).set( "function_parameters", tCapacityTheta );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPrescTheta" );
        aParameterLists( 0 ).set( "function_parameters", tPrescTheta );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // properties for phi problem

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensityPhi" );
        aParameterLists( 0 ).set( "function_parameters", tDensityPhi );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropCapacityPhi" );
        aParameterLists( 0 ).set( "function_parameters", tCapacityPhi );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPrescPhi" );
        aParameterLists( 0 ).set( "function_parameters", tPrescPhi );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropEigenStrainPhi" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLevelSetConst" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLevelSetGradxConst" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );
        aParameterLists( 0 ).set( "function_parameters", tLevelSetGradxConstant );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLevelSet" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "tLevelSetFunc" );
        aParameterLists( 0 ).set( "dv_derivative_functions", "tDerLevelSetFunc" );
        aParameterLists( 0 ).set( "dv_dependencies", "LS1" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLevelSetReal" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "tLevelSetFuncReal" );
        aParameterLists( 0 ).set( "dv_derivative_functions", "tDerLevelSetFunc" );
        aParameterLists( 0 ).set( "dv_dependencies", "LS1" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLevelSetGradx" );
        aParameterLists( 0 ).set( "function_parameters", "1.0" );
        aParameterLists( 0 ).set( "value_function", "tLevelSetGradxFunc" );
        aParameterLists( 0 ).set( "dv_derivative_functions", "tDerLevelSetGradxFunc" );
        aParameterLists( 0 ).set( "dv_dependencies", "LS1" );

        // time continuity weights
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropWeightCurrent" );
        aParameterLists( 0 ).set( "function_parameters", "10.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropWeightPrevious" );
        aParameterLists( 0 ).set( "function_parameters", "10.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // initial condition
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropInitialCondition" );
        aParameterLists( 0 ).set( "function_parameters", "0.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso_Frame" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        aParameterLists( 1 ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso_Interior" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        aParameterLists( 1 ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        // create parameter list for constitutive model - Theta problem
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusionTheta" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "THETA", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity      , Conductivity;"
                "PropDensityTheta      , Density;"
                "PropCapacityTheta     , HeatCapacity" );

        // create parameter list for constitutive model - Phi problem
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusionPhi" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "PHID;THETA", "Temperature,Theta" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity    , Conductivity;"
                "PropDensityPhi      , Density;"
                "PropCapacityPhi     , HeatCapacity;"
                "PropEigenStrainPhi  , EigenStrain" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheDirichletBC" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungs,Material" );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", std::string( "SPNitscheFrameInteriorInterface" ) );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists( 2 ).set( "function_parameters", std::string( "100.0" ) );
        aParameterLists( 2 ).set( "leader_properties", std::string( "PropYoungs,Material" ) );
        aParameterLists( 2 ).set( "follower_properties", std::string( "PropYoungs,Material" ) );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", std::string( "SPGhost_Frame" ) );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", std::string( "0.005" ) );
        aParameterLists( 2 ).set( "leader_properties", std::string( "PropYoungs,Material" ) );
        aParameterLists( 2 ).set( "follower_properties", std::string( "PropYoungs,Material" ) );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", std::string( "SPGhost_Interior" ) );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", std::string( "0.005" ) );
        aParameterLists( 2 ).set( "leader_properties", std::string( "PropYoungs,Material" ) );
        aParameterLists( 2 ).set( "follower_properties", std::string( "PropYoungs,Material" ) );

        // create parameter list for ghost stabilization parameter for theta and phi problems
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPTemp" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity,Material" );
        aParameterLists( 2 ).set( "follower_properties", "PropConductivity,Material" );

        // create parameter list for DBC on interface for theta problem
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", tNitschePenalty );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // create IWG  - bulk diffusion
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionThetaBulkInterior" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "THETA" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists( 3 ).set( "mesh_set_names", tInteriorSets );

        // create IWG  - bulk diffusion
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionThetaBulkVoid" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "THETA" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tHeatDofs );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists( 3 ).set( "mesh_set_names", tVoidSets );

        // create parameter list for single side interface condition
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGSurfaceInnerTheta" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "THETA" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists( 3 ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tInterfaceVoidSSets3 );

        // create parameter list for single side interface condition
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGSurfaceOuterTheta" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "THETA" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tHeatDofs );
        aParameterLists( 3 ).set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tInterfaceVoidSSets4 );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGTimeContinuityThetaVoid" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        aParameterLists( 3 ).set( "dof_residual", "THETA" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tHeatDofs );
        aParameterLists( 3 ).set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        aParameterLists( 3 ).set( "mesh_set_names", tVoidSets );
        aParameterLists( 3 ).set( "time_continuity", true );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGTimeContinuityTheta" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        aParameterLists( 3 ).set( "dof_residual", "THETA" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists( 3 ).set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        aParameterLists( 3 ).set( "mesh_set_names", tInteriorSets );
        aParameterLists( 3 ).set( "time_continuity", true );

        // create IWG - bulk diffusion
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionOuterBulkVoid" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "PHID" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tHeatDofs );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists( 3 ).set( "mesh_set_names", tVoidSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionOuterBulkInterior" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "PHID" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists( 3 ).set( "mesh_set_names", tInteriorSets );

        // create parameter list for single side interface condition
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGSurfaceInnerPhi" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "PHID" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists( 3 ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tInterfaceVoidSSets3 );

        // create parameter list for single side interface condition
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGSurfaceOuterPhi" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "PHID" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tHeatDofs );
        aParameterLists( 3 ).set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tInterfaceVoidSSets4 );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGBulkU_Frame" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", tDofStrg );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofStrg );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        aParameterLists( 3 ).set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists( 3 ).set( "mesh_set_names", tFrameSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGBulkU_FrameInterior" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", tDofStrg );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso_Interior,ElastLinIso" );
        aParameterLists( 3 ).set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists( 3 ).set( "mesh_set_names", tInteriorSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletU" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", tDofStrg );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofStrg );
        aParameterLists( 3 ).set( "leader_properties", "PropDirichletU,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", tFrameSupportSSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGTraction" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
        aParameterLists( 3 ).set( "dof_residual", tDofStrg );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofStrg );
        aParameterLists( 3 ).set( "leader_properties", "PropTraction,Traction" );
        aParameterLists( 3 ).set( "mesh_set_names", tFrameLoadSSets );

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", std::string( "IWGFrameInteriorInterface" ) );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", tDofStrg );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists( 3 ).set( "follower_dof_dependencies", tDofStrg );
        aParameterLists( 3 ).set( "leader_constitutive_models", std::string( "CMStrucLinIso_Frame,ElastLinIso" ) );
        aParameterLists( 3 ).set( "follower_constitutive_models", std::string( "CMStrucLinIso_Interior,ElastLinIso" ) );
        aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPNitscheFrameInteriorInterface,NitscheInterface" ) );
        aParameterLists( 3 ).set( "mesh_set_names", tFrameInteriorDSets );

        if ( tUseGhost )
        {
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", std::string( "IWGGhostFrame" ) );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( 3 ).set( "dof_residual", tDofStrg );
            aParameterLists( 3 ).set( "leader_dof_dependencies", tDofStrg );
            aParameterLists( 3 ).set( "follower_dof_dependencies", tDofStrg );
            aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPGhost_Frame,GhostSP" ) );
            aParameterLists( 3 ).set( "ghost_order", (uint)tDispOrder );
            aParameterLists( 3 ).set( "mesh_set_names", tFrameGhost );

            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", std::string( "IWGGhostInterior" ) );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists( 3 ).set( "dof_residual", tDofStrg );
            aParameterLists( 3 ).set( "leader_dof_dependencies", tDofStrgAll );
            aParameterLists( 3 ).set( "follower_dof_dependencies", tDofStrgAll );
            aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPGhost_Interior,GhostSP" ) );
            aParameterLists( 3 ).set( "ghost_order", (uint)tDispOrder );
            aParameterLists( 3 ).set( "mesh_set_names", tInteriorGhost );
        }

        //------------------------------------------------------------------------------
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkUX" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", tDofStrg );
        aParameterLists( 4 ).set( "leader_dof_dependencies", tDofStrg );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomainSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkUY" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", tDofStrg );
        aParameterLists( 4 ).set( "dof_quantity", tDofStrg );
        aParameterLists( 4 ).set( "vectorial_field_index", 1 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomainSets );

        if ( tIs3D )
        {
            aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
            aParameterLists( 4 ).set( "IQI_name", "IQIBulkUZ" );
            aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
            aParameterLists( 4 ).set( "leader_dof_dependencies", tDofStrg );
            aParameterLists( 4 ).set( "dof_quantity", tDofStrg );
            aParameterLists( 4 ).set( "vectorial_field_index", 2 );
            aParameterLists( 4 ).set( "mesh_set_names", tTotalDomainSets );
            }

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkStrainEnergy_Frame" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", tDofStrg );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso_Frame,Elast" );
        aParameterLists( 4 ).set( "mesh_set_names", tFrameSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkStrainEnergy_Interior" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", tDofStrg );
        aParameterLists( 4 ).set( "leader_constitutive_models", "CMStrucLinIso_Interior,Elast" );
        aParameterLists( 4 ).set( "mesh_set_names", tInteriorSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkVolume_Frame" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists( 4 ).set( "leader_properties", "PropDensity,Density" );
        aParameterLists( 4 ).set( "mesh_set_names", tFrameSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkVolume_Interior" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists( 4 ).set( "leader_properties", "PropDensity,Density" );
        aParameterLists( 4 ).set( "mesh_set_names", tInteriorSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIPerimeter_InterfaceVoid" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", tDofStrg );
        aParameterLists( 4 ).set( "mesh_set_names", tInterfaceVoidSSets );

        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkTHETA" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "THETA" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "THETA" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // Nodal PHID IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkPHID" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "PHID" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "PHID" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // H1 Error if reference is constant
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIH1ErrorConst" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::H1_ERROR ) ;
        aParameterLists( 4 ).set( "dof_quantity", "PHID" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "leader_properties", "PropLevelSetConst,L2_Reference;PropLevelSetGradxConst,H1S_Reference" );
        aParameterLists( 4 ).set( "function_parameters", "1.0 / 1.0 / 1.0" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // H1 Error if reference is design dependent
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIH1Error" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::H1_ERROR ) ;
        aParameterLists( 4 ).set( "dof_quantity", "PHID" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "leader_properties", "PropLevelSet,L2_Reference;PropLevelSetGradx,H1S_Reference" );
        aParameterLists( 4 ).set( "function_parameters", "1.0 / 1.0 / 1.0" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // H1 Error if reference is design dependent
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIHeatMethodPenalty" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::HEAT_METHOD_PENALTY ) ;
        aParameterLists( 4 ).set( "dof_quantity", "PHID" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "THETA;PHID" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "leader_properties", "PropLevelSet,L2_Reference;PropLevelSetGradx,H1S_Reference" );
        aParameterLists( 4 ).set( "function_parameters", "1.6 / 1.0 / 36.0 / 0.2 / 0.2 / 0.333 / 0.333 " );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain );

        // Design level set function
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQILevelSet" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::PROPERTY ) ;
        aParameterLists( 4 ).set( "leader_properties", "PropLevelSetReal,Property" );
        aParameterLists( 4 ).set( "mesh_set_names", tTotalDomain1 );

        // create computation  parameter list
        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
        aParameterLists( 0 ).set( "Solver_Type", "Amesos_Umfpack" );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 0
        aParameterLists( 2 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 0
        aParameterLists( 2 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );    // nonlinear algorithm index 1
        aParameterLists( 2 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1.0e-9 );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );

        //------------------------------------------------------------------------------

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );    // nonlinear solver index 0
        aParameterLists( 3 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0" );             // set nonlinear algorithm with index 0
        aParameterLists( 3 ).set( "NLA_DofTypes", "THETA" );
        aParameterLists( 3 ).set( "NLA_Secondary_DofTypes", "" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );    // nonlinear solver index 1
        aParameterLists( 3 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "1" );             // set nonlinear algorithm with index 0
        aParameterLists( 3 ).set( "NLA_DofTypes", "PHID" );
        aParameterLists( 3 ).set( "NLA_Secondary_DofTypes", "" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) ;
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "3" );    // set nonlinear algorithm with index 0
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );    // nonlinear solver index 2
        aParameterLists( 3 ).set( "NLA_Solver_Implementation",  moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) ;
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "2" );             // set nonlinear algorithm with index 1.
        aParameterLists( 3 ).set( "NLA_Sub_Nonlinear_Solver", "0,1,2" );                // set sub nonlinear solvers with index 0 and 1
        aParameterLists( 3 ).set( "NLA_DofTypes", "THETA;PHID;UX,UY" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Nonlinear_Solver", 3 );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "THETA;PHID;UX,UY" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists( 0 ).set( "UX", 0 );
        aParameterLists( 0 ).set( "UY", 0 );
        if ( tIs3D )
        {
            aParameterLists( 0 ).set( "UZ", 0 );
        }

        aParameterLists( 0 ).set( "THETA", 1 );
        aParameterLists( 0 ).set( "PHID", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        // aParameterLists( 0 ).set( "Set_Names"  , tTotalDomain );

        if ( tIs3D )
        {
            aParameterLists( 0 ).set( "Field_Names", std::string( "UX,UY,UZ,StrainEnergyFrame,StrainEnergyInterior,VolumeInterior,PerimeterInteriorVoid" ) );
            aParameterLists( 0 ).set( "Field_Type", std::string( "NODAL,NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL,GLOBAL" ) );
            aParameterLists( 0 ).set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkUZ,IQIBulkStrainEnergy_Frame,IQIBulkStrainEnergy_Interior,"
                                                                    "IQIBulkVolume_Interior,IQIPerimeter_InterfaceVoid" ) );
        }
        else
        {

            aParameterLists( 0 ).set( "Set_Names", tTotalDomain1 );
            aParameterLists( 0 ).set( "Field_Names", std::string( "UX,UY,PHID,THETA,LEVELSET" ) );
            aParameterLists( 0 ).set( "Field_Type", std::string( "NODAL,NODAL,NODAL,NODAL,NODAL" ) );
            aParameterLists( 0 ).set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkPHID,IQIBulkTHETA,IQILevelSet" ) );

            // aParameterLists( 0 ).set( "Field_Names", std::string("THETA,PHID") );
            // aParameterLists( 0 ).set( "Field_Type" , std::string("NODAL,NODAL") );
            // aParameterLists( 0 ).set( "IQI_Names"  , std::string("IQIBulkTHETA,IQIBulkPHID") );
        }

        aParameterLists( 0 ).set( "Save_Frequency", 1 );
        aParameterLists( 0 ).set( "Time_Offset", 10.0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( Parameter_List( "" ) );
        prm::create_remeshing_parameterlist( aParameterLists( 0 )( 0 ) );
        aParameterLists( 0 ).set( "mode", "ab_initio" );
        aParameterLists( 0 ).set( "remeshing_field_names", "Box,Level_Set_Field" );
        aParameterLists( 0 ).set( "remeshing_levels_of_refinement", "2,1;2,1" );
        aParameterLists( 0 ).set( "remeshing_refinement_pattern", "0,1;0,1" );
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif

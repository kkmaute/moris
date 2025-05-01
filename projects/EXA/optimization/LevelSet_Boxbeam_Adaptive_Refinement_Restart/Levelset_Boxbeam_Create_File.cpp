/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Levelset_Boxbeam_Create_File.cpp
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

    std::string tName = "Levelset_Boxbeam_Create_File";

    bool tIs3D     = false;
    bool tIsOpt    = true;
    bool tUseGhost = true;

    bool tUseAbsolutValue = true;

    // wall thickness
    real tWallThickness = 0.05;

    // background mesh parameters
    Vector< uint > tNumElementsPerDir = tIs3D ? Vector< uint >{ 45, 15, 15 } : Vector< uint >{ 30, 10 };
    Vector< real > tDimensions        = tIs3D ? Vector< real >{ 3, 1, 1 } : Vector< real >{ 3, 1 };

    int tDispOrder = 1;

    // Hole Seeding parameters
    sint tNumHolesY = 2;
    sint tNumHolesZ = 1 * tNumHolesY;
    sint tNumHolesX = tNumHolesY + 1;

    real tHoleRadius   = tIs3D ? 0.4771 / tNumHolesY : 0.4771 / tNumHolesY;
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
    int  tMMAMaxIter  = 12;

    real tBsplineLimit = tHoleRadius;

    // other mesh depedendent parameters
    real tElementEdgeLength = 1.0 / tNumHolesX / pow( 2, 3 );
    real tLoadLimitY        = std::floor( 0.2 / tElementEdgeLength ) * tElementEdgeLength;

    //------------------------------------------------------------------------------
    // Derived problem paramters

    std::string tOutputFileName = tName + ".exo";
    std::string tLibraryName    = "Levelset_Boxbeam_Create_File.so";
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
            std::string( "SideSet_1_n_p2,SideSet_1_c_p2" ) +    //
            std::string( "SideSet_3_n_p2,SideSet_3_c_p2" );

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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        // return absolute value of level set function
        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::LS1 )->val();
    }

    // Level set function defining property in FEM
    void
    tLevelSetFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
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
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
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
            const Vector< real >&   aGeometryParameters )
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
            const Vector< real >&   aGeometryParameters )
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
        aParameterLists.set( "is_optimization_problem", tIsOpt );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tLibraryName );
        aParameterLists.set( "restart_file", "" );
        aParameterLists.set( "reinitialize_interface_iter", 10 );

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
        aParameterLists.set( "number_of_elements_per_dimension", tNumElementsPerDir );
        aParameterLists.set( "domain_dimensions", tDimensions );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "1,1" );
        aParameterLists.set( "lagrange_pattern", "0,1" );

        aParameterLists.set( "bspline_orders", "1,1" );
        aParameterLists.set( "bspline_pattern", "0,1" );

        aParameterLists.set( "pattern_initial_refinement", 0, 0 );

        aParameterLists.set( "lagrange_to_bspline", "0,1;-1" );

        aParameterLists.set( "refinement_buffer", 1 );
        aParameterLists.set( "staircase_buffer", 1 );
        aParameterLists.set( "use_refine_low_level_elements", false );
        aParameterLists.set( "write_refinement_pattern_file", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
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

        aParameterLists.set( "IQI_types", "IQIBulkStrainEnergy_Frame", "IQIBulkStrainEnergy_Interior", "IQIBulkVolume_Interior", "IQIPerimeter_InterfaceVoid", "IQIHeatMethodPenalty" );
        // aParameterLists.set("IQI_types"       , "IQIBulkStrainEnergy_Frame,IQIBulkStrainEnergy_Interior,IQIBulkVolume_Interior,IQIPerimeter_InterfaceVoid,IQIH1ErrorConst,IQIH1Error");
        aParameterLists.set( "output_mesh_file", tGENOutputFile );
        aParameterLists.set( "time_offset", 10.0 );

        Matrix< DDUMat > tPhaseMap( 4, 1, 0 );
        tPhaseMap( 1 ) = 1;
        tPhaseMap( 2 ) = 2;
        tPhaseMap( 3 ) = 2;
        aParameterLists.set( "phase_table", moris::ios::stringify( tPhaseMap ) );

        aParameterLists.set( "print_phase_table", true );

        // outer frame
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "isocontour_tolerance", 10e-14 );
        aParameterLists.set( "field_function_name", "Box_2D3D" );
        aParameterLists.set( "number_of_refinements", 1, 0 );
        aParameterLists.set( "refinement_mesh_index", 0, 1 );
        aParameterLists.set( "name", "Box" );

        // initialize fins as swiss cheese geometry
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Hole_Pattern_2D3D" );
        aParameterLists.set( "name", "Level_Set_Field" );
        aParameterLists.set( "number_of_refinements", 1, 0 );
        aParameterLists.set( "refinement_mesh_index", 0, 1 );

        if ( tIsOpt )
        {
            aParameterLists.set( "discretization_mesh_index", 1 );
            aParameterLists.set( "discretization_lower_bound", -tBsplineLimit );
            aParameterLists.set( "discretization_upper_bound", tBsplineLimit );
        }

        aParameterLists( GEN::PROPERTIES ).add_parameter_list( gen::Field_Type::SCALED_FIELD );
        uint tParamCounter = 0;

        aParameterLists.set( "name", "LevelsetField" );
        aParameterLists.set( "dependencies", "Level_Set_Field" );
        aParameterLists.set( "scaling_factor", 1.0 );
        aParameterLists.set( "pdv_type", "LS1" );
        aParameterLists.set( "pdv_mesh_set_names", "HMR_dummy_n_p0", "HMR_dummy_c_p0", "HMR_dummy_n_p1", "HMR_dummy_c_p1", "HMR_dummy_n_p2", "HMR_dummy_c_p2" );
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
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropBedding" );
        aParameterLists.set( "function_parameters", "1.0e-6" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 5
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFlux" );
        aParameterLists.set( "function_parameters", "10.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletU" );
        aParameterLists.set( "function_parameters", tDirichletStr );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 10
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropTraction" );
        aParameterLists.set( "function_parameters", "1.0" );
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
        aParameterLists.set( "function_parameters", tLevelSetGradxConstant );

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
        aParameterLists.set( "constitutive_name", "CMStrucLinIso_Frame" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso_Interior" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
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

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheDirichletBC" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPNitscheFrameInteriorInterface" ) );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists.set( "function_parameters", std::string( "100.0" ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungs,Material" ) );
        aParameterLists.set( "follower_properties", std::string( "PropYoungs,Material" ) );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPGhost_Frame" ) );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", std::string( "0.005" ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungs,Material" ) );
        aParameterLists.set( "follower_properties", std::string( "PropYoungs,Material" ) );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", std::string( "SPGhost_Interior" ) );
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
        aParameterLists.set( "follower_properties", "PropConductivity,Material" );

        // create parameter list for DBC on interface for theta problem
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", tNitschePenalty );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // create IWG  - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionThetaBulkInterior" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists.set( "mesh_set_names", tInteriorSets );

        // create IWG  - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionThetaBulkVoid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", tHeatDofs );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists.set( "mesh_set_names", tVoidSets );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceInnerTheta" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists.set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tInterfaceVoidSSets3 );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceOuterTheta" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", tHeatDofs );
        aParameterLists.set( "leader_properties", "PropPrescTheta,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionTheta,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tInterfaceVoidSSets4 );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTimeContinuityThetaVoid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", tHeatDofs );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        aParameterLists.set( "mesh_set_names", tVoidSets );
        aParameterLists.set( "time_continuity", true );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTimeContinuityTheta" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists.set( "leader_properties",
                "PropWeightCurrent,       WeightCurrent;"
                "PropWeightPrevious,      WeightPrevious;"
                "PropInitialCondition,    InitialCondition" );
        aParameterLists.set( "mesh_set_names", tInteriorSets );
        aParameterLists.set( "time_continuity", true );

        // create IWG - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionOuterBulkVoid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", tHeatDofs );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists.set( "mesh_set_names", tVoidSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionOuterBulkInterior" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists.set( "mesh_set_names", tInteriorSets );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceInnerPhi" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists.set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tInterfaceVoidSSets3 );

        // create parameter list for single side interface condition
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGSurfaceOuterPhi" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "PHID" );
        aParameterLists.set( "leader_dof_dependencies", tHeatDofs );
        aParameterLists.set( "leader_properties", "PropPrescPhi,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionPhi,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tInterfaceVoidSSets4 );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkU_Frame" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tFrameSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkU_FrameInterior" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Interior,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tInteriorSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichletU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_properties", "PropDirichletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tFrameSupportSSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTraction" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_properties", "PropTraction,Traction" );
        aParameterLists.set( "mesh_set_names", tFrameLoadSSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", std::string( "IWGFrameInteriorInterface" ) );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrgAll );
        aParameterLists.set( "follower_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", std::string( "CMStrucLinIso_Frame,ElastLinIso" ) );
        aParameterLists.set( "follower_constitutive_models", std::string( "CMStrucLinIso_Interior,ElastLinIso" ) );
        aParameterLists.set( "stabilization_parameters", std::string( "SPNitscheFrameInteriorInterface,NitscheInterface" ) );
        aParameterLists.set( "mesh_set_names", tFrameInteriorDSets );

        if ( tUseGhost )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGGhostFrame" ) );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", tDofStrg );
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "follower_dof_dependencies", tDofStrg );
            aParameterLists.set( "stabilization_parameters", std::string( "SPGhost_Frame,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tFrameGhost );

            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGGhostInterior" ) );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", tDofStrg );
            aParameterLists.set( "leader_dof_dependencies", tDofStrgAll );
            aParameterLists.set( "follower_dof_dependencies", tDofStrgAll );
            aParameterLists.set( "stabilization_parameters", std::string( "SPGhost_Interior,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tInteriorGhost );
        }

        //------------------------------------------------------------------------------
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkUX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomainSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkUY" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "dof_quantity", tDofStrg );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tTotalDomainSets );

        if ( tIs3D )
        {
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "IQIBulkUZ" );
            aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "dof_quantity", tDofStrg );
            aParameterLists.set( "vectorial_field_index", 2 );
            aParameterLists.set( "mesh_set_names", tTotalDomainSets );
            }

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy_Frame" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Frame,Elast" );
        aParameterLists.set( "mesh_set_names", tFrameSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy_Interior" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Interior,Elast" );
        aParameterLists.set( "mesh_set_names", tInteriorSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume_Frame" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tFrameSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume_Interior" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tInteriorSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIPerimeter_InterfaceVoid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "mesh_set_names", tInterfaceVoidSSets );

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
        aParameterLists.set( "function_parameters", "1.6 / 1.0 / 36.0 / 0.2 / 0.2 / 0.333 / 0.333 " );
        aParameterLists.set( "mesh_set_names", tTotalDomain );

        // Design level set function
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQILevelSet" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropLevelSetReal,Property" );
        aParameterLists.set( "mesh_set_names", tTotalDomain1 );

        // create computation  parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
        aParameterLists.set( "Solver_Type", "Amesos_Umfpack" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-9 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 0
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-9 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();    // nonlinear algorithm index 1
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NLBGS_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-9 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", moris::NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-9 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 1 );

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

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "UX", 0 );
        aParameterLists.set( "UY", 0 );
        if ( tIs3D )
        {
            aParameterLists.set( "UZ", 0 );
        }

        aParameterLists.set( "THETA", 1 );
        aParameterLists.set( "PHID", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        // aParameterLists.set( "Set_Names"  , tTotalDomain );

        if ( tIs3D )
        {
            aParameterLists.set( "Field_Names", std::string( "UX,UY,UZ,StrainEnergyFrame,StrainEnergyInterior,VolumeInterior,PerimeterInteriorVoid" ) );
            aParameterLists.set( "Field_Type", std::string( "NODAL,NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL,GLOBAL" ) );
            aParameterLists.set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkUZ,IQIBulkStrainEnergy_Frame,IQIBulkStrainEnergy_Interior,"
                                                                    "IQIBulkVolume_Interior,IQIPerimeter_InterfaceVoid" ) );
        }
        else
        {

            aParameterLists.set( "Set_Names", tTotalDomain1 );
            aParameterLists.set( "Field_Names", std::string( "UX,UY,PHID,THETA,LEVELSET" ) );
            aParameterLists.set( "Field_Type", std::string( "NODAL,NODAL,NODAL,NODAL,NODAL" ) );
            aParameterLists.set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkPHID,IQIBulkTHETA,IQILevelSet" ) );

            // aParameterLists.set( "Field_Names", std::string("THETA,PHID") );
            // aParameterLists.set( "Field_Type" , std::string("NODAL,NODAL") );
            // aParameterLists.set( "IQI_Names"  , std::string("IQIBulkTHETA,IQIBulkPHID") );
        }

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
        aParameterLists.set( "remeshing_levels_of_refinement", "2,1;2,1" );
        aParameterLists.set( "remeshing_refinement_pattern", "0,1;0,1" );
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif

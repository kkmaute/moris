/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Levelset_Boxbeam_Lin_Quad_LevelSet.cpp
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
    //------------------------------------------------------------------------------

    // Main problem parameters

    std::string tName = "Levelset_Boxbeam_Lin_Quad_LevelSet";

    bool tIs3D     = false;
    bool tIsOpt    = true;
    bool tUseGhost = true;

    // wall thickness
    real tWallThickness = 0.05;

    // background mesh parameters
    std::string tNumElementsPerDir = tIs3D ? "45,15,15" : "15,5";
    std::string tDimensions        = tIs3D ? "3,1,1" : "3,1";
    std::string tOffSet            = tIs3D ? "0.0,0.0,0.0" : "0.0,0.0";

    int tLevelsetOrder = 2;
    int tDispOrder     = 1;

    int tLevelsetInitialRef = 2;
    int tDispInitialRef     = 3;
    int tRefineBuffer       = 0;

    // note: pattern 0 - Levelset field  pattern 1 - displacement field
    std::string tLagrangeOrder   = "1";    // std::to_string(std::max(tLevelsetOrder,tDispOrder));
    std::string tBsplineOrder    = std::to_string( tLevelsetOrder ) + "," + std::to_string( tDispOrder );
    std::string tInitialRef      = std::to_string( tLevelsetInitialRef ) + "," + std::to_string( tDispInitialRef );
    std::string tLagrangePattern = tLevelsetInitialRef > tDispInitialRef ? "0" : "1";

    // Hole Seeding parameters
    sint tNumHolesY = 3;
    sint tNumHolesZ = 1 * tNumHolesY;
    sint tNumHolesX = 3 * tNumHolesY;

    real tHoleRadius   = tIs3D ? 0.355 / tNumHolesY : 0.395 / tNumHolesY;
    real tHoleExponent = 6.0;
    real tHoleScaling  = 1.0;

    real tHoleXdim = 3.0;
    real tHoleYdim = 1.0;
    real tHoleZdim = 1.0;

    real tHoleXOrg = 0.0;
    real tHoleYOrg = 0.0;
    real tHoleZOrg = 0.0;

    bool tHoleOffsetRow = false;

    // optimization parameters
    real tInitialStrainEnergy = tIs3D ? 1.49721 + 1.47892 : 4.306294e+00 + 4.817279e+00;
    real tInitialPerimeter    = tIs3D ? 27.3111 : 41.5976;
    real tPerimeterPenalty    = 0.2;

    real tMaxMass = tIs3D ? 1.5 : 1.0;

    real tMMAPenalty  = 5.0;
    real tMMAStepSize = 0.05;
    int  tMMAMaxIter  = 3;

    real tBsplineLimit = tHoleRadius;

    // other mesh depedendent parameters
    real tElementEdgeLength = 1.0 / tNumHolesX / pow( 2, tDispInitialRef );
    real tLoadLimitY        = std::floor( 0.2 / tElementEdgeLength ) * tElementEdgeLength;

    //------------------------------------------------------------------------------
    // Derived problem paramters

    std::string tOutputFileName = tName + ".exo";
    std::string tLibraryName    = tName + "so";
    std::string tGENOutputFile  = "GEN_" + tName + ".exo";

    std::string tFrameSets    = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tInteriorSets = "HMR_dummy_n_p1,HMR_dummy_c_p1";

    std::string tTotalDomainSets = tFrameSets + "," + tInteriorSets;

    std::string tFrameGhost    = "ghost_p2";
    std::string tInteriorGhost = "ghost_p1";

    std::string tFrameLoadSSets    = "SideSet_2_n_p2,SideSet_2_c_p2";
    std::string tFrameSupportSSets = "SideSet_4_n_p2,SideSet_4_c_p2";
    std::string tFrameFreeSSets =
            std::string( "SideSet_1_n_p2,SideSet_1_c_p2" ) + std::string( "SideSet_3_n_p2,SideSet_3_c_p2" );

    std::string tInterfaceVoidSSets = "iside_b0_2_b1_0,iside_b0_1_b1_0";

    std::string tFrameInteriorDSets = "dbl_iside_p0_1_p1_2";

    std::string tDofStrg = tIs3D ? "UX,UY,UZ" : "UX,UY";

    std::string tDirichletStr = tIs3D ? "0.0;0.0;0.0" : "0.0;0.0";

    //------------------------------------------------------------------------------

    // Hole pattern
    real Box_2D3D(
            const Matrix< DDRMat >  & aCoordinates,
            const Vector< real >     & aGeometryParameters )
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
    real Hole_Pattern_2D3D(
            const Matrix< DDRMat >  & aCoordinates,
            const Vector< real >     & aGeometryParameters )
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

        tObjectives( 0, 0 ) = obj1 + obj2 + obj3;

        std::cout << "% --------------------------------- % \n";
        std::cout << "Objective                = " << tObjectives( 0, 0 ) << " \n";
        std::cout << "Strain Energy (Frame)    = " << aCriteria( 0 ) << " ( " << obj1 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Strain Energy (Interior) = " << aCriteria( 1 ) << " ( " << obj2 / tObjectives( 0, 0 ) << " )\n";
        std::cout << "Perimeter                = " << aCriteria( 3 ) << " ( " << obj3 / tObjectives( 0, 0 ) << " )\n";
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
    OPTParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( moris::prm::create_opt_problem_parameter_list() );
        aParameterLists.set( "is_optimization_problem", tIsOpt );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tLibraryName );
        aParameterLists.set( "restart_file", "" );

        aParameterLists( 2 ).push_back( opt::Optimization_Algorithm_Type::GCMMA );
        aParameterLists.set( "step_size", tMMAStepSize );
        aParameterLists.set( "penalty", tMMAPenalty );
        aParameterLists.set( "max_its", tMMAMaxIter );    // Maximum number of iterations
        aParameterLists.set( "restart_index", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_hmr_parameter_list() );

        aParameterLists.set( "number_of_elements_per_dimension", tNumElementsPerDir );
        aParameterLists.set( "domain_dimensions", tDimensions );
        aParameterLists.set( "domain_offset", tOffSet );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", tLagrangeOrder );
        aParameterLists.set( "lagrange_pattern", tLagrangePattern );

        aParameterLists.set( "bspline_orders", tBsplineOrder );
        aParameterLists.set( "bspline_pattern", "0,1" );

        aParameterLists.set( "initial_refinement", tInitialRef );
        aParameterLists.set( "initial_refinement_pattern", "0,1" );

        aParameterLists.set( "lagrange_to_bspline", "0,1" );

        aParameterLists.set( "refinement_buffer", tRefineBuffer );
        aParameterLists.set( "staircase_buffer", tRefineBuffer );

        aParameterLists.set( "use_number_aura", 1 );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_xtk_parameter_list() );
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
    GENParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {

        aParameterLists( 0 ).push_back( moris::prm::create_gen_parameter_list() );
        aParameterLists.set( "IQI_types", "IQIBulkStrainEnergy_Frame", "IQIBulkStrainEnergy_Interior", "IQIBulkVolume_Interior", "IQIPerimeter_InterfaceVoid" );
        aParameterLists.set( "output_mesh_file", tGENOutputFile );
        aParameterLists.set( "time_offset", 10.0 );

        Matrix< DDUMat > tPhaseMap( 4, 1, 0 );
        tPhaseMap( 1 ) = 1;
        tPhaseMap( 2 ) = 2;
        tPhaseMap( 3 ) = 2;
        aParameterLists.set( "phase_table", moris::ios::stringify( tPhaseMap ) );

        aParameterLists.set( "print_phase_table", true );

        // outer frame
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Box_2D3D" );

        // initialize fins as swiss cheese geometry
        aParameterLists( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Hole_Pattern_2D3D" );

        if ( tIsOpt )
        {
            aParameterLists.set( "discretization_mesh_index", 0 );
            aParameterLists.set( "discretization_lower_bound", -tBsplineLimit );
            aParameterLists.set( "discretization_upper_bound", tBsplineLimit );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // create parameter list for property 1
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropBedding" );
        aParameterLists.set( "function_parameters", "1.0e-6" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 5
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropFlux" );
        aParameterLists.set( "function_parameters", "10.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropDirichletU" );
        aParameterLists.set( "function_parameters", tDirichletStr );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 10
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropTraction" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Neumann_U" );

        // create parameter list for property 7
        aParameterLists( 0 ).push_back( prm::create_property_parameter_list() );
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMStrucLinIso_Frame" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
        aParameterLists.set( "constitutive_name", "CMStrucLinIso_Interior" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( tDofStrg, "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", "SPNitscheDirichletBC" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", std::string( "SPNitscheFrameInteriorInterface" ) );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::NITSCHE_INTERFACE ) ;
        aParameterLists.set( "function_parameters", std::string( "100.0" ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungs,Material" ) );
        aParameterLists.set( "follower_properties", std::string( "PropYoungs,Material" ) );

        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", std::string( "SPGhost_Frame" ) );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", std::string( "0.005" ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungs,Material" ) );

        aParameterLists( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists.set( "stabilization_name", std::string( "SPGhost_Interior" ) );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists.set( "function_parameters", std::string( "0.005" ) );
        aParameterLists.set( "leader_properties", std::string( "PropYoungs,Material" ) );

        //------------------------------------------------------------------------------
        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGBulkU_Frame" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tFrameSets );

        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGBulkU_Frame" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Interior,ElastLinIso" );
        aParameterLists.set( "leader_properties", "PropBedding,Bedding" );
        aParameterLists.set( "mesh_set_names", tInteriorSets );

        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGDirichletU" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_properties", "PropDirichletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tFrameSupportSSets );

        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", "IWGTraction" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_properties", "PropTraction,Traction" );
        aParameterLists.set( "mesh_set_names", tFrameLoadSSets );

        aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
        aParameterLists.set( "IWG_name", std::string( "IWGFrameInteriorInterface" ) );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "follower_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", std::string( "CMStrucLinIso_Frame,ElastLinIso" ) );
        aParameterLists.set( "follower_constitutive_models", std::string( "CMStrucLinIso_Interior,ElastLinIso" ) );
        aParameterLists.set( "stabilization_parameters", std::string( "SPNitscheFrameInteriorInterface,NitscheInterface" ) );
        aParameterLists.set( "mesh_set_names", tFrameInteriorDSets );

        if ( tUseGhost )
        {
            aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", std::string( "IWGGhostFrame" ) );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual", tDofStrg );
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "follower_dof_dependencies", tDofStrg );
            aParameterLists.set( "stabilization_parameters", std::string( "SPGhost_Frame,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tFrameGhost );

            aParameterLists( 3 ).push_back( prm::create_IWG_parameter_list() );
            aParameterLists.set( "IWG_name", std::string( "IWGGhostInterior" ) );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual", tDofStrg );
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "follower_dof_dependencies", tDofStrg );
            aParameterLists.set( "stabilization_parameters", std::string( "SPGhost_Interior,GhostSP" ) );
            aParameterLists.set( "ghost_order", (uint)tDispOrder );
            aParameterLists.set( "mesh_set_names", tInteriorGhost );
        }

        //------------------------------------------------------------------------------
        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkUX" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", tDofStrg );
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTotalDomainSets );

        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkUY" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "dof_quantity", tDofStrg );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tTotalDomainSets );

        if ( tIs3D )
        {
            aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
            aParameterLists.set( "IQI_name", "IQIBulkUZ" );
            aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
            aParameterLists.set( "leader_dof_dependencies", tDofStrg );
            aParameterLists.set( "dof_quantity", tDofStrg );
            aParameterLists.set( "vectorial_field_index", 2 );
            aParameterLists.set( "mesh_set_names", tTotalDomainSets );
            }

        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy_Frame" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Frame,Elast" );
        aParameterLists.set( "mesh_set_names", tFrameSets );

        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy_Interior" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::STRAIN_ENERGY ) ;
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso_Interior,Elast" );
        aParameterLists.set( "mesh_set_names", tInteriorSets );

        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkVolume_Frame" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tFrameSets );

        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIBulkVolume_Interior" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tInteriorSets );

        aParameterLists( 4 ).push_back( prm::create_IQI_parameter_list() );
        aParameterLists.set( "IQI_name", "IQIPerimeter_InterfaceVoid" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::VOLUME ) ;
        aParameterLists.set( "leader_dof_dependencies", tDofStrg );
        aParameterLists.set( "mesh_set_names", tInterfaceVoidSSets );

        // create computation  parameter list
        aParameterLists( 5 ).push_back( prm::create_computation_parameter_list() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {

        aParameterLists( 0 ).push_back( add_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( 1 ).push_back( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).push_back( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists.set( "NLA_combined_res_jac_assembly", true );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.00 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.00 );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( 3 ).push_back( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists.set( "NLA_DofTypes", tDofStrg );

        aParameterLists( 4 ).push_back( moris::prm::create_time_solver_algorithm_parameter_list() );

        aParameterLists( 5 ).push_back( moris::prm::create_time_solver_parameter_list() );
        aParameterLists.set( "TSA_DofTypes", tDofStrg );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( 6 ).push_back( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).push_back(  sol::PreconditionerType::NONE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_msi_parameter_list() );
        aParameterLists.set( "UX", 1 );
        aParameterLists.set( "UY", 1 );
        if ( tIs3D )
        {
            aParameterLists.set( "UZ", 1 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
        aParameterLists( 0 ).push_back( prm::create_vis_parameter_list() );
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", tTotalDomainSets );

        if ( tIs3D )
        {
            aParameterLists.set( "Field_Names", std::string( "UX,UY,UZ,StrainEnergyFrame,StrainEnergyInterior,VolumeInterior,PerimeterInteriorVoid" ) );
            aParameterLists.set( "Field_Type", std::string( "NODAL,NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL,GLOBAL" ) );
            aParameterLists.set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkUZ,IQIBulkStrainEnergy_Frame,IQIBulkStrainEnergy_Interior,"
                                                                    "IQIBulkVolume_Interior,IQIPerimeter_InterfaceVoid" ) );
        }
        else
        {
            aParameterLists.set( "Field_Names", std::string( "UX,UY,StrainEnergyFrame,StrainEnergyInterior,VolumeInterior,PerimeterInteriorVoid" ) );
            aParameterLists.set( "Field_Type", std::string( "NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL,GLOBAL" ) );
            aParameterLists.set( "IQI_Names", std::string( "IQIBulkUX,IQIBulkUY,IQIBulkStrainEnergy_Frame,IQIBulkStrainEnergy_Interior,"
                                                                    "IQIBulkVolume_Interior,IQIPerimeter_InterfaceVoid" ) );
        }

        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 10.0 );
    }

    void
    MORISGENERALParameterList( Vector< Vector< ParameterList > >& aParameterLists )
    {
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif

/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 */
#include "cl_Bitset.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_Matrix.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "parameters.hpp"
#include "moris_typedefs.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>

#ifdef MORIS_COULOMB_SHEAR_DRIVER
extern moris::uint gShearElements;
extern moris::uint gShearSubsteps;
extern moris::real gShearFriction;
extern moris::real gShearYoungsModulus;
extern moris::real gShearCompression;
extern moris::real gShearPenalty;
#else
// Standalone input defaults do not depend on EXA driver symbols.
moris::uint gShearElements      = 8;
moris::uint gShearSubsteps      = 1;
moris::real gShearFriction      = 0.3;
moris::real gShearYoungsModulus = 1000.0;
moris::real gShearCompression   = 1.0e-4;
moris::real gShearPenalty       = 500.0;
#endif

namespace
{
    constexpr moris::real tInitialGap                     = 0.02;
    constexpr moris::real tSolidHeight                    = 1.0 - tInitialGap;
    constexpr moris::uint tStages                         = 20;
    const moris::real     tDisplacementFactors[ tStages ] = {
        0,
        0.25,
        0.5,
        1,
        1.5,
        2,
        2,
        1.5,
        1,
        0.5,
        0,
        -0.5,
        -2,
        -2,
        -1.5,
        -1,
        -0.5,
        0,
        0.5,
        2
    };

    struct Loading_State
    {
        moris::real displacement;
        moris::real shear;
    };

    Loading_State loading_state( moris::real aTime )
    {
        const moris::real tG                    = 0.5 * gShearYoungsModulus;
        const moris::real tLambda               = 1.0 - gShearCompression / tSolidHeight;
        const moris::real tLimit                = gShearFriction * tG * ( 1.0 / tLambda - tLambda );
        const moris::real tCriticalDisplacement = tLimit * tSolidHeight / tG;
        const moris::uint tStep                 = static_cast< moris::uint >( std::round( aTime ) );
        const moris::uint tSteps                = tStages * gShearSubsteps;
        moris::real       tSlip                 = 0.0;
        Loading_State     tState{ 0.0, 0.0 };

        // This analytic history supplies the manufactured lateral traction only.
        // The EXA checker uses independent stored stage references.
        for ( moris::uint i = 1; i <= std::min( tStep, tSteps ); ++i )
        {
            const moris::uint tStage    = ( i - 1 ) / gShearSubsteps;
            const moris::real tFraction = static_cast< moris::real >( ( i - 1 ) % gShearSubsteps + 1 ) / gShearSubsteps;
            const moris::real tPrevious = tStage == 0 ? 0.0 : tDisplacementFactors[ tStage - 1 ];
            tState.displacement         = ( tPrevious + tFraction * ( tDisplacementFactors[ tStage ] - tPrevious ) ) * tCriticalDisplacement;
            tState.shear                = std::max( -tLimit, std::min( tLimit, tG * ( tState.displacement - tSlip ) ) );
            tSlip                       = tState.displacement - tState.shear / tG;
        }
        return tState;
    }

    moris::real evaluation_time( moris::fem::Field_Interpolator_Manager* aManager )
    {
        return aManager->get_IP_geometry_interpolator()->get_time_coeff()( 1 );
    }
}    // namespace

extern "C" {
namespace moris
{
    uint Shear_Phase( const Bitset< 2 >& aSigns )
    {
        return aSigns.test( 1 ) ? ( aSigns.test( 0 ) ? 1 : 0 ) : 2;
    }

    void Shear_Constant( Matrix< DDRMat >& aValue, Vector< Matrix< DDRMat > >& aParameters, fem::Field_Interpolator_Manager* )
    {
        aValue = aParameters( 0 );
    }

    void Shear_Top_Displacement( Matrix< DDRMat >& aValue, Vector< Matrix< DDRMat > >&, fem::Field_Interpolator_Manager* aManager )
    {
        const auto tState = loading_state( evaluation_time( aManager ) );
        aValue            = { { tState.displacement }, { -( tInitialGap + gShearCompression ) } };
    }

    void Shear_Lateral_Traction( Matrix< DDRMat >& aValue, Vector< Matrix< DDRMat > >& aParameters, fem::Field_Interpolator_Manager* aManager )
    {
        const auto tState = loading_state( evaluation_time( aManager ) );
        // For F=[[1,k],[0,lambda]], nu=0 gives P_yx=G*k/lambda.
        // Applying this side traction removes finite-width end effects in stick.
        aValue = { { 0.0 }, { aParameters( 0 )( 0 ) * tState.shear / ( 1.0 - gShearCompression / tSolidHeight ) } };
    }

    bool Shear_Output( tsa::Time_Solver* ) { return true; }

    void OPTParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "is_optimization_problem", false );
    }

    void HMRParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "number_of_elements_per_dimension", gShearElements, gShearElements );
        aLists.set( "domain_dimensions", 1.0, 1.0 );
        aLists.set( "lagrange_pattern", "0" );
        aLists.set( "lagrange_orders", "1" );
        aLists.set( "lagrange_output_meshes", "0" );
        aLists.set( "bspline_pattern", "0" );
        aLists.set( "bspline_orders", "1" );
        aLists.set( "lagrange_to_bspline", "0" );
        aLists.set( "refinement_buffer", 0 );
        aLists.set( "staircase_buffer", 0 );
    }

    void XTKParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "decompose", true );
        aLists.set( "decomposition_type", "conformal" );
        aLists.set( "enrich_mesh_indices", "0" );
        aLists.set( "ghost_stab", false );
        aLists.set( "multigrid", false );
        aLists.set( "high_to_low_dbl_side_sets", true );
    }

    void GENParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "number_of_phases", 3 );
        aLists.set( "phase_function_name", "Shear_Phase" );
        for ( real tY : { 0.5 + 0.5 * tInitialGap, 0.5 - 0.5 * tInitialGap } )
        {
            aLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            aLists.set( "center_x", 0.0 );
            aLists.set( "center_y", tY );
            aLists.set( "normal_x", 0.0 );
            aLists.set( "normal_y", 1.0 );
            aLists.set( "number_of_refinements", 0 );
        }
    }

    void FEMParameterList( Module_Parameter_Lists& aLists )
    {
        auto tVoid = prm::create_phase_parameter_list();
        tVoid.set( "phase_name", "Void" );
        tVoid.set( "phase_indices", "0" );
        aLists( FEM::PHASES ).add_parameter_list( tVoid );

        for ( const std::string tSide : { "Bottom", "Top" } )
        {
            const bool tTop   = tSide == "Top";
            auto       tPhase = prm::create_phase_parameter_list();
            tPhase.set( "phase_name", tSide );
            tPhase.set( "phase_indices", tTop ? "1" : "2" );
            aLists( FEM::PHASES ).add_parameter_list( tPhase );
        }

        const auto tProperty = [ & ]( const std::string& aName, const std::string& aParameters, const std::string& aFunction ) {
            auto tList = prm::create_property_parameter_list();
            tList.set( "property_name", aName );
            tList.set( "function_parameters", aParameters );
            tList.set( "value_function", aFunction );
            aLists( FEM::PROPERTIES ).add_parameter_list( tList );
        };
        tProperty( "Youngs", std::to_string( gShearYoungsModulus ), "Shear_Constant" );
        tProperty( "Poisson", "0.0", "Shear_Constant" );
        tProperty( "Fixed", "0.0;0.0", "Shear_Constant" );
        tProperty( "Driven", "0.0", "Shear_Top_Displacement" );
        tProperty( "WestLoad", "-1.0", "Shear_Lateral_Traction" );
        tProperty( "EastLoad", "1.0", "Shear_Lateral_Traction" );

        for ( const std::string tSide : { "Bottom", "Top" } )
        {
            auto tCM = prm::create_constitutive_model_parameter_list();
            tCM.set( "constitutive_name", "Elastic" + tSide );
            tCM.set( "phase_name", tSide );
            tCM.set( "constitutive_type", fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS );
            tCM.set( "model_type", fem::Model_Type::PLANE_STRAIN );
            tCM.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
            tCM.set( "properties", "Youngs,YoungsModulus;Poisson,PoissonRatio" );
            aLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( tCM );

            for ( const std::string tKind : { "Boundary", "Contact" } )
            {
                auto tSP = prm::create_stabilization_parameter_parameter_list();
                tSP.set( "stabilization_name", tKind + tSide );
                tSP.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
                tSP.set( "function_parameters", tKind == "Boundary" ? "100.0" : "10.0" );
                tSP.set( "leader_phase_name", tSide );
                tSP.set( "leader_properties", "Youngs,Material" );
                aLists( FEM::STABILIZATION ).add_parameter_list( tSP );
            }

            auto tBulk = prm::create_IWG_parameter_list();
            tBulk.set( "IWG_name", "Bulk" + tSide );
            tBulk.set( "IWG_type", fem::IWG_Type::STRUC_NON_LINEAR_BULK_PF );
            tBulk.set( "IWG_bulk_type", fem::Element_Type::BULK );
            tBulk.set( "dof_residual", "UX,UY" );
            tBulk.set( "leader_phase_name", tSide );
            tBulk.set( "leader_dof_dependencies", "UX,UY" );
            tBulk.set( "leader_constitutive_models", "Elastic" + tSide + ",ElastLinIso" );
            aLists( FEM::IWG ).add_parameter_list( tBulk );

            auto tBoundary = prm::create_IWG_parameter_list();
            tBoundary.set( "IWG_name", "Dirichlet" + tSide );
            tBoundary.set( "IWG_type", fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE_PF );
            tBoundary.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
            tBoundary.set( "dof_residual", "UX,UY" );
            tBoundary.set( "leader_phase_name", tSide );
            tBoundary.set( "side_ordinals", tSide == "Top" ? "3" : "1" );
            tBoundary.set( "leader_dof_dependencies", "UX,UY" );
            tBoundary.set( "leader_properties", tSide == "Top" ? "Driven,Dirichlet" : "Fixed,Dirichlet" );
            tBoundary.set( "leader_constitutive_models", "Elastic" + tSide + ",ElastLinIso" );
            tBoundary.set( "stabilization_parameters", "Boundary" + tSide + ",DirichletNitsche" );
            aLists( FEM::IWG ).add_parameter_list( tBoundary );

            for ( const std::string tEdge : { "West", "East" } )
            {
                auto tLoad = prm::create_IWG_parameter_list();
                tLoad.set( "IWG_name", tEdge + tSide );
                tLoad.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
                tLoad.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
                tLoad.set( "dof_residual", "UX,UY" );
                tLoad.set( "leader_phase_name", tSide );
                tLoad.set( "side_ordinals", tEdge == "West" ? "4" : "2" );
                tLoad.set( "leader_dof_dependencies", "UX,UY" );
                tLoad.set( "leader_properties", tEdge + "Load,Traction" );
                aLists( FEM::IWG ).add_parameter_list( tLoad );
            }

            const std::string tOther   = tSide == "Top" ? "Bottom" : "Top";
            auto              tContact = prm::create_IWG_parameter_list();
            tContact.set( "IWG_name", "Friction" + tSide );
            tContact.set( "IWG_type", fem::IWG_Type::STRUC_NONLINEAR_CONTACT_SEITZ_UNBIASED_UNSYMMETRIC );
            tContact.set( "IWG_bulk_type", fem::Element_Type::NONCONFORMAL_SIDESET );
            tContact.set( "dof_residual", "UX,UY" );
            tContact.set( "leader_phase_name", tSide );
            tContact.set( "neighbor_phases", "Void" );
            tContact.set( "follower_phase_name", tOther );
            tContact.set( "leader_dof_dependencies", "UX,UY" );
            tContact.set( "follower_dof_dependencies", "UX,UY" );
            tContact.set( "leader_constitutive_models", "Elastic" + tSide + ",ElastLinIso" );
            tContact.set( "follower_constitutive_models", "Elastic" + tOther + ",ElastLinIso" );
            tContact.set( "stabilization_parameters", "Contact" + tSide + ",NitscheInterface" );
            tContact.set( "function_parameters", std::to_string( gShearFriction ) + "/" + std::to_string( gShearPenalty ) + "/0/0" );
            tContact.set( "analytical_jacobian", true );
            aLists( FEM::IWG ).add_parameter_list( tContact );

            for ( uint tComponent = 0; tComponent < 2; ++tComponent )
            {
                auto tQI = prm::create_IQI_parameter_list();
                tQI.set( "IQI_name", tSide + ( tComponent == 0 ? "UX" : "UY" ) );
                tQI.set( "IQI_type", fem::IQI_Type::DOF );
                tQI.set( "IQI_bulk_type", fem::Element_Type::BULK );
                tQI.set( "leader_phase_name", tSide );
                tQI.set( "dof_quantity", "UX,UY" );
                tQI.set( "leader_dof_dependencies", "UX,UY" );
                tQI.set( "vectorial_field_index", static_cast< sint >( tComponent ) );
                aLists( FEM::IQI ).add_parameter_list( tQI );
            }
        }

        for ( uint tComponent = 0; tComponent < 2; ++tComponent )
        {
            auto tQI = prm::create_IQI_parameter_list();
            tQI.set( "IQI_name", tComponent == 0 ? "ShearForce" : "NormalForce" );
            tQI.set( "IQI_type", fem::IQI_Type::TRACTION );
            tQI.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
            tQI.set( "leader_phase_name", "Top" );
            tQI.set( "side_ordinals", "3" );
            tQI.set( "leader_dof_dependencies", "UX,UY" );
            tQI.set( "leader_constitutive_models", "ElasticTop,TractionCM" );
            tQI.set( "vectorial_field_index", static_cast< sint >( tComponent ) );
            aLists( FEM::IQI ).add_parameter_list( tQI );
        }

        aLists( FEM::COMPUTATION ).set( "nonconformal_integration_order", static_cast< uint >( mtk::Integration_Order::BAR_4 ) );
        aLists.set( "nonconformal_max_negative_ray_length", -0.01 );
        aLists.set( "nonconformal_max_positive_ray_length", 0.05 );
    }

    void SOLParameterList( Module_Parameter_Lists& aLists )
    {
        auto tLinear = prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
        tLinear.set( "Solver_Type", "Amesos_Umfpack" );
        aLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( tLinear );
        aLists( SOL::LINEAR_SOLVERS ).add_parameter_list( prm::create_linear_solver_parameter_list() );

        auto tNewton = prm::create_nonlinear_algorithm_parameter_list();
        tNewton.set( "NLA_combined_res_jac_assembly", true );
        tNewton.set( "NLA_max_iter", 60 );
        tNewton.set( "NLA_hard_break", true );
        tNewton.set( "NLA_rel_res_norm_drop", 1.0e-11 );
        tNewton.set( "NLA_tot_res_norm", 1.0e-11 );
        // Sliding is much smaller than a facet: retain the initial candidate
        // facets and update the actual ray coordinates locally in every IWG call.
        // This patch isolates the friction law from global search/facet switching.
        tNewton.set( "NLA_remap_strategy", static_cast< uint >( sol::SolverRaytracingStrategy::None ) );
        tNewton.set( "NLA_remap_iteration_frequency", 1 );
        aLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list( tNewton );
        auto tNonlinear = prm::create_nonlinear_solver_parameter_list();
        tNonlinear.set( "NLA_DofTypes", "UX,UY" );
        aLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( tNonlinear );

        auto tTime = prm::create_time_solver_algorithm_parameter_list();
        tTime.set( "TSA_Num_Time_Steps", static_cast< sint >( tStages * gShearSubsteps ) );
        tTime.set( "TSA_Time_Frame", static_cast< real >( tStages * gShearSubsteps ) );
        aLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list( tTime );
        auto tSolver = prm::create_time_solver_parameter_list();
        tSolver.set( "TSA_DofTypes", "UX,UY" );
        tSolver.set( "TSA_Output_Indices", "0" );
        tSolver.set( "TSA_Output_Criteria", "Shear_Output" );
        aLists( SOL::TIME_SOLVERS ).add_parameter_list( tSolver );
        aLists( SOL::PRECONDITIONERS ).add_parameter_list( prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void MSIParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "UX", 0 );
        aLists.set( "UY", 0 );
    }

    void VISParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Coulomb_Shear.exo" ) );
        aLists.set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        aLists.set( "Set_Names", "HMR_dummy_n_p2,HMR_dummy_c_p2,HMR_dummy_n_p1,HMR_dummy_c_p1,SideSet_3_n_p1,SideSet_3_c_p1" );
        aLists.set( "Field_Names", "TOP_UX,TOP_UY,BOTTOM_UX,BOTTOM_UY,SHEAR_FORCE,NORMAL_FORCE" );
        aLists.set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,GLOBAL,GLOBAL" );
        aLists.set( "IQI_Names", "TopUX,TopUY,BottomUX,BottomUY,ShearForce,NormalForce" );
        aLists.set( "Save_Frequency", 1 );
        aLists.set( "Time_Offset", 0.0 );
    }
}    // namespace moris
}

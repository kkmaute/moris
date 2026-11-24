/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#include <cmath>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"
#include <string>
#include "parameters.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "fn_trans.hpp"

#ifdef MORIS_DG_CANTILEVER_DRIVER
extern "C" {
extern bool        gDGNonlinear;
extern moris::uint gDGTimeSteps;
extern moris::real gDGFinalTime;
extern moris::real gDGBeamHeight;
extern moris::real gDGDensity;
extern moris::real gDGAmplitude;
}
#else
// Standalone external inputs supply defaults without EXA driver symbols.
bool        gDGNonlinear  = false;
moris::uint gDGTimeSteps  = 6;
moris::real gDGFinalTime  = 0.3;
moris::real gDGBeamHeight = 0.43;
moris::real gDGDensity    = 1.0;
moris::real gDGAmplitude  = 0.025;
#endif

extern "C" {
namespace moris
{
    constexpr real    tLength   = 3.0;
    constexpr real    tMu       = 3.0;
    constexpr real    tLambda   = 5.0;
    const std::string tBulk     = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    const std::string tClamp    = "SideSet_4_n_p0,SideSet_4_c_p0";
    const std::string tTraction = "SideSet_1_n_p0,SideSet_1_c_p0,SideSet_2_n_p0,SideSet_2_c_p0,iside_b0_0_b1_1";

    void Constant( Matrix< DDRMat >& aValue, Vector< Matrix< DDRMat > >& aParameters, fem::Field_Interpolator_Manager* )
    {
        aValue = aParameters( 0 );
    }

    void InitialVelocity( Matrix< DDRMat >& aValue, Vector< Matrix< DDRMat > >&, fem::Field_Interpolator_Manager* aManager )
    {
        const real x = aManager->get_IP_geometry_interpolator()->valx()( 0 );
        aValue       = { { 0.0 }, { -gDGAmplitude * x * x } };
    }

    void BodyLoad( Matrix< DDRMat >& aValue, Vector< Matrix< DDRMat > >&, fem::Field_Interpolator_Manager* aManager )
    {
        const real t = aManager->get_IP_geometry_interpolator()->valt()( 0 );
        // Manufactured u=(0,-a*x^2*t), with zero acceleration and det(F)=1.
        const real x = aManager->get_IP_geometry_interpolator()->valx()( 0 );
        const real k = -2.0 * gDGAmplitude * x * t;
        const real c = 0.5 * tLambda + tMu;
        aValue       = { { gDGNonlinear ? 4.0 * c * gDGAmplitude * t * k : 0.0 },
                  { 2.0 * gDGAmplitude * t * ( tMu + ( gDGNonlinear ? 3.0 * c * k * k : 0.0 ) ) } };
    }

    void ManufacturedTraction( Matrix< DDRMat >& aValue, Vector< Matrix< DDRMat > >&, fem::Field_Interpolator_Manager* aManager )
    {
        const real       x = aManager->get_IP_geometry_interpolator()->valx()( 0 );
        const real       t = aManager->get_IP_geometry_interpolator()->valt()( 0 );
        const real       k = -2.0 * gDGAmplitude * x * t;
        const real       c = 0.5 * tLambda + tMu;
        Matrix< DDRMat > P = { { 0.0, tMu * k }, { tMu * k, 0.0 } };
        if ( gDGNonlinear )
        {
            P( 0, 0 ) = c * k * k;
            P( 1, 0 ) += c * k * k * k;
            P( 1, 1 ) = c * k * k;
        }
        aValue = P * aManager->get_IG_geometry_interpolator()->get_normal();
    }

    bool OutputCriterion( tsa::Time_Solver* ) { return true; }

    void OPTParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "is_optimization_problem", false );
    }

    void HMRParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "number_of_elements_per_dimension", Vector< uint >{ 12, 4 } );
        aLists.set( "domain_dimensions", Vector< real >{ tLength, 0.8 } );
        aLists.set( "domain_offset", Vector< real >{ 0.0, 0.0 } );
        aLists.set( "lagrange_output_meshes", "0" );
        aLists.set( "lagrange_orders", "2" );
        aLists.set( "lagrange_pattern", "0" );
        aLists.set( "bspline_orders", "2" );
        aLists.set( "bspline_pattern", "0" );
        aLists.set( "lagrange_to_bspline", "0" );
    }

    void XTKParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "decompose", true );
        aLists.set( "decomposition_type", "conformal" );
        aLists.set( "enrich_mesh_indices", "0" );
        aLists.set( "ghost_stab", true );
        aLists.set( "high_to_low_dbl_side_sets", true );
    }

    void GENParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "IQI_types", "FinalUY", "FinalVolume" );
        aLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aLists.set( "center_x", 0.0 );
        aLists.set( "center_y", gDGBeamHeight );

        aLists.set( "normal_x", 0.0 );
        aLists.set( "normal_y", 1.0 );
    }

    void FEMParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.hack_for_legacy_fem();
        auto property = [ & ]( const std::string& aName, const std::string& aValue, const std::string& aFunction = "Constant" ) {
            aLists( FEM::PROPERTIES ).add_parameter_list();
            aLists.set( "property_name", aName );
            aLists.set( "function_parameters", aValue );
            aLists.set( "value_function", aFunction );
        };
        property( "Youngs", "7.875" );
        property( "Poisson", "0.3125" );
        property( "Zero", "0.0;0.0" );
        property( "InitialV", "0.0;0.0", "InitialVelocity" );
        property( "Body", "0.0;0.0", "BodyLoad" );
        property( "Traction", "0.0;0.0", "ManufacturedTraction" );
        property( "Rho", std::to_string( gDGDensity ), "Constant" );

        aLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aLists.set( "constitutive_name", "Elasticity" );
        aLists.set( "constitutive_type", gDGNonlinear ? fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF : fem::Constitutive_Type::STRUC_LIN_ISO );
        aLists.set( "model_type", fem::Model_Type::PLANE_STRAIN );
        aLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aLists.set( "properties", "Youngs,YoungsModulus;Poisson,PoissonRatio" );

        aLists( FEM::STABILIZATION ).add_parameter_list();
        aLists.set( "stabilization_name", "ClampPenalty" );
        aLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aLists.set( "function_parameters", "100.0" );
        aLists.set( "leader_properties", "Youngs,Material" );
        aLists( FEM::STABILIZATION ).add_parameter_list();
        aLists.set( "stabilization_name", "GhostPenalty" );
        aLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aLists.set( "function_parameters", "0.01" );
        aLists.set( "leader_properties", "Youngs,Material" );

        auto iwg = [ & ]( const std::string& aName, fem::IWG_Type aType, const std::string& aSets, bool aSpatial ) {
            aLists( FEM::IWG ).add_parameter_list();
            aLists.set( "IWG_name", aName );
            aLists.set( "IWG_type", aType );
            aLists.set( "dof_residual", "UX,UY" );
            aLists.set( "leader_dof_dependencies", "UX,UY" );
            aLists.set( "mesh_set_names", aSets );
            aLists.set( "time_dg", aSpatial );
        };
        iwg( "Bulk", gDGNonlinear ? fem::IWG_Type::STRUC_NON_LINEAR_BULK_PF : fem::IWG_Type::STRUC_LINEAR_BULK, tBulk, true );
        aLists.set( "leader_constitutive_models", "Elasticity,ElastLinIso" );
        aLists.set( "leader_properties", "Body,Load" );
        iwg( "Clamp", gDGNonlinear ? fem::IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF : fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE, tClamp, true );
        aLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aLists.set( "leader_constitutive_models", "Elasticity,ElastLinIso" );
        aLists.set( "leader_properties", "Zero,Dirichlet" );
        aLists.set( "stabilization_parameters", "ClampPenalty,DirichletNitsche" );
        iwg( "Tractions", fem::IWG_Type::STRUC_LINEAR_NEUMANN, tTraction, true );
        aLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aLists.set( "leader_properties", "Traction,Traction" );

        iwg( "Ghost", fem::IWG_Type::GHOST_NORMAL_FIELD, "ghost_p0", true );
        aLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
        aLists.set( "follower_dof_dependencies", "UX,UY" );
        aLists.set( "stabilization_parameters", "GhostPenalty,GhostSP" );
        aLists.set( "ghost_order", 2 );

        iwg( "Inertia", fem::IWG_Type::TIME_DG_MASS, tBulk, false );
        aLists.set( "time_continuity", true );
        aLists.set( "leader_properties", "Rho,Density;Zero,InitialDisplacement;InitialV,InitialVelocity" );

        for ( uint c = 0; c < 2; ++c )
        {
            aLists( FEM::IQI ).add_parameter_list();
            aLists.set( "IQI_name", c == 0 ? "UX" : "UY" );
            aLists.set( "IQI_type", fem::IQI_Type::DOF );
            aLists.set( "dof_quantity", "UX,UY" );
            aLists.set( "leader_dof_dependencies", "UX,UY" );
            aLists.set( "vectorial_field_index", static_cast< sint >( c ) );
            aLists.set( "mesh_set_names", tBulk );
        }
        aLists( FEM::IQI ).add_parameter_list();
        aLists.set( "IQI_name", "FinalUY" );
        aLists.set( "IQI_type", fem::IQI_Type::DOF );
        aLists.set( "dof_quantity", "UX,UY" );
        aLists.set( "leader_dof_dependencies", "UX,UY" );
        aLists.set( "vectorial_field_index", 1 );
        aLists.set( "time_boundary", true );
        aLists.set( "time_final", gDGFinalTime );
        aLists.set( "mesh_set_names", tBulk );
        aLists( FEM::IQI ).add_parameter_list();
        aLists.set( "IQI_name", "FinalVolume" );
        aLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aLists.set( "time_boundary", true );
        aLists.set( "time_final", gDGFinalTime );
        aLists.set( "mesh_set_names", tBulk );
        aLists( FEM::COMPUTATION );
        aLists.set( "is_analytical_forward", true );
    }

    void SOLParameterList( Module_Parameter_Lists& aLists )
    {
        aLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
        aLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aLists.set( "NLA_combined_res_jac_assembly", false );
        aLists.set( "NLA_rel_res_norm_drop", 1e-11 );
        aLists.set( "NLA_tot_res_norm", 1e-12 );
        aLists.set( "NLA_max_iter", 30 );
        aLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aLists.set( "NLA_DofTypes", "UX,UY" );
        aLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aLists.set( "TSA_Num_Time_Steps", static_cast< sint >( gDGTimeSteps ) );
        aLists.set( "TSA_Time_Frame", gDGFinalTime );
        aLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aLists.set( "TSA_DofTypes", "UX,UY" );
        aLists.set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0" );
        aLists.set( "TSA_time_level_per_type", "UX,2;UY,2" );
        aLists.set( "TSA_Output_Indices", "0" );
        aLists.set( "TSA_Output_Criteria", "OutputCriterion" );
        aLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    void MSIParameterList( Module_Parameter_Lists& aLists ) { aLists.set( "order_adofs_by_host", false ); }
    void VISParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Time_DG_Cantilever_Forward.exo" ) );
        aLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aLists.set( "Set_Names", tBulk );
        aLists.set( "Field_Names", "UX,UY,FinalUY,FinalVolume" );
        aLists.set( "Field_Type", "NODAL,NODAL,GLOBAL,GLOBAL" );
        aLists.set( "IQI_Names", "UX,UY,FinalUY,FinalVolume" );
        aLists.set( "Save_Frequency", 1 );
    }
    void MORISGENERALParameterList( Module_Parameter_Lists& ) {}
}    // namespace moris
}

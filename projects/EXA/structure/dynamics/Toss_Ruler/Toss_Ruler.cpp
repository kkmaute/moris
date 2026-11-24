/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#include <algorithm>
#include <cmath>
#include <string>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"
#include "parameters.hpp"
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "fn_trans.hpp"

#ifdef MORIS_TOSS_DRIVER
extern "C" {
extern moris::uint gTossElementsX;
extern moris::uint gTossElementsY;
extern moris::uint gTossTimeSteps;
extern moris::real gTossFinalTime;
}
#else
// The standalone input runs the full flight without EXA driver symbols.
moris::uint gTossElementsX = 30;
moris::uint gTossElementsY = 4;
moris::uint gTossTimeSteps = 800;
moris::real gTossFinalTime = 0.04;
#endif

extern "C" {
namespace moris
{
    const std::string tBulk          = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    constexpr real    tLength        = 0.3;
    constexpr real    tHeight        = 0.06;
    constexpr real    tThickness     = 0.002;
    constexpr real    tPulseDuration = 0.004;
    constexpr real    tPeakLineLoad  = 4.0e4;

    void Constant( Matrix< DDRMat >& aValue, Vector< Matrix< DDRMat > >& aParameters, fem::Field_Interpolator_Manager* )
    {
        aValue = aParameters( 0 );
    }

    // Loads are dead nominal tractions on the reference end faces. The line
    // resultant R(t) [N/m of extrusion] is converted to traction by dividing by
    // the in-plane height; the Neumann IWG applies the extrusion thickness.
    void TossTraction( Matrix< DDRMat >& aValue, Vector< Matrix< DDRMat > >& aParameters, fem::Field_Interpolator_Manager* aManager )
    {
        const real t          = aManager->get_IP_geometry_interpolator()->valt()( 0 );
        const real tAmplitude = std::max( 0.0, 1.0 - std::abs( 2.0 * t / tPulseDuration - 1.0 ) );
        aValue                = { { 0.0 }, { aParameters( 0 )( 0 ) * tPeakLineLoad * tAmplitude / tHeight } };
    }

    // Export the incoming displacement trace alongside the outgoing one.
    // Together with the preceding outgoing trace it permits reconstruction of
    // the committed endpoint velocity; gradt(u) is not that velocity.
    void IncomingDisplacement( Matrix< DDRMat >& aValue, Vector< Matrix< DDRMat > >&, fem::Field_Interpolator_Manager* aManager )
    {
        auto*                  tFI    = aManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
        const uint             tSpace = tFI->get_number_of_space_bases();
        const Matrix< DDRMat > tShape = tFI->NBuild()( { 0, 0 }, { 0, tSpace - 1 } )
                                      + tFI->NBuild()( { 0, 0 }, { tSpace, 2 * tSpace - 1 } );
        aValue = trans( tShape * tFI->get_coeff()( { 0, tSpace - 1 }, { 0, 1 } ) );
    }

    bool OutputCriterion( tsa::Time_Solver* ) { return true; }

    void OPTParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "is_optimization_problem", false );
    }

    void HMRParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "number_of_elements_per_dimension", Vector< uint >{ gTossElementsX, gTossElementsY } );
        aLists.set( "domain_dimensions", Vector< real >{ tLength, tHeight } );
        aLists.set( "domain_offset", Vector< real >{ 0.0, -0.5 * tHeight } );
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
        aLists.set( "ghost_stab", false );
    }

    void GENParameterList( Module_Parameter_Lists& aLists )
    {
        // A line outside the fitted beam leaves a single active phase. No cut
        // cells, ghost penalty, supports or artificial springs affect flight.
        aLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aLists.set( "center_x", 0.0 );
        aLists.set( "center_y", 2.0 * tHeight );
        aLists.set( "normal_x", 0.0 );
        aLists.set( "normal_y", 1.0 );
    }

    void FEMParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.hack_for_legacy_fem();
        auto tProperty = [ & ]( const std::string& aName, const std::string& aParameters, const std::string& aFunction = "Constant" ) {
            aLists( FEM::PROPERTIES ).add_parameter_list();
            aLists.set( "property_name", aName );
            aLists.set( "function_parameters", aParameters );
            aLists.set( "value_function", aFunction );
        };
        tProperty( "Youngs", "2.06e11" );
        tProperty( "Poisson", "0.0" );
        tProperty( "Density", "7800.0" );
        tProperty( "Thickness", std::to_string( tThickness ) );
        tProperty( "Zero", "0.0;0.0" );
        tProperty( "LeftPulse", "1.0", "TossTraction" );
        tProperty( "RightPulse", "-0.5", "TossTraction" );
        tProperty( "Incoming", "0.0", "IncomingDisplacement" );

        aLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aLists.set( "constitutive_name", "Elasticity" );
        aLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF );
        aLists.set( "model_type", fem::Model_Type::PLANE_STRESS );
        aLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aLists.set( "properties", "Youngs,YoungsModulus;Poisson,PoissonRatio" );

        auto tIWG = [ & ]( const std::string& aName, fem::IWG_Type aType, const std::string& aSets, bool aSpatial ) {
            aLists( FEM::IWG ).add_parameter_list();
            aLists.set( "IWG_name", aName );
            aLists.set( "IWG_type", aType );
            aLists.set( "dof_residual", "UX,UY" );
            aLists.set( "leader_dof_dependencies", "UX,UY" );
            aLists.set( "mesh_set_names", aSets );
            aLists.set( "time_dg", aSpatial );
        };
        tIWG( "Bulk", fem::IWG_Type::STRUC_NON_LINEAR_BULK_PF, tBulk, true );
        aLists.set( "leader_constitutive_models", "Elasticity,ElastLinIso" );
        aLists.set( "leader_properties", "Thickness,Thickness" );

        tIWG( "LeftEnd", fem::IWG_Type::STRUC_LINEAR_NEUMANN, "SideSet_4_n_p0,SideSet_4_c_p0", true );
        aLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aLists.set( "leader_properties", "LeftPulse,Traction;Thickness,Thickness" );
        tIWG( "RightEnd", fem::IWG_Type::STRUC_LINEAR_NEUMANN, "SideSet_2_n_p0,SideSet_2_c_p0", true );
        aLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aLists.set( "leader_properties", "RightPulse,Traction;Thickness,Thickness" );

        // Consistent mass makes the free-body transient operator invertible.
        tIWG( "Inertia", fem::IWG_Type::TIME_DG_MASS, tBulk, false );
        aLists.set( "time_continuity", true );
        aLists.set( "leader_properties", "Density,Density;Thickness,Thickness;Zero,InitialDisplacement;Zero,InitialVelocity" );

        for ( uint c = 0; c < 2; ++c )
        {
            aLists( FEM::IQI ).add_parameter_list();
            aLists.set( "IQI_name", c == 0 ? "UX" : "UY" );
            aLists.set( "IQI_type", fem::IQI_Type::DOF );
            aLists.set( "dof_quantity", "UX,UY" );
            aLists.set( "leader_dof_dependencies", "UX,UY" );
            aLists.set( "vectorial_field_index", static_cast< sint >( c ) );
            aLists.set( "mesh_set_names", tBulk );

            aLists( FEM::IQI ).add_parameter_list();
            aLists.set( "IQI_name", c == 0 ? "IncomingUX" : "IncomingUY" );
            aLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
            aLists.set( "leader_properties", "Incoming,Property" );
            aLists.set( "leader_dof_dependencies", "UX,UY" );
            aLists.set( "vectorial_field_index", static_cast< sint >( c ) );
            aLists.set( "mesh_set_names", tBulk );
        }
        aLists( FEM::IQI ).add_parameter_list();
        aLists.set( "IQI_name", "FinalArea" );
        aLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        // The endpoint-only set still needs the displacement interpolation map
        // when the solver traverses all equation sets during assembly.
        aLists.set( "leader_dof_dependencies", "UX,UY" );
        aLists.set( "time_boundary", true );
        aLists.set( "time_final", gTossFinalTime );
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
        aLists.set( "NLA_rel_res_norm_drop", 1.0e-9 );
        aLists.set( "NLA_tot_res_norm", 1.0e-10 );
        aLists.set( "NLA_max_iter", 30 );
        aLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aLists.set( "NLA_DofTypes", "UX,UY" );
        aLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aLists.set( "TSA_Num_Time_Steps", static_cast< sint >( gTossTimeSteps ) );
        aLists.set( "TSA_Time_Frame", gTossFinalTime );
        aLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aLists.set( "TSA_DofTypes", "UX,UY" );
        aLists.set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0" );
        aLists.set( "TSA_time_level_per_type", "UX,2;UY,2" );
        aLists.set( "TSA_Output_Indices", "0" );
        aLists.set( "TSA_Output_Criteria", "OutputCriterion" );
        aLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    void MSIParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "order_adofs_by_host", false );
    }

    void VISParameterList( Module_Parameter_Lists& aLists )
    {
        aLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Toss_Ruler.exo" ) );
        aLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aLists.set( "Set_Names", tBulk );
        aLists.set( "Field_Names", "UX,UY,IncomingUX,IncomingUY,FinalArea" );
        aLists.set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,GLOBAL" );
        aLists.set( "IQI_Names", "UX,UY,IncomingUX,IncomingUY,FinalArea" );
        // Every slab is needed for reconstructed endpoint velocity diagnostics.
        aLists.set( "Save_Frequency", 1 );
    }
}    // namespace moris
}    // extern "C"

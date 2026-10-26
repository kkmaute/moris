/* Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. */

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

#include "parameters.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_TSA_Time_Solver.hpp"

// Test and benchmark settings are defined by the example driver, following the
// other EXA input libraries. The driver sets them before invoking the workflow.
extern moris::uint gElementsPerSubdomain;
extern bool        gUseMaterialJump;
extern moris::real gYoungsModulusRatio;
extern moris::real gLinearTolerance;
extern std::string gLinearSolver;
extern std::string gOutputFileName;
extern bool        gUseEdges;
extern std::string gFactorSolverType;
extern std::string gPetscOptions;

// Parameter callbacks exported for the standard HMR/XTK/FEM workflow.
extern "C" {
namespace moris
{
    constexpr uint tDimension = 2;

    // ----------------------------------------------------------------------------

    void Func_Const(
            Matrix< DDRMat >&           aValue,
            Vector< Matrix< DDRMat > >& aParameters,
            fem::Field_Interpolator_Manager* )
    {
        aValue = aParameters( 0 );
    }

    // ----------------------------------------------------------------------------

    void Func_Youngs(
            Matrix< DDRMat >&                aValue,
            Vector< Matrix< DDRMat > >&      aParameters,
            fem::Field_Interpolator_Manager* aFIManager )
    {
        const auto tCoordinates = aFIManager->get_IP_geometry_interpolator()->valx();
        bool       tInside      = true;

        // Material jumps occupy the central square/cube of the unit domain.
        for ( uint d = 0; d < tDimension; ++d )
        {
            tInside = tInside && tCoordinates( d ) > 0.25 && tCoordinates( d ) < 0.75;
        }

        aValue.set_size( 1, 1, tInside ? aParameters( 0 )( 1 ) : aParameters( 0 )( 0 ) );
    }

    // ----------------------------------------------------------------------------

    bool Output_Criterion( tsa::Time_Solver* )
    {
        return true;
    }

    // ----------------------------------------------------------------------------

    void OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
        aParameterLists.set( "workflow", "HMR_XTK" );
    }

    // ----------------------------------------------------------------------------

    void HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // Match HMR's Cartesian rank decomposition; the global controls the
        // number of Q1 elements along each subdomain direction.
        int tGrid[ 3 ] = { 0, 0, 0 };
        MPI_Dims_create( par_size(), tDimension, tGrid );

        Vector< uint > tNumElements( tDimension );

        for ( uint d = 0; d < tDimension; ++d )
        {
            MORIS_ERROR( gElementsPerSubdomain > 0
                                 && gElementsPerSubdomain <= uint( std::numeric_limits< sint >::max() ) / uint( tGrid[ d ] ),
                    "Invalid number of elements per subdomain." );
            tNumElements( d ) = gElementsPerSubdomain * tGrid[ d ];

            MORIS_ERROR( !gUseMaterialJump || tNumElements( d ) % 4 == 0,
                    "Material inclusion interfaces must coincide with HMR element boundaries." );
        }

        aParameterLists.set( "number_of_elements_per_dimension", tNumElements );
        aParameterLists.set( "domain_dimensions", Vector< real >( tDimension, 1.0 ) );
        aParameterLists.set( "processor_decomposition_method", 1 );

        // Match first-order interpolation and B-spline meshes without refinement.
        aParameterLists.set( "lagrange_output_meshes", "0" );
        aParameterLists.set( "lagrange_orders", "1" );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", "1" );
        aParameterLists.set( "bspline_pattern", "0" );
        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "refinement_buffer", 0 );
        aParameterLists.set( "staircase_buffer", 0 );
    }

    // ----------------------------------------------------------------------------

    void XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // A geometry outside the domain gives a single uncut phase. XTK supplies
        // the usual enriched mesh/set interfaces while preserving QUAD4/HEX8 cells.
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "triangulate_all", false );
        aParameterLists.set( "ghost_stab", false );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", false );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", false );
    }

    // ----------------------------------------------------------------------------

    real Outside_Domain( const Matrix< DDRMat >&, const Vector< real >& )
    {
        return 1.0;
    }

    // ----------------------------------------------------------------------------

    void GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Outside_Domain" );
    }

    // ----------------------------------------------------------------------------

    void FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        const std::string tDofs  = tDimension == 2 ? "UX,UY" : "UX,UY,UZ";
        const std::string tBulk  = "HMR_dummy_n_p1";
        const std::string tLeft  = "SideSet_4_n_p1";
        const std::string tRight = "SideSet_2_n_p1";

        aParameterLists.hack_for_legacy_fem();

        MORIS_ERROR( std::isfinite( gYoungsModulusRatio ) && gYoungsModulusRatio > 0.0,
                "The inclusion modulus must be positive and finite." );

        // Material properties are evaluated by FEM on the HMR interpolation mesh.
        std::ostringstream tMaterial;
        tMaterial << std::setprecision( 17 ) << "1.0," << gYoungsModulusRatio;

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", gUseMaterialJump ? tMaterial.str() : "3.0e7,3.0e7" );
        aParameterLists.set( "value_function", "Func_Youngs" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", "0.3" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichlet" );
        aParameterLists.set( "function_parameters", tDimension == 2 ? "0.0;0.0" : "0.0;0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // Apply unit horizontal traction, independent of mesh resolution.
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropTraction" );
        aParameterLists.set( "function_parameters", tDimension == 2 ? "1.0;0.0" : "1.0;0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // Integrate the traction magnitude to verify the total applied force.
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropLoadMagnitude" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // Linear isotropic plane stress / full 3D elasticity.
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMElasticity" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "model_type", tDimension == 2 ? fem::Model_Type::PLANE_STRESS : fem::Model_Type::FULL );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( tDofs, "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        // Use the standard symmetric Nitsche clamp. The MSI workflow currently
        // exposes no strongly constrained DOF list; this distinction from the
        // paper's strong clamp is documented with the example's benchmark results.
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPClamp" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        // Assemble interior elastic forces and tangent contributions.
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGElasticity" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", tDofs );
        aParameterLists.set( "leader_dof_dependencies", tDofs );
        aParameterLists.set( "leader_constitutive_models", "CMElasticity,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tBulk );

        // Impose zero displacement on the left boundary through symmetric Nitsche.
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGClamp" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", tDofs );
        aParameterLists.set( "leader_dof_dependencies", tDofs );
        aParameterLists.set( "leader_properties", "PropDirichlet,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMElasticity,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPClamp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tLeft );

        // Apply horizontal loading on the opposite boundary.
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTraction" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", tDofs );
        aParameterLists.set( "leader_dof_dependencies", tDofs );
        aParameterLists.set( "leader_properties", "PropTraction,Traction" );
        aParameterLists.set( "mesh_set_names", tRight );

        // Nodal fields and global IQIs are written through the normal VIS model.
        const std::string tNames[ 3 ] = { "UX", "UY", "UZ" };

        for ( uint d = 0; d < tDimension; ++d )
        {
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "IQI" + tNames[ d ] );
            aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
            aParameterLists.set( "dof_quantity", tDofs );
            aParameterLists.set( "leader_dof_dependencies", tDofs );
            aParameterLists.set( "vectorial_field_index", sint( d ) );
            aParameterLists.set( "mesh_set_names", tBulk );
        }

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIEnergy" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", tDofs );
        aParameterLists.set( "leader_constitutive_models", "CMElasticity,Elast" );
        aParameterLists.set( "mesh_set_names", tBulk );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIVolume" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "mesh_set_names", tBulk );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQILoad" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropLoadMagnitude,Property" );
        aParameterLists.set( "mesh_set_names", tRight );

        aParameterLists( FEM::COMPUTATION );
    }

    // ----------------------------------------------------------------------------

    void SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {
        const std::string tDofs = tDimension == 2 ? "UX,UY" : "UX,UY,UZ";

        MORIS_ERROR( gLinearSolver == "bddc" || gLinearSolver == "aij" || gLinearSolver == "amesos",
                "Unknown BDDC example solver." );
        MORIS_ERROR( std::isfinite( gLinearTolerance ) && gLinearTolerance > 0.0 && gLinearTolerance < 1.0,
                "The linear relative tolerance must be between zero and one." );

        // Use the same FEM model for the direct reference and PETSc comparisons.
        if ( gLinearSolver == "amesos" )
        {
            aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
        }
        else
        {
            aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::PETSC );
            aParameterLists.set( "KSPType", gLinearSolver == "bddc" ? "fgmres" : "preonly" );
            aParameterLists.set( "KSPTol", gLinearTolerance );
            aParameterLists.set( "KSPMaxits", 2000 );
            aParameterLists.set( "KSPMGMRESRestart", 1000 );
            aParameterLists.set( "preconditioners", "0" );
            aParameterLists.set( "KSPErrorIfNotConverged", true );
            aParameterLists.set( "KSPMonitorTrueResidual", true );
            aParameterLists.set( "KSPConvergedReason", true );

            // BDDC needs subdomain element matrices (MATIS); retain AIJ for LU.
            aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_TPL_Type", sol::MapType::Petsc );
            aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_PETSc_Log_File", "petsc.log" );
            aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_PETSc_Options", gPetscOptions );
            aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_PETSc_Mat_Type", gLinearSolver == "bddc" ? "is" : "aij" );
        }

        // Wrap the linear solve in the standard nonlinear/time-solver workflow.
        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_combined_res_jac_assembly", true );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.e-10 );
        aParameterLists.set( "NLA_max_iter", 5 );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", tDofs );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Num_Time_Steps", 1 );
        aParameterLists.set( "TSA_Time_Frame", 1.0 );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", tDofs );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( gLinearSolver == "amesos" ? sol::PreconditionerType::NONE : sol::PreconditionerType::PETSC );
        if ( gLinearSolver == "bddc" )
        {
            aParameterLists.set( "PCType", "bddc" );
            aParameterLists.set( "use_bddc_defaults", false );
            aParameterLists.set( "bddc_symmetric", true );
            aParameterLists.set( "bddc_detect_disconnected", true );
            aParameterLists.set( "bddc_use_vertices", true );
            aParameterLists.set( "bddc_use_edges", gUseEdges );
            aParameterLists.set( "bddc_use_faces", gUseEdges );
            aParameterLists.set( "bddc_use_stiffness_scaling", true );

            // Identify displacement fields and physical subdomain corners for
            // robust primal constraints on floating elasticity subdomains.
            aParameterLists.set( "bddc_use_dof_splitting", true );
            aParameterLists.set( "bddc_use_subdomain_corners", true );

            aParameterLists.set( "bddc_dirichlet_solver", gFactorSolverType );
            aParameterLists.set( "bddc_neumann_solver", gFactorSolverType );
            if ( par_size() > 1 )
            {
                // The serial coarse problem can be empty; select a parallel
                // factor package only for the nontrivial parallel coarse solve.
                aParameterLists.set( "bddc_coarse_solver", gFactorSolverType );
            }
        }
        else if ( gLinearSolver == "aij" )
        {
            aParameterLists.set( "PCType", "lu" );
            aParameterLists.set( "factor_solver_type", gFactorSolverType );
        }
    }

    // ----------------------------------------------------------------------------

    void MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "order_adofs_by_host", true );
    }

    // ----------------------------------------------------------------------------

    void VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", gOutputFileName ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );

        // Include the loaded sideset so VIS evaluates the total-force IQI.
        aParameterLists.set( "Set_Names", "HMR_dummy_n_p1,SideSet_2_n_p1" );
        aParameterLists.set( "Field_Names", tDimension == 2 ? "UX,UY,STRAIN_ENERGY,VOLUME,LOAD" : "UX,UY,UZ,STRAIN_ENERGY,VOLUME,LOAD" );
        aParameterLists.set( "Field_Type", tDimension == 2 ? "NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL" : "NODAL,NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL" );
        aParameterLists.set( "IQI_Names", tDimension == 2 ? "IQIUX,IQIUY,IQIEnergy,IQIVolume,IQILoad" : "IQIUX,IQIUY,IQIUZ,IQIEnergy,IQIVolume,IQILoad" );
        aParameterLists.set( "Save_Frequency", 1 );
    }

    // ----------------------------------------------------------------------------

    void MORISGENERALParameterList( Module_Parameter_Lists& )
    {
    }
}    // namespace moris
}    // extern "C"

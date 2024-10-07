/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * single_element.cpp
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
// #include "cl_DLA_Eigen_Solver.hpp"
#include "parameters.hpp"
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"
#include "fn_stringify_matrix.hpp"

#include "AztecOO.h"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

//------------------------------------------------------------------------------

// interpolation order
extern std::string tOrder;

// spatial dimensions
extern uint tDim;

// stress used in geometric stiffness (PK2 or Linear)
extern std::string tStressType;

// output file
extern std::string tOutputFileName;

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------------

namespace moris
{
    void create_trilinos_parameter_list( Module_Parameter_Lists& );
    void create_petsc_parameter_list( Module_Parameter_Lists& );

    /* ------------------------------------------------------------------------ */
    // General

    // use Petsc or Trilinos eigensolver (Trilinos does not work for this problem yet)
    bool tUsePetsc = true;

    // Traction load
    std::string tTraction = tDim == 3 ? "0.0;-1.0;0.0" : "0.0;-1.0";

    // Prescribed displacements
    std::string tDirichlet = tDim == 3 ? "0.0;0.0;0.0" : "0.0;0.0";

    // Number of eigenvectors
    sint tNumEigenVectors = 5;

    std::string tDofString = tDim == 3 ? "UX,UY,UZ" : "UX,UY";

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    std::string tNumElemsPerDim = tDim == 3 ? "10,55,10" : "10,55";
    std::string tDomainDims     = tDim == 3 ? "2,11,2" : "2,11";
    std::string tDomainOffset   = tDim == 3 ? "-0.5,-0.5,-0.5" : "-0.5,-0.5";
    std::string tDomainSidesets = tDim == 3 ? "1,2,3,4,5,6" : "1,2,3,4";

    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    // Bulk Phases
    std::string tBulk = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    std::string tBulkGhost = "ghost_p0";

    // boundaries
    std::string tBottom = "iside_b0_0_b1_1";
    std::string tTop    = "iside_b0_0_b1_2";

    /* ------------------------------------------------------------------------ */
    // material parameters, kg is scaled with a factor 1e-6

    // Copper material properties
    std::string tYoungsModulus = "100";
    std::string tPoissonRatio  = "0.3";

    /* ------------------------------------------------------------------------ */
    // Output Config
    bool tOutputCriterion = true;

    /* ------------------------------------------------------------------------ */
    // Solver config

    moris::real tNLA_rel_res_norm_drop    = 1.0e-10;
    moris::real tNLA_relaxation_parameter = 1.0;

    /* ------------------------------------------------------------------------ */
    // DUMMY FUNCTIONS
    /* ------------------------------------------------------------------------ */

    // Output criterion for VIS mesh
    bool Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return tOutputCriterion;
    }

    /* ------------------------------------------------------------------------ */
    // PARAMETER LISTS
    /* ------------------------------------------------------------------------ */

    void OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", false );
    }

    /* ------------------------------------------------------------------------ */

    void HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists( 0 ).set( "domain_dimensions", tDomainDims );
        aParameterLists( 0 ).set( "domain_offset", tDomainOffset );
        aParameterLists( 0 ).set( "domain_sidesets", tDomainSidesets );
        aParameterLists( 0 ).set( "lagrange_output_meshes", std::string( "0" ) );

        aParameterLists( 0 ).set( "lagrange_orders", tOrder );
        aParameterLists( 0 ).set( "lagrange_pattern", "0" );
        aParameterLists( 0 ).set( "bspline_orders", tOrder );
        aParameterLists( 0 ).set( "bspline_pattern", "0" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );

        aParameterLists( 0 ).set( "use_number_aura", 1 );

        aParameterLists( 0 ).set( "initial_refinement", "0" );

        aParameterLists( 0 ).set( "initial_refinement_pattern", "0" );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );
    }

    /* ------------------------------------------------------------------------ */

    void XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", std::string( "conformal" ) );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", std::string( "bspline" ) );
        aParameterLists( 0 ).set( "enrich_mesh_indices", std::string( "0" ) );
        aParameterLists( 0 ).set( "ghost_stab", true );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", false );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    uint get_phase_index_2d( const Bitset< 4 >& aGeometrySigns )
    {
        // left of left or right of right
        if ( !aGeometrySigns.test( 2 ) || !aGeometrySigns.test( 3 ) )
        {
            return 3;
        }

        // below bottom
        if ( !aGeometrySigns.test( 0 ) )
        {
            return 1;
        }

        // above top
        if ( !aGeometrySigns.test( 1 ) )
        {
            return 2;
        }

        // in solid
        return 0;
    }

    uint get_phase_index_3d( const Bitset< 6 >& aGeometrySigns )
    {
        // left of left or right of right
        if ( !aGeometrySigns.test( 2 ) || !aGeometrySigns.test( 3 ) || !aGeometrySigns.test( 4 ) || !aGeometrySigns.test( 5 ) )
        {
            return 3;
        }

        // below bottom
        if ( !aGeometrySigns.test( 0 ) )
        {
            return 1;
        }

        // above top
        if ( !aGeometrySigns.test( 1 ) )
        {
            return 2;
        }

        // in solid
        return 0;
    }

    /* ------------------------------------------------------------------------ */

    void GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );

        aParameterLists( 0 ).set( "number_of_phases", 4 );
        aParameterLists( 0 ).set( "phase_function_name", tDim == 3 ? "get_phase_index_3d" : "get_phase_index_2d" );

        gen::Field_Type tFtype = tDim == 3 ? gen::Field_Type::PLANE : gen::Field_Type::LINE;

        // Bottom
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( tFtype ) );
        aParameterLists( 1 ).set( "normal_x", 0.0 );
        aParameterLists( 1 ).set( "normal_y", 1.0 );

        // Top
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( tFtype ) );
        aParameterLists( 1 ).set( "normal_x", 0.0 );
        aParameterLists( 1 ).set( "center_y", 10.0 );
        aParameterLists( 1 ).set( "normal_y", -1.0 );

        // Left
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( tFtype ) );
        aParameterLists( 1 ).set( "normal_x", 1.0 );

        // Right
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( tFtype ) );
        aParameterLists( 1 ).set( "center_x", 1.0 );
        aParameterLists( 1 ).set( "normal_x", -1.0 );

        if ( tDim == 3 )
        {
            // front
            aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( tFtype ) );
            aParameterLists( 1 ).set( "normal_x", 0.0 );
            aParameterLists( 1 ).set( "normal_z", 1.0 );

            // back
            aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( tFtype ) );
            aParameterLists( 1 ).set( "normal_x", 0.0 );
            aParameterLists( 1 ).set( "center_z", 1.0 );
            aParameterLists( 1 ).set( "normal_z", -1.0 );
            }
    }

    /* ------------------------------------------------------------------------ */

    void FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // create a cell of cell of parameter list for fem

        uint tPhaseIndex   = 7;

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", std::string( "PhaseBulk" ) );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "0" );

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", std::string( "PhaseBottom" ) );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "1" );

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", std::string( "PhaseTop" ) );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "2" );

        aParameterLists( tPhaseIndex ).add_parameter_list( prm::create_phase_parameter_list() );
        aParameterLists( tPhaseIndex ).set( "phase_name", std::string( "PhaseSides" ) );
        aParameterLists( tPhaseIndex ).set( "phase_indices", "3" );

        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - STRUCTURE (ni-w-alloy?)
        //------------------------------------------------------------------------------

        // Youngs Modulus Shell
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropYoungsModulus" );
        aParameterLists( 0 ).set( "function_parameters", tYoungsModulus );

        // Poisson Ratio Shell
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", std::string( "PropPoissonRatio" ) );
        aParameterLists( 0 ).set( "function_parameters", tPoissonRatio );

        //------------------------------------------------------------------------------
        // BOUNDARY CONDITIONS
        //------------------------------------------------------------------------------

        // Compressive load
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", std::string( "PropTraction" ) );
        aParameterLists( 0 ).set( "function_parameters", tTraction );

        // Dirichlet structure
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", std::string( "PropDirichletStruct" ) );
        aParameterLists( 0 ).set( "function_parameters", tDirichlet );

        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // LINEAR ELASTICITY
        //------------------------------------------------------------------------------

        // linear elasticity - shell - 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso" );
        aParameterLists( 1 ).set( "model_type", static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );
        aParameterLists( 1 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( tDofString, "Displacement" ) );
        aParameterLists( 1 ).set( "properties", "PropYoungsModulus, YoungsModulus; PropPoissonRatio, PoissonRatio" );
        aParameterLists( 1 ).set( "phase_name", "PhaseBulk" );

        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMStrucNonLinIso" );
        aParameterLists( 1 ).set( "model_type", static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );
        aParameterLists( 1 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF ) );
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( tDofString, "Displacement" ) );
        aParameterLists( 1 ).set( "properties", "PropYoungsModulus, YoungsModulus; PropPoissonRatio, PoissonRatio" );
        aParameterLists( 1 ).set( "phase_name", "PhaseBulk" );

        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // NITSCHE DIRICHLET
        //------------------------------------------------------------------------------

        // Displacements - Shell - back wall
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheStruc" );
        aParameterLists( 2 ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        aParameterLists( 2 ).set( "function_parameters", "100.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropYoungsModulus,Material" );
        aParameterLists( 2 ).set( "leader_phase_name", "PhaseBulk" );

        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", std::string( "SPGhost" ) );
        aParameterLists( 2 ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::GHOST_DISPL ) );
        aParameterLists( 2 ).set( "function_parameters", std::string( "0.001" ) );
        aParameterLists( 2 ).set( "leader_properties", std::string( "PropYoungsModulus,Material" ) );
        aParameterLists( 2 ).set( "leader_phase_name", "PhaseBulk" );
        aParameterLists( 2 ).set( "follower_phase_name", "PhaseBulk" );

        ////////////////////////////////////////////////////////////////////////////////
        //------------------------------------------------------------------------------
        // BULK IWGs

        //------------------------------------------------------------------------------

        // linear elasticity
        ///*
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGStructBulk" );
        aParameterLists( 3 ).set( "IWG_bulk_type", (uint)fem::Element_Type::BULK );
        aParameterLists( 3 ).set( "leader_phase_name", "PhaseBulk" );
        aParameterLists( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        aParameterLists( 3 ).set( "dof_residual", tDofString );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofString );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        //*/

        //------------------------------------------------------------------------------

        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGStructBulkLBA" );
        aParameterLists( 3 ).set( "IWG_bulk_type", (uint)fem::Element_Type::BULK );
        aParameterLists( 3 ).set( "leader_phase_name", "PhaseBulk" );
        aParameterLists( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_NON_LINEAR_GEOMETRIC_STIFFNESS ) );
        aParameterLists( 3 ).set( "dof_residual", tDofString );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofString );
        if ( tStressType == "PK2" )
        {
            aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucNonLinIso,ElastLinIso" );
        }
        else
        {
            aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        }

        //------------------------------------------------------------------------------
        // NEUMANN BCs - IWGs
        //------------------------------------------------------------------------------

        // Traction BC ( Compressive force at top )
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGNeumannTraction" );
        aParameterLists( 3 ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        aParameterLists( 3 ).set( "leader_phase_name", "PhaseBulk" );
        aParameterLists( 3 ).set( "neighbor_phases", "PhaseTop" );
        aParameterLists( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
        aParameterLists( 3 ).set( "dof_residual", tDofString );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofString );
        aParameterLists( 3 ).set( "leader_properties", "PropTraction,Traction" );

        //------------------------------------------------------------------------------
        // DIRICHLET BCS - IWGs
        //------------------------------------------------------------------------------
        // Fixed bottom edge
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDirichletStruct" );
        aParameterLists( 3 ).set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
        aParameterLists( 3 ).set( "leader_phase_name", "PhaseBulk" );
        aParameterLists( 3 ).set( "neighbor_phases", "PhaseBottom" );
        aParameterLists( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) );
        aParameterLists( 3 ).set( "dof_residual", tDofString );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofString );
        aParameterLists( 3 ).set( "leader_properties", "PropDirichletStruct,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );

        //------------------------------------------------------------------------------
        // GHOST - IWGs
        //------------------------------------------------------------------------------
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", std::string( "IWGGhostMaterial" ) );
        aParameterLists( 3 ).set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
        aParameterLists( 3 ).set( "leader_phase_name", "PhaseBulk" );
        aParameterLists( 3 ).set( "follower_phase_name", "PhaseBulk" );
        aParameterLists( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::GHOST_NORMAL_FIELD ) );
        aParameterLists( 3 ).set( "dof_residual", tDofString );
        aParameterLists( 3 ).set( "leader_dof_dependencies", tDofString );
        aParameterLists( 3 ).set( "follower_dof_dependencies", tDofString );
        aParameterLists( 3 ).set( "stabilization_parameters", std::string( "SPGhost,GhostSP" ) );
        aParameterLists( 3 ).set( "ghost_order", (uint)std::stoi( tOrder ) );

        ////////////////////////////////////////////////////////////////////////////////
        // Volume IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQITotalVolume" );
        aParameterLists( 4 ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        aParameterLists( 4 ).set( "leader_dof_dependencies", tDofString );
        aParameterLists( 4 ).set( "leader_phase_name", "PhaseBulk" );

        // Eigen values and Eigenvectors

        for ( sint i = 0; i < tNumEigenVectors; ++i )
        {
            // eigen value
            aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
            aParameterLists( 4 ).set( "IQI_name", "IQIEIGENVALUE" + std::to_string( i + 1 ) );
            aParameterLists( 4 ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::EIGEN_VALUE ) );
            aParameterLists( 4 ).set( "leader_dof_dependencies", tDofString );
            aParameterLists( 4 ).set( "vectorial_field_index", i );
            aParameterLists( 4 ).set( "function_parameters", "0" );
            aParameterLists( 4 ).set( "leader_phase_name", "PhaseBulk" );

            // X-displacement
            aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
            aParameterLists( 4 ).set( "IQI_name", "IQIEIGENVEC" + std::to_string( i + 1 ) + "X" );
            aParameterLists( 4 ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::EIGEN_VECTOR ) );
            aParameterLists( 4 ).set( "dof_quantity", tDofString );
            aParameterLists( 4 ).set( "leader_dof_dependencies", tDofString );
            aParameterLists( 4 ).set( "vectorial_field_index", 0 );
            aParameterLists( 4 ).set( "function_parameters", std::to_string( i ) );
            aParameterLists( 4 ).set( "leader_phase_name", "PhaseBulk" );

            // Y-displacement
            aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
            aParameterLists( 4 ).set( "IQI_name", "IQIEIGENVEC" + std::to_string( i + 1 ) + "Y" );
            aParameterLists( 4 ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::EIGEN_VECTOR ) );
            aParameterLists( 4 ).set( "dof_quantity", tDofString );
            aParameterLists( 4 ).set( "leader_dof_dependencies", tDofString );
            aParameterLists( 4 ).set( "vectorial_field_index", 1 );
            aParameterLists( 4 ).set( "function_parameters", std::to_string( i ) );
            aParameterLists( 4 ).set( "leader_phase_name", "PhaseBulk" );

            if ( tDim == 3 )
            {
                // Z-displacement
                aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
                aParameterLists( 4 ).set( "IQI_name", "IQIEIGENVEC" + std::to_string( i + 1 ) + "Z" );
                aParameterLists( 4 ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::EIGEN_VECTOR ) );
                aParameterLists( 4 ).set( "dof_quantity", tDofString );
                aParameterLists( 4 ).set( "leader_dof_dependencies", tDofString );
                aParameterLists( 4 ).set( "vectorial_field_index", 2 );
                aParameterLists( 4 ).set( "function_parameters", std::to_string( i ) );
                aParameterLists( 4 ).set( "leader_phase_name", "PhaseBulk" );
                    }
        }

        // create computation parameter list
        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    /* ------------------------------------------------------------------------ */

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        if ( tUsePetsc )
        {
            create_petsc_parameter_list( aParameterLists );
        }
        else
        {
            create_trilinos_parameter_list( aParameterLists );
        }
    }

    //---------------------------------------------------------------------------------------------------------------------------------------

    void
    create_trilinos_parameter_list( Module_Parameter_Lists& aParameterLists )
    {
        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::EIGEN_SOLVER ) );
        aParameterLists( 0 ).set( "Eigen_Algorithm", "EIGALG_BLOCK_DAVIDSON" );
        aParameterLists( 0 ).set( "Verbosity", false );
        aParameterLists( 0 ).set( "Which", "SM" );
        aParameterLists( 0 ).set( "Block_Size", 5 );          // Block Size should be same as Number of Eigen values
        aParameterLists( 0 ).set( "Num_Eig_Vals", 5 );        // Number of Eigen values should be same as Block Size
        aParameterLists( 0 ).set( "Num_Blocks", 2 );          // Number of Blocks should satisfy : Num_Blocks*Block_Size < InitVec Length
        aParameterLists( 0 ).set( "MaxSubSpaceDims", 75 );    // Max Subspace Dimension = 3*Block_Size*Num_Eig_Vals
        aParameterLists( 0 ).set( "Initial_Guess", 0 );
        aParameterLists( 0 ).set( "MaxRestarts", 20 );
        aParameterLists( 0 ).set( "Convergence_Tolerance", 1e-01 );
        aParameterLists( 0 ).set( "Relative_Convergence_Tolerance", true );
        aParameterLists( 0 ).set( "preconditioners", "0" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "0" );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "1" );
        aParameterLists( 1 ).set( "RHS_Matrix_Type", "GeomStiffMat" );    // MassMat or IdentityMat

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", 2 );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", true );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 1 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", false );
        aParameterLists( 2 ).set( "NLA_is_eigen_problem", true );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0,1" );
        aParameterLists( 3 ).set( "NLA_DofTypes", tDofString );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", tDofString );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        // aParameterLists( 6 ).set( "SOL_save_operator_to_matlab", "LBAMat" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK ) );
        // Ifpack Preconditioner parameters
        aParameterLists( 7 ).set( "ifpack_prec_type", "Amesos" );
        aParameterLists( 7 ).set( "amesos: solver type", "Amesos_Pardiso" );    // Amesos_Umfpack or Amesos_Pardiso

        // Preconditioner parameters
        aParameterLists( 7 ).set( "overlap-level", 0 );
        aParameterLists( 7 ).set( "schwarz: combine mode", "add" );    // for Amesos_Umfpack and Amesos_Pardiso provide this parameter with "add" mode
    }

    void
    create_petsc_parameter_list( Module_Parameter_Lists& aParameterLists )
    {
        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::PETSC ) );
        aParameterLists( 0 ).set( "KSPType", "gmres" );
        aParameterLists( 0 ).set( "preconditioners", "0" );

        // find max eigen value
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::SLEPC_SOLVER ) );
        aParameterLists( 0 ).set( "Eigen_Algorithm", "krylovschur" );
        aParameterLists( 0 ).set( "Which", std::string( "LM" ) );
        aParameterLists( 0 ).set( "Num_Eig_Vals", 5 );
        aParameterLists( 0 ).set( "STType", "shift_invert" );
        aParameterLists( 0 ).set( "sub_linear_solver", "0" );    // 10 shift_invert
        aParameterLists( 0 ).set( "is_symmetric", false );       // 10 shift_invert
        aParameterLists( 0 ).set( "Update_Flag", true );         // 10 shift_invert
        aParameterLists( 0 ).set( "Verbosity", false );

        // precondioerr
        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::PETSC ) );
        aParameterLists( 7 ).set( "PCType", "mumps" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "0" );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "1" );
        aParameterLists( 1 ).set( "RHS_Matrix_Type", "GeomStiffMat" );    // MassMat or IdentityMat

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", 2 );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", true );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_Linear_solver", 1 );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", 1 );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", false );
        aParameterLists( 2 ).set( "NLA_is_eigen_problem", true );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0,1" );
        aParameterLists( 3 ).set( "NLA_DofTypes", tDofString );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", tDofString );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        aParameterLists( 6 ).set( "SOL_TPL_Type", sol::MapType::Petsc );
    }

    /* ------------------------------------------------------------------------ */

    void MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists( 0 ).set( "number_eigen_vectors", tNumEigenVectors );
    }

    /* ------------------------------------------------------------------------ */

    void VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        aParameterLists( 0 ).set( "Set_Names", tBulk );

        std::string Field_Names = "VOLUME";
        std::string Field_Type  = "GLOBAL";
        std::string IQI_Names   = "IQITotalVolume";

        for ( sint i = 0; i < tNumEigenVectors; ++i )
        {
            Field_Names += ",";
            Field_Type += ",";
            IQI_Names += ",";

            Field_Names += ",EVAL" + std::to_string( i + 1 ) + ",EVEC" + std::to_string( i + 1 ) + "X,EVEC" + std::to_string( i + 1 ) + "Y";
            Field_Type += ",GLOBAL,NODAL,NODAL";
            IQI_Names += ",IQIEIGENVALUE" + std::to_string( i + 1 ) + ",IQIEIGENVEC" + std::to_string( i + 1 ) + "X,IQIEIGENVEC" + std::to_string( i + 1 ) + "Y";
            if ( tDim == 3 )
            {
                Field_Names += ",EVEC" + std::to_string( i + 1 ) + "Z";
                Field_Type += ",NODAL";
                IQI_Names += ",IQIEIGENVEC" + std::to_string( i + 1 ) + "Z";
            }
        }

        aParameterLists( 0 ).set( "Field_Names", Field_Names );
        aParameterLists( 0 ).set( "Field_Type", Field_Type );
        aParameterLists( 0 ).set( "IQI_Names", IQI_Names );

        aParameterLists( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif

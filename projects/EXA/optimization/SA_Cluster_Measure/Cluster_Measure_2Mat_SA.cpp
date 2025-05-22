/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Cluster_Measure_2Mat_SA.cpp
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

#include "AztecOO.h"

//---------------------------------------------------------------

// global variable for interpolation order
extern uint gInterpolationOrder;

// problem dimension: 2D or 3D
extern uint gDim;

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    // ghost
    bool tUseGhost = true;

    // ghost penalty value
    std::string tGhostPenalty = "0.001";

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName = "Cluster_Measure_2Mat_SA" + std::to_string( gDim )
                                + ".exo";
    std::string tGENOutputFileName = "Cluster_Measure_2Mat_SA" + std::to_string( gDim )
                                   + "_GEN.exo";
    std::string tSoFileName   = "Cluster_Measure_2Mat_SA.so";
    std::string tHdf5FileName = "Cluster_Measure_2Mat_SA" + std::to_string( gDim )
                              + ".hdf5";

    /* ------------------------------------------------------------------------ */
    // Mesh Set Information
    std::string tPhase1      = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tPhase0      = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    std::string tTotalDomain = tPhase0 + "," + tPhase1;

    // FD
    std::string tFDEpsilon = "1e-6";
    std::string tWhichFD   = "all";

    /* ------------------------------------------------------------------------ */
    // boundary condition

    // prescribed temperature
    moris::real sTpre = 0.0;

    // flux at tip fac
    moris::real sP2 = 5.0;

    /* ------------------------------------------------------------------------ */
    // material parameters

    // capacity
    std::string sCap1 = "0.0";

    // density
    std::string tDens1 = "0.0";

    // conductivity
    moris::real sK1 = 1.0;

    // body flux
    moris::real sQ1 = 1.0;

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    Vector< uint > tNumElemsPerDim( gDim, 1 );
    Vector< real > tDomainDims( gDim, 1.0 );

    int tRefineBuffer = 1;

    /* ------------------------------------------------------------------------ */
    // Phase assignment

    std::string tGetPhaseIndex = gDim == 2 ? "get_phase_index_2d" : "get_phase_index_3d";

    uint get_phase_index_2d( const Bitset< 1 >& aGeometrySigns )
    {
        // by default in void
        uint tPhaseIndex = 0;

        // Phase solid 1
        if ( !aGeometrySigns.test( 0 ) )
        {
            tPhaseIndex = 1;
        }

        return tPhaseIndex;
    }

    uint get_phase_index_3d( const Bitset< 1 >& aGeometrySigns )
    {
        // by default in void
        uint tPhaseIndex = 0;

        // Phase solid 1
        if ( !aGeometrySigns.test( 0 ) )
        {
            tPhaseIndex = 1;
        }

        return tPhaseIndex;
    }

    /* ------------------------------------------------------------------------ */
    // Solver config

    moris::real tNLA_rel_res_norm_drop    = 1.0e-08;
    moris::real tNLA_relaxation_parameter = 1.0;
    int         tNLA_max_iter             = 2;

    int         tTSA_Num_Time_Steps = 1;
    moris::real tTSA_Time_Frame     = 1.0e0;

    /* ------------------------------------------------------------------------ */

    bool Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */

    uint             tNumConstraints = 8;
    Matrix< DDSMat > get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, tNumConstraints, 1 );

        return tConstraintTypes;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat > compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );
        tObjectives( 0 ) = aCriteria( 0 );

        return tObjectives;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat > compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, tNumConstraints );
        tConstraints( 0 ) = aCriteria( 1 );
        tConstraints( 1 ) = aCriteria( 2 );
        tConstraints( 2 ) = aCriteria( 3 );
        tConstraints( 3 ) = aCriteria( 4 );
        tConstraints( 4 ) = aCriteria( 5 );
        tConstraints( 5 ) = aCriteria( 6 );
        tConstraints( 6 ) = aCriteria( 7 );
        tConstraints( 7 ) = aCriteria( 8 );

        return tConstraints;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat > compute_dobjective_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );

        return tDObjectiveDADV;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat > compute_dobjective_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.size(), 0.0 );
        tDObjectiveDCriteria( 0 ) = 1;

        return tDObjectiveDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat > compute_dconstraint_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( tNumConstraints, aADVs.size(), 0.0 );
        return tDConstraintDADV;
    }

    /* ------------------------------------------------------------------------ */

    Matrix< DDRMat > compute_dconstraint_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( tNumConstraints, aCriteria.size(), 0.0 );
        tDConstraintDCriteria( 0, 1 ) = 1.0;
        tDConstraintDCriteria( 1, 2 ) = 1.0;
        tDConstraintDCriteria( 2, 3 ) = 1.0;
        tDConstraintDCriteria( 3, 4 ) = 1.0;
        tDConstraintDCriteria( 4, 5 ) = 1.0;
        tDConstraintDCriteria( 5, 6 ) = 1.0;
        tDConstraintDCriteria( 6, 7 ) = 1.0;
        tDConstraintDCriteria( 7, 8 ) = 1.0;

        return tDConstraintDCriteria;
    }

    /* ------------------------------------------------------------------------ */

    void OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", true );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tSoFileName );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::SWEEP );
        aParameterLists.set( "hdf5_path", tHdf5FileName );
        aParameterLists.set( "evaluate_objective_gradients", true );
        aParameterLists.set( "evaluate_constraint_gradients", true );
        aParameterLists.set( "num_evaluations_per_adv", "1" );
        aParameterLists.set( "include_bounds", false );
        aParameterLists.set( "finite_difference_type", tWhichFD );
        aParameterLists.set( "finite_difference_epsilons", tFDEpsilon );
        aParameterLists.set( "print", true );
    }

    void HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists.set( "domain_dimensions", tDomainDims );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", std::to_string( gInterpolationOrder ) );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "refinement_buffer", tRefineBuffer );
        aParameterLists.set( "staircase_buffer", tRefineBuffer );
    }

    void XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    void GENParameterList( Module_Parameter_Lists& aParameterLists )
    {

        // main gen parameter list
        aParameterLists.set( "number_of_phases", 2 );
        aParameterLists.set( "phase_function_name", tGetPhaseIndex );
        aParameterLists.set( "IQI_types",
                "IQIVolumeInterface1",
                "IQIVolumeMat0",
                "IQIVolumeMat1",
                "IQISidesetMeasure",
                "IQISidesetLengthMeasure",
                "IQIBulkMeasure",
                "IQIBulkLengthMeasure",
                "IQIBulkStrainEnergy0",
                "IQIBulkStrainEnergy1" );

        // interface plane
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        aParameterLists.set( "center_x", 1.35, 1.35, 1.35 );
        aParameterLists.set( "center_y", 0.0, 0.0, 0.0 );
        aParameterLists.set( "normal_x", 1.0, 1.0, 1.0 );
        aParameterLists.set( "normal_y", 1.0, 1.0, 1.0 );
    }

    void FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
         //------------------------------------------------------------------------------

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseMaterial0" );
        aParameterLists.set( "phase_indices", "0" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseMaterial1" );
        aParameterLists.set( "phase_indices", "1" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseAll" );
        aParameterLists.set( "phase_indices", "0,1" );

        //------------------------------------------------------------------------------

        // properties for material 0
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity0" );
        aParameterLists.set( "function_parameters", tDens1 );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacity0" );
        aParameterLists.set( "function_parameters", sCap1 );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity0" );
        aParameterLists.set( "function_parameters", std::to_string( sK1 ) );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropHeatLoad0" );
        aParameterLists.set( "function_parameters", std::to_string( sQ1 ) );

        // properties for material 1
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity1" );
        aParameterLists.set( "function_parameters", tDens1 );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacity1" );
        aParameterLists.set( "function_parameters", sCap1 );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity1" );
        aParameterLists.set( "function_parameters", std::to_string( sK1 ) );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropHeatLoad1" );
        aParameterLists.set( "function_parameters", std::to_string( sQ1 ) );

        // surface flux
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropSurfaceFlux" );
        aParameterLists.set( "function_parameters", std::to_string( sP2 ) );

        // temperature at back surface
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropImposedTemperature" );
        aParameterLists.set( "function_parameters", "2.0" );

        //------------------------------------------------------------------------------

        // create diffusion CM for material 0
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion0" );
        aParameterLists.set( "phase_name", "PhaseMaterial0" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity0 , Conductivity;"
                "PropDensity0      , Density;"
                "PropCapacity0     , HeatCapacity" );

        // create diffusion CM for material 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion1" );
        aParameterLists.set( "phase_name", "PhaseMaterial1" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity1 , Conductivity;"
                "PropDensity1      , Density;"
                "PropCapacity1     , HeatCapacity" );

        //------------------------------------------------------------------------------

        // Nitsche for Dirichlet on material 0
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheTemp0" );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial0" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity0,Material" );

        // Nitsche for Dirichlet on material 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheTemp1" );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity1,Material" );

        // Nitsche stabilization parameter for mat0 - mat1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPInterfaceNitsche01" );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial0" );
        aParameterLists.set( "follower_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity0,Material" );
        aParameterLists.set( "follower_properties", "PropConductivity1,Material" );

        // cluster measure for sideset
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPSidesetMeasure" );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::MEASURE );
        aParameterLists.set( "cluster_measures", std::pair< std::string, std::string >( "CELL_SIDE_MEASURE,PRIMARY,LEADER", "ElementSize" ) );

        // cluster measure for sideset length
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPSidesetLengthMeasure" );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::MEASURE );
        aParameterLists.set( "cluster_measures", std::pair< std::string, std::string >( "CELL_LENGTH_MEASURE,PRIMARY,LEADER", "ElementSize" ) );

        // cluster measure for bulk
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPBulkMeasure" );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::MEASURE );
        aParameterLists.set( "cluster_measures", std::pair< std::string, std::string >( "CELL_MEASURE,PRIMARY,LEADER", "ElementSize" ) );

        // cluster measure for bulk length
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPBulkLengthMeasure" );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::MEASURE );
        aParameterLists.set( "cluster_measures", std::pair< std::string, std::string >( "CELL_LENGTH_MEASURE,PRIMARY,LEADER", "ElementSize" ) );

        if ( tUseGhost )
        {
            // ghost penalty for material 0
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGPTemp0" );
            aParameterLists.set( "leader_phase_name", "PhaseMaterial0" );
            aParameterLists.set( "follower_phase_name", "PhaseMaterial0" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            aParameterLists.set( "function_parameters", tGhostPenalty );
            aParameterLists.set( "leader_properties", "PropConductivity0,Material" );
            aParameterLists.set( "follower_properties", "PropConductivity0,Material" );

            // ghost penalty for material 1
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGPTemp1" );
            aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
            aParameterLists.set( "follower_phase_name", "PhaseMaterial1" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            aParameterLists.set( "function_parameters", tGhostPenalty );
            aParameterLists.set( "leader_properties", "PropConductivity1,Material" );
            aParameterLists.set( "follower_properties", "PropConductivity1,Material" );
        }

        //------------------------------------------------------------------------------
        // create IWG for material 0 - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusion0Bulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial0" );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion0,Diffusion" );
        aParameterLists.set( "leader_properties", "PropHeatLoad0,Load" );

        // create IWG for material 1 - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusion1Bulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        aParameterLists.set( "leader_properties", "PropHeatLoad1,Load" );

        // create IWG for interface conditions mat0 - mat1
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInterface01" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial0" );
        aParameterLists.set( "follower_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion0,Diffusion" );
        aParameterLists.set( "follower_constitutive_models", "CMDiffusion1,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPInterfaceNitsche01,NitscheInterface" );

        // create IWG for Neumann boundary conditions
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInletFlux" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_NEUMANN );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "side_ordinals", "4" );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropSurfaceFlux,Neumann" );

        // create IWG for Dirichlet boundary conditions on mat 0
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGOutletTemp0" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial0" );
        aParameterLists.set( "side_ordinals", "4" );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropImposedTemperature,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion0,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp0,DirichletNitsche" );

        // create IWG for Dirichlet boundary conditions on mat 1
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGOutletTemp1" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "side_ordinals", "4" );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropImposedTemperature,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion1,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp1,DirichletNitsche" );

        if ( tUseGhost )
        {
            // create IWG for ghost on material 0
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGP0Temp" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "leader_phase_name", "PhaseMaterial0" );
            aParameterLists.set( "follower_phase_name", "PhaseMaterial0" );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTemp0,GhostSP" );
                // create IWG for ghost on material 1
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGP1Temp" );
            aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
            aParameterLists.set( "follower_phase_name", "PhaseMaterial1" );
            aParameterLists.set( "dof_residual", "TEMP" );
            aParameterLists.set( "stabilization_parameters", "SPGPTemp1,GhostSP" );
        }

        //------------------------------------------------------------------------------
        // Nodal Temperature IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "leader_phase_name", "PhaseAll" );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // Volume mat0
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIVolumeMat0" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial0" );

        // Volume mat1
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIVolumeMat1" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );

        // Interface length
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIVolumeInterface1" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "neighbor_phases", "PhaseMaterial0" );

        // strain energy mat0
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy0" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial0" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion0,Elast" );

        // strain energy mat1
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy1" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion1,Elast" );

        // IQI with sideset measure
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQISidesetMeasure" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STABILIZATION );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "neighbor_phases", "PhaseMaterial0" );
        aParameterLists.set( "stabilization_parameters", "SPSidesetMeasure,Stabilization" );

        // IQI with sideset length measure
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQISidesetLengthMeasure" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STABILIZATION );
        aParameterLists.set( "IQI_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "neighbor_phases", "PhaseMaterial0" );
        aParameterLists.set( "stabilization_parameters", "SPSidesetLengthMeasure,Stabilization" );

        // IQI with bulk measure
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkMeasure" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STABILIZATION );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "stabilization_parameters", "SPBulkMeasure,Stabilization" );

        // IQI with bulk length measure
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkLengthMeasure" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STABILIZATION );
        aParameterLists.set( "leader_phase_name", "PhaseMaterial1" );
        aParameterLists.set( "stabilization_parameters", "SPBulkLengthMeasure,Stabilization" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
        aParameterLists.set( "print_physics_model", false );
        // sensitivity analysis
        aParameterLists.set( "is_analytical_sensitivity", false );
        aParameterLists.set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        aParameterLists.set( "finite_difference_perturbation_size", stod( tFDEpsilon ) );
        // forward analysis
        aParameterLists.set( "is_analytical_forward", true );
        aParameterLists.set( "finite_difference_scheme_forward", fem::FDScheme_Type::POINT_3_CENTRAL );
        aParameterLists.set( "finite_difference_perturbation_size_forward", stod( tFDEpsilon ) );
        // perturbation strategy
        aParameterLists.set( "finite_difference_perturbation_strategy", fem::Perturbation_Type::ABSOLUTE );
    }

    void SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "TEMP" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        aParameterLists.set( "TSA_Time_Frame", tTSA_Time_Frame );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    void MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "order_adofs_by_host", false );
    }

    void VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tTotalDomain );
        aParameterLists.set( "Field_Names",
                "TEMP,VOLUME_MAT1,VOLUME_MAT0,VOLUME_INT1" );
        aParameterLists.set( "Field_Type", "NODAL,GLOBAL,GLOBAL,GLOBAL" );
        aParameterLists.set( "IQI_Names",
                "IQIBulkTEMP,IQIVolumeMat1,IQIVolumeMat0,IQIVolumeInterface1" );
        aParameterLists.set( "Save_Frequency", 1 );
    }

    void MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif

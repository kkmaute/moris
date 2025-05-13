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
 #include "cl_Bitset.hpp"
 #include "linalg_typedefs.hpp"
 
 #include "cl_FEM_Field_Interpolator_Manager.hpp"
 #include "cl_MSI_Equation_Object.hpp"
 #include "cl_TSA_Time_Solver.hpp"
 #include "cl_DLA_Solver_Interface.hpp"
 #include "cl_DLA_Linear_Solver_Aztec.hpp"
 #include "fn_PRM_FEM_Parameters.hpp"
 #include "fn_PRM_MSI_Parameters.hpp"
 #include "fn_PRM_SOL_Parameters.hpp"
 #include "fn_PRM_VIS_Parameters.hpp"
 #include "fn_PRM_HMR_Parameters.hpp"
 #include "fn_PRM_GEN_Parameters.hpp"
 #include "fn_PRM_XTK_Parameters.hpp"
 #include "fn_PRM_OPT_Parameters.hpp"
 #include "fn_equal_to.hpp"
 #include "fn_stringify_matrix.hpp"
 
 #include "AztecOO.h"
 #include "fn_dot.hpp"
 
 //------------------------------------------------------------------------------
 
 // interpolation order
 std::string tOrder = "1";
 
 // spatial dimensions
 uint tDim = 2;
 
 // Mesh refinement - first entry for FEM background mesh, second entry for level set background mesh
 std::string tInitialRefinement = "3,1";
 
 int  tNLA_max_iter  = 10;
 
 // file name
 std::string tName = "Interpolation_based_immersed_initialize";
 
 
 #ifdef __cplusplus
 extern "C" {
 #endif
 
 //------------------------------------------------------------------------------
 
 namespace moris
 {
     void create_trilinos_parameter_list( Vector< Vector< Parameter_List > >& );
     void create_petsc_parameter_list( Vector< Vector< Parameter_List > >& );
 
     /* ------------------------------------------------------------------------ */
     // General
     bool tIsOpt = true;
 
     // Prescribed displacements
     std::string tDirichlet = "0.0;0.0";
     
     // TSA arguments
     int tTSA_Num_Time_Steps = 3;
     moris::real tTSA_Time_Frame = 0.2;
 
 
     // Define library name for OPT and output file name.
     std::string tOutputFileName  = tName + ".exo";
     std::string tLibraryName     = tName + ".so";
    
     std::string tDofString = "UX,UY";
     int tNumMaxGcmmaIts = 2;
 
     // relative step size for MMA
     moris::real tStepSize = 0.1;
 
     
     real tInitialEnergy = 1.76839e-05;
     real tInitialVolume = 3.80329;
     real tVolumeLimit = 1.0;
     real tPenalty = 0.;
 
 
     // pi
     real tPi = std::acos( -1.0 );
 
     /* ------------------------------------------------------------------------ */
     // HMR parameters
     // Done
 
     std::string tNumElemsPerDim = tDim == 3 ? "110,10,10" : "8,8";
     std::string tDomainDims     = tDim == 3 ? "11,11,2" : "4,4";
     std::string tDomainOffset   = tDim == 3 ? "0.0,0.0,0.0" : "-2.0,-2.0";
     std::string tDomainSidesets = tDim == 3 ? "1,2,3,4,5,6" : "1,2,3,4";
     
     /* ------------------------------------------------------------------------ */
     // Mesh Set Information - Not needed since BCs on phase.
 
     // Bulk Phases
     std::string tSolidDomain        = "HMR_dummy_n_p1,HMR_dummy_c_p1";
 
     // Side sets
     std::string tDirichletSets      = "SideSet_1_n_p1"; 
     std::string tLoadSetsY          = "SideSet_3_n_p1"; 
     std::string tLoadSetsXR         = "SideSet_2_n_p1";
     std::string tLoadSetsXL         = "SideSet_4_n_p1";
 
     // Ghost 
     std::string tSolidGhost         = "ghost_p1";
     /* ------------------------------------------------------------------------ */
     // material parameters, kg is scaled with a factor 1e-6
     std::string tYoungsModulus = "130000";
     std::string tPoissonRatio  = "0.34";
 
     /* ------------------------------------------------------------------------ */
     // Output Config
     bool tOutputCriterion = true;
 
     /* ------------------------------------------------------------------------ */
     // DUMMY FUNCTIONS
     /* ------------------------------------------------------------------------ */
 
     // Output criterion for VIS mesh
     bool Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
     {
         return tOutputCriterion;
     }
 
     /* ------------------------------------------------------------------------ */
     // Shape Optimization FUNCTIONS - Required for MORIS to exit cleanly out of the workflow.
     /* ------------------------------------------------------------------------ */
    Matrix< DDSMat >
     get_constraint_types()
     {
         Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );
 
         return tConstraintTypes;
     }
 
     /* ------------------------------------------------------------------------ */
 
     Matrix< DDRMat >
     compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
     {
         Matrix< DDRMat > tObjectives( 1, 1 );
         
         tObjectives( 0 ) = aCriteria( 0 ) / tInitialEnergy;
         
         real constraint = aCriteria( 1 ) / tInitialVolume - tVolumeLimit;
         
         if (constraint > 0.0)        
         {
              tObjectives( 0 ) += tPenalty * std::pow( constraint, 2.0);
         }
         
         std::cout << "aCriteria( 0 ) : " << aCriteria( 0 )   << std::endl;
         std::cout << "Objective      : " << tObjectives( 0 ) << std::endl;
 
         return tObjectives;
     }
 
     /* ------------------------------------------------------------------------ */
 
 
 
     Matrix< DDRMat >
     compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
     {
         Matrix< DDRMat > tConstraints( 1, 1 );
         
         tConstraints( 0 ) = aCriteria( 1 ) / tInitialVolume - tVolumeLimit;
 
         std::cout << "aCriteria( 1 ) : " << aCriteria( 1 )    << std::endl;
         std::cout << "Constraint     : " << tConstraints( 0 ) << std::endl;
 
         return tConstraints;
     }
 
     /* ------------------------------------------------------------------------ */
 
     Matrix< DDRMat >
     compute_dobjective_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
     {
         Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );
 
         return tDObjectiveDADV;
     }
 
     /* ------------------------------------------------------------------------ */
 
     Matrix< DDRMat >
     compute_dobjective_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
     {
         Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.size(), 0.0 );
        
         tDObjectiveDCriteria( 0 ) = 1.0 / tInitialEnergy;
 
         real constraint = aCriteria( 1 ) / tInitialVolume - tVolumeLimit;
         
         if (constraint > 0.0)        
         {
              tDObjectiveDCriteria( 1 ) = 2.0 * tPenalty * constraint / tInitialVolume;
         }
 
         return tDObjectiveDCriteria;
     }
 
     /* ------------------------------------------------------------------------ */
 
     Matrix< DDRMat >
     compute_dconstraint_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
     {
         Matrix< DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );
 
         return tDConstraintDADV;
     }
 
     /* ------------------------------------------------------------------------ */
 
     Matrix< DDRMat >
     compute_dconstraint_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
     {
         Matrix< DDRMat > tDConstraintDCriteria( 1, aCriteria.size(), 0.0 );
         
         tDConstraintDCriteria( 0, 1 ) = 1.0 / tInitialVolume;
 
         return tDConstraintDCriteria;
     }
 
    
     /* ------------------------------------------------------------------------ */
     // NEUMANN BC FUNCTION
     /* ------------------------------------------------------------------------ */
 
     /* ------------------------------------------------------------------------ */
     // PARAMETER LISTS
     /* ------------------------------------------------------------------------ */
 
     void OPTParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
     {
         tParameterlist.resize( 3 );
         tParameterlist( 0 ).resize( 1 );
         tParameterlist( 1 ).resize( 0 );
         tParameterlist( 2 ).resize( 1 );
 
         // General Config.
 
         tParameterlist( 0 )( 0 ) = moris::prm::create_opt_problem_parameter_list();
         tParameterlist( 0 )( 0 ).set( "is_optimization_problem", true );
         tParameterlist( 0 )( 0 ).set( "problem", "user_defined" );
         tParameterlist( 0 )( 0 ).set( "library", tLibraryName );
         // Uncomment below line for the design iteration input file.
         //tParameterlist( 0 )( 0 ).set( "restart_file", "ADV_Alg_0_IterDisc2.hdf5" );
          

         tParameterlist( 2 )( 0 ) = moris::prm::create_sweep_parameter_list();
         tParameterlist( 2 )( 0 ).set( "hdf5_path", "shape_opt_test.hdf5" );
         tParameterlist( 2 )( 0 ).set( "num_evaluations_per_adv", "1" );
         tParameterlist( 2 )( 0 ).set( "finite_difference_type", "all" );
         tParameterlist( 2 )( 0 ).set( "finite_difference_epsilons", "1e-6" );


     
     }
 
     /* ------------------------------------------------------------------------ */
 
     void HMRParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
     {
         tParameterlist.resize( 1 );
         tParameterlist( 0 ).resize( 1 );
 
         tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();
 
         tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
         tParameterlist( 0 )( 0 ).set( "domain_dimensions", tDomainDims );
         tParameterlist( 0 )( 0 ).set( "domain_offset", tDomainOffset );
         tParameterlist( 0 )( 0 ).set( "domain_sidesets", tDomainSidesets );
         tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", std::string( "0" ) );
 
         tParameterlist( 0 )( 0 ).set( "lagrange_orders", tOrder );
         tParameterlist( 0 )( 0 ).set( "lagrange_pattern", "0" );
         tParameterlist( 0 )( 0 ).set( "bspline_orders", "1" );
         tParameterlist( 0 )( 0 ).set( "bspline_pattern", "0" );
         tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );
 
         tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );
 
         tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );
 
         tParameterlist( 0 )( 0 ).set( "initial_refinement", "1" );
 
         tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );
 
         tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
         tParameterlist( 0 )( 0 ).set( "severity_level", 0 );
         tParameterlist( 0 )( 0 ).set( "use_advanced_T_matrix_scheme", 0 );
 
     }
 
     /* ------------------------------------------------------------------------ */
 
     void XTKParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
     {
         tParameterlist.resize( 1 );
         tParameterlist( 0 ).resize( 1 );
 
         tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
         tParameterlist( 0 )( 0 ).set( "decompose", true );
         tParameterlist( 0 )( 0 ).set( "decomposition_type", std::string( "conformal" ) );
         tParameterlist( 0 )( 0 ).set( "enrich", true );
         tParameterlist( 0 )( 0 ).set( "basis_rank", std::string( "bspline" ) );
         tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", std::string( "0" ) );
         tParameterlist( 0 )( 0 ).set( "ghost_stab", true );
         tParameterlist( 0 )( 0 ).set( "multigrid", false );
         tParameterlist( 0 )( 0 ).set( "verbose", true );
         tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", false );
         tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
         tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
         tParameterlist( 0 )( 0 ).set( "write_cluster_measures_to_exo", false );
         tParameterlist( 0 )( 0 ).set( "triangulate_all", true );
         tParameterlist( 0 )( 0 ).set( "use_SPG_based_enrichment", true );
         tParameterlist( 0 )( 0 ).set( "global_T_matrix_output_file", std::string("Global_Extraction_Operator") );
         tParameterlist( 0 )( 0 ).set( "activate_basis_agglomeration", true );
 
 
     }
 
     /* ------------------------------------------------------------------------ */
 
     
 
     /* ------------------------------------------------------------------------ */

 
     /* ------------------------------------------------------------------------ */
 
     void GENParameterList( Vector< Vector< Parameter_List > >& tParameterlist )    // Done
     {
         
         tParameterlist.resize( 3 );
         tParameterlist( 0 ).resize( 1 );
 
 
         // Main GEN parameter list
         tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();
 
         // Phase table assignment - please change as required.
         Matrix< DDUMat > tPhaseMap( 32, 1, 1 );
         tPhaseMap( 21 )             = 1;    
         tPhaseMap( 23 )             = 1;    
         tPhaseMap( 22 )             = 1;    
         tPhaseMap( 31 )             = 0;   
         tPhaseMap( 30 )             = 2;    
         tPhaseMap( 26 )             = 3;    
         tPhaseMap( 27 )             = 3;    
         tPhaseMap( 25 )             = 3; 
         tPhaseMap( 29 )             = 4;    
         tPhaseMap( 15 )             = 5;
         std::string tPhaseMapString = moris::ios::stringify( tPhaseMap );

         // MORIS will print the phase table as part of the output log with this argument set to true
         tParameterlist( 0 )( 0 ).set( "print_phase_table", false );

         // Assign the created phase table.
         tParameterlist( 0 )( 0 ).set( "phase_table", tPhaseMapString );
         
         //Argument to write design extraction operator to file
         tParameterlist( 0 )( 0 ).set( "write_design_extraction_operator_to_file", true );
 
 
         // init geometry counter
         uint tGeoCounter = 0;
 
         // This will create a field array of super ellipses.
         tParameterlist( 1 ).push_back( prm::create_field_array_parameter_list( gen::Field_Type::SUPERELLIPSE ) );
         tParameterlist( 1 )( tGeoCounter ).set( "semidiameter_x", 0.1875 );
         tParameterlist( 1 )( tGeoCounter ).set( "semidiameter_y", 0.1875 );
         tParameterlist( 1 )( tGeoCounter ).set( "exponent", 2.0 );
         tParameterlist( 1 )( tGeoCounter ).set( "lower_bound_x",  -1.5 + (0.5/6.0) + 0.25 );           // Left-most hole center -1.5 + (0.5/6.0) + 0.25 -1
         tParameterlist( 1 )( tGeoCounter ).set( "upper_bound_x",   1.5 - (0.5/6.0) - 0.25 );           // Right-most hole center 1.5 - (0.5/6.0) - 0.25 1
         tParameterlist( 1 )( tGeoCounter ).set( "lower_bound_y",  -1.0 + 0.375 );           // Bottom-most hole center -1.0 + 0.375 -0.5625
         tParameterlist( 1 )( tGeoCounter ).set( "upper_bound_y",   1.0 - 0.375 );           // Top-most hole center 1.0 - 0.375 0.5625
 
         tParameterlist( 1 )( tGeoCounter ).set( "number_of_fields_x", 5 );    // Number of holes in the x direction
         tParameterlist( 1 )( tGeoCounter ).set( "number_of_fields_y", 3);    // Number of holes in the y direction
 
         tParameterlist( 1 )( tGeoCounter ).set( "discretization_mesh_index", 1 );           // Index of B-spline mesh to create level set field on (-1 = none)
         tParameterlist( 1 )( tGeoCounter ).set( "discretization_lower_bound", -10.0 );    // Lower bound of level set field (if bspline_mesh_index >= 0)
         tParameterlist( 1 )( tGeoCounter ).set( "discretization_upper_bound",  10.0 );     // Upper bound of level set field (if bspline_mesh_index >= 0)
         tGeoCounter++;
         
         // Uncomment below if you wish to have a single hole instead of a field array of holes.
         /*tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( gen::Field_Type::SUPERELLIPSE ) );
         tParameterlist( 1 )( tGeoCounter ).set( "center_x", 0.0 );
         tParameterlist( 1 )( tGeoCounter ).set( "center_y", 0.0 );
         tParameterlist( 1 )( tGeoCounter ).set( "semidiameter_x", 0.25);
         tParameterlist( 1 )( tGeoCounter ).set( "semidiameter_y", 0.25);
         tParameterlist( 1 )( tGeoCounter ).set( "exponent", 2.0 );
         tParameterlist( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", true );
         tParameterlist( 1 )( tGeoCounter ).set( "discretization_mesh_index",     0);
         tParameterlist( 1 )( tGeoCounter ).set( "discretization_lower_bound", -20.0  );
         tParameterlist( 1 )( tGeoCounter ).set( "discretization_upper_bound",  20.0  );
         tGeoCounter++;*/
 
         gen::Field_Type tFtype = tDim == 3 ? gen::Field_Type::PLANE : gen::Field_Type::LINE;
 
         // Bottom
         tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( tFtype ) );
         tParameterlist( 1 )( tGeoCounter ).set( "center_y", -1.0 );
         tParameterlist( 1 )( tGeoCounter ).set( "normal_x", 0.0 );
         tParameterlist( 1 )( tGeoCounter ).set( "normal_y", 1.0 );
         tGeoCounter++;
 
         // Top
         tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( tFtype ) );
         tParameterlist( 1 )( tGeoCounter ).set( "normal_x", 0.0 );
         tParameterlist( 1 )( tGeoCounter ).set( "center_y", 1.0 );
         tParameterlist( 1 )( tGeoCounter ).set( "normal_y", -1.0 );
         tGeoCounter++;
 
         // Left
         tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( tFtype ) );
         tParameterlist( 1 )( tGeoCounter ).set( "normal_x", 1.0 );
         tParameterlist( 1 )( tGeoCounter ).set( "center_x", -1.5 );
         tGeoCounter++;
 
         // Right
         tParameterlist( 1 ).push_back( prm::create_level_set_geometry_parameter_list( tFtype ) );
         tParameterlist( 1 )( tGeoCounter ).set( "center_x", 1.5 );
         tParameterlist( 1 )( tGeoCounter ).set( "normal_x", -1.0 );
         tGeoCounter++;
         
     
     
     }
 
     /* ------------------------------------------------------------------------ */ 
 
     void
     MORISGENERALParameterList( Vector< Vector< Parameter_List > >& tParameterlist )
     {
     }
 
     /* ------------------------------------------------------------------------ */
 }    // namespace moris
 
 //------------------------------------------------------------------------------
 #ifdef __cplusplus
 }
 #endif
 

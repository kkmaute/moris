/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Library_IO_Meshgen.cpp
 *
 */

#include "cl_Library_IO_Meshgen.hpp"

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Meshgen::Library_IO_Meshgen()
            : Library_IO()    // initialize base class data as usual
    {
        // set the type of this library
        mLibraryType = Library_Type::MESHGEN;

        // list of supported parameter list types
        mSupportedParamListTypes = { Parameter_List_Type::HMR, Parameter_List_Type::XTK, Parameter_List_Type::GEN };
    }

    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Meshgen::~Library_IO_Meshgen()
    {
        // do nothing extra
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::finalize()
    {
        // FIXME: uncomment once complete
        // check that an .xml input file has been specified
        // MORIS_ERROR( mXmlParserIsInitialized,
        //         "Library_IO_Meshgen::finalize() - No .xml input file has been specified. "
        //         "This is required for the mesh generation workflow." );

        // load the standard parameters into the member variables
        this->load_all_standard_parameters();

        // if an .so file has been parsed, first use its parameters (if any were defined in it) to overwrite or add to the standard parameters
        if( mSoLibIsInitialized )
        {
            this->load_parameters_from_shared_object_library();
        }

        // load parameters from xml, overwrites parameters specified in either the standard parameters or an .so file if parsed
        if( mXmlParserIsInitialized )
        {
            this->load_parameters_from_xml();
        }

        // check the parameters for validity
        // TODO: this->check_parameters();

        // mark this library as finalized and lock it from modification
        mLibraryIsFinalized = true;

        // print receipt of the finalized library
        this->print_parameter_receipt( "./Parameter_Receipt.xml" ); // TODO: the file name and location should be user defineable
    }

    //------------------------------------------------------------------------------------------------------------------
    // STANDARD PARAMETER LIST FUNCTIONS
    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::load_all_standard_parameters()
    {
        // go over all modules and check that
        for( uint iModule = 0; iModule < (uint)( Parameter_List_Type::END_ENUM ); iModule++ )
        {
            // get the current module
            Parameter_List_Type tParamListType = (Parameter_List_Type)( iModule );

            // fill the parameterlist entry with the standard parameters
            this->create_standard_parameter_list_for_module( tParamListType, mParameterLists( iModule ) );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_parameter_list_for_module( 
            Parameter_List_Type aParamListType,
            ModuleParameterList& aParameterList )
    {
        switch( aParamListType )
        {
            case Parameter_List_Type::OPT : 
                this->create_standard_OPT_parameter_list( aParameterList );
                break;

            case Parameter_List_Type::HMR :
                this->create_standard_HMR_parameter_list( aParameterList );
                break;

            case Parameter_List_Type::XTK :
                this->create_standard_XTK_parameter_list( aParameterList );
                break;

            case Parameter_List_Type::GEN :
                this->create_standard_GEN_parameter_list( aParameterList );
                break;

            case Parameter_List_Type::END_ENUM :
                MORIS_ERROR( false, "Library_IO_Meshgen::create_standard_parameter_list_for_module() - "
                        "No standard library defined for module END_ENUM" );
                break;

            // create an empty parameter list for modules that are not needed
            default : 
                aParameterList = ModuleParameterList();
                break;
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_OPT_parameter_list( ModuleParameterList& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList.resize( 1 );
        aParameterList( 0 ).resize( 1 );
        aParameterList( 0 )( 0 ) = ParameterList();    //  prm::create_opt_problem_parameter_list();

        //----------------------------------------------------------------------------
        // FIXME: use standard parameters from PRM
        {
            aParameterList( 0 )( 0 ).insert( "is_optimization_problem", false );    
            aParameterList( 0 )( 0 ).insert( "workflow", "HMR_XTK" );               
            aParameterList( 0 )( 0 ).insert( "problem", "user_defined" );           
            aParameterList( 0 )( 0 ).insert( "restart_file", "" );                  
            aParameterList( 0 )( 0 ).insert( "finite_difference_type", "none" );    
            aParameterList( 0 )( 0 ).insert( "finite_difference_epsilons", "1E-8" );
            aParameterList( 0 )( 0 ).insert( "library", "" );                       
            aParameterList( 0 )( 0 ).insert( "reinitialize_interface_iter", INT_MAX );       
            aParameterList( 0 )( 0 ).insert( "first_reinitialize_interface_iter", INT_MAX ); 
        }

        //  nothing else
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_XTK_parameter_list( ModuleParameterList& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList.resize( 1 );
        aParameterList( 0 ).resize( 1 );
        aParameterList( 0 )( 0 ) = ParameterList();    //  prm::create_xtk_parameter_list();

        //----------------------------------------------------------------------------
        // FIXME: use standard parameters from PRM
        {
            aParameterList( 0 )( 0 ).insert( "decompose", true );
            aParameterList( 0 )( 0 ).insert( "decomposition_type", "conformal" );
            aParameterList( 0 )( 0 ).insert( "octree_refinement_level", "-1" );
            aParameterList( 0 )( 0 ).insert( "triangulate_all", false );
            aParameterList( 0 )( 0 ).insert( "ig_element_order", moris::uint( 1 ) );
            aParameterList( 0 )( 0 ).insert( "cleanup_cut_mesh", false );
            aParameterList( 0 )( 0 ).insert( "enrich", false );
            aParameterList( 0 )( 0 ).insert( "use_SPG_based_enrichment", false );
            aParameterList( 0 )( 0 ).insert( "basis_rank", "node" );
            aParameterList( 0 )( 0 ).insert( "enrich_mesh_indices", "0" );
            aParameterList( 0 )( 0 ).insert( "sort_basis_enrichment_levels", false );
            aParameterList( 0 )( 0 ).insert( "unenriched_mesh_indices", "" );
            aParameterList( 0 )( 0 ).insert( "ghost_stab", false );
            aParameterList( 0 )( 0 ).insert( "visualize_ghost", false );
            aParameterList( 0 )( 0 ).insert( "multigrid", false );
            aParameterList( 0 )( 0 ).insert( "contact_sandbox", false );
            aParameterList( 0 )( 0 ).insert( "potential_phases_in_contact", "" );
            aParameterList( 0 )( 0 ).insert( "bb_epsilon", 0.1 );
            aParameterList( 0 )( 0 ).insert( "verbose", false );
            aParameterList( 0 )( 0 ).insert( "verbose_level", moris::uint( 0 ) );
            aParameterList( 0 )( 0 ).insert( "diagnostics", false );
            aParameterList( 0 )( 0 ).insert( "diagnostics_id", "" );
            aParameterList( 0 )( 0 ).insert( "diagnostics_path", "" );
            aParameterList( 0 )( 0 ).insert( "deactivate_empty_sets", false );
            aParameterList( 0 )( 0 ).insert( "deactivate_all_but_blocks", "" );
            aParameterList( 0 )( 0 ).insert( "deactivate_all_but_side_sets", "" );
            aParameterList( 0 )( 0 ).insert( "write_enrichment_fields", false );
            aParameterList( 0 )( 0 ).insert( "write_enrichment_fields_probe_spheres", "" );
            aParameterList( 0 )( 0 ).insert( "global_T_matrix_output_file", "" );
            aParameterList( 0 )( 0 ).insert( "nodal_T_matrix_output_file", "" );
            aParameterList( 0 )( 0 ).insert( "MPC_output_file", "" );
            aParameterList( 0 )( 0 ).insert( "triangulate_all_in_post", false );
            aParameterList( 0 )( 0 ).insert( "output_path", "./" );
            aParameterList( 0 )( 0 ).insert( "output_file", "xtk_temp.exo" );
            aParameterList( 0 )( 0 ).insert( "print_enriched_ig_mesh", false );
            aParameterList( 0 )( 0 ).insert( "keep_all_opt_iters", false );
            aParameterList( 0 )( 0 ).insert( "exodus_output_XTK_ig_mesh", false );
            aParameterList( 0 )( 0 ).insert( "output_cut_ig_mesh", false );
            aParameterList( 0 )( 0 ).insert( "output_intersection_mesh", false );
            aParameterList( 0 )( 0 ).insert( "high_to_low_dbl_side_sets", false );
            aParameterList( 0 )( 0 ).insert( "print_memory", false );
            aParameterList( 0 )( 0 ).insert( "low_memory", true );
            aParameterList( 0 )( 0 ).insert( "probe_bg_cells", "" );
            aParameterList( 0 )( 0 ).insert( "union_blocks", "" );
            aParameterList( 0 )( 0 ).insert( "union_block_names", "" );
            aParameterList( 0 )( 0 ).insert( "union_block_colors", "" );
            aParameterList( 0 )( 0 ).insert( "union_side_sets", "" );
            aParameterList( 0 )( 0 ).insert( "union_side_set_names", "" );
            aParameterList( 0 )( 0 ).insert( "union_side_set_colors", "" );
            aParameterList( 0 )( 0 ).insert( "identify_hanging_nodes", false );
            aParameterList( 0 )( 0 ).insert( "delete_xtk_after_generation", true );
        }
        // FIXME: use standard parameters from PRM
        //----------------------------------------------------------------------------

        // enrichment
        aParameterList( 0 )( 0 ).set( "enrich", true );
        aParameterList( 0 )( 0 ).set( "basis_rank", "bspline" );

        // output
        aParameterList( 0 )( 0 ).set( "output_file", "foreground_mesh.exo" );
        aParameterList( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterList( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
        aParameterList( 0 )( 0 ).set( "print_enriched_ig_mesh", true );
        aParameterList( 0 )( 0 ).set( "global_T_matrix_output_file", "Global_Extraction_Operator" );
        aParameterList( 0 )( 0 ).set( "nodal_T_matrix_output_file", "Nodal_Extraction_Operators" );

        // TODO: ...
        // aParameterList( 0 )( 0 ).set( "enrich_mesh_indices",       "0" ) ; // TODO
        // aParameterList( 0 )( 0 ).set( "triangulate_all",           false ); // TODO
        // aParameterList( 0 )( 0 ).set( "ig_element_order",          tIgElementOrder ); // TODO
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_HMR_parameter_list( ModuleParameterList& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList.resize( 1 );
        aParameterList( 0 ).resize( 1 );
        aParameterList( 0 )( 0 ) = ParameterList();    //  prm::create_xtk_parameter_list();

        { // FIXME: use standard parameters from PRM
            aParameterList( 0 )( 0 ).insert( "number_of_elements_per_dimension", "2, 2" );
            aParameterList( 0 )( 0 ).insert( "processor_decomposition_method", 1 );
            aParameterList( 0 )( 0 ).insert( "processor_dimensions", "2, 2" );
            aParameterList( 0 )( 0 ).insert( "domain_dimensions", "1, 1" );
            aParameterList( 0 )( 0 ).insert( "domain_offset", "0, 0 " );
            aParameterList( 0 )( 0 ).insert( "domain_sidesets", "" );
            aParameterList( 0 )( 0 ).insert( "lagrange_output_meshes", "" );
            aParameterList( 0 )( 0 ).insert( "lagrange_output_meshe_names", "" );
            aParameterList( 0 )( 0 ).insert( "lagrange_input_meshes", "" );
            aParameterList( 0 )( 0 ).insert( "refinement_buffer", 0 );
            aParameterList( 0 )( 0 ).insert( "staircase_buffer", 0 );
            aParameterList( 0 )( 0 ).insert( "lagrange_orders", "1" );
            aParameterList( 0 )( 0 ).insert( "lagrange_pattern", "0" );
            aParameterList( 0 )( 0 ).insert( "bspline_orders", "1" );
            aParameterList( 0 )( 0 ).insert( "bspline_pattern", "0" );
            aParameterList( 0 )( 0 ).insert( "union_pattern", 6 );
            aParameterList( 0 )( 0 ).insert( "working_pattern", 7 );
            aParameterList( 0 )( 0 ).insert( "lagrange_to_bspline", "0" );
            aParameterList( 0 )( 0 ).insert( "severity_level", 0 );
            aParameterList( 0 )( 0 ).insert( "truncate_bsplines", 1 );
            aParameterList( 0 )( 0 ).insert( "use_multigrid", 0 );
            aParameterList( 0 )( 0 ).insert( "use_number_aura", 1 );
            aParameterList( 0 )( 0 ).insert( "initial_refinement", "0" );
            aParameterList( 0 )( 0 ).insert( "initial_refinement_pattern", "0" );
            aParameterList( 0 )( 0 ).insert( "write_background_mesh", "" );
            aParameterList( 0 )( 0 ).insert( "write_lagrange_output_mesh", "" );
            aParameterList( 0 )( 0 ).insert( "write_lagrange_output_mesh_to_exodus", "" );
            aParameterList( 0 )( 0 ).insert( "write_refinement_pattern_file", false );
            aParameterList( 0 )( 0 ).insert( "restart_refinement_pattern_file", "" );
            aParameterList( 0 )( 0 ).insert( "basis_function_vtk_file", "" );
            aParameterList( 0 )( 0 ).insert( "max_refinement_level", -1 );
            aParameterList( 0 )( 0 ).insert( "additional_lagrange_refinement", 0 );
            aParameterList( 0 )( 0 ).insert( "adaptive_refinement_level", 0 );
            aParameterList( 0 )( 0 ).insert( "use_refinement_interrelation", 0 );
            aParameterList( 0 )( 0 ).insert( "renumber_lagrange_nodes", 0 );
            aParameterList( 0 )( 0 ).insert( "use_advanced_T_matrix_scheme", 0 );
            aParameterList( 0 )( 0 ).insert( "refinement_function_names", "" );
            aParameterList( 0 )( 0 ).insert( "use_refine_low_level_elements", false );
        }

        // nothing else
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_GEN_parameter_list( ModuleParameterList& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList.resize( 3 );
        aParameterList( 0 ).resize( 1 );
        aParameterList( 0 )( 0 ) = ParameterList();

        { // FIXME: use standard parameters from PRM
            aParameterList( 0 )( 0 ).insert( "intersection_mode", "LEVEL_SET" );
            aParameterList( 0 )( 0 ).insert( "isocontour_threshold", 0.0 );
            aParameterList( 0 )( 0 ).insert( "isocontour_tolerance", 1e-12 );
            aParameterList( 0 )( 0 ).insert( "intersection_tolerance", 1e-12 );
            aParameterList( 0 )( 0 ).insert( "evaluate_new_pts_as_linear", false );
            aParameterList( 0 )( 0 ).insert( "output_mesh_file", "" );
            aParameterList( 0 )( 0 ).insert( "geometry_field_file", "" );
            aParameterList( 0 )( 0 ).insert( "time_offset", 0.0 );
            aParameterList( 0 )( 0 ).insert( "initial_advs", "" );
            aParameterList( 0 )( 0 ).insert( "advs_size", 0 );
            aParameterList( 0 )( 0 ).insert( "initial_advs_fill", 0.0 );
            aParameterList( 0 )( 0 ).insert( "lower_bounds", "" );
            aParameterList( 0 )( 0 ).insert( "lower_bounds_fill", 0.0 );
            aParameterList( 0 )( 0 ).insert( "upper_bounds", "" );
            aParameterList( 0 )( 0 ).insert( "upper_bounds_fill", 0.0 );
            aParameterList( 0 )( 0 ).insert( "IQI_types", "" );
            aParameterList( 0 )( 0 ).insert( "PDV_types", "" );
            aParameterList( 0 )( 0 ).insert( "phase_table", "" );
            aParameterList( 0 )( 0 ).insert( "phase_function_name", "" );
            aParameterList( 0 )( 0 ).insert( "number_of_phases", 0 );
            aParameterList( 0 )( 0 ).insert( "print_phase_table", false );
            aParameterList( 0 )( 0 ).insert( "diagnostics", false );
            aParameterList( 0 )( 0 ).insert( "diagnostics_id", "" );
            aParameterList( 0 )( 0 ).insert( "diagnostics_path", "" );
        }

        // nothing else
    }

    //------------------------------------------------------------------------------------------------------------------

    // FIXME: use standard parameters from PRM
    void
    Library_IO_Meshgen::create_standard_geometry_parameter_list( ParameterList & aParameterList )
    {
        aParameterList = ParameterList();
        aParameterList.insert( "type", "" );                          
        aParameterList.insert( "name", "" );                          
        aParameterList.insert( "field_variable_indices", "" );        
        aParameterList.insert( "adv_indices", "" );                   
        aParameterList.insert( "constant_parameters", "" );           
        aParameterList.insert( "number_of_refinements", "" );         
        aParameterList.insert( "refinement_mesh_index", "" );         
        aParameterList.insert( "refinement_function_index", -1 );     
        aParameterList.insert( "discretization_mesh_index", -2 );     
        aParameterList.insert( "discretization_lower_bound", -1.0 );  
        aParameterList.insert( "discretization_upper_bound", 1.0 );   
        aParameterList.insert( "multilinear_intersections", false );
    }

    //------------------------------------------------------------------------------------------------------------------

    // FIXME: use standard parameters from PRM
    void
    Library_IO_Meshgen::create_standard_user_defined_geometry_parameter_list( ParameterList & aParameterList )
    {
        this->create_standard_geometry_parameter_list( aParameterList );
        aParameterList.set( "type", "user_defined" );            
        aParameterList.insert( "field_function_name", "" );      
        aParameterList.insert( "sensitivity_function_name", "" );
    }

    //------------------------------------------------------------------------------------------------------------------

    // FIXME: use standard parameters from PRM
    void
    Library_IO_Meshgen::create_standard_voxel_field_parameter_list( ParameterList & aParameterList )
    {
        this->create_standard_geometry_parameter_list( aParameterList );
        aParameterList.set( "type", "voxel" );            
        aParameterList.insert( "voxel_field_file", "" );  
        aParameterList.insert( "domain_dimensions", "" ); 
        aParameterList.insert( "domain_offset", "" );     
        aParameterList.insert( "grain_id_value_map", "" );
    }

    //------------------------------------------------------------------------------------------------------------------

    // FIXME: use standard parameters from PRM
    void
    Library_IO_Meshgen::create_standard_image_sdf_field_parameter_list( ParameterList & aParameterList )
    {
        this->create_standard_geometry_parameter_list( aParameterList );
        aParameterList.set( "type", "image_sdf" );              
        aParameterList.insert( "image_file", "" );              
        aParameterList.insert( "image_dimensions", "" );        
        aParameterList.insert( "image_offset", "" );            
        aParameterList.insert( "image_sdf_scaling", 0.0 );      
        aParameterList.insert( "image_sdf_shift", 0.0 );        
        aParameterList.insert( "image_sdf_default", -1.0 );     
        aParameterList.insert( "image_sdf_interpolate", false );
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris

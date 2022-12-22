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
#include "enums.hpp"
#include "parameters.hpp"

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
        this->print_parameter_receipt( "./Parameter_Receipt.xml" ); // TODO: the file name and location should be user define-able
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

            // fill the parameter list entry with the standard parameters
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
        aParameterList( 0 )( 0 ) = prm::create_opt_problem_parameter_list(); // ParameterList();
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_XTK_parameter_list( ModuleParameterList& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList.resize( 1 );
        aParameterList( 0 ).resize( 1 );
        aParameterList( 0 )( 0 ) = prm::create_xtk_parameter_list(); // ParameterList();

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
        aParameterList( 0 )( 0 ) = prm::create_hmr_parameter_list(); // ParameterList(); 
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Meshgen::create_standard_GEN_parameter_list( ModuleParameterList& aParameterList )
    {
        // resize and initialize with standard parameters
        aParameterList.resize( 3 );
        aParameterList( 0 ).resize( 1 );
        aParameterList( 0 )( 0 ) = prm::create_gen_parameter_list(); // ParameterList(); 
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris

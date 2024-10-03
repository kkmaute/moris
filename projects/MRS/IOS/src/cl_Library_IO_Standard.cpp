/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Library_IO_Standard.cpp
 *
 */

#include "cl_Library_IO_Standard.hpp"
#include "enums.hpp"
#include "parameters.hpp"

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Standard::Library_IO_Standard()
            : Library_IO()    // initialize base class data as usual
    {
        // set the type of this library
        mLibraryType = Library_Type::STANDARD;

        // list of supported parameter list types
        mSupportedParamListTypes = {
            Parameter_List_Type::OPT,
            Parameter_List_Type::HMR,
            Parameter_List_Type::STK,
            Parameter_List_Type::XTK,
            Parameter_List_Type::GEN,
            Parameter_List_Type::FEM,
            Parameter_List_Type::SOL,
            Parameter_List_Type::MSI,
            Parameter_List_Type::VIS,
            Parameter_List_Type::MIG,
            Parameter_List_Type::WRK,
            Parameter_List_Type::MORISGENERAL
        };
    }

    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Standard::~Library_IO_Standard()
    {
        // do nothing extra
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Standard::finalize( std::string aFilePath )
    {
        // check that an .xml input file has been specified
        MORIS_ERROR( mSoLibIsInitialized || mXmlParserIsInitialized,
                "Library_IO_Standard::finalize() - Neither an .xml nor a .so input file has been specified. "
                "At least one input file is required." );

        // load the standard parameters into the member variables
        this->load_all_standard_parameters();

        // if an .so file has been parsed, first use its parameters (if any were defined in it) to overwrite or add to the standard parameters
        if ( mSoLibIsInitialized )
        {
            this->load_parameters_from_shared_object_library();
        }

        // load parameters from xml, overwrites parameters specified in either the standard parameters or an .so file if parsed
        if ( mXmlParserIsInitialized )
        {
            // this->load_parameters_from_xml();
        }

        // check the parameters for validity
        this->check_parameters();

        // mark this library as finalized and lock it from modification
        mLibraryIsFinalized = true;

        // print receipt of the finalized library
        // this->print_parameter_receipt( "./Parameter_Receipt.xml" ); // TODO: the file name and location should be user defineable
        this->print_parameter_receipt( aFilePath );    // TODO: the file name and location should be user defineable
    }

    //------------------------------------------------------------------------------------------------------------------
    // STANDARD PARAMETER LIST FUNCTIONS
    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Standard::load_all_standard_parameters()
    {
        // go over all modules and check that
        for ( uint iModule = 0; iModule < (uint)( Parameter_List_Type::END_ENUM ); iModule++ )
        {
            // get the current module
            Parameter_List_Type tParamListType = (Parameter_List_Type)( iModule );

            // get access to this module's parameter list
            ModuleParameterList tModuleParamList = mParameterLists( iModule );

            // most modules must use a 1x1 parameter list, so resize here to this default
            tModuleParamList.resize( 1 );

            // for each parameter list type, initialize it with a default
            switch ( tParamListType )
            {
                case Parameter_List_Type::OPT:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );
                    break;

                case Parameter_List_Type::HMR:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );
                    break;

                case Parameter_List_Type::STK:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_stk_parameter_list() );
                    break;

                case Parameter_List_Type::XTK:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
                    break;

                case Parameter_List_Type::GEN:
                    tModuleParamList.resize( 3 );
                    tModuleParamList( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
                    break;

                case Parameter_List_Type::FEM:
                    tModuleParamList.resize( 9 );
                    tModuleParamList( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
                    break;

                case Parameter_List_Type::SOL:
                    tModuleParamList.resize( 8 );
                    tModuleParamList( 0 ).add_parameter_list( prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
                    tModuleParamList( 1 ).add_parameter_list( prm::create_linear_solver_parameter_list() );
                    tModuleParamList( 2 ).add_parameter_list( prm::create_nonlinear_algorithm_parameter_list() );
                    tModuleParamList( 3 ).add_parameter_list( prm::create_nonlinear_solver_parameter_list() );
                    tModuleParamList( 4 ).add_parameter_list( prm::create_time_solver_algorithm_parameter_list() );
                    tModuleParamList( 5 ).add_parameter_list( prm::create_time_solver_parameter_list() );
                    tModuleParamList( 6 ).add_parameter_list( prm::create_solver_warehouse_parameterlist() );
                    tModuleParamList( 7 ).add_parameter_list( prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
                    break;

                case Parameter_List_Type::MSI:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
                    break;

                case Parameter_List_Type::VIS:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
                    tModuleParamList( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
                    tModuleParamList( 0 )( 0 ).set( "Save_Frequency", 1 );
                    break;

                case Parameter_List_Type::MIG:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_mig_parameter_list() );
                    break;

                case Parameter_List_Type::WRK:
                    tModuleParamList( 0 ).add_parameter_list( prm::create_wrk_parameter_list() );
                    break;

                case Parameter_List_Type::MORISGENERAL:
                    tModuleParamList.resize( 3 );
                    break;

                case Parameter_List_Type::END_ENUM:
                    MORIS_ERROR( false,
                            "Library_IO_Standard::load_all_standard_parameters() - "
                            "Parameter_List_Type is UNDEFINED. This loop shouldn't get here." );
                    break;

                default:
                    MORIS_ERROR( false,
                            "Library_IO_Standard::load_all_standard_parameters() - "
                            "No standard library defined for module parameter list type #%i. "
                            "Most likely a new enum has been defined that isn't used in this function.",
                            (uint)( tParamListType ) );
                    break;
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris

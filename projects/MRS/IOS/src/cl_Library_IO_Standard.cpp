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

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Standard::Library_IO_Standard()
            : Library_IO() // initialize base class data as usual
    {
        // set the type of this library
        mLibraryType = Library_Type::STANDARD;
    }

    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Standard::~Library_IO_Standard()
    {
        // do nothing extra
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO_Standard::finalize()
    {
        // check that an .xml input file has been specified
        MORIS_ERROR( mSoLibIsInitialized,  // || mXmlParserIsInitialized , // TODO
                "Library_IO_Standard::finalize() - Neither an .xml nor a .so input file has been specified. "
                "At least one input file is required." );

        // load the standard parameters into the member variables
        this->load_all_standard_parameters();

        // if an .so file has been parsed, first use its parameters (if any were defined in it) to overwrite or add to the standard parameters
        if( mSoLibIsInitialized )
        {
            // go through the various parameter list names and see if they exist in the provided .so file
            for( uint iParamListType = 0; iParamListType < (uint)( Parameter_List_Type::END_ENUM ); iParamListType++ )
            {
                // get the current enum
                Parameter_List_Type tParamListType = (Parameter_List_Type)( iParamListType );
                
                // get the name of the parameter list function
                std::string tParamListFuncName = get_name_for_parameter_list_type( tParamListType );

                // see if a function for this parameter list function exists in the provide input file
                Parameter_Function tUserDefinedParmListFunc = reinterpret_cast< Parameter_Function >( dlsym( mLibraryHandle, tParamListFuncName.c_str() ) );
                bool tParamListFuncExists = ( tUserDefinedParmListFunc != nullptr );

                // if the parameter list function exists, use it to overwrite and add to the standard parameters
                if( tParamListFuncExists )
                {
                    // log that the parameter list has been recognized
                    MORIS_LOG( "Parameters provided for %s in .so file.", convert_parameter_list_enum_to_string( tParamListType ).c_str() );

                    // throw out a warning if unknown parameter list types are used
                    if( mSupportedParamListTypes.find( tParamListType ) == mSupportedParamListTypes.end() )
                    {
                        MORIS_LOG( "These parameters are irrelevant for mesh generation and will be ignored." );
                    }
                    else // otherwise, if parameter list is supported, overwrite and add parameters to standard parameters
                    {
                        // create a copy of the parameter list the .so file provided
                        ModuleParameterList tUserDefinedParamList;
                        tUserDefinedParmListFunc( tUserDefinedParamList );

                        // get the parameter list to the currently 
                        ModuleParameterList& tCurrentModuleStandardParamList = mParameterLists( iParamListType );

                        // superseed standard parameters with user-defined parameters
                        this->overwrite_and_add_parameters( tCurrentModuleStandardParamList, tUserDefinedParamList );
                    }
                } 

            } // end for: parameter list types that could be specified

        } // end if: SO library is initialized

        // load parameters from xml
        // TODO: this->load_parameters_from_xml();

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
    Library_IO_Standard::load_all_standard_parameters()
    {
        // FIXME: use PRM functions once everything is moved
        // do nothing for now
    }

    //------------------------------------------------------------------------------------------------------------------

} // namespace moris

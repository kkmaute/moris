/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Library_IO.cpp
 *
 */

#include "cl_Library_IO.hpp"
#include "cl_XML_Parser.hpp"

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    Library_IO::Library_IO()
            : mSoFilePath( "" )
            , mLibraryHandle( nullptr )
            , mSoLibIsInitialized( false )
            , mXmlFilePath( "" )
            , mXmlReader( std::make_unique< XML_Parser >() )
            , mXmlParserIsInitialized( false )
            , mLibraryIsFinalized( false )
            , mLibraryType( Library_Type::UNDEFINED )                       // base class library-type is undefined
            , mParameterLists( (uint)( Parameter_List_Type::END_ENUM ) )    // list of module parameter lists sized to the number of modules that exist
            , mXmlWriter( std::make_unique< XML_Parser >() )
            , mSupportedParamListTypes()
    {
        // do nothing else
    }

    //------------------------------------------------------------------------------------------------------------------

    Library_IO::~Library_IO()
    {
        // close handle to shared object library if it has been opened
        if ( mSoLibIsInitialized )
        {
            dlclose( mLibraryHandle );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    void*
    Library_IO::get_shared_object_library_handle()
    {
        MORIS_ASSERT( mSoLibIsInitialized,
                "Library_IO::get_shared_object_library_handle() - "
                "A shared object library has not been initialized." );
        return mLibraryHandle;
    }

    //------------------------------------------------------------------------------------------------------------------

    std::string
    Library_IO::convert_to_absolute_file_path( const std::string& aFilePath )
    {
        // check the first letter of file path
        if ( aFilePath.at( 0 ) == '/' )    // this is already an absolute path
        {
            // just return the same file path
            return aFilePath;
        }
        else    // this is a relative path
        {
            // get the current absolute working directory path
            std::string tCurrentDir = std::string( std::getenv( "PWD" ) );

            // add the directory to the file path
            return tCurrentDir + "/" + aFilePath;
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    bool
    Library_IO::check_if_parameter_list_function_name( const std::string& aFunctionName )
    {
        std::string tParamListEnding  = "ParameterList";
        size_t      tNumCharsInEnding = tParamListEnding.length();

        bool tIsParameterFunction = false;

        if ( aFunctionName.length() > tNumCharsInEnding )
        {
            size_t tPos          = aFunctionName.length() - tNumCharsInEnding;
            int    tCompare      = aFunctionName.compare( tPos, tNumCharsInEnding, tParamListEnding );
            tIsParameterFunction = ( tCompare == 0 );
        }

        return tIsParameterFunction;
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO::overwrite_and_add_parameters(
            ModuleParameterList& aParamListToModify,
            ModuleParameterList& aParamListToAdd )
    {
        // get the sizes of the different parameter lists
        uint tOuterSizeToMod = aParamListToModify.size();
        uint tOuterSizeToAdd = aParamListToAdd.size();

        // resize the outer cell, if the parameter list to add is longer
        if ( tOuterSizeToAdd > tOuterSizeToMod )
        {
            aParamListToModify.resize( tOuterSizeToAdd );
        }

        // resize the inner cells, if the parameter lists to add are longer than the original ones
        for ( uint iOuterCell = 0; iOuterCell < tOuterSizeToAdd; iOuterCell++ )
        {
            uint tInnerSizeToMod = aParamListToModify( iOuterCell ).size();
            uint tInnerSizeToAdd = aParamListToAdd( iOuterCell ).size();

            if ( tInnerSizeToAdd > tInnerSizeToMod )
            {
                aParamListToModify( iOuterCell ).resize( tInnerSizeToAdd );
            }
        }

        // go over the various pieces of the parameter list and add non-existing parameters,
        // or overwrite them with the parameters provided
        for ( uint iOuterCell = 0; iOuterCell < tOuterSizeToAdd; iOuterCell++ )
        {
            uint tInnerSizeToAdd = aParamListToAdd( iOuterCell ).size();

            for ( uint iInnerCell = 0; iInnerCell < tInnerSizeToAdd; iInnerCell++ )
            {
                // get access to the current parameter lists
                Parameter_List& tParamsToMod = aParamListToModify( iOuterCell )( iInnerCell );
                Parameter_List& tParamsToAdd = aParamListToAdd( iOuterCell )( iInnerCell );

                // if the existing parameter list is empty, just replace it with whatever is in the one to add
                if ( tParamsToMod.is_empty() )    // FIXME, potentially
                {
                    tParamsToMod = tParamsToAdd;
                }
                else    // otherwise compare and add/overwrite values
                {
                    tParamsToMod.copy_parameters( tParamsToAdd );
                }
            }    // end for: inner cells
        }    // end for: outer cells
    }

    //------------------------------------------------------------------------------------------------------------------

    std::string
    Library_IO::get_path( File_Type aFileType ) const
    {
        // switch between the file types that could be requested
        switch ( aFileType )
        {
            // shared object library file
            case File_Type::SO_FILE:

                MORIS_ASSERT( mSoLibIsInitialized,
                        "Library_IO::get_path() - "
                        "Trying to get the path to the .so file used, but no shared object library has been initialized." );
                return mSoFilePath;

            // xml input file
            case File_Type::XML_FILE:

                MORIS_ASSERT( mXmlParserIsInitialized,
                        "Library_IO::get_path() - "
                        "Trying to get the path to the .xml file used, but no XML parser has been initialized." );
                return mXmlFilePath;

            // unknown file type to the base class
            default:
                MORIS_ERROR( false, "Library_IO_MeshGen::get_path() - File type unknown." );
                return "";
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO::load_parameter_list(
            const std::string& aFileName,
            File_Type          aFileType )
    {
        // check that this library has not been fully initialized yet and isn't locked
        MORIS_ERROR( !mLibraryIsFinalized,
                "Library_IO::load_parameter_list() - "
                "This Library has already been finalized and cannot load any additional parameters." );

        // initialization procedure for various input file types
        switch ( aFileType )
        {
            /* -------------------------- */
            case File_Type::SO_FILE:

                // check that no shared object library has been initialized yet
                MORIS_ASSERT( !mSoLibIsInitialized,
                        "Library_IO::load_parameter_list() - "
                        "Trying to initialize a shared object library, but one has already been initialized." );

                // get and store the absolute file path to the .so input file
                mSoFilePath = this->convert_to_absolute_file_path( aFileName );

                // try to open library file
                mLibraryHandle = dlopen( mSoFilePath.c_str(), RTLD_NOW );

                // test if loading succeeded
                if ( !mLibraryHandle )
                {
                    // get error string
                    std::string tError = dlerror();

                    // throw error
                    MORIS_ERROR( mLibraryHandle, "%s", tError.c_str() );
                }

                // if loading succeeded set the shared object library to initialized
                mSoLibIsInitialized = true;

                // stop switch case here
                break;

            /* -------------------------- */
            case File_Type::XML_FILE:

                // check that no other xml file has been initialized yet
                MORIS_ASSERT( !mXmlParserIsInitialized,
                        "Library_IO::load_parameter_list() - "
                        "An XML file has already been initialized. Only one xml file accepted as an input argument." );

                // get and store the absolute file path to the .xml input file
                mXmlFilePath = this->convert_to_absolute_file_path( aFileName );

                // load the xml file in the parser
                mXmlReader->initialize_read( mXmlFilePath );

                // mark the xml-parser as initialized
                mXmlParserIsInitialized = true;

                // stop switch case here
                break;

            /* -------------------------- */
            default:
                MORIS_ERROR( false, "Library_IO::load_parameter_list() - File type unknown." );
                break;

        }    // end: switch( aFileType )
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO::load_parameters_from_shared_object_library()
    {
        // check that an .so file has been initialized
        MORIS_ERROR( mSoLibIsInitialized, "Library_IO::load_parameters_from_shared_object_library() - No .so file has been loaded." );

        // go through the various parameter list names and see if they exist in the provided .so file
        for ( uint iParamListType = 0; iParamListType < (uint)( Parameter_List_Type::END_ENUM ); iParamListType++ )
        {
            // get the current enum
            Parameter_List_Type tParamListType = (Parameter_List_Type)( iParamListType );

            // get the name of the parameter list function
            std::string tParamListFuncName = get_name_for_parameter_list_type( tParamListType );

            // see if a function for this parameter list function exists in the provide input file
            Parameter_Function tUserDefinedParmListFunc = reinterpret_cast< Parameter_Function >( dlsym( mLibraryHandle, tParamListFuncName.c_str() ) );
            bool               tParamListFuncExists     = ( tUserDefinedParmListFunc != nullptr );

            // if the parameter list function exists, use it to overwrite and add to the standard parameters
            if ( tParamListFuncExists )
            {
                // log that the parameter list has been recognized
                MORIS_LOG( "Parameters for %s provided in .so file.", convert_parameter_list_enum_to_string( tParamListType ).c_str() );

                // throw out a warning if unknown parameter list types are used
                if ( mSupportedParamListTypes.find( tParamListType ) == mSupportedParamListTypes.end() )
                {
                    MORIS_LOG( "These parameters are irrelevant for chosen workflow and will be ignored." );
                }
                else    // otherwise, if parameter list is supported, overwrite and add parameters to standard parameters
                {
                    // create a copy of the parameter list the .so file provided
                    ModuleParameterList tUserDefinedParamList;
                    tUserDefinedParmListFunc( tUserDefinedParamList );

                    // get the parameter list to the currently
                    ModuleParameterList& tCurrentModuleStandardParamList = mParameterLists( iParamListType );

                    // supersede standard parameters with user-defined parameters
                    this->overwrite_and_add_parameters( tCurrentModuleStandardParamList, tUserDefinedParamList );
                }
            }

        }    // end for: parameter list types that could be specified
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO::load_parameters_from_xml()
    {
        // TODO ...
        //         // check that an XML file has been initialized
        //         MORIS_ERROR( mXmlParserIsInitialized, "Library_IO::load_parameters_from_xml() - No XML file has been loaded." );
        //
        //         // go through the various parameter list names and see if they exist in the provided xml file
        //         for( uint iParamListType = 0; iParamListType < (uint)( Parameter_List_Type::UNDEFINED ); iParamListType++ )
        //         {
        //             // get the current enum
        //             Parameter_List_Type tParamListType = (Parameter_List_Type)( iParamListType );
        //
        //             // get the name of the module, e.g. "HMR"
        //             std::string tModuleName = convert_parameter_list_enum_to_string( tParamListType );
        //
        //             // find this module name in the XML-tree
        //             size_t tCount = mXmlReader->count_keys_in_subtree( XML_PARAMETER_FILE_ROOT, tModuleName );
        //
        //             // Make sure that a module is not listed twice
        //             MORIS_ERROR( tCount < 2, "Library_IO::load_parameters_from_xml() - "
        //                     "Module '%s' has been specified %i times in XML input file. Only one definition allowed.",
        //                     tModuleName,
        //                     tCount );
        //
        //             // check if an entry for this Module exists in the XML file, if not, skip to next module
        //             if( tCount == 0 )
        //             {
        //                 continue;
        //             }
        //
        //             // get the number of sub-parameter lists that could be specified
        //             uint tMaxNumSubParamLists = get_number_of_sub_parameter_lists_in_module( tParamListType );
        //
        //             // get the root of the current Module parameter list
        //             std::string tModuleRoot = XML_PARAMETER_FILE_ROOT + "." + tModuleName;
        //
        //             // go over each of the sub-parameter lists
        //             for( uint iSubParamList = 0; iSubParamList < tMaxNumSubParamLists; iSubParamList++ )
        //             {
        //                 // get the name for this sub-parameter list
        //                 std::string tOuterSubParamListName = get_outer_sub_parameter_list_name( tParamListType, iSubParamList );
        //
        //                 // check if the sub-parameter list is found
        //                 size_t tSubParamListCount = mXmlReader->count_keys_in_subtree( tModuleRoot, tOuterSubParamListName );
        //                 MORIS_ERROR( tSubParamListCount < 2, "Library_IO::load_parameters_from_xml() - "
        //                         "Sub-parameter list '%s' in module '%s' has been specified %i times in XML input file. Only one definition allowed.",
        //                         tOuterSubParamListName,
        //                         tModuleName,
        //                         tSubParamListCount );
        //
        //                 // if the sub-parameter list is missing skip everything here after
        //                 if ( tSubParamListCount == 0 )
        //                 {
        //                     continue;
        //                 }
        //
        //                 // get the names of the inner parameter lists for this sub-parameter list
        //                 std::string tInnerSubParamListName = get_inner_sub_parameter_list_name( tParamListType, iSubParamList );
        //
        //                 // if there is no name for the inner parameter list, this is just a general parameter list without inner children
        //
        //                 std::string tInnerSubParamListRoot = tModuleRoot + "." + tOuterSubParamListName;
        //                 if( tInnerSubParamListName != "" )
        //                 {
        //                     tInnerSubParamListRoot = tInnerSubParamListRoot + "." + tInnerSubParamListName;
        //                 }
        //
        //                 // get the name of the
        //
        //                 // count the inner sub-parameter lists if there are any
        //
        //             }
        //
        //         }
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO::print_parameter_receipt( const std::string& aOutputFileName )
    {
        // initialize the xml writer by defining the root
        mXmlWriter->initialize_write( aOutputFileName );

        // write root of the tree
        // mXmlWriter->flush_buffer_to_tree( XML_PARAMETER_FILE_ROOT );

        // go through the modules and print their parameters to file
        for ( uint iModule = 0; iModule < (uint)( Parameter_List_Type::END_ENUM ); iModule++ )
        {
            // get the enum and name of the current module
            Parameter_List_Type tModule = (Parameter_List_Type)( iModule );

            // write this module to the xml tree
            this->write_module_parameter_list_to_xml_tree( tModule );
        }

        // write the xml file
        mXmlWriter->save();
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO::write_module_parameter_list_to_xml_tree( const Parameter_List_Type aModule )
    {
        // get the location of the current
        uint tModuleIndex = (uint)( aModule );

        // get this module's parameter list
        ModuleParameterList& tModuleParamList    = mParameterLists( tModuleIndex );
        uint                 tOuterParamListSize = tModuleParamList.size();

        // go through the individual sub-parameter lists and write them to the file
        for ( uint iOuterSubParamList = 0; iOuterSubParamList < tOuterParamListSize; iOuterSubParamList++ )
        {
            // get the number of inner parameter sub lists
            uint tInnerParamListSize = tModuleParamList( iOuterSubParamList ).size();

            for ( uint iInnerSubParamList = 0; iInnerSubParamList < tInnerParamListSize; iInnerSubParamList++ )
            {
                // add an index if there are multiple
                if ( tInnerParamListSize > 1 )
                {
                    mXmlWriter->set_attribute_in_buffer( "ind", iInnerSubParamList );
                }

                // store the current parameter list in xml writer buffer
                this->write_parameter_list_to_xml_buffer( tModuleParamList( iOuterSubParamList )( iInnerSubParamList ) );

                // get the path where the buffered parameter list should be stored
                std::string tParamListLocation = this->get_sub_parameter_list_location_in_xml_tree( aModule, iOuterSubParamList, true );

                // store the parameter list to the xml tree
                mXmlWriter->flush_buffer_to_tree( tParamListLocation );
            }
        }
    }

    // -----------------------------------------------------------------------------

    void
    Library_IO::write_parameter_list_to_xml_buffer( Parameter_List& aParameterList )
    {
        // go over the entries of the parameter list ...
        for ( const Parameter_Iterator& iParamToAdd : aParameterList )
        {
            // ... and add or modify them
            mXmlWriter->set_in_buffer( iParamToAdd.get_name(), iParamToAdd.get_parameter().get_string() );
        }
    }

    // -----------------------------------------------------------------------------

    std::string
    Library_IO::get_sub_parameter_list_location_in_xml_tree(
            const Parameter_List_Type aModule,
            const uint                aSubParamListIndex,
            const bool                aIsInnerParamList )
    {
        // initialize the location with the root of the xml tree
        std::string tLocation = XML_PARAMETER_FILE_ROOT;

        // get the name of the module and add it to the location
        std::string tModuleName = convert_parameter_list_enum_to_string( aModule );
        tLocation               = tLocation + "." + tModuleName;

        // if a sub-parameter list index has been provided
        if ( aSubParamListIndex < MORIS_UINT_MAX )
        {
            // get the name of the sub-parameter list
            std::string tSubParamListName = get_outer_sub_parameter_list_name( aModule, aSubParamListIndex );

            // add it to the location if not empty, otherwise don't add anything
            if ( tSubParamListName != "" )
            {
                tLocation = tLocation + "." + tSubParamListName;
            }
        }

        if ( aIsInnerParamList )
        {
            // get the name of the sub-parameter list
            std::string tInnerSubParamListName = get_inner_sub_parameter_list_name( aModule, aSubParamListIndex );

            // add it to the location if not empty, otherwise don't add anything
            if ( tInnerSubParamListName != "" )
            {
                tLocation = tLocation + "." + tInnerSubParamListName;
            }
        }

        // return the location path
        return tLocation;
    }

    //------------------------------------------------------------------------------------------------------------------

    ModuleParameterList
    Library_IO::get_parameters_for_module( Parameter_List_Type aParamListType ) const
    {
        // check that the parameter lists are complete
        MORIS_ERROR( mLibraryIsFinalized,
                "Library_IO::get_parameters_for_module() - "
                "Library has not been fully initialized. "
                "The Library needs to be finalized before parameters can be loaded." );

        // get the parameter list for the module and return it
        uint tParamListIndex = (uint)( aParamListType );
        return mParameterLists( tParamListIndex );
    }

    //------------------------------------------------------------------------------------------------------------------
    void Library_IO::check_parameters()
    {
        // Loop over all modules
        for ( uint iModuleIndex = 0; iModuleIndex < mParameterLists.size(); iModuleIndex++ )
        {
            // Loop over module sub-vectors
            for ( uint iOuterIndex = 0; iOuterIndex < mParameterLists( iModuleIndex ).size(); iOuterIndex++ )
            {
                // Loop over parameter lists in sub-vector
                for ( uint iInnerIndex = 0; iInnerIndex < mParameterLists( iModuleIndex )( iOuterIndex ).size(); iInnerIndex++ )
                {
                    // Get parameter list
                    const Parameter_List& tParameterList = mParameterLists( iModuleIndex )( iOuterIndex )( iInnerIndex );

                    // Loop over mapped parameters
                    for ( const Parameter_Iterator& iParameterPair : tParameterList )
                    {
                        // Get external validator
                        const External_Validator& tExternalValidator = iParameterPair.get_parameter().get_external_validator();

                        // Go through validation cases
                        switch ( tExternalValidator.mValidationType )
                        {
                            case Validation_Type::NONE:
                            {
                                // Do nothing
                                break;
                            }
                            case Validation_Type::SELECTION:
                            {
                                // Build internal variants that need to be checked
                                Vector< Variant > iInternalVariants = split_variant( iParameterPair.get_parameter().get_value() );

                                // Get valid external variants
                                Vector< Variant > tExternalVariants = this->get_external_variants( tExternalValidator, tParameterList );

                                // Loop over all internal variants
                                bool        tAllMatchesFound = true;
                                std::string tInternalVariantNotFound;
                                for ( const Variant& iInternalOption : iInternalVariants )
                                {
                                    // Check options for a match
                                    bool tMatchFound = false;
                                    for ( const Variant& iExternalOption : tExternalVariants )
                                    {
                                        if ( iInternalOption == iExternalOption )
                                        {
                                            tMatchFound = true;
                                            break;
                                        }
                                    }

                                    // If match not found, exit
                                    if ( not tMatchFound )
                                    {
                                        tInternalVariantNotFound = convert_variant_to_string( iInternalOption );
                                        tAllMatchesFound         = false;
                                        break;
                                    }
                                }

                                // Error if no match was found
                                if ( not tAllMatchesFound )
                                {
                                    // Comma-separated list of options
                                    std::string tExternalOptionList;
                                    std::string tDelimiter;
                                    for ( const Variant& iExternalOption : tExternalVariants )
                                    {
                                        tExternalOptionList += tDelimiter + convert_variant_to_string( iExternalOption );
                                        tDelimiter = ", ";
                                    }

                                    // Additional info about where the checked options came from
                                    std::string tExternalValidatorString =
                                            "These selections are taken from parameter " + tExternalValidator.mParameterName + ", located in the ";
                                    if ( tExternalValidator.mParameterListType == Parameter_List_Type::END_ENUM )
                                    {
                                        tExternalValidatorString +=
                                                convert_parameter_list_enum_to_string( (Parameter_List_Type)iModuleIndex )
                                                + " parameters with outer vector index " + std::to_string( iOuterIndex ) + " and inner vector index " + std::to_string( iInnerIndex );
                                    }
                                    else
                                    {
                                        tExternalValidatorString +=
                                                convert_parameter_list_enum_to_string( tExternalValidator.mParameterListType )
                                                + " parameters with outer vector index " + std::to_string( tExternalValidator.mParameterListIndex );
                                    }

                                    // Execute error
                                    MORIS_ERROR( tAllMatchesFound,
                                            "Parameter %s was set with an invalid value, %s. It requires one of the following selections:\n\t%s\n(%s)",
                                            iParameterPair.get_name().c_str(),
                                            tInternalVariantNotFound.c_str(),
                                            tExternalOptionList.c_str(),
                                            tExternalValidatorString.c_str() );
                                }
                                break;
                            }
                            case Validation_Type::SIZE:
                            {
                                // Get valid external variants
                                Vector< Variant > tExternalVariants = this->get_external_variants( tExternalValidator, tParameterList );

                                // Check that we have single external variant (for now, this can be changed in the future)
                                MORIS_ERROR( tExternalVariants.size() == 1,
                                        "%lu variants were found for the external size validation of parameter %s.",
                                        tExternalVariants.size(),
                                        iParameterPair.get_name().c_str() );

                                bool tSizeMatch = get_size( iParameterPair.get_parameter().get_value() ) == get_size( tExternalVariants( 0 ) );
                                if ( not tSizeMatch )
                                {
                                    // Additional info about where the checked options came from
                                    std::string tExternalValidatorString =
                                            "This size is based on parameter " + tExternalValidator.mParameterName + ", located in the ";
                                    if ( tExternalValidator.mParameterListType == Parameter_List_Type::END_ENUM )
                                    {
                                        tExternalValidatorString +=
                                                convert_parameter_list_enum_to_string( (Parameter_List_Type)iModuleIndex )
                                                + " parameters with outer vector index " + std::to_string( iOuterIndex ) + " and inner vector index " + std::to_string( iInnerIndex );
                                    }
                                    else
                                    {
                                        tExternalValidatorString +=
                                                convert_parameter_list_enum_to_string( tExternalValidator.mParameterListType )
                                                + " parameters with outer vector index " + std::to_string( tExternalValidator.mParameterListIndex );
                                    }

                                    // Execute error
                                    MORIS_ERROR( tSizeMatch,
                                            "Parameter %s was set with an invalid value, %s. It requires a size of %u.\n(%s)",
                                            iParameterPair.get_name().c_str(),
                                            iParameterPair.get_parameter().get_string().c_str(),
                                            get_size( tExternalVariants( 0 ) ),
                                            tExternalValidatorString.c_str() );
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    Vector< Variant > Library_IO::get_external_variants(
            const External_Validator& aExternalValidator,
            const Parameter_List&     aContainingParameterList )
    {
        // Start with empty vector
        Vector< Variant > tExternalOptions;

        // Check if external validation is required
        if ( aExternalValidator.mParameterListType == Parameter_List_Type::END_ENUM )
        {
            // Return single parameter
            tExternalOptions = { aContainingParameterList.get_variant( aExternalValidator.mParameterName ) };
        }
        else
        {
            // Check if parameter list exists
            if ( mParameterLists( (uint)aExternalValidator.mParameterListType ).size() > 0 )
            {
                // Loop over parameter lists in external module sub-vector
                for ( const auto& iExternalParameterList : mParameterLists( (uint)aExternalValidator.mParameterListType )( aExternalValidator.mParameterListIndex ) )
                {
                    // Check if external parameter exists here
                    if ( iExternalParameterList.exists( aExternalValidator.mParameterName ) )
                    {
                        tExternalOptions.push_back( iExternalParameterList.get_variant( aExternalValidator.mParameterName ) );
                    }
                }
            }
        }

        // Return options
        return tExternalOptions;
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris

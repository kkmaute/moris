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
#include "cl_Library_Enums.hpp"
#include "cl_XML_Parser.hpp"
#include <cctype>
#include <cstddef>
#include <iterator>

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    // Declare helper function for reading a module's parameter lists
    Module_Parameter_Lists read_module( uint aRoot );

    //------------------------------------------------------------------------------------------------------------------

    Library_IO::Library_IO()
            : mSoFilePath( "" )
            , mLibraryHandle( nullptr )
            , mSoLibIsInitialized( false )
            , mXmlFilePath( "" )
            , mXmlReader( std::make_unique< XML_Parser >() )
            , mXmlParserIsInitialized( false )
            , mLibraryIsFinalized( false )
            , mXmlWriter( std::make_unique< XML_Parser >() )
            , mSupportedParamListTypes()
    {
        for ( uint iTypeIndex = 0; iTypeIndex < static_cast< uint >( Module_Type::END_ENUM ); iTypeIndex++ )
        {
            mParameterLists.push_back( Module_Parameter_Lists( static_cast< Module_Type >( iTypeIndex ) ) );
        }
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
            Module_Parameter_Lists& aParamListToModify,
            Module_Parameter_Lists& aParamListToAdd )
    {
        // get the sizes of the different parameter lists
        uint tOuterSizeToMod = aParamListToModify.size();
        uint tOuterSizeToAdd = aParamListToAdd.size();

        // resize the outer cell, if the parameter list to add is longer
        if ( tOuterSizeToAdd > tOuterSizeToMod )
        {
            // TODO if we want to handle this again, give to the module parameter list for copying
            MORIS_ERROR( false, "Parameter lists are not the same size" );
            // aParamListToModify.resize( tOuterSizeToAdd );
        }

        // resize the inner cells, if the parameter lists to add are longer than the original ones
        for ( uint iOuterCell = 0; iOuterCell < tOuterSizeToAdd; iOuterCell++ )
        {
            uint tInnerSizeToMod = aParamListToModify( iOuterCell ).size();
            uint tInnerSizeToAdd = aParamListToAdd( iOuterCell ).size();

            // Modify existing parameters
            for ( uint iInnerCell = 0; iInnerCell < tInnerSizeToMod; iInnerCell++ )
            {
                // get access to the current parameter lists
                Parameter_List& tParamsToMod = aParamListToModify( iOuterCell )( iInnerCell );
                Parameter_List& tParamsToAdd = aParamListToAdd( iOuterCell )( iInnerCell );

                tParamsToMod.copy_parameters( tParamsToAdd );
            }

            // Add parameters that don't exist yet
            for ( uint iInnerIndex = tInnerSizeToMod; iInnerIndex < tInnerSizeToAdd; iInnerIndex++ )
            {
                aParamListToModify( iOuterCell ).add_parameter_list( aParamListToAdd( iOuterCell )( iInnerIndex ) );
            }
        }
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

    void Library_IO::finalize( const std::string& aFilePath )
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
            this->load_parameters_from_xml();
        }

        // check the parameters for validity
        this->check_parameters();

        // mark this library as finalized and lock it from modification
        mLibraryIsFinalized = true;

        // print receipt of the finalized library
        this->print_parameter_receipt( aFilePath );
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO::load_parameters_from_shared_object_library()
    {
        // check that an .so file has been initialized
        MORIS_ERROR( mSoLibIsInitialized, "Library_IO::load_parameters_from_shared_object_library() - No .so file has been loaded." );

        // go through the various parameter list names and see if they exist in the provided .so file
        for ( uint iParamListType = 0; iParamListType < (uint)( Module_Type::END_ENUM ); iParamListType++ )
        {
            // get the current enum
            Module_Type tParamListType = (Module_Type)( iParamListType );

            // get the name of the parameter list function
            std::string tParamListFuncName = get_name_for_parameter_list_type( tParamListType );

            // see if a function for this parameter list function exists in the provide input file
            Parameter_Function tUserDefinedParamListFunc = reinterpret_cast< Parameter_Function >( dlsym( mLibraryHandle, tParamListFuncName.c_str() ) );
            bool               tParamListFuncExists      = ( tUserDefinedParamListFunc != nullptr );

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
                    tUserDefinedParamListFunc( mParameterLists( iParamListType ) );
                }
            }
            else
            {
                mParameterLists( iParamListType ).clear();
            }
        }    // end for: parameter list types that could be specified
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO::load_parameters_from_xml()
    {
        // check that an XML file has been initialized
        MORIS_ERROR( mXmlParserIsInitialized, "Library_IO::load_parameters_from_xml() - No XML file has been loaded." );

        // go through the various parameter list names and see if they exist in the provided xml file
        for ( uint iParamListType = 0; iParamListType < (uint)( Module_Type::END_ENUM ); iParamListType++ )
        {
            // get the current enum
            Module_Type tParamListType = (Module_Type)( iParamListType );

            // get the name of the module, e.g. "HMR"
            std::string tModuleName = convert_parameter_list_enum_to_string( tParamListType );

            // find this module name in the XML-tree
            size_t tCount = mXmlReader->count_keys_in_subtree( XML_PARAMETER_FILE_ROOT, tModuleName );

            // Make sure that a module is not listed twice
            // COMPILATION ERROR HERE (CHECK WITH ADAM)
            // MORIS_ERROR( tCount < 2,
            //         "Library_IO::load_parameters_from_xml() - "
            //         "Module '%s' has been specified %i times in XML input file. Only one definition allowed.",
            //         tModuleName,
            //         tCount );

            // check if an entry for this Module exists in the XML file, if not, skip to next module


            // get the number of sub-parameter lists that could be specified
            uint tMaxNumSubParamLists = get_number_of_sub_parameter_lists_in_module( tParamListType );
            if ( tParamListType == Module_Type::FEM )
            {
                tMaxNumSubParamLists = 8;
            }
            // temporary storage for the parameter lists, later addded to mParamterLists
            Module_Parameter_Lists tParameterList( tParamListType );

            // If there are no modules of this type, create a default parameter list
            if ( tCount == 0 )
            {
                // Loop over all the sub-modules in this module and create the respective parameter lists with the defaults
                // for ( uint iSubModule = 0; iSubModule < tMaxNumSubParamLists; iSubModule++ )
                // {
                //     tParameterList( iSubModule ).push_back( create_parameter_list( tParamListType, iSubModule, 0 ) );
                // }
                mParameterLists( iParamListType ) = tParameterList;
                continue;
            }

            // get the root of the current Module parameter list
            std::string tModuleRoot = XML_PARAMETER_FILE_ROOT + "." + tModuleName;

            // get all submodule names
            Vector< std::string > tSubmoduleNames = get_submodule_names( tParamListType );

            // go over each of the sub-parameter lists
            for ( uint iSubParamList = 0; iSubParamList < tMaxNumSubParamLists; iSubParamList++ )
            {
                // get the name for this sub-parameter list
                std::string tOuterSubParamListName = tSubmoduleNames( iSubParamList );

                // check if the sub-parameter list is found
                size_t tSubParamListCount = mXmlReader->count_keys_in_subtree( tModuleRoot, tOuterSubParamListName );

                // COMPILATION ERROR HERE (CHECK WITH ADAM)
                // MORIS_ERROR( tSubParamListCount < 2,
                //         "Library_IO::load_parameters_from_xml() - "
                //         "Sub-parameter list '%s' in module '%s' has been specified %i times in XML input file. Only one definition allowed.",
                //         tOuterSubParamListName,
                //         tModuleName,
                //         tSubParamListCount );

                // if the sub-parameter list is missing skip everything here after
                if ( tSubParamListCount == 0 )
                {
                    // if ( tOuterSubParamListName != "Geometries" && tOuterSubParamListName != "Algorithms" && tOuterSubParamListName != "Linear_Algorithm" && tOuterSubParamListName != "Interface" )
                    // {
                    //     tParameterList( iSubParamList ).push_back( create_parameter_list( tParamListType, iSubParamList, 0 ) );
                    // }
                    continue;
                }

                // get the names of the inner parameter lists for this sub-parameter list
                std::string tInnerSubParamListName = get_inner_sub_parameter_list_name( tParamListType, iSubParamList );

                // if there is no name for the inner parameter list, this is just a general parameter list without inner children

                std::string           tInnerSubParamListRoot = tModuleRoot + "." + tOuterSubParamListName;
                Vector< std::string > tKeys;
                Vector< std::string > tValues;
                if ( tInnerSubParamListName != "" )
                {
                    size_t tInnerSubParamListCount = mXmlReader->count_keys_in_subtree( tInnerSubParamListRoot, tInnerSubParamListName );
                    // tInnerSubParamListRoot         = tInnerSubParamListRoot + "." + tInnerSubParamListName;

                    if ( tInnerSubParamListCount == 1 )
                    {
                        // tKeys and tValues are filled with the keys and values of the inner sub-parameter list
                        mXmlReader->get_keys_from_subtree( tInnerSubParamListRoot, tInnerSubParamListName, 0, tKeys, tValues );

                        // getting the index of the inner sub-module type by reading from the XML file parameter list (used for special forms like "GEN/Geometry", "OPT/Algorithm" and "SOL/Linear_Algorithm")
                        uint tIndex = get_subchild_index_from_xml_list( tParamListType, tInnerSubParamListName, tKeys, tValues );
                        // Adding the parameter list with set values from the XML file to tParameterList (that is added to mParameterLists)
                        // If tParameterList is empty, create and set a new parameter list, else set the parameter list 
                        if ( tParameterList( iSubParamList ).empty() )
                        {
                            tParameterList( iSubParamList ).add_parameter_list( create_and_set_parameter_list( tParamListType, iSubParamList, tIndex, tKeys, tValues ) );
                        }
                        else
                        {
                            set_parameter_list( tParameterList( iSubParamList )( 0 ), tKeys, tValues );
                        }
                    }
                    else
                    {
                        // If tParameterList(iSubParamList) is empty, create and set a new parameter list for each inner sub-parameter list
                        // Else, set the parameter list for the first element of tParameterList(iSubParamList) and add the rest of the inner sub-parameter lists
                        if ( tParameterList( iSubParamList ).empty() )
                        {
                            for ( uint iInnerSubParamList = 0; iInnerSubParamList < tInnerSubParamListCount; iInnerSubParamList++ )
                        {
                            tKeys.clear();
                            tValues.clear();
                            // tKeys and tValues are filled with the keys and values of the inner sub-parameter list
                            mXmlReader->get_keys_from_subtree( tInnerSubParamListRoot, tInnerSubParamListName, iInnerSubParamList, tKeys, tValues );

                            // getting the index of the inner sub-module type by reading from the XML file parameter list (used for special forms like "GEN/Geometry", "OPT/Algorithm" and "SOL/Linear_Algorithm")
                            uint tIndex = get_subchild_index_from_xml_list( tParamListType, tInnerSubParamListName, tKeys, tValues );

                            // Adding the parameter list with set values from the XML file to tParameterList (that is added to mParameterLists)
                            tParameterList( iSubParamList ).add_parameter_list( create_and_set_parameter_list( tParamListType, iSubParamList, tIndex, tKeys, tValues ) );
                        }
                        }
                        else
                        {
                            mXmlReader->get_keys_from_subtree( tInnerSubParamListRoot, tInnerSubParamListName, 0, tKeys, tValues );
                            set_parameter_list( tParameterList( iSubParamList )( 0 ), tKeys, tValues );
                            for ( uint iInnerSubParamList = 1; iInnerSubParamList < tInnerSubParamListCount; iInnerSubParamList++ )
                        {
                            tKeys.clear();
                            tValues.clear();
                            // tKeys and tValues are filled with the keys and values of the inner sub-parameter list
                            mXmlReader->get_keys_from_subtree( tInnerSubParamListRoot, tInnerSubParamListName, iInnerSubParamList, tKeys, tValues );

                            // getting the index of the inner sub-module type by reading from the XML file parameter list (used for special forms like "GEN/Geometry", "OPT/Algorithm" and "SOL/Linear_Algorithm")
                            uint tIndex = get_subchild_index_from_xml_list( tParamListType, tInnerSubParamListName, tKeys, tValues );

                            // Adding the parameter list with set values from the XML file to tParameterList (that is added to mParameterLists)
                            tParameterList( iSubParamList ).add_parameter_list( create_and_set_parameter_list( tParamListType, iSubParamList, tIndex, tKeys, tValues ) );
                        }
                        }
                    }
                }
                else
                {
                    // If there is no inner sub-parameter list, just adding the set parameter list from the xml file to the Module_Parameter_Lists
                    // mXmlReader->get_keys_from_subtree( tInnerSubParamListRoot, tInnerSubParamListName, 0, tKeys, tValues );
                    mXmlReader->get_keys_from_subtree( tModuleRoot, tOuterSubParamListName, 0, tKeys, tValues );
                    if ( tParameterList( iSubParamList ).empty() )
                    {
                        tParameterList( iSubParamList ).add_parameter_list( create_and_set_parameter_list( tParamListType, iSubParamList, 0, tKeys, tValues ) );
                    }
                    else
                    {
                        set_parameter_list( tParameterList( iSubParamList )( 0 ), tKeys, tValues );
                    }
                }
            }
            // adding the Module_Parameter_Lists to the mParameterLists (which is a vector of Module_Parameter_Listss)
            mParameterLists( iParamListType ) = tParameterList;
            // this->overwrite_and_add_parameters(mParameterLists(iParamListType),tParameterList);
            MORIS_LOG( "Parameters for %s provided in .xml file.", convert_parameter_list_enum_to_string( (Module_Type)iParamListType ).c_str() );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO::create_new_module_parameterlist()
    {
        Vector< Module_Parameter_Lists > tParameterList;
        for ( uint iParamListType = 0; iParamListType < (uint)( Module_Type::END_ENUM ); iParamListType++ )
        {
            mParameterLists( iParamListType ) = read_module( iParamListType );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO::print_parameter_receipt( const std::string& aOutputFileName )
    {
        if ( not aOutputFileName.empty() )
        {
            // initialize the xml writer by defining the root
            mXmlWriter->initialize_write( aOutputFileName );

            // write root of the tree
            // mXmlWriter->flush_buffer_to_tree( XML_PARAMETER_FILE_ROOT );

            // go through the modules and print their parameters to file
            for ( uint iModule = 0; iModule < (uint)( Module_Type::END_ENUM ); iModule++ )
            {
                // get the enum and name of the current module
                Module_Type tModule = (Module_Type)( iModule );

                // write this module to the xml tree
                this->write_module_parameter_list_to_xml_tree( tModule );
            }

            // write the xml file
            mXmlWriter->save();
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    void
    Library_IO::write_module_parameter_list_to_xml_tree( const Module_Type aModule )
    {
        // get the location of the current
        uint tModuleIndex = (uint)( aModule );

        // get this module's parameter list
        Module_Parameter_Lists& tModuleParamList    = mParameterLists( tModuleIndex );
        uint                    tOuterParamListSize = tModuleParamList.size();

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
        for ( const auto& iParamToAdd : aParameterList )
        {
            // ... and add or modify them
            // check if the get_name() has whitespaces, if so, replace them with $
            std::string tName = iParamToAdd.get_name();
            std::replace( tName.begin(), tName.end(), ' ', '$' );

            mXmlWriter->set_in_buffer( tName, iParamToAdd.get_parameter().get_string() );
        }
    }

    // -----------------------------------------------------------------------------

    std::string
    Library_IO::get_sub_parameter_list_location_in_xml_tree(
            const Module_Type aModule,
            const uint        aSubParamListIndex,
            const bool        aIsInnerParamList )
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
            std::string tSubParamListName = get_submodule_names( aModule )( aSubParamListIndex );

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

    Module_Parameter_Lists
    Library_IO::get_parameters_for_module( Module_Type aParamListType ) const
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
                    for ( auto iParameterPair : tParameterList )
                    {
                        // Get external validator
                        const External_Validator& tExternalValidator = iParameterPair.get_parameter().get_external_validator();

                        // Check if parameter needs linking (has not been cross-validated yet)
                        if ( iParameterPair.get_parameter().needs_linking() )
                        {
                            // Go through entry types
                            switch ( iParameterPair.get_parameter().get_entry_type() )
                            {
                                case Entry_Type::FREE:
                                {
                                    // Do nothing
                                    break;
                                }
                                case Entry_Type::SELECTION:
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
                                        if ( tExternalValidator.mParameterListType == Module_Type::END_ENUM )
                                        {
                                            tExternalValidatorString +=
                                                    convert_parameter_list_enum_to_string( (Module_Type)iModuleIndex )
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
                                case Entry_Type::LINKED_SIZE_VECTOR:
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
                                        if ( tExternalValidator.mParameterListType == Module_Type::END_ENUM )
                                        {
                                            tExternalValidatorString +=
                                                    convert_parameter_list_enum_to_string( (Module_Type)iModuleIndex )
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
    }

    //------------------------------------------------------------------------------------------------------------------

    Vector< Variant > Library_IO::get_external_variants(
            const External_Validator& aExternalValidator,
            const Parameter_List&     aContainingParameterList )
    {
        // Start with empty vector
        Vector< Variant > tExternalOptions;

        // Check if external validation is required
        if ( aExternalValidator.mParameterListType == Module_Type::END_ENUM )
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

    // FREE FUNCTIONS

    /**
     * @brief get_subchild_index_from_xml_list - Get the index of the sub-module type from the XML file
     * @param tInnerSubParamListName - The name of the inner sub-parameter list
     * @param tKeys - The keys of the XML file parameter list
     * @param tValues - The values of the XML file parameter list
     * @return uint - The index of the sub-module type for special forms like "GEN/Geometry", "OPT/Algorithm" and "SOL/Linear_Algorithm", if not these forms, returns 0
     */
    uint get_subchild_index_from_xml_list( Module_Type iParameterListType, std::string tInnerSubParamListName, Vector< std::string >& tKeys, Vector< std::string >& tValues )
    {
        uint tIndex = 0;
        // Check if the paramter_list_type is GEN and the inner sub-parameter list name is Geometry or Field
        if ( iParameterListType == Module_Type::GEN && ( tInnerSubParamListName == "Geometry" || tInnerSubParamListName == "Field" ) )
        {
            // Finding field_type in the keys of the XML file parameter list
            auto it = std::find( tKeys.begin(), tKeys.end(), "field_type" );
            if ( it != tKeys.end() )
            {
                // Index is the value of field_type
                tIndex = (uint)std::stoi( tValues( std::distance( tKeys.begin(), it ) ) );
            }

            // tParameterList( iSubParamList ).push_back( create_and_set_parameter_list( tParamListType, iSubParamList, tFieldType, tKeys, tValues ) );
        }
        else if ( tInnerSubParamListName == "Algorithms" )
        {
            std::string tAlgorithm;

            // Finding algorithm in the keys of the XML file parameter list
            auto it = std::find( tKeys.begin(), tKeys.end(), "algorithm" );
            if ( it != tKeys.end() )
            {
                tAlgorithm = tValues( std::distance( tKeys.begin(), it ) );
                // Need to create cl_OPT_Enums and create the algorithms enum to check against that
                // For now, just using the string with a else if statement
                if ( tAlgorithm.find( "gcmma" ) != std::string::npos )
                {
                    tIndex = 0;
                }
                else if ( tAlgorithm.find( "lbfgs" ) != std::string::npos )
                {
                    tIndex = 1;
                }
                else if ( tAlgorithm.find( "sql" ) != std::string::npos )
                {
                    tIndex = 2;
                }
                else
                {
                    tIndex = 3;
                }
            }
        }
        else if ( tInnerSubParamListName == "LinearAlgorithm" )
        {
            // Finding solver_type in the keys of the XML file parameter list
            auto it = std::find( tKeys.begin(), tKeys.end(), "Solver_Implementation" );
            if ( it != tKeys.end() )
            {
                // Index is the value of Solver_Implementation

                // solver_type returns a string name of the solver type
                tIndex = convert_parameter_from_string_to_type< uint >( tValues( std::distance( tKeys.begin(), it ) ) );


                // // Transforming the string to uppercase
                // transform( tSolverType.begin(), tSolverType.end(), tSolverType.begin(), ::toupper );
                // //remove any \ or " from the string
                // tSolverType.erase( std::remove( tSolverType.begin(), tSolverType.end(), '\"' ), tSolverType.end() );
                // tSolverType.erase( std::remove( tSolverType.begin(), tSolverType.end(), '\\' ), tSolverType.end() );

                // if ( tSolverType == "AMESOS_MUMPS") {
                //     tIndex = (uint) (sol::SolverType::AMESOS_IMPL);
                // }
                // else if ( tSolverType == "GMRES")
                // else {

                // // Finding the index of the solver type in SolverType ENUM
                // auto iFind = std::find( sol::SolverType_String::values.begin(), sol::SolverType_String::values.end(), tSolverType );
                // if ( iFind != sol::SolverType_String::values.end() )
                // {
                //     tIndex = std::distance( sol::SolverType_String::values.begin(), iFind );
                // }
                // }
            }
        }
        else if ( tInnerSubParamListName == "Preconditioner" )
        {
            // Finding preconditioner_type in the keys of the XML file parameter list
            auto it = std::find( tKeys.begin(), tKeys.end(), "Preconditioner_Implementation" );
            if ( it != tKeys.end() )
            {
                // Index is the value of preconditioner_type
                tIndex = (uint)std::stoi( tValues( std::distance( tKeys.begin(), it ) ) );
            }
        }
        else
        {
            tIndex = 0;
        }

        return tIndex;
    }

    /**
     * @brief create_and_set_parameter_list - Calls the create_parameter_list function and sets the parameter list with the values from the XML file in the correct data type
     * @param aModule - Module in Module_Type enum type
     * @param aChild - The index of the sub-module
     * @param aSubChild - The index of the sub-module type for special forms like "GEN/Geometry", "OPT/Algorithm" and "SOL/Linear_Algorithm", if not these forms, then 0
     * @param tKeys - The keys of the XML file parameter list
     * @param tValues - The values of the XML file parameter list
     * @return Parameter_List - The parameter list with the set values from the XML file
     */

    Parameter_List create_and_set_parameter_list( Module_Type aModule,
            uint                                              aChild,
            uint                                              aSubChild,
            Vector< std::string >&                            aKeys,
            Vector< std::string >&                            aValues )
    {
        // Create the parameter list with default values
        Parameter_List tParameterList = create_parameter_list( aModule, aChild, aSubChild );

        // Loop through the default parameter list

        // Instead of looping through iElements, the following for loop should loop over the keys and check that against the parameter list
        // If the key is found in the parameter list, then set the value of the parameter list with the value from the XML file

        set_parameter_list( tParameterList, aKeys, aValues );
        
        return tParameterList;
    }

    void set_parameter_list( Parameter_List& tParameterList,
            Vector< std::string >& aKeys,
            Vector< std::string >& aValues )
    {
        for ( uint iIndex = 0; iIndex < aKeys.size(); iIndex++ )
        {
            // Find the key in the parameter list
            bool tCheck = false;
            for ( Parameter_List::iterator iElements : tParameterList )
            {
                if ( iElements.get_name() == aKeys( iIndex ) )
                {
                    Parameter& aParameter = iElements.get_parameter();
                    tCheck                = true;
                    set_new_parameter( aParameter, aKeys( iIndex ), aValues( iIndex ) );
                    break;
                }
            }

            // Skip if the key is not found in the parameter list
            if ( not tCheck )
            {
                Design_Variable tDesign_Variable = convert_parameter_from_string_to_type< real >( aValues( iIndex ) );

                tParameterList.insert( aKeys( iIndex ), tDesign_Variable );
                continue;
            }
        }
    }

    void set_new_parameter( Parameter& aParameter,
            std::string&               aKey,
            const std::string&         aValue )
    {
        if ( aParameter.is_locked() )
        {
            return;
        }

        // Replace $ with whitespace in key
        std::replace( aKey.begin(), aKey.end(), '$', ' ' );


        // Set the value based on the type of parameter
        if ( aParameter.index() == variant_index< bool >() )
        {
            bool tBool = convert_parameter_from_string_to_type< bool >( aValue );
            aParameter.set_value( aKey, tBool, false );
        }
        else if ( aParameter.index() == variant_index< uint >() )
        {
            uint tUint = convert_parameter_from_string_to_type< uint >( aValue );
            aParameter.set_value( aKey, tUint, false );
        }
        else if ( aParameter.index() == variant_index< sint >() )
        {
            sint tSint = convert_parameter_from_string_to_type< sint >( aValue );
            aParameter.set_value( aKey, tSint, false );
        }
        else if ( aParameter.index() == variant_index< real >() )
        {
            real tReal = convert_parameter_from_string_to_type< real >( aValue );
            aParameter.set_value( aKey, tReal, false );
        }
        else if ( aParameter.index() == variant_index< std::string >() )
        {
            std::string tString = aValue;
            // Strip the string of leading and trailing quotation marks
            tString.erase( std::remove( tString.begin(), tString.end(), '\"' ), tString.end() );
            aParameter.set_value( aKey, tString, false );
        }
        else if ( aParameter.index() == variant_index< std::pair< std::string, std::string > >() )
        {
            std::string tPairString = aValue;
            // Process the pair string
            tPairString.erase( std::remove( tPairString.begin(), tPairString.end(), '{' ), tPairString.end() );
            tPairString.erase( std::remove( tPairString.begin(), tPairString.end(), '}' ), tPairString.end() );
            tPairString.erase( std::remove( tPairString.begin(), tPairString.end(), ' ' ), tPairString.end() );
            size_t                                tPos  = tPairString.find( '#' );
            std::pair< std::string, std::string > tPair = { tPairString.substr( 0, tPos ), tPairString.substr( tPos + 1 ) };
            aParameter.set_value( aKey, tPair, false );
        }
        else if ( aParameter.index() == variant_index< Vector< uint > >() )
        {
            Vector< uint > tVec = string_to_vector< uint >( aValue );
            aParameter.set_value( aKey, tVec, false );
        }
        else if ( aParameter.index() == variant_index< Vector< sint > >() )
        {
            Vector< sint > tVec = string_to_vector< sint >( aValue );
            aParameter.set_value( aKey, tVec, false );
        }
        else if ( aParameter.index() == variant_index< Vector< real > >() )
        {
            Vector< real > tVec = string_to_vector< real >( aValue );
            aParameter.set_value( aKey, tVec, false );
        }
        else if ( aParameter.index() == variant_index< Vector< std::string > >() )
        {
            std::string tString = aValue;
            tString.erase( std::remove( tString.begin(), tString.end(), '{' ), tString.end() );
            tString.erase( std::remove( tString.begin(), tString.end(), '}' ), tString.end() );
            tString.erase( std::remove( tString.begin(), tString.end(), ' ' ), tString.end() );
            Vector< std::string > tVec = split_string( tString, "," );
            aParameter.set_value( aKey, tVec, false );
        }
        else
        {
            // For other types, process accordingly
            Design_Variable tDesignVariable = convert_parameter_from_string_to_type< real >( aValue );
            aParameter.set_value( aKey, tDesignVariable, false );
        }
    }

    /**
     * @brief convert_parameter_from_string_to_type - Converts the string value from the XML file to the correct data type
     * @tparam T - The data type of the parameter
     * @param aValue - The string value of the parameter from the XML file
     * @return T - The value of the parameter in the correct data type
     */

    template< typename T >
    T convert_parameter_from_string_to_type( const std::string& aValue )
    {
        return (T)std::stod( aValue );
    }

    // Specialization for bool
    template<>
    bool convert_parameter_from_string_to_type< bool >( const std::string& aValue )
    {
        if ( aValue == "0" || aValue == "false" || aValue == "False" || aValue == "FALSE" )
        {
            return false;
        }
        else if ( aValue == "1" || aValue == "true" || aValue == "True" || aValue == "TRUE" )
        {
            return true;
        }
        else
        {
            throw std::invalid_argument( "Invalid value for conversion to bool" );
        }
    }

    /**
     * @brief create_parameter_list - This function creates a parameter list for a given module, child, and sub-child
     * @param aModule - The project index
     * @param aChild - The sub-module index
     * @param aSubChild - Should be 0 unless a sub-module has inner types (for instance GEN/Geometries, OPT/Algorithms and SOL/LinearAlgorithms)
     * @return Parameter_List
     */
    Parameter_List create_parameter_list( Module_Type aModule, uint aChild, uint aSubChild )
    {
        /*
        function name: create_parameter_list
        parameters:
          Module_Type aModule (ENUM) -> this gives the project name
          uint aChild -> gives the child index
          uint aSubChild -> gives the Sub-Child (inner sub-module) index
        returns:
            QList <QStringList>
                the create_function returns a ParameterList object that is a type of map
                The 0th index of the QList gives the "keys" of the map
                The 1st index of the QList gives the default "values" of the map
        */
        switch ( aModule )
        {
            case Module_Type::OPT:
                switch ( aChild )
                {
                    case 0:
                        return prm::create_opt_problem_parameter_list();

                    case 1:
                        return prm::create_opt_interface_parameter_list();

                        // Commented out the Interface manager for now

                        // switch ( aSubChild )
                        // {
                        //     case 0:
                        //         return prm::create_opt_interface_parameter_list();

                        //         break;

                        //     case 1:
                        //         return prm::create_opt_interface_manager_parameter_list();

                        //         break;

                        //     default:
                        //         break;
                        // }

                        break;

                    case 2:
                        switch ( aSubChild )
                        {

                            // Eventually create an enum to check this
                            case 0:
                                return prm::create_gcmma_parameter_list();

                            case 1:
                                return prm::create_lbfgs_parameter_list();

                            case 2:
                                return prm::create_sqp_parameter_list();

                            case 3:
                                return prm::create_sweep_parameter_list();
                            default:
                                break;
                        }
                        break;

                    default:
                        break;
                }
                // Free
                break;

            case Module_Type::HMR:
                return prm::create_hmr_parameter_list();

            case Module_Type::STK:
                return prm::create_stk_parameter_list();

            case Module_Type::XTK:
                return prm::create_xtk_parameter_list();

            case Module_Type::GEN:
                switch ( aChild )
                {
                    case 0:
                    {
                        return prm::create_gen_parameter_list();
                    }

                    case 1:
                    {
                        if ( aSubChild <= (uint)gen::Field_Type::USER_DEFINED )
                        {
                            return prm::create_level_set_geometry_parameter_list( (gen::Field_Type)aSubChild );
                        }
                        else if ( aSubChild == (uint)gen::Field_Type::USER_DEFINED + 1 )
                        {
                            return prm::create_surface_mesh_geometry_parameter_list();
                        }
                        else
                        {
                            return prm::create_voxel_geometry_parameter_list();
                        }
                        break;
                    }
                    case 2:
                    {
                        return prm::create_gen_property_parameter_list( (gen::Field_Type)aSubChild );
                    }
                    default:
                    {
                        break;
                    }
                }
                break;

            case Module_Type::FEM:
                /*
                 * Set of Dropdowns for tParameterList[0] (property_name in FEM)
                 * //Dropdown
                 * PropDensity, PropYoungs, PropPoisson,
                 * PropCTE, PropRefTemp, PropConductivity,
                 * PropCapacity, PropDirichlet, PropSelectX,
                 * PropSelectY, PropSelectZ, PropInnerPressureLoad,
                 * Pro#include <QApplication>
                 * pOuterPressureLoad, PropOuterTemperature
                 */
                switch ( aChild )
                {
                    case 0:
                        return prm::create_property_parameter_list();

                    case 1:
                        return prm::create_constitutive_model_parameter_list();

                    case 2:
                        return prm::create_stabilization_parameter_parameter_list();

                    case 3:
                        return prm::create_IWG_parameter_list();

                    case 4:
                        return prm::create_IQI_parameter_list();

                    case 5:
                        return prm::create_computation_parameter_list();

                    case 6:
                        return prm::create_fem_field_parameter_list();

                    case 7:
                        return prm::create_phase_parameter_list();

                    case 8:
                        return prm::create_material_model_parameter_list();

                    default:
                        break;
                }


                break;

            case Module_Type::SOL:
                //
                switch ( aChild )
                {
                    case 0:

                        return prm::create_linear_algorithm_parameter_list( (sol::SolverType)aSubChild );

                    case 1:
                        return prm::create_linear_solver_parameter_list();

                    case 2:
                        return prm::create_nonlinear_algorithm_parameter_list();

                    case 3:
                        return prm::create_nonlinear_solver_parameter_list();

                    case 4:
                        return prm::create_time_solver_algorithm_parameter_list();

                    case 5:
                        return prm::create_time_solver_parameter_list();

                    case 6:
                        return prm::create_solver_warehouse_parameterlist();

                    case 7:
                        // Need to add Preconditioners
                        return moris::prm::create_preconditioner_parameter_list( (moris::sol::PreconditionerType)aSubChild );
                        break;

                    default:
                        break;
                }

                break;

            case Module_Type::MSI:
                return prm::create_msi_parameter_list();

            case Module_Type::VIS:
                return prm::create_vis_parameter_list();    //

            case Module_Type::MIG:
                return prm::create_mig_parameter_list();

            case Module_Type::WRK:
                return prm::create_wrk_parameter_list();

            case Module_Type::MORISGENERAL:
                switch ( aChild )
                {
                    case 0:
                    {
                        Parameter_List tParameterList = ( moris::prm::create_moris_general_parameter_list() );
                        moris::prm::create_remeshing_parameterlist( tParameterList );
                        tParameterList.set( "mode", "none" );
                        return tParameterList;
                    }
                    break;

                    case 1:
                    {
                        Parameter_List tParameterList = ( moris::prm::create_moris_general_parameter_list() );
                        moris::prm::create_refinement_parameterlist( tParameterList );
                        return tParameterList;
                    }
                    break;

                    case 2:
                    {
                        Parameter_List tParameterList = ( moris::prm::create_moris_general_parameter_list() );
                        moris::prm::create_mapping_parameterlist( tParameterList );
                        return tParameterList;
                    }
                    break;

                    default:
                        break;
                }
                break;

            default:
                break;
        }

        MORIS_ERROR( false, "Library_Enums::get_number_of_sub_parameter_lists_in_module() - Parameter list type enum unknown." );
        return Parameter_List( "" );
    }

    Module_Parameter_Lists read_module( uint aRoot )
    {
        // Create the 3d vector
        Module_Parameter_Lists tParameterList( static_cast< Module_Type >( aRoot ) );
        for ( uint iChild = 0; iChild < get_number_of_sub_parameter_lists_in_module( (Module_Type)aRoot ); iChild++ )
        {
            if ( ( aRoot == (uint)( Module_Type::OPT ) && iChild == (uint)( OPT_Submodule::ALGORITHMS ) )
                    || ( aRoot == (uint)( Module_Type::GEN ) && iChild == (uint)( GEN_Submodule::GEOMETRIES ) )
                    || ( aRoot == (uint)( Module_Type::SOL ) && iChild == (uint)( SOL_Submodule::LINEAR_ALGORITHMS ) ) )
            {
            }
            else
            {
                tParameterList( iChild ).add_parameter_list( create_parameter_list( (Module_Type)aRoot, iChild, 0 ) );
            }
        }

        // Resize based on the projects/sub-projects
        // Populate based on the parameter_list
        return tParameterList;
    }

}    // namespace moris

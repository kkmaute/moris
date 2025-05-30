/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Library_IO.hpp
 *
 */

#pragma once

#include <string>
#include <set>
#include <algorithm>
#include <cctype>
#include <sys/types.h>
#include "dlfcn.h"
#include "assert.hpp"

#include "cl_Library_Enums.hpp"
#include "cl_Vector.hpp"

#include "parameters.hpp"


namespace moris
{
    // -----------------------------------------------------------------------------

    // Forward declare the XML-parser
    class XML_Parser;

    // Define what a parameter function is
    typedef void ( *Parameter_Function )( Module_Parameter_Lists& aParameterList );

    // -----------------------------------------------------------------------------

    /**
     * Wrapper class for any inputs that may be provided
     */
    class Library_IO
    {
        // -----------------------------------------------------------------------------

      protected:
        // -----------------------------------------------------------------------------

        // data for shared object input file
        std::string mSoFilePath;            // path to library file
        void*       mLibraryHandle;         // handle to shared object library
        bool        mSoLibIsInitialized;    // flag whether initialization of the

        // data for xml input file
        std::string                   mXmlFilePath;               // path to the .xml input file
        std::unique_ptr< XML_Parser > mXmlReader;                 // XML_Parser
        bool                          mXmlParserIsInitialized;    // flag whether the xml parser has its file

        // flag indicating whether the library is complete and initialized
        bool mLibraryIsFinalized;

        // storage for the parameters for the various Modules
        Vector< Module_Parameter_Lists > mParameterLists;

        // XML parser for output
        std::unique_ptr< XML_Parser > mXmlWriter;

        // -----------------------------------------------------------------------------

        /**
         * @brief Get handle to any external shared libraries that the library might have
         *
         * @return void* external shared library handle
         */
        virtual void*
        get_shared_object_library_handle();

        // -----------------------------------------------------------------------------

        /**
         * @brief takes any type of file path and converts it to an absolute path
         *
         * @param aFilePath file path & name to be converted
         * @return std::string absolute file path
         */
        std::string
        convert_to_absolute_file_path( const std::string& aFilePath );

        // -----------------------------------------------------------------------------

        bool
        check_if_parameter_list_function_name( const std::string& aFunctionName );

        // -----------------------------------------------------------------------------

        /**
         * @brief Take a parameter list and modify it by adding or overwriting its parameter with those of a second
         *
         * @param aParamListToModify parameter list that will be modified with the next one
         * @param aParamsToAdd parameters to overwrite or add to the previous
         */
        void
        overwrite_and_add_parameters(
                Module_Parameter_Lists& aParamListToModify,
                Module_Parameter_Lists& aParamsToAdd );

      public:
        // -----------------------------------------------------------------------------

        /**
         * Default constructor doing nothing
         */
        Library_IO();

        // -----------------------------------------------------------------------------

        /**
         * Default destructor
         */
        virtual ~Library_IO();

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the path and filename to any input file that might have been specified
         *
         * @param aFileType type of input file to the the filename for
         * @return std::string path/filename
         */
        virtual std::string
        get_path( File_Type aFileType ) const;

        // -----------------------------------------------------------------------------

        Vector< Module_Parameter_Lists >&
        get_parameter_lists()
        {
            return mParameterLists;
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief give a parameter list file to the library and read it
         *
         * @param aFileName Name/path to the file to be used
         * @param aFileType type of the file to be used
         */
        virtual void
        load_parameter_list( const std::string& aFileName, File_Type aFileType );

        /**
         * Finalizes this library and locks it from modification.
         *
         * @param aFilePath Optional file path/name for printing out an XML parameter receipt
         */
        void finalize( const std::string& aFilePath = "" );

        /**
         * @brief loads parameters from an shared object library and overwrites any previously specified parameters by it
         */
        virtual void
        load_parameters_from_shared_object_library();

        // -----------------------------------------------------------------------------

        /**
         * @brief loads parameters from an xml file and overwrites any previously specified parameters by it
         */
        virtual void
        load_parameters_from_xml();

        // -----------------------------------------------------------------------------

        /**
         * @brief Uses the read() function to create a new module parameter list from scratch when no XML file is given
         */
        virtual void
        create_new_module_parameterlist();

        // -----------------------------------------------------------------------------

        /**
         * @brief print a summary of all parameters after finalizing the library
         *
         * @param aOutputFileName name of the xml file to print the parameters to
         */
        void
        print_parameter_receipt( const std::string& aOutputFileName );

        // -----------------------------------------------------------------------------

        void
        write_module_parameter_list_to_xml_tree( const Module_Type aModule );

        // -----------------------------------------------------------------------------

        void
        write_parameter_list_to_xml_buffer( Parameter_List& aParameterList );

        // -----------------------------------------------------------------------------

        std::string
        get_sub_parameter_list_location_in_xml_tree(
                const Module_Type aModule,
                const uint                aSubParamListIndex = MORIS_UINT_MAX,
                const bool                aIsInnerParamList  = false );

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the parameters for module object
         *
         * @param aParamListType
         * @return Module_Parameter_Lists
         */
        Module_Parameter_Lists
        get_parameters_for_module( Module_Type aParamListType ) const;

        // -----------------------------------------------------------------------------

        /**
         * Loads a function from a shared object file if it exists.
         * This function is in the base class as it is templated.
         *
         * @tparam Function_Type Function pointer type
         * @param aFunctionName Function name to look for in the file
         * @param aThrowError parameter to check if the list exists
         * @return Function pointer
         */
        template< typename Function_Type >
        Function_Type
        load_function( const std::string& aFunctionName, bool aThrowError = true )
        {
            // make sure the library is fully initialized before any functions are loaded from it
            MORIS_ERROR( mLibraryIsFinalized,
                    "Library_IO::load_function() - "
                    "Accessing a function from the Library before it has been finalized." );

            // make sure that parameter lists are not loaded through this function
            MORIS_ERROR( !( this->check_if_parameter_list_function_name( aFunctionName ) ),
                    "Library_IO::load_function() - "
                    "Trying to load parameters for %s directly from shared object library. Use 'get_parameters_for_module()' instead.",
                    aFunctionName.c_str() );

            // get a handle to the library handle
            void* tLibraryHandle = this->get_shared_object_library_handle();

            // get the pointer to the user-defined function
            Function_Type aUserFunction = reinterpret_cast< Function_Type >( dlsym( tLibraryHandle, aFunctionName.c_str() ) );

            // depending on the flag throw an error
            if ( aThrowError )
            {
                // make sure that loading succeeded
                MORIS_ERROR( aUserFunction != 0,
                        "Could not find symbol %s  within file %s",
                        aFunctionName.c_str(),
                        this->get_path( File_Type::SO_FILE ).c_str() );
            }

            // return function handle
            return aUserFunction;
        }

        /**
         * Checks all parameters by running external validations. Called during finalize().
         */
        void check_parameters();

        /**
         * Gets all of the external options for a parameter based on its external validator.
         *
         * @param aExternalValidator External validator of a parameter
         * @param aContainingParameterList Parameter list that contains the parameter with the external validator
         * @return Valid external options: If size zero, all options are valid.
         */
        Vector< Variant > get_external_variants(
                const External_Validator& aExternalValidator,
                const Parameter_List&     aContainingParameterList );

      private:

        /**
         * Gets if a given module is being supported by the current library.
         *
         * @return If module is supported
         */
        virtual bool is_module_supported( Module_Type aModuleType ) = 0;

    };    // class Library_IO

    // -----------------------------------------------------------------------------

    // FREE FUNCTIONS
    // These are for reading parameter lists from xml files and set to the correct type

    /**
     * @brief Checks whether a string has a certain ending. Good for checking file types.
     *
     * @param aString string to check the ending of
     * @param aEnding ending to check for
     * @return true/false whether this string has this ending
     */
    bool
    string_ends_with(
            std::string const & aString,
            std::string const & aEnding );

    /**
     * @brief get_subchild_index_from_xml_list - Get the index of the sub-module type from the XML file
     * @param tInnerSubParamListName - The name of the inner sub-parameter list
     * @param tKeys - The keys of the XML file parameter list
     * @param tValues - The values of the XML file parameter list
     * @return uint - The index of the sub-module type for special forms like "GEN/Geometry", "OPT/Algorithm" and "SOL/Linear_Algorithm", if not these forms, returns 0
     */
    uint get_subchild_index_from_xml_list( std::string tInnerSubParamListName, Vector< std::string >& aKeys, Vector< std::string >& aValues );

    /**
     * @brief convert_parameter_from_string_to_type - Converts the string value from the XML file to the correct data type
     * @tparam T - The data type of the parameter
     * @param aValue - The string value of the parameter from the XML file
     * @return T - The value of the parameter in the correct data type
     */
    template< typename T >
    T convert_parameter_from_string_to_type( const std::string& aString );

    template<>
    bool convert_parameter_from_string_to_type< bool >( const std::string& aString );

    /**
     * @brief create_and_set_parameter_list - Calls the create_parameter_list function and sets the parameter list with the values from the XML file in the correct data type
     * @param aModule - Module in Parameter_List_Type enum type
     * @param aChild - The index of the sub-module
     * @param aSubChild - The index of the sub-module type for special forms like "GEN/Geometry", "OPT/Algorithm" and "SOL/Linear_Algorithm", if not these forms, then 0
     * @param tKeys - The keys of the XML file parameter list
     * @param tValues - The values of the XML file parameter list
     * @return Parameter_List - The parameter list with the set values from the XML file
     */

    Parameter_List create_and_set_parameter_list( Module_Type aModule,
            uint                                                      aChild,
            uint                                                      aSubChild,
            const Vector< std::string >&                              aKeys,
            const Vector< std::string >&                              aValues );

    /**
     * @brief Create a parameter list for a given module, child, and sub-child
     
     * @param aModule module to create the parameter list for
     * @param aChild child to create the parameter list for
     * @param aSubChild sub-child to create the parameter list for
     * @return Parameter_List
     */

    Parameter_List create_parameter_list( Module_Type aModule, uint aChild, uint aSubChild );

}    // namespace moris

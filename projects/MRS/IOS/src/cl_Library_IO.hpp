/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Library_IO.hpp
 *
 */

#ifndef MORIS_CL_LIBRARY_IO_HPP
#define MORIS_CL_LIBRARY_IO_HPP

#include <string>
#include <set>
#include "dlfcn.h"
#include "assert.hpp"

#include "cl_Library_Enums.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_Vector.hpp"

namespace moris
{
    // -----------------------------------------------------------------------------

    // Forward declare the XML-parser
    class XML_Parser;

    // Define what a module parameter list is
    typedef Vector< Vector< moris::Parameter_List > > ModuleParameterList;

    // Define what a parameter function is
    typedef void ( *Parameter_Function )( ModuleParameterList& aParameterList );

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

        // Library type
        Library_Type mLibraryType;

        // storage for the parameters for the various Modules
        Vector< ModuleParameterList > mParameterLists;

        // XML parser for output
        std::unique_ptr< XML_Parser > mXmlWriter;

        // list of parameter lists supported by the particular workflow
        std::set< Parameter_List_Type > mSupportedParamListTypes;

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
        convert_to_absolute_file_path( const std::string aFilePath );

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
                ModuleParameterList& aParamListToModify,
                ModuleParameterList& aParamsToAdd );

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

        /**
         * @brief give a parameter list file to the library and read it
         *
         * @param aFileName Name/path to the file to be used
         * @param aFileType type of the file to be used
         */
        virtual void
        load_parameter_list( std::string aFileName, File_Type aFileType );

        // -----------------------------------------------------------------------------

        /**
         * @brief finishes the initialization of the library and locks it from modification
         */
        virtual void
        finalize()
        {
            MORIS_ERROR( false, "Library_IO::finalize() - Function not implemented in this base class." );
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief fills the member parameter lists with the standard parameters for all modules
         */
        virtual void
        load_all_standard_parameters()
        {
            MORIS_ERROR( false, "Library_IO::load_all_standard_parameters() - Function not implemented in this base class." );
        }

        // -----------------------------------------------------------------------------

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
         * @brief print a summary of all parameters after finalizing the library
         *
         * @param aOutputFileName name of the xml file to print the parameters to
         */
        void
        print_parameter_receipt( const std::string aOutputFileName );

        // -----------------------------------------------------------------------------

        void
        write_module_parameter_list_to_xml_tree( const Parameter_List_Type aModule );

        // -----------------------------------------------------------------------------

        void
        write_parameter_list_to_xml_buffer( Parameter_List& aParameterList );

        // -----------------------------------------------------------------------------

        std::string
        get_sub_parameter_list_location_in_xml_tree(
                const Parameter_List_Type aModule,
                const uint                aSubParamListIndex = MORIS_UINT_MAX,
                const bool                aIsInnerParamList  = false );

        // -----------------------------------------------------------------------------

        /**
         * @brief Get the parameters for module object
         *
         * @param aParamListType
         * @return ModuleParameterList
         */
        ModuleParameterList
        get_parameters_for_module( Parameter_List_Type aParamListType ) const;

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
        load_function( std::string aFunctionName, bool aThrowError = true )
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

        // -----------------------------------------------------------------------------

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
                std::string const & aEnding )
        {
            return aString.substr( aString.length() - aEnding.length() ) == aEnding;
        }

        // -----------------------------------------------------------------------------

    };    // class Library_IO

    // -----------------------------------------------------------------------------

}    // namespace moris

#endif    // MORIS_CL_LIBRARY_IO_HPP

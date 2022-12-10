/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XML_Parser.hpp
 *
 */

#ifndef PROJECTS_MRS_IOS_SRC_CL_XML_PARSER_HPP_
#define PROJECTS_MRS_IOS_SRC_CL_XML_PARSER_HPP_

#include <string>
#include <set>
#include <exception>
#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Library_Enums.hpp"

namespace moris
{
    // -----------------------------------------------------------------------------

    class XML_Parser
    {
        // flag whether this has been initialized
        XML_Mode mMode;

        // path to the xml file
        std::string mFilePath;

        // tree object containing the xml information
        boost::property_tree::ptree mTree;

        // subtree that can be handled/modified independent of the main tree
        boost::property_tree::ptree mBuffer;

        // -----------------------------------------------------------------------------

      public:
        // -----------------------------------------------------------------------------

        /**
         * opens a path and writes the content into mTree
         */
        XML_Parser()
            : mMode( XML_Mode::UNINITIALIZED )
            , mFilePath( "" )
            , mTree()
            , mBuffer()
        {
            // nothing else
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief constructor for both creating and initializing the XML_Parser at the same time
         * 
         * @param aFilePath file path to open
         * @param aMode whether the XML parser is initialized as a reader or writer
         */
        XML_Parser( const std::string& aFilePath, XML_Mode aMode = XML_Mode::READ );

        // -----------------------------------------------------------------------------

        void
        initialize_read( const std::string& aFilePath );

        // -----------------------------------------------------------------------------

        void
        initialize_write( const std::string& aFilePath );

        // -----------------------------------------------------------------------------

        /**
         * trivial destructor
         */
        ~XML_Parser()
        {
        }

        // -----------------------------------------------------------------------------

        /**
         * gets a key from the tree and writes the value into aValue
         */
        template< typename T >
        void
        get( const std::string aKey, T& aValue ) const
        {
            aValue = mTree.get< T >( aKey );
        }

        // -----------------------------------------------------------------------------

        /**
         * sets a value in the tree
         */
        template< typename T >
        void
        set( const std::string aKey, const T& aValue )
        {
            mTree.put( aKey, aValue );
        }

        // -----------------------------------------------------------------------------

        /**
         * sets a value in the tree
         */
        template< typename T >
        void
        set_in_buffer( const std::string aKey, const T& aValue )
        {
            mBuffer.put( aKey, aValue );
        }

        // -----------------------------------------------------------------------------

        /**
         * sets a value in the tree
         */
        template< typename T >
        void
        set_attribute_in_buffer(
                const std::string& aAttributeName,
                const T&           aAttributeValue,
                const std::string  aLocation = "" )
        {
            std::string tKey = "<xmlattr>." + aAttributeName;
            if( aLocation != "" )
            { 
                tKey = aLocation + "." + tKey; 
            }
            mBuffer.put( tKey, aAttributeValue );
        }

        // -----------------------------------------------------------------------------

        void
        flush_buffer_to_tree( const std::string& aLocation );

        // -----------------------------------------------------------------------------

        /**
         * writes the data into a file
         */
        void save( const std::string& aFilePath )
        {
            boost::property_tree::xml_writer_settings< std::string > tWriteSettings( '\t', 1 );
            boost::property_tree::write_xml( aFilePath, mTree, std::locale(), tWriteSettings );
        }

        // -----------------------------------------------------------------------------

        void
        save()
        {
            this->save( mFilePath );
        }

        // -----------------------------------------------------------------------------

        /**
         * counts number of entries in a subtree
         */
        moris::size_t
        count_keys_in_subtree(
                const std::string& aSubTree,
                const std::string& aLabel );
        // -----------------------------------------------------------------------------

        /**
         * returns entries from a subtree, assuming that the subtree is flat
         */
        void
        get_keys_from_subtree(
                const std::string&   aSubTree,
                const std::string&   aLabel,
                const moris::size_t& aIndex,
                Cell< std::string >& aFirst,
                Cell< std::string >& aSecond );

        // -----------------------------------------------------------------------------

        bool
        to_bool( std::string const & aStr );

        // -----------------------------------------------------------------------------

    };    // class XML_Parser

}    // namespace moris

#endif /* PROJECTS_MRS_IOS_SRC_CL_XML_PARSER_HPP_ */

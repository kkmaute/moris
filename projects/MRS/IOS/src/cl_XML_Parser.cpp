/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XML_Parser.cpp
 *
 */

#include "cl_XML_Parser.hpp"

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
#include "fn_assert.hpp"

namespace moris
{

    // -----------------------------------------------------------------------------

    XML_Parser::XML_Parser( const std::string& aFilePath, XML_Mode aMode )
    {
        // initialize the parser
        if( aMode == XML_Mode::READ )
        {
            this->initialize_read( aFilePath );
        }
        else if( aMode == XML_Mode::WRITE )
        {
            this->initialize_write( aFilePath );
        }
        else
        {
            MORIS_ERROR( false, "XML_Parser::XML_Parser() - Initializing an XML parser with in unknown mode." );
        }
    }

    // -----------------------------------------------------------------------------

    void
    XML_Parser::initialize_read( const std::string& aFilePath )
    {
        // check that this parser has not been initialized before
        if( mMode != XML_Mode::UNINITIALIZED )
        {
            std::cerr << "XML_Parser::initialize_read() - " <<
                    "Trying to initialize an XML-parser that has already been initialized." << aFilePath << std::endl;
            throw;
        }

        // test if the file specified exists
        if ( boost::filesystem::exists( aFilePath ) )
        {
            mFilePath = aFilePath;
            boost::property_tree::read_xml( mFilePath, mTree );
        }
        else
        {
            std::cerr << "XML_Parser::initialize_read() - " << 
                    "Something went wrong while trying to load from XML file: " << aFilePath << std::endl;
            throw;
        }

        // mark the parser as initialized
        mMode = XML_Mode::READ;
    }

    // -----------------------------------------------------------------------------

    void
    XML_Parser::initialize_write( const std::string& aFilePath )
    {
        // check that this parser has not been initialized before
        if( mMode != XML_Mode::UNINITIALIZED )
        {
            std::cerr << "XML_Parser::initialize_write() - " <<
                    "Trying to initialize an XML-parser that has already been initialized." << aFilePath << std::endl;
            throw;
        }

        // test if the file specified exists
        if ( boost::filesystem::exists( aFilePath ) )
        {
            MORIS_LOG( "An XML writer is initialized with path '%s'." , aFilePath.c_str() ); 
            MORIS_LOG( "This file already exists and will be overwritten." ); 
        }
        
        // set the file path
        mFilePath = aFilePath;

        // mark the parser as initialized
        mMode = XML_Mode::WRITE;
    }

    // -----------------------------------------------------------------------------

    void
    XML_Parser::flush_buffer_to_tree( const std::string& aLocation )
    {
        // check that the xml-parser is in writing mode
        MORIS_ASSERT( mMode == XML_Mode::WRITE, "XML_Parser::flush_buffer_to_tree() - xml-parser must be in write mode to use this function." );

        // add the subtree stored in the buffer to the desired location in the main tree
        mTree.add_child( aLocation, mBuffer );

        // clear the buffer
        mBuffer.clear();
    }

    // -----------------------------------------------------------------------------

    /**
     * counts number of entries in a subtree
     */
    moris::size_t
    XML_Parser::count_keys_in_subtree( 
            const std::string& aSubTree,
            const std::string& aLabel )
    {
        // initialize counter
        moris::size_t aCount = 0;

        // loop over all entries in this tag
        BOOST_FOREACH ( boost::property_tree::ptree::value_type &v, mTree.get_child( aSubTree ) )
        {
            if ( v.first.data() == aLabel )
            {
                // increment counter
                ++aCount;
            }
        }

        // return counter
        return aCount;
    }

    // -----------------------------------------------------------------------------

    /**
     * returns entries from a subtree, assuming that the subtree is flat
     */
    void
    XML_Parser::get_keys_from_subtree( 
            const std::string&   aSubTree,
            const std::string&   aLabel,
            const moris::size_t& aIndex,
            Cell< std::string >& aFirst,
            Cell< std::string >& aSecond )
    {
        // tidy up output data
        aFirst.clear();
        aSecond.clear();

        // initialize counter
        moris::size_t tCount = 0;

        // loop over all entries in this tag
        BOOST_FOREACH ( boost::property_tree::ptree::value_type &v,
                mTree.get_child( aSubTree ) )
        {
            if ( v.first.data() == aLabel )
            {
                if ( tCount == aIndex )
                {
                    if ( !v.second.empty() )
                    {
                        BOOST_FOREACH ( boost::property_tree::ptree::value_type &w,
                                v.second )
                        {
                            if ( w.second.empty() )
                            {
                                aFirst.push_back( w.first.data() );
                                aSecond.push_back( w.second.data() );
                            }
                        }
                    }
                    break;
                }
                ++tCount;
            }
        }
    }

    // -----------------------------------------------------------------------------

    bool
    XML_Parser::to_bool( std::string const &aStr )
    {
        if ( aStr == "true" || aStr == "1" )
        {
            return true;
        }
        else if ( aStr == "false" || aStr == "0" )
        {
            return false;
        }
        else
        {
            MORIS_ERROR( false, 
                    "XML_Parser::to_bool() - Unrecognized string passed into to_bool. "
                    "Needs to be true or 1 for bool = true or false or 0 for bool = false, "
                    "be sure to check for extraneous spaces" );
            return false;
        };
    }

    // -----------------------------------------------------------------------------

}    // namespace moris

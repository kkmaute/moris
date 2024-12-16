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

#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
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
            std::cerr << "XML_Parser::initialize_read() - " << "Trying to initialize an XML-parser that has already been initialized." << aFilePath << '\n';
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
            std::cerr << "XML_Parser::initialize_read() - " << "Something went wrong while trying to load from XML file: " << aFilePath << '\n';
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
            std::cerr << "XML_Parser::initialize_write() - " << "Trying to initialize an XML-parser that has already been initialized." << aFilePath << '\n';
            throw;
        }

        // test if the file specified exists
        MORIS_LOG( "An XML writer is initialized with path '%s'." , aFilePath.c_str() ); 
        if ( boost::filesystem::exists( aFilePath ) )
        {
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
            Vector< std::string >& aFirst,
            Vector< std::string >& aSecond )
    {
        // tidy up output data
        aFirst.clear();
        aSecond.clear();

        // initialize counter
        moris::size_t tCount = 0;

        // loop over all entries in this tag
        BOOST_FOREACH ( boost::property_tree::ptree::value_type &v, mTree.get_child( aSubTree ) )
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
                                // replace all $ in aFirst with whitespace
                                std::string tFirst = w.first.data();
                                std::replace( tFirst.begin(), tFirst.end(), '$', ' ' );
                                aFirst.push_back( tFirst );
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

    // A function that checks if a label exists in a subtree
    bool XML_Parser::label_exists_in_subtree( const std::string& aSubTree, const std::string& aLabel )
    {
        // loop over all entries in this tag
        BOOST_FOREACH ( boost::property_tree::ptree::value_type &v, mTree.get_child( aSubTree ) )
        {
            if ( v.first.data() == aLabel )
            {
                return true;
            }
        }

        return false;
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

    bool
    XML_Parser::copy_subtree_into_buffer(
            const std::string&  aNodeNestedUnder,
            const std::string&  aSubTreeLabel,
            const moris::size_t aIndex )
    {
        // check that the buffer is a read buffer
        MORIS_ASSERT( mMode == XML_Mode::READ, "XML_Parser::copy_subtree_into_buffer() - "
                "This operation should only be performed in READ mode. XML-parser is currently in mode #%i.",
                (uint)( mMode ) );

        // initialize counter to keep 
        moris::size_t tCount = 0;

        // loop over all entries in this tag
        BOOST_FOREACH ( boost::property_tree::ptree::value_type & iNode, mTree.get_child( aNodeNestedUnder ) )
        {
            // if this is the instance of type of subtree this function is looking for
            if ( iNode.first.data() == aSubTreeLabel )
            {
                // if this is the particular instance we're looking for, copy it into the buffer
                if( aIndex == tCount )
                {
                    // store subtree in buffer
                    mBuffer = iNode.second;
                    
                    // stop and return that the subtree has been buffered successfully
                    return true;
                }

                // otherwise skip this instance and look for the next one
                tCount++;
            }
        }

        // return that the subtree has NOT been buffered successfully
        return false;
    }

    // -----------------------------------------------------------------------------

}    // namespace moris

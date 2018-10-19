/*
 * cl_XML_Parser.hpp
 *
 *  Created on: Sep 10, 2018
 *      Author: messe
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

namespace moris
{
// -----------------------------------------------------------------------------

    class XML_Parser
    {
        //! path to the xml file
        std::string                 mFilePath;

        //! tree object containing the xml information
        boost::property_tree::ptree mTree;

// -----------------------------------------------------------------------------
    public:
// -----------------------------------------------------------------------------

        /**
         * opens a path and writes the content into mTree
         */
        XML_Parser( const std::string & aFilePath ) : mFilePath( aFilePath )
        {
            // test if file exists
            if ( boost::filesystem::exists( aFilePath ) )
            {
                boost::property_tree::read_xml( aFilePath, mTree );
            }
            else
            {
                std::cerr << "Something went wrong while trying to load from XML file " <<
                                  aFilePath << "." << std::endl;
                exit(-1);
            }
        }

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
        template < typename T >
        void
        get( const std::string aKey, T & aValue ) const
        {
            aValue = mTree.get< T >( aKey );
        }

// -----------------------------------------------------------------------------

        /**
         * sets a value in the tree
         */
        template < typename T >
        void
        set( const std::string aKey, const T & aValue )
        {
            mTree.put( aKey, aValue );
        }

// -----------------------------------------------------------------------------

        /**
         * writes the data into a file
         */
        void
        save( const std::string & aFilePath )
        {
            boost::property_tree::write_xml( aFilePath, mTree );
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
                const std::string & aSubTree,
                const std::string & aLabel )
        {
            // initialize counter
            moris::size_t aCount = 0;

            // loop over all entries in this tag
            BOOST_FOREACH( boost::property_tree::ptree::value_type &v,
                    mTree.get_child( aSubTree ) )
            {
                if( v.first.data() == aLabel )
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
        get_keys_from_subtree(
                const std::string   & aSubTree,
                const std::string   & aLabel,
                const moris::size_t & aIndex,
                Cell< std::string > & aFirst,
                Cell< std::string > & aSecond )
        {
            // tidy up output data
            aFirst.clear();
            aSecond.clear();

            // initialize counter
            moris::size_t tCount = 0;

            // loop over all entries in this tag
            BOOST_FOREACH( boost::property_tree::ptree::value_type &v,
                    mTree.get_child( aSubTree ) )
            {
                if( v.first.data() == aLabel )
                {
                    if ( tCount == aIndex )
                    {

                        if( ! v.second.empty() )
                        {
                            BOOST_FOREACH( boost::property_tree::ptree::value_type &w,
                                    v.second )
                            {
                                if( w.second.empty() )
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
    };
}

#endif /* PROJECTS_MRS_IOS_SRC_CL_XML_PARSER_HPP_ */

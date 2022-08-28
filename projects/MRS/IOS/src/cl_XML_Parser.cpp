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
/*
 * cl_XML_Parser.cpp
 *
 *  Created on: Sep 10, 2018
 *      Author: messe
 */

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

        /**
         * counts number of entries in a subtree
         */
        moris::size_t XML_Parser::count_keys_in_subtree( const std::string & aSubTree,
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
        void XML_Parser::get_keys_from_subtree( const std::string         & aSubTree,
                                                const std::string         & aLabel,
                                                const moris::size_t       & aIndex,
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

        bool XML_Parser::to_bool(std::string const & aStr)
        {
            if(aStr == "true" || aStr == "1")
            {
                return true;
            }
            else if(aStr == "false" || aStr == "0")
            {
                 return false;
            }
            else
            {
                MORIS_ERROR(0,"Unrecognized string passed into to_bool. Needs to be true or 1 for bool = true or false or 0 for bool = false, be sure to check for extraneous spaces");
                return false;
            };
        }

}


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
                std::cerr << "Something went wrong while trying to load from " <<
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
    };
}

#endif /* PROJECTS_MRS_IOS_SRC_CL_XML_PARSER_HPP_ */

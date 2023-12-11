//
// Created by frank on 12/7/23.
//

#ifndef MORIS_CL_JSON_OBJECT_HPP
#define MORIS_CL_JSON_OBJECT_HPP

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "cl_Matrix.hpp"
#include "typedefs.hpp"

namespace moris
{
    /**
     * @brief The Json object is a simple wrapper around boost::property_tree::ptree. This prevents boost from being
     * included in moris classes.
     */
    using Json = boost::property_tree::ptree;

    void write_json( std::string const &aFileName, Json const &aTree );

    // to_json implementation for all integral types
    template< typename T,
            std::enable_if_t< std::is_arithmetic< T >::value, bool > = true >
    Json to_json( T const &aValue )
    {
        Json tElement;
        tElement.put_value( aValue );
        return tElement;
    }

    template< typename T >
    Json to_json( moris::Cell< T > const &aVector )
    {
        Json tList;
        for ( auto const &tValue : aVector )
        {
            tList.push_back( { "", to_json( tValue ) } );
        }
        return tList;
    }

    Json to_json( std::string const &aString );

    Json to_json( Json const &aTree );

    Json to_json( Matrix< DDRMat > const &aMatrix );

    std::string to_string( Json const &aTree );
}    // namespace moris


#endif    // MORIS_CL_JSON_OBJECT_HPP

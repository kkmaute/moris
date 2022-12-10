/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Param_List.hpp
 *
 */

#ifndef MORIS_CONTAINERS_CL_PARAM_LIST_HPP_
#define MORIS_CONTAINERS_CL_PARAM_LIST_HPP_

// c++ header files.
#include <cstring>
#include <map>
#include <string>
#include <sstream>

// Third-party header files.
#include <boost/variant.hpp>

// MORIS library header files.
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "assert.hpp"
#include "core.hpp"
#include "ios.hpp"
#include "fn_Parsing_Tools.hpp"

namespace moris
{
    namespace containers
    {
        struct strcmp
        {
            bool
            operator()( const std::string& a, const std::string& b ) const
            {
                return std::strcmp( a.c_str(), b.c_str() ) < 0;
            }
        };

    }    // namespace containers
}    // namespace moris

namespace moris
{

    /**
     * The parameter list class.
     *
     * A parameter list can be created as follows:
     * @include CON/src/cl_param_list/cl_param_list_insert.inc
     */
    template< typename Variant >
    class Param_List
    {
        //------------------------------------------------------------------------------
      private:

        std::map< std::string, Variant, moris::containers::strcmp > mParamMap;

        //------------------------------------------------------------------------------

      public:

        //------------------------------------------------------------------------------

        /**
         * Constructor
         */
        Param_List() = default;

        //------------------------------------------------------------------------------

        /**
         * Destructor
         */
        ~Param_List() = default;

        //------------------------------------------------------------------------------

        /**
         * @brief Extends an existing container by the element inserted.
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed; should not contain trailing or leading whitespaces
         * @param[in] aVal Value corresponding to aKey; if value is a string,
         *            pair of strings, or character array leading and trailing whitespaces
         *            are trimmed off
         */
        void
        insert( const std::string& aKey, Variant aVal )
        {
            // check for leading and trailing whitespaces in key
            std::string tKey = aKey;
            split_trim_string( tKey, "" );
            MORIS_ERROR( !aKey.compare( tKey ),
                    "Param_List::insert - key contains whitespaces" );

            // if value is a string remove leading and trailing whitespaces
            if ( aVal.type() == typeid( std::string ) )
            {
                // extract string from variant variable
                std::string tVal( boost::get< std::string >( aVal ) );

                // trim off leading and trailing whitespaces
                split_trim_string( tVal, ",;" );

                // overwrite input parameter
                aVal = tVal;
            }

            // check for whitespaces in value if value is a pair of strings
            if ( aVal.type() == typeid( std::pair< std::string, std::string > ) )
            {
                // extract pair of strings from variant variable
                std::pair< std::string, std::string > tVal( boost::get< std::pair< std::string, std::string > >( aVal ) );

                // extract elements of pair
                std::string tFirst  = tVal.first;
                std::string tSecond = tVal.second;

                // trim off leading and trailing whitespaces
                split_trim_string( tFirst, ",;" );
                split_trim_string( tSecond, ",;" );

                // overwrite input parameter
                aVal = std::make_pair( tFirst, tSecond );
            }

            // if variant type is const char* convert it into a string
            // and check for whitespaces
            if ( aVal.type() == typeid( const char* ) )
            {
                // extract string from variant variable
                std::string tVal( boost::get< const char* >( aVal ) );

                // trim off leading and trailing whitespaces
                split_trim_string( tVal, ",;" );

                // overwrite input parameter
                aVal = tVal;
            }
            mParamMap.insert( { aKey, aVal } );
        }

        //------------------------------------------------------------------------------

        /**
         * Removes the specified key and its associated value from the map
         *
         * @param aKey the key to be erased
         */
        void
        erase( const std::string& aKey )
        {
            mParamMap.erase( aKey );
        }

        //------------------------------------------------------------------------------

        /**
         * @brief Sets an element to a value if it exists, otherwise an error is thrown
         *        whitespaces in key and value string will be removed
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed; spaces will be removed
         * @param[in] aVal Value corresponding to aKey; if values is a string, pair of
         *            strings or character array, whitespaces will be removed
         */
        void
        set( const std::string& aKey, Variant aVal )
        {
            // create copy of key string such that it can be manipulated
            std::string tKey = aKey;

            // remove spurious whitespaces from key string
            split_trim_string( tKey, "" );

            // find key in map
            auto it = mParamMap.find( tKey );

            // if key does not exist in map
            if ( it == mParamMap.end() )
            {
                // create error message
                std::string tError = "The requested parameter '" + tKey + "' can not be set because it does not exist.\n";

                // throw error
                MORIS_ERROR( false, tError.c_str() );
            }
            else
            {
                // if variant type is string remove whitespaces
                if ( aVal.type() == typeid( std::string ) )
                {
                    std::string tVal( boost::get< std::string >( aVal ) );

                    // remove whitespaces from string
                    split_trim_string( tVal, ",;" );

                    // overwrite input argument
                    aVal = tVal;
                }

                // if variant type is pair of strings remove whitespaces
                if ( aVal.type() == typeid( std::pair< std::string, std::string > ) )
                {
                    std::pair< std::string, std::string > tVal( boost::get< std::pair< std::string, std::string > >( aVal ) );

                    // extract elements of pair
                    std::string tFirst  = tVal.first;
                    std::string tSecond = tVal.second;

                    // remove whitespaces from strings
                    split_trim_string( tFirst, ",;" );
                    split_trim_string( tSecond, ",;" );

                    // overwrite input argument with new pair
                    aVal = make_pair( tFirst, tSecond );
                }

                // if variant type is const char* convert it into a string and remove whitespaces
                if ( aVal.type() == typeid( const char* ) )
                {
                    std::string tVal( boost::get< const char* >( aVal ) );

                    // remove whitespaces from string
                    split_trim_string( tVal, ",;" );

                    // overwrite input argument string
                    aVal = tVal;
                }
                it->second = aVal;
            }
        }

        //------------------------------------------------------------------------------

        /**
         * @brief Sets an element to a value if it exists, if it doesn't it is added to the map
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed; spaces will be removed
         * @param[in] aVal Value corresponding to aKey; if values is a string, pair of
         *            strings or character array, whitespaces will be removed
         * @return true / false whether the key already existed
         */
        bool
        set_or_insert( const std::string& aKey, Variant aVal )
        {
            // create copy of key string such that it can be manipulated
            std::string tKey = aKey;

            // remove spurious whitespaces from key string
            split_trim_string( tKey, "" );

            // see if the key is already in the map
            auto it = mParamMap.find( tKey );

            // if key does not exist in map
            if ( it == mParamMap.end() )
            {
                // add it to the map
                this->insert( tKey, aVal );

                // return that the key has NOT already existed
                return false;
            }
            // if the key already exists
            else
            {
                // set the key to its new value
                this->set( tKey, aVal );

                // return that the key has already existed
                return true;
            }
        }

        //------------------------------------------------------------------------------

        /**
         * @brief Access operator.
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         *
         * The mapped values can be accessed directly using this
         * operator. If the input parameter matches the key of an
         * element in the container, the function returns a reference to
         * the mapped value else an error is encountered.
         *
         * @include CON/src/cl_param_list/cl_param_list_access.inc
         */
        Variant&
        operator()( const std::string& aKey )
        {
            auto it = mParamMap.find( aKey );

            // check if parameter exists
            if ( it == mParamMap.end() )
            {
                // create error message
                std::string tError = "The requested parameter '" + aKey + "' does not exist.\n";

                // throw error
                MORIS_ERROR( false, tError.c_str() );
            }

            return it->second;
        }

        //------------------------------------------------------------------------------

        /**
         * @brief Determine the underlying type of the entry.
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         *
         * @return The index of the type as ordered in the variant entry
         *         of mParamMap
         */
        moris::sint
        which( const std::string& aKey )
        {
            auto it = mParamMap.find( aKey );

            // check if parameter exists
            if ( it == mParamMap.end() )
            {
                // create error message
                std::string tError = "The requested parameter '" + aKey + "' does not exist.\n";

                // throw error
                MORIS_ERROR( false, tError.c_str() );
            }

            return ( it->second ).which();
        }

        //------------------------------------------------------------------------------

        /**
         * @brief Get the value corresponding to a given key.
         *
         * @param[in] aKey Key corresponding to the mapped value that
         *            needs to be accessed
         *
         * @return The value corresponding to aKey.
         */
        template< typename Key >
        const Key&
        get( const std::string& aKey )
        {
            auto it = mParamMap.find( aKey );

            // check if parameter exists
            if ( it == mParamMap.end() )
            {
                // create error message
                std::string tError = "The requested parameter '" + aKey + "' does not exist.\n";

                // throw error
                MORIS_ERROR( false, tError.c_str() );
            }

            // check if the requested type is correct
            if ( boost::get< Key >( &( it->second ) ) == nullptr )
            {
                // create error message
                std::string tError = "The parameter '" + aKey + "' was requested with an incorrect type.\n";
                MORIS_ERROR( false, tError.c_str() );
            }

            return boost::get< Key >( it->second );
        }

        //------------------------------------------------------------------------------

        /**
         * @brief Get the beginning iterator of the underlying map.
         *
         * @return An iterator pointing to the beginning of mParamMap.
         */
        auto
        begin() -> decltype( mParamMap.begin() )
        {
            return mParamMap.begin();
        }

        //------------------------------------------------------------------------------

        /**
         * @brief Get the end iterator of the underlying map.
         *
         * @return An iterator pointing to the end of mParamMap.
         */
        auto
        end() -> decltype( mParamMap.end() )
        {
            return mParamMap.end();
        }

        //------------------------------------------------------------------------------

        /**
         * @brief check if parameter list is empty
         * 
         * @return true if parameter list doesn't have any entries
         */
        bool 
        isempty()
        {
            return mParamMap.empty();
        }

        //------------------------------------------------------------------------------

        /**
         * @brief get the size of the parameter list
         * 
         * @return size_t number of entries in the parameter list
         */
        size_t 
        size()
        {
            return mParamMap.size();
        }

        //------------------------------------------------------------------------------

    }; // class ParameterList

    //------------------------------------------------------------------------------

    template< typename Variant >
    std::string
    convert_param_value_to_string( Variant aParameterValue )
    {
        std::stringstream tStringStream;

        if ( boost::get< bool >( & aParameterValue ) != nullptr )
        {
            tStringStream << boost::get< bool >( aParameterValue );
        }

        else if ( boost::get< sint >( & aParameterValue ) != nullptr )
        {
            tStringStream << boost::get< sint >( aParameterValue );
        }

        else if ( boost::get< real >( & aParameterValue ) != nullptr )
        {
            tStringStream << boost::get< real >( aParameterValue );
        }

        else if ( boost::get< std::string >( & aParameterValue ) != nullptr )
        {
            tStringStream << boost::get< std::string >( aParameterValue );
        }

        else if ( boost::get< const char* >( & aParameterValue ) != nullptr )
        {
            tStringStream << boost::get< const char* >( aParameterValue );
        }

        else if ( boost::get< uint >( & aParameterValue ) != nullptr )
        {
            tStringStream << boost::get< uint >( aParameterValue );
        }

        else
        {
            tStringStream << "<print_error>" << std::endl;
        }

        return tStringStream.str();
    }

    //------------------------------------------------------------------------------

    // datatype for parameter lists
    typedef boost::variant< bool, sint, real, const char*, std::string, uint, std::pair< std::string, std::string > > ParameterListTypes;
    typedef Param_List< ParameterListTypes > ParameterList;

    //------------------------------------------------------------------------------

}    // namespace moris

#endif /* MORIS_CONTAINERS_CL_PARAM_LIST_HPP_ */

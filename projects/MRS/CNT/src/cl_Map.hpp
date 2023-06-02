/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Map.hpp
 *
 */

#ifndef MORIS_CONTAINERS_CL_MAP_HPP_
#define MORIS_CONTAINERS_CL_MAP_HPP_

// C++ header files.
#include <exception>
#include <map>
#include <iostream>

// MORIS library header files.
#include "typedefs.hpp" // COR/src
#include "assert.hpp"

namespace moris
{
    template< typename T1, typename T2 >
    class map
    {
    private:

        /**
         * MORIS map
         */
        std::map< T1, T2 > mMap;

    public:

    /**
     * empty container constructor aMap. Map provides a way to organize containers by keyword identifiers.
     * Each element associates a key to a mapped value: Keys are meant to identify the elements whose main
     * content is the mapped value.
     */
        explicit map() = default;

        /**
         * container constructor: creates a container copy of each of the elements in aMap.
         */
        map( moris::map< T1, T2 > const & aMap ) = default;

        /*
         * Default destructor
         */
        ~map() = default;

        /**
         * The function insert extends an existing container by the number of elements inserted. During this operation,
         * a check is performed to ensure all key values are distinct. If the inserted key value is non-unique, the element
         * is not inserted. An alternative way to insert elements is to use operator[].
         *
         * @return Non-const reference to underlying standard pair
         */
        template<class pair>
        std::pair< T2, bool >
        insert( const pair & apair )
        {
            auto tIterator = mMap.insert(apair);
            return std::make_pair( tIterator->first.second, tIterator->second );
        }

        /**
         * Access operator. The mapped values can be accessed directly or defined using this operator. if the input parameter aK
         * matches the key of an element in the container, the function returns a reference to the mapped value. If the input
         * parameter aK does not match a key of any element in the container, the function inserts a new element with that key and
         * returns a reference to its mapped value. Notice that this always increases the container size by one, even if no mapped
         * value is assigned to the element.
         *
         * @param[in] aK Key value of the element whose mapped value is accessed
         *
         * @include CON/src/cl_Map/cl_Map.inc
         */
        T2 &
        operator[]( const T1 & aK )
        {
            return mMap[aK];
        }

        /**
         * Checks if key exists in map
         *
         * @param[in] aK Key to be searched for
         */
        bool
        key_exists( const T1 & aK ) const
        {
            auto tIterator = mMap.find(aK);

            bool tReturn = true;

            if ( tIterator == mMap.end() )
            {
                tReturn = false;
            }

            return tReturn;
        }

        /**
         * Searches container for an element with a key equivalent to the input parameter aK. If the key is found,
         * it returns the value of the element. If the key is not found, this function currently throws a logic error.
         * This first entry corresponds to constant maps
         *
         * @param[in] aK Key to be searched for
         */
        const T2 &
        find( const T1 & aK ) const
        {
            auto tIterator = mMap.find(aK);

            MORIS_ERROR( tIterator != mMap.end(), "moris::map.find - Key not found" );

            return tIterator->second;
        }
        /**
         * This overloaded find functionality calls the first and recasts the constant element to a non-constant type.
         *
         * @param[in] aK a Key to be searched for
         *
         * @return
         */
        T2 &
        find( const T1 & aK )
        {
            return const_cast< T2 & >( static_cast< const moris::map< T1, T2 > * >( this )->find( aK ) );
        }

        /**
         * Checks if map is empty. This function returns true if the container size is 0, otherwise it returns false.
         *
         * @include CON/src/cl_Map/cl_Map_empty.inc
         */
        bool
        empty()
        {
            return mMap.empty();
        }

        /**
         * Returns size of map.
         */

        /**
         * Removes all elements from the map container, leaving the container with a size of 0. It is important to note
         * that the map type (i.e. double, int, etc) remains unchanged.
         */
        void
        clear()
        {
            mMap.clear();
        }

        /**
         * Returns size of map which is an unsigned integral type.
         *
         * @include CON/src/cl_Map/cl_Map_clear.inc
         */
        moris::uint
        size() const
        {
            return mMap.size();
        }

        /**
         * Returns number of elements that match the incoming key.
         *
         * @param[in] aK a Key to be counted in the map.
         */
        moris::uint
        count( const T1 & aK )
        {
            return mMap.count( aK );
        }

        /**
         * @brief Returns an iterator to the first element
         */
        auto
        begin()
            -> decltype( mMap.begin() )
        {
            return mMap.begin();
        }

        /**
         * @brief  Returns an iterator pointing to the past-the-end element.
         */
        auto
        end()
            -> decltype( mMap.end() )
        {
            return mMap.end();
        }

        /**
         * @brief  Removes an element from the map
         */
        void
        erase( const T1 & aK )
        {
            auto tIterator = mMap.find(aK);

            mMap.erase( tIterator );
        }

        /**
         * @brief Print moris map on screen
         */
        void print( const std::string & aVarName = std::string() )
        {
            std::cout << "\n-------------------------------------------------\n\n";

            std::string tVarName;

            if ( aVarName.empty() )
            {
                tVarName = aVarName;
            }
            else
            {
                tVarName = "morisMap";
            }

            for (const auto &p : mMap)
            {
                std::cout << tVarName << "[" << p.first << "] = " << p.second << '\n';
            }
        }

        //--------------------------------------------------------------------------------

        /**
         * @brief returns the underlying std::map such that other std libraries can called
         *
         * @return std::map< T1, T2 >&  underlying std::map
         */
        std::map< T1, T2 > &
        data()
        {
            return mMap;
        }

        //--------------------------------------------------------------------------------
    };
}

#endif /* CL_MAP_HPP_ */


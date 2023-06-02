/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Array.hpp
 *
 */

#ifndef MORIS_CONTAINERS_CL_ARRAY_HPP_
#define MORIS_CONTAINERS_CL_ARRAY_HPP_

// C++ header files.
#include <array>

// MORIS library header files.
#include "typedefs.hpp" // COR/src

namespace moris
{
    template < typename T, moris::size_t N >
    class Array
    {
    private:

        /**
         * Underlying std::array
         */
        std::array< T, N > mArray;

    public:

        /**
         * moris::Array constructor.
         */
        Array(
                const std::array< T, N > & aArray )
            : mArray( aArray )
        {
        }

        /**
         * moris::Array destructor.
         */
        ~Array() = default; // 'default' tells the compiler to automatically delete the underlying array

        /**
         * @brief Returns the number of elements in the array.
         *
         * @return Number of elements in the array.
         */
        moris::size_t
        size()
        {
            return mArray.size();
        }

        /**
         * @brief Checks whether the array container is empty, i.e. whether its size is 0.
         *
         * @return bool value indicating whether the array container is empty.
         */
        bool
        empty()
        {
            return mArray.empty();
        }

        /**
         * @brief Operator for accessing an array element
         *
         * @param[in] i_index Position of an element in the array.The function .at()
         *            automatically checks whether an element is within the
         *            bounds of valid elements in the container, throwing an out_of_range
         *            exception if it is not. This is in contrast with member operator[],
         *            that does not check against bounds and hence is faster.
         */
        auto
        operator[](
                const moris::size_t i_index )
        -> decltype(
#ifdef MORIS_HAVE_DEBUG
                ( mArray.at( i_index ) )
#else
                ( mArray[ i_index ] )
#endif
        )
        {
#ifdef MORIS_HAVE_DEBUG
            return( mArray.at( i_index ) );
#else
            return( mArray[ i_index ] );
#endif
        }

        /**
         * @brief Returns the pointer to the beginning of the array.
         *
         * @return Pointer to the first element of the array.
         */
        auto
        begin()->decltype( mArray.begin() )
        {
            return mArray.begin();
        }

        /**
         * @brief Returns the reference to the beginning of the array.
         *
         * @return Direct reference to the first element of the array.
         */
        auto
        front()->decltype( mArray.front() )
        {
            return mArray.front();
        }

        /**
         * @brief Returns the pointer to the end of the array.
         *
         * @return Pointer to the last+1 element of the array.
         */
        auto
        end()->decltype( mArray.end() )
        {
            return mArray.end();
        }

        /**
         * @brief Returns the reference to the end of the array.
         *
         * @return Direct reference to the last element of the array.
         */
        auto
        back()->decltype( mArray.back() )
        {
            return mArray.back();
        }

        /**
         * @brief Sets val as the value for all the elements in the array object.
         *
         * @return Nothing.
         */
        void
        fill(
                const T val )
        {
            mArray.fill( val );
            return;
        }

    };
}

#endif /* MORIS_CONTAINERS_CL_ARRAY_HPP_ */


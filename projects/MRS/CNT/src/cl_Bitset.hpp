/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Bitset.hpp
 *
 */

#ifndef SRC_CONTAINERS_CL_BITSET_HPP_
#define SRC_CONTAINERS_CL_BITSET_HPP_

// C++ header files.
#include <bitset>
#include <string>

// MORIS library header files.
#include "moris_typedefs.hpp" //COR/src

namespace moris
{
    template < moris::size_t N >
    class Bitset
    {
    private:
        /**
         * MORIS Bitset
         */
        std::bitset<N> mBitset;

    public:

        /**
         * moris::Bitset constructor
         */
        Bitset() = default;

        // -----------------------------------------------------------------------------

        /**
         * moris::Bitset constructor
         *
         * @param[in] aBitset A Bitset
         */
        Bitset( const std::bitset< N > & aBitset )
                : mBitset( aBitset )
        {
        }

        // -----------------------------------------------------------------------------

        /**
         * Constructor using an unsigned integer
         *
         * @param[in] aUint Unsigned integer
         */
        Bitset( luint aUint )
                : mBitset( std::bitset<N>( aUint ) )
        {
        }

        // -----------------------------------------------------------------------------

        /**
         * moris::Bitset destructor
         */
        ~Bitset() = default; // 'default' tells the compiler to automatically
                             // delete the underlying Bitset

        // -----------------------------------------------------------------------------

        /**
         * @brief Returns the number of elements in the Bitset.
         *
         * @return Number of elements in the Bitset.
         */
        moris::size_t
        size() const
        {
            return mBitset.size();
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief Returns the number of bits set to true.
         *
         * @return Returns the number of bits set to true.
         */
        moris::size_t
        count() const
        {
            return mBitset.count();
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief Returns the bitset in a string
         *
         * @return Returns a string of the bitset
         */
        std::string
        to_string() const
        {
            return mBitset.to_string();
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief Returns the number of bits set to true in an uinsigned long int.
         *
         * @return Returns the number of bits set to true in an uinsigned long int.
         */
        moris::uint
        to_ulong() const
        {
            return mBitset.to_ulong();
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief Returns the number of bits set to true in an uinsigned long long int.
         *
         * @return Returns the number of bits set to true in an uinsigned long long int.
         */
        unsigned long long int
        to_ullong() const
        {
            return mBitset.to_ullong();
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief Sets the position index in the bitset to 1.
         *
         * @return Bitset.
         */
        auto
        set( moris::size_t aIndex, bool aValue = true )
        -> decltype( mBitset.set( aIndex, aValue ) )
        {
            return( mBitset.set( aIndex, aValue ) );
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief Sets the position "index" in the bitset to 0.
         *
         * @return Bitset.
         */
        auto
        reset(
                const moris::size_t index )
        -> decltype( mBitset.reset( index ) )
        {
            return( mBitset.reset( index ) );
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief Sets all bits to 0.
         *
         * @return Bitset.
         */
        auto
        reset()
            -> decltype( mBitset.reset(  ) )
        {
            return( mBitset.reset(  ) );
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief Toggles the values of bits .
         *
         * @return toggles the values of bits.
         */
        auto
        flip( const moris::size_t index )
            -> decltype( mBitset.flip( index ) )
        {
            return( mBitset.flip( index ) );
        }

        // -----------------------------------------------------------------------------

        /**
         * @brief Accesses specific bit .
         *
         * @return Specific bit .
         */
        auto
        test( const moris::size_t index ) const
            -> decltype( mBitset.test( index ) )
        {
            return( mBitset.test( index ) );
        }

        // -----------------------------------------------------------------------------

        auto
        operator&() const
            -> decltype( mBitset & mBitset )
        {
            return( mBitset & mBitset );
        }

        // -----------------------------------------------------------------------------
    };
}

#endif /* SRC_CONTAINERS_CL_BITSET_HPP_ */


/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_BoostBitset.hpp
 *
 */

#ifndef SRC_CONTAINERS_CL_BOOSTBITSET_HPP_
#define SRC_CONTAINERS_CL_BOOSTBITSET_HPP_

// C++ header files.
#include <boost/dynamic_bitset.hpp>
#include <string>
// MORIS library header files.
#include "typedefs.hpp" // COR/src

namespace moris
{

    class BoostBitset
    {
    private:

        /**
         * MORIS Bitset
         */
        boost::dynamic_bitset<> mBitset;

    public:

        /**
         * moris::Bitset constructor
         */
        BoostBitset() = default;

        /**
         * moris::Bitset constructor
         *
         * @param[in] aBitset A Bitset
         */
        BoostBitset(
                const boost::dynamic_bitset<> & aBitset )
            : mBitset( aBitset )
        {
        }

        /**
         * moris::Bitbool constructor
         *
         * @param[in] aSize Size of the Cell to be initialized
         */
        BoostBitset(
                moris::uint const aSize )
        : mBitset( aSize )
        {
        }

//        Bitset(N & aNumBits )
//            : mBitset( aNumBits )
//        {
//        }

        /**
         * moris::Bitset destructor
         */
        ~BoostBitset() = default; // 'default' tells the compiler to automatically
                             // delete the underlying Bitset

        /**
         * @brief Returns the number of elements in the Bitset.
         *
         * @return Number of elements in the Bitset.
         */
        moris::size_t
        size()
        {
            return mBitset.size();
        }

        /**
         * @brief Returns the number of elements in the Bitset. (const version)
         *
         * @return Number of elements in the Bitset.
         */
        moris::size_t
        size()
        const {
            return mBitset.size();
        }

        /**
         * @brief Returns the number of bits set to true.
         *
         * @return Returns the number of bits set to true.
         */
        moris::size_t
        count()
        {
            return mBitset.count();
        }

        /**
        * @brief Returns the number of bits set to true (const version)
        *
        * @return Returns the number of bits set to true.
        */
        moris::size_t
        count() const
        {
            return mBitset.count();
        }

        /**
         * @brief Sets the position index in the bitset to 1.
         *
         * @return Bitset.
         */
        auto
        resize(
                const moris::size_t index )
        -> decltype( mBitset.resize( index ) )
        {
            return( mBitset.resize( index ) );
        }

        /**
         * @brief Sets the position index in the bitset to 1.
         *
         * @return Bitset.
         */
//        void
//        to_string(
//                const boost::dynamic_bitset<> & mBitset,
//                 std::basic_string const & aString )
//        {
//            boost::to_string( mBitset, aString );
//        }

        /**
         * @brief Sets the position index in the bitset to 1.
         *
         * @return Bitset.
         */
        auto
        set(
                const moris::size_t index )
        -> decltype( mBitset.set( index ) )
        {
            return( mBitset.set( index ) );
        }

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

        /**
         * @brief Clear the Bitset
         *
         * @return Cleared Bitset.
         */
        auto
        clear()
        -> decltype( mBitset.clear(  ) )
        {
            return( mBitset.clear(  ) );
        }

        /**
         * @brief Toggles the values of bits .
         *
         * @return toggles the values of bits.
         */
        auto
        flip(
                const moris::size_t index )
        -> decltype( mBitset.flip( index ) )
        {
            return( mBitset.flip( index ) );
        }

        /**
         * @brief Accesses specific bit .
         *
         * @return Specific bit .
         */
        auto
        test(
                const moris::size_t index )
        -> decltype( mBitset.test( index ) )
        {
            return( mBitset.test( index ) );
        }

        /**
         * @brief Accesses specific bit . (const variant)
         *
         * @return Specific bit .
         */
        auto
        test(
                const moris::size_t index ) const
        -> decltype( mBitset.test( index ) )
         {
            return( mBitset.test( index ) );
        }

        /**
         * @brief Accesses specific bit .
         *
         * @Find first bit, which is on.
         */
        auto
        find_first()
        -> decltype( mBitset.find_first() )
        {
            return( mBitset.find_first() );
        }
//
//
//        /**
//         * @brief Accesses specific bit .
//         *
//         * @Find next bit, which is on.
//         */
//        auto
  //                const moris::size_t index)
//        -> decltype( mBitset.find_next( index ) )
//        {
//            return( mBitset.find_next( index ) );
//        }

        auto
        operator&()
        -> decltype( mBitset & mBitset )
        {
            return( mBitset & mBitset );
        }

    };
}

#endif /* SRC_CONTAINERS_CL_BOOSTBITSET_HPP_ */


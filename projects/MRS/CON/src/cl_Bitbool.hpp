/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Bitbool.hpp
 *
 */

#ifndef SRC_CONTAINERS_CL_BITBOOL_HPP_
#define SRC_CONTAINERS_CL_BITBOOL_HPP_

// C++ header files.
#include <string>
// MORIS library header files.
#include "typedefs.hpp" // COR/src

namespace moris
{

    class Bitbool
    {
    private:

        /**
         * MORIS Bitbool
         */
        std::vector< bool > mBitbool;

    public:

        /**
         * moris::Bitbool constructor
         */
        Bitbool() = default;

        /**
         * moris::Bitbool constructor
         *
         * @param[in] aBitbool      A Bitbool
         */
        Bitbool(
                const std::vector< bool > & aBitbool )
            : mBitbool( aBitbool )
        {
        }

        /**
         * moris::Bitbool constructor
         *
         * @param[in] aSize Size of the Cell to be initialized
         */
        Bitbool(
                moris::uint const aSize )
        : mBitbool( aSize )
        {
        }

//        Bitset(N & aNumBits )
//            : mBitset( aNumBits )
//        {
//        }

        /**
         * moris::Bitset destructor
         */
        ~Bitbool() = default; // 'default' tells the compiler to automatically
                             // delete the underlying Bitset

        /**
         * @brief Returns the number of elements in the Bitbool.
         *
         * @return Number of elements in the Bitbool.
         */
        moris::size_t
        size()
        {
            return mBitbool.size();
        }

        /**
         * const version of above
         */
        moris::size_t
        size() const
        {
            return mBitbool.size();
        }

        /**
         * Clears the contents of the Bitbool.
         */
        void
        clear()
        {
            mBitbool.clear();
        }

        /**
         * @brief Sets the position index in the Bitbool to 1.
         *
         * @return Bitset.
         */
        void
        set(
                const moris::size_t index )
        {
            mBitbool[ index ] = true;
        }

        /**
         * @brief Sets the position "index" in the Bitbool to 0.
         *
         * @return Bitset.
         */
        void
        reset(
                const moris::size_t index )
        {
            mBitbool[ index ] = false;
        }

        /**
         * @brief Accesses specific bit .
         *
         * @return Specific bit .
         */
        auto
        test(
                const moris::size_t index )
        -> decltype( mBitbool[ index ] )
        {
            return( mBitbool[ index ] );
        }
//
//
//        auto
//        operator&()
//        -> decltype( mBitset & mBitset )
//        {
//            return( mBitset & mBitset );
//        }

    };
}

#endif /* SRC_CONTAINERS_CL_BITBOOL_HPP_ */


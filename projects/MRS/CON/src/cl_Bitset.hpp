#ifndef SRC_CONTAINERS_CL_BITSET_HPP_
#define SRC_CONTAINERS_CL_BITSET_HPP_

// C++ header files.
#include <bitset>
#include <string>

// MORIS library header files.
#include "typedefs.hpp" //COR/src

namespace moris
{

    template < moris::size_t N >
    class Bitset
    {
// -----------------------------------------------------------------------------
    private:
// -----------------------------------------------------------------------------
        /**
         * MORIS Bitset
         */
        std::bitset<N> mBitset;
// -----------------------------------------------------------------------------
    public:
// -----------------------------------------------------------------------------

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
        Bitset(
                const std::bitset< N > & aBitset )
            : mBitset( aBitset )
        {
        }

// -----------------------------------------------------------------------------

        /**
         * moris::Bitset constructor
         *
         * @param[in] aChar an unsigned char
         */
        Bitset( const unsigned char & aChar )
        : mBitset( std::bitset<8>( aChar ) )
        {
        }

// -----------------------------------------------------------------------------

        /**
         * moris::Bitset constructor
         *
         * @param[in] aChar an unsigned char
         */
        Bitset( const unsigned int & aUint )
        : mBitset( std::bitset<32>( aUint ) )
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
        set(
                const moris::size_t index )
        -> decltype( mBitset.set( index ) )
        {
            return( mBitset.set( index ) );
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

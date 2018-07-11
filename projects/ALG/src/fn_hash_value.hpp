#ifndef SRC_ALGORITHMS_FN_HASH_VALUE_HPP_
#define SRC_ALGORITHMS_FN_HASH_VALUE_HPP_

// Third-party header files.
#include <boost/functional/hash.hpp>

// MORIS header files.
#include "typedefs.hpp" // COR/src
#include "cl_Mat.hpp" // LNA/src

namespace moris
{
    /**
     * hash_value function for moris matrices
     *
     * Overloads boost::hash_value to be used in conjunction with boost::hash.
     * For matrices, hash_value adds the hashed values of each component of the
     * matrix and returns that sum.
     *
     * @param[in] aMat Matrix to be hashed
     *
     * @return hash value for given matrix
     */
    template< typename T >
    moris::size_t
    hash_value(
            moris::Mat< T > const & aMat )
    {
        boost::hash<moris::real > hasher;
        moris::size_t             tResult = 0;

        for( moris::uint il = 0; il < aMat.numel() ; il++ )
        {
            tResult += hasher( aMat(il) );
        }

        return tResult;
    }
}

#endif /* SRC_ALGORITHMS_FN_HASH_VALUE_HPP_ */

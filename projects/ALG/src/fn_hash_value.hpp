/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_hash_value.hpp
 *
 */

#ifndef SRC_ALGORITHMS_FN_HASH_VALUE_HPP_
#define SRC_ALGORITHMS_FN_HASH_VALUE_HPP_

// Third-party header files.
#include <boost/functional/hash.hpp>

// MORIS header files.
#include "moris_typedefs.hpp" // COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

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
            Matrix < T > const & aMat )
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


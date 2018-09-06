/*
 * op_less.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_OP_LESS_HPP_
#define PROJECTS_LINALG_SRC_OP_LESS_HPP_


#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/op_less_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/op_less_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Element-wise check if one matrix is less than another..
     *
     * @param[in] aA Matrix.
     * @param[in] aB Matrix.
     *
     * @return Checks if: @f$ \mathbf{A}_{ij} <  \mathbf{B}_{ij} @f$ \n
     *         Element-wise equality evaluation of two objects; generates a
     *         matrix with entries that indicate whether at a given position
     *         the two elements from the two objects are less (1) or not (0).
     */
    template< typename Type, typename Matrix_Type >
    auto
    operator<(
            Matrix< Type, Matrix_Type > const & aA,
            Matrix< Type, Matrix_Type > const & aB )
    ->decltype( operator<( aA.matrix_data(), aB.matrix_data() ) )
    {
        return operator<( aA.matrix_data(), aB.matrix_data() );
    }

    template< typename Type, typename Matrix_Type >
    auto
    operator<(
            Matrix< Type, Matrix_Type > & aA,
            Matrix< Type, Matrix_Type > & aB )
    ->decltype( operator<( aA.matrix_data(), aB.matrix_data() ) )
    {
        return operator<( aA.matrix_data(), aB.matrix_data() );
    }

    template< typename Type, typename Matrix_Type >
    auto
    operator<( Type aA,
                Matrix< Type, Matrix_Type > & aB )
    ->decltype( operator<( aA, aB.matrix_data() ) )
    {
        return operator<( aA, aB.matrix_data() );
    }

    template< typename Type, typename Matrix_Type >
    auto
    operator<( Matrix< Type, Matrix_Type > & aA,
                Type                          aB)
    ->decltype( operator<( aA.matrix_data(), aB ) )
    {

        return operator<( aA.matrix_data(), aB );
    }

}



#endif /* PROJECTS_LINALG_SRC_OP_LESS_HPP_ */

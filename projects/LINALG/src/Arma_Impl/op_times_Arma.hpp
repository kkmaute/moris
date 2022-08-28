/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_times_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_OP_TIMES_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_OP_TIMES_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename Matrix_Type, typename ET >
    auto
    operator*(
            const ET                    & aA,
            Matrix< Matrix_Type > const & aB )
    ->decltype( aA * aB.matrix_data() )
    {
        return  aA * aB.matrix_data();
    }

    template< typename Matrix_Type, typename ET >
    auto
    operator*(
            Matrix< Matrix_Type > const & aA,
            const ET                    & aB)
    ->decltype( aA.matrix_data() * aB )
    {
        return  aA.matrix_data() * aB;
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_OP_TIMES_ARMA_HPP_ */


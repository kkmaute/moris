/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_div.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_OP_DIV_HPP_
#define PROJECTS_LINALG_SRC_OP_DIV_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename Matrix_Type >
    auto
    operator/(
            const  Matrix< Matrix_Type > & aMatrix,
            const typename Matrix< Matrix_Type >::Data_Type & aScalar )
    ->decltype ( aMatrix.matrix_data() / aScalar )
    {
        return aMatrix.matrix_data() / aScalar;
    }
}
#endif /* PROJECTS_LINALG_SRC_OP_DIV_HPP_ */


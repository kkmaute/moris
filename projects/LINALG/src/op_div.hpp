/*
 * op_div.hpp
 *
 *  Created on: Sep 5, 2018
 *      Author: messe
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

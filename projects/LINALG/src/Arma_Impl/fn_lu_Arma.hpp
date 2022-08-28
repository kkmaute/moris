/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_lu_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_LU_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_LU_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ETL, typename ETU , typename ETP, typename ETA  >
    void
    lu(     ETL       & aL,
            ETU       & aU,
            ETP       & aP,
            ETA const & aA )
    {
        if ( ! arma::lu( aL, aU, aP, aA ) )
        {
            throw std::runtime_error( "LU decomposition with partial pivoting failed." );
        }
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_LU_ARMA_HPP_ */


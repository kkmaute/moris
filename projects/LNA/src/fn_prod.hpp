#ifndef SRC_LINALG_FN_RPROD_HPP_
#define SRC_LINALG_FN_PROD_HPP_

#include "op_times.hpp" // LNA/src

namespace moris
{

    /**
     * Wrapper function for operator* of tensor products.
     *
     * Typically, the higher order tensor is given first.
     *
     * @param[in] aTensor1 First tensor, higher order
     * @param[in] aTensor2 Second tensor
     *
     */
    template< int Order1, int Order2, typename T1, typename T2, bool Sym1, bool Sym2, int Dim >
    auto
    prod(
            moris::Tensor< T1, Order1, Dim, Sym1 > const & aTensor1 ,
            moris::Tensor< T2, Order2, Dim, Sym2 > const & aTensor2 )
    -> decltype( aTensor1 * aTensor2 )
    {
        return aTensor1 * aTensor2;
    }

}

#endif /* SRC_LINALG_FN_PROD_HPP_ */

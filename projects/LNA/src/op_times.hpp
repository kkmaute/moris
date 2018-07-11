#ifndef MORIS_LINALG_OP_TIMES_HPP_
#define MORIS_LINALG_OP_TIMES_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"
#include "cl_Sp_Mat.hpp"
#include "cl_Tensor.hpp"
#include "fn_trans.hpp"

namespace moris
{
    /**
     * @brief Matrix multiplication
     *
     * @param[in] aA Input matrix
     * @param[in] aB Input matrix
     *
     * This function returns a matrix such that @f$ \mathbf{C}_{ij} = \mathbf{A}_{ik} \mathbf{B}_{kj} @f$
     * where C is the product of A and B. If A is an m-by-p and B is a p-by-n matrix, then C is an m-by-n matrix.
     *
     * Example:
     * @include LNA/src/op_times.inc
     *
     */
    template< typename T1, typename T2 >
    auto
    operator*(
            moris::Base_Mat< T1 > const & aA,
            moris::Base_Mat< T2 > const & aB )
    ->decltype( aA.data() * aB.data() )
    {
        return  aA.data() * aB.data();
    }

    template< typename T1, typename T2,
    bool M = ( moris::is_Mat< T2 >::value || moris::is_Sp_Mat< T2 >::value ),
    typename std::enable_if< ! M >::type* = nullptr >
    auto
    operator*(
            moris::Base_Mat< T1 > const & aA,
            T2                    const & aB )
    ->decltype( aA.data() * aB )
    {
        return aA.data() * aB;
    }

    template< typename T1, typename T2,
    bool M = ( moris::is_Mat< T1 >::value || moris::is_Sp_Mat< T1 >::value ),
    typename std::enable_if< ! M >::type* = nullptr >
    auto
    operator*(
            T1                    const & aA,
            moris::Base_Mat< T2 > const & aB )
    ->decltype( aA * aB.data() )
    {
        return aA * aB.data();
    }

    /**
     * @brief Tensor double-dot (:) product
     *
     * Performs the tensor product defined as @f$ A_{ij}:B_{ij} @f$.
     *
     * Currently, both tensors must be of the same type.
     *
     */
    template< typename T1, int Dim >
    T1
    operator*(
            moris::Tensor< T1, 2, Dim, true > const & aTensor1,
            moris::Tensor< T1, 2, Dim, true > const & aTensor2 )
    {
        moris::Mat< T1> result1 = moris::trans(aTensor1.data())*aTensor2.data()
                                + moris::trans(aTensor1.data().rows(aTensor1.dim(),aTensor1.data().n_rows()-1))
                                * aTensor2.data().rows(aTensor2.dim(),aTensor2.data().n_rows()-1);
        return result1(0, 0);

    }

    /**
     * @brief Tensor double-dot (:) product
     *
     * Performs the tensor for 2nd order, unsymmetric tensors.
     */
    template< typename T1, int Dim >
    T1
    operator*(
            moris::Tensor< T1, 2, Dim, false > const & aTensor1,
            moris::Tensor< T1, 2, Dim, false > const & aTensor2 )
    {
        moris::Mat< T1> result1 = moris::trans(aTensor1.data())*aTensor2.data();
        return result1(0, 0);
    }

    /**
     * @brief Tensor dot product
     *
     * Performs the tensor for 1st order tensors.
     */
    template< typename T1, int Dim , bool Sym1, bool Sym2 >
    T1
    operator*(
            moris::Tensor< T1, 1, Dim, Sym1 > const & aTensor1,
            moris::Tensor< T1, 1, Dim, Sym2 > const & aTensor2 )
    {
        moris::Mat< T1> result1 = moris::trans(aTensor1.data())*aTensor2.data();
        return result1(0, 0);
    }


    /**
     * @brief Tensor product of a 4th and 2nd order Tensors
     *
     * Performs the tensor product defined as @f$ s_{ij} = A_{ijkl}:B_{kl} @f$.
     *
     * Currently, both tensors must be of the same type.
     */
    template< typename T1, int Dim >
    moris::Tensor< T1, 2, Dim, true >
    operator*(
            moris::Tensor< T1, 4, Dim, true > const & aTensor1,
            moris::Tensor< T1, 2, Dim, true > const & aTensor2 )
//    -> decltype( aTensor1.data()*aTensor2.data()  )
    {

        moris::Mat<T1> tResult = aTensor1.data()*aTensor2.data();

        tResult.rows(0,aTensor1.data().n_rows()-1) += aTensor1.data().cols(aTensor1.dim(),aTensor1.data().n_rows()-1)
                                                    * aTensor2.data().rows(aTensor1.dim(),aTensor1.data().n_rows()-1);

        return moris::Tensor< T1, 2, Dim, true >(tResult);
    }

    /**
     * @brief Tensor product of a 4th and 2nd order Tensors
     *
     * Performs the tensor product defined as @f$ s_{ij} = A_{ijkl}:B_{kl} @f$.
     *
     */
    template< typename T1, int Dim >
    moris::Tensor< T1, 2, Dim, false >
    operator*(
            moris::Tensor< T1, 4, Dim, false > const & aTensor1,
            moris::Tensor< T1, 2, Dim, false > const & aTensor2 )
    {
        return moris::Tensor<T1, 2, Dim, false >(aTensor1.data()*aTensor2.data());
    }

    /**
     * @brief Tensor product of a 3rd and 1st order Tensors
     *
     * Performs the tensor product defined as @f$ s_{ij} = A_{ijk}B_{k} @f$.
     * Function is templated for symmetry, since the symmetry that A may have
     * in the first two indices does not affect the summation over its third index.
     * Thus, the product "inherits" the symmetry of A.
     *
     * Note that this is different that the product @f$ s_{ik} = A_{ijk}B_{j} @f$,
     * which is not supported by the * operator.
     *
     */
    template< typename T1 , int Dim, bool Sym1, bool Sym2 >
//    moris::Tensor< T1, 2, Dim, Sym1 >
    auto
    operator*(
            moris::Tensor< T1, 3, Dim, Sym1 > const & aTensor1,
            moris::Tensor< T1, 1, Dim, Sym2 > const & aTensor2 )
    -> decltype( aTensor1.data()*aTensor2.data() )
    {
        return aTensor1.data()*aTensor2.data();
    }

// -------------------------------------------------------------------------------------

    /**
     * @brief Tensor product of a 2nd and 1st order 2D Tensors
     *
     * Performs the tensor product defined as @f$ s_{i} = T_{ij}v_{j} @f$.
     *
     * Because this product will be written out explicitly, the symmetric of the 2nd order tensor
     * does not matter. This function is templated to allow any combination of symmetric/unsymmetric
     * 2nd and 1st order tensors. (Though, for a 1st order tensor, the use of symmetry is lost).
     */
    template< typename T1, bool Sym >
    moris::Mat< T1 >
    operator*(
            moris::Tensor< T1, 2, 2, true > const & aTensor1,
            moris::Tensor< T1, 1, 2, Sym >  const & aTensor2 )
    {

        moris::Mat< T1 > tVector(2, 1);

        tVector(0) = aTensor1[0]*aTensor2[0] + aTensor1[2]*aTensor2[1];
        tVector(1) = aTensor1[2]*aTensor2[0] + aTensor1[1]*aTensor2[1];

        return tVector;
    }

    template< typename T1, bool Sym >
    moris::Mat< T1 >
    operator*(
            moris::Tensor< T1, 2, 2, false > const & aTensor1,
            moris::Tensor< T1, 1, 2, Sym >   const & aTensor2 )
    {

        moris::Mat< T1 > tVector(2, 1);

        tVector(0) = aTensor1[0]*aTensor2[0] + aTensor1[2]*aTensor2[1];
        tVector(1) = aTensor1[3]*aTensor2[0] + aTensor1[1]*aTensor2[1];

        return tVector;
    }

    /**
     * @brief Tensor product of a 2nd and 1st order 3D Tensors
     *
     * Performs the tensor product defined as @f$ s_{i} = T_{ij}v_{j} @f$.
     *
     * Because this product will be written out explicitly, the symmetric of the 2nd order tensor
     * does not matter. This function is templated to allow any combination of symmetric/unsymmetric
     * 2nd and 1st order tensors. (Though, for a 1st order tensor, the use of symmetry is lost).
     */
    template< typename T1 , bool Sym >
    moris::Mat< T1 >
    operator*(
            moris::Tensor< T1, 2, 3, true > const & aTensor1,
            moris::Tensor< T1, 1, 3, Sym >  const & aTensor2 )
    {

        moris::Mat< T1 > tVector(3, 1);

        tVector(0) = aTensor1[0]*aTensor2[0] + aTensor1[5]*aTensor2[1] + aTensor1[4]*aTensor2[2];
        tVector(1) = aTensor1[5]*aTensor2[0] + aTensor1[1]*aTensor2[1] + aTensor1[3]*aTensor2[2];
        tVector(2) = aTensor1[4]*aTensor2[0] + aTensor1[3]*aTensor2[1] + aTensor1[2]*aTensor2[2];

        return tVector; // moris::Tensor<T1, 1, 3, true>(tVector);
    }

    template< typename T1, bool Sym >
    moris::Mat< T1 >
    operator*(
            moris::Tensor< T1, 2, 3, false > const & aTensor1,
            moris::Tensor< T1, 1, 3, Sym >   const & aTensor2 )
    {

        moris::Mat< T1 > tVector(3, 1);

        tVector(0) = aTensor1[0]*aTensor2[0] + aTensor1[5]*aTensor2[1] + aTensor1[4]*aTensor2[2];
        tVector(1) = aTensor1[8]*aTensor2[0] + aTensor1[1]*aTensor2[1] + aTensor1[3]*aTensor2[2];
        tVector(2) = aTensor1[7]*aTensor2[0] + aTensor1[6]*aTensor2[1] + aTensor1[2]*aTensor2[2];

        return tVector;
    }

}

#endif  /* MORIS_LINALG_OP_TIMES_HPP_ */


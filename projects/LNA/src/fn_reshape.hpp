/*
 * fn_reshape.hpp
 *
 *  Created on: Jul 10, 2017
 *      Author: gleim
 */

#ifndef SRC_LINALG_FN_RESHAPE_HPP_
#define SRC_LINALG_FN_RESHAPE_HPP_

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template<  typename T1, typename T2, typename T3 >
    moris::Mat<T1>
    reshape(
            moris::Mat< T1 > const & aA,
            T2               const & aB,
            T3               const & aC)
            {
        return arma::conv_to< arma::Mat<T1>>::from(reshape((arma::Mat<T1>)aA.data(),aB,aC));
            }

    template<  typename T1, typename T2, typename T3 >
    moris::Mat<T1>
    reshape(
            moris::Sp_Mat< T1 > const & aA,
            T2               const & aB,
            T3               const & aC)
            {
        return arma::conv_to< arma::Mat<T1>>::from(reshape((arma::Mat<T1>)aA.data(),aB,aC));
            }
}

#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template<  typename T1, typename T2, typename T3 >
    moris::Mat< T1 >
    reshape(
            moris::Mat< T1 > const & aA,
            T2               const & aB,
            T3               const & aC)
            {
        MORIS_LOG_INFO << "Eigen doesn't support reshape, an explicit loop is implemented.\n";
        moris::Mat<T1> tMat(aB,aC);
        moris::uint tvar = 0;
        moris::uint ncol = aA.n_cols();
        moris::uint nrow = aA.n_rows();
        for(moris::uint i = 0; i<(moris::uint)aB; i++)
        {
            tvar = 0;
            for(moris::uint j = 0; j<(moris::uint)aC; j++)
            {
                tMat(i,j) = aA(tvar*aB+i);
                tvar++;
            }
        }
        return tMat;
            }

    template<  typename T1, typename T2, typename T3 >
    moris::Mat< T1 >
    reshape(
            moris::Sp_Mat< T1 > const & aA,
            T2               const & aB,
            T3               const & aC)
            {
        MORIS_LOG_INFO << "Eigen doesn't support reshape, an explicit loop is implemented.";
        moris::Mat<T1> tMat(aB,aC);
        moris::uint tvar = 0;
        moris::uint ncol = aA.n_cols();
        moris::uint nrow = aA.n_rows();
        for(moris::uint i = 0; i<(moris::uint)aB; i++)
        {
            tvar = 0;
            for(moris::uint j = 0; j<(moris::uint)aC; j++)
            {
                tMat(i,j) = aA(tvar*aB+i);
                tvar++;
            }
        }
        return tMat;
            }

}
#endif
// ----------------------------------------------------------------------------

namespace moris
{

    /**
     * @brief find vector.
     *
     * @param[in] aMat     Vector.
     * @param[in] WhichValue    Value, which is in the vector.
     *
     * @return  Vector of items found
     */
    template<  typename T1, typename T2, typename T3 >
    auto
    reshape(
            moris::Mat< T1 > const & aA,
            T2               const & aB,
            T3               const & aC)
    -> decltype( moris::Math::reshape( aA,aB,aC) )
    {
        return moris::Math::reshape( aA,aB,aC );
    }

    template<  typename T1, typename T2, typename T3 >
    auto
    reshape(
            moris::Sp_Mat< T1 > const & aA,
            T2               const & aB,
            T3               const & aC)
    -> decltype( moris::Math::reshape( aA,aB,aC) )
    {
        return moris::Math::reshape( aA,aB,aC );
    }


}


#endif /* SRC_LINALG_FN_RESHAPE_HPP_ */

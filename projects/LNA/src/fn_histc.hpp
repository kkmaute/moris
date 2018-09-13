/*
 * fn_histc.hpp
 *
 *  Created on: Jul 10, 2017
 *      Author: gleim
 */

#ifndef SRC_LINALG_FN_HISTC_HPP_
#define SRC_LINALG_FN_HISTC_HPP_

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template<  typename T1, typename T2 >
    moris::Mat<T1>
    histc(
            moris::Mat< T1 > const & aA,
            moris::Mat< T2 > const & aB)
            {
        return arma::conv_to< arma::Mat<T1>>::from(histc((arma::Mat<T1>)aA.data(),(arma::Mat<T1>)aB.data()));
            }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template<  typename T1, typename T2 >
    moris::Mat<T1>
    histc(
            moris::Mat< T1 > const & aA,
            moris::Mat< T2 > const & aB)
            {
        // Eigen does not have an internal histc function
//        MORIS_LOG_INFO << "Eigen doesn't support histc, an explicit loop is implemented.";
        moris::Mat<moris::uint> CountNumber(aB.numel(),1,0);
        moris::uint tvar=0; // temporary variable for a loop
        for(moris::uint i=0; i<aB.numel(); i++)
        {
            tvar = 0;
            for(moris::uint j=0; j<aA.numel(); j++)
            {
                if(aB(i) == aA(j))
                {
                    tvar++;
                }
            }
            CountNumber(i)=tvar;
        }
        return CountNumber;
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
    template<  typename T1, typename T2 >
    auto
    histc(
            moris::Mat< T1 > const & aA,
            moris::Mat< T2 > const & aB)
    -> decltype( moris::Math::histc( aA, aB ) )
    {
        return moris::Math::histc( aA,aB );
    }


}


#endif /* SRC_LINALG_FN_HISTC_HPP_ */

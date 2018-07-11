/*
 * fn_sort.hpp
 *
 *  Created on: Feb 20, 2017
 *      Author: doble
 */

#ifndef SRC_LINALG_FN_SORT_HPP_
#define SRC_LINALG_FN_SORT_HPP_

// ----------------------------------------------------------------------------
// FIXME: Behavior of sort function is not clear; in eigen the input vector content is changed; what's about arma

// FIXME: There is no check of dimension of moris::mat; should sort also work on matrices?

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template<typename T>
    moris::Mat<T>
    sort( moris::Mat<T> &aMat )
    {
        return arma::sort(aMat.data());
    }

    // aDir means sort direction, ascend by default
    template<typename T>
    moris::Mat<T>
    sort( moris::Mat<T> & aMat,
          T               aDir)
    {
        return arma::sort(aMat.data(),aDir);
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template<typename T>
    moris::Mat<T>
    sort( moris::Mat<T> &aMat )
    {
        // Eigen does not have an internal sort function
        // 1. Create copy  of input  matrix
        // 2. Sort matrix using std::sort
        moris::Mat<T> tMat = aMat;
        T* tData = (tMat.data()).data();
        moris::uint tLen = tMat.numel();
        std::sort(tData,tData+tLen);
        return tMat;
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{

    /**
     * @brief Sort vector.
     *
     * @param[in] aMat Vector.
     *
     * @return  sorted vector
     */
    template< typename T>
    moris::Mat<T>
    sort( moris::Mat<T> &aMat )
    {
        return moris::Math::sort( aMat );
    }
}


#endif /* SRC_LINALG_FN_SORT_HPP_ */

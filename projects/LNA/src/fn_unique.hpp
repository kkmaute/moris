/*
 * fn_unique.hpp
 *
 *  Created on: Apr 18, 2017
 *      Author: gleim
 */

#ifndef SRC_LINALG_FN_UNIQUE_HPP_
#define SRC_LINALG_FN_UNIQUE_HPP_

#include "fn_isrow.hpp"

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template<typename T>
    moris::Mat<T>
    unique( moris::Mat<T> &aMat )
    {
        return arma::unique(aMat.data());
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template<typename T>
    moris::Mat<T>
    unique( moris::Mat<T> &aMat )
    {
        // Eigen does not have an internal unique function
        // 1. Create copy  of input  matrix
        // 2. Sort matrix using std::sort
        // 3. Taking only the unique part and saving in ttMat
        // 4. Finding position and resize tMat
        moris::Mat<T> tMat = aMat;
        T* tData = (tMat.data()).data();
        moris::uint tLen = tMat.numel();
        std::sort(tData,tData+tLen);
        auto last  = std::unique(tData, tData+tLen);
        auto pos = std::distance(tData, last);
        if (moris::Math::isrow( aMat ) == 1)
        {
            tMat.resize(1,pos);
        }
        else
        {
            tMat.resize(pos,1);
        }
        return tMat;
    }
}
#endif
// ----------------------------------------------------------------------------

namespace moris
{

    /**
     * @brief unique vector.
     *
     * @param[in] aMat Vector.
     *
     * @return  sorted unique vector
     */
    template< typename T>
    moris::Mat<T>
    unique( moris::Mat<T> &aMat )
    {
        return moris::Math::unique( aMat );
    }
}


#endif /* SRC_LINALG_FN_UNIQUE_HPP_ */


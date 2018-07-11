/*
 * fn_find_unique.hpp
 *
 *  Created on: Sep 8, 2017
 *      Author: gleim
 */

#ifndef SRC_LINALG_FN_FIND_UNIQUE_HPP_
#define SRC_LINALG_FN_FIND_UNIQUE_HPP_

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    moris::Mat< uint >
    find_unique(moris::Mat< T > const & aA)
    {
        return arma::conv_to< arma::Mat<moris::uint>>::from(find_unique((arma::Mat<T>)aA.data(),false));
    }

}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    moris::Mat<uint>
    find_unique(
            moris::Mat< T > const & aA)
            {
        // Eigen does not have an internal find function
        //        MORIS_LOG_INFO << "Eigen doesn't support find, an explicit loop is implemented.";
        //Create first a unique list
        moris::Mat<T> tMat = aA;
        T* tData = (tMat.data()).data();
        moris::uint tLen = tMat.numel();
        std::sort(tData,tData+tLen);
        auto last  = std::unique(tData, tData+tLen);
        auto pos = std::distance(tData, last);
        if (moris::Math::isrow( aA ) == 1)
        {
            tMat.resize(1,pos);
        }
        else
        {
            tMat.resize(pos,1);
        }
        //Use the unique list to find the unique entries in the original matrix
        moris::Mat<T> tMatFind = aA;
        moris::Mat<moris::uint> WhereInMat(tLen,1,0); // Wector with the specific value
        uint tvar=0; // temporary variable for a loop
        for(uint i=0; i<pos; i++)
        {
            for(uint j=0; j<tLen; j++)
            {
                if( tMat( i ) == tMatFind( j ) )
                {
                    WhereInMat(tvar)= j;
                    tvar++;
                    break;
                }
            }
        }
        WhereInMat.resize(tvar,1);
        return WhereInMat;
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
    template< typename T >
    auto
    find_unique(
            moris::Mat< T > const & aA)
    -> decltype( moris::Math::find_unique( aA) )
    {
        return moris::Math::find_unique( aA );
    }

}


#endif /* SRC_LINALG_FN_FIND_UNIQUE_HPP_ */

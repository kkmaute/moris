/*
 * fn_find.hpp
 *
 *  Created on: Jul 10, 2017
 *      Author: gleim
 */

#ifndef SRC_LINALG_FN_FIND_HPP_
#define SRC_LINALG_FN_FIND_HPP_

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    //Find all indices
    template< typename T >
    moris::Mat< uint >
    find(moris::Mat< T > const & aA)
    {
        return arma::conv_to< arma::Mat<moris::uint>>::from(find((arma::Mat<T>)aA.data()));
    }

    template< typename T >
    moris::Mat< uint >
    find(moris::Sp_Mat< T > const & aA)
    {
        return arma::conv_to< arma::Mat<moris::uint>>::from(find((arma::Mat<T>)aA.data()));
    }

    // Find a specific number (ab) of indices
    template< typename T >
    moris::Mat< uint >
    find(moris::Mat< T > const & aA,
            moris::uint const & ab)
    {
        return arma::conv_to< arma::Mat<moris::uint>>::from(find((arma::Mat<T>)aA.data(),ab));
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    moris::Mat<uint>
    find(
            moris::Mat< T > const & aA)
    {
        // Eigen does not have an internal find function
//        MORIS_LOG_INFO << "Eigen doesn't support find, an explicit loop is implemented.";
        moris::Mat<T> tMat = aA;
        moris::uint tLen = tMat.numel();
        moris::Mat<moris::uint> WhereInMat(tLen,1,0); // Wector with the specific value
        uint tvar=0; // temporary variable for a loop
        for(uint i=0; i<tLen; i++)
        {
            if(tMat(i) == 1)
            {
                WhereInMat(tvar)= i;
                tvar++;
            }
        }
        WhereInMat.resize(tvar,1);
        return WhereInMat;
    }

    template< typename T >
    moris::Mat<uint>
    find(
            moris::Sp_Mat< T > const & aA)
    {
        // Eigen does not have an internal find function
        MORIS_LOG_INFO << "Eigen doesn't support find, an explicit loop is implemented.";
        moris::Sp_Mat<T> tMat = aA;
        moris::uint tLen = tMat.numel();
        moris::Mat<moris::uint> WhereInMat(tLen,1,0); // Vector with the specific value
        uint tvar=0; // temporary variable for a loop
        for(uint i=0; i<tLen; i++)
        {
            if(tMat(i) == 1)
            {
                WhereInMat(tvar)= i;
                tvar++;
            }
        }
        WhereInMat.resize(tvar,1);
        return WhereInMat;
    }

// Find a a number of indices
    template< typename T >
    moris::Mat<uint>
    find(
            moris::Mat< T > const & aA,
            moris::uint const & ab)
    {
        // Eigen does not have an internal find function
//        MORIS_LOG_INFO << "Eigen doesn't support find, an explicit loop is implemented.";
        moris::Mat<T> tMat = aA;
        moris::uint tLen = tMat.numel();
        moris::Mat<moris::uint> WhereInMat(tLen,1,0); // Wector with the specific value
        uint tvar=0; // temporary variable for a loop
        for(uint i=0; i<tLen; i++)
        {
            if(tMat(i) == 1)
            {
                WhereInMat(tvar)= i;
                tvar++;
                break;
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
     *
     * @return  Vector of items found
     */
    template< typename T >
    auto
    find(
            moris::Mat< T > const & aA)
    -> decltype( moris::Math::find( aA) )
    {
        return moris::Math::find( aA );
    }

    template< typename T >
    auto
    find(
            moris::Sp_Mat< T > const & aA)
    -> decltype( moris::Math::find( aA) )
    {
        return moris::Math::find( aA );
    }

    /**
     * @brief find vector.
     *
     * @param[in] aMat     Vector.
     * @param[in] ab    Find a specific number of indices.
     *
     * @return  Vector of items found
     */
    template< typename T >
    auto
    find(
            moris::Mat< T > const & aA,
            moris::uint const & ab)
    -> decltype( moris::Math::find( aA,ab) )
    {
        return moris::Math::find( aA,ab );
    }

}


#endif /* SRC_LINALG_FN_FIND_HPP_ */

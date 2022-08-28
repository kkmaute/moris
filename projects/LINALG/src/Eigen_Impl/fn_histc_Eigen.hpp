/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_histc_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_HISTC_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_HISTC_EIGEN_HPP_
#include <Eigen/Dense>

namespace moris
{
    template< typename T1 >
    Eigen::Matrix<T1,Eigen::Dynamic,Eigen::Dynamic>
    histc( const Eigen::Matrix<T1,Eigen::Dynamic,Eigen::Dynamic> & aA,
           const Eigen::Matrix<T1,Eigen::Dynamic,Eigen::Dynamic> & aB)
    {
        // Eigen does not have an internal histc function
        //typedef Eigen::Matrix< uint, aB.size(), 1 > aCountNumber;
    	Eigen::Matrix<T1,Eigen::Dynamic,Eigen::Dynamic> tCountNumber(aB.size(),1);
    	tCountNumber.fill(0);
        moris::uint tVar = 0;

        for( moris::uint i = 0; i < aB.size(); i++ )
        {
            tVar = 0;
            for( moris::uint j = 0; j < aA.size(); j++ )
            {
                if( aB( i, 0 ) == aA( j, 0 ) )
                {
                    tVar++;
                }
            }
            tCountNumber( i, 0 ) = tVar;
        }
        return tCountNumber;
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_HISTC_EIGEN_HPP_ */


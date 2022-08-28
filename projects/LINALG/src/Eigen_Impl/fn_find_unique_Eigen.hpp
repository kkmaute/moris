/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_find_unique_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_FIND_UNIQUE_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_FIND_UNIQUE_EIGEN_HPP_
#include <Eigen/Dense>

namespace moris
{
    template< typename T1 >
    Eigen::Matrix< uint, Eigen::Dynamic, Eigen::Dynamic >
    find_unique( const Eigen::Matrix<T1,Eigen::Dynamic,Eigen::Dynamic> & aA )
    {
        // make sure that matrix is row or col matrix
        MORIS_ASSERT( aA.rows() == 1 || aA.cols() == 1, "unique() can only be performed on a vector" );

        // Eigen does not have an internal find unique function
        moris::uint tLength = aA.size();

        Eigen::Matrix< T1, Eigen::Dynamic,Eigen::Dynamic> aUniqueMatrix = aA;

        // get pointer to raw data
        T1 * tData = aUniqueMatrix.data();

        // sort data
        std::sort( tData, tData+tLength );

        // find positions
        auto tLast  = std::unique( tData, tData+tLength );
        auto tPos   = std::distance( tData, tLast );

        // resize matrix
        if ( aUniqueMatrix.rows() == 1)
        {
            aUniqueMatrix.conservativeResize( 1, tPos );
        }
        else
        {
            aUniqueMatrix.conservativeResize( tPos, 1 );
        }

        // Create position matrix
        Eigen::Matrix< uint, Eigen::Dynamic,Eigen::Dynamic> tPosition( tPos, 1 );

        tPosition.fill( 0 );

        moris::uint tVar = 0;

        // Loop over all values in unique vector
        for( uint Ik = 0; Ik < tPos; Ik++ )
        {
            // Loop over all values in original vector
            for( uint Ij = 0; Ij < tLength; Ij++ )
            {
                // If value in unique vector and value in original vector are the same, save position and break loop
                if( equal_to( aUniqueMatrix( Ik ), aA( Ij ) ) )
                {
                    tPosition( tVar++ )= Ij;
                    break;
                }
            }
        }

        return tPosition;
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_FIND_UNIQUE_EIGEN_HPP_ */


/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_find_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_FIND_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_FIND_EIGEN_HPP_
#include <Eigen/Dense>
#include "fn_equal_to.hpp"

namespace moris
{
    template< typename T1 >
    Eigen::Matrix< uint, Eigen::Dynamic, Eigen::Dynamic >
    find( const Eigen::Matrix<T1,Eigen::Dynamic,Eigen::Dynamic> & aA )
    {
        // Eigen does not have an internal find function
        moris::uint tLength = aA.size();

        Eigen::Matrix< uint, Eigen::Dynamic,Eigen::Dynamic> tPosition( tLength, 1 );

        tPosition.fill( 0 );

        moris::uint tVar = 0;

        for( uint i = 0; i < tLength; i++ )
        {
            if( !equal_to( aA( i ), 0) )
            {
                tPosition( tVar )= i;
                tVar++;
            }
        }

        tPosition.conservativeResize(tVar,1);

        return tPosition;
    }

    template< typename T1 >
    Eigen::Matrix< uint, Eigen::Dynamic, Eigen::Dynamic >
    find( const Eigen::Matrix<T1,Eigen::Dynamic,Eigen::Dynamic> & aA,
          const moris::uint                                     & ab )
    {
        // Eigen does not have an internal find function
        moris::uint tLength = aA.size();

        Eigen::Matrix< uint, Eigen::Dynamic, Eigen::Dynamic > tPosition( tLength, 1 );

        tPosition.fill( 0 );

        moris::uint tVar = 0;

        for( uint i = 0; i < tLength; i++ )
        {
            if( !equal_to( aA( i ), 0) )
            {
                tPosition( tVar++ )= i;

                if ( tVar == ab )
                {
                    break;
                }
            }
        }
        tPosition.conservativeResize( tVar, 1 );
        return tPosition;
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_FIND_EIGEN_HPP_ */


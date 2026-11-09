/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cluster.cpp
 *
 */

#include "cl_MTK_Cluster.hpp"
#include <cmath>

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Cluster::compute_cluster_cell_measure_derivatives(
            const Matrix< DDRMat >& aVertexCoordinates,
            const Matrix< DDSMat >& aRequested,
            Primary_Void            aPrimaryOrVoid,
            Leader_Follower         aIsLeader ) const
    {
        MORIS_ERROR( aRequested.n_rows() == aVertexCoordinates.n_rows()
                             && aRequested.n_cols() == aVertexCoordinates.n_cols(),
                "Cluster::compute_cluster_cell_measure_derivatives - Invalid request mask." );

        Matrix< DDRMat > tResult( aVertexCoordinates.n_rows(), aVertexCoordinates.n_cols(), 0.0 );

        for ( uint iVertex = 0; iVertex < aVertexCoordinates.n_rows(); ++iVertex )
        {
            const Matrix< DDRMat > tCoordinates = aVertexCoordinates.get_row( iVertex );

            for ( uint iDirection = 0; iDirection < aVertexCoordinates.n_cols(); ++iDirection )
            {
                if ( aRequested( iVertex, iDirection ) != 0 )
                {
                    // Unknown cluster implementations retain their scalar virtual behavior.
                    tResult( iVertex, iDirection ) = this->compute_cluster_cell_measure_derivative(
                            tCoordinates, iDirection, aPrimaryOrVoid, aIsLeader );
                }
            }
        }

        return tResult;
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Cluster::compute_cell_measure_derivatives_from_cells(
            const Vector< Cell const * >& aCells,
            const Matrix< DDRMat >&       aVertexCoordinates,
            const Matrix< DDSMat >&       aRequested )
    {
        MORIS_ERROR( aRequested.n_rows() == aVertexCoordinates.n_rows()
                             && aRequested.n_cols() == aVertexCoordinates.n_cols(),
                "Cluster::compute_cell_measure_derivatives_from_cells - Invalid request mask." );

        const uint       tDimension = aVertexCoordinates.n_cols();
        Matrix< DDRMat > tResult( aVertexCoordinates.n_rows(), tDimension, 0.0 );
        bool             tAnyRequested = false;

        for ( uint i = 0; i < aRequested.numel(); ++i )
        {
            tAnyRequested = tAnyRequested || aRequested.data()[ i ] != 0;
        }

        if ( !tAnyRequested )
        {
            return tResult;
        }

        for ( Cell const * tCell : aCells )
        {
            const Matrix< DDRMat > tCellCoordinates = tCell->get_vertex_coords();
            MORIS_ERROR( tCellCoordinates.n_cols() == tDimension,
                    "Cluster::compute_cell_measure_derivatives_from_cells - Coordinate dimensions must match." );

            Vector< moris_index > tMatchedNodes( aVertexCoordinates.n_rows(), -1 );
            Matrix< DDSMat >      tCellRequested( tCellCoordinates.n_rows(), tDimension, 0 );
            uint                  tNumRequested = 0;

            for ( uint iVertex = 0; iVertex < aVertexCoordinates.n_rows(); ++iVertex )
            {
                bool tActive = false;

                for ( uint iDirection = 0; iDirection < tDimension; ++iDirection )
                {
                    tActive = tActive || aRequested( iVertex, iDirection ) != 0;
                }

                if ( !tActive )
                {
                    continue;
                }

                for ( uint iNode = 0; iNode < tCellCoordinates.n_rows(); ++iNode )
                {
                    // Preserve the old coordinate-based correspondence, including
                    // coincident enriched vertices and the strict 1e-12 tolerance.
                    real tSquaredDistance = 0.0;

                    for ( uint iDirection = 0; iDirection < tDimension; ++iDirection )
                    {
                        const real tDelta = tCellCoordinates( iNode, iDirection ) - aVertexCoordinates( iVertex, iDirection );

                        if ( std::abs( tDelta ) >= 1e-12 )
                        {
                            tSquaredDistance = 1e-24;
                            break;
                        }

                        tSquaredDistance += tDelta * tDelta;
                    }

                    if ( tSquaredDistance < 1e-24 )
                    {
                        tMatchedNodes( iVertex ) = iNode;

                        for ( uint iDirection = 0; iDirection < tDimension; ++iDirection )
                        {
                            if ( aRequested( iVertex, iDirection ) != 0 && tCellRequested( iNode, iDirection ) == 0 )
                            {
                                tCellRequested( iNode, iDirection ) = 1;
                                ++tNumRequested;
                            }
                        }

                        break;
                    }
                }
            }

            if ( tNumRequested == 0 )
            {
                continue;
            }

            Matrix< DDRMat > tCellDerivatives;

            if ( tNumRequested == 1 )
            {
                // A single query retains the inexpensive scalar path.
                tCellDerivatives.set_size( tCellCoordinates.n_rows(), tDimension, 0.0 );

                for ( uint iNode = 0; iNode < tCellCoordinates.n_rows(); ++iNode )
                {
                    for ( uint iDirection = 0; iDirection < tDimension; ++iDirection )
                    {
                        if ( tCellRequested( iNode, iDirection ) != 0 )
                        {
                            tCellDerivatives( iNode, iDirection ) = tCell->compute_cell_measure_deriv( iNode, iDirection );
                        }
                    }
                }
            }
            else
            {
                tCell->compute_cell_measure_derivatives( tCellDerivatives, tCellRequested );
            }

            for ( uint iVertex = 0; iVertex < tMatchedNodes.size(); ++iVertex )
            {
                if ( tMatchedNodes( iVertex ) >= 0 )
                {
                    for ( uint iDirection = 0; iDirection < tDimension; ++iDirection )
                    {
                        if ( aRequested( iVertex, iDirection ) != 0 )
                        {
                            tResult( iVertex, iDirection ) += tCellDerivatives( tMatchedNodes( iVertex ), iDirection );
                        }
                    }
                }
            }
        }

        return tResult;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::mtk

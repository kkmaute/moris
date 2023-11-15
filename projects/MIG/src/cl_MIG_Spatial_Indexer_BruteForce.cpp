//
// Created by frank on 11/13/23.
//

#include "cl_MIG_Spatial_Indexer_BruteForce.h"

namespace moris::mig
{
    moris::map< moris_index, moris_index > Spatial_Indexer_BruteForce::perform( real epsilon )
    {
        map< moris_index, moris_index > tClosestVertices;
        auto                            tCoordinates = get_deformed_coordinates();
        auto                            tNumVertices = static_cast< moris_index >( mCoordinates.n_rows() );

        // outer loop over all coordinates
        for ( moris_index first_vertex = 0; first_vertex < tNumVertices; first_vertex++ )
        {
            moris_index tClosestVertex = MORIS_INDEX_MAX;
            real        tMinDistance   = MORIS_REAL_MAX;

            // inner loop to compare all other coordinates to the first one
            for ( moris_index second_vertex = 0; second_vertex < tNumVertices; second_vertex++ )
            {

                // skip the same vertex
                if ( first_vertex == second_vertex )
                {
                    continue;
                }

                // check if the second vertex is a neighbor of the first one
                bool tIsNeighbor = false;
                for ( auto& tNeighbor : mNeighbors( first_vertex ) )
                {
                    if ( tNeighbor == second_vertex )
                    {
                        tIsNeighbor = true;
                        break;
                    }
                }

                // if the second vertex is a neighbor, skip it
                if ( tIsNeighbor ) continue;

                // check if the distance between the two vertices is smaller than the current distance
                auto tDistance = norm( tCoordinates.get_row( first_vertex ) - tCoordinates.get_row( second_vertex ) );
                if ( tDistance < tMinDistance )
                {
                    tMinDistance   = tDistance;
                    tClosestVertex = second_vertex;
                }
            }    // end inner loop over all coordinates

            if ( tMinDistance < epsilon )
            {
                tClosestVertices[ first_vertex ] = tClosestVertex;
            }
        }    // end outer loop over all coordinates

        return tClosestVertices;
    }
}    // namespace moris::mig

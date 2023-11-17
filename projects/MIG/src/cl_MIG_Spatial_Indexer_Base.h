//
// Created by frank on 11/13/23.
//

#ifndef MORIS_CL_MIG_SPATIAL_INDEXER_BASE_H
#define MORIS_CL_MIG_SPATIAL_INDEXER_BASE_H

#include "cl_Matrix.hpp"
#include "cl_Map.hpp"
#include "cl_MTK_Surface_Mesh.hpp"
namespace moris::mig
{
    class Spatial_Indexer_Base
    {
      protected:
        /**
         * @brief List of coordinates that are used to query the spatial index
         */
        Matrix< DDRMat > mCoordinates;

        /**
         * @brief List of neighbors for each coordinate to sort out neighbors as "false" closest points
         * First Cell is for all coordinates, second is for all neighbors of the respective coordinate
         */
        Cell< Cell< moris_index > > mNeighbors;

        /**
         * @brief List of displacements for each coordinate to be applied before the spatial index.
         * @details This is required to use the spatial indexing class in every step of the Newton algorithm. It will
         * update the new solution vector that will be then be used to query the new deformed mesh.
         */
        Matrix< DDRMat > mDisplacements;

        Matrix< DDRMat > mVertexNormals;

      public:
        Spatial_Indexer_Base(
                const Matrix< DDRMat >            &aCoordinates,
                const Cell< Cell< moris_index > > &aNeighbors,
                const Matrix< DDRMat >            &aDisplacements,
                const Matrix< DDRMat >            &aVertexNormals )
                : mCoordinates( aCoordinates )
                , mNeighbors( aNeighbors )
                , mDisplacements( aDisplacements )
                , mVertexNormals( aVertexNormals )
        {
        }

        Spatial_Indexer_Base(
                moris::mtk::Surface_Mesh       aSurfaceMesh,
                moris::Matrix< DDRMat > const &aDisplacements )
                : mCoordinates( aSurfaceMesh.get_vertex_coordinates() )
                , mNeighbors( aSurfaceMesh.get_vertex_neighbors() )
                , mDisplacements( aDisplacements )
                , mVertexNormals( aSurfaceMesh.get_vertex_normals() )
        {
        }

        virtual ~Spatial_Indexer_Base() = default;

        virtual moris::map< moris_index, moris_index > perform( real epsilon ) = 0;

      protected:
        Matrix< DDRMat > get_deformed_coordinates()
        {
            MORIS_ASSERT( mCoordinates.numel() == mDisplacements.numel(),
                    "Spatial_Index: The number of coordinates and displacements must be the same." );
            auto tDeformed = mCoordinates.copy();
            tDeformed += mDisplacements;
            return tDeformed;
        }
    };


}    // namespace moris::mig

#endif    // MORIS_CL_MIG_SPATIAL_INDEXER_BASE_H

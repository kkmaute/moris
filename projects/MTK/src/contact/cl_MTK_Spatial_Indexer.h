//
// Created by frank on 11/13/23.
//

#ifndef MORIS_CL_MTK_SPATIAL_INDEXER_H
#define MORIS_CL_MTK_SPATIAL_INDEXER_H

#include "cl_Matrix.hpp"
#include "cl_Map.hpp"
#include "cl_MTK_Surface_Mesh.hpp"

namespace moris::mtk
{
    /**
     * @brief Holds the result of a spatial indexing operation (i.e. closest point search) for a single source surface mesh and (potentially) multiple target surface meshes.
     */
    struct Spatial_Indexing_Result
    {
        /**
         * @brief Contains all the information about the closest target vertex (index, mesh index and distance).
         */
        struct Target_Info
        {
            moris_index vertex;
            moris_index mesh_index;
            real        distance;
        };

        /**
         * @brief Contains the info about the closest target vertex for each source vertex.
         */
        moris::map< moris_index, Target_Info > mTargetMap;

        Target_Info& operator[]( moris_index aIndex )
        {
            return mTargetMap[ aIndex ];
        }

        Target_Info const & operator[]( moris_index aIndex ) const
        {
            return mTargetMap.at( aIndex );
        }

        /**
         * @brief Merges the result of two spatial indexing operations.
         * @param aOther
         * @return
         */
        Spatial_Indexing_Result merge( Spatial_Indexing_Result& aOther )
        {
            for ( auto& tPair : mTargetMap )
            {
                auto tVertex = tPair.first;
                auto tTarget = tPair.second;
                if ( aOther.mTargetMap.key_exists( tVertex ) )
                {
                    auto tOtherTarget = aOther.mTargetMap[ tVertex ];
                    if ( tTarget.distance > tOtherTarget.distance )
                    {
                        mTargetMap[ tVertex ] = tOtherTarget;
                    }
                }
            }

            // add the vertices that are not in the current map
            for ( auto& tPair : aOther.mTargetMap )
            {
                auto tVertex = tPair.first;
                auto tTarget = tPair.second;
                if ( !mTargetMap.key_exists( tVertex ) )
                {
                    mTargetMap[ tVertex ] = tTarget;
                }
            }
            return *this;
        }
    };

    class Spatial_Indexer
    {
      public:    // constructors
        Spatial_Indexer(
                moris::Cell< mtk::Surface_Mesh > const &                     aSurfaceMeshes,
                moris::Cell< std::pair< moris_index, moris_index > > const & aCandidatePairs )
                : mSurfaceMeshes( aSurfaceMeshes )
                , mCandidatePairs( aCandidatePairs )
        {
        }

        virtual ~Spatial_Indexer() = default;

      public:    // methods


        virtual Spatial_Indexing_Result perform( moris_index aSourceMeshIndex, real epsilon ) = 0;

      protected:    // methods
        /**
         * @brief List of surface mesh entities.
         */
        moris::Cell< mtk::Surface_Mesh > const mSurfaceMeshes;

        /**
         * @brief List of pairs that define which surface meshes should be compared. The first index defines the source mesh, the second index defines the target mesh.
         * To get a bidirectional comparison, the pair (i,j) and (j,i) must be added to the list.
         * @example If the list contains the pair (0,1), then the first surface mesh will be compared to the second one.
         * The pair (0, 0) means that the first surface mesh will be compared to itself (thus allowing a self-intersection).
         */
        moris::Cell< std::pair< moris_index, moris_index > > const mCandidatePairs;
    };


}    // namespace moris::mtk

#endif    // MORIS_CL_MTK_SPATIAL_INDEXER_H

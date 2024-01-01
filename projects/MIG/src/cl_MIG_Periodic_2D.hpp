/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MIG_Periodic_2D.hpp
 *
 */

#ifndef SRC_cl_MIG_Periodic_2D
#define SRC_cl_MIG_Periodic_2D

#include "cl_MTK_Mesh_Manager.hpp"
//#include "cl_MTK_Side_Cluster_ISC_Impl.hpp"
#include "cl_Param_List.hpp"
#include "cl_MTK_Intersec_Mesh_Data.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"

namespace moris::mig
{
    class Periodic_2D
    {
      private:
        std::shared_ptr< moris::mtk::Mesh_Manager > mMeshManager;
        moris::moris_index                          mMeshIndex;
        Vector< Vector< std::string > >   mMeshSideSetPairs;
        moris::uint                                 mNumBulkPhases;

        // Index of double sided cluster ( relevant indices of each cluster )
        // Each index shows a distinct interaction between leader phases and salve phases
        Vector< moris_index > mDoubleSidedClustersIndex;

        //the outer cell is the side cluster number
        //the inner cell vertex indices
        Vector< Vector< moris_index > > mSideClusterToVertexIndices ;

        //the outer cell is the side cluster number
        Vector< Vector< moris_index > > mSideClusterToCells;

        //the outer cell
        Vector< moris_index >                mSideClusterToIPCell;

        //cell to vertex connectivity
        //outer cell is the side cluster number
        Vector< Vector<moris_index> >  mCellToVertexIndices;

        Matrix< DDRMat > mVerticesCoords;
        Matrix< DDRMat > mVertexParametricCoords;

        uint mNumVertices = 0 ;
        uint mNumCells = 0 ;
        uint mNumSideClusters = 0;
        uint mNumDblSideCluster = 0;
        uint mNumParamCoords = 0 ;

        //maps, outer cell is the pair number
        //inside map, key: ijk of the IP cell, value: side clusters within the same background cell
        Vector< std::unordered_map< moris::moris_index, Vector< moris_index > > > mBackgroundCellToSideClusterMap1;
        Vector< std::unordered_map< moris::moris_index, Vector< moris_index > > > mBackgroundCellToSideClusterMap2;

      public:
        // ----------------------------------------------------------------------------
        Periodic_2D() = default ;
        /*
         *Default Constructor
         */
        Periodic_2D( std::shared_ptr< moris::mtk::Mesh_Manager > aMeshManager,
            moris::moris_index                                   aMeshIndex,
            moris::ParameterList                                &aParameterList,
            moris::uint                                          aNumBulkPhases );

        // ----------------------------------------------------------------------------
        /*
         *Default deconstructor
         */
        ~Periodic_2D();

        // ----------------------------------------------------------------------------
        /*
         * performs intersection of side set pairs and adds double sided sets to the mesh
         */
        void
        perform();

        // ----------------------------------------------------------------------------
        /**
         * @ makes new pairs of side cluster and associated double sided cluster
         * @param[in] tP cell containing vertices of lines created
         * @param[in] aRightInterpCell follower interpolation cell
         * @param[in] aLeftInterpCell leader interpolation cell
         * @param[in] aPairCount number of the pair in the periodic side set pair
         * @param[in] tPhaseToPhaseIndex a index showing interaction of leader-side and follower-side phases
         */
        void create_dbl_sided_cluster(
            Vector< Matrix< DDRMat > > tP,
            Vector< moris_index >     &aIndicesInCutCell,
            moris::mtk::Cell const         &aRightInterpCell,
            moris::mtk::Cell const         &aLeftInterpCell,
            uint                            aPairCount,
            moris_index                     aPhaseToPhase );

        // ----------------------------------------------------------------------------
        /**
         * @ creates a leader Integration cell
         * @param [ in ] tLeaderVertices all the leader vertices created from intersection of two side clusters
         * @param [ in ] aLeaderInterpCell interpolation cell of the leader side cluster
         * @param[in] aPairCount number of the pair in the periodic side set pair
         */
        void
        create_leader_ig_cell(
            moris::mtk::Cell const             &aLeaderInterpCell,
            uint                                aPairCount );

        // ----------------------------------------------------------------------------
        /**
         * @ creates a leader Integration cell
         * @param [ in ] tFollowerVertices all the leader vertices created from intersection of two side clusters
         * @param [ in ] aFollowerInterpCell interpolation cell of the leader side cluster
         * @param[in] aPairCount number of the pair in the periodic side set pair
         */
        void
        create_follower_ig_cell(
            moris::mtk::Cell const             &aFollowerInterpCell,
            uint                                aPairCount );

        // ----------------------------------------------------------------------------
        /*
         * constructs and add all the double sided sets to the mesh
         * @param[ in ] tPhaseInteractionTable a matrix specifying an index to each phase to phase interaction
         */
        void construct_add_dbl_sided_set( moris::Matrix< IndexMat > tPhaseInteractionTable );

        // ----------------------------------------------------------------------------
        /**
         * Calculate the offset vector between two surfaces and all the side sets attached to two surfaces
         * @param[ out ] tOffsetVector A row vector specifying offset for each pair
         * @param[ out ] tFirstSideSetNames Side set names of the first side in the pair
         * @param[ out ] tSecondSideSetNames Side set names of the second side in the pair
         * @param[ in ] aPairCount Number of the pair in the periodic side set pair
         */
        void
        offset_vector( moris::Matrix< DDRMat > &tOffsetVector,
            Vector< std::string >         &tFirstSideSetNames,
            Vector< std::string >         &tSecondSideSetNames,
            uint                                aPairCount );

        // ----------------------------------------------------------------------------
        /**
         * construct vertices based on the given coordinates
         * @param[ in ] tUniqueIntersectedPoints Unique coordinates obtained from intersecting two side clusters
         * @param[ in ] aLeaderInterpCell Leader interpolation cell
         * @param[ in ] aPairCount Number of the pair in the periodic side set pair
         */
        void
        create_leader_vertices( Matrix< DDRMat > tUniqueIntersectedPoints,
            moris::mtk::Cell const             &aFollowerInterpCell,
            uint                                aPairCount);

        // ----------------------------------------------------------------------------
        /**
         * construct vertices based on the given coordinates
         * @param[ in ] tUniqueIntersectedPoints Unique coordinates obtained from intersecting two side clusters
         * @param[ in ] aFollowerInterpCell Follower interpolation cell
         * @param[ in ] aPairCount Number of the pair in the periodic side set pair
         */
        void
        create_follower_vertices( Matrix< DDRMat > tUniqueIntersectedPoints,
            moris::mtk::Cell const             &aFollowerInterpCell,
            uint                                aPairCount );

        // ----------------------------------------------------------------------------
        /*
         * Gives the rotation matrix and inverse rotation matrix
         * InverseRotation matrix converts 3d coordinates to sudo-2d
         * Rotation matrix converts sudo-2d coordinates 3d
         * @param[ in ] aPairCount Number of the pair in the periodic side set pair
         * @param[ out ] aRotation Rotation matrix
         * @param[ out ] aInverseRotation Inverse rotation matrix
         */
        void
        rotation_matrix(
            moris::Matrix< DDRMat > &aRotation,
            moris::Matrix< DDRMat > &aInverseRotation,
            uint                     aPairCount );

        // ----------------------------------------------------------------------------
        /**
         * @brief this function generates a unique identifier( a number) based on the ijk position of the background cell
         *  All the clusters that are contained in the background cell will have the same identifier
         * @param[ in ]  aSideClusters all the clusters on 1 side
         * @param[ in ]  aPairCount Number of the pair in the periodic side set pair
         * @param[ out ] aMap map relating the background cells to the side clusters
         */
        void
        generate_identifier( Vector< mtk::Cluster const * >                 &aSideClusters,
            uint                                                                 &aPairCount,
            std::unordered_map< moris::moris_index, Vector< moris_index > > &aMap ) const;

        // ----------------------------------------------------------------------------
        /*
         * @ the main function in the intersection algorithm which takes 2 surfaces and return
         * @ their intersection
         * @param[ in ]  aFirstTRICoords coordinates of first triangle
         * @param[ in ]  aSecondTRICoords coordinates of second triangle
         * @param[ in ]  aFirstTRINodeIndex Local indices of the first TRI mesh
         * @param[ in ]  aSecondTRINodeIndex Local indices of the second TRI mesh
         */

        void
        elementwise_bruteforce_search(
            Vector< moris::Matrix< DDRMat > > const &tParamCoordsCell,
            moris::Matrix< moris::IndexMat > const       &tIGCellToSideClusterMap,
            Vector< moris::Matrix< DDRMat > > const &tParamCoordsCell2,
            moris::Matrix< moris::IndexMat > const       &tIGCellToSideClusterMap2,
            Vector< moris::Matrix< DDRMat > >       &tCutTriangles,
            moris::Matrix< moris::IndexMat >             &tCutTrianglesIdentifier ) const;

        // ----------------------------------------------------------------------------
        /*
         * computes intersection of two line segments
         * @param[ in ]  aFirstTRICoords coordinates of first triangle
         * @param[ in ]  aSecondTRICoords coordinates of second triangle
         * @param[ out ] aIntersectedPoints intersection points and
         * @param[ out ] aIntersectVecintersection vector showing which neighbors of first
         *  are also intersecting the second triangle
         */
        void
        Intersect(
            moris::Matrix< moris::DDRMat > const &aFirstTRICoords,
            moris::Matrix< moris::DDRMat > const &aSecondTRICoords,
            moris::Matrix< moris::DDRMat >       &aIntersectedPoints ) const;

        // ----------------------------------------------------------------------------

        /**
         * This function returns a pair specifying which indices should be converted from 3d coordinates to 2d coordinates
         * In order to generate surfaces.
         * @param[ in ] aPairCount Number of the pair in the periodic side set pair
         * @param[ out ] indices which we should use
         */
        void
        group_cut_cells( moris::Matrix< IndexMat > const                  &aCutCellIdentifierMatrix,
            std::unordered_map< moris_index, Vector< moris_index > > &aCutCellIdentifierToCutCellIndex ) const;

        // ----------------------------------------------------------------------------

        /**
         * This function returns a pair specifying which indices should be converted from 3d coordinates to 2d coordinates
         * In order to generate surfaces.
         * @param[ in ] aPairCount Number of the pair in the periodic side set pair
         * @param[ out ] indices which we should use
         */
        uint
        permutation_order( uint const &aPairCount ) const;

        // ----------------------------------------------------------------------------

        /**
         * @brief determine max size of the entities
         *
         */
        void
        determine_max_size();

        friend class Periodic_Mesh_Editor;
    };
}// namespace moris::mtk

#endif /* cl_MTK_Periodic_2D.hpp */

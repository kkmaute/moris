/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MIG_Periodic_3D.hpp
 *
 */

#ifndef SRC_cl_MIG_Periodic_3D
#define SRC_cl_MIG_Periodic_3D

#include "cl_MTK_Mesh_Manager.hpp"
// #include "cl_MTK_Side_Cluster_ISC_Impl.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_MTK_Intersec_Mesh_Data.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"

namespace moris::mig
{
    class Periodic_3D
    {
      private:
        std::shared_ptr< moris::mtk::Mesh_Manager > mMeshManager;
        moris::moris_index                          mMeshIndex;
        Vector< Vector< std::string > >             mMeshSideSetPairs;
        moris::uint                                 mNumBulkPhases;

        // Index of double sided cluster ( relevant indices of each cluster )
        // Each index shows a distinct interaction between leader phases and salve phases
        Vector< moris_index > mDoubleSidedClustersIndex;

        // the outer cell is the side cluster number
        // the inner cell vertex indices
        Vector< Vector< moris_index > > mSideClusterToVertexIndices;

        // the outer cell is the side cluster number
        Vector< Vector< moris_index > > mSideClusterToCells;

        // the outer cell
        Vector< moris_index > mSideClusterToIPCell;

        // cell to vertex connectivity
        // outer cell is the side cluster number
        Vector< Vector< moris_index > > mCellToVertexIndices;

        // coordinates and parametric coordinates of the added vertices
        Matrix< DDRMat > mVerticesCoords;
        Matrix< DDRMat > mVertexParametricCoords;

        // counter member data used to adjust the size at the end
        uint mNumVertices       = 0;
        uint mNumCells          = 0;
        uint mNumSideClusters   = 0;
        uint mNumDblSideCluster = 0;
        uint mNumParamCoords    = 0;

        // maps, outer cell is the pair number
        // inside map, key: ijk of the IP cell, value: side clusters within the same background cell
        Vector< std::unordered_map< moris::moris_index, Vector< moris_index > > > mBackgroundCellToSideClusterMap1;
        Vector< std::unordered_map< moris::moris_index, Vector< moris_index > > > mBackgroundCellToSideClusterMap2;

      public:
        // ----------------------------------------------------------------------------
        Periodic_3D() = default;
        /*
         *Default Constructor
         */
        Periodic_3D( std::shared_ptr< moris::mtk::Mesh_Manager > aMeshManager,
                moris::moris_index                               aMeshIndex,
                moris::Parameter_List                           &aParameterList,
                moris::uint                                      aNumBulkPhases );

        // ----------------------------------------------------------------------------
        /*
         *Default deconstructor
         */
        ~Periodic_3D();

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
                Vector< moris_index >     &aIndicesinCutCell,
                moris::mtk::Cell const    &aRightInterpCell,
                moris::mtk::Cell const    &aLeftInterpCell,
                uint                       aPairCount,
                moris_index                aPhaseToPhase );

        // ----------------------------------------------------------------------------
        /**
         * @ creates a leader Integration cell
         * @param [ in ] tLeaderVertices all the leader vertices created from intersection of two side clusters
         * @param [ in ] aLeaderInterpCell interpolation cell of the leader side cluster
         * @param[in] aPairCount number of the pair in the periodic side set pair
         */
        void
        create_leader_ig_cell( Matrix< IndexMat > const &tVertexIndex,
                moris::mtk::Cell const                  &aLeaderInterpCell,
                uint                                     aPairCount );

        // ----------------------------------------------------------------------------
        /**
         * @ creates a leader Integration cell
         * @param [ in ] tFollowerVertices all the leader vertices created from intersection of two side clusters
         * @param [ in ] aFollowerInterpCell interpolation cell of the leader side cluster
         * @param[in] aPairCount number of the pair in the periodic side set pair
         */
        void
        create_follower_ig_cell( Matrix< IndexMat > const &tVertexIndex,
                moris::mtk::Cell const                    &aFollowerInterpCell,
                uint                                       aPairCount );

        // ----------------------------------------------------------------------------
        /*
         * constrcuts and add all the double sided sets to the mesh
         * @param[ in ] tPhaseInteractionTable a matrix specifying an index to each phase to phase interaction
         */
        void construct_add_dbl_sided_set( const Matrix< IndexMat > &tPhaseInteractionTable );

        // ----------------------------------------------------------------------------
        /**
         * Calculate the offset vector between two surfaces and all the side sets attached to two surfaces
         * @param[ out ] tOffsetVector A row vector specifying offset for each pair
         * @param[ out ] tFirstSideSetNames Side set names of the first side in the pair
         * @param[ out ] tSecondSideSetNames Side set names of the second side in the pair
         * @param[ in ] aPairCount Number of the pair in the periodic side set pair
         */
        void
        offset_vector( Matrix< DDRMat > &tOffsetVector,
                Vector< std::string >   &tFirstSideSetNames,
                Vector< std::string >   &tSecondSideSetNames,
                uint                     aPairCount );

        // ----------------------------------------------------------------------------
        /**
         * construct vertices based on the given coordinates
         * @param[ in ] tUniqueIntersectedPoints Unique coordinates obtained from intersecting two side clusters
         * @param[ in ] aLeaderInterpCell Leader interpolation cell
         * @param[ in ] aPairCount Number of the pair in the periodic side set pair
         */
        void
        create_leader_vertices( Matrix< DDRMat > &tPhysicalCoordinates3D, moris::mtk::Cell const &aLeaderInterpCell, uint aPairCount, Matrix< DDRMat > &tNewNodeCoordinates );

        // ----------------------------------------------------------------------------
        /**
         * construct vertices based on the given coordinates
         * @param[ in ] aSurfaceNodesCoordinates Unique coordinates obtained from intersecting two side clusters
         * @param[ in ] aFollowerInterpCell Follower interpolation cell
         * @param[ in ] aPairCount Number of the pair in the periodic side set pair
         * @param[ in ] tTopNodeCoordinates head node of tetrahedran physical coordinates
         */
        void
        create_follower_vertices( Matrix< DDRMat > &tPhysicalCoordinates3D, moris::mtk::Cell const &aFollowerInterpCell, uint aPairCount, Matrix< DDRMat > &tNewNodeCoordinates );

        // ----------------------------------------------------------------------------
        /*
         * Gives the rotation matrix and inverse roation matrix
         * InverseRotation matrix converts 3d coordinates to sudo-2d
         * Rotation matrix converts sudo-2d coordinates 3d
         * @param[ in ] aPairCount Number of the pair in the periodic side set pair
         * @param[ out ] aRotation Rotation matrix
         * @param[ out ] aInverseRotation Inverse rotation matrix
         */
        void
        rotation_matrix(
                Matrix< DDRMat > &aRotation,
                Matrix< DDRMat > &aInverseRotation,
                uint              aPairCount );

        // ----------------------------------------------------------------------------
        /**
         * @brief this function generates a unique identifier( a number) based on the ijk position of the background cell
         *  All the clusters that are contained in the background cell will have the same identifier
         * @param[ in ]  aSideClusters all the clusters on 1 side
         * @param[ in ]  aPairCount Number of the pair in the periodic side set pair
         * @param[ out ] aMap map relating the background cells to the side clusters
         */
        void
        generate_identifier( Vector< mtk::Cluster const * >                     &aSideClusters,
                uint                                                            &aPairCount,
                std::unordered_map< moris::moris_index, Vector< moris_index > > &aBackgroundCellToSideClusterMap ) const;

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
                Vector< Matrix< DDRMat > > const &tParamCoordsCell,
                Vector< moris_index > const      &tIGCellToSideClusterMap,
                Vector< Matrix< DDRMat > > const &tParamCoordsCell2,
                Vector< moris_index > const      &tIGCellToSideClusterMap2,
                Vector< Matrix< DDRMat > >       &tCutTriangles,
                Vector< moris_index >            &tCutTrianglesIdentifier ) const;

        // ----------------------------------------------------------------------------

        /**
         * This function returns a pair specifying which indices should be converted from 3d coordinates to 2d coordinates
         * In order to generate surfaces.
         * @param[ in ] aPairCount Number of the pair in the periodic side set pair
         * @param[ out ] indices which we should use
         */
        void
        group_cut_cells( Vector< moris_index > const                     &aCutCellIndetiferMatrix,
                std::unordered_map< moris_index, Vector< moris_index > > &aCutCellIdentifierToCutCellIndex ) const;

        // ----------------------------------------------------------------------------

        /**
         * This function returns a pair specifying which indices should be converted from 3d coordinates to 2d coordinates
         * In order to generate surfaces.
         * @param[ in ] aPairCount Number of the pair in the periodic side set pair
         * @param[ out ] indices which we should use
         */
        std::pair< moris_id, moris_id >
        permutation_pair( uint const &aPairCount ) const;

        // ----------------------------------------------------------------------------

        /**
         * @brief determine max size of the entities
         *
         */
        void
        determine_max_size();

        friend class Periodic_Mesh_Editor;
    };
}    // namespace moris::mig

#endif /* cl_MTK_Periodic_3D.hpp */

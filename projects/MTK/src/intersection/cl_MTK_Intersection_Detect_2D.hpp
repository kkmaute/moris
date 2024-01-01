/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Intersection_Detect_2D.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_INTERSECTION_DETECT_2D_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_INTERSECTION_DETECT_2D_HPP_

#include "cl_MTK_Mesh_Manager.hpp"
//#include "cl_MTK_Side_Cluster_ISC_Impl.hpp"
#include "cl_Param_List.hpp"
#include "cl_MTK_Intersec_Mesh_Data.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"

namespace moris
{
    namespace mtk
    {
        class Intersection_Detect_2D
        {
            private:

                std::shared_ptr<moris::mtk::Mesh_Manager>             mMeshManager;
                moris::moris_index                                    mMeshIndex;
                Vector< Vector< std::string > >             mMeshSideSetPairs;
                moris::uint                                           mNumBulkPhases;

                // To keep track of id and index of added data
                Mesh_Intersection_Data                                mIntersectedMeshData;

                //All the double sided cluster
                Vector<moris::mtk::Cluster const* >              mDoubleSidedClusters;
                Vector<moris::mtk::Cluster const* >              mFollowerSidedClusters;
                Vector<moris::mtk::Cluster const* >              mLeaderSidedClusters;

                // All the leader side cells created in the intersection process
                Vector<moris::mtk::Cell const *>                         mLeaderSideCells;
                Vector<moris::mtk::Cell const *>                         mFollowerSideCells;

                // All the vertices created in the intersection process
                Vector<moris::mtk::Vertex const *>                         mLeaderVertices;
                Vector<moris::mtk::Vertex const *>                         mFollowerVertices;

                // Index of double sided cluster ( relevant indices of each cluster )
                // Each index shows a distinct interaction between leader phases and salve phases
                Vector<moris::moris_index >                      mDoubleSidedClustersIndex;

            public:

                // ----------------------------------------------------------------------------
                /*
                 *Default Constructor
                 */
                Intersection_Detect_2D( std::shared_ptr<moris::mtk::Mesh_Manager> aMeshManager,
                        moris::moris_index                        aMeshIndex,
                        moris::ParameterList &                    aParameterList,
                        moris::uint                               aNumBulkPhases);

                // ----------------------------------------------------------------------------
                /*
                 *Default deconstructor
                 */
                ~Intersection_Detect_2D();

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
                        Vector< Matrix < DDRMat> > tP,
                        Vector<moris_index> & aIndicesinCutCell,
                        moris::mtk::Cell const & aRightInterpCell,
                        moris::mtk::Cell const & aLeftInterpCell,
                        uint aPairCount,
                        moris_index aPhaseToPhase);

                // ----------------------------------------------------------------------------
                /**
                 * @ creates a leader Integration cell
                 * @param [ in ] tLeaderVertices all the leader vertices created from intersection of two side clusters
                 * @param [ in ] aLeaderInterpCell interpolation cell of the leader side cluster
                 * @param[in] aPairCount number of the pair in the periodic side set pair
                 */
                moris::mtk::Cell const * create_leader_ig_cell(
                        Vector<moris::mtk::Vertex *> tLeaderVertices ,
                        moris::mtk::Cell const & aLeaderInterpCell,
                        uint aPairCount);

                // ----------------------------------------------------------------------------
                /**
                 * @ creates a leader Integration cell
                 * @param [ in ] tFollowerVertices all the leader vertices created from intersection of two side clusters
                 * @param [ in ] aFollowerInterpCell interpolation cell of the leader side cluster
                 * @param[in] aPairCount number of the pair in the periodic side set pair
                 */
                moris::mtk::Cell const *  create_follower_ig_cell(
                        Vector<moris::mtk::Vertex *> tFollowerVertices ,
                        moris::mtk::Cell const & aFollowerInterpCell,
                        uint aPairCount );

                // ----------------------------------------------------------------------------
                /*
                 * constrcuts and add all the double sided sets to the mesh
                 * @param[ in ] tPhaseInteractionTable a matrix specifying an index to each phase to phase interaction
                 */
                void constrcuct_add_dbl_sided_set ( moris::Matrix < IndexMat > tPhaseInteractionTable);

                // ----------------------------------------------------------------------------
                /**
                 * Calculate the offset vector between two surfaces and all the side sets attached to two surfaces
                 * @param[ out ] tOffsetVector A row vector specifying offset for each pair
                 * @param[ out ] tFirstSideSetNames Side set names of the first side in the pair
                 * @param[ out ] tSecondSideSetNames Side set names of the second side in the pair
                 * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                 */
                void
                offset_vector(moris::Matrix<DDRMat > &tOffsetVector,
                        Vector< std::string > &tFirstSideSetNames,
                        Vector< std::string > &tSecondSideSetNames,
                        uint aPairCount);

                // ----------------------------------------------------------------------------
                /**
                 * construct vertices based on the given coordinates
                 * @param[ in ] tUniqueIntersectedPoints Unique coordinates obtained from intersecting two side clusters
                 * @param[ in ] aLeaderInterpCell Leader interpolation cell
                 * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                 */
                Vector< moris::mtk::Vertex *>
                create_leader_vertices( Matrix < DDRMat> tUniqueIntersectedPoints,
                        moris::mtk::Cell const & aLeaderInterpCell,
                        uint aPairCount);

                // ----------------------------------------------------------------------------
                /**
                 * construct vertices based on the given coordinates
                 * @param[ in ] tUniqueIntersectedPoints Unique coordinates obtained from intersecting two side clusters
                 * @param[ in ] aFollowerInterpCell Follower interpolation cell
                 * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                 */
                Vector< moris::mtk::Vertex *>
                create_follower_vertices( Matrix < DDRMat> tUniqueIntersectedPoints,
                        moris::mtk::Cell const & aFollowerInterpCell,
                        uint aPairCount);

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
                        moris::Matrix< DDRMat > & aRotation,
                        moris::Matrix< DDRMat > & aInverseRotation,
                        uint aPairCount);

                // ----------------------------------------------------------------------------
                /**
                 * @brief this function generates a unique identifier( a number) based on the ijk position of the background cell
                 *  All the clusters that are contained in the background cell will have the same identifier
                 * @param[ in ]  aSideClusters all the clusters on 1 side
                 * @param[ in ]  aPairCount Number of the pair in the periodic side set pair
                 * @param[ out ] aMap map relating the background cells to the side clusters
                 */
                void
                generate_identifier( Vector< mtk::Cluster const * > & aSideClusters,
                        uint &                                  aPairCount,
                        std::unordered_map< moris::moris_index, Vector<moris_index> > & aMap) const ;

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
                elementwise_bruteforce_search (
                        Vector < moris::Matrix <DDRMat> >  const &  tParamCoordsCell,
                        moris::Matrix< moris::IndexMat>  const &         tIGCellToSideClusterMap,
                        Vector < moris::Matrix <DDRMat> >  const &  tParamCoordsCell2,
                        moris::Matrix< moris::IndexMat>  const &         tIGCellToSideClusterMap2,
                        Vector < moris::Matrix <DDRMat> > &         tCutTriangles,
                        moris::Matrix< moris::IndexMat> &                tCutTrianglesIdentifier) const;

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
                        moris::Matrix < moris::DDRMat >  const &  aFirstTRICoords,
                        moris::Matrix < moris::DDRMat  > const &  aSecondTRICoords,
                        moris::Matrix < moris::DDRMat  > &        aIntersectedPoints) const;

                // ----------------------------------------------------------------------------

                /**
                 * This function returns a pair specifying which indices should be converted from 3d coordinates to 2d coordinates
                 * In order to generate surfaces.
                 * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                 * @param[ out ] indices which we should use
                 */
                void
                group_cut_cells(moris::Matrix< IndexMat> const & aCutCellIndetiferMatrix,
                        std::unordered_map< moris_index , Vector< moris_index > > & aCutCellIdentifierToCutCellIndex) const;

                // ----------------------------------------------------------------------------

                /**
                 * This function returns a pair specifying which indices should be converted from 3d coordinates to 2d coordinates
                 * In order to generate surfaces.
                 * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                 * @param[ out ] indices which we should use
                 */
                uint
                permutation_order(uint const & aPairCount) const;

        };
    }
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_INTERSECTION_DETECT_HPP_ */


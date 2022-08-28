/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Intersection_Detect.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_INTERSECTION_DETECT_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_INTERSECTION_DETECT_HPP_

#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Side_Cluster_ISC_Impl.hpp"
#include "cl_Param_List.hpp"
#include "cl_MTK_Intersec_Mesh_Data.hpp"

namespace moris
{
    namespace mtk
    {

        class Intersection_Detect
        {

            private:

                std::shared_ptr<moris::mtk::Mesh_Manager>             mMeshManager;
                moris::moris_index                                    mMeshIndex;
                moris::Cell< moris::Cell< std::string > >             mMeshSideSetPairs;
                moris::uint                                           mNumBulkPhases;

                // To keep track of id and index of added data
                Mesh_Intersection_Data                                mIntersectedMeshData;

                // All the double sided cluster
                moris::Cell<moris::mtk::Cluster const* >              mDoubleSidedClusters;
                moris::Cell<moris::mtk::Cluster const* >              mMasterSidedClusters;
                moris::Cell<moris::mtk::Cluster const* >              mSlaveSidedClusters;

                // Index of double sided cluster ( relevant indices of each cluster )
                // Each index shows a distinct interaction between master phases and salve phases
                moris::Cell<moris::moris_index >                      mDoubleSidedClustersIndex;

                // Double side sets
                moris::Cell<moris::Cell<mtk::Cluster const*> >        mDoubleSideSets;
                moris::Cell<moris::Cell<mtk::Cluster const*> >        mMasterSideSets;

                // Local to global Id Entity Matrix
                // The outer cell is the entity rank
                // The inner cell contains Ids
                moris::Cell<moris::Cell<moris::moris_index>>           mEntityLocaltoGlobalMap;
                //moris::Cell<std::unordered_map< moris::moris_index, moris::moris_id >> mEntityLocaltoGlobalMap;

                // Vertex constructed by the decomposition process
                std::unordered_map< moris_id, moris_index>             mVertexGlbToLocalMap;

                //coordinates of new nodes
                moris::Cell< Matrix<DDRMat> >                          mNewNodeCoords;

                // Side Set labels and a corresponding map
                std::unordered_map<std::string, moris_index>           mSideSideSetLabelToOrd;
                moris::Cell<std::string>                               mSideSetLabels;

                // block set labels and a corresponding map
                moris::Cell<std::string>                                mBlockSetLabels;
                std::unordered_map<std::string, moris_index>            mBlockSetLabelToOrd;

                // All the master side cells created in the intersection process
                moris::Cell<moris::mtk::Cell const *>                   mMasterSideCells;
                moris::Cell<moris::mtk::Cell const *>                   mSlaveSideCells;

                moris::Cell<moris::mtk::Vertex const *>                 mMasterVertices;
                moris::Cell<moris::mtk::Vertex const *>                 mSlaveVertices;

                // A map from cell index to index of the cell containing all the mtk::cells ( mMasterSideCells );
                std::unordered_map<moris_index, moris_index >           mMasterCellIndextoCellMap;

                //Visualization flag for the intersection mesh
                bool                                                    mVisFlag=true;

            public:

                // ----------------------------------------------------------------------------
                /*
                 * Default constructor
                 */
                Intersection_Detect( std::shared_ptr<moris::mtk::Mesh_Manager> aMeshManager,
                        moris::moris_index                                     aMeshIndex,
                        moris::ParameterList &                                 aParameterList,
                        moris::uint                                            aNumBulkPhases);

                // ----------------------------------------------------------------------------
                /*
                 * Default deconstructor
                 */
                ~Intersection_Detect();
                // ----------------------------------------------------------------------------
                /*
                 * the main performer which creates cluster and double sided side sets
                 */
                void
                perform();

                // ----------------------------------------------------------------------------
                /*
                 * computes the edge intersection of two triangles
                 * @param[ in ]  aFirstTRICoords coordinates of first triangle
                 * @param[ in ]  aSecondTRICoords coordinates of second triangle
                 * @param[ out ] aIntersectedPoints intersection points and
                 * @param[ out ] aIntersectVecintersection vector showing which neighbors of first
                 *  are also intersecting the second triangle
                 */
                void
                EdgeIntersect(
                        moris::Matrix < moris::DDRMat >  const &  aFirstTRICoords,
                        moris::Matrix < moris::DDRMat  > const &  aSecondTRICoords,
                        moris::Matrix < moris::DDRMat  > &        aIntersectedPoints,
                        moris::Matrix < moris::DDUMat  > &        aIntersectVec) const;

                // ----------------------------------------------------------------------------
                /*
                 * finds corners of one triangle within another one
                 * computes for the two given triangles
                 * @param[ in ]  aFirstTRICoords coordinates of first triangle
                 * @param[ in ]  aSecondTRICoords coordinates of second triangle
                 * @param[ out ] aIntersectedPoints intersection points and
                 * @note (point coordinates are stored column-wise, in counter clock
                 *order) the corners of first which lie in the interior of second.
                 */
                void  PointsXInY(
                        moris::Matrix < moris::DDRMat >  const &  aFirstTRICoords,
                        moris::Matrix < moris::DDRMat  > const &  aSecondTRICoords,
                        moris::Matrix < moris::DDRMat  > &        aIntersectedPoints) const;

                // ----------------------------------------------------------------------------
                /*
                 * sort points and remove duplicates
                 * orders polygon corners in counterclock wise and removes duplicates
                 * @param[ in ] aIntersectedPoints polygon points, not ordered
                 * @param[ out ] aIntersectedPoints polygon points, ordered
                 */
                void
                SortandRemove ( moris::Matrix < moris::DDRMat  > & aIntersectedPoints ) const;

                // ----------------------------------------------------------------------------
                /*
                 * All intersection points of two triangles
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
                        moris::Matrix < moris::DDRMat  > &        aIntersectedPoints,
                        moris::Matrix < moris::DDUMat  > &        aIntersectVec) const;

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
                        moris::Cell < moris::Matrix <DDRMat> >  const &  tParamCoordsCell,
                        moris::Matrix< moris::IndexMat>  const &         tIGCellToSideClusterMap,
                        moris::Cell < moris::Matrix <DDRMat> >  const &  tParamCoordsCell2,
                        moris::Matrix< moris::IndexMat>  const &         tIGCellToSideClusterMap2,
                        moris::Cell < moris::Matrix <DDRMat> > &         tCutTriangles,
                        moris::Matrix< moris::IndexMat> &                tCutTrianglesIdentifier) const;

                // ----------------------------------------------------------------------------
                /*
                 * @ makes new pairs of side cluster and associated double sided cluster
                 * @param[in] tP cell containing vertices of triangles/polygons created
                 * @param[in] aRightInterpCell slave interpolation cell
                 * @param[in] aLeftInterpCell master interpolation cell
                 * @param[in] aPairCount number of the pair in the periodic side set pair
                 * @param[in] tPhaseToPhaseIndex a index showing interaction of master-side and slave-side phases
                 */
                void
                create_dbl_sided_cluster(
                        moris::Cell< Matrix < DDRMat> > &   tP,
                        moris::Cell<moris_index> &          tIndicesinCutCell,
                        moris::mtk::Cell const &            aRightInterpCell,
                        moris::mtk::Cell const &            aLeftInterpCell,
                        uint &                              aPairCount,
                        moris_index  &                      aPhaseToPhaseIndex);

                // ----------------------------------------------------------------------------
                /**
                 * @ creates a master Integration cell
                 * @param [ in ] tMasterVertices all the master vertices created from intersection of two side clusters
                 * @param [ in ] tP indices of the tMasterVertices which create a cell
                 * @param [ in ] aMasterInterpCell interpolation cell of the master side clsuter
                 */
                moris::mtk::Cell const *
                create_master_ig_cell(moris::Cell<moris::mtk::Vertex *> aMasterVertices,
                        Matrix<IndexMat>                                aP ,
                        moris::mtk::Cell const &                        aMasterInterpCell,
                        uint                                            aPairCount);

                // ----------------------------------------------------------------------------
                /*
                 * @ creates a master Integration cell
                 * @param [ in ] tSlaveVertices all the master vertices created from intersection of two side clusters
                 * @param [ in ] tP indices of the tSlaveVertices which create a cell
                 * @param [ in ] aMasterInterpCell interpolation cell of the master side clsuter
                 */
                moris::mtk::Cell const *
                create_slave_ig_cell(moris::Cell<moris::mtk::Vertex *> aSlaveVertices,
                        Matrix<IndexMat>                               aP ,
                        moris::mtk::Cell const &                       aSlaveInterpCell,
                        uint                                           aPairCount );

                // ----------------------------------------------------------------------------
                /*
                 * @ construct and add double sided set
                 * @param[ in ] tPhaseInteractionTable a matrix specifying an index to each phase to phase interaction
                 */
                void
                construct_add_dbl_sided_set (moris::Matrix < IndexMat > const & aPhaseInteractionTable);

                // ----------------------------------------------------------------------------
                /**
                 * Calculate the offset vector between two surfaces and all the side sets attached to two surfaces
                 * @param[ out ] tOffsetVector A row vector specifying offset for each pair
                 * @param[ out ] tFirstSideSetNames Side set names of the first side in the pair
                 * @param[ out ] tSecondSideSetNames Side set names of the second side in the pair
                 * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                 */
                void
                offset_vector(
                        moris::Matrix<DDRMat > &      tOffsetVector,
                        moris::Cell< std::string > &  tFirstSideSetNames,
                        moris::Cell< std::string > &  tSecondSideSetNames,
                        uint                           aPairCount);

                // ----------------------------------------------------------------------------
                /**
                 * construct vertices based on the given coordinates
                 * @param[ in ] aSurfaceNodesCoordinates 3D physical coordinates on the surface
                 * @param[ in ] aMasterInterpCell Master interpolation cell
                 * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                 * @param[ in ] aTopNodeCoordinates head node of tetrahedran physical coordinates
                 */
                 void
                create_master_vertices(moris::Cell< moris::mtk::Vertex *> & aMasterVertices,
                        Matrix < DDRMat> &                                  aSurfaceNodesCoordinates,
                        moris::mtk::Cell const &                            aMasterInterpCell,
                        uint                                                aPairCount,
                        Matrix < DDRMat> &                                  aTopNodeCoordinates);

                // ----------------------------------------------------------------------------
                /**
                 * construct vertices based on the given coordinates
                 * @param[ in ] aSurfaceNodesCoordinates Unique coordinates obtained from intersecting two side clusters
                 * @param[ in ] aSlaveInterpCell Slave interpolation cell
                 * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                 * @param[ in ] tTopNodeCoordinates head node of tetrahedran physical coordinates
                 */
                 void
                 create_slave_vertices(moris::Cell< moris::mtk::Vertex *> & aSlaveVertices,
                         Matrix < DDRMat> &                                 aSurfaceNodesCoordinates,
                         moris::mtk::Cell const &                           aSlaveInterpCell,
                         uint                                               aPairCount,
                         Matrix < DDRMat> &                                 aTopNodeCoordinates);

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
                        uint                      aPairCount) const;

                // ----------------------------------------------------------------------------
                /*
                 * Builds the necessary data for intersection algorithm
                 * @param[ in ] aSideCluster Side Cluster
                 * @param[ in ] aInverseRotation Inverse rotation matrix
                 * @param[ out ] aUniqueParamCoords List of all coordinates
                 * @param[ out ] aLocalElemToNode Element to node number map
                 * @param[ out ] aElemToElemLocal Element to element connectivity
                 */
                void
                build_input_data_for_intersection(
                        mtk::Cluster const *                          aSideCluster,
                        moris::Matrix< DDRMat >   &                   aInverseRotation,
                        moris::Matrix<DDRMat>  &                      aUniqueParamCoords,
                        moris::Matrix<moris::IdMat> &                 aLocalElemToNode);

                // ----------------------------------------------------------------------------
                /**
                 * Check if two clusters align
                 * @param[ in ] aFirstSideCluster First side cluster
                 * @param[ in ] aSecondSideCluster Second side cluste
                 * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                 */
                bool
                clusters_align(
                        mtk::Cluster const *                          aFirstSideCluster,
                        mtk::Cluster const *                          aSecondSideCluster,
                        uint                                          aPaircount,
                        moris::Matrix< DDRMat > const &               aOffsetMatrix) const;

                // ----------------------------------------------------------------------------
                /**
                 * @ the main function in the intersection algorithm which takes 2 surfaces and return
                 * @ their intersection
                 * @param[ in ]  aFirstTRICoords coordinates of first triangle
                 * @param[ in ]  aSecondTRICoords coordinates of second triangle
                 * @param[ in ]  aFirstTRINodeIndex Local indices of the first TRI mesh
                 * @param[ in ]  aSecondTRINodeIndex Local indices of the second TRI mesh
                 * @param[ in ]  aFirstTRIConnect Connection map of the first side cluster
                 * @param[ in ]  aSecondTRIConnect Connection map of the second side cluster
                 */
                moris::Cell< moris::Matrix < moris::DDRMat > >
                elementwise_intersection_search(
                        moris::Matrix < moris::DDRMat >   const &  aFirstTRICoords,
                        moris::Matrix < moris::DDRMat >   const &  aSecondTRICoords,
                        moris::Matrix < moris::IdMat >    const &  aFirstTRINodeIndex,
                        moris::Matrix < moris::IdMat >    const &  aSecondTRINodeIndex,
                        moris::Matrix < moris::IndexMat > const &  aFirstTRIConnect,
                        moris::Matrix < moris::IndexMat > const &  aSecondTRIConnect);

                // ----------------------------------------------------------------------------
                /**
                 * @brief this function generates a unique identifier( a number) based on the ijk position of the background cell
                 *  All the clusters that are contained in the background cell will have the same identifier
                 * @param[ in ]  aSideClusters all the clusters on 1 side
                 * @param[ in ]  aPairCount Number of the pair in the periodic side set pair
                 * @param[ out ] aIdentifierMatrix an identifying number for every cluster calculated based on the background cell
                 * @param[ out ] aMap map relating the background cells to the side clusters
                 */
                void
                generate_identifier( moris::Cell< mtk::Cluster const * > & aSideClusters,
                        uint &                                  aPairCount,
                        Matrix<IndexMat>               &        aIdentifierMatrix,
                        std::unordered_map< moris::moris_index, moris::Cell<moris_index> > & aMap) const ;

                // ----------------------------------------------------------------------------

                /**
                 * This function returns a pair specifying which indices should be converted from 3d coordinates to 2d coordinates
                 * In order to generate surfaces.
                 * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                 * @param[ out ] indices which we should use
                 */
                std::pair< moris_id, moris_id >
                permutation_pair(uint const & aPairCount) const;

                 // ----------------------------------------------------------------------------

                /**
                 * This function returns a pair specifying which indices should be converted from 3d coordinates to 2d coordinates
                 * In order to generate surfaces.
                 * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                 * @param[ out ] indices which we should use
                 */
                void
                group_cut_cells(moris::Matrix< IndexMat> const & aCutCellIndetiferMatrix,
                std::unordered_map< moris_index , moris::Cell< moris_index > > & aCutCellIdentifierToCutCellIndex) const;

                // ----------------------------------------------------------------------------

                /**
                 * This function returns a pair specifying which indices should be converted from 3d coordinates to 2d coordinates
                 * In order to generate surfaces.
                 * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                 * @param[ out ] indices which we should use
                 */

                void set_vis_flag(bool aVisRequested)
                {
                mVisFlag =aVisRequested ;
                }

                friend class Intersection_Mesh;
        };
    }
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_INTERSECTION_DETECT_HPP_ */


/*
 * cl_MTK_Intersection_Detect.hpp
 *
 *  Created on: Jun 7, 2021
 *      Author: momo
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
                moris::Cell< moris::Cell< std::string > >             mMeshSideSetPairs;
                moris::uint                                           mNumBulkPhases;

                // To keep track of id and index of added data
                Mesh_Intersection_Data                                mIntersectedMeshData;

                //All the double sided cluster
                moris::Cell<moris::mtk::Cluster const* >              mDoubleSidedClusters;
                moris::Cell<moris::mtk::Cluster const* >              mSlaveSidedClusters;
                moris::Cell<moris::mtk::Cluster const* >              mMasterSidedClusters;

                // Index of double sided cluster ( relevant indices of each cluster )
                // Each index shows a distinct interaction between master phases and salve phases
                moris::Cell<moris::moris_index >                      mDoubleSidedClustersIndex;

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
                 * @param[in] aRightInterpCell slave interpolation cell
                 * @param[in] aLeftInterpCell master interpolation cell
                 * @param[in] aPairCount number of the pair in the periodic side set pair
                 * @param[in] tPhaseToPhaseIndex a index showing interaction of master-side and slave-side phases
                 */
                void make_new_pairs(
                        moris::Cell< Matrix < DDRMat> > tP,
                        moris::mtk::Cell const & aRightInterpCell,
                        moris::mtk::Cell const & aLeftInterpCell,
                        uint aPairCount,
                        moris_index aPhaseToPhase);

                // ----------------------------------------------------------------------------
                /**
                 * @ creates a master Integration cell
                 * @param [ in ] tMasterVertices all the master vertices created from intersection of two side clusters
                 * @param [ in ] aMasterInterpCell interpolation cell of the master side cluster
                 * @param[in] aPairCount number of the pair in the periodic side set pair
                 */
                moris::mtk::Cell const * create_master_ig_cell(
                        moris::Cell<moris::mtk::Vertex *> tMasterVertices ,
                        moris::mtk::Cell const & aMasterInterpCell,
                        uint aPairCount);

                // ----------------------------------------------------------------------------
                /**
                 * @ creates a master Integration cell
                 * @param [ in ] tSlaveVertices all the master vertices created from intersection of two side clusters
                 * @param [ in ] aSlaveInterpCell interpolation cell of the master side cluster
                 * @param[in] aPairCount number of the pair in the periodic side set pair
                 */
                moris::mtk::Cell const *  create_slave_ig_cell(
                        moris::Cell<moris::mtk::Vertex *> tSlaveVertices ,
                        moris::mtk::Cell const & aSlaveInterpCell,
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
                        moris::Cell< std::string > &tFirstSideSetNames,
                        moris::Cell< std::string > &tSecondSideSetNames,
                        uint aPairCount);

                // ----------------------------------------------------------------------------
                /**
                 * construct vertices based on the given coordinates
                 * @param[ in ] tUniqueIntersectedPoints Unique coordinates obtained from intersecting two side clusters
                 * @param[ in ] aMasterInterpCell Master interpolation cell
                 * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                 */
                moris::Cell< moris::mtk::Vertex *>
                create_master_vertices( Matrix < DDRMat> tUniqueIntersectedPoints,
                        moris::mtk::Cell const & aMasterInterpCell,
                        uint aPairCount);

                // ----------------------------------------------------------------------------
                /**
                  * construct vertices based on the given coordinates
                  * @param[ in ] tUniqueIntersectedPoints Unique coordinates obtained from intersecting two side clusters
                  * @param[ in ] aSlaveInterpCell Slave interpolation cell
                  * @param[ in ] aPairCount Number of the pair in the periodic side set pair
                  */
                moris::Cell< moris::mtk::Vertex *>
                create_slave_vertices( Matrix < DDRMat> tUniqueIntersectedPoints,
                        moris::mtk::Cell const & aSlaveInterpCell,
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

        };
    }
}



#endif /* PROJECTS_MTK_SRC_CL_MTK_INTERSECTION_DETECT_HPP_ */

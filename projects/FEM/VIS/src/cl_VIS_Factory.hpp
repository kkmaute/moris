/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_VIS_Factory.hpp
 *
 */

#ifndef SRC_FEM_CL_VIS_FACTORY_HPP_
#define SRC_FEM_CL_VIS_FACTORY_HPP_

#include "cl_Vector.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#include "cl_VIS_Vertex_Visualization.hpp"
#include "cl_VIS_Cell_Visualization.hpp"
#include "cl_VIS_Cell_Cluster_Visualization.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Output_Enums.hpp"

namespace moris
{
    namespace mtk
    {
        class Set;
        class Cluster;
        class Cell;
        class Vertex;
        class Interpolation_Mesh;
        class Integration_Mesh;
    }

    namespace vis
    {
        struct Output_Data;

        class VIS_Factory
        {
            //-----------------------------------------------------------------------------------------------------------

          private:
            //-----------------------------------------------------------------------------------------------------------

            /// @brief VIS-mesh to be constructed by this factory
            vis::Visualization_Mesh* mVisMesh = nullptr;

            /// @brief Names of sets requested for output
            Vector< std::string > mAllRequestedSetNames;
            Vector< std::string > mRequestedBlockSetNames;
            Vector< std::string > mRequestedSideSetNames;
            Vector< std::string > mRequestedDoubleSideSetNames;

            /// @brief mtk/fem mesh sets corresponding to the sets listed above
            Vector< moris::mtk::Set* > mFemBlockSets;
            Vector< moris::mtk::Set* > mFemSideSets;
            Vector< moris::mtk::Set* > mFemDoubleSideSets;

            /// @brief map relating the FEM/MTK cells in each block to the VIS cells (to be) created
            // || input: (1) index of blockset in vis mesh (2) index of cell in the FEM/MTK IG mesh 
            // || output: index of VIS cell created from this FEM/MTK cell 
            Vector< Vector< moris_index > > mBlockAndFemCellIndexToVisCellIndex;

            /// @brief map relating the FEM/MTK cells in each block to the VIS vertices (to be) created
            // || input: (1) index of blockset in vis mesh (2) index of cell in the FEM/MTK IG mesh 
            // || output: list of indices of VIS vertices in element local order
            Vector< Vector< Vector< moris_index > > > mBlockAndFemCellIndexToVisVertexIndices;

            /// @brief map relating the fem cell index to the vis cell index
            // || input: index of IG cell from FEM mesh 
            // || output: index of the corresponding VIS cell which is primary material
            // (Note: for overlapping meshes there may be multiple VIS cells created on the same FEM/MTK cell 
            // but only one is primary wrt. to one of the material phases)
            Vector< moris_index > mPrimaryFemCellIndexToVisCellIndex;

            /// @brief map relating the fem cell index to the vis cell index
            // || input: index of IG cell from FEM mesh 
            // || output: index of the VIS block set in which this cell sits
            // (Note: for overlapping meshes there may be multiple VIS cells created on the same FEM/MTK cell 
            // but only one is primary wrt. to one of the material phases)
            Vector< moris_index > mPrimaryFemCellIndexToBlockIndex;

            /// @brief which cells to output depending on VIS mesh type used
            bool mOnlyPrimaryCells = false;

            /// @brief flag indicating whether the the mesh is to be constructed using fully discontinuous elements (each element gets their own nodes)
            bool mConstructDiscontinuousMesh = false;

            /// @brief access to the mtk::mesh
            mtk::Interpolation_Mesh*             mInterpolationMesh = nullptr;
            mtk::Integration_Mesh*               mIntegrationMesh   = nullptr;

          public:

            /**
             * Constructor
             *
             * @param aMesh Mesh manager
             * @param aMeshPairIndex Mesh pair index
             */
            VIS_Factory(
                    std::shared_ptr< mtk::Mesh_Manager > aMesh,
                    uint                                 aMeshPairIndex );

            /**
             * Destructor, checks for memory leak
             */
            ~VIS_Factory();

            /**
             * @brief initialize member data and sort requested mesh sets by type
             *
             * @param aOutputData input parameters specifying output request
             */
            void
            initialize( moris::vis::Output_Data const & aOutputData );

            //-----------------------------------------------------------------------------------------------------------

            vis::Visualization_Mesh*
            hand_off_VIS_mesh()
            {
                MORIS_ERROR( mVisMesh != nullptr,
                        "VIS_Factory::hand_off_VIS_mesh() - "
                        "VIS mesh is nullptr. The mesh has either already been handed off to another owner, or not been constructed." );

                MORIS_ERROR( mVisMesh->mMeshIsFinalized,
                        "VIS_Factory::hand_off_VIS_mesh() - "
                        "The VIS mesh has not been finalized yet and is not ready to be handed off." );

                // copy the pointer for passing off
                vis::Visualization_Mesh* tVisMeshPtr = mVisMesh;

                // delete pointer from member data
                mVisMesh = nullptr;

                // pass on pointer to VIS mesh
                return tVisMeshPtr;
            }

            //-----------------------------------------------------------------------------------------------------------

            mtk::Mesh*
            create_visualization_mesh( moris::vis::Output_Data& aOutputData );

            //-----------------------------------------------------------------------------------------------------------

            /**
             * @brief Generate the vertices on the VIS mesh
             *
             */
            void
            create_visualization_vertices();

            //-----------------------------------------------------------------------------------------------------------
            
            /**
             * @brief Generate the vertices assuming continuous fields within blocks
             */
            void
            create_visualization_vertices_standard();

            //-----------------------------------------------------------------------------------------------------------
            
            /**
             * @brief Generate the vertices such that fields discontinuous between clusters can be computed
             */
            void
            create_visualization_vertices_full_discontinuous();

            //-----------------------------------------------------------------------------------------------------------

            /**
             * @brief Generated the IG cells on the VIS mesh
             *
             */
            void
            create_visualization_cells();

            //-----------------------------------------------------------------------------------------------------------

            /**
             * @brief Generate the clusters in the VIS mesh to interface with FEM for evaluation
             *
             */
            void
            create_visualization_clusters();

            //-----------------------------------------------------------------------------------------------------------

            /**
             * @brief Generate the side clusters in the VIS mesh to interface with FEM for evaluation
             *
             */
            void
            create_visualization_side_clusters();

            //-----------------------------------------------------------------------------------------------------------

            /**
             * @brief Generate the double sided side clusters in the VIS mesh to interface with FEM for evaluation
             *
             */
            void
            create_visualization_double_side_clusters();

            //-----------------------------------------------------------------------------------------------------------

            /**
             * @brief Generate the VIS block sets
             *
             */
            void
            create_visualization_blocks();

            //-----------------------------------------------------------------------------------------------------------

            /**
             * @brief generate the VIS side sets
             *
             */
            void
            create_visualization_side_sets();

            //-----------------------------------------------------------------------------------------------------------

            /**
             * @brief generate the VIS double sided side sets
             *
             */
            void
            create_visualization_double_side_sets();

            //-----------------------------------------------------------------------------------------------------------
        };
    }    // namespace vis
} /* namespace moris */

#endif /* SRC_FEM_CL_VIS_FACTORY_HPP_ */

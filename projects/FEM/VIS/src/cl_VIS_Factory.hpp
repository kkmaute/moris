/*
 * cl_VIS_Factory.hpp
 *
 *  Created on: Dez 02, 2019
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_VIS_FACTORY_HPP_
#define SRC_FEM_CL_VIS_FACTORY_HPP_

#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#include "../../../MRS/CON/src/cl_Cell.hpp"
#include "../../../MTK/src/cl_MTK_Cell.hpp"
#include "../../../MTK/src/cl_MTK_Integration_Mesh.hpp"
#include "../../../MTK/src/cl_MTK_Mesh_Manager.hpp"       //MTK/src
#include "../../../MTK/src/cl_MTK_Set.hpp"
#include "../../../MTK/src/cl_MTK_Vertex.hpp"

#include "cl_VIS_Vertex_Visualization.hpp"
#include "cl_VSI_Cell_Visualization.hpp"
#include "cl_VIS_Cell_Cluster_Visualization.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"

namespace moris
{
    namespace vis
    {
        class Factory
        {
        private:
            moris::Cell< moris::mtk::Set * > mListofBlocks;

            moris::Cell< moris::Cell< mtk::Cluster * > >        mClustersOnBlock;   //FIXME delete can be used temporary
            moris::Cell< moris::Cell< const mtk::Cluster * > >  mClustersOnBlock_1;
            moris::Cell< moris::Cell< mtk::Cell * > >   mCellsOnBlock;
            moris::Cell< moris::Cell< mtk::Vertex * > > mVerticesOnBlock;

            moris::Cell< Matrix< DDSMat > >             mVertexMapOnBlock;
            moris::Cell< Matrix< DDSMat > >             mCellMapOnBlock;

            uint mRequestedBlocks;
			bool mOnlyPrimary = false;

            mtk::Mesh_Manager* mMesh = nullptr;
            const uint         mMeshPairIndex;

            Matrix< IndexMat > mListOfRequestedBlocks;


        public:
            Factory( mtk::Mesh_Manager* aMesh,
                     const uint         aMeshPairIndex ) : mMesh( aMesh ),
                                                           mMeshPairIndex( aMeshPairIndex )
            {
                mListOfRequestedBlocks = { { 0, 1, 2, 3} };
				mOnlyPrimary = false;

                mRequestedBlocks = mListOfRequestedBlocks.numel();

                this->create_visualization_vertices();

                this->create_visualization_cells();

                this->create_visualization_clusters();

                this->create_visualization_blocks();
            };

//-----------------------------------------------------------------------------------------------------------

            ~Factory(){};

//-----------------------------------------------------------------------------------------------------------

            void create_visualization_vertices()
            {
                mVerticesOnBlock.resize( mRequestedBlocks );
                mVertexMapOnBlock.resize( mRequestedBlocks );

                mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
                mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
                mMesh->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

                moris_index tVertexIndexCounter = 0;
                moris_id tVertexIdCounter = 0;

                for( uint Ij = 0; Ij <mRequestedBlocks; Ij++ )
                {
                    moris::mtk::Set * tMeshSet = tIntegrationMesh->get_block_by_index( mListOfRequestedBlocks( Ij ) );

                    uint tNumVerticesOnSet = tMeshSet->get_num_vertieces_on_set( mOnlyPrimary );

                    moris::Matrix< DDSMat > tVertexIndOnBlock = tMeshSet->get_vertieces_inds_on_block( mOnlyPrimary );

                    mVerticesOnBlock( Ij ).resize( tNumVerticesOnSet, nullptr );
                    mVertexMapOnBlock( Ij ).set_size( tVertexIndOnBlock.max() + 1, 1, -1 );

                    for( uint Ik = 0; Ik < tNumVerticesOnSet; Ik++ )
                    {
                        mVertexMapOnBlock( Ij )( tVertexIndOnBlock( Ik ) ) = Ik;

                        mVerticesOnBlock( Ij )( Ik ) = new Vertex_Visualization( tVertexIdCounter,
                                                                                 tVertexIndexCounter,
                                                                                 &tIntegrationMesh->get_mtk_vertex( tVertexIndOnBlock( Ik ) ) );

                        tVertexIndexCounter++;
                        tVertexIdCounter++;
                    }
                }
            };

//-----------------------------------------------------------------------------------------------------------

            void create_visualization_cells()
            {
                mCellsOnBlock  .resize( mRequestedBlocks );
                mCellMapOnBlock.resize( mRequestedBlocks );

                mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
                mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
                mMesh->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

                moris_index tCellIndexCounter = 0;
                moris_id tCellIdCounter = 0;

                for( uint Ij = 0; Ij <mRequestedBlocks; Ij++ )
                {
                    moris::mtk::Set * tMeshSet = tIntegrationMesh->get_block_by_index( mListOfRequestedBlocks( Ij ) );

                    uint tNumCellsOnSet = tMeshSet->get_num_cells_on_set( mOnlyPrimary );

                    moris::Matrix< DDSMat > tCellIndOnBlock = tMeshSet->get_cell_inds_on_block( mOnlyPrimary );

                    mCellsOnBlock( Ij ).resize( tNumCellsOnSet, nullptr );

                    mCellMapOnBlock( Ij ).set_size( tCellIndOnBlock.max() + 1, 1, -1 );

                    for( uint Ik = 0; Ik < tNumCellsOnSet; Ik++ )
                    {
                        Matrix< IndexMat > tIntegrationVertexInd = tIntegrationMesh->get_mtk_cell( tCellIndOnBlock( Ik ) ).get_vertex_inds();

                        moris::Cell< mtk::Vertex* > tCellVerices( tIntegrationVertexInd.numel(), nullptr );

                        for( uint Ii = 0; Ii < tIntegrationVertexInd.numel(); Ii++ )
                        {
                            tCellVerices( Ii ) = mVerticesOnBlock( Ij )( mVertexMapOnBlock( Ij )( tIntegrationVertexInd( Ii ) ) );
                        }

                        mCellMapOnBlock( Ij )( tCellIndOnBlock( Ik ) ) = Ik;

                        mCellsOnBlock( Ij )( Ik ) = new Cell_Visualization( tCellIdCounter,
                                                                            tCellIndexCounter,
                                                                            tCellVerices,
                                                                            &tIntegrationMesh->get_mtk_cell( tCellIndOnBlock( Ik ) ) );

                        tCellIndexCounter++;
                        tCellIdCounter++;
                    }
                }
            }

//-----------------------------------------------------------------------------------------------------------

            void create_visualization_clusters()
            {
                mClustersOnBlock.resize( mRequestedBlocks );
                mClustersOnBlock_1.resize( mRequestedBlocks );

                mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
                mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
                mMesh->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

                for( uint Ij = 0; Ij <mRequestedBlocks; Ij++ )
                {
                    moris::mtk::Set * tMeshSet = tIntegrationMesh->get_block_by_index( mListOfRequestedBlocks( Ij ) );

                    uint tNumClustersOnSet = tMeshSet->get_num_clusters_on_set();

                    moris::Cell< mtk::Cluster const * > tClustersOnBlock = tMeshSet->get_clusters_on_set();

                    mClustersOnBlock  ( Ij ).resize( tNumClustersOnSet, nullptr );
                    mClustersOnBlock_1( Ij ).resize( tNumClustersOnSet, nullptr );

                    for( uint Ik = 0; Ik < tNumClustersOnSet; Ik++ )
                    {
                        moris::Cell<moris::mtk::Cell const *> tPrimaryCells = tClustersOnBlock( Ik )->get_primary_cells_in_cluster();

                        moris::Cell< mtk::Cell * > tClusterPrimaryCells( tPrimaryCells.size(), nullptr );

                        for( uint Ii = 0; Ii < tPrimaryCells.size(); Ii++ )
                        {
                            tClusterPrimaryCells( Ii ) = mCellsOnBlock( Ij )( mCellMapOnBlock( Ij )( tPrimaryCells( Ii )->get_index() ) );
                        }

                        moris::Cell< moris::mtk::Cell const * > tVoidCells = tClustersOnBlock( Ik )->get_void_cells_in_cluster();

                        if( !mOnlyPrimary )
                        {
                            moris::Cell< mtk::Cell * > tClusterVoidCells( tVoidCells.size(), nullptr );

                            for( uint Ii = 0; Ii < tVoidCells.size(); Ii++ )
                            {
                                tClusterVoidCells( Ii ) = mCellsOnBlock( Ij )( mCellMapOnBlock( Ij )( tVoidCells( Ii )->get_index() ) );
                             }
                        }

                        mClustersOnBlock( Ij )( Ik ) = new Cell_Cluster_Visualization;


                        if( tVoidCells.size() > 0 && !mOnlyPrimary )
                        {
                             mClustersOnBlock( Ij )( Ik )->mark_as_nontrivial();

                             mClustersOnBlock( Ij )( Ik )->add_void_integration_cell( tVoidCells );
                        }

                        mClustersOnBlock( Ij )( Ik )->add_primary_integration_cell( tPrimaryCells );
                        mClustersOnBlock( Ij )( Ik )->set_interpolation_cell( &tClustersOnBlock( Ik )->get_interpolation_cell() );
                    }

                    for( uint Ik = 0; Ik < tNumClustersOnSet; Ik++ )
                    {
                        mClustersOnBlock_1( Ij )( Ik ) = mClustersOnBlock( Ij )( Ik );
                    }
                }
            }

//-----------------------------------------------------------------------------------------------------------

            void create_visualization_blocks()
            {
                mListofBlocks.resize( mRequestedBlocks );

                mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
                mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
                mMesh->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

                for( uint Ij = 0; Ij <mRequestedBlocks; Ij++ )
                {
                    moris::mtk::Set * tMeshSet = tIntegrationMesh->get_block_by_index( mListOfRequestedBlocks( Ij ) );

                    mListofBlocks( Ij ) = new moris::mtk::Block( tMeshSet->get_set_name(),
                                                                 mClustersOnBlock_1( Ij ),
                                                                 tMeshSet->get_spatial_dim() );
                }
            };

//-----------------------------------------------------------------------------------------------------------

            mtk::Mesh * create_visualization_mesh()
            {
                mtk::Mesh * tMesh = new Visualization_Mesh( mListofBlocks, mCellsOnBlock, mVerticesOnBlock, mOnlyPrimary );

                return tMesh;
            };
        };
    } /* namespace VIS */
} /* namespace moris */

#endif /* SRC_FEM_CL_VIS_FACTORY_HPP_ */

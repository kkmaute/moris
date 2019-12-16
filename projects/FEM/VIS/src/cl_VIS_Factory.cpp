
#include "cl_VIS_Output_Manager.hpp"

#include "cl_VIS_Factory.hpp"

#include "../../../MTK/src/cl_MTK_Cell.hpp"
#include "../../../MTK/src/cl_MTK_Integration_Mesh.hpp"
#include "../../../MTK/src/cl_MTK_Mesh_Manager.hpp"       //MTK/src
#include "../../../MTK/src/cl_MTK_Set.hpp"
#include "../../../MTK/src/cl_MTK_Vertex.hpp"

extern moris::Comm_Manager gMorisComm;

namespace moris
{
    namespace vis
    {

    mtk::Mesh * VIS_Factory::create_visualization_mesh( moris::vis::Output_Data & aOutputData )
    {
        switch( aOutputData.mMeshType )
        {
            case ( vis::VIS_Mesh_Type::STANDARD ):
                mOnlyPrimaryCells = true ;
                break;

            case ( vis::VIS_Mesh_Type::OVERLAPPING_INTERFACE ):
                mOnlyPrimaryCells = false;
                break;

            case ( vis::VIS_Mesh_Type::FULL_DISCONTINOUS ):
                 MORIS_ERROR( false, "create_visualization_mesh() - Mesh type FULL_DISCONTINOUS not implemented yet. " );
                 break;

            default:
                MORIS_ERROR( false, "create_visualization_mesh() - Mesh type not specified. " );
                break;
        }

        tRequestedSetNames = aOutputData.mSetNames;

        mNumRequestedBlocks = tRequestedSetNames.size();

        this->create_visualization_vertices();

        this->create_visualization_cells();

        this->create_visualization_clusters();

        this->create_visualization_blocks();

        return new Visualization_Mesh( mListofBlocks, mCellsOnBlock, mVerticesOnBlock, mOnlyPrimaryCells );
    }

//-----------------------------------------------------------------------------------------------------------

    void VIS_Factory::create_visualization_vertices()
    {
        mVerticesOnBlock.resize( mNumRequestedBlocks );
        mVertexMapOnBlock.resize( mNumRequestedBlocks );

        mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
        mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
        mMesh->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

        moris_index tVertexIndexCounter = 0;
        moris_id tVertexIdCounter = 0;

        for( uint Ij = 0; Ij < mNumRequestedBlocks; Ij++ )
        {
            moris_index tSetIndex = tIntegrationMesh->get_block_set_index( tRequestedSetNames( Ij ) );

            moris::mtk::Set * tMeshSet = tIntegrationMesh->get_block_by_index( tSetIndex );

            uint tNumVerticesOnSet = tMeshSet->get_num_vertieces_on_set( mOnlyPrimaryCells );

            moris::Matrix< DDSMat > tVertexIndOnBlock = tMeshSet->get_vertieces_inds_on_block( mOnlyPrimaryCells );

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
    }

//-----------------------------------------------------------------------------------------------------------

    void VIS_Factory::create_visualization_cells()
    {
        mCellsOnBlock  .resize( mNumRequestedBlocks );
        mCellMapOnBlock.resize( mNumRequestedBlocks );

        mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
        mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
        mMesh->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

        moris_index tCellIndexCounter = 0;
        moris_id    tCellIdCounter    = 0;

        for( uint Ij = 0; Ij <mNumRequestedBlocks; Ij++ )
        {
            moris_index tSetIndex = tIntegrationMesh->get_block_set_index( tRequestedSetNames( Ij ) );
            moris::mtk::Set * tMeshSet = tIntegrationMesh->get_block_by_index( tSetIndex );

            uint tNumCellsOnSet = tMeshSet->get_num_cells_on_set( mOnlyPrimaryCells );

            moris::Matrix< DDSMat > tCellIndOnBlock = tMeshSet->get_cell_inds_on_block( mOnlyPrimaryCells );

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

    void VIS_Factory::create_visualization_clusters()
    {
        mClustersOnBlock.resize( mNumRequestedBlocks );
        mClustersOnBlock_1.resize( mNumRequestedBlocks );

        mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
        mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
        mMesh->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

        for( uint Ij = 0; Ij <mNumRequestedBlocks; Ij++ )
        {
            moris_index tSetIndex = tIntegrationMesh->get_block_set_index( tRequestedSetNames( Ij ) );
            moris::mtk::Set * tMeshSet = tIntegrationMesh->get_block_by_index( tSetIndex );

            uint tNumClustersOnSet = tMeshSet->get_num_clusters_on_set();

            moris::Cell< mtk::Cluster const * > tClustersOnBlock = tMeshSet->get_clusters_on_set();

            mClustersOnBlock  ( Ij ).resize( tNumClustersOnSet, nullptr );
            mClustersOnBlock_1( Ij ).resize( tNumClustersOnSet, nullptr );

            for( uint Ik = 0; Ik < tNumClustersOnSet; Ik++ )
            {
                moris::Cell<moris::mtk::Cell const *> tPrimaryCells = tClustersOnBlock( Ik )->get_primary_cells_in_cluster();

                moris::Cell< moris::mtk::Cell const * > tClusterPrimaryCells( tPrimaryCells.size(), nullptr );

                for( uint Ii = 0; Ii < tPrimaryCells.size(); Ii++ )
                {
                    tClusterPrimaryCells( Ii ) = mCellsOnBlock( Ij )( mCellMapOnBlock( Ij )( tPrimaryCells( Ii )->get_index() ) );
                }

                moris::Cell< moris::mtk::Cell const * > tVoidCells = tClustersOnBlock( Ik )->get_void_cells_in_cluster();

                moris::Cell< moris::mtk::Cell const * > tClusterVoidCells( tVoidCells.size(), nullptr );

                if( !mOnlyPrimaryCells )
                {
                    for( uint Ii = 0; Ii < tVoidCells.size(); Ii++ )
                    {
                        tClusterVoidCells( Ii ) = mCellsOnBlock( Ij )( mCellMapOnBlock( Ij )( tVoidCells( Ii )->get_index() ) );
                     }
                }

                mClustersOnBlock( Ij )( Ik ) = new Cell_Cluster_Visualization;


                if( tVoidCells.size() > 0 && !mOnlyPrimaryCells )
                {
                     mClustersOnBlock( Ij )( Ik )->mark_as_nontrivial();

                     mClustersOnBlock( Ij )( Ik )->add_void_integration_cell( tClusterVoidCells );
                }

                mClustersOnBlock( Ij )( Ik )->add_primary_integration_cell( tClusterPrimaryCells );

                mClustersOnBlock( Ij )( Ik )->set_interpolation_cell( &tClustersOnBlock( Ik )->get_interpolation_cell() );
            }

            for( uint Ik = 0; Ik < tNumClustersOnSet; Ik++ )
            {
                mClustersOnBlock_1( Ij )( Ik ) = mClustersOnBlock( Ij )( Ik );
            }
        }
    }

//-----------------------------------------------------------------------------------------------------------

    void VIS_Factory::create_visualization_blocks()
    {
        mListofBlocks.resize( mNumRequestedBlocks );

        mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
        mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
        mMesh->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

        for( uint Ij = 0; Ij <mNumRequestedBlocks; Ij++ )
        {
            moris_index tSetIndex = tIntegrationMesh->get_block_set_index( tRequestedSetNames( Ij ) );
            moris::mtk::Set * tMeshSet = tIntegrationMesh->get_block_by_index( tSetIndex );

            mListofBlocks( Ij ) = new moris::mtk::Block( tMeshSet->get_set_name(),
                                                         mClustersOnBlock_1( Ij ),
                                                         tMeshSet->get_spatial_dim() );
        }
    }

//-----------------------------------------------------------------------------------------------------------


    }
}


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
            Tracer tTracer( "VisFactory", "VisMesh", "CreateVisMesh" );

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

            // Get set names
            tRequestedSetNames = aOutputData.mSetNames;

            // get num requested sets
            mNumRequestedSets = tRequestedSetNames.size();

            this->create_visualization_vertices();

            this->create_visualization_cells();

            this->create_visualization_clusters( aOutputData);

            this->create_visualization_blocks();

            // Create vis mesh
            return new Visualization_Mesh( mListofBlocks,mClustersOnBlock, mCellsOnSet, mVerticesOnSet, mOnlyPrimaryCells );
        }

        //-----------------------------------------------------------------------------------------------------------

        void VIS_Factory::create_visualization_vertices()
        {
            Tracer tTracer( "VisFactory", "VisMesh", "CreateVertices" );

            // resize list of vertices on set
            mVerticesOnSet.resize( mNumRequestedSets );
            mVertexMapOnSet.resize( mNumRequestedSets );

            // get integration mesh
            mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
            mMesh->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

            // define vertex id and index counter
            moris_index tVertexIndexCounter = 0;
            moris_id tVertexIdCounter = 0;

            uint tMaxVertexInd = tIntegrationMesh->get_num_entities( EntityRank::NODE );

            // loop over all requested sets
            for( uint Ij = 0; Ij < mNumRequestedSets; Ij++ )
            {
                // get set index for label
                moris_index tSetIndex = tIntegrationMesh->get_set_index_by_name( tRequestedSetNames( Ij ) );

                // get set for index
                moris::mtk::Set * tMeshSet = tIntegrationMesh->get_set_by_index( tSetIndex );

                // get all vertices on set
                uint tNumVerticesOnSet = tMeshSet->get_num_vertices_on_set( mOnlyPrimaryCells );

                if( tNumVerticesOnSet > 0 )
                {
                    // get vertex indices on set
                    moris::Matrix< DDSMat > tVertexIndOnBlock = tMeshSet->get_ig_vertices_inds_on_block( mOnlyPrimaryCells );

                    // resize list of vertices for this set
                    mVerticesOnSet( Ij ).resize( tNumVerticesOnSet, nullptr );
                    mVertexMapOnSet( Ij ).set_size( tMaxVertexInd + 1, 1, -1 );

                    // Loop over all vertices on this set and create new vis vertices
                    for( uint Ik = 0; Ik < tNumVerticesOnSet; Ik++ )
                    {
                        // create old index to set local index map
                        mVertexMapOnSet( Ij )( tVertexIndOnBlock( Ik ) ) = Ik;

                        // create vis vertex and renumber id and index
                        mVerticesOnSet( Ij )( Ik ) = new Vertex_Visualization( tVertexIdCounter++,
                                tVertexIndexCounter++,
                                &tIntegrationMesh->get_mtk_vertex( tVertexIndOnBlock( Ik ) ) );
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void VIS_Factory::create_visualization_cells()
        {
            Tracer tTracer( "VisFactory", "VisMesh", "CreateCells" );

            // resize list of cells on set
            mCellsOnSet  .resize( mNumRequestedSets );
            mCellMapOnSet.resize( mNumRequestedSets );

            // get integration mesh
            mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
            mMesh->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

            // define cell id and index counter
            moris_index tCellIndexCounter = 0;
            moris_id    tCellIdCounter    = 0;

            // loop over all requested sets
            for( uint Ij = 0; Ij <mNumRequestedSets; Ij++ )
            {
                // get set index for label
                moris_index tSetIndex = tIntegrationMesh->get_set_index_by_name( tRequestedSetNames( Ij ) );

                // get set for index
                moris::mtk::Set * tMeshSet = tIntegrationMesh->get_set_by_index( tSetIndex );

                // get number of cells on set
                uint tNumCellsOnSet = tMeshSet->get_num_cells_on_set( mOnlyPrimaryCells );

                if( tNumCellsOnSet > 0 )
                {
                    // get cell indices on set
                    moris::Matrix< DDSMat > tCellIndOnSet = tMeshSet->get_cell_inds_on_block( mOnlyPrimaryCells );

                    // resize list of cells for this set
                    mCellsOnSet( Ij ).resize( tNumCellsOnSet, nullptr );
                    mCellMapOnSet( Ij ).set_size( tCellIndOnSet.max() + 1, 1, -1 );

                    // Loop over all cells on this set and create new vis cells
                    for( uint Ik = 0; Ik < tNumCellsOnSet; Ik++ )
                    {
                        // get list of integration vertex indices
                        Matrix< IndexMat > tIntegrationVertexInd = tIntegrationMesh->get_mtk_cell( tCellIndOnSet( Ik ) ).get_vertex_inds();

                        // Create List of vis vertex pointers for this cell
                        moris::Cell< mtk::Vertex* > tCellVerices( tIntegrationVertexInd.numel(), nullptr );

                        // loop over integration vertices and get the corresponding vis vertices for this set
                        for( uint Ii = 0; Ii < tIntegrationVertexInd.numel(); Ii++ )
                        {
                            tCellVerices( Ii ) = mVerticesOnSet( Ij )( mVertexMapOnSet( Ij )( tIntegrationVertexInd( Ii ) ) );
                        }

                        // create old index to set local index map
                        mCellMapOnSet( Ij )( tCellIndOnSet( Ik ) ) = Ik;

                        // create vis cells and renumber id and index
                        mCellsOnSet( Ij )( Ik ) = new Cell_Visualization( tCellIdCounter++,
                                tCellIndexCounter++,
                                tCellVerices,
                                &tIntegrationMesh->get_mtk_cell( tCellIndOnSet( Ik ) ) );
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void VIS_Factory::create_visualization_clusters( moris::vis::Output_Data & aOutputData )
        {
            Tracer tTracer( "VisFactory", "VisMesh", "CreateClusters" );

            // resize list of cells on set
            mClustersOnBlock.resize( mNumRequestedSets );

            // get integration mesh
            mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
            mMesh->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

            // loop over requested sets
            for( uint Ij = 0; Ij <mNumRequestedSets; Ij++ )
            {
                // get set index for label
                moris_index tSetIndex = tIntegrationMesh->get_set_index_by_name( tRequestedSetNames( Ij ) );

                // get set for index
                moris::mtk::Set * tMeshSet = tIntegrationMesh->get_set_by_index( tSetIndex );

                // get number of clusters on set
                uint tNumClustersOnSet = tMeshSet->get_num_clusters_on_set();

                if( tNumClustersOnSet > 0 )
                {
                    // get list of old clusters on set
                    moris::Cell< mtk::Cluster const * > tClustersOnSet = tMeshSet->get_clusters_on_set();

                    // resize list of clustersfor this set
                    mClustersOnBlock  ( Ij ).resize( tNumClustersOnSet, nullptr );

                    Matrix< DDSMat > tIdentifierMat( mVertexMapOnSet( Ij ).numel(), 1, -1);

                    // loop over clusters on set
                    for( uint Ik = 0; Ik < tNumClustersOnSet; Ik++ )
                    {
                        // get primary/void cells on old cluster
                        const moris::Cell< moris::mtk::Cell const * > & tPrimaryCells = tClustersOnSet( Ik )->get_primary_cells_in_cluster();
                        const moris::Cell< moris::mtk::Cell const * > & tVoidCells    = tClustersOnSet( Ik )->get_void_cells_in_cluster();

                        // resize primary/void cell list for vis cluster
                        moris::Cell< moris::mtk::Cell const * > tClusterPrimaryCells( tPrimaryCells.size(), nullptr );
                        moris::Cell< moris::mtk::Cell const * > tClusterVoidCells   ( tVoidCells.size()   , nullptr );

                        // find vis primary cells corresponding to integration primary cells
                        for( uint Ii = 0; Ii < tPrimaryCells.size(); Ii++ )
                        {
                            tClusterPrimaryCells( Ii ) = mCellsOnSet( Ij )( mCellMapOnSet( Ij )( tPrimaryCells( Ii )->get_index() ) );
                        }

                        // find vis void cells corresponding to integration void cells
                        if( !mOnlyPrimaryCells )
                        {
                            for( uint Ii = 0; Ii < tVoidCells.size(); Ii++ )
                            {
                                tClusterVoidCells( Ii ) = mCellsOnSet( Ij )( mCellMapOnSet( Ij )( tVoidCells( Ii )->get_index() ) );
                            }
                        }

                        // create vis cluster
                        Cell_Cluster_Visualization * tVisCellCluster = new Cell_Cluster_Visualization;

                        // add vis primary cells to vis cluster
                        tVisCellCluster->add_primary_integration_cell( tClusterPrimaryCells );

                        // add interpolation cell to vis cluster
                        tVisCellCluster->set_interpolation_cell( &tClustersOnSet( Ik )->get_interpolation_cell() );

                        // mark as non trivial if old cluster was trivial
                        if( !tClustersOnSet( Ik )->is_trivial() )
                        {
                            tVisCellCluster->mark_as_nontrivial();
                        }

                        // add void cells to cluster if requested by mesh type
                        if( tVoidCells.size() > 0 && !mOnlyPrimaryCells )
                        {
                            tVisCellCluster->add_void_integration_cell( tClusterVoidCells );
                        }

                        // add vertices and local coordinates to vis cluster if non-trivial
                        // if( !tVisCellCluster->is_trivial() )
                        // {
                        // get vertices from old cluster
                        moris::Cell<moris::mtk::Vertex const *> tAllVertices = tClustersOnSet( Ik )->get_vertices_in_cluster();
                        moris::Cell<moris::mtk::Vertex const *> tVertices;

                        switch( aOutputData.mMeshType )
                        {
                            case ( vis::VIS_Mesh_Type::STANDARD ):
                            {
                                 moris::Cell<moris::mtk::Vertex *> tPrimVertices = tClustersOnSet( Ik )->get_primary_vertices_in_cluster();

                                 tVertices.resize(tPrimVertices.size());

                                 for( uint Ib = 0; Ib < tPrimVertices.size(); Ib ++ )
                                 {
                                      tVertices( Ib ) = tPrimVertices( Ib );
                                 }
                            }

                            break;

                            case ( vis::VIS_Mesh_Type::OVERLAPPING_INTERFACE ):
                                tVertices = tClustersOnSet( Ik )->get_vertices_in_cluster();
                            break;

                            case ( vis::VIS_Mesh_Type::FULL_DISCONTINOUS ):
                                  MORIS_ERROR( false, "create_visualization_mesh() - Mesh type FULL_DISCONTINOUS not implemented yet. " );
                            break;

                            default:
                                MORIS_ERROR( false, "create_visualization_mesh() - Mesh type not specified. " );
                                break;
                        }

                        moris::Cell< mtk::Vertex const * > tVisClusterVertices( tVertices.size() );

                        // identifier and pos map
                        tIdentifierMat.fill( -1 );

                        uint tCounter = 0;

                        for( uint Ii = 0; Ii < tVertices.size(); Ii++ )
                        {
                            // get old vertex index
                            moris_index tIndex = tVertices( Ii )->get_index();

                            if( tIdentifierMat(tIndex) == -1)
                            {
                                tIdentifierMat(tIndex) = tCounter;

                                if( mVertexMapOnSet( Ij )( tIndex ) != -1 )
                                {
                                    tVisClusterVertices( tCounter++ ) = mVerticesOnSet( Ij )( mVertexMapOnSet( Ij )( tIndex ) );
                                }
                            }
                        }

                        tVisClusterVertices.resize( tCounter );

                        tVisCellCluster->add_vertex_to_cluster( tVisClusterVertices );

                        uint tNumCols = tClustersOnSet( Ik )->get_vertices_local_coordinates_wrt_interp_cell().n_cols();

                        moris::Matrix< DDRMat > tVisClusterVerticesLocalCoords( tCounter, tNumCols, 0.0 );

                        tCounter = 0;

                        for( uint Ii = 0; Ii < tAllVertices.size(); Ii++ )
                        {
                            // get old vertex index
                            moris_index tIndex = tAllVertices( Ii )->get_index();

                            if( tIdentifierMat(tIndex) != -1)
                            {
                                sint tPos = tIdentifierMat(tIndex);
                                tIdentifierMat(tIndex) = -1;

                                if( mVertexMapOnSet( Ij )( tIndex ) != -1 )
                                {
                                    tVisClusterVerticesLocalCoords.set_row( tPos, tClustersOnSet( Ik )->get_vertices_local_coordinates_wrt_interp_cell().get_row( Ii ) );
                                }
                            }
                        }

                        // add local coordinates to vis cluster
                        tVisCellCluster->add_vertex_local_coordinates_wrt_interp_cell( tVisClusterVerticesLocalCoords );
                        // }

                        mClustersOnBlock( Ij )( Ik ) = tVisCellCluster;
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void VIS_Factory::create_visualization_blocks()
        {
            Tracer tTracer( "VisFactory", "VisMesh", "CreateBlocks" );

            mListofBlocks.resize( mNumRequestedSets );

            mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;

            mMesh->get_mesh_pair( mMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

            for( uint Ij = 0; Ij <mNumRequestedSets; Ij++ )
            {
                moris_index tSetIndex = tIntegrationMesh->get_set_index_by_name( tRequestedSetNames( Ij ) );

                moris::mtk::Set * tMeshSet = tIntegrationMesh->get_set_by_index( tSetIndex );

                mListofBlocks( Ij ) = new moris::mtk::Block( tMeshSet->get_set_name(),
                        mClustersOnBlock( Ij ),
                        tMeshSet->get_set_colors(),
                        tMeshSet->get_spatial_dim() );

                mListofBlocks( Ij )->set_cell_topology( tMeshSet->get_cell_topology() );
                mListofBlocks( Ij )->set_IG_cell_shape( tMeshSet->get_IG_cell_shape() );
            }
        }

        //-----------------------------------------------------------------------------------------------------------
    }
}

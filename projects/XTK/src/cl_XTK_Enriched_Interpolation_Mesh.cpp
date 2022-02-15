/*
 * cl_XTK_Enriched_Interpolation_Mesh.cpp
 *
 *  Created on: Jul 10, 2019
 *      Author: doble
 */

#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_Map.hpp"

#include "cl_XTK_Multigrid.hpp"
#include "fn_TOL_Capacities.hpp"

namespace xtk
{
    // ----------------------------------------------------------------------------

    Enriched_Interpolation_Mesh::Enriched_Interpolation_Mesh(
            Model* aXTKModel )
            : mXTKModel( aXTKModel )
            , mNumVerts( 0 )
            , mNumVertsPerInterpCell( MORIS_UINT_MAX )
            , mCellInfo( nullptr )
    {
    }

    // ----------------------------------------------------------------------------

    Enriched_Interpolation_Mesh::~Enriched_Interpolation_Mesh()
    {
        for ( moris::uint i = 0; i < mInterpVertEnrichment.size(); i++ )
        {
            for ( auto it : mInterpVertEnrichment( i ) )
            {
                delete it;
            }

            mInterpVertEnrichment( i ).clear();
        }

        mInterpVertEnrichment.clear();

        for ( auto it : mEnrichedInterpCells )
        {
            delete it;
        }

        mEnrichedInterpCells.clear();

        for ( auto it : mEnrichedInterpVerts )
        {
            delete it;
        }

        mEnrichedInterpVerts.clear();
    }

    // ----------------------------------------------------------------------------

    MeshType
    Enriched_Interpolation_Mesh::get_mesh_type() const
    {
        return MeshType::XTK;
    }

    // ----------------------------------------------------------------------------

    moris::uint
    Enriched_Interpolation_Mesh::get_spatial_dim() const
    {
        return mXTKModel->get_spatial_dim();
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_num_entities(
            enum EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        switch ( aEntityRank )
        {
            case EntityRank::NODE:
            {
                return mEnrichedInterpVerts.size();
                break;
            }
            case EntityRank::ELEMENT:
            {
                return mEnrichedInterpCells.size();
                break;
            }
            default:
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
                return 0;
        }
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_max_num_coeffs_on_proc( const uint aBSplineMeshIndex ) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index( aBSplineMeshIndex );

        return mEnrichCoeffLocToGlob( tLocalMeshIndex ).numel();
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Interpolation_Mesh::get_entity_connected_to_entity_loc_inds(
            moris_index       aEntityIndex,
            enum EntityRank   aInputEntityRank,
            enum EntityRank   aOutputEntityRank,
            const moris_index aIndex ) const
    {
        MORIS_ERROR( aInputEntityRank == EntityRank::ELEMENT && aOutputEntityRank == EntityRank::NODE,
                "Only support element to node connectivity" );

        MORIS_ASSERT( aEntityIndex < (moris_index)mEnrichedInterpCells.size(),
                "Element index out of bounds" );

        return mEnrichedInterpCells( aEntityIndex )->get_vertex_inds();
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Interpolation_Mesh::get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const
    {
        MORIS_ERROR( 0, "XTK ENRICHED MESH ERROR: get_elements_connected_to_element_and_face_ind_loc_inds no implemented" );

        return Matrix< IndexMat >( 0, 0 );
    }

    // ----------------------------------------------------------------------------

    Cell< mtk::Vertex const * >
    Enriched_Interpolation_Mesh::get_all_vertices() const
    {
        moris::uint tNumNodes = this->get_num_entities( EntityRank::NODE );

        Cell< mtk::Vertex const * > tVertices( tNumNodes );

        for ( moris::uint i = 0; i < tNumNodes; i++ )
        {
            tVertices( i ) = mEnrichedInterpVerts( i );
        }

        return tVertices;
    }

    // ----------------------------------------------------------------------------

    moris_id
    Enriched_Interpolation_Mesh::get_glb_entity_id_from_entity_loc_index(
            moris_index       aEntityIndex,
            enum EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        moris::uint tMapIndex = (uint)aEntityRank;

        MORIS_ASSERT( aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT,
                "XTK ENRICHED MESH ERROR: Only support glb to loc conversion for vertices and cells" );

        MORIS_ASSERT( aEntityIndex < (moris_index)mLocalToGlobalMaps( tMapIndex ).numel(),
                "Entityindex out of bounds" );

        return mLocalToGlobalMaps( tMapIndex )( aEntityIndex );
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_loc_entity_ind_from_entity_glb_id(
            moris_id          aEntityId,
            enum EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        moris::uint tMapIndex = (uint)aEntityRank;

        auto tIter = mGlobaltoLobalMaps( tMapIndex ).find( aEntityId );

        MORIS_ASSERT( aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT,
                "XTK ENRICHED MESH ERROR: Only support glb to loc conversion for vertices and cells" );

        if ( tIter == mGlobaltoLobalMaps( tMapIndex ).end() )
        {
            std::cout << "Not Found  Entity Id = " << aEntityId << " | par_rank = " << par_rank() << std::endl;
        }

        MORIS_ASSERT( tIter != mGlobaltoLobalMaps( tMapIndex ).end(), "Id does not appear in map" );

        return tIter->second;
    }

    // ----------------------------------------------------------------------------

    std::unordered_map< moris_id, moris_index >
    Enriched_Interpolation_Mesh::get_vertex_glb_id_to_loc_vertex_ind_map() const
    {
        moris::uint tMapIndex = (uint)EntityRank::NODE;

        return mGlobaltoLobalMaps( tMapIndex );
    }

    // ----------------------------------------------------------------------------

    Matrix< IdMat >
    Enriched_Interpolation_Mesh::get_entity_connected_to_entity_glob_ids(
            moris_id          aEntityId,
            enum EntityRank   aInputEntityRank,
            enum EntityRank   aOutputEntityRank,
            const moris_index aIndex ) const
    {
        moris_index tEntityIndex =
                get_loc_entity_ind_from_entity_glb_id(
                        aEntityId,
                        aInputEntityRank );

        Matrix< IndexMat > tEntityToEntityLoc =
                this->get_entity_connected_to_entity_loc_inds(
                        tEntityIndex,
                        aInputEntityRank,
                        aOutputEntityRank );

        return convert_indices_to_ids( tEntityToEntityLoc, aOutputEntityRank );
    }

    // ----------------------------------------------------------------------------

    Matrix< DDRMat >
    Enriched_Interpolation_Mesh::get_node_coordinate( moris_index aNodeIndex ) const
    {
        mtk::Vertex const & tVertex = get_mtk_vertex( aNodeIndex );

        return tVertex.get_coords();
    }

    // ----------------------------------------------------------------------------

    Matrix< DDRMat >
    Enriched_Interpolation_Mesh::get_base_node_coordinate( moris_index aBaseNodeIndex ) const
    {
        return mXTKModel->get_background_mesh().get_mtk_vertex( aBaseNodeIndex ).get_coords();
    }

    // ----------------------------------------------------------------------------

    mtk::Vertex&
    Enriched_Interpolation_Mesh::get_mtk_vertex( moris_index aVertexIndex )
    {
        MORIS_ASSERT( aVertexIndex < (moris_index)mEnrichedInterpVerts.size(),
                "Vertex index out of bounds" );

        return *mEnrichedInterpVerts( aVertexIndex );
    }

    // ----------------------------------------------------------------------------

    mtk::Vertex const &
    Enriched_Interpolation_Mesh::get_mtk_vertex( moris_index aVertexIndex ) const
    {
        MORIS_ASSERT( aVertexIndex < (moris_index)mEnrichedInterpVerts.size(),
                "Vertex index out of bounds" );

        return *mEnrichedInterpVerts( aVertexIndex );
    }

    // ----------------------------------------------------------------------------

    mtk::Cell const &
    Enriched_Interpolation_Mesh::get_mtk_cell( moris_index aElementIndex ) const
    {
        MORIS_ASSERT( aElementIndex < (moris_index)mEnrichedInterpCells.size(),
                "Cell index out of bounds" );

        return *mEnrichedInterpCells( aElementIndex );
    }

    // ----------------------------------------------------------------------------------

    mtk::Cell&
    Enriched_Interpolation_Mesh::get_mtk_cell( moris_index aElementIndex )
    {
        MORIS_ASSERT( aElementIndex < (moris_index)mEnrichedInterpCells.size(),
                "Cell index out of bounds" );

        return *mEnrichedInterpCells( aElementIndex );
    }
    // ----------------------------------------------------------------------------

    mtk::Cell&
    Enriched_Interpolation_Mesh::get_writable_mtk_cell( moris_index aElementIndex )
    {
        MORIS_ASSERT( aElementIndex < (moris_index)mEnrichedInterpCells.size(),
                "Cell index out of bounds" );

        return *mEnrichedInterpCells( aElementIndex );
    }

    // ----------------------------------------------------------------------------

    Matrix< IdMat >
    Enriched_Interpolation_Mesh::get_communication_table() const
    {
        return mXTKModel->get_cut_integration_mesh()->get_communication_table();
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_num_elements()
    {
        return mEnrichedInterpCells.size();
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Interpolation_Mesh::get_element_indices_in_block_set( uint aSetIndex )
    {
        Matrix< IndexMat > tElementIndices( mEnrichedInterpCells.size(), 1 );
        if ( aSetIndex == 0 )
        {
            for ( uint tElementIndex = 0; tElementIndex < mEnrichedInterpCells.size(); tElementIndex++ )
            {
                tElementIndices( tElementIndex ) = tElementIndex;
            }
        }
        else
        {
            tElementIndices.set_size( 0, 0 );
        }

        return tElementIndices;
    }

    // ----------------------------------------------------------------------------

    enum CellTopology
    Enriched_Interpolation_Mesh::get_blockset_topology( const std::string& aSetName )
    {
        uint tNumberOfDimensions = this->get_spatial_dim();

        moris::mtk::Interpolation_Order tOrder = mEnrichedInterpCells( 0 )->get_interpolation_order();

        enum CellTopology tCellTopology = CellTopology::END_ENUM;

        switch ( tNumberOfDimensions )
        {
            case 1:
            {
                MORIS_ERROR( false, "1D not implemented" );
                break;
            }
            case 2:
            {
                switch ( tOrder )
                {
                    case moris::mtk::Interpolation_Order::LINEAR:
                    {
                        tCellTopology = CellTopology::QUAD4;
                        break;
                    }
                    case moris::mtk::Interpolation_Order::QUADRATIC:
                    {
                        tCellTopology = CellTopology::QUAD9;
                        break;
                    }
                    case moris::mtk::Interpolation_Order::CUBIC:
                    {
                        tCellTopology = CellTopology::QUAD16;
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, " Order not implemented" );
                        break;
                    }
                }
                break;
            }
            case 3:
            {
                switch ( tOrder )
                {
                    case moris::mtk::Interpolation_Order::LINEAR:
                    {
                        tCellTopology = CellTopology::HEX8;
                        break;
                    }
                    case moris::mtk::Interpolation_Order::QUADRATIC:
                    {
                        tCellTopology = CellTopology::HEX27;
                        break;
                    }
                    case moris::mtk::Interpolation_Order::CUBIC:
                    {
                        tCellTopology = CellTopology::HEX64;
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, " Order not implemented" );
                        break;
                    }
                }
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Number of dimensions not implemented" );
                break;
            }
        }
        return tCellTopology;
    }

    // ----------------------------------------------------------------------------

    enum CellShape
    Enriched_Interpolation_Mesh::get_IG_blockset_shape( const std::string& aSetName )
    {
        // get the clusters in the set
        moris::Cell< mtk::Cluster const * > tSetClusters = this->get_set_by_name( aSetName )->get_clusters_on_set();

        // init cell shape
        CellShape tCellShape = CellShape::EMPTY;

        // if the set isn't empty exist
        if ( tSetClusters.size() > 0 )
        {
            // get the cells in the first cluster
            moris::Cell< moris::mtk::Cell const * > tClusterCells = tSetClusters( 0 )->get_primary_cells_in_cluster();

            // compute the cell shape of the first cell
            tCellShape = tClusterCells( 0 )->get_cell_info()->compute_cell_shape( tClusterCells( 0 ) );
        }

        // within debug, checking all cells to make sure that they are the same Cell Shape
        // if cells exist
        // looping through the clusters
        for ( uint iCluster = 0; iCluster < tSetClusters.size(); iCluster++ )
        {
            // get cell of cells in the cluster
            moris::Cell< moris::mtk::Cell const * > tClusterCellsCheck = tSetClusters( iCluster )->get_primary_cells_in_cluster();

            // looping through the cells in the cluster
            for ( uint iCheckCell = 0; iCheckCell < tClusterCellsCheck.size(); iCheckCell++ )
            {
                MORIS_ASSERT( tClusterCellsCheck( iCheckCell )->get_cell_info()->compute_cell_shape( tClusterCellsCheck( iCheckCell ) ) == tCellShape,
                        "Mesh_Core_STK::get_IG_blockset_shape - cell shape is not consistent in the block" );
            }
        }

        return tCellShape;
    }

    // ----------------------------------------------------------------------------

    enum CellShape
    Enriched_Interpolation_Mesh::get_IP_blockset_shape( const std::string& aSetName )
    {
        // get the clusters in the set
        moris::Cell< mtk::Cluster const * > tSetClusters = this->get_set_by_name( aSetName )->get_clusters_on_set();

        // init cell shape
        CellShape tCellShape = CellShape::EMPTY;

        // if the set isn't empty exist
        if ( tSetClusters.size() > 0 )
        {
            // get the cells in the first cluster
            mtk::Cell const & tClusterCell = tSetClusters( 0 )->get_interpolation_cell();

            // compute the cell shape of the first cell
            tCellShape = tClusterCell.get_cell_info()->compute_cell_shape( &tClusterCell );
        }

        // within debug, checking all cells to make sure that they are the same Cell Shape
        // if cells exist
        if ( tSetClusters.size() > 0 )
        {
            // looping through the clusters
            for ( uint iCluster = 1; iCluster < tSetClusters.size(); iCluster++ )
            {
                MORIS_ASSERT( tSetClusters( iCluster )->get_interpolation_cell().get_cell_info()->compute_cell_shape( &tSetClusters( iCluster )->get_interpolation_cell() ) == tCellShape,
                        "Enriched_Interpolation_Mesh::get_IP_blockset_shape - cell shape is not consistent in the block" );
            }
        }

        return tCellShape;
    }

    // ----------------------------------------------------------------------------

    moris_id
    Enriched_Interpolation_Mesh::get_max_entity_id(
            enum EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        MORIS_ASSERT( aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT,
                "Only Elements or Nodes have max id" );

        moris::uint tNumEntities = this->get_num_entities( aEntityRank );

        moris_id tLocalMaxId = 0;

        for ( moris::uint i = 0; i < tNumEntities; i++ )
        {
            moris_id tId = this->get_glb_entity_id_from_entity_loc_index( i, aEntityRank );

            if ( tId > tLocalMaxId )
            {
                tLocalMaxId = tId;
            }
        }

        moris_id tGlobalMaxId = moris::max_all( tLocalMaxId );
        return tGlobalMaxId;
    }

    //------------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_node_owner( moris_index aNodeIndex ) const
    {
        return this->get_mtk_vertex( aNodeIndex ).get_owner();
    }

    //------------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_element_owner( moris_index aElementIndex ) const
    {
        return this->get_mtk_cell( aElementIndex ).get_owner();
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::get_adof_map(
            const uint                    aBSplineIndex,
            map< moris_id, moris_index >& aAdofMap ) const
    {
        aAdofMap.clear();

        moris_index tLocalMeshIndex = this->get_local_mesh_index( aBSplineIndex );

        for ( moris::uint iB = 0; iB < mEnrichCoeffLocToGlob( tLocalMeshIndex ).numel(); iB++ )
        {
            MORIS_ASSERT( !aAdofMap.key_exists( mEnrichCoeffLocToGlob( tLocalMeshIndex )( iB ) ),
                    "Duplicate id in the basis map detected" );

            aAdofMap[mEnrichCoeffLocToGlob( tLocalMeshIndex )( iB )] = (moris_index)iB;
        }
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat > const &
    Enriched_Interpolation_Mesh::get_enriched_coefficients_at_background_coefficient(
            moris_index const & aMeshIndex,
            moris_index         aBackgroundCoeffIndex ) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index( aMeshIndex );

        MORIS_ASSERT( aBackgroundCoeffIndex < (moris_index)mCoeffToEnrichCoeffs( tLocalMeshIndex ).size(),
                "Background coefficient index out of bounds. Be sure this is not an enriched coefficient index passed in." );

        return mCoeffToEnrichCoeffs( tLocalMeshIndex )( aBackgroundCoeffIndex );
    }

    // ----------------------------------------------------------------------------

    Cell< Matrix< IndexMat > > const &
    Enriched_Interpolation_Mesh::get_enriched_coefficients_to_background_coefficients( moris_index const & aMeshIndex ) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index( aMeshIndex );
        return mCoeffToEnrichCoeffs( tLocalMeshIndex );
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat > const &
    Enriched_Interpolation_Mesh::get_enriched_coefficient_local_to_global_map( moris_index const & aMeshIndex ) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index( aMeshIndex );

        return mEnrichCoeffLocToGlob( tLocalMeshIndex );
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Interpolation_Mesh::get_background_coefficient_local_to_global_map() const
    {
        moris::mtk::Mesh& tBackgroundMeshData = mXTKModel->get_background_mesh();

        moris::uint tNumBackgroundCoeffs = tBackgroundMeshData.get_num_entities( mBasisRank );

        Matrix< IndexMat > tCoefficientLocalToGlobal( 1, tNumBackgroundCoeffs );

        for ( moris::uint i = 0; i < tNumBackgroundCoeffs; i++ )
        {
            tCoefficientLocalToGlobal( i ) =
                    tBackgroundMeshData.get_glb_entity_id_from_entity_loc_index( (moris_index)i, mBasisRank );
        }

        return tCoefficientLocalToGlobal;
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_num_verts_per_interp_cell()
    {
        MORIS_ASSERT( mNumVertsPerInterpCell != MORIS_UINT_MAX,
                "Number of verts per interpolation cell not set" );

        return mNumVertsPerInterpCell;
    }

    // ----------------------------------------------------------------------------

    Interpolation_Vertex_Unzipped*
    Enriched_Interpolation_Mesh::get_unzipped_vertex_pointer( moris_index aVertexIndex )
    {
        MORIS_ASSERT( aVertexIndex < (moris_index)mEnrichedInterpVerts.size(),
                "Provided vertex index is out of bounds" );

        return mEnrichedInterpVerts( aVertexIndex );
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_num_background_coefficients( moris_index const & aMeshIndex ) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index( aMeshIndex );
        return mCoeffToEnrichCoeffs( tLocalMeshIndex ).size();
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::add_vertex_enrichment(
            moris_index const & aMeshIndex,
            mtk::Vertex*        aBaseInterpVertex,
            Vertex_Enrichment&  aVertexEnrichment,
            bool&               aNewVertex )
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index( aMeshIndex );

        // vertex index of the base interpolation vertex
        moris_index tBaseVertIndex = aBaseInterpVertex->get_index();

        // Number of enriched vertices related to the base vertex
        moris::uint tNumVertsEnrOnBaseVert = mBaseInterpVertToVertEnrichmentIndex( tLocalMeshIndex )( tBaseVertIndex ).size();

        // not new until we make it to the end
        aNewVertex = false;

        // if the vertex has interpolation then there is a potential that the t-matrices
        // are the same. If not we always construct a new interpolation vertex here
        if ( aBaseInterpVertex->has_interpolation( aMeshIndex ) )
        {
            // iterate through the enriched vertices related to the base vertex and see if any are equal
            for ( moris::uint i = 0; i < tNumVertsEnrOnBaseVert; i++ )
            {
                // std::cout<<" i = "<<i<<" | Base basis ind = "<<<<std::endl;
                moris_index tVertEnrIndex = mBaseInterpVertToVertEnrichmentIndex( tLocalMeshIndex )( tBaseVertIndex )( i );

                if ( aVertexEnrichment == *mInterpVertEnrichment( tLocalMeshIndex )( tVertEnrIndex ) )
                {
                    return tVertEnrIndex;
                }
            }
        }

        // if we make it through the loop without finding an enrichment vertex
        // make a new one
        aNewVertex = true;

        // index of the vertex enrichment
        moris_index tVertexEnrichmentIndex = mInterpVertEnrichment( tLocalMeshIndex ).size();

        // add to member data
        mInterpVertEnrichment( tLocalMeshIndex ).push_back( new Vertex_Enrichment( aVertexEnrichment ) );

        // add a dummy value to the parent vertex index of a vertex interpolation
        mVertexEnrichmentParentVertexIndex( tLocalMeshIndex ).push_back( tVertexEnrichmentIndex );

        mBaseInterpVertToVertEnrichmentIndex( tLocalMeshIndex )( tBaseVertIndex ).push_back( tVertexEnrichmentIndex );

        return tVertexEnrichmentIndex;
    }
    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::merge_duplicate_interpolation_vertices()
    {
        // TODO: only owning processor should merge and the non-owning procs should then request merged information

        Tracer tTracer( "XTK", "Enriched_Interpolation_Mesh", "merge_duplicate_interpolation_vertices", mXTKModel->mVerboseLevel, 0 );
        MORIS_LOG_SPEC( "Num Enriched Interpolation Vertices Pre Merge ", mEnrichedInterpVerts.size() );

        // this is how many vertices I start with
        moris_index tStartNumVerts = mEnrichedInterpVerts.size();

        // collect the vertices
        moris::Cell< moris::Cell< Interpolation_Vertex_Unzipped* > > tBaseVertexToEnrichedVertex;
        this->collect_base_vertex_to_enriched_vertex_connectivity( tBaseVertexToEnrichedVertex );

        // new vertex index
        Cell< moris_index > tNewIndex( mEnrichedInterpVerts.size(), MORIS_INDEX_MAX );
        Cell< moris_index > tNodesToDelete;
        tNodesToDelete.reserve( mEnrichedInterpVerts.size() );

        for ( moris::uint iBV = 0; iBV < tBaseVertexToEnrichedVertex.size(); iBV++ )
        {
            // keep track of the first vertex that a given vertex wants to merge with
            Cell< moris_index > tMergeWithVertex( tBaseVertexToEnrichedVertex( iBV ).size(), MORIS_INDEX_MAX );

            for ( moris::uint iEV = 0; iEV < tBaseVertexToEnrichedVertex( iBV ).size(); iEV++ )
            {
                for ( moris::uint iOtherEV = 0; iOtherEV < iEV + 1; iOtherEV++ )
                {
                    if ( iEV != iOtherEV )
                    {
                        bool tMerge = true;
                        for ( moris::uint iMeshIndex = 0; iMeshIndex < mInterpVertEnrichment.size(); iMeshIndex++ )
                        {
                            Vertex_Enrichment* tMyVertexEnrichment    = tBaseVertexToEnrichedVertex( iBV )( iEV )->get_xtk_interpolation( iMeshIndex );
                            Vertex_Enrichment* tOtherVertexEnrichment = tBaseVertexToEnrichedVertex( iBV )( iOtherEV )->get_xtk_interpolation( iMeshIndex );

                            if ( tMyVertexEnrichment != tOtherVertexEnrichment )
                            {
                                tMerge = false;
                                break;
                            }
                        }

                        if ( tMerge )
                        {
                            if ( iEV == 0 )
                            {
                                tMergeWithVertex( iEV ) = iEV;
                            }
                            else
                            {
                                tMergeWithVertex( iEV ) = iOtherEV;
                            }
                            break;
                        }
                    }
                }
            }

            // at this point, entries in tMergeWithVertex should be such that nodes merged to are MORIS_INDEX_MAX everyone else has their new vertex index.
            // iterate through these and change the new vertex index list
            for ( moris::uint iMerge = 0; iMerge < tMergeWithVertex.size(); iMerge++ )
            {
                moris_index tOldIndex = tBaseVertexToEnrichedVertex( iBV )( iMerge )->get_index();
                if ( tMergeWithVertex( iMerge ) == MORIS_INDEX_MAX )
                {
                    MORIS_ASSERT( tNewIndex( tOldIndex ) == MORIS_INDEX_MAX, "OVERWRITTING A NEW NODE INDEX" );
                    tNewIndex( tOldIndex ) = tOldIndex;
                }
                else
                {
                    MORIS_ASSERT( tNewIndex( tOldIndex ) == MORIS_INDEX_MAX, "OVERWRITTING A NEW NODE INDEX" );
                    tNewIndex( tOldIndex ) = tBaseVertexToEnrichedVertex( iBV )( tMergeWithVertex( iMerge ) )->get_index();
                    tNodesToDelete.push_back( tOldIndex );
                }
            }
        }

        // update cell to vertex connectivity (change the pointers)
        for ( moris::uint iCell = 0; iCell < mEnrichedInterpCells.size(); iCell++ )
        {
            for ( moris::uint iV = 0; iV < mEnrichedInterpCells( iCell )->mVertices.size(); iV++ )
            {
                moris_index tNewVertexIndex                    = tNewIndex( mEnrichedInterpCells( iCell )->mVertices( iV )->get_index() );
                mEnrichedInterpCells( iCell )->mVertices( iV ) = mEnrichedInterpVerts( tNewVertexIndex );
            }
        }

        // delete the removed nodes
        for ( moris::uint iDelete = 0; iDelete < tNodesToDelete.size(); iDelete++ )
        {
            delete mEnrichedInterpVerts( tNodesToDelete( iDelete ) );
            mEnrichedInterpVerts( tNodesToDelete( iDelete ) ) = nullptr;
        }

        mEnrichedInterpVerts.data().erase( std::remove( std::begin( mEnrichedInterpVerts.data() ), std::end( mEnrichedInterpVerts.data() ), nullptr ), std::end( mEnrichedInterpVerts.data() ) );

        MORIS_ERROR( mEnrichedInterpVerts.size() == tStartNumVerts - tNodesToDelete.size(), "Something went wrong merging nodes." );

        // reindex the vertices
        for ( moris::uint iV = 0; iV < mEnrichedInterpVerts.size(); iV++ )
        {
            mEnrichedInterpVerts( iV )->mVertexIndex = (moris_index)iV;
        }

        for ( moris::uint iCell = 0; iCell < mEnrichedInterpCells.size(); iCell++ )
        {
            for ( moris::uint iV = 0; iV < mEnrichedInterpCells( iCell )->mVertices.size(); iV++ )
            {
                MORIS_ERROR( mEnrichedInterpCells( iCell )->mVertices( iV ) != nullptr, "Nullptr in cell vertices" );
            }
        }

        mNumVerts = mEnrichedInterpVerts.size();

        MORIS_LOG_SPEC( "Num Enriched Interpolation Vertices Post Merge", mEnrichedInterpVerts.size() );
    }

    // ----------------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::collect_base_vertex_to_enriched_vertex_connectivity( moris::Cell< moris::Cell< Interpolation_Vertex_Unzipped* > >& aBaseVertexToEnrichedVertex )
    {
        // allocate space in the cell of cell
        aBaseVertexToEnrichedVertex.resize( mXTKModel->get_background_mesh().get_num_nodes() );

        for ( moris::uint iN = 0; iN < mEnrichedInterpVerts.size(); iN++ )
        {
            moris_index tBaseVertIndex = mEnrichedInterpVerts( iN )->get_base_vertex()->get_index();
            aBaseVertexToEnrichedVertex( tBaseVertIndex ).push_back( mEnrichedInterpVerts( iN ) );
        }
    }

    // ----------------------------------------------------------------------------

    Vertex_Enrichment*
    Enriched_Interpolation_Mesh::get_vertex_enrichment(
            moris_index const & aMeshIndex,
            moris_index const & aVertexEnrichmentIndex )
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index( aMeshIndex );

        MORIS_ASSERT( aVertexEnrichmentIndex < (moris_index)mInterpVertEnrichment( tLocalMeshIndex ).size(),
                "Provided vertex enrichment index out of bounds" );

        return mInterpVertEnrichment( tLocalMeshIndex )( aVertexEnrichmentIndex );
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_vertex_related_to_vertex_enrichment(
            moris_index const & aMeshIndex,
            moris_index         aVertexEnrichmentIndex ) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index( aMeshIndex );

        MORIS_ASSERT( aVertexEnrichmentIndex < (moris_index)mVertexEnrichmentParentVertexIndex( tLocalMeshIndex ).size(),
                "Provided vertex enrichment index out of bounds" );

        return mVertexEnrichmentParentVertexIndex( tLocalMeshIndex )( aVertexEnrichmentIndex );
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_local_mesh_index( moris_index const & aMeshIndex ) const
    {
        auto tIter = mMeshIndexToLocMeshIndex.find( aMeshIndex );

        MORIS_ASSERT( tIter != mMeshIndexToLocMeshIndex.end(), "Mesh index not in map" );

        return tIter->second;
    }

    // ----------------------------------------------------------------------------

    bool
    Enriched_Interpolation_Mesh::basis_exists_on_partition(
            moris_index const & aMeshIndex,
            moris_index const & aBasisId )
    {
        moris_index tMeshIndex = this->get_local_mesh_index( aMeshIndex );

        if ( mGlobaltoLocalBasisMaps( tMeshIndex ).find( aBasisId ) == mGlobaltoLocalBasisMaps( tMeshIndex ).end() )
        {
            return false;
        }

        return true;
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::add_basis_function(
            moris_index const & aMeshIndex,
            moris_index const & aBasisIdToAdd,
            moris_index const & aBasisOwner,
            moris_index const & aBasisBulkPhase )
    {
        MORIS_ASSERT( !this->basis_exists_on_partition( aMeshIndex, aBasisIdToAdd ),
                "Basis that you are trying to add already exists in this mesh" );

        moris_index tLocMesh  = this->get_local_mesh_index( aMeshIndex );
        moris_index tNewIndex = mEnrichCoeffLocToGlob( tLocMesh ).numel();

        // add a size of 1
        mEnrichCoeffLocToGlob( tLocMesh ).resize( 1, tNewIndex + 1 );
        mEnrichCoeffOwnership( tLocMesh ).resize( 1, tNewIndex + 1 );
        mEnrichCoeffBulkPhase( tLocMesh ).resize( 1, tNewIndex + 1 );

        // add the local to glb map
        mEnrichCoeffLocToGlob( tLocMesh )( tNewIndex ) = aBasisIdToAdd;
        mEnrichCoeffOwnership( tLocMesh )( tNewIndex ) = aBasisOwner;
        mEnrichCoeffBulkPhase( tLocMesh )( tNewIndex ) = aBasisBulkPhase;

        // add to glb to local map
        mGlobaltoLocalBasisMaps( tLocMesh )[aBasisIdToAdd] = tNewIndex;

        return tNewIndex;
    }

    // ----------------------------------------------------------------------------

    moris::Cell< Interpolation_Cell_Unzipped const * >
    Enriched_Interpolation_Mesh::get_enriched_cells_from_base_cells(
            moris::Cell< moris::mtk::Cell const * > const & aBaseCells ) const
    {
        uint tNumBaseCells = aBaseCells.size();

        moris::Cell< Interpolation_Cell_Unzipped const * > tEnrichedCells;

        for ( moris::uint i = 0; i < tNumBaseCells; i++ )
        {
            tEnrichedCells.append( this->get_enriched_cells_from_base_cell( aBaseCells( i ) ) );
        }

        return tEnrichedCells;
    }

    // ----------------------------------------------------------------------------

    Cell< Interpolation_Cell_Unzipped* > const &
    Enriched_Interpolation_Mesh::get_enriched_interpolation_cells() const
    {
        return mEnrichedInterpCells;
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_num_interpolation_types() const
    {
        return mMeshIndices.numel();
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_interpolation_index( moris_index const & aLocalInterpIndex ) const
    {
        MORIS_ASSERT( aLocalInterpIndex < (moris_index)mMeshIndices.numel(),
                "Local interpolation index out of bounds" );

        return mMeshIndices( aLocalInterpIndex );
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_basis_owner(
            moris_index aBasisIndex,
            moris_index aMeshIndex )
    {
        moris_index tLocMesh = this->get_local_mesh_index( aMeshIndex );

        return mEnrichCoeffOwnership( tLocMesh )( aBasisIndex );
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_basis_bulk_phase( moris_index const & aBasisIndex,
            moris_index const &                                            aMeshIndex ) const
    {
        moris_index tLocMesh = this->get_local_mesh_index( aMeshIndex );

        return mEnrichCoeffBulkPhase( tLocMesh )( aBasisIndex );
    }

    // ----------------------------------------------------------------------------

    Cell< Interpolation_Cell_Unzipped* >&
    Enriched_Interpolation_Mesh::get_enriched_interpolation_cells()
    {
        return mEnrichedInterpCells;
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::get_owned_and_not_owned_enriched_interpolation_cells(
            Cell< Interpolation_Cell_Unzipped* >&         aOwnedInterpCells,
            Cell< Cell< Interpolation_Cell_Unzipped* > >& aNotOwnedInterpCells,
            Cell< uint >&                                 aProcRanks )
    {
        // get all interp cells
        Cell< Interpolation_Cell_Unzipped* >& tEnrInterpCells = this->get_enriched_interpolation_cells();

        // reserve space
        aOwnedInterpCells.resize( 0 );
        aOwnedInterpCells.reserve( tEnrInterpCells.size() );
        aNotOwnedInterpCells.resize( 0 );

        // proc rank
        moris_index tParRank = par_rank();

        // counter
        moris::uint         tOwnerCount = 0;
        Cell< moris::uint > tCounts( 0 );

        // map
        std::unordered_map< moris_id, moris_id > tProcRankToDataIndex;

        // access the communication table
        Matrix< IdMat > tCommTable = this->get_communication_table();

        // resize proc ranks and setup map to comm table
        aProcRanks.resize( tCommTable.numel() );
        for ( moris::uint i = 0; i < tCommTable.numel(); i++ )
        {
            tProcRankToDataIndex[tCommTable( i )] = i;
            aProcRanks( i )                       = ( tCommTable( i ) );
            aNotOwnedInterpCells.push_back( Cell< Interpolation_Cell_Unzipped* >( 0 ) );
        }

        for ( moris::uint i = 0; i < tEnrInterpCells.size(); i++ )
        {
            moris_index tOwnerProc = tEnrInterpCells( i )->get_owner();

            if ( tParRank == tOwnerProc )
            {
                aOwnedInterpCells.push_back( tEnrInterpCells( i ) );
                tOwnerCount++;
            }
            else
            {
                moris_index tProcDataIndex = tProcRankToDataIndex[tOwnerProc];
                aNotOwnedInterpCells( tProcDataIndex ).push_back( tEnrInterpCells( i ) );
            }
        }
    }

    // ----------------------------------------------------------------------------

    Interpolation_Vertex_Unzipped&
    Enriched_Interpolation_Mesh::get_xtk_interp_vertex( moris::uint aVertexIndex )
    {
        return *mEnrichedInterpVerts( aVertexIndex );
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::add_proc_to_comm_table( moris_index aProcRank )
    {
        mXTKModel->get_cut_integration_mesh()->add_proc_to_comm_table( aProcRank );
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_enr_basis_index_from_enr_basis_id(
            moris_index const & aMeshIndex,
            moris_index const & aBasisId ) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index( aMeshIndex );

        auto tIter = mGlobaltoLocalBasisMaps( tLocalMeshIndex ).find( aBasisId );

        MORIS_ASSERT( tIter != mGlobaltoLocalBasisMaps( tLocalMeshIndex ).end(),
                "Basis id not in map" );

        return tIter->second;
    }

    // ----------------------------------------------------------------------------

    moris::Cell< Interpolation_Cell_Unzipped const * >
    Enriched_Interpolation_Mesh::get_enriched_cells_from_base_cell( moris::mtk::Cell const * aBaseCells ) const
    {
        moris_index tBaseIndex = aBaseCells->get_index();

        MORIS_ASSERT( tBaseIndex < (moris_index)mBaseCelltoEnrichedCell.size(),
                "Base Cell index is out of bounds. This index is related to the non-enriched interpolation mesh. Make sure enriched cell is not passed into this function" );

        uint tNumEnrichedCells = mBaseCelltoEnrichedCell( tBaseIndex ).size();

        moris::Cell< Interpolation_Cell_Unzipped const * > tEnrichedCellPtrs( tNumEnrichedCells );

        for ( moris::uint i = 0; i < tNumEnrichedCells; i++ )
        {
            tEnrichedCellPtrs( i ) = ( mBaseCelltoEnrichedCell( tBaseIndex )( i ) );
        }

        return tEnrichedCellPtrs;
    }

    // ----------------------------------------------------------------------------

    Interpolation_Vertex_Unzipped const &
    Enriched_Interpolation_Mesh::get_xtk_interp_vertex( moris::uint aVertexIndex ) const
    {
        return *mEnrichedInterpVerts( aVertexIndex );
    }

    // ----------------------------------------------------------------------------

    Matrix< IdMat >
    Enriched_Interpolation_Mesh::convert_indices_to_ids(
            Matrix< IndexMat > const & aIndices,
            enum EntityRank            aEntityRank ) const
    {
        moris::uint tNRow = aIndices.n_rows();
        moris::uint tNCol = aIndices.n_cols();

        Matrix< IdMat > tIds( tNRow, tNCol );

        for ( moris::uint i = 0; i < tNRow; i++ )
        {
            for ( moris::uint j = 0; j < tNCol; j++ )
            {
                tIds( i, j ) = this->get_glb_entity_id_from_entity_loc_index( aIndices( i, j ), aEntityRank );
            }
        }

        return tIds;
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Interpolation_Mesh::convert_ids_to_indices(
            Matrix< IdMat > const & aIds,
            enum EntityRank         aEntityRank ) const
    {
        moris::uint tNRow = aIds.n_rows();
        moris::uint tNCol = aIds.n_cols();

        Matrix< IdMat > tIndices( tNRow, tNCol );

        for ( moris::uint i = 0; i < tNRow; i++ )
        {
            for ( moris::uint j = 0; j < tNCol; j++ )
            {
                tIndices( i, j ) = this->get_loc_entity_ind_from_entity_glb_id( aIds( i, j ), aEntityRank );
            }
        }

        return tIndices;
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::convert_enriched_basis_indices_to_ids(
            moris_index const &        aMeshIndex,
            Matrix< IndexMat > const & aEnrichedIndices,
            Matrix< IdMat >&           aEnrichedIds ) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index( aMeshIndex );

        aEnrichedIds.resize( aEnrichedIndices.n_rows(), aEnrichedIndices.n_cols() );

        for ( moris::uint i = 0; i < aEnrichedIndices.n_rows(); i++ )
        {
            for ( moris::uint j = 0; j < aEnrichedIndices.n_cols(); j++ )
            {
                aEnrichedIds( i, j ) = mEnrichCoeffLocToGlob( tLocalMeshIndex )( aEnrichedIndices( i, j ) );
            }
        }
    }

    // ----------------------------------------------------------------------------

    moris::Memory_Map
    Enriched_Interpolation_Mesh::get_memory_usage()
    {
        // memory map of ig mesh
        moris::Memory_Map tMM;
        tMM.mMemoryMapData["mXTKModel ptr"]  = sizeof( mXTKModel );
        tMM.mMemoryMapData["mBasisRank ptr"] = sizeof( mBasisRank );
        tMM.mMemoryMapData["mMeshIndices"]   = mMeshIndices.capacity();
        // FIXME: add mMeshIndexToLocMeshIndex
        tMM.mMemoryMapData["mNumVerts"]                            = sizeof( mNumVerts );
        tMM.mMemoryMapData["mNumVertsPerInterpCell"]               = sizeof( mNumVertsPerInterpCell );
        tMM.mMemoryMapData["mBaseInterpVertToVertEnrichmentIndex"] = moris::internal_capacity_nested( mBaseInterpVertToVertEnrichmentIndex );
        tMM.mMemoryMapData["mInterpVertEnrichment"]                = moris::internal_capacity_nested_ptr( mInterpVertEnrichment );
        tMM.mMemoryMapData["mVertexEnrichmentParentVertexIndex"]   = moris::internal_capacity( mVertexEnrichmentParentVertexIndex );
        tMM.mMemoryMapData["mVertexBulkPhase"]                     = mVertexBulkPhase.capacity();
        tMM.mMemoryMapData["mCoeffToEnrichCoeffs"]                 = moris::internal_capacity_nested( mCoeffToEnrichCoeffs );
        tMM.mMemoryMapData["mEnrichCoeffLocToGlob"]                = moris::internal_capacity( mEnrichCoeffLocToGlob );
        // fixme: add me mGlobaltoLocalBasisMaps
        tMM.mMemoryMapData["mEnrichCoeffOwnership"] = moris::internal_capacity( mEnrichCoeffOwnership );
        tMM.mMemoryMapData["mLocalToGlobalMaps"]    = moris::internal_capacity( mLocalToGlobalMaps );
        // fixme: add mGlobaltoLobalMaps
        tMM.mMemoryMapData["mBaseCelltoEnrichedCell"] = moris::internal_capacity( mBaseCelltoEnrichedCell );
        tMM.mMemoryMapData["mCellInfo"]               = sizeof( mCellInfo );
        tMM.mMemoryMapData["mNotOwnedVerts"]          = mNotOwnedVerts.capacity();
        tMM.mMemoryMapData["mNotOwnedBasis"]          = mNotOwnedBasis.capacity();
        tMM.mMemoryMapData["mOwnedBasis"]             = mOwnedBasis.capacity();
        return tMM;
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::convert_enriched_basis_indices_to_ids(
            moris_index const &                aMeshIndex,
            Cell< Matrix< IndexMat > > const & aEnrichedIndices,
            Cell< Matrix< IdMat > >&           aEnrichedIds ) const
    {
        aEnrichedIds.resize( aEnrichedIndices.size() );

        for ( moris::uint i = 0; i < aEnrichedIndices.size(); i++ )
        {
            this->convert_enriched_basis_indices_to_ids(
                    aMeshIndex,
                    aEnrichedIndices( i ),
                    aEnrichedIds( i ) );
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::write_diagnostics()
    {
        std::string tCellDiagFile   = mXTKModel->get_diagnostic_file_name( std::string( "Enr_IP_Cells" ) );
        std::string tVertexDiagFile = mXTKModel->get_diagnostic_file_name( std::string( "Enr_IP_Verts" ) );
        // this->print_vertex_maps();
        // this->print_enriched_cell_maps();
        this->print_enriched_cells( tCellDiagFile );
        this->print_enriched_verts( tVertexDiagFile );

        for ( moris::size_t iBasisType = 0; iBasisType < mMeshIndices.numel(); iBasisType++ )
        {
            // get the mesh index
            moris_index tMeshIndex = mMeshIndices( iBasisType );

            std::string tEnrVertexInterpolationFile = mXTKModel->get_diagnostic_file_name( std::string( "Ip_Vertex_Interp" + std::to_string( tMeshIndex ) ) );

            this->print_enriched_verts_interpolation( tMeshIndex, tEnrVertexInterpolationFile );
        }
        // this->print_basis_to_enriched_basis();
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::print_enriched_cells( std::string aFile )
    {
        Cell< Interpolation_Cell_Unzipped* > const & tEnrIPCells = this->get_enriched_interpolation_cells();

        std::ostringstream tStringStream;
        // max num verts to cells
        uint tMaxVertsToCell = 0;
        for ( moris::uint i = 0; i < this->get_num_entities( EntityRank::ELEMENT, 0 ); i++ )
        {
            Interpolation_Cell_Unzipped const * tCell = tEnrIPCells( i );
            if ( tCell->get_number_of_vertices() > tMaxVertsToCell )
            {
                tMaxVertsToCell = tCell->get_number_of_vertices();
            }
        }

        tStringStream << "Cell_Id,";
        tStringStream << "Cell_Ind,";
        tStringStream << "Owner,";
        tStringStream << "PRank,";
        tStringStream << "Base_Cell_Id,";
        tStringStream << "Bulk_Phase,";
        tStringStream << "Sub_Phase_Id,";
        tStringStream << "Measure,";
        for ( moris::uint iVH = 0; iVH < tMaxVertsToCell; iVH++ )
        {
            tStringStream << "Vert_" + std::to_string( iVH );

            if ( iVH != tMaxVertsToCell - 1 )
            {
                tStringStream << ",";
            }
        }
        tStringStream << "\n";

        for ( moris::uint i = 0; i < this->get_num_entities( EntityRank::ELEMENT, 0 ); i++ )
        {
            Interpolation_Cell_Unzipped const * tCell     = tEnrIPCells( (moris_index)i );
            moris::Cell< moris::mtk::Vertex* >  tVertices = tCell->get_vertex_pointers();

            tStringStream << tCell->get_id() << ",";
            tStringStream << tCell->get_index() << ",";
            tStringStream << std::to_string( tCell->get_owner() ) << ",";
            tStringStream << std::to_string( par_rank() ) << ",";
            tStringStream << tCell->get_base_cell()->get_id() << ",";
            tStringStream << tCell->get_bulkphase_index() << ",";
            tStringStream << mXTKModel->get_cut_integration_mesh()->get_subphase_id( (moris_index)tCell->get_subphase_index() ) << ",";
            tStringStream << std::scientific << tCell->compute_cell_measure() << ",";

            for ( moris::uint j = 0; j < tMaxVertsToCell; j++ )
            {
                if ( j < tVertices.size() )
                {
                    tStringStream << std::to_string( tVertices( j )->get_id() );
                }
                else
                {
                    tStringStream << std::to_string( MORIS_INDEX_MAX );
                }

                if ( j != tMaxVertsToCell - 1 )
                {
                    tStringStream << ",";
                }
            }
            tStringStream << "\n";
        }

        if ( aFile.empty() == false )
        {
            std::ofstream tOutputFile( aFile );
            tOutputFile << tStringStream.str() << std::endl;
            tOutputFile.close();
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::print_enriched_verts( std::string aFile )
    {
        std::ostringstream tStringStream;
        tStringStream.clear();
        tStringStream.str( "" );

        tStringStream << "Vert_Id,";
        tStringStream << "Base_Vert_Id,";
        tStringStream << "Vert Ind,";
        tStringStream << "Owner,";
        tStringStream << "Prank,";
        tStringStream << "Bulk_Phase,";
        tStringStream << "Subphase,";

        for ( moris::uint iVH = 0; iVH < this->get_spatial_dim(); iVH++ )
        {
            tStringStream << "Coords_" + std::to_string( iVH );

            if ( iVH != this->get_spatial_dim() - 1 )
            {
                tStringStream << ",";
            }
        }

        tStringStream << std::endl;

        for ( moris::uint i = 0; i < this->get_num_entities( EntityRank::NODE, 0 ); i++ )
        {
            Interpolation_Vertex_Unzipped const & tVertex = this->get_xtk_interp_vertex( (moris_index)i );
            tStringStream.precision( 16 );

            tStringStream << tVertex.get_id() << ",";
            tStringStream << tVertex.get_base_vertex()->get_id() << ",";
            tStringStream << tVertex.get_index() << ",";
            tStringStream << tVertex.get_owner() << ",";
            tStringStream << par_rank() << ",";
            tStringStream << mVertexBulkPhase( i ) << ",";
            tStringStream << mVertexMaxSubphase( i ) << ",";

            moris::Matrix< moris::DDRMat > tCoords = tVertex.get_coords();

            for ( moris::uint iSp = 0; iSp < this->get_spatial_dim(); iSp++ )
            {
                tStringStream << std::scientific << tCoords( iSp );

                if ( iSp != this->get_spatial_dim() - 1 )
                {
                    tStringStream << ",";
                }
            }

            //
            tStringStream << std::endl;
        }
        if ( aFile.empty() == false )
        {
            std::ofstream tOutputFile( aFile );
            tOutputFile << tStringStream.str() << std::endl;
            tOutputFile.close();
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::print_enriched_verts_interpolation(
            const moris_index& aMeshIndex,
            std::string        aFileName )
    {
        std::ostringstream tStringStream;
        tStringStream.precision( 16 );
        tStringStream << "Vert_Id,";
        tStringStream << "Num_Coeffs,";

        // global max size of
        moris_index tLocalTMatrixSize = 0;
        for ( moris::uint iV = 0; iV < this->get_num_nodes(); iV++ )
        {
            mtk::Vertex_Interpolation* tVertexInterp = this->get_mtk_vertex( (moris_index)iV ).get_interpolation( aMeshIndex );
            Matrix< IdMat >            tBasisIds     = tVertexInterp->get_ids();

            if ( (moris_index)tBasisIds.numel() > tLocalTMatrixSize )
            {
                tLocalTMatrixSize = (moris_index)tBasisIds.numel();
            }
        }

        moris_index tGlbMaxTMatrixSize = moris::max_all( tLocalTMatrixSize );

        for ( moris_index iCH = 0; iCH < tGlbMaxTMatrixSize; iCH++ )
        {
            tStringStream << "Basis_ID" + std::to_string( iCH ) << ",";
            tStringStream << "Basis_Weight" + std::to_string( iCH );

            if ( iCH != tGlbMaxTMatrixSize - 1 )
            {
                tStringStream << ",";
            }
        }

        tStringStream << "\n";

        for ( moris::uint iV = 0; iV < this->get_num_nodes(); iV++ )
        {
            tStringStream << this->get_mtk_vertex( (moris_index)iV ).get_id() << ",";
            mtk::Vertex_Interpolation* tVertexInterp = this->get_mtk_vertex( (moris_index)iV ).get_interpolation( aMeshIndex );
            Matrix< IdMat >            tBasisIds     = tVertexInterp->get_ids();
            const Matrix< DDRMat >*    tBasisWeights = tVertexInterp->get_weights();

            tStringStream << tBasisIds.numel() << ",";
            MORIS_ASSERT( tBasisIds.numel() == tBasisWeights->numel(), "Size mismatch" );

            for ( moris::uint iB = 0; iB < tBasisIds.numel(); iB++ )
            {
                tStringStream << std::to_string( tBasisIds( iB ) ) << ",";

                tStringStream << ( *tBasisWeights )( iB );
                if ( iB != tBasisIds.numel() - 1 )
                {
                    tStringStream << ",";
                }
            }

            // maybe a t-matrix isn't needed on this proc in this case we have no data in the vertex enrichment data
            if ( tBasisIds.numel() == 0 )
            {
                tStringStream << std::to_string( MORIS_INDEX_MAX ) << ",";
                tStringStream << std::numeric_limits< moris::real >::quiet_NaN();
            }

            tStringStream << "\n";
        }

        if ( aFileName.empty() == false )
        {
            std::ofstream tOutputFile( aFileName );
            tOutputFile << tStringStream.str() << std::endl;
            tOutputFile.close();
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::print_vertex_maps() const
    {
        moris::uint tNumNodes = this->get_num_entities( EntityRank::NODE );

        moris::uint tMapIndex = (uint)( EntityRank::NODE );

        MORIS_ASSERT( mLocalToGlobalMaps( tMapIndex ).numel() == tNumNodes,
                "Enriched_Interpolation_Mesh::print_vertex_maps: number of nodes and size of map do not match." );

        std::cout << "\nVertex Map:" << std::endl;

        for ( moris::uint i = 0; i < tNumNodes; i++ )
        {
            std::cout << "    Vertex Index: " << std::setw( 9 ) << i << " | Vertex Id: " << std::setw( 9 ) << mLocalToGlobalMaps( tMapIndex )( i ) << std::endl;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::print_enriched_cell_maps() const
    {
        moris::uint tNumCells = this->get_num_entities( EntityRank::ELEMENT );

        moris::uint tMapIndex = (uint)( EntityRank::ELEMENT );

        MORIS_ASSERT( mLocalToGlobalMaps( tMapIndex ).numel() == tNumCells,
                "Enriched_Interpolation_Mesh::print_enriched_cell_maps: number of elements and size of map do not match." );

        std::cout << "\nCell Map:" << std::endl;

        for ( moris::uint i = 0; i < tNumCells; i++ )
        {
            std::cout << "    Cell Index: " << std::setw( 9 ) << i << " | Cell Id: " << std::setw( 9 ) << mLocalToGlobalMaps( tMapIndex )( i ) << std::endl;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::print_basis_to_enriched_basis() const
    {
        for ( moris::uint iM = 0; iM < mMeshIndices.numel(); iM++ )
        {
            moris::uint tNumBasis = mCoeffToEnrichCoeffs( iM ).size();

            std::cout << "\nBackground Basis to Enriched Basis Indices For Mesh: " << mMeshIndices( iM ) << std::endl;

            for ( moris::uint iB = 0; iB < tNumBasis; iB++ )
            {
                std::cout << "    Basis Index: " << std::setw( 9 ) << iB << " | Enriched Indices";

                for ( moris::uint iEB = 0; iEB < mCoeffToEnrichCoeffs( iM )( iB ).numel(); iEB++ )
                {
                    std::cout << std::setw( 9 ) << mCoeffToEnrichCoeffs( iM )( iB )( iEB );
                }
                std::cout << std::endl;
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::print_vertex_interpolation() const
    {
        moris::uint tNumVerts = this->get_num_entities( EntityRank::NODE );

        std::cout << "\nVertex Interpolation:" << std::endl;

        for ( moris::moris_index i = 0; i < (moris_index)tNumVerts; i++ )
        {
            Interpolation_Vertex_Unzipped const & tVertex = this->get_xtk_interp_vertex( i );

            std::cout << "\nVertex Id: " << std::setw( 9 ) << tVertex.get_id() << std::endl;
            std::cout << *tVertex.get_xtk_interpolation( 0 ) << std::endl;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::print_basis_information() const
    {
        std::cout << "\nBasis Information on proc " << par_rank() << ":" << std::endl;

        for ( moris::moris_index iM = 0; iM < (moris_index)mMeshIndices.numel(); iM++ )
        {
            std::cout << " Mesh Index: " << mMeshIndices( iM ) << std::endl;

            for ( moris::moris_index i = 0; i < (moris_index)mEnrichCoeffLocToGlob( iM ).numel(); i++ )
            {
                moris_id tId = mEnrichCoeffLocToGlob( iM )( i );

                moris_index tIndex = this->get_enr_basis_index_from_enr_basis_id( mMeshIndices( iM ), tId );

                std::cout << "    Basis Id: " << std::setw( 9 ) << tId << " | Basis Index: " << std::setw( 9 ) << tIndex << std::endl;
            }
        }
    }

    // ----------------------------------------------------------------------------

    bool
    Enriched_Interpolation_Mesh::verify_basis_interpolating_into_cluster(
            mtk::Cluster const &    aCluster,
            moris_index const &     aMeshIndex,
            const mtk::Master_Slave aIsMaster )
    {
        bool tDiagnosticFlag = true;

        // interpolation cell
        moris::mtk::Cell const & tIpCell = aCluster.get_interpolation_cell( aIsMaster );

        // get the xtk interpolation cell
        Interpolation_Cell_Unzipped* tEnrichedIpCell = mEnrichedInterpCells( tIpCell.get_index() );

        MORIS_ASSERT( tIpCell.get_id() == tEnrichedIpCell->get_id(), "Id mismatch" );

        // get the bulkphase index
        moris_index tBulkPhase = tEnrichedIpCell->get_bulkphase_index();

        // vertices attached to the interpolation cells
        moris::Cell< mtk::Vertex* > tVertexPointers = tIpCell.get_vertex_pointers();

        // iterate through vertex pointers
        for ( moris::uint i = 0; i < tVertexPointers.size(); i++ )
        {
            // todo:remove
            if ( !tVertexPointers( i )->has_interpolation( aMeshIndex ) )
            {
                std::cout << " Id = " << tVertexPointers( i )->get_id()
                          << " | back vertex id = " << mEnrichedInterpVerts( tVertexPointers( i )->get_index() )->get_id()
                          << " | owner = " << tVertexPointers( i )->get_owner()
                          << " | my rank = " << par_rank()
                          << " | IP Cell Owner = " << tEnrichedIpCell->get_owner()
                          << " | IP Cell ID = " << tEnrichedIpCell->get_id() << std::endl;

                moris::print( tVertexPointers( i )->get_coords(), "Coordinate of vertex" );
            }

            //            MORIS_ASSERT(tVertexPointers(i)->has_interpolation(aMeshIndex),"Vertex does not have interpolation.");

            mtk::Vertex_Interpolation* tVertexInterp = tVertexPointers( i )->get_interpolation( aMeshIndex );

            // get the basis indices
            Matrix< IndexMat > tBasisIndices = tVertexInterp->get_indices();
            Matrix< IndexMat > tBasisIds     = tVertexInterp->get_ids();

            if ( tBasisIds.numel() == 0 )
            {
                std::cout << "Vertex id = " << tVertexPointers( i )->get_id()
                          << " | back vertex id = " << mEnrichedInterpVerts( tVertexPointers( i )->get_index() )->get_id() << std::endl;
                MORIS_ASSERT( 0, "Interpolation Vertex in cluster does not have at least 1 basis interpolating into it" );
                tDiagnosticFlag = false;
            }

            //  iterate through basis indices and check bulk phases
            for ( moris::uint iB = 0; iB < tBasisIndices.numel(); iB++ )
            {
                if ( this->get_basis_bulk_phase( tBasisIndices( iB ), aMeshIndex ) != tBulkPhase )
                {
                    std::cout << "Vertex id = " << tVertexPointers( i )->get_id()
                              << " | back vertex id = " << mEnrichedInterpVerts( tVertexPointers( i )->get_index() )->get_id()
                              << " | tBasisIds(iB) = " << tBasisIds( iB )
                              << " | MeshBP" << this->get_basis_bulk_phase( tBasisIndices( iB ), aMeshIndex )
                              << " | tBulkPhase = " << tBulkPhase << std::endl;
                    MORIS_ASSERT( 0, "Basis in cluster does not interpolate into same bulk phase as interpolation cell" );
                    tDiagnosticFlag = false;
                }
            }
        }

        return tDiagnosticFlag;
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::finalize_setup()
    {
        this->setup_local_to_global_maps();

        this->setup_vertex_to_bulk_phase();

        this->setup_basis_ownership();

        this->setup_basis_to_bulk_phase();

        MORIS_ASSERT( this->verify_basis_support(), "Issue detected in basis support." );
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_basis_to_bulk_phase()
    {
        // size member data
        mEnrichCoeffBulkPhase.resize( mMeshIndices.max() + 1 );

        // verify all the subphases in a enriched basis support are the same bulk phase
        for ( moris::size_t iMesh = 0; iMesh < mMeshIndices.numel(); iMesh++ )
        {
            // Mesh index
            moris_index tMeshIndex = mMeshIndices( iMesh );

            // Number of enriched functions
            moris::size_t tNumEnrBasis = mEnrichCoeffLocToGlob( tMeshIndex ).numel();

            // allocate interpolation cells in basis support
            moris::Cell< moris::Cell< Interpolation_Cell_Unzipped* > > tCellsInEnrSupports( tNumEnrBasis );

            Cell< Interpolation_Cell_Unzipped* > const & tEnrIpCells = this->get_enriched_interpolation_cells();

            // number of cells
            moris_index tNumCells = this->get_num_entities( EntityRank::ELEMENT );

            MORIS_ASSERT( tNumCells == (moris_index)tEnrIpCells.size(), "Inconsistent num cells information." );

            // create the enriched interpolation basis to interpolation cell interpolation
            for ( moris_index i = 0; i < tNumCells; i++ )
            {
                moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const & tVertices = tEnrIpCells( i )->get_xtk_interpolation_vertices();

                for ( moris::size_t iV = 0; iV < tVertices.size(); iV++ )
                {
                    if ( tVertices( iV )->has_interpolation( tMeshIndex ) )
                    {
                        Vertex_Enrichment* tVertInterp = tVertices( iV )->get_xtk_interpolation( tMeshIndex );

                        Matrix< IndexMat > tBasisIndices = tVertInterp->get_indices();

                        // iterate through vertices and add the ip cell to the support
                        for ( moris::uint iB = 0; iB < tBasisIndices.numel(); iB++ )
                        {
                            tCellsInEnrSupports( tBasisIndices( iB ) ).push_back( tEnrIpCells( i ) );
                        }
                    }
                }
            }

            mEnrichCoeffBulkPhase( tMeshIndex ).resize( 1, tNumEnrBasis );
            mEnrichCoeffBulkPhase( tMeshIndex ).fill( MORIS_INDEX_MAX );

            // iterate through enriched basis functions
            for ( moris::size_t i = 0; i < tNumEnrBasis; i++ )
            {
                // iterate through enriched interpolation cells in the support
                moris_index tNumSubphaseInSupport = tCellsInEnrSupports( i ).size();

                for ( moris::moris_index iSP = 0; iSP < tNumSubphaseInSupport; iSP++ )
                {
                    Interpolation_Cell_Unzipped* tIpCell = tCellsInEnrSupports( i )( iSP );

                    moris_index tBulkPhase = tIpCell->get_bulkphase_index();

                    MORIS_ASSERT( tBulkPhase != MORIS_INDEX_MAX, "Bulk phase index not set." );

                    if ( iSP == 0 )
                    {
                        mEnrichCoeffBulkPhase( tMeshIndex )( i ) = tBulkPhase;
                    }

                    if ( tBulkPhase != mEnrichCoeffBulkPhase( tMeshIndex )( i ) )
                    {
                        std::cout << "tExpectedBulkPhase=  " << mEnrichCoeffBulkPhase( tMeshIndex )( iSP ) << " | tBulkPhase = " << tBulkPhase << std::endl;
                        MORIS_ERROR( 0, "Subphase in enriched basis function support not consistent bulkphase" );
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_vertex_to_bulk_phase()
    {
        moris::size_t tNumVertices = this->get_num_nodes();

        // size member data
        mVertexBulkPhase.resize( 1, tNumVertices );
        mVertexBulkPhase.fill( MORIS_INDEX_MAX );
        mVertexMaxSubphase.resize( 1, tNumVertices );
        mVertexMaxSubphase.fill( -1 );

        Cell< Interpolation_Cell_Unzipped* > const & tEnrIpCells = this->get_enriched_interpolation_cells();

        // number of cells
        moris_index tNumCells = this->get_num_entities( EntityRank::ELEMENT );

        // create the enriched interpolation basis to interpolation cell interpolation
        for ( moris::moris_index i = 0; i < tNumCells; i++ )
        {
            moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const & tVertices = tEnrIpCells( i )->get_xtk_interpolation_vertices();

            moris_index const tSubphaseIndex = tEnrIpCells( i )->get_subphase_index();
            moris_index const tSubphaseId    = mXTKModel->get_cut_integration_mesh()->get_subphase_id( tSubphaseIndex );

            for ( moris::size_t iV = 0; iV < tVertices.size(); iV++ )
            {
                moris_index tVertexIndex = tVertices( iV )->get_index();

                if ( mVertexBulkPhase( tVertexIndex ) == MORIS_INDEX_MAX )
                {
                    mVertexBulkPhase( tVertexIndex ) = tEnrIpCells( i )->get_bulkphase_index();
                }

                if ( tSubphaseId > mVertexMaxSubphase( tVertexIndex ) )
                {
                    mVertexMaxSubphase( tVertexIndex ) = tSubphaseId;
                }

                else
                {
                    if ( mVertexBulkPhase( tVertexIndex ) != tEnrIpCells( i )->get_bulkphase_index() )
                    {

                        std::cout << " Vert Id = " << tVertices( iV )->get_id()
                                  << " | Vert Owner = " << tVertices( iV )->get_owner()
                                  << " | mVertexBulkPhase(tVertexIndex) = " << mVertexBulkPhase( tVertexIndex )
                                  << " | tEnrIpCells(i)->get_bulkphase_index() =" << tEnrIpCells( i )->get_bulkphase_index()
                                  << std::endl;
                    }
                    MORIS_ASSERT( mVertexBulkPhase( tVertexIndex ) == tEnrIpCells( i )->get_bulkphase_index(), "Inconsistent vertex bulk phase" );
                }
            }
        }
    }

    // ----------------------------------------------------------------------------

    bool
    Enriched_Interpolation_Mesh::verify_basis_support()
    {
        bool tSubphaseBulkPhasesInSupportDiag = true;

        // verify all the subphases in a enriched basis support are the same bulk phase
        for ( moris::size_t iMesh = 0; iMesh < mMeshIndices.numel(); iMesh++ )
        {
            moris_index tMeshIndex = mMeshIndices( iMesh );

            // Number of enriched functions
            moris::size_t tNumEnrBasis = mEnrichCoeffLocToGlob( tMeshIndex ).numel();

            // allocate interpolation cells in basis support
            moris::Cell< moris::Cell< Interpolation_Cell_Unzipped* > > tCellsInEnrSupports( tNumEnrBasis );

            Cell< Interpolation_Cell_Unzipped* > const & tEnrIpCells = this->get_enriched_interpolation_cells();

            // number of cells
            moris_index tNumCells = this->get_num_entities( EntityRank::ELEMENT );

            MORIS_ASSERT( tNumCells == (moris_index)tEnrIpCells.size(), "Inconsistent num cells information." );

            // create the enriched interpolation basis to interpolation cell interpolation
            for ( moris::sint i = 0; i < tNumCells; i++ )
            {
                moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const & tVertices = tEnrIpCells( i )->get_xtk_interpolation_vertices();

                for ( moris::size_t iV = 0; iV < tVertices.size(); iV++ )
                {
                    if ( tVertices( iV )->has_interpolation( tMeshIndex ) )
                    {
                        Vertex_Enrichment* tVertInterp = tVertices( iV )->get_xtk_interpolation( tMeshIndex );

                        Matrix< IndexMat > tBasisIndices = tVertInterp->get_indices();

                        // iterate through vertices and add the ip cell to the support
                        for ( moris::uint iB = 0; iB < tBasisIndices.numel(); iB++ )
                        {
                            tCellsInEnrSupports( tBasisIndices( iB ) ).push_back( tEnrIpCells( i ) );
                        }
                    }
                }
            }

            // iterate through enriched basis functions
            for ( moris::size_t i = 0; i < tNumEnrBasis; i++ )
            {
                // iterate through enriched interpolation cells in the support
                moris_index tNumSubphaseInSupport = tCellsInEnrSupports( i ).size();

                moris_index tExpectedBulkPhase = MORIS_INDEX_MAX;

                for ( moris::moris_index iSP = 0; iSP < tNumSubphaseInSupport; iSP++ )
                {
                    Interpolation_Cell_Unzipped* tIpCell = tCellsInEnrSupports( i )( iSP );

                    moris_index tBulkPhase = tIpCell->get_bulkphase_index();

                    MORIS_ASSERT( tBulkPhase != MORIS_INDEX_MAX, "Bulk phase index not set." );

                    if ( iSP == 0 )
                    {
                        tExpectedBulkPhase = tBulkPhase;
                    }

                    if ( tBulkPhase != tExpectedBulkPhase )
                    {
                        std::cout << "tExpectedBulkPhase=  " << tExpectedBulkPhase << " | tBulkPhase = " << tBulkPhase << std::endl;
                        tSubphaseBulkPhasesInSupportDiag = false;
                    }
                }
            }
        }

        MORIS_ERROR( tSubphaseBulkPhasesInSupportDiag, "All bulk phases in support of enriched basis  function do not match" );
        return true;
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_local_to_global_maps()
    {
        // initialize local to global maps
        mLocalToGlobalMaps = Cell< Matrix< IdMat > >( 4 );
        mGlobaltoLobalMaps = Cell< std::unordered_map< moris_id, moris_index > >( 4 );

        this->setup_cell_maps();

        this->setup_basis_maps();

        this->assign_ip_vertex_ids();

        this->setup_vertex_maps();
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::assign_ip_vertex_ids()
    {
        // mpitag
        moris_index tTag = 600001;

        // owned requests and shared requests sorted by owning proc
        Cell< uint >         tOwnedVertices;
        Cell< Cell< uint > > tNotOwnedVertices;
        Cell< Cell< uint > > tNotOwnedIpCells;
        Cell< uint >         tProcRanks;

        std::unordered_map< moris_id, moris_id > tProcRankToDataIndex;
        this->sort_ip_vertices_by_owned_and_not_owned(
                tOwnedVertices,
                tNotOwnedVertices,
                tNotOwnedIpCells,
                tProcRanks,
                tProcRankToDataIndex );

        moris::moris_id tVertexId = this->allocate_entity_ids( tOwnedVertices.size(), EntityRank::NODE, true );

        // Assign owned request identifiers
        this->assign_owned_ip_vertex_ids( tOwnedVertices, tVertexId );

        // prepare node information request data
        Cell< Matrix< IndexMat > > tOutwardBaseVertexIds;
        Cell< Matrix< IndexMat > > tOutwardIpCellIds;
        this->setup_outward_ip_vertex_requests(
                tNotOwnedVertices,
                tNotOwnedIpCells,
                tProcRanks,
                tProcRankToDataIndex,
                tOutwardBaseVertexIds,
                tOutwardIpCellIds );

        // send requests to owning processor
        mXTKModel->send_outward_requests( tTag, tProcRanks, tOutwardBaseVertexIds );
        mXTKModel->send_outward_requests( tTag + 1, tProcRanks, tOutwardIpCellIds );

        // hold on to make sure everyone has sent all their information
        barrier();

        // receive the requests
        Cell< Matrix< IndexMat > > tReceivedBaseVertexIds;
        Cell< Matrix< IndexMat > > tReceivedIpCellIds;

        Cell< uint > tProcsReceivedFrom1;
        Cell< uint > tProcsReceivedFrom2;
        mXTKModel->inward_receive_requests(
                tTag, 1, tReceivedBaseVertexIds, tProcsReceivedFrom1 );

        mXTKModel->inward_receive_requests(
                tTag + 1, 1, tReceivedIpCellIds, tProcsReceivedFrom2 );

        MORIS_ASSERT( tProcsReceivedFrom1.size() == tProcsReceivedFrom2.size(),
                "Size mismatch between procs received from child cell ids and number of child cells" );

        Cell< Matrix< IndexMat > > tVertexIdsAnswers;
        this->prepare_ip_vertex_id_answers( tReceivedBaseVertexIds, tReceivedIpCellIds, tVertexIdsAnswers );

        // return information
        mXTKModel->return_request_answers( tTag + 2, tVertexIdsAnswers, tProcsReceivedFrom1 );

        // receive the information
        barrier();

        // receive the answers
        Cell< Matrix< IndexMat > > tReceivedVertexIds;

        mXTKModel->inward_receive_request_answers( tTag + 2, 1, tProcRanks, tReceivedVertexIds );

        // add child cell ids to not owned child meshes
        this->handle_received_ip_vertex_ids( tNotOwnedVertices, tReceivedVertexIds );

        barrier();
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::sort_ip_vertices_by_owned_and_not_owned(
            Cell< uint >&                             aOwnedVertices,
            Cell< Cell< uint > >&                     aNotOwnedVertices,
            Cell< Cell< uint > >&                     aNotOwnedIPCells,
            Cell< uint >&                             aProcRanks,
            std::unordered_map< moris_id, moris_id >& aProcRankToIndexInData )
    {
        // access the mesh data behind the background mesh
        // access the communication
        Matrix< IdMat > tCommTable = this->get_communication_table();

        // Par rank
        moris::moris_index tParRank = par_rank();

        // resize proc ranks and setup map to comm table
        aProcRanks.resize( tCommTable.numel() );

        for ( moris::uint i = 0; i < tCommTable.numel(); i++ )
        {
            aProcRankToIndexInData[tCommTable( i )] = i;

            aProcRanks( i ) = ( tCommTable( i ) );

            aNotOwnedVertices.push_back( Cell< uint >( 0 ) );
            aNotOwnedIPCells.push_back( Cell< uint >( 0 ) );
        }

        // Number of cells
        uint tNumCells = this->get_num_entities( EntityRank::ELEMENT );

        // number of nodes
        uint tNumVerts = this->get_num_entities( EntityRank::NODE );

        // Keep track of nodes which have ids
        Cell< uint > tNodeTracker( tNumVerts, 0 );

        // Reserve memory for indices of owned vertices
        aOwnedVertices.reserve( tNumVerts );

        // iterate through cells
        for ( moris::uint iC = 0; iC < tNumCells; iC++ )
        {
            Interpolation_Cell_Unzipped const & tCell = *mEnrichedInterpCells( (moris_index)iC );

            moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const & tVerts =
                    tCell.get_xtk_interpolation_vertices();

            // iterate through vertices
            for ( moris::uint iV = 0; iV < tVerts.size(); iV++ )
            {
                moris_index tVertexIndex = tVerts( iV )->get_index();
                moris_index tOwner       = tVerts( iV )->get_owner();

                if ( tNodeTracker( tVertexIndex ) == 0 )
                {
                    if ( tOwner == tParRank )
                    {
                        aOwnedVertices.push_back( tVertexIndex );
                    }
                    else
                    {
                        moris_index tIndex = aProcRankToIndexInData[tOwner];

                        aNotOwnedVertices( tIndex ).push_back( tVertexIndex );
                        aNotOwnedIPCells( tIndex ).push_back( tCell.get_index() );
                    }
                    tNodeTracker( tVertexIndex ) = 1;
                }
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::assign_owned_ip_vertex_ids(
            Cell< uint > const & aOwnedIpVerts,
            moris::moris_id&     aNodeId )
    {
        // iterate through vertices that i own and assign a node id to them
        for ( moris::uint i = 0; i < aOwnedIpVerts.size(); i++ )
        {
            mEnrichedInterpVerts( aOwnedIpVerts( i ) )->set_vertex_id( aNodeId );
            aNodeId++;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_outward_ip_vertex_requests(
            Cell< Cell< uint > > const &              aNotOwnedIpVerts,
            Cell< Cell< uint > > const &              aNotOwnedIpCells,
            Cell< uint > const &                      aProcRanks,
            std::unordered_map< moris_id, moris_id >& aProcRankToIndexInData,
            Cell< Matrix< IndexMat > >&               aOutwardBaseVertexIds,
            Cell< Matrix< IndexMat > >&               aOutwardIpCellIds )
    {
        MORIS_ASSERT( aNotOwnedIpVerts.size() == aNotOwnedIpCells.size(), "Dimension mismatch" );

        aOutwardBaseVertexIds.resize( aNotOwnedIpVerts.size() );
        aOutwardIpCellIds.resize( aNotOwnedIpVerts.size() );

        // iterate through the verts and collect them for communication
        for ( moris::uint iP = 0; iP < aNotOwnedIpVerts.size(); iP++ )
        {

            aOutwardBaseVertexIds( iP ).resize( 1, aNotOwnedIpVerts( iP ).size() );
            aOutwardIpCellIds( iP ).resize( 1, aNotOwnedIpVerts( iP ).size() );

            for ( moris::uint iV = 0; iV < aNotOwnedIpVerts( iP ).size(); iV++ )
            {
                moris_index tVertexIndex = aNotOwnedIpVerts( iP )( iV );
                moris_index tCellIndex   = aNotOwnedIpCells( iP )( iV );

                aOutwardBaseVertexIds( iP )( iV ) = mEnrichedInterpVerts( tVertexIndex )->get_base_vertex()->get_id();
                aOutwardIpCellIds( iP )( iV )     = mEnrichedInterpCells( tCellIndex )->get_id();
            }

            if ( aNotOwnedIpVerts( iP ).size() == 0 )
            {
                aOutwardBaseVertexIds( iP ).resize( 1, 1 );
                aOutwardIpCellIds( iP ).resize( 1, 1 );

                aOutwardBaseVertexIds( iP )( 0 ) = MORIS_INDEX_MAX;
                aOutwardIpCellIds( iP )( 0 )     = MORIS_INDEX_MAX;
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::prepare_ip_vertex_id_answers(
            Cell< Matrix< IndexMat > >& aReceivedBaseVertexIds,
            Cell< Matrix< IndexMat > >& aReceivedIpCellIds,
            Cell< Matrix< IndexMat > >& aVertexIdAnswer )
    {
        MORIS_ASSERT( aReceivedBaseVertexIds.size() == aReceivedIpCellIds.size(),
                "Mismatch in received base vertex ids and received ip cell" );

        // allocate answer size
        aVertexIdAnswer.resize( aReceivedBaseVertexIds.size() );

        // iterate through received data
        for ( moris::uint i = 0; i < aReceivedBaseVertexIds.size(); i++ )
        {
            uint tNumReceivedReqs = aReceivedBaseVertexIds( i ).n_cols();

            aVertexIdAnswer( i ).resize( 1, tNumReceivedReqs );

            if ( aReceivedBaseVertexIds( i )( 0 ) != MORIS_INDEX_MAX )
            {
                // iterate through received requests
                for ( moris::uint j = 0; j < tNumReceivedReqs; j++ )
                {
                    bool tSuccess = false;

                    // base vertex id
                    moris_index tBaseVertexId = aReceivedBaseVertexIds( i )( j );

                    // ip cell id
                    moris_index tIpCellId = aReceivedIpCellIds( i )( j );

                    // ip cell index
                    moris_index tIpCellInd = this->get_loc_entity_ind_from_entity_glb_id( tIpCellId, EntityRank::ELEMENT );

                    // ip cell
                    Interpolation_Cell_Unzipped* tIpCell = mEnrichedInterpCells( tIpCellInd );

                    // vertices attached to the ip cell
                    moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const & tVerts = tIpCell->get_xtk_interpolation_vertices();

                    for ( moris::uint iV = 0; iV < tVerts.size(); iV++ )
                    {
                        if ( tVerts( iV )->get_base_vertex()->get_id() == tBaseVertexId )
                        {
                            tSuccess                  = true;
                            aVertexIdAnswer( i )( j ) = tVerts( iV )->get_id();
                        }
                    }

                    if ( tSuccess == false )
                    {
                        MORIS_ASSERT( false, "Vertex Not Found" );
                    }
                }
            }
            else
            {
                aVertexIdAnswer( i )( 0 ) = MORIS_INDEX_MAX;
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::handle_received_ip_vertex_ids(
            Cell< Cell< uint > > const &       aNotOwnedVertices,
            Cell< Matrix< IndexMat > > const & aReceivedVertexIds )
    {
        // iterate through received data
        for ( moris::uint i = 0; i < aNotOwnedVertices.size(); i++ )
        {
            uint tNumReceivedReqs = aNotOwnedVertices( i ).size();

            // iterate through received requests
            for ( moris::uint j = 0; j < tNumReceivedReqs; j++ )
            {
                mEnrichedInterpVerts( aNotOwnedVertices( i )( j ) )->set_vertex_id( aReceivedVertexIds( i )( j ) );
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::communicate_select_vertex_interpolation(
            moris::Cell< mtk::Vertex* >& aVerticesToCommunicate )
    {
        Tracer tTracer( "XTK", "Enriched_Interpolation_Mesh", "communicate_select_vertex_interpolation", mXTKModel->mVerboseLevel, 0 );

        // get communication table
        Matrix< IndexMat > tCommTable = this->get_communication_table();

        // build maps for communication
        std::unordered_map< moris_index, moris_index > tCommMap;
        Cell< uint >                                   tProcs( tCommTable.numel() );

        for ( moris::uint i = 0; i < tCommTable.numel(); i++ )
        {
            tCommMap[tCommTable( i )] = i;
            tProcs( i )               = tCommTable( i );
        }

        // count vertices to be communicated with other processors
        Cell< moris::uint > tCounts( tCommTable.numel() );

        for ( const auto& iVert : aVerticesToCommunicate )
        {
            MORIS_ASSERT( iVert->get_owner() != moris::par_rank(),
                    "Enriched_Interpolation_Mesh::communicate_select_vertex_interpolation - Requested communication of vertex that current proc owns" );

            tCounts( tCommMap[iVert->get_owner()] )++;
        }

        // set size for data to be communicated
        Cell< Matrix< IndexMat > > tVertexIdsToProcs( tCommTable.numel() );

        for ( moris::uint iProcs = 0; iProcs < tCommTable.numel(); iProcs++ )
        {
            if ( tCounts( iProcs ) == 0 )
            {
                tVertexIdsToProcs( iProcs ).resize( 1, 1 );
            }
            else
            {
                tVertexIdsToProcs( iProcs ).resize( 1, tCounts( iProcs ) );
            }

            tVertexIdsToProcs( iProcs ).fill( MORIS_INDEX_MAX );
        }

        // populate data to be communicated with vertex ids
        Cell< moris::uint > tCurrentCounts( tCommTable.numel(), 0 );

        for ( const auto& iVert : aVerticesToCommunicate )
        {
            // get processor index
            const moris_index tProcOrd = tCommMap[iVert->get_owner()];

            // get current number of vertices for this processor
            const moris::uint tCount = tCurrentCounts( tProcOrd );

            // set vertex id
            tVertexIdsToProcs( tProcOrd )( tCount ) = iVert->get_id();

            // increase vertex counter for this processor
            tCurrentCounts( tProcOrd )++;
        }

        // iterate through all bspline meshes
        for ( moris::uint iMeshIndex = 0; iMeshIndex < mMeshIndices.numel(); iMeshIndex++ )
        {
            // current mesh index
            moris_index tMeshIndex = mMeshIndices( iMeshIndex );

            // set communication tag
            moris::uint tMPITag = 5001;

            // FIXME: the following communication needs to be refactored to properly use non-blocking sends and receives

            // send list with vertex ids to communicating processors
            mXTKModel->send_outward_requests( tMPITag, tProcs, tVertexIdsToProcs );

            Cell< Matrix< IndexMat > > tReceivedVertexIds;
            Cell< uint >               tProcsReceivedFrom;
            mXTKModel->inward_receive_requests( tMPITag, 1, tReceivedVertexIds, tProcsReceivedFrom );

            // prepare the t-matrices for sending
            Cell< Matrix< DDRMat > >   tTMatrixWeights;
            Cell< Matrix< IndexMat > > tTMatrixIndices;
            Cell< Matrix< IndexMat > > tTMatrixOwners;
            Cell< Matrix< IndexMat > > tTMatrixOffsets;

            this->prepare_t_matrix_request_answers(
                    tMeshIndex,
                    tReceivedVertexIds,
                    tTMatrixWeights,
                    tTMatrixIndices,
                    tTMatrixOwners,
                    tTMatrixOffsets );

            // send information
            mXTKModel->return_request_answers_reals( tMPITag + 3, tTMatrixWeights, tProcsReceivedFrom );
            mXTKModel->return_request_answers( tMPITag + 4, tTMatrixIndices, tProcsReceivedFrom );
            mXTKModel->return_request_answers( tMPITag + 5, tTMatrixOwners, tProcsReceivedFrom );
            mXTKModel->return_request_answers( tMPITag + 6, tTMatrixOffsets, tProcsReceivedFrom );

            // wait
            barrier();

            // receive
            Cell< Matrix< DDRMat > >   tRequestedTMatrixWeights;
            Cell< Matrix< IndexMat > > tRequestedTMatrixIndices;
            Cell< Matrix< IndexMat > > tRequestedTMatrixOwners;
            Cell< Matrix< IndexMat > > tRequestedTMatrixOffsets;

            // receive the answers
            mXTKModel->inward_receive_request_answers_reals( tMPITag + 3, 1, tProcs, tRequestedTMatrixWeights );
            mXTKModel->inward_receive_request_answers( tMPITag + 4, 1, tProcs, tRequestedTMatrixIndices );
            mXTKModel->inward_receive_request_answers( tMPITag + 5, 1, tProcs, tRequestedTMatrixOwners );
            mXTKModel->inward_receive_request_answers( tMPITag + 6, 1, tProcs, tRequestedTMatrixOffsets );

            barrier();

            // commit it to my data
            this->handle_received_interpolation_data(
                    tMeshIndex,
                    tVertexIdsToProcs,
                    tRequestedTMatrixWeights,
                    tRequestedTMatrixIndices,
                    tRequestedTMatrixOwners,
                    tRequestedTMatrixOffsets );

            // wait
            barrier();
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::prepare_t_matrix_request_answers(
            moris_index const &                aMeshIndex,
            Cell< Matrix< IndexMat > > const & aRequestedEnrIPVertexIds,
            Cell< Matrix< DDRMat > >&          aTMatrixWeights,
            Cell< Matrix< IndexMat > >&        aTMatrixIndices,
            Cell< Matrix< IndexMat > >&        aBasisOwners,
            Cell< Matrix< IndexMat > >&        aTMatrixOffsets )
    {

        // information about size of interpolation mats
        Cell< uint > tSizes( aRequestedEnrIPVertexIds.size(), 0 );

        // resize the input data
        aTMatrixWeights.resize( aRequestedEnrIPVertexIds.size() );
        aTMatrixIndices.resize( aRequestedEnrIPVertexIds.size() );
        aBasisOwners.resize( aRequestedEnrIPVertexIds.size() );
        aTMatrixOffsets.resize( aRequestedEnrIPVertexIds.size() );

        // collect the vertex interpolations
        Cell< Cell< Vertex_Enrichment* > > tVertexInterpolations( aRequestedEnrIPVertexIds.size() );

        // collect size information throughout loop
        Cell< moris_index > tDataSizes( aRequestedEnrIPVertexIds.size(), 0 );

        // iterate through and figure out how big to make the weights and indices mats
        // also collect vertex interpolations
        for ( moris::uint iP = 0; iP < aRequestedEnrIPVertexIds.size(); iP++ )
        {
            // no information requested prepare a dummy response
            if ( aRequestedEnrIPVertexIds( iP ).numel() == 1 and aRequestedEnrIPVertexIds( iP )( 0 ) == MORIS_INDEX_MAX )
            {
                aTMatrixWeights( iP ).resize( 1, 1 );
                aTMatrixIndices( iP ).resize( 1, 1 );
                aBasisOwners( iP ).resize( 1, 1 );

                aTMatrixWeights( iP )( 0 ) = MORIS_REAL_MAX;
                aTMatrixIndices( iP )( 0 ) = MORIS_INDEX_MAX;
                aBasisOwners( iP )( 0 )    = MORIS_INDEX_MAX;
                continue;
            }

            // size the tmatrix offset for each vertex requested (num verts +1)
            aTMatrixOffsets( iP ).resize( 1, aRequestedEnrIPVertexIds( iP ).numel() + 1 );

            // set the first one to 0
            aTMatrixOffsets( iP )( 0 ) = 0;

            // iterate through the vertices and get their interpolations and figure out
            // how big it is
            for ( moris::uint iV = 0; iV < aRequestedEnrIPVertexIds( iP ).numel(); iV++ )
            {

                moris_index tVertexIndex = this->get_loc_entity_ind_from_entity_glb_id( aRequestedEnrIPVertexIds( iP )( iV ), EntityRank::NODE );

                Interpolation_Vertex_Unzipped* tVertex = this->get_unzipped_vertex_pointer( tVertexIndex );

                MORIS_ERROR( tVertex->get_id() == aRequestedEnrIPVertexIds( iP )( iV ), "ID Mismatch" );
                MORIS_ERROR( tVertex->get_owner() == par_rank(), "Must be a vertex ownership issue." );

                // get the vertex interpolation
                Vertex_Enrichment* tVertexInterp = tVertex->get_xtk_interpolation( aMeshIndex );

                MORIS_ASSERT( tVertexInterp->get_base_vertex_interpolation() != nullptr, "Owning proc has a nullptr for the vertex interpolation." );

                tVertexInterpolations( iP ).push_back( tVertexInterp );

                // number of basis functions interpolating into this vertex
                moris_index tNumBasis = tVertexInterp->get_basis_indices().numel();

                // offsets
                aTMatrixOffsets( iP )( iV + 1 ) = aTMatrixOffsets( iP )( iV ) + tNumBasis;

                // add to size
                tDataSizes( iP ) = tDataSizes( iP ) + tNumBasis;
            }
        }

        //  iterate through and size data
        for ( moris::uint iP = 0; iP < aRequestedEnrIPVertexIds.size(); iP++ )
        {
            if ( aRequestedEnrIPVertexIds( iP ).numel() == 1 and aRequestedEnrIPVertexIds( iP )( 0 ) == MORIS_INDEX_MAX )
            {
                continue;
            }

            aTMatrixWeights( iP ).resize( 1, tDataSizes( iP ) );
            aTMatrixIndices( iP ).resize( 1, tDataSizes( iP ) );
            aBasisOwners( iP ).resize( 1, tDataSizes( iP ) );
        }

        // populate the data
        for ( moris::uint iP = 0; iP < aRequestedEnrIPVertexIds.size(); iP++ )
        {

            if ( aRequestedEnrIPVertexIds( iP ).numel() == 1 and aRequestedEnrIPVertexIds( iP )( 0 ) == MORIS_INDEX_MAX )
            {
                aTMatrixOffsets( iP ).resize( 1, 1 );
                aTMatrixOffsets( iP )( 0 ) = MORIS_INDEX_MAX;
                continue;
            }

            moris::uint tCount = 0;

            for ( moris::uint iV = 0; iV < tVertexInterpolations( iP ).size(); iV++ )
            {
                this->add_vertex_interpolation_to_communication_data(
                        tCount,
                        tVertexInterpolations( iP )( iV ),
                        aTMatrixWeights( iP ),
                        aTMatrixIndices( iP ),
                        aBasisOwners( iP ),
                        aTMatrixOffsets( iP ) );
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::add_vertex_interpolation_to_communication_data(
            moris::uint&        aCount,
            Vertex_Enrichment*  aInterpolation,
            Matrix< DDRMat >&   aTMatrixWeights,
            Matrix< IndexMat >& aTMatrixIndices,
            Matrix< IndexMat >& aTMatrixOwners,
            Matrix< IndexMat >& aTMatrixOffsets )
    {
        // access the basis indices and weights
        moris::Matrix< moris::IndexMat > const & tBasisIndices = aInterpolation->get_basis_ids();
        moris::Matrix< moris::DDRMat > const &   tBasisWeights = aInterpolation->get_basis_weights();
        moris::Matrix< moris::IndexMat >         tBasisOwners  = aInterpolation->get_owners();

        for ( moris::uint i = 0; i < tBasisIndices.numel(); i++ )
        {
            aTMatrixIndices( aCount ) = tBasisIndices( i );
            aTMatrixWeights( aCount ) = tBasisWeights( i );

            MORIS_ASSERT( tBasisOwners( i ) < par_size() || tBasisOwners( i ) > 0, "Bad ownership for basis function." );

            aTMatrixOwners( aCount ) = tBasisOwners( i );
            aCount++;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::extract_vertex_interpolation_from_communication_data(
            moris::uint const &         aNumVerts,
            Matrix< DDRMat > const &    aTMatrixWeights,
            Matrix< IndexMat > const &  aTMatrixIndices,
            Matrix< IndexMat > const &  aTMatrixOwners,
            Matrix< IndexMat > const &  aTMatrixOffsets,
            Cell< Matrix< DDRMat > >&   aExtractedTMatrixWeights,
            Cell< Matrix< IndexMat > >& aExtractedTMatrixIndices,
            Cell< Matrix< IndexMat > >& aExtractedBasisOwners )
    {
        // size output data
        aExtractedTMatrixWeights.resize( aNumVerts );
        aExtractedTMatrixIndices.resize( aNumVerts );
        aExtractedBasisOwners.resize( aNumVerts );

        // current starting index
        moris_index tStart = 0;

        // extract the data into the cells
        for ( moris::uint iV = 0; iV < aNumVerts; iV++ )
        {
            // number of basis interpolating into the vertex
            moris::moris_index tNumBasis = aTMatrixOffsets( iV + 1 ) - tStart;

            aExtractedTMatrixWeights( iV ).resize( tNumBasis, 1 );
            aExtractedTMatrixIndices( iV ).resize( tNumBasis, 1 );
            aExtractedBasisOwners( iV ).resize( 1, tNumBasis );

            // itere and grab  data
            for ( moris::moris_index iIp = 0; iIp < tNumBasis; iIp++ )
            {
                aExtractedTMatrixWeights( iV )( iIp ) = aTMatrixWeights( tStart + iIp );
                aExtractedTMatrixIndices( iV )( iIp ) = aTMatrixIndices( tStart + iIp );
                aExtractedBasisOwners( iV )( iIp )    = aTMatrixOwners( tStart + iIp );
            }

            tStart = aTMatrixOffsets( iV + 1 );
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_basis_ownership()
    {
        // size data
        mEnrichCoeffOwnership.resize( mMeshIndices.max() + 1 );

        // iterate through meshes
        for ( moris::uint iM = 0; iM < mMeshIndices.numel(); iM++ )
        {
            moris_index tMeshIndex = mMeshIndices( iM );

            mEnrichCoeffOwnership( tMeshIndex ).resize( 1, mEnrichCoeffLocToGlob( tMeshIndex ).numel() );

            mEnrichCoeffOwnership( tMeshIndex ).fill( MORIS_INDEX_MAX );

            // iterate through basis functions
            for ( moris::uint iB = 0; iB < mCoeffToEnrichCoeffs( tMeshIndex ).size(); iB++ )
            {
                moris_index tOwner = mXTKModel->get_background_mesh().get_entity_owner( (moris_index)iB, mBasisRank, tMeshIndex );

                for ( moris::uint iEB = 0; iEB < mCoeffToEnrichCoeffs( tMeshIndex )( iB ).numel(); iEB++ )
                {
                    moris_index tEnrIndex                            = mCoeffToEnrichCoeffs( tMeshIndex )( iB )( iEB );
                    mEnrichCoeffOwnership( tMeshIndex )( tEnrIndex ) = tOwner;

                    if ( tOwner != par_rank() )
                    {
                        mNotOwnedBasis.push_back( mEnrichCoeffLocToGlob( tMeshIndex )( tEnrIndex ) );
                    }
                    else
                    {
                        mOwnedBasis.push_back( mEnrichCoeffLocToGlob( tMeshIndex )( tEnrIndex ) );
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::handle_received_interpolation_data(
            moris_index const &                aMeshIndex,
            Cell< Matrix< IndexMat > > const & aVertexIdsToProc,
            Cell< Matrix< DDRMat > > const &   aRequestedTMatrixWeights,
            Cell< Matrix< IndexMat > > const & aRequestedTMatrixIndices,
            Cell< Matrix< IndexMat > > const & aRequestedBasisOwners,
            Cell< Matrix< IndexMat > > const & aRequestedTMatrixOffsets )
    {
        // access the communication
        Matrix< IdMat > tCommTable = this->get_communication_table();

        std::unordered_map< moris_id, moris_id > tProcRankToIndexInData;

        moris::uint tCount = tCommTable.numel();

        // resize proc ranks and setup map to comm table
        for ( moris::uint i = 0; i < tCommTable.numel(); i++ )
        {
            tProcRankToIndexInData[tCommTable( i )] = i;
        }

        // iterate through returned information
        for ( moris::uint iP = 0; iP < aVertexIdsToProc.size(); iP++ )
        {
            // MORIS_LOG_SPEC("aVertexIdsToProc(iP).numel() on " + std::to_string(par_rank()), aVertexIdsToProc(iP).numel());

            if ( aVertexIdsToProc( iP ).numel() == 1 and aVertexIdsToProc( iP )( 0 ) == MORIS_INDEX_MAX )
            {
                // do nothing for this iP
            }
            // standard case
            else
            {
                // extract the t-matrices and basis ids/owners for the proc ip
                Cell< Matrix< DDRMat > >   tExtractedTMatrixWeights;
                Cell< Matrix< IndexMat > > tExtractedTMatrixIds;
                Cell< Matrix< IndexMat > > tExtractedTBasisOwners;

                this->extract_vertex_interpolation_from_communication_data(
                        aVertexIdsToProc( iP ).numel(),
                        aRequestedTMatrixWeights( iP ),
                        aRequestedTMatrixIndices( iP ),
                        aRequestedBasisOwners( iP ),
                        aRequestedTMatrixOffsets( iP ),
                        tExtractedTMatrixWeights,
                        tExtractedTMatrixIds,
                        tExtractedTBasisOwners );

                // verify consistent sizes
                MORIS_ASSERT( aVertexIdsToProc( iP ).numel() == tExtractedTMatrixWeights.size(),
                        "Size mismatch in t-matrix weights." );
                MORIS_ASSERT( aVertexIdsToProc( iP ).numel() == tExtractedTMatrixIds.size(),
                        "Size mismatch in t-matrix ids." );
                MORIS_ASSERT( aVertexIdsToProc( iP ).numel() == tExtractedTBasisOwners.size(),
                        "Size mismatch in basis owners." );

                // iterate through vertices and set their interpolation weights and basis ids
                for ( moris::uint iV = 0; iV < aVertexIdsToProc( iP ).numel(); iV++ )
                {
                    // get the vertex
                    moris_index tVertexIndex = this->get_loc_entity_ind_from_entity_glb_id( aVertexIdsToProc( iP )( iV ), EntityRank::NODE );

                    Interpolation_Vertex_Unzipped& tVertex = this->get_xtk_interp_vertex( tVertexIndex );

                    // get the enriched vertex interpolation
                    Vertex_Enrichment* tVertexInterp = tVertex.get_xtk_interpolation( aMeshIndex );

                    // iterate through basis functions and find local indices
                    moris::Matrix< IndexMat > tBasisIndices( tExtractedTMatrixIds( iV ).n_rows(), tExtractedTMatrixIds( iV ).n_cols() );

                    for ( moris::uint iBs = 0; iBs < tExtractedTMatrixIds( iV ).numel(); iBs++ )
                    {
                        // basis id
                        moris_id tId = tExtractedTMatrixIds( iV )( iBs );

                        // add this basis to the mesh if it doesnt exists on the current partition
                        if ( !this->basis_exists_on_partition( aMeshIndex, tId ) )
                        {
                            MORIS_ASSERT( tExtractedTBasisOwners( iV )( iBs ) != par_rank(), "Owned basis should already exist on partition." );

                            this->add_basis_function( aMeshIndex, tId, tExtractedTBasisOwners( iV )( iBs ), MORIS_INDEX_MAX );
                        }

                        tBasisIndices( iBs ) = this->get_enr_basis_index_from_enr_basis_id( aMeshIndex, tId );

                        moris_id tBasisOwner = tExtractedTBasisOwners( iV )( iBs );

                        MORIS_ASSERT( this->get_basis_owner( tBasisIndices( iBs ), aMeshIndex ) == tBasisOwner, "Ownership discrepancy." );

                        // if the basis has an owning proc that is not in the comm table, add it to the comm table
                        if ( tProcRankToIndexInData.find( tBasisOwner ) == tProcRankToIndexInData.end() && tBasisOwner != par_rank() )
                        {
                            this->add_proc_to_comm_table( tBasisOwner );
                            tProcRankToIndexInData[tBasisOwner] = tCount;
                            tCount++;
                        }
                    }

                    // iterate through basis in the base vertex interpolation
                    moris::uint tNumCoeffs = tExtractedTMatrixIds( iV ).numel();

                    // Setup the map in the basis function
                    std::unordered_map< moris::moris_index, moris::moris_index >& tVertEnrichMap = tVertexInterp->get_basis_map();

                    for ( moris::uint iB = 0; iB < tNumCoeffs; iB++ )
                    {
                        moris::moris_index tBasisIndex = tBasisIndices( iB );

                        tVertEnrichMap[tBasisIndex] = iB;
                    }

                    // get the basis indices from the basis ids
                    tVertexInterp->add_basis_information( tBasisIndices, tExtractedTMatrixIds( iV ) );
                    tVertexInterp->add_basis_weights( tBasisIndices, tExtractedTMatrixWeights( iV ) );
                    tVertexInterp->add_basis_owners( tBasisIndices, tExtractedTBasisOwners( iV ) );
                    tVertexInterp->add_base_vertex_interpolation( nullptr );// base vertex interpolation does not exists (other  proc)
                }
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_vertex_maps()
    {
        moris::uint tNumNodes = this->get_num_entities( EntityRank::NODE );

        mLocalToGlobalMaps( 0 ) = Matrix< IdMat >( tNumNodes, 1 );

        for ( moris::uint i = 0; i < tNumNodes; i++ )
        {
            mLocalToGlobalMaps( 0 )( mEnrichedInterpVerts( i )->get_index() ) = mEnrichedInterpVerts( i )->get_id();

            MORIS_ASSERT( mEnrichedInterpVerts( i )->get_index() == (moris_index)i, "Index alignment issue in vertices" );

            MORIS_ASSERT( mGlobaltoLobalMaps( 0 ).find( mEnrichedInterpVerts( i )->get_id() ) == mGlobaltoLobalMaps( 0 ).end(),
                    "Duplicate id in the vertex map detected" );

            mGlobaltoLobalMaps( 0 )[mEnrichedInterpVerts( i )->get_id()] = mEnrichedInterpVerts( i )->get_index();
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_basis_maps()
    {
        mGlobaltoLocalBasisMaps.resize( mMeshIndices.max() + 1 );

        // iterate through meshes
        for ( moris::uint iM = 0; iM < mEnrichCoeffLocToGlob.size(); iM++ )
        {
            for ( moris::uint iB = 0; iB < mEnrichCoeffLocToGlob( iM ).numel(); iB++ )
            {
                // MORIS_LOG_SPEC("mEnrichCoeffLocToGlob(iM)(iB)",mEnrichCoeffLocToGlob(iM)(iB));
                // MORIS_ASSERT(mGlobaltoLocalBasisMaps(iM).find(mEnrichCoeffLocToGlob(iM)(iB)) == mGlobaltoLocalBasisMaps(iM).end(),
                //         "Duplicate id in the basis map detected");

                if ( mGlobaltoLocalBasisMaps( iM ).find( mEnrichCoeffLocToGlob( iM )( iB ) ) == mGlobaltoLocalBasisMaps( iM ).end() )
                {
                    mGlobaltoLocalBasisMaps( iM )[mEnrichCoeffLocToGlob( iM )( iB )] = (moris_index)iB;
                    MORIS_ASSERT( this->get_enr_basis_index_from_enr_basis_id( iM, mEnrichCoeffLocToGlob( iM )( iB ) ) == (moris_index)iB,
                            "Issue setting up the basis map" );
                }
                else
                {
                    // MORIS_ASSERT(this->get_basis_owner(mEnrichCoeffLocToGlob(iM)(iB),iM) != par_rank(),"Merging basis required for owned basis.");
                }
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_cell_maps()
    {
        moris::uint tNumCells = this->get_num_entities( EntityRank::ELEMENT );

        mLocalToGlobalMaps( 3 ) = Matrix< IdMat >( tNumCells, 1 );

        //        mGlobaltoLobalMaps(3).clear();

        for ( moris::uint i = 0; i < tNumCells; i++ )
        {
            mLocalToGlobalMaps( 3 )( mEnrichedInterpCells( i )->get_index() ) = mEnrichedInterpCells( i )->get_id();

            MORIS_ASSERT( mEnrichedInterpCells( i )->get_index() == (moris_index)i, "Index alignment issue in cells" );

            MORIS_ASSERT( mGlobaltoLobalMaps( 3 ).find( mEnrichedInterpCells( i )->get_id() ) == mGlobaltoLobalMaps( 3 ).end(),
                    "Duplicate id in the cell map detected" );

            mGlobaltoLobalMaps( 3 )[mEnrichedInterpCells( i )->get_id()] = mEnrichedInterpCells( i )->get_index();
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_mesh_index_map()
    {
        for ( moris::uint i = 0; i < (uint)mMeshIndices.max() + 1; i++ )
        {

            mMeshIndexToLocMeshIndex[i] = i;
        }
    }

    // ----------------------------------------------------------------------------

    moris_id
    Enriched_Interpolation_Mesh::allocate_entity_ids(
            moris::size_t   aNumReqs,
            enum EntityRank aEntityRank,
            bool            aStartFresh )
    {
        MORIS_ASSERT( aEntityRank == EntityRank::NODE || aEntityRank == EntityRank::ELEMENT,
                "Only Elements or Nodes have ids" );

        moris_id tGlobalMax = 1;
        if ( !aStartFresh )
        {
            this->get_max_entity_id( aEntityRank );
        }

        int tProcRank = par_rank();
        int tProcSize = par_size();

        moris::Cell< moris::moris_id > aGatheredInfo;
        moris::Cell< moris::moris_id > tFirstId( 1 );
        moris::Cell< moris::moris_id > tNumIdsRequested( 1 );

        tNumIdsRequested( 0 ) = (moris::moris_id)aNumReqs;

        moris::gather( tNumIdsRequested, aGatheredInfo );

        moris::Cell< moris::moris_id > tProcFirstID( tProcSize );

        if ( tProcRank == 0 )
        {
            // Loop over entities print the number of entities requested by each processor
            for ( int iProc = 0; iProc < tProcSize; ++iProc )
            {
                // Give each processor their desired amount of IDs
                tProcFirstID( iProc ) = tGlobalMax;

                // Increment the first available node ID
                tGlobalMax = tGlobalMax + aGatheredInfo( iProc );
            }
        }

        moris::scatter( tProcFirstID, tFirstId );

        return tFirstId( 0 );
    }

    // ----------------------------------------------------------------------------
    // multigrid accessor functions
    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_num_interpolations()
    {
        return mMeshIndices.max() + 1;
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_max_level( const moris_index aInterpolationIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_max_level( aInterpolationIndex );
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_num_basis( const moris_index aInterpolationIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_num_basis( aInterpolationIndex );
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_basis_level(
            const moris_index aInterpolationIndex,
            const moris_index aBasisIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_basis_level(
                aInterpolationIndex,
                aBasisIndex );
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_num_coarse_basis_of_basis(
            const moris_index aInterpolationIndex,
            const moris_index aBasisIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_num_coarse_basis_of_basis(
                aInterpolationIndex,
                aBasisIndex );
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_coarse_basis_index_of_basis(
            const moris_index aInterpolationIndex,
            const moris_index aBasisIndex,
            const moris_index aCoarseParentIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_coarse_basis_index_of_basis(
                aInterpolationIndex,
                aBasisIndex,
                aCoarseParentIndex );
    }

    // ----------------------------------------------------------------------------

    moris::Matrix< DDSMat >
    Enriched_Interpolation_Mesh::get_fine_basis_inds_of_basis(
            const moris_index aInterpolationIndex,
            const moris_index aBasisIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_fine_basis_inds_of_basis(
                aInterpolationIndex,
                aBasisIndex );
    }

    // ----------------------------------------------------------------------------

    moris::Matrix< DDRMat >
    Enriched_Interpolation_Mesh::get_fine_basis_weights_of_basis(
            const moris_index aInterpolationIndex,
            const moris_index aBasisIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_fine_basis_weights_of_basis(
                aInterpolationIndex,
                aBasisIndex );
    }

    // ----------------------------------------------------------------------------

}// namespace xtk

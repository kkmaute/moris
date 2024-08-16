/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Enriched_Interpolation_Mesh.cpp
 *
 */

#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_Map.hpp"

#include "cl_XTK_Multigrid.hpp"
#include "fn_TOL_Capacities.hpp"
#include "fn_stringify_matrix.hpp"
#include "cl_XTK_Field.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Field_Discrete.hpp"

// debug
#include "cl_Stopwatch.hpp"    //CHR/src

namespace moris::xtk
{
    // ----------------------------------------------------------------------------

    Enriched_Interpolation_Mesh::Enriched_Interpolation_Mesh(
            Model* aXTKModel )
            : mXTKModel( aXTKModel )
            , mNumVerts( 0 )
            , mNumVertsPerInterpCell( MORIS_UINT_MAX )
            , mCellInfo( nullptr )
            , mFieldLabelToIndex( 2 )
    {
    }

    // ----------------------------------------------------------------------------

    Enriched_Interpolation_Mesh::~Enriched_Interpolation_Mesh()
    {
        for ( uint i = 0; i < mInterpVertEnrichment.size(); i++ )
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

    mtk::MeshType
    Enriched_Interpolation_Mesh::get_mesh_type() const
    {
        return mtk::MeshType::XTK;
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_spatial_dim() const
    {
        return mXTKModel->get_spatial_dim();
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_num_entities(
            mtk::EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        switch ( aEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                return mEnrichedInterpVerts.size();
                break;
            }
            case mtk::EntityRank::ELEMENT:
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
        moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( aBSplineMeshIndex );

        return mEnrichCoeffLocToGlob( tLocalMeshIndex ).numel();
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Interpolation_Mesh::get_entity_connected_to_entity_loc_inds(
            moris_index       aEntityIndex,
            mtk::EntityRank   aInputEntityRank,
            mtk::EntityRank   aOutputEntityRank,
            const moris_index aIndex ) const
    {
        MORIS_ERROR( aInputEntityRank == mtk::EntityRank::ELEMENT && aOutputEntityRank == mtk::EntityRank::NODE,
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

    Vector< mtk::Vertex const * >
    Enriched_Interpolation_Mesh::get_all_vertices() const
    {
        uint tNumNodes = this->get_num_entities( mtk::EntityRank::NODE );

        Vector< mtk::Vertex const * > tVertices( tNumNodes );

        for ( uint i = 0; i < tNumNodes; i++ )
        {
            tVertices( i ) = mEnrichedInterpVerts( i );
        }

        return tVertices;
    }

    // ----------------------------------------------------------------------------

    moris_id
    Enriched_Interpolation_Mesh::get_glb_entity_id_from_entity_loc_index(
            moris_index       aEntityIndex,
            mtk::EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        uint tMapIndex = (uint)aEntityRank;

        MORIS_ASSERT( aEntityRank == mtk::EntityRank::NODE || aEntityRank == mtk::EntityRank::ELEMENT,
                "XTK ENRICHED MESH ERROR: Only support glb to loc conversion for vertices and cells" );

        MORIS_ASSERT( aEntityIndex < (moris_index)mLocalToGlobalMaps( tMapIndex ).numel(),
                "Entity index out of bounds" );

        return mLocalToGlobalMaps( tMapIndex )( aEntityIndex );
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_loc_entity_ind_from_entity_glb_id(
            moris_id          aEntityId,
            mtk::EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        uint tMapIndex = (uint)aEntityRank;

        auto tIter = mGlobalToLocalMaps( tMapIndex ).find( aEntityId );

        MORIS_ASSERT( aEntityRank == mtk::EntityRank::NODE || aEntityRank == mtk::EntityRank::ELEMENT,
                "XTK ENRICHED MESH ERROR: Only support glb to loc conversion for vertices and cells" );

        if ( tIter == mGlobalToLocalMaps( tMapIndex ).end() )
        {
            std::cout << "Not Found  Entity Id = " << aEntityId << " | par_rank = " << par_rank() << std::endl;
        }

        MORIS_ASSERT( tIter != mGlobalToLocalMaps( tMapIndex ).end(), "Id does not appear in map" );

        return tIter->second;
    }

    // ----------------------------------------------------------------------------

    const mtk::Interpolation_Mesh&
    Enriched_Interpolation_Mesh::get_background_mesh()
    {
        return mXTKModel->get_background_mesh();
    }

    // ----------------------------------------------------------------------------

    std::unordered_map< moris_id, moris_index >
    Enriched_Interpolation_Mesh::get_vertex_glb_id_to_loc_vertex_ind_map() const
    {
        uint tMapIndex = (uint)mtk::EntityRank::NODE;

        return mGlobalToLocalMaps( tMapIndex );
    }

    // ----------------------------------------------------------------------------

    Matrix< IdMat >
    Enriched_Interpolation_Mesh::get_entity_connected_to_entity_glob_ids(
            moris_id          aEntityId,
            mtk::EntityRank   aInputEntityRank,
            mtk::EntityRank   aOutputEntityRank,
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
        Matrix< IndexMat > tElementIndices( mElementIndicesInBlock( aSetIndex ).size(), 1 );

        std::copy( mElementIndicesInBlock( aSetIndex ).begin(), mElementIndicesInBlock( aSetIndex ).end(), tElementIndices.begin() );

        return tElementIndices;
    }

    // ----------------------------------------------------------------------------

    mtk::CellTopology
    Enriched_Interpolation_Mesh::get_blockset_topology( const std::string& aSetName )
    {
        uint tNumberOfDimensions = this->get_spatial_dim();

        moris::mtk::Interpolation_Order tOrder = mEnrichedInterpCells( 0 )->get_interpolation_order();

        mtk::CellTopology tCellTopology = mtk::CellTopology::UNDEFINED;

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
                        tCellTopology = mtk::CellTopology::QUAD4;
                        break;
                    }
                    case moris::mtk::Interpolation_Order::QUADRATIC:
                    {
                        tCellTopology = mtk::CellTopology::QUAD9;
                        break;
                    }
                    case moris::mtk::Interpolation_Order::CUBIC:
                    {
                        tCellTopology = mtk::CellTopology::QUAD16;
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
                        tCellTopology = mtk::CellTopology::HEX8;
                        break;
                    }
                    case moris::mtk::Interpolation_Order::QUADRATIC:
                    {
                        tCellTopology = mtk::CellTopology::HEX27;
                        break;
                    }
                    case moris::mtk::Interpolation_Order::CUBIC:
                    {
                        tCellTopology = mtk::CellTopology::HEX64;
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

    mtk::CellShape
    Enriched_Interpolation_Mesh::get_IG_blockset_shape( const std::string& aSetName )
    {
        // get the clusters in the set
        Vector< mtk::Cluster const * > tSetClusters = this->get_set_by_name( aSetName )->get_clusters_on_set();

        // init cell shape
        mtk::CellShape tCellShape = mtk::CellShape::EMPTY;

        // if the set isn't empty exist
        if ( tSetClusters.size() > 0 )
        {
            // get the cells in the first cluster
            Vector< moris::mtk::Cell const * > tClusterCells = tSetClusters( 0 )->get_primary_cells_in_cluster();

            // compute the cell shape of the first cell
            tCellShape = tClusterCells( 0 )->get_cell_info()->compute_cell_shape( tClusterCells( 0 ) );
        }

        // within debug, checking all cells to make sure that they are the same Cell Shape
        // if cells exist
        // looping through the clusters
        for ( uint iCluster = 0; iCluster < tSetClusters.size(); iCluster++ )
        {
            // get cell of cells in the cluster
            Vector< moris::mtk::Cell const * > tClusterCellsCheck = tSetClusters( iCluster )->get_primary_cells_in_cluster();

            // looping through the cells in the cluster
            for ( uint iCheckCell = 0; iCheckCell < tClusterCellsCheck.size(); iCheckCell++ )
            {
                MORIS_ASSERT( tClusterCellsCheck( iCheckCell )->get_cell_info()->compute_cell_shape( tClusterCellsCheck( iCheckCell ) ) == tCellShape,
                        "Enriched_Interpolation_Mesh::get_IG_blockset_shape() - cell shape is not consistent in the block" );
            }
        }

        return tCellShape;
    }

    // ----------------------------------------------------------------------------

    mtk::CellShape
    Enriched_Interpolation_Mesh::get_IP_blockset_shape( const std::string& aSetName )
    {
        // get the clusters in the set
        Vector< mtk::Cluster const * > tSetClusters = this->get_set_by_name( aSetName )->get_clusters_on_set();

        // init cell shape
        mtk::CellShape tCellShape = mtk::CellShape::EMPTY;

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
            mtk::EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        MORIS_ASSERT( aEntityRank == mtk::EntityRank::NODE || aEntityRank == mtk::EntityRank::ELEMENT,
                "Only Elements or Nodes have max id" );

        uint tNumEntities = this->get_num_entities( aEntityRank );

        moris_id tLocalMaxId = 0;

        for ( uint i = 0; i < tNumEntities; i++ )
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
        // reset the ADoF map before constructing it
        aAdofMap.clear();

        // get the index of the current mesh in the list of B-spline mesh indices
        moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( aBSplineIndex );

        // get the array mapping proc local enr. basis indices to their IDs
        Matrix< IdMat > const & tEnrBfIndToIdMap = mEnrichCoeffLocToGlob( tLocalMeshIndex );
        uint                    tNumEnrBFs       = tEnrBfIndToIdMap.numel();

        // go over basis functions and build the reverse ID to index map
        for ( uint iBF = 0; iBF < tNumEnrBFs; iBF++ )
        {
            // get the current BF's ID
            moris_id tEnrBfId = tEnrBfIndToIdMap( iBF );

            // make sure there is no other enr. BF that has been assigned this ID
            MORIS_ASSERT( !aAdofMap.key_exists( tEnrBfId ),
                    "Enriched_Interpolation_Mesh::get_adof_map() - "
                    "Duplicate enriched basis function ID in the basis map detected. "
                    "Trying to assign ID %i to index #%i, but it has already been assigned to #%i",
                    tEnrBfId,
                    iBF,
                    aAdofMap.find( tEnrBfId ) );

            // populate the ID to index map
            aAdofMap[ tEnrBfId ] = (moris_index)iBF;
        }
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat > const &
    Enriched_Interpolation_Mesh::get_enriched_coefficients_at_background_coefficient(
            moris_index const & aMeshIndex,
            moris_index         aBackgroundCoeffIndex ) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( aMeshIndex );

        MORIS_ASSERT( aBackgroundCoeffIndex < (moris_index)mCoeffToEnrichCoeffs( tLocalMeshIndex ).size(),
                "Background coefficient index out of bounds. Be sure this is not an enriched coefficient index passed in." );

        return mCoeffToEnrichCoeffs( tLocalMeshIndex )( aBackgroundCoeffIndex );
    }

    // ----------------------------------------------------------------------------

    Vector< Matrix< IndexMat > > const &
    Enriched_Interpolation_Mesh::get_enriched_coefficients_to_background_coefficients( moris_index const & aMeshIndex ) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( aMeshIndex );
        return mCoeffToEnrichCoeffs( tLocalMeshIndex );
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat > const &
    Enriched_Interpolation_Mesh::get_enriched_coefficient_local_to_global_map( moris_index const & aMeshIndex ) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( aMeshIndex );

        return mEnrichCoeffLocToGlob( tLocalMeshIndex );
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Interpolation_Mesh::get_background_coefficient_local_to_global_map() const
    {
        moris::mtk::Mesh& tBackgroundMeshData = mXTKModel->get_background_mesh();

        uint tNumBackgroundCoeffs = tBackgroundMeshData.get_num_entities( mBasisRank );

        Matrix< IndexMat > tCoefficientLocalToGlobal( 1, tNumBackgroundCoeffs );

        for ( uint i = 0; i < tNumBackgroundCoeffs; i++ )
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
        moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( aMeshIndex );
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
        // get index of B-spline mesh index in local list of associated B-spline meshes
        moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( aMeshIndex );

        // vertex index of the base interpolation vertex
        moris_index tBaseVertIndex = aBaseInterpVertex->get_index();

        // Number of enriched vertices related to the base vertex
        uint tNumVertsEnrOnBaseVert = mBaseInterpVertToVertEnrichmentIndex( tLocalMeshIndex )( tBaseVertIndex ).size();

        // not new until we make it to the end
        aNewVertex = false;

        // check whether the base vertex already has any kind of T-matrix on it
        bool tBaseVertexHasAnyInterpolation = aBaseInterpVertex->has_interpolation( aMeshIndex );

        // if the vertex has an interpolation then there is a potential that the t-matrices
        // are the same. If not, we always construct a new interpolation vertex here
        if ( tBaseVertexHasAnyInterpolation )
        {
            // iterate through the enriched vertices related to the base vertex and see if any are equal
            for ( uint iVertEnrLvl = 0; iVertEnrLvl < tNumVertsEnrOnBaseVert; iVertEnrLvl++ )
            {
                // std::cout<<" i = "<<i<<" | Base basis ind = "<<<<std::endl;

                // get the index of the enriched vertex
                moris_index tVertEnrIndex = mBaseInterpVertToVertEnrichmentIndex( tLocalMeshIndex )( tBaseVertIndex )( iVertEnrLvl );

                // get the current Vertex Enrichment to check against
                Vertex_Enrichment* tCurrentVertEnrichment = mInterpVertEnrichment( tLocalMeshIndex )( tVertEnrIndex );

                // check if the Vertex-Enrichment (VE) (i.e. enriched T-matrix) passed into this fnct. is the same as the VE already associated with this particular enr. vertex
                if ( tCurrentVertEnrichment != nullptr )
                {
                    // check whether the two T-matrices are equal
                    bool tVertEnrichmentsAreEqual = ( aVertexEnrichment == *tCurrentVertEnrichment );

                    // if the two T-matrices are the same, return the index of the already existing one and stop the function here
                    if ( tVertEnrichmentsAreEqual )
                    {
                        return tVertEnrIndex;
                    }
                }
            }
        }

        // if we make it through the loop without finding an enrichment vertex
        // make a new one
        aNewVertex = true;

        // index of the vertex enrichment (just check how many VEs there already are, next one up is the new index)
        moris_index tVertexEnrichmentIndex = mInterpVertEnrichment( tLocalMeshIndex ).size();

        // add new VE to member data
        mInterpVertEnrichment( tLocalMeshIndex ).push_back( new Vertex_Enrichment( aVertexEnrichment ) );

        // add a dummy value to the parent vertex index of a vertex interpolation
        mVertexEnrichmentParentVertexIndex( tLocalMeshIndex ).push_back( aBaseInterpVertex->get_index() );

        // update list of VEs living on current base vertex associated with the current DMI
        mBaseInterpVertToVertEnrichmentIndex( tLocalMeshIndex )( tBaseVertIndex ).push_back( tVertexEnrichmentIndex );

        // return the index of the new VE
        return tVertexEnrichmentIndex;
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::add_empty_vertex_enrichment(
            moris_index const & aMeshIndex,
            mtk::Vertex*        aBaseInterpVertex,
            bool&               aNewVertex )
    {
        // get index of B-spline mesh index in local list of associated B-spline meshes
        moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( aMeshIndex );

        // vertex index of the base interpolation vertex
        moris_index tBaseVertIndex = aBaseInterpVertex->get_index();

        // Number of enriched vertices related to the base vertex
        uint tNumVertsEnrOnBaseVert = mBaseInterpVertToVertEnrichmentIndex( tLocalMeshIndex )( tBaseVertIndex ).size();

        // not new until we make it to the end
        aNewVertex = false;

        // if the vertex has an interpolation then there is a potential that the t-matrices
        // are the same. If not we always construct a new interpolation vertex here
        if ( aBaseInterpVertex->has_interpolation( aMeshIndex ) )
        {
            // iterate through the enriched vertices related to the base vertex and see if any are equal
            for ( uint iVertEnrLvl = 0; iVertEnrLvl < tNumVertsEnrOnBaseVert; iVertEnrLvl++ )
            {
                // std::cout<<" i = "<<i<<" | Base basis ind = "<<<<std::endl;

                // get the index of the enriched vertex
                moris_index tVertEnrIndex = mBaseInterpVertToVertEnrichmentIndex( tLocalMeshIndex )( tBaseVertIndex )( iVertEnrLvl );

                // check if the Vertex-Enrichment (VE) passed into this fnct. is the same as the VE already associated with this particular enr. vertex
                if ( mInterpVertEnrichment( tLocalMeshIndex )( tVertEnrIndex ) == nullptr )
                {
                    return tVertEnrIndex;
                }
            }
        }

        // if we make it through the loop without finding an enrichment vertex
        // make a new one
        aNewVertex = true;

        // index of the vertex enrichment (just check how many VEs there already are, next one up is the new index)
        moris_index tVertexEnrichmentIndex = mInterpVertEnrichment( tLocalMeshIndex ).size();

        // add new VE to member data
        mInterpVertEnrichment( tLocalMeshIndex ).push_back( nullptr );

        // add a dummy value to the parent vertex index of a vertex interpolation
        mVertexEnrichmentParentVertexIndex( tLocalMeshIndex ).push_back( tVertexEnrichmentIndex );

        // update list of VEs living on current base vertex associated with the current DMI
        mBaseInterpVertToVertEnrichmentIndex( tLocalMeshIndex )( tBaseVertIndex ).push_back( tVertexEnrichmentIndex );

        // return the index of the new VE
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
        Vector< Vector< Interpolation_Vertex_Unzipped* > > tBaseVertexToEnrichedVertex;
        this->collect_base_vertex_to_enriched_vertex_connectivity( tBaseVertexToEnrichedVertex );

        // new vertex index
        Vector< moris_index > tNewIndex( mEnrichedInterpVerts.size(), MORIS_INDEX_MAX );
        Vector< moris_index > tNodesToDelete;
        tNodesToDelete.reserve( mEnrichedInterpVerts.size() );

        for ( uint iBV = 0; iBV < tBaseVertexToEnrichedVertex.size(); iBV++ )
        {
            // keep track of the first vertex that a given vertex wants to merge with
            Vector< moris_index > tMergeWithVertex( tBaseVertexToEnrichedVertex( iBV ).size(), MORIS_INDEX_MAX );

            for ( uint iEV = 0; iEV < tBaseVertexToEnrichedVertex( iBV ).size(); iEV++ )
            {
                for ( uint iOtherEV = 0; iOtherEV < iEV + 1; iOtherEV++ )
                {
                    if ( iEV != iOtherEV )
                    {
                        bool tMerge = true;
                        for ( uint iMeshIndex = 0; iMeshIndex < mInterpVertEnrichment.size(); iMeshIndex++ )
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
            for ( uint iMerge = 0; iMerge < tMergeWithVertex.size(); iMerge++ )
            {
                moris_index tOldIndex = tBaseVertexToEnrichedVertex( iBV )( iMerge )->get_index();
                if ( tMergeWithVertex( iMerge ) == MORIS_INDEX_MAX )
                {
                    MORIS_ASSERT( tNewIndex( tOldIndex ) == MORIS_INDEX_MAX, "OVERWRITING A NEW NODE INDEX" );
                    tNewIndex( tOldIndex ) = tOldIndex;
                }
                else
                {
                    MORIS_ASSERT( tNewIndex( tOldIndex ) == MORIS_INDEX_MAX, "OVERWRITING A NEW NODE INDEX" );
                    tNewIndex( tOldIndex ) = tBaseVertexToEnrichedVertex( iBV )( tMergeWithVertex( iMerge ) )->get_index();
                    tNodesToDelete.push_back( tOldIndex );
                }
            }
        }

        // update cell to vertex connectivity (change the pointers)
        for ( uint iCell = 0; iCell < mEnrichedInterpCells.size(); iCell++ )
        {
            for ( uint iV = 0; iV < mEnrichedInterpCells( iCell )->mVertices.size(); iV++ )
            {
                moris_index tNewVertexIndex                    = tNewIndex( mEnrichedInterpCells( iCell )->mVertices( iV )->get_index() );
                mEnrichedInterpCells( iCell )->mVertices( iV ) = mEnrichedInterpVerts( tNewVertexIndex );
            }
        }

        // delete the removed nodes
        for ( uint iDelete = 0; iDelete < tNodesToDelete.size(); iDelete++ )
        {
            delete mEnrichedInterpVerts( tNodesToDelete( iDelete ) );
            mEnrichedInterpVerts( tNodesToDelete( iDelete ) ) = nullptr;
        }

        mEnrichedInterpVerts.data().erase( std::remove( std::begin( mEnrichedInterpVerts.data() ), std::end( mEnrichedInterpVerts.data() ), nullptr ), std::end( mEnrichedInterpVerts.data() ) );

        MORIS_ERROR( mEnrichedInterpVerts.size() == tStartNumVerts - tNodesToDelete.size(), "Something went wrong merging nodes." );

        // reindex the vertices
        for ( uint iV = 0; iV < mEnrichedInterpVerts.size(); iV++ )
        {
            mEnrichedInterpVerts( iV )->mVertexIndex = (moris_index)iV;
        }

        for ( uint iCell = 0; iCell < mEnrichedInterpCells.size(); iCell++ )
        {
            for ( uint iV = 0; iV < mEnrichedInterpCells( iCell )->mVertices.size(); iV++ )
            {
                MORIS_ERROR( mEnrichedInterpCells( iCell )->mVertices( iV ) != nullptr, "Nullptr in cell vertices" );
            }
        }

        mNumVerts = mEnrichedInterpVerts.size();

        MORIS_LOG_SPEC( "Num Enriched Interpolation Vertices Post Merge", mEnrichedInterpVerts.size() );
    }

    // ----------------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::collect_base_vertex_to_enriched_vertex_connectivity( Vector< Vector< Interpolation_Vertex_Unzipped* > >& aBaseVertexToEnrichedVertex )
    {
        // allocate space in the cell of cell
        aBaseVertexToEnrichedVertex.resize( mXTKModel->get_background_mesh().get_num_nodes() );

        for ( uint iN = 0; iN < mEnrichedInterpVerts.size(); iN++ )
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
        moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( aMeshIndex );

        MORIS_ASSERT( aVertexEnrichmentIndex < (moris_index)mInterpVertEnrichment( tLocalMeshIndex ).size(),
                "Provided vertex enrichment index out of bounds" );

        xtk::Vertex_Enrichment* tVertexEnrichment = mInterpVertEnrichment( tLocalMeshIndex )( aVertexEnrichmentIndex );

        return tVertexEnrichment;
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_vertex_related_to_vertex_enrichment(
            moris_index const & aMeshIndex,
            moris_index         aVertexEnrichmentIndex ) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( aMeshIndex );

        MORIS_ASSERT( aVertexEnrichmentIndex < (moris_index)mVertexEnrichmentParentVertexIndex( tLocalMeshIndex ).size(),
                "Provided vertex enrichment index out of bounds" );

        return mVertexEnrichmentParentVertexIndex( tLocalMeshIndex )( aVertexEnrichmentIndex );
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_local_mesh_index_xtk( moris_index const & aMeshIndex ) const
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
        moris_index tMeshIndex = this->get_local_mesh_index_xtk( aMeshIndex );

        if ( mGlobalToLocalBasisMaps( tMeshIndex ).find( aBasisId ) == mGlobalToLocalBasisMaps( tMeshIndex ).end() )
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

        moris_index tLocMesh  = this->get_local_mesh_index_xtk( aMeshIndex );
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
        mGlobalToLocalBasisMaps( tLocMesh )[ aBasisIdToAdd ] = tNewIndex;

        return tNewIndex;
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::add_basis_functions(
            moris_index const &           aMeshIndex,
            Vector< moris_id > const &    aBfIdsToAdd,
            Vector< moris_id > const &    aBfOwners,
            Vector< moris_index > const & aBfBulkPhases )
    {
        // perform some checks on the inputs
        uint tNumBfs = aBfIdsToAdd.size();
        MORIS_ASSERT( aBfOwners.size() == tNumBfs && aBfBulkPhases.size() == tNumBfs,
                "xtk::Enriched_Interpolation_Mesh::add_basis_functions() - "
                "Input arrays not of equal size" );

        // get the discretization mesh index
        moris_index tLocMesh       = this->get_local_mesh_index_xtk( aMeshIndex );
        moris_index tFirstNewIndex = mEnrichCoeffLocToGlob( tLocMesh ).numel();

        // add a size of 1
        mEnrichCoeffLocToGlob( tLocMesh ).resize( 1, tFirstNewIndex + tNumBfs );
        mEnrichCoeffOwnership( tLocMesh ).resize( 1, tFirstNewIndex + tNumBfs );
        mEnrichCoeffBulkPhase( tLocMesh ).resize( 1, tFirstNewIndex + tNumBfs );

        for ( uint iBf = 0; iBf < tNumBfs; iBf++ )
        {
            // get the current BF's index
            moris_index tBfIndex = tFirstNewIndex + iBf;

            // get the info
            moris_id    tBfId        = aBfIdsToAdd( iBf );
            moris_id    tBfOwner     = aBfOwners( iBf );
            moris_index tBfBulkPhase = aBfBulkPhases( iBf );

            // check that the BF doesn't already exist
            MORIS_ASSERT(
                    !this->basis_exists_on_partition( aMeshIndex, tBfId ),
                    "Enriched_Interpolation_Mesh::add_basis_functions() - "
                    "The basis function (ID: %i) to be added already exists in this mesh.",
                    tBfId );

            // add the local to glb map
            mEnrichCoeffLocToGlob( tLocMesh )( tBfIndex ) = tBfId;
            mEnrichCoeffOwnership( tLocMesh )( tBfIndex ) = tBfOwner;
            mEnrichCoeffBulkPhase( tLocMesh )( tBfIndex ) = tBfBulkPhase;

            // add to glb to local map
            mGlobalToLocalBasisMaps( tLocMesh )[ tBfId ] = tBfIndex;
        }

    }    // end function: Enriched_Interpolation_Mesh::add_basis_functions()

    // ----------------------------------------------------------------------------

    Vector< Interpolation_Cell_Unzipped const * >
    Enriched_Interpolation_Mesh::get_enriched_cells_from_base_cells(
            Vector< moris::mtk::Cell const * > const & aBaseCells ) const
    {
        uint tNumBaseCells = aBaseCells.size();

        Vector< Interpolation_Cell_Unzipped const * > tEnrichedCells;

        for ( uint i = 0; i < tNumBaseCells; i++ )
        {
            tEnrichedCells.append( this->get_enriched_cells_from_base_cell( aBaseCells( i ) ) );
        }

        return tEnrichedCells;
    }

    // ----------------------------------------------------------------------------

    Vector< Interpolation_Cell_Unzipped* > const &
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

    Matrix< IndexMat >
    Enriched_Interpolation_Mesh::get_enriched_mesh_indices() const
    {
        return mMeshIndices;
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
        moris_index tLocMesh = this->get_local_mesh_index_xtk( aMeshIndex );

        return mEnrichCoeffOwnership( tLocMesh )( aBasisIndex );
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_basis_bulk_phase( moris_index const & aBasisIndex,
            moris_index const &                                            aMeshIndex ) const
    {
        moris_index tLocMesh = this->get_local_mesh_index_xtk( aMeshIndex );

        return mEnrichCoeffBulkPhase( tLocMesh )( aBasisIndex );
    }

    // ----------------------------------------------------------------------------

    Vector< Interpolation_Cell_Unzipped* >&
    Enriched_Interpolation_Mesh::get_enriched_interpolation_cells()
    {
        return mEnrichedInterpCells;
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::get_owned_and_not_owned_enriched_interpolation_cells(
            Vector< Interpolation_Cell_Unzipped* >&           aOwnedInterpCells,
            Vector< Vector< Interpolation_Cell_Unzipped* > >& aNotOwnedInterpCells,
            Vector< uint >&                                   aProcRanks )
    {
        // get all interp cells
        Vector< Interpolation_Cell_Unzipped* >& tEnrInterpCells = this->get_enriched_interpolation_cells();

        // reserve space
        aOwnedInterpCells.resize( 0 );
        aOwnedInterpCells.reserve( tEnrInterpCells.size() );
        aNotOwnedInterpCells.resize( 0 );

        // proc rank
        moris_index tParRank = par_rank();

        // counter
        uint           tOwnerCount = 0;
        Vector< uint > tCounts( 0 );

        // map
        std::unordered_map< moris_id, moris_id > tProcRankToDataIndex;

        // access the communication table
        Matrix< IdMat > tCommTable = this->get_communication_table();

        // resize proc ranks and setup map to comm table
        aProcRanks.resize( tCommTable.numel() );
        for ( uint i = 0; i < tCommTable.numel(); i++ )
        {
            tProcRankToDataIndex[ tCommTable( i ) ] = i;
            aProcRanks( i )                         = ( tCommTable( i ) );
            aNotOwnedInterpCells.push_back( Vector< Interpolation_Cell_Unzipped* >( 0 ) );
        }

        for ( uint i = 0; i < tEnrInterpCells.size(); i++ )
        {
            moris_index tOwnerProc = tEnrInterpCells( i )->get_owner();

            if ( tParRank == tOwnerProc )
            {
                aOwnedInterpCells.push_back( tEnrInterpCells( i ) );
                tOwnerCount++;
            }
            else
            {
                moris_index tProcDataIndex = tProcRankToDataIndex[ tOwnerProc ];
                aNotOwnedInterpCells( tProcDataIndex ).push_back( tEnrInterpCells( i ) );
            }
        }
    }

    // ----------------------------------------------------------------------------

    Interpolation_Vertex_Unzipped&
    Enriched_Interpolation_Mesh::get_xtk_interp_vertex( uint aVertexIndex )
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
        moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( aMeshIndex );

        auto tIter = mGlobalToLocalBasisMaps( tLocalMeshIndex ).find( aBasisId );

        MORIS_ASSERT( tIter != mGlobalToLocalBasisMaps( tLocalMeshIndex ).end(),
                "Basis id not in map" );

        return tIter->second;
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_enr_basis_id_from_enr_basis_index(
            moris_index const & aMeshIndex,
            moris_index const & aBasisIndex ) const
    {
        moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( aMeshIndex );

        return mEnrichCoeffLocToGlob( tLocalMeshIndex )( aBasisIndex );
    }

    // ----------------------------------------------------------------------------

    Vector< Interpolation_Cell_Unzipped const * >
    Enriched_Interpolation_Mesh::get_enriched_cells_from_base_cell( moris::mtk::Cell const * aBaseCells ) const
    {
        moris_index tBaseIndex = aBaseCells->get_index();

        MORIS_ASSERT( tBaseIndex < (moris_index)mBaseCellToEnrichedCell.size(),
                "Base Cell index is out of bounds. This index is related to the non-enriched interpolation mesh. Make sure enriched cell is not passed into this function" );

        uint tNumEnrichedCells = mBaseCellToEnrichedCell( tBaseIndex ).size();

        Vector< Interpolation_Cell_Unzipped const * > tEnrichedCellPtrs( tNumEnrichedCells );

        for ( uint i = 0; i < tNumEnrichedCells; i++ )
        {
            tEnrichedCellPtrs( i ) = ( mBaseCellToEnrichedCell( tBaseIndex )( i ) );
        }

        return tEnrichedCellPtrs;
    }

    // ----------------------------------------------------------------------------

    Interpolation_Vertex_Unzipped const &
    Enriched_Interpolation_Mesh::get_xtk_interp_vertex( uint aVertexIndex ) const
    {
        return *mEnrichedInterpVerts( aVertexIndex );
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::set_enriched_basis_to_bulkphase_map(
            const moris_index aMeshIndex,
            Matrix< IdMat >   aBulkPhaseInEnrichedBasis )
    {
        // set the enr. BF index to bulk-phase index map
        mEnrichCoeffBulkPhase( aMeshIndex ) = aBulkPhaseInEnrichedBasis;
    }

    // ----------------------------------------------------------------------------

    Matrix< IdMat >
    Enriched_Interpolation_Mesh::convert_indices_to_ids(
            Matrix< IndexMat > const & aIndices,
            mtk::EntityRank            aEntityRank ) const
    {
        uint tNRow = aIndices.n_rows();
        uint tNCol = aIndices.n_cols();

        Matrix< IdMat > tIds( tNRow, tNCol );

        for ( uint i = 0; i < tNRow; i++ )
        {
            for ( uint j = 0; j < tNCol; j++ )
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
            mtk::EntityRank         aEntityRank ) const
    {
        uint tNRow = aIds.n_rows();
        uint tNCol = aIds.n_cols();

        Matrix< IdMat > tIndices( tNRow, tNCol );

        for ( uint i = 0; i < tNRow; i++ )
        {
            for ( uint j = 0; j < tNCol; j++ )
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
        moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( aMeshIndex );

        aEnrichedIds.resize( aEnrichedIndices.n_rows(), aEnrichedIndices.n_cols() );

        for ( uint i = 0; i < aEnrichedIndices.n_rows(); i++ )
        {
            for ( uint j = 0; j < aEnrichedIndices.n_cols(); j++ )
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
        tMM.mMemoryMapData[ "mXTKModel ptr" ]  = sizeof( mXTKModel );
        tMM.mMemoryMapData[ "mBasisRank ptr" ] = sizeof( mBasisRank );
        tMM.mMemoryMapData[ "mMeshIndices" ]   = mMeshIndices.capacity();
        // FIXME: add mMeshIndexToLocMeshIndex
        tMM.mMemoryMapData[ "mNumVerts" ]                            = sizeof( mNumVerts );
        tMM.mMemoryMapData[ "mNumVertsPerInterpCell" ]               = sizeof( mNumVertsPerInterpCell );
        tMM.mMemoryMapData[ "mBaseInterpVertToVertEnrichmentIndex" ] = moris::internal_capacity_nested( mBaseInterpVertToVertEnrichmentIndex );
        tMM.mMemoryMapData[ "mInterpVertEnrichment" ]                = moris::internal_capacity_nested_ptr( mInterpVertEnrichment );
        tMM.mMemoryMapData[ "mVertexEnrichmentParentVertexIndex" ]   = moris::internal_capacity( mVertexEnrichmentParentVertexIndex );
        tMM.mMemoryMapData[ "mVertexBulkPhase" ]                     = mVertexBulkPhase.capacity();
        tMM.mMemoryMapData[ "mCoeffToEnrichCoeffs" ]                 = moris::internal_capacity_nested( mCoeffToEnrichCoeffs );
        tMM.mMemoryMapData[ "mEnrichCoeffLocToGlob" ]                = moris::internal_capacity( mEnrichCoeffLocToGlob );
        // fixme: add me mGlobalToLocalBasisMaps
        tMM.mMemoryMapData[ "mEnrichCoeffOwnership" ] = moris::internal_capacity( mEnrichCoeffOwnership );
        tMM.mMemoryMapData[ "mLocalToGlobalMaps" ]    = moris::internal_capacity( mLocalToGlobalMaps );
        // fixme: add mGlobalToLocalMaps
        tMM.mMemoryMapData[ "mBaseCellToEnrichedCell" ] = moris::internal_capacity( mBaseCellToEnrichedCell );
        tMM.mMemoryMapData[ "mCellInfo" ]               = sizeof( mCellInfo );
        tMM.mMemoryMapData[ "mNotOwnedBasis" ]          = mNotOwnedBasis.capacity();
        tMM.mMemoryMapData[ "mOwnedBasis" ]             = mOwnedBasis.capacity();
        return tMM;
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::convert_enriched_basis_indices_to_ids(
            moris_index const &                  aMeshIndex,
            Vector< Matrix< IndexMat > > const & aEnrichedIndices,
            Vector< Matrix< IdMat > >&           aEnrichedIds ) const
    {
        aEnrichedIds.resize( aEnrichedIndices.size() );

        for ( uint i = 0; i < aEnrichedIndices.size(); i++ )
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
        Vector< Interpolation_Cell_Unzipped* > const & tEnrIPCells = this->get_enriched_interpolation_cells();

        std::ostringstream tStringStream;
        // max num verts to cells
        uint tMaxVertsToCell = 0;
        for ( uint i = 0; i < this->get_num_entities( mtk::EntityRank::ELEMENT, 0 ); i++ )
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
        for ( uint iVH = 0; iVH < tMaxVertsToCell; iVH++ )
        {
            tStringStream << "Vert_" + std::to_string( iVH );

            if ( iVH != tMaxVertsToCell - 1 )
            {
                tStringStream << ",";
            }
        }
        tStringStream << "\n";

        for ( uint i = 0; i < this->get_num_entities( mtk::EntityRank::ELEMENT, 0 ); i++ )
        {
            Interpolation_Cell_Unzipped const * tCell     = tEnrIPCells( (moris_index)i );
            Vector< moris::mtk::Vertex* >       tVertices = tCell->get_vertex_pointers();

            tStringStream << tCell->get_id() << ",";
            tStringStream << tCell->get_index() << ",";
            tStringStream << std::to_string( tCell->get_owner() ) << ",";
            tStringStream << std::to_string( par_rank() ) << ",";
            tStringStream << tCell->get_base_cell()->get_id() << ",";
            tStringStream << tCell->get_bulkphase_index() << ",";
            tStringStream << mXTKModel->get_cut_integration_mesh()->get_subphase_id( (moris_index)tCell->get_subphase_index() ) << ",";
            tStringStream << std::scientific << tCell->compute_cell_measure() << ",";

            for ( uint j = 0; j < tMaxVertsToCell; j++ )
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

        for ( uint iVH = 0; iVH < this->get_spatial_dim(); iVH++ )
        {
            tStringStream << "Coords_" + std::to_string( iVH );

            if ( iVH != this->get_spatial_dim() - 1 )
            {
                tStringStream << ",";
            }
        }

        tStringStream << std::endl;

        for ( uint i = 0; i < this->get_num_entities( mtk::EntityRank::NODE, 0 ); i++ )
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

            Matrix< DDRMat > tCoords = tVertex.get_coords();

            for ( uint iSp = 0; iSp < this->get_spatial_dim(); iSp++ )
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
        for ( uint iV = 0; iV < this->get_num_nodes(); iV++ )
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

        for ( uint iV = 0; iV < this->get_num_nodes(); iV++ )
        {
            tStringStream << this->get_mtk_vertex( (moris_index)iV ).get_id() << ",";
            mtk::Vertex_Interpolation* tVertexInterp = this->get_mtk_vertex( (moris_index)iV ).get_interpolation( aMeshIndex );
            Matrix< IdMat >            tBasisIds     = tVertexInterp->get_ids();
            const Matrix< DDRMat >*    tBasisWeights = tVertexInterp->get_weights();

            tStringStream << tBasisIds.numel() << ",";
            MORIS_ASSERT( tBasisIds.numel() == tBasisWeights->numel(), "Size mismatch" );

            for ( uint iB = 0; iB < tBasisIds.numel(); iB++ )
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
        uint tNumNodes = this->get_num_entities( mtk::EntityRank::NODE );

        uint tMapIndex = (uint)( mtk::EntityRank::NODE );

        MORIS_ASSERT( mLocalToGlobalMaps( tMapIndex ).numel() == tNumNodes,
                "Enriched_Interpolation_Mesh::print_vertex_maps: number of nodes and size of map do not match." );

        std::cout << "\nVertex Map:" << std::endl;

        for ( uint i = 0; i < tNumNodes; i++ )
        {
            std::cout << "    Vertex Index: " << std::setw( 9 ) << i << " | Vertex Id: " << std::setw( 9 ) << mLocalToGlobalMaps( tMapIndex )( i ) << std::endl;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::print_enriched_cell_maps() const
    {
        uint tNumCells = this->get_num_entities( mtk::EntityRank::ELEMENT );

        uint tMapIndex = (uint)( mtk::EntityRank::ELEMENT );

        MORIS_ASSERT( mLocalToGlobalMaps( tMapIndex ).numel() == tNumCells,
                "Enriched_Interpolation_Mesh::print_enriched_cell_maps: number of elements and size of map do not match." );

        std::cout << "\nCell Map:" << std::endl;

        for ( uint i = 0; i < tNumCells; i++ )
        {
            std::cout << "    Cell Index: " << std::setw( 9 ) << i << " | Cell Id: " << std::setw( 9 ) << mLocalToGlobalMaps( tMapIndex )( i ) << std::endl;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::print_basis_to_enriched_basis() const
    {
        for ( uint iM = 0; iM < mMeshIndices.numel(); iM++ )
        {
            uint tNumBasis = mCoeffToEnrichCoeffs( iM ).size();

            std::cout << "\nBackground Basis to Enriched Basis Indices For Mesh: " << mMeshIndices( iM ) << std::endl;

            for ( uint iB = 0; iB < tNumBasis; iB++ )
            {
                std::cout << "    Basis Index: " << std::setw( 9 ) << iB << " | Enriched Indices";

                for ( uint iEB = 0; iEB < mCoeffToEnrichCoeffs( iM )( iB ).numel(); iEB++ )
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
        uint tNumVerts = this->get_num_entities( mtk::EntityRank::NODE );

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
            mtk::Cluster const &       aCluster,
            moris_index const &        aMeshIndex,
            const mtk::Leader_Follower aIsLeader )
    {
        bool tDiagnosticFlag = true;

        // interpolation cell
        moris::mtk::Cell const & tIpCell = aCluster.get_interpolation_cell( aIsLeader );

        // get the xtk interpolation cell
        Interpolation_Cell_Unzipped* tEnrichedIpCell = mEnrichedInterpCells( tIpCell.get_index() );

        MORIS_ASSERT( tIpCell.get_id() == tEnrichedIpCell->get_id(), "Id mismatch" );

        // get the bulkphase index
        moris_index tBulkPhase = tEnrichedIpCell->get_bulkphase_index();

        // vertices attached to the interpolation cells
        Vector< mtk::Vertex* > tVertexPointers = tIpCell.get_vertex_pointers();

        // iterate through vertex pointers
        for ( uint i = 0; i < tVertexPointers.size(); i++ )
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
            for ( uint iB = 0; iB < tBasisIndices.numel(); iB++ )
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
    Enriched_Interpolation_Mesh::finalize_setup_new()
    {
        this->setup_local_to_global_maps();

        this->setup_basis_ownership();

        // note: this information is feed to the enr. IP mesh from the outside, as it is already constructed in the enrichment data
        // this->setup_basis_to_bulk_phase(); // not needed anymore
        mEnrichCoeffBulkPhase.resize( mMeshIndices.max() + 1 );

        // TODO: this check needs to be re-written without the assumption one IP cell == one particular bulk phase
        // MORIS_ASSERT( this->verify_basis_support(), "Issue detected in basis support." );
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_basis_to_bulk_phase()
    {
        // size member data
        mEnrichCoeffBulkPhase.resize( mMeshIndices.max() + 1 );

        // verify all the subphases in a enriched basis support are the same bulk phase
        for ( moris::size_t iBspMesh = 0; iBspMesh < mMeshIndices.numel(); iBspMesh++ )
        {
            // Mesh index
            moris_index tMeshIndex = mMeshIndices( iBspMesh );

            // Number of enriched functions
            moris::size_t tNumEnrBasis = mEnrichCoeffLocToGlob( tMeshIndex ).numel();

            // allocate interpolation cells in basis support // input: enr. BF index || output: list of UIPCs in the BF's support
            Vector< Vector< Interpolation_Cell_Unzipped* > > tCellsInEnrSupports( tNumEnrBasis );

            Vector< Interpolation_Cell_Unzipped* > const & tEnrIpCells = this->get_enriched_interpolation_cells();

            // number of cells
            moris_index tNumCells = this->get_num_entities( mtk::EntityRank::ELEMENT );

            MORIS_ASSERT( tNumCells == (moris_index)tEnrIpCells.size(), "Inconsistent num cells information." );

            // collect all UIPCs that a given enriched basis function interpolates into
            for ( moris_index iUIPC = 0; iUIPC < tNumCells; iUIPC++ )
            {
                Vector< xtk::Interpolation_Vertex_Unzipped* > const & tVertices =
                        tEnrIpCells( iUIPC )->get_xtk_interpolation_vertices();

                for ( moris::size_t iVert = 0; iVert < tVertices.size(); iVert++ )
                {
                    if ( tVertices( iVert )->has_interpolation( tMeshIndex ) )
                    {
                        Vertex_Enrichment* tVertInterp = tVertices( iVert )->get_xtk_interpolation( tMeshIndex );

                        Matrix< IndexMat > tBasisIndices = tVertInterp->get_indices();

                        // iterate through vertices and add the UIPCs cell to the list of cells supported by the basis
                        for ( uint iBF = 0; iBF < tBasisIndices.numel(); iBF++ )
                        {
                            tCellsInEnrSupports( tBasisIndices( iBF ) ).push_back( tEnrIpCells( iUIPC ) );
                        }
                    }
                }
            }

            // initialize list associating enr. BFs and the bulk-phase they interpolate into
            mEnrichCoeffBulkPhase( tMeshIndex ).resize( 1, tNumEnrBasis );
            mEnrichCoeffBulkPhase( tMeshIndex ).fill( MORIS_INDEX_MAX );

            // set error flag
            bool tIsConsistent = true;

            // iterate through enriched basis functions
            for ( moris::size_t iBF = 0; iBF < tNumEnrBasis; iBF++ )
            {
                // iterate through enriched interpolation cells in the support
                moris_index tNumSubphaseInSupport = tCellsInEnrSupports( iBF ).size();

                for ( moris::moris_index iSP = 0; iSP < tNumSubphaseInSupport; iSP++ )
                {
                    Interpolation_Cell_Unzipped* tIpCell = tCellsInEnrSupports( iBF )( iSP );

                    moris_index tBulkPhase = tIpCell->get_bulkphase_index();

                    MORIS_ASSERT( tBulkPhase != MORIS_INDEX_MAX, "Bulk phase index not set." );

                    // check that all subphases in basis support have same bulk phase as first subphase
                    if ( iSP == 0 )
                    {
                        mEnrichCoeffBulkPhase( tMeshIndex )( iBF ) = tBulkPhase;
                    }

                    if ( tBulkPhase != mEnrichCoeffBulkPhase( tMeshIndex )( iBF ) )
                    {
                        tIsConsistent = false;

                        std::cout << "enriched basis = " << iBF                                                  //
                                  << "  subphase = " << iSP                                                      //
                                  << "  expected bulk phase = " << mEnrichCoeffBulkPhase( tMeshIndex )( iBF )    //
                                  << "  current bulk phase  = " << tBulkPhase                                    //
                                  << std::endl;

                        // print( tCellsInEnrSupports( iBF )( 0 )->get_vertex_coords(), "Vertex coordinates of subphase 0" );
                        // print( tIpCell->get_vertex_coords(), "Vertex coordinates of subphase i" );
                    }
                }
            }
            MORIS_ERROR( tIsConsistent,
                    "Subphase in enriched basis function support not consistent bulk phase" );
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

        Vector< Interpolation_Cell_Unzipped* > const & tEnrIpCells = this->get_enriched_interpolation_cells();

        // number of cells
        moris_index tNumEnrIpCells = this->get_num_entities( mtk::EntityRank::ELEMENT );

        // create the enriched interpolation basis to interpolation cell interpolation
        for ( moris::moris_index iEnrIpCell = 0; iEnrIpCell < tNumEnrIpCells; iEnrIpCell++ )
        {
            Vector< xtk::Interpolation_Vertex_Unzipped* > const & tVertices = tEnrIpCells( iEnrIpCell )->get_xtk_interpolation_vertices();

            moris_index const tSubphaseIndex = tEnrIpCells( iEnrIpCell )->get_subphase_index();
            moris_index const tSubphaseId    = mXTKModel->get_cut_integration_mesh()->get_subphase_id( tSubphaseIndex );

            for ( moris::size_t iVertex = 0; iVertex < tVertices.size(); iVertex++ )
            {
                moris_index tVertexIndex = tVertices( iVertex )->get_index();

                if ( mVertexBulkPhase( tVertexIndex ) == MORIS_INDEX_MAX )
                {
                    mVertexBulkPhase( tVertexIndex ) = tEnrIpCells( iEnrIpCell )->get_bulkphase_index();
                }

                if ( tSubphaseId > mVertexMaxSubphase( tVertexIndex ) )
                {
                    mVertexMaxSubphase( tVertexIndex ) = tSubphaseId;
                }

                else
                {
                    if ( mVertexBulkPhase( tVertexIndex ) != tEnrIpCells( iEnrIpCell )->get_bulkphase_index() )
                    {

                        std::cout << " Vert Id = " << tVertices( iVertex )->get_id()
                                  << " | Vert Owner = " << tVertices( iVertex )->get_owner()
                                  << " | mVertexBulkPhase(tVertexIndex) = " << mVertexBulkPhase( tVertexIndex )
                                  << " | tEnrIpCells(i)->get_bulkphase_index() =" << tEnrIpCells( iEnrIpCell )->get_bulkphase_index()
                                  << std::endl;
                    }

                    // check that vertices on current Enr. IP cell have same Bulk phase as Enr. IP cell
                    MORIS_ASSERT( mVertexBulkPhase( tVertexIndex ) == tEnrIpCells( iEnrIpCell )->get_bulkphase_index(),
                            "Enriched_Interpolation_Mesh::setup_vertex_to_bulk_phase() - Inconsistent vertex bulk phase" );
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
            Vector< Vector< Interpolation_Cell_Unzipped* > > tCellsInEnrSupports( tNumEnrBasis );

            Vector< Interpolation_Cell_Unzipped* > const & tEnrIpCells = this->get_enriched_interpolation_cells();

            // number of cells
            moris_index tNumCells = this->get_num_entities( mtk::EntityRank::ELEMENT );

            MORIS_ASSERT( tNumCells == (moris_index)tEnrIpCells.size(), "Inconsistent num cells information." );

            // create the enriched interpolation basis to interpolation cell interpolation
            for ( moris::sint iUIPC = 0; iUIPC < tNumCells; iUIPC++ )
            {
                Vector< xtk::Interpolation_Vertex_Unzipped* > const & tVertices = tEnrIpCells( iUIPC )->get_xtk_interpolation_vertices();

                for ( moris::size_t iVert = 0; iVert < tVertices.size(); iVert++ )
                {
                    if ( tVertices( iVert )->has_interpolation( tMeshIndex ) )
                    {
                        Vertex_Enrichment* tVertInterp = tVertices( iVert )->get_xtk_interpolation( tMeshIndex );

                        Matrix< IndexMat > tBasisIndices = tVertInterp->get_indices();

                        // iterate through vertices and add the ip cell to the support
                        for ( uint iBF = 0; iBF < tBasisIndices.numel(); iBF++ )
                        {
                            tCellsInEnrSupports( tBasisIndices( iBF ) ).push_back( tEnrIpCells( iUIPC ) );
                        }
                    }
                }
            }

            // iterate through enriched basis functions
            for ( moris::size_t iBF = 0; iBF < tNumEnrBasis; iBF++ )
            {
                // iterate through enriched interpolation cells in the support
                moris_index tNumSubphaseInSupport = tCellsInEnrSupports( iBF ).size();

                moris_index tExpectedBulkPhase = MORIS_INDEX_MAX;

                for ( moris::moris_index iSP = 0; iSP < tNumSubphaseInSupport; iSP++ )
                {
                    Interpolation_Cell_Unzipped* tIpCell = tCellsInEnrSupports( iBF )( iSP );

                    moris_index tBulkPhase = tIpCell->get_bulkphase_index();

                    MORIS_ASSERT( tBulkPhase != MORIS_INDEX_MAX, "Bulk phase index not set." );

                    if ( iSP == 0 )
                    {
                        tExpectedBulkPhase = tBulkPhase;
                    }

                    if ( tBulkPhase != tExpectedBulkPhase )
                    {
                        MORIS_LOG_WARNING(
                                "Enr. BF index = %ld | SP index = %d  | tExpectedBulkPhase = %d | tBulkPhase = %d",
                                iBF,
                                tCellsInEnrSupports( iBF )( iSP )->get_index(),
                                tExpectedBulkPhase,
                                tBulkPhase );

                        tSubphaseBulkPhasesInSupportDiag = false;
                    }
                }
            }
        }

        // print error message if check fails
        MORIS_ERROR( tSubphaseBulkPhasesInSupportDiag,
                "Enriched_Interpolation_Mesh::verify_basis_support() - "
                "Bulk phases of IG cells in enriched basis support do not match" );

        // return true if check does not fail before
        return true;
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_local_to_global_maps()
    {
        // initialize local to global maps
        mLocalToGlobalMaps = Vector< Matrix< IdMat > >( 4 );
        mGlobalToLocalMaps = Vector< std::unordered_map< moris_id, moris_index > >( 4 );

        this->setup_cell_maps();

        this->setup_basis_maps();

        this->assign_ip_vertex_ids();

        this->setup_vertex_maps();
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::assign_ip_vertex_ids()
    {
        // log this function when verbose output is requested
        Tracer tTracer( "XTK", "Enriched Interpolation Mesh", "assign unzipped vertex IDs", mXTKModel->mVerboseLevel, 1 );

// debug
tic tTimer;

        /* ---------------------------------------------------------------------------------------- */
        /* Step 0: Sort into owned and not owned entities */

        // initialize a temporary list storing the UIPCs associated with the UIPVs
        Vector< moris_index > tUipcsAssociatedWithNotOwnedUipvs;

        // determine which UIPVs can be assigned an ID and which need to be communicated
        this->sort_unzipped_vertices_into_owned_and_not_owned( tUipcsAssociatedWithNotOwnedUipvs );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 0.5: Get the communication table */

        // get the communication table and map
        Matrix< IdMat > tCommTable     = mXTKModel->get_communication_table();
        uint            tCommTableSize = tCommTable.numel();

        std::map< moris_id, moris_index > tProcIdToCommTableIndex = mXTKModel->get_communication_map();

        /* ---------------------------------------------------------------------------------------- */
        /* Step 1: Let each proc decide how many entity IDs it needs & communicate ID ranges */

        // reserve IDs for his proc
        moris_id tMyFirstId = get_processor_offset( mOwnedUnzippedVertices.size() ) + 1;

        /* ---------------------------------------------------------------------------------------- */
        /* Step 2: Assign IDs to owned entities */

        // iterate through vertices that the current proc owns and assign a node id to them
        for ( uint iVert = 0; iVert < mOwnedUnzippedVertices.size(); iVert++ )
        {
            mEnrichedInterpVerts( mOwnedUnzippedVertices( iVert ) )->set_vertex_id( tMyFirstId );
            tMyFirstId++;
        }

        /* ---------------------------------------------------------------------------------------- */
        /* The following steps are only necessary if code runs in parallel */

        if ( par_size() == 1 )    // serial
        {
            // check that all entities are owned in serial
            MORIS_ASSERT( mNotOwnedUnzippedVertices.size() == 0,
                    "Enriched_Interpolation_Mesh::assign_ip_vertex_ids() - "
                    "Code running in serial, but not all UIPVs are owned by proc 0." );
        }
        else    // parallel
        {
            // check that NOT all entities are owned in parallel
            MORIS_ASSERT( mNotOwnedUnzippedVertices.size() > 0,
                    "Enriched_Interpolation_Mesh::assign_ip_vertex_ids() - "
                    "Code running in parallel, but all UIPVs are owned by current proc #%i.",
                    par_rank() );

            /* ---------------------------------------------------------------------------------------- */
            /* Step 3: Prepare requests for non-owned entities */

            // initialize lists of information that identifies entities (on other procs)
            Vector< Vector< moris_index > > tNotOwnedUIPVsToProcs;    // UIPV indices for communication (local to current proc, just used for construction of arrays)
            Vector< Matrix< IdMat > >       tBaseVertexIds;           // base vertex's ID the UIPVs live on
            Vector< Matrix< IdMat > >       tUnzippedIpCellIds;       // UIPC IDs the UIPVs belong to

            // fill identifying information
            this->prepare_requests_for_not_owned_unzipped_vertex_IDs(
                    tUipcsAssociatedWithNotOwnedUipvs,
                    tNotOwnedUIPVsToProcs,
                    tBaseVertexIds,
                    tUnzippedIpCellIds );

            /* ---------------------------------------------------------------------------------------- */
            /* Step 4: Send and Receive requests about non-owned entities to and from other procs */

            // initialize arrays for receiving
            Vector< Matrix< IdMat > > tReceivedBaseVertexIds;
            Vector< Matrix< IdMat > > tReceivedUnzippedIpCellIds;

            // communicate information
            moris::communicate_mats( tCommTable, tBaseVertexIds, tReceivedBaseVertexIds );
            moris::communicate_mats( tCommTable, tUnzippedIpCellIds, tReceivedUnzippedIpCellIds );

            // clear memory not needed anymore
            tBaseVertexIds.clear();
            tUnzippedIpCellIds.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 5: Find answers to the requests */

            // initialize lists of ID answers to other procs
            Vector< Matrix< IdMat > > tVertIds( tCommTableSize );

            this->prepare_answers_for_owned_unzipped_vertex_IDs( tVertIds, tReceivedBaseVertexIds, tReceivedUnzippedIpCellIds );

            // clear memory from requests (the answers to which have been found)
            tReceivedBaseVertexIds.clear();
            tReceivedUnzippedIpCellIds.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 6: Send and receive answers to and from other procs */

            // initialize arrays for receiving
            Vector< Matrix< IdMat > > tReceivedVertIds;

            // communicate answers
            moris::communicate_mats( tCommTable, tVertIds, tReceivedVertIds );

            // clear unused memory
            tVertIds.clear();

            /* ---------------------------------------------------------------------------------------- */
            /* Step 7: Use answers to assign IDs to non-owned entities */

            this->handle_requested_unzipped_vertex_ID_answers( tNotOwnedUIPVsToProcs, tReceivedVertIds );

        }    // end if: parallel

// debug
real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;
std::cout << "Proc #" << par_rank() << ": Enriched_Interpolation_Mesh::assign_ip_vertex_ids took " << tElapsedTime / 1000.0 << " seconds." << std::endl;

    }    // end function: Enriched_Interpolation_Mesh::assign_ip_vertex_ids

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_basis_ownership()
    {
        // size data
        mEnrichCoeffOwnership.resize( mMeshIndices.max() + 1 );

        // iterate through meshes
        for ( uint iM = 0; iM < mMeshIndices.numel(); iM++ )
        {
            moris_index tMeshIndex = mMeshIndices( iM );

            mEnrichCoeffOwnership( tMeshIndex ).resize( 1, mEnrichCoeffLocToGlob( tMeshIndex ).numel() );

            mEnrichCoeffOwnership( tMeshIndex ).fill( MORIS_INDEX_MAX );

            // iterate through basis functions
            for ( uint iB = 0; iB < mCoeffToEnrichCoeffs( tMeshIndex ).size(); iB++ )
            {
                moris_index tOwner = mXTKModel->get_background_mesh().get_entity_owner( (moris_index)iB, mBasisRank, tMeshIndex );

                for ( uint iEB = 0; iEB < mCoeffToEnrichCoeffs( tMeshIndex )( iB ).numel(); iEB++ )
                {
                    moris_index tEnrIndex = mCoeffToEnrichCoeffs( tMeshIndex )( iB )( iEB );

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
    Enriched_Interpolation_Mesh::setup_vertex_maps()
    {
        uint tNumNodes = this->get_num_entities( mtk::EntityRank::NODE );

        mLocalToGlobalMaps( 0 ) = Matrix< IdMat >( tNumNodes, 1 );

        for ( uint i = 0; i < tNumNodes; i++ )
        {
            mLocalToGlobalMaps( 0 )( mEnrichedInterpVerts( i )->get_index() ) = mEnrichedInterpVerts( i )->get_id();

            MORIS_ASSERT( mEnrichedInterpVerts( i )->get_index() == (moris_index)i, "Index alignment issue in vertices" );

            MORIS_ASSERT( mGlobalToLocalMaps( 0 ).find( mEnrichedInterpVerts( i )->get_id() ) == mGlobalToLocalMaps( 0 ).end(),
                    "Duplicate id in the vertex map detected" );

            mGlobalToLocalMaps( 0 )[ mEnrichedInterpVerts( i )->get_id() ] = i;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_basis_maps()
    {
        mGlobalToLocalBasisMaps.resize( mMeshIndices.max() + 1 );

        // iterate through meshes
        for ( uint iM = 0; iM < mEnrichCoeffLocToGlob.size(); iM++ )
        {
            for ( uint iB = 0; iB < mEnrichCoeffLocToGlob( iM ).numel(); iB++ )
            {
                // MORIS_LOG_SPEC("mEnrichCoeffLocToGlob(iM)(iB)",mEnrichCoeffLocToGlob(iM)(iB));
                // MORIS_ASSERT(mGlobalToLocalBasisMaps(iM).find(mEnrichCoeffLocToGlob(iM)(iB)) == mGlobalToLocalBasisMaps(iM).end(),
                //         "Duplicate id in the basis map detected");

                if ( mGlobalToLocalBasisMaps( iM ).find( mEnrichCoeffLocToGlob( iM )( iB ) ) == mGlobalToLocalBasisMaps( iM ).end() )
                {
                    mGlobalToLocalBasisMaps( iM )[ mEnrichCoeffLocToGlob( iM )( iB ) ] = (moris_index)iB;
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
    Enriched_Interpolation_Mesh::sort_unzipped_vertices_into_owned_and_not_owned(
            Vector< moris_index >& aUipcsAssociatedWithNotOwnedUipvs )
    {
        // get the number of enriched interpolation cells
        uint tNumUIPCs = this->get_num_entities( mtk::EntityRank::ELEMENT );
        uint tNumUIPVs = this->get_num_entities( mtk::EntityRank::NODE );

        // reserve memory
        mOwnedUnzippedVertices.reserve( tNumUIPVs );
        mNotOwnedUnzippedVertices.reserve( tNumUIPVs );
        aUipcsAssociatedWithNotOwnedUipvs.reserve( tNumUIPVs );

        // Keep track of which vertices have been treated
        Vector< bool > tVertexTracker( tNumUIPVs, true );

        // go over all enriched cells and their respective vertices to collect all owned and non-owned UIPVs
        for ( uint iCell = 0; iCell < tNumUIPCs; iCell++ )
        {
            // get access to the UIPC
            Interpolation_Cell_Unzipped const * tUIPC      = mEnrichedInterpCells( (moris_index)iCell );
            moris_index                         tUipcIndex = tUIPC->get_index();

            // get access to the UIPVs living on this UIPC
            Vector< xtk::Interpolation_Vertex_Unzipped* > const & tUIPVs          = tUIPC->get_xtk_interpolation_vertices();
            uint                                                  tNumUIPVsOnCell = tUIPVs.size();

            // go over the vertices and sort them into owned and not owned
            for ( uint iVert = 0; iVert < tNumUIPVsOnCell; iVert++ )
            {
                // get the index and owner of the current UIPV
                moris_index tVertexIndex = tUIPVs( iVert )->get_index();
                moris_index tOwner       = tUIPVs( iVert )->get_owner();

                // make sure vertices are not re-numbered
                MORIS_ASSERT( tVertexTracker( tVertexIndex ),
                        "Enriched_Interpolation_Mesh::assign_ip_vertex_ids() - "
                        "Trying to assign an index to a UIPV that has already been assigned one." );
                tVertexTracker( tVertexIndex ) = false;

                // sort into owned and not owned
                if ( tOwner == par_rank() )    // owned
                {
                    mOwnedUnzippedVertices.push_back( tVertexIndex );
                }
                else    // not owned
                {
                    mNotOwnedUnzippedVertices.push_back( tVertexIndex );
                    aUipcsAssociatedWithNotOwnedUipvs.push_back( tUipcIndex );
                }
            }
        }

        // size out unused memory
        mOwnedUnzippedVertices.shrink_to_fit();
        mNotOwnedUnzippedVertices.shrink_to_fit();
        aUipcsAssociatedWithNotOwnedUipvs.shrink_to_fit();

    }    // end function: Enriched_Interpolation_Mesh::sort_unzipped_vertices_into_owned_and_not_owned()

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::prepare_requests_for_not_owned_unzipped_vertex_IDs(
            Vector< moris_index > const &    aUipcsAssociatedWithNotOwnedUipvs,
            Vector< Vector< moris_index > >& aNotOwnedUIPVsToProcs,
            Vector< Matrix< IdMat > >&       aBaseVertexIds,
            Vector< Matrix< IdMat > >&       aUnzippedIpCellIds )
    {
        // get the communication table and map
        Matrix< IdMat > tCommTable     = mXTKModel->get_communication_table();
        uint            tCommTableSize = tCommTable.numel();

        std::map< moris_id, moris_index > tProcIdToCommTableIndex = mXTKModel->get_communication_map();

        // initialize lists of identifying information
        aNotOwnedUIPVsToProcs.resize( tCommTableSize );
        aBaseVertexIds.resize( tCommTableSize );
        aUnzippedIpCellIds.resize( tCommTableSize );

        // get the number of non-owned entities on the executing processor
        uint tNumNotOwnedUIPVs = mNotOwnedUnzippedVertices.size();

        // prepare list that give the position of the requested UIPV in the array of non-owned UIPVs
        Vector< Vector< moris_index > > tUipvPositionInNotOwnedList( tCommTableSize );

        // go through SPGs that executing proc knows about, but doesn't own, ...
        for ( uint iNotOwnedVert = 0; iNotOwnedVert < tNumNotOwnedUIPVs; iNotOwnedVert++ )
        {
            // ... get their index ...
            moris_index tVertIndex = mNotOwnedUnzippedVertices( iNotOwnedVert );

            // ... get their respective owners, and position in the comm table ...
            moris_index tOwnerProc = mEnrichedInterpVerts( tVertIndex )->get_owner();
            auto        tIter      = tProcIdToCommTableIndex.find( tOwnerProc );
            MORIS_ASSERT(
                    tIter != tProcIdToCommTableIndex.end(),
                    "Enriched_Interpolation_Mesh::prepare_requests_for_not_owned_unzipped_vertex_IDs() - "
                    "Entity owner (Proc #%i) not found in communication table of current proc #%i which is: %s",
                    tOwnerProc,
                    par_rank(),
                    ios::stringify_log( tCommTable ).c_str() );

            moris_index tProcDataIndex = tIter->second;

            // ... and finally add the non-owned SPGs in the list of SPs to be requested from that owning proc
            aNotOwnedUIPVsToProcs( tProcDataIndex ).push_back( tVertIndex );

            // store the where the current UIPV can be found in the list of not owned UIPVs
            tUipvPositionInNotOwnedList( tProcDataIndex ).push_back( iNotOwnedVert );
        }

        // size out unused memory
        aNotOwnedUIPVsToProcs.shrink_to_fit();
        tUipvPositionInNotOwnedList.shrink_to_fit();

        // assemble identifying information for every processor communicated with
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of non-owned entities to be sent to each processor processor
            uint tNumNotOwnedEntitiesOnProc = aNotOwnedUIPVsToProcs( iProc ).size();

            // allocate matrix
            aBaseVertexIds( iProc ).resize( tNumNotOwnedEntitiesOnProc, 1 );
            aUnzippedIpCellIds( iProc ).resize( tNumNotOwnedEntitiesOnProc, 1 );

            // go through the Subphase groups for which IDs will be requested by the other processor
            for ( uint iVert = 0; iVert < tNumNotOwnedEntitiesOnProc; iVert++ )
            {
                // get the index of the UIPV on the executing proc and its ID
                moris_index tVertIndex = aNotOwnedUIPVsToProcs( iProc )( iVert );
                moris_id    tVertId    = mEnrichedInterpVerts( tVertIndex )->get_base_vertex()->get_id();

                // get the index and ID of the UIPC the current UIPV is attached to
                moris_index tIndexInNotOwnedList = tUipvPositionInNotOwnedList( iProc )( iVert );
                moris_index tCellIndex           = aUipcsAssociatedWithNotOwnedUipvs( tIndexInNotOwnedList );
                moris_id    tCellId              = mEnrichedInterpCells( tCellIndex )->get_id();

                // store the identifying information in the output arrays
                aBaseVertexIds( iProc )( iVert )     = tVertId;
                aUnzippedIpCellIds( iProc )( iVert ) = tCellId;
            }

        }    // end for: each proc communicated with

    }    // end function: Enriched_Interpolation_Mesh::prepare_requests_for_not_owned_unzipped_vertex_IDs()

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::prepare_answers_for_owned_unzipped_vertex_IDs(
            Vector< Matrix< IdMat > >&        aVertIds,
            Vector< Matrix< IdMat > > const & aReceivedBaseVertexIds,
            Vector< Matrix< IdMat > > const & aReceivedUnzippedIpCellIds )
    {
        // get the communication table
        Matrix< IdMat > tCommTable     = mXTKModel->get_communication_table();
        uint            tCommTableSize = tCommTable.numel();

        // initialize answer array with correct size
        aVertIds.resize( tCommTableSize );

        // check that the received data is complete
        MORIS_ASSERT(
                aReceivedBaseVertexIds.size() == tCommTableSize && aReceivedUnzippedIpCellIds.size() == tCommTableSize,
                "Enriched_Interpolation_Mesh::prepare_answers_for_owned_unzipped_vertex_IDs() - "
                "Received information incomplete." );

        // go through the list of processors in the array of ID requests
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of entity IDs requested from the current proc position
            uint tNumReceivedReqs = aReceivedBaseVertexIds( iProc ).numel();

            // size the list of answers / IDs accordingly
            aVertIds( iProc ).resize( 1, tNumReceivedReqs );

            // iterate through the entities for which the IDs are requested
            for ( uint iVert = 0; iVert < tNumReceivedReqs; iVert++ )
            {
                // get the ID of the received base vertex
                moris_id tBaseVertexId = aReceivedBaseVertexIds( iProc )( iVert );

                // get the the UIPC
                moris_id    tUipcId    = aReceivedUnzippedIpCellIds( iProc )( iVert );
                moris_index tUipcIndex = this->get_loc_entity_ind_from_entity_glb_id( tUipcId, mtk::EntityRank::ELEMENT );

                Interpolation_Cell_Unzipped* tIpCell = mEnrichedInterpCells( tUipcIndex );

                // ge the vertices that are attached to the unzipped IP cell
                Vector< xtk::Interpolation_Vertex_Unzipped* > const & tVertsOnCell =
                        tIpCell->get_xtk_interpolation_vertices();

                uint tNumVertsOnCell = tVertsOnCell.size();

                // check which of the vertices is the one requested
                bool tFound = false;
                for ( uint iVertOnCell = 0; iVertOnCell < tNumVertsOnCell; iVertOnCell++ )
                {
                    // get access to the the base vertex and its id to test against
                    xtk::Interpolation_Vertex_Unzipped const * tUIPV = tVertsOnCell( iVertOnCell );

                    moris_id tBaseVertexOnCellId = tUIPV->get_base_vertex()->get_id();

                    // check if this is the one we're looking for
                    if ( tBaseVertexOnCellId == tBaseVertexId )
                    {
                        // store ID answer
                        aVertIds( iProc )( iVert ) = tUIPV->get_id();

                        // mark as found
                        tFound = true;

                        // stop loop over vertices on IP cell to save time
                        break;
                    }
                }

                // make sure that an answer has been found
                MORIS_ERROR( tFound,
                        "Enriched_Interpolation_Mesh::prepare_answers_for_owned_unzipped_vertex_IDs() - "
                        "No unzipped vertex with base vertex ID %i, requested by proc #%i, has been found on this UIPC (ID: %i, #%i).",
                        tBaseVertexId,
                        tCommTable( iProc ),
                        tUipcId,
                        tUipcIndex );

            }    // end for: communication for each entity with current processor

        }    // end for: communication list for each processor

    }    // end function: Enriched_Interpolation_Mesh::prepare_answers_for_owned_unzipped_vertex_IDs()

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::handle_requested_unzipped_vertex_ID_answers(
            Vector< Vector< moris_index > > const & tNotOwnedUIPVsToProcs,
            Vector< Matrix< IdMat > > const &       tReceivedVertIds )
    {
        // get the communication table
        Matrix< IdMat > tCommTable     = mXTKModel->get_communication_table();
        uint            tCommTableSize = tCommTable.numel();

        // process answers from each proc communicated with
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of requests and answers from the current proc
            uint tNumReceivedEntityIds = tReceivedVertIds( iProc ).numel();

            // make sure everything has been answered
            MORIS_ASSERT( tNumReceivedEntityIds == tNotOwnedUIPVsToProcs( iProc ).size(),
                    "Enriched_Interpolation_Mesh::handle_requested_unzipped_vertex_ID_answers() - "
                    "Request arrays sent to and answers received from proc #%i have different size.",
                    tCommTable( iProc ) );

            // assign IDs to each communicated entity
            for ( uint iVert = 0; iVert < tNumReceivedEntityIds; iVert++ )
            {
                // get the current SPG index and ID from the data provided
                moris_index tUipvIndex = tNotOwnedUIPVsToProcs( iProc )( iVert );
                moris_id    tUipvId    = tReceivedVertIds( iProc )( iVert );

                // store the received entity ID
                mEnrichedInterpVerts( tUipvIndex )->set_vertex_id( tUipvId );
            }

        }    // end for: each processor communicated with

    }    // end function: Enriched_Interpolation_Mesh::handle_requested_unzipped_vertex_ID_answers()

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_cell_maps()
    {
        uint tNumCells = this->get_num_entities( mtk::EntityRank::ELEMENT );

        mLocalToGlobalMaps( 3 ) = Matrix< IdMat >( tNumCells, 1 );

        //        mGlobalToLocalMaps(3).clear();

        for ( uint iCell = 0; iCell < tNumCells; iCell++ )
        {
            mLocalToGlobalMaps( 3 )( mEnrichedInterpCells( iCell )->get_index() ) = mEnrichedInterpCells( iCell )->get_id();

            MORIS_ASSERT( mEnrichedInterpCells( iCell )->get_index() == (moris_index)iCell,
                    "Enriched_Interpolation_Mesh::setup_cell_maps() - Index alignment issue in cells" );

            MORIS_ASSERT( mGlobalToLocalMaps( 3 ).find( mEnrichedInterpCells( iCell )->get_id() ) == mGlobalToLocalMaps( 3 ).end(),
                    "Enriched_Interpolation_Mesh::setup_cell_maps() - Duplicate id in the cell map detected" );

            mGlobalToLocalMaps( 3 )[ mEnrichedInterpCells( iCell )->get_id() ] = mEnrichedInterpCells( iCell )->get_index();
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::setup_mesh_index_map()
    {
        for ( uint iMeshIndex = 0; iMeshIndex < (uint)mMeshIndices.max() + 1; iMeshIndex++ )
        {
            mMeshIndexToLocMeshIndex[ iMeshIndex ] = iMeshIndex;
        }
    }

    // ----------------------------------------------------------------------------

    moris_id
    Enriched_Interpolation_Mesh::allocate_entity_ids(
            moris::size_t   aNumReqs,
            mtk::EntityRank aEntityRank,
            bool            aStartFresh )
    {
        MORIS_ASSERT( aEntityRank == mtk::EntityRank::NODE || aEntityRank == mtk::EntityRank::ELEMENT,
                "Only Elements or Nodes have ids" );

        moris_id tGlobalMax = 1;
        if ( !aStartFresh )
        {
            this->get_max_entity_id( aEntityRank );
        }

        int tProcRank = par_rank();
        int tProcSize = par_size();

        Vector< moris::moris_id > aGatheredInfo;
        Vector< moris::moris_id > tFirstId( 1 );
        Vector< moris::moris_id > tNumIdsRequested( 1 );

        tNumIdsRequested( 0 ) = (moris::moris_id)aNumReqs;

        moris::gather_vector( tNumIdsRequested, aGatheredInfo );

        Vector< moris::moris_id > tProcFirstID( tProcSize );

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

        moris::scatter_vector( tProcFirstID, tFirstId );

        return tFirstId( 0 );
    }

    // ----------------------------------------------------------------------------
    // multi-grid accessor functions
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

    Matrix< DDSMat >
    Enriched_Interpolation_Mesh::get_fine_basis_inds_of_basis(
            const moris_index aInterpolationIndex,
            const moris_index aBasisIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_fine_basis_inds_of_basis(
                aInterpolationIndex,
                aBasisIndex );
    }

    // ----------------------------------------------------------------------------

    Matrix< DDRMat >
    Enriched_Interpolation_Mesh::get_fine_basis_weights_of_basis(
            const moris_index aInterpolationIndex,
            const moris_index aBasisIndex )
    {
        return mXTKModel->get_multigrid_ptr()->get_fine_basis_weights_of_basis(
                aInterpolationIndex,
                aBasisIndex );
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::determine_unenriched_meshes_are_enriched_beforehand() const
    {
        // loop over the unenriched meshes to see if they have been enriched before
        for ( auto const & iUnenrichedMeshIndex : mUnenrichedMeshIndices )
        {
            // check if unenriched mesh index, exits in the mesh indices ( enriched meshes)
            bool tMeshIsUnenriched =
                    std::any_of( mMeshIndices.cbegin(), mMeshIndices.cend(), [ &iUnenrichedMeshIndex ]( moris_index aMeshIndex )    //
                            { return iUnenrichedMeshIndex == aMeshIndex; } );

            // throw an error specifying which mesh number is not enriched before
            MORIS_ERROR( tMeshIsUnenriched, "Mesh %u is not enriched beforehand", iUnenrichedMeshIndex );
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::override_maps()
    {
        // loop over the meshes that will be unenriched
        for ( auto const & iMeshIndex : mUnenrichedMeshIndices )
        {
            // get the local mesh index
            moris::moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( iMeshIndex );

            // get the global to local basis map from HMR that corresponds to the unenriched version
            map< moris_id, moris_index > tGlobalToLocalHMRBasisMap;
            mXTKModel->mBackgroundMesh->get_adof_map( iMeshIndex, tGlobalToLocalHMRBasisMap );

            // clear the maps and cells that need to be overwritten
            mGlobalToLocalBasisMaps( tLocalMeshIndex ).clear();
            mEnrichCoeffLocToGlob( tLocalMeshIndex ).set_size( 1, tGlobalToLocalHMRBasisMap.size(), MORIS_INDEX_MAX );
            mEnrichCoeffOwnership( tLocalMeshIndex ).set_size( 1, tGlobalToLocalHMRBasisMap.size(), MORIS_INDEX_MAX );

            // loop over the global to local hmr map and fill out the xtk maps
            for ( auto const & iGlobalToLocal : tGlobalToLocalHMRBasisMap )
            {
                // set the local to global and global to local maps
                mGlobalToLocalBasisMaps( tLocalMeshIndex )[ iGlobalToLocal.first ] = iGlobalToLocal.second;
                mEnrichCoeffLocToGlob( tLocalMeshIndex )( iGlobalToLocal.second )  = iGlobalToLocal.first;

                // get the owner of the non-enriched basis owner and store it
                moris_index tOwner = mXTKModel->get_background_mesh().get_entity_owner( iGlobalToLocal.second, mBasisRank, iMeshIndex );

                mEnrichCoeffOwnership( tLocalMeshIndex )( iGlobalToLocal.second ) = tOwner;
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::set_unenriched_mesh_indices( Matrix< IndexMat > const & aMeshIndices )
    {
        mUnenrichedMeshIndices = aMeshIndices;

        // check if the mesh index is compatible
        this->determine_unenriched_meshes_are_enriched_beforehand();
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::override_vertex_enrichment_id_index()
    {
        // loop over the unenriched mesh indices
        for ( auto const & iMeshIndex : mUnenrichedMeshIndices )
        {
            moris_index tLocalMeshIndex = this->get_local_mesh_index_xtk( iMeshIndex );

            // loop over the vertex enrichments to change their id and index
            for ( Vertex_Enrichment* iVertexEnrichment : mInterpVertEnrichment( tLocalMeshIndex ) )
            {
                // if it is not empty
                if ( iVertexEnrichment != nullptr )
                {
                    // and if it has interpolation basis , it is assumed that it has a base vertex
                    if ( iVertexEnrichment->has_interpolation() )
                    {
                        // get the base vertex interpolation
                        mtk::Vertex_Interpolation const * tBaseVertexInterpolation = iVertexEnrichment->get_base_vertex_interpolation();

                        // extract basis indices of the base one
                        Matrix< IndexMat > tBaseCoeffInds = tBaseVertexInterpolation->get_indices();

                        // get access to the basis to local index map of the vertex enrichment for modification
                        IndexMap& tVertEnrichMap = iVertexEnrichment->get_basis_map();

                        // clear the map and populate it with the new indices
                        tVertEnrichMap.clear();
                        for ( uint iBC = 0; iBC < tBaseCoeffInds.numel(); iBC++ )
                        {
                            moris::moris_index tBasisIndex = tBaseCoeffInds( iBC );
                            tVertEnrichMap[ tBasisIndex ]  = iBC;
                        }

                        // replace the index and the id of the vertex interpolations with the base one
                        iVertexEnrichment->add_basis_information( tBaseCoeffInds, tBaseVertexInterpolation->get_ids() );
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------

    Vector< std::string >
    Enriched_Interpolation_Mesh::get_set_names( mtk::EntityRank aSetEntityRank ) const
    {
        switch ( aSetEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                return Vector< std::string >( 0 );
                break;
            }
            case mtk::EntityRank::EDGE:
            {
                return Vector< std::string >( 0 );
                break;
            }
            case mtk::EntityRank::FACE:
            {
                return Vector< std::string >( 0 );
                ;
                break;
            }
            case mtk::EntityRank::ELEMENT:
            {
                return mBlockSetNames;
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Currently only supporting block, node and side sets in XTK enriched integration meshes" );
            }
                return Vector< std::string >( 0 );
                break;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::create_set_names()
    {
        // get number of phases and backgroud basis
        uint tNumPhases = mXTKModel->mGeometryEngine->get_num_bulk_phase();

        Vector< std::string > tBlockSetsNames = mXTKModel->mBackgroundMesh->get_set_names( mtk::EntityRank::ELEMENT );

        // populate the names of the sets
        mBlockSetNames.resize( tNumPhases );
        for ( uint i = 0; i < tNumPhases; i++ )
        {
            mBlockSetNames( i ) = tBlockSetsNames( 0 ) + "_p" + std::to_string( i );
        }

        // this data strcuture is needed to create sets in exodus writer
        mElementIndicesInBlock.resize( tNumPhases );

        // reserve space for each block set
        for ( int iBlockIndex = 0; iBlockIndex < (int)tNumPhases; iBlockIndex++ )
        {
            mElementIndicesInBlock( iBlockIndex ).reserve( mEnrichedInterpCells.size() );
        }

        // loop over the UIPC and put them them in the corret set based on the bulk phase index
        for ( const auto& iEnrIPCell : mEnrichedInterpCells )
        {
            // get bulk phase index and index of the UIPC
            moris_index tBulkPhaseIndex = iEnrIPCell->get_bulkphase_index();
            moris_index tIPCellIndex    = iEnrIPCell->get_index();

            // insert the UPIC at the corret set
            mElementIndicesInBlock( tBulkPhaseIndex ).push_back( tIPCellIndex );
        }

        // shrink all the internal cells
        shrink_to_fit_all( mElementIndicesInBlock );
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::create_basis_support_fields( Matrix< DDRMat > const & aProbeSpheres )
    {

// DEBUG because basis coordinates is only defined on debug mode
#ifdef MORIS_HAVE_DEBUG

        MORIS_ASSERT( aProbeSpheres.n_cols() == 4, "Probe sphere should be r, xc, yc, zc" );
        moris_index tNumSpheres = aProbeSpheres.n_rows();

        // background mesh data
        moris::mtk::Interpolation_Mesh& tMeshData = mXTKModel->get_background_mesh();

        // base string of field
        std::string tBaseStr = "weights";

        // determine which basis functions we are visualizing
        Vector< Vector< moris_index > > tActiveBasis( mMeshIndices.numel() );

        Vector< std::unordered_map< moris_index, moris_index > > tEnrCoeffActiveIndexFieldIndex( mMeshIndices.numel() );

        moris_index tFieldIndex = 0;

        for ( uint iBT = 0; iBT < mMeshIndices.numel(); iBT++ )
        {
            moris_index tMeshIndex = iBT;

            // iterate through background basis functions
            for ( uint iBackBasisIndex = 0; iBackBasisIndex < this->get_num_background_coefficients( tMeshIndex ); iBackBasisIndex++ )
            {
                // get the basis coordinate of the background basis function
                Matrix< DDRMat > tBasisCoords = tMeshData.get_basis_coords( tMeshIndex, (moris_index)iBackBasisIndex );

                // iterate through circles, see if the basis is active
                for ( moris_index iSp = 0; iSp < tNumSpheres; iSp++ )
                {
                    // initialize  the level set value
                    moris::real tLSVal = 0;

                    // determine the value of level set at the basis coordinates based on the dimension
                    if ( this->get_spatial_dim() == 3 )
                    {
                        tLSVal = sqrt( pow( tBasisCoords( 0 ) - aProbeSpheres( iSp, 1 ), 2 ) + pow( tBasisCoords( 1 ) - aProbeSpheres( iSp, 2 ), 2 ) + pow( tBasisCoords( 2 ) - aProbeSpheres( iSp, 3 ), 2 ) ) - aProbeSpheres( iSp, 0 );
                    }
                    else
                    {
                        tLSVal = sqrt( pow( tBasisCoords( 0 ) - aProbeSpheres( iSp, 1 ), 2 ) + pow( tBasisCoords( 1 ) - aProbeSpheres( iSp, 2 ), 2 ) ) - aProbeSpheres( iSp, 0 );
                    }

                    if ( tLSVal < 0.0 )
                    {
                        // iterate through enriched interpolation coeffs
                        Matrix< IndexMat > const & tEnrCoeffs = this->get_enriched_coefficients_at_background_coefficient( tMeshIndex, (moris_index)iBackBasisIndex );

                        for ( uint iEnrBasisOrd = 0; iEnrBasisOrd < tEnrCoeffs.numel(); iEnrBasisOrd++ )
                        {
                            const moris_index tEnrIndex = tEnrCoeffs( iEnrBasisOrd );

                            tActiveBasis( tMeshIndex ).push_back( tEnrIndex );
                            tEnrCoeffActiveIndexFieldIndex( tMeshIndex )[ tEnrIndex ] = tFieldIndex;
                            tFieldIndex++;
                        }
                    }
                }
            }
        }

        // field information for internal use
        Vector< std::string >      tFieldNames( tFieldIndex );
        Vector< moris_index >      tFieldIndices( tFieldIndex );
        Vector< Matrix< DDRMat > > tFieldData( tFieldIndex, Matrix< DDRMat >( 1, this->get_num_nodes(), 0.0 ) );

        // iterate through interpolation types and for each basis declare the field in mesh
        for ( uint iBT = 0; iBT < this->get_num_interpolation_types(); iBT++ )
        {

            moris_index tMeshIndex     = iBT;
            std::string tInterpTypeStr = "_mi_" + std::to_string( tMeshIndex );

            // iterate through basis functions
            for ( uint iB = 0; iB < tActiveBasis( tMeshIndex ).size(); iB++ )
            {
                MORIS_ASSERT( tEnrCoeffActiveIndexFieldIndex( tMeshIndex ).find( tActiveBasis( tMeshIndex )( iB ) ) != tEnrCoeffActiveIndexFieldIndex( tMeshIndex ).end(), "Not in map" );
                tFieldIndex = tEnrCoeffActiveIndexFieldIndex( tMeshIndex ).find( tActiveBasis( tMeshIndex )( iB ) )->second;

                tFieldNames( tFieldIndex ) = tBaseStr + tInterpTypeStr + "_ind_" + std::to_string( tActiveBasis( tMeshIndex )( iB ) );

                // declare the field in this mesh
                tFieldIndices( tFieldIndex ) = this->create_field( tFieldNames( tFieldIndex ), mtk::EntityRank::NODE, 0 );
            }
        }

        // iterate through primary cells
        for ( const auto& iUnzippedIPCell : mEnrichedInterpCells )
        {
            // get vertices attached to primary cells
            Vector< moris::mtk::Vertex* > tVertices = iUnzippedIPCell->get_vertex_pointers();

            // iterate through vertices and mark them as in support of all coefficients in tCoeffsIPIntoCluster
            for ( uint iV = 0; iV < tVertices.size(); iV++ )
            {
                for ( uint iBT = 0; iBT < this->get_num_interpolation_types(); iBT++ )
                {
                    Matrix< IndexMat >      tBasisIndices = tVertices( iV )->get_interpolation( iBT )->get_indices();
                    const Matrix< DDRMat >* tBasisWeights = tVertices( iV )->get_interpolation( iBT )->get_weights();

                    for ( uint iBasisOrd = 0; iBasisOrd < tBasisIndices.numel(); iBasisOrd++ )
                    {
                        auto tFieldIndIter = tEnrCoeffActiveIndexFieldIndex( iBT ).find( tBasisIndices( iBasisOrd ) );
                        if ( tFieldIndIter != tEnrCoeffActiveIndexFieldIndex( iBT ).end() )
                        {
                            tFieldData( tFieldIndIter->second )( tVertices( iV )->get_index() ) = ( *tBasisWeights )( iBasisOrd );
                        }
                    }
                }
            }
        }

        // add field data to mesh
        // iterate through interpolation
        for ( uint iField = 0; iField < tFieldIndices.size(); iField++ )
        {
            this->add_field_data( tFieldIndices( iField ), mtk::EntityRank::NODE, tFieldData( iField ) );
        }
#endif
    }

    // ----------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::create_basis_function_fields( Matrix< DDRMat > const & aProbeSpheres )
    {

// DEBUG because basis coordinates is only defined on debug mode
#ifdef MORIS_HAVE_DEBUG

        MORIS_ASSERT( aProbeSpheres.n_cols() == 4, "Probe sphere should be r, xc, yc, zc" );
        moris_index tNumSpheres = aProbeSpheres.n_rows();

        // background mesh data
        moris::mtk::Interpolation_Mesh& tMeshData = mXTKModel->get_background_mesh();

        // base string of field
        std::string tBaseStr = "Basis";

        // determine which basis functions we are visualizing
        Vector< Vector< moris_index > >                          tActiveBasis( mMeshIndices.numel() );
        Vector< std::unordered_map< moris_index, moris_index > > tEnrCoeffActiveIndexFieldIndex( mMeshIndices.numel() );

        moris_index tFieldIndex = 0;

        for ( uint iBT = 0; iBT < mMeshIndices.numel(); iBT++ )
        {
            moris_index tMeshIndex = iBT;

            // iterate through background basis functions
            for ( uint iBackBasisIndex = 0; iBackBasisIndex < this->get_num_background_coefficients( tMeshIndex ); iBackBasisIndex++ )
            {
                // get the basis coordinate of the background basis function
                Matrix< DDRMat > tBasisCoords = tMeshData.get_basis_coords( tMeshIndex, (moris_index)iBackBasisIndex );

                // iterate through circles, see if the basis is active
                for ( moris_index iSp = 0; iSp < tNumSpheres; iSp++ )
                {
                    // initialize  the level set value
                    moris::real tLSVal = 0;

                    // determine the value of level set at the basis coordinates based on the dimension
                    if ( this->get_spatial_dim() == 3 )
                    {
                        tLSVal = sqrt( pow( tBasisCoords( 0 ) - aProbeSpheres( iSp, 1 ), 2 ) + pow( tBasisCoords( 1 ) - aProbeSpheres( iSp, 2 ), 2 ) + pow( tBasisCoords( 2 ) - aProbeSpheres( iSp, 3 ), 2 ) ) - aProbeSpheres( iSp, 0 );
                    }
                    else
                    {
                        tLSVal = sqrt( pow( tBasisCoords( 0 ) - aProbeSpheres( iSp, 1 ), 2 ) + pow( tBasisCoords( 1 ) - aProbeSpheres( iSp, 2 ), 2 ) ) - aProbeSpheres( iSp, 0 );
                    }

                    if ( tLSVal < 0.0 )
                    {
                        // iterate through enriched interpolation coeffs
                        Matrix< IndexMat > const & tEnrCoeffs = this->get_enriched_coefficients_at_background_coefficient( tMeshIndex, (moris_index)iBackBasisIndex );

                        for ( uint iEnrBasisOrd = 0; iEnrBasisOrd < tEnrCoeffs.numel(); iEnrBasisOrd++ )
                        {
                            const moris_index tEnrIndex = tEnrCoeffs( iEnrBasisOrd );

                            tActiveBasis( tMeshIndex ).push_back( tEnrIndex );
                            tEnrCoeffActiveIndexFieldIndex( tMeshIndex )[ tEnrIndex ] = tFieldIndex;
                            tFieldIndex++;
                        }
                    }
                }
            }
        }

        // field information for internal use
        Vector< std::string >      tFieldNames( tFieldIndex );
        Vector< moris_index >      tFieldIndices( tFieldIndex );
        Vector< Matrix< DDRMat > > tFieldData( tFieldIndex, Matrix< DDRMat >( 1, this->get_num_nodes(), 0.0 ) );

        mtk::Mesh_Pair      tMeshPair( this, nullptr );
        mtk::Field_Discrete tFieldDiscrete( tMeshPair );

        // iterate through interpolation types and for each basis declare the field in mesh
        for ( uint iBT = 0; iBT < this->get_num_interpolation_types(); iBT++ )
        {

            moris_index tMeshIndex     = iBT;
            std::string tInterpTypeStr = "_mi_" + std::to_string( tMeshIndex );

            // iterate through basis functions
            for ( uint iB = 0; iB < tActiveBasis( tMeshIndex ).size(); iB++ )
            {
                MORIS_ASSERT( tEnrCoeffActiveIndexFieldIndex( tMeshIndex ).find( tActiveBasis( tMeshIndex )( iB ) ) != tEnrCoeffActiveIndexFieldIndex( tMeshIndex ).end(), "Not in map" );
                tFieldIndex = tEnrCoeffActiveIndexFieldIndex( tMeshIndex ).find( tActiveBasis( tMeshIndex )( iB ) )->second;

                tFieldNames( tFieldIndex ) = tBaseStr + tInterpTypeStr + "_ind_" + std::to_string( tActiveBasis( tMeshIndex )( iB ) );

                // declare the field in this mesh
                tFieldIndices( tFieldIndex ) = this->create_field( tFieldNames( tFieldIndex ), mtk::EntityRank::NODE, 0 );

                Matrix< DDRMat > tCoeffMatrix( this->get_max_num_coeffs_on_proc( tMeshIndex ), 1, 0.0 );
                tCoeffMatrix( tActiveBasis( tMeshIndex )( iB ) ) = 1.0;

                tFieldDiscrete.unlock_field();
                tFieldDiscrete.set_coefficients( tCoeffMatrix );
                tFieldDiscrete.compute_nodal_values();
                Matrix< DDRMat > const & tVals = tFieldDiscrete.get_values();

                tFieldData( tFieldIndex ) = tVals;
            }
        }

        // add field data to mesh
        // iterate through interpolation
        for ( uint iField = 0; iField < tFieldIndices.size(); iField++ )
        {
            this->add_field_data( tFieldIndices( iField ), mtk::EntityRank::NODE, tFieldData( iField ) );
        }
#endif
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Interpolation_Mesh::get_coefficient_indices_of_node(
            uint aNodeIndex,
            uint aDiscretizationMeshIndex )
    {
        mtk::Vertex_Interpolation* tVertexInterp = this->get_mtk_vertex( aNodeIndex ).get_interpolation( aDiscretizationMeshIndex );
        // std::cout << "Index Node: " << this->get_mtk_vertex( aNodeIndex ).get_index() << std::endl;
        // print_as_row_vector( tVertexInterp->get_indices(), "tVertexInterp->get_indices()" );
        return tVertexInterp->get_indices();
    }

    // ----------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Enriched_Interpolation_Mesh::get_t_matrix_of_node_loc_ind(
            uint aNodeIndex,
            uint aDiscretizationMeshIndex )
    {
        mtk::Vertex_Interpolation* tVertexInterp = this->get_mtk_vertex( aNodeIndex ).get_interpolation( aDiscretizationMeshIndex );
        // std::cout << "Index Node: " << this->get_mtk_vertex( aNodeIndex ).get_index() << std::endl;
        // print_as_row_vector( *tVertexInterp->get_weights(), "*tVertexInterp->get_weights()" );
        return *tVertexInterp->get_weights();
    }

    // ----------------------------------------------------------------------------

    Matrix< IdMat >
    Enriched_Interpolation_Mesh::get_coefficient_owners_of_node(
            uint aNodeIndex,
            uint aBSplineMeshIndex )
    {
        mtk::Vertex_Interpolation* tVertexInterp = this->get_mtk_vertex( aNodeIndex ).get_interpolation( aBSplineMeshIndex );
        // std::cout << "Index Node: " << this->get_mtk_vertex( aNodeIndex ).get_index() << std::endl;
        // print_as_row_vector( tVertexInterp->get_owners(), "tVertexInterp->get_owners()" );
        return tVertexInterp->get_owners();
    }

    // ----------------------------------------------------------------------------

    Matrix< IdMat >
    Enriched_Interpolation_Mesh::get_coefficient_IDs_of_node(
            uint aNodeIndex,
            uint aDiscretizationIndex )
    {
        mtk::Vertex_Interpolation* tVertexInterp = this->get_mtk_vertex( aNodeIndex ).get_interpolation( aDiscretizationIndex );
        // std::cout << "Index Node: " << this->get_mtk_vertex( aNodeIndex ).get_index() << std::endl;
        // print_as_row_vector( tVertexInterp->get_ids(), "tVertexInterp->get_ids()" );
        return tVertexInterp->get_ids();
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Interpolation_Mesh::get_entity_owner(
            moris_index       aEntityIndex,
            mtk::EntityRank   aEntityRank,
            const moris_index aDiscretizationMeshIndex ) const
    {
        switch ( aEntityRank )
        {
            case mtk::EntityRank::BSPLINE:
            {
                return mEnrichCoeffOwnership( aDiscretizationMeshIndex )( aEntityIndex );
                break;
            }

            case mtk::EntityRank::NODE:
            {
                return this->get_mtk_vertex( aEntityIndex ).get_owner();
            }

            default:
                return 0;
                break;
        }
    }

    // ----------------------------------------------------------------------------

    moris::moris_index
    Enriched_Interpolation_Mesh::get_field_index(
            std::string     aLabel,
            mtk::EntityRank aEntityRank )
    {
        MORIS_ASSERT( field_exists( aLabel, aEntityRank ), "Field does not exist in mesh" );

        moris_index tIndex = get_entity_rank_field_index( aEntityRank );
        auto        tIter  = mFieldLabelToIndex( tIndex ).find( aLabel );
        return tIter->second;
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::add_field_data(
            moris::moris_index       aFieldIndex,
            mtk::EntityRank          aEntityRank,
            Matrix< DDRMat > const & aFieldData )
    {
        mFields( aFieldIndex ).mFieldData = aFieldData.copy();
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > const &
    Enriched_Interpolation_Mesh::get_field_data(
            moris::moris_index aFieldIndex,
            mtk::EntityRank    aEntityRank ) const
    {
        return mFields( aFieldIndex ).mFieldData;
    }

    //------------------------------------------------------------------------------

    Vector< std::string >
    Enriched_Interpolation_Mesh::get_field_names( mtk::EntityRank aEntityRank )
    {
        Vector< std::string > tOutputFieldNames;

        moris_index tRankFieldIndex = this->get_entity_rank_field_index( aEntityRank );

        for ( auto const & iter : mFieldLabelToIndex( tRankFieldIndex ) )
        {
            tOutputFieldNames.push_back( iter.first );
        }

        return tOutputFieldNames;
    }

    //------------------------------------------------------------------------------

    moris::moris_index
    Enriched_Interpolation_Mesh::create_field(
            std::string        aLabel,
            mtk::EntityRank    aEntityRank,
            moris::moris_index aBulkPhaseIndex )
    {
        MORIS_ASSERT( !field_exists( aLabel, aEntityRank ), "Field already created" );

        moris::moris_index tFieldIndex                                                   = mFields.size();
        mFieldLabelToIndex( this->get_entity_rank_field_index( aEntityRank ) )[ aLabel ] = tFieldIndex;
        mFields.push_back( Field( aLabel, aBulkPhaseIndex ) );

        return tFieldIndex;
    }

    //------------------------------------------------------------------------------

    bool
    Enriched_Interpolation_Mesh::field_exists(
            std::string     aLabel,
            mtk::EntityRank aEntityRank )
    {
        moris::moris_index tIndex = this->get_entity_rank_field_index( aEntityRank );

        return mFieldLabelToIndex( tIndex ).find( aLabel ) != mFieldLabelToIndex( tIndex ).end();
    }

    //------------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_entity_rank_field_index( mtk::EntityRank aEntityRank )
    {
        MORIS_ERROR( aEntityRank == mtk::EntityRank::NODE || aEntityRank == mtk::EntityRank::ELEMENT,
                "Only node and cell fields are supported" );

        moris_index tIndex = MORIS_INDEX_MAX;

        if ( aEntityRank == mtk::EntityRank::NODE )
        {
            tIndex = 0;
        }

        else if ( aEntityRank == mtk::EntityRank::ELEMENT )
        {
            tIndex = 1;
        }

        return tIndex;
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::write_mesh( moris::Parameter_List* aParamList )
    {
        this->create_set_names();

        // get path to output XTK files to
        std::string tOutputPath = aParamList->get< std::string >( "output_path" );
        std::string tOutputFile = aParamList->get< std::string >( "output_file" );
        std::string tOutputBase = tOutputFile.substr( 0, tOutputFile.find( "." ) ) + "_ip";
        std::string tOutputExt  = tOutputFile.substr( tOutputFile.find( "." ), tOutputFile.length() );

        MORIS_ASSERT( tOutputExt == ".exo" || tOutputExt == ".e", "Invalid file extension, needs to be .exo or .e" );

        // Write mesh
        moris::mtk::Writer_Exodus writer( this );

        // if user requests to keep XTK output for all iterations, add iteration count to output file name
        if ( aParamList->get< bool >( "keep_all_opt_iters" ) )
        {
            // get optimization iteration ( function returns zero if no optimization )
            uint tOptIter = gLogger.get_opt_iteration();

            writer.write_mesh(
                    "", tOutputPath + tOutputBase + "." + std::to_string( tOptIter ) + tOutputExt, "", tOutputPath + "xtk_temp2." + std::to_string( tOptIter ) + tOutputExt );
        }
        // else: proceed as usual and overwrite xtk_temp.exo each iteration
        else
        {
            writer.write_mesh(
                    "", tOutputPath + tOutputBase + tOutputExt, "", tOutputPath + "xtk_temp3.exo" );
        }

        if ( aParamList->get< bool >( "write_enrichment_fields" ) )
        {
            std::string tProbeSpheresStr = aParamList->get< std::string >( "write_enrichment_fields_probe_spheres" );

            Vector< std::string > tNodeFields;

            if ( !tProbeSpheresStr.empty() )
            {
                Matrix< DDRMat > tProbeSpheres = string_to_mat< DDRMat >( tProbeSpheresStr );

                // set up the nodal fields for basis support
                this->create_basis_support_fields( tProbeSpheres );
            }
        }

        if ( aParamList->get< bool >( "write_basis_functions" ) )
        {
            std::string           tProbeSpheresStr = aParamList->get< std::string >( "write_enrichment_fields_probe_spheres" );
            Vector< std::string > tNodeFields;

            if ( !tProbeSpheresStr.empty() )
            {
                Matrix< DDRMat > tProbeSpheres = string_to_mat< DDRMat >( tProbeSpheresStr );

                // set up the nodal fields for basis support
                this->create_basis_function_fields( tProbeSpheres );
            }
        }

        Vector< std::string > tNodeFields = this->get_field_names( mtk::EntityRank::NODE );
        writer.set_nodal_fields( tNodeFields );

        for ( uint iF = 0; iF < tNodeFields.size(); iF++ )
        {
            moris::moris_index tFieldIndex = this->get_field_index( tNodeFields( iF ), mtk::EntityRank::NODE );
            writer.write_nodal_field( tNodeFields( iF ), this->get_field_data( tFieldIndex, mtk::EntityRank::NODE ) );
        }

        // create element id field
        this->create_cell_id_fields();

        // iterate through blocks
        Vector< std::string > tCellFields = this->get_field_names( mtk::EntityRank::ELEMENT );

        writer.set_elemental_fields( tCellFields );

        Vector< std::string > tBlockNames = this->get_set_names( mtk::EntityRank::ELEMENT );

        for ( uint iField = 0; iField < tCellFields.size(); iField++ )
        {

            moris::moris_index       tFieldIndex = this->get_field_index( tCellFields( iField ), mtk::EntityRank::ELEMENT );
            Matrix< DDRMat > const & tFieldData  = this->get_field_data( tFieldIndex, mtk::EntityRank::ELEMENT );

            for ( uint iBlock = 0; iBlock < this->get_num_blocks(); iBlock++ )
            {
                std::string tBlockName  = tBlockNames( iBlock );
                moris_index tBlockIndex = this->get_block_set_index( tBlockName );

                Matrix< IndexMat > tCellIndices = this->get_element_indices_in_block_set( tBlockIndex );

                Matrix< DDRMat > tBlockFieldData( 1, tCellIndices.numel(), -10.0 );

                for ( uint iCell = 0; iCell < tCellIndices.numel(); iCell++ )
                {
                    tBlockFieldData( iCell ) = tFieldData( tCellIndices( iCell ) );
                }

                if ( tBlockFieldData.numel() > 0 )
                {

                    writer.write_elemental_field( tBlockName, tCellFields( iField ), tBlockFieldData );
                }
            }
        }

        // Write the fields
        writer.set_time( 0.0 );
        writer.close_file();
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::create_cell_id_fields()
    {
        // Fields constructed here
        Vector< std::string > tCellFields = { "cell_id" };

        moris_index tFieldIndex = this->create_field( tCellFields( 0 ), mtk::EntityRank::ELEMENT, 0 );

        Matrix< DDRMat > tCellIdField( 1, this->get_num_elems() );

        for ( uint iCell = 0; iCell < this->get_num_elems(); iCell++ )
        {
            tCellIdField( iCell ) = (moris::real)this->get_mtk_cell( iCell ).get_id();
        }

        this->add_field_data( tFieldIndex, mtk::EntityRank::ELEMENT, tCellIdField );
    }

    //------------------------------------------------------------------------------

    moris_index
    Enriched_Interpolation_Mesh::get_block_set_index( std::string aBlockSetLabel ) const
    {
        // find the iterator
        auto it = std::find_if( mBlockSetNames.begin(), mBlockSetNames.end(), [ & ]( const std::string& aStr ) {
            return aStr == aBlockSetLabel;
        } );

        // find the index
        moris_index tIndex = std::distance( mBlockSetNames.begin(), it );

        return tIndex;
    }

    //------------------------------------------------------------------------------

    moris::uint
    Enriched_Interpolation_Mesh::get_num_blocks() const
    {
        return mBlockSetNames.size();
    }

    //------------------------------------------------------------------------------

    Vector< mtk::Cell const * >
    Enriched_Interpolation_Mesh::get_set_cells( std::string aSetLabel ) const
    {
        // get the index of the set
        moris_index tSetIndex = this->get_block_set_index( aSetLabel );

        // initailize the output
        Vector< mtk::Cell const * > tSetCells;
        tSetCells.reserve( ( mElementIndicesInBlock( tSetIndex ).size() ) );

        // fill out the cell based on indices
        for ( const auto& iCellIndex : mElementIndicesInBlock( tSetIndex ) )
        {
            mtk::Cell const & tIPCell = this->get_mtk_cell( iCellIndex );

            tSetCells.push_back( &tIPCell );
        }

        return tSetCells;
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Interpolation_Mesh::update_communication_table( Vector< moris_id > const & aNewCommunicationTable )
    {
        mXTKModel->get_cut_integration_mesh()->update_communication_table( aNewCommunicationTable );
    }
}    // namespace moris::xtk

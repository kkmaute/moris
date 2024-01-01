/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Background_Mesh.cpp
 *
 */

#include "cl_XTK_Background_Mesh.hpp"
#include "fn_unique.hpp"
#include "cl_MTK_Enums.hpp"

namespace xtk
{
    // ----------------------------------------------------------------------------------
    // Constructor/Deconstructor Source code
    // ----------------------------------------------------------------------------------

    Background_Mesh::Background_Mesh(
            moris::mtk::Interpolation_Mesh* aMeshData )
            : mMeshData( aMeshData )
            , mEntityLocaltoGlobalMap( 4 )
            , mChildMtkCells( 0 )
            , mXtkMtkVertices( 0 )
            , mNodeIndexToChildMeshIndex( 0, 0 )

    {
        // initialize_downward_inheritance();
        // mExternalMeshData.set_up_external_entity_data(mMeshData);
        // initialize_background_mesh_vertices();
        // setup_local_to_global_maps();
    }

    // ----------------------------------------------------------------------------------

    Background_Mesh::Background_Mesh(
            moris::mtk::Interpolation_Mesh* aMeshData,
            moris::ge::Geometry_Engine*     aGeometryEngine )
            : mMeshData( aMeshData )
            , mEntityLocaltoGlobalMap( 4 )
            , mChildMtkCells( 0 )
            , mXtkMtkVertices( 0 )
            , mNodeIndexToChildMeshIndex( 0, 0 )
    {
        // initialize_downward_inheritance();
        // mExternalMeshData.set_up_external_entity_data(mMeshData);
        // initialize_background_mesh_vertices();
        // setup_local_to_global_maps();
    }

    // ----------------------------------------------------------------------------------

    Background_Mesh::~Background_Mesh()
    {
        for ( auto i : mChildMtkCells )
        {
            delete i;
        }
        mChildMtkCells.clear();
    }

    // ----------------------------------------------------------------------------------
    // Entity related functions
    // ----------------------------------------------------------------------------------

    moris::size_t
    Background_Mesh::get_num_entities( mtk::EntityRank aEntityRank ) const
    {
        // Initialize
        if ( aEntityRank == mtk::EntityRank::NODE )
        {
            moris::size_t tNumBackgroundEntities = mMeshData->get_num_entities( (mtk::EntityRank)aEntityRank );
            moris::size_t tExternalEntities      = mExternalMeshData.get_num_entities_external_data( aEntityRank );
            return tNumBackgroundEntities + tExternalEntities;
        }
        else if ( aEntityRank == mtk::EntityRank::ELEMENT )
        {
            moris::size_t tNumBackgroundEntities = mMeshData->get_num_entities( (mtk::EntityRank)aEntityRank );
            return tNumBackgroundEntities + mChildMtkCells.size();
        }
        else
        {
            MORIS_ERROR( 0, "Only cells and verts supported." );
            return 0;
        }
    }

    // ----------------------------------------------------------------------------------

    moris::size_t
    Background_Mesh::get_num_entities_background( mtk::EntityRank aEntityRank ) const
    {
        return mMeshData->get_num_entities( (mtk::EntityRank)aEntityRank );
    }

    // ----------------------------------------------------------------------------------

    moris::moris_index
    Background_Mesh::get_vertex_owner( moris::moris_index aVertexIndex ) const

    {
        if ( mExternalMeshData.is_external_entity( aVertexIndex, mtk::EntityRank::NODE ) )
        {
            return mExternalMeshData.get_external_vertex_owner( aVertexIndex );
        }
        else
        {
            return mMeshData->get_entity_owner( aVertexIndex, mtk::EntityRank::NODE );
        }
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Background_Mesh::get_vertices_owner( moris::Matrix< moris::IndexMat > const & aVertexIndices ) const
    {
        moris::uint tNumVerts = aVertexIndices.numel();

        moris::Matrix< moris::IndexMat > tVertexOwners( 1, tNumVerts );

        for ( moris::uint i = 0; i < tNumVerts; i++ )
        {
            tVertexOwners( i ) = this->get_vertex_owner( aVertexIndices( i ) );
        }

        return tVertexOwners;
    }

    // ----------------------------------------------------------------------------------

    Vector< moris::mtk::Vertex const * >
    Background_Mesh::get_mtk_vertices( Matrix< IndexMat > const & aVertexIndices )
    {
        Vector< moris::mtk::Vertex const * > tVertices( aVertexIndices.numel() );

        for ( moris::uint i = 0; i < aVertexIndices.numel(); i++ )
        {
            tVertices( i ) = &this->get_mtk_vertex( aVertexIndices( i ) );
        }
        return tVertices;
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Vertex&
    Background_Mesh::get_mtk_vertex( moris::moris_index aVertexIndex )
    {
        MORIS_ASSERT( aVertexIndex < (moris::moris_index)mXtkMtkVertices.size(),
                "Vertex index is out of bounds" );

        return mXtkMtkVertices( aVertexIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Vertex_XTK&
    Background_Mesh::get_mtk_vertex_xtk( moris::moris_index aVertexIndex )
    {
        MORIS_ASSERT( aVertexIndex < (moris::moris_index)mXtkMtkVertices.size(),
                "Vertex index is out of bounds" );

        return mXtkMtkVertices( aVertexIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Cell&
    Background_Mesh::get_mtk_cell( moris::moris_index aCellIndex )
    {

        if ( !this->is_background_cell( aCellIndex ) )
        {
            return get_child_element_mtk_cell( aCellIndex );
        }
        // otherwise the cell comes from the background mesh
        else
        {
            return mMeshData->get_mtk_cell( aCellIndex );
        }
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Cell const &
    Background_Mesh::get_mtk_cell( moris::moris_index aCellIndex ) const
    {
        if ( !this->is_background_cell( aCellIndex ) )
        {
            return get_child_element_mtk_cell( aCellIndex );
        }
        // otherwise the cell comes from the background mesh
        else
        {
            return mMeshData->get_mtk_cell( aCellIndex );
        }
    }

    // ----------------------------------------------------------------------------------

    moris::moris_id
    Background_Mesh::allocate_entity_ids(
            moris::size_t   aNumReqs,
            mtk::EntityRank aChildEntityRank )
    {
        return mExternalMeshData.allocate_entity_ids_external_entity_data( aNumReqs, aChildEntityRank );
    }

    // ----------------------------------------------------------------------------------

    moris::moris_index
    Background_Mesh::get_first_available_index( mtk::EntityRank aEntityRank ) const
    {
        return mExternalMeshData.get_first_available_index_external_data( aEntityRank );
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::update_first_available_index(
            moris::size_t   aNewFirstAvailableIndex,
            mtk::EntityRank aEntityRank )
    {
        mExternalMeshData.update_first_available_index_external_data( aNewFirstAvailableIndex, aEntityRank );
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::batch_create_new_nodes(
            Vector< moris_index > const &                    aNewNodeIds,
            Vector< moris_index > const &                    aNewNodeIndices,
            Vector< moris_index > const &                    aNewNodeOwningProc,
            Vector< moris::Matrix< moris::DDRMat > > const & aNewNodeCoordinates )
    {
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::batch_create_new_nodes_as_copy_of_other_nodes(
            moris::Matrix< moris::IndexMat > const & aExistingNodeIndices,
            moris::Matrix< moris::IndexMat > const & aNewNodeIds,
            moris::Matrix< moris::IndexMat > const & aNewNodeIndices )
    {
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::allocate_external_node_to_child_mesh_associations()
    {
        // Hard coded to 4 which is the case when a node is created on a parent edge
        // This should always be the max for XTK.
        const moris::size_t tNumCMPerNode = 4;

        // Number of external nodes
        moris::size_t tExtNumNodes = mExternalMeshData.get_num_entities_external_data( mtk::EntityRank::NODE );

        // Allocate matrix filled with moris_index max
        mNodeIndexToChildMeshIndex.resize( tExtNumNodes, tNumCMPerNode );
        mNodeIndexToChildMeshIndex.fill( std::numeric_limits< moris::moris_index >::max() );
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::associate_external_nodes_to_child_mesh(
            moris::moris_index                       aChildMeshIndex,
            moris::Matrix< moris::IndexMat > const & aNodeIndices )
    {
        // Make sure the node to child mesh matrix has been allocated
        MORIS_ASSERT( mNodeIndexToChildMeshIndex.n_rows() == mExternalMeshData.get_num_entities_external_data( mtk::EntityRank::NODE ),
                "mNodeIndexToChildMeshIndex has not been allocated" );

        MORIS_ASSERT( moris::isvector( aNodeIndices ), "Provided node indices need to be a vector" );

        // Number of columns in mNodeIndexToChildMeshIndex
        size_t tNumCols = mNodeIndexToChildMeshIndex.n_cols();

        // Iterate over nodes and create associated between node index and child mesh index
        for ( size_t i = 0; i < aNodeIndices.numel(); i++ )
        {
            if ( mExternalMeshData.is_external_entity( aNodeIndices( i ), mtk::EntityRank::NODE ) )
            {
                moris::size_t tExtIndex = mExternalMeshData.get_external_entity_index( aNodeIndices( i ), mtk::EntityRank::NODE );
                for ( size_t j = 0; j < tNumCols; j++ )
                {
                    if ( mNodeIndexToChildMeshIndex( tExtIndex, j ) == std::numeric_limits< moris::moris_index >::max() )
                    {
                        mNodeIndexToChildMeshIndex( tExtIndex, j ) = aChildMeshIndex;
                        break;
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    associate_external_nodes_to_child_mesh(
            moris::moris_index                       aChildMeshIndex,
            moris::Matrix< moris::IndexMat > const & aNodeIndices );

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Background_Mesh::get_node_child_mesh_assocation( moris::moris_index aNodeIndex ) const
    {
        MORIS_ASSERT( mExternalMeshData.is_external_entity( aNodeIndex, mtk::EntityRank::NODE ),
                "Provided node index needs to be one created during the decomposition process" );

        // External index
        moris::size_t tExtIndex = mExternalMeshData.get_external_entity_index( aNodeIndex, mtk::EntityRank::NODE );

        moris::Matrix< moris::IndexMat > tChildMeshIndices = mNodeIndexToChildMeshIndex.get_row( tExtIndex );

        size_t tCount = 0;
        for ( size_t i = 0; i < tChildMeshIndices.n_cols(); i++ )
        {
            if ( tChildMeshIndices( i ) == std::numeric_limits< moris::moris_index >::max() )
            {
                break;
            }

            tCount++;
        }

        tChildMeshIndices.resize( 1, tCount );
        return tChildMeshIndices;
    }

    // ----------------------------------------------------------------------------------

    moris::size_t
    Background_Mesh::get_glb_entity_id_from_entity_loc_index(
            moris::size_t       aEntityIndex,
            mtk::EntityRank     aEntityRank,
            moris_index const & aMeshIndex ) const
    {

        if ( aEntityRank == mtk::EntityRank::ELEMENT || aEntityRank == mtk::EntityRank::NODE )
        {
            MORIS_ASSERT( aEntityIndex < mEntityLocaltoGlobalMap( (uint)aEntityRank ).size(),
                    "Entity Index out of bounds" );

            return mEntityLocaltoGlobalMap( (uint)aEntityRank )( aEntityIndex );
        }
        else
        {
            return mMeshData->get_glb_entity_id_from_entity_loc_index( aEntityIndex, aEntityRank, aMeshIndex );
        }
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IdMat >
    Background_Mesh::get_glb_entity_id_from_entity_loc_index_range(
            moris::Matrix< moris::IndexMat > const & tEntityIndices,
            mtk::EntityRank                          aEntityRank ) const
    {
        MORIS_ASSERT( moris::isvector( tEntityIndices ), "Entity indices are not provided in a vector" );
        moris::size_t                 tNumEntities = tEntityIndices.numel();
        moris::Matrix< moris::IdMat > tEntityIds( 1, tNumEntities );

        for ( moris::size_t i = 0; i < tNumEntities; i++ )
        {
            tEntityIds( 0, i ) = this->get_glb_entity_id_from_entity_loc_index( tEntityIndices( i ), aEntityRank );
        }
        return tEntityIds;
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::convert_loc_entity_ind_to_glb_entity_ids(
            mtk::EntityRank                   aEntityRank,
            moris::Matrix< moris::IndexMat >& aEntityIndices ) const
    {
        moris::uint tNumRows = aEntityIndices.n_rows();
        moris::uint tNumCols = aEntityIndices.n_cols();
        for ( moris::uint j = 0; j < tNumCols; j++ )
        {
            for ( moris::uint i = 0; i < tNumRows; i++ )
            {
                aEntityIndices( i, j ) =
                        this->get_glb_entity_id_from_entity_loc_index( aEntityIndices( i, j ), aEntityRank );
            }
        }
    }

    // ----------------------------------------------------------------------------------

    bool
    Background_Mesh::is_background_node( moris::moris_index aNodeIndex )
    {
        return !mExternalMeshData.is_external_entity( aNodeIndex, mtk::EntityRank::NODE );
    }

    // ----------------------------------------------------------------------------------

    bool
    Background_Mesh::is_background_cell( moris::moris_index aCellIndex ) const
    {
        if ( aCellIndex < (moris_index)mMeshData->get_num_elems() )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Background_Mesh::get_all_node_coordinates_loc_inds() const
    {
        // Get counts inside of STK and in External Data (both happen in this function call)
        moris::size_t tNumNodes   = this->get_num_entities( mtk::EntityRank::NODE );
        moris::size_t tNumBGNodes = mMeshData->get_num_entities( (mtk::EntityRank) mtk::EntityRank::NODE );

        moris::Matrix< moris::DDRMat > tAllNodeCoordinates = this->get_all_node_coordinates_loc_inds_background_mesh();

        tAllNodeCoordinates.resize( tNumNodes, mMeshData->get_spatial_dim() );

        // Get node coordinates from external entities
        mExternalMeshData.get_all_node_coordinates_loc_inds_external_data( tNumBGNodes, tAllNodeCoordinates );

        return tAllNodeCoordinates;
    }

    // ---------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Background_Mesh::get_selected_node_coordinates_loc_inds(
            moris::Matrix< moris::IndexMat > const & aNodeIndices ) const
    {
        MORIS_ERROR( moris::isvector( aNodeIndices ), "Provided Node indices need to be a vector" );
        // TODO: Add external entity check to see if xtk has the coordinate field or stk has it
        // Number of spatial dimensions
        moris::size_t tSpatialDimension = mMeshData->get_spatial_dim();

        // Get number of nodes provided
        moris::size_t tNumNodes = aNodeIndices.numel();

        // Initialize output matrix
        moris::Matrix< moris::DDRMat > tSelectedNodesCoords( tNumNodes, tSpatialDimension );

        // Loop over all nodes
        mtk::EntityRank tEntityRank = mtk::EntityRank::NODE;

        for ( moris::size_t n = 0; n < tNumNodes; ++n )
        {
            if ( mExternalMeshData.is_external_entity( aNodeIndices( n ), tEntityRank ) )
            {
                moris::Matrix< moris::DDRMat > const & tNodeCoords =
                        mExternalMeshData.get_selected_node_coordinates_loc_inds_external_data( aNodeIndices( n ) );

                tSelectedNodesCoords.set_row( n, tNodeCoords );
            }
            else
            {
                Matrix< DDRMat > tNodeCoord = mMeshData->get_node_coordinate( (moris_index)aNodeIndices( n ) );

                for ( moris::uint iDim = 0; iDim < tSpatialDimension; iDim++ )
                {
                    tSelectedNodesCoords( n, iDim ) = tNodeCoord( iDim );
                }
            }
        }

        return tSelectedNodesCoords;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IdMat >
    Background_Mesh::get_local_to_global_map( mtk::EntityRank aEntityRank ) const
    {
        MORIS_ERROR( aEntityRank == mtk::EntityRank::NODE, " This function is only implemented for node maps" );

        moris::uint tNumNodes = mMeshData->get_num_entities( mtk::EntityRank::NODE );

        moris::Matrix< moris::IndexMat > tLocalToGlobalBM( 1, tNumNodes );

        for ( size_t i = 0; i < tNumNodes; i++ )
        {
            tLocalToGlobalBM( i ) = mMeshData->get_glb_entity_id_from_entity_loc_index( i, mtk::EntityRank::NODE );
        }

        moris::Matrix< moris::IndexMat > const & tLocalToGlobalExt = mExternalMeshData.get_local_to_global_node_map();

        tNumNodes               = tLocalToGlobalBM.numel() + tLocalToGlobalExt.numel();
        size_t tFirstExtNodeInd = tLocalToGlobalBM.numel();

        MORIS_ERROR( tNumNodes = tLocalToGlobalBM.numel() + tLocalToGlobalExt.numel(),
                "Number of nodes returned does not match the number in the map" );

        // combine the two maps
        moris::Matrix< moris::IdMat > tLocalToGlobal( 1, tNumNodes );

        for ( moris::size_t i = 0; i < tLocalToGlobalBM.numel(); i++ )
        {
            tLocalToGlobal( i ) = tLocalToGlobalBM( i );
        }

        for ( moris::size_t i = 0; i < tLocalToGlobalExt.numel(); i++ )
        {
            tLocalToGlobal( i + tFirstExtNodeInd ) = tLocalToGlobalExt( i );
        }

        return tLocalToGlobal;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IdMat >
    Background_Mesh::get_full_non_intersected_node_to_element_glob_ids() const
    {
        moris::size_t     tNumElementsBG = this->get_num_entities( mtk::EntityRank::ELEMENT );
        mtk::CellTopology tElemTopo      = this->get_parent_cell_topology();

        moris::size_t tNumNodesPerElem = 0;
        if ( tElemTopo == mtk::CellTopology::TET4 )
        {
            tNumNodesPerElem = 4;
        }
        else if ( tElemTopo == mtk::CellTopology::HEX8 )
        {
            tNumNodesPerElem = 8;
        }
        else
        {
            tNumNodesPerElem = 8;
        }

        moris::Matrix< moris::IdMat > tElementToNode( tNumElementsBG, tNumNodesPerElem );
        moris::size_t                 tCount = 0;
        for ( moris::size_t i = 0; i < tElementToNode.n_rows(); i++ )
        {
            if ( !this->entity_has_children( i, mtk::EntityRank::ELEMENT ) )
            {
                moris::Matrix< moris::IndexMat > tElementToNodeInd = mMeshData->get_entity_connected_to_entity_loc_inds( (moris::moris_index)i,
                        mtk::EntityRank::ELEMENT,
                        mtk::EntityRank::NODE );

                moris::Matrix< moris::IdMat > tElementToNodeId =
                        moris::mtk::convert_entity_indices_to_ids( tElementToNodeInd, mtk::EntityRank::NODE, mMeshData );

                tElementToNode.set_row( tCount, tElementToNodeId );
                tCount++;
            }
        }

        tElementToNode.resize( tCount, tNumNodesPerElem );
        return tElementToNode;
    }

    // ----------------------------------------------------------------------------------

    Vector< moris::Matrix< moris::IdMat > >
    Background_Mesh::get_non_intersected_element_to_node_by_phase( moris::uint aNumPhases )
    {
        uint              tNumElementsBG = mMeshData->get_num_entities( mtk::EntityRank::ELEMENT );
        mtk::CellTopology tElemTopo      = this->get_parent_cell_topology();

        moris::size_t tNumNodesPerElem = 0;
        if ( tElemTopo == mtk::CellTopology::TET4 )
        {
            tNumNodesPerElem = 4;
        }
        else if ( tElemTopo == mtk::CellTopology::HEX8 )
        {
            tNumNodesPerElem = 8;
        }
        else if ( tElemTopo == mtk::CellTopology::QUAD4 )
        {
            tNumNodesPerElem = 4;
        }
        else
        {
            tNumNodesPerElem = 8;
        }

        Vector< moris::Matrix< moris::IdMat > > tElementToNodeByPhase(
                aNumPhases,
                moris::Matrix< moris::IdMat >( tNumElementsBG, tNumNodesPerElem ) );

        moris::Matrix< moris::DDUMat > tPhaseCount( 1, aNumPhases, 0 );

        for ( moris::size_t i = 0; i < tNumElementsBG; i++ )
        {
            if ( !this->entity_has_children( i, mtk::EntityRank::ELEMENT ) )
            {
                moris::moris_id tId = mMeshData->get_glb_entity_id_from_entity_loc_index( i, mtk::EntityRank::ELEMENT );

                moris::Matrix< moris::IdMat > tElementToNodeId =
                        mMeshData->get_entity_connected_to_entity_glob_ids( tId, mtk::EntityRank::ELEMENT, mtk::EntityRank::NODE );

                moris::moris_index tPhaseIndex = this->get_element_phase_index( i );

                for ( moris::uint j = 0; j < tNumNodesPerElem; j++ )
                {
                    tElementToNodeByPhase( tPhaseIndex )( tPhaseCount( tPhaseIndex ), j ) = tElementToNodeId( j );
                }

                tPhaseCount( tPhaseIndex )++;
            }
        }

        for ( moris::uint i = 0; i < aNumPhases; i++ )
        {
            tElementToNodeByPhase( i ).resize( tPhaseCount( i ), tNumNodesPerElem );
        }

        return tElementToNodeByPhase;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IdMat >
    Background_Mesh::get_all_non_intersected_elements() const
    {
        moris::size_t tNumElementsBG = this->get_num_entities( mtk::EntityRank::ELEMENT );

        moris::Matrix< moris::IdMat > tElementIds( tNumElementsBG, 1 );
        moris::size_t                 tCount = 0;
        for ( moris::size_t i = 0; i < tElementIds.numel(); i++ )
        {
            if ( !this->entity_has_children( i, mtk::EntityRank::ELEMENT ) )
            {
                moris::moris_id tElementToNodeId =
                        mMeshData->get_glb_entity_id_from_entity_loc_index( (moris::moris_index)i, mtk::EntityRank::ELEMENT );

                tElementIds( tCount ) = tElementToNodeId;
                tCount++;
            }
        }

        tElementIds.resize( tCount, 1 );
        return tElementIds;
    }

    // ----------------------------------------------------------------------------------

    Vector< moris::Matrix< moris::IdMat > >
    Background_Mesh::get_all_non_intersected_elements_by_phase( uint aNumPhases ) const
    {
        uint tNumElementsBG = mMeshData->get_num_entities( mtk::EntityRank::ELEMENT );

        // Initialize output
        Vector< moris::Matrix< moris::IdMat > > tElementsByPhase( aNumPhases );
        moris::Matrix< moris::DDUMat >        tPhaseCount( 1, aNumPhases, 0 );
        for ( uint i = 0; i < aNumPhases; i++ )
        {
            tElementsByPhase( i ) = moris::Matrix< moris::IdMat >( 1, tNumElementsBG );
        }

        for ( uint iE = 0; iE < tNumElementsBG; iE++ )
        {
            if ( !this->entity_has_children( iE, mtk::EntityRank::ELEMENT ) )
            {
                moris::moris_index tPhaseInd = this->get_element_phase_index( iE );
                uint               tCount    = tPhaseCount( tPhaseInd );

                tElementsByPhase( tPhaseInd )( tCount ) = mMeshData->get_glb_entity_id_from_entity_loc_index(
                        (moris::moris_index)iE,
                        mtk::EntityRank::ELEMENT );

                tPhaseCount( tPhaseInd )++;
            }
        }

        // Resize
        for ( uint i = 0; i < aNumPhases; i++ )
        {
            tElementsByPhase( i ).resize( 1, tPhaseCount( i ) );
        }

        return tElementsByPhase;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Background_Mesh::get_all_non_intersected_elements_loc_inds() const
    {
        moris::size_t tNumElementsBG = mMeshData->get_num_entities( mtk::EntityRank::ELEMENT );

        moris::Matrix< moris::IdMat > tElementInds( tNumElementsBG, 1 );
        moris::size_t                 tCount = 0;
        for ( moris::size_t i = 0; i < tElementInds.numel(); i++ )
        {
            if ( !this->entity_has_children( i, mtk::EntityRank::ELEMENT ) )
            {
                tElementInds( tCount ) = i;
                tCount++;
            }
        }
        tElementInds.resize( tCount, 1 );
        return tElementInds;
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::register_new_downward_inheritance(
            Vector< std::pair< moris::moris_index, moris::moris_index > > const & aNewElementToChildMeshPairs )
    {
        moris::size_t tNumNewPairs = aNewElementToChildMeshPairs.size();
        for ( moris::size_t i = 0; i < tNumNewPairs; i++ )
        {
            mElementDownwardInheritance.register_new_inheritance_pair(
                    aNewElementToChildMeshPairs( i ).first,
                    aNewElementToChildMeshPairs( i ).second );
        }
    }

    void
    Background_Mesh::setup_downward_inheritance( Cut_Mesh& aCutMesh )
    {
        // reset the downward inheritance
        mElementDownwardInheritance = Downward_Inheritance< moris::moris_index, moris::moris_index >();

        for ( moris::size_t i = 0; i < aCutMesh.get_num_child_meshes(); i++ )
        {
            moris::uint tParentElementIndex = aCutMesh.get_parent_element_index( i );

            mElementDownwardInheritance.register_new_inheritance_pair(
                    (moris::moris_index)tParentElementIndex,
                    (moris::moris_index)i );
        }
    }

    // ----------------------------------------------------------------------------------

    bool
    Background_Mesh::entity_has_children(
            moris::size_t   aEntityIndex,
            mtk::EntityRank aEntityRank ) const
    {
        MORIS_ASSERT( aEntityRank == mtk::EntityRank::ELEMENT, "ONLY ELEMENT DOWNWARD INHERITANCE SUPPORTED" );

        return mElementDownwardInheritance.has_inheritance( aEntityIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::moris_index const &
    Background_Mesh::child_mesh_index(
            moris::size_t   aEntityIndex,
            mtk::EntityRank aEntityRank )
    {
        MORIS_ASSERT( aEntityRank == mtk::EntityRank::ELEMENT, "ONLY ELEMENT DOWNWARD INHERITANCE SUPPORTED" );

        return mElementDownwardInheritance.get_inheritance( aEntityIndex );
    }

    // ---------------------------------------------------------------------------------

    void
    Background_Mesh::initialize_interface_node_flags(
            moris::size_t const & aNumNodes,
            moris::size_t const & aNumGeometry )
    {
        mInterfaceNodeFlag = moris::Matrix< moris::IndexMat >( aNumNodes, aNumGeometry, 0 );
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::allocate_space_in_interface_node_flags(
            moris::size_t aNumNodes,
            moris::size_t aNumGeometry )
    {
        moris::size_t tCurrentSize = mInterfaceNodeFlag.n_rows();
        mInterfaceNodeFlag.resize( tCurrentSize + aNumNodes, aNumGeometry );

        for ( moris::size_t i = tCurrentSize; i < tCurrentSize + aNumNodes; i++ )
        {
            for ( moris::size_t j = 0; j < aNumGeometry; j++ )
            {
                mInterfaceNodeFlag( i, j ) = 0;
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::mark_node_as_interface_node( moris::moris_index aNodeIndex,
            moris::size_t                                            aGeomIndex )
    {
        mInterfaceNodeFlag( aNodeIndex, aGeomIndex ) = 1;
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::mark_nodes_as_interface_node_loc_inds( moris::Matrix< moris::IndexMat > aNodeIndices,
            moris::size_t                                                                    aGeomIndex )
    {
        for ( uint i = 0; i < aNodeIndices.numel(); i++ )
        {
            mInterfaceNodeFlag( aNodeIndices( i ), aGeomIndex ) = 1;
        }
    }

    // ----------------------------------------------------------------------------------

    bool
    Background_Mesh::is_interface_node(
            moris::moris_index aNodeIndex,
            moris::size_t      aGeomIndex ) const
    {
        MORIS_ASSERT( aNodeIndex < (moris::moris_index)mInterfaceNodeFlag.n_rows(),
                "Attempting to access interface node flag for node index out of bounds. Have you called allocate_space_in_interface_node_flags?" );

        if ( mInterfaceNodeFlag( aNodeIndex, aGeomIndex ) == 1 )
        {
            return true;
        }

        else
        {
            return false;
        }
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Background_Mesh::get_interface_nodes_loc_inds( moris::moris_index aGeometryIndex ) const
    {
        // initialize output
        moris::size_t                    tNumNodes = this->get_num_entities( mtk::EntityRank::NODE );
        moris::Matrix< moris::IndexMat > tInterfaceNodes( 1, tNumNodes );

        // keep track of how many interface nodes
        moris::size_t tCount = 0;

        for ( moris::size_t i = 0; i < tNumNodes; i++ )
        {
            if ( is_interface_node( i, aGeometryIndex ) )
            {
                tInterfaceNodes( 0, tCount ) = i;
                tCount++;
            }
        }

        tInterfaceNodes.resize( 1, tCount );
        return tInterfaceNodes;
    }

    // ----------------------------------------------------------------------------------

    Vector< moris::Matrix< moris::IdMat > >
    Background_Mesh::get_interface_nodes_loc_inds() const
    {
        // initialize output
        moris::size_t tNumGeoms = mInterfaceNodeFlag.n_cols();

        Vector< moris::Matrix< moris::IdMat > > tInterfaceNodes( tNumGeoms );

        for ( moris::size_t i = 0; i < tNumGeoms; i++ )
        {
            tInterfaceNodes( i ) = get_interface_nodes_loc_inds( i );
        }

        return tInterfaceNodes;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IdMat >
    Background_Mesh::get_interface_nodes_glob_ids( moris::moris_index aGeometryIndex ) const
    {
        // initialize output
        moris::size_t                 tNumNodes = this->get_num_entities( mtk::EntityRank::NODE );
        moris::Matrix< moris::IdMat > tInterfaceNodes( 1, tNumNodes );

        // keep track of how many interface nodes
        moris::size_t tCount = 0;

        for ( moris::size_t i = 0; i < tNumNodes; i++ )
        {
            if ( is_interface_node( i, aGeometryIndex ) )
            {
                tInterfaceNodes( 0, tCount ) = get_glb_entity_id_from_entity_loc_index( i, mtk::EntityRank::NODE );
                tCount++;
            }
        }

        tInterfaceNodes.resize( 1, tCount );
        return tInterfaceNodes;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Background_Mesh::restrict_vertex_list_to_owned_by_this_proc_loc_inds( moris::Matrix< moris::IndexMat > const & aNodeIndexList ) const
    {
        // proc rank
        moris_index tMyProcRank = par_rank();

        //  initialize list of owned nodes
        moris::Matrix< moris::IndexMat > tOwnedNodeList( aNodeIndexList.numel(), 1 );

        // keep track of owned nodes
        moris::uint tCount = 0;

        for ( moris::uint i = 0; i < aNodeIndexList.numel(); i++ )
        {
            if ( this->get_vertex_owner( i ) == tMyProcRank )
            {
                tOwnedNodeList( tCount ) = aNodeIndexList( i );
            }
        }

        // size out extra space
        tOwnedNodeList.resize( tCount, 1 );

        return tOwnedNodeList;
    }

    // ----------------------------------------------------------------------------------

    Vector< moris::Matrix< moris::IdMat > >
    Background_Mesh::get_interface_nodes_glob_ids()
    {
        // initialize output
        moris::size_t tNumGeoms = mInterfaceNodeFlag.n_cols();

        Vector< moris::Matrix< moris::IdMat > > tInterfaceNodes( tNumGeoms );

        for ( moris::size_t i = 0; i < tNumGeoms; i++ )
        {
            tInterfaceNodes( i ) = get_interface_nodes_glob_ids( i );
        }

        return tInterfaceNodes;
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::print_interface_node_flags()
    {
        for ( moris::size_t i = 0; i < mInterfaceNodeFlag.n_rows(); i++ )
        {
            std::cout << "Vertex Id: " << std::setw( 8 ) << this->get_glb_entity_id_from_entity_loc_index( i, mtk::EntityRank::NODE ) << " | ";
            for ( moris::size_t j = 0; j < mInterfaceNodeFlag.n_cols(); j++ )
            {
                std::cout << mInterfaceNodeFlag( i, j ) << " ";
            }
            std::cout << std::endl;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::print_vertex_map()
    {
        std::cout << "Background Mesh Vertex Map:" << std::endl;
        std::cout << "   Id    |   Index " << std::endl;
        std::cout << "-------------------" << std::endl;
        for ( auto it : mVertexGlbToLocalMap )
        {
            std::cout << std::setw( 8 ) << it.first << "|" << std::setw( 8 ) << it.second << std::endl;
        }
    }

    // ----------------------------------------------------------------------------------

    moris::Memory_Map
    Background_Mesh::get_memory_usage()
    {
        // memory map
        moris::Memory_Map tMemoryMap;

        tMemoryMap.mMemoryMapData[ "mMeshData (ptr)" ]         = sizeof( mMeshData );
        tMemoryMap.mMemoryMapData[ "mExternalMeshData" ]       = mExternalMeshData.capacity();
        tMemoryMap.mMemoryMapData[ "mEntityLocaltoGlobalMap" ] = moris::internal_capacity( mEntityLocaltoGlobalMap );
        tMemoryMap.mMemoryMapData[ "mCommunicationMap" ]       = mCommunicationMap.capacity();
        tMemoryMap.mMemoryMapData[ "mChildMtkCells" ]          = moris::internal_capacity_ptr( mChildMtkCells );
        // fixme: add mVertexGlbToLocalMap
        //  tMemoryMap.mMemoryMapData["mVertexGlbToLocalMap"] =moris::internal_capacity(mVertexGlbToLocalMap);
        tMemoryMap.mMemoryMapData[ "mXtkMtkVertices" ]            = mXtkMtkVertices.capacity();
        tMemoryMap.mMemoryMapData[ "mNodeIndexToChildMeshIndex" ] = mNodeIndexToChildMeshIndex.capacity();
        tMemoryMap.mMemoryMapData[ "mElementPhaseIndex" ]         = mElementPhaseIndex.capacity();
        tMemoryMap.mMemoryMapData[ "mInterfaceNodeFlag" ]         = mInterfaceNodeFlag.capacity();

        return tMemoryMap;
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::initialize_element_phase_indices( moris::size_t const & aNumElements )
    {
        mElementPhaseIndex = moris::Matrix< moris::IndexMat >( aNumElements, 1, -1 );
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::set_element_phase_index(
            moris::size_t const & aElementIndex,
            moris::size_t const & aElementPhaseIndex )
    {
        mElementPhaseIndex( aElementIndex, 0 ) = aElementPhaseIndex;
    }

    // ----------------------------------------------------------------------------------

    moris::moris_index const &
    Background_Mesh::get_element_phase_index( moris::size_t const & aElementIndex ) const
    {
        return (mElementPhaseIndex)( aElementIndex, 0 );
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Background_Mesh::get_element_phase_inds( moris::Matrix< moris::DDSTMat > const & aElementInds )
    {
        moris::Matrix< moris::IndexMat > tElementPhasesInds( 1, aElementInds.n_cols() );

        for ( moris::size_t i = 0; i < aElementInds.n_cols(); i++ )
        {
            tElementPhasesInds( 0, i ) = this->get_element_phase_index( aElementInds( i ) );
        }

        return tElementPhasesInds;
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Background_Mesh::get_loc_entity_ind_from_entity_glb_id(
            moris_id        aEntityId,
            mtk::EntityRank aEntityRank ) const
    {

        MORIS_ERROR( aEntityRank == mtk::EntityRank::NODE,
                "Background_Mesh::get_loc_entity_ind_from_entity_glb_id() - Only a node map is implemented in XTK" );

        auto tIter = mVertexGlbToLocalMap.find( aEntityId );

        MORIS_ERROR( tIter != mVertexGlbToLocalMap.end(),
                "Background_Mesh::get_loc_entity_ind_from_entity_glb_id() - "
                "Provided Entity Id is not in the map, Has the map been initialized?: aEntityId =%u mtk::EntityRank = %u on process %u",
                aEntityId,
                (uint)aEntityRank,
                par_rank() );

        return tIter->second;
    }

    std::unordered_map< moris_id, moris_index >
    Background_Mesh::get_vertex_glb_id_to_loc_vertex_ind_map() const
    {
        return mVertexGlbToLocalMap;
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::add_child_element_to_mtk_cells(
            moris::moris_index aElementIndex,
            moris::moris_index aElementId,
            moris::moris_index aElementOwner,
            moris::moris_index aCMElementIndex,
            Child_Mesh*        aChildMeshPtr )
    {
        mChildMtkCells.push_back( new Cell_XTK_CM(
                aElementId,
                aElementIndex,
                aElementOwner,
                aCMElementIndex,
                aChildMeshPtr,
                this ) );

        MORIS_ASSERT( mChildMtkCellMap.find( aElementIndex ) == mChildMtkCellMap.end(),
                "Background_Mesh::add_child_element_to_mtk_cells() - Element index already has an mtk cell associated with it" );

        mChildMtkCellMap[ aElementIndex ] = mChildMtkCells.size() - 1;
    }

    //-------------------------------------------------------------------------------

    void
    Background_Mesh::add_new_cell_to_mesh( moris::mtk::Cell* aCell )
    {
        mChildMtkCells.push_back( aCell );

        MORIS_ASSERT( mChildMtkCellMap.find( aCell->get_index() ) == mChildMtkCellMap.end(),
                "Background_Mesh::add_new_cell_to_mesh() - Element index already has an mtk cell associated with it" );

        mChildMtkCellMap[ aCell->get_index() ] = mChildMtkCells.size() - 1;

        MORIS_ASSERT( (moris_index)mEntityLocaltoGlobalMap( 3 ).size() == aCell->get_index(),
                "Background_Mesh::add_new_cell_to_mesh() - Index mismatch." );

        mEntityLocaltoGlobalMap( 3 ).push_back( aCell->get_id() );
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::add_cells_to_global_to_local_map(
            Matrix< IndexMat > const & aCellIndices,
            Matrix< IdMat > const &    aCellIds )
    {
#ifdef MORIS_HAVE_DEBUG

        if ( aCellIndices.numel() > 0 )
        {
            MORIS_ASSERT( (moris_index)mEntityLocaltoGlobalMap( 3 ).size() == aCellIndices.min(),
                    "Background_Mesh::add_cells_to_global_to_local_map() - "
                    "The minimum indices calling this function needs to correspond with the size of the number of elements already in the map" );
        }
#endif
        MORIS_ASSERT( aCellIndices.numel() == aCellIds.numel(), "Dimension mismatch between indices and ids" );
        // number of cells
        moris::uint tNumCells    = aCellIndices.numel();
        moris::uint tCurrentSize = mEntityLocaltoGlobalMap( 3 ).size();

        // allocate space
        mEntityLocaltoGlobalMap( 3 ).resize( tNumCells + tCurrentSize );

        // add to map
        for ( moris::uint i = 0; i < tNumCells; i++ )
        {
            mEntityLocaltoGlobalMap( 3 )( aCellIndices( i ) ) = aCellIds( i );
        }

#ifdef MORIS_HAVE_DEBUG
        // since I can't write these functions in one line, need to have ifdef
        moris::Vector< moris::moris_index > tUniqueCellIds = moris::unique_index( mEntityLocaltoGlobalMap( 3 ) );

        MORIS_ASSERT( mEntityLocaltoGlobalMap( 3 ).size() == tUniqueCellIds.size(),
                "Background_Mesh::add_cells_to_global_to_local_map() - Duplicate cell id detected." );
#endif
    }

    // ----------------------------------------------------------------------------------

    const moris::mtk::Cell*
    Background_Mesh::get_child_element_mtk_cell_ptr( moris::moris_index aElementIndex ) const
    {
        auto tIter = mChildMtkCellMap.find( aElementIndex );

        MORIS_ASSERT( mChildMtkCellMap.find( aElementIndex ) != mChildMtkCellMap.end(), "Element index is not in the map" );

        moris::moris_index tIndex = tIter->second;
        return mChildMtkCells( tIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Cell const &
    Background_Mesh::get_child_element_mtk_cell( moris::moris_index aElementIndex ) const
    {
        auto tIter = mChildMtkCellMap.find( aElementIndex );

        MORIS_ASSERT( mChildMtkCellMap.find( aElementIndex ) != mChildMtkCellMap.end(), "Element index is not in the map" );

        moris::moris_index tIndex = tIter->second;
        return *mChildMtkCells( tIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Interpolation_Mesh&
    Background_Mesh::get_mesh_data()
    {
        return *mMeshData;
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Interpolation_Mesh const &
    Background_Mesh::get_mesh_data() const
    {
        return *mMeshData;
    }

    // ----------------------------------------------------------------------------------

    mtk::CellTopology
    Background_Mesh::get_parent_cell_topology() const
    {
        mtk::Cell&                     tCell     = mMeshData->get_mtk_cell( 0 );
        enum moris::mtk::Geometry_Type tGeomType = tCell.get_geometry_type();

        mtk::CellTopology tTopo = mtk::CellTopology::UNDEFINED;
        if ( tGeomType == moris::mtk::Geometry_Type::HEX )
        {
            tTopo = mtk::CellTopology::HEX8;
        }
        else if ( tGeomType == moris::mtk::Geometry_Type::TET )
        {
            tTopo = mtk::CellTopology::TET4;
        }
        else if ( tGeomType == moris::mtk::Geometry_Type::QUAD )
        {
            tTopo = mtk::CellTopology::QUAD4;
        }
        else
        {
            MORIS_ERROR( 0, "Provided parent cell topology not implemented." );
        }

        return tTopo;
    }

    // ----------------------------------------------------------------------------------

    Matrix< IdMat > const &
    Background_Mesh::get_communication_table() const
    {
        MORIS_ERROR( 0, "Depracated" );
        return mCommunicationMap;
    }

    //-------------------------------------------------------------------------------

    void
    Background_Mesh::add_proc_to_comm_table( moris_index aProcRank )
    {
        moris_index tIndex = mCommunicationMap.numel();

        for ( moris::uint i = 0; i < mCommunicationMap.numel(); i++ )
        {
            MORIS_ERROR( mCommunicationMap( i ) != aProcRank, "Processor rank already in communication table" );
        }

        mCommunicationMap.resize( 1, mCommunicationMap.numel() + 1 );
        mCommunicationMap( tIndex ) = aProcRank;
    }
    // ----------------------------------------------------------------------------------
    void
    Background_Mesh::remove_cells_from_mesh( Vector< moris_index > const & aCellsToRemove,
            Vector< moris_index >&                                         aOldIndexToNewCellIndex )
    {
        aOldIndexToNewCellIndex = Vector< moris_index >( this->get_num_entities( mtk::EntityRank::ELEMENT ), MORIS_INDEX_MAX );

        // fill the background cell with their original indices
        for ( moris::uint iCell = 0; iCell < mMeshData->get_num_entities( mtk::EntityRank::ELEMENT ); iCell++ )
        {
            aOldIndexToNewCellIndex( iCell ) = iCell;
        }

        // iterate through the cells and remove
        // delete the pointers
        for ( moris::uint iCell = 0; iCell < aCellsToRemove.size(); iCell++ )
        {
            MORIS_ASSERT( !this->is_background_cell( aCellsToRemove( iCell ) ), "Cannot remove background cells" );

            auto tIter = mChildMtkCellMap.find( aCellsToRemove( iCell ) );
            MORIS_ASSERT( mChildMtkCellMap.find( aCellsToRemove( iCell ) ) != mChildMtkCellMap.end(), "Element index is not in the map" );
            moris::moris_index tIndex = tIter->second;

            delete mChildMtkCells( tIndex );

            mChildMtkCells( tIndex ) = nullptr;
        }

        // clear out the map
        mChildMtkCellMap.clear();

        // resize the member data and update the map
        uint        tCount = 0;
        moris_index tIndex = mMeshData->get_num_elems();
        for ( moris::uint iCell = 0; iCell < mChildMtkCells.size(); iCell++ )
        {
            if ( mChildMtkCells( iCell ) != nullptr )
            {
                mChildMtkCells( tCount ) = mChildMtkCells( iCell );

                MORIS_ASSERT( mChildMtkCellMap.find( mChildMtkCells( iCell )->get_index() ) == mChildMtkCellMap.end(),
                        "Cell index already has an child mtk cell associated with it" );

                // update the output
                aOldIndexToNewCellIndex( mChildMtkCells( iCell )->get_index() ) = tIndex;

                // change the cell index
                mChildMtkCells( tCount )->set_index( tIndex );
                tIndex++;

                mChildMtkCellMap[ mChildMtkCells( iCell )->get_index() ] = tCount;
                tCount++;
            }
        }

        mChildMtkCells.resize( tCount );
        this->setup_local_to_global_maps();
    }
    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Background_Mesh::get_all_node_coordinates_loc_inds_background_mesh() const
    {
        size_t                         tNumNodes = mMeshData->get_num_entities( (mtk::EntityRank) mtk::EntityRank::NODE );
        moris::Matrix< moris::DDRMat > tNodeCoords( tNumNodes, mMeshData->get_spatial_dim() );
        for ( size_t i = 0; i < tNumNodes; i++ )
        {
            tNodeCoords.set_row( i, mMeshData->get_node_coordinate( i ) );
        }
        return tNodeCoords;
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::intialize_downward_inheritance()
    {
        moris::size_t tNumElements = mMeshData->get_num_entities( (mtk::EntityRank) mtk::EntityRank::ELEMENT );
        MORIS_ASSERT( tNumElements != 0, "Empty Mesh Given to XTK Mesh" );

        mElementDownwardInheritance = Downward_Inheritance< moris::moris_index, moris::moris_index >( tNumElements );
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::initialize_background_mesh_vertices()
    {
        // moris::uint tNumNodes = mMeshData->get_num_entities(EntityRank::NODE);

        // for(moris::uint  i = 0; i <tNumNodes; i++)
        // {
        //     moris::mtk::Vertex * tVertex = &mMeshData->get_mtk_vertex(i);

        //     MORIS_ASSERT(mVertexGlbToLocalMap.find(tVertex->get_id()) == mVertexGlbToLocalMap.end(),"Vertex already in map");

        //     mVertexGlbToLocalMap[tVertex->get_id()] = tVertex->get_index();
        //     mXtkMtkVertices.push_back(tVertex);
        // }
    }

    // ----------------------------------------------------------------------------------

    void
    Background_Mesh::setup_local_to_global_maps()
    {
        // this only creates node local to global maps and cell local to global maps

        // setup nodes
        uint tNumNodes               = this->get_num_entities( mtk::EntityRank::NODE );
        mEntityLocaltoGlobalMap( 0 ) = Vector< moris_index >( tNumNodes );

        for ( moris::uint iN = 0; iN < tNumNodes; iN++ )
        {
            mtk::Vertex& tVertex = this->get_mtk_vertex( (moris_index)iN );
            MORIS_ASSERT( tVertex.get_index() == (moris_index)iN, "Vertex index mismatch" );
            mEntityLocaltoGlobalMap( 0 )( tVertex.get_index() ) = tVertex.get_id();
        }

        // setup cells
        uint tNumCells               = this->get_num_entities( mtk::EntityRank::ELEMENT );
        mEntityLocaltoGlobalMap( 3 ) = Vector< moris_index >( tNumCells );
        for ( moris::uint iC = 0; iC < tNumCells; iC++ )
        {
            mtk::Cell& tCell = this->get_mtk_cell( (moris_index)iC );
            MORIS_ASSERT( tCell.get_index() == (moris_index)iC, "Cell index mismatch" );
            mEntityLocaltoGlobalMap( 3 )( tCell.get_index() ) = tCell.get_id();
        }
    }

    // ----------------------------------------------------------------------------------

}    // namespace xtk

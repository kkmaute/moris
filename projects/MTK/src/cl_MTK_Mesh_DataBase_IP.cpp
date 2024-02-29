/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_DataBase_IP.cpp
 *
 */

#include "cl_MTK_Mesh_DataBase_IP.hpp"
#include "cl_MTK_Vertex_DataBase.hpp"
#include "cl_MTK_Vertex_Interpolation_DataBase.hpp"
#include "cl_MTK_Cell_DataBase.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_Tracer.hpp"

namespace moris::mtk
{

    // uint Vertex_DataBase::mDim;
    //-----------------------------------------------------------------------------
    Interpolation_Mesh_DataBase_IP::Interpolation_Mesh_DataBase_IP( const Interpolation_Mesh& aBackgroundMesh )
            : mBackgroundMesh( aBackgroundMesh )
    {
    }

    //-----------------------------------------------------------------------------

    Interpolation_Mesh_DataBase_IP::~Interpolation_Mesh_DataBase_IP()
    {
    }

    //-----------------------------------------------------------------------------

    uint
    Interpolation_Mesh_DataBase_IP::get_spatial_dim() const
    {
        return mVertexCoordinates.n_rows();
    }

    //-----------------------------------------------------------------------------

    uint
    Interpolation_Mesh_DataBase_IP::get_num_entities(
            enum EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        switch ( aEntityRank )
        {
            case EntityRank::NODE:
            {
                return mVertices.size();
                break;
            }
            case EntityRank::ELEMENT:
            {
                return mCells.size();
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }
                return 0;
        }
    }

    //-----------------------------------------------------------------------------

    MeshType
    Interpolation_Mesh_DataBase_IP::get_mesh_type() const
    {
        return MeshType::MTK;
    }

    // ----------------------------------------------------------------------------

    const Interpolation_Mesh& Interpolation_Mesh_DataBase_IP::get_background_mesh()
    {
        return mBackgroundMesh;
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Interpolation_Mesh_DataBase_IP::get_entity_connected_to_entity_loc_inds(
            moris_index       aEntityIndex,
            enum EntityRank   aInputEntityRank,
            enum EntityRank   aOutputEntityRank,
            const moris_index aDiscretizationIndex ) const
    {
        MORIS_ERROR( 0, "get_entity_connected_to_entity_loc_inds not implemented for Interpolation_Mesh_DataBase_IP" );
        return { {} };
        // return mIPMesh.get_entity_connected_to_entity_loc_inds( aEntityIndex, aInputEntityRank, aOutputEntityRank, aDiscretizationIndex );
    }

    // ----------------------------------------------------------------------------

    Matrix< DDRMat >
    Interpolation_Mesh_DataBase_IP::get_node_coordinate( moris_index aNodeIndex ) const
    {
        return mVertexCoordinates.get_column( aNodeIndex );
    }

    // ----------------------------------------------------------------------------

    moris_id
    Interpolation_Mesh_DataBase_IP::get_glb_entity_id_from_entity_loc_index(
            moris_index aEntityIndex,
            EntityRank  aEntityRank,
            moris_index aDiscretizationIndex ) const
    {
        switch ( aEntityRank )
        {
            case ( EntityRank::NODE ):
                return mVertexIdList( aEntityIndex );
            case ( EntityRank::ELEMENT ):
                return mCellIdList( aEntityIndex );
            default:
                MORIS_ERROR( false, "Interpolation_Mesh_DataBase_IP::get_glb_entity_id_from_entity_loc_index() does not support given entity type" );
                return 0;
        }
    }

    // ----------------------------------------------------------------------------

    moris_index
    Interpolation_Mesh_DataBase_IP::get_loc_entity_ind_from_entity_glb_id(
            moris_id    aEntityId,
            EntityRank  aEntityRank,
            moris_index aDiscretizationIndex ) const
    {
        switch ( aEntityRank )
        {
            case ( EntityRank::NODE ):
            {
                auto tSearch = mVertexGlobalIdToLocalIndex.find( aEntityId );
                if ( tSearch not_eq mVertexGlobalIdToLocalIndex.end() )
                {
                    return tSearch->second;
                }
                else
                {
                    return -1;
                }
            }
            default:
                MORIS_ERROR( false, "Interpolation_Mesh_DataBase_IP::get_loc_entity_ind_from_entity_glb_id() does not support given entity type" );
                return -1;
        }
    }

    // ----------------------------------------------------------------------------

    uint
    Interpolation_Mesh_DataBase_IP::get_node_owner( moris_index aNodeIndex ) const
    {
        MORIS_ERROR( 0, "get_node_owner not implemented for Interpolation_Mesh_DataBase_IP" );
        return 0;
        // return mIPMesh.get_node_owner( aNodeIndex );
    }

    // ----------------------------------------------------------------------------

    Matrix< IdMat >
    Interpolation_Mesh_DataBase_IP::get_communication_table() const
    {
        return mCommunicationTable;
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Interpolation_Mesh_DataBase_IP::get_element_indices_in_block_set( uint aSetIndex )
    {
        MORIS_ERROR( 0, "get_element_indices_in_block_set not implemented for Interpolation_Mesh_DataBase_IP" );
        return { {} };
        // return mIPMesh.get_element_indices_in_block_set( aSetIndex );
    }

    // ----------------------------------------------------------------------------

    enum CellTopology
    Interpolation_Mesh_DataBase_IP::get_blockset_topology( const std::string& aSetName )
    {
        MORIS_ERROR( 0, "get_blockset_topology not implemented for Interpolation_Mesh_DataBase_IP" );
        return CellTopology::UNDEFINED;
        // return mIPMesh.get_blockset_topology( aSetName );
    }

    // ----------------------------------------------------------------------------
    enum CellShape
    Interpolation_Mesh_DataBase_IP::get_IG_blockset_shape( const std::string& aSetName )
    {
        MORIS_ERROR( 0, "get_IG_blockset_shape not implemented for Interpolation_Mesh_DataBase_IP" );
        return CellShape::UNDEFINED;
        // return mIPMesh.get_IG_blockset_shape( aSetName );
    }

    // ----------------------------------------------------------------------------

    enum CellShape
    Interpolation_Mesh_DataBase_IP::get_IP_blockset_shape( const std::string& aSetName )
    {
        MORIS_ERROR( 0, "get_IP_blockset_shape not implemented for Interpolation_Mesh_DataBase_IP" );
        return CellShape::UNDEFINED;
        // return mIPMesh.get_IP_blockset_shape( aSetName );
    }

    // ----------------------------------------------------------------------------

    uint
    Interpolation_Mesh_DataBase_IP::get_element_owner( moris_index aElementIndex ) const
    {
        MORIS_ERROR( 0, "get_element_owner not implemented for Interpolation_Mesh_DataBase_IP" );
        return 0;
        // return mIPMesh.get_element_owner( aElementIndex );
    }

    // ----------------------------------------------------------------------------
    std::unordered_map< moris_id, moris_index >
    Interpolation_Mesh_DataBase_IP::get_vertex_glb_id_to_loc_vertex_ind_map() const
    {
        return mVertexGlobalIdToLocalIndex;
    }

    // ----------------------------------------------------------------------------
    Vertex&
    Interpolation_Mesh_DataBase_IP::get_mtk_vertex( moris_index aVertexIndex )
    {
        MORIS_ASSERT( aVertexIndex < (moris_index)mVertices.size(), "index of the vertex specified exceeds the bounds" );
        return mVertices( aVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vertex const &
    Interpolation_Mesh_DataBase_IP::get_mtk_vertex( moris_index aVertexIndex ) const
    {
        MORIS_ASSERT( aVertexIndex < (moris_index)mVertices.size(), "index of the vertex specified exceeds the bounds" );
        return mVertices( aVertexIndex );
    }
    //--------------------------------------------------------------------------------------------------------------

    void
    Interpolation_Mesh_DataBase_IP::get_adof_map(
            const uint                    aBSplineIndex,
            map< moris_id, moris_index >& aAdofMap ) const
    {
        aAdofMap = mAdofMap( aBSplineIndex );
    }

    //--------------------------------------------------------------------------------------------------------------
    uint
    Interpolation_Mesh_DataBase_IP::get_num_interpolations()
    {
        return mMeshIndices.max() + 1;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Interpolation_Mesh_DataBase_IP::get_max_num_coeffs_on_proc( uint aDiscretizationIndex ) const
    {
        return mAdofMap( aDiscretizationIndex ).size();
    }

    //--------------------------------------------------------------------------------------------------------------

    mtk::Cell&
    Interpolation_Mesh_DataBase_IP::get_mtk_cell( moris_index aCellIndex )
    {
        MORIS_ASSERT( aCellIndex < (moris_index)mCells.size(), "index of the vertex specified exceeds the bounds" );
        return mCells( aCellIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    mtk::Cell const &
    Interpolation_Mesh_DataBase_IP::get_mtk_cell( moris_index aCellIndex ) const
    {
        MORIS_ASSERT( aCellIndex < (moris_index)mCells.size(), "index of the vertex specified exceeds the bounds" );
        return mCells( aCellIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id
    Interpolation_Mesh_DataBase_IP::get_entity_id( enum EntityRank aEntityRank, moris_index aEntityIndex ) const
    {
        switch ( aEntityRank )
        {
            case EntityRank::NODE:
            {
                MORIS_ASSERT( aEntityIndex < (moris_index)mVertices.size(), "index of the vertex specified exceeds the bounds" );
                return mVertexIdList( aEntityIndex );
                break;
            }
            case EntityRank::ELEMENT:
            {
                MORIS_ASSERT( aEntityIndex < (moris_index)mCells.size(), "index of the vertex specified exceeds the bounds" );
                return mCellIdList( aEntityIndex );
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }
                return 0;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id
    Interpolation_Mesh_DataBase_IP::get_entity_owner( enum EntityRank aEntityRank, moris_index aEntityIndex ) const
    {
        switch ( aEntityRank )
        {
            case EntityRank::NODE:
            {
                MORIS_ASSERT( aEntityIndex < (moris_index)mVertices.size(), "index of the vertex specified exceeds the bounds" );
                return mVertexOwnerList( aEntityIndex );
                break;
            }
            case EntityRank::ELEMENT:
            {
                MORIS_ASSERT( aEntityIndex < (moris_index)mCells.size(), "index of the vertex specified exceeds the bounds" );
                return mCellOwnerList( aEntityIndex );
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }
                return 0;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    moris::real*
    Interpolation_Mesh_DataBase_IP::get_vertex_coords_ptr( moris_index aVertexIndex )
    {
        moris::real* tCoordPointer = const_cast< real* >( mVertexCoordinates.colptr( aVertexIndex ) );

        return tCoordPointer;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vertex**
    Interpolation_Mesh_DataBase_IP::get_cell_vertices( moris_index aCellIndex )
    {
        return mCellToVertices.memptr() + mCellToVertexOffSet( aCellIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vertex_Interpolation**
    Interpolation_Mesh_DataBase_IP::get_vertex_interpolation( moris_index aVertexIndex )
    {
        return mVertexInterpoltionsPtrs.memptr() + mMeshIndices.numel() * aVertexIndex;
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Interpolation_Mesh_DataBase_IP::get_local_mesh_index( const uint aBsplineMeshIndex )
    {
        std::unordered_map< moris_index, moris_index >::const_iterator tLoc =
                mGlobalMeshIndexToLocalMeshIndex.find( aBsplineMeshIndex );

        MORIS_ASSERT( tLoc != mGlobalMeshIndexToLocalMeshIndex.end(),
                "Interpolation_Mesh_DataBase_IP::get_local_mesh_index - Bspline mesh index not found: %d",
                aBsplineMeshIndex );

        return tLoc->second;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id*
    Interpolation_Mesh_DataBase_IP::get_basis_ids( moris_index aVertexIndex, moris_index aOrder )
    {
        // moris_id* tBasisIdPointer = const_cast< moris_id* >( mBasisIds( aOrder ).memptr() + mOffSetTMatrix( aOrder )( aVertexIndex ) );
        moris_id* tBasisIdPointer = mBasisIds( aOrder ).memptr() + mOffSetTMatrix( aOrder )( aVertexIndex );
        return tBasisIdPointer;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index*
    Interpolation_Mesh_DataBase_IP::get_basis_indicies( moris_index aVertexIndex, moris_index aOrder )
    {
        moris_id* tBasisIndPointer = mBasisIndices( aOrder ).memptr() + mOffSetTMatrix( aOrder )( aVertexIndex );
        return tBasisIndPointer;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris::real*
    Interpolation_Mesh_DataBase_IP::get_basis_weights( moris_index aVertexIndex, moris_index aOrder )
    {
        real* tBasisWeightsPointer = const_cast< real* >( mWeights( aOrder ).memptr() + mOffSetTMatrix( aOrder )( aVertexIndex ) );

        return tBasisWeightsPointer;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id*
    Interpolation_Mesh_DataBase_IP::get_basis_owners( moris_index aVertexIndex, moris_index aOrder )
    {
        return mBasisOwners( aOrder ).memptr() + mOffSetTMatrix( aOrder )( aVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index
    Interpolation_Mesh_DataBase_IP::get_basis_length( moris_index aVertexIndex, moris_index aOrder )
    {
        return mOffSetTMatrix( aOrder )( aVertexIndex + 1 ) - mOffSetTMatrix( aOrder )( aVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< mtk::Cell_Info >
    Interpolation_Mesh_DataBase_IP::get_cell_info_sp( moris_index aEntityIndex ) const
    {
        return mCellInfo;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris::Memory_Map
    Interpolation_Mesh_DataBase_IP::get_memory_usage()
    {
        moris::Memory_Map tMemoryMap;

        tMemoryMap.mMemoryMapData.reserve( 6 );

        // tMemoryMap.mMemoryMapData["misc"] = sizeof( mSpatilDim ) + sizeof( mIPDataBase ) + sizeof( mCellInfo );

        // tMemoryMap.mMemoryMapData["vertices"] = moris::internal_capacity( mVertices )
        //                                         + moris::internal_capacity( mVertexInterpolations )
        //                                         + mVertexInterpoltionsPtrs.capacity() * ( 1 + sizeof( void* ) );

        // tMemoryMap.mMemoryMapData["cells"] = mCellToVertices.capacity() * ( 1 + sizeof( moris_index ) )
        //                                      + moris::internal_capacity( mCells );

        // tMemoryMap.mMemoryMapData["id-owner-list"] = mVertexIdList.capacity() * ( 1 + sizeof( moris_id ) )
        //                                              + mCellIdList.capacity() * ( 1 + sizeof( moris_id ) )
        //                                              + mVertexOwnerList.capacity() * ( 1 + sizeof( moris_id ) )
        //                                              + mCellOwnerList.capacity() * ( 1 + sizeof( moris_id ) );

        // // Determine the map size for the
        // size_t tVertexGlbToLocalMapSize = 0;

        // uint  tBucketCount = mVertexGlobalIdToLocalIndex.bucket_count();
        // float tLoadFactor  = mVertexGlobalIdToLocalIndex.max_load_factor();
        // if ( tLoadFactor > 1.0 )
        // {
        //     tVertexGlbToLocalMapSize = (size_t)tBucketCount * tLoadFactor;
        // }
        // else
        // {
        //     tVertexGlbToLocalMapSize = (size_t)tBucketCount;
        // }

        // size_t tDofMapSize = 0;
        // for ( const auto& iMap : mAdofMap )
        // {
        //     tDofMapSize += ( sizeof( moris_id ) + sizeof( moris_id ) + sizeof( std::_Rb_tree_node_base ) ) * iMap.size();//+ sizeof(std::_Rb_tree);
        // }

        // tMemoryMap.mMemoryMapData["maps"] = tVertexGlbToLocalMapSize + +sizeof( mCommunicationTable ) + mCommunicationTable.capacity()
        //                                     + tDofMapSize;

        // // moris::Memory_Map tIPDataBaseMM;

        // // tIPDataBaseMM = mIPDataBase->get_memory_usage();

        // // tMemoryMap.mMemoryMapData["database"] = tIPDataBaseMM.sum();

        // tMemoryMap.par_print( "IP Mesh" );

        return tMemoryMap;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Interpolation_Mesh_DataBase_IP::free_memory()
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< IndexMat >
    Interpolation_Mesh_DataBase_IP::get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const
    {
        MORIS_ERROR( 0, " get_elements_connected_to_element_and_face_ind_loc_inds not implemnted" );
        return { {} };
    }

}    // namespace moris::mtk

#ifndef MORIS_MESH_CL_MESH_HPP_
#define MORIS_MESH_CL_MESH_HPP_

// MORIS project header files.
//#include "cl_Mesh_Mtk.hpp"

#include <climits>
#include <string>
#include <vector>

#include "algorithms.hpp"
#include "fn_assert.hpp" // ASR/src
#include "cl_Cell.hpp" // CON/src
#include "typedefs.hpp" // COR/src
#include "cl_Logger.hpp" // IOS/src
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_Interpolation.hpp" // TOL/src
#include "cl_Pairing.hpp" // TOL/src
#include "cl_Entity_Tracker.hpp" // May move this to the mesh library // XTK/src
#include "cl_Mesh_Enums.hpp" // MTK/src
#include "cl_Database.hpp" // MTK/src
#include "cl_STK_Implementation.hpp" // STK/src
#include "cl_Debug.hpp" // TOL/src

namespace moris
{

class mesh
{
protected:

public:
    database* mDatabase;

    /**
     * Mesh constructor (mesh generated internally or obtained from an Exodus file )
     *
     * @param[in] aFileName  .................    String with mesh file name.
     */
    mesh(
            enum MeshType   aMeshType,
            std::string     aFileName,
            MtkSetsInfo*    aSetsInfo=NULL )
    {
        switch (aMeshType)
        {
        case(MeshType::MTK):
        {
            mDatabase = new STK_Implementation( aFileName, aSetsInfo );
            break;
        }
        default:
        {
            MORIS_ASSERT( 0, "Specified mesh type not supported by MORIS" );
        }
        }
    }

    /**
     * Mesh destructor.
     */
    virtual
    ~mesh()
    {
        delete mDatabase;
    }

    ///////////////////////////////////////////
    ///  Beginning of function declarations  ///
    ///////////////////////////////////////////

    /**
     * Create stk mesh with information given by the user.
     *
     * @param[in] aFileName   .................   String with mesh file name.
     *
     */
    void
    build_mesh(
            std::string    aFileName ,
            MtkSetsInfo*   aSetsInfo )
    {
        return mDatabase->build_mesh( aFileName, aSetsInfo );
    }

    mesh(
            enum MeshType   aMeshType,
            MtkMeshData     aMeshData )

    {
        switch ( aMeshType )
        {
        case ( MeshType::MTK ):
        {
            mDatabase = new STK_Implementation( aMeshData );
            break;
        }
        default:
        {
            MORIS_ASSERT( 0, "Specified mesh type not supported by MORIS" );
        }
        }
    }

    /**
     * Create a MORIS HMR mesh .
     *
     *
     */
    mesh( database* aMeshDatabase ):
        mDatabase( aMeshDatabase )
    {

    }

    /**
     * Create a MORIS mesh with information given by the user.
     *
     * @param[in] aMeshType   .................   string specifying which underlying mesh database should be used
     * @param[in] aMeshData   .................   struct with the following mesh information
     *
     */
    void
    build_mesh(
            MtkMeshData   aMeshData )
    {
        return mDatabase->build_mesh( aMeshData );
    }

    /*
     * Create an exodus mesh database with the specified
     * filename.
     *
     * @param[in] filename The full pathname to the file which will be
     *   created and the mesh data written to. If the file already
     *   exists, it will be overwritten.
     * @param[out] output_file_index
     *
     *   Description from create_output_mesh() in StkMeshIoBroker.hpp
     */
    void
    create_output_mesh(
            std::string   &aFileName )
    {
        return mDatabase->create_output_mesh( aFileName );
    }

    /**
     * Return number of elements in mesh.
     *
     * @return Number of elements.
     */
    uint
    get_num_entities_universal(
            enum EntityRank   aEntityRank ) const
    {
        return mDatabase->get_num_entities_universal( aEntityRank );
    }

    uint
    get_num_elems( ) const
    {
        return mDatabase->get_num_elems( );
    }

    uint
    get_num_faces( ) const
    {
        return mDatabase->get_num_faces( );
    }

    uint
    get_num_edges( ) const
    {
        return mDatabase->get_num_edges( );
    }

    uint
    get_num_nodes( ) const
    {
        return mDatabase->get_num_nodes( );
    }

    uint
    get_num_spatial_dims( ) const
    {
        return mDatabase->get_num_spatial_dims( );
    }

    uint
    get_num_elems_current_proc( ) const
    {
        return mDatabase->get_num_elems_current_proc( );
    }

    uint
    get_num_faces_current_proc( ) const
    {
        return mDatabase->get_num_faces_current_proc( );
    }

    uint
    get_num_edges_current_proc( ) const
    {
        return mDatabase->get_num_edges_current_proc( );
    }

    uint
    get_num_nodes_current_proc( ) const
    {
        return mDatabase->get_num_nodes_current_proc( );
    }
    uint
    get_num_entities_aura(
            enum EntityRank   aEntityRank ) const
    {
        return mDatabase->get_num_entities_aura( aEntityRank );
    }

    uint
    get_num_entities_locally_owned_globally_shared(
            enum EntityRank   aEntityRank ) const
    {
        return mDatabase->get_num_entities_locally_owned_globally_shared( aEntityRank );
    }

    /**
     * Return number of nodes in current processor.
     *
     * @return Number of nodes.
     */
    Mat< uint >
    get_entities_universal(
            enum EntityRank   aEntityRank ) const
    {
        return mDatabase->get_entities_universal( aEntityRank );
    }

    Mat< uint >
    get_entities_glb_shared_current_proc(
            EntityRank   aEntityRank ) const
    {
        return mDatabase->get_entities_glb_shared_current_proc( aEntityRank );
    }

    Mat< uint >
    get_entities_owned_current_proc(
            EntityRank   aEntityRank ) const
    {
        return mDatabase->get_entities_owned_current_proc( aEntityRank );
    }

    Mat<uint>
    get_entities_in_aura(
            EntityRank   aEntityRank ) const
    {
        return mDatabase->get_entities_in_aura( aEntityRank );
    }

    Mat< uint >
    get_entities_owned_and_shared_by_current_proc(
            EntityRank   aEntityRank ) const
    {
        return mDatabase->get_entities_owned_and_shared_by_current_proc( aEntityRank );
    }

    /**
     * Return number of nodes in current processor.
     *
     * @return Number of nodes.
     */
    uint
    parallel_owner_rank_by_entity_id(
            uint              aEntityIndex,
            enum EntityRank   aEntityRank ) const
    {
        return mDatabase->parallel_owner_rank_by_entity_id( aEntityIndex, aEntityRank );
    }

    uint
    parallel_owner_rank_by_entity_index(
            uint              aEntityIndex,
            enum EntityRank   aEntityRank ) const
    {
        return mDatabase->parallel_owner_rank_by_entity_index( aEntityIndex, aEntityRank );
    }

    Mat< uint >
    get_procs_sharing_entity_by_id(
            uint              aEntityID,
            enum EntityRank   aEntityRank ) const
    {
        return mDatabase->get_procs_sharing_entity_by_id( aEntityID, aEntityRank );
    }

    Mat<uint>
    get_procs_sharing_entity_by_index(
            uint              aEntityIndex,
            enum EntityRank   aEntityRank ) const
    {
        return mDatabase->get_procs_sharing_entity_by_index( aEntityIndex, aEntityRank );
    }

    uint
    get_parallel_size() const
    {
        return mDatabase->get_parallel_size( );
    }

    uint
    get_parallel_rank() const
    {
        return mDatabase->get_parallel_rank( );
    }

    /**
     * Return number of faces in element of a particular topology.
     *
     * @return Number of face.
     */
    uint
    get_elem_topology_num_faces(
            uint aElemId ) const
    {
        return mDatabase->get_elem_topology_num_faces( aElemId );
    }

    uint
    get_elem_topology_num_edges(
            uint aElemId ) const
    {
        return mDatabase->get_elem_topology_num_edges( aElemId );
    }

    uint
    get_elem_topology_num_nodes(
            uint aElemId ) const
    {
        return mDatabase->get_elem_topology_num_nodes( aElemId );
    }

    uint
    get_face_topology_num_edges(
            uint aFaceId ) const
    {
        return mDatabase->get_face_topology_num_edges( aFaceId );
    }

    uint
    get_face_topology_num_nodes(
            uint aFaceId ) const
    {
        return mDatabase->get_face_topology_num_nodes( aFaceId );
    }

    uint
    get_edge_topology_num_nodes(
            uint aEdgeId ) const
    {
        return mDatabase->get_edge_topology_num_nodes( aEdgeId );
    }

    /**
     * Get elements connected to a face
     *
     * @param[in]  aFaceId    ............................   Entity Id
     * @param[out] aElementsConnectedToFace   ............   Connected entities
     *
     */
    Mat< uint >
    get_elements_connected_to_element(
            uint   const aElemId ) const
    {
        return mDatabase->get_elements_connected_to_element( aElemId );
    }

    Mat< uint >
    get_elements_connected_to_face(
            uint   const aFaceId ) const
    {
        return mDatabase->get_elements_connected_to_face( aFaceId );
    }

    Mat< uint >
    get_elements_connected_to_edge(
            uint   const aEdgeId ) const
    {
        return mDatabase->get_elements_connected_to_edge( aEdgeId );
    }

    Mat< uint >
    get_elements_connected_to_node(
            uint   const aNodeId ) const
    {
        return mDatabase->get_elements_connected_to_node( aNodeId );
    }

    Mat< uint >
    get_faces_connected_to_element(
            uint   const aElementId ) const
    {
        return mDatabase->get_faces_connected_to_element( aElementId );
    }

    Mat< uint >
    get_faces_connected_to_edge(
            uint   const aEdgeId ) const
    {
        return mDatabase->get_faces_connected_to_edge( aEdgeId );
    }

    Mat< uint >
    get_faces_connected_to_node(
            uint   const aNodeId ) const
    {
        return mDatabase->get_faces_connected_to_node( aNodeId );
    }

    Mat< uint >
    get_edges_connected_to_element(
            uint   const aElementId ) const
    {
        return mDatabase->get_edges_connected_to_element( aElementId );
    }

    Mat< uint >
    get_edges_connected_to_face(
            uint   const aFaceId ) const
    {
        return mDatabase->get_edges_connected_to_face( aFaceId );
    }

    Mat< uint >
    get_edges_connected_to_node(
            uint   const aNodeId ) const
    {
        return mDatabase->get_edges_connected_to_node( aNodeId );
    }

    Mat< uint >
    get_nodes_connected_to_element(
            uint   const aElementId ) const
    {
        return mDatabase->get_nodes_connected_to_element( aElementId );
    }

    Mat< uint >
    get_nodes_connected_to_face(
            uint   const aFaceId ) const
    {
        return mDatabase->get_nodes_connected_to_face( aFaceId );
    }

    Mat< uint >
    get_nodes_connected_to_edge(
            uint   const aEdgeId ) const
    {
        return mDatabase->get_nodes_connected_to_edge( aEdgeId );
    }

    /**
     * Get nodes connected to an edge
     *
     * @param[in]  aNodeId       ......................   Entity Id
     * @param[in]  aNodeId       ......................   Entity Id
     * @param[in]  aNodeId       ......................   Entity Id
     * @param[out] aEntitiesConnectedToEdge   ........   Connected entities
     *
     */
    Mat< uint >
    entities_connected_to_given_entity(
            uint         const aEntityId ,
            EntityRank   const aInputEntityRank,
            EntityRank   const aOutputEntityRank ) const
    {
        return mDatabase->entities_connected_to_given_entity( aEntityId, aInputEntityRank, aOutputEntityRank );
    }

    Mat< uint >
    get_nodes_in_node_set(
            uint   const aNodeSetId ) const
    {
        return mDatabase->get_nodes_in_node_set( aNodeSetId );
    }

    Mat< uint >
    get_nodes_in_side_set(
            uint   const aSideSetId ) const
    {
        return mDatabase->get_nodes_in_side_set( aSideSetId );
    }

    Mat< uint >
    get_nodes_in_block_set(
            uint   const aBlockSetId ) const
    {
        return mDatabase->get_nodes_in_block_set( aBlockSetId );
    }

    Mat< uint >
    get_edges_in_side_set(
            uint   const aSideSetId ) const
    {
        return mDatabase->get_edges_in_side_set( aSideSetId );
    }

    Mat< uint >
    get_edges_in_block_set(
            uint   const aSideSetId ) const
    {
        return mDatabase->get_edges_in_block_set( aSideSetId );
    }

    Mat< uint >
    get_faces_in_side_set(
            uint   const aSideSetId ) const
    {
        return mDatabase->get_faces_in_side_set( aSideSetId );
    }

    Mat< uint >
    get_faces_in_block_set(
            uint   const aBlockSetId ) const
    {
        return mDatabase->get_faces_in_block_set( aBlockSetId );
    }

    uint
    get_entity_index(
            enum EntityRank   aEntityRank,
            uint              aEntityID ) const
    {
        return mDatabase->get_entity_index( aEntityRank, aEntityID );
    }

    Mat< uint >
    get_entity_local_ids_connected_to_entity(
            uint         const aEntityId,
            EntityRank   const aInputEntityRank,
            EntityRank   const aOutputEntityRank ) const
    {
        return mDatabase->get_entity_local_ids_connected_to_entity( aEntityId, aInputEntityRank, aOutputEntityRank );
    }

    /**
     * Return coordinates of all nodes in mesh
     *
     * @return coordinates.
     */
    Mat< real >
    get_all_nodes_coords( ) const
    {
        return mDatabase->get_all_nodes_coords( );
    }

    Mat< real >
    get_all_nodes_coords_aura( ) const
    {
        return mDatabase->get_all_nodes_coords_aura( );
    }

    Mat< real >
    get_selected_nodes_coords(
            Mat< uint > aNodeIds ) const
    {
        return mDatabase->get_selected_nodes_coords( aNodeIds );
    }

    Mat< real >
    get_selected_nodes_coords_lcl_ind(
            Mat< uint > aNodeIds ) const
    {
        return mDatabase->get_selected_nodes_coords_lcl_ind( aNodeIds );
    }

    Mat< uint >
    get_node_ids_from_local_map(
            Mat< uint >   aLocalInds ) const
    {
        return mDatabase->get_node_ids_from_local_map( aLocalInds );
    }

    /*
     * Using a local node index, return a global node ID
     * @param[in]  aNodeInds - row vector of node index (processor unique)
     * @Return row vector of corresponding global node ID
     */
    Mat< uint >
    get_nodal_local_map( )
    {
        return mDatabase->get_nodal_local_map( );
    }

    Mat< uint >
    get_elemental_local_map( )
    {
        return mDatabase->get_elemental_local_map( );
    }

    Mat< uint >
    get_edge_local_map( )
    {
        return mDatabase->get_edge_local_map( );
    }

    Mat< uint >
    get_face_local_map( )
    {
        return mDatabase->get_face_local_map( );
    }

    /*
     * Using a local node index, return a global node ID
     * @param[in]  aNodeInds - row vector of node index (processor unique)
     * @Return row vector of corresponding global node ID
     */
    Mat< uint >
    get_nodal_owner_proc_map( )
    {
        return mDatabase->get_nodal_owner_proc_map( );
    }

    Mat< uint >
    get_elemental_owner_proc_map( )
    {
        return mDatabase->get_elemental_owner_proc_map( );
    }

    Mat< uint >
    get_edge_owner_proc_map( )
    {
        return mDatabase->get_edge_owner_proc_map( );
    }

    Mat< uint >
    get_face_owner_proc_map( )
    {
        return mDatabase->get_face_owner_proc_map( );
    }

    /*
     * Using a local node index, return a global node ID
     * @param[in]  aNodeInds - row vector of node index (processor unique)
     * @Return row vector of corresponding global node ID
     */
    Cell < Cell < uint > >
    get_nodes_shared_processors( )
    {
        return mDatabase->get_nodes_shared_processors( );
    };

    Cell < Cell < uint > >
    get_elements_shared_processors( )
    {
        return mDatabase->get_elements_shared_processors( );
    };

    Cell < Cell < uint > >
    get_edges_shared_processors( )
    {
        return mDatabase->get_edges_shared_processors( );
    };

    Cell < Cell < uint > >
    get_faces_shared_processors( )
    {
        return mDatabase->get_faces_shared_processors( );
    };

    /*
     * Returns a list of globally unique element ids
     * @param[in]  aNumNodes         - number of element ids requested
     * @param[out] aAvailableNodeIDs - list of globally unique element IDs
     */
    Mat< uint >
    generate_unique_elem_ids(
            uint aNumElems ) const
    {
        return mDatabase->generate_unique_elem_ids( aNumElems );
    }

    Mat< uint >
    generate_unique_face_ids(
            uint aNumFaces ) const
    {
        return mDatabase->generate_unique_face_ids( aNumFaces );
    }

    Mat< uint >
    generate_unique_edge_ids(
            uint aNumEdges ) const
    {
        return mDatabase->generate_unique_edge_ids( aNumEdges );
    }

    Mat< uint >
    generate_unique_node_ids(
            uint aNumNodes ) const
    {
        return mDatabase->generate_unique_node_ids( aNumNodes );
    }

    /*
     *
     */
    Mat< uint >
    get_field_entities(
            enum EntityRank   aNewEntityRank,
            std::string       aNewFieldName )
    {
        return mDatabase->get_field_entities( aNewEntityRank, aNewFieldName );
    }

    Mat< real >
    get_field_values(
            enum EntityRank aNewEntityRank,
            std::string aNewFieldName )
    {
        return mDatabase->get_field_values( aNewEntityRank, aNewFieldName );
    }

    uint
    get_entity_id_from_entity_key(
            enum EntityRank aEntityRank,
            uint aEntityKeyNumber ) const
    {
        return mDatabase->get_entity_id_from_entity_key( aEntityRank, aEntityKeyNumber );
    }

    uint
    get_entity_key_from_entity_id(
            enum EntityRank aEntityRank,
            uint aEntityIdNumber ) const
    {
        return mDatabase->get_entity_key_from_entity_id( aEntityRank, aEntityIdNumber );
    }

    Mat< uint >
    get_set_entity_ids(
            enum EntityRank aNewEntityRank,
            std::string aNewFieldName )const
    {
        return mDatabase->get_set_entity_ids( aNewEntityRank, aNewFieldName );
    }
    
    /*
     *
     */
    Mat< uint >
    get_intersected_entities_field_set(
            enum EntityRank   aEntityRank,
            std::string       aFieldName,
            std::string       aSetName )
    {
        return mDatabase->get_intersected_entities_field_set( aEntityRank, aFieldName, aSetName);
    }
    
    Mat< real >
    get_intersected_data_field_set(
            enum EntityRank   aEntityRank,
            std::string       aFieldName,
            std::string       aSetName )
    {
        return mDatabase->get_intersected_data_field_set( aEntityRank, aFieldName, aSetName);
    }

    /**
     * Get duplicates of a coordinate and id list
     *
     * @param[in]  aCoord         .... Coordinate list with x,y,z
     * @param[in]  aId            .... Id list of the coordinates
     * @param[out] duplicate_list .... Shows the duplicates [Position(i) Position(j)]
     *
     */
    Mat< uint >
    duplicate_node_coord_check()
    {
        return mDatabase->duplicate_node_coord_check( );
    }

    Mat< uint >
    duplicate_node_coord_check(
            Mat< real >&   aCoord )
    {
        return mDatabase->duplicate_node_coord_check( aCoord );
    }

    Mat< uint >
    duplicate_node_coord_and_id_check(
            Mat< real >&   aCoord,
            Mat< uint >&   aId )
    {
        return mDatabase->duplicate_node_coord_and_id_check( aCoord, aId );
    }

    Mat< uint >
    duplicate_node_coord_and_id_check(
            Cell< Mat< real > >&   aCoord,
            Cell< Mat< uint > >&   aId )
    {
        return mDatabase->duplicate_node_coord_and_id_check( aCoord, aId );
    }

    Mat< uint >
    duplicate_node_coord_and_id_check_problems(
            Mat< real >&   aCoord,
            Mat< uint >&   aId )
    {
        return mDatabase->duplicate_node_coord_and_id_check_problems( aCoord, aId );
    }

    Mat< uint >
    duplicate_node_coord_and_id_check_problems(
            Cell< Mat< real > >&   aCoord,
            Cell< Mat< uint > >&   aId )
    {
        return mDatabase->duplicate_node_coord_and_id_check_problems( aCoord, aId );
    }

    //////////////////////////
    // Deprecated functions //
    //////////////////////////

    Mat< real >
    interpolate_to_location_on_entity(
            enum EntityRank   aParentEntityRank,
            uint              aParentEntityIndex,
            Mat<real>         aLclCoord )
    {
        return mDatabase->interpolate_to_location_on_entity( aParentEntityRank, aParentEntityIndex, aLclCoord );
    }

    /////////////////////////////
};
}   // namespace moris

#endif /* MORIS_MESH_CL_MESH_HPP_ */

#ifndef MORIS_MESH_CL_DATABASE_HPP_
#define MORIS_MESH_CL_DATABASE_HPP_

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
#include "cl_Mesh_Enums.hpp" // MTK/src
#include "cl_Debug.hpp" // TOL/src

namespace moris
{
    // Arguments list for for creating a mesh form data

    //////////////////////////
    // STRUC FOR BLOCK SET  //
    //////////////////////////
    struct MtkBlockSetsInfo
    {
        Mat< uint >*          BSetInds;
        Cell< std::string >   BSetNames;

        MtkBlockSetsInfo():
            BSetInds(),
            BSetNames(){}
    };

    /////////////////////////
    // STRUC FOR SIDE SET  //
    /////////////////////////
    struct MtkSideSetsInfo
    {
        Cell< Mat< uint > >*    ElemIdsAndSideOrds;
        Cell< std::string >     SSetNames;

        MtkSideSetsInfo():
            ElemIdsAndSideOrds(),
            SSetNames(){}
    };

    //////////////////////////
    // STRUC FOR NODE SET  //
    //////////////////////////
    struct MtkNodeSetsInfo
    {
        Cell< Mat< uint > >*    EntIds;
        Cell< std::string >     NSetNames;

        MtkNodeSetsInfo():
            EntIds(),
            NSetNames(){}
    };

    ///////////////////////////////
    // STRUC FOR SETS CONTAINER  //
    ///////////////////////////////
    struct MtkSetsInfo
    {
        MtkNodeSetsInfo* NodeSetsInfo;
        MtkSideSetsInfo* SideSetsInfo;
        MtkBlockSetsInfo* BlockSetsInfo;

        MtkSetsInfo():
            NodeSetsInfo(),
            SideSetsInfo(),
            BlockSetsInfo(){}
    };

    ///////////////////////
    // STRUC FOR FIELDS  //
    ///////////////////////
    //template <typename T>
    //template< typename Variant = boost::variant< bool, sint, real, const char* > >
    struct MtkFieldsInfo
    {
        Cell< Mat< real > >*    FieldsData;
        Cell< std::string >     FieldsName;
        Cell< enum EntityRank > FieldsRank;
        Cell< std::string >*    SetsOwner;

        MtkFieldsInfo():
            FieldsData(),
            FieldsName(),
            FieldsRank(),
            SetsOwner() {}
    };

    //////////////////////////
    // STRUC FOR BLOCK SET  //
    //////////////////////////
    struct MtkMeshData
    {
        uint*        SpatialDim ;
        Mat< uint >* ElemConn;
        Mat< uint >* EntProcOwner;
        Mat< real >* NodeCoords;
        Mat< uint >* LocaltoGlobalElemMap;
        Mat< uint >* LocaltoGlobalNodeMap;
        bool         CreateAllEdgesAndFaces;
        MtkFieldsInfo* FieldsInfo;
        MtkSetsInfo* SetsInfo;

        MtkMeshData():
            SpatialDim(),
            ElemConn(),
            EntProcOwner(),
            NodeCoords(),
            LocaltoGlobalElemMap(),
            LocaltoGlobalNodeMap(),
            CreateAllEdgesAndFaces(false),
            FieldsInfo(),
            SetsInfo(){}
    };

    class database
    {
    protected:

    public:

        /**
         * Mesh constructor (mesh generated internally or obtained from an Exodus file )
         *
         * @param[in] aFileName  .................    String with mesh file name.
         */
        database()
    {
    }

        /**
         * Mesh destructor.
         */
        virtual ~database() = default;
        /// Begining of function declarations

        /**
         * Create stk mesh with information given by the user.
         *
         * @param[in] aFileName   .................   String with mesh file name.
         *
         */
        virtual void build_mesh( std::string   aFileName ,
                MtkSetsInfo*   aSetsInfo)
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
        };

        /**
         * Create a MORIS mesh with information given by the user.
         *
         * @param[in] aMeshType   .................   string specifying which underlying mesh database should be used
         * @param[in] aMeshData   .................   struct with the following mesh information
         * @param[in] SpatialDim   ...............   problem dimensions (1D = 1, 2D = 2 , or 3D = 3).
         * @param[in] ElemConn     ...............   table containing element connectivity.
         * @param[in] NodeCoords         .........   node coordinates
         * @param[in] LocaltoGlobalNodeMap   .....   node local ind to global id map
         * @param[in] LocaltoGlobalElemMap   .....   element local ind to global id map
         * @param[in] EntProcOwner................   number of processors required.
         * @param[in] PartNames   ................   information of parts where elements will be stored.
         * @param[in] FieldsData  ................   vectors with field data for all fields
         * @param[in] FieldsRank  ................   entity types on which the fields are acting
         * @param[in] FieldsName  ................   names for all fields
         *
         */
        virtual void
        build_mesh(
                MtkMeshData   aMeshData )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
        };

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
        virtual void
        create_output_mesh(
                std::string&   aFileName )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
        };

        /**
         * Return number of elements in mesh.
         *
         * @return Number of elements.
         */
        virtual uint
        get_num_entities_universal(
                enum EntityRank   aEntityRank ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_num_elems( ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_num_faces( ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_num_edges( ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_num_nodes( ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_num_spatial_dims( ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_num_elems_current_proc( ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_num_faces_current_proc( ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_num_edges_current_proc( ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_num_nodes_current_proc( ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_num_entities_aura(
                enum EntityRank   aEntityRank ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_num_entities_locally_owned_globally_shared(
                enum EntityRank   aEntityRank) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };
        /**
         * Return number of nodes in current processor.
         *
         * @return Number of nodes.
         */
        virtual Mat< uint >
        get_entities_universal(
                enum EntityRank   aEntityRank ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_entities_glb_shared_current_proc(
                EntityRank   aEntityRank ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_entities_owned_current_proc(
                EntityRank   aEntityRank ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_entities_in_aura(
                EntityRank   aEntityRank ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_entities_owned_and_shared_by_current_proc(
                EntityRank   aEntityRank ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        /**
         * Return number of nodes in current processor.
         *
         * @return Number of nodes.
         */

        virtual uint
        parallel_owner_rank_by_entity_id(
                uint              aEntityIndex,
                enum EntityRank   aEntityRank ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        parallel_owner_rank_by_entity_index(
                uint              aEntityIndex,
                enum EntityRank   aEntityRank ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual Mat< uint >
        get_procs_sharing_entity_by_id(
                uint              aEntityID,
                enum EntityRank   aEntityRank ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_procs_sharing_entity_by_index(
                uint              aEntityIndex,
                enum EntityRank   aEntityRank ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
         };

        virtual uint
        get_parallel_size() const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_parallel_rank() const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };
        /**
         * Return number of faces in element of a particular topology.
         *
         * @return Number of face.
         */
        virtual uint
        get_elem_topology_num_faces(
                uint   aElemId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_elem_topology_num_edges(
                uint   aElemId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_elem_topology_num_nodes(
                uint   aElemId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_face_topology_num_edges(
                uint   aFaceId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_face_topology_num_nodes(
                uint   aFaceId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_edge_topology_num_nodes(
                uint   aEdgeId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        /**
         * Get elements connected to a face
         *
         * @param[in]  aFaceId    ............................   Entity Id
         * @param[out] aElementsConnectedToFace   ............   Connected entities
         *
         */
        virtual Mat< uint >
        get_elements_connected_to_element(
                uint const   aElemId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_elements_connected_to_face(
                uint const   aFaceId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_elements_connected_to_edge(
                uint const   aEdgeId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_elements_connected_to_node(
                uint const   aNodeId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_faces_connected_to_element(
                uint const   aElementId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_faces_connected_to_edge(
                uint const   aEdgeId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_faces_connected_to_node(
                uint const   aNodeId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_edges_connected_to_element(
                uint const   aElementId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_edges_connected_to_face(
                uint const   aFaceId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_edges_connected_to_node(
                uint const   aNodeId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_nodes_connected_to_element(
                uint const   aElementId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_nodes_connected_to_face(
                uint const   aFaceId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_nodes_connected_to_edge(
                uint const   aEdgeId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        /**
         * Get nodes connected to an edge
         *
         * @param[in]  aNodeId       ......................   Entity Id
         * @param[in]  aNodeId       ......................   Entity Id
         * @param[in]  aNodeId       ......................   Entity Id
         * @param[out] aEntitiesConnectedToEdge   ........   Connected entities
         *
         */
        virtual Mat< uint >
        entities_connected_to_given_entity(
                uint const         aEntityId,
                EntityRank const   aInputEntityRank,
                EntityRank const   aOutputEntityRank ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_nodes_in_node_set(
                uint const   aNodeSetId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_nodes_in_side_set(
                uint const   aSideSetId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_nodes_in_block_set(
                uint const   aBlockSetId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_edges_in_side_set(
                uint const   aSideSetId ) const
       {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
       };

        virtual Mat< uint >
        get_edges_in_block_set(
                uint const   aSideSetId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_faces_in_side_set(
                uint const   aSideSetId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_faces_in_block_set (
                uint const   aBlockSetId ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual uint
        get_entity_index(
                enum EntityRank   aEntityRank,
                uint              aEntityID) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        /**
         * Get nodes connected to an edge
         *
         * @param[in]  aNodeId       ......................   Entity Id
         * @param[in]  aNodeId       ......................   Entity Id
         * @param[in]  aNodeId       ......................   Entity Id
         * @param[out] aEntitiesConnectedToEdge   ........   Connected entities
         *
         */
        virtual Mat< uint >
        get_entity_local_ids_connected_to_entity(
                uint const         aEntityId,
                EntityRank const   aInputEntityRank,
                EntityRank const   aOutputEntityRank) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        /**
         * Return coordinates of all nodes in mesh
         *
         * @return coordinates.
         */
        virtual Mat< real >
        get_all_nodes_coords( ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< real >( 0, 0 );
        };

        virtual Mat< real >
        get_all_nodes_coords_aura( ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< real >( 0, 0 );
        };

        virtual Mat< real >
        get_selected_nodes_coords(
                Mat< uint > aNodeIds ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< real >( 0, 0 );
        };

        virtual Mat< real >
        get_selected_nodes_coords_lcl_ind(
                Mat< uint > aNodeIds ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< real >( 0, 0 );
        };

        virtual Mat < uint >
        get_node_ids_from_local_map(
                Mat< uint >   aLocalInds ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
         };

        /*
         * Using a local node index, return a global node ID
         * @param[in]  aNodeInds - row vector of node index (processor unique)
         * @Return row vector of corresponding global node ID
         */
        virtual Mat< uint >
        get_nodal_local_map( )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_elemental_local_map( )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_edge_local_map( )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_face_local_map( )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        /*
         * Using a local node index, return a global node ID
         * @param[in]  aNodeInds - row vector of node index (processor unique)
         * @Return row vector of corresponding global node ID
         */
        virtual Mat< uint >
        get_nodal_owner_proc_map( )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_elemental_owner_proc_map( )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_edge_owner_proc_map( )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        get_face_owner_proc_map( )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        /*
         * Using a local node index, return a global node ID
         * @param[in]  aNodeInds - row vector of node index (processor unique)
         * @Return row vector of corresponding global node ID
         */
        virtual Cell < Cell < uint > >
        get_nodes_shared_processors( )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Cell < Cell < uint > >( 1 );
        };

        virtual Cell < Cell < uint > >
        get_elements_shared_processors( )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Cell < Cell < uint > >( 1 );
        };

        virtual Cell < Cell < uint > >
        get_edges_shared_processors( )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Cell < Cell < uint > >( 1 );
        };

        virtual Cell < Cell < uint > >
        get_faces_shared_processors( )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Cell < Cell < uint > >( 1 );
        };
        /*
         * Returns a list of globally unique element ids
         * @param[in]  aNumNodes         - number of element ids requested
         * @param[out] aAvailableNodeIDs - list of globally unique element IDs
         */

        virtual Mat< uint >
        generate_unique_elem_ids(
                uint   aNumElems ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        generate_unique_face_ids(
                uint   aNumFaces ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        generate_unique_edge_ids(
                uint   aNumEdges ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        generate_unique_node_ids(
                uint   aNumNodes ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        /*
         *
         */

        virtual Mat< uint >
        get_field_entities(
                enum EntityRank   aNewEntityRank,
                std::string       aNewFieldName)
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< real >
        get_field_values(
                enum EntityRank   aNewEntityRank,
                std::string       aNewFieldName)
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< real >( 0, 0 );
        };

        virtual uint
        get_entity_id_from_entity_key(
                enum EntityRank   aEntityRank,
                uint              aEntityKeyNumber) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual uint
        get_entity_key_from_entity_id(
                enum EntityRank   aEntityRank,
                uint              aEntityIdNumber ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return 0;
        };

        virtual Mat< uint >
        get_set_entity_ids(
                enum EntityRank   aNewEntityRank,
                std::string       aNewFieldName ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        /*
         *
         */

        virtual Mat< uint >
        get_intersected_entities_field_set(
                enum EntityRank   aEntityRank,
                std::string       aFieldName,
                std::string       aSetName ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< real > get_intersected_data_field_set(
                enum EntityRank   aEntityRank,
                std::string       aFieldName,
                std::string       aSetName ) const
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< real >( 0, 0 );
        };

        /**
         * Get duplicates of a coordinate and id list
         *
         * @param[in]  aCoord         .... Coordinate list with x,y,z
         * @param[in]  aId            .... Id list of the coordinates
         * @param[out] duplicate_list .... Shows the duplicates [Position(i) Position(j)]
         *
         */

        virtual Mat< uint >
        duplicate_node_coord_check()
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        duplicate_node_coord_check(
                Mat< real >&   aCoord )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        duplicate_node_coord_and_id_check(
                Mat< real >&   aCoord,
                Mat< uint >&   aId )
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        duplicate_node_coord_and_id_check(
                Cell< Mat< real > >&   aCoord,
                Cell< Mat< uint > >&   aId)
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        duplicate_node_coord_and_id_check_problems(
                Mat< real >&   aCoord,
                Mat< uint >&   aId)
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        virtual Mat< uint >
        duplicate_node_coord_and_id_check_problems(
                Cell< Mat< real > >&   aCoord,
                Cell< Mat< uint > >&   aId)
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< uint >( 0, 0 );
        };

        //////////////////////////
        // Deprecated functions //
        //////////////////////////

        virtual Mat< real >
        interpolate_to_location_on_entity(
                enum EntityRank   aParentEntityRank,
                uint              aParentEntityIndex,
                Mat< real >         aLclCoord)
        {
            MORIS_ASSERT( 0, "Current Database does not support this functionality");
            return Mat< real >( 0, 0 );
        };

        /////////////////////////////
    };
}   // namespace moris

#endif /* MORIS_MESH_CL_MESH_HPP_ */

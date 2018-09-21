#ifndef MORIS_MESH_CL_STK_IMPLEMENTATION_HPP_
#define MORIS_MESH_CL_STK_IMPLEMENTATION_HPP_

// Third-party header files.
#include <stk_io/StkMeshIoBroker.hpp>     // for StkMeshIoBroker
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <stk_mesh/base/MetaData.hpp>     // for MetaData
#include <stk_mesh/base/BulkData.hpp>     // for BulkData
#include <stk_mesh/base/Selector.hpp>     // for Selector
#include <stk_mesh/base/FEMHelpers.hpp>   // for Selector
#include "stk_io/DatabasePurpose.hpp"     // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/CoordinateSystems.hpp" // for Cartesian
#include "stk_mesh/base/CreateFaces.hpp"  // for handling faces
#include "stk_mesh/base/CreateEdges.hpp"  // for handling faces
#include "stk_mesh/base/Bucket.hpp"       // for buckets
#include "stk_mesh/base/Field.hpp"    // for coordinates
#include "stk_mesh/base/GetEntities.hpp"    // for coordinates
#include "stk_mesh/base/FieldParallel.hpp"  // for handling parallel fields

// MORIS project header files.
#include "cl_Cell.hpp" // CON/src
#include "typedefs.hpp" // COR/src
//#include "cl_Mat.hpp" // LNA/src
#include "algorithms.hpp"
#include "cl_Mesh_Enums.hpp" // MTK/src
#include "fn_sort.hpp"             // For use of sort function // LNA/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_Database.hpp" // MTK/src

// Class forward declarations.
namespace moris
{

// timestep global for exodus output
//double gStkTimeStep = 0.0;

/**
 *  Mtk wrapper class around Trilinos STK
 */
class STK_Implementation : public database
{
protected:
    // local to global entity communication lists (for parallel runs)
    Matrix< DDUMat >  mLocalToGlobalElemMap;   // Local to global element map
    Matrix< DDUMat >  mLocalToGlobalNodeMap;   // Local to global node map
    Matrix< DDUMat >  mLocalToGlobalEdgeMap;   // Local to global edge map
    Matrix< DDUMat >  mLocalToGlobalFaceMap;   // Local to global face map

    // Entity processor owners
    Matrix< DDUMat >  mElemMapToOwnerProc;   // node map to owner proc
    Matrix< DDUMat >  mNodeMapToOwnerProc;   // node map to owner proc
    Matrix< DDUMat >  mEdgeMapToOwnerProc;   // edge map to owner proc
    Matrix< DDUMat >  mFaceMapToOwnerProc;   // face map to owner proc

    // Entity processor shared. Note: Elements are owned by current processor
    Cell < Cell < uint > > mElemMapToSharingProcs;   // node map to sharing procs
    Cell < Cell < uint > > mNodeMapToSharingProcs;   // node map to sharing procs
    Cell < Cell < uint > > mEdgeMapToSharingProcs;   // edge map to sharing procs
    Cell < Cell < uint > > mFaceMapToSharingProcs;   // face map to sharing procs

    // Lists of sharing processors
    Matrix< DDUMat >  mElemSharedProcsList;   // node map to owner proc
    Matrix< DDUMat >  mNodeSharedProcsList;   // node map to owner proc
    Matrix< DDUMat >  mEdgeSharedProcsList;   // edge map to owner proc
    Matrix< DDUMat >  mFaceSharedProcsList;   // face map to owner proc

    std::map < uint, uint > mProcsSharedToIndex;

public:
      // member variables for mtk
      stk::io::StkMeshIoBroker*    mMeshReader;
      stk::mesh::MetaData*         mMtkMeshMetaData;
      stk::mesh::BulkData*         mMtkMeshBulkData;

      Matrix< DDUMat >  mEntProcOwner;

      bool mDataGeneratedMesh = false;
      bool mFieldInDataGiven  = false;
      uint mMaxNumFields = 20;
      uint mNumDims = 0;

      typedef stk::mesh::Field<double>    Field1Comp;
      typedef stk::mesh::Field<double,stk::mesh::Cartesian2d>  Field2Comp;
      typedef stk::mesh::Field<double,stk::mesh::Cartesian3d>  Field3Comp;
      typedef stk::mesh::Field<double,stk::mesh::FullTensor22> Field4Comp;
      typedef stk::mesh::Field<double,stk::mesh::FullTensor>   Field9Comp;

      std::vector<Field1Comp*> mField1CompVec;
      std::vector<Field2Comp*> mField2CompVec;
      std::vector<Field3Comp*> mField3CompVec;
      std::vector<Field4Comp*> mField4CompVec;
      std::vector<Field9Comp*> mField9CompVec;
      std::vector < bool > mSetRankFlags;   // Flags for user-defined node [0], side [1], and block [2] sets.
      std::vector < std::vector < std::string > > mSetNames;  // User-defined names for node [0], side [1], and block [2] sets.
    
      /**
      * mtk destructor.
      */
      ~STK_Implementation();

      /**
      * mtk constructor (mesh generated internally or obtained from an Exodus file )
      *
      * @param[in] aFileName  .................    String with mesh file name.
      */
      STK_Implementation(
              std::string    aFileName,
              MtkSetsInfo*   aSetsInfo );

      /**
       * Create stk mesh with information given by the user.
       *
       * @param[in] aFileName   .................   String with mesh file name.
       *
       */
      void
      build_mesh(
              std::string    aFileName,
              MtkSetsInfo*   aSetsInfo );

      /**
       * Create a MORIS mesh with information given by the user.
       *
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
      STK_Implementation(
              MtkMeshData   aMeshData );

      /**
       * Create a MORIS mesh with information given by the user.
       *
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
      void
      build_mesh(
              MtkMeshData   aMeshData );

      /**
       * Verifies that the information given by the user makes sense.
       *
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
    void
    check_and_update_input_data(
            MtkMeshData&   aMeshData );

    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    check_and_update_fields_data(
            MtkMeshData&   aMeshData );
    
    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    check_and_update_sets_data(
            MtkMeshData&   aMeshData );
    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    declare_mesh_fields(
            MtkMeshData   aMeshData );
    /*
     * Returns
     * @param[in]  aMeshData
     * @param[in]  aFieldsInfo
     */
    void
    internal_declare_mesh_field(
            MtkMeshData   aMeshData,
            uint          iField );
    /*
     * Returns
     * @param[in]  aNumModelDims
     * @param[in]  aNumNodesPerElem
     * @param[in]  aPartNames
     */
    void
    declare_mesh_parts(
            MtkMeshData   aMeshData );

    /*
     * Returns
     * @param[in]  aModelDim
     * @param[in]  aNumNodesInElem
     */
    stk::topology::topology_t
    get_mesh_topology(
            uint   aModelDim,
            uint   aNumNodesInElem );

    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    populate_mesh_database(
            MtkMeshData   aMeshData );
    
    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    process_block_sets(
            MtkMeshData   aMeshData );

    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    process_node_sets(
            MtkMeshData   aMeshData );

    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    process_side_sets(
            MtkMeshData   aMeshData );

    /*
     * Returns
     * @param[in]  aStkElemConn
     * @param[in]  aElemParts
     * @param[in]  aMeshData
     */
    void
    populate_mesh_database_serial(
            MtkMeshData                            aMeshData,
            std::vector< stk::mesh::PartVector >   aElemParts,
            Matrix< DDUMat >                             aOwnerPartInds);
    /*
     * Returns
     * @param[in]  aStkElemConn
     * @param[in]  aElemParts
     * @param[in]  aMeshData
     */
    void
    populate_mesh_database_parallel(
            MtkMeshData                            aMeshData,
            std::vector< stk::mesh::PartVector >   aElemParts,
            Matrix< DDUMat >                             aOwnerPartInds);

    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    populate_mesh_fields(
            MtkMeshData   aMeshData );
    /**
     * Accessor to member variable mMtkMeshMetaData.
     *
     * @return mMtkMeshMetaData.
     */
    stk::mesh::MetaData*
    get_mtk_mesh_meta_data( )
    {
        return mMtkMeshMetaData;
    }

    /**
     * Accessor to member variable mMtkMeshBulkData.
     *
     * @return mMtkMeshBulkData.
     */
    stk::mesh::BulkData*
    get_mtk_mesh_bulk_data( )
    {
        return mMtkMeshBulkData;
    }
    /*
     * Create an exodus mesh database with the specified
     * filename.
     *
     * @param[in] filename The full pathname to the file which will be
     *   created and the mesh data written to. If the file already
     *   exists, it will be overwritten.
     *
     *   Description from create_output_mesh() in StkMeshIoBroker.hpp
     */
    void
    create_output_mesh(
            std::string  &aFileName );

    Matrix< DDUMat >
    get_nodal_local_map( )
    {
        return mLocalToGlobalNodeMap;
    }

    Matrix< DDUMat >
    get_elemental_local_map( )
    {
        return mLocalToGlobalElemMap;
    }

    Matrix< DDUMat >
    get_edge_local_map( )
    {
        return mLocalToGlobalEdgeMap;
    }

    Matrix< DDUMat >
    get_face_local_map( )
    {
        return mLocalToGlobalFaceMap;
    }
    Matrix< DDUMat >
    get_nodal_owner_proc_map( )
    {
        return mNodeMapToOwnerProc;
    }

    Matrix< DDUMat >
    get_elemental_owner_proc_map( )
    {
        return mElemMapToOwnerProc;
    }

    Matrix< DDUMat >
    get_edge_owner_proc_map( )
    {
        return mEdgeMapToOwnerProc;
    }

    Matrix< DDUMat >
    get_face_owner_proc_map( )
    {
        return mFaceMapToOwnerProc;
    }

    Cell < Cell < uint > >
    get_nodes_shared_processors( )
    {
        return mNodeMapToSharingProcs;
    }

    Cell < Cell < uint > >
    get_elements_shared_processors( )
    {
        return mElemMapToSharingProcs;
    }

    Cell < Cell < uint > >
    get_edges_shared_processors( )
    {
        return mEdgeMapToSharingProcs;
    }

    Cell < Cell < uint > >
    get_faces_shared_processors( )
    {
        return mFaceMapToSharingProcs;
    }
    Matrix< DDUMat >
    get_node_ids_from_local_map(
            Matrix< DDUMat >    aLocalInds ) const;

    Matrix< DDUMat >
    get_element_ids_from_local_map_interface(
            Matrix< DDUMat >    aLocalInds ) const;

    uint
    get_num_entities_universal(
            enum EntityRank  aEntityRank ) const
    {
        return stk::mesh::count_selected_entities( mMtkMeshMetaData->universal_part(), mMtkMeshBulkData->buckets(this->get_stk_entity_rank(aEntityRank)) );
    }
    /**
     * Return number of elements in mesh.
     *
     * @return Number of elements.
     */
    uint
    get_num_elems( ) const
    {
        // Return number of elements extracted from entity counts vector
        return this->get_num_entities_universal( EntityRank::ELEMENT );
    }
    /**
     * Return number of faces in mesh.
     *
     * @return Number of faces.
     */
    uint
    get_num_faces( ) const
    {
        // Return number of faces extracted from entity counts vector
        return this->get_num_entities_universal( EntityRank::FACE );
    }

    /**
     * Return number of edges in mesh.
     *
     * @return Number of edges.
     */
    uint
    get_num_edges( ) const
    {
        // Return number of edges extracted from entity counts vector
        return this->get_num_entities_universal( EntityRank::EDGE );
    }

    /**
     * Return number of nodes in mesh.
     *
     * @return Number of nodes.
     */
    uint
    get_num_nodes( ) const
    {
        // Return number of nodes extracted from entity counts vector
        return this->get_num_entities_universal( EntityRank::NODE );
    }

    /**
     * Return number of nodes in mesh.
     *
     * @return Number of nodes.
     */
    uint
    get_num_spatial_dims( ) const
    {
        // Return number of nodes extracted from entity counts vector
        return mNumDims;
    }

    /**
     * Return number of nodes in mesh.
     *
     * @return Number of nodes.
     */
    EntityRank
    get_side_rank( ) const
    {
        // Return number of nodes extracted from entity counts vector
        return this->get_mtk_entity_rank( mMtkMeshMetaData->side_rank() );
    }

    /**
     * Return number of entities shared by all processors given an entity rank.
     *
     * @return Number of entities requested.
     */
    uint
    get_num_entities_globally_shared_interface(
            enum EntityRank aRequestedRank ) const
    {
        return stk::mesh::count_selected_entities(mMtkMeshMetaData->globally_shared_part(), mMtkMeshBulkData->buckets(this->get_stk_entity_rank(aRequestedRank)) );
    }

    uint
    get_num_entities_aura(
            enum EntityRank  aRequestedRank ) const
    {
        return stk::mesh::count_selected_entities( mMtkMeshMetaData->aura_part(), mMtkMeshBulkData->buckets(this->get_stk_entity_rank(aRequestedRank)) );
    }

    /**
     * Return number of elements in mesh.
     *
     * @return Number of elements.
     */
    uint
    get_num_elems_current_proc( ) const
    {
        // Return number of elements extracted from entity counts vector
        return this->get_num_entities_current_proc( EntityRank::ELEMENT );
    }

    /**
     * Return number of faces in mesh.
     *
     * @return Number of faces.
     */
    uint
    get_num_faces_current_proc( ) const
    {
        // Return number of faces extracted from entity counts vector
        return this->get_num_entities_current_proc( EntityRank::FACE );
    }

    /**
     * Return number of edges in mesh.
     *
     * @return Number of edges.
     */
    uint
    get_num_edges_current_proc( ) const
    {
        // Return number of edges extracted from entity counts vector
        return this->get_num_entities_current_proc( EntityRank::EDGE );
    }

    /**
     * Return number of nodes in mesh.
     *
     * @return Number of nodes.
     */
    uint
    get_num_nodes_current_proc( ) const
    {
        // Return number of nodes extracted from entity counts vector
        return this->get_num_entities_current_proc( EntityRank::NODE );
    }

    /**
     * Return number of entities in current processor given the appropriate rank.
     *
     * @return Number of entities requested.
     */
    uint
    get_num_entities_current_proc(
            enum EntityRank   aRequestedRank ) const
    {
        return this->get_num_entities_locally_owned_globally_shared( aRequestedRank );
    }

    /*
     * Return the number of a entities
     */
    uint
    get_num_entities_locally_owned_globally_shared(
            enum EntityRank   aRequestedRank ) const
    {
        stk::mesh::Selector aLocOwnedAndGlobSharedEntities = mMtkMeshMetaData->globally_shared_part() | mMtkMeshMetaData->locally_owned_part();
        return stk::mesh::count_selected_entities( aLocOwnedAndGlobSharedEntities, mMtkMeshBulkData->buckets(this->get_stk_entity_rank(aRequestedRank)) );
    }

    /*
     * Return the number of a entities
     */
    Matrix< DDUMat >
    get_entities_universal(
            EntityRank   aEntityRank ) const
    {
        return this->get_entities_in_selector_interface( aEntityRank, mMtkMeshMetaData->universal_part() );
    }
    /*
    * Returns the local entity index of entities globally shared by current process
    *
    */
    Matrix< DDUMat >
    get_entities_glb_shared_current_proc(
            EntityRank   aEntityRank ) const
    {
        return this->get_entities_in_selector_interface( aEntityRank, mMtkMeshMetaData->globally_shared_part() );
    }
    /*
    * Returns the local entity index of entities owned by current proc
    *  @return
    */
    Matrix< DDUMat >
    get_entities_owned_current_proc(
            EntityRank   aEntityRank ) const
    {
        return this->get_entities_in_selector_interface( aEntityRank, mMtkMeshMetaData->locally_owned_part() );
    }
    /*
    * Returns the local entity index of entities in aura part by current proc
    *  @return
    */
    Matrix< DDUMat >
    get_entities_in_aura(
            EntityRank   aEntityRank ) const
    {
        return this->get_entities_in_selector_interface( aEntityRank, mMtkMeshMetaData->aura_part() );
    }
    /*
    * Returns the local entity index of entities owned and shared by current proc
    *  @return
    */
    Matrix< DDUMat >
    get_entities_owned_and_shared_by_current_proc(
            EntityRank   aEntityRank ) const;
    /*
    * Returns the local entity index of entities owned and shared by current proc
    *  @return
    */
    Matrix< DDUMat >
    get_entities_in_selector_interface(
            EntityRank            aEntityRank,
            stk::mesh::Selector   aSelectedEntities ) const;

    /**
         * Return number of in element of a particular topology.
         *
         * @return Number of face.
    */
    uint
    get_elem_topology_num_faces(
            uint   aElemId ) const
    {
        return this->get_entity_topology_num_entities( stk::topology::ELEMENT_RANK, stk::topology::FACE_RANK, aElemId );
    }
    /**
         * Return number of in element of a particular topology.
         *
         * @return Number of face.
         */
    uint
    get_elem_topology_num_edges(
            uint   aElemId ) const
    {
        return this->get_entity_topology_num_entities( stk::topology::ELEMENT_RANK, stk::topology::EDGE_RANK, aElemId );
    }
    /**
         * Return number of in element of a particular topology.
         *
         * @return Number of face.
         */
    uint
    get_elem_topology_num_nodes(
            uint   aElemId ) const
    {
        return this->get_entity_topology_num_entities( stk::topology::ELEMENT_RANK, stk::topology::NODE_RANK, aElemId );
    }
    /**
         * Return number of in element of a particular topology.
         *
         * @return Number of face.
         */
    uint
    get_face_topology_num_edges(
            uint   aFaceId ) const
    {
        return this->get_entity_topology_num_entities( stk::topology::FACE_RANK, stk::topology::EDGE_RANK, aFaceId );
    }
    /**
         * Return number of in element of a particular topology.
         *
         * @return Number of face.
         */
    uint
    get_face_topology_num_nodes(
            uint   aFaceId ) const
    {
        return this->get_entity_topology_num_entities( stk::topology::FACE_RANK, stk::topology::NODE_RANK, aFaceId );
    }
    /**
         * Return number of in element of a particular topology.
         *
         * @return Number of face.
         */
    uint
    get_edge_topology_num_nodes(
            uint   aEdgeId ) const
    {
        return this->get_entity_topology_num_entities( stk::topology::EDGE_RANK, stk::topology::NODE_RANK, aEdgeId );
    }
    /**
     * Get elements connected to a face
     *
     * @param[in]  aFaceId    ............................   Entity Id
     * @param[out] aElementsConnectedToFace   ............   Connected entities
     *
     */

    uint
    get_entity_topology_num_entities(
            stk::mesh::EntityRank   aEntityRank,
            stk::mesh::EntityRank   aRequestedRank,
            uint                    aEntityId ) const;

    Matrix< DDUMat >
    get_elements_connected_to_element(
            uint   const aElemId ) const;

    Matrix< DDUMat >
    get_elements_connected_to_face(
            uint   const aFaceId ) const
    {
        const stk::mesh::Entity aFaceEntity = mMtkMeshBulkData->get_entity(stk::topology::FACE_RANK, aFaceId );
        return this->entities_connected_to_entity_generic( &aFaceEntity , stk::topology::ELEMENT_RANK );
    }
    /**
     * Get elements connected to a edge
     *
     * @param[in]  aEdgeId    ............................   Entity Id
     * @param[out] aElementsConnectedToEdge   ............   Connected entities
     *
     */
    Matrix< DDUMat >
    get_elements_connected_to_edge(
            uint   const aEdgeId ) const
    {
        const stk::mesh::Entity aEdgeEntity = mMtkMeshBulkData->get_entity(stk::topology::EDGE_RANK, aEdgeId );
        return this->entities_connected_to_entity_generic( &aEdgeEntity , stk::topology::ELEMENT_RANK );
    }
    /**
     * Get elements connected to a node
     *
     * @param[in]  aNodeId    ............................   Entity Id
     * @param[out] aElementsConnectedToNode   ............   Connected entities
     *
     */
    Matrix< DDUMat >
    get_elements_connected_to_node(
            uint   const aNodeId ) const
    {
        const stk::mesh::Entity aNodeEntity = mMtkMeshBulkData->get_entity(stk::topology::NODE_RANK, aNodeId );
        return this->entities_connected_to_entity_generic( &aNodeEntity , stk::topology::ELEMENT_RANK );
    }

    /**
     * Get faces connected to an element
     *
     * @param[in]  aElemId       ......................   Entity Id
     * @param[out] aFacesConnectedToElem   ............   Connected entities
     *
     */
    Matrix< DDUMat >
    get_faces_connected_to_element(
            uint   const aElemId  ) const
    {
        const stk::mesh::Entity aElemEntity = mMtkMeshBulkData->get_entity(stk::topology::ELEMENT_RANK, aElemId );
        return this->entities_connected_to_entity_generic( &aElemEntity , stk::topology::FACE_RANK );
    }
    /**
     * Get faces connected to an edge
     *
     * @param[in]  aEdgeId       ......................   Entity Id
     * @param[out] aFacesConnectedToEdge   ............   Connected entities
     *
     */
    Matrix< DDUMat >
    get_faces_connected_to_edge(
            uint   const aEdgeId ) const
    {
        const stk::mesh::Entity aEdgeEntity = mMtkMeshBulkData->get_entity(stk::topology::EDGE_RANK, aEdgeId );
        return this->entities_connected_to_entity_generic( &aEdgeEntity , stk::topology::FACE_RANK );
    }
    /**
     * Get faces connected to a node
     *
     * @param[in]  aNodeId       ......................   Entity Id
     * @param[out] aFacesConnectedToNode   ............   Connected entities
     *
     */
    Matrix< DDUMat >
    get_faces_connected_to_node(
            uint   const aNodeId ) const
    {
        const stk::mesh::Entity aNodeEntity = mMtkMeshBulkData->get_entity(stk::topology::NODE_RANK, aNodeId );
        return this->entities_connected_to_entity_generic( &aNodeEntity , stk::topology::FACE_RANK );
    }

    /**
     * Get edges connected to an element
     *
     * @param[in]  aElemId       ......................   Entity Id
     * @param[out] aEdgesConnectedToElem   ............   Connected entities
     *
     */
    Matrix< DDUMat >
    get_edges_connected_to_element(
            uint   const aElemId ) const
    {
        const stk::mesh::Entity aElemEntity = mMtkMeshBulkData->get_entity(stk::topology::ELEMENT_RANK, aElemId );
        return this->entities_connected_to_entity_generic( &aElemEntity , stk::topology::EDGE_RANK );
    }
    /**
     * Get edges connected to a face
     *
     * @param[in]  aFaceId       ......................   Entity Id
     * @param[out] aEdgesConnectedToFace   ............   Connected entities
     *
     */
    Matrix< DDUMat >
    get_edges_connected_to_face(
            uint   const aFaceId ) const
    {
        const stk::mesh::Entity aFaceEntity = mMtkMeshBulkData->get_entity(stk::topology::FACE_RANK, aFaceId );
        return this->entities_connected_to_entity_generic( &aFaceEntity , stk::topology::EDGE_RANK );
    }
    /**
     * Get edges connected to a node
     *
     * @param[in]  aNodeId       ......................   Entity Id
     * @param[out] aEdgesConnectedToNode   ............   Connected entities
     *
     */
    Matrix< DDUMat >
    get_edges_connected_to_node(
            uint   const aNodeId ) const
    {
        const stk::mesh::Entity aNodeEntity = mMtkMeshBulkData->get_entity(stk::topology::NODE_RANK, aNodeId );
        return this->entities_connected_to_entity_generic( &aNodeEntity , stk::topology::EDGE_RANK );
    }

    /**
     * Get nodes connected to an element
     *
     * @param[in]  aElemId       ......................   Entity Id
     * @param[out] aNodesConnectedToElem   ............   Connected entities
     *
     */
    Matrix< DDUMat >
    get_nodes_connected_to_element(
            uint   const aElemId ) const
    {
        const stk::mesh::Entity aElemEntity = mMtkMeshBulkData->get_entity(stk::topology::ELEMENT_RANK, aElemId );
        return this->entities_connected_to_entity_generic( &aElemEntity , stk::topology::NODE_RANK );
    }
    /**
     * Get nodes connected to a face
     *
     * @param[in]  aFaceId       ......................   Entity Id
     * @param[out] aNodesConnectedToFace   ............   Connected entities
     *
     */
    Matrix< DDUMat >
    get_nodes_connected_to_face(
            uint   const aFaceId ) const
    {
        const stk::mesh::Entity aFaceEntity = mMtkMeshBulkData->get_entity(stk::topology::FACE_RANK, aFaceId );
        return this->entities_connected_to_entity_generic( &aFaceEntity , stk::topology::NODE_RANK );
    }
    /**
     * Get nodes connected to an edge
     *
     * @param[in]  aNodeId       ......................   Entity Id
     * @param[out] aNodesConnectedToEdge   ............   Connected entities
     *
     */
    Matrix< DDUMat >
    get_nodes_connected_to_edge(
            uint   const aEdgeId ) const
    {
        const stk::mesh::Entity aEdgeEntity = mMtkMeshBulkData->get_entity(stk::topology::EDGE_RANK, aEdgeId );
        return this->entities_connected_to_entity_generic( &aEdgeEntity , stk::topology::NODE_RANK );
    }

    /**
     * Get upward connectivity (general)
     *
     * @param[in]  entity given as input   ............................   Id of the reference entity
     * @param[in]  rank of input entity    ............................
     * @param[in]  rank of output entity   ............................
     * @param[out] entities connected to given entity      ............   Ids of entities of specified rank to given entity
     *
     */
    Matrix< DDUMat >
    entities_connected_to_given_entity(
            uint         const aEntityId,
            EntityRank   const aInputEntityRank,
            EntityRank   const aOutputEntityRank ) const
    {
        const stk::mesh::Entity aStkEntity = mMtkMeshBulkData->get_entity( this->get_stk_entity_rank( aInputEntityRank), aEntityId );
        return this->entities_connected_to_entity_generic( &aStkEntity , this->get_stk_entity_rank( aOutputEntityRank) );
    }

    /**
     * Get upward connectivity (general)
     *
     * @param[in]  entity given as input   ............................   Id of the reference entity
     * @param[out] entities connected to given entity      ............   Ids of entities of specified rank to given entity
     *
     */
    Matrix< DDUMat >
    entities_connected_to_entity_generic(
            const stk::mesh::Entity*      aInputEntity,
            stk::mesh::EntityRank   needConnectivityOfType ) const;
    /**
     * Get upward connectivity (general)
     *
     * @param[in]  entity given as input   ............................   Id of the reference entity
     * @param[in]  rank of input entity    ............................
     * @param[in]  rank of output entity   ............................
     * @param[out] entities connected to given entity      ............   Ids of entities of specified rank to given entity
     *
     */
    Matrix< DDUMat >
    entity_ordinals_connected_to_given_entity(
            uint         const aEntityId,
            EntityRank   const aInputEntityRank,
            EntityRank   const aOutputEntityRank ) const;

    Matrix< DDUMat >
    entity_ordinals_connected_to_entity_generic(
            const stk::mesh::Entity*   aInputEntity ,
            stk::mesh::EntityRank      aNeedConnectivityOfType ) const;

    Matrix< DDUMat >
    get_nodes_in_node_set(
            uint    const aNodeSetId ) const
    {
        return this->get_entities_in_set_generic( stk::topology::NODE_RANK, stk::topology::NODE_RANK, aNodeSetId );
    }

    Matrix< DDUMat >
    get_nodes_in_side_set(
            uint    const aSideSetId ) const
    {
        return this->get_entities_in_set_generic( stk::topology::NODE_RANK, mMtkMeshMetaData->side_rank(), aSideSetId );
    }

    Matrix< DDUMat >
    get_nodes_in_block_set(
            uint    const aBlockSetId ) const
    {
        return this->get_entities_in_set_generic( stk::topology::NODE_RANK, stk::topology::ELEMENT_RANK, aBlockSetId );
    }

    Matrix< DDUMat >
    get_edges_in_side_set(
            uint    const aSideSetId ) const
    {
        return this->get_entities_in_set_generic( stk::topology::EDGE_RANK, mMtkMeshMetaData->side_rank(), aSideSetId );
    }

    Matrix< DDUMat >
    get_edges_in_block_set(
            uint    const aBlockSetId ) const
    {
        return this->get_entities_in_set_generic( stk::topology::EDGE_RANK, stk::topology::ELEMENT_RANK, aBlockSetId );
    }

    Matrix< DDUMat >
    get_faces_in_side_set(
            uint    const aSideSetId ) const
    {
        return this->get_entities_in_set_generic( stk::topology::FACE_RANK, mMtkMeshMetaData->side_rank(), aSideSetId );
    }

    Matrix< DDUMat >
    get_faces_in_block_set(
            uint    const aBlockSetId ) const
    {
        return this->get_entities_in_set_generic( stk::topology::FACE_RANK, stk::topology::ELEMENT_RANK, aBlockSetId );
    }

    Matrix< DDUMat >
    get_entities_in_set_generic(
            stk::mesh::EntityRank   aRequiredEntityRank,
            stk::mesh::EntityRank   aSetRank,
            uint                    const aNodeSetId ) const;

    Matrix< DDUMat >
    get_nodes_in_set(
            std::string    const aSetName ) const
    {
        return this->get_set_entity_ids( stk::topology::NODE_RANK, aSetName );
    }

    Matrix< DDUMat >
    get_edges_in_set(
            std::string    const aSetName ) const
    {
        return this->get_set_entity_ids( stk::topology::EDGE_RANK, aSetName );
    }

    Matrix< DDUMat >
    get_faces_in_set(
            std::string    const aSetName ) const
    {
        return this->get_set_entity_ids( stk::topology::FACE_RANK, aSetName );
    }

    Matrix< DDUMat >
    get_elements_in_set(
            std::string    const aSetName ) const
    {
        return this->get_set_entity_ids( stk::topology::ELEMENT_RANK, aSetName );
    }
    /*
     * Returns
     * @param[in]  aEntityRank
     * @param[in]  aFieldName
     * @param[out] field values
     */
    Matrix< DDUMat >
    get_set_entity_ids(
            stk::mesh::EntityRank   aEntityRank,
            std::string             aFieldName ) const;

    /*
     * Returns
     * @param[in]  aEntityRank
     * @param[in]  aFieldName
     * @param[out] field values
     */
    Matrix< DDUMat >
    get_intersected_entities_field_set(
            enum EntityRank   aEntityRank,
            std::string        aFieldName,
            std::string        aSetName ) const;

    Matrix< DDRMat >
    get_intersected_data_field_set(
            enum EntityRank   aEntityRank,
            std::string        aFieldName,
            std::string        aSetName ) const;
    
    /*
     * Returns
     * @param[in]  aEntityRank
     * @param[in]  aFieldName
     * @param[out] field values
     */
    Matrix< DDUMat >
    get_set_entity_ids(
            enum EntityRank   aEntityRank,
            std::string       aFieldName ) const
    {
        return this->get_set_entity_ids( this->get_stk_entity_rank(aEntityRank), aFieldName );
    }

    uint
    get_entity_index(
            enum EntityRank   aEntityRank,
            uint   const aEntityId ) const
    {
        return mMtkMeshBulkData->local_id( mMtkMeshBulkData->get_entity( get_stk_entity_rank(aEntityRank), aEntityId ) );
    }
    /**
     * Get upward connectivity (general)
     *
     * @param[in]  entity given as input   ............................   Id of the reference entity
     * @param[in]  rank of input entity    ............................
     * @param[in]  rank of output entity   ............................
     * @param[out] entities connected to given entity      ............   Ids of entities of specified rank to given entity
     *
     */
    Matrix< DDUMat >
    get_entity_local_ids_connected_to_entity(
            uint         const aEntityId ,
            EntityRank   const aInputEntityRank,
            EntityRank   const aOutputEntityRank ) const;

    /**
         * Return coordinates of the node with the corresponding Ids.
         *
         * @return coordinates.
         */
    Matrix< DDRMat >
    get_selected_nodes_coords(
            Matrix< DDUMat >  aNodeIds ) const
    {
        return this->get_nodes_coords_generic( aNodeIds );
    }

    Matrix< DDRMat >
    get_selected_nodes_coords_lcl_ind(
            Matrix< DDUMat >  aNodeInds ) const;
    /**
         * Return coordinates of the node with the corresponding Ids.
         *
         * @return coordinates.
         */
    Matrix< DDRMat >
    get_all_nodes_coords( ) const
    {
        Matrix< DDUMat >  aUniversalNodeIds = this->get_entities_universal( EntityRank::NODE );
        return this->get_nodes_coords_generic( aUniversalNodeIds );
    }
    /*
     * Returns coordinates of nodes found in aura
     */
    Matrix< DDRMat >
    get_all_nodes_coords_aura( ) const
    {
        Matrix< DDUMat >  aAuraNodeIds = this->get_entities_in_aura( EntityRank::NODE );
        return this->get_nodes_coords_generic( aAuraNodeIds );
    }

    Matrix< DDRMat >
    get_nodes_coords_generic(
            Matrix< DDUMat >  aNodeIds ) const;

    /*
     * Returns a list of globally unique element ids
     * @param[in]  aNumNodes - number of element ids requested
     * @param[out] aAvailableNodeIDs - list of globally unique element IDs
     */
    Matrix< DDUMat >
    generate_unique_elem_ids(
            uint   aNumEntities ) const
    {
        return this->generate_unique_entity_ids( aNumEntities, stk::topology::ELEMENT_RANK );
    }
    /*
     * Returns a list of globally unique face ids
     * @param[in]  aNumNodes - number of face ids requested
     * @param[out] aAvailableNodeIDs - list of globally unique face IDs
     */
    Matrix< DDUMat >
    generate_unique_face_ids(
            uint   aNumEntities ) const
    {
        return this->generate_unique_entity_ids( aNumEntities, stk::topology::FACE_RANK );
    }

    /*
     * Returns a list of globally unique edge ids
     * @param[in]  aNumNodes - number of edge ids requested
     * @param[out] aAvailableNodeIDs - list of globally unique edge IDs
     */
    Matrix< DDUMat >
    generate_unique_edge_ids(
            uint   aNumEntities ) const
    {
        return this->generate_unique_entity_ids( aNumEntities, stk::topology::EDGE_RANK );
    }

    /*
     * Returns a list of globally unique node ids
     * @param[in]  aNumNodes - number of node ids requested
     * @param[out] aAvailableNodeIDs - list of globally unique node IDs
     */
    Matrix< DDUMat >
    generate_unique_node_ids(
            uint   aNumEntities ) const
    {
        return this->generate_unique_entity_ids( aNumEntities, stk::topology::NODE_RANK );
    }

    /*
     * Returns a list of globally unique entity ids
     * @param[in]  aNumNodes - number of entity ids requested
     * @param[out] aAvailableNodeIDs - list of globally unique entity IDs
     */
    Matrix< DDUMat >
    generate_unique_entity_ids(
            uint                    aNumEntities,
            stk::mesh::EntityRank   aNeedConnectivityOfType ) const;

    /*
     * Prints information in parallel runs (for debugging)
     */
    void add_shared_field_interface( ) const;
    /*
     * Interface between stk entity rank and MORIS enum
     */
    stk::mesh::EntityRank
    get_stk_entity_rank(
            enum EntityRank   aEntityRank ) const ;

    /*
     * Interface between stk entity rank and MORIS enum
     */
    EntityRank
    get_mtk_entity_rank(
            stk::mesh::EntityRank   aEntityRank ) const ;

    /*
     * Returns
     * @param[in]  aEntityIndex
     * @param[in]  aEntityRank
     * @param[out] rank of owner processor
     */
    uint
    parallel_owner_rank_by_entity_id(
            uint            aEntityIndex,
            enum EntityRank aEntityRank) const;

    /*
     * Returns
     * @param[in]  aEntityIndex
     * @param[in]  aEntityRank
     * @param[out] rank of owner processor
     */
    uint
    parallel_owner_rank_by_entity_index(
            uint            aEntityIndex,
            enum EntityRank aEntityRank) const;

    /*
     * Returns
     * @param[in]  aEntityRank       - entity rank
     * @param[in]  aEntityRank       - entity ID
     * @param[out] list of entities shared
     */
    Matrix< DDUMat >
    get_procs_sharing_entity_by_id(
            uint              aEntityID,
            enum EntityRank   aEntityRank ) const;

    /*
     * Returns
     * @param[in]  aEntityIndex
     * @param[in]  aEntityRank
     * @param[out] list of shared processors
     */
    Matrix< DDUMat >
    get_procs_sharing_entity_by_index(
            uint              aEntityIndex,
            enum EntityRank   aEntityRank ) const;

    /*
     * Returns
     * @param[out] number of processors
     */
    uint
    get_parallel_size() const
    {
        return mMtkMeshBulkData->parallel_size();
    }

    /*
     * Returns
     * @param[out] current processor rank
     */
    uint
    get_parallel_rank() const
    {
        return mMtkMeshBulkData->parallel_rank();
    }

    /*
     * Returns
     * @param[in]  aEntityRank
     * @param[in]  aFieldName
     * @param[out] field values
     */
    Matrix< DDUMat >
    get_field_entities(
            enum EntityRank   aEntityRank,
            std::string       aFieldName )
    {
        stk::mesh::FieldBase const * aField_base = mMtkMeshMetaData->get_field( this->get_stk_entity_rank( aEntityRank ), aFieldName );
        return this->get_field_entities( aField_base );
    }

    Matrix< DDUMat >
    get_field_entities(
            stk::mesh::FieldBase   const * aField ) ;

    Matrix< DDRMat >
    get_field_values(
            enum EntityRank   aEntityRank,
            std::string       aFieldName) ;

    Matrix< DDUMat >
    get_part_entities(
            enum EntityRank   aEntityRank,
            std::string       aPartName ) ;

    /*
     * Returns
     * @param[out] entity id
     */
    uint
    get_entity_key_from_entity_id(
            enum EntityRank   aEntityRank,
            uint              aEntityIdNumber ) const
    {
        stk::mesh::Entity aEnitity = mMtkMeshBulkData->get_entity(this->get_stk_entity_rank( aEntityRank), aEntityIdNumber );
        return stk::mesh::EntityKey::entity_key_t( mMtkMeshBulkData->entity_key( aEnitity ));
    }

    /*
     * Returns
     * @param[out] entity id
     */
    uint
    get_entity_id_from_entity_key(
            enum EntityRank   aEntityRank,
            uint              aEntityKeyNumber ) const
    {
        stk::mesh::EntityKey aEntKey( static_cast<stk::mesh::EntityKey::entity_key_t>(aEntityKeyNumber) );
        return (uint) aEntKey.id();
    }

    /*
     *
     * Create communication list (only in parallel runs)
     */
    void
    create_communication_lists_from_string();

    void
    create_additional_communication_lists_from_data();

    void
    create_facets_communication_lists();

    void
    create_owners_communication_lists();

    void
    create_shared_communication_lists();

    Cell < Cell < uint > >
    get_shared_info_by_entity( uint aNumActiveSharedProcs, enum EntityRank  aEntityRank );
    /*
     *
     * Print mesh info
     */
    void
    print_parts_to_screen();

    void
    print_fields_to_screen();

    /**
     * Get duplicates of a coordinate and id list
     *
     * @param[in]  aCoord         .... Coordinate list with x,y,z
     * @param[in]  aId            .... Id list of the coordinates
     * @param[out] duplicate_list .... Shows the duplicates [Position(i) Position(j)]
     *
     */
    Matrix< DDUMat >
    duplicate_node_coord_check();

    Matrix< DDUMat >
    duplicate_node_coord_check(
            Matrix< DDRMat > &   aCoord );

    Matrix< DDUMat >
    duplicate_node_coord_and_id_check(
            Matrix< DDRMat > &   aCoord,
            Matrix< DDUMat > &   aId );

    Matrix< DDUMat >
    duplicate_node_coord_and_id_check(
            Cell< Matrix< DDRMat >  >&   aCoord,
            Cell< Matrix< DDUMat >  >&   aId );

    Matrix< DDUMat >
    duplicate_node_coord_and_id_check_problems(
            Matrix< DDRMat > &   aCoord,
            Matrix< DDUMat > &   aId );

    Matrix< DDUMat >
    duplicate_node_coord_and_id_check_problems(
            Cell< Matrix< DDRMat >  >&   aCoord,
            Cell< Matrix< DDUMat >  >&   aId );

    //////////////////////////
    // Deprecated functions //
    //////////////////////////

    Matrix< DDRMat >
    interpolate_to_location_on_entity(
            enum EntityRank   aParentEntityRank,
            uint              aParentEntityIndex,
            Matrix< DDRMat >          aLclCoord);
};

}   // namespace moris
#endif /* MORIS_MESH_CL_MESH_MTK_HPP_ */

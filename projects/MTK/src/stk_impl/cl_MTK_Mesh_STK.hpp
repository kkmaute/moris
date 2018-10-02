/*
 * cl_MTK_Mesh_STK_New.hpp
 *
 *  Created on: Sep 18, 2018
 *      Author: barrera/doble
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_MESH_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_MESH_STK_HPP_

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


#include "../cl_MTK_Sets_Info.hpp"
#include "../cl_MTK_Mesh_Data_Input.hpp"
#include "../cl_MTK_Block.hpp"
#include "../cl_Mesh_Enums.hpp"
#include "../cl_MTK_Mesh.hpp"

#include "cl_Cell.hpp"

// For moris::Cell and vertex APi
#include "cl_MTK_Cell_STK.hpp"
#include "cl_MTK_Vertex_STK.hpp"


namespace moris
{
namespace mtk
{
class Mesh_STK: public Mesh
{
    //! timestamp for stk output. Set in cosnstructor over MtkMeshData
    double mTimeStamp = 0.0;

public:
    //##############################################
    // Build mesh functionality
    //##############################################
    /**
     * STK destructor.
     */
    ~Mesh_STK();

    /**
     * STK constructor (mesh generated internally or obtained from an Exodus file )
     *
     * @param[in] aFileName  .................    String with mesh file name.
     */
    Mesh_STK(
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
    Mesh_STK(
            MtkMeshData   aMeshData );


    //##############################################
    // General mesh information access
    //##############################################

    /*
     * Get spatial dimension of the mesh
     */
    uint
    get_spatial_dim() const
    {
        return (uint)mMtkMeshMetaData->spatial_dimension();
    }

    /*
     * Get number of entities for specified rank in the STK universal part
     */
    uint
    get_num_entities(
            enum EntityRank aEntityRank) const;


    //##############################################
    // Access Mesh Data by index Functions
    //##############################################
    /*
     * Generic get local index of entities connected to
     * entity using an entities local index
     */
    Matrix<IndexMat>
    get_entity_connected_to_entity_loc_inds(moris_index     aEntityIndex,
                                            enum EntityRank aInputEntityRank,
                                            enum EntityRank aOutputEntityRank) const;

    /*
     * Since the connectivity between entities of the same rank are considered
     * invalid by STK standards, we need a seperate function for element to element
     * specifically
     *      *
     * @param[in]  aElementId - element id
     * @param[out] Element to element connectivity and face ordinal shared
     *                   (where elements are all by index)
     */
    Matrix< IndexMat >
    get_elements_connected_to_element_loc_inds(moris_index aElementIndex) const;




    //##############################################
    // global id functions
    //##############################################

    /*
     * Get global identifier of an entity from a local index and entity rank
     */
    moris_id
    get_glb_entity_id_from_entity_loc_index(moris_index     aEntityIndex,
                                            enum EntityRank aEntityRank) const;

    /*
     * Get local indec of an entity from a global index
     */
    moris_index
    get_loc_entity_ind_from_entity_glb_id(moris_id        aEntityId,
                                          enum EntityRank aEntityRank) const;

    /*
     * Generic get global id of entities connected to
     * entity using an entities global id
     */
    Matrix< IdMat >
    get_entity_connected_to_entity_glob_ids( moris_id     aEntityId,
                                             enum EntityRank aInputEntityRank,
                                             enum EntityRank aOutputEntityRank) const;

    /*
     * Since the connectivity between entities of the same rank are considered
     * invalid by STK standards, we need a seperate function for element to element
     * specifically
     *
     * @param[in]  aElementId - element id
     * @param[out] Element to element connectivity and face ordinal shared
     */
    Matrix< IdMat >
    get_element_connected_to_element_glob_ids(moris_index aElementId) const;

    /*
     * Returns a list of globally unique entity ids for entities
     * of the provided rank
     * @param[in]  aNumNodes - number of node ids requested
     * @param[in]  aEntityRank - Entity rank to assign ids for
     * @param[out] aAvailableNodeIDs - list of globally unique node IDs
     */
    Matrix< IdMat >
    generate_unique_entity_ids( uint            aNumEntities,
                                enum EntityRank aEntityRank) const;

    //##############################################
    // Coordinate Field Functions
    //##############################################
    Matrix< DDRMat >
    get_node_coordinate( moris_index aNodeIndex ) const;


    //##############################################
    // Entity Ownership Functions
    //##############################################
    moris_id
    get_entity_owner(  moris_index     aEntityIndex,
                       enum EntityRank aEntityRank ) const;


    //##############################################
    // moris::Cell and Vertex Pointer Functions
    //##############################################

    /*
     * Get an mtk cell by index
     */
    const mtk::Cell &
    get_mtk_cell(moris_index aCellIndex);

    /*
     * get an mtk vertex by index
     */
    const mtk::Vertex &
    get_mtk_vertex(moris_index aVertexIndex);

    //------------------------------------------------------------------------------

    //fixme: this function needs to go
    /**
     * populates the member variables of the relevant nodes
     * with their T-Matrices
     */
    void
    finalize()
    {
        MORIS_ERROR(0,"Not implemented in STK");
    }

    //------------------------------------------------------------------------------

    /**
     * provides a moris::Mat<uint> containing the IDs this mesh has
     * to communicate with
     */
    Matrix< IdMat >
    get_communication_table() const
    {
        MORIS_ERROR(0,"Not implemented in STK");
        return mEntityLocaltoGlobalMap(0);
    }

    //##############################################
    //  Output Mesh To a File
    //##############################################
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



private:

    // STK specific Member variables
    stk::io::StkMeshIoBroker*    mMeshReader;
    stk::mesh::MetaData*         mMtkMeshMetaData;
    stk::mesh::BulkData*         mMtkMeshBulkData;

    // General mesh trait member variables
    bool mDataGeneratedMesh = false;

    // Local to Global c
    // moris::Cell(0) - Node Local to Global
    // moris::Cell(1) - Edge Local to Global
    // moris::Cell(2) - Face Local to Global
    // moris::Cell(3) - Edge Local to Global
    moris::Cell<moris::Matrix< IdMat >> mEntityLocaltoGlobalMap;

    //Exterior moris::Cell Entity Rank (Same structure as local to global
    // Interior moris::Cell (processor rank)
    moris::Cell<moris::Cell<moris::Matrix< IndexMat >>> mEntitySendList;
    moris::Cell<moris::Cell<moris::Matrix< IndexMat >>> mEntityReceiveList;

    //    // moris::Cell and Vertex
    moris::Cell<mtk::Cell_STK>   mMtkCells;
    moris::Cell<mtk::Vertex_STK> mMtkVertices;

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
    // TODO: ADD edge sets (not sure why these were neglected in previous implementation
    std::vector < std::vector < std::string > > mSetNames;  // User-defined names for node [0], side [1], and block [2] sets.

    // Entity processor shared. Note: Elements are owned by current processor
    moris::Cell < moris::Cell < uint > > mElemMapToSharingProcs;   // node map to sharing procs
    moris::Cell < moris::Cell < uint > > mNodeMapToSharingProcs;   // node map to sharing procs
    moris::Cell < moris::Cell < uint > > mEdgeMapToSharingProcs;   // edge map to sharing procs
    moris::Cell < moris::Cell < uint > > mFaceMapToSharingProcs;   // face map to sharing procs

    // Entity processor owners
    Matrix< DDUMat >  mElemMapToOwnerProc;   // node map to owner proc
    Matrix< DDUMat >  mNodeMapToOwnerProc;   // node map to owner proc
    Matrix< DDUMat >  mEdgeMapToOwnerProc;   // edge map to owner proc
    Matrix< DDUMat >  mFaceMapToOwnerProc;   // face map to owner proc

    std::map < uint, uint > mProcsSharedToIndex;

    // Dummy Block
    Block* mDummyBlock;

    //##############################################
    // Private functions to build mesh
    //##############################################
    /*
     *
     * Create communication list
     */
    void
    create_communication_lists_and_local_to_global_map(enum EntityRank aEntityRank);
    //------------------------------------------------------------------------------

    /*
     * sets up the vertex and cell api
     */
    void
    set_up_vertices_and_cell();
    //------------------------------------------------------------------------------


    //##############################################
    // Private functions to access mesh information
    //##############################################

    /*
     * Return the STK rank of an entity given the moris
     * entity rank
     */
    stk::mesh::EntityRank
    get_stk_entity_rank(enum EntityRank aMRSEntityRank) const;
    //------------------------------------------------------------------------------

    /*
    * Returns the local entity index of entities owned and shared by current proc
    *  @return
    */
    Matrix< IdMat >
    get_entities_owned_and_shared_by_current_proc(
            EntityRank   aEntityRank ) const;
    //------------------------------------------------------------------------------

    /*
     * Return the number of a entities
     */
    Matrix< DDUMat >
    get_entities_universal(
            EntityRank   aEntityRank ) const
    {
        return this->get_entities_in_selector_interface( aEntityRank, mMtkMeshMetaData->universal_part() );
    }
    //------------------------------------------------------------------------------

    /*
    * Returns the local entity index of entities owned and shared by current proc
    *  @return
    */
    Matrix< DDUMat >
    get_entities_in_selector_interface(
            EntityRank            aEntityRank,
            stk::mesh::Selector   aSelectedEntities ) const;
    //------------------------------------------------------------------------------

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
    //------------------------------------------------------------------------------
    moris::Cell < moris::Cell < uint > >
    get_shared_info_by_entity( uint aNumActiveSharedProcs, enum EntityRank  aEntityRank );
    //------------------------------------------------------------------------------


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

    //##############################################
    // internal id functions
    //##############################################
    std::vector<stk::mesh::Entity>
    entities_connected_to_entity_stk(stk::mesh::Entity*    const aInputEntity,
                                     const stk::mesh::EntityRank aInputEntityRank,
                                     const stk::mesh::EntityRank aOutputEntityRank) const;


    //##############################################
    // Build mesh from data functions internal
    //##############################################
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
     * @param[in]  aMeshData
     */
    void
    check_and_update_fields_data( MtkMeshData&   aMeshData );

    /*
     * @param[in]  aMeshData
     */
    void
    check_and_update_sets_data(
            MtkMeshData&   aMeshData );


    //##############################################
    // internal build functions
    //##############################################
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

//------------------------------------------------------------------------------

    /*
     * Returns
     * @param[in]  aNumModelDims
     * @param[in]  aNumNodesPerElem
     * @param[in]  aPartNames
     */
    void
    declare_mesh_parts(
            MtkMeshData   aMeshData );

//------------------------------------------------------------------------------

    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    declare_mesh_fields(
            MtkMeshData   aMeshData );

//------------------------------------------------------------------------------

    /*
     * Returns
     * @param[in]  aMeshData
     * @param[in]  aFieldsInfo
     */
    void
    internal_declare_mesh_field(
            MtkMeshData   aMeshData,
            uint          iField );

//------------------------------------------------------------------------------

    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    populate_mesh_database(
            MtkMeshData   aMeshData );
//------------------------------------------------------------------------------
    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    process_block_sets(
            MtkMeshData   aMeshData );
//------------------------------------------------------------------------------
    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    populate_mesh_fields(
            MtkMeshData   aMeshData );
//------------------------------------------------------------------------------
    /*
     * Returns
     * @param[in]  aModelDim
     * @param[in]  aNumNodesInElem
     */
    stk::topology::topology_t
    get_mesh_topology(
            uint   aModelDim,
            uint   aNumNodesInElem );
//------------------------------------------------------------------------------
    void
    create_additional_communication_lists_from_data();
//------------------------------------------------------------------------------
    void
    create_facets_communication_lists();
//------------------------------------------------------------------------------
    void
    create_owners_communication_lists();
//------------------------------------------------------------------------------
    void
    create_shared_communication_lists();
//------------------------------------------------------------------------------
    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    process_node_sets(
            MtkMeshData   aMeshData );
//------------------------------------------------------------------------------
    /*
     * Returns
     * @param[in]  aMeshData
     */
    void
    process_side_sets(
            MtkMeshData   aMeshData );
//------------------------------------------------------------------------------
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
            Matrix< DDUMat >                       aOwnerPartInds);
//------------------------------------------------------------------------------
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
            Matrix< DDUMat >                       aOwnerPartInds);
//------------------------------------------------------------------------------
    /*
     * Returns
     * @param[in]  aEntityRank
     * @param[in]  aFieldName
     * @param[out] field values
     */
    Matrix< IdMat >
    get_set_entity_ids(
            stk::mesh::EntityRank   aEntityRank,
            std::string             aFieldName ) const;
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

};
}
}


#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_MESH_STK_HPP_ */

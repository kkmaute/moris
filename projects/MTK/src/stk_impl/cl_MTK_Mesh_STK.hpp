/*
 * cl_MTK_Mesh_STK_New.hpp
 *
 *  Created on: Sep 18, 2018
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_MESH_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_MESH_STK_HPP_

#include <stk_io/StkMeshIoBroker.hpp>     // for StkMeshIoBroker
#include <stk_mesh/base/MetaData.hpp>     // for MetaData
#include <stk_mesh/base/BulkData.hpp>     // for BulkData

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
    get_element_connected_to_element_loc_inds(moris_index aElementIndex) const;




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
     * Generic get global id of entities connected to
     * entity using an entities global id
     */
    Matrix< IdMat >
    get_entity_connected_to_entity_glob_ids( moris_id     aEntityId,
                                            enum EntityRank aInputEntityRank,
                                            enum EntityRank aOutputEntityRank);

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
    mtk::Cell const &
    get_mtk_cell(moris_index aCellIndex);

    /*
     * get an mtk vertex by index
     */
    mtk::Vertex const &
    get_mtk_vertex(moris_index aVertexIndex);

    //##############################################
    //  Access block information
    //##############################################


   /**
    * returns the number of blocks on this mesh
    */
   uint
   get_number_of_blocks() const
   {
       MORIS_ERROR(0,"Not implemented in STK");
       return 0;
   }

//------------------------------------------------------------------------------

   /**
    * returns a pointer to a block
    */
   Block *
   get_block_by_index( const moris_index & aIndex )
   {
       MORIS_ERROR(0,"Not implemented in STK");
       return mDummyBlock;
   }

//------------------------------------------------------------------------------

   /**
    * returns a pointer to a block ( const version )
    */
   const Block *
   get_block_by_index( const moris_index & aIndex ) const
   {
       MORIS_ERROR(0,"Not implemented in STK");
       return mDummyBlock;
   }

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


private:
    // Set names
    // TODO: ADD edge sets (not sure why these were neglected in previous implementation
    std::vector < std::vector < std::string > > mSetNames;  // User-defined names for node [0], side [1], and block [2] sets.

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

    /*
     * sets up the vertex and cell api
     */
    void
    set_up_vertices_and_cell();

    //##############################################
    // Private functions to access mesh information
    //##############################################

    /*
     * Return the STK rank of an entity given the moris
     * entity rank
     */
    stk::mesh::EntityRank
    get_stk_entity_rank(enum EntityRank aMRSEntityRank) const;


    //##############################################
    // internal id functions
    //##############################################
    std::vector<stk::mesh::Entity>
    entities_connected_to_entity_stk(stk::mesh::Entity*    const aInputEntity,
                                            stk::mesh::EntityRank const aInputEntityRank,
                                            stk::mesh::EntityRank const aOutputEntityRank) const;

};
}
}


#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_MESH_STK_HPP_ */

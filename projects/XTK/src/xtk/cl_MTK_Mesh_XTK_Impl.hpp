/*
 * cl_MTK_Mesh_XTK_Impl.hpp
 *
 *  Created on: Feb 20, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_MTK_MESH_XTK_IMPL_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_MTK_MESH_XTK_IMPL_HPP_

#include "cl_MTK_Mesh.hpp"

namespace xtk
{
class Model;
}

namespace moris
{

namespace mtk
{
class XTK_Impl: public Mesh
{
public :
    XTK_Impl(){}

    XTK_Impl(xtk::Model* aModelPtr);

    ~XTK_Impl()
    {
        delete mOutputMeshPtr;
    }


    //##############################################
    // General mesh information access
    //##############################################

    //------------------------------------------------------------------------------

    /**
     * returns the type enum of this mesh
     */
    MeshType
    get_mesh_type() const
    {
        return MeshType::XTK;

    }

    //------------------------------------------------------------------------------

    uint
    get_spatial_dim() const
    {
        return 3;
    }

    //------------------------------------------------------------------------------
    /*
     * Get number of entities for specified rank
     */
    uint
    get_num_entities(
            enum EntityRank aEntityRank) const;
    //------------------------------------------------------------------------------


    const Matrix< DDRMat > &
    get_t_matrix_of_node_loc_ind(
            const moris_index aNodeIndex,
            const EntityRank  aBSplineRank );

//------------------------------------------------------------------------------

    Matrix< IndexMat >
    get_bspline_inds_of_node_loc_ind(
            const moris_index aNodeIndex,
            const EntityRank  aBSplineRank );

    //------------------------------------------------------------------------------
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
    //------------------------------------------------------------------------------
    /*
     * Since the connectivity between entities of the same rank are considered
     * invalid by STK standards, we need a separate function for element to element
     * specifically
     *      *
     * @param[in]  aElementId - element id
     * @param[out] A 2 row matrix where the first row it the neighbor elements index and the
     *             second row is the shared face ordinal corresponding to the neighbor
     */
    
    Matrix< IndexMat >
    get_elements_connected_to_element_and_face_ord_loc_inds(moris_index aElementIndex) const
    {
        MORIS_ERROR(0,"Entered  function in Mesh base class, (function is not implemented)");
        return Matrix<IndexMat>(0,0);
    }

    /*
     * Since the connectivity between entities of the same rank are considered
     * invalid by STK standards, we need a seperate function for element to element
     * specifically
     *
     * @param[in]  aElementId - element id
     * @param[out] Element to element connectivity and face index shared
     *                   (where elements are all by index)
     */

    
    Matrix< IndexMat >
    get_elements_connected_to_element_and_face_ind_loc_inds(moris_index aElementIndex) const
    {
        MORIS_ERROR(0,"Entered  function in Mesh base class, (function is not implemented)");
        return Matrix<IndexMat>(0,0);
    }

    //------------------------------------------------------------------------------
    /*
     * Get elements connected to node
     */
    
    Matrix < IndexMat >
    get_elements_connected_to_node_loc_inds( moris_index aNodeIndex ) const
    {
        return get_entity_connected_to_entity_loc_inds(aNodeIndex,EntityRank::NODE, EntityRank::ELEMENT);
    }
    //------------------------------------------------------------------------------
    /*
     * Get faces connected to node
     */
    
    Matrix < IndexMat >
    get_faces_connected_to_node_loc_inds( moris_index aNodeIndex ) const
    {
        return get_entity_connected_to_entity_loc_inds(aNodeIndex,EntityRank::NODE, EntityRank::FACE);
    }
    //------------------------------------------------------------------------------
    /*
     * Get edges connected to node
     */
    
    Matrix < IndexMat >
    get_edges_connected_to_node_loc_inds( moris_index aNodeIndex ) const
    {
        return get_entity_connected_to_entity_loc_inds(aNodeIndex,EntityRank::NODE, EntityRank::EDGE);
    }
    //------------------------------------------------------------------------------
    /*
     * Get elements connected to edge
     */
    
    Matrix < IndexMat >
    get_elements_connected_to_edge_loc_inds( moris_index aEdgeIndex ) const
    {
        return get_entity_connected_to_entity_loc_inds(aEdgeIndex,EntityRank::EDGE, EntityRank::ELEMENT);
    }
    //------------------------------------------------------------------------------
    /*
     * Get faces connected to edge
     */
    
    Matrix < IndexMat >
    get_faces_connected_to_edge_loc_inds( moris_index aEdgeIndex ) const
    {
        return get_entity_connected_to_entity_loc_inds(aEdgeIndex,EntityRank::EDGE, EntityRank::FACE);
    }
    //------------------------------------------------------------------------------
    
    Matrix< IndexMat >
    get_elements_connected_to_face_loc_inds( moris_index aFaceIndex ) const
    {
        return get_entity_connected_to_entity_loc_inds(aFaceIndex,EntityRank::FACE, EntityRank::ELEMENT);
    }
    //------------------------------------------------------------------------------
    /*
      * Get faces connected to an element
      */
     
     Matrix< IndexMat >
     get_faces_connected_to_element_loc_inds(moris_index aElementIndex) const
     {
         return get_entity_connected_to_entity_loc_inds(aElementIndex,EntityRank::ELEMENT, EntityRank::FACE);
     }
     //------------------------------------------------------------------------------
     /*
      * Get edges connected to an element
      */
     
     Matrix< IndexMat >
     get_edges_connected_to_element_loc_inds(moris_index aElementIndex) const
     {
         return get_entity_connected_to_entity_loc_inds(aElementIndex,EntityRank::ELEMENT, EntityRank::EDGE);
     }

     //------------------------------------------------------------------------------
     /*
      * Get nodes connected to an element
      */
     
     Matrix< IndexMat >
     get_nodes_connected_to_element_loc_inds(moris_index aElementIndex) const
     {
         return get_entity_connected_to_entity_loc_inds(aElementIndex,EntityRank::ELEMENT, EntityRank::NODE);
     }

     //------------------------------------------------------------------------------
     /*
      * Get number of basis functions. For Lagrange meshes, the number of basis functions and the number of nodes
      * are equivalent. Therefore, a default implementation using get_num_nodes() is used here.
      */

     
     uint
     get_num_basis_functions()
     {
         return get_num_nodes();
     }

     /*
      * Get elements interpolated into by a basis function. For a Lagrange mesh,
      * the elements in support of basis is equivalent to the elements connected
      * to a node. Therefore, a call to get_elements
      */
     
     Matrix< IndexMat >
     get_elements_in_support_of_basis(moris_index aBasisIndex)
     {
         return get_entity_connected_to_entity_loc_inds(aBasisIndex, EntityRank::NODE, EntityRank::ELEMENT);
     }

     //##############################################
     // global id functions
     //##############################################
     /*
      * Get global identifier of an entity from a local index and entity rank
      */
     
     moris_id
     get_glb_entity_id_from_entity_loc_index(moris_index     aEntityIndex,
                                             enum EntityRank aEntityRank) const
     {
         MORIS_ERROR(0,"Entered  function in Mesh base class, (function is not implemented)");
         return 0;
     }
     //------------------------------------------------------------------------------
     /*
      * Get global identifier of an entity from a local index and entity rank
      */
     
     moris_index
     get_loc_entity_ind_from_entity_glb_id(moris_id        aEntityId,
                                           enum EntityRank aEntityRank) const
     {
         MORIS_ERROR(0,"Entered  function in Mesh base class, (function is not implemented)");
         return 0;
     }
     //------------------------------------------------------------------------------
     /*
      * Generic get global id of entities connected to
      * entity using an entities global id
      */
     
     Matrix<IdMat>
     get_entity_connected_to_entity_glob_ids( moris_id     aEntityId,
                                             enum EntityRank aInputEntityRank,
                                             enum EntityRank aOutputEntityRank) const
      {
         MORIS_ERROR(0,"Entered  function in Mesh base class, (function is not implemented)");
         return Matrix<IdMat>(0,0);
      }
     //------------------------------------------------------------------------------
     /*
      * Since the connectivity between entities of the same rank are considered
      * invalid by STK standards, we need a seperate function for element to element
      * specifically
      *
      * @param[in]  aElementId - element id
      * @param[out] Element to element connectivity and face ordinal shared
      */
     
     Matrix< IdMat >
     get_element_connected_to_element_glob_ids(moris_id aElementId) const
     {
         MORIS_ERROR(0,"Entered  function in Mesh base class, (function is not implemented)");
         return Matrix<IdMat>(0,0);
     }
     //------------------------------------------------------------------------------
     /*
      * Returns a list of globally unique entity ids for entities
      * of the provided rank
      * @param[in]  aNumNodes - number of node ids requested
      * @param[in]  aEntityRank - Entity rank to assign ids for
      * @param[out] aAvailableNodeIDs - list of globally unique node IDs
      */
     
     Matrix< IdMat >
     generate_unique_entity_ids( uint            aNumEntities,
                                 enum EntityRank aEntityRank) const
      {
         return Matrix<IdMat>(1,1,this->get_num_entities(aEntityRank)+1);
      }

     //------------------------------------------------------------------------------
     /*
      * Returns a matrix of globally unique node ids
      * @aNumNodes -- number of node ids needed
      */
     
     Matrix < IdMat >
     generate_unique_node_ids(uint aNumNodes)
     {
         return generate_unique_entity_ids(aNumNodes,EntityRank::NODE);
     }

     //------------------------------------------------------------------------------
     /*
      * Get elements connected to node
      */
     
     Matrix < IdMat >
     get_elements_connected_to_node_glob_ids( moris_id aNodeId )
     {
         return get_entity_connected_to_entity_glob_ids(aNodeId,EntityRank::NODE, EntityRank::ELEMENT);
     }
     //------------------------------------------------------------------------------
     /*
      * Get faces connected to node
      */
     
     Matrix < IdMat >
     get_faces_connected_to_node_glob_ids( moris_id aNodeId )
     {
         return get_entity_connected_to_entity_glob_ids(aNodeId,EntityRank::NODE, EntityRank::FACE);
     }
     //------------------------------------------------------------------------------
     /*
      * Get edges connected to node
      */
     
     Matrix < IdMat >
     get_edges_connected_to_node_glob_ids( moris_id aNodeId )
     {
         return get_entity_connected_to_entity_glob_ids(aNodeId,EntityRank::NODE, EntityRank::EDGE);
     }

     //------------------------------------------------------------------------------
     /*
      * Get elements connected to edge
      */
     
     Matrix < IdMat >
     get_elements_connected_to_edge_glob_ids( moris_id aEdgeId )
     {
         return get_entity_connected_to_entity_glob_ids(aEdgeId,EntityRank::EDGE, EntityRank::ELEMENT);
     }
     //------------------------------------------------------------------------------
     /*
      * Get faces connected to edge
      */
     
     Matrix < IdMat >
     get_faces_connected_to_edge_glob_ids( moris_id aEdgeId )
     {
         return get_entity_connected_to_entity_glob_ids(aEdgeId,EntityRank::EDGE, EntityRank::FACE);
     }
     //------------------------------------------------------------------------------
     /*
      * Get elements connected to face
      */
     
     Matrix< IdMat >
     get_elements_connected_to_face_glob_ids( moris_id aFaceId )
     {
         return get_entity_connected_to_entity_glob_ids(aFaceId,EntityRank::FACE, EntityRank::ELEMENT);
     }
     //------------------------------------------------------------------------------
     /*
      * Get faces connected to an element
      */
     
     Matrix< IdMat >
     get_faces_connected_to_element_glob_ids(moris_id aElementId)
     {
         return get_entity_connected_to_entity_glob_ids(aElementId,EntityRank::ELEMENT, EntityRank::FACE);
     }
     //------------------------------------------------------------------------------
     /*
      * Get edges connected to an element
      */
     
     Matrix< IdMat >
     get_edges_connected_to_element_glob_ids(moris_id aElementId)
     {
         return get_entity_connected_to_entity_glob_ids(aElementId,EntityRank::ELEMENT, EntityRank::EDGE);
     }
     //------------------------------------------------------------------------------
     /*
      * Get nodes connected to an element
      */
     
     Matrix< IdMat >
     get_nodes_connected_to_element_glob_ids(moris_id aElementId)
     {
         return get_entity_connected_to_entity_glob_ids(aElementId,EntityRank::ELEMENT, EntityRank::NODE);
     }
     //------------------------------------------------------------------------------
     //##############################################
     // Coordinate Field Functions
     //##############################################
     /*
      * Get coordinate of a node
      */
     
     Matrix< DDRMat >
     get_node_coordinate( moris_index aNodeIndex ) const
     {
         MORIS_ERROR(0,"Entered  function in Mesh base class, (function is not implemented)");
         return Matrix<DDRMat>(0,0);
     }

     //------------------------------------------------------------------------------
     //##############################################
     // Entity Ownership Functions
     //##############################################

     /*
      * Get the entity owner
      */
     
     moris_id
     get_entity_owner(  moris_index     aEntityIndex,
                        enum EntityRank aEntityRank ) const
     {
         MORIS_ERROR(0," get entity owner has no base implementation");
         return 0;
     }

     /*
      * Processors whom share a given entity
      * @param[in]  - Entity Index
      * @param[in]  - Entity Rank
      * @param[out] - Processors whom share an entity vector
      */
     
     void
     get_processors_whom_share_entity(moris_index       aEntityIndex,
                                      enum EntityRank   aEntityRank,
                                      Matrix< IdMat > & aProcsWhomShareEntity) const
     {
         MORIS_ERROR(0," get_processors_whom_share_entity has no base implementation");
     }

     
     uint
     get_num_of_entities_shared_with_processor(moris_id        aProcessorRank,
                                               enum EntityRank aEntityRank,
                                               bool aSendFlag) const
     {
         MORIS_ERROR(0," get_num_of_entities_shared_with_processor has no base implementation");
         return 0;
     }

     //##############################################
     // Mesh Sets Access
     //##############################################

     
     Matrix< IndexMat >
     get_set_entity_loc_inds( enum EntityRank aSetEntityRank,
                              std::string     aSetName) const
     {
         MORIS_ERROR(0," get_set_entity_ids has no base implementation");
         return Matrix< IndexMat >(0,0);
     }

     //##############################################
     // Field Access
     //##############################################

     /*
      * Access an entity
      *
      */
     //TODO: introduce a concept of field indices to prevent accessing via a name which
     //TODO: involves a string comparison
     
     Matrix< DDRMat >
     get_entity_field_value_real_scalar(const Matrix<IndexMat> & aEntityIndices,
                                                      const std::string & aFieldName,
                                        enum EntityRank            aFieldEntityRank) const
     {
         MORIS_ERROR(0,"Entered  function in Mesh base class, (get_entity_field_value_real_scalar is not implemented)");
         return Matrix< DDRMat >(0,0);
     }


     /*
      * Given a field name and rank associated with field, add the field data
      * For now, this is just for real type single component fields
      *
      */
     
     void
     add_mesh_field_real_scalar_data_loc_inds(const std::string & aFieldName,
                                          const enum EntityRank & aFieldEntityRank,
                                          const Matrix<DDRMat> & aFieldData)
     {
         MORIS_ERROR(0,"Entered  function in Mesh base class, (add_mesh_field_real_scalar_data_loc_inds is not implemented)");

     }


     //------------------------------------------------------------------------------
     //##############################################
     // Face Cluster Access
     //##############################################

     /*
      * Is this face a member of a face cluster?
      */
     
     bool
     has_face_cluster_membership(moris_index aFaceIndex) const
     {
         MORIS_ERROR(0,"Entered  function in Mesh base class, (has_face_cluster_membership is not implemented)");
         return 0;
     }

     //------------------------------------------------------------------------------

     /*
      * Is this face a parent of a face cluster?
      */
     
     bool
     is_face_cluster_a_parent(moris_index aFaceIndex)
     {
         MORIS_ERROR(0,"Entered  function in Mesh base class, (has_face_cluster_membership is not implemented)");
         return 0;
     }

     //------------------------------------------------------------------------------

     /*
      * Get the a face cluster (provided a parent face index)
      */
     
     const Facet_Cluster &
     get_face_cluster(moris_index aParentFaceIndex)
     {
         MORIS_ERROR(0,"Entered  function in Mesh base class, (register_face_cluster is not implemented)");
         return mDummyFaceCluster;
     }

     //------------------------------------------------------------------------------

     /*
      * Manually add a face cluster to mesh.
      */
     
     void
     register_face_cluster(Facet_Cluster & aFaceCluster)
     {
         MORIS_ERROR(0,"Entered  function in Mesh base class, (register_face_cluster is not implemented)");
     }


     //------------------------------------------------------------------------------
     //##############################################
     // Cell and Vertex Pointer Functions
     //##############################################
     /*
      * Returns a reference to a cell in the mesh
      */
     
     mtk::Cell  &
     get_mtk_cell( moris_index aElementIndex);

     mtk::Cell const &
     get_mtk_cell( moris_index aElementIndex) const;

     /*
      * Returns a reference to a vertex in the mesh
      */
     mtk::Vertex &
     get_mtk_vertex( moris_index aVertexIndex );

     //##############################################
     // For FEM
     //##############################################

      mtk::Cell  &
     get_writable_mtk_cell( moris_index aElementIndex )
     {
         MORIS_ERROR(0,"Entered  function in Mesh base class, (function is not implemented)");
         return *mDummyCells;
     }

     void
      get_adof_map( const uint aOrder,
                           map< moris_id, moris_index > & aAdofMap ) const
     {
         MORIS_ERROR(0,"Entered  function in Mesh base class, (function is not implemented)");
     }
     //--------------------------------------------------------------
     
     moris_id
     get_max_entity_id( enum EntityRank aEntityRank ) const;
     //------------------------------------------------------------------------------


    //FIXME: THIS FUNCTION DESCRIPTION NEEDS TO BE IMPROVED
    //FIXME: Also, a unit test (not clear what needs to be provided)
    /**
     * provides a moris::Mat<uint> containing the IDs this mesh has
     * to communicate with
     *
     */
     Matrix< IdMat >
    get_communication_table() const
    {
         MORIS_ERROR(0,"Communication table not implemented in XTK");
         return Matrix< IdMat >(0,0);
    }

    //------------------------------------------------------------------------------

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
            std::string  &aFileName,
            bool          aAddElemCmap = false)
    {
        MORIS_ERROR(0,"Create output mesh not implemented");
    }
    //------------------------------------------------------------------------------
//------------------------------------------------------------------------------

    //##############################################
    //  Field Functions
    //##############################################

    /**
     * return the number of fields that are connected to this field
     */
     uint
    get_num_fields(  const enum EntityRank aEntityRank ) const
    {
        MORIS_ERROR( false ,"get_num_fields() not implemented" );
        return 0;
    }

//------------------------------------------------------------------------------

    /**
     * return the index of the field of this label
     * return gNoIndex if not found
     */
     moris_index
    get_field_ind(
            const std::string & aFieldLabel,
            const enum EntityRank aEntityRank ) const
    {
        MORIS_ERROR( false ,"get_field_ind() not implemented" );
        return gNoIndex;
    }

//------------------------------------------------------------------------------

    /**
     * return the interpolation order of this field
     */
     uint
    get_order_of_field(
            const moris_index     aFieldIndex,
            const enum EntityRank aEntityRank )
    {
        MORIS_ERROR( false ,"get_order_of_field() not implemented" );
        return 0;
    }

    //------------------------------------------------------------------------------

    /**
     * add a scalar field to the database
     *
     * fixme: how to make sure that field does not exist ?
     */
     moris_index
    create_scalar_field(
            const std::string   & aFieldLabel,
            const enum EntityRank aEntityRank )
    {
        MORIS_ERROR( false ,"create_scalar_field() not implemented" );
        return gNoIndex;
    }

    //------------------------------------------------------------------------------

    /**
     * add a vector field to the database
     */
     moris_index
    create_vector_field(
            const std::string   & aFieldLabel,
            const enum EntityRank aEntityRank,
            const uint            aDimension )
    {
        MORIS_ERROR( false ,"create_vector_field() not implemented" );
        return gNoIndex;
    }

    //------------------------------------------------------------------------------

    /**
     * get value of entity
     */
     real &
    get_value_of_scalar_field(
            const moris_index     aFieldIndex,
            const enum EntityRank aEntityRank,
            const uint            aEntityIndex )
    {
        MORIS_ERROR( false ,"get_value_of_scalar_field() not implemented" );
        return mDummyReal;
    }

    //------------------------------------------------------------------------------

    /**
     * get value of entity ( const version )
     */
     const real &
    get_value_of_scalar_field(
            const moris_index     aFieldIndex,
            const enum EntityRank aEntityRank,
            const uint            aEntityIndex ) const
    {
        MORIS_ERROR( false ,"get_value_of_scalar_field() const not implemented" );
        return mDummyReal;
    }

    //------------------------------------------------------------------------------

    /**
     * fixme: need opinion: sould we always return a DDRMat?
     *        should this be a row or column vector?
     */
     Matrix<DDRMat> &
    get_value_of_vector_field(
            const moris_index     aFieldIndex,
            const enum EntityRank aEntityRank,
            const uint            aEntityIndex )
            {
        MORIS_ERROR( false ,"get_value_of_vector_field() not implemented" );
        return  mDummyMatrix;
            }

    //------------------------------------------------------------------------------

    /**
     * return the entry of a vector field ( const version )
     */
     const Matrix<DDRMat> &
    get_value_of_vector_field(
            const moris_index     aFieldIndex,
            const enum EntityRank aEntityRank,
            const uint            aEntityIndex ) const
            {
        MORIS_ERROR( false ,"get_value_of_vector_field() not implemented" );
        return mDummyMatrix;
            }

    //------------------------------------------------------------------------------

    /**
     * returns a moris::Matrix with the field
     * This function is specific to HMR, and called by the mapper
     * if HMR is used.
     */
     Matrix<DDRMat> &
    get_field( const moris_index     aFieldIndex,
               const enum EntityRank aEntityRank )
               {
        MORIS_ERROR( false ,"get_field() not implemented" );
        return mDummyMatrix;
               }

//------------------------------------------------------------------------------

    //##############################################
    //  Multigrid
    //##############################################

    /**
     * returns the number of levels
     */
     uint
    get_num_level( const enum EntityRank aEntityRank )
    {
        // no error is thrown here
        return 0;
    }

     //------------------------------------------------------------------------------

    /**
     * needed for multigrid and HMR
     */
     uint
    get_max_level_of_entity( const enum EntityRank aEntityRank )
    {
        // no error is thrown here
        return 0;
    }

//------------------------------------------------------------------------------

    /**
     * returns the level of an entity. Makes only sense for HMR
     */
    uint
    get_level_of_entity_loc_ind( const enum EntityRank aEntityRank,
                                 const uint            aEntityIndex )
    {
        // no error is thrown here
        return 0;
    }


protected:
    xtk::Model*       mXTKModelPtr;
    moris::mtk::Mesh* mOutputMeshPtr;

};
}
}


#endif /* PROJECTS_XTK_SRC_XTK_CL_MTK_MESH_XTK_IMPL_HPP_ */

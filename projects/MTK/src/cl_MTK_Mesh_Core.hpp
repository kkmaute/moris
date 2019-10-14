/*
 * cl_MTK_Mesh_Core.hpp
 *
 *  Created on: Apr 15, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_MESH_CORE_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_MESH_CORE_HPP_

#include "assert.hpp"
#include "cl_Matrix.hpp"

#include "cl_Mesh_Enums.hpp"
#include "MTK_Tools.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_Map.hpp"
#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Cell.hpp" //MTK/src
#include "cl_MTK_Facet.hpp"

namespace moris
{
namespace hmr
{
class Database;
class Lagrange_Mesh_Base;
}

namespace mtk
{
class Mesh :  public std::enable_shared_from_this< Mesh >
{
protected:
    // Note these members are here only to allow for throwing in
    // get_mtk_cell and get_mtk_vertex function
    mtk::Vertex*     mDummyVertex;
    mtk::Cell*       mDummyCells;
    real             mDummyReal = 0.0;
    Matrix<DDRMat>   mDummyMatrix;

    //------------------------------------------------------------------------------
    //! ref to hmr object for multigrid
    std::shared_ptr< hmr::Database > mDatabase;

    hmr::Lagrange_Mesh_Base * mMesh = nullptr;
    //------------------------------------------------------------------------------

public:
    // Verbose flag
    bool mVerbose = false;

    /**
     * trivial constructor
     */
    Mesh(){};

    //------------------------------------------------------------------------------

    /**
     * virtual destructor
     */
    virtual
    ~Mesh(){};

    //##############################################
    // 1.) General mesh information access
    //##############################################

    //------------------------------------------------------------------------------

    /**
     * returns the type enum of this mesh
     */
    virtual MeshType
    get_mesh_type() const = 0;

    //------------------------------------------------------------------------------

    virtual
    uint
    get_spatial_dim() const = 0;

    //------------------------------------------------------------------------------
    /*
     * Get number of entities for specified rank
     */
    virtual
    uint
    get_num_entities(
            enum EntityRank aEntityRank) const = 0;

    //------------------------------------------------------------------------------
    // end of pure virtual functions in section 1
    // all functions below this line need to be able to have a default implementation
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    /*
     * Get number of nodes
     */
    virtual
    uint
    get_num_nodes() const
    {
        return get_num_entities(EntityRank::NODE);
    }
    //------------------------------------------------------------------------------
    /*
     * Get number of edges
     */
    virtual
    uint
    get_num_edges() const
    {
        return get_num_entities(EntityRank::EDGE);
    }
    //------------------------------------------------------------------------------
    /*
     * Get number of faces
     */
    virtual
    uint
    get_num_faces() const
    {
        return get_num_entities(EntityRank::FACE);
    }
    //------------------------------------------------------------------------------
    /*
     * Get number of elements
     */
    virtual
    uint
    get_num_elems() const
    {
        return get_num_entities(EntityRank::ELEMENT);
    }

    //------------------------------------------------------------------------------
    //##############################################
    // 2.) Access Mesh Data by index Functions
    //##############################################
    //------------------------------------------------------------------------------
    //##############################################
    // 2.a.) Access standard mesh data
    //##############################################
    /*
     * Generic get local index of entities connected to
     * entity using an entities local index
     */
    virtual
    Matrix<IndexMat>
    get_entity_connected_to_entity_loc_inds(moris_index     aEntityIndex,
                                            enum EntityRank aInputEntityRank,
                                            enum EntityRank aOutputEntityRank) const = 0;
    //------------------------------------------------------------------------------
    /*
     * Since the connectivity between entities of the same rank are considered
     * invalid by STK standards, we need a separate function for element to element
     * specifically
     *
     * @param[in]  aElementId - element id
     * @param[out] A 2 row matrix where the first row it the neighbor elements index and the
     *             second row is the shared face ordinal corresponding to the neighbor
     */

    virtual
    Matrix< IndexMat >
    get_elements_connected_to_element_and_face_ord_loc_inds(moris_index aElementIndex) const
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
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
    virtual
    Matrix< IndexMat >
    get_elements_connected_to_element_and_face_ind_loc_inds(moris_index aElementIndex) const = 0;


    /*
     *  Returns all the vertices
     */
    virtual
    moris::Cell<moris::mtk::Vertex const *>
    get_all_vertices() const
    {
        MORIS_ERROR(0,"No default implementation of get_all_vertices_no_aura");

        return moris::Cell<moris::mtk::Vertex const *>(0,nullptr);
    }


    //------------------------------------------------------------------------------
    // end of pure virtual functions in section 2.1
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    /*
     * Get elements connected to node
     */
    virtual
    Matrix < IndexMat >
    get_elements_connected_to_node_loc_inds( moris_index aNodeIndex ) const
    {
        return get_entity_connected_to_entity_loc_inds(aNodeIndex,EntityRank::NODE, EntityRank::ELEMENT);
    }
    //------------------------------------------------------------------------------
    /*
     * Get faces connected to node
     */
    virtual
    Matrix < IndexMat >
    get_faces_connected_to_node_loc_inds( moris_index aNodeIndex ) const
    {
        return get_entity_connected_to_entity_loc_inds(aNodeIndex,EntityRank::NODE, EntityRank::FACE);
    }
    //------------------------------------------------------------------------------
    /*
     * Get edges connected to node
     */
    virtual
    Matrix < IndexMat >
    get_edges_connected_to_node_loc_inds( moris_index aNodeIndex ) const
    {
        return get_entity_connected_to_entity_loc_inds(aNodeIndex,EntityRank::NODE, EntityRank::EDGE);
    }
    //------------------------------------------------------------------------------
    /*
     * Get elements connected to edge
     */
    virtual
    Matrix < IndexMat >
    get_elements_connected_to_edge_loc_inds( moris_index aEdgeIndex ) const
    {
        return get_entity_connected_to_entity_loc_inds(aEdgeIndex,EntityRank::EDGE, EntityRank::ELEMENT);
    }
    //------------------------------------------------------------------------------
    /*
     * Get faces connected to edge
     */
    virtual
    Matrix < IndexMat >
    get_faces_connected_to_edge_loc_inds( moris_index aEdgeIndex ) const
    {
        return get_entity_connected_to_entity_loc_inds(aEdgeIndex,EntityRank::EDGE, EntityRank::FACE);
    }
    //------------------------------------------------------------------------------
    virtual
    Matrix< IndexMat >
    get_elements_connected_to_face_loc_inds( moris_index aFaceIndex ) const
    {
        return get_entity_connected_to_entity_loc_inds(aFaceIndex,EntityRank::FACE, EntityRank::ELEMENT);
    }
    //------------------------------------------------------------------------------
    /*
     * Get faces connected to an element
     */
    virtual
    Matrix< IndexMat >
    get_faces_connected_to_element_loc_inds(moris_index aElementIndex) const
    {
        return get_entity_connected_to_entity_loc_inds(aElementIndex,EntityRank::ELEMENT, EntityRank::FACE);
    }
    //------------------------------------------------------------------------------
    /*
     * Get edges connected to an element
     */
    virtual
    Matrix< IndexMat >
    get_edges_connected_to_element_loc_inds(moris_index aElementIndex) const
    {
        return get_entity_connected_to_entity_loc_inds(aElementIndex,EntityRank::ELEMENT, EntityRank::EDGE);
    }

    //------------------------------------------------------------------------------
    /*
     * Get nodes connected to an element
     */
    virtual
    Matrix< IndexMat >
    get_nodes_connected_to_element_loc_inds(moris_index aElementIndex) const
    {
        return get_entity_connected_to_entity_loc_inds(aElementIndex,EntityRank::ELEMENT, EntityRank::NODE);
    }

    //------------------------------------------------------------------------------

    //##############################################
    // 2.a.) Access mesh data from ids
    //##############################################

    //##############################################
    // global id functions
    //##############################################
    /*
     * Get global identifier of an entity from a local index and entity rank
     */
    virtual
    moris_id
    get_glb_entity_id_from_entity_loc_index(moris_index     aEntityIndex,
                                            enum EntityRank aEntityRank) const
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return 0;
    }
    //------------------------------------------------------------------------------
    /*
     * Get global identifier of an entity from a local index and entity rank
     */
    virtual
    moris_index
    get_loc_entity_ind_from_entity_glb_id(moris_id        aEntityId,
                                          enum EntityRank aEntityRank) const
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return 0;
    }
    //------------------------------------------------------------------------------
    /*
     * Generic get global id of entities connected to
     * entity using an entities global id
     */
    virtual
    Matrix<IdMat>
    get_entity_connected_to_entity_glob_ids( moris_id     aEntityId,
                                             enum EntityRank aInputEntityRank,
                                             enum EntityRank aOutputEntityRank) const
                                             {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
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
    virtual
    Matrix< IdMat >
    get_element_connected_to_element_glob_ids(moris_id aElementId) const
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return Matrix<IdMat>(0,0);
    }

    virtual
    moris::moris_index
    get_facet_ordinal_from_cell_and_facet_loc_inds(moris::moris_index aFaceIndex,
                                                   moris::moris_index aCellIndex) const
    {
        Matrix<IdMat> tElementFaces = get_entity_connected_to_entity_loc_inds(aCellIndex,EntityRank::ELEMENT, EntityRank::FACE);

        moris_index tOrdinal = MORIS_INDEX_MAX;
        for(moris_index iOrd = 0; iOrd<(moris_index)tElementFaces.numel(); iOrd++)
        {
            if(tElementFaces(iOrd) == aFaceIndex)
            {
                tOrdinal = iOrd;
                return tOrdinal;
            }
        }
        MORIS_ERROR(tOrdinal!=MORIS_INDEX_MAX," Facet ordinal not found");
        return tOrdinal;
    }

    //------------------------------------------------------------------------------
    /*
     * Returns a list of globally unique entity ids for entities
     * of the provided rank
     * @param[in]  aNumNodes - number of node ids requested
     * @param[in]  aEntityRank - Entity rank to assign ids for
     * @param[out] aAvailableNodeIDs - list of globally unique node IDs
     */
    virtual
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
    virtual
    Matrix < IdMat >
    generate_unique_node_ids(uint aNumNodes)
    {
        return generate_unique_entity_ids(aNumNodes,EntityRank::NODE);
    }

    //------------------------------------------------------------------------------
    /*
     * Get elements connected to node
     */
    virtual
    Matrix < IdMat >
    get_elements_connected_to_node_glob_ids( moris_id aNodeId )
    {
        return get_entity_connected_to_entity_glob_ids(aNodeId,EntityRank::NODE, EntityRank::ELEMENT);
    }
    //------------------------------------------------------------------------------
    /*
     * Get faces connected to node
     */
    virtual
    Matrix < IdMat >
    get_faces_connected_to_node_glob_ids( moris_id aNodeId )
    {
        return get_entity_connected_to_entity_glob_ids(aNodeId,EntityRank::NODE, EntityRank::FACE);
    }
    //------------------------------------------------------------------------------
    /*
     * Get edges connected to node
     */
    virtual
    Matrix < IdMat >
    get_edges_connected_to_node_glob_ids( moris_id aNodeId )
    {
        return get_entity_connected_to_entity_glob_ids(aNodeId,EntityRank::NODE, EntityRank::EDGE);
    }

    //------------------------------------------------------------------------------
    /*
     * Get elements connected to edge
     */
    virtual
    Matrix < IdMat >
    get_elements_connected_to_edge_glob_ids( moris_id aEdgeId )
    {
        return get_entity_connected_to_entity_glob_ids(aEdgeId,EntityRank::EDGE, EntityRank::ELEMENT);
    }
    //------------------------------------------------------------------------------
    /*
     * Get faces connected to edge
     */
    virtual
    Matrix < IdMat >
    get_faces_connected_to_edge_glob_ids( moris_id aEdgeId )
    {
        return get_entity_connected_to_entity_glob_ids(aEdgeId,EntityRank::EDGE, EntityRank::FACE);
    }
    //------------------------------------------------------------------------------
    /*
     * Get elements connected to face
     */
    virtual
    Matrix< IdMat >
    get_elements_connected_to_face_glob_ids( moris_id aFaceId )
    {
        return get_entity_connected_to_entity_glob_ids(aFaceId,EntityRank::FACE, EntityRank::ELEMENT);
    }
    //------------------------------------------------------------------------------
    /*
     * Get faces connected to an element
     */
    virtual
    Matrix< IdMat >
    get_faces_connected_to_element_glob_ids(moris_id aElementId)
    {
        return get_entity_connected_to_entity_glob_ids(aElementId,EntityRank::ELEMENT, EntityRank::FACE);
    }
    //------------------------------------------------------------------------------
    /*
     * Get edges connected to an element
     */
    virtual
    Matrix< IdMat >
    get_edges_connected_to_element_glob_ids(moris_id aElementId)
    {
        return get_entity_connected_to_entity_glob_ids(aElementId,EntityRank::ELEMENT, EntityRank::EDGE);
    }
    //------------------------------------------------------------------------------
    /*
     * Get nodes connected to an element
     */
    virtual
    Matrix< IdMat >
    get_nodes_connected_to_element_glob_ids(moris_id aElementId)
    {
        return get_entity_connected_to_entity_glob_ids(aElementId,EntityRank::ELEMENT, EntityRank::NODE);
    }

    /*
     * Get elements interpolated into by a basis function. For a Lagrange mesh,
     * the elements in support of basis is equivalent to the elements connected
     * to a node. Therefore, a call to get_elements
     */
    virtual
    void
    get_elements_in_support_of_basis(const uint           aMeshIndex,
                                     const uint           aBasisIndex,
                                     Matrix< IndexMat > & aElementIndices )
    {
        MORIS_ERROR(0,"get_elements_in_support_of_basis not implemented");
    }

    //------------------------------------------------------------------------------
    //##############################################
    // Coordinate Field Functions
    //##############################################
    /*
     * Get coordinate of a node
     */
    virtual
    Matrix< DDRMat >
    get_node_coordinate( moris_index aNodeIndex ) const
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return Matrix<DDRMat>(0,0);
    }

    //------------------------------------------------------------------------------

    //##############################################
    // Field Access
    //##############################################

    /*
     * Access an entity
     *
     */
    //TODO: introduce a concept of field indices to prevent accessing via a name which
    //TODO: involves a string comparison
    virtual
    Matrix< DDRMat >
    get_entity_field_value_real_scalar(const Matrix<IndexMat> & aEntityIndices,
                                       const std::string & aFieldName,
                                       enum EntityRank            aFieldEntityRank) const
                                       {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (get_entity_field_value_real_scalar is not implemented)");
        return Matrix< DDRMat >(0,0);
                                       }


    /*
     * Given a field name and rank associated with field, add the field data
     * For now, this is just for real type single component fields
     *
     */
    virtual
    void
    add_mesh_field_real_scalar_data_loc_inds(const std::string & aFieldName,
                                             const enum EntityRank & aFieldEntityRank,
                                             const Matrix<DDRMat> & aFieldData)
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (add_mesh_field_real_scalar_data_loc_inds is not implemented)");

    }


    //------------------------------------------------------------------------------
    //##############################################
    // Facet Access
    //##############################################

    virtual
    moris::mtk::Facet*
    get_facet(moris_index)
    {
        MORIS_ERROR(0,"get facet not implemented");
        return nullptr;
    }

    //------------------------------------------------------------------------------
    //##############################################
    // Cell and Vertex Pointer Functions
    //##############################################
    /*
     * Returns a reference to a cell in the mesh
     */
    virtual
    mtk::Cell  &
    get_mtk_cell( moris_index aElementIndex)
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return *mDummyCells;
    }

    virtual
    mtk::Cell const &
    get_mtk_cell( moris_index aElementIndex) const
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return *mDummyCells;
    }

    /*
     * Returns a reference to a vertex in the mesh
     */
    virtual
    mtk::Vertex &
    get_mtk_vertex( moris_index aVertexIndex )
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return *mDummyVertex;
    }

    virtual
    mtk::Vertex const &
    get_mtk_vertex( moris_index aVertexIndex ) const
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return *mDummyVertex;
    }

    //##############################################
    // For FEM
    //##############################################

    virtual mtk::Cell  &
    get_writable_mtk_cell( moris_index aElementIndex )
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return *mDummyCells;
    }

    //--------------------------------------------------------------
    // FIXME: REMOVE SINCE NOT USED
    virtual
    moris_id
    get_max_entity_id( enum EntityRank aEntityRank ) const
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return 0;
    }

    //------------------------------------------------------------------------------

    //##############################################
    // Entity Ownership Functions
    //##############################################

    /*
     * Get the entity owner
     */
    virtual
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
    virtual
    void
    get_processors_whom_share_entity(moris_index       aEntityIndex,
                                     enum EntityRank   aEntityRank,
                                     Matrix< IdMat > & aProcsWhomShareEntity) const
    {
        MORIS_ERROR(0," get_processors_whom_share_entity has no base implementation");
    }

    virtual
    uint
    get_num_of_entities_shared_with_processor(moris_id        aProcessorRank,
                                              enum EntityRank aEntityRank,
                                              bool aSendFlag) const
    {
        MORIS_ERROR(0," get_num_of_entities_shared_with_processor has no base implementation");
        return 0;
    }

    //------------------------------------------------------------------------------

    //##############################################
    // Communication Tables
    //##############################################

    //FIXME: THIS FUNCTION DESCRIPTION NEEDS TO BE IMPROVED
    //FIXME: Also, a unit test (not clear what needs to be provided)
    /**
     * provides a moris::Mat<uint> containing the IDs this mesh has
     * to communicate with
     *
     */
    virtual Matrix< IdMat >
    get_communication_table() const = 0;

    virtual
    Matrix< IdMat >
    get_communication_proc_ranks() const
    {
        MORIS_ERROR(0,"get_communication_proc_ranks not implemented");
        return Matrix< IdMat >(0,0);
    }

    virtual
    moris::Cell<Matrix< IdMat >>
    get_communication_vertex_pairing() const
    {
        MORIS_ERROR(0,"get_communication_vertex_pairing not implemented");
        return moris::Cell<Matrix< IdMat >>(0);
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
    virtual
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
    virtual uint
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
    virtual moris_index get_field_ind( const std::string     & aFieldLabel,
                                       const enum EntityRank   aEntityRank ) const
    {
        MORIS_ERROR( false ,"get_field_ind() not implemented" );
        return gNoIndex;
    }

    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------

    /**
     * add a scalar field to the database
     *
     * fixme: how to make sure that field does not exist ?
     */
    virtual moris_index
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
    virtual moris_index
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
    virtual real &
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
    virtual const real &
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
    virtual Matrix<DDRMat> &
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
    virtual const Matrix<DDRMat> &
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
    virtual Matrix<DDRMat> &
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
    virtual uint
    get_num_level( const enum EntityRank aEntityRank )
    {
        // no error is thrown here
        return 0;
    }

    //------------------------------------------------------------------------------

    /**
     * needed for multigrid and HMR
     */
    virtual uint
    get_max_level_of_entity( const enum EntityRank aEntityRank )
    {
        // no error is thrown here
        return 0;
    }

    //------------------------------------------------------------------------------

    /**
     * returns the level of an entity. Makes only sense for HMR
     */
    virtual uint
    get_level_of_entity_loc_ind( const enum EntityRank aEntityRank,
                                 const uint            aEntityIndex )
    {
        // no error is thrown here
        return 0;
    }

    //------------------------------------------------------------------------------
    /**
     * returns HMR database pointer if MTK is build with HMR
     */
    std::shared_ptr< hmr::Database > get_HMR_database( )
    {
        MORIS_ERROR( this->get_mesh_type() == MeshType::HMR ,"Not HMR" );
        return mDatabase;
    }

    hmr::Lagrange_Mesh_Base * get_HMR_lagrange_mesh( )
    {
        MORIS_ERROR( this->get_mesh_type() == MeshType::HMR ,"Not HMR" );
        MORIS_ERROR( mMesh != nullptr ,"get_HMR_lagrange_mesh(), Lagrange mesh is nullptr" );

        return mMesh;
    }

    /*
     * Get number of B-Spline coefficients
     */
    virtual uint
    get_num_coeffs(const uint aBSplineMeshIndex) const
    {
        MORIS_ERROR( false, "get_num_coeffs() not implemented for this mesh" );
        return 0;
    }

    //------------------------------------------------------------------------------
    //FIXME: IMPLEMENT THIS FUNCTION IN STK
    virtual const Matrix< DDRMat > &
    get_t_matrix_of_node_loc_ind(
            const moris_index aNodeIndex,
            const EntityRank  aBSplineRank )
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return mDummyMatrix;
    }

    //FIXME: IMPLEMENT THIS FUNCTION IN XTK
    virtual Matrix< IndexMat >
    get_bspline_inds_of_node_loc_ind(
            const moris_index aNodeIndex,
            const EntityRank  aBSplineRank )
            {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
        return Matrix<IndexMat>(0,0);
            }

    /*
     * Get number of basis functions. For Lagrange meshes, the number of basis functions and the number of nodes
     * are equivalent. Therefore, a default implementation using get_num_nodes() is used here.
     */
    virtual
    uint
    get_num_basis_functions(const uint aMeshIndex = 0)
    {
        return this->get_num_nodes();
    }

    //FIXME: Rename or use get loc entity id from global entity id
    void
    virtual
    get_adof_map( const uint                     aBSplineIndex,
                  map< moris_id, moris_index > & aAdofMap ) const
    {
        MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
    }

    /**
     * return the interpolation order of this field
     */
    virtual uint
    get_order_of_field(
            const moris_index     aFieldIndex,
            const enum EntityRank aEntityRank )
    {
        MORIS_ERROR( false ,"get_order_of_field() not implemented" );
        return 0;
    }

    // fixme: move these functions to integration base class

    /*
     * Returns all the set names ordered by set index for a provided entity rank
     */
    virtual
    moris::Cell<std::string>
    get_set_names(enum EntityRank aSetEntityRank) const
    {
        MORIS_ERROR(0," get_set_names has no base implementation");
        return moris::Cell<std::string>(0);
    }

    /*
     * Returns the local indices of the entities in a set. This includes aura entities
     */
    virtual
    Matrix< IndexMat >
    get_set_entity_loc_inds( enum EntityRank aSetEntityRank,
                             std::string     aSetName) const
                             {
        MORIS_ERROR(0," get_set_entity_loc_inds has no base implementation");
        return Matrix< IndexMat >(0,0);
                             }


    /*
     * Topology of cells in block set
     */
    virtual
    enum CellTopology
    get_blockset_topology(const  std::string & aSetName)
    {
        MORIS_ERROR(0," get_blockset_topology has no base implementation");
        return CellTopology::INVALID;
    }

    /*
     * Topology of sides in side set
     */
    virtual
    enum CellTopology
    get_sideset_topology(const  std::string & aSetName)
    {
        MORIS_ERROR(0," get_sideset_topology has no base implementation");
        return CellTopology::INVALID;
    }


    /*
     * Returns the mtk cells in a block set. Contains the aura entities
     */
    virtual
    moris::Cell<mtk::Cell const *>
    get_block_set_cells( std::string     aSetName) const
    {
        Matrix< IndexMat > tBlockSetElementInd = this->get_set_entity_loc_inds( EntityRank::ELEMENT, aSetName );

        moris::Cell<mtk::Cell const *> tBlockSetCells(tBlockSetElementInd.numel());

        for( luint k=0; k < tBlockSetElementInd.numel(); ++k )
        {
            tBlockSetCells( k ) = & this->get_mtk_cell( tBlockSetElementInd(k) );
        }


        return tBlockSetCells;
    }


    /*
     * Returns the cell index and side ordinals in a provided side set name
     */
    virtual
    void
    get_sideset_elems_loc_inds_and_ords(
            const  std::string     & aSetName,
            Matrix< IndexMat >     & aElemIndices,
            Matrix< IndexMat >     & aSidesetOrdinals ) const
    {
        MORIS_ERROR(0," get_sideset_elems_loc_inds_and_ords has no base implementation");
    }

    /*
     * Returns the mtk cell and side ordinals in a provided side set name
     */
    virtual
    void
    get_sideset_cells_and_ords(
            const  std::string & aSetName,
            moris::Cell< mtk::Cell const * > & aCells,
            Matrix< IndexMat > &       aSidesetOrdinals ) const
    {
        moris::Matrix<moris::IndexMat> tElemIndices;
        this->get_sideset_elems_loc_inds_and_ords(aSetName,tElemIndices,aSidesetOrdinals);

        // convert element indices to cell pointers
        moris::uint tNumCellsInSet = tElemIndices.numel();
        aCells = moris::Cell< mtk::Cell const * >(tNumCellsInSet);

        for(moris::uint i = 0 ; i < tNumCellsInSet; i++)
        {
            aCells(i) = &this->get_mtk_cell(tElemIndices(i));
        }
    }


    /*
     * returns the number of faces in a side set.
     */
    virtual
    uint get_sidesets_num_faces( moris::Cell< moris_index > aSideSetIndex ) const
    {
        moris::uint tNumSideSetFaces = 0;

        moris::Cell<std::string> tSideSetsNames = this->get_set_names( this->get_facet_rank() );

        for( luint Ik=0; Ik < aSideSetIndex.size(); ++Ik )
        {
            // get the treated sideset name
            std::string tTreatedSideset = tSideSetsNames( aSideSetIndex ( Ik ) );

            // get the sideset face indices
            Matrix< IndexMat > tSideSetElementInd = this->get_set_entity_loc_inds( this->get_facet_rank(), tTreatedSideset );

            // add up the sideset number of faces
            tNumSideSetFaces = tNumSideSetFaces + tSideSetElementInd.numel();
        }

        return tNumSideSetFaces;
    }

    /*
     * Returns the vertices in a node set. Does not include the aura.
     */
    virtual
    moris::Cell<moris::mtk::Vertex const *>
    get_vertices_in_vertex_set_no_aura(std::string aSetName) const
    {
        MORIS_ERROR(0,"No default implementation of get_vertices_in_vertex_set");
        return moris::Cell<moris::mtk::Vertex const *> (0);
    }


    virtual
    enum EntityRank
    get_facet_rank() const
    {
        uint tSpatialDim = this->get_spatial_dim();
        if(tSpatialDim  == 1)
        {
            return EntityRank::NODE;
        }
        else if(tSpatialDim == 2)
        {
            return EntityRank::EDGE;
        }
        else if(tSpatialDim == 3)
        {
            return EntityRank::FACE;
        }
        else
        {
            MORIS_ASSERT(0,"Invalid Mesh dimension detected in get_facet_rank ");
            return EntityRank::INVALID;
        }

    }


    void
    get_mtk_cells( Matrix< IndexMat > aCellInds,
                   moris::Cell<moris::mtk::Cell const *> & aCells) const
    {
        aCells = moris::Cell< mtk::Cell const * >(aCellInds.numel());

        for(moris::uint i = 0 ; i < aCellInds.numel(); i++)
        {
            aCells(i) = &this->get_mtk_cell(aCellInds(i));
        }
    }
};
}
}


#endif /* PROJECTS_MTK_SRC_CL_MTK_MESH_CORE_HPP_ */

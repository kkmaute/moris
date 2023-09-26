/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh.hpp
 *
 */

#ifndef SRC_MESH_CL_MTK_MESH_HPP_
#define SRC_MESH_CL_MTK_MESH_HPP_

#include <memory>
#include "typedefs.hpp" //MRS/COR/src
#include "cl_MTK_Enums.hpp"
#include "MTK_Tools.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_Map.hpp"
#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Cell.hpp" //MTK/src

#include "cl_MTK_Mesh_Core.hpp"

namespace moris
{

    namespace hmr
    {
        class Database;
    }

    namespace mtk
    {
//------------------------------------------------------------------------------

        /*!
         * The MTK mesh class is virtual but not pure virtual to allow for
         * partial implementation of the class depending on restrictions/abilities
         * of various libraries. All pure virtual functions appear at the top of the
         * sections. The functions are loosely grouped in this header file as follows:
         *  -  1.) General mesh information access
         *  -  2.) Access mesh connectivity
         *         2.a) standard connectivity access by index
         *         2.b) geometric connectivity access by index
         *  -  3.) Access mesh data by global ids functions
         *  -  4.) Coordinate Field Functions
         *  -  5.) Entity Ownership Functions
         *  -  6.) Mesh Sets Functions
         *  -  7.) Field Functions
         *  -  8.) Face Cluster Functions
         *  -  9.) Cell and Vertex Pointer Functions
         *  - 10.) Outputting Mesh
         *  Use this->shared_from_this() to create a shared
         * pointer of this ( the Mesh )
         */
        class Mesh_OLD : public std::enable_shared_from_this< Mesh_OLD >
        {
        public :
            // Verbose flag
            bool mVerbose = false;

            /**
             * trivial constructor
             */
            Mesh_OLD(){};

            //------------------------------------------------------------------------------

            /**
             * virtual destructor
             */
            virtual
            ~Mesh_OLD(){};

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
            //FIXME: IMPLEMENT THIS FUNCTION IN STK,XTK
            /*
             * Get number of B-Spline coefficients
             */
            virtual uint
            get_max_num_coeffs_on_proc(const uint aOrder) const
            {
                MORIS_ERROR( false, "get_max_num_coeffs_on_proc() not implemented for this mesh" );
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
            get_coefficient_indices_of_node(
                    const moris_index aNodeIndex,
                    const EntityRank  aBSplineRank )
                    {
                MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
                return Matrix<IndexMat>(0,0);
                    }

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

            uint get_sidesets_num_faces( moris::Cell< moris_index > aSidesetOrdinalList ) const
            {
                moris::uint tNumSideSetFaces = 0;

                moris::Cell<std::string> tSideSetsNames = this->get_set_names( EntityRank::FACE );

                for( luint Ik=0; Ik < aSidesetOrdinalList.size(); ++Ik )
                {
                    // get the treated sideset name
                    std::string tTreatedSideset = tSideSetsNames( aSidesetOrdinalList ( Ik ) );

                    // get the sideset face indices
                    Matrix< IndexMat > tSideSetElementInd = this->get_set_entity_loc_inds( EntityRank::FACE, tTreatedSideset );

                    // add up the sideset number of faces
                    tNumSideSetFaces = tNumSideSetFaces + tSideSetElementInd.numel();
                }

                return tNumSideSetFaces;
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

            /*
             * Get number of basis functions. For Lagrange meshes, the number of basis functions and the number of nodes
             * are equivalent. Therefore, a default implementation using get_num_nodes() is used here.
             */
            virtual
            uint
            get_num_basis_functions()
            {
                return get_num_nodes();
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

            //FIXME: change to pure virtual
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

            //------------------------------------------------------------------------------

            /**
             * @brief Get the neighboring elements connected to a given element through a given side ordinal
             *
             * @param aElementIndex element wrt. which the neighbors are to be determined
             * @param aSideOrdinal index of side ordinal relative to the given element
             * @param aMyRefineLevel
             * @param aMyOctreePosition
             * @param aNeighborElements output: list of neighboring elements connected through the side ordinal
             * @param aNeighborSideOrdinals output: list of neighboring element's side ordinals
             * @param aTransitionLocations output: list of the transition locations for the facet connectivity
             * @param aNeighborRefinementLevels
             */
            virtual void
            get_elements_connected_to_element_through_face_ord(
                    moris_index                 aElementIndex,
                    moris_index                 aSideOrdinal,
                    moris_index&                aMyRefineLevel,
                    moris_index&                aMyOctreePosition,
                    moris::Cell< moris_index >& aNeighborElements,
                    moris::Cell< moris_index >& aNeighborSideOrdinals,
                    moris::Cell< moris_index >& aTransitionLocations,
                    moris::Cell< moris_index >& aNeighborRefinementLevels ) const
            {
                MORIS_ERROR( false,
                    "mtk::Mesh::get_elements_connected_to_element_through_face_ord() - "
                    "Entered virtual function in Mesh base class; function is not implemented." );
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
             /*
              * Get elements interpolated into by a basis function. For a Lagrange mesh,
              * the elements in support of basis is equivalent to the elements connected
              * to a node. Therefore, a call to get_elements
              */
             virtual
             Matrix< IndexMat >
             get_elements_in_support_of_basis(moris_index aBasisIndex)
             {
                 return get_entity_connected_to_entity_loc_inds(aBasisIndex, EntityRank::NODE, EntityRank::ELEMENT);
             }

             //##############################################
             // 2.a.) Access geometric mesh data
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
                 MORIS_ERROR(0,"Entered virtual function in Mesh base class, (get_facet_ordinal_from_cell_and_facet_id_loc_inds is not implemented)");
                 return 0;
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

             /*
              * Get coordinate of a node
              */
             virtual
             Matrix< DDRMat >
             get_base_node_coordinate( moris_index aBaseNodeIndex ) const
             {
                 return this->get_node_coordinate(aBaseNodeIndex);
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

             //##############################################
             // Mesh Sets Access
             //##############################################

             virtual
             moris::Cell<std::string>
             get_set_names(enum EntityRank aSetEntityRank) const
             {
                 MORIS_ERROR(0," get_set_names has no base implementation");
                 return moris::Cell<std::string>(0);
             }

             virtual
             Matrix< IndexMat >
             get_set_entity_loc_inds( enum EntityRank aSetEntityRank,
                                      std::string     aSetName) const
             {
                 MORIS_ERROR(0," get_set_entity_ids has no base implementation");
                 return Matrix< IndexMat >(0,0);
             }

             virtual
             void
             get_sideset_elems_loc_inds_and_ords(
                     const  std::string     & aSetName,
                     Matrix< IndexMat >     & aElemIndices,
                     Matrix< IndexMat >     & aSidesetOrdinals ) const
             {
                 MORIS_ERROR(0," get_sideset_elems_loc_inds_and_ords has no base implementation");
             }

             virtual
             void
             get_sideset_cells_and_ords(
                     const  std::string & aSetName,
                     moris::Cell< mtk::Cell * > & aCells,
                     Matrix< IndexMat > &       aSidesetOrdinals )
             {
                 Matrix<IndexMat> tCellInds;

                 this->get_sideset_elems_loc_inds_and_ords(aSetName, tCellInds,aSidesetOrdinals);

                 aCells.resize(tCellInds.numel());

                 // iterate through cell inds and get cell ptrs
                 for(moris::uint i = 0; i <tCellInds.numel(); i++)
                 {
                     aCells(i) = &this->get_mtk_cell(tCellInds(i));
                 }
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
             // Face Cluster Access
             //##############################################

             /*
              * Is this face a member of a face cluster?
              */
             virtual
             bool
             has_face_cluster_membership(moris_index aFaceIndex) const
             {
                 MORIS_ERROR(0,"Entered virtual function in Mesh base class, (has_face_cluster_membership is not implemented)");
                 return 0;
             }

             //------------------------------------------------------------------------------

             /*
              * Is this face a parent of a face cluster?
              */
             virtual
             bool
             is_face_cluster_a_parent(moris_index aFaceIndex)
             {
                 MORIS_ERROR(0,"Entered virtual function in Mesh base class, (has_face_cluster_membership is not implemented)");
                 return 0;
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

             //##############################################
             // For FEM
             //##############################################

             virtual mtk::Cell  &
             get_writable_mtk_cell( moris_index aElementIndex )
             {
                 MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
                 return *mDummyCells;
             }

             void
             virtual get_adof_map( const uint aOrder,
                                   map< moris_id, moris_index > & aAdofMap ) const
             {
                 MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
             }
             //--------------------------------------------------------------

             virtual
             moris_id
             get_max_entity_id( enum EntityRank aEntityRank ) const
             {
                 MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
                 return 0;
             }

             //------------------------------------------------------------------------------

            //FIXME: THIS FUNCTION DESCRIPTION NEEDS TO BE IMPROVED
            //FIXME: Also, a unit test (not clear what needs to be provided)
            /**
             * provides a moris::Mat<uint> containing the IDs this mesh has
             * to communicate with
             *
             */
            virtual Matrix< IdMat >
            get_communication_table() const = 0;

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
            virtual moris_index
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
            virtual uint
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
            std::shared_ptr< hmr::Database >
            get_HMR_database( )
            {
                MORIS_ERROR( this->get_mesh_type() == MeshType::HMR ,"Not HMR " );
                return mDatabase;
            }

        protected:
            // Note these members are here only to allow for throwing in
            // get_mtk_cell and get_mtk_vertex function
            mtk::Vertex*     mDummyVertex;
            mtk::Cell*       mDummyCells;
            real             mDummyReal = 0.0;
            Matrix<DDRMat>   mDummyMatrix;

            //! ref to hmr object
            std::shared_ptr< hmr::Database > mDatabase;

        };
    } /* namespace mtk */
} /* namespace moris */

#endif /* SRC_MESH_CL_MTK_MESH_HPP_ */


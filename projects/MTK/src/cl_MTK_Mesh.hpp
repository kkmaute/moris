/*
 * cl_MTK_Mesh.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: messe
 */

#ifndef SRC_MESH_CL_MTK_MESH_HPP_
#define SRC_MESH_CL_MTK_MESH_HPP_

#include <memory>
#include "typedefs.hpp" //MRS/COR/src
#include "cl_MTK_Block.hpp" //MTK/src
#include "cl_Mesh_Enums.hpp"

namespace moris
{

    namespace mtk
    {
//------------------------------------------------------------------------------

        /**
         * the Mesh class. Use this->shared_from_this() to create a shared
         * pointer of this ( the Mesh )
         */
        class Mesh : public std::enable_shared_from_this< Mesh >
        {
//------------------------------------------------------------------------------
        public :
//------------------------------------------------------------------------------

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
            // General mesh information access
            //##############################################

            virtual
            uint
            get_spatial_dim() const
            {
                MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
                return 0;
            }

//------------------------------------------------------------------------------
            /*
             * Get number of entities for specified rank
             */
            virtual
            uint
            get_num_entities(
                    enum EntityRank aEntityRank) const
            {
                MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
                return 0;
            }
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

            /*
             * Get number of B-Spline coefficients
             */
            virtual uint
            get_num_coeffs() const
            {
                MORIS_ERROR( false, "get_num_coeffs() not implemented for this mesh" );
                return 0;
            }

//------------------------------------------------------------------------------
            //##############################################
            // Access Mesh Data by index Functions
            //##############################################
            /*
             * Generic get local index of entities connected to
             * entity using an entities local index
             */
            virtual
            Matrix<IndexMat>
            get_entity_connected_to_entity_loc_inds(moris_index     aEntityIndex,
                                                    enum EntityRank aInputEntityRank,
                                                    enum EntityRank aOutputEntityRank) const
            {
                MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
                return Matrix<IndexMat>(0,0);
            }
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
            get_elements_connected_to_element_and_face_ind_loc_inds(moris_index aElementIndex) const
            {
                MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
                return Matrix<IndexMat>(0,0);
            }

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
              * Get number of basis functions. For Lagrange meshes, the number of basis functions and the number of nodes
              * are equivalent. Therefore, a default implementation using get_num_nodes() is used here.
              */

             virtual
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
             virtual
             Matrix< IndexMat >
             get_elements_in_support_of_basis(moris_index aBasisIndex)
             {
                 return get_entity_connected_to_entity_loc_inds(aBasisIndex, EntityRank::NODE, EntityRank::ELEMENT);
             }




             //------------------------------------------------------------------------------
             //##############################################
             // global id functions
             //##############################################
             //------------------------------------------------------------------------------
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
//                 MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
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
             virtual
             Matrix< DDRMat >
             get_entity_field_value_real_scalar(Matrix< IndexMat > const & aEntityIndices,
                                                std::string        const & aFieldName,
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
             add_mesh_field_real_scalar_data_loc_inds(std::string     const & aFieldName,
                                                      enum EntityRank const & aFieldEntityRank,
                                                      Matrix<DDRMat>  const & aFieldData)
             {
                 MORIS_ERROR(0,"Entered virtual function in Mesh base class, (add_mesh_field_real_scalar_data_loc_inds is not implemented)");

             }

//------------------------------------------------------------------------------
             //##############################################
             // Cell and Vertex Pointer Functions
             //##############################################
             /*
              * Returns a reference to a cell in the mesh
              */
             virtual
             const mtk::Cell  &
             get_mtk_cell( moris_index aElementIndex)
             {
                 MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
                 return *mDummyCells;
             }

             /*
              * Returns a reference to a vertex in the mesh
              */
             virtual
             mtk::Vertex const &
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
             virtual get_adof_map( map< moris_id, moris_index > & aAdofMap ) const
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
            //FIXME: Also, a unit test (not clear what STK needs to provide)
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
                    std::string  &aFileName )
            {
                MORIS_ERROR(0,"Create output mesh not implemented");
            }
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------



        private:
            // Note these members are here only to allow for throwing in
            // get_mtk_cell and get_mtk_vertex function
            mtk::Vertex* mDummyVertex;
            mtk::Cell* mDummyCells;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */



#endif /* SRC_MESH_CL_MTK_MESH_HPP_ */

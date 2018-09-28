/*
 * cl_MTK_Mesh.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: messe
 */

#ifndef SRC_MESH_CL_MTK_MESH_HPP_
#define SRC_MESH_CL_MTK_MESH_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_MTK_Block.hpp" //MTK/src
#include "cl_Mesh_Enums.hpp"

namespace moris
{
// timestep global for exodus output
double gStkTimeStep = 0.0;

    namespace mtk
    {
//------------------------------------------------------------------------------
        class Mesh
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
             * Get number of faces
             */
            virtual
            uint
            get_num_elems() const
            {
                return get_num_entities(EntityRank::ELEMENT);
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
             * @param[out] Element to element connectivity and face ordinal shared
             *                   (where elements are all by index)
             */
            virtual
            Matrix< IndexMat >
            get_elements_connected_to_element_loc_inds(moris_index aElementIndex) const
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
                 MORIS_ERROR(0,"Entered virtual function in Mesh base class, (function is not implemented)");
                 return Matrix<IdMat>(0,0);
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
//------------------------------------------------------------------------------
             //##############################################
             // Cell and Vertex Pointer Functions
             //##############################################
             /*
              * Returns a reference to a cell in the mesh
              */
             virtual
             mtk::Cell const &
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
//------------------------------------------------------------------------------
             //##############################################
             //  Access block information
             //##############################################


            /**
             * returns the number of blocks on this mesh
             */
            virtual uint
            get_number_of_blocks() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns a pointer to a block
             */
            virtual Block *
            get_block_by_index( const moris_index & aIndex ) = 0;

//------------------------------------------------------------------------------

            /**
             * returns a pointer to a block ( const version )
             */
            const virtual Block *
            get_block_by_index( const moris_index & aIndex ) const = 0;

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

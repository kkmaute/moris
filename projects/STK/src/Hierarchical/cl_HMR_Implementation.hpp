/*
 * cl_HMR_Implementation.hpp
 *
 *  Created on: Feb 9, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HMR_IMPLEMENTATION_HPP_
#define SRC_MESH_CL_HMR_IMPLEMENTATION_HPP_
#include "cl_Hierarchical_Mesh_Main.hpp" // STK/src/Hierarchical

namespace moris
{
    class HMR_Implementation : public database
    {
    public:
        HMR_Implementation( Hierarchical_Mesh_Main & aHMRMesh ):
            mHMRMesh( aHMRMesh )
    {
    }

        ~HMR_Implementation()
        {

        }

        //Data from current proc ---------------------------------------------------------------------------------------------------------------------------------------------
        /**
         * Computes the number of nodes,edges,faces or elements  on the current processor, which are owned by this proc
         *
         * @param[out] tNumOwnedEntityrank   ............  Number of entities, which are owned by current proc
         *
         */
        uint
        get_num_nodes_current_proc() const
        {
            return mHMRMesh.get_num_entities_owned_current_proc( EntityRank::NODE );
        }

        uint
        get_num_edges_current_proc() const
        {
            return mHMRMesh.get_num_entities_owned_current_proc( EntityRank::EDGE );
        }

        uint
        get_num_faces_current_proc() const
        {
            return mHMRMesh.get_num_entities_owned_current_proc( EntityRank::FACE );
        }

        uint
        get_num_elems_current_proc() const
        {
            return mHMRMesh.get_num_entities_owned_current_proc( EntityRank::ELEMENT );
        }

        /**
         * Provides a list of entities, which are owned by the current proc
         *
         * @param[out] tOwnedEntityrank   ............  A list in a vector with entities, which are owned by the current processor
         *
         */
        Matrix< DDUMat >
        get_entities_owned_current_proc(
                enum EntityRank  aEntityRank) const
                {
            return mHMRMesh.get_entities_owned_current_proc( aEntityRank );
                }

        //Universal data (All data, that the current proc knows) --------------------------------------------------------------------------
        /**
         * Computes the number of nodes,edges,faces or elements  on the current processor and the aura (independent of ownership)
         *
         * @param[out] tNumEntityrank   ............  Number of entities, which are on the current processor and their aura
         *
         */
        uint
        get_num_nodes() const
        {
            return mHMRMesh.get_num_entities_universal( EntityRank::NODE );
        }

        uint
        get_num_edges() const
        {
            return mHMRMesh.get_num_entities_universal( EntityRank::EDGE );
        }

        uint
        get_num_faces() const
        {
            return mHMRMesh.get_num_entities_universal( EntityRank::FACE );
        }

        uint
        get_num_elems() const
        {
            return mHMRMesh.get_num_entities_universal( EntityRank::ELEMENT );
        }

        //Data from current proc with aura------------------------------------------------------------------------------------------
        /**
         * Provides a list of entities, which are known by the current proc with the aura
         *
         * @param[out] tUniveralEntityrank   ............  A list in a vector with entities, which are known from the current processor and the aura
         *
         */
        uint
        get_num_entities_universal(
                enum EntityRank  aEntityRank) const
        {
            return mHMRMesh.get_num_entities_universal( aEntityRank );
        }

        /**
         * Provides a list of entities, which are known by the current proc with the aura
         *
         * @param[out] tUniveralEntityrank   ............  A list in a vector with entities, which are known from the current processor and the aura
         *
         */
        Matrix< DDUMat >
        get_entities_universal(
                enum EntityRank  aEntityRank) const
                {
            return mHMRMesh.get_entities_universal( aEntityRank );
                }

        //Data from current proc, but independent of ownership------------------------------------------------------------------------------------------
        /**
         * Computes the number of nodes,edges,faces or elements  on the current processor, without the aura (inedpendent of ownership)
         *
         * @param[out] tNumEntityrank   ............  Number of entities, which are in the domain of the owned elements, (independen of ownership and without aura)
         *
         */
        uint
        get_num_entities_locally_owned_globally_shared(
                enum EntityRank  aEntityRank) const
        {
            return mHMRMesh.get_num_entities_owned_and_shared_by_current_proc( aEntityRank );
        }

        /**
         * Provides a list of entities, which are on the current processor, without the aura (independent of ownership)
         *
         * @param[out] tEntityrank   ............  A list in a vector with entities, which are on the current processor, without aura
         *
         *
         */
        Matrix< DDUMat >
        get_entities_owned_and_shared_by_current_proc(
                enum EntityRank  aEntityRank) const
                {
            return mHMRMesh.get_entities_owned_and_shared_by_current_proc( aEntityRank );
                }

        //Data of aura from current proc------------------------------------------------------------------------------------------------------------------
        /**
         * Computes the number of nodes,edges,faces or elements  on the current processor, which are living in the aura ( Universal_entities minus entities_owned_and_shared )
         *
         * @param[out] tNumEntityrank   ............  Number of entities, which are in the aura of the current processor ( Universal_entities minus entities_owned_and_shared )
         *
         */
        uint
        get_num_entities_aura(
                enum EntityRank  aEntityRank) const
        {
            return mHMRMesh.get_num_entities_in_aura( aEntityRank );
        }

        /**
         * Provides a list of entities, which are in the aura on the current processor ( Universal_entities minus entities_owned_and_shared )
         *
         * @param[out] tEntityrank   ............  A list in a vector with entities, which are in the aura on the current processor ( Universal_entities minus entities_owned_and_shared )
         *
         */
        Matrix< DDUMat >
        get_entities_in_aura(
                enum EntityRank  aEntityRank) const
                {
            return mHMRMesh.get_entities_in_aura( aEntityRank );
                }

        //Data information about edges, faces, nodes, (Mesh database wants to provide an element, face, edge ID, but this is independent for HMR--------------------------------------------
        /**
         * Computes the number of faces. This information is independent for the HMR mesh, but can be dependent for an arbitrary mesh
         *
         * @param[out] tNumTopology   ............  Number of topology relations between element and face
         *
         */
        uint
        get_elem_topology_num_faces(
                uint aElemId ) const
        {
            MORIS_ASSERT( mHMRMesh.give_modeldim() >= 2 && mHMRMesh.give_modeldim() <=3, " Faces are only in 2D and 3D available" );
            return 2 * mHMRMesh.give_modeldim();
        };

        /**
         * Computes the number of edges. This information is independent for the HMR mesh, but can be dependent for an arbitrary mesh
         *
         * @param[out] tNumTopology   ............  Number of topology relations between element and edges
         *
         */
        uint
        get_elem_topology_num_edges(
                uint aElemId ) const;

        /**
         * Computes the number of nodes. This information is independent for the HMR mesh, but can be dependent for an arbitrary mesh
         *
         * @param[out] tNumTopology   ............  Number of topology relations between element and nodes
         *
         */
        uint
        get_elem_topology_num_nodes(
                uint aElemId ) const
        {
            return pow( mHMRMesh.give_polynomial() + 1, mHMRMesh.give_modeldim() );
        };

        /**
         * Computes the number of edges of a face. This information is independent for the HMR mesh, but can be dependent for an arbitrary mesh
         *
         * @param[out] tNumTopology   ............  Number of topology relations between face and edges
         *
         */
        uint
        get_face_topology_num_edges(
                uint aFaceId ) const;

        /**
         * Computes the number of faces. This information is independent for the HMR mesh, but can be dependent for an arbitrary mesh
         *
         * @param[out] tNumTopology   ............  Number of topology relations between faces and nodes
         *
         */
        uint
        get_face_topology_num_nodes(
                uint aFaceId ) const
        {
            MORIS_ASSERT( mHMRMesh.give_modeldim() >= 2 && mHMRMesh.give_modeldim() <=3, " Edges are only in 2D and 3D available" );
            return pow( mHMRMesh.give_polynomial() + 1, mHMRMesh.give_modeldim() - 1 );
        };

        /**
         * Computes the number of nodes. This information is independent for the HMR mesh, but can be dependent for an arbitrary mesh
         *
         * @param[out] tNumTopology   ............  Number of topology relations between edges and nodes
         *
         */
        uint
        get_edge_topology_num_nodes(
                uint aEdgeId ) const
        {
            MORIS_ASSERT( mHMRMesh.give_modeldim() >= 2 && mHMRMesh.give_modeldim() <=3, " Edges are only in 2D and 3D available" );
            return mHMRMesh.give_polynomial() + 1;
        };
        //Topology information between elements, edges, faces, nodes ----------------------------------------------------
        /**
         * Computes the neighbor elements, which are connected between faces (no diagonal neighbors are possible)
         *
         * @param[in] aElemId ......................... Element Id
         *
         * @param[out] tNeighborElements   ............  A list of elements in a matrix, which are connected to the current element "aElemId", First column: Element ID, second column: ordinal
         *
         */
        Matrix< DDUMat >
        get_elements_connected_to_element(
                uint const aElemId ) const
                {
            return mHMRMesh.give_active_face_neighbor_of_element( aElemId );
                };

        /**
         * Computes the elements, which are connected to a face
         *
         * @param[in] aFaceId ......................... Face Id
         *
         * @param[out] tElements   ............  A list of elements, which are connected to a face
         *
         */
        Matrix< DDUMat > get_elements_connected_to_face(
                uint const aFaceId ) const
                        {
            return mHMRMesh.give_active_elements_of_face( aFaceId );
                        }

        /**
         * Computes the elements, which are connected to a face
         *
         * @param[in] aEdgeId .................  Edgecd Id
         *
         * @param[out] tElements   ............  A list of elements, which are connected to a face
         *
         */
        Matrix< DDUMat >
        get_elements_connected_to_edge(
                uint const aEdgeId ) const
                {
            return mHMRMesh.give_active_elements_of_edge( aEdgeId );
                };

        /**
         * Computes the elements, which are connected to a face
         *
         * @param[in] aEdgeId .................  Edgecd Id
         *
         * @param[out] tElements   ............  A list of elements, which are connected to a face
         *
         */
        Matrix< DDUMat >
        get_elements_connected_to_node(
                uint const aNodeId ) const
                {
            return mHMRMesh.give_element_of_lagrange_basis( aNodeId );
                };

        /**
         * Computes the faces, which are connected to an element
         *
         * @param[in] aElementId .................  Element Id
         *
         * @param[out] tFaces   ............  A list of faces, which are connected to an element
         *
         */
        Matrix< DDUMat >
        get_faces_connected_to_element(
                uint const aElementId ) const;

        /**
         * Computes the edges, which are connected to an element
         *
         * @param[in] aElementId .................  Element Id
         *
         * @param[out] tEdges   ............  A list of edges, which are connected to an element
         *
         */
        Matrix< DDUMat >
        get_edges_connected_to_element(
                uint const aElementId ) const;

        /**
         * Computes the nodes, which are connected to an element
         *
         * @param[in] aElementId .................  Element Id
         *
         * @param[out] tNodes   ............  A list of nodes, which are connected to an element
         *
         */
        Matrix< DDUMat >
        get_nodes_connected_to_element(
                uint const aElementId ) const
                {
            return mHMRMesh.give_active_lagrange_basis_of_element( aElementId );
                };

        /**
         * Computes the nodes, which are connected to a face
         *
         * @param[in] aFaceId .................  Face Id
         *
         * @param[out] tNodes   ............  A list of nodes, which are connected to a face
         *
         */
        Matrix< DDUMat >
        get_nodes_connected_to_face(
                uint const aFaceId ) const
                {
            return mHMRMesh.give_face_active_lagrange_basis( aFaceId );
                };

        /**
         * Computes the nodes, which are connected to an edge
         *
         * @param[in] aEdgeId .................  Edge Id
         *
         * @param[out] tNodes   ............  A list of nodes, which are connected to an edge
         *
         */
        Matrix< DDUMat >
        get_nodes_connected_to_edge(
                uint const aEdgeId ) const
                {
            return mHMRMesh.give_edge_active_lagrange_nodes( aEdgeId );
                };

        /**
         * Provides the owner of an entity ID
         *
         * @param[in] aEntityID ......................... aEntity Id
         *
         * @param[out] tOwner   ............  Owner of the entity Id with respective entity rank
         *
         */
        uint
        parallel_owner_rank_by_entity_id(
                uint aEntityID,
                enum EntityRank aEntityRank) const
        {
            return mHMRMesh.parallel_owner_rank_by_entity_id( aEntityID, aEntityRank );
        };

        /**
         * Provides a list of entities, which share a specific ID
         *
         * @param[in] aEntityID ......................... aEntity Id
         *
         * @param[out] tShare   ............  List of vector, which share this ID
         *
         */
        Matrix< DDUMat >
        get_procs_sharing_entity_by_id(
                uint aEntityID,
                enum EntityRank aEntityRank) const
                {
            return mHMRMesh.get_procs_sharing_entity_by_id( aEntityID, aEntityRank );
                };

        /**
         * Provides a (3xn) matrix containing node coordinates
         *
         * @param[in] aNodeIds .............. vector containing node indices
         *
         *
         */
        Matrix< DDRMat >
        get_selected_nodes_coords(
                Matrix< DDUMat > aNodeIds ) const
        {
            Mat < real > aOutput(3, aNodeIds.length());
            Mat < real > tNodeCoords(3,1);

            for(moris::uint k=0; k<aNodeIds.length(); ++k)
            {
                tNodeCoords = mHMRMesh.give_coordinate_from_lagrange_basis(aNodeIds(k));
                for(moris::uint i=0; i<3; ++i)
                {
                    aOutput(i,k) = tNodeCoords(i);
                }
            }

            return aOutput;
        }

        /* @brief: converts an entity ID to a local one on the proc
         *
         * @param[in] aEntityRank   type of entity
         * @param[in] aID           gobal ID of the entity
         */
        /* uint
        get_local_id_on_proc_from_global_id( const EntityRank   aEntityRank , const moris::uint aId )
        const {

            switch ( aEntityRank )
            {
                case (EntityRank::ELEMENT):
                {
                    return mHMRMesh.give_local_element_id(aId);
                }
                case (EntityRank::NODE):
                {
                    return mHMRMesh.give_local_lagrange_basis_id(aId);
                }
                default:
                {
                    MORIS_LOG_ERROR << "Current Database does not support this functionality";
                    return 0;
                }
            }
        } */

        /**
         * Computes the nodes, which are connected to an element
         *
         * @param[in] aElementId .................  Element Id
         *
         * @param[out] tNodes   ............  A list of nodes, which are connected to an element
         *
         */
       /*  Matrix< DDUMat >
        get_local_nodes_on_proc_connected_to_local_element_on_proc(
                uint const aLocalElementId ) const
        {
            return mHMRMesh.give_local_active_lagrange_basis_of_element( aLocalElementId );
        }; */
    private:
        /* @TODO Discuss this with Keenan */
        const Hierarchical_Mesh_Main& mHMRMesh;
    };
}

#endif /* SRC_MESH_CL_HMR_IMPLEMENTATION_HPP_ */

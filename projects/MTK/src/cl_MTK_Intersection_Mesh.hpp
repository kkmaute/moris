/*
 * cl_MTK_Intersection_Mesh.hpp
 *
 *  Created on: Jul 30, 2021
 *      Author: momo
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_INTERSECTION_MESH_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_INTERSECTION_MESH_HPP_

#include "cl_MTK_Intersection_Detect.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

namespace moris
{
    namespace mtk
    {
        class Intersection_Mesh : public mtk::Integration_Mesh
        {
            private:
                moris::mtk::Integration_Mesh*    mBackGroundMesh;  /**< Background mesh that is not intersected */
                moris::mtk::Intersection_Detect* mIntersectionDetect; /**< Holds all the information and data created */

            public:
                // ----------------------------------------------------------------------------
                /*
                 * Default constructor
                 */
                Intersection_Mesh( moris::mtk::Integration_Mesh* aBackGroundMesh,
                        moris::mtk::Intersection_Detect* aIntersectionDetect );

                // ----------------------------------------------------------------------------
                /**
                 * Returns the spatial dimension of the mesh
                 */
                uint
                get_spatial_dim() const
                {
                    return mBackGroundMesh->get_spatial_dim();
                }

                // ----------------------------------------------------------------------------
                /*
                 * Get number of nodes/elements in the mesh
                 * @param[ in ] aEntityRank Rank of the entity
                 * @param[ in ] aIndex A mesh index
                 */
                moris::uint
                get_num_entities(enum EntityRank aEntityRank, const moris_index aIndex = 0) const;

                // ----------------------------------------------------------------------------
                /*
                 * adds side sets and block sets created in the process of intersection to the mesh
                 * @param[ in ] aEntityRank Rank of the entity
                 */
                moris::Cell< std::string >
                get_set_names(enum EntityRank aSetEntityRank) const;

                // ----------------------------------------------------------------------------
                /*
                 * Determines if a given side set name belongs to intersection data or background mesh
                 * @param[ in ] aSetName Name of the set
                 */
                bool
                is_intersection_side_set(const  std::string & aSetName) const;

                // ----------------------------------------------------------------------------

                /**
                 *Determines if a given block set name belongs to intersection data or background mesh
                 * @param[ in ] aSetName Name of the set
                 */
                bool
                is_intersection_block_set(const  std::string & aSetName) const;

                // ----------------------------------------------------------------------------
                /*
                 * get side set elements indices and their respective ordinal
                 * @param[ in ] aSetName Name of the set
                 * @param[ out ] aElemIndices Indices of the element in the set
                 * @param[ out ] aSidesetOrdinals Side ordinals of the element that make the side set
                 */
                void
                get_sideset_elems_loc_inds_and_ords(
                        const  std::string & aSetName,
                        Matrix< IndexMat > & aElemIndices,
                        Matrix< IndexMat > & aSidesetOrdinals ) const;

                // ----------------------------------------------------------------------------
                /*
                 * Get element indices in a block set
                 * @param[ in ] aSetIndex Index of the set
                 */
                Matrix<IndexMat>
                get_element_indices_in_block_set( uint aSetIndex );

                // ----------------------------------------------------------------------------
                /*
                 * Get element ids in a block set
                 * @param[ in ] aSetIndex Index of the set
                 */
                Matrix<IndexMat>
                get_element_ids_in_block_set( uint aSetIndex ) const;

                // ----------------------------------------------------------------------------
                /*
                 * Get the node coordinates for an index
                 * @param[ in ] aNodeIndex Index of the node
                 */
                moris::Matrix< moris::DDRMat >
                get_node_coordinate(moris_index aNodeIndex) const;

                // ----------------------------------------------------------------------------
                /*
                 * Get global entity id from entity local index
                 * @param[ in ] aEntityIndex Index of the entity
                 * @param[ in ] aEntityRank Rank of the entity
                 * @param[ in ] aMeshIndex Index of the mesh
                 */
                moris::moris_id
                get_glb_entity_id_from_entity_loc_index(moris_index   aEntityIndex,
                        enum EntityRank aEntityRank,
                        const moris_index aMeshIndex ) const;

                // ----------------------------------------------------------------------------
                /**
                 * Returns block set toplogy
                 * @param[ in ] aSetName Name of the set
                 */
                enum CellTopology get_blockset_topology(const std::string & aSetName);

                // ----------------------------------------------------------------------------
                /**
                 * Get nodes connected to an element
                 * @param[ in ] aElementIndex An element index
                 */
                Matrix< IndexMat >
                get_nodes_connected_to_element_loc_inds(moris_index aElementIndex) const;

                //##############################################
                // Pure virtual function ( not implemented )
                //##############################################

                // ----------------------------------------------------------------------------
                /**
                 * Get nodes connected to an element
                 * @param[ in ] aInterpCell An Interpolation cell
                 */
                Cell_Cluster const &
                get_cell_cluster(Cell const & aInterpCell) const
                {
                    MORIS_ERROR( false, "Not Implemented for Intersection Mesh");
                    mtk::Cell_Cluster*       tDummyCells = nullptr;
                    return *tDummyCells;
                }

                // ----------------------------------------------------------------------------

                /**
                 * Get block set names
                 */
                moris::Cell<std::string>
                get_block_set_names() const
                {
                    MORIS_ERROR( false, "Not Implemented for Intersection Mesh");
                    return moris::Cell< std::string > (0);
                }

                // ----------------------------------------------------------------------------

                /**
                 * Get cell clusters in set
                 * @param[ in ] aBlockSetOrdinal Ordinal of the block set
                 */
                moris::Cell<Cluster const *>
                get_cell_clusters_in_set(moris_index aBlockSetOrdinal) const
                {
                    MORIS_ERROR( false, "Not Implemented for Intersection Mesh");
                    return moris::Cell< Cluster const * >(0);
                }

                // ----------------------------------------------------------------------------
                /**
                 * Get side clusters within a side set
                 * @param[ in ] aSideSetOrdinal Ordinal of the side set
                 */
                moris::Cell<Cluster const *>
                get_side_set_cluster(moris_index aSideSetOrdinal) const
                {
                    MORIS_ERROR( false, "Not Implemented for Intersection Mesh");
                    return moris::Cell< Cluster const *>(0);
                }

                // ----------------------------------------------------------------------------
                /**
                 * get number of side sets
                 */
                uint
                get_num_side_sets() const
                {
                    MORIS_ERROR( false, "get_num_side_sets(), Not Implemented for Intersection Mesh");
                    return 0;
                }

                // ----------------------------------------------------------------------------
                /**
                 * Returns the label
                 * @param[ in ] aSideSetOrdinal Ordinal of the side set
                 */
                std::string
                get_side_set_label(moris_index aSideSetOrdinal) const
                {
                    MORIS_ERROR( false, " get_side_set_label, Not Implemented for Intersection Mesh");
                    return std::string();
                }

                // ----------------------------------------------------------------------------
                /*!
                 * Returns the index given a label
                 * @param[ in ] aSideSetLabel Label of the side set
                 */
                moris_index
                get_side_set_index(std::string aSideSetLabel) const
                {
                    MORIS_ERROR( false, "get_side_set_index, Not Implemented for Intersection Mesh");
                    return 0;
                }

                // ----------------------------------------------------------------------------
                /*!
                 * Returns the number of double sided side sets in the mesh
                 */
                uint
                get_num_double_sided_sets() const
                {
                    MORIS_ERROR( false, "get_num_double_sided_sets, Not Implemented for Intersection Mesh");
                    return 0;
                }

                // ----------------------------------------------------------------------------
                /*!
                 * Returns the label
                 * @param[ in ] aSideSetOrdinal Ordinal of the side set
                 */
                std::string
                get_double_sided_set_label(moris_index aSideSetOrdinal) const
                {
                    MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh");
                    return std::string();
                }

                // ----------------------------------------------------------------------------
                /*!
                 * Returns the label
                 */
                MeshType
                get_mesh_type() const
                {
                    MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh");
                    return  MeshType::END_ENUM;
                }


                // ----------------------------------------------------------------------------
                /*!
                 * Get local indices of the entities connected to a specific entity
                 * @param[ in ] aEntityIndex Entity index
                 * @param[ in ] aInputEntityRank Entity rank
                 * @param[ in ] aOutputEntityRank Entity rank of the output
                 * @param[ in ] aDiscretizationIndex A discretization index
                 */
                Matrix<IndexMat>
                get_entity_connected_to_entity_loc_inds(
                        moris_index        aEntityIndex,
                        enum EntityRank    aInputEntityRank,
                        enum EntityRank    aOutputEntityRank,
                        const moris_index  aDiscretizationIndex = 0) const
                {
                    MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh");
                    return Matrix<IndexMat>(0,0);
                }

                // ----------------------------------------------------------------------------
                /*!
                 * Returns node owner
                 * @param[ in ] aNodeIndex Index of the node
                 */
                uint get_node_owner(moris_index aNodeIndex) const
                {
                    MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh");
                    return 0;
                }

                // ----------------------------------------------------------------------------
                /*!
                 * Returns the element owner
                 * @param[ in ] aNodeIndex Index of the node
                 */
                uint get_element_owner(moris_index aNodeIndex) const
                {
                    MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh");
                    return 0;
                }

                // ----------------------------------------------------------------------------
                /*!
                 * get communication table
                 */
                Matrix< IdMat >
                get_communication_table() const
                {
                    MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh");
                    return Matrix<IndexMat>(0,0);
                }

                // ----------------------------------------------------------------------------
                /*!
                 * Returns integration block set shape
                 * @param[ in ] aSetName Name of the set
                 */
                enum CellShape
                get_IG_blockset_shape(const std::string & aSetName)
                {
                    MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh");
                    return CellShape::END_ENUM;
                }

                // ----------------------------------------------------------------------------
                /*!
                 * Returns interpolation block set shape
                 * @param[ in ] aSetName Name of the set
                 */
                enum CellShape
                get_IP_blockset_shape(const std::string & aSetName)
                {
                    MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh");
                    return CellShape::END_ENUM;
                }

                // ----------------------------------------------------------------------------
                /*!
                 * Get double sided cluster
                 * @param[ in ] aSideSetOrdinal Ordinal of the side set
                 */
                moris::Cell<Cluster const*>
                get_double_side_set_cluster(moris_index aSideSetOrdinal) const
                {
                    MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh");
                    return moris::Cell<Cluster const*>(0);
                }

                // ----------------------------------------------------------------------------
                /*!
                 * Returns the label
                 * @param[ in ] aElementIndex Index of the element
                 */
                Matrix< IndexMat >
                get_elements_connected_to_element_and_face_ind_loc_inds(moris_index aElementIndex) const
                {
                    MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh");
                    return Matrix< IndexMat >(0,0);
                }

                // ----------------------------------------------------------------------------
                /*!
                 * Returns the label
                 * @param[ in ]
                 */
                Matrix< IndexMat >
                get_set_entity_loc_inds(
                        enum EntityRank aSetEntityRank,
                        std::string     aSetName) const
                {
                    return Matrix< IndexMat >(0,0);
                }

        };
    }
}


#endif /* PROJECTS_MTK_SRC_CL_MTK_INTERSECTION_MESH_HPP_ */

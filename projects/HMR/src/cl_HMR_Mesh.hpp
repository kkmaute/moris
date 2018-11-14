/*
 * cl_HMR_Mesh.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_INTERFACE_HPP_
#define SRC_HMR_CL_HMR_INTERFACE_HPP_

#include "cl_Cell.hpp" //CON/src
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Mesh.hpp" //MTK/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp"

namespace moris
{

    namespace hmr
    {
//-------------------------------------------------------------------------------

        // forward declaration of HMR Field object
        class Field;

        class Database;

//-------------------------------------------------------------------------------
        /**
         * \brief mesh interface class
         */
        class Mesh : public mtk::Mesh
        {

            //! ref to hmr object
            std::shared_ptr< Database > mDatabase;

            //! describing label
            std::string mLabel;

            Lagrange_Mesh_Base * mMesh;
//-------------------------------------------------------------------------------
        public:
//-------------------------------------------------------------------------------

            /**
             * mesh constructor, to be called from HMR
             */
            Mesh( std::shared_ptr< Database > aDatabase,
                    const uint & aLagrangeOrder,
                    const uint & aActivationPattern );

//-------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Mesh();

//-------------------------------------------------------------------------------

            /**
             * return the type of this mesh
             */
            MeshType
            get_mesh_type() const
            {
                return MeshType::HMR;
            }

//-------------------------------------------------------------------------------

            /**
             * provides a moris::Matrix< DDUMat > containing the IDs this mesh has
             * to communicate with
             */
            Matrix< IdMat >
            get_communication_table() const ;

//-------------------------------------------------------------------------------

            /**
             * creates a new field pointer that is linked to this mesh
             */
            std::shared_ptr< Field >
            create_field( const std::string & aLabel, const uint & aBSplineOrder );

//-------------------------------------------------------------------------------

            /**
             * returns a pointer to the underlying lagrange mesh
             */
            Lagrange_Mesh_Base *
            get_lagrange_mesh()
            {
                return mMesh;
            }

//-------------------------------------------------------------------------------

            /**
             * return a pointer to the database
             */
            std::shared_ptr< Database >
            get_database()
            {
                return mDatabase;
            }

//-------------------------------------------------------------------------------
// Functions for MTK
//-------------------------------------------------------------------------------

            /**
             * returns the number of dimensions of this mesh
             */
            uint
            get_spatial_dim() const;

//-------------------------------------------------------------------------------

            uint
            get_num_entities(
                    enum EntityRank aEntityRank) const;

//-------------------------------------------------------------------------------

            uint
            get_num_nodes() const;

//-------------------------------------------------------------------------------

            uint
            get_num_edges() const;

//-------------------------------------------------------------------------------

            uint
            get_num_faces() const;

//-------------------------------------------------------------------------------

            uint
            get_num_elems() const;

//-------------------------------------------------------------------------------

            uint
            get_num_coeffs( const uint aOrder ) const;

//-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_bspline_inds_of_node_loc_ind(
                    const moris_index aNodeIndex,
                    const EntityRank  aBSplineRank );

//-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_entity_connected_to_entity_loc_inds(
                                       moris_index     aEntityIndex,
                                       enum EntityRank aInputEntityRank,
                                       enum EntityRank aOutputEntityRank) const;

//-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_nodes_connected_to_node_loc_inds( moris_index aNodeIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_nodes_connected_to_edge_loc_inds( moris_index aEdgeIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_nodes_connected_to_face_loc_inds( moris_index aFaceIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_nodes_connected_to_element_loc_inds( moris_index aElementIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix < IndexMat >
            get_edges_connected_to_node_loc_inds( moris_index aNodeIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_edges_connected_to_element_loc_inds( moris_index aElementIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix < IndexMat >
            get_faces_connected_to_node_loc_inds( moris_index aNodeIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_faces_connected_to_element_loc_inds( moris_index aElementIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix < IndexMat >
            get_elements_connected_to_node_loc_inds( moris_index aNodeIndex ) const;

//-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_elements_connected_to_face_loc_inds( moris_index aFaceIndex ) const ;

//-------------------------------------------------------------------------------

            /*
             * Since the connectivity between entities of the same rank are considered
             * invalid by STK standards, we need a separate function for element to element
             * specifically
             *      *
             * @param[in]  aElementId - element id
             * @param[out] Element to element connectivity and face ordinal shared
             *                   (where elements are all by index)
             */
            Matrix< IndexMat >
            get_elements_connected_to_element_loc_inds( moris_index aElementIndex ) const;

//-------------------------------------------------------------------------------
//          Global ID Functions
//-------------------------------------------------------------------------------

            moris_id
            get_glb_entity_id_from_entity_loc_index(
                    moris_index     aEntityIndex,
                    enum EntityRank aEntityRank) const ;
//-------------------------------------------------------------------------------

            moris_id
            get_max_entity_id( enum EntityRank aEntityRank ) const ;

//-------------------------------------------------------------------------------
//          Coordinate Field Functions
//-------------------------------------------------------------------------------

            Matrix< DDRMat >
            get_node_coordinate( moris_index aNodeIndex ) const;

//-------------------------------------------------------------------------------
//           Entity Ownership Functions
//-------------------------------------------------------------------------------

            moris_id
            get_entity_owner(  moris_index     aEntityIndex,
                    enum EntityRank aEntityRank ) const;

//-------------------------------------------------------------------------------
//           Set Functions
//-------------------------------------------------------------------------------

            void
            get_sideset_elems_loc_inds_and_ords(
                    const  std::string     & aSetName,
                    Matrix< IndexMat >     & aElemIndices,
                    Matrix< IndexMat >     & aSidesetOrdinals );

//-------------------------------------------------------------------------------
//           Pointer Functions for FEM
//-------------------------------------------------------------------------------

            const  mtk::Vertex &
            get_mtk_vertex( moris_index aVertexIndex )
            {
                return *mMesh->get_node_by_index( aVertexIndex );
            }

//-------------------------------------------------------------------------------

            const  mtk::Cell &
            get_mtk_cell( moris_index aElementIndex )
            {
                return *mMesh->get_element( aElementIndex );
            }

//-------------------------------------------------------------------------------

            mtk::Cell  &
            get_writable_mtk_cell( moris_index aElementIndex )
            {
                return *mMesh->get_element( aElementIndex );
            }

//-------------------------------------------------------------------------------

            void
            get_adof_map(
                    const uint aOrder,
                    map< moris_id, moris_index > & aAdofMap ) const;

//-------------------------------------------------------------------------------

            moris_index
            get_field_ind(
                    const std::string     & aFieldLabel,
                    const enum EntityRank   aEntityRank  ) const;

//-------------------------------------------------------------------------------

            uint
            get_num_fields( const enum EntityRank aEntityRank ) const;

//-------------------------------------------------------------------------------

            real &
            get_value_of_scalar_field(
                    const      moris_index  aFieldIndex,
                    const enum EntityRank   aEntityRank,
                    const uint              aEntityIndex );

//-------------------------------------------------------------------------------

            const real &
            get_value_of_scalar_field(
                    const      moris_index  aFieldIndex,
                    const enum EntityRank   aEntityRank,
                    const uint              aEntityIndex ) const;

//-------------------------------------------------------------------------------

            Matrix<DDRMat> &
            get_field( const moris_index     aFieldIndex,
                       const enum EntityRank aEntityRank );

//-------------------------------------------------------------------------------
            private:
//-------------------------------------------------------------------------------

            void
            get_element_indices_from_memory_indices(
                    const Matrix< DDLUMat>      & aMemoryIndices,
                          Matrix< IndexMat >    & aIndices ) const;

//-------------------------------------------------------------------------------
            /**
             * subroutine for get_elements_connected_to_element_loc_inds(
             */
            void
            collect_memory_indices_of_active_element_neighbors(
                    const moris_index  aElementIndex,
                    Matrix< DDLUMat> & aMemoryIndices,
                    luint            & aCounter ) const;

//-------------------------------------------------------------------------------

            /**
             * subroutine for get_elements_connected_to_node_loc_inds
             */
            void
            collect_memory_indices_of_active_elements_connected_to_node(
                    const moris_index  aNodeIndex,
                    Matrix< DDLUMat> & aMemoryIndices ) const;

//-------------------------------------------------------------------------------

            const Matrix< DDRMat > &
            get_t_matrix_of_node_loc_ind(
                    const moris_index aNodeIndex,
                    const EntityRank  aBSplineRank )
            {
                return *mMesh->get_node_by_index(
                            aNodeIndex )->get_interpolation(
                                entity_rank_to_order( aBSplineRank ) )
                                ->get_weights();
            }

//-------------------------------------------------------------------------------
        };

    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_INTERFACE_HPP_ */

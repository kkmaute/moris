/*
 * cl_HMR_Mesh.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_CL_HMR_INTERFACE_HPP_
#define SRC_HMR_CL_HMR_INTERFACE_HPP_

#include "cl_HMR_Lagrange_Mesh_Base.hpp"
#include "cl_Cell.hpp" //CON/src
#include "cl_Mesh_Enums.hpp"
#include "MTK_Tools.hpp"
#include "cl_MTK_Mesh_Core.hpp" //MTK/src


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
        class Mesh : public virtual mtk::Mesh
        {
            //! describing label
            std::string mLabel;

//            Lagrange_Mesh_Base * mMesh = nullptr;

            moris::Cell< moris::hmr::BSpline_Mesh_Base* > mDummyBSplineMeshes;

            // give access to the member data to stk based interpolation and integration meshes
            friend class Interpolation_Mesh_HMR;
            friend class Integration_Mesh_HMR;

            // Global to Local Entity Map
            // moris::Cell(0) - Node Global to Local
            // moris::Cell(1) - Edge Global to Local
            // moris::Cell(2) - Face Global to Local
            // moris::Cell(3) - Cell Global to Local
            // Keep in mind not all of these are created
            moris::Cell<std::unordered_map<moris_id,moris_index>> mEntityGlobaltoLocalMap;

//-------------------------------------------------------------------------------
        public:
//-------------------------------------------------------------------------------

            /**
             * mesh constructor, to be called from HMR
             */
            Mesh(       std::shared_ptr< Database >   aDatabase,
                  const uint                        & aLagrangeOrder,
                  const uint                        & aLagrangePattern );


            Mesh(       std::shared_ptr< Database >   aDatabase,
                  const uint                        & aOrder,
                  const uint                        & aLagrangePattern,
                  const uint                        & aBsplinePattern);

            Mesh(       std::shared_ptr< Database >   aDatabase,
                  const uint                        & aLagrangeMeshIndex );

//-------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Mesh();

//-------------------------------------------------------------------------------

            /**
             * return the type of this mesh
             */
            MeshType get_mesh_type() const
            {
                return MeshType::HMR;
            }

//-------------------------------------------------------------------------------

            /**
             * provides a moris::Matrix< DDUMat > containing the IDs this mesh has
             * to communicate with
             */
            Matrix< IdMat > get_communication_table() const ;

//-------------------------------------------------------------------------------

            /**
             * creates a new field pointer that is linked to this mesh
             */
            std::shared_ptr< Field > create_field( const std::string & aLabel,
                                                   const uint        & aBSplineOrder );

//-------------------------------------------------------------------------------

            /**
             * returns a pointer to the underlying lagrange mesh
             */
            Lagrange_Mesh_Base * get_lagrange_mesh()
            {
                return mMesh;
            }

//-------------------------------------------------------------------------------

            /**
             * return a pointer to the database
             */
            std::shared_ptr< Database > get_database()
            {
                return mDatabase;
            }

//-------------------------------------------------------------------------------
// Functions for MTK
//-------------------------------------------------------------------------------

            /**
             * returns the number of dimensions of this mesh
             */
            uint get_spatial_dim() const;

//-------------------------------------------------------------------------------

            uint get_num_entities( enum EntityRank aEntityRank) const;

//-------------------------------------------------------------------------------

            uint get_num_elemens_including_aura() const ;

//-------------------------------------------------------------------------------

            uint get_num_nodes() const;

//-------------------------------------------------------------------------------

            uint get_num_nodes_including_aura() const;

//-------------------------------------------------------------------------------

            uint get_num_edges() const;

//-------------------------------------------------------------------------------

            uint get_num_faces() const;

//-------------------------------------------------------------------------------

            uint get_num_elems() const;

//-------------------------------------------------------------------------------

            uint get_num_coeffs( const uint aBSplineMeshIndex ) const;

//-------------------------------------------------------------------------------

            Matrix< IndexMat > get_bspline_inds_of_node_loc_ind( const moris_index      aNodeIndex,
                                                                 const enum EntityRank  aBSplineRank );

//-------------------------------------------------------------------------------

            Matrix< IndexMat > get_entity_connected_to_entity_loc_inds(
                                       moris_index     aEntityIndex,
                                       enum EntityRank aInputEntityRank,
                                       enum EntityRank aOutputEntityRank) const;

//-------------------------------------------------------------------------------

            Matrix< IndexMat > get_entity_connected_to_entity_glob_ids(
                                       moris_index     aEntityId,
                                       enum EntityRank aInputEntityRank,
                                       enum EntityRank aOutputEntityRank) const;

//-------------------------------------------------------------------------------

            Matrix< IndexMat > get_nodes_connected_to_node_loc_inds( moris_index aNodeIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix< IndexMat > get_nodes_connected_to_edge_loc_inds( moris_index aEdgeIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix< IndexMat > get_nodes_connected_to_face_loc_inds( moris_index aFaceIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix< IndexMat > get_nodes_connected_to_element_loc_inds( moris_index aElementIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix < IndexMat > get_edges_connected_to_node_loc_inds( moris_index aNodeIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix< IndexMat > get_edges_connected_to_element_loc_inds( moris_index aElementIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix < IndexMat > get_faces_connected_to_node_loc_inds( moris_index aNodeIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix< IndexMat > get_faces_connected_to_element_loc_inds( moris_index aElementIndex ) const ;

//-------------------------------------------------------------------------------

            Matrix < IndexMat > get_elements_connected_to_node_loc_inds( moris_index aNodeIndex ) const;

//-------------------------------------------------------------------------------

            Matrix< IndexMat > get_elements_connected_to_face_loc_inds( moris_index aFaceIndex ) const ;

//-------------------------------------------------------------------------------

            void get_elements_in_support_of_basis(const uint           aMeshIndex,
                                                  const uint           aBasisIndex,
                                                  Matrix< IndexMat > & aElementIndices);

//-------------------------------------------------------------------------------

            void get_nodes_indices_in_bounding_box( const moris::Matrix< DDRMat >   & aPoint,
                                                    const moris::Matrix< DDRMat >   & aBoundingBoxSize,
                                                          moris::Matrix< IndexMat > & aNodeIndices );

//-------------------------------------------------------------------------------
            uint get_num_basis_functions(const uint aMeshIndex);
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
            Matrix< IndexMat > get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const;

//-------------------------------------------------------------------------------
//          Global ID Functions
//-------------------------------------------------------------------------------

            moris_id get_glb_entity_id_from_entity_loc_index(
                    moris_index     aEntityIndex,
                    enum EntityRank aEntityRank) const ;

            moris_id get_glb_element_id_from_element_loc_index( moris_index aEntityIndex ) const;

//-------------------------------------------------------------------------------
            moris_index get_loc_entity_ind_from_entity_glb_id(
                    moris_id        aEntityId,
                    enum EntityRank aEntityRank) const;

//-------------------------------------------------------------------------------

            moris_id get_max_entity_id( enum EntityRank aEntityRank ) const ;

//-------------------------------------------------------------------------------
//          Coordinate Field Functions
//-------------------------------------------------------------------------------

            Matrix< DDRMat > get_node_coordinate( moris_index aNodeIndex ) const;

//-------------------------------------------------------------------------------
//           Entity Ownership Functions
//-------------------------------------------------------------------------------

            moris_id get_entity_owner( moris_index     aEntityIndex,
                                       enum EntityRank aEntityRank ) const;

            //FIXME Needs parallel implementation
            void
            get_processors_whom_share_entity(moris_index       aEntityIndex,
                                             enum EntityRank   aEntityRank,
                                             Matrix< IdMat > & aProcsWhomShareEntity) const;


            enum EntityRank
            get_facet_rank() const;

//-------------------------------------------------------------------------------
//           Set Functions
//-------------------------------------------------------------------------------

            void get_sideset_elems_loc_inds_and_ords(
                    const  std::string     & aSetName,
                    Matrix< IndexMat >     & aElemIndices,
                    Matrix< IndexMat >     & aSidesetOrdinals ) const;

//-------------------------------------------------------------------------------

            moris::Cell<std::string> get_set_names(enum EntityRank aSetEntityRank) const;

//-------------------------------------------------------------------------------

            Matrix< IndexMat > get_set_entity_loc_inds( enum EntityRank aSetEntityRank,
                                                        std::string     aSetName ) const;

//-------------------------------------------------------------------------------
//           Pointer Functions for FEM
//-------------------------------------------------------------------------------

            mtk::Vertex & get_mtk_vertex( moris_index aVertexIndex )
            {
                return *mMesh->get_node_by_index( aVertexIndex );
            }

//-------------------------------------------------------------------------------

            mtk::Vertex const & get_mtk_vertex( moris_index aVertexIndex ) const
            {
                return *mMesh->get_node_by_index( aVertexIndex );
            }

//-------------------------------------------------------------------------------

            mtk::Vertex & get_mtk_vertex_including_aura( moris_index aVertexIndex )
            {
                return *mMesh->get_node_by_index_including_aura( aVertexIndex );
            }

//-------------------------------------------------------------------------------

            mtk::Vertex const & get_mtk_vertex_including_aura( moris_index aVertexIndex ) const
            {
                return *mMesh->get_node_by_index_including_aura( aVertexIndex );
            }

//-------------------------------------------------------------------------------

            mtk::Cell & get_mtk_cell( moris_index aElementIndex )
            {
//                return *mMesh->get_element( aElementIndex );
                return *mMesh->get_element_including_aura( aElementIndex );
            }

//-------------------------------------------------------------------------------

            mtk::Cell const & get_mtk_cell( moris_index aElementIndex ) const
            {
//                return *mMesh->get_element( aElementIndex );
                return *mMesh->get_element_including_aura( aElementIndex );
            }
//-------------------------------------------------------------------------------

            moris::mtk::Facet*
            get_facet(moris_index aFacetIndex)
            {
                return mMesh->get_facet(aFacetIndex);
            }

//-------------------------------------------------------------------------------


            moris::Cell<moris::mtk::Vertex const *>
            get_all_vertices() const;

//-------------------------------------------------------------------------------

            moris::Cell<moris::mtk::Vertex const *>
            get_all_vertices_including_aura() const;

//-------------------------------------------------------------------------------

            void get_adof_map( const uint                           aBSplineIndex,
                                     map< moris_id, moris_index > & aAdofMap ) const;

//-------------------------------------------------------------------------------

            moris_index get_field_ind( const std::string     & aFieldLabel,
                                       const enum EntityRank   aEntityRank  ) const;

//-------------------------------------------------------------------------------

            uint get_num_fields( const enum EntityRank aEntityRank ) const;

//-------------------------------------------------------------------------------

            real & get_value_of_scalar_field(
                    const      moris_index  aFieldIndex,
                    const enum EntityRank   aEntityRank,
                    const uint              aEntityIndex );

//-------------------------------------------------------------------------------

            const real & get_value_of_scalar_field(
                    const      moris_index  aFieldIndex,
                    const enum EntityRank   aEntityRank,
                    const uint              aEntityIndex ) const;

//-------------------------------------------------------------------------------

            Matrix<DDRMat> & get_field( const moris_index     aFieldIndex,
                                        const enum EntityRank aEntityRank );

//-------------------------------------------------------------------------------
            private:
//-------------------------------------------------------------------------------

            void get_element_indices_from_memory_indices(
                    const Matrix< DDLUMat>      & aMemoryIndices,
                          Matrix< IndexMat >    & aIndices ) const;

//-------------------------------------------------------------------------------
            /**
             * subroutine for get_elements_connected_to_element_loc_inds(
             */
            void collect_memory_indices_of_active_element_neighbors(
                    const moris_index  aElementIndex,
                    Matrix< DDLUMat> & aMemoryIndices,
                    Matrix< DDLUMat> & aThisCellFacetOrds,
                    luint            & aCounter ) const;

//-------------------------------------------------------------------------------

            /**
             * subroutine for get_elements_connected_to_node_loc_inds
             */
            void collect_memory_indices_of_active_elements_connected_to_node(
                    const moris_index  aNodeIndex,
                    Matrix< DDLUMat> & aMemoryIndices ) const;

//-------------------------------------------------------------------------------
public:
            const Matrix< DDRMat > & get_t_matrix_of_node_loc_ind(
                    const moris_index aNodeIndex,
                    const enum EntityRank  aBSplineRank )
            {
                return *mMesh->get_node_by_index( aNodeIndex )->get_interpolation( 0)
                                                              ->get_weights();
            }
private:

//-------------------------------------------------------------------------------

            uint get_level_of_entity_loc_ind(
                    const enum EntityRank aEntityRank,
                    const uint            aEntityIndex );


//-------------------------------------------------------------------------------

            uint get_max_level_of_entity( const enum EntityRank aEntityRank );

//-------------------------------------------------------------------------------

            Matrix< IndexMat > get_inds_of_active_elements_connected_to_basis( const Basis * aBasis ) const;

//-------------------------------------------------------------------------------

            void setup_glb_to_local_maps();

//-------------------------------------------------------------------------------

            void setup_entity_global_to_local_map(enum EntityRank aEntityRank);

//-------------------------------------------------------------------------------

            void setup_element_global_to_local_map();


//-------------------------------------------------------------------------------
        };

    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_INTERFACE_HPP_ */

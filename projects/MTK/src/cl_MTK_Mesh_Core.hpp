/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_Core.hpp
 *
 */

#ifndef MORIS_CL_MTK_MESH_CORE_HPP_
#define MORIS_CL_MTK_MESH_CORE_HPP_

#include "assert.hpp"
#include "cl_Matrix.hpp"

#include "cl_MTK_Enums.hpp"
#include "MTK_Tools.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_Map.hpp"
#include "cl_MTK_Vertex.hpp"    //MTK/src
#include "cl_MTK_Cell.hpp"      //MTK/src
#include "cl_MTK_Facet.hpp"
#include "cl_MTK_Set.hpp"

#include <unordered_map>

namespace moris
{
    namespace hmr
    {
        class Database;
    }

    namespace mtk
    {
        class Set;
        class Vertex;
        class Vertex_Interpolation;
        class Interpolation_Mesh;

        class Mesh : public std::enable_shared_from_this< Mesh >
        {
          protected:

            moris::Cell<moris_index> mMesh_GEN_map;

            // FIXME these members are here only to allow for throwing, should be removed later
            mtk::Vertex* mDummyVertex = nullptr;
            mtk::Cell*   mDummyCells  = nullptr;
            real         mDummyReal   = 0.0;

            Matrix< DDRMat > mDummyMatrix;
            Matrix< DDSMat > mDummyMatrix2;

            //! ref to hmr object for multigrid FIXME
            std::shared_ptr< hmr::Database > mDatabase;

            // ----------------------------------------------------------------------------

          public:

            // ----------------------------------------------------------------------------
            // Verbose flag
            bool mVerbose = false;

            // ----------------------------------------------------------------------------

            /**
             * Constructor
             */
            Mesh();

            // ----------------------------------------------------------------------------

            /**
             * Destructor
             */
            virtual ~Mesh();

            //##############################################
            // 1.) General mesh information access
            //##############################################

            virtual void get_Mesh_GEN_map(moris::Cell<moris_index> &aMesh_GEN_map);

            // ----------------------------------------------------------------------------

            /**
             * Returns the type enum for this mesh.
             *
             * @return Mesh type
             */
            virtual MeshType get_mesh_type() const = 0;

            // ----------------------------------------------------------------------------

            /**
             * Returns the spatial dimension of this mesh.
             *
             * @return Spatial dimension
             */
            virtual uint get_spatial_dim() const = 0;

            // ----------------------------------------------------------------------------

            /**
             * Gets the polynomial order of this mesh
             *
             * @return Polynomial degree
             */
            virtual uint get_order();

            // ----------------------------------------------------------------------------

            /**
             * Gets the polynomial order of an underlying discretization of this mesh
             *
             * @param aDiscretizationIndex Discretization index, if applicable
             * @return Polynomial degree
             */
            virtual uint get_discretization_order( uint aDiscretizationIndex = 0 );

            // ----------------------------------------------------------------------------

            virtual luint get_num_active_bg_elements_on_discretization_mesh_index_including_aura( moris_index const aDiscretizationMeshIndex );

            // ----------------------------------------------------------------------------

            virtual void get_active_bg_element_indices_on_discretization_mesh_index_including_aura(
                    moris_index const  aDiscretizationMeshIndex,
                    Matrix< DDLUMat >& aElementIDs );

            // -----------------------------------------------------------------------------

            /**
             * @brief Get the lagrange elements within one bspline element
             *
             * @param aBspElementIndex index of the b-spline element
             * @param aDiscretizationMeshIndex discretization mesh index
             * @param aCells output: list of Lagrange elements as mtk::cells that sit inside the B-spline element
             */
            virtual void get_elements_in_bspline_element(
                    moris_index const          aBspElementIndex,
                    moris_index const          aDiscretizationMeshIndex,
                    moris::Cell< mtk::Cell* >& aCells );

            // -----------------------------------------------------------------------------

            /**
             * @brief Get the lagrange elements in the bspline elements for the whole mesh
             *
             * @param aDiscretizationMeshIndex discretization mesh index
             * @param aCells list of lists of Lagrange elements (mtk::cells) inside each B-spline element
             * @param aCellIndices list of lists of Lagrange elements (indices) inside each B-spline element
             */
            virtual void
            get_lagrange_elements_in_bspline_elements(
                    moris_index const                          aDiscretizationMeshIndex,
                    moris::Cell< moris::Cell< mtk::Cell* > >&  aCells,
                    moris::Cell< moris::Cell< moris_index > >& aCellIndices,
                    moris::Cell< moris_index >&                aLagToBspCellIndices,
                    moris::Cell< uint >&                       aBspCellRefineLevels,
                    moris::Cell< mtk::Cell* >&                 aBspCells );

            // -----------------------------------------------------------------------------

            /**
             * @brief Get the bspline element ijk and level
             *
             * @param aDiscretizationMeshIndex discretization mesh index
             * @param aBsplineElementIndex 
             * @param aIJK
             * @param aLevel
             */

            virtual const luint*
            get_bspline_element_ijk_level(
                    moris_index      aDiscretizationMeshIndex,
                    const mtk::Cell* aBsplineElement,
                    uint&            aLevel );

            // -----------------------------------------------------------------------------

            /**
             * @brief Get the bspline element ijk level object
             *
             * @param aDiscretizationMeshIndex
             * @param aBsplineElementIndex
             * @param aIJK
             * @param aLevel
             */

            virtual void
            get_extended_t_matrix(
                    moris_index                                 aDiscretizationMeshIndex,
                    moris_index                                 aBSplineCellIndex,
                    moris::mtk::Cell&                           aLagrangeCell,
                    moris::Cell< moris::Cell< mtk::Vertex* > >& tBsplineBasis,
                    moris::Cell< Matrix< DDRMat > >&            tWeights );

            // -----------------------------------------------------------------------------

            /**
             * @brief Get the L2 projection matrix object
             *
             * @param aDiscretizationMeshIndex discretization mesh index
             * @param aRootBSplineCellIndex 
             * @param aExtendedBSplineCellIndex
             * @param tRootBsplineBasis
             * @param tExtendedBsplineBasis
             * @param tWeights coefficients to relate extended basis to root basis
             */

            virtual void
            get_L2_projection_matrix(
                    moris_index                                 aDiscretizationMeshIndex,
                    const mtk::Cell*                            aRootBSplineCell,
                    const mtk::Cell*                            aExtendedBSplineCell,
                    moris::Cell< moris::Cell< const mtk::Vertex* > >& tRootBsplineBasis,
                    moris::Cell< const mtk::Vertex* >&                tExtendedBsplineBasis,
                    moris::Cell< Matrix< DDRMat > >&            tWeights );

            // ----------------------------------------------------------------------------

            /**
             * Returns lagrange elements inside the same B-Spline elements as the input lagrange element
             *
             * @param aElementIndex Lagrange element Index
             * @param aDiscretizationMeshIndex dicretization mesh index
             *
             * @return Lagrange elements
             */
            virtual void get_elements_in_interpolation_cluster(
                    moris_index                aElementIndex,
                    moris_index                aDiscretizationMeshIndex,
                    moris::Cell< mtk::Cell* >& tCells );

            // ----------------------------------------------------------------------------

            /**
             * Returns lagrange elements on given side ordinal inside a given B-Spline element
             *
             * @param aElementIndex Lagrange element Index
             * @param aDiscretizationMeshIndex dicretization mesh index
             * @param aSideOrdinal Side Ordinal
             *
             * @return Lagrange elements
             */
            virtual void
            get_elements_in_bspline_element_and_side_ordinal(
                    moris_index const          aBsplineElementIndex,
                    moris_index const          aDiscretizationMeshIndex,
                    moris_index const          aSideOrdinal,
                    moris::Cell< mtk::Cell* >& aCells );

            // ----------------------------------------------------------------------------

            /**
             * Returns lagrange elements inside the same B-Spline elements as the input lagrange element and side ordinal
             *
             * @param aElementIndex Lagrange element Index
             * @param aDiscretizationMeshIndex dicretization mesh index
             * @param aSideOrdinal Side Ordinal
             *
             * @return Lagrange elements
             */
            virtual void get_elements_in_interpolation_cluster_and_side_ordinal(
                    moris_index const          aElementIndex,
                    moris_index const          aDiscretizationMeshIndex,
                    moris_index const          aSideOrdinal,
                    moris::Cell< mtk::Cell* >& aCells );

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
            virtual bool
            get_elements_connected_to_element_through_face_ord(
                    moris_index                 aBaseElementIndex,
                    moris_index                 aMySideOrdinal,
                    moris_index&                aMyRefineLevel,
                    moris::Cell< moris_index >& aNeighborElements,
                    moris::Cell< moris_index >& aNeighborElementSideOrdinals,
                    moris::Cell< moris_index >& aTransitionLocations,
                    moris::Cell< moris_index >& aNeighborRefinementLevels ) const
            {
                MORIS_ERROR( false,
                        "mtk::Mesh::get_elements_connected_to_element_through_face_ord() - "
                        "Entered virtual function in Mesh base class; function is not implemented." );

                return false;
            }

            //------------------------------------------------------------------------------

            // FIXME This should be default, individual calls should be virtual
            /**
             * Gets the number of entities for a specified entity rank.
             *
             * @param aEntityRank Entity rank (node, edge, etc.)
             * @param aIndex Entity index
             * @return Number of entities of this rank
             */
            virtual uint get_num_entities(
                    enum EntityRank   aEntityRank,
                    const moris_index aIndex = 0 ) const = 0;

            // ----------------------------------------------------------------------------
            
            /**
             * Get the number of sets on this mesh.
             *
             * @return Number of sets
             */
            virtual moris::uint get_num_sets() const;

            // ----------------------------------------------------------------------------
            
            virtual moris::mtk::Set* get_set_by_index( moris::uint aSetIndex ) const;

            // ----------------------------------------------------------------------------
            
            virtual moris::mtk::Set* get_set_by_name( std::string aSetLabel ) const;

            // end of pure virtual functions in section 1
            // all functions below this line need to be able to have a default implementation

            // ----------------------------------------------------------------------------
            
            // FIXME pure virtual
            /**
             * Gets the number of nodes on this mesh.
             *
             * @return Number of nodes
             */
            virtual uint get_num_nodes() const;

            // ----------------------------------------------------------------------------
            
            // FIXME pure virtual
            /**
             * Gets the number of edges on this mesh.
             *
             * @return Number of edges
             */
            virtual uint get_num_edges() const;

            // ----------------------------------------------------------------------------
            
            // FIXME pure virtual
            /**
             * Gets the number of faces on this mesh.
             *
             * @return Number of faces
             */
            virtual uint get_num_faces() const;

            // ----------------------------------------------------------------------------
            
            // FIXME pure virtual
            /**
             * Gets the number of elements on this mesh.
             *
             * @return Number of elements
             */
            virtual uint get_num_elems() const;

            //##############################################
            // 2.) Access Mesh Data by index Functions
            //##############################################
            //##############################################
            // 2.a.) Access standard mesh data
            //##############################################

            // FIXME this should be default, not pure virtual. Individual calls should be pure virtual.
            /**
             * Get all entity indices of a given output entity rank which are connected to an entity of a given index
             * and input rank.
             *
             * @param aEntityIndex Input entity index
             * @param aInputEntityRank Input entity rank
             * @param aOutputEntityRank Output entity rank
             * @param aDiscretizationIndex discretization mesh index
             * @return Output entity indices
             */
            virtual Matrix< IndexMat >
            get_entity_connected_to_entity_loc_inds(
                    moris_index       aEntityIndex,
                    enum EntityRank   aInputEntityRank,
                    enum EntityRank   aOutputEntityRank,
                    const moris_index aDiscretizationIndex = 0 ) const = 0;

            // ----------------------------------------------------------------------------
            
            /**
             * Since the connectivity between entities of the same rank are considered
             * invalid by STK standards, we need a separate function for element to element
             * specifically.
             *
             * @param[in]  aElementId - element id
             * @return Matrix< IndexMat > matrix containing connected elements and corresponding side ordinals
             *         Format:  1st row are indices of other connected elements/cells
             *                  2nd row are the side ordinal indices wrt. treated Bg cell
             *                  3rd row are the side ordinal indices wrt. respective neighbor Bg cell
             *                  4th row ? // TODO: need to understand this ~Nils
             */
            virtual Matrix< IndexMat >
            get_elements_connected_to_element_and_face_ord_loc_inds( moris_index aElementIndex ) const;

            // ----------------------------------------------------------------------------
            
            /**
             * @brief Since the connectivity between entities of the same rank are considered
             *        invalid by STK standards, we need a separate function for element to element
             *        specifically
             *
             * @param aElementIndex - index of Bg-cell for which connectivity is to be retrieved
             * @return Matrix< IndexMat > matrix containing connected elements and corresponding facets
             *         Format:  1st row are indices of other connected elements/cells
             *                  2nd row are the indices of the Bg-cell's facets through which they are connected
             */
            virtual Matrix< IndexMat >
            get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const
            {
                MORIS_ERROR( 0, "No default implementation" );
                return Matrix< IndexMat >( 0, 0 );
            }

            // ----------------------------------------------------------------------------
            
            // FIXME Remove access to vertex
            /**
             * Deprecated
             */
            virtual moris::Cell< moris::mtk::Vertex const * >
            get_all_vertices() const;

            // end of pure virtual functions in section 2.1

            // ----------------------------------------------------------------------------
            
            // FIXME pure virtual
            /**
             * Get all element indices connected to a node.
             *
             * @param aNodeIndex Node index
             * @return Element indices
             */
            virtual Matrix< IndexMat >
            get_elements_connected_to_node_loc_inds( moris_index aNodeIndex ) const;

            // ----------------------------------------------------------------------------
            
            // FIXME pure virtual
            /**
             * Get all face indices connected to a node.
             *
             * @param aNodeIndex Node index
             * @return Face indices
             */
            virtual Matrix< IndexMat >
            get_faces_connected_to_node_loc_inds( moris_index aNodeIndex ) const;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get all edge indices connected to a node.
             *
             * @param aNodeIndex Node index
             * @return Edge indices
             */
            virtual Matrix< IndexMat >
            get_edges_connected_to_node_loc_inds( moris_index aNodeIndex ) const;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get all element indices connected to an edge.
             *
             * @param aEdgeIndex
             * @return Element indices
             */
            virtual Matrix< IndexMat >
            get_elements_connected_to_edge_loc_inds( moris_index aEdgeIndex ) const;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get all face indices connected to an edge.
             *
             * @param aEdgeIndex Edge index
             * @return Face indices
             */
            virtual Matrix< IndexMat >
            get_faces_connected_to_edge_loc_inds( moris_index aEdgeIndex ) const;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get all element indices connected to a face
             *
             * @param aFaceIndex Face index
             * @return Element indices
             */
            virtual Matrix< IndexMat >
            get_elements_connected_to_face_loc_inds( moris_index aFaceIndex ) const;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get all face indices connected to an element
             *
             * @param aElementIndex Element index
             * @return Face indices
             */
            virtual Matrix< IndexMat >
            get_faces_connected_to_element_loc_inds( moris_index aElementIndex ) const;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get all edge indices connected to an element
             *
             * @param aElementIndex Element index
             * @return Edge indices
             */
            virtual Matrix< IndexMat >
            get_edges_connected_to_element_loc_inds( moris_index aElementIndex ) const;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get all node indices connected to an element
             *
             * @param aElementIndex Element index
             * @return Node indices
             */
            virtual Matrix< IndexMat >
            get_nodes_connected_to_element_loc_inds( moris_index aElementIndex ) const;

            //##############################################
            // 2.a.) Access mesh data from ids
            //##############################################

            //##############################################
            // global id functions
            //##############################################

            /**
             * Get a global entity ID from an entity rank and local index.
             *
             * @param aEntityIndex Local entity index
             * @param aEntityRank Entity rank
             * @param aDiscretizationIndex discretization mesh Index
             * @return Global entity ID
             */
            virtual moris_id
            get_glb_entity_id_from_entity_loc_index(
                    moris_index aEntityIndex,
                    EntityRank  aEntityRank,
                    moris_index aDiscretizationIndex = 0 ) const = 0;

            // ----------------------------------------------------------------------------

            virtual std::unordered_map< moris_id, moris_index >
            get_vertex_glb_id_to_loc_vertex_ind_map() const
            {
                MORIS_ERROR( false, "get_vertex_glb_id_to_loc_vertex_ind_map(), Not implemented for this mesh type" );
                std::unordered_map< moris_id, moris_index > tMap;

                return tMap;
            };

            // ----------------------------------------------------------------------------

            // FIXME pure virtual or default implementation without error
            /**
             * Get a local entity index from a global ID.
             *
             * @param aEntityId Global entity ID
             * @param aEntityRank Entity rank
             * @param aDiscretizationIndex discretization mesh index
             * @return Local entity index
             */
            virtual moris_index
            get_loc_entity_ind_from_entity_glb_id(
                    moris_id    aEntityId,
                    EntityRank  aEntityRank,
                    moris_index aDiscretizationIndex = 0 ) const = 0;

            // ----------------------------------------------------------------------------

            /**
             * Get a facet ordinal from a face index and cell index.
             *
             * @param aFaceIndex Face index
             * @param aCellIndex Cell index
             * @return Facet ordinal
             */
            virtual moris::moris_index
            get_facet_ordinal_from_cell_and_facet_loc_inds( moris::moris_index aFaceIndex,
                    moris::moris_index                                         aCellIndex ) const;

            // ----------------------------------------------------------------------------

            // FIXME split into individual calls
            /**
             * Generates unique entity IDs for a given entity rank.
             *
             * @param aNumEntities Number of entities to generate IDs for
             * @param aEntityRank Entity rank
             * @return Unique IDs
             */
            virtual Matrix< IdMat >
            generate_unique_entity_ids(
                    uint            aNumEntities,
                    enum EntityRank aEntityRank ) const;

            // ----------------------------------------------------------------------------

            /**
             * Generate unique node IDs.
             *
             * @param aNumNodes Number of nodes to generate IDs for
             * @return Unique IDs
             */
            virtual Matrix< IdMat >
            generate_unique_node_ids( uint aNumNodes );

            // ----------------------------------------------------------------------------

            // FIXME should have default instead of pure virtual, individual functions should be pure virtual instead
            /**
             * Generic get global id of entities connected to
             * entity using an entities global id
             */
            virtual Matrix< IdMat >
            get_entity_connected_to_entity_glob_ids(
                    moris_id          aEntityId,
                    enum EntityRank   aInputEntityRank,
                    enum EntityRank   aOutputEntityRank,
                    const moris_index aDiscretizationMeshIndex = 0 ) const;

            // ----------------------------------------------------------------------------

            /**
             * Since the connectivity between entities of the same rank are considered
             * invalid by STK standards, we need a seperate function for element to element
             * specifically
             *
             * @param[in]  aElementId - element id
             * @param[out] Element to element connectivity and face ordinal shared
             */
            virtual Matrix< IdMat >
            get_element_connected_to_element_glob_ids( moris_id aElementId ) const;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get element IDs connected to a node.
             *
             * @param aNodeId Node ID
             * @return Element IDs
             */
            virtual Matrix< IdMat >
            get_elements_connected_to_node_glob_ids( moris_id aNodeId );

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get face IDs connected to a node.
             *
             * @param aNodeId Node ID
             * @return Face IDs
             */
            virtual Matrix< IdMat >
            get_faces_connected_to_node_glob_ids( moris_id aNodeId );

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get edge IDs connected to a node.
             *
             * @param aNodeId Node ID
             * @return Edge IDs
             */
            virtual Matrix< IdMat >
            get_edges_connected_to_node_glob_ids( moris_id aNodeId );

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get element IDs connected to an edge.
             *
             * @param aEdgeId Edge ID
             * @return Element IDs
             */
            virtual Matrix< IdMat >
            get_elements_connected_to_edge_glob_ids( moris_id aEdgeId );

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get face IDs connected to an edge.
             *
             * @param aEdgeId Edge ID
             * @return Face IDs
             */
            virtual Matrix< IdMat >
            get_faces_connected_to_edge_glob_ids( moris_id aEdgeId );

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get element IDs connected to a face.
             *
             * @param aFaceId Face ID
             * @return Element IDs
             */
            virtual Matrix< IdMat >
            get_elements_connected_to_face_glob_ids( moris_id aFaceId );

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get face IDs connected to an element.
             *
             * @param aElementId Element ID
             * @return Face IDs
             */
            virtual Matrix< IdMat >
            get_faces_connected_to_element_glob_ids( moris_id aElementId );

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get edge IDs connected to an element.
             *
             * @param aElementId Element ID
             * @return Edge IDs
             */
            virtual Matrix< IdMat >
            get_edges_connected_to_element_glob_ids( moris_id aElementId );

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Get node IDs connected to an element.
             *
             * @param aElementId Element ID
             * @return Node IDs
             */
            virtual Matrix< IdMat >
            get_nodes_connected_to_element_glob_ids( moris_id aElementId );

            // ----------------------------------------------------------------------------

            // FIXME default implemenation with no error
            /**
             * Get elements interpolated into by a basis function. For a Lagrange mesh,
             * the elements in support of basis is equivalent to the elements connected
             * to a node. Therefore, a call to get_elements
             */
            virtual void
            get_elements_in_support_of_basis(
                    const uint          aMeshIndex,
                    const uint          aBasisIndex,
                    Matrix< IndexMat >& aElementIndices );

            // ----------------------------------------------------------------------------

            // TODO determine if we can remove this
            /**
             * Get the node indices in a bounding box.
             *
             * @param aPoint Point to evaluate
             * @param aBoundingBoxSize Bounding box size
             * @param aNodeIndices Returned node indices
             */
            virtual void
            get_nodes_indices_in_bounding_box(
                    const moris::Matrix< DDRMat >& aPoint,
                    const moris::Matrix< DDRMat >& aBoundingBoxSize,
                    moris::Matrix< IndexMat >&     aNodeIndices );

            // ----------------------------------------------------------------------------

            /**
             * Get the spatial coordinates of a node.
             *
             * @param aNodeIndex Node index
             * @return Node coordinates
             */
            virtual Matrix< DDRMat >
            get_node_coordinate( moris_index aNodeIndex ) const = 0;

            //##############################################
            // Field Access
            //##############################################

            // TODO: introduce a concept of field indices to prevent accessing via a name which
            // TODO: involves a string comparison
            //  FIXME remove from base class
            virtual Matrix< DDRMat >
            get_entity_field_value_real_scalar(
                    const Matrix< IndexMat >& aEntityIndices,
                    const std::string&        aFieldName,
                    enum EntityRank           aFieldEntityRank ) const;

            // ----------------------------------------------------------------------------

            // FIXME remove from base class
            virtual void
            add_mesh_field_real_scalar_data_loc_inds(
                    const std::string&      aFieldName,
                    const enum EntityRank&  aFieldEntityRank,
                    const Matrix< DDRMat >& aFieldData );

            //##############################################
            // Facet Access
            //##############################################

            // FIXME remove access to facet
            virtual moris::mtk::Facet*
                    get_facet( moris_index );

            //##############################################
            // Cell and Vertex Pointer Functions
            //##############################################

            // FIXME remove access to cell
            virtual mtk::Cell&
            get_mtk_cell( moris_index aElementIndex );

            // ----------------------------------------------------------------------------

            // FIXME remove access to cell
            virtual mtk::Cell const &
            get_mtk_cell( moris_index aElementIndex ) const;

            // ----------------------------------------------------------------------------

            // FIXME remove access to vertex
            virtual mtk::Vertex&
            get_mtk_vertex( moris_index aVertexIndex );

            // ----------------------------------------------------------------------------

            // FIXME remove access to vertex
            virtual mtk::Vertex const &
            get_mtk_vertex( moris_index aVertexIndex ) const;

            // ----------------------------------------------------------------------------

            /**
             * Gets the index of a node's base node. Can be different than the node index itself if the node was created
             * based on a different node.
             *
             * @param aNodeIndex Node index
             */
            uint get_base_node_index( uint aNodeIndex );

            //##############################################
            // For FEM
            //##############################################

            // FIXME remove access to cell, add set functions for cell instead
            virtual mtk::Cell&
            get_writable_mtk_cell( moris_index aElementIndex );

            // ----------------------------------------------------------------------------

            // FIXME split into only the needed calls (node from what I can tell), make pure virtual
            /**
             * Gets the max entity ID for a given entity rank.
             *
             * @param aEntityRank Entity rank
             * @param aDiscretizationIndex discretization mesh index
             * @return Max entity ID
             */
            virtual moris_id
            get_max_entity_id( enum EntityRank aEntityRank,
                    const moris_index          aDiscretizationIndex = 0 ) const;

            //##############################################
            // Entity Ownership Functions
            //##############################################

            /**
             * Gets the owner of a node.
             *
             * @param aNodeIndex Node index
             * @return Node owner
             */
            virtual uint get_node_owner( moris_index aNodeIndex ) const = 0;

            // ----------------------------------------------------------------------------

            /**
             * Gets the owner of an element.
             *
             * @param aElementIndex Element index
             * @return Element owner
             */
            virtual uint get_element_owner( moris_index aElementIndex ) const = 0;

            // ----------------------------------------------------------------------------

            /**
             * Gets the owner of a given entity.
             *
             * @param aEntityIndex Entity index
             * @param aEntityRank Entity rank
             * @param aDiscretizationMeshIndex discretization mesh index
             * @return Entity owner
             */
            virtual uint get_entity_owner(
                    moris_index       aEntityIndex,
                    enum EntityRank   aEntityRank,
                    const moris_index aDiscretizationMeshIndex = 0 ) const;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual or default implementation
            /**
             * Processors whom share a given entity
             * @param[in]  - Entity Index
             * @param[in]  - Entity Rank
             * @param[out] - Processors whom share an entity vector
             */
            virtual void
            get_processors_whom_share_entity(
                    moris_index      aEntityIndex,
                    enum EntityRank  aEntityRank,
                    Matrix< IdMat >& aProcsWhomShareEntity ) const;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual or default implementation
            /**
             * Gets the number of entities shared with a processor.
             *
             * @param aProcessorRank Processor rank
             * @param aEntityRank Entity rank
             * @param aSendFlag Send flag
             * @return Number of entities
             */
            virtual uint
            get_num_of_entities_shared_with_processor(
                    moris_id        aProcessorRank,
                    enum EntityRank aEntityRank,
                    bool            aSendFlag ) const;

            //##############################################
            // Communication Tables
            //##############################################

            // FIXME: THIS FUNCTION DESCRIPTION NEEDS TO BE IMPROVED
            // FIXME: Also, a unit test (not clear what needs to be provided)
            /**
             * provides a moris::Mat<uint> containing the IDs this mesh has
             * to communicate with
             *
             */
            virtual Matrix< IdMat >
            get_communication_table() const = 0;

            // ----------------------------------------------------------------------------

            virtual Matrix< IdMat >
            get_communication_proc_ranks() const;

            // ----------------------------------------------------------------------------

            virtual moris::Cell< Matrix< IdMat > >
            get_communication_vertex_pairing() const;

            //##############################################
            //  Output Mesh To a File
            //##############################################

            // FIXME default implementation with exodus writer, or remove
            /**
             * Create an exodus mesh database with the specified
             * filename.
             *
             * @param[in] filename The full pathname to the file which will be
             *   created and the mesh data written to. If the file already
             *   exists, it will be overwritten.
             *
             *   Description from create_output_mesh() in StkMeshIoBroker.hpp
             */
            virtual void
            create_output_mesh(
                    std::string& aFileName,
                    bool         aAddElemCmap = false );

            //##############################################
            //  Field Functions TODO sort out which of these need to be in the base
            //##############################################

            /**
             * Return the number of fields on this mesh.
             *
             * @param aEntityRank Entity rank
             * @param aDiscretizationMeshIndex discretization mesh index
             * @return Number of fields
             */
            virtual uint
            get_num_fields(
                    const enum EntityRank aEntityRank,
                    const moris_index     aDiscretizationMeshIndex = 0 ) const;

            // ----------------------------------------------------------------------------

            /**
             * return the index of the field of this label
             * return gNoIndex if not found
             */
            virtual moris_index get_field_ind(
                    const std::string&    aFieldLabel,
                    const enum EntityRank aEntityRank ) const;

            // ----------------------------------------------------------------------------

            /**
             * add a scalar field to the database
             *
             * fixme: how to make sure that field does not exist ?
             */
            virtual moris_index
            create_scalar_field(
                    const std::string&    aFieldLabel,
                    const enum EntityRank aEntityRank );

            // ----------------------------------------------------------------------------

            /**
             * add a vector field to the database
             */
            virtual moris_index
            create_vector_field(
                    const std::string&    aFieldLabel,
                    const enum EntityRank aEntityRank,
                    const uint            aDimension );

            // ----------------------------------------------------------------------------

            /**
             * get value of entity
             */
            virtual real&
            get_value_of_scalar_field(
                    const moris_index     aFieldIndex,
                    const enum EntityRank aEntityRank,
                    const uint            aEntityIndex,
                    const moris_index     aDiscretizationMeshIndex = 0 );

            // ----------------------------------------------------------------------------

            /**
             * get value of entity ( const version )
             */
            virtual const real&
            get_value_of_scalar_field(
                    const moris_index     aFieldIndex,
                    const enum EntityRank aEntityRank,
                    const uint            aEntityIndex,
                    const moris_index     aDiscretizationMeshIndex = 0 ) const;

            // ----------------------------------------------------------------------------

            /**
             * fixme: need opinion: sould we always return a DDRMat?
             *        should this be a row or column vector?
             */
            virtual Matrix< DDRMat >&
            get_value_of_vector_field(
                    const moris_index     aFieldIndex,
                    const enum EntityRank aEntityRank,
                    const uint            aEntityIndex );

            // ----------------------------------------------------------------------------

            /**
             * return the entry of a vector field ( const version )
             */
            virtual const Matrix< DDRMat >&
            get_value_of_vector_field(
                    const moris_index     aFieldIndex,
                    const enum EntityRank aEntityRank,
                    const uint            aEntityIndex ) const;

            // ----------------------------------------------------------------------------

            // FIXME not required by mapper anymore, can/should remove from here and mapper
            /**
             * returns a moris::Matrix with the field
             * This function is specific to HMR, and called by the mapper
             * if HMR is used.
             */
            virtual Matrix< DDRMat >&
            get_field( const moris_index  aFieldIndex,
                    const enum EntityRank aEntityRank,
                    const moris_index     aDiscretizationMeshIndex = 0 );

            //##############################################
            //  Multigrid
            //##############################################

            /**
             * Gets the max level of an entity.
             *
             * @param aEntityRank Entity rank
             * @param aDiscretizationMeshIndex discretization mesh index
             * @return Max level
             */
            virtual uint
            get_max_level_of_entity(
                    const enum EntityRank aEntityRank,
                    const moris_index     aDiscretizationMeshIndex = 0 );

            // ----------------------------------------------------------------------------

            /**
             * Gets the level of an entity.
             *
             * @param aEntityRank Entity rank
             * @param aEntityIndex Entity index
             * @param aDiscretizationMeshIndex discretization mesh index
             * @return
             */
            virtual uint
            get_level_of_entity_loc_ind(
                    const enum EntityRank aEntityRank,
                    const uint            aEntityIndex,
                    const moris_index     aDiscretizationMeshIndex = 0 );

            // ----------------------------------------------------------------------------

            // FIXME breaks inheritance
            std::shared_ptr< hmr::Database > get_HMR_database();

            // ----------------------------------------------------------------------------

            /**
             * Gets the shared IDs of the coefficients which define a discretization.
             *
             * Note. This function is pretty much useless and should not be implemented in mtk.
             * I strongly recommend to not use this function.
             * I will not take liability for possible foreseeable problems arising from the use of this function
             *
             * @param aNodeIndices Node indices to check for coefficients
             * @param aDiscretizationIndex Index of the specific discretization
             * @return Shared coefficient IDs
             */
            virtual Matrix< DDUMat > get_shared_discretization_coefficient_indices(
                    const Matrix< DDUMat >& aNodeIndices,
                    uint                    aDiscretizationIndex );

            // ----------------------------------------------------------------------------

            /**
             * Gets the owned IDs of the coefficients which define a discretization.
             *
             * Note. This function is pretty much useless and should not be implemented in mtk.
             * I strongly recommend to not use this function.
             * I will not take liability for possible foreseeable problems arising from the use of this function
             *
             * @param aNodeIndices Node indices to check for coefficients
             * @param aDiscretizationIndex Index of the specific discretization
             * @return Owned coefficient IDs
             */
            virtual Matrix< DDUMat > get_owned_discretization_coefficient_indices(
                    const Matrix< DDUMat >& aNodeIndices,
                    uint                    aDiscretizationIndex );

            // ----------------------------------------------------------------------------

            /**
             * Gets the number of discretization coefficients on a discretization mesh.
             *
             * @param aDiscretizationIndex discretization mesh index
             * @return Number of discretization coefficients
             */
            virtual uint
            get_max_num_coeffs_on_proc( uint aDiscretizationIndex ) const;

            // ----------------------------------------------------------------------------

            /**
             * Get the T-matrix of a node.
             *
             * @param aNodeIndex Node index
             * @param aDiscretizationIndex discretization mesh index
             * @return T-matrix
             */
            virtual const Matrix< DDRMat >& get_t_matrix_of_node_loc_ind(
                    uint aNodeIndex,
                    uint aDiscretizationIndex );

            // ----------------------------------------------------------------------------

            /**
             * Get the indices of the discretization coefficients of a node.
             *
             * @param aNodeIndex Node index
             * @param aDiscretizationIndex discretization mesh index
             * @return discretization coefficient indices
             */
            virtual Matrix< IndexMat > get_coefficient_indices_of_node(
                    uint aNodeIndex,
                    uint aDiscretizationIndex );

            // ----------------------------------------------------------------------------

            virtual Matrix< IdMat > get_coefficient_owners_of_node(
                    uint aNodeIndex,
                    uint aBSplineMeshIndex );

            // ----------------------------------------------------------------------------

            virtual Matrix< IdMat > get_coefficient_ijkl_IDs_of_node(
                    uint aNodeIndex,
                    uint aBSplineMeshIndex );

            // ----------------------------------------------------------------------------

            /**
             * Get the IDs of the discretization coefficients of a node.
             *
             * @param aNodeIndex Node index
             * @param aDiscretizationIndex discretization mesh index
             * @return discretization coefficient indices
             */
            virtual Matrix< IdMat > get_coefficient_IDs_of_node(
                    uint aNodeIndex,
                    uint aDiscretizationIndex );

            // ----------------------------------------------------------------------------

            /**
             * Gets the number of basis functions. For Lagrange meshes, the number of basis functions and the number of
             * nodes are equivalent. Therefore, a default implementation using get_num_nodes() is used here.
             *
             * @param aMeshIndex Mesh index
             * @return Number of basis functions
             */
            virtual uint
            get_num_basis_functions( const uint aMeshIndex = 0 );

            // ----------------------------------------------------------------------------

            // FIXME: Rename or use get loc entity id from global entity id
            // FIXME pure virtual
            /**
             * Get the adof map.
             *
             * @param aBSplineIndex discretization index
             * @param aAdofMap Adof map
             */
            void virtual get_adof_map(
                    const uint                    aBSplineIndex,
                    map< moris_id, moris_index >& aAdofMap ) const;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Gets the interpolation order of a field on this mesh.
             *
             * @param aFieldIndex Field index
             * @param aEntityRank Entity rank
             * @return Interpolation order
             */
            virtual uint
            get_order_of_field(
                    const moris_index     aFieldIndex,
                    const enum EntityRank aEntityRank );

            // ----------------------------------------------------------------------------

            // FIXME pure virtual or remove
            /**
             * Gets all set names for a given entity rank.
             *
             * @param aSetEntityRank Entity rank of the set
             * @return All set names
             */
            virtual moris::Cell< std::string >
            get_set_names( enum EntityRank aSetEntityRank ) const;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual or remove
            /**
             * Gets the indices of a set entity.
             *
             * @param aSetEntityRank Entity rank of the set
             * @param aSetName Set name
             * @return Entity indices in the set
             */
            virtual Matrix< IndexMat >
            get_set_entity_loc_inds(
                    enum EntityRank aSetEntityRank,
                    std::string     aSetName ) const;

            // ----------------------------------------------------------------------------

            /**
             * Gets element indices in a block set.
             *
             * @param aSetIndex Block set index
             * @return Element indices in the set
             */
            virtual Matrix< IndexMat > get_element_indices_in_block_set( uint aSetIndex ) = 0;

            // ----------------------------------------------------------------------------

            /**
             * Gets the element IDs in a block set, in order by index. Default implementation is to get the indices
             * and then transform them, but this can be overridden.
             *
             * @param aSetIndex Block set index
             * @return Element IDs in the set
             */
            virtual Matrix< IdMat > get_element_ids_in_block_set( uint aSetIndex );

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Gets the cell topology of a block set.
             *
             * @param aSetName Set name
             * @return Cell topology type
             */
            virtual enum CellTopology
            get_blockset_topology( const std::string& aSetName ) = 0;

            // ----------------------------------------------------------------------------

            /**
             * Gets the cell shape of a block set.
             *
             * @param aSetName Set name
             * @return Cell shape type
             */
            virtual enum CellShape
            get_IG_blockset_shape( const std::string& aSetName ) = 0;

            // ----------------------------------------------------------------------------

            /**
             * Gets the IP cell shape of a block set.
             *
             * @param aSetName Set name
             * @return Cell shape type
             */
            virtual enum CellShape
            get_IP_blockset_shape( const std::string& aSetName ) = 0;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Gets the cell topology of a sideset.
             *
             * @param aSetName Set name
             * @return Cell topology type
             */
            virtual enum CellTopology
            get_sideset_topology( const std::string& aSetName );

            // ----------------------------------------------------------------------------

            // FIXME remove access to cell
            /**
             * Deprecated
             */
            virtual moris::Cell< mtk::Cell const * > get_set_cells( std::string aSetLabel ) const;

            // ----------------------------------------------------------------------------

            // FIXME remove access to cell
            /**
             * Deprecated
             */
            virtual moris::Cell< mtk::Cell const * >
            get_block_set_cells( std::string aSetName ) const;

            // ----------------------------------------------------------------------------

            // FIXME pure virtual
            /**
             * Gets the element indices and their ordinals which lie on a given sideset.
             *
             * @param aSetName Sideset name
             * @param aElemIndices Element indices
             * @param aSidesetOrdinals Sideset ordinals
             */
            virtual void get_sideset_elems_loc_inds_and_ords(
                    const std::string&  aSetName,
                    Matrix< IndexMat >& aElemIndices,
                    Matrix< IndexMat >& aSidesetOrdinals ) const;

            // ----------------------------------------------------------------------------

            // FIXME remove access to cell
            /**
             * Deprecated
             */
            virtual void get_sideset_cells_and_ords(
                    const std::string&                aSetName,
                    moris::Cell< mtk::Cell const * >& aCells,
                    Matrix< IndexMat >&               aSidesetOrdinals ) const;

            // ----------------------------------------------------------------------------

            /**
             * Get the number of faces on the given sideset.
             *
             * @param aSideSetIndex Sideset indices
             * @return Number of faces
             */
            virtual uint get_sidesets_num_faces( moris::Cell< moris_index > aSideSetIndex ) const;

            // ----------------------------------------------------------------------------

            // FIXME remove access to vertex
            virtual moris::Cell< moris::mtk::Vertex const * >
            get_vertices_in_vertex_set_no_aura( std::string aSetName ) const;

            // ----------------------------------------------------------------------------

            /**
             * Gets the entity rank of a "facet" on this mesh.
             *
             * @return Facet rank
             */
            virtual enum EntityRank
            get_facet_rank() const;

            // ----------------------------------------------------------------------------

            // FIXME remove access to cell
            void
            get_mtk_cells(
                    Matrix< IndexMat >                       aCellInds,
                    moris::Cell< moris::mtk::Cell const * >& aCells );

            //##############################################
            // Multigrid acessor functions FIXME default implementation for non-multigrid + documentation
            //##############################################

            virtual uint get_num_interpolations();

            // ----------------------------------------------------------------------------

            virtual uint get_max_level( const moris_index aInterpolationIndex );

            // ----------------------------------------------------------------------------

            virtual uint get_num_basis( const moris_index aInterpolationIndex );

            // ----------------------------------------------------------------------------

            virtual uint get_basis_level( const moris_index aInterpolationIndex,
                    const moris_index                       aBasisIndex );

            // ----------------------------------------------------------------------------

            virtual uint get_num_coarse_basis_of_basis(
                    const moris_index aInterpolationIndex,
                    const moris_index aBasisIndex );

            // ----------------------------------------------------------------------------

            virtual uint get_coarse_basis_index_of_basis(
                    const moris_index aInterpolationIndex,
                    const moris_index aBasisIndex,
                    const moris_index aCoarseParentIndex );

            // ----------------------------------------------------------------------------

            virtual moris::Matrix< DDSMat > get_fine_basis_inds_of_basis(
                    const moris_index aInterpolationIndex,
                    const moris_index aBasisIndex );

            // ----------------------------------------------------------------------------

            virtual moris::Matrix< DDRMat > get_fine_basis_weights_of_basis(
                    const moris_index aInterpolationIndex,
                    const moris_index aBasisIndex );

#ifdef MORIS_HAVE_DEBUG
            virtual Matrix< DDRMat > get_basis_coords(
                    const moris_index aInterpolationIndex,
                    const moris_index aBasisIndex );

            virtual sint get_basis_status(
                    const moris_index aInterpolationIndex,
                    const moris_index aBasisIndex );
#endif

            /////////////////////////
            // Accessor functions for the data base entities
            /////////////////////////

            //--------------------------------------------------------------------------------------------------------------

            virtual moris_id
            get_entity_id( enum EntityRank aEntityRank,
                    moris_index            aEntityIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

            virtual moris_id
            get_entity_owner( enum EntityRank aEntityRank,
                    moris_index               aEntityIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

            virtual uint
            get_local_mesh_index( const uint aBsplineMeshIndex );

            //--------------------------------------------------------------------------------------------------------------

            virtual moris_id*
            get_basis_ids( moris_index aVertexIndex,
                    moris_index        aOrder );

            //--------------------------------------------------------------------------------------------------------------

            virtual moris_index*
            get_basis_indicies( moris_index aVertexIndex,
                    moris_index             aOrder );

            //--------------------------------------------------------------------------------------------------------------

            virtual real*
            get_basis_weights( moris_index aVertexIndex,
                    moris_index            aOrder );

            //--------------------------------------------------------------------------------------------------------------

            virtual moris_id*
            get_basis_owners( moris_index aVertexIndex,
                    moris_index           aOrder );

            //--------------------------------------------------------------------------------------------------------------

            virtual moris_index
            get_basis_length( moris_index aVertexIndex,
                    moris_index           aOrder );

            //--------------------------------------------------------------------------------------------------------------

            virtual moris::real*
            get_vertex_coords_ptr( moris_index aVertexIndex );

            //--------------------------------------------------------------------------------------------------------------

            virtual Vertex**
            get_cell_vertices( moris_index aCellIndex );

            //--------------------------------------------------------------------------------------------------------------

            virtual Vertex_Interpolation**
            get_vertex_interpolation( moris_index aVertexIndex );

            //--------------------------------------------------------------------------------------------------------------

            virtual mtk::Cell*
            get_ip_cell_in_cluster( enum ClusterType aClusterType,
                    moris_index                      aClusterIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

            virtual mtk::Cell* const *
            get_ig_cells_in_cluster( enum ClusterType aClusterType,
                    Primary_Void            aPrimaryOrVoid,
                    moris_index                       aClusterIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

            virtual uint
            get_num_cells_in_cluster( enum ClusterType aClusterType,
                    Primary_Void             aPrimaryOrVoid,
                    moris_index                        aClusterIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

            virtual moris_index*
            get_side_ordinals_in_cluster( enum ClusterType aClusterType,
                    moris_index                            aClusterIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

            virtual bool
            cluster_is_trivial( enum ClusterType aClusterType,
                    moris_index                  aClusterIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

            virtual Vertex* const *
            get_vertices_in_cluster( enum ClusterType aClusterType,
                    moris_index                       aClusterIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

            virtual uint
            get_num_vertices_in_cluster( enum ClusterType aClusterType,
                    moris_index                           aClusterIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

            virtual Matrix< DDRMat >*
            get_local_coord_matrix_ptr( enum ClusterType aClusterType,
                    moris_index                          aClusterIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

            virtual uint
            get_row_number_local_coords_matrix( enum ClusterType aClusterType,
                    moris_index                                  aClusterIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

            virtual mtk::Cell_Cluster const *
            get_associated_cell_cluster( moris_index aClusterIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

            virtual Matrix< IndexMat >
            get_enriched_mesh_indices() const;

            //--------------------------------------------------------------------------------------------------------------

            virtual bool
            is_secondary_cluster( moris_index aClusterIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

            virtual std::shared_ptr< mtk::Cell_Info >
            get_cell_info_sp( moris_index aEntityIndex ) const;

            //--------------------------------------------------------------------------------------------------------------

        }; // class mtk::Mesh

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace mtk
}    // namespace moris

//--------------------------------------------------------------------------------------------------------------

#endif /* MORIS_CL_MTK_MESH_CORE_HPP_ */

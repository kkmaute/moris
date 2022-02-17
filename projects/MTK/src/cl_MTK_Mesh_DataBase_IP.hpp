/*
 * cl_MTK_Mesh_DataBase_IP.hpp
 *
 *  Created on: Jan  24, 2022
 *      Author: momo
 */
#ifndef SRC_cl_MTK_Mesh_DataBase_IP
#define SRC_cl_MTK_Mesh_DataBase_IP

#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_TOL_Memory_Map.hpp"

namespace moris
{
    namespace mtk
    {
        // Forward declratios to
        class Interpolation_Mesh_DataBase;
        class Vertex_DataBase;
        class Vertex_Interpolation_DataBase;
        class Vertex_Interpolation;
        class Vertex;
        class Cell_DataBase;


        class Interpolation_Mesh_DataBase_IP : public mtk::Interpolation_Mesh
        {
          private:
            // Interpolation_Mesh_DataBase* mIPDataBase;
            // mtk::Interpolation_Mesh&     mIPMesh;
            uint                         mSpatilDim;

            // Vertex Information
            moris::Cell< Vertex_DataBase > mVertices;

            // coordinates matrix to store all the coordinates
            Matrix< moris::DDRMat > mVertexCoordinates;

            // vertex ip info
            moris::Cell< Vertex_Interpolation_DataBase > mVertexInterpoltions;
            moris::Cell< Vertex_Interpolation* >         mVertexInterpoltionsPtrs;

            // outer cell is the bspline mesh index
            // inner cell is the combined nodal t-matrices
            moris::Cell< moris::Cell< size_t > >             mOffSetTMatrix;
            moris::Cell< moris::Cell< moris::real > >        mWeights;
            moris::Cell< moris::Cell< moris::moris_id > >    mBasisIds;
            moris::Cell< moris::Cell< moris::moris_index > > mBasisOwners;
            moris::Cell< moris::Cell< moris::moris_index > > mBasisIndices;

            // Cell to vertex connectivity
            moris::Cell< mtk::Vertex* > mCellToVertices;
            moris::Cell< moris_index >  mCellToVertexOffSet;

            // Cell Information
            moris::Cell< Cell_DataBase >             mCells;
            std::shared_ptr< moris::mtk::Cell_Info > mCellInfo;

            // Communication table
            Matrix< IdMat > mCommunicationTable;

            // Adof Map (global to local map for enriched vertices)
            moris::Cell< map< moris_id, moris_index > > mAdofMap;

            // vertex map ( used in GEN)
            std::unordered_map< moris_id, moris_index > mVertexGlobalIdToLocalIndex;

            // id and owner information
            moris::Cell< moris_id > mVertexIdList;
            moris::Cell< moris_id > mCellIdList;
            moris::Cell< moris_id > mVertexOwnerList;
            moris::Cell< moris_id > mCellOwnerList;

            Matrix< IndexMat > mMeshIndices;

          public:
            // ----------------------------------------------------------------------------

            /**
             * @brief Construct a new Interpolation_Mesh_Analysis object
             *
             * @param aIPDataBase an IP mesh data based containing the raw data
             * @param aIPMesh an already existing IP mesh
             */
            Interpolation_Mesh_DataBase_IP();

            // ----------------------------------------------------------------------------

            /**
             * @brief Destroy the Interpolation_Mesh_Analysis object
             *
             */

            ~Interpolation_Mesh_DataBase_IP();

             // ----------------------------------------------------------------------------

            /**
             * @brief Get the communication table object
             *
             * @return Matrix< IdMat >  a coomunication table
             */

            virtual Matrix< IdMat >
            get_communication_table() const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the spatial dim
             *
             * @return uint the spatial dimension of the mesh
             */

            virtual uint get_spatial_dim() const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the num entities object
             *
             * @param aEntityRank the entity rank (node, element)
             * @param aIndex the index of the mesh
             * @return uint
             */

            virtual uint get_num_entities(
                enum EntityRank   aEntityRank,
                const moris_index aIndex = 0 ) const override;
            // ----------------------------------------------------------------------------

            /**
             * @brief Get the mtk vertex object / non-writable
             *
             * @param aVertexIndex index of the vertex
             * @return Vertex const&
             */

            virtual Vertex const& get_mtk_vertex( moris_index aVertexIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the mtk vertex object / writable
             *
             * @param aVertexIndex index of the vertex
             * @return Vertex& a constant vertex object
             */

            virtual Vertex& get_mtk_vertex( moris_index aVertexIndex ) override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the mtk vertex object / non-writable
             *
             * @param aVertexIndex index of the vertex
             * @return Vertex const&
             */

            virtual mtk::Cell& get_mtk_cell( moris_index aCellIndex ) override;


            // moris::Cell<mtk::Cell_DataBase> const & get_mtk_cells() const;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the mtk vertex object / non-writable
             *
             * @param aVertexIndex index of the vertex
             * @return Vertex const&
             */

            virtual mtk::Cell const& get_mtk_cell( moris_index aCellIndex ) const override;


            // ----------------------------------------------------------------------------

            /**
             * @brief Get the adof map for specific mesh
             *
             * @param aBSplineIndex mesh index
             * @param aAdofMap the global to local adof map
             */
            void virtual get_adof_map(
                const uint                    aBSplineIndex,
                map< moris_id, moris_index >& aAdofMap ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the vertex glb id to loc vertex ind map object , used in gen
             *
             * @return std::unordered_map< moris_id, moris_index >
             */

            virtual std::unordered_map< moris_id, moris_index >
            get_vertex_glb_id_to_loc_vertex_ind_map() const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the num interpolations
             *
             * @return uint number of bspline meshes
             */
            virtual uint get_num_interpolations() override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the max num coeffs on proc
             *
             * @param aDiscretizationIndex bspline mesh index
             * @return uint max number of adof on a processor
             */
            virtual uint
            get_max_num_coeffs_on_proc( uint aDiscretizationIndex ) const override;


            // /////////////////------------------------------------------------------------
            // Non Used Mesh Function (so far)
            // /////////////////------------------------------------------------------------

            /**
             * @brief Get the mesh type
             *
             * @return MeshType
             */

            virtual MeshType get_mesh_type() const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the entity connected to entity loc inds
             *
             * @param aEntityIndex
             * @param aInputEntityRank
             * @param aOutputEntityRank
             * @param aDiscretizationIndex
             * @return Matrix< IndexMat >
             */

            Matrix< IndexMat >
            get_entity_connected_to_entity_loc_inds(
                moris_index       aEntityIndex,
                enum EntityRank   aInputEntityRank,
                enum EntityRank   aOutputEntityRank,
                const moris_index aDiscretizationIndex = 0 ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the node coordinate
             *
             * @param aNodeIndex
             * @return Matrix< DDRMat >
             */

            Matrix< DDRMat >
            get_node_coordinate( moris_index aNodeIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the glb entity id from entity loc index object
             *
             * @param aEntityIndex
             * @param aEntityRank
             * @param aDiscretizationIndex
             * @return moris_id
             */

            virtual moris_id
            get_glb_entity_id_from_entity_loc_index(
                moris_index       aEntityIndex,
                enum EntityRank   aEntityRank,
                const moris_index aDiscretizationIndex = 0 ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the node owner
             *
             * @param aNodeIndex
             * @return uint
             */

            virtual uint get_node_owner( moris_index aNodeIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the element indices in block set
             *
             * @param aSetIndex
             * @return Matrix< IndexMat >
             */

            virtual Matrix< IndexMat > get_element_indices_in_block_set( uint aSetIndex ) override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the blockset topology
             *
             * @param aSetName
             * @return enum CellTopology
             */

            virtual enum CellTopology
            get_blockset_topology( const std::string& aSetName ) override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the IG blockset shape
             *
             * @param aSetName
             * @return enum CellShape
             */

            virtual enum CellShape
            get_IG_blockset_shape( const std::string& aSetName ) override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the IP blockset shape
             *
             * @param aSetName
             * @return enum CellShape
             */

            virtual enum CellShape
            get_IP_blockset_shape( const std::string& aSetName ) override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the element owner
             *
             * @param aElementIndex
             * @return uint
             */

            virtual uint get_element_owner( moris_index aElementIndex ) const override;

            // ----------------------------------------------------------------------------

            /**
             * @brief Get the elements connected to element and face ind loc inds
             *
             * @param aElementIndex
             * @return Matrix< IndexMat >
             */

            virtual Matrix< IndexMat >
            get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const override;

            // -------------------------------------------------------------------------------

            /////////////////////////
            // Accessor functions for mesh entities
            /////////////////////////

            /**
             * @brief Get the entity id
             *
             * @param aEntityRank
             * @param aEntityIndex
             * @return moris_id const
             */

            virtual moris_id
            get_entity_id( enum EntityRank aEntityRank, moris_index aEntityIndex ) const override;

            // -------------------------------------------------------------------------------

            /**
             * @brief Get the entity owner
             *
             * @param aEntityRank
             * @param aEntityIndex
             * @return moris_id const
             */

            virtual moris_id
            get_entity_owner( enum EntityRank aEntityRank, moris_index aEntityIndex ) const override;

            // -------------------------------------------------------------------------------
            virtual std::shared_ptr< mtk::Cell_Info >
            get_cell_info_sp( moris_index aEntityIndex ) const override;

            /**
             * @brief Get the basis ids
             *
             * @param aVertexIndex
             * @param aOrder
             * @return moris_id* const
             */

            virtual moris_id* 
            get_basis_ids( moris_index aVertexIndex, moris_index aOrder ) override;

            // -------------------------------------------------------------------------------

            /**
             * @brief Get the basis indicies
             *
             * @param aVertexIndex
             * @param aOrder
             * @return moris_index* const
             */

            virtual moris_index*
            get_basis_indicies( moris_index aVertexIndex, moris_index aOrder ) override;

            // -------------------------------------------------------------------------------

            /**
             * @brief Get the basis weights
             *
             * @param aVertexIndex
             * @param aOrder
             * @return real* const
             */

            virtual real* 
            get_basis_weights( moris_index aVertexIndex, moris_index aOrder ) override;

            // -------------------------------------------------------------------------------

            /**
             * @brief Get the basis owners
             *
             * @param aVertexIndex
             * @param aOrder
             * @return moris_id* const
             */

            virtual moris_id*
            get_basis_owners( moris_index aVertexIndex, moris_index aOrder )  override;

            // -------------------------------------------------------------------------------

            /**
             * @brief Get the basis length
             *
             * @param aVertexIndex
             * @param aOrder
             * @return moris_index const
             */

            virtual moris_index
            get_basis_length( moris_index aVertexIndex, moris_index aOrder ) override;

            // -------------------------------------------------------------------------------

            /**
             * @brief Get the vertex coords ptr
             *
             * @param aVertexIndex
             * @return moris::real* const
             */

            virtual moris::real*
            get_vertex_coords_ptr( moris_index aVertexIndex );

            // -------------------------------------------------------------------------------

            /**
             * @brief Get the cell vertices object
             *
             * @param aCellIndex
             * @return Vertex**
             */

            virtual Vertex**
            get_cell_vertices( moris_index aCellIndex ) override;

            // -------------------------------------------------------------------------------

            /**
             * @brief Get the vertex interpolation object
             *
             * @param aVertexIndex
             * @return Vertex_Interpolation**
             */

            virtual Vertex_Interpolation**
            get_vertex_interpolation( moris_index aVertexIndex ) override;

            // -------------------------------------------------------------------------------

            /**
             * @brief Get the memory usage
             *
             * @return moris::Memory_Map
             */
            moris::Memory_Map
            get_memory_usage();

            // ----------------------------------------------------------------------------

            /**
             * @brief delete unused data
             *
             */

            void
            free_memory();

            // ----------------------------------------------------------------------------

            friend class Interpolation_Mesh_Editor;
        };
    }// namespace mtk
}// namespace moris


#endif /* cl_MTK_Mesh_DataBase_IP.hpp */
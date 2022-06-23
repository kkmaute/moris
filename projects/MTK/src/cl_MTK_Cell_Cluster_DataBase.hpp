/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_MTK_Cell_Cluster_DataBase.hpp  
 * 
 */
#ifndef SRC_cl_MTK_Cell_Cluster_DataBase
#define SRC_cl_MTK_Cell_Cluster_DataBase

#include "cl_MTK_Cell_Cluster.hpp"

namespace moris
{
    namespace mtk
    {
        // forward declare the classes
        class Cell_DataBase;
        class Mesh;

        class Cell_Cluster_DataBase : public Cell_Cluster
        {
        private:
            // FIXME: old data member that needs to be deleted
            moris::Cell< moris::mtk::Cell const* > mPrimaryIntegrationCells;
            moris::Cell< moris::mtk::Cell const* > mVoidIntegrationCells;

            moris_index mCellClusterIndex;// cell cluster index
            mtk::Mesh*  mMesh;// mesh pointer

        public:
            //------------------------------------------------------------------------------

            /**
             * @brief Construct a new Cell_Cluster_DataBase object
             *
             */

            Cell_Cluster_DataBase() = default;

            //------------------------------------------------------------------------------

            /**
             * @brief Construct a new Cell_Cluster_DataBase object
             *
             * @param aCellClusterIndex an index identifying the cluster
             * @param aMesh the mesh pointer
             */

            Cell_Cluster_DataBase( 
                    moris_index aCellClusterIndex,
                    mtk::Mesh*  aMesh );


            //------------------------------------------------------------------------------

            /**
             * @brief Destroy the Cell_Cluster_DataBase object
             *
             */

            ~Cell_Cluster_DataBase() = default;


            //------------------------------------------------------------------------------

            /**
             * @brief
             *
             * @param aIsMaster is master or slave , typically not used in cell cluster
             * @return true master side
             * @return false slave side
             */

            virtual 
            bool
            is_trivial( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const override;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the primary cells in cluster object
             *
             * @param aIsMaster master or slave side , typically not used in cell cluster
             * @return moris::Cell< moris::mtk::Cell const* > const&  cell containing primary cells
             */

            virtual 
            moris::Cell< moris::mtk::Cell const* > const&
            get_primary_cells_in_cluster( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const override;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the void cells in cluster object
             *
             * @return moris::Cell< moris::mtk::Cell const* > const& cell containing void cells
             */

            virtual 
            moris::Cell< moris::mtk::Cell const* > const&
            get_void_cells_in_cluster() const override;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the interpolation cell object
             *
             * @param aIsMaster master or slave side , typically not used in cell cluster
             * @return moris::mtk::Cell const& interpolation cell of the cluster(this is unique)
             */

            virtual 
            moris::mtk::Cell const&
            get_interpolation_cell( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const override;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the vertices in cluster object
             *
             * @param aIsMaster master or slave side , typically not used in cell cluster
             * @return moris::Cell< moris::mtk::Vertex const* > get vertices in the cluster
             */

            virtual 
            moris::Cell< moris::mtk::Vertex const* >
            get_vertices_in_cluster( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const override;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the vertices local coordinates wrt interp cell object
             *
             * @param aIsMaster master or slave side , typically not used in cell cluster
             * @return moris::Matrix< moris::DDRMat > coordinates of the vertices (nVert*nDim)
             */

            virtual 
            moris::Matrix< moris::DDRMat >
            get_vertices_local_coordinates_wrt_interp_cell( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const override;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the vertex local coordinate wrt interp cell object
             *
             * @param aVertex a vertex that is in the cluster
             * @param aIsMaster master or slave side , typically not used in cell cluster
             * @return moris::Matrix< moris::DDRMat > coordinates of the vertex (1*nDim)
             */

            virtual 
            moris::Matrix< moris::DDRMat >
            get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const* aVertex,
                const mtk::Master_Slave                                            aIsMaster = mtk::Master_Slave::MASTER ) const override;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the dim of param coord object
             *
             * @param aIsMaster master or slave side , typically not used in cell cluster
             * @return moris_index spatial dimension of the probelm
             */

            virtual 
            moris_index
            get_dim_of_param_coord( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const override;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the primary cell local coords on side wrt interp cell object
             *
             * @param aPrimaryCellClusterIndex a local cell index in the cluster
             * @return moris::Matrix< moris::DDRMat > coordinates of the cell in the cluster (nVert*nDim)
             */

            virtual 
            moris::Matrix< moris::DDRMat >
            get_primary_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aPrimaryCellClusterIndex ) const override;

            //------------------------------------------------------------------------------

            virtual
            moris::real
            compute_cluster_group_cell_measure(
                    const moris_index       aDiscretizationMeshIndex,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER) const override
            {
                MORIS_ERROR( false, "mtk::Side_Cluster_DataBase::compute_cluster_group_cell_measure() - Only implemented in child mtk::Cell_Cluster_DataBase" );
                return 0.0;
            }

            //------------------------------------------------------------------------------

            virtual
            moris::real
            compute_cluster_group_cell_measure_derivative(
                    const moris_index       aDiscretizationMeshIndex,
                    const Matrix< DDRMat > & aPerturbedVertexCoords,
                    uint aDirection,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const override
            {
                MORIS_ERROR( false, "mtk::Side_Cluster_DataBase::compute_cluster_group_cell_measure_derivative() - Only implemented in child mtk::Cell_Cluster_DataBase" );
                return 0.0;
            }

            //------------------------------------------------------------------------------

            virtual
            moris::real
            compute_cluster_group_cell_side_measure(
                    const moris_index       aDiscretizationMeshIndex,
                    const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                    const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const override
            {
                MORIS_ERROR( false, "mtk::Cell_Cluster_DataBase::compute_cluster_group_cell_side_measure() - Only valid on mtk::Side_Cluster_DataBase" );
                return 0.0;
            }

            //------------------------------------------------------------------------------

            virtual
            moris::real
            compute_cluster_group_cell_side_measure_derivative(
                    const moris_index       aDiscretizationMeshIndex,
                    const Matrix< DDRMat >& aPerturbedVertexCoords,
                    uint aDirection,
                    const mtk::Primary_Void aPrimaryOrVoid,
                    const mtk::Master_Slave aIsMaster ) const override
            {
                MORIS_ERROR( false, "mtk::Cell_Cluster_DataBase::compute_cluster_group_cell_side_measure_derivative() - Only valid on mtk::Side_Cluster_DataBase" );
                return 0.0;
            }

            //------------------------------------------------------------------------------

            /**
             * @brief populate the data that is returned by refernce
             *
             */

            void
            set_outward_data();

            //------------------------------------------------------------------------------

            /**
             * @brief return the memory of the cell cluster
             *
             * @return size_t
             */

            size_t
            capacity();

            //------------------------------------------------------------------------------
        
        }; // class Cell_Cluster_DataBase
        
    } // namespace mtk
} // namespace moris


#endif /* cl_MTK_Cell_Cluster_DataBase.hpp */
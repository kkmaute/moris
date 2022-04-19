/* 
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_MTK_Double_Side_Cluster.hpp  
 * 
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_DOUBLE_SIDE_CLUSTER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_DOUBLE_SIDE_CLUSTER_HPP_

#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Side_Cluster.hpp"

#include "cl_MTK_Cluster.hpp"

namespace moris
{
    namespace mtk
    {
        class Double_Side_Cluster : public Cluster
        {
                /*
                 * An assumption here is made that the first side appearing in the master side cluster is paired
                 * with the first side appearing in the slave side cluster.
                 */

                // Master side cluster
                moris::mtk::Cluster const * mMasterSideCluster;

                // Slave side cluster
                moris::mtk::Cluster const * mSlaveSideCluster;


                /*!
                 * A one way pairing from master vertices to slave vertices
                 */
                moris::Cell<moris::mtk::Vertex const *> mMasterToSlaveVertexPairs;

            public:
                // ----------------------------------------------------------------------------------

                /**
                 * Default constructor, both master and slave side clusters are intialized
                 * as nullptr.
                 */
                Double_Side_Cluster();

                // ----------------------------------------------------------------------------------


                /**
                 * Default virtual destructor
                 */
                virtual
                ~Double_Side_Cluster();

                // ----------------------------------------------------------------------------------

                /**
                 * Constructor which should be used if to construct operational double side cluster
                 * @param[in] aMasterSideCluster Master side cluster
                 * @param[in] aSlaveSideCluster Slave side cluster
                 * @param[in] aMasterToSlaveVertexPair) Vertices on the sides of the clusters
                 */
                Double_Side_Cluster(
                        moris::mtk::Cluster const *                           aMasterSideCluster,
                        moris::mtk::Cluster const *                           aSlaveSideCluster,
                        moris::Cell<moris::mtk::Vertex const *> const & aMasterToSlaveVertexPair);

                //##############################################
                // Side Cluster traits access
                //##############################################

                /*!
                 * Ask a cluster in the double side cluster whether it is trivial or not. A trivial
                 * cluster has a 1-to-1 relationship between integration entities and interpolation entities.
                 * This does not mean there is the same id for each entity.
                 * @param[in] aIsMaster enum specifying which cluster the question is for
                 * @return True if trivial cluster
                 */
                //virtual
                bool
                is_trivial( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;

                /*!
                 * @return true if master cluster is trivial
                 */
                bool
                is_master_trivial() const;

                /*!
                 * @return true if slave cluster is trivial
                 */
                bool
                is_slave_trivial() const;

                //##############################################
                // Single Side Cluster Access
                //##############################################

                /*!
                 * @return const master side cluster
                 */
                moris::mtk::Cluster const &
                get_master_side_cluster() const;

                //----------------------------------------------------------------

                /*!
                 * @return const slave side cluster
                 */
                moris::mtk::Cluster const &
                get_slave_side_cluster() const;

                //----------------------------------------------------------------

                moris::mtk::Cluster const &
                get_cluster(
                        const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;


                //##############################################
                // Vertex Pair Access
                //##############################################
                /*!
                 * @param[in] aMasterVertex A vertex on the master side of double side cluster
                 * @return Corresponding slave vertex
                 */
                moris::mtk::Vertex const *
                get_master_vertex_pair(moris::mtk::Vertex const * aMasterVertex) const;


                //----------------------------------------------------------------

                //##############################################
                // Cell Side Ordinals/Vertex Access
                //##############################################

                /*
                 * @param[in] aVertex Vertex in mesh
                 * @param[in] aIsMaster Master/Slave selector enum
                 * @return Vertex local cluster index wrt master/slave
                 */
                moris_index
                get_vertex_cluster_index(
                        const Vertex *          aVertex,
                        const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;

                //----------------------------------------------------------------

                /*
                 * @param[in] aIsMaster Master/Slave selector enum
                 * @return Interpolation cell of the side cluster.
                 */
                moris::mtk::Cell const &
                get_interpolation_cell( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;


                /*!
                 * @return Interpolation cell of master side cluster
                 */
                moris::mtk::Cell const &
                get_master_interpolation_cell() const;
                //----------------------------------------------------------------

                /*!
                 * @return Interpolation cell of slave side cluster
                 */
                moris::mtk::Cell const &
                get_slave_interpolation_cell() const;

                //----------------------------------------------------------------

                /*!
                 * @param[in] aIsMaster Master/Slave selector enum
                 * @return Primary integration cells in the side cluster
                 */
                moris::Cell<mtk::Cell const *> const &
                get_primary_cells_in_cluster( const mtk::Master_Slave aIsMaster ) const;

                /*!
                 * @return Primary integration cells in the master side cluster
                 */
                moris::Cell<mtk::Cell const *> const &
                get_master_integration_cells() const;

                //----------------------------------------------------------------

                /*!
                 * @return Primary integration cells in the master side cluster
                 */
                moris::Cell<mtk::Cell const *> const &
                get_slave_integration_cells() const;

                //----------------------------------------------------------------
                /*
                 * @param[in] aIsMaster Master/Slave selector enum
                 * @return All integration cell side ordinals
                 */
                moris::Matrix<moris::IndexMat>
                get_cell_side_ordinals( const mtk::Master_Slave aIsMaster ) const;

                //----------------------------------------------------------------
                /*
                 *
                 */
                moris_index
                get_cell_side_ordinal(
                        moris::moris_index aCellIndexInCluster,
                        const mtk::Master_Slave aIsMaster ) const;

                /*!
                 * @return all integration cell side ordinals on master side of the
                 * double sided side cluster
                 */
                moris::Matrix<moris::IndexMat>
                get_master_integration_cell_side_ordinals() const;

                //----------------------------------------------------------------

                /*!
                 * @return all integration cell side ordinals on master side of the
                 * double sided side cluster
                 */
                moris_index
                get_master_cell_side_ordinal(moris::moris_index aMasterCellIndexInCluster) const;
                //----------------------------------------------------------------

                /*!
                 * @return all integration cell side ordinals on slave side of the
                 * double sided side cluster
                 */
                moris::Matrix<moris::IndexMat>
                get_slave_integration_cell_side_ordinals() const;

                //----------------------------------------------------------------

                /*!
                 * Single side ordinal version of above
                 */

                moris_index
                get_slave_cell_side_ordinal(moris::moris_index aSlaveCellIndexInCluster) const;

                //----------------------------------------------------------------

                /*!
                 * @return all the master vertices in this cluster
                 */
                moris::Cell<moris::mtk::Vertex const *>
                get_vertices_in_cluster( const mtk::Master_Slave aIsMaster ) const;



                moris::Cell<moris::mtk::Vertex const *>
                get_master_vertices_in_cluster() const;

                //----------------------------------------------------------------

                /*!
                 * Returns all the slave vertices in this cluster
                 */

                moris::Cell<moris::mtk::Vertex const *>
                get_slave_vertices_in_cluster() const;

                //----------------------------------------------------------------

                moris_index
                get_slave_vertex_ord_on_facet( moris_index  aCellClusterIndex,
                        moris::mtk::Vertex const * aSlaveVertex) const;

                //----------------------------------------------------------------

                moris::Matrix<moris::IndexMat>
                get_master_vertex_indices_in_cluster() const;

                //----------------------------------------------------------------

                moris::Matrix<moris::IndexMat>
                get_slave_vertex_indices_in_cluster() const;

                //----------------------------------------------------------------

                //virtual
                moris::Matrix<moris::IndexMat>
                get_vertex_indices_in_cluster() const;

                //----------------------------------------------------------------

                //##############################################
                // Local Coordinate Access
                //##############################################


                //----------------------------------------------------------------

                /*
                 * Access the full array of local coordinates on the master
                 */

                moris::Matrix<moris::DDRMat>
                get_vertices_local_coordinates_wrt_interp_cell( const mtk::Master_Slave aIsMaster ) const;

                //----------------------------------------------------------------

                moris::Matrix<moris::DDRMat>
                get_master_vertices_local_coordinates_wrt_interp_cell() const;

                //----------------------------------------------------------------

                /*
                 * Access the full array of local coordinates on the slave
                 */
                moris::Matrix<moris::DDRMat>
                get_slave_vertices_local_coordinates_wrt_interp_cell() const;

                //----------------------------------------------------------------

                moris::Matrix<moris::DDRMat>
                get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex,
                        const mtk::Master_Slave aIsMaster ) const;
                //----------------------------------------------------------------

                /*
                 * Access a single local coordinate of a vertex on the master
                 */
                moris::Matrix<moris::DDRMat>
                get_master_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex ) const;

                //----------------------------------------------------------------

                /*
                 * Access a single local coordinate of a vertex on the slave
                 */
                moris::Matrix<moris::DDRMat>
                get_slave_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex ) const;

                //----------------------------------------------------------------

                moris::Matrix<moris::DDRMat>
                get_cell_local_coords_on_side_wrt_interp_cell(moris::moris_index aClusterLocalIndex,
                        const mtk::Master_Slave aIsMaster ) const;

                //----------------------------------------------------------------

                /*!
                 * Access an integration cells parametric coordinates on a side master side
                 * @param[in] - Local integration cell index with respect to the cluster (not proc local index)
                 */
                moris::Matrix<moris::DDRMat>
                get_master_cell_local_coords_on_side_wrt_interp_cell(moris::moris_index aMasterClusterLocalIndex) const;

                //----------------------------------------------------------------

                /*!
                 * Access an integration cells parametric coordinates on a side slave side
                 * @param[in] - Local integration cell index with respect to the cluster (not proc local index)
                 */
                moris::Matrix<moris::DDRMat>
                get_slave_cell_local_coords_on_side_wrt_interp_cell(moris::moris_index aSlaveClusterLocalIndex) const;

                //----------------------------------------------------------------

                //##############################################
                // Size Access
                //##############################################
                /*!
                 * Size of the xsi vector of master
                 */
                moris_index
                get_dim_of_param_coord( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;

                //----------------------------------------------------------------

                moris::real
                compute_cluster_cell_measure(
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER) const;

                //----------------------------------------------------------------

                Matrix<DDRMat>
                compute_cluster_ig_cell_measures(
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER) const;

                //----------------------------------------------------------------

                moris::real
                compute_cluster_cell_measure_derivative(
                        const Matrix< DDRMat > & aPerturbedVertexCoords,
                        uint aDirection,
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const;

                //----------------------------------------------------------------

                moris::real
                compute_cluster_cell_side_measure(
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER) const;

                // ----------------------------------------------------------------------------------

                Matrix<DDRMat>
                compute_cluster_ig_cell_side_measures(
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER) const;

                //----------------------------------------------------------------

                moris::real
                compute_cluster_cell_side_measure_derivative(
                        const Matrix< DDRMat > & aPerturbedVertexCoords,
                        uint aDirection,
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const;

                //----------------------------------------------------------------

                moris_index
                get_master_dim_of_param_coord() const;

                //----------------------------------------------------------------

                /*!
                 * Size of the xsi vector of slave
                 */
                moris_index
                get_slave_dim_of_param_coord() const;
                //----------------------------------------------------------------

                moris::uint
                get_master_num_vertices_in_cluster() const;

                //----------------------------------------------------------------

                moris::uint
                get_slave_num_vertices_in_cluster() const;
                
                //----------------------------------------------------------------
                
                /**
                 * @brief Get the master vertex pairs cell
                 * 
                 * @return moris::Cell<moris::mtk::Vertex const *> const& a cell of the master side vertices
                 */
                 moris::Cell<moris::mtk::Vertex const *> const &
                get_master_vertex_pairs() const;

                //----------------------------------------------------------------
                
                /**
                 * @brief memory usage of the double sided cluster
                 * 
                 * @return size_t 
                 */

                 size_t
                 capacity();
                 
        };

        inline
        std::ostream &
        operator<<(std::ostream & os, const Double_Side_Cluster & dt)
        {
            os<<"\n  Master Interpolation Cell: "<<std::setw(9)<<dt.get_master_side_cluster().get_interpolation_cell().get_id();
            os<<" | Slave Interpolation Cell: "<<std::setw(9)<<dt.get_slave_side_cluster().get_interpolation_cell().get_id();


            moris::Cell<mtk::Cell const *> const & tMasterIGCells =dt.get_master_side_cluster().get_primary_cells_in_cluster( );
            moris::Matrix<moris::IndexMat>         tMasterIGCellSideOrds = dt.get_master_side_cluster().get_cell_side_ordinals();
            moris::Cell<mtk::Cell const *> const & tSlaveIGCells = dt.get_slave_side_cluster().get_primary_cells_in_cluster( );
            moris::Matrix<moris::IndexMat>         tSlaveIGCellSideOrds = dt.get_slave_side_cluster().get_cell_side_ordinals();

            os<<"       Cell Pairs: "<<std::endl;

            for(moris::uint  i = 0; i < tMasterIGCells.size(); i++)
            {
                std::cout<<"  Master Cell ID/Ord: "<<std::setw(9)<<tMasterIGCells(i)->get_id()<<std::setw(9)<<tMasterIGCellSideOrds(i);
                std::cout<<"  Slave Cell ID/Ord: "<<std::setw(9)<<tSlaveIGCells(i)->get_id()<<std::setw(9)<<tSlaveIGCellSideOrds(i)<<std::endl;
            }

            moris::print(dt.get_vertices_local_coordinates_wrt_interp_cell(mtk::Master_Slave::MASTER),"Master Local Coords");
            moris::print(dt.get_vertices_local_coordinates_wrt_interp_cell(mtk::Master_Slave::SLAVE),"Slave Local Coords");

            return os;
        }

        inline
        std::ostream &
        operator<<(std::ostream & os, Double_Side_Cluster const * const & dt)
        {
            os<<*dt;

            return os;
        }

    }
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_DOUBLE_SIDE_CLUSTER_HPP_ */

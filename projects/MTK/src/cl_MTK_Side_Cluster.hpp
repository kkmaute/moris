/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_MTK_Side_Cluster.hpp  
 * 
 */
#ifndef PROJECTS_MTK_SRC_CL_MTK_SIDE_CLUSTER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_SIDE_CLUSTER_HPP_

#include "typedefs.hpp"
#include "fn_assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Cell.hpp"

#include "cl_MTK_Cluster.hpp"

namespace moris
{
    namespace mtk
    {
        class Cluster_Group;

        class Side_Cluster : public Cluster
        {
            public:

                // ----------------------------------------------------------------------------------

                /*!
                 * Default constructor
                 */
                Side_Cluster();

                // ----------------------------------------------------------------------------------

                virtual
                ~Side_Cluster(){};


                // ----------------------------------------------------------------------------------
                // Cell Side Ordinals/Vertex Access
                // (Pure Virtual)
                // ----------------------------------------------------------------------------------
                /*!
                 * Indicates there is a 1 to 1 relationship between
                 * integration cell and interpolation cells in this cluster
                 */
                virtual
                bool
                is_trivial( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const  = 0;

                // ----------------------------------------------------------------------------------

                /*!
                 * Get interpolation cell interpolating into this side cluster
                 * @param[in] aIsMaster  Master or Slave Selector Enum (for Double side clusters only)
                 * @return Interpolation cell related to cluster
                 */
                virtual
                moris::mtk::Cell const &
                get_interpolation_cell( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const = 0;

                // ----------------------------------------------------------------------------------

                /*!
                 * Get all integration cells in this side cluster
                 */
                virtual
                moris::Cell<mtk::Cell const *> const &
                get_cells_in_side_cluster() const = 0;

                // ----------------------------------------------------------------------------------
                /*!
                 * @param[in] aIsMaster  Master or Slave Selector Enum (for Double side clusters only)
                 */
                moris::Cell<mtk::Cell const *> const &
                get_primary_cells_in_cluster( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;

                // ----------------------------------------------------------------------------------

                /*!
                 * @param[in] aIsMaster  Master or Slave Selector Enum (for Double side clusters only)
                 * @return all integration cell side ordinals
                 */
                virtual
                moris::Matrix<moris::IndexMat>
                get_cell_side_ordinals( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const  = 0;

                // ----------------------------------------------------------------------------------

                /*!
                 * Single side ordinal version of above
                 * @param[in] aCellIndexInCluster Integration cell cluster index
                 * @return single integration cell side ordinal
                 *
                 */
                virtual
                moris_index
                get_cell_side_ordinal(moris::moris_index aCellIndexInCluster,
                        const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const = 0;

                // ----------------------------------------------------------------------------------

                /*!
                 * @return all the vertices in this cluster
                 */
                virtual
                moris::Cell<moris::mtk::Vertex const *>
                get_vertices_in_cluster( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const = 0;


                // ----------------------------------------------------------------------------------
                // Local Coordinate Access
                // ----------------------------------------------------------------------------------
                /*!
                 * Access the full array of local coordinates
                 * @param[in] aIsMaster - Master or Slave Selector Enum (for Double side clusters only)
                 * @return All vertex local coordinates
                 */
                virtual
                moris::Matrix<moris::DDRMat>
                get_vertices_local_coordinates_wrt_interp_cell( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const = 0;


                /*!
                 * Access a single local coordinate of a vertex
                 * @param[in] aIsMaster - Master or Slave Selector Enum (for Double side clusters only)
                 * @return Single vertex local coordinates
                 */
                virtual
                moris::Matrix<moris::DDRMat>
                get_vertex_local_coordinate_wrt_interp_cell(
                        moris::mtk::Vertex const * aVertex,
                        const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER) const  = 0;


                //----------------------------------------------------------------

                virtual
                bool
                has_cluster_group( const moris_index aBsplineMeshListIndex ) const override
                {
                    MORIS_ERROR( false, "mtk::Side_Cluster::has_cluster_group() - not implemented for this class" );
                    return false;
                }

                //----------------------------------------------------------------

                virtual
                std::shared_ptr< Cluster_Group >
                get_cluster_group( const moris_index aBsplineMeshListIndex ) const override
                {
                    MORIS_ERROR( false, "mtk::Side_Cluster::get_cluster_group() - not implemented for this class" );
                    return nullptr;
                }

                // ----------------------------------------------------------------------------------

                virtual
                void
                set_cluster_group( 
                        const moris_index aBsplineMeshListIndex,
                        std::shared_ptr< Cluster_Group > aClusterGroupPtr ) override
                {
                    MORIS_ERROR( false, "mtk::Side_Cluster::set_cluster_group() - only implemented for child xtk::Side_Cluster" );
                }

                // ----------------------------------------------------------------------------------

                //##############################################
                // Size Access
                //##############################################
                // ----------------------------------------------------------------------------------
                /*!
                 * @return Size of the xsi vector in this side cluster
                 */
                virtual
                moris_index
                get_dim_of_param_coord( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const  = 0;

                // ----------------------------------------------------------------------------------
                /*!
                 * @param[in] aPrimaryOrVoid Primary or Void Integration Cell Selector Enum
                 * @param[in] aIsMaster Master or Slave Selector Enum (for Double side clusters only)
                 * @return Cell Measure
                 */
                moris::real
                compute_cluster_cell_measure(
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER) const;

                //------------------------------------------------------------------------------

                virtual
                moris::real
                compute_cluster_group_cell_measure(
                        const moris_index       aBsplineMeshListIndex,
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER) const
                {
                    MORIS_ERROR( false, "mtk::Side_Cluster::compute_cluster_group_cell_measure() - Only implemented in child xtk::Side_Cluster" );
                    return 0.0;
                }

                // ----------------------------------------------------------------------------------
                /*!
                 * @param[in] aPrimaryOrVoid Primary or Void Integration Cell Selector Enum
                 * @param[in] aIsMaster Master or Slave Selector Enum (for Double side clusters only)
                 * @return vector of IG cell Measures
                 */
                Matrix<DDRMat>
                compute_cluster_ig_cell_measures(
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER) const;

                // ----------------------------------------------------------------------------------
                /*!
                 * @param[in] aPerturbedVertexCoords coordinate of perturbed vertex
                 * @param[in] aDirection spatial direction for perturbation
                 * @param[in] aPrimaryOrVoid Primary or Void Integration Cell Selector Enum
                 * @param[in] aIsMaster Master or Slave Selector Enum (for Double side clusters only)
                 * @return Derivative of Cell Measure
                 */
                moris::real
                compute_cluster_cell_measure_derivative(
                        const Matrix< DDRMat > & aPerturbedVertexCoords,
                        uint aDirection,
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const;

                //------------------------------------------------------------------------------

                virtual
                moris::real
                compute_cluster_group_cell_measure_derivative(
                        const moris_index       aBsplineMeshListIndex,
                        const Matrix< DDRMat > & aPerturbedVertexCoords,
                        uint aDirection,
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
                {
                    MORIS_ERROR( false, "mtk::Side_Cluster::compute_cluster_group_cell_measure_derivative() - Only implemented in child xtk::Side_Cluster" );
                    return 0.0;
                }

                // ----------------------------------------------------------------------------------
                /*!
                 * @param[in] aPrimaryOrVoid Primary or Void Integration Cell Selector Enum
                 * @param[in] aIsMaster Master or Slave Selector Enum (for Double side clusters only)
                 * @return Cell Side Measure
                 */
                virtual
                moris::real
                compute_cluster_cell_side_measure(
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER) const;

                //------------------------------------------------------------------------------

                virtual
                moris::real
                compute_cluster_group_cell_side_measure(
                        const moris_index       aBsplineMeshListIndex,
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER) const
                {
                    MORIS_ERROR( false, "mtk::Side_Cluster::compute_cluster_group_cell_side_measure() - Only implemented in child xtk::Side_Cluster" );
                    return 0.0;
                }

                // ----------------------------------------------------------------------------------
                /*!
                 * @param[in] aPrimaryOrVoid Primary or Void Integration Cell Selector Enum
                 * @param[in] aIsMaster Master or Slave Selector Enum (for Double side clusters only)
                 * @return vector of individual side measures for each IG cell
                 */
                virtual
                Matrix<DDRMat>
                compute_cluster_ig_cell_side_measures(
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER) const;

                // ----------------------------------------------------------------------------------
                /*!
                 * @param[in] aPerturbedVertexCoords coordinate of perturbed vertex
                 * @param[in] aDirection spatial direction for perturbation
                 * @param[in] aPrimaryOrVoid Primary or Void Integration Cell Selector Enum
                 * @param[in] aIsMaster Master or Slave Selector Enum (for Double side clusters only)
                 * @return Derivative of Cell Side Measure
                 */
                virtual
                moris::real
                compute_cluster_cell_side_measure_derivative(
                        const Matrix< DDRMat > & aPerturbedVertexCoords,
                        uint aDirection,
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const;

                //------------------------------------------------------------------------------

                virtual
                moris::real
                compute_cluster_group_cell_side_measure_derivative(
                        const moris_index       aBsplineMeshListIndex,
                        const Matrix< DDRMat > & aPerturbedVertexCoords,
                        uint aDirection,
                        const mtk::Primary_Void aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                        const mtk::Master_Slave aIsMaster      = mtk::Master_Slave::MASTER ) const
                {
                    MORIS_ERROR( false, "mtk::Side_Cluster::compute_cluster_group_cell_side_measure_derivative() - Only implemented in child xtk::Side_Cluster" );
                    return 0.0;
                }

                // ----------------------------------------------------------------------------------
                /*!
                 * @return Matrix of integration cell indices in cluster
                 */
                virtual
                moris::Matrix<moris::IndexMat>
                get_cell_indices_in_cluster() const;

                // ----------------------------------------------------------------------------------

                /*!
                 * @return Matrix of primary integration cell indices in cluster
                 */
                moris::Matrix<moris::IndexMat>
                get_primary_cell_indices_in_cluster() const;

                // ----------------------------------------------------------------------------------
                /*!
                 * @return Interpolation cell index
                 */
                virtual
                moris::moris_index
                get_interpolation_cell_index() const;

                // ----------------------------------------------------------------------------------
                /*!
                 * @return Matrix of vertex indices in cluster
                 */
                virtual
                moris::Matrix<moris::IndexMat>
                get_vertex_indices_in_cluster() const;
                // ----------------------------------------------------------------------------------
                // Cell/Vertex Id Access
                // ----------------------------------------------------------------------------------
                /*!
                 * @return Matrix of integration cell ids in cluster
                 */
                virtual
                moris::Matrix<moris::IdMat>
                get_cell_ids_in_cluster() const;

                // ----------------------------------------------------------------------------------
                /*!
                 * @return Matrix of primary integration cell ids in cluster
                 */
                moris::Matrix<moris::IdMat>
                get_primary_cell_ids_in_cluster() const;

                // ----------------------------------------------------------------------------------
                /*!
                 * @return Matrix of vertex ids in cluster
                 */
                virtual
                moris::Matrix<moris::IndexMat>
                get_vertex_ids_in_cluster() const;

                // ----------------------------------------------------------------------------------
                // Local Coordinate access
                // ----------------------------------------------------------------------------------

                /*!
                 * @brief Access an integration cells parametric coordinates on a side.
                 * @param[in] aClusterLocalIndex Local integration cell index with respect to the cluster (not proc local index)
                 * @param[in] aIsMaster Master or Slave Selector Enum (for Double side clusters only)
                 * @return Local coordinates wrt interpolation cell of vertices on integration cell side
                 */
                virtual
                moris::Matrix<moris::DDRMat>
                get_cell_local_coords_on_side_wrt_interp_cell(
                        moris::moris_index      aClusterLocalIndex,
                        const mtk::Master_Slave aIsMaster =         mtk::Master_Slave::MASTER ) const;

                // ----------------------------------------------------------------------------------
                /*!
                 * @brief Access an integration cells parametric coordinates on a side.
                 * @param[in] aClusterLocalIndex Local integration cell index with respect to the cluster (not proc local index)
                 * @param[in] aIsMaster Master or Slave Selector Enum (for Double side clusters only)
                 * @return Local coordinates wrt interpolation cell of vertices on integration cell side (primary)
                 */
                virtual
                moris::Matrix<moris::DDRMat>
                get_primary_cell_local_coords_on_side_wrt_interp_cell(
                        moris::moris_index      aPrimaryCellClusterIndex,
                        const mtk::Master_Slave aIsMaster                = mtk::Master_Slave::MASTER ) const;


                // ----------------------------------------------------------------------------------
                // Size Access
                // ----------------------------------------------------------------------------------
                /*!
                 * @return Number of sides in side cluster
                 */
                virtual
                moris::uint
                get_num_sides_in_cluster() const;

                // ----------------------------------------------------------------------------------
                /*!
                 * @return Number of primary cells in cluster
                 */
                moris::uint
                get_num_primary_cells() const;

                // ----------------------------------------------------------------------------------
                /*!
                 * @return number of vertices in cluster
                 */
                virtual
                moris::uint
                get_num_vertices_in_cluster() const;

                // ----------------------------------------------------------------------------------

                // ----------------------------------------------------------------------------------
                /*!
                 * @return Matrix of vertex ids in cluster
                 */
                virtual
                moris::Matrix<moris::DDRMat>
                get_vertex_coords_in_cluster() const;

                // ----------------------------------------------------------------------------------

        };
    }
}



#endif /* PROJECTS_MTK_SRC_CL_MTK_SIDE_CLUSTER_HPP_ */

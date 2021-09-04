/*
 * cl_VIS_Side_Cluster_Visualization.hpp
 *
 *  Created on: Jul 27, 2021
 *      Author: momo
 */

#ifndef PROJECTS_FEM_VIS_SRC_CL_VIS_SIDE_CLUSTER_VISUALIZATION_HPP_
#define PROJECTS_FEM_VIS_SRC_CL_VIS_SIDE_CLUSTER_VISUALIZATION_HPP_

#include <unordered_map>


#include "cl_MTK_Side_Cluster.hpp"

namespace moris
{
    namespace mtk
    {
        class Cell;
        class Vertex;
    }

    namespace vis
    {
        class Side_Cluster_Visualization: public mtk::Side_Cluster
        {
            private:
                bool                                    mTrivial;
                moris::mtk::Cell const *                mInterpolationCell;
                moris::Cell<moris::mtk::Cell const *>   mPrimaryIntegrationCells;
                moris::Matrix<moris::IndexMat>          mIntegrationCellSideOrdinals;
                //moris::Cell<moris::mtk::Cell const *>   mVoidIntegrationCells;
                moris::Cell<moris::mtk::Vertex const *> mVerticesInCluster;
                moris::Matrix<moris::DDRMat>            mVertexParamCoords;

                // map from vertex id to local index
                std::unordered_map<moris_index,moris_index> mVertexIdToLocalIndex;


            public:
                //----------------------------------------------------------------

                Side_Cluster_Visualization():mTrivial( true ),
                mInterpolationCell( nullptr ),
                mPrimaryIntegrationCells( 0, nullptr ),
                mIntegrationCellSideOrdinals( 0, 0 ),
                mVerticesInCluster( 0, nullptr ),
                mVertexParamCoords( 0, 0 )
                {};
                //----------------------------------------------------------------

                ~Side_Cluster_Visualization(){}
                //----------------------------------------------------------------

                bool
                is_trivial( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;

                //##############################################
                // Add and setup of cluster
                //##############################################
                void mark_as_nontrivial();

                //----------------------------------------------------------------

                void set_interpolation_cell(moris::mtk::Cell const * aInterpCell);

                //----------------------------------------------------------------

                void add_primary_integration_cell(moris::mtk::Cell  const * aIntegrationCell);

                //----------------------------------------------------------------

                void add_primary_integration_cell(moris::Cell<moris::mtk::Cell  const *> const & aIntegrationCell);

                //----------------------------------------------------------------

                void add_integration_cell_side_ordinals( moris::Matrix<moris::IndexMat> aIntegrationCellSideOrdinals );

                //----------------------------------------------------------------

                void add_integration_cell_side_ordinals( moris_index aSideOrdinal );

                //----------------------------------------------------------------

                void add_vertex_to_cluster(moris::Cell<moris::mtk::Vertex const *> const & aVertex);

                //----------------------------------------------------------------

                void add_vertex_local_coordinates_wrt_interp_cell(moris::Matrix<moris::DDRMat> const & aLocalCoords);

                //----------------------------------------------------------------

                //##############################################
                // Required Access Functions
                //##############################################

                moris::mtk::Cell const &
                get_interpolation_cell( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER) const;

                //----------------------------------------------------------------

                moris::Cell<moris::mtk::Cell const *> const &
                get_cells_in_side_cluster() const;

                //----------------------------------------------------------------

                moris::Matrix<moris::IndexMat>
                get_cell_side_ordinals( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;

                //----------------------------------------------------------------

//                moris_index
//                get_cell_side_ordinal(moris::moris_index aCellIndexInCluster,
//                        const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;

                //----------------------------------------------------------------

                moris::Cell<moris::mtk::Vertex const *>
                get_vertices_in_cluster( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;

                //----------------------------------------------------------------

                moris::Matrix<moris::DDRMat>
                get_vertices_local_coordinates_wrt_interp_cell(const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER) const;

                //----------------------------------------------------------------

//                moris_index
//                get_vertex_cluster_index( moris::mtk::Vertex const * aVertex,
//                        const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER  ) const;

                //----------------------------------------------------------------

//                moris_index
//                get_vertex_ordinal_on_facet( moris_index aCellIndexInCluster, moris::mtk::Vertex const * aVertex ) const;

                //----------------------------------------------------------------

                moris::Matrix<moris::DDRMat>
                get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex,
                        const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER) const;

                //----------------------------------------------------------------

                moris_index
                get_dim_of_param_coord( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;

                //----------------------------------------------------------------

                moris_index
                get_vertex_cluster_local_index(moris_id aVertexId) const;

                //----------------------------------------------------------------

            private:
                void
                add_vertex_to_map(moris_id aVertexId,
                        moris_index aVertexLocalIndex);
        };

    }


}
#endif /* PROJECTS_FEM_VIS_SRC_CL_VIS_SIDE_CLUSTER_VISUALIZATION_HPP_ */

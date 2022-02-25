/*
 * cl_MTK_Double_Side_Cluster.cpp
 *
 *  Created on: Aug 3, 2020
 *      Author: kedo3694
 */
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "fn_TOL_Capacities.hpp"

//----------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
    {
        Double_Side_Cluster::Double_Side_Cluster()
        :mMasterSideCluster(nullptr),
         mSlaveSideCluster(nullptr)
        {

        }

        //----------------------------------------------------------------------------------

        Double_Side_Cluster::~Double_Side_Cluster(){}

        //----------------------------------------------------------------------------------

        Double_Side_Cluster::Double_Side_Cluster(
                moris::mtk::Cluster const *                      aMasterSideCluster,
                moris::mtk::Cluster const *                      aSlaveSideCluster,
                moris::Cell<moris::mtk::Vertex const *> const & aLeftToRightVertexPair )
        : mMasterSideCluster(aMasterSideCluster),
          mSlaveSideCluster(aSlaveSideCluster)
        {
            // This check prohibits the construction of double side interfaces between child meshes therefore it is being removed.
            // if(!this->is_master_trivial())
            // {
                // MORIS_ASSERT(this->get_master_num_vertices_in_cluster() == this->get_slave_num_vertices_in_cluster(),"Number of vertices mismatch in double cluster");
            // }

            mMasterToSlaveVertexPairs.append(aLeftToRightVertexPair);
        }

        //----------------------------------------------------------------------------------

        bool
        Double_Side_Cluster::is_trivial(
                const mtk::Master_Slave aIsMaster ) const
        {
            if ( aIsMaster == mtk::Master_Slave::MASTER )
            {
                return this->get_master_side_cluster().is_trivial();
            }
            else if( aIsMaster == mtk::Master_Slave::SLAVE )
            {
                return this->get_slave_side_cluster().is_trivial();
            }
            else
            {
                MORIS_ERROR(false, "is_trivial(): can only be MASTER or SLAVE");
                return false;
            }
        }

        //----------------------------------------------------------------------------------

        bool
        Double_Side_Cluster::is_master_trivial() const
        {
            return this->get_master_side_cluster().is_trivial();
        }

        //----------------------------------------------------------------------------------

        bool
        Double_Side_Cluster::is_slave_trivial() const
        {
            return this->get_slave_side_cluster().is_trivial();

        }

        //----------------------------------------------------------------------------------

        moris::mtk::Cluster const &
        Double_Side_Cluster::get_master_side_cluster() const
        {
            return *mMasterSideCluster;
        }

        //----------------------------------------------------------------------------------

        moris::mtk::Cluster const &
        Double_Side_Cluster::get_slave_side_cluster() const
        {
            return *mSlaveSideCluster;
        }

        //----------------------------------------------------------------------------------

        moris::mtk::Cluster const &
        Double_Side_Cluster::get_cluster(
                const mtk::Master_Slave aIsMaster ) const
        {
            if(aIsMaster == mtk::Master_Slave::MASTER)
            {
                return this->get_master_side_cluster();
            }
            else
            {
                return this->get_slave_side_cluster();
            }
        }


        //----------------------------------------------------------------------------------

        moris::mtk::Vertex const *
        Double_Side_Cluster::get_master_vertex_pair(
                moris::mtk::Vertex const * aMasterVertex) const
        {
            moris_index tMasterClusterIndex = this->get_master_side_cluster().get_vertex_cluster_index(aMasterVertex);

            MORIS_ASSERT(tMasterClusterIndex < (moris_index)mMasterToSlaveVertexPairs.size(),"Vertex index out of bounds in pairing.");

            return mMasterToSlaveVertexPairs(tMasterClusterIndex);
        }

        //----------------------------------------------------------------------------------

        moris::Cell<moris::mtk::Vertex const *> const &
        Double_Side_Cluster::get_master_vertex_pairs() const
        {
             return mMasterToSlaveVertexPairs;
        }


        //----------------------------------------------------------------------------------


        moris_index
        Double_Side_Cluster::get_vertex_cluster_index(
                const Vertex * aVertex,
                const mtk::Master_Slave aIsMaster ) const
        {
            if ( aIsMaster == mtk::Master_Slave::MASTER )
            {
                return this->get_master_side_cluster().get_vertex_cluster_index( aVertex );
            }
            else if( aIsMaster == mtk::Master_Slave::SLAVE )
            {
                return this->get_slave_side_cluster().get_vertex_cluster_index( aVertex );
            }
            else
            {
                MORIS_ERROR(false, "get_vertex_cluster_index(): can only be MASTER and SLAVE");
                return 0;
            }
        }


        //----------------------------------------------------------------------------------

        moris::mtk::Cell const &
        Double_Side_Cluster::get_interpolation_cell(
                const mtk::Master_Slave aIsMaster ) const
        {
            if ( aIsMaster == mtk::Master_Slave::MASTER )
            {
                return this->get_master_side_cluster().get_interpolation_cell();
            }
            else if( aIsMaster == mtk::Master_Slave::SLAVE )
            {
                return this->get_slave_side_cluster().get_interpolation_cell();
            }
            else
            {
                MORIS_ERROR(false, "get_interpolation_cell(): can only be 0 and 1");
                
                //This function will be never be used
                return this->get_master_side_cluster().get_interpolation_cell();
            }
        }

        //----------------------------------------------------------------------------------

        moris::mtk::Cell const &
        Double_Side_Cluster::get_master_interpolation_cell() const
        {
            return this->get_master_side_cluster().get_interpolation_cell();
        }


        //----------------------------------------------------------------------------------

        moris::mtk::Cell const &
        Double_Side_Cluster::get_slave_interpolation_cell() const
        {
            return this->get_slave_side_cluster().get_interpolation_cell();
        }

        //----------------------------------------------------------------------------------

        moris::Cell<mtk::Cell const *> const &
        Double_Side_Cluster::get_primary_cells_in_cluster(
                const mtk::Master_Slave aIsMaster ) const
        {
            if ( aIsMaster == mtk::Master_Slave::MASTER )
            {
                return this->get_master_side_cluster().get_primary_cells_in_cluster();
            }
            else if( aIsMaster ==  mtk::Master_Slave::SLAVE )
            {
                return this->get_slave_side_cluster().get_primary_cells_in_cluster();
            }
            else
            {
                MORIS_ERROR(false, "get_primary_cells_in_cluster(): can only be MASTER and SLAVE");
                
                //create a dummy cell that never will be generated
                moris::Cell< mtk::Cell const * > *tDummyCell = new moris::Cell< mtk::Cell const * >( 0 );
                return *tDummyCell;
            }
        }

        //----------------------------------------------------------------------------------

        moris::Cell<mtk::Cell const *> const &
        Double_Side_Cluster::get_master_integration_cells() const
        {
            return this->get_master_side_cluster().get_primary_cells_in_cluster();
        }

        //----------------------------------------------------------------------------------

        moris::Cell<mtk::Cell const *> const &
        Double_Side_Cluster::get_slave_integration_cells() const
        {
            return this->get_slave_side_cluster().get_primary_cells_in_cluster();
        }

        //----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Double_Side_Cluster::get_cell_side_ordinals(
                const mtk::Master_Slave aIsMaster ) const
        {
            if ( aIsMaster == mtk::Master_Slave::MASTER )
            {
                return this->get_master_side_cluster().get_cell_side_ordinals();
            }
            else if( aIsMaster == mtk::Master_Slave::SLAVE )
            {
                return this->get_slave_side_cluster().get_cell_side_ordinals();
            }
            else
            {
                MORIS_ERROR(false, "get_cell_side_ordinals(): can only be MASTER and SLAVE");
                return moris::Matrix<moris::IndexMat>(0,0);
            }
        }

        //----------------------------------------------------------------------------------

        moris_index
        Double_Side_Cluster::get_cell_side_ordinal(
                moris::moris_index aCellIndexInCluster,
                const mtk::Master_Slave aIsMaster ) const
        {
            if ( aIsMaster == mtk::Master_Slave::MASTER )
            {
                return this->get_master_side_cluster().get_cell_side_ordinal(aCellIndexInCluster);
            }
            else if( aIsMaster == mtk::Master_Slave::SLAVE )
            {
                return this->get_slave_side_cluster().get_cell_side_ordinal(aCellIndexInCluster);
            }
            else
            {
                MORIS_ERROR(false, "get_cell_side_ordinal(): can only be MASTER and SLAVE");
                return 0;
            }
        }

        //----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Double_Side_Cluster::get_master_integration_cell_side_ordinals() const
        {
            return this->get_master_side_cluster().get_cell_side_ordinals();
        }

        //----------------------------------------------------------------------------------

        moris_index
        Double_Side_Cluster::get_master_cell_side_ordinal(
                moris::moris_index aMasterCellIndexInCluster ) const
        {
            return this->get_master_side_cluster().get_cell_side_ordinal(aMasterCellIndexInCluster);
        }

        //----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Double_Side_Cluster::get_slave_integration_cell_side_ordinals() const
        {
            return this->get_slave_side_cluster().get_cell_side_ordinals();
        }


        //----------------------------------------------------------------------------------

        moris_index
        Double_Side_Cluster::get_slave_cell_side_ordinal(
                moris::moris_index aSlaveCellIndexInCluster ) const
        {
            return this->get_slave_side_cluster().get_cell_side_ordinal(aSlaveCellIndexInCluster);
        }

        //----------------------------------------------------------------------------------

        moris::Cell<moris::mtk::Vertex const *>
        Double_Side_Cluster::get_vertices_in_cluster(
                const mtk::Master_Slave aIsMaster ) const
        {
            if ( aIsMaster == mtk::Master_Slave::MASTER )
            {
                return this->get_master_side_cluster().get_vertices_in_cluster();
            }
            else if( aIsMaster == mtk::Master_Slave::SLAVE )
            {
                return this->get_slave_side_cluster().get_vertices_in_cluster();
            }
            else
            {
                MORIS_ERROR(false, "get_vertices_in_cluster(): can only be MASTER and SLAVE");
                
                return moris::Cell<moris::mtk::Vertex const *>(0);
            }
        }

        //----------------------------------------------------------------------------------

        moris::Cell<moris::mtk::Vertex const *>
        Double_Side_Cluster::get_master_vertices_in_cluster() const
        {
            return this->get_master_side_cluster().get_vertices_in_cluster();
        }

        //----------------------------------------------------------------------------------

        moris::Cell<moris::mtk::Vertex const *>
        Double_Side_Cluster::get_slave_vertices_in_cluster() const
        {
            return this->get_slave_side_cluster().get_vertices_in_cluster();
        }

        //----------------------------------------------------------------------------------

        moris_index
        Double_Side_Cluster::get_slave_vertex_ord_on_facet(
                moris_index                aCellClusterIndex,
                moris::mtk::Vertex const * aSlaveVertex       ) const
        {
            return mSlaveSideCluster->get_vertex_ordinal_on_facet(aCellClusterIndex,aSlaveVertex);
        }

        //----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Double_Side_Cluster::get_master_vertex_indices_in_cluster() const
        {
            return this->get_master_side_cluster().get_vertex_indices_in_cluster();
        }

        //----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Double_Side_Cluster::get_slave_vertex_indices_in_cluster() const
        {
            return this->get_slave_side_cluster().get_vertex_indices_in_cluster();
        }

        //----------------------------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Double_Side_Cluster::get_vertex_indices_in_cluster() const
        {

            // number of cells in cluster
            moris::uint tMasterNumVertices = mMasterSideCluster->get_num_vertices_in_cluster();
            moris::uint tSlaveNumVertices = mSlaveSideCluster->get_num_vertices_in_cluster();

            // access the vertices in a side cluster
            moris::Cell<moris::mtk::Vertex const *> tMasterVertices = mMasterSideCluster->get_vertices_in_cluster();
            moris::Cell<moris::mtk::Vertex const *> const & tSlaveVertices = mSlaveSideCluster->get_vertices_in_cluster();

            // initialize output
            moris::Matrix<moris::IndexMat> tVertexIndices(1,tMasterNumVertices+tSlaveNumVertices);

            // get cell indices and store
            for(moris::uint i = 0 ; i < tMasterNumVertices; i++)
            {
                tVertexIndices(i) = tMasterVertices(i)->get_index();
            }
            for(moris::uint i = 0 ; i < tSlaveNumVertices; i++)
            {
                tVertexIndices(tMasterNumVertices + i) = tSlaveVertices(i)->get_index();
            }

            return tVertexIndices;
        }


        //----------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        Double_Side_Cluster::get_vertices_local_coordinates_wrt_interp_cell( const mtk::Master_Slave aIsMaster ) const
        {
            if ( aIsMaster == mtk::Master_Slave::MASTER )
            {
                return this->get_master_side_cluster().get_vertices_local_coordinates_wrt_interp_cell();
            }
            else if( aIsMaster == mtk::Master_Slave::SLAVE )
            {
                return this->get_slave_side_cluster().get_vertices_local_coordinates_wrt_interp_cell();
            }
            else
            {
                MORIS_ERROR(false, "get_vertices_local_coordinates_wrt_interp_cell(): can only be MASTER and SLAVE");
                
                return {{}};
            }
        }


        //----------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        Double_Side_Cluster::get_master_vertices_local_coordinates_wrt_interp_cell() const
        {
            return this->get_master_side_cluster().get_vertices_local_coordinates_wrt_interp_cell();
        }

        //----------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        Double_Side_Cluster::get_slave_vertices_local_coordinates_wrt_interp_cell() const
        {
            return this->get_slave_side_cluster().get_vertices_local_coordinates_wrt_interp_cell();
        }


        //----------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        Double_Side_Cluster::get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex,
                const mtk::Master_Slave aIsMaster ) const
        {
            if ( aIsMaster == mtk::Master_Slave::MASTER )
            {
                return this->get_master_side_cluster().get_vertex_local_coordinate_wrt_interp_cell(aVertex);
            }
            else if( aIsMaster == mtk::Master_Slave::SLAVE )
            {
                return this->get_slave_side_cluster().get_vertex_local_coordinate_wrt_interp_cell(aVertex);
            }
            else
            {
                MORIS_ERROR(false, "get_vertex_local_coordinate_wrt_interp_cell(): can only be MASTER and SLAVE");
                return moris::Matrix<moris::DDRMat>(0,0);
            }
        }

        //----------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        Double_Side_Cluster::get_master_vertex_local_coordinate_wrt_interp_cell(
                moris::mtk::Vertex const * aVertex ) const
        {
            return this->get_master_side_cluster().get_vertex_local_coordinate_wrt_interp_cell(aVertex);
        }

        //----------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        Double_Side_Cluster::get_slave_vertex_local_coordinate_wrt_interp_cell(
                moris::mtk::Vertex const * aVertex ) const
        {
            return this->get_slave_side_cluster().get_vertex_local_coordinate_wrt_interp_cell(aVertex);
        }

        //----------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        Double_Side_Cluster::get_cell_local_coords_on_side_wrt_interp_cell(
                moris::moris_index      aClusterLocalIndex,
                const mtk::Master_Slave aIsMaster           ) const
        {
            if ( aIsMaster == mtk::Master_Slave::MASTER )
            {
                return this->get_master_side_cluster().get_cell_local_coords_on_side_wrt_interp_cell(aClusterLocalIndex);
            }
            else if( aIsMaster == mtk::Master_Slave::SLAVE )
            {
                return this->get_slave_side_cluster().get_cell_local_coords_on_side_wrt_interp_cell(aClusterLocalIndex);
            }
            else
            {
                MORIS_ERROR(false, "get_cell_local_coords_on_side_wrt_interp_cell(): can only be MASTER and SLAVE");
                return moris::Matrix<moris::DDRMat>(0,0);
            }
        }

        //----------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        Double_Side_Cluster::get_master_cell_local_coords_on_side_wrt_interp_cell(
                moris::moris_index aMasterClusterLocalIndex) const
        {
            return this->get_master_side_cluster().get_cell_local_coords_on_side_wrt_interp_cell(aMasterClusterLocalIndex);
        }

        //----------------------------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        Double_Side_Cluster::get_slave_cell_local_coords_on_side_wrt_interp_cell(
                moris::moris_index aSlaveClusterLocalIndex) const
        {
            return this->get_slave_side_cluster().get_cell_local_coords_on_side_wrt_interp_cell(aSlaveClusterLocalIndex);
        }

        //----------------------------------------------------------------------------------

        moris_index
        Double_Side_Cluster::get_dim_of_param_coord(
                const mtk::Master_Slave aIsMaster ) const
        {
            if ( aIsMaster == mtk::Master_Slave::MASTER )
            {
                return this->get_master_side_cluster().get_dim_of_param_coord();
            }
            else if( aIsMaster == mtk::Master_Slave::SLAVE )
            {
                return this->get_slave_side_cluster().get_dim_of_param_coord();
            }
            else
            {
                MORIS_ERROR(false, "get_dim_of_param_coord(): can only be MASTER and SLAVE");
                return 0;
            }
        }

        //----------------------------------------------------------------------------------

        moris::real
        Double_Side_Cluster::compute_cluster_cell_measure(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster       ) const
        {
            moris::mtk::Cluster const & tCluster = this->get_cluster(aIsMaster);

            return tCluster.compute_cluster_cell_measure(aPrimaryOrVoid,aIsMaster);
        }

        //----------------------------------------------------------------------------------

        Matrix<DDRMat>
        Double_Side_Cluster::compute_cluster_ig_cell_measures(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster       ) const
        {
            moris::mtk::Cluster const & tCluster = this->get_cluster(aIsMaster);

            return tCluster.compute_cluster_ig_cell_measures(aPrimaryOrVoid,aIsMaster);
        }

        //----------------------------------------------------------------------------------

        moris::real
        Double_Side_Cluster::compute_cluster_cell_measure_derivative(
                const Matrix< DDRMat > & aPerturbedVertexCoords,
                uint aDirection,
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster ) const
        {
            moris::mtk::Cluster const & tCluster = this->get_cluster(aIsMaster);
            return tCluster.compute_cluster_cell_measure_derivative(
                    aPerturbedVertexCoords,
                    aDirection,
                    aPrimaryOrVoid,
                    aIsMaster );
        }

        //----------------------------------------------------------------------------------

        moris::real
        Double_Side_Cluster::compute_cluster_cell_side_measure(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster       ) const
        {
            moris::mtk::Cluster const & tCluster = this->get_cluster(aIsMaster);

            return tCluster.compute_cluster_cell_side_measure(aPrimaryOrVoid,aIsMaster);
        }

        //----------------------------------------------------------------------------------

        Matrix<DDRMat>
        Double_Side_Cluster::compute_cluster_ig_cell_side_measures(
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster       ) const
        {
            moris::mtk::Cluster const & tCluster = this->get_cluster(aIsMaster);

            return tCluster.compute_cluster_ig_cell_side_measures(aPrimaryOrVoid,aIsMaster);
        }


        //----------------------------------------------------------------

        moris::real
        Double_Side_Cluster::compute_cluster_cell_side_measure_derivative(
                const Matrix< DDRMat > & aPerturbedVertexCoords,
                uint aDirection,
                const mtk::Primary_Void aPrimaryOrVoid,
                const mtk::Master_Slave aIsMaster ) const
        {
            moris::mtk::Cluster const & tCluster = this->get_cluster(aIsMaster);
            return tCluster.compute_cluster_cell_side_measure_derivative(
                    aPerturbedVertexCoords,
                    aDirection,
                    aPrimaryOrVoid,
                    aIsMaster );
        }

        //----------------------------------------------------------------------------------

        moris_index
        Double_Side_Cluster::get_master_dim_of_param_coord() const
        {
            return this->get_master_side_cluster().get_dim_of_param_coord();
        }

        //----------------------------------------------------------------------------------

        moris_index
        Double_Side_Cluster::get_slave_dim_of_param_coord() const
        {
            return this->get_slave_side_cluster().get_dim_of_param_coord();
        }


        //----------------------------------------------------------------------------------


        moris::uint
        Double_Side_Cluster::get_master_num_vertices_in_cluster() const
        {
            return this->get_master_side_cluster().get_num_vertices_in_cluster();
        }

        //----------------------------------------------------------------------------------

        moris::uint
        Double_Side_Cluster::get_slave_num_vertices_in_cluster() const
        {
            return this->get_slave_side_cluster().get_num_vertices_in_cluster();
        }

        //----------------------------------------------------------------------------------

        size_t
        Double_Side_Cluster::capacity()
        {
            size_t tCapacity = 0;

            //sum up the member data size
            tCapacity += sizeof( mMasterSideCluster );
            tCapacity += sizeof( mSlaveSideCluster );
            tCapacity += mMasterToSlaveVertexPairs.capacity() * ( ( sizeof( void * ) ) + 1 );

            return tCapacity;
        }
    }
}

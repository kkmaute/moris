#include <string>
#include <iomanip>
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Block.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Cell_Info.cpp"
#include "fn_norm.hpp"

#include "HDF5_Tools.hpp"

namespace moris
{
    namespace mtk
    {
        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::save_IG_node_TMatrices_to_file()
        {
            // get total number of integration vertices on IG mesh
            uint tNumVertices = this->get_num_nodes();
            MORIS_LOG_SPEC( "Num of IG vertices to save T-matrices for", tNumVertices );

            // initialize maps mapping (Background-BSp)-(Background-Lag) and (Background-Lag)-(Integration-Lag)
            Matrix< IdMat >  tBSToIPMap( this->get_num_basis_functions( 0 ), 1, gNoID );
            Matrix< IdMat >  tIPToIGMap( tNumVertices, 1, gNoID );

            // create/compute the above maps
            this->create_TMatrix_maps( tBSToIPMap, tIPToIGMap );
            
            // initialize cells containing info for T-matrices mapping (Background-BSp)-(Background-Lag)
            moris::Cell< Matrix< IdMat > >   tBSToIGIds( tNumVertices );
            moris::Cell< Matrix< DDRMat > >  tBSToIGWeights( tNumVertices );

            // build the B-spline BM - Lagrange BM relationship
            this->build_TMatrix_maps( tBSToIGIds, tBSToIGWeights, tBSToIPMap, tIPToIGMap );

            // initialize cells containing info for T-matrices mapping (Background-Lag)-(IG-Vertices)
            moris::Cell< Matrix< IdMat > >   tIPtoVertIds( tNumVertices );
            moris::Cell< Matrix< DDRMat > >  tIPtoVertWeights( tNumVertices );

            // go through all bulk sets
            for( uint iSet = 0; iSet < this->get_num_blocks(); iSet++ )
            {
                // get number of clusters on current set
                uint tNumClustersOnSet = this->get_set_by_index( iSet )->get_num_clusters_on_set();

                // proceed with saving, if set is non-empty
                if( tNumClustersOnSet > 0 )
                {
                    // go through clusters in current set
                    for( uint iCluster = 0; iCluster < tNumClustersOnSet; iCluster++ )
                    {
                        // get current cluster
                        const Cluster * tCluster = this->get_set_by_index( iSet )->get_clusters_by_index( iCluster );

                        // get info for IP cell associated with cluster currently handled
                        Cell_Info const * tCellInfo = tCluster->get_interpolation_cell().get_cell_info();

                        // creating interpolation rule
                        Interpolation_Rule tInterpolationRule(
                                tCellInfo->get_cell_geometry(),             // basic geometry of IP element (i.e. QUAD or HEX)
                                Interpolation_Type::LAGRANGE,               // space IP type
                                tCellInfo->get_cell_interpolation_order(),  // space IP order
                                Interpolation_Type::LAGRANGE,               // time IP type
                                Interpolation_Order::LINEAR);               // time IP order

                        // create a space interpolator
                        Space_Interpolator tSpaceInterpolator(tInterpolationRule);

                        // get list of list of coordinate of vertices on current cluster
                        moris::Matrix< moris::DDRMat > tLocalCoords = tCluster->get_vertices_local_coordinates_wrt_interp_cell();

                        // get list of (pointers to) vertices on current cluster
                        moris::Cell< moris::mtk::Vertex const * > tIGVertices = tCluster->get_vertices_in_cluster();

                        // get number of vertices on the treated mesh cluster
                        uint tNumVerticesOnCluster = tLocalCoords.n_rows();

                        // loop over the vertices on the treated mesh cluster
                        for( uint iVert = 0; iVert < tNumVerticesOnCluster; iVert++ )
                        {
                            // set interpolation point in space to current vertex
                            tSpaceInterpolator.set_space_time( trans( tLocalCoords.get_row( iVert ) ) );

                            // evaluate vector of shape functions at current vertex and retrieve it
                            const Matrix< DDRMat > & tN = tSpaceInterpolator.NXi();

                            // get current vertex ID and index
                            moris_id tIGVertexId = tIGVertices( iVert )->get_id();

                            // get number of shape functions interpolating into current vertex
                            uint tNumSFs = tN.n_cols();

                            // initialize size of T-Matrix for current vertex
                            tIPtoVertIds( tIGVertexId - 1 ).set_size( tNumSFs, 1, gNoID );
                            tIPtoVertWeights( tIGVertexId - 1 ).set_size( tNumSFs, 1, -1.0 );

                            // initialize counter
                            uint tCount = 0;

                            // loop over all T-Matrix entries
                            for( uint iSF = 0; iSF < tNumSFs; iSF++ )
                            {
                                // ignore T-matrix entries which are zero close to machine precision
                                if ( std::abs( tN( iSF ) ) > 10.0 * MORIS_REAL_EPS )
                                {
                                    // copy pointer of dof and convert to mtk::Vertex
                                    tIPtoVertIds( tIGVertexId - 1 )( tCount ) = tIGVertices( iSF )->get_id();

                                    // copy entry of T-Matrix
                                    tIPtoVertWeights( tIGVertexId - 1 )( tCount ) = tN( iSF );

                                    // increment counter
                                    tCount++;
                                }
                            }
                        } // end: loop over vertices on cluster
                    } // end: loop over clusters
                } // end: if cluster is empty
            } // end:: loop over bulk sets

            // initialize 
            Matrix< IdMat > tVertexIDs( 1, tNumVertices );
            Matrix< IdMat > tTmatrixIDs( tBSToIGWeights( 0 ).numel(), tNumVertices );
            Matrix< DDRMat > tTmatrixWeights( tBSToIGWeights( 0 ).numel(), tNumVertices );

            for( uint iVert = 0; iVert < tNumVertices; iVert++ )
            {
                tVertexIDs( iVert ) = iVert + 1;

                for( uint iTmat = 0; iTmat < tBSToIGWeights( 0 ).numel(); iTmat++ )
                {
                    tTmatrixIDs( iTmat, iVert )     = tBSToIGIds( iVert )( iTmat );
                    tTmatrixWeights( iTmat, iVert ) = tBSToIGWeights( iVert )( iTmat );
                }
            }

            // create file
            std::string tFilePath = "T_Matrices.hdf5";

            // make path parallel
            tFilePath = parallelize_path( tFilePath );

            // Create a new file using default properties
            hid_t tFileID = H5Fcreate( tFilePath.c_str(),
                    H5F_ACC_TRUNC,
                    H5P_DEFAULT,
                    H5P_DEFAULT);

            // error handler
            herr_t tStatus;

            // save ids to file
            save_matrix_to_hdf5_file(
                    tFileID,
                    "VertexIDs",
                    tVertexIDs,
                    tStatus );

            // save ids to file
            save_matrix_to_hdf5_file(
                    tFileID,
                    "IDs",
                    tTmatrixIDs,
                    tStatus );

            // save weights to file
            save_matrix_to_hdf5_file(
                    tFileID,
                    "Weights",
                    tTmatrixWeights,
                    tStatus );

            // close file
            tStatus = H5Fclose( tFileID );

            MORIS_ERROR( false, "Stop after writing T-Matrices." );

        }

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::create_TMatrix_maps(
                Matrix< IdMat > & tBSToIPMap,
                Matrix< IdMat > & tIPToIGMap)
        {
            // matrices that mark which ID pairs have been found
            Matrix< DDSMat > tIdentifierMatBS( tBSToIPMap.numel(),1 ,0 );
            Matrix< DDSMat > tIdentifierMatIP( tIPToIGMap.numel(),1 ,0 );

            // loop over all bulk sets
            for( uint iSet = 0; iSet < this->get_num_blocks(); iSet ++ )
            {
                // get number of clusters on current set
                uint tNumClustersOnSet = this->get_set_by_index( iSet )->get_num_clusters_on_set();

                // check if cluster is empty, if so skip procedure
                if( tNumClustersOnSet > 0 )
                {
                    // loop over all clusters in current bulk set
                    for( uint iCluster = 0; iCluster < tNumClustersOnSet; iCluster++ )
                    {
                        // get current cluster
                        const Cluster * tCluster = this->get_set_by_index( iSet )->get_clusters_by_index( iCluster );

                        // get IP cell current cluster is on
                        const moris::mtk::Cell & tInterpolationCell = tCluster->get_interpolation_cell();

                        // get list of all vertices on IP cell associated with current cluster
                        moris::Cell< Vertex *> tIPVertices = tInterpolationCell.get_vertex_pointers();

                        // get list of all vertices on IG cells in current cluster
                        moris::Cell< moris::mtk::Vertex const * > tIGVertices = tCluster->get_vertices_in_cluster();

                        // debug - ouput number of IG and IP vertices
                        //std::cout<<" IGVert: "<<tIGVertices.size()<<" IPVert: "<<tIPVertices.size()<<std::endl;

                        // go through all IP vertices, build BSp - Lag map
                        for( uint iIpVert = 0; iIpVert < tIPVertices.size(); iIpVert++)
                        {
                            // check if current IP vertex is empty
                            MORIS_ASSERT( tIPVertices( iIpVert )->has_interpolation( 0 ), "IP Vertex does not have a B-spline interpolation 0.");

                            // get IP vertex ID
                            uint tIP_ID = tIPVertices( iIpVert )->get_id();

                            // get list of B-spline IDs associated with current IP vertex
                            Matrix< IdMat > tBSpIDs = tIPVertices( iIpVert )->get_interpolation( 0 )->get_ids();

                            // build BSpline to IP map
                            if( tBSpIDs.numel() == 1 )
                            {
                                // get ID of first B-spline associated with IP vertex
                                moris_id tBSId = tBSpIDs( 0 );

                                // check if BSp-Lag pair has not already been found
                                if( tIdentifierMatBS( tBSId - 1 ) == 0 )
                                {
                                    // associating B-spline BF ID with Lagrange BF ID
                                    tBSToIPMap( tBSId - 1 ) = tIP_ID;

                                    // tick off that the BSp-Lag pair has been found
                                    tIdentifierMatBS( tBSId - 1 ) = 1;
                                }
                            }

                            // build IP to IG map (only fill IG vertices that fall on IP vertices)
                            // assume corner vertices come first

                            // get ID of Lag. BF associated with current IG vertex (IG vert. = current IP vert.)
                            moris_id tIG_ID = tIGVertices( iIpVert )->get_id();

                            // check if Lag-IGvert pair has not already been found
                            if( tIdentifierMatIP( tIP_ID - 1 ) == 0 )
                            {
                                // associating Lag ID with IG vertex ID
                                tIPToIGMap( tIP_ID - 1 ) = tIG_ID;

                                // tick off that the Lag-IGvert pair has been found
                                tIdentifierMatIP( tIP_ID - 1 ) = 1;
                            }
                        } // end: loop over IP vertices
                    } // end: loop over clusters
                } // skip: check if cluster is empty
            } // end:: loop over bulk sets
        }

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::build_TMatrix_maps(
                moris::Cell< Matrix< IdMat > > & tBSToIGIds,
                moris::Cell< Matrix< DDRMat > > & tBSToIGWeights,
                const Matrix< IdMat > & tBSToIPMap,
                const Matrix< IdMat > & tIPToIGMap)
        {
            // loop over all bulk sets
            for( uint iSet = 0; iSet < this->get_num_blocks(); iSet ++ )
            {
                // get number of clusters on set
                uint tNumClustersOnSet = this->get_set_by_index( iSet )->get_num_clusters_on_set();

                // check if cluster is empty, if so skip procedure
                if( tNumClustersOnSet > 0 )
                {
                    // loop over all clusters in set
                    for( uint iCluster = 0; iCluster < tNumClustersOnSet; iCluster++ )
                    {
                        // get current cluster
                        const Cluster * tCluster = this->get_set_by_index( iSet )->get_clusters_by_index( iCluster );

                        // IP cell associated with cluster
                        const moris::mtk::Cell & tInterpolationCell = tCluster->get_interpolation_cell();

                        // get list of vertices on IP cell of current cluster
                        moris::Cell< Vertex *> tIPVertices = tInterpolationCell.get_vertex_pointers(); 

                        // get list of IG vertices of IG cells on current cluster
                        moris::Cell< moris::mtk::Vertex const * > tIGVertices = tCluster->get_vertices_in_cluster();

                        // go through all IP vertices
                        for( uint iIpVert = 0; iIpVert < tIPVertices.size(); iIpVert++)
                        {
                            // build BSpline to iP map
                            if( tIPVertices( iIpVert )->get_interpolation( 0 )->get_ids().numel() > 1 )
                            {
                                // get ID of Lag. BF associated with current IG corner vertex (IG vert. = current IP vert.)
                                uint tIG_ID = tIGVertices( iIpVert )->get_id();

                                // get IDs and weights of B-splines associated with current IP vertex
                                Matrix< IdMat > tBSId = tIPVertices( iIpVert )->get_interpolation( 0 )->get_ids();
                                Matrix< DDRMat > tBSWeights= *tIPVertices( iIpVert )->get_interpolation( 0 )->get_weights();

                                // initialize size of ID and weight lists for T-Matrices
                                tBSToIGIds( tIG_ID - 1 ).resize( tBSId.numel(), 1 );
                                tBSToIGWeights( tIG_ID - 1 ).resize( tBSId.numel(), 1 );

                                // go through list of B-spline IDs associated with current IP vertex
                                for( uint iBSp = 0; iBSp < tBSId.numel(); iBSp++ )
                                {
                                    // copy weights and IDs onto T-matrices
                                    tBSToIGIds    ( tIG_ID - 1 )( iBSp ) = tIPToIGMap( tBSToIPMap( tBSId( iBSp ) - 1 ) );
                                    tBSToIGWeights( tIG_ID - 1 )( iBSp ) = tBSWeights( iBSp );
                                }
                            } // check: IP vertex has B-Spline interpolation
                        } // end: loop over IP vertices
                    } // end: loop over clusters
                } // check: cluster is non-empty
            } // end: loop over bulk sets
        }

        // ----------------------------------------------------------------------------

    }
}
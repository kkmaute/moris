/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Mesh_write_T_matrices.cpp
 *
 */

#include <string>
#include <iomanip>
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Block_Set.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Cell_Info.cpp"
#include "cl_Tracer.hpp"

#include "fn_norm.hpp"
#include "fn_stringify_matrix.hpp"
#include "fn_unique.hpp"
#include "fn_join_horiz.hpp"

#include "HDF5_Tools.hpp"

namespace moris
{
    namespace mtk
    {
        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::save_elemental_T_matrices_to_file(
                std::string aFileName,
                const uint  aNumBsplineMeshes )
        {
            // initialize arrays of T-matrix information
            uint                                           tNumIgCells = 0;
            Matrix< IdMat >                                tIgCellIds;
            Matrix< IndexMat >                             tIgCellIndices;
            moris::Vector< moris::Vector< Matrix< IdMat > > >  tIGtoBSIds;        // outer cell: B-spline mesh index | inner cell: IG cell index | matrix/vector: assembly map
            moris::Vector< moris::Vector< Matrix< DDRMat > > > tIGtoBSWeights;    // outer cell: B-spline mesh index | inner cell: IG cell index | matrix: extraction operator

            // -------------------------------------
            // get T-Matrices
            {
                // trace these operations
                Tracer tTracer( "MTK", "Compute elemental T-matrices" );

                // compute T-matrices from integration mesh to Lagrange mesh
                moris::Vector< moris_id >           tIgCellIdList;
                moris::Vector< Matrix< IndexMat > > tIgToIpIndices;
                moris::Vector< Matrix< DDRMat > >   tIgToIpTmatrices;
                this->get_IG_to_IP_elemental_T_matrices( tIgCellIdList, tIgToIpIndices, tIgToIpTmatrices );

                // compute T-matrices from Lagrange mesh to the different B-spline meshes
                moris::Vector< moris::Vector< Matrix< IdMat > > >  tIPtoBSIds;
                moris::Vector< moris::Vector< Matrix< DDRMat > > > tIPtoBSWeights;
                this->get_IP_to_BS_nodal_T_matrices( tIPtoBSIds, tIPtoBSWeights, aNumBsplineMeshes );

                // combine the above information to obtain the T-matrices from the IG mesh to the B-spline meshes
                this->get_elemental_IG_to_BS_T_matrices( tIgToIpIndices, tIgToIpTmatrices, tIPtoBSIds, tIPtoBSWeights, tIGtoBSIds, tIGtoBSWeights );

                // count all used cell IDs in the output mesh
                for ( const moris_id iIgCellId : tIgCellIdList )
                {
                    // do not count unused indices
                    if ( iIgCellId > -1 )
                    {
                        tNumIgCells++;
                    }
                }

                // convert the index to ID map to a matrix list of all used IDs
                tIgCellIds.set_size( tNumIgCells, 1, MORIS_ID_MAX );
                tIgCellIndices.set_size( tNumIgCells, 1, MORIS_INDEX_MAX );
                uint tIdCounter    = 0;
                uint tIndexCounter = 0;
                for ( const moris_id iIgCellId : tIgCellIdList )
                {
                    // only list used cells
                    if ( iIgCellId > -1 )
                    {
                        tIgCellIds( tIdCounter )     = iIgCellId;
                        tIgCellIndices( tIdCounter ) = tIndexCounter;
                        tIdCounter++;
                    }

                    // count indices for all cells though
                    tIndexCounter++;
                }

            }    // end: operations to get extraction operators (curly brackets to delete arrays which are no longer needed to free up memory)

            // -------------------------------------
            // write to file

            // trace these operations
            Tracer tTracer( "MTK", "Save elemental T-matrices" );

            // get the number of B-spline meshes
            uint tNumBsplineMeshes = tIGtoBSIds.size();

            // save extraction operators for all B-spline meshes
            for ( uint iBspMesh = 0; iBspMesh < tNumBsplineMeshes; iBspMesh++ )
            {
                // construct string from output file name
                std::string tStrOutputFile = aFileName + "_B" + std::to_string( iBspMesh ) + ".hdf5";

                // log/print that the extraction operator is output to
                MORIS_LOG_INFO( "Saving elemental extraction operators for B-spline mesh #%i to file: %s", iBspMesh, tStrOutputFile.c_str() );

                // initialize hdf5 file
                hid_t  tFileID = create_hdf5_file( tStrOutputFile );
                herr_t tStatus = 0;

                // save the index to ID map to the .hdf5
                save_matrix_to_hdf5_file( tFileID, "Cell_IDs", tIgCellIds, tStatus );

                // write elemental information for all IG cells
                for ( uint iIgCell = 0; iIgCell < tNumIgCells; iIgCell++ )
                {
                    // get the IG cell's index
                    moris_index tIgCellIndex = tIgCellIndices( iIgCell );

                    // get the current IG cell's ID
                    moris_id tID = tIgCellIds( iIgCell );

                    // assemble names
                    std::string tIDsName     = "IDs_" + std::to_string( tID );
                    std::string tWeightsName = "Weights_" + std::to_string( tID );

                    // write the IDs
                    save_matrix_to_hdf5_file( tFileID, tIDsName, tIGtoBSIds( iBspMesh )( tIgCellIndex ), tStatus );
                    MORIS_ERROR( tStatus == 0,
                            "Integration_Mesh::save_elemental_T_matrices_to_file() - "
                            "HDF5 writer returned status %i writing the basis IDs for IG cell %i.",
                            tStatus,
                            tID );

                    // write the IDs
                    save_matrix_to_hdf5_file( tFileID, tWeightsName, tIGtoBSWeights( iBspMesh )( tIgCellIndex ), tStatus );
                    MORIS_ERROR( tStatus == 0,
                            "Integration_Mesh::save_elemental_T_matrices_to_file() - "
                            "HDF5 writer returned status %i writing the weights for IG cell %i.",
                            tStatus,
                            tID );
                }

                // close the hdf5 file
                close_hdf5_file( tFileID );
            }
        }

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::save_IG_global_T_matrix_to_file( std::string aFileName )
        {
            // trace this function
            Tracer tTracer( "MTK", "Save Nodal T-Matrices to File" );

            // get total number of integration vertices on IG mesh
            uint tNumVertices = this->get_num_nodes();
            MORIS_LOG_SPEC( "Num of IG vertices to save T-matrices for", tNumVertices );

            // get Number of Block Sets
            uint tNumSets = this->get_num_blocks();

            // go through all bulk sets
            for ( uint iSet = 0; iSet < tNumSets; iSet++ )
            {
                // -------------------------------------
                // get T-Matrices

                // initialize cells containing info for T-matrices mapping (Background-BSp)-(Background-Lag) for each IP vertex
                moris::Vector< Matrix< IdMat > >  tIPtoBSIds;
                moris::Vector< Matrix< DDRMat > > tIPtoBSWeights;

                // get T-matrices mapping from IP nodes to B-splines
                this->get_IP_to_BS_nodal_T_matrices( tIPtoBSIds, tIPtoBSWeights, iSet );

                // initialize cells containing info for T-matrices mapping (Background-Lag)-(IG-Vertices) for each IG vertex
                moris::Vector< Matrix< IdMat > >  tIGtoIPIds;
                moris::Vector< Matrix< DDRMat > > tIGtoIPWeights;

                // get T-matrices mapping from IG nodes to IP nodes
                this->get_IG_to_IP_nodal_T_matrices( tIGtoIPIds, tIGtoIPWeights, iSet );

                // initialize cells containing info for T-matrices mapping (Background-BSp)-(IG-Vertices) for each IG vertex
                moris::Vector< Matrix< IdMat > >  tIGtoBSIds;
                moris::Vector< Matrix< DDRMat > > tIGtoBSWeights;

                // combine T-Matrices to get
                this->get_IG_to_BS_nodal_T_matrices( tIGtoBSIds, tIGtoBSWeights, tIPtoBSIds, tIPtoBSWeights, tIGtoIPIds, tIGtoIPWeights );

                // -------------------------------------
                // combine everything into one matrix

                // initialize matrices for sparse mat
                Matrix< DDUMat > tSparseIndices;
                Matrix< DDRMat > tWeights;

                // combine everything
                this->build_sparse_extraction_operator( tIGtoBSIds, tIGtoBSWeights, tSparseIndices, tWeights );

                // -------------------------------------
                // write to file

                // construct string from output file name
                std::string tStrOutputFile = aFileName + "." + ios::stringify( iSet ) + ".hdf5";

                // log/print that the extraction operator is output to
                MORIS_LOG_INFO( "Saving Extraction Operator for Set-# %i to file: %s", iSet, tStrOutputFile.c_str() );

                // initialize hdf5 file
                hid_t  tFileID = create_hdf5_file( tStrOutputFile );
                herr_t tStatus = 0;

                // write to file
                save_matrix_to_hdf5_file( tFileID, "IDs", tSparseIndices, tStatus );
                save_matrix_to_hdf5_file( tFileID, "Weights", tWeights, tStatus );

                // close file
                close_hdf5_file( tFileID );

            }    // end: loop over all sets
        }

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::get_IG_to_IP_elemental_T_matrices(
                moris::Vector< moris_id >           &aIgCellIds,
                moris::Vector< Matrix< IndexMat > > &aIgToIpIndices,
                moris::Vector< Matrix< DDRMat > >   &aIgToIpTmatrices )
        {
            // trace this function
            Tracer tTracer( "MTK", "Compute elemental IG to IP T-Matrices" );

            // get the number of IG elements in the mesh
            uint tNumIgCells = this->get_num_elems();

            // initialize the size of the lists
            aIgCellIds.resize( tNumIgCells, gNoID );
            aIgToIpIndices.resize( tNumIgCells );
            aIgToIpTmatrices.resize( tNumIgCells );

            // get the number of sets in the mesh
            uint tNumSetsInMesh = this->get_num_blocks();

            // loop over all bulk clusters in the mesh by first looping over the sets and then clusters within them
            for ( uint iSet = 0; iSet < tNumSetsInMesh; iSet++ )
            {
                // get access to the current set
                const mtk::Set *tSet              = this->get_set_by_index( iSet );
                uint            tNumClustersOnSet = tSet->get_num_clusters_on_set();

                // skip this set if it doesn't contain any clusters
                if ( tNumClustersOnSet == 0 )
                {
                    continue;
                }

                // go over the individual clusters
                for ( uint iClusterInSet = 0; iClusterInSet < tNumClustersOnSet; iClusterInSet++ )
                {
                    // get the currently treated cluster
                    const Cluster *tCluster = tSet->get_clusters_by_index( iClusterInSet );

                    // get IP cell (UIPC) the current cluster is constructed from
                    const moris::mtk::Cell &tInterpolationCell = tCluster->get_interpolation_cell();

                    // get info for IP cell associated with cluster currently handled
                    Cell_Info const *tCellInfo = tCluster->get_interpolation_cell().get_cell_info();

                    // create an interpolation rule
                    Interpolation_Rule tInterpolationRule(
                            tCellInfo->get_cell_geometry(),               // basic geometry of IP element (i.e. QUAD or HEX)
                            Interpolation_Type::LAGRANGE,                 // space IP type
                            tCellInfo->get_cell_interpolation_order(),    // space IP order
                            Interpolation_Type::LAGRANGE,                 // time IP type
                            Interpolation_Order::LINEAR );                // time IP order

                    // create a space interpolator
                    Space_Interpolator tSpaceInterpolator( tInterpolationRule );

                    // get list of indices and coordinates of vertices on current cluster
                    const Matrix< IndexMat > tVertexIndsInCluster      = tCluster->get_vertex_indices_in_cluster();
                    const Matrix< DDRMat >   tVertexLocalCoords        = tCluster->get_vertices_local_coordinates_wrt_interp_cell();
                    uint                     tNumPrimaryVertsInCluster = tVertexIndsInCluster.numel();

                    // build map stating where in the list of all IG vertices on the current cluster the primary IG vertices sit
                    std::unordered_map< moris_index, uint > tIgVertexMap;
                    for ( uint iIgVert = 0; iIgVert < tNumPrimaryVertsInCluster; iIgVert++ )
                    {
                        tIgVertexMap[ tVertexIndsInCluster( iIgVert ) ] = iIgVert;
                    }

                    // initialize list of nodal T-matrices
                    // input: IG vertex index local to current cluster || output: nodal T-matrix (weights)
                    moris::Vector< Matrix< DDRMat > > tNodalTmatWeights( tNumPrimaryVertsInCluster );

                    // loop over the vertices on the treated mesh cluster
                    for ( uint iPrimaryIgVert = 0; iPrimaryIgVert < tNumPrimaryVertsInCluster; iPrimaryIgVert++ )
                    {
                        // set interpolation point in space to current vertex
                        tSpaceInterpolator.set_space_time( trans( tVertexLocalCoords.get_row( iPrimaryIgVert ) ) );

                        // evaluate vector of shape functions at current vertex and retrieve it
                        const Matrix< DDRMat > &tN = tSpaceInterpolator.NXi();

                        // copy into the weights
                        tNodalTmatWeights( iPrimaryIgVert ) = tN;

                    }    // end for: vertices on cluster

                    // get list of (unzipped ?) vertex indices on this cell, corresponding to the basis functions evaluated
                    Matrix< IndexMat > tUnzippedIpVertIndices       = tInterpolationCell.get_vertex_inds();
                    uint               tNumBfsInterpolatingIntoCell = tUnzippedIpVertIndices.numel();

                    // get a list of primary IG cells
                    const moris::Vector< const moris::mtk::Cell * > &tCellsInCluster  = tCluster->get_primary_cells_in_cluster();
                    uint                                           tNumPrimaryCells = tCellsInCluster.size();

                    // go over the primary elements
                    for ( uint iPrimaryCell = 0; iPrimaryCell < tNumPrimaryCells; iPrimaryCell++ )
                    {
                        // quick access to this IG cell
                        const moris::mtk::Cell *tIgCell = tCellsInCluster( iPrimaryCell );

                        // get this cell's index and ID
                        moris_index tCellIndex = tIgCell->get_index();
                        moris_id    tCellId    = tIgCell->get_id();

                        // store ID in index to ID map
                        MORIS_ERROR( aIgCellIds( tCellIndex ) == gNoID,
                                "Integration_Mesh::get_IG_element_to_IP_T_matrices() - "
                                "IG cell (#%i) already has a T-matrix assigned to it. It should not be overwritten.",
                                tCellIndex );
                        aIgCellIds( tCellIndex ) = tCellId;

                        // store the IP vertex indices in the T-matrices
                        aIgToIpIndices( tCellIndex ) = tUnzippedIpVertIndices;

                        // get the vertices on this element
                        Matrix< IndexMat > tVertIndsOnIgCell = tIgCell->get_vertex_inds();
                        uint               tNumIgVertsOnCell = tVertIndsOnIgCell.numel();

                        // get the vertex indices in this cell
                        aIgToIpTmatrices( tCellIndex ).set_size( tNumIgVertsOnCell, tNumBfsInterpolatingIntoCell );

                        // go over
                        for ( uint iVertOnIgCell = 0; iVertOnIgCell < tNumIgVertsOnCell; iVertOnIgCell++ )
                        {
                            // get the index of the vertex treated
                            moris_index tVertexIndex = tVertIndsOnIgCell( iVertOnIgCell );

                            // get the vertex index inside the cluster
                            auto tIter                 = tIgVertexMap.find( tVertexIndex );
                            uint tVertexIndexInCluster = tIter->second;
                            MORIS_ASSERT( tIter != tIgVertexMap.end(),
                                    "Integration_Mesh::get_IG_element_to_IP_T_matrices() - "
                                    "Vertex index not found in list of IG vertices on the current IG cell." );

                            Matrix< DDRMat > &tWeightsOnCurrentIgVert     = tNodalTmatWeights( tVertexIndexInCluster );
                            Matrix< DDRMat > &tIgToIpWeightsOnCurrentCell = aIgToIpTmatrices( tCellIndex );

                            // add the T-matrix weights of the current IG vertex weights into the elemental T-matrix
                            tIgToIpWeightsOnCurrentCell(
                                    { iVertOnIgCell, iVertOnIgCell },
                                    { 0, tNumBfsInterpolatingIntoCell - 1 } ) =
                                    tWeightsOnCurrentIgVert.matrix_data();

                        }    // end for: vertices in primary IG cell
                    }        // end for: primary IG cells in cluster
                }            // end for: clusters in set
            }                // end for: sets in mesh
        }                    // Integration_Mesh::get_IG_element_to_IP_T_matrices()

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::get_IG_to_IP_nodal_T_matrices(
                moris::Vector< Matrix< IdMat > >  &aIGtoIPIds,
                moris::Vector< Matrix< DDRMat > > &aIGtoIPWeights,
                uint                             aSetIndex )
        {
            // trace this function
            Tracer tTracer( "MTK", "Compute nodal IG to IP T-Matrices" );

            // get number of IG vertices on mesh
            uint tNumVertices = this->get_num_nodes();

            // initialize cells containing info for T-matrices mapping (Background-Lag)-(IG-Vertices) for each IG vertex
            aIGtoIPIds.resize( tNumVertices );
            aIGtoIPWeights.resize( tNumVertices );

            // get number of clusters on current set
            uint tNumClustersOnSet = this->get_set_by_index( aSetIndex )->get_num_clusters_on_set();

            // proceed with saving, if set is non-empty
            if ( tNumClustersOnSet > 0 )
            {
                // go through clusters in current set
                for ( uint iCluster = 0; iCluster < tNumClustersOnSet; iCluster++ )
                {
                    // get current cluster
                    const Cluster *tCluster = this->get_set_by_index( aSetIndex )->get_clusters_by_index( iCluster );

                    // get IP cell current cluster is on
                    const moris::mtk::Cell &tInterpolationCell = tCluster->get_interpolation_cell();

                    // get list of all vertices on IP cell associated with current cluster
                    moris::Vector< Vertex * > tIPVertices = tInterpolationCell.get_vertex_pointers();

                    // get info for IP cell associated with cluster currently handled
                    Cell_Info const *tCellInfo = tCluster->get_interpolation_cell().get_cell_info();

                    // creating interpolation rule
                    Interpolation_Rule tInterpolationRule(
                            tCellInfo->get_cell_geometry(),               // basic geometry of IP element (i.e. QUAD or HEX)
                            Interpolation_Type::LAGRANGE,                 // space IP type
                            tCellInfo->get_cell_interpolation_order(),    // space IP order
                            Interpolation_Type::LAGRANGE,                 // time IP type
                            Interpolation_Order::LINEAR );                // time IP order

                    // create a space interpolator
                    Space_Interpolator tSpaceInterpolator( tInterpolationRule );

                    // get list of list of coordinate of vertices on current cluster
                    moris::Matrix< moris::DDRMat > tLocalCoords = tCluster->get_vertices_local_coordinates_wrt_interp_cell();

                    // get list of (pointers to) all vertices on current cluster
                    moris::Vector< moris::mtk::Vertex const * > tIGVertices = tCluster->get_vertices_in_cluster();

                    // build map stating where in the list of all IG vertices on the current cluster the primary IG vertices sit
                    std::unordered_map< moris_id, uint > tIgVertexMap;
                    for ( uint iIgVert = 0; iIgVert < tIGVertices.size(); iIgVert++ )
                    {
                        tIgVertexMap[ tIGVertices( iIgVert )->get_id() ] = iIgVert;
                    }

                    // get list of (pointers to) primary vertices on current cluster
                    moris::Vector< moris::mtk::Vertex * > tPrimaryIGVertices = tCluster->get_primary_vertices_in_cluster();

                    // get number of primary vertices on the treated mesh cluster
                    uint tNumPrimaryVerticesOnCluster = tPrimaryIGVertices.size();

                    // loop over the vertices on the treated mesh cluster
                    for ( uint iPrimaryIgVert = 0; iPrimaryIgVert < tNumPrimaryVerticesOnCluster; iPrimaryIgVert++ )
                    {
                        // get the current primary IG vertex ID
                        moris_id tPrimaryVertexID = tPrimaryIGVertices( iPrimaryIgVert )->get_id();

                        // look for position of primary vertex in list of all vertices on cluster
                        auto tIter = tIgVertexMap.find( tPrimaryVertexID );

                        // check that the primary vertex is actually in the list of
                        MORIS_ERROR( tIter != tIgVertexMap.end(),
                                "Integration_Mesh::get_IG_to_IP_nodal_T_matrices() - Primary vertex ID on cluster not found in list of all vertex IDs on cluster" );

                        // get the position of the current primary IG vertex in the list of all IG vertices on the cluster
                        uint tIgVertListIndex = tIter->second;

                        // set interpolation point in space to current vertex
                        tSpaceInterpolator.set_space_time( trans( tLocalCoords.get_row( tIgVertListIndex ) ) );

                        // evaluate vector of shape functions at current vertex and retrieve it
                        const Matrix< DDRMat > &tN = tSpaceInterpolator.NXi();

                        // get number of shape functions interpolating into current vertex
                        uint tNumSFs = tN.n_cols();

                        // initialize size of T-Matrix for current vertex
                        aIGtoIPIds( tPrimaryVertexID - 1 ).set_size( 1, tNumSFs, gNoID );
                        aIGtoIPWeights( tPrimaryVertexID - 1 ).set_size( 1, tNumSFs, -1.0 );

                        // initialize counter
                        uint tCount = 0;

                        // loop over all T-Matrix entries
                        for ( uint iSF = 0; iSF < tNumSFs; iSF++ )
                        {
                            // ignore T-matrix entries which are zero close to machine precision
                            if ( std::abs( tN( iSF ) ) > 10.0 * MORIS_REAL_EPS )
                            {
                                // copy pointer of dof and convert to mtk::Vertex
                                aIGtoIPIds( tPrimaryVertexID - 1 )( tCount ) = tIPVertices( iSF )->get_id();

                                // copy entry of T-Matrix
                                aIGtoIPWeights( tPrimaryVertexID - 1 )( tCount ) = tN( iSF );

                                // increment counter
                                tCount++;
                            }
                        }
                    }    // end: loop over vertices on cluster
                }        // end: loop over clusters
            }            // end: if cluster is empty
        }

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::get_IP_to_BS_nodal_T_matrices(
                moris::Vector< moris::Vector< Matrix< IdMat > > >  &aIPtoBSIds,        // input: B-spline mesh index, UIPV index || output: B-spline BF IDs interpolating into UIPV
                moris::Vector< moris::Vector< Matrix< DDRMat > > > &aIPtoBSWeights,    // input: B-spline mesh index, UIPV index || output: T-matrix weights
                const uint                                      aNumBsplineMeshes )
        {
            // trace this function
            Tracer tTracer( "MTK", "Compute IP to B-Spline T-Matrices" );

            // get number of IP nodes on mesh
            uint tNumLagBfs = this->get_max_IP_index_in_mesh() + 1;

            // initialize storage for T-matrices
            aIPtoBSIds.resize( aNumBsplineMeshes );
            aIPtoBSWeights.resize( aNumBsplineMeshes );
            for ( uint iBspMesh = 0; iBspMesh < aNumBsplineMeshes; iBspMesh++ )
            {
                aIPtoBSIds( iBspMesh ).resize( tNumLagBfs );
                aIPtoBSWeights( iBspMesh ).resize( tNumLagBfs );
            }

            // punch card marking which UIPVs have already been treated
            moris::Vector< bool > tUipvIndFound( tNumLagBfs, false );

            // get the number of sets in the mesh
            uint tNumSetsInMesh = this->get_num_blocks();

            // loop over all bulk clusters in the mesh by first looping over the sets and then clusters within them
            for ( uint iSet = 0; iSet < tNumSetsInMesh; iSet++ )
            {
                // get access to the current set
                mtk::Set const *tSet = this->get_set_by_index( iSet );

                // get number of clusters on current set
                uint tNumClustersOnSet = tSet->get_num_clusters_on_set();

                // check if cluster is empty, if so skip procedure
                if ( tNumClustersOnSet == 0 )
                {
                    continue;
                }

                // loop over all clusters in current bulk set
                for ( uint iClusterInSet = 0; iClusterInSet < tNumClustersOnSet; iClusterInSet++ )
                {
                    // get current cluster
                    const Cluster *tCluster = this->get_set_by_index( iSet )->get_clusters_by_index( iClusterInSet );

                    // get IP cell current cluster is on
                    const moris::mtk::Cell &tInterpolationCell = tCluster->get_interpolation_cell();

                    // get list of all vertices on IP cell associated with current cluster
                    moris::Vector< Vertex * > tIPVertices          = tInterpolationCell.get_vertex_pointers();
                    uint                    tNumIpVertsOnCluster = tIPVertices.size();

                    // go through all IP vertices, build BSp - Lag map
                    for ( uint iIpVert = 0; iIpVert < tNumIpVertsOnCluster; iIpVert++ )
                    {
                        // FIXME: this check needs to go back in when we have a method for checking whether a vertex has an associated interpolation
                        // check if current IP vertex is empty
                        // MORIS_ASSERT( tIPVertices( iIpVert )->has_interpolation( 0 ), "IP Vertex does not have a B-spline interpolation 0.");

                        // get unzipped IP vertex index
                        uint tUipvIndex = tIPVertices( iIpVert )->get_index();

                        // if this Basis function already has its T-matrices assigned, skip it
                        if ( tUipvIndFound( tUipvIndex ) )
                        {
                            continue;
                        }

                        // mark this UIPV as found
                        tUipvIndFound( tUipvIndex ) = true;

                        //
                        for ( uint iBspMesh = 0; iBspMesh < aNumBsplineMeshes; iBspMesh++ )
                        {
                            // get list of B-spline IDs associated with current IP vertex
                            Matrix< IdMat >  tBSpIDs    = tIPVertices( iIpVert )->get_interpolation( iBspMesh )->get_ids();
                            Matrix< DDRMat > tBSWeights = *tIPVertices( iIpVert )->get_interpolation( iBspMesh )->get_weights();

                            // copy the weights and basis IDs onto the output T-matrices
                            aIPtoBSIds( iBspMesh )( tUipvIndex )     = trans( tBSpIDs );
                            aIPtoBSWeights( iBspMesh )( tUipvIndex ) = trans( tBSWeights );

                        }    // end for: B-spline meshes for which the T-matrix should be retrieved
                    }        // end for: unzipped IP vertices on cluster
                }            // end for: clusters in set
            }                // end for: sets in mesh
        }                    // end function: Integration_Mesh::get_IP_to_BS_nodal_T_matrices()

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::get_IP_to_BS_nodal_T_matrices(
                moris::Vector< Matrix< IdMat > >  &aIPtoBSIds,
                moris::Vector< Matrix< DDRMat > > &aIPtoBSWeights,
                uint                             aSetIndex )
        {
            // trace this function
            Tracer tTracer( "MTK", "Compute IP to B-Spline T-Matrices" );

            // get number of IP nodes on mesh
            uint tNumIpNodes = this->get_max_IP_ID_on_set( aSetIndex );

            // initialize map size
            aIPtoBSIds.resize( tNumIpNodes );
            aIPtoBSWeights.resize( tNumIpNodes );

            // matrices that mark which ID pairs have been found
            Matrix< DDSMat > tIdentifierMatIP( tNumIpNodes, 1, 0 );

            // get number of clusters on current set
            uint tNumClustersOnSet = this->get_set_by_index( aSetIndex )->get_num_clusters_on_set();

            // check if cluster is empty, if so skip procedure
            if ( tNumClustersOnSet > 0 )
            {
                // loop over all clusters in current bulk set
                for ( uint iCluster = 0; iCluster < tNumClustersOnSet; iCluster++ )
                {
                    // get current cluster
                    const Cluster *tCluster = this->get_set_by_index( aSetIndex )->get_clusters_by_index( iCluster );

                    // get IP cell current cluster is on
                    const moris::mtk::Cell &tInterpolationCell = tCluster->get_interpolation_cell();

                    // get list of all vertices on IP cell associated with current cluster
                    moris::Vector< Vertex * > tIPVertices = tInterpolationCell.get_vertex_pointers();

                    // go through all IP vertices, build BSp - Lag map
                    for ( uint iIpVert = 0; iIpVert < tIPVertices.size(); iIpVert++ )
                    {
                        // FIXME: this check needs to go back in when we have a method for checking whether a vertex has an associated interpolation
                        // check if current IP vertex is empty
                        // MORIS_ASSERT( tIPVertices( iIpVert )->has_interpolation( 0 ), "IP Vertex does not have a B-spline interpolation 0.");

                        // get IP vertex ID
                        uint tIP_ID = tIPVertices( iIpVert )->get_id();

                        MORIS_ERROR( tIP_ID > 0, "Integration_Mesh::get_IP_to_BS_nodal_T_matrices() - ID expected to be greater equal 1." );

                        // check if IP vertex has already been treated
                        if ( tIdentifierMatIP( tIP_ID - 1 ) == 0 )
                        {
                            // get list of B-spline IDs associated with current IP vertex
                            Matrix< IdMat >  tBSpIDs    = tIPVertices( iIpVert )->get_interpolation( 0 )->get_ids();
                            Matrix< DDRMat > tBSWeights = *tIPVertices( iIpVert )->get_interpolation( 0 )->get_weights();

                            // get number of B-splines associated with current IP vertex
                            uint tNumBspOnIpNode = tBSpIDs.numel();

                            // resize nodal T-matrix to accommodate all B-spline IDs and weights
                            aIPtoBSIds( tIP_ID - 1 ).set_size( 1, tNumBspOnIpNode, gNoID );
                            aIPtoBSWeights( tIP_ID - 1 ).set_size( 1, tNumBspOnIpNode, -1 );

                            // create T-matrices - copy list of associated B-spline IDs and their weights
                            for ( uint iBsp = 0; iBsp < tNumBspOnIpNode; iBsp++ )
                            {
                                aIPtoBSIds( tIP_ID - 1 )( iBsp )     = tBSpIDs( iBsp );
                                aIPtoBSWeights( tIP_ID - 1 )( iBsp ) = tBSWeights( iBsp );
                            }

                            // tick off IP vertex as treated
                            tIdentifierMatIP( tIP_ID - 1 ) = 1;
                        }
                    }    // end for: IP vertices
                }        // end for: clusters
            }            // end if: cluster is non-empty
        }

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::get_elemental_IG_to_BS_T_matrices(
                moris::Vector< Matrix< IndexMat > > const              &aIgToIpIndices,      // outer cell: IG cell index
                moris::Vector< Matrix< DDRMat > > const                &aIgToIpTmatrices,    // outer cell: IG cell index
                moris::Vector< moris::Vector< Matrix< IdMat > > > const  &aIPtoBSIds,          // outer cell: B-spline mesh index | inner cell: enr. Lagrange BF index
                moris::Vector< moris::Vector< Matrix< DDRMat > > > const &aIPtoBSWeights,      // outer cell: B-spline mesh index | inner cell: enr. Lagrange BF index
                moris::Vector< moris::Vector< Matrix< IdMat > > >        &aIGtoBSIds,          // outer cell: B-spline mesh index | inner cell: IG cell index
                moris::Vector< moris::Vector< Matrix< DDRMat > > >       &aIGtoBSWeights )           // outer cell: B-spline mesh index | inner cell: IG cell index
        {
            // trace this function
            Tracer tTracer( "MTK", "Compute elemental IG to IP T-Matrices" );

            // get the number of B-spline meshes
            uint tNumBspMeshes = aIPtoBSIds.size();

            // resize the output arrays
            aIGtoBSIds.resize( tNumBspMeshes );
            aIGtoBSWeights.resize( tNumBspMeshes );

            // compute extraction operators for all B-spline meshes
            for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
            {
                // get the number of IG cells
                uint tNumIgCells = aIgToIpIndices.size();

                // resize output arrays
                aIGtoBSIds( iBspMesh ).resize( tNumIgCells );
                aIGtoBSWeights( iBspMesh ).resize( tNumIgCells );

                // build T-matrices for each IG cell
                for ( uint iIgCell = 0; iIgCell < tNumIgCells; iIgCell++ )
                {
                    // initialize map that relates a B-spline basis ID to its position in the element local extraction operator
                    std::map< moris_id, uint > tBsplineIdToLocalIndexMap;

                    // initialize counter that tracks how many B-spline basis functions interpolate into the current IG cell
                    uint tNumBspBfsInElem = 0;

                    // get the number of Lagrange basis functions that interpolate into the current IG element
                    uint tNumLagBfs = aIgToIpIndices( iIgCell ).numel();

                    // go over the individual vertices on the IG cell and collect the related B-spline basis IDs
                    for ( uint iLagBf = 0; iLagBf < tNumLagBfs; iLagBf++ )
                    {
                        // get the current Lagrange BF's index
                        moris_index tLagBfIndex = aIgToIpIndices( iIgCell )( iLagBf );

                        // get the B-spline basis functions interpolating into this Lagrange basis function
                        Matrix< IdMat > const &tBspIdsInLagBf    = aIPtoBSIds( iBspMesh )( tLagBfIndex );
                        uint                   tNumBspBfsInLagBf = tBspIdsInLagBf.numel();

                        for ( uint iBspBf = 0; iBspBf < tNumBspBfsInLagBf; iBspBf++ )
                        {
                            // get the ID of the Bspline basis function
                            moris_id tBspBfId = tBspIdsInLagBf( iBspBf );

                            // check if this B-spline BF's ID is already in the map and add it if not
                            auto tIter = tBsplineIdToLocalIndexMap.find( tBspBfId );
                            if ( tIter == tBsplineIdToLocalIndexMap.end() )
                            {
                                tBsplineIdToLocalIndexMap[ tBspBfId ] = tNumBspBfsInElem;
                                tNumBspBfsInElem++;
                            }
                        }    // end for: B-spline basis functions interpolating into Lagrange basis function
                    }        // end for: Lagrange basis functions interpolating into IG element

                    // get the number of nodes on the current IG cell
                    uint tNumNodesOnIgElem = aIgToIpTmatrices( iIgCell ).n_rows();

                    // get access to where the current IG element's T-matrix will be stored
                    Matrix< IdMat >  &tCurrentIgCellToBSIds     = aIGtoBSIds( iBspMesh )( iIgCell );
                    Matrix< DDRMat > &tCurrentIgCellToBSWeights = aIGtoBSWeights( iBspMesh )( iIgCell );

                    // correctly size the T-matrix information
                    tCurrentIgCellToBSIds.set_size( 1, tNumBspBfsInElem );
                    tCurrentIgCellToBSWeights.set_size( tNumNodesOnIgElem, tNumBspBfsInElem, 0.0 );

                    // loop over map entries and dump them into the T-matrix ID information
                    for ( auto const &iBspBF : tBsplineIdToLocalIndexMap )
                    {
                        tCurrentIgCellToBSIds( iBspBF.second ) = iBspBF.first;
                    }

                    // go over the weights for each Lagrange basis and multiply them by the weights for the B-spline basis functions
                    for ( uint iLagBf = 0; iLagBf < tNumLagBfs; iLagBf++ )
                    {
                        // get the current Lagrange BF's index
                        moris_index tLagBfIndex = aIgToIpIndices( iIgCell )( iLagBf );

                        // get the B-spline basis function IDs interpolating into this Lagrange basis function
                        Matrix< IdMat > const &tBspIdsInLagBf    = aIPtoBSIds( iBspMesh )( tLagBfIndex );
                        uint                   tNumBspBfsInLagBf = tBspIdsInLagBf.numel();

                        // get the IP T-matrix entries for this Lagrange BF
                        auto tIpWeights = aIgToIpTmatrices( iIgCell )( { 0, tNumNodesOnIgElem - 1 }, { iLagBf, iLagBf } );

                        // get the weights for this B-spline basis
                        Matrix< DDRMat > const &tBspWeightsInLagBf = aIPtoBSWeights( iBspMesh )( tLagBfIndex );

                        // multiply weights
                        Matrix< DDRMat > tBspWeights = tIpWeights * tBspWeightsInLagBf;

                        // put the B-spline weights for every B-spline ID back into the elemental T-matrix
                        for ( uint iBspBf = 0; iBspBf < tNumBspBfsInLagBf; iBspBf++ )
                        {
                            // get the ID of the Bspline basis function
                            moris_id tBspBfId = tBspIdsInLagBf( iBspBf );

                            // get the local index of this B-spline basis ID
                            auto tIter = tBsplineIdToLocalIndexMap.find( tBspBfId );
                            MORIS_ERROR( tIter != tBsplineIdToLocalIndexMap.end(),
                                    "Integration_Mesh::get_IG_to_BS_nodal_T_matrices() - "
                                    "B-spline ID not found in IG cell local map." );
                            uint tLocalIndex = tIter->second;

                            // grab the part of the elemental T-matrix which the weights should be written into
                            auto tCopyWeightsIntoThis = tCurrentIgCellToBSWeights( { 0, tNumNodesOnIgElem - 1 }, { tLocalIndex, tLocalIndex } );
                            auto tCopyWeightsFromThis = tBspWeights( { 0, tNumNodesOnIgElem - 1 }, { iBspBf, iBspBf } );

                            // copy the weights into this
                            tCopyWeightsIntoThis += tCopyWeightsFromThis;

                        }    // end for: B-spline basis functions interpolating into Lagrange basis function
                    }        // end for: Lagrange basis functions interpolating into IG element
                }            // end for: IG elements
            }                // end for: B-spline meshes
        }                    // Integration_Mesh::get_IG_to_BS_nodal_T_matrices()

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::get_IG_to_BS_nodal_T_matrices(
                moris::Vector< Matrix< IdMat > >  &aIGtoBSIds,
                moris::Vector< Matrix< DDRMat > > &aIGtoBSWeights,
                moris::Vector< Matrix< IdMat > >  &aIPtoBSIds,
                moris::Vector< Matrix< DDRMat > > &aIPtoBSWeights,
                moris::Vector< Matrix< IdMat > >  &aIGtoIPIds,
                moris::Vector< Matrix< DDRMat > > &aIGtoIPWeights )
        {
            // trace this function
            Tracer tTracer( "MTK", "Compute IG to IP T-Matrices" );

            // get number of IG vertices
            uint tNumIgNodes = aIGtoIPIds.size();
            MORIS_ERROR( tNumIgNodes == aIGtoIPWeights.size(),
                    "Integration_Mesh::get_IG_to_BS_nodal_T_matrices() - Size of ID/weight maps for IG-IP T-Matrices don't match." );

            // initialize size of T-Matrix cells
            aIGtoBSIds.resize( tNumIgNodes );
            aIGtoBSWeights.resize( tNumIgNodes );

            // loop over all IG vertices
            for ( uint iIgNode = 0; iIgNode < tNumIgNodes; iIgNode++ )
            {
                // initialize temporary vectors containing all
                Matrix< IdMat >  tNodalIGtoBSIds( 0, 0 );
                Matrix< DDRMat > tNodalIGtoBSWeights( 0, 0 );

                // number of IP nodes associated with current IG node
                uint tNumIpBfs = aIGtoIPIds( iIgNode ).numel();

                // loop over all IP nodes current IG node is related to
                for ( uint jIpNode = 0; jIpNode < tNumIpBfs; jIpNode++ )
                {
                    // get current IP ID
                    moris_id tIpId = aIGtoIPIds( iIgNode )( jIpNode );

                    // skip if IP-ID doesn't exist
                    if ( tIpId > gNoID )
                    {
                        // get what B-spline IDs and Weights are associated with current Node
                        Matrix< IdMat >  tAddIdList     = aIPtoBSIds( tIpId - 1 );
                        Matrix< DDRMat > tAddWeightList = aIGtoIPWeights( iIgNode )( jIpNode ) * aIPtoBSWeights( tIpId - 1 );

                        // append vectors of IDs and weights
                        tNodalIGtoBSIds     = join_horiz( tNodalIGtoBSIds, tAddIdList );
                        tNodalIGtoBSWeights = join_horiz( tNodalIGtoBSWeights, tAddWeightList );
                    }
                }    // end: loop over all IP Nodes

                // get unique list of B-Spline IDs associated with each IG node
                moris::unique( tNodalIGtoBSIds, aIGtoBSIds( iIgNode ) );

                // get number of unique IDs in list
                uint tNumUniqueIDs = aIGtoBSIds( iIgNode ).numel();

                // initialize list of weights associated with unique list of IDs
                aIGtoBSWeights( iIgNode ).set_size( tNumUniqueIDs, 1, 0.0 );

                // create list of weights associated with unique list of IDs
                for ( uint iUnique = 0; iUnique < tNumUniqueIDs; iUnique++ )
                {
                    for ( uint jNonUnique = 0; jNonUnique < tNodalIGtoBSIds.numel(); jNonUnique++ )
                    {
                        if ( aIGtoBSIds( iIgNode )( iUnique ) == tNodalIGtoBSIds( jNonUnique ) )
                        {
                            aIGtoBSWeights( iIgNode )( iUnique ) = aIGtoBSWeights( iIgNode )( iUnique ) + tNodalIGtoBSWeights( jNonUnique );
                        }
                    }
                }
            }    // end: loop over all IG Nodes
        }

        // ----------------------------------------------------------------------------

        uint
        Integration_Mesh::get_max_IP_ID_on_set( uint aSetIndex )
        {
            // initialize counter
            uint tMaxID = 0;

            // get number of clusters on current set
            uint tNumClustersOnSet = this->get_set_by_index( aSetIndex )->get_num_clusters_on_set();

            // check if cluster is empty, if so skip procedure
            if ( tNumClustersOnSet > 0 )
            {
                // loop over all clusters in current bulk set
                for ( uint iCluster = 0; iCluster < tNumClustersOnSet; iCluster++ )
                {
                    // get current cluster
                    const Cluster *tCluster = this->get_set_by_index( aSetIndex )->get_clusters_by_index( iCluster );

                    // get IP cell current cluster is on
                    const moris::mtk::Cell &tInterpolationCell = tCluster->get_interpolation_cell();

                    // get list of all vertices on IP cell associated with current cluster
                    moris::Vector< Vertex * > tIPVertices = tInterpolationCell.get_vertex_pointers();

                    // go through all IP vertices, build BSp - Lag map
                    for ( uint iIpVert = 0; iIpVert < tIPVertices.size(); iIpVert++ )
                    {
                        // get IP vertex ID
                        uint tIP_ID = tIPVertices( iIpVert )->get_id();

                        tMaxID = std::max( tMaxID, tIP_ID );
                    }
                }
            }

            // return max ID
            return tMaxID;
        }

        // ----------------------------------------------------------------------------

        uint
        Integration_Mesh::get_max_IP_index_in_mesh()
        {
            // initialize counter
            moris_index tMaxIndex = 0;

            // get the number of sets in the mesh
            uint tNumSetsInMesh = this->get_num_blocks();

            // loop over all bulk clusters in the mesh by first looping over the sets and then clusters within them
            for ( uint iSet = 0; iSet < tNumSetsInMesh; iSet++ )
            {
                // get access to the current set
                mtk::Set const *tSet = this->get_set_by_index( iSet );

                // get number of clusters on current set
                uint tNumClustersOnSet = tSet->get_num_clusters_on_set();

                // check if cluster is empty, if so skip procedure
                if ( tNumClustersOnSet == 0 )
                {
                    continue;
                }

                // loop over all clusters in current bulk set
                for ( uint iCluster = 0; iCluster < tNumClustersOnSet; iCluster++ )
                {
                    // get current cluster
                    const Cluster *tCluster = tSet->get_clusters_by_index( iCluster );

                    // get IP cell current cluster is on
                    const moris::mtk::Cell &tInterpolationCell = tCluster->get_interpolation_cell();

                    // get list of all vertices on IP cell associated with current cluster
                    Matrix< IndexMat > tIpVertexIndices = tInterpolationCell.get_vertex_inds();

                    // get and store the maximum index
                    moris_index tMaxIpVertIndOnIpElem = tIpVertexIndices.max();
                    tMaxIndex                         = std::max( tMaxIndex, tMaxIpVertIndOnIpElem );
                }
            }

            // return max ID
            return tMaxIndex;
        }

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh::build_sparse_extraction_operator(
                moris::Vector< Matrix< IdMat > >  &aIGtoBSIds,
                moris::Vector< Matrix< DDRMat > > &aIGtoBSWeights,
                Matrix< DDUMat >                &aSparseIndices,
                Matrix< DDRMat >                &aWeights )
        {
            // trace this function
            Tracer tTracer( "MTK", "Build Sparse Extraction Operator Matrix" );

            // get number of IG nodes
            uint tNumIgNodes = aIGtoBSIds.size();

            // initialize length of global sparse operator
            uint tNumIdWeightPairs = 0;

            // loop over all IG Nodes to figure out total number of weights
            for ( uint iG = 0; iG < tNumIgNodes; iG++ )
            {
                tNumIdWeightPairs += aIGtoBSIds( iG ).numel();
            }

            // initialize sparse matrix
            aSparseIndices.set_size( tNumIdWeightPairs, 2 );
            aWeights.set_size( tNumIdWeightPairs, 1 );

            // initialize index counter
            uint tIndex = 0;

            // loop over all IG Nodes
            for ( uint iG = 0; iG < tNumIgNodes; iG++ )
            {
                // check if IG node has interpolation
                if ( aIGtoBSIds( iG ).numel() > 0 )
                {
                    // loop over B-splines interpolating into current IG node
                    for ( uint iBsp = 0; iBsp < aIGtoBSIds( iG ).numel(); iBsp++ )
                    {
                        // get
                        uint tBspId = aIGtoBSIds( iG )( iBsp );

                        // write IG/BS indices and weights to list
                        aSparseIndices( tIndex, 0 ) = iG + 1;
                        aSparseIndices( tIndex, 1 ) = tBspId;
                        aWeights( tIndex )          = aIGtoBSWeights( iG )( iBsp );

                        // increment index
                        tIndex++;
                    }
                }
            }    // end: loop over all IG nodes
        }

        // ----------------------------------------------------------------------------
    }    // namespace mtk
}    // namespace moris

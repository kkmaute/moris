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
#include <numeric>
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
#include "cl_MTK_Cell_Info.hpp"
#include "cl_Tracer.hpp"

#include "fn_norm.hpp"
#include "fn_stringify_matrix.hpp"
#include "fn_unique.hpp"
#include "fn_join_horiz.hpp"
#include "fn_trans.hpp"
#include "fn_isempty.hpp"
#include "cl_Communication_Tools.hpp"

#include "HDF5_Tools.hpp"

namespace moris::mtk
{
    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::save_elemental_T_matrices_to_file(
            const std::string &aFileName,
            const uint         aNumBsplineMeshes )
    {
        // initialize arrays of T-matrix information
        uint                                 tNumIgCells = 0;
        Matrix< IdMat >                      tIgCellIds;
        Matrix< IndexMat >                   tIgCellIndices;
        Vector< Vector< Matrix< IdMat > > >  tIGtoBSIds;        // outer cell: B-spline mesh index | inner cell: IG cell index | matrix/vector: assembly map
        Vector< Vector< Matrix< DDRMat > > > tIGtoBSWeights;    // outer cell: B-spline mesh index | inner cell: IG cell index | matrix: extraction operator

        // -------------------------------------
        // get T-Matrices
        {
            // trace these operations
            Tracer tTracer( "MTK", "Compute elemental T-matrices" );

            // compute T-matrices from integration mesh to Lagrange mesh
            Vector< moris_id >           tIgCellIdList;
            Vector< Matrix< IndexMat > > tIgToIpIndices;
            Vector< Matrix< DDRMat > >   tIgToIpTmatrices;
            this->get_IG_to_IP_elemental_T_matrices( tIgCellIdList, tIgToIpIndices, tIgToIpTmatrices );

            // compute T-matrices from Lagrange mesh to the different B-spline meshes
            Vector< Vector< Matrix< IdMat > > >  tIPtoBSIds;
            Vector< Vector< Matrix< DDRMat > > > tIPtoBSWeights;
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
    Integration_Mesh::save_IG_global_T_matrix_to_file( const std::string &aFileName )
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

            // Only write the t-matrices to file for set 0
            if (this->get_set_by_index(iSet)->get_set_name() != "HMR_dummy_c_p0")
            {
                continue;
            }

            // initialize cells containing info for T-matrices mapping (Background-BSp)-(Background-Lag) for each IP vertex
            Vector< Matrix< IdMat > >  tIPtoBSIds;
            Vector< Matrix< DDRMat > > tIPtoBSWeights;
           
            // get T-matrices mapping from IP nodes to B-splines
            this->get_IP_to_BS_nodal_T_matrices( tIPtoBSIds, tIPtoBSWeights, iSet );

            // initialize cells containing info for T-matrices mapping (Background-Lag)-(IG-Vertices) for each IG vertex
            Vector< Matrix< IdMat > >  tIGtoIPIds;
            Vector< Matrix< IdMat > >  tIGtoIPIdsX;
            Vector< Matrix< IdMat > >  tIGtoIPIdsY;
            Vector< Matrix< IdMat > >  tIGtoIPIdsZ;

            Vector< Matrix< DDRMat > > tIGtoIPWeights;
            Vector< Matrix< DDRMat > > tIGtoIPGradientWeightsX;
            Vector< Matrix< DDRMat > > tIGtoIPGradientWeightsY;
            Vector< Matrix< DDRMat > > tIGtoIPGradientWeightsZ;

            std::vector< int >         tIGNodeIDs;
            std::unordered_map< moris_index, moris_id > tIGNodeIDMap;

            // get T-matrices mapping from IG nodes to IP nodes
            this->get_IG_to_IP_nodal_T_matrices( tIGtoIPIds, tIGtoIPIdsX, tIGtoIPIdsY, tIGtoIPIdsZ, tIGtoIPWeights, tIGtoIPGradientWeightsX, tIGtoIPGradientWeightsY, tIGtoIPGradientWeightsZ, tIGNodeIDs, tIGNodeIDMap, iSet );

            // initialize cells containing info for T-matrices mapping (Background-BSp)-(IG-Vertices) for each IG vertex
            Vector< Matrix< IdMat > >  tIGtoBSIds;
            Vector< Matrix< IdMat > >  tIGtoBSIdsX;
            Vector< Matrix< IdMat > >  tIGtoBSIdsY;
            Vector< Matrix< IdMat > >  tIGtoBSIdsZ;
            
            Vector< Matrix< DDRMat > > tIGtoBSWeights;
            Vector< Matrix< DDRMat > > tIGtoBSGradientWeightsX;
            Vector< Matrix< DDRMat > > tIGtoBSGradientWeightsY;
            Vector< Matrix< DDRMat > > tIGtoBSGradientWeightsZ;


            // combine T-Matrices to get
            this->get_IG_to_BS_nodal_T_matrices( tIGtoBSIds, tIGtoBSIdsX, tIGtoBSIdsY, tIGtoBSIdsZ, tIGtoBSWeights, tIGtoBSGradientWeightsX, tIGtoBSGradientWeightsY, tIGtoBSGradientWeightsZ, tIPtoBSIds, tIPtoBSWeights, tIGtoIPIds, tIGtoIPIdsX, tIGtoIPIdsY, 
                  tIGtoIPIdsZ, tIGtoIPWeights, tIGtoIPGradientWeightsX, tIGtoIPGradientWeightsY, tIGtoIPGradientWeightsZ );

            // -------------------------------------
            // combine everything into one matrix

            // initialize matrices for sparse mat
            Matrix< DDUMat > tSparseIndices;
            Matrix< DDUMat > tSparseIndicesX;
            Matrix< DDUMat > tSparseIndicesY;
            Matrix< DDUMat > tSparseIndicesZ;
            
            Matrix< DDRMat > tWeights;
            Matrix< DDRMat > tGradientWeightsX;
            Matrix< DDRMat > tGradientWeightsY;
            Matrix< DDRMat > tGradientWeightsZ;
            
            

            // combine everything
            this->build_sparse_extraction_operator( tIGtoBSIds, tIGtoBSIdsX, tIGtoBSIdsY, tIGtoBSIdsZ, tIGtoBSWeights, tIGtoBSGradientWeightsX, tIGtoBSGradientWeightsY, tIGtoBSGradientWeightsZ, tSparseIndices, tSparseIndicesX, tSparseIndicesY, tSparseIndicesZ, tWeights, tGradientWeightsX, tGradientWeightsY, tGradientWeightsZ, tIGNodeIDMap );

            // -------------------------------------
            // write to file

            // construct string from output file name
            std::string tStrOutputFile = aFileName + "." + ios::stringify( iSet ) + ".hdf5";
            std::string tStrOutputFileGradientX = aFileName + "Gradient" + "X" + "." + ios::stringify( iSet ) + ".hdf5";
            std::string tStrOutputFileGradientY = aFileName + "Gradient" + "Y" + "." + ios::stringify( iSet ) + ".hdf5";
            std::string tStrOutputFileGradientZ = aFileName + "Gradient" + "Z" + "." + ios::stringify( iSet ) + ".hdf5";
            std::string tStrOutputFileNodeIDs = aFileName + "_IG_Node_IDs" + "." + ios::stringify( iSet ) + ".hdf5";


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

            hid_t  tFileID1 = create_hdf5_file( tStrOutputFileGradientX );
            herr_t tStatus1 = 0;

            // write to file
            save_matrix_to_hdf5_file( tFileID1, "IDs", tSparseIndicesX, tStatus1 );
            save_matrix_to_hdf5_file( tFileID1, "Weights", tGradientWeightsX, tStatus1 );

            // close file
            close_hdf5_file( tFileID1 );

            hid_t  tFileID2 = create_hdf5_file( tStrOutputFileGradientY );
            herr_t tStatus2 = 0;

            // write to file
            save_matrix_to_hdf5_file( tFileID2, "IDs", tSparseIndicesY, tStatus2 );
            save_matrix_to_hdf5_file( tFileID2, "Weights", tGradientWeightsY, tStatus2 );

            // close file
            close_hdf5_file( tFileID2 );

            // initialize hdf5 file - This time for the IG Node IDs
            hid_t  tFileID3 = create_hdf5_file( tStrOutputFileNodeIDs );
            herr_t tStatus3 = 0;

            // write to file
            save_vector_to_hdf5_file( tFileID3, "IDs", tIGNodeIDs, tStatus3 );

            // close file
            close_hdf5_file( tFileID3 );

            if (this->get_spatial_dim() == 3)
            {
                hid_t  tFileID4 = create_hdf5_file( tStrOutputFileGradientZ );
                herr_t tStatus4 = 0;

                // write to file
                save_matrix_to_hdf5_file( tFileID4, "IDs", tSparseIndicesZ, tStatus4 );
                save_matrix_to_hdf5_file( tFileID4, "Weights", tGradientWeightsZ, tStatus4 );

                // close file
                close_hdf5_file( tFileID4 );
            }

        }    // end: loop over all sets
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::get_IG_to_IP_elemental_T_matrices(
            Vector< moris_id >           &aIgCellIds,
            Vector< Matrix< IndexMat > > &aIgToIpIndices,
            Vector< Matrix< DDRMat > >   &aIgToIpTmatrices )
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
                Vector< Matrix< DDRMat > > tNodalTmatWeights( tNumPrimaryVertsInCluster );

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
                const Vector< const moris::mtk::Cell * > &tCellsInCluster  = tCluster->get_primary_cells_in_cluster();
                uint                                      tNumPrimaryCells = tCellsInCluster.size();

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
                        MORIS_ERROR( tIter != tIgVertexMap.end(),
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
                }    // end for: primary IG cells in cluster
            }    // end for: clusters in set
        }    // end for: sets in mesh
    }    // Integration_Mesh::get_IG_element_to_IP_T_matrices()

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::get_IG_to_IP_nodal_T_matrices(
            Vector< Matrix< IdMat > >  &aIGtoIPIds,
            Vector< Matrix< IdMat > >  &aIGtoIPIdsX,
            Vector< Matrix< IdMat > >  &aIGtoIPIdsY,
            Vector< Matrix< IdMat > >  &aIGtoIPIdsZ,
            Vector< Matrix< DDRMat > > &aIGtoIPWeights,
            Vector< Matrix< DDRMat > > &aIGtoIPGradientWeightsX,
            Vector< Matrix< DDRMat > > &aIGtoIPGradientWeightsY,
            Vector< Matrix< DDRMat > > &aIGtoIPGradientWeightsZ,
            std::vector< int >             &aIGNodeIDs,
            std::unordered_map< moris_index, moris_id > &aIGNodeIDMap,
            uint                        aSetIndex )
    {
        // trace this function
        Tracer tTracer( "MTK", "Compute nodal IG to IP T-Matrices" );

        // get number of IG vertices on mesh
        uint tNumVertices = this->get_num_nodes();

        // initialize cells containing info for T-matrices mapping (Background-Lag)-(IG-Vertices) for each IG vertex
        aIGtoIPIds.resize( tNumVertices );
        aIGtoIPIdsX.resize( tNumVertices );
        aIGtoIPIdsY.resize( tNumVertices );
        aIGtoIPWeights.resize( tNumVertices );
        aIGtoIPGradientWeightsX.resize( tNumVertices );
        aIGtoIPGradientWeightsY.resize( tNumVertices );
        if ( this->get_spatial_dim() == 3 )
        {
            aIGtoIPIdsZ.resize( tNumVertices );
            aIGtoIPGradientWeightsZ.resize( tNumVertices );
        }
        

        // To ensure node IDs are not repeatedly added into aIGNodeIDs
        std::unordered_map< moris_id, uint > tIgVertexCount;

        // initialize counter for number of IG vertices
        uint tIgCount = 0;

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
                Vector< Vertex * > tIPVertices = tInterpolationCell.get_vertex_pointers();

                // get base IP cell 
                const moris::mtk::Cell* tBaseCell = tInterpolationCell.get_base_cell();

                // Get vertices associated with base IP cell
                Vector< Vertex * > tBaseIPVertices = tBaseCell->get_vertex_pointers();
                
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

                // Get physical coordinates from base cell
                Matrix< DDRMat > tPhysicalIPCoords = tBaseCell->get_vertex_coords();

                // Set physical coordinates
                tSpaceInterpolator.set_space_coeff( tPhysicalIPCoords );

                // get list of list of coordinate of vertices on current cluster
                moris::Matrix< moris::DDRMat > tLocalCoords = tCluster->get_vertices_local_coordinates_wrt_interp_cell();

                // get list of (pointers to) all vertices on current cluster
                Vector< moris::mtk::Vertex const * > tIGVertices = tCluster->get_vertices_in_cluster();

                // build map stating where in the list of all IG vertices on the current cluster the primary IG vertices sit
                std::unordered_map< moris_id, uint > tIgVertexMap;
                for ( uint iIgVert = 0; iIgVert < tIGVertices.size(); iIgVert++ )
                {
                    tIgVertexMap[ tIGVertices( iIgVert )->get_id() ] = iIgVert;
                }

                // get list of (pointers to) primary vertices on current cluster
                Vector< moris::mtk::Vertex * > tPrimaryIGVertices = tCluster->get_primary_vertices_in_cluster();

                // get number of primary vertices on the treated mesh cluster
                uint tNumPrimaryVerticesOnCluster = tPrimaryIGVertices.size();

                // Get coordinates of base IP cell for checking if primary IG cell lies on a cell edge.
                //Matrix< DDRMat > tVert0 = tBaseIPVertices( 0 )->get_coords();
                //Matrix< DDRMat > tVert1 = tBaseIPVertices( 1 )->get_coords();
                //Matrix< DDRMat > tVert2 = tBaseIPVertices( 2 )->get_coords();

                // loop over the vertices on the treated mesh cluster
                for ( uint iPrimaryIgVert = 0; iPrimaryIgVert < tNumPrimaryVerticesOnCluster; iPrimaryIgVert++ )
                {                    
                    // get the current primary IG vertex ID
                    moris_id tPrimaryVertexID = tPrimaryIGVertices( iPrimaryIgVert )->get_id();   

                    // get the current primary IG vertex index
                    moris_index tPrimaryVertexIndex = tPrimaryIGVertices( iPrimaryIgVert )->get_index();

                    // Get ownership rank of current IG vertex
                    int tOwningProc = tPrimaryIGVertices( iPrimaryIgVert )->get_owner();

                    // Only proceed if current processor owns the vertex
                    if ( tOwningProc != par_rank() )
                    {
                        continue;
                    }

                    // Check if it has already been looped over
                    if ( tIgVertexCount.count( tPrimaryVertexID ) == 0 )
                    {
                        // If not looped over, add to the list of IG node IDs, increment the count and set the vertex as treated
                        
                        aIGNodeIDs.push_back( tPrimaryVertexID );
                        tIgVertexCount[ tPrimaryVertexID ] = 1;
                        aIGNodeIDMap[ tPrimaryVertexIndex ] = tPrimaryVertexID;
                        tIgCount++;
                    
                        
                    }

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

                    // evaluate vector of gradient of shape functions at the current vertex and retrieve it.
                    const Matrix< DDRMat >& tdNdXi = tSpaceInterpolator.dNdXi();

                    // Get the geometric Jacobian for the gradients
                    // evaluate the space Jacobian from the geometry interpolator
                    const Matrix< DDRMat >& tInvJGeot = tSpaceInterpolator.inverse_space_jacobian();

                    // compute first derivative of the space shape function wrt x
                    Matrix< DDRMat > tdNdX = tInvJGeot * tdNdXi;

                    // get number of shape functions interpolating into current vertex
                    uint tNumSFs = tN.n_cols();

                    // initialize size of T-Matrix for current vertex
                    // aIGtoIPIds( tPrimaryVertexID - 1 ).set_size( 1, tNumSFs, gNoID );
                    // aIGtoIPIdsX( tPrimaryVertexID - 1 ).set_size( 1, tNumSFs, gNoID );
                    // aIGtoIPIdsY( tPrimaryVertexID - 1 ).set_size( 1, tNumSFs, gNoID );
                    // if ( this->get_spatial_dim() == 3 )
                    // {
                    //     aIGtoIPIdsZ( tPrimaryVertexID - 1 ).set_size( 1, tNumSFs, gNoID );
                    // }

                    aIGtoIPIds( tPrimaryVertexIndex ).set_size( 1, tNumSFs, gNoID );
                    aIGtoIPIdsX( tPrimaryVertexIndex ).set_size( 1, tNumSFs, gNoID );
                    aIGtoIPIdsY( tPrimaryVertexIndex ).set_size( 1, tNumSFs, gNoID );
                    if ( this->get_spatial_dim() == 3 )
                    {
                        aIGtoIPIdsZ( tPrimaryVertexIndex ).set_size( 1, tNumSFs, gNoID );
                    }

                    // aIGtoIPWeights( tPrimaryVertexID - 1 ).set_size( 1, tNumSFs, -1.0 );
                    // aIGtoIPGradientWeightsX( tPrimaryVertexID - 1 ).set_size( 1, tNumSFs, -1.0 );
                    // aIGtoIPGradientWeightsY( tPrimaryVertexID - 1 ).set_size( 1, tNumSFs, -1.0 );
                    // if ( this->get_spatial_dim() == 3 )
                    // {
                    //     aIGtoIPGradientWeightsZ( tPrimaryVertexID - 1 ).set_size( 1, tNumSFs, -1.0 );
                    // }

                    aIGtoIPWeights( tPrimaryVertexIndex ).set_size( 1, tNumSFs, -1.0 );
                    aIGtoIPGradientWeightsX( tPrimaryVertexIndex ).set_size( 1, tNumSFs, -1.0 );
                    aIGtoIPGradientWeightsY( tPrimaryVertexIndex ).set_size( 1, tNumSFs, -1.0 );
                    if ( this->get_spatial_dim() == 3 )
                    {
                        aIGtoIPGradientWeightsZ( tPrimaryVertexIndex ).set_size( 1, tNumSFs, -1.0 );
                    }
                    // initialize counter
                    uint tCount = 0;

                    // initialize counter x IDs
                    uint tCountX = 0;

                    // initialize counter y IDs
                    uint tCountY = 0;

                    // initialize counter z IDs
                    uint tCountZ = 0;

                    // loop over all T-Matrix entries
                    for ( uint iSF = 0; iSF < tNumSFs; iSF++ )
                    {
                        
                        // ignore T-matrix entries which are zero close to machine precision
                        if ( std::abs( tN( iSF ) ) > 10.0 * MORIS_REAL_EPS )
                        {
                        
                            // copy pointer of dof and convert to mtk::Vertex
                            //aIGtoIPIds( tPrimaryVertexID - 1 )( tCount ) = tIPVertices( iSF )->get_id();
                            aIGtoIPIds( tPrimaryVertexIndex )( tCount ) = tIPVertices( iSF )->get_id();

                            // copy entry of T-Matrix
                            //aIGtoIPWeights( tPrimaryVertexID - 1 )( tCount ) = tN( iSF );
                            aIGtoIPWeights( tPrimaryVertexIndex )( tCount ) = tN( iSF );

                            // increment counter
                            tCount++;
                        }
                        // copy entry of gradient
                        if ( std::abs( tdNdX( 0 , iSF ) ) > 10.0 * MORIS_REAL_EPS )
                        {
                            // copy pointer of dof and convert to mtk::Vertex
                            //aIGtoIPIdsX( tPrimaryVertexID - 1 )( tCountX ) = tIPVertices( iSF )->get_id();
                            aIGtoIPIdsX( tPrimaryVertexIndex )( tCountX ) = tIPVertices( iSF )->get_id();

                            //aIGtoIPGradientWeightsX( tPrimaryVertexID - 1 )( tCountX ) = tdNdX( 0 , iSF );
                            aIGtoIPGradientWeightsX( tPrimaryVertexIndex )( tCountX ) = tdNdX( 0, iSF );

                            // increment counter
                            tCountX++;
                            
                        }
                        if ( std::abs( tdNdX( 1 , iSF ) ) > 10.0 * MORIS_REAL_EPS )
                        {
                            // copy pointer of dof and convert to mtk::Vertex
                            //aIGtoIPIdsY( tPrimaryVertexID - 1 )( tCountY ) = tIPVertices( iSF )->get_id();
                            aIGtoIPIdsY( tPrimaryVertexIndex )( tCountY ) = tIPVertices( iSF )->get_id();

                            //aIGtoIPGradientWeightsY( tPrimaryVertexID - 1 )( tCountY ) = tdNdX( 1 , iSF );
                            aIGtoIPGradientWeightsY( tPrimaryVertexIndex )( tCountY ) = tdNdX( 1, iSF ); 

                            // increment counter
                            tCountY++;
                            
                        } 
                        if ( this->get_spatial_dim() == 3 )
                        {
                            if ( std::abs( tdNdX( 2 , iSF ) ) > 10.0 * MORIS_REAL_EPS )
                            {
                                // copy pointer of dof and convert to mtk::Vertex
                                //aIGtoIPIdsZ( tPrimaryVertexID - 1 )( tCountZ ) = tIPVertices( iSF )->get_id();
                                aIGtoIPIdsZ( tPrimaryVertexIndex )( tCountZ ) = tIPVertices( iSF )->get_id();
                                
                                //aIGtoIPGradientWeightsZ( tPrimaryVertexID - 1 )( tCountZ ) = tdNdX( 2 , iSF );
                                aIGtoIPGradientWeightsZ( tPrimaryVertexIndex )( tCountZ ) = tdNdX( 2, iSF ); 

                                // increment counter
                                tCountZ++;
                                
                            } 
                        }                            
                    }
                }    // end: loop over vertices on cluster
            }    // end: loop over clusters
        }    // end: if cluster is empty
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::get_IP_to_BS_nodal_T_matrices(
            Vector< Vector< Matrix< IdMat > > >  &aIPtoBSIds,        // input: B-spline mesh index, UIPV index || output: B-spline BF IDs interpolating into UIPV
            Vector< Vector< Matrix< DDRMat > > > &aIPtoBSWeights,    // input: B-spline mesh index, UIPV index || output: T-matrix weights
            const uint                            aNumBsplineMeshes )
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
        Vector< bool > tUipvIndFound( tNumLagBfs, false );

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
                Vector< Vertex * > tIPVertices          = tInterpolationCell.get_vertex_pointers();
                uint               tNumIpVertsOnCluster = tIPVertices.size();

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
                        Matrix< DDRMat >  tCoords    = tIPVertices( iIpVert )->get_coords();

                        // copy the weights and basis IDs onto the output T-matrices
                        aIPtoBSIds( iBspMesh )( tUipvIndex )     = trans( tBSpIDs );
                        aIPtoBSWeights( iBspMesh )( tUipvIndex ) = trans( tBSWeights );

                    }    // end for: B-spline meshes for which the T-matrix should be retrieved
                }    // end for: unzipped IP vertices on cluster
            }    // end for: clusters in set
        }    // end for: sets in mesh
    }    // end function: Integration_Mesh::get_IP_to_BS_nodal_T_matrices()

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::get_IP_to_BS_nodal_T_matrices(
            Vector< Matrix< IdMat > >  &aIPtoBSIds,
            Vector< Matrix< DDRMat > > &aIPtoBSWeights,
            uint                        aSetIndex )
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
                Vector< Vertex * > tIPVertices = tInterpolationCell.get_vertex_pointers();

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
                        Matrix< DDRMat > tCoords = tIPVertices( iIpVert )->get_coords();

                        /*if( tCoords( 0 ) == 0.0 )
                        {
                            if( tCoords( 1 ) == 0.375 )
                            {
                                // print B-spline IDs and weights
                                std::cout << "B-spline IDs: " << tBSpIDs( 0 ) << std::endl;
                                std::cout << "B-spline Weights: " << tBSWeights( 0 ) << std::endl;
                                std::cout << "B-spline IDs: " << tBSpIDs( 1 ) << std::endl;
                                std::cout << "B-spline Weights: " << tBSWeights( 1 ) << std::endl;
                                std::cout << "B-spline IDs: " << tBSpIDs( 2 ) << std::endl;
                                std::cout << "B-spline Weights: " << tBSWeights( 2 ) << std::endl;
                                std::cout << "B-spline IDs: " << tBSpIDs( 3 ) << std::endl;
                                std::cout << "B-spline Weights: " << tBSWeights( 3 ) << std::endl;
                                //Print the coordinates
                                std::cout << "Coordinates: " << tCoords( 0 ) << std::endl;
                                std::cout << "Coordinates: " << tCoords( 1 ) << std::endl;
                            
                            }
                        }*/

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
            }    // end for: clusters
        }    // end if: cluster is non-empty
    }

    //----------------------------------------------------------------------------

    void
    Integration_Mesh::get_IP_to_BS_nodal_T_matrix_sensitivities(
        Vector< Matrix< IdMat > >  &aIPtoBSIds,
        Vector< Vector< Vector< Matrix< DDRMat > > > > &aIPtoBSSensitivityWeights,
        Vector< Vector< Vector< moris_id > > > &aSensitivityIGIDs,
        uint                        aSetIndex )
    {
        // aIPtoBSIds is same as in get_IP_to_BS_nodal_T_matrices, but aIPtoBSWeights is now a vector of matrices, one vector per BS vertex in IP to BS IDs.
        // aSensitivityIGIDs is a vector of matrices, one per BS vertex in IP to BS IDs. 
        // Per BS vertex, essentially we have all the IG nodes w.r.t. which it has sensitivities and their values.
        
        // trace this function
        Tracer tTracer( "MTK", "Compute IP to B-Spline T-Matrix Sensitivities" );

        // get number of IP nodes on mesh
        uint tNumIpNodes = this->get_max_IP_ID_on_set( aSetIndex );

        // initialize map size
        aIPtoBSIds.resize( tNumIpNodes );
        aIPtoBSSensitivityWeights.resize( tNumIpNodes );
        aSensitivityIGIDs.resize( tNumIpNodes );

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
                Vector< Vertex * > tIPVertices = tInterpolationCell.get_vertex_pointers();

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
                        Vector< Vector< real >  > tBSWeights = tIPVertices( iIpVert )->get_interpolation( 0 )->get_basis_weight_per_leader_basis();

                        // get list of numerator and denominator SPs associated with the B-Spline basis associated with the current IP vertex
                        Vector< Vector< moris_index > > tNumeratorSPs = tIPVertices( iIpVert )->get_interpolation( 0 )->get_basis_weight_numerators();
                        Vector< Vector< moris_index > > tDenominatorSPs = tIPVertices( iIpVert )->get_interpolation( 0 )->get_basis_weight_denominators();
                        Vector< Vector< uint > > tNumNumeratorSPs = tIPVertices( iIpVert )->get_interpolation( 0 )->get_basis_weight_numerator_count();
                        Vector< Vector< uint > > tNumDenominatorSPs = tIPVertices( iIpVert )->get_interpolation( 0 )->get_basis_weight_denominator_count();

                        // get number of B-splines associated with current IP vertex
                        uint tNumBspOnIpNode = tBSpIDs.numel();

                        // resize nodal T-matrix sensitivity IDs to accommodate all B-spline IDs, weights and associated SP IDs
                        aIPtoBSIds( tIP_ID - 1 ).set_size( 1, tNumBspOnIpNode, gNoID );
                        aIPtoBSSensitivityWeights( tIP_ID - 1 ).resize( tNumBspOnIpNode );
                        aSensitivityIGIDs( tIP_ID - 1 ).resize( tNumBspOnIpNode );
                        

                        // create T-matrix derivatives - copy list of associated B-spline IDs and their weights
                        for ( uint iBsp = 0; iBsp < tNumBspOnIpNode; iBsp++ )
                        {
                            // clear before proceeding to ensure that the sensitivity weights and IG IDs are empty before filling them in
                            aSensitivityIGIDs( tIP_ID - 1 )( iBsp ).clear();
                            aIPtoBSSensitivityWeights( tIP_ID - 1 )( iBsp ).clear();

                            std::unordered_map< moris_id, Matrix< DDRMat > > tVolumeFractionDerivMap;
                            aIPtoBSIds( tIP_ID - 1 )( iBsp )     = tBSpIDs( iBsp );
                            Vector< real > tBSpWtsPerLeaderBasis = tBSWeights( iBsp );
                            
                            // figure out which IG IDs are associated with each numerator/denominator pair for current B-Spline basis' volume fraction
                            // if numerator/denominator indices are empty/ only have a single -1 entry, then the sensitivity is zero and we don't need to store any IG IDs or weights
                            if ( tNumNumeratorSPs( iBsp ).size() >= 1 && tNumDenominatorSPs( iBsp ).size() >= 1 )
                            {
                                uint tNumeratorOffset   = 0;
                                uint tDenominatorOffset = 0;
                                for ( uint iNumSP = 0; iNumSP < tNumNumeratorSPs( iBsp ).size(); iNumSP++ )
                                {
                                    // get numerator and denominator IG cell indices
                                    const auto &tNumeratorIndices = tNumeratorSPs( iBsp );
                                    uint        tNumeratorCount   = tNumNumeratorSPs( iBsp )( iNumSP );

                                    Vector< moris_index > tNumIGIndex;
                                    tNumIGIndex.insert(
                                            0,
                                            tNumeratorIndices.begin() + tNumeratorOffset,
                                            tNumeratorIndices.begin() + tNumeratorOffset + tNumeratorCount );
                                    tNumeratorOffset += tNumeratorCount;

                                    const auto &tDenominatorIndices = tDenominatorSPs( iBsp );
                                    uint        tDenominatorCount   = tNumDenominatorSPs( iBsp )( iNumSP );

                                    Vector< moris_index > tDenIGIndex;
                                    tDenIGIndex.insert(
                                            0,
                                            tDenominatorIndices.begin() + tDenominatorOffset,
                                            tDenominatorIndices.begin() + tDenominatorOffset + tDenominatorCount );
                                    tDenominatorOffset += tDenominatorCount;

                                    // Compute measure of each IG cell in the numerator and denominator and store them in a vector
                                    // also compute sensitivities of cell measure w.r.t the IG vertices
                                    Vector< real > tNumIGAreas( tNumIGIndex.size() );
                                    std::unordered_map< moris_id, Matrix< DDRMat > > tNumIGVertexMap; 
                                    for ( uint iNumIG = 0; iNumIG < tNumIGIndex.size(); iNumIG++ )
                                    {
                                        // skip if the IG cell index is -1, which indicates that the there is a single extension operator for the current basis - not multiple (volume fraction 1)
                                        if ( tNumIGIndex( iNumIG ) == -1 )
                                        {
                                            continue;
                                        }
                                        // compute cell measure
                                        tNumIGAreas( iNumIG ) = this->get_mtk_cell( tNumIGIndex( iNumIG ) ).compute_cell_measure();

                                        // compute the sensitivities of the cell measure w.r.t the IG vertices
                                        for (uint iIGNodeID = 0; iIGNodeID < this->get_mtk_cell( tNumIGIndex( iNumIG ) ).get_vertex_ids().numel(); iIGNodeID++)
                                        {
                                            moris_id tIGNodeID = this->get_mtk_cell( tNumIGIndex( iNumIG ) ).get_vertex_ids()( iIGNodeID );
                                            Matrix< DDRMat > tCellSensitivities = this->get_mtk_cell( tNumIGIndex( iNumIG ) ).compute_cell_measure_deriv_general( iIGNodeID, -1 );
                                            if (tNumIGVertexMap.find( tIGNodeID ) == tNumIGVertexMap.end())
                                            {
                                                tNumIGVertexMap[tIGNodeID] = tCellSensitivities;
                                            }
                                            else
                                            {
                                                tNumIGVertexMap[tIGNodeID] += tCellSensitivities;
                                            }
                                        }
                                    }
                                    
                                    // Repeat the same process for the denominator IG cells
                                    Vector< real > tDenIGAreas( tDenIGIndex.size() );
                                    std::unordered_map< moris_id, Matrix< DDRMat > > tDenIGVertexMap;
                                    for ( uint iDenIG = 0; iDenIG < tDenIGIndex.size(); iDenIG++ )
                                    {
                                        // skip if the IG cell index is -1, which indicates that the there is a single extension operator for the current basis - not multiple, so no derivatives (volume fraction 1)
                                        if ( tDenIGIndex( iDenIG ) == -1 )
                                        {
                                            continue;
                                        }
                                        // compute cell measure
                                        tDenIGAreas( iDenIG ) = this->get_mtk_cell( tDenIGIndex( iDenIG ) ).compute_cell_measure();

                                        // compute the sensitivities of the cell measure w.r.t the IG vertices
                                        for ( uint iIGNodeID = 0; iIGNodeID < this->get_mtk_cell( tDenIGIndex( iDenIG ) ).get_vertex_ids().numel(); iIGNodeID++ )
                                        {
                                            moris_id         tIGNodeID          = this->get_mtk_cell( tDenIGIndex( iDenIG ) ).get_vertex_ids()( iIGNodeID );
                                            Matrix< DDRMat > tCellSensitivities = this->get_mtk_cell( tDenIGIndex( iDenIG ) ).compute_cell_measure_deriv_general( iIGNodeID, -1 );
                                            if ( tDenIGVertexMap.find( tIGNodeID ) == tDenIGVertexMap.end() )
                                            {
                                                tDenIGVertexMap[ tIGNodeID ] = tCellSensitivities;
                                            }
                                            else
                                            {
                                                tDenIGVertexMap[ tIGNodeID ] += tCellSensitivities;
                                            }
                                        }
                                    }

                                    // compute the derivative of the volume fraction w.r.t the IG vertices using the quotient rule
                                    for ( const auto &tPair : tDenIGVertexMap )
                                    {
                                        // get the vertex ID
                                        moris_id tIGNodeID = tPair.first;

                                        // get the derivative of denominator IG cells w.r.t. to this vertex
                                        Matrix< DDRMat > tDenDeriv = tPair.second;

                                        // check if this vertex is also in the numerator IG cells
                                        if ( tNumIGVertexMap.find( tIGNodeID ) != tNumIGVertexMap.end() )
                                        {
                                            // get the derivative of numerator IG cells w.r.t. to this vertex
                                            Matrix< DDRMat > tNumDeriv = tNumIGVertexMap[ tIGNodeID ];

                                            const real tNumArea =
                                                    std::accumulate( tNumIGAreas.begin(), tNumIGAreas.end(), real( 0 ) );
                                            const real tDenArea =
                                                    std::accumulate( tDenIGAreas.begin(), tDenIGAreas.end(), real( 0 ) );

                                            // compute the derivative of the volume fraction w.r.t this vertex using the quotient rule
                                            Matrix< DDRMat > tVolumeFractionDeriv = ( tNumDeriv * tDenArea - tDenDeriv * tNumArea ) / ( tDenArea * tDenArea );

                                            // compute derivative of b-spline weight this numerator-denominator pairing.
                                            if ( tVolumeFractionDerivMap.find( tIGNodeID ) == tVolumeFractionDerivMap.end() )
                                            {
                                                tVolumeFractionDerivMap[ tIGNodeID ] = tBSpWtsPerLeaderBasis( iNumSP ) * tVolumeFractionDeriv;
                                            }
                                            else
                                            {
                                                tVolumeFractionDerivMap[ tIGNodeID ] += tBSpWtsPerLeaderBasis( iNumSP ) * tVolumeFractionDeriv;
                                            }
                                        }
                                        else
                                        {
                                            // No contribution from numerator, so derivative is just the negative of the denominator contribution
                                            const real tNumArea =
                                                    std::accumulate( tNumIGAreas.begin(), tNumIGAreas.end(), real( 0 ) );
                                            const real tDenArea =
                                                    std::accumulate( tDenIGAreas.begin(), tDenIGAreas.end(), real( 0 ) );
                                            Matrix< DDRMat > tVolumeFractionDeriv = ( -tDenDeriv * tNumArea ) / ( tDenArea * tDenArea );

                                            // compute derivative of b-spline weight this numerator-denominator pairing.
                                            if ( tVolumeFractionDerivMap.find( tIGNodeID ) == tVolumeFractionDerivMap.end() )
                                            {
                                                tVolumeFractionDerivMap[ tIGNodeID ] = tBSpWtsPerLeaderBasis( iNumSP ) * tVolumeFractionDeriv;
                                            }
                                            else
                                            {
                                                tVolumeFractionDerivMap[ tIGNodeID ] += tBSpWtsPerLeaderBasis( iNumSP ) * tVolumeFractionDeriv;
                                            }

                                        }
                                        
                                    }


                                }
                            }
                            //  populate the output vectors with the IG IDs and their corresponding sensitivities for the current B-Spline basis
                            for ( const auto &tPair : tVolumeFractionDerivMap )
                            {
                                moris_id tIGNodeID = tPair.first;
                                Matrix< DDRMat > tVolumeFractionDeriv = tPair.second;
                                aSensitivityIGIDs( tIP_ID - 1 )( iBsp ).push_back( tIGNodeID );
                                aIPtoBSSensitivityWeights( tIP_ID - 1 )( iBsp ).push_back( tVolumeFractionDeriv );
                            }                         
                        }

                        // tick off IP vertex as treated
                        tIdentifierMatIP( tIP_ID - 1 ) = 1;
                    }
                }    // end for: IP vertices
            }    // end for: clusters
        }    // end if: cluster is non-empty
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::get_elemental_IG_to_BS_T_matrices(
            Vector< Matrix< IndexMat > > const         &aIgToIpIndices,      // outer cell: IG cell index
            Vector< Matrix< DDRMat > > const           &aIgToIpTmatrices,    // outer cell: IG cell index
            Vector< Vector< Matrix< IdMat > > > const  &aIPtoBSIds,          // outer cell: B-spline mesh index | inner cell: enr. Lagrange BF index
            Vector< Vector< Matrix< DDRMat > > > const &aIPtoBSWeights,      // outer cell: B-spline mesh index | inner cell: enr. Lagrange BF index
            Vector< Vector< Matrix< IdMat > > >        &aIGtoBSIds,          // outer cell: B-spline mesh index | inner cell: IG cell index
            Vector< Vector< Matrix< DDRMat > > >       &aIGtoBSWeights )           // outer cell: B-spline mesh index | inner cell: IG cell index
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
                }    // end for: Lagrange basis functions interpolating into IG element

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
                        // uint tLocalIndex = tIter->second;

                        // grab the part of the elemental T-matrix which the weights should be written into
                        // auto tCopyWeightsIntoThis = tCurrentIgCellToBSWeights( { 0, tNumNodesOnIgElem - 1 }, { tLocalIndex, tLocalIndex } );
                        // auto tCopyWeightsFromThis = tBspWeights( { 0, tNumNodesOnIgElem - 1 }, { iBspBf, iBspBf } );

                        // copy the weights into this
                        // tCopyWeightsIntoThis += tCopyWeightsFromThis;

                    }    // end for: B-spline basis functions interpolating into Lagrange basis function
                }    // end for: Lagrange basis functions interpolating into IG element
            }    // end for: IG elements
        }    // end for: B-spline meshes
    }    // Integration_Mesh::get_IG_to_BS_nodal_T_matrices()

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::get_IG_to_BS_nodal_T_matrices(
            Vector< Matrix< IdMat > >  &aIGtoBSIds,
            Vector< Matrix< IdMat > >  &aIGtoBSIdsX,
            Vector< Matrix< IdMat > >  &aIGtoBSIdsY,
            Vector< Matrix< IdMat > >  &aIGtoBSIdsZ,
            Vector< Matrix< DDRMat > > &aIGtoBSWeights,
            Vector< Matrix< DDRMat > > &aIGtoBSGradientWeightsX,
            Vector< Matrix< DDRMat > > &aIGtoBSGradientWeightsY,
            Vector< Matrix< DDRMat > > &aIGtoBSGradientWeightsZ,
            Vector< Matrix< IdMat > >  &aIPtoBSIds,
            Vector< Matrix< DDRMat > > &aIPtoBSWeights,
            Vector< Matrix< IdMat > >  &aIGtoIPIds,
            Vector< Matrix< IdMat > >  &aIGtoIPIdsX,
            Vector< Matrix< IdMat > >  &aIGtoIPIdsY,
            Vector< Matrix< IdMat > >  &aIGtoIPIdsZ,
            Vector< Matrix< DDRMat > > &aIGtoIPWeights,
            Vector< Matrix< DDRMat > > &aIGtoIPGradientWeightsX,
            Vector< Matrix< DDRMat > > &aIGtoIPGradientWeightsY,
            Vector< Matrix< DDRMat > > &aIGtoIPGradientWeightsZ )
    {
        // trace this function
        Tracer tTracer( "MTK", "Compute IG to IP T-Matrices" );

        // get number of IG vertices
        uint tNumIgNodes = aIGtoIPIds.size();
        MORIS_ERROR( tNumIgNodes == aIGtoIPWeights.size(),
                "Integration_Mesh::get_IG_to_BS_nodal_T_matrices() - Size of ID/weight maps for IG-IP T-Matrices don't match." );

        // initialize size of T-Matrix cells
        aIGtoBSIds.resize( tNumIgNodes );
        aIGtoBSIdsX.resize( tNumIgNodes );
        aIGtoBSIdsY.resize( tNumIgNodes );
        
        aIGtoBSWeights.resize( tNumIgNodes );
        aIGtoBSGradientWeightsX.resize( tNumIgNodes );
        aIGtoBSGradientWeightsY.resize( tNumIgNodes );
        
        if ( this->get_spatial_dim() == 3 )
        {
            aIGtoBSIdsZ.resize( tNumIgNodes );
            aIGtoBSGradientWeightsZ.resize( tNumIgNodes );
        }

        // loop over all IG vertices
        for ( uint iIgNode = 0; iIgNode < tNumIgNodes; iIgNode++ )
        {
            // initialize temporary vectors containing all
            Matrix< IdMat >  tNodalIGtoBSIds( 0, 0 );
            Matrix< IdMat >  tNodalIGtoBSIdsX( 0, 0 );
            Matrix< IdMat >  tNodalIGtoBSIdsY( 0, 0 );
            Matrix< IdMat >  tNodalIGtoBSIdsZ( 0, 0 );

            Matrix< DDRMat > tNodalIGtoBSWeights( 0, 0 );
            Matrix< DDRMat > tNodalIGtoBSGradientWeightsX( 0 , 0 );
            Matrix< DDRMat > tNodalIGtoBSGradientWeightsY( 0 , 0 );
            Matrix< DDRMat > tNodalIGtoBSGradientWeightsZ( 0 , 0 );

            // number of IP nodes associated with current IG node
            uint tNumIpBfs = aIGtoIPIds( iIgNode ).numel();

            // number of IP nodes associated with current IG node for grad x
            uint tNumIpBfsX = aIGtoIPIdsX( iIgNode ).numel();

            // number of IP nodes associated with current IG node for grad y
            uint tNumIpBfsY = aIGtoIPIdsY( iIgNode ).numel();

            // number of IP nodes associated with current IG node for grad z
            uint tNumIpBfsZ = 0;
            if ( this->get_spatial_dim() == 3 )
            {
                tNumIpBfsZ = aIGtoIPIdsZ( iIgNode ).numel();
            }

            // loop over all IP nodes current IG node is related to
            for ( uint jIpNode = 0; jIpNode < tNumIpBfs; jIpNode++ )
            {
                // get current IP ID
                moris_id tIpId = aIGtoIPIds( iIgNode )( jIpNode );

                // skip if IP-ID doesn't exist
                if ( tIpId > gNoID )
                {
                    // get what B-spline IDs and Weights are associated with current Node
                    const Matrix< IdMat > &tAddIdList     = aIPtoBSIds( tIpId - 1 );
                    Matrix< DDRMat >       tAddWeightList = aIGtoIPWeights( iIgNode )( jIpNode ) * aIPtoBSWeights( tIpId - 1 );
                    //Matrix< DDRMat >       tAddGradientWeightListX = aIGtoIPGradientWeightsX( iIgNode )( jIpNode ) * aIPtoBSWeights( tIpId - 1 );
                    //Matrix< DDRMat >       tAddGradientWeightListY = aIGtoIPGradientWeightsY( iIgNode )( jIpNode ) * aIPtoBSWeights( tIpId - 1 );

                    // append vectors of IDs and weights
                    tNodalIGtoBSIds     = join_horiz( tNodalIGtoBSIds, tAddIdList );
                    tNodalIGtoBSWeights = join_horiz( tNodalIGtoBSWeights, tAddWeightList );
                    //tNodalIGtoBSGradientWeightsX = join_horiz( tNodalIGtoBSGradientWeightsX, tAddGradientWeightListX );
                    //tNodalIGtoBSGradientWeightsY = join_horiz( tNodalIGtoBSGradientWeightsY, tAddGradientWeightListY );
                }
            }    // end: loop over all IP Nodes

            // loop over all IP nodes current IG node is related to for X
            for ( uint jIpNode = 0; jIpNode < tNumIpBfsX; jIpNode++ )
            {
                // get current IP ID
                moris_id tIpId = aIGtoIPIdsX( iIgNode )( jIpNode );

                if ( tIpId > gNoID )
                {
                    const Matrix< IdMat > &tAddIdList  = aIPtoBSIds( tIpId - 1 );
                    Matrix< DDRMat >       tAddWeightList = aIGtoIPGradientWeightsX( iIgNode )( jIpNode ) * aIPtoBSWeights( tIpId - 1 );

                    // append vectors of IDs and weights
                    tNodalIGtoBSIdsX     = join_horiz( tNodalIGtoBSIdsX, tAddIdList );
                    tNodalIGtoBSGradientWeightsX = join_horiz( tNodalIGtoBSGradientWeightsX, tAddWeightList );
                }
            }

            // loop over all IP nodes current IG node is related to for Y
            for ( uint jIpNode = 0; jIpNode < tNumIpBfsY; jIpNode++ )
            {
                // get current IP ID
                moris_id tIpId = aIGtoIPIdsY( iIgNode )( jIpNode );

                if ( tIpId > gNoID )
                {
                    const Matrix< IdMat > &tAddIdList  = aIPtoBSIds( tIpId - 1 );
                    Matrix< DDRMat >       tAddWeightList = aIGtoIPGradientWeightsY( iIgNode )( jIpNode ) * aIPtoBSWeights( tIpId - 1 );

                    // append vectors of IDs and weights
                    tNodalIGtoBSIdsY     = join_horiz( tNodalIGtoBSIdsY, tAddIdList );
                    tNodalIGtoBSGradientWeightsY = join_horiz( tNodalIGtoBSGradientWeightsY, tAddWeightList );
                }
            }

            // if spatial dimension is 3, loop over all IP nodes current IG node is related to for Z
            if ( this->get_spatial_dim() == 3 )
            {
                for ( uint jIpNode = 0; jIpNode < tNumIpBfsZ; jIpNode++ )
                {
                    // get current IP ID
                    moris_id tIpId = aIGtoIPIdsZ( iIgNode )( jIpNode );

                    if ( tIpId > gNoID )
                    {
                        const Matrix< IdMat > &tAddIdList     = aIPtoBSIds( tIpId - 1 );
                        Matrix< DDRMat >       tAddWeightList = aIGtoIPGradientWeightsZ( iIgNode )( jIpNode ) * aIPtoBSWeights( tIpId - 1 );

                        // append vectors of IDs and weights
                        tNodalIGtoBSIdsZ             = join_horiz( tNodalIGtoBSIdsZ, tAddIdList );
                        tNodalIGtoBSGradientWeightsZ = join_horiz( tNodalIGtoBSGradientWeightsZ, tAddWeightList );
                    }
                }
            }

            // get unique list of B-Spline IDs associated with each IG node
            moris::unique( tNodalIGtoBSIds, aIGtoBSIds( iIgNode ) );

            // get unique list of B-Spline IDs associated with each IG node for x
            moris::unique( tNodalIGtoBSIdsX, aIGtoBSIdsX( iIgNode ) );

            // get unique list of B-Spline IDs associated with each IG node for y
            moris::unique( tNodalIGtoBSIdsY, aIGtoBSIdsY( iIgNode ) );

            // get unique list of B-Spline IDs associated with each IG node for z
            if (this->get_spatial_dim() == 3 )
            {
                moris::unique( tNodalIGtoBSIdsZ, aIGtoBSIdsZ( iIgNode ) );
            }

            // get number of unique IDs in list
            uint tNumUniqueIDs = aIGtoBSIds( iIgNode ).numel();

            // get number of unique IDs in list X
            uint tNumUniqueIDsX = aIGtoBSIdsX( iIgNode ).numel();

            // get number of unique IDs in list Y
            uint tNumUniqueIDsY = aIGtoBSIdsY( iIgNode ).numel();

            // get number of unique IDs in list Z
            uint tNumUniqueIDsZ = 0;
            if ( this->get_spatial_dim() == 3 )
            {
                tNumUniqueIDsZ = aIGtoBSIdsZ( iIgNode ).numel();
            }

            // initialize list of weights associated with unique list of IDs
            aIGtoBSWeights( iIgNode ).set_size( tNumUniqueIDs, 1, 0.0 );
            aIGtoBSGradientWeightsX( iIgNode ).set_size( tNumUniqueIDsX, 1, 0.0 );
            aIGtoBSGradientWeightsY( iIgNode ).set_size( tNumUniqueIDsY, 1, 0.0 );
            if ( this->get_spatial_dim() == 3 )
            {
                aIGtoBSGradientWeightsZ( iIgNode ).set_size( tNumUniqueIDsZ, 1, 0.0 );
            }

            // create list of weights associated with unique list of IDs
            for ( uint iUnique = 0; iUnique < tNumUniqueIDs; iUnique++ )
            {
                for ( uint jNonUnique = 0; jNonUnique < tNodalIGtoBSIds.numel(); jNonUnique++ )
                {
                    if ( aIGtoBSIds( iIgNode )( iUnique ) == tNodalIGtoBSIds( jNonUnique ) )
                    {
                        aIGtoBSWeights( iIgNode )( iUnique ) = aIGtoBSWeights( iIgNode )( iUnique ) + tNodalIGtoBSWeights( jNonUnique );
                        //aIGtoBSGradientWeightsX( iIgNode )( iUnique ) = aIGtoBSGradientWeightsX( iIgNode )( iUnique ) + tNodalIGtoBSGradientWeightsX( jNonUnique );
                        //aIGtoBSGradientWeightsY( iIgNode )( iUnique ) = aIGtoBSGradientWeightsY( iIgNode )( iUnique ) + tNodalIGtoBSGradientWeightsY( jNonUnique );
                    }
                }
            }

            // create list of weights associated with unique list of IDs x
            for ( uint iUnique = 0; iUnique < tNumUniqueIDsX; iUnique++ )
            {
                for ( uint jNonUnique = 0; jNonUnique < tNodalIGtoBSIdsX.numel(); jNonUnique++ )
                {
                    if ( aIGtoBSIdsX( iIgNode )( iUnique ) == tNodalIGtoBSIdsX( jNonUnique ) )
                    {
                        aIGtoBSGradientWeightsX( iIgNode )( iUnique ) = aIGtoBSGradientWeightsX( iIgNode )( iUnique ) + tNodalIGtoBSGradientWeightsX( jNonUnique );
                        //aIGtoBSGradientWeightsX( iIgNode )( iUnique ) = aIGtoBSGradientWeightsX( iIgNode )( iUnique ) + tNodalIGtoBSGradientWeightsX( jNonUnique );
                        //aIGtoBSGradientWeightsY( iIgNode )( iUnique ) = aIGtoBSGradientWeightsY( iIgNode )( iUnique ) + tNodalIGtoBSGradientWeightsY( jNonUnique );
                    }
                }
            }

            // create list of weights associated with unique list of IDs y
            for ( uint iUnique = 0; iUnique < tNumUniqueIDsY; iUnique++ )
            {
                for ( uint jNonUnique = 0; jNonUnique < tNodalIGtoBSIdsY.numel(); jNonUnique++ )
                {
                    if ( aIGtoBSIdsY( iIgNode )( iUnique ) == tNodalIGtoBSIdsY( jNonUnique ) )
                    {
                        aIGtoBSGradientWeightsY( iIgNode )( iUnique ) = aIGtoBSGradientWeightsY( iIgNode )( iUnique ) + tNodalIGtoBSGradientWeightsY( jNonUnique );
                        //aIGtoBSGradientWeightsX( iIgNode )( iUnique ) = aIGtoBSGradientWeightsX( iIgNode )( iUnique ) + tNodalIGtoBSGradientWeightsX( jNonUnique );
                        //aIGtoBSGradientWeightsY( iIgNode )( iUnique ) = aIGtoBSGradientWeightsY( iIgNode )( iUnique ) + tNodalIGtoBSGradientWeightsY( jNonUnique );
                    }
                }
            }
            // create list of weights associated with unique list of IDs z
            if ( this->get_spatial_dim() == 3 )
            {
                for ( uint iUnique = 0; iUnique < tNumUniqueIDsZ; iUnique++ )
                {
                    for ( uint jNonUnique = 0; jNonUnique < tNodalIGtoBSIdsZ.numel(); jNonUnique++ )
                    {
                        if ( aIGtoBSIdsZ( iIgNode )( iUnique ) == tNodalIGtoBSIdsZ( jNonUnique ) )
                        {
                            aIGtoBSGradientWeightsZ( iIgNode )( iUnique ) = aIGtoBSGradientWeightsZ( iIgNode )( iUnique ) + tNodalIGtoBSGradientWeightsZ( jNonUnique );
                            //aIGtoBSGradientWeightsX( iIgNode )( iUnique ) = aIGtoBSGradientWeightsX( iIgNode )( iUnique ) + tNodalIGtoBSGradientWeightsX( jNonUnique );
                            //aIGtoBSGradientWeightsY( iIgNode )( iUnique ) = aIGtoBSGradientWeightsY( iIgNode )( iUnique ) + tNodalIGtoBSGradientWeightsY( jNonUnique );
                        }
                    }
                }
            }
        }    // end: loop over all IG Nodes
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh::get_IG_to_BS_nodal_T_matrices(
            Vector< Matrix< IdMat > >                      &aIGtoBSIds,
            Vector< Matrix< IdMat > >                      &aIGtoBSIdsForSensitivityX,
            Vector< Matrix< IdMat > >                      &aIGtoBSIdsForSensitivityY,
            Vector< Matrix< IdMat > >                      &aIGtoBSIdsForSensitivityZ,
            Vector< Vector< Vector< moris_id > > >         &aIGtoBSSensitivityIdsX,
            Vector< Vector< Vector< moris_id > > >         &aIGtoBSSensitivityIdsY,
            Vector< Vector< Vector< moris_id > > >         &aIGtoBSSensitivityIdsZ,
            Vector< Matrix< DDRMat > >                     &aIGtoBSWeights,
            Vector< Vector< Vector< real > > >             &aIGtoBSGradientWeightsX,
            Vector< Vector< Vector< real > > >             &aIGtoBSGradientWeightsY,
            Vector< Vector< Vector< real > > >             &aIGtoBSGradientWeightsZ,
            Vector< Matrix< IdMat > >                      &aIPtoBSIds,
            Vector< Matrix< DDRMat > >                     &aIPtoBSWeights,
            Vector< Vector< Vector< Matrix< DDRMat > > > > &aIPtoBSSensitivityWeights,
            Vector< Vector< Vector< moris_id > > >         &aSensitivityIGIDs,
            Vector< Matrix< IdMat > >                      &aIGtoIPIds,
            Vector< Matrix< IdMat > >                      &aIGtoIPIdsX,
            Vector< Matrix< IdMat > >                      &aIGtoIPIdsY,
            Vector< Matrix< IdMat > >                      &aIGtoIPIdsZ,
            Vector< Matrix< DDRMat > >                     &aIGtoIPWeights,
            Vector< Matrix< DDRMat > >                     &aIGtoIPGradientWeightsX,
            Vector< Matrix< DDRMat > >                     &aIGtoIPGradientWeightsY,
            Vector< Matrix< DDRMat > >                     &aIGtoIPGradientWeightsZ )
    {
        // trace this function
        Tracer tTracer( "MTK", "Compute IG to IP T-Matrices" );

        // get number of IG vertices
        uint tNumIgNodes = aIGtoIPIds.size();
        MORIS_ERROR( tNumIgNodes == aIGtoIPWeights.size(),
                "Integration_Mesh::get_IG_to_BS_nodal_T_matrices() - Size of ID/weight maps for IG-IP T-Matrices don't match." );

        // initialize size of T-Matrix cells
        aIGtoBSIds.resize( tNumIgNodes );
        aIGtoBSIdsForSensitivityX.resize( tNumIgNodes );
        aIGtoBSIdsForSensitivityY.resize( tNumIgNodes );
        
        aIGtoBSSensitivityIdsX.resize( tNumIgNodes );
        aIGtoBSSensitivityIdsY.resize( tNumIgNodes );
        
        if ( this->get_spatial_dim() == 3 )
        {
            aIGtoBSIdsForSensitivityZ.resize( tNumIgNodes );
            aIGtoBSSensitivityIdsZ.resize( tNumIgNodes );
        }

        aIGtoBSWeights.resize( tNumIgNodes );
        aIGtoBSGradientWeightsX.resize( tNumIgNodes );
        aIGtoBSGradientWeightsY.resize( tNumIgNodes );

        if ( this->get_spatial_dim() == 3 )
        {
            aIGtoBSGradientWeightsZ.resize( tNumIgNodes );
        }

        // loop over all IG vertices
        for ( uint iIgNode = 0; iIgNode < tNumIgNodes; iIgNode++ )
        {
            // initialize temporary vectors containing all
            Matrix< IdMat >              tNodalIGtoBSIds( 0, 0 );
            Vector< Vector< moris_id > > tNodalIGtoBSSensitivityIdsX( 0 );
            Matrix< IdMat >              tNodalIGtoBSIdsForSensitivityX( 0 );
            Vector< Vector< moris_id > > tNodalIGtoBSSensitivityIdsY( 0 );
            Matrix< IdMat >              tNodalIGtoBSIdsForSensitivityY( 0 );
            Vector< Vector< moris_id > > tNodalIGtoBSSensitivityIdsZ( 0 );
            Matrix< IdMat >              tNodalIGtoBSIdsForSensitivityZ( 0 );
            

            Matrix< DDRMat >          tNodalIGtoBSWeights( 0, 0 );
            Vector< Vector< real > >  tNodalIGtoBSGradientWeightsX( 0 );
            Vector< Vector< real > >  tNodalIGtoBSGradientWeightsY( 0 );
            Vector< Vector< real > >  tNodalIGtoBSGradientWeightsZ( 0 );

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
                    const Matrix< IdMat > &tAddIdList     = aIPtoBSIds( tIpId - 1 );
                    Matrix< DDRMat >       tAddWeightList = aIGtoIPWeights( iIgNode )( jIpNode ) * aIPtoBSWeights( tIpId - 1 );

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

            // There are two terms here - one is derivative of FG Lagrange to BG Lagrange times Lagrange extrarction and second term is 
            // derivative of Lagrange extraction times FG Lagrange to BG Lagrange. 

            // Computing the first term - derivative of FG Lagrange to BG Lagrange times Lagrange extraction
            uint tNumIpBfsX = aIGtoIPIdsX( iIgNode ).numel();
            for ( uint jIpNode = 0; jIpNode < tNumIpBfsX; jIpNode++ )
            {
                // get current IP ID
                moris_id tIpId = aIGtoIPIdsX( iIgNode )( jIpNode );
                if ( tIpId > gNoID )
                {
                    const Matrix< IdMat > &tAddIdList  = aIPtoBSIds( tIpId - 1 );
                    Matrix< DDRMat >       tAddWeightList = aIGtoIPGradientWeightsX( iIgNode )( jIpNode ) * aIPtoBSWeights( tIpId - 1 );

                    // append vectors of IDs and weights
                    tNodalIGtoBSIdsForSensitivityX = join_horiz( tNodalIGtoBSIdsForSensitivityX, tAddIdList );

                    for ( uint iBf = 0; iBf < tAddWeightList.numel(); iBf++ )
                    {
                        Vector< real > tAddGradientWeightList( 1 );
                        tAddGradientWeightList( 0 ) = tAddWeightList( iBf );
                        tNodalIGtoBSGradientWeightsX.push_back( tAddGradientWeightList );
                    }

                    // append IG element IDs for sensitivity
                    for ( uint iBf = 0; iBf < tAddIdList.numel(); iBf++ )
                    {
                        Vector< moris_id > tAddGradientIdList( 1 );
                        const moris_id     tIGNodeId = this->get_mtk_vertex( iIgNode ).get_id();
                        tAddGradientIdList( 0 )      = tIGNodeId;
                        tNodalIGtoBSSensitivityIdsX.push_back( tAddGradientIdList );
                    }
                }
            }
            // computing second term - derivative of Lagrange extraction times FG Lagrange to BG Lagrange
            for ( uint jIpNode = 0; jIpNode < tNumIpBfs; jIpNode++ )
            {
                // get current IP ID
                moris_id tIpId = aIGtoIPIds( iIgNode )( jIpNode );
                if ( tIpId > gNoID )
                {
                    const Matrix< IdMat > &tAddIdList  = aIPtoBSIds( tIpId - 1 );
                    tNodalIGtoBSIdsForSensitivityX = join_horiz( tNodalIGtoBSIdsForSensitivityX, tAddIdList );

                    // loop over the sensitivities and multiply with IG to IP weights for all 
                    // B spline basis functions associated with this node.
                    for ( uint iBf = 0; iBf < tAddIdList.numel(); iBf++ )
                    {
                        // get basis function ID
                        moris_id tBfId = tAddIdList( iBf );
                        
                        // get all the sensitivities corresponding to this basis function ID
                        Vector< Matrix< DDRMat > > tSensitivityWeights = aIPtoBSSensitivityWeights( tIpId - 1 )( iBf );

                        // if tSensitivityWeights is empty the loop below won't run.

                        // declare vectors for holding sensitivity values and IDs
                        Vector< real > tAddGradientWeightList;
                        Vector< moris_id > tAddGradientIdList;
                        const auto        &tSensitivityIds = aSensitivityIGIDs( tIpId - 1 )( iBf );

                        MORIS_ERROR(
                                tSensitivityWeights.size() == tSensitivityIds.size(),
                                "Sensitivity weights and IG node IDs have different lengths." );

                        // we are assembling the x gradient here, so we need to multiply x component of tSensitivityWeights with the IG to IP weights
                        for ( uint iSens = 0; iSens < tSensitivityWeights.size(); iSens++ )
                        {
                            // get the sensitivity weight
                            Matrix< DDRMat > tSensWeight = tSensitivityWeights( iSens );

                            MORIS_ERROR(
                                    tSensWeight.n_rows() == this->get_spatial_dim()
                                            && tSensWeight.n_cols() == 1,
                                    "Sensitivity gradient has an unexpected shape." );

                            // multiply with IG to IP weight
                            real tWeight = aIGtoIPWeights( iIgNode )( jIpNode ) * tSensWeight( 0 );

                            tAddGradientWeightList.push_back( tWeight );
                            tAddGradientIdList.push_back( tSensitivityIds( iSens ) );
                        }

                        // Append these vectors into the main vectors
                        tNodalIGtoBSGradientWeightsX.push_back( tAddGradientWeightList );
                        tNodalIGtoBSSensitivityIdsX.push_back( tAddGradientIdList );
                    }

                    MORIS_ASSERT(
                            tNodalIGtoBSIdsForSensitivityX.numel()
                                            == tNodalIGtoBSGradientWeightsX.size()
                                    && tNodalIGtoBSGradientWeightsX.size()
                                               == tNodalIGtoBSSensitivityIdsX.size(),
                            "Temporary X sensitivity arrays are misaligned." );
                }
            }

            // repeating for Y components
            // Computing the first term - derivative of FG Lagrange to BG Lagrange times Lagrange extraction
            uint tNumIpBfsY = aIGtoIPIdsY( iIgNode ).numel();
            for ( uint jIpNode = 0; jIpNode < tNumIpBfsY; jIpNode++ )
            {
                // get current IP ID
                moris_id tIpId = aIGtoIPIdsY( iIgNode )( jIpNode );
                if ( tIpId > gNoID )
                {
                    const Matrix< IdMat > &tAddIdList     = aIPtoBSIds( tIpId - 1 );
                    Matrix< DDRMat >       tAddWeightList = aIGtoIPGradientWeightsY( iIgNode )( jIpNode ) * aIPtoBSWeights( tIpId - 1 );

                    // append vectors of IDs and weights
                    tNodalIGtoBSIdsForSensitivityY = join_horiz( tNodalIGtoBSIdsForSensitivityY, tAddIdList );

                    for ( uint iBf = 0; iBf < tAddWeightList.numel(); iBf++ )
                    {
                        Vector< real > tAddGradientWeightList( 1 );
                        tAddGradientWeightList( 0 ) = tAddWeightList( iBf );
                        tNodalIGtoBSGradientWeightsY.push_back( tAddGradientWeightList );
                    }

                    // append IG element IDs for sensitivity
                    for ( uint iBf = 0; iBf < tAddIdList.numel(); iBf++ )
                    {
                        Vector< moris_id > tAddGradientIdList( 1 );
                        const moris_id     tIGNodeId = this->get_mtk_vertex( iIgNode ).get_id();
                        tAddGradientIdList( 0 )      = tIGNodeId;
                        tNodalIGtoBSSensitivityIdsY.push_back( tAddGradientIdList );
                    }
                }
            }
            // computing second term - derivative of Lagrange extraction times FG Lagrange to BG Lagrange
            for ( uint jIpNode = 0; jIpNode < tNumIpBfs; jIpNode++ )
            {
                // get current IP ID
                moris_id tIpId = aIGtoIPIds( iIgNode )( jIpNode );
                if ( tIpId > gNoID )
                {
                    const Matrix< IdMat > &tAddIdList = aIPtoBSIds( tIpId - 1 );
                    tNodalIGtoBSIdsForSensitivityY    = join_horiz( tNodalIGtoBSIdsForSensitivityY, tAddIdList );

                    // loop over the sensitivities and multiply with IG to IP weights for all
                    // B spline basis functions associated with this node.
                    for ( uint iBf = 0; iBf < tAddIdList.numel(); iBf++ )
                    {
                        // get basis function ID
                        moris_id tBfId = tAddIdList( iBf );

                        // get all the sensitivities corresponding to this basis function ID
                        Vector< Matrix< DDRMat > > tSensitivityWeights = aIPtoBSSensitivityWeights( tIpId - 1 )( iBf );

                        // declare vectors for holding sensitivity values and IDs
                        Vector< real >     tAddGradientWeightList;
                        Vector< moris_id > tAddGradientIdList;
                        const auto        &tSensitivityIds = aSensitivityIGIDs( tIpId - 1 )( iBf );

                        MORIS_ERROR(
                                tSensitivityWeights.size() == tSensitivityIds.size(),
                                "Sensitivity weights and IG node IDs have different lengths." );

                        // we are assembling the x gradient here, so we need to multiply x component of tSensitivityWeights with the IG to IP weights
                        for ( uint iSens = 0; iSens < tSensitivityWeights.size(); iSens++ )
                        {
                            // get the sensitivity weight
                            Matrix< DDRMat > tSensWeight = tSensitivityWeights( iSens );

                            MORIS_ERROR(
                                    tSensWeight.n_rows() == this->get_spatial_dim()
                                            && tSensWeight.n_cols() == 1,
                                    "Sensitivity gradient has an unexpected shape." );

                            // multiply with IG to IP weight
                            real tWeight = aIGtoIPWeights( iIgNode )( jIpNode ) * tSensWeight( 1 );

                            tAddGradientWeightList.push_back( tWeight );
                            tAddGradientIdList.push_back( tSensitivityIds( iSens ) );
                        }

                        // Append these vectors into the main vectors
                        tNodalIGtoBSGradientWeightsY.push_back( tAddGradientWeightList );
                        tNodalIGtoBSSensitivityIdsY.push_back( tAddGradientIdList );
                    }

                    MORIS_ASSERT(
                            tNodalIGtoBSIdsForSensitivityY.numel()
                                            == tNodalIGtoBSGradientWeightsY.size()
                                    && tNodalIGtoBSGradientWeightsY.size()
                                               == tNodalIGtoBSSensitivityIdsY.size(),
                            "Temporary X sensitivity arrays are misaligned." );
                }
            }

            // if problem is 3D, repeat for Z components
            if ( this->get_spatial_dim() == 3 )
            {
                // Computing the first term - derivative of FG Lagrange to BG Lagrange times Lagrange extraction
                uint tNumIpBfsZ = aIGtoIPIdsZ( iIgNode ).numel();
                for ( uint jIpNode = 0; jIpNode < tNumIpBfsZ; jIpNode++ )
                {
                    // get current IP ID
                    moris_id tIpId = aIGtoIPIdsZ( iIgNode )( jIpNode );
                    if ( tIpId > gNoID )
                    {
                        const Matrix< IdMat > &tAddIdList     = aIPtoBSIds( tIpId - 1 );
                        Matrix< DDRMat >       tAddWeightList = aIGtoIPGradientWeightsZ( iIgNode )( jIpNode ) * aIPtoBSWeights( tIpId - 1 );

                        // append vectors of IDs and weights
                        tNodalIGtoBSIdsForSensitivityZ = join_horiz( tNodalIGtoBSIdsForSensitivityZ, tAddIdList );

                        for ( uint iBf = 0; iBf < tAddWeightList.numel(); iBf++ )
                        {
                            Vector< real > tAddGradientWeightList( 1 );
                            tAddGradientWeightList( 0 ) = tAddWeightList( iBf );
                            tNodalIGtoBSGradientWeightsZ.push_back( tAddGradientWeightList );
                        }

                        // append IG element IDs for sensitivity
                        for ( uint iBf = 0; iBf < tAddIdList.numel(); iBf++ )
                        {
                            Vector< moris_id > tAddGradientIdList( 1 );
                            const moris_id     tIGNodeId = this->get_mtk_vertex( iIgNode ).get_id();
                            tAddGradientIdList( 0 )      = tIGNodeId;
                            tNodalIGtoBSSensitivityIdsZ.push_back( tAddGradientIdList );
                        }
                    }
                }
                // computing second term - derivative of Lagrange extraction times FG Lagrange to BG Lagrange
                for ( uint jIpNode = 0; jIpNode < tNumIpBfs; jIpNode++ )
                {
                    // get current IP ID
                    moris_id tIpId = aIGtoIPIds( iIgNode )( jIpNode );
                    if ( tIpId > gNoID )
                    {
                        const Matrix< IdMat > &tAddIdList = aIPtoBSIds( tIpId - 1 );
                        tNodalIGtoBSIdsForSensitivityZ    = join_horiz( tNodalIGtoBSIdsForSensitivityZ, tAddIdList );

                        // loop over the sensitivities and multiply with IG to IP weights for all
                        // B spline basis functions associated with this node.
                        for ( uint iBf = 0; iBf < tAddIdList.numel(); iBf++ )
                        {
                            // get basis function ID
                            moris_id tBfId = tAddIdList( iBf );

                            // get all the sensitivities corresponding to this basis function ID
                            Vector< Matrix< DDRMat > > tSensitivityWeights = aIPtoBSSensitivityWeights( tIpId - 1 )( iBf );

                            // declare vectors for holding sensitivity values and IDs
                            Vector< real >     tAddGradientWeightList;
                            Vector< moris_id > tAddGradientIdList;
                            const auto        &tSensitivityIds = aSensitivityIGIDs( tIpId - 1 )( iBf );

                            MORIS_ERROR(
                                    tSensitivityWeights.size() == tSensitivityIds.size(),
                                    "Sensitivity weights and IG node IDs have different lengths." );

                            // we are assembling the x gradient here, so we need to multiply x component of tSensitivityWeights with the IG to IP weights
                            for ( uint iSens = 0; iSens < tSensitivityWeights.size(); iSens++ )
                            {
                                // get the sensitivity weight
                                Matrix< DDRMat > tSensWeight = tSensitivityWeights( iSens );

                                MORIS_ERROR(
                                        tSensWeight.n_rows() == this->get_spatial_dim()
                                                && tSensWeight.n_cols() == 1,
                                        "Sensitivity gradient has an unexpected shape." );

                                // multiply with IG to IP weight
                                real tWeight = aIGtoIPWeights( iIgNode )( jIpNode ) * tSensWeight( 2 );

                                tAddGradientWeightList.push_back( tWeight );
                                tAddGradientIdList.push_back( tSensitivityIds( iSens ) );
                            }

                            // Append these vectors into the main vectors
                            tNodalIGtoBSGradientWeightsZ.push_back( tAddGradientWeightList );
                            tNodalIGtoBSSensitivityIdsZ.push_back( tAddGradientIdList );

                        }

                        MORIS_ASSERT(
                                tNodalIGtoBSIdsForSensitivityZ.numel()
                                                == tNodalIGtoBSGradientWeightsZ.size()
                                        && tNodalIGtoBSGradientWeightsZ.size()
                                                   == tNodalIGtoBSSensitivityIdsZ.size(),
                                "Temporary X sensitivity arrays are misaligned." );
                    }
                }
            }

            // create a list of unique B Spline IDs associated with each IG node for sensitivity X
            moris::unique( tNodalIGtoBSIdsForSensitivityX, aIGtoBSIdsForSensitivityX( iIgNode ) );

            // get number of unique IDs in list
            uint tNumUniqueIDsX = aIGtoBSIdsForSensitivityX( iIgNode ).numel();

            // initialize list of weights associated with unique list of IDs
            aIGtoBSGradientWeightsX( iIgNode ).resize( tNumUniqueIDsX );

            aIGtoBSSensitivityIdsX( iIgNode ).resize( tNumUniqueIDsX );

            // create a list of unique B Spline IDs associated with each IG node for sensitivity Y
            moris::unique( tNodalIGtoBSIdsForSensitivityY, aIGtoBSIdsForSensitivityY( iIgNode ) );

            // get number of unique IDs in list
            uint tNumUniqueIDsY = aIGtoBSIdsForSensitivityY( iIgNode ).numel();

            // initialize list of weights associated with unique list of IDs
            aIGtoBSGradientWeightsY( iIgNode ).resize( tNumUniqueIDsY );

            aIGtoBSSensitivityIdsY( iIgNode ).resize( tNumUniqueIDsY );

            uint tNumUniqueIDsZ = 0;
            if ( this->get_spatial_dim() == 3 )
            {
                // create a list of unique B Spline IDs associated with each IG node for sensitivity Z
                moris::unique( tNodalIGtoBSIdsForSensitivityZ, aIGtoBSIdsForSensitivityZ( iIgNode ) );

                // get number of unique IDs in list
                tNumUniqueIDsZ = aIGtoBSIdsForSensitivityZ( iIgNode ).numel();

                // initialize list of weights associated with unique list of IDs
                aIGtoBSGradientWeightsZ( iIgNode ).resize( tNumUniqueIDsZ );

                aIGtoBSSensitivityIdsZ( iIgNode ).resize( tNumUniqueIDsZ );
            }

            // Loop over the sensitivity weights and see if there are repeated contributions from the same B-spline functions.
            // if that is the case, if the IG node with respect to which there are derivatives also match, we add contributions.
            // otherwise, the vector of sensitivities associated with each B-spline function will get an append, so will the vector of IG nodes with respect to which the sensitivities are computed.

            for ( uint iUnique = 0; iUnique < tNumUniqueIDsX; iUnique++ )
            {
                aIGtoBSGradientWeightsX( iIgNode )( iUnique ).clear();
                aIGtoBSSensitivityIdsX( iIgNode )( iUnique ).clear();
                for ( uint jNonUnique = 0; jNonUnique < tNodalIGtoBSIdsForSensitivityX.numel(); jNonUnique++ )
                {
                    if ( aIGtoBSIdsForSensitivityX( iIgNode )( iUnique ) == tNodalIGtoBSIdsForSensitivityX( jNonUnique ) )
                    {
                        // get the vector of weights from the non-unique list
                        Vector< real > tWeightVector = tNodalIGtoBSGradientWeightsX( jNonUnique );

                        // get the vector of IG nodes from the non-unique list
                        Vector< moris_id > tIGNodeVector = tNodalIGtoBSSensitivityIdsX( jNonUnique );

                        // get the vector of weights from the unique list
                        Vector< real > &tUniqueWeightVector = aIGtoBSGradientWeightsX( iIgNode )( iUnique );

                        // get the vector of IG nodes from the unique list
                        Vector< moris_id > &tUniqueIGNodeVector = aIGtoBSSensitivityIdsX( iIgNode )( iUnique );

                        // skip if it is first instance, because unique list will be empty and we can just append the non-unique list to it.
                        if ( tUniqueWeightVector.size() == 0 )
                        {
                            tUniqueWeightVector = tWeightVector;
                            tUniqueIGNodeVector = tIGNodeVector;    
                        }
                        else
                        {
                            // loop over the non-unique list and see if there are any IG nodes that match with the unique list. If so, add contributions.
                            for ( uint iNonUnique = 0; iNonUnique < tIGNodeVector.size(); iNonUnique++ )
                            {
                                moris_id tNonUniqueIGNode = tIGNodeVector( iNonUnique );
                                real tNonUniqueWeight = tWeightVector( iNonUnique );

                                // check if this IG node is already in the unique list
                                auto it = std::find( tUniqueIGNodeVector.begin(), tUniqueIGNodeVector.end(), tNonUniqueIGNode );
                                if ( it != tUniqueIGNodeVector.end() )
                                {
                                    // if found, add the weight to the corresponding index in the unique weight vector
                                    uint index = std::distance( tUniqueIGNodeVector.begin(), it );
                                    tUniqueWeightVector( index ) += tNonUniqueWeight;
                                }
                                else
                                {
                                    // if not found, append the IG node and weight to the unique lists
                                    tUniqueIGNodeVector.push_back( tNonUniqueIGNode );
                                    tUniqueWeightVector.push_back( tNonUniqueWeight );
                                }
                            }
                        }
                    }
                }
            }

            // repeat for Y and Z components
            for ( uint iUnique = 0; iUnique < tNumUniqueIDsY; iUnique++ )
            {
                aIGtoBSGradientWeightsY( iIgNode )( iUnique ).clear();
                aIGtoBSSensitivityIdsY( iIgNode )( iUnique ).clear();
                for ( uint jNonUnique = 0; jNonUnique < tNodalIGtoBSIdsForSensitivityY.numel(); jNonUnique++ )
                {
                    if ( aIGtoBSIdsForSensitivityY( iIgNode )( iUnique ) == tNodalIGtoBSIdsForSensitivityY( jNonUnique ) )
                    {
                        // get the vector of weights from the non-unique list
                        Vector< real > tWeightVector = tNodalIGtoBSGradientWeightsY( jNonUnique );

                        // get the vector of IG nodes from the non-unique list
                        Vector< moris_id > tIGNodeVector = tNodalIGtoBSSensitivityIdsY( jNonUnique );

                        // get the vector of weights from the unique list
                        Vector< real > &tUniqueWeightVector = aIGtoBSGradientWeightsY( iIgNode )( iUnique );

                        // get the vector of IG nodes from the unique list
                        Vector< moris_id > &tUniqueIGNodeVector = aIGtoBSSensitivityIdsY( iIgNode )( iUnique );

                        // skip if it is first instance, because unique list will be empty and we can just append the non-unique list to it.
                        if ( tUniqueWeightVector.size() == 0 )
                        {
                            tUniqueWeightVector = tWeightVector;
                            tUniqueIGNodeVector = tIGNodeVector;
                        }
                        else
                        {
                            // loop over the non-unique list and see if there are any IG nodes that match with the unique list. If so, add contributions.
                            for ( uint iNonUnique = 0; iNonUnique < tIGNodeVector.size(); iNonUnique++ )
                            {
                                moris_id tNonUniqueIGNode = tIGNodeVector( iNonUnique );
                                real     tNonUniqueWeight = tWeightVector( iNonUnique );

                                // check if this IG node is already in the unique list
                                auto it = std::find( tUniqueIGNodeVector.begin(), tUniqueIGNodeVector.end(), tNonUniqueIGNode );
                                if ( it != tUniqueIGNodeVector.end() )
                                {
                                    // if found, add the weight to the corresponding index in the unique weight vector
                                    uint index = std::distance( tUniqueIGNodeVector.begin(), it );
                                    tUniqueWeightVector( index ) += tNonUniqueWeight;
                                }
                                else
                                {
                                    // if not found, append the IG node and weight to the unique lists
                                    tUniqueIGNodeVector.push_back( tNonUniqueIGNode );
                                    tUniqueWeightVector.push_back( tNonUniqueWeight );
                                }
                            }
                        }
                    }
                }
            }

            // repeat for Z component if 3D
            if ( this->get_spatial_dim() == 3 )
            {
                for ( uint iUnique = 0; iUnique < tNumUniqueIDsZ; iUnique++ )
                {
                    aIGtoBSGradientWeightsZ( iIgNode )( iUnique ).clear();
                    aIGtoBSSensitivityIdsZ( iIgNode )( iUnique ).clear();

                    for ( uint jNonUnique = 0; jNonUnique < tNodalIGtoBSIdsForSensitivityZ.numel(); jNonUnique++ )
                    {
                        if ( aIGtoBSIdsForSensitivityZ( iIgNode )( iUnique ) == tNodalIGtoBSIdsForSensitivityZ( jNonUnique ) )
                        {
                            // get the vector of weights from the non-unique list
                            Vector< real > tWeightVector = tNodalIGtoBSGradientWeightsZ( jNonUnique );

                            // get the vector of IG nodes from the non-unique list
                            Vector< moris_id > tIGNodeVector = tNodalIGtoBSSensitivityIdsZ( jNonUnique );

                            // get the vector of weights from the unique list
                            Vector< real > &tUniqueWeightVector = aIGtoBSGradientWeightsZ( iIgNode )( iUnique );

                            // get the vector of IG nodes from the unique list
                            Vector< moris_id > &tUniqueIGNodeVector = aIGtoBSSensitivityIdsZ( iIgNode )( iUnique );

                            // skip if it is first instance, because unique list will be empty and we can just append the non-unique list to it.
                            if ( tUniqueWeightVector.size() == 0 )
                            {
                                tUniqueWeightVector = tWeightVector;
                                tUniqueIGNodeVector = tIGNodeVector;
                            }
                            else
                            {
                                // loop over the non-unique list and see if there are any IG nodes that match with the unique list. If so, add contributions.
                                for ( uint iNonUnique = 0; iNonUnique < tIGNodeVector.size(); iNonUnique++ )
                                {
                                    moris_id tNonUniqueIGNode = tIGNodeVector( iNonUnique );
                                    real     tNonUniqueWeight = tWeightVector( iNonUnique );

                                    // check if this IG node is already in the unique list
                                    auto it = std::find( tUniqueIGNodeVector.begin(), tUniqueIGNodeVector.end(), tNonUniqueIGNode );
                                    if ( it != tUniqueIGNodeVector.end() )
                                    {
                                        // if found, add the weight to the corresponding index in the unique weight vector
                                        uint index = std::distance( tUniqueIGNodeVector.begin(), it );
                                        tUniqueWeightVector( index ) += tNonUniqueWeight;
                                    }
                                    else
                                    {
                                        // if not found, append the IG node and weight to the unique lists
                                        tUniqueIGNodeVector.push_back( tNonUniqueIGNode );
                                        tUniqueWeightVector.push_back( tNonUniqueWeight );
                                    }
                                }
                            }
                        }
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
                Vector< Vertex * > tIPVertices = tInterpolationCell.get_vertex_pointers();

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
            Vector< Matrix< IdMat > >  &aIGtoBSIds,
            Vector< Matrix< IdMat > >  &aIGtoBSIdsX,
            Vector< Matrix< IdMat > >  &aIGtoBSIdsY,
            Vector< Matrix< IdMat > >  &aIGtoBSIdsZ,
            Vector< Matrix< DDRMat > > &aIGtoBSWeights,
            Vector< Matrix< DDRMat > > &aIGtoBSGradientWeightsX,
            Vector< Matrix< DDRMat > > &aIGtoBSGradientWeightsY,
            Vector< Matrix< DDRMat > > &aIGtoBSGradientWeightsZ,
            Matrix< DDUMat >           &aSparseIndices,
            Matrix< DDUMat >           &aSparseIndicesX,
            Matrix< DDUMat >           &aSparseIndicesY,
            Matrix< DDUMat >           &aSparseIndicesZ,
            Matrix< DDRMat >           &aWeights,
            Matrix< DDRMat >           &aGradientWeightsX,
            Matrix< DDRMat >           &aGradientWeightsY,
            Matrix< DDRMat >           &aGradientWeightsZ,
            std::unordered_map< moris_index, moris_id >     &aIGNodeIDMap )
    {
        // trace this function
        Tracer tTracer( "MTK", "Build Sparse Extraction Operator Matrix" );

        // get number of IG nodes
        uint tNumIgNodes = aIGtoBSIds.size();

        // initialize length of global sparse operator
        uint tNumIdWeightPairs = 0;

        // initialize length of global sparse operator X
        uint tNumIdWeightPairsX = 0;

        // initialize length of global sparse operator Y
        uint tNumIdWeightPairsY = 0;

        // initialize length of global sparse operator Z
        uint tNumIdWeightPairsZ = 0;

        // loop over all IG Nodes to figure out total number of weights
        for ( uint iG = 0; iG < tNumIgNodes; iG++ )
        {
            tNumIdWeightPairs += aIGtoBSIds( iG ).numel();
            tNumIdWeightPairsX += aIGtoBSIdsX( iG ).numel();
            tNumIdWeightPairsY += aIGtoBSIdsY( iG ).numel();
            if (this->get_spatial_dim() == 3 )
            {
                tNumIdWeightPairsZ += aIGtoBSIdsZ( iG ).numel();
            }
        }

        // initialize sparse matrix
        aSparseIndices.set_size( tNumIdWeightPairs, 2 );
        aSparseIndicesX.set_size( tNumIdWeightPairsX, 2 );
        aSparseIndicesY.set_size( tNumIdWeightPairsY, 2 );
        if (this->get_spatial_dim() == 3 )
        {
            aSparseIndicesZ.set_size( tNumIdWeightPairsZ, 2 );
        }
        
        aWeights.set_size( tNumIdWeightPairs, 1 );
        aGradientWeightsX.set_size( tNumIdWeightPairsX, 1 );
        aGradientWeightsY.set_size( tNumIdWeightPairsY, 1 );

        if (this->get_spatial_dim() == 3 )
        {
            aGradientWeightsZ.set_size( tNumIdWeightPairsZ, 1 );
        }   

        // initialize index counter
        uint tIndex = 0;

        // initialize index counter x
        uint tIndexX = 0;
        
        // initialize index counter y
        uint tIndexY = 0;

        // initialize index counter z
        uint tIndexZ = 0;

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
                    //aSparseIndices( tIndex, 0 ) = iG + 1;
                    aSparseIndices( tIndex, 0 ) = aIGNodeIDMap.at( iG );
                    aSparseIndices( tIndex, 1 ) = tBspId;
                    aWeights( tIndex )          = aIGtoBSWeights( iG )( iBsp );
                    //aGradientWeightsX( tIndex ) = aIGtoBSGradientWeightsX( iG )( iBsp );
                    //aGradientWeightsY( tIndex ) = aIGtoBSGradientWeightsY( iG )( iBsp );

                    // increment index
                    tIndex++;
                }
            }
            
            // check if IG node has interpolation x
            if ( aIGtoBSIdsX( iG ).numel() > 0 )
            {
                // loop over B-splines interpolating into current IG node x
                for ( uint iBsp = 0; iBsp < aIGtoBSIdsX( iG ).numel(); iBsp++ )
                {
                    // get
                    uint tBspId = aIGtoBSIdsX( iG )( iBsp );

                    // write IG/BS indices and weights to list
                    //aSparseIndicesX( tIndexX, 0 ) = iG + 1;
                    aSparseIndicesX( tIndexX, 0 ) = aIGNodeIDMap.at( iG );
                    aSparseIndicesX( tIndexX, 1 ) = tBspId;
                    aGradientWeightsX( tIndexX ) = aIGtoBSGradientWeightsX( iG )( iBsp );

                    // increment index
                    tIndexX++;

                }
            }

            // check if IG node has interpolation Y
            if ( aIGtoBSIdsY( iG ).numel() > 0 )
            {
                // loop over B-splines interpolating into current IG node
                for ( uint iBsp = 0; iBsp < aIGtoBSIdsY( iG ).numel(); iBsp++ )
                {
                    // get
                    uint tBspId = aIGtoBSIdsY( iG )( iBsp );

                    // write IG/BS indices and weights to list
                    //aSparseIndicesY( tIndexY, 0 ) = iG + 1;
                    aSparseIndicesY( tIndexY, 0 ) = aIGNodeIDMap.at( iG );
                    aSparseIndicesY( tIndexY, 1 ) = tBspId;
                    aGradientWeightsY( tIndexY ) = aIGtoBSGradientWeightsY( iG )( iBsp );

                    // increment index
                    tIndexY++;
                }
            }

            // check if IG node has interpolation Z
            if ( this->get_spatial_dim() == 3 )
            {
                if ( aIGtoBSIdsZ( iG ).numel() > 0 )
                {
                    // loop over B-splines interpolating into current IG node
                    for ( uint iBsp = 0; iBsp < aIGtoBSIdsZ( iG ).numel(); iBsp++ )
                    {
                        // get
                        uint tBspId = aIGtoBSIdsZ( iG )( iBsp );

                        // write IG/BS indices and weights to list
                        //aSparseIndicesZ( tIndexZ, 0 ) = iG + 1;
                        aSparseIndicesZ( tIndexZ, 0 ) = aIGNodeIDMap.at( iG );
                        aSparseIndicesZ( tIndexZ, 1 ) = tBspId;
                        aGradientWeightsZ( tIndexZ )  = aIGtoBSGradientWeightsZ( iG )( iBsp );

                        // increment index
                        tIndexZ++;
                    }
                }
            }
        }    // end: loop over all IG nodes
    }

    //----------------------------------------------------------------------------

    // new function to incorporate IG-IP sensitivities and IP-BS sensitivities into the IG-BS sensitivities.
    // IP-BS sensitivities due to basis extension.
    void
    Integration_Mesh::build_sparse_extraction_operator(
            Vector< Matrix< IdMat > >                   &aIGtoBSIds,
            Vector< Matrix< IdMat > >                   &aIGtoBSIdsForSensitivitiesX,
            Vector< Matrix< IdMat > >                   &aIGtoBSIdsForSensitivitiesY,
            Vector< Matrix< IdMat > >                   &aIGtoBSIdsForSensitivitiesZ,
            Vector< Vector< Vector< moris_id > > >      &aIGtoBSSensitivityIdsX,
            Vector< Vector< Vector< moris_id > > >      &aIGtoBSSensitivityIdsY,
            Vector< Vector< Vector< moris_id > > >      &aIGtoBSSensitivityIdsZ,
            Vector< Matrix< DDRMat > >                  &aIGtoBSWeights,
            Vector< Vector< Vector< real > > >          &aIGtoBSGradientWeightsX,
            Vector< Vector< Vector< real > > >          &aIGtoBSGradientWeightsY,
            Vector< Vector< Vector< real > > >          &aIGtoBSGradientWeightsZ,
            Matrix< DDUMat >                            &aSparseIndices,
            Matrix< DDUMat >                            &aSparseIndicesX,
            Matrix< DDUMat >                            &aSparseIndicesY,
            Matrix< DDUMat >                            &aSparseIndicesZ,
            Matrix< DDRMat >                            &aWeights,
            Matrix< DDRMat >                            &aGradientWeightsX,
            Matrix< DDRMat >                            &aGradientWeightsY,
            Matrix< DDRMat >                            &aGradientWeightsZ,
            std::unordered_map< moris_index, moris_id > &aIGNodeIDMap )
    {
        // the intention is to store the sensitivities in a format like thus -
        // (IG node ID, B-Spline Node ID, IG sensitivity node ID, weight) 
        // for each of the X, Y and Z components. 
        // The IG sensitivity node ID is the ID of the IG node with respect to which the derivative is being taken. 
        // The B-Spline Node ID is the ID of the B-Spline basis function with respect to which the derivative is being taken. 
        // The weight is the value of the derivative of the B-Spline basis function with respect to the IG node.
        
        // trace this function
        Tracer tTracer( "MTK", "Build Sparse Extraction Operator Matrix" );

        // get number of IG nodes
        uint tNumIgNodes = aIGtoBSIds.size();

        // initialize length of global sparse operator
        uint tNumIdWeightPairs = 0;

        // initialize length of global sparse operator X
        uint tNumIdWeightPairsX = 0;

        // initialize length of global sparse operator Y
        uint tNumIdWeightPairsY = 0;

        // initialize length of global sparse operator Z
        uint tNumIdWeightPairsZ = 0;

        // loop over all IG Nodes to figure out total number of weights
        for ( uint iG = 0; iG < tNumIgNodes; iG++ )
        {
            tNumIdWeightPairs += aIGtoBSIds( iG ).numel();

            // for the sensitivities we have multiple weights corresponding to same B-Spline function. Need to count them all
            //tNumIdWeightPairsX += aIGtoBSIdsForSensitivitiesX( iG ).numel();
            for ( uint iBf = 0; iBf < aIGtoBSIdsForSensitivitiesX( iG ).numel(); iBf++ )
            {
                tNumIdWeightPairsX += aIGtoBSGradientWeightsX( iG )( iBf ).size();
            }

            //tNumIdWeightPairsY += aIGtoBSIdsForSensitivitiesY( iG ).numel();
            for ( uint iBf = 0; iBf < aIGtoBSIdsForSensitivitiesY( iG ).numel(); iBf++ )
            {
                tNumIdWeightPairsY += aIGtoBSGradientWeightsY( iG )( iBf ).size();
            }


            if ( this->get_spatial_dim() == 3 )
            {
                //tNumIdWeightPairsZ += aIGtoBSIdsForSensitivitiesZ( iG ).numel();
                for ( uint iBf = 0; iBf < aIGtoBSIdsForSensitivitiesZ( iG ).numel(); iBf++ )
                {
                    tNumIdWeightPairsZ += aIGtoBSGradientWeightsZ( iG )( iBf ).size();
                }
            }
        }

        // initialize sparse matrix
        aSparseIndices.set_size( tNumIdWeightPairs, 2 );
        aSparseIndicesX.set_size( tNumIdWeightPairsX, 3 );
        aSparseIndicesY.set_size( tNumIdWeightPairsY, 3 );
        if ( this->get_spatial_dim() == 3 )
        {
            aSparseIndicesZ.set_size( tNumIdWeightPairsZ, 3 );
        }

        aWeights.set_size( tNumIdWeightPairs, 1 );
        aGradientWeightsX.set_size( tNumIdWeightPairsX, 1 );
        aGradientWeightsY.set_size( tNumIdWeightPairsY, 1 );

        if ( this->get_spatial_dim() == 3 )
        {
            aGradientWeightsZ.set_size( tNumIdWeightPairsZ, 1 );
        }

        // initialize index counter
        uint tIndex = 0;

        // initialize index counter x
        uint tIndexX = 0;

        // initialize index counter y
        uint tIndexY = 0;

        // initialize index counter z
        uint tIndexZ = 0;

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
                    // aSparseIndices( tIndex, 0 ) = iG + 1;
                    aSparseIndices( tIndex, 0 ) = aIGNodeIDMap.at( iG );
                    aSparseIndices( tIndex, 1 ) = tBspId;
                    aWeights( tIndex )          = aIGtoBSWeights( iG )( iBsp );
                    // aGradientWeightsX( tIndex ) = aIGtoBSGradientWeightsX( iG )( iBsp );
                    // aGradientWeightsY( tIndex ) = aIGtoBSGradientWeightsY( iG )( iBsp );

                    // increment index
                    tIndex++;
                }
            }
            // Populate the derivative contributions to the IG-BS sensitivities. The IG-BS sensitivities are stored in a format like thus -
            // (IG node ID, B-Spline Node ID, IG sensitivity node ID, weight) for each of the X, Y and Z components.  
            // The IG sensitivity node ID is the ID of the IG node with respect to which the derivative is being taken.
            // check if IG node has interpolation x
            if ( aIGtoBSIdsForSensitivitiesX( iG ).numel() > 0 )
            {
                // loop over B-splines interpolating into current IG node x
                for ( uint iBsp = 0; iBsp < aIGtoBSIdsForSensitivitiesX( iG ).numel(); iBsp++ )
                {
                    // get B Spline ID for the current sensitivity.
                    uint tBspId = aIGtoBSIdsForSensitivitiesX( iG )( iBsp );

                    for ( uint iSens = 0; iSens < aIGtoBSGradientWeightsX( iG )( iBsp ).size(); iSens++ )
                    {
                        // write IG/BS indices and weights to list

                        // get node ID of the IG node with respect to which the derivative is being taken.
                        uint tIGSensID = aIGtoBSSensitivityIdsX( iG )( iBsp )( iSens );
                        // aSparseIndicesX( tIndexX, 0 ) = iG + 1;
                        aSparseIndicesX( tIndexX, 0 ) = aIGNodeIDMap.at( iG );
                        aSparseIndicesX( tIndexX, 1 ) = tBspId;
                        aSparseIndicesX( tIndexX, 2 ) = tIGSensID;
                        aGradientWeightsX( tIndexX )  = aIGtoBSGradientWeightsX( iG )( iBsp )( iSens );
                        
                        // increment index
                        tIndexX++;
                    }
                }
            }

            if ( aIGtoBSIdsForSensitivitiesY( iG ).numel() > 0 )
            {
                // loop over B-splines interpolating into current IG node x
                for ( uint iBsp = 0; iBsp < aIGtoBSIdsForSensitivitiesY( iG ).numel(); iBsp++ )
                {
                    // get B Spline ID for the current sensitivity.
                    uint tBspId = aIGtoBSIdsForSensitivitiesY( iG )( iBsp );

                    for ( uint iSens = 0; iSens < aIGtoBSGradientWeightsY( iG )( iBsp ).size(); iSens++ )
                    {
                        // write IG/BS indices and weights to list

                        // get node ID of the IG node with respect to which the derivative is being taken.
                        uint tIGSensID = aIGtoBSSensitivityIdsY( iG )( iBsp )( iSens );
                        // aSparseIndicesX( tIndexX, 0 ) = iG + 1;
                        aSparseIndicesY( tIndexY, 0 ) = aIGNodeIDMap.at( iG );
                        aSparseIndicesY( tIndexY, 1 ) = tBspId;
                        aSparseIndicesY( tIndexY, 2 ) = tIGSensID;
                        aGradientWeightsY( tIndexY )  = aIGtoBSGradientWeightsY( iG )( iBsp )( iSens );

                        // increment index
                        tIndexY++;
                    }
                }
            }

            // check if IG node has interpolation Z
            if ( this->get_spatial_dim() == 3 )
            {
                if ( aIGtoBSIdsForSensitivitiesZ( iG ).numel() > 0 )
                {
                    // loop over B-splines interpolating into current IG node x
                    for ( uint iBsp = 0; iBsp < aIGtoBSIdsForSensitivitiesZ( iG ).numel(); iBsp++ )
                    {
                        // get B Spline ID for the current sensitivity.
                        uint tBspId = aIGtoBSIdsForSensitivitiesZ( iG )( iBsp );

                        for ( uint iSens = 0; iSens < aIGtoBSGradientWeightsZ( iG )( iBsp ).size(); iSens++ )
                        {
                            // write IG/BS indices and weights to list

                            // get node ID of the IG node with respect to which the derivative is being taken.
                            uint tIGSensID = aIGtoBSSensitivityIdsZ( iG )( iBsp )( iSens );
                            // aSparseIndicesX( tIndexX, 0 ) = iG + 1;
                            aSparseIndicesZ( tIndexZ, 0 ) = aIGNodeIDMap.at( iG );
                            aSparseIndicesZ( tIndexZ, 1 ) = tBspId;
                            aSparseIndicesZ( tIndexZ, 2 ) = tIGSensID;
                            aGradientWeightsZ( tIndexZ )  = aIGtoBSGradientWeightsZ( iG )( iBsp )( iSens );

                            // increment index
                            tIndexZ++;
                        }
                    }
                }
            }
        }    // end: loop over all IG nodes
    }

    // ----------------------------------------------------------------------------
}    // namespace moris::mtk

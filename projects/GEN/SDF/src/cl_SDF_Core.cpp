/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Core.cpp
 *
 */

#include <fstream>

#include "cl_Stopwatch.hpp"
#include "cl_Communication_Tools.hpp"
#include "SDF_Tools.hpp"
#include "cl_SDF_Core.hpp"
#include "cl_Tracer.hpp"
#include "fn_sort.hpp"
#include "fn_print.hpp"

namespace moris::sdf
{
    //-------------------------------------------------------------------------------

    Core::Core( Mesh& aMesh, Object& aObject, bool aVerbose )
            : mMesh( aMesh )
            , mObject( aObject )
            , mVerbose( aVerbose )
    {
        // fill unsure nodes list
        // uint tNumberOfNodes = aMesh.get_num_nodes();
    }
    //-------------------------------------------------------------------------------

    void
    Core::calculate_raycast(
            Matrix< IndexMat >& aElementsAtSurface,
            Matrix< IndexMat >& aElementsInVolume )
    {

        // call private routine
        this->raycast_mesh();

        // assign element containers
        aElementsAtSurface.set_size( mSurfaceElements, 1 );

        aElementsInVolume.set_size( mVolumeElements, 1 );

        // counters
        uint tSurfaceCount = 0;
        uint tVolumeCount  = 0;

        // get number of elements
        uint tNumberOfElements = mMesh.get_num_elems();

        // loop over all elements
        for ( uint k = 0; k < tNumberOfElements; ++k )
        {
            // get pointer to element
            Cell* tElement = mMesh.get_cell( k );

            if ( tElement->is_on_surface() )
            {
                aElementsAtSurface( tSurfaceCount++ ) = k;
            }
            else if ( tElement->is_in_volume() )
            {
                aElementsInVolume( tVolumeCount++ ) = k;
            }
        }

        // make sure that everything is OK
        MORIS_ASSERT( tSurfaceCount == mSurfaceElements,
                "Number of surface elements does not match. Surface Count after raycast: %d\t Struct Surface Elements: %d.",
                tSurfaceCount,
                mSurfaceElements );

        MORIS_ASSERT( tVolumeCount == mVolumeElements,
                "Number of volume elements does not match. Volume Count after raycast: %d\t Struct Volume Elements: %d.",
                tVolumeCount,
                mVolumeElements );
    }

    //-------------------------------------------------------------------------------

    void
    Core::calculate_raycast(
            Matrix< IndexMat >& aElementsAtSurface )
    {

        // call private routine
        this->raycast_mesh();

        // assign element containers
        aElementsAtSurface.set_size( mSurfaceElements, 1 );

        // counters
        uint tSurfaceCount = 0;

        // get number of elements
        uint tNumberOfElements = mMesh.get_num_elems();

        // loop over all elements
        for ( uint k = 0; k < tNumberOfElements; ++k )
        {
            // get pointer to element
            Cell* tElement = mMesh.get_cell( k );

            if ( tElement->is_on_surface() )
            {
                aElementsAtSurface( tSurfaceCount++ ) = k;
            }
        }

        // make sure that everything is OK
        MORIS_ASSERT( tSurfaceCount = mSurfaceElements,
                "Number of surface elements does not match" );
    }

    //-------------------------------------------------------------------------------

    void
    Core::raycast_mesh()
    {
        // time this function
        Tracer tTracer( "SDF", "Perform Ray-Cast" );

        // set unsure flag of all nodes to true
        uint tNumberOfNodes = mMesh.get_num_nodes();

        for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfNodes; ++iNodeIndex )
        {
            mMesh.get_vertex( iNodeIndex )->reset();
        }

        for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfNodes; ++iNodeIndex )
        {
            if ( mMesh.get_vertex( iNodeIndex )->is_flagged() )
            {
                // get node coordinate
                const Matrix< DDRMat >& tPoint = mMesh.get_node_coordinate( iNodeIndex );

                // raycast on this point until the point is determined
                Object_Region tPointRegion = raycast_point( mObject, tPoint );

                switch ( tPointRegion )
                {
                    case OUTSIDE:
                    {
                        mMesh.get_vertex( iNodeIndex )->unset_inside_flag();
                        mMesh.get_vertex( iNodeIndex )->unflag();
                        break;
                    }
                    case INSIDE:
                    {
                        mMesh.get_vertex( iNodeIndex )->set_inside_flag();
                        mMesh.get_vertex( iNodeIndex )->unflag();
                        break;
                    }
                    case UNSURE:
                    {
                        mMesh.get_vertex( iNodeIndex )->flag();
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "SDF_Core - raycast_mesh(): Unexpected inside condition returned from Raycast class. Inside condition = %u, should be [0,2]", tPointRegion );
                    }
                }
            }
        }

        this->calculate_candidate_points_and_buffer_diagonal();
    }

    //-------------------------------------------------------------------------------

    void
    Core::calculate_raycast_and_sdf( Matrix< DDRMat >& aSDF )
    {
        this->raycast_mesh();

        Vector< Vertex* > tCandidateList;               //========================================
        tCandidateList = this->set_candidate_list();    //===================================

        this->calculate_udf( tCandidateList );
        this->sweep();
        this->fill_sdf_with_values( aSDF );
    }

    //-------------------------------------------------------------------------------

    void
    Core::calculate_raycast_and_sdf(
            Matrix< DDRMat >&   aSDF,
            Matrix< IndexMat >& aElementsAtSurface,
            Matrix< IndexMat >& aElementsInVolume )
    {
        this->calculate_raycast( aElementsAtSurface, aElementsInVolume );

        Vector< Vertex* > tCandidateList;               //========================================
        tCandidateList = this->set_candidate_list();    //===================================

        this->calculate_udf( tCandidateList );
        this->sweep();
        this->fill_sdf_with_values( aSDF );
    }

    //-------------------------------------------------------------------------------

    void
    Core::calculate_udf( Vector< Vertex* >& aCandidateList )
    {
        Tracer tTracer( "SDF", "Compute UDF" );

        // get number of triangles
        uint tNumberOfFacets = mObject.get_num_facets();

        // loop over all triangles
        for ( uint k = 0; k < tNumberOfFacets; ++k )
        {
            // get pointer to triangle
            Facet& tFacet = mObject.get_facet( k );

            // get nodes within triangle
            Vector< Vertex* > tNodes;

            this->get_nodes_withing_bounding_box_of_triangle(
                    tFacet, tNodes, aCandidateList );

            // get number of nodes
            uint tNumberOfNodes = tNodes.size();

            // calculate distance of this point to the triangle
            // and update udf value if it is smaller
            for ( uint i = 0; i < tNumberOfNodes; ++i )
            {
                // update UDF of this node
                tNodes( i )->update_udf( tFacet );
            }

        }    // end loop over all triangles
    }

    //-------------------------------------------------------------------------------

    void
    Core::calculate_candidate_points_and_buffer_diagonal()
    {
        // get number of elements
        uint tNumberOfElements = mMesh.get_num_elems();

        // counter for elements near surface
        mSurfaceElements = 0;

        // counter for elements in volume
        mVolumeElements = 0;

        // reset buffer diagonal
        mBufferDiagonal = 0;

        // search all elements for sign change
        for ( uint e = 0; e < tNumberOfElements; ++e )
        {
            // unflag this element

            Cell* tElement = mMesh.get_cell( e );

            // reset flags of this element
            tElement->unflag();
            tElement->unset_surface_flag();
            tElement->unset_volume_flag();

            // get pointer to nodes
            const Vector< Vertex* > tNodes = tElement->get_vertices();

            // get number of nodes
            uint tNumberOfNodes = tNodes.size();

            // get first sign
            bool tIsInside = tNodes( 0 )->is_inside();

            // assume element is not intersected
            bool tIsIntersected = false;

            // loop over all other nodes
            for ( uint k = 1; k < tNumberOfNodes; ++k )
            {
                // check of sign is the same
                if ( tNodes( k )->is_inside() != tIsInside )
                {
                    // sign is not same
                    tIsIntersected = true;

                    // cancel loop
                    break;
                }
            }

            // test if there is a sign change
            if ( tIsIntersected )
            {
                // flag this element as surface element
                tElement->set_surface_flag();
                tElement->unset_volume_flag();

                // increment counter
                ++mSurfaceElements;

                // update buffer diagonal
                mBufferDiagonal = std::max(
                        mBufferDiagonal,
                        tElement->get_buffer_diagonal() );

                // flag to indicate that the buffer of this element
                // has been calculated
                tElement->flag();

                // flag all nodes of this element as candidates
                for ( uint k = 0; k < tNumberOfNodes; ++k )
                {
                    tNodes( k )->set_candidate_flag();
                }
            }
            else if ( tIsInside )
            {
                // flag this element as volume element
                tElement->unset_surface_flag();
                tElement->set_volume_flag();

                // increment counter
                ++mVolumeElements;
            }
            else
            {
                // unflag element
                tElement->unset_surface_flag();
                tElement->unset_volume_flag();
            }
        }

        // add additional search depth
        for ( uint d = 1; d < mCandidateSearchDepth; ++d )
        {
            // loop over all elements
            for ( uint e = 0; e < tNumberOfElements; ++e )
            {
                // get pointer to element
                Cell* tElement = mMesh.get_cell( e );

                // test if element is not flagged
                if ( !tElement->is_flagged() )
                {
                    // get pointer to nodes
                    const Vector< Vertex* > tNodes = tElement->get_vertices();

                    // get number of nodes
                    uint tNumberOfNodes = tNodes.size();

                    bool tIsCandidate = false;

                    // test if at least one node of this element is flagged
                    // as candidate
                    for ( uint k = 0; k < tNumberOfNodes; ++k )
                    {
                        if ( tNodes( k )->is_candidate() )
                        {
                            tIsCandidate = true;
                            break;
                        }
                    }

                    // test if candidtae flag is set
                    if ( tIsCandidate )
                    {
                        // update buffer diagonal
                        mBufferDiagonal = std::max(
                                mBufferDiagonal,
                                tElement->get_buffer_diagonal() );

                        // flag this element
                        tElement->flag();

                        // flag all nodes of this element
                        for ( uint k = 0; k < tNumberOfNodes; ++k )
                        {
                            tNodes( k )->set_candidate_flag();
                        }
                    }
                }    // end loop over all elements
            }
        }    // end candidate search depth loop
    }

    //-------------------------------------------------------------------------------

    Vector< Vertex* >
    Core::set_candidate_list()
    {
        uint tNumberOfNodes = mMesh.get_num_nodes();
        //        	std::cout<<"number of nodes in mesh   : "<<tNumberOfNodes<<std::endl;
        Vector< Vertex* > tCandidateVertices;

        for ( uint k = 0; k < tNumberOfNodes; k++ )
        {
            Vertex* tNode = mMesh.get_vertex( k );

            if ( tNode->is_candidate() )
            {
                tCandidateVertices.push_back( tNode );
            }
            else
            {
                continue;
            }
        }
        //        	std::cout<<"number of candidate nodes : "<<tCandidateVertices.size()<<std::endl;
        return tCandidateVertices;
    }

    //-------------------------------------------------------------------------------

    void
    Core::get_nodes_withing_bounding_box_of_triangle(
            Facet&             aFacet,
            Vector< Vertex* >& aNodes,
            Vector< Vertex* >& aCandList )    //===========================================
    {
        // calculate minimum and maximum coordinate

        Matrix< F31RMat > tMinCoord( 3, 1 );
        Matrix< F31RMat > tMaxCoord( 3, 1 );

        for ( uint i = 0; i < 3; ++i )
        {
            tMinCoord( i ) = aFacet.get_min_coord( i ) - mBufferDiagonal;
            tMaxCoord( i ) = aFacet.get_max_coord( i ) + mBufferDiagonal;
        }

        // why is this necessary?

        for ( uint i = 0; i < 3; ++i )
        {
            tMinCoord( i ) = std::max( tMinCoord( i ), mMesh.get_min_coord( i ) );
            tMaxCoord( i ) = std::min( tMaxCoord( i ), mMesh.get_max_coord( i ) );
        }

        //            // get number of nodes on this mesh
        //            uint tNumberOfNodes = mMesh.get_num_nodes();

        // number of candidate nodes
        uint tNumberOfCandidates = aCandList.size();    //========================================

        // node counter
        uint tCount = 0;

        // loop over all nodes of this mesh
        //            for( uint k=0; k<tNumberOfNodes; ++k )
        //            {

        // loop over only the candidate nodes
        for ( uint k = 0; k < tNumberOfCandidates; k++ )    //========================================
        {
            // get pointer to node
            //                Vertex * tNode = mMesh.get_vertex( k );
            Vertex* tNode = aCandList( k );    //=================================================

            // unflag this node
            tNode->unflag();

            // test if node is a candidate
            //                if( tNode->is_candidate() )
            //                {
            // get coords of this node
            const Matrix< F31RMat >& tPoint = tNode->get_coords();

            // assume that node is in triangle
            bool tNodeIsWithinTriangle = true;

            for ( uint i = 0; i < 3; ++i )
            {
                if ( tPoint( i ) < tMinCoord( i ) || tPoint( i ) > tMaxCoord( i ) )
                {
                    tNodeIsWithinTriangle = false;
                    break;
                }
            }

            // if node is in triangle
            if ( tNodeIsWithinTriangle )
            {
                // flag this node
                tNode->flag();

                // increment counter
                ++tCount;
            }
            //                } // end node is candidate
        }    // end loop over all nodes

        // reset output array
        aNodes.resize( tCount, nullptr );

        // reset counter
        tCount = 0;

        // loop over all nodes of this mesh
        //            for( uint k=0; k<tNumberOfNodes; ++k )
        //            {

        // loop over only the candidate nodes+
        for ( uint k = 0; k < tNumberOfCandidates; k++ )    //========================================
        {
            // get pointer to node
            //                Vertex * tNode = mMesh.get_vertex( k );
            Vertex* tNode = aCandList( k );    //==================================================

            // test if node is flagged
            if ( tNode->is_flagged() )
            {
                aNodes( tCount++ ) = tNode;
            }
        }
    }

    //-------------------------------------------------------------------------------

    void
    Core::sweep()
    {
        Tracer tTracer( "SDF", "Sweep" );

        uint tNumberOfVertices = mMesh.get_num_nodes();

        uint tSweepCount = 1;

        while ( tSweepCount != 0 )
        {
            tic tTimer;
            tSweepCount = 0;

            // loop over all vertices
            for ( uint k = 0; k < tNumberOfVertices; ++k )
            {
                // get vertex
                Vertex* tVertex = mMesh.get_vertex( k );

                // test if node has sdf
                if ( tVertex->has_sdf() )
                {
                    // sweep this vertex
                    tSweepCount += tVertex->sweep();
                }
            }

            if ( mVerbose )
            {
                // stop the timer
                real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;

                // print elapsed time
                if ( par_size() == 1 )
                {
                    std::fprintf( stdout, "Time for sweeping              : %5.3f [sec]\nSwept %i nodes\n", tElapsedTime / 1000.0, (int)tSweepCount );
                }
                else
                {
                    std::fprintf( stdout, "Proc % i - Time for sweeping              : %5.3f [sec]\nSwept %i nodes\n", (int)par_rank(), tElapsedTime / 1000.0, (int)tSweepCount );
                }
            }
        }
    }

    // -----------------------------------------------------------------------------

    void
    Core::fill_sdf_with_values( Matrix< DDRMat >& aSDF )
    {
        // get number of vertices
        uint tNumberOfVertices = mMesh.get_num_nodes();

        // min and max value
        real tMinSDF = -1e-12;    // std::numeric_limits<real>::max();
        real tMaxSDF = 1e-12;     // std::numeric_limits<real>::min();

        // allocate matrix
        aSDF.set_size( tNumberOfVertices, 1 );

        // loop over all nodes and write real values
        for ( uint k = 0; k < tNumberOfVertices; ++k )
        {
            // get pointer to vertex
            Vertex* tVertex = mMesh.get_vertex( k );

            // test if vertex has SDF
            if ( tVertex->has_sdf() )
            {
                real tSDF = tVertex->get_sdf();
                // write value
                aSDF( tVertex->get_index() ) = tSDF;

                if ( tVertex->is_inside() )
                {
                    tMinSDF = std::min( tMinSDF, tSDF );
                }
                else
                {
                    tMaxSDF = std::max( tMaxSDF, tSDF );
                }
            }
        }

        // if parallel, synchronize min and max values for SDF
        if ( par_size() > 1 )
        {
            // container for min and max values
            Matrix< DDRMat > tValues;

            // communicate minimal value
            allgather_scalar( tMinSDF, tValues );
            tMinSDF = tValues.min();

            // communicate maximal value
            allgather_scalar( tMaxSDF, tValues );
            tMaxSDF = tValues.max();
        }

        // loop over all nodes and write fake values
        for ( uint k = 0; k < tNumberOfVertices; ++k )
        {
            // get pointer to vertex
            Vertex* tVertex = mMesh.get_vertex( k );

            // test if vertex does not have an SDF
            if ( !tVertex->has_sdf() )
            {
                if ( tVertex->is_inside() )
                {
                    aSDF( tVertex->get_index() ) = tMinSDF;
                }
                else
                {
                    aSDF( tVertex->get_index() ) = tMaxSDF;
                }
            }
        }
    }

    // -----------------------------------------------------------------------------

    void
    Core::save_to_vtk( const std::string& aFilePath )
    {
        // open the file
        std::ofstream tFile( aFilePath, std::ios::binary );

        // containers
        float tFChar = 0;
        int   tIChar = 0;

        Matrix< DDRMat > tSDF;

        this->fill_sdf_with_values( tSDF );

        tFile << "# vtk DataFile Version 3.0" << '\n';
        tFile << "GO BUFFS!" << '\n';
        tFile << "BINARY" << '\n';
        // tFile << "ASCII" << std::endl;
        uint tNumberOfNodes = mMesh.get_num_nodes();

        // write node data
        tFile << "DATASET UNSTRUCTURED_GRID" << '\n';

        tFile << "POINTS " << tNumberOfNodes << " float" << '\n';

        // loop over all nodes
        for ( luint k = 0; k < tNumberOfNodes; ++k )
        {
            // get coordinate from node
            const Matrix< F31RMat >& tCoords = mMesh.get_vertex( k )->get_coords();

            // write coordinates to mesh
            tFChar = swap_byte_endian( (float)tCoords( 0 ) );
            tFile.write( (char*)&tFChar, sizeof( float ) );
            tFChar = swap_byte_endian( (float)tCoords( 1 ) );
            tFile.write( (char*)&tFChar, sizeof( float ) );
            tFChar = swap_byte_endian( (float)tCoords( 2 ) );
            tFile.write( (char*)&tFChar, sizeof( float ) );
            // tFile << tCoords( 0 ) << " " << tCoords( 1 ) << " " << tCoords( 2 ) << std::endl;
        }
        tFile << '\n';

        uint tNumberOfElements = mMesh.get_num_elems();

        // write header for cells
        tFile << "CELLS " << tNumberOfElements << " "
              << 9 * tNumberOfElements << '\n';

        // cell types
        Matrix< DDUMat > tCellTypes( tNumberOfElements, 1 );

        for ( luint k = 0; k < tNumberOfElements; ++k )
        {
            // get pointet to cell
            Cell* tCell = mMesh.get_cell( k );

            // get number of vertices
            uint tNumberOfCellVerts = tCell->get_number_of_vertices();

            // get vertex indices
            Matrix< DDUMat > tIndices( tNumberOfCellVerts, 1 );

            // VTK cell type
            uint tCellType = 0;

            switch ( tNumberOfCellVerts )
            {
                case ( 4 ):    // TET 4
                {
                    tCellType = 10;
                    break;
                }
                case ( 10 ):    // TET 10
                {
                    tCellType = 24;
                    break;
                }
                case ( 8 ):    // HEX 8
                {
                    tCellType = 12;
                    break;
                }
                case ( 27 ):    // HEX27
                {
                    tCellType = 29;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "unknown cell type" );
                }
            }

            // remember cell type
            tCellTypes( k ) = tCellType;

            if ( tCellType == 29 )
            {
                // special case for HEX27
                tIndices( 0 )  = tCell->get_vertex( 0 )->get_index();
                tIndices( 1 )  = tCell->get_vertex( 1 )->get_index();
                tIndices( 2 )  = tCell->get_vertex( 2 )->get_index();
                tIndices( 3 )  = tCell->get_vertex( 3 )->get_index();
                tIndices( 4 )  = tCell->get_vertex( 4 )->get_index();
                tIndices( 5 )  = tCell->get_vertex( 5 )->get_index();
                tIndices( 6 )  = tCell->get_vertex( 6 )->get_index();
                tIndices( 7 )  = tCell->get_vertex( 7 )->get_index();
                tIndices( 8 )  = tCell->get_vertex( 8 )->get_index();
                tIndices( 9 )  = tCell->get_vertex( 9 )->get_index();
                tIndices( 10 ) = tCell->get_vertex( 10 )->get_index();
                tIndices( 11 ) = tCell->get_vertex( 11 )->get_index();
                tIndices( 12 ) = tCell->get_vertex( 16 )->get_index();
                tIndices( 13 ) = tCell->get_vertex( 17 )->get_index();
                tIndices( 14 ) = tCell->get_vertex( 18 )->get_index();
                tIndices( 15 ) = tCell->get_vertex( 19 )->get_index();
                tIndices( 16 ) = tCell->get_vertex( 12 )->get_index();
                tIndices( 17 ) = tCell->get_vertex( 13 )->get_index();
                tIndices( 18 ) = tCell->get_vertex( 14 )->get_index();
                tIndices( 19 ) = tCell->get_vertex( 15 )->get_index();
                tIndices( 20 ) = tCell->get_vertex( 23 )->get_index();
                tIndices( 21 ) = tCell->get_vertex( 24 )->get_index();
                tIndices( 22 ) = tCell->get_vertex( 25 )->get_index();
                tIndices( 23 ) = tCell->get_vertex( 26 )->get_index();
                tIndices( 24 ) = tCell->get_vertex( 21 )->get_index();
                tIndices( 25 ) = tCell->get_vertex( 22 )->get_index();
                tIndices( 26 ) = tCell->get_vertex( 20 )->get_index();
            }
            else
            {
                for ( uint i = 0; i < tNumberOfCellVerts; ++i )
                {
                    tIndices( i ) = tCell->get_vertex( i )->get_index();
                }
            }

            tIChar = swap_byte_endian( (int)tNumberOfCellVerts );
            tFile.write( (char*)&tIChar, sizeof( int ) );

            // write indices to file
            for ( uint i = 0; i < tNumberOfCellVerts; ++i )
            {
                tIChar = swap_byte_endian( (int)tIndices( i ) );
                tFile.write( (char*)&tIChar, sizeof( int ) );
            }
        }
        tFile << '\n';
        // write cell types
        tFile << "CELL_TYPES " << tNumberOfElements << '\n';

        for ( luint k = 0; k < tNumberOfElements; ++k )
        {
            tIChar = swap_byte_endian( tCellTypes( k ) );
            tFile.write( (char*)&tIChar, sizeof( int ) );
        }
        tFile << '\n';

        tFile << "POINT_DATA " << tNumberOfNodes << '\n';

        tFile << "SCALARS SDF float" << '\n';
        tFile << "LOOKUP_TABLE default" << '\n';
        for ( uint k = 0; k < tNumberOfNodes; ++k )
        {

            tFChar = swap_byte_endian( (float)tSDF( k ) );
            tFile.write( (char*)&tFChar, sizeof( float ) );
        }

        tFile << "SCALARS HAS_SDF int" << '\n';
        tFile << "LOOKUP_TABLE default" << '\n';
        for ( uint k = 0; k < tNumberOfNodes; ++k )
        {
            if ( mMesh.get_vertex( k )->has_sdf() )
            {
                tIChar = swap_byte_endian( 1 );
            }
            else
            {
                tIChar = swap_byte_endian( 0 );
            }
            tFile.write( (char*)&tIChar, sizeof( int ) );
        }

        tFile << "SCALARS IS_INSIDE int" << '\n';
        tFile << "LOOKUP_TABLE default" << '\n';
        for ( uint k = 0; k < tNumberOfNodes; ++k )
        {
            if ( mMesh.get_vertex( k )->is_inside() )
            {
                tIChar = swap_byte_endian( 1 );
            }
            else
            {
                tIChar = swap_byte_endian( 0 );
            }
            tFile.write( (char*)&tIChar, sizeof( int ) );
        }

        tFile << "SCALARS VERTEX_INDEX int" << '\n';
        tFile << "LOOKUP_TABLE default" << '\n';
        for ( uint k = 0; k < tNumberOfNodes; ++k )
        {
            tIChar = swap_byte_endian( (int)mMesh.get_vertex( k )->get_index() );
            tFile.write( (char*)&tIChar, sizeof( int ) );
        }

        // close the output file
        tFile.close();
    }
    // -----------------------------------------------------------------------------

    void
    Core::save_unsure_to_vtk( const std::string& aFilePath )
    {
        // open the file
        std::ofstream tFile( aFilePath, std::ios::binary );

        // containers
        float tFChar = 0;
        int   tIChar = 0;

        Matrix< DDRMat > tSDF;

        this->fill_sdf_with_values( tSDF );

        tFile << "# vtk DataFile Version 3.0" << '\n';
        tFile << "GO BUFFS!" << '\n';
        tFile << "BINARY" << '\n';
        // tFile << "ASCII" << std::endl;
        uint tNumberOfNodes = mMesh.get_num_nodes();

        // write node data
        tFile << "DATASET UNSTRUCTURED_GRID" << '\n';

        tFile << "POINTS " << tNumberOfNodes << " float" << '\n';

        // loop over all nodes
        for ( luint k = 0; k < tNumberOfNodes; ++k )
        {
            // get coordinate from node
            const Matrix< F31RMat >& tCoords = mMesh.get_vertex( k )->get_coords();

            // write coordinates to mesh
            tFChar = swap_byte_endian( (float)tCoords( 0 ) );
            tFile.write( (char*)&tFChar, sizeof( float ) );
            tFChar = swap_byte_endian( (float)tCoords( 1 ) );
            tFile.write( (char*)&tFChar, sizeof( float ) );
            tFChar = swap_byte_endian( (float)tCoords( 2 ) );
            tFile.write( (char*)&tFChar, sizeof( float ) );
            // tFile << tCoords( 0 ) << " " << tCoords( 1 ) << " " << tCoords( 2 ) << std::endl;
        }
        tFile << '\n';

        uint tNumberOfElements = mMesh.get_num_elems();

        // write header for cells
        tFile << "CELLS " << tNumberOfElements << " "
              << 9 * tNumberOfElements << '\n';

        // cell types
        Matrix< DDUMat > tCellTypes( tNumberOfElements, 1 );

        for ( luint k = 0; k < tNumberOfElements; ++k )
        {
            // get pointet to cell
            Cell* tCell = mMesh.get_cell( k );

            // get number of vertices
            uint tNumberOfCellVerts = tCell->get_number_of_vertices();

            // get vertex indices
            Matrix< DDUMat > tIndices( tNumberOfCellVerts, 1 );

            // VTK cell type
            uint tCellType = 0;

            switch ( tNumberOfCellVerts )
            {
                case ( 4 ):    // TET 4
                {
                    tCellType = 10;
                    break;
                }
                case ( 10 ):    // TET 10
                {
                    tCellType = 24;
                    break;
                }
                case ( 8 ):    // HEX 8
                {
                    tCellType = 12;
                    break;
                }
                case ( 27 ):    // HEX27
                {
                    tCellType = 29;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "unknown cell type" );
                }
            }

            // remember cell type
            tCellTypes( k ) = tCellType;

            if ( tCellType == 29 )
            {
                // special case for HEX27
                tIndices( 0 )  = tCell->get_vertex( 0 )->get_index();
                tIndices( 1 )  = tCell->get_vertex( 1 )->get_index();
                tIndices( 2 )  = tCell->get_vertex( 2 )->get_index();
                tIndices( 3 )  = tCell->get_vertex( 3 )->get_index();
                tIndices( 4 )  = tCell->get_vertex( 4 )->get_index();
                tIndices( 5 )  = tCell->get_vertex( 5 )->get_index();
                tIndices( 6 )  = tCell->get_vertex( 6 )->get_index();
                tIndices( 7 )  = tCell->get_vertex( 7 )->get_index();
                tIndices( 8 )  = tCell->get_vertex( 8 )->get_index();
                tIndices( 9 )  = tCell->get_vertex( 9 )->get_index();
                tIndices( 10 ) = tCell->get_vertex( 10 )->get_index();
                tIndices( 11 ) = tCell->get_vertex( 11 )->get_index();
                tIndices( 12 ) = tCell->get_vertex( 16 )->get_index();
                tIndices( 13 ) = tCell->get_vertex( 17 )->get_index();
                tIndices( 14 ) = tCell->get_vertex( 18 )->get_index();
                tIndices( 15 ) = tCell->get_vertex( 19 )->get_index();
                tIndices( 16 ) = tCell->get_vertex( 12 )->get_index();
                tIndices( 17 ) = tCell->get_vertex( 13 )->get_index();
                tIndices( 18 ) = tCell->get_vertex( 14 )->get_index();
                tIndices( 19 ) = tCell->get_vertex( 15 )->get_index();
                tIndices( 20 ) = tCell->get_vertex( 23 )->get_index();
                tIndices( 21 ) = tCell->get_vertex( 24 )->get_index();
                tIndices( 22 ) = tCell->get_vertex( 25 )->get_index();
                tIndices( 23 ) = tCell->get_vertex( 26 )->get_index();
                tIndices( 24 ) = tCell->get_vertex( 21 )->get_index();
                tIndices( 25 ) = tCell->get_vertex( 22 )->get_index();
                tIndices( 26 ) = tCell->get_vertex( 20 )->get_index();
            }
            else
            {
                for ( uint i = 0; i < tNumberOfCellVerts; ++i )
                {
                    tIndices( i ) = tCell->get_vertex( i )->get_index();
                }
            }

            tIChar = swap_byte_endian( (int)tNumberOfCellVerts );
            tFile.write( (char*)&tIChar, sizeof( int ) );

            // write indices to file
            for ( uint i = 0; i < tNumberOfCellVerts; ++i )
            {
                tIChar = swap_byte_endian( (int)tIndices( i ) );
                tFile.write( (char*)&tIChar, sizeof( int ) );
                // tFile << " " << tIndices( i );
            }
            // tFile << std::endl;
        }
        tFile << '\n';
        // write cell types
        tFile << "CELL_TYPES " << tNumberOfElements << '\n';

        for ( luint k = 0; k < tNumberOfElements; ++k )
        {
            tIChar = swap_byte_endian( tCellTypes( k ) );
            tFile.write( (char*)&tIChar, sizeof( int ) );
        }
        tFile << '\n';

        tFile << "POINT_DATA " << tNumberOfNodes << '\n';

        tFile << "SCALARS RAYCAST int" << '\n';
        tFile << "LOOKUP_TABLE default" << '\n';
        for ( uint k = 0; k < tNumberOfNodes; ++k )
        {
            // test if vertex is determined
            if ( mMesh.get_vertex( k )->is_flagged() )
            {
                tIChar = swap_byte_endian( 0 );
            }
            else if ( mMesh.get_vertex( k )->is_inside() )
            {
                tIChar = swap_byte_endian( ( -1 ) );
            }
            else
            {
                tIChar = swap_byte_endian( 1 );
            }
            tFile.write( (char*)&tIChar, sizeof( int ) );
        }

        tFile << "SCALARS VERTEX_INDEX int" << '\n';
        tFile << "LOOKUP_TABLE default" << '\n';
        for ( uint k = 0; k < tNumberOfNodes; ++k )
        {
            tIChar = swap_byte_endian( (int)mMesh.get_vertex( k )->get_index() );
            tFile.write( (char*)&tIChar, sizeof( int ) );
        }

        // close the output file
        tFile.close();
    }

    // -----------------------------------------------------------------------------

}    // namespace moris::sdf

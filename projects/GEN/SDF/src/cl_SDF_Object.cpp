/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Object.cpp
 *
 */

#include <cl_SDF_Object.hpp>
#include <fstream>
#include <iostream>

#include "cl_Cell.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "SDF_Tools.hpp"
#include "cl_SDF_Facet_Vertex.hpp"
#include "fn_print.hpp"

namespace moris
{
    namespace sdf
    {
        //-------------------------------------------------------------------------------

        Object::Object( const std::string& aFilePath,
                const Cell< real >&        aOffsets,
                const Cell< real >&        aScale )
                : mNumberOfFacets( 0 )
        {
            MORIS_ERROR( aOffsets.size() > 0, "SDF - Object(): Null offset matrix provided. If no offset is needed, use the default value" );

            // check the file extension
            auto tFileExt = aFilePath.substr( aFilePath.find_last_of( "." ) + 1, aFilePath.length() );

            if ( tFileExt == "obj" )
            {
                this->load_from_object_file( aFilePath, aOffsets, aScale );
            }
            else if ( tFileExt == "stl" )
            {
                this->load_from_stl_file( aFilePath );
            }
            else
            {
                MORIS_ERROR( false, "Object(), file type is not supported" );
            }
            MORIS_ASSERT( mNumberOfFacets == mFacets.size(), "SDF - Object(): number of facets not consistent" );
        }

        //-------------------------------------------------------------------------------

        void
        Object::load_from_object_file( const std::string& aFilePath, const moris::Cell< real >& aOffsets, const moris::Cell< real >& aScale )
        {
            // copy file into buffer
            moris::Cell< std::string > tBuffer;
            this->load_ascii_to_buffer( aFilePath, tBuffer );

            // step 1: count number of dimensions, vertices and facets in file

            // flag to determine dimensionality only once
            mDimension           = 0;
            bool tDimsDetermined = false;

            // reset counter for vertices
            uint tNumberOfVertices = 0;

            // reset counter for facets
            uint tNumberOfFacets = 0;

            // get length of buffer
            uint tBufferLength = tBuffer.size();

            // loop over all lines
            for ( uint k = 0; k < tBufferLength; ++k )
            {
                if ( tBuffer( k ).substr( 0, 2 ) == "v " )
                {
                    ++tNumberOfVertices;
                    // for the first line only, determine the number of columns
                    if ( !tDimsDetermined )
                    {
                        for ( char c : tBuffer( k ) )
                        {
                            if ( c == ' ' or c == '\t' )
                            {
                                mDimension++;
                            }
                        }
                        tDimsDetermined = true;
                    }
                }
                else if ( tBuffer( k ).substr( 0, 2 ) == "f " )
                {
                    ++tNumberOfFacets;
                }
            }

            // step 2: create vertices
            moris::Cell< std::shared_ptr< Facet_Vertex > > tVertices( tNumberOfVertices );

            // reset counter

            uint tCount = 0;
            // loop over all lines
            for ( uint k = 0; k < tBufferLength; ++k )
            {
                if ( tBuffer( k ).substr( 0, 2 ) == "v " )
                {
                    // create matrix with coordinates
                    Matrix< DDRMat > tNodeCoords( mDimension, 1 );

                    float tX[ 3 ];

                    if ( mDimension == 3 )
                    {
                        // read ascii data into coordinates
                        std::sscanf( tBuffer( k ).substr( 2, tBuffer( k ).length() ).c_str(),
                                "%f %f %f",
                                &tX[ 0 ],
                                &tX[ 1 ],
                                &tX[ 2 ] );
                    }
                    else if ( mDimension == 2 )
                    {
                        // read ascii data into coordinates
                        std::sscanf( tBuffer( k ).substr( 2, tBuffer( k ).length() ).c_str(),
                                "%f %f",
                                &tX[ 0 ],
                                &tX[ 1 ] );
                    }
                    else
                    {
                        MORIS_ERROR( false, "SDF Object() - .obj files with %d dimensions not supported.", mDimension );
                    }

                    // test coordinates for highpass
                    for ( uint i = 0; i < mDimension; ++i )
                    {
                        if ( std::abs( tX[ i ] ) > mMeshHighPass )
                        {
                            tNodeCoords( i ) = tX[ i ] * aScale( i ) + aOffsets( i );
                        }
                        else
                        {
                            // use zero value
                            tNodeCoords( i ) = 0.0 + aOffsets( i );
                        }
                    }

                    // create vertex
                    tVertices( tCount ) = std::make_shared< Facet_Vertex >( tCount, tNodeCoords );

                    // increment counter
                    ++tCount;
                }
            }

            // step 3
            // create facets

            // reset counter
            tCount = 0;

            // reserve memory
            mFacets.resize( tNumberOfFacets );

            // temporary one-based Ids for facet nodes 1, 2 and 3

            // loop over all lines
            for ( uint k = 0; k < tBufferLength; ++k )
            {
                if ( tBuffer( k ).substr( 0, 2 ) == "f " )
                {
                    // temporary container for vertices
                    moris::Cell< std::shared_ptr< Facet_Vertex > > tNodes( mDimension, nullptr );
                    Matrix< DDUMat >                        tNodeIndices( 3, 1 );
                    // read facet topology
                    if ( mDimension == 3 )
                    {
                        std::sscanf( tBuffer( k ).substr( 2, tBuffer( k ).length() ).c_str(),
                                "%u %u %u",
                                &tNodeIndices( 0 ),
                                &tNodeIndices( 1 ),
                                &tNodeIndices( 2 ) );
                    }
                    else if ( mDimension == 2 )
                    {
                        std::sscanf( tBuffer( k ).substr( 2, tBuffer( k ).length() ).c_str(),
                                "%u %u",
                                &tNodeIndices( 0 ),
                                &tNodeIndices( 1 ) );
                    }
                    else
                    {
                        MORIS_ERROR( false, "SDF Object() - .obj files with %d dimensions not supported.", mDimension );
                    }

                    // assign vertices with facet
                    for ( uint i = 0; i < mDimension; ++i )
                    {
                        // make sure that file is sane
                        MORIS_ERROR( 0 < tNodeIndices( i ) && tNodeIndices( i ) <= tNumberOfVertices,
                                "Invalid vertex ID in object file" );

                        // copy vertex into cell
                        tNodes( i ) = tVertices( tNodeIndices( i ) - 1 );
                    }

                    // create facet pointer
                    switch ( mDimension )
                    {
                        case 2:
                        {
                            mFacets( tCount ) = std::make_shared< Line >( tCount, tNodes );
                            mNumberOfFacets++;
                            break;
                        }
                        case 3:
                        {
                            mFacets( tCount ) = std::make_shared< Triangle >( tCount, tNodes );
                            mNumberOfFacets++;
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "Facet implementation for %dD facets not implemented yet.", mDimension );
                        }
                    }


                    // increment counter
                    ++tCount;
                }
            }
        }

        //-------------------------------------------------------------------------------
        void
        Object::load_ascii_to_buffer( const std::string& aFilePath,
                moris::Cell< std::string >&              aBuffer )
        {
            // try to open ascii file
            std::ifstream tAsciiFile( aFilePath );
            std::string   tLine;

            // load file into buffer, otherwise throw error
            if ( tAsciiFile )
            {

                // count number of lines
                uint tBufferLength = 0;
                while ( !tAsciiFile.eof() )
                {
                    std::getline( tAsciiFile, tLine );
                    ++tBufferLength;
                }
                tAsciiFile.close();
                tAsciiFile.open( aFilePath );

                // load file into buffer
                aBuffer.reserve( tBufferLength );
                while ( !tAsciiFile.eof() )
                {
                    std::getline( tAsciiFile, tLine );
                    aBuffer.push_back( tLine );
                }
                tAsciiFile.close();
            }
            else
            {
                MORIS_ERROR( false, "Something went wrong while trying to load from %s.", aFilePath.c_str() );
                std::cerr << "Something went wrong while trying to load from " << aFilePath << "." << std::endl;
            }
        }

        //-------------------------------------------------------------------------------

        // note: this routines reads the ascii stl files
        // for binary stl see https://en.wikipedia.org/wiki/STL_(file_format)
        void
        Object::load_from_stl_file( const std::string& aFilePath )
        {
            // copy file into buffer
            moris::Cell< std::string > tBuffer;
            this->load_ascii_to_buffer( aFilePath, tBuffer );

            // get length of buffer
            uint tBufferLength = tBuffer.size();

            // - - - - - - - - - - - - -
            // step 1: count triangles
            // - - - - - - - - - - - - -

            // initialize counter
            uint tTriangleCount = 0;

            // loop over all lines
            for ( uint k = 0; k < tBufferLength; ++k )
            {
                // extract first word from string
                std::string tWord = clean( tBuffer( k ) );

                if ( tWord.substr( 0, 5 ) == "facet" )
                {
                    ++tTriangleCount;
                }
            }

            MORIS_ERROR( tTriangleCount > 0, "Could not find any facets in this file. Maybe not an ASCII STL?" );

            // remember number of triangles
            uint tNumberOfTriangles = tTriangleCount;

            // - - - - - - - - - - - - -
            // step 2: create vertices
            // - - - - - - - - - - - - -

            moris::Cell< std::shared_ptr< Facet_Vertex > > tVertices( 3 * tTriangleCount );

            // initialize vertex counter
            uint tVertexCount = 0;

            // create matrix with coordinates
            Matrix< DDRMat > tNodeCoords( 3, 1 );

            // loop over all lines
            for ( uint k = 0; k < tBufferLength; ++k )
            {
                // extract first word from string
                moris::Cell< std::string > tWords = string_to_words( tBuffer( k ) );

                if ( tWords.size() > 0 )
                {
                    if ( tWords( 0 ) == "vertex" )
                    {
                        // parse words to coords
                        tNodeCoords( 0 ) = stod( tWords( 1 ) );
                        tNodeCoords( 1 ) = stod( tWords( 2 ) );
                        tNodeCoords( 2 ) = stod( tWords( 3 ) );

                        // create vertex
                        tVertices( tVertexCount ) = std::make_shared< Facet_Vertex >( tVertexCount, tNodeCoords );

                        // increment vertex counter
                        ++tVertexCount;
                    }
                }
            }

            // make sure that number of triangles is correct
            MORIS_ERROR( tVertexCount == tVertices.size(), "Number of vertices does not match" );

            // - - - - - - - - - - - - -
            // step 3: create triangles
            // - - - - - - - - - - - - -

            // allocate memory
            mFacets.resize( tNumberOfTriangles );

            // reset triangle counter
            tTriangleCount = 0;

            // temporary container for vertices
            moris::Cell< std::shared_ptr< Facet_Vertex > > tNodes( 3, nullptr );

            // create triangles
            for ( uint iTriangle = 0; iTriangle < tNumberOfTriangles; ++iTriangle )
            {
                tNodes( 0 ) = tVertices( tTriangleCount++ );
                tNodes( 1 ) = tVertices( tTriangleCount++ );
                tNodes( 2 ) = tVertices( tTriangleCount++ );

                // create triangle pointer
                mFacets( iTriangle ) = std::make_shared< Triangle >( iTriangle, tNodes );
                mNumberOfFacets++;
            }
        }

        //-------------------------------------------------------------------------------

        Matrix< IndexMat >
        Object::get_nodes_connected_to_element_loc_inds( moris_index aElementIndex ) const
        {
            // get pointer to facet
            return mFacets( aElementIndex )->get_vertex_inds();
        }

        //-------------------------------------------------------------------------------

        void
        Object::rotate( const Matrix< DDRMat >& aRotationMatrix )
        {
            // go through each facet and rotate its vertices
            for ( uint iFacet = 0; iFacet < mFacets.size(); iFacet++ )
            {
                mFacets( iFacet )->rotate( aRotationMatrix );

                // recompute information about the facet (normal, center, etc.)
                mFacets( iFacet )->update_data();
            }

            // reset the transformed flags for all the vertices
            for ( uint iFacet = 0; iFacet < mFacets.size(); iFacet++ )
            {
                mFacets( iFacet )->reset_vertex_transformed_flags();
            }
        }

        //-------------------------------------------------------------------------------

        void
        Object::scale( const moris::Cell< real >& aScaling )
        {
            MORIS_ASSERT( aScaling.size() == mDimension, "SDF_Object: scale_object() - Scaling factors must be equal to object dimension." );

            // go through each facet and scale its vertices
            for ( uint iFacet = 0; iFacet < mFacets.size(); iFacet++ )
            {
                mFacets( iFacet )->scale( aScaling );

                // recompute information about the facet (normal, center, etc.)
                mFacets( iFacet )->update_data();
            }

            // reset the transformed flags for all the vertices
            for ( uint iFacet = 0; iFacet < mFacets.size(); iFacet++ )
            {
                mFacets( iFacet )->reset_vertex_transformed_flags();
            }
        }

        //-------------------------------------------------------------------------------

        void
        Object::shift( const moris::Cell< real >& aShift )
        {
            MORIS_ASSERT( aShift.size() == mDimension, "SDF_Object::shift_object() - Shift must be equal to object dimension." );

            // go through each facet and shift its vertices
            for ( uint iFacet = 0; iFacet < mFacets.size(); iFacet++ )
            {
                mFacets( iFacet )->shift( aShift );

                // recompute information about the facet (normal, center, etc.)
                mFacets( iFacet )->update_data();
            }

            // reset the transformed flags for all the vertices
            for ( uint iFacet = 0; iFacet < mFacets.size(); iFacet++ )
            {
                mFacets( iFacet )->reset_vertex_transformed_flags();
            }
        }

        //-------------------------------------------------------------------------------


        void
        Object::reset_coordinates()
        {
            // go through each facet and reset its coordinates
            for ( uint iFacet = 0; iFacet < mFacets.size(); iFacet++ )
            {
                mFacets( iFacet )->reset_coordinates();

                // recompute information about the facet (normal, center, etc.)
                mFacets( iFacet )->update_data();
            }
        }

        //-------------------------------------------------------------------------------

        real
        Object::get_facet_min_coord(
                uint aFacetIndex,
                uint aAxis )
        {
            MORIS_ASSERT( aFacetIndex >= 0, "SDF_Object:get_facet_min_coord() - aFacetIndex must be >= 0. Current index: %u", aFacetIndex );
            MORIS_ASSERT( aFacetIndex >= 0, "SDF_Object:get_facet_min_coord() - aAxis must be >= 0. Current index: %u", aAxis );

            return mFacets( aFacetIndex )->get_min_coord( aAxis );
        }

        //-------------------------------------------------------------------------------

        real
        Object::get_facet_max_coord(
                uint aFacetIndex,
                uint aAxis )

        {
            MORIS_ASSERT( aFacetIndex >= 0, "SDF_Object:get_facet_max_coord() - aFacetIndex must be >= 0. Current index: %u", aFacetIndex );
            MORIS_ASSERT( aFacetIndex >= 0, "SDF_Object:get_facet_max_coord() - aAxis must be >= 0. Current index: %u", aAxis );

            return mFacets( aFacetIndex )->get_max_coord( aAxis );
        }

        //-------------------------------------------------------------------------------

    } /* namespace sdf */
} /* namespace moris */

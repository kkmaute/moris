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
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "SDF_Tools.hpp"
#include "cl_SDF_Triangle_Vertex.hpp"
#include "cl_SDF_Triangle.hpp"
#include "fn_print.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        Object::Object ( const std::string & aFilePath,
                         Matrix< DDRMat >    aOffsets)
        {
            if( aOffsets.numel()>0)
            {
                mOffsets = aOffsets;
            }

            // check the file extension
            auto tFileExt = aFilePath.substr(aFilePath.find_last_of(".")+1,aFilePath.length());

            if (tFileExt == "obj")
            {
                this->load_from_object_file( aFilePath );
            }
            else if(tFileExt == "stl")
            {
                this->load_from_stl_file( aFilePath );
            }
            else
            {
                MORIS_ERROR(false, "Object(), file type is not supported");
            }
        }

//-------------------------------------------------------------------------------

        void
        Object::load_from_object_file( const std::string& aFilePath )
        {
            // copy file into buffer
            moris::Cell<std::string> tBuffer;
            this->load_ascii_to_buffer( aFilePath, tBuffer );

            // step 1: count number of vertices and triangles in file

            // reset counter for vertices
            uint tNumberOfVertices = 0;

            // reset counter for triangles
            uint tNumberOfTriangles = 0;

            // get length of buffer
            uint tBufferLength = tBuffer.size();

            // loop over all lines
            for ( uint k=0; k<tBufferLength; ++k )
            {
                if( tBuffer( k ).substr( 0, 2 ) == "v ")
                {
                    ++tNumberOfVertices;
                }
                else if ( tBuffer( k ).substr( 0, 2 ) == "f ")
                {
                    ++tNumberOfTriangles;
                }
            }

            // step 2: create vertices
            mVertices.resize( tNumberOfVertices, nullptr );

            // reset counter

            uint tCount = 0;
            // loop over all lines
            for ( uint k=0; k<tBufferLength; ++k )
            {
                if( tBuffer( k ).substr( 0, 2 ) == "v ")
                {
                    // create matrix with coordinates
                    Matrix< DDRMat > tNodeCoords( 3, 1 );

                    float tX[ 3 ];

                    // read ascii data into coordinates
                    std::sscanf(tBuffer(k).substr(
                            2,
                            tBuffer(k).length()).c_str(),
                            "%f %f %f",
                            &tX[ 0 ],
                            &tX[ 1 ],
                            &tX[ 2 ] );

                    // test coordinates for highpass
                    for( uint i=0; i<3; ++i )
                    {
                        if( std::abs( tX[ i ] ) > mMeshHighPass )
                        {
                            tNodeCoords( i ) = tX[ i ] + mOffsets( i );
                        }
                        else
                        {
                            // use zero value
                            tNodeCoords( i ) = 0.0 + mOffsets( i ) ;
                        }
                    }

                    // create vertex
                    mVertices( tCount ) = new Triangle_Vertex( tCount, tNodeCoords );

                    // increment counter
                    ++tCount;
                }
            }

            // step 3
            // create triangles

            // reset counter
            tCount = 0;

            // reserve memory
            mTriangles.resize( tNumberOfTriangles, nullptr );

            // temporary one-based Ids for triangle nodes 1, 2 and 3

            // loop over all lines
            for ( uint k=0; k<tBufferLength; ++k )
            {
                if( tBuffer( k ).substr( 0, 2 ) == "f " )
                {
                    // temporary container for vertices
                    Cell< Triangle_Vertex * > tNodes( 3, nullptr );
                    Matrix< DDUMat > tNodeIndices( 3, 1 );
                    // read triangle topology
                    std::sscanf(tBuffer(k).substr(2,tBuffer(k).length()).c_str(),
                            "%u %u %u",
                            &tNodeIndices( 0 ),
                            &tNodeIndices( 1 ),
                            &tNodeIndices( 2 ));

                    // assign vertices with triangle
                    for( uint i=0; i<3; ++i )
                    {
                        // make sure that file is sane
                        MORIS_ERROR( 0< tNodeIndices( i ) && tNodeIndices( i ) <= tNumberOfVertices,
                                "Invalid vertex ID in object file" );

                        // copy vertex into cell
                        tNodes( i ) = mVertices( tNodeIndices( i ) - 1 );
                    }

                    // create triangle pointer
                    mTriangles( tCount ) = new Triangle( tCount, tNodes );

                    // increment counter
                    ++tCount;
                }
            }
        }

//-------------------------------------------------------------------------------

        Object::~Object()
        {
            for( auto tTriangle : mTriangles )
            {
                delete tTriangle;
            }

            for( auto tVertex : mVertices )
            {
                delete tVertex;
            }
        }

//-------------------------------------------------------------------------------
        void
        Object::load_ascii_to_buffer( const std::string& aFilePath,
                moris::Cell<std::string>& aBuffer)
        {
            // try to open ascii file
            std::ifstream tAsciiFile( aFilePath );
            std::string tLine;

            // load file into buffer, otherwise throw error
            if( tAsciiFile )
            {

                // count number of lines
                uint tBufferLength = 0;
                while ( !tAsciiFile.eof() )
                {
                    std::getline(tAsciiFile, tLine);
                    ++tBufferLength;
                }
                tAsciiFile.close();
                tAsciiFile.open(aFilePath);

                // load file into buffer
                aBuffer.reserve(tBufferLength);
                while (!tAsciiFile.eof())
                {
                    std::getline(tAsciiFile, tLine);
                    aBuffer.push_back(tLine);
                }
                tAsciiFile.close();
            }
            else
            {
                std::cerr << "Something went wrong while trying to load from " <<
                        aFilePath << "." << std::endl;
            }
        }

//-------------------------------------------------------------------------------

        // note: this routines reads the ascii stl files
        // for binary stl see https://en.wikipedia.org/wiki/STL_(file_format)
        void
        Object::load_from_stl_file( const std::string& aFilePath )
        {
            // copy file into buffer
            moris::Cell<std::string> tBuffer;
            this->load_ascii_to_buffer( aFilePath, tBuffer );

            // get length of buffer
            uint tBufferLength = tBuffer.size();

            // - - - - - - - - - - - - -
            // step 1: count triangles
            // - - - - - - - - - - - - -

            // initialize counter
            uint tCount = 0;

            // loop over all lines
            for( uint k=0; k<tBufferLength; ++k )
            {
                // extract first word from string
                std::string tWord = clean( tBuffer ( k ) );

                if( tWord.substr( 0, 5 ) == "facet" )
                {
                    ++tCount;
                }
            }

            MORIS_ERROR( tCount > 0, "Could not find any facets in this file. Maybe not an ASCII STL?" );

            // remember number of triangles
            uint tNumberOfTriangles = tCount;

            // - - - - - - - - - - - - -
            // step 2: create vertices
            // - - - - - - - - - - - - -

            mVertices.resize( 3*tCount, nullptr );

            // reset counter
            tCount = 0;

            // create matrix with coordinates
            Matrix< DDRMat > tNodeCoords( 3, 1 );

            // loop over all lines
            for( uint k=0; k<tBufferLength; ++k )
            {
                // extract first word from string
                Cell< std::string > tWords = string_to_words( tBuffer ( k ) );

                if( tWords.size() > 0 )
                {
                    if( tWords( 0 ) == "vertex" )
                    {
                        // parse words to coords
                        tNodeCoords( 0 ) = stod( tWords( 1 ) );
                        tNodeCoords( 1 ) = stod( tWords( 2 ) );
                        tNodeCoords( 2 ) = stod( tWords( 3 ) );

                        // create vertex
                       mVertices( tCount ) = new Triangle_Vertex( tCount, tNodeCoords );

                       // increment vertex counter
                        ++tCount;
                    }
                }
            }

            // make sure that number of triangles is correct
            MORIS_ERROR( tCount == mVertices.size(), "Number of vertices does not match" );

            // - - - - - - - - - - - - -
            // step 3: create triangles
            // - - - - - - - - - - - - -

            // allocate memory
            mTriangles.resize( tNumberOfTriangles, nullptr );

            // reset counter
            tCount = 0;

            // temporary container for vertices
            Cell< Triangle_Vertex * > tNodes( 3, nullptr );

            // create triangles
            for( uint k=0; k<tNumberOfTriangles; ++k )
            {
                tNodes( 0 ) = mVertices( tCount++ );
                tNodes( 1 ) = mVertices( tCount++ );
                tNodes( 2 ) = mVertices( tCount++ );

                // create triangle pointer
                mTriangles( k ) = new Triangle( k, tNodes );
            }
        }

//-------------------------------------------------------------------------------

        Matrix< IndexMat >
        Object::get_nodes_connected_to_element_loc_inds
                        ( moris_index aElementIndex ) const
        {
            // get pointer to triangle
            return mTriangles( aElementIndex )->get_vertex_inds();
        }

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */


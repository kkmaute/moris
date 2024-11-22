/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_MTK_Load_External_Surface_Mesh.cpp
 *
 */

#include "fn_MTK_Load_External_Surface_Mesh.hpp"

namespace moris::mtk
{
    void
    load_ascii_to_buffer( const std::string& aFilePath,
            Vector< std::string >&           aBuffer )
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

    Matrix< DDRMat >
    load_vertices_from_object_file( const std::string& aFilePath, const Vector< real >& aOffsets, const Vector< real >& aScale )
    {
        // copy file into buffer
        Vector< std::string > tBuffer;
        load_ascii_to_buffer( aFilePath, tBuffer );

        // step 1: count number of dimensions, vertices and facets in file

        // flag to determine dimensionality only once
        uint tDimension      = 0;
        bool tDimsDetermined = false;

        // reset counter for vertices
        uint tNumberOfVertices = 0;

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
                    for ( uint iCharacter = 0; iCharacter < tBuffer( k ).length(); iCharacter++ )
                    {
                        if ( ( tBuffer( k )[ iCharacter ] == ' ' or tBuffer( k )[ iCharacter ] == '\t' ) && iCharacter != tBuffer( k ).length() - 1 )
                        {
                            tDimension++;
                        }
                    }
                    tDimsDetermined = true;
                }
            }
        }

        // intialize output matrix
        Matrix< DDRMat > tVertexCoordinates( tDimension, tNumberOfVertices );

        // reset counter
        uint tVertexCount = 0;

        // loop over all lines
        for ( uint k = 0; k < tBufferLength; ++k )
        {
            if ( tBuffer( k ).substr( 0, 2 ) == "v " )
            {
                real tX[ 3 ];

#if MORIS_USE_32BIT
                // read ascii data into coordinates
                std::sscanf( tBuffer( k ).substr( 2, tBuffer( k ).length() ).c_str(),
                        "%lf %lf %lf",
                        &tX[ 0 ],
                        &tX[ 1 ],
                        &tX[ 2 ] );
#else
                // read ascii data into coordinates
                std::sscanf( tBuffer( k ).substr( 2, tBuffer( k ).length() ).c_str(),
                        "%Lf %Lf %Lf",
                        &tX[ 0 ],
                        &tX[ 1 ],
                        &tX[ 2 ] );
#endif

                // test coordinates for highpass
                for ( uint i = 0; i < tDimension; ++i )
                {
                    if ( std::abs( tX[ i ] ) > 1e-9 )
                    {
                        tVertexCoordinates( i, tVertexCount ) = tX[ i ] * aScale( i ) + aOffsets( i );
                    }
                    else
                    {
                        // use zero value
                        tVertexCoordinates( i, tVertexCount ) = 0.0 + aOffsets( i );
                    }
                }
                // increment counter
                ++tVertexCount;
            }
        }

        return tVertexCoordinates;
    }

    //-------------------------------------------------------------------------------

    Vector< Vector< moris_index > >
    load_facets_from_object_file( const std::string& aFilePath )
    {
        // initialize output vector
        Vector< Vector< moris_index > > tFacetConnectivity;

        // copy file into buffer
        Vector< std::string > tBuffer;
        load_ascii_to_buffer( aFilePath, tBuffer );

        // flag to determine dimensionality only once
        uint tDimension      = 0;
        bool tDimsDetermined = false;

        // get length of buffer
        uint tBufferLength = tBuffer.size();

        // loop over all lines
        for ( uint k = 0; k < tBufferLength; ++k )
        {
            if ( tBuffer( k ).substr( 0, 2 ) == "f " )
            {
                if ( !tDimsDetermined )
                {
                    for ( uint iCharacter = 0; iCharacter < tBuffer( k ).length(); iCharacter++ )
                    {
                        if ( ( tBuffer( k )[ iCharacter ] == ' ' or tBuffer( k )[ iCharacter ] == '\t' ) && iCharacter != tBuffer( k ).length() - 1 )
                        {
                            tDimension++;
                        }
                    }
                    tDimsDetermined = true;
                }

                // temporary container for vertices
                Vector< moris_index > tNodeIndices( tDimension );

                // get the facet indices excluding any texture and normal vector information. Facets are formatted as "f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3"
                // note: the indices are one-based
                moris::size_t tStart = 2;
                moris::size_t tEnd   = tBuffer( k ).find( "/" ) == std::string::npos ? tBuffer( k ).find( " ", tStart ) : tBuffer( k ).find( "/", tStart );

                for ( uint iVertex = 0; iVertex < tDimension; iVertex++ )
                {
                    // Parse the v1 part of vertex (subtract 1 as object files are 1 indexed)
                    tNodeIndices( iVertex ) = std::stoi( tBuffer( k ).substr( tStart, tEnd - tStart ) ) - 1;

                    // Move indices to the next vertex
                    tStart = tBuffer( k ).find( " ", tEnd );
                    if ( tStart == std::string::npos )
                    {
                        break;
                    }
                    tStart++;
                    tEnd = tBuffer( k ).find( "/", tStart ) == std::string::npos ? tBuffer( k ).find( " ", tStart ) : tBuffer( k ).find( "/", tStart );
                }

                // add the facet to the list
                tFacetConnectivity.push_back( tNodeIndices );
            }
        }

        return tFacetConnectivity;
    }
}    // namespace moris::mtk

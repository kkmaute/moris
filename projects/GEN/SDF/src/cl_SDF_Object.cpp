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
#include <filesystem>

#include "cl_Vector.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_inv.hpp"
#include "fn_eye.hpp"
#include "op_elemwise_div.hpp"

#include "SDF_Tools.hpp"
#include "cl_SDF_Facet_Vertex.hpp"
#include "fn_print.hpp"

namespace moris::sdf
{
    //-------------------------------------------------------------------------------

    Object::Object( const std::string& aFilePath,
            real                       aIntersectionTolerance,
            const Vector< real >&      aOffsets,
            const Vector< real >&      aScale )
            : mIntersectionTolerance( aIntersectionTolerance )
    {
        MORIS_ERROR( aOffsets.size() > 0, "SDF - Object(): Null offset matrix provided. If no offset is needed, use the default value" );

        // check the file extension
        std::string tFileExt = aFilePath.substr( aFilePath.find_last_of( "." ) + 1, aFilePath.length() );

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
            MORIS_ERROR( false, "sdf::Object(), file type of %s is not supported", tFileExt.c_str() );
        }

        mRotation = eye( mDimension, mDimension );
        mShift.set_size( mDimension, 1, 0.0 );
        mScale.set_size( mDimension, 1, 1.0 );
    }

    //-------------------------------------------------------------------------------

    void
    Object::load_from_object_file( const std::string& aFilePath, const Vector< real >& aOffsets, const Vector< real >& aScale )
    {
        // copy file into buffer
        Vector< std::string > tBuffer;
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
                    for ( uint iCharacter = 0; iCharacter < tBuffer( k ).length(); iCharacter++ )
                    {
                        if ( ( tBuffer( k )[ iCharacter ] == ' ' or tBuffer( k )[ iCharacter ] == '\t' ) && iCharacter != tBuffer( k ).length() - 1 )
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
        mVertices.resize( tNumberOfVertices );

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
                mVertices( tCount ) = std::make_shared< Facet_Vertex >( tCount, tNodeCoords );

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
                Vector< std::shared_ptr< Facet_Vertex > > tNodes( mDimension, nullptr );
                Vector< uint >                            tNodeIndices( mDimension );

                // get the facet indices excluding any texture and normal vector information. Facets are formatted as "f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3"
                // note: the indices are one-based
                moris::size_t tStart = 2;
                moris::size_t tEnd   = tBuffer( k ).find( "/" ) == std::string::npos ? tBuffer( k ).find( " ", tStart ) : tBuffer( k ).find( "/", tStart );

                for ( uint iVertex = 0; iVertex < mDimension; iVertex++ )
                {
                    // Parse the v1 part of vertex
                    tNodeIndices( iVertex ) = std::stoi( tBuffer( k ).substr( tStart, tEnd - tStart ) );

                    // Move indices to the next vertex
                    tStart = tBuffer( k ).find( " ", tEnd );
                    if ( tStart == std::string::npos )
                    {
                        break;
                    }
                    tStart++;
                    tEnd = tBuffer( k ).find( "/", tStart ) == std::string::npos ? tBuffer( k ).find( " ", tStart ) : tBuffer( k ).find( "/", tStart );
                }

                // assign vertices with facet
                for ( uint i = 0; i < mDimension; ++i )
                {
                    // make sure that file is sane
                    MORIS_ERROR( 0 < tNodeIndices( i ) && tNodeIndices( i ) <= tNumberOfVertices,
                            "Invalid vertex ID in object file" );

                    // copy vertex into cell
                    tNodes( i ) = mVertices( tNodeIndices( i ) - 1 );
                }

                // create facet pointer
                switch ( mDimension )
                {
                    case 2:
                    {
                        // FIXME: The intersection tolerance is the same for every facet, therefore it should be necessary member data in Line
                        mFacets( tCount ) = std::make_shared< Line >( tCount, tNodes, mIntersectionTolerance );
                        break;
                    }
                    case 3:
                    {
                        // FIXME: The intersection tolerance is the same for every facet, therefore it should be necessary member data in Line
                        mFacets( tCount ) = std::make_shared< Triangle >( tCount, tNodes, mIntersectionTolerance );
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
            Vector< std::string >&                   aBuffer )
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
        Vector< std::string > tBuffer;
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

        mVertices.resize( 3 * tTriangleCount );

        // initialize vertex counter
        uint tVertexCount = 0;

        // create matrix with coordinates
        Matrix< DDRMat > tNodeCoords( 3, 1 );

        // loop over all lines
        for ( uint k = 0; k < tBufferLength; ++k )
        {
            // extract first word from string
            Vector< std::string > tWords = string_to_words( tBuffer( k ) );

            if ( tWords.size() > 0 )
            {
                if ( tWords( 0 ) == "vertex" )
                {
                    // parse words to coords
                    tNodeCoords( 0 ) = stod( tWords( 1 ) );
                    tNodeCoords( 1 ) = stod( tWords( 2 ) );
                    tNodeCoords( 2 ) = stod( tWords( 3 ) );

                    // create vertex
                    mVertices( tVertexCount ) = std::make_shared< Facet_Vertex >( tVertexCount, tNodeCoords );

                    // increment vertex counter
                    ++tVertexCount;
                }
            }
        }

        // make sure that number of triangles is correct
        MORIS_ERROR( tVertexCount == mVertices.size(), "Number of vertices does not match" );

        // - - - - - - - - - - - - -
        // step 3: create triangles
        // - - - - - - - - - - - - -

        // allocate memory
        mFacets.resize( tNumberOfTriangles );

        // reset triangle counter
        tTriangleCount = 0;

        // temporary container for vertices
        Vector< std::shared_ptr< Facet_Vertex > > tNodes( 3, nullptr );

        // create triangles
        for ( uint iTriangle = 0; iTriangle < tNumberOfTriangles; ++iTriangle )
        {
            tNodes( 0 ) = mVertices( tTriangleCount++ );
            tNodes( 1 ) = mVertices( tTriangleCount++ );
            tNodes( 2 ) = mVertices( tTriangleCount++ );

            // create triangle pointer
            mFacets( iTriangle ) = std::make_shared< Triangle >( iTriangle, tNodes );
        }
    }

    //-------------------------------------------------------------------------------

    void
    Object::write_to_file( std::string aFilePath )
    {
        // Open file for writing
        std::ofstream tFile;
        tFile.open( aFilePath );

        // Write vertices
        for ( auto iVertex : mVertices )
        {
            tFile << "v ";
            for ( uint iDimension = 0; iDimension < mDimension; iDimension++ )
            {
                tFile << iVertex->get_coord( iDimension ) << " ";
            }
            tFile << std::endl;
        }

        // Write facets
        for ( auto iFacet : mFacets )
        {
            tFile << "f ";
            Matrix< IdMat > tIndices = iFacet->get_vertex_inds();
            for ( uint iDimension = 0; iDimension < mDimension; iDimension++ )
            {
                tFile << tIndices( iDimension ) + 1 << " ";
            }
            tFile << std::endl;
        }

        // close file
        tFile.close();
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
        // rotate each vertex
        for ( uint iVertexIndex = 0; iVertexIndex < mVertices.size(); iVertexIndex++ )
        {
            mVertices( iVertexIndex )->rotate_node_coords( aRotationMatrix );
        }

        // update the object's rotation matrix to include the new rotation
        mRotation = aRotationMatrix * mRotation;

        // update the facet normals/centers/hesse/etc.
        this->update_all_facets();
    }

    //-------------------------------------------------------------------------------

    void
    Object::scale( const Vector< real >& aScaling )
    {
        // scale each facet vertex
        for ( uint iVertexIndex = 0; iVertexIndex < mVertices.size(); iVertexIndex++ )
        {
            mVertices( iVertexIndex )->scale_node_coords( aScaling );
        }

        // update the object's scale vector by multiplying the new scale
        std::transform( mScale.begin(), mScale.end(), aScaling.cbegin(), mScale.begin(), std::multiplies< real >() );

        // update the facet normals/centers/hesse/etc.
        this->update_all_facets();
    }

    //-------------------------------------------------------------------------------

    void
    Object::shift( const Vector< real >& aShift )
    {
        // move each facet vertex
        for ( uint iVertexIndex = 0; iVertexIndex < mVertices.size(); iVertexIndex++ )
        {
            mVertices( iVertexIndex )->shift_node_coords_from_current( aShift );
        }

        // update the object's total shift vector by adding the new shift
        std::transform( mShift.begin(), mShift.end(), aShift.cbegin(), mShift.begin(), std::plus< real >() );

        // update the facet normals/centers/hesse/etc.
        this->update_all_facets();
    }

    //-------------------------------------------------------------------------------

    void Object::reset_coordinates()
    {
        // Get the inverse of the rotation matrix
        Matrix< DDRMat > tInverseRotation = inv( mRotation );

        for ( uint iVertex = 0; iVertex < mVertices.size(); iVertex++ )
        {
            mVertices( iVertex )->set_node_coords( tInverseRotation * ( mVertices( iVertex )->get_coords_reference() - mShift ) / mScale );
        }

        // reset the object's rotation matrix
        mRotation = eye( mDimension, mDimension );

        // reset the object's shift
        std::fill( mShift.begin(), mShift.end(), 0.0 );

        // reset the object's scale
        std::fill( mScale.begin(), mScale.end(), 1.0 );

        // update the facet normals/centers/hesse/etc.
        this->update_all_facets();
    }

    //-------------------------------------------------------------------------------

    void
    Object::update_all_facets()
    {
        // update information about the facets (normal, center, etc.)
        for ( uint iFacet = 0; iFacet < mFacets.size(); iFacet++ )
        {
            mFacets( iFacet )->update_data();
        }
    }

    //-------------------------------------------------------------------------------

    real
    Object::get_facet_min_coord(
            uint aFacetIndex,
            uint aAxis )
    {
        return mFacets( aFacetIndex )->get_min_coord( aAxis );
    }

    //-------------------------------------------------------------------------------

    real
    Object::get_facet_max_coord(
            uint aFacetIndex,
            uint aAxis )

    {
        return mFacets( aFacetIndex )->get_max_coord( aAxis );
    }

    //-------------------------------------------------------------------------------

}    // namespace moris::sdf

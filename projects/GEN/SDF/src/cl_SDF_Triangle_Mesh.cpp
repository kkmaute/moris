
#include <GEN/SDF/src/cl_SDF_Triangle_Mesh.hpp>
#include <fstream>
#include <iostream>

#include "cl_Cell.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_SDF_Triangle_Vertex.hpp"
#include "cl_SDF_Triangle.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        Triangle_Mesh::Triangle_Mesh ( const std::string & aFilePath )
        {
            this->load_from_object_file( aFilePath );
        }

//-------------------------------------------------------------------------------

        void
        Triangle_Mesh::load_from_object_file( const std::string& aFilePath )
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
                            tNodeCoords( i ) = tX[ i ];
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
            uint tNode0 = 0;
            uint tNode1 = 0;
            uint tNode2 = 0;

            // loop over all lines
            for ( uint k=0; k<tBufferLength; ++k )
            {
                if( tBuffer( k ).substr( 0, 2 ) == "f " )
                {
                    // temporary container for vertices
                    Cell< Triangle_Vertex * > tNodes( 3, nullptr );

                    // read triangle topology
                    std::sscanf(tBuffer(k).substr(2,tBuffer(k).length()).c_str(),
                            "%u %u %u", &tNode0, &tNode1, &tNode2);

                    // populate container
                    tNodes( 0 ) = mVertices( tNode0 - 1 );
                    tNodes( 1 ) = mVertices( tNode1 - 1 );
                    tNodes( 2 ) = mVertices( tNode2 - 1 );

                    // create triangle pointer
                    mTriangles( tCount ) = new Triangle( tCount, tNodes );

                    // increment counter
                    ++tCount;
                }
            }
            std::cout << " v " << tNumberOfVertices << " f " << tNumberOfTriangles  << std::endl;
        }

//-------------------------------------------------------------------------------

        Triangle_Mesh::~Triangle_Mesh()
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
        Triangle_Mesh::load_ascii_to_buffer( const std::string& aFilePath,
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

        Matrix< IndexMat >
        Triangle_Mesh::get_nodes_connected_to_element_loc_inds
                        ( moris_index aElementIndex ) const
        {
            // get pointer to triangle
            return mTriangles( aElementIndex )->get_vertex_inds();

        }

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

#include "cl_Stopwatch.hpp"
#include "cl_Communication_Tools.hpp"

#include "fn_unique.hpp"

#include "cl_MTK_Enums.hpp"
#include "cl_SDF_Mesh.hpp"


namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        Mesh::Mesh( mtk::Mesh * aMesh, bool aVerbose ) :
                mMesh( aMesh ),
                mVerbose( aVerbose ),
                mMinCoord( 3, 1 ),
                mMaxCoord( 3, 1 )
        {
            // determine interpolation order of mesh
            // pick first element
            Matrix< IndexMat > tNodeIndices
                = aMesh->get_nodes_connected_to_element_loc_inds( 0 );
            // determine interpolation order of mesh
            switch( tNodeIndices.length() )
            {
                case( 4 ) : // tet4
                case( 6 ) : // penta 6
                case( 8 ) : // hex8
                {
                    mOrder = 1;
                    break;
                }
                case( 10 ) : // tet10
                case( 20 ) : // hex20
                case( 27 ) : // hex27
                {
                    mOrder = 2;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "Can't determine order of 3D cell");
                }
            }


            // get number of nodes
            uint tNumberOfNodes = aMesh->get_num_nodes();

            // reserve memory
            mVertices.resize( tNumberOfNodes, nullptr );

            // populate container
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                mVertices( k ) = new Vertex( k, aMesh->get_node_coordinate( k ) );
            }

            // get number of elements
            uint tNumberOfElements = aMesh->get_num_elems();

            // reserve memory
            mCells.resize( tNumberOfElements, nullptr );

            // populate container
            for( uint k=0; k<tNumberOfElements; ++k )
            {
                // get cell indices from mesh

                mCells( k ) = new Cell(
                        k,
                        aMesh->get_glb_entity_id_from_entity_loc_index( k, EntityRank::ELEMENT),
                        aMesh->get_nodes_connected_to_element_loc_inds( k ),
                        mVertices );
            }

            this->link_vertex_cells();
            this->link_vertex_neighbors();

            // find min and max coordinate
            for( uint i=0; i<3; ++i )
            {
                mMinCoord( i ) = std::numeric_limits<real>::max();
            }
            for( uint i=0; i<3; ++i )
            {
                mMaxCoord( i ) = std::numeric_limits<real>::min();
            }

            for( Vertex * tVertex : mVertices )
            {
                Matrix< F31RMat > tPoint = tVertex->get_coords();
                for( uint i=0; i<3; ++i )
                {
                    mMinCoord( i ) = std::min( tPoint( i ), mMinCoord( i ) );
                    mMaxCoord( i ) = std::max( tPoint( i ), mMaxCoord( i ) );
                }

            }

            for( uint i=0; i<3; ++i )
            {
                real tMin = mMinCoord( i );
                real tMax = mMaxCoord( i );

                real tMinGlobCoord = min_all( tMin );
                real tMaxGlobCoord = max_all( tMax );

                mMinCoord( i ) = tMinGlobCoord;
                mMaxCoord( i ) = tMaxGlobCoord;
            }

            // matrix that contains node IDs
            mNodeIDs.set_size( tNumberOfNodes, 1 );
            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                mNodeIDs( k )
                        = aMesh->get_glb_entity_id_from_entity_loc_index( k, EntityRank::NODE );
            }
        }

//-------------------------------------------------------------------------------

        Mesh::Mesh( std::shared_ptr< mtk::Mesh > aMesh, bool aVerbose ) : Mesh ( aMesh.get(), aVerbose )
        {

        }

//-------------------------------------------------------------------------------

        Mesh::~Mesh()
        {
            // delete cell pointers
            for( Cell * tCell : mCells )
            {
                delete tCell;
            }

            // delete vertex pointers
            for( Vertex * tVertex : mVertices )
            {
                delete tVertex;
            }
        }

//-------------------------------------------------------------------------------

        void
        Mesh::link_vertex_cells()
        {
            // count elements per node
            for( Cell * tCell : mCells )
            {
                // get vertex pointers
                moris::Cell< Vertex* > tVertices = tCell->get_vertices();

                // count elements connected to this nodes
                for( Vertex * tVertex : tVertices )
                {
                    tVertex->increment_cell_counter();
                }
            }

            // initialize cell containers
            for( Vertex * tVertex : mVertices )
            {
                tVertex->init_cell_container();
            }

            // insert cells
            for( Cell * tCell : mCells )
            {
                // get vertex pointers
                moris::Cell< Vertex* > tVertices = tCell->get_vertices();

                // count elements connected to this nodes
                for( Vertex * tVertex : tVertices )
                {
                    tVertex->insert_cell( tCell );
                }
            }
        }

//-------------------------------------------------------------------------------

        void
        Mesh::link_vertex_neighbors()
        {
            // loop over all vertices
            for( Vertex * tVertex : mVertices )
            {
                // initialize counter
                uint tCount = 0;

                // get nuber of cells
                uint tNumberOfCells = tVertex->get_number_of_cells();

                // count number of cells
                for( uint k=0; k<tNumberOfCells; ++k )
                {
                    tCount += tVertex->get_cell( k )->get_number_of_vertices();
                }

                // allocate array with indices
                Matrix< DDUMat > tIndices( tCount, 1 );

                // reset counter
                tCount = 0;
                for( uint k=0; k<tNumberOfCells; ++k )
                {
                    // get pointer to cell
                    Cell * tCell = tVertex->get_cell( k );

                    // get number of vertices per element
                    uint tNumberOfVertices = tCell->get_number_of_vertices();

                    // loop over all vertices
                    for( uint i=0; i<tNumberOfVertices; ++i )
                    {
                        tIndices( tCount++ ) = tCell->get_vertex( i )->get_index();
                    }
                }

                // make indices unique
                Matrix< DDUMat > tUniqueIndices;
                unique( tIndices, tUniqueIndices );

                // get my index
                uint tIndex = tVertex->get_index();

                // reset counter
                tCount = 0;

                uint tNumberOfVertices = tUniqueIndices.length();

                // init container
                tVertex->init_neighbor_container( tNumberOfVertices-1 );

                // loop over all indices
                for( uint i=0; i<tNumberOfVertices; ++i )
                {
                    // test that this is not me
                    if( tUniqueIndices( i ) != tIndex )
                    {
                        // insert neighbor
                        tVertex->insert_neighbor( mVertices( tUniqueIndices( i ) ), tCount++ );
                    }
                }
            }
        }

//-------------------------------------------------------------------------------
    } /* namespace sdf */
} /* namespace moris */

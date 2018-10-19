#include "cl_SDF_Vertex.hpp"
#include "cl_SDF_Cell.hpp"
#include "cl_SDF_Triangle.hpp"
namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------
        /**
         * constructor
         */
        Vertex::Vertex( const moris_index aIndex, const Matrix< DDRMat > & aNodeCoords ) :
            mIndex( aIndex ),
            mNodeCoords( 3, 1 )
        {
            // convert dynamic array to fixed array
            for( uint k=0; k<3; ++k )
            {
                mNodeCoords( k ) = aNodeCoords( k );
            }
        }

//-------------------------------------------------------------------------------

        void
        Vertex::update_udf( Triangle *  aTriangle )
        {
            // calculate distance to this point
            real tDistance = aTriangle->get_distance_to_point( mNodeCoords );

            if( tDistance < mSDF )
            {
                // remember value
                mSDF = tDistance;

                // set sdf flag
                mHasSDF = true;

                // remember triangle
                mClosestTriangle = aTriangle;
            }
        }

//-------------------------------------------------------------------------------


        void
        Vertex::insert_cell( Cell * aCell )
        {
            mCells( mCellCounter++ ) = aCell;
        }

// -----------------------------------------------------------------------------

        uint
        Vertex::sweep()
        {
            bool tSwept = false;

            // loop over all neighbors
            for( Vertex * tNeighbor : mNeighbors )
            {
                // get pointer to triangle
                Triangle * tTriangle = tNeighbor->get_closest_triangle();

                if( tTriangle != NULL )
                {
                    // get distance to triangle of neighbor
                    real tDistance
                        = tTriangle->get_distance_to_point( mNodeCoords );

                    if( tDistance < mSDF )
                    {
                        tSwept = true;
                        mSDF = tDistance;
                        mClosestTriangle = tTriangle;
                    }
                }
            }

            if( tSwept )
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }

//-------------------------------------------------------------------------------

    } /* namespace sdf */
} /* namespace moris */

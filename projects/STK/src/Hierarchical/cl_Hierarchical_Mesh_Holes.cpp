#include "cl_Hierarchical_Mesh_Holes.hpp" // STK/src/Hierarchical
namespace moris
{
    void
    create_holes(
            Hierarchical_Mesh_Main & aHMR,
            const moris::database  & aMesh )
    {


        // get radius from settings
        real tRadius = aHMR.mSettings.HoleRadius;

        // get exponent from settings
        real tExpo   = aHMR.mSettings.HoleExponent;

        // distance parameters for SDF 1 and SDF 2
        real tDelta1 = 0.0; //0.1;
        real tDelta2 = 0.0; //0.1;
        
        // flag to indicate where holes should be generated
        bool tHolesEverywhere = true;

        // relative distance between hole centers with respect to radius
        real tRelativeDistance = 4;

        // only do something if radius is set and there are at least three SDFs
        if ( tRadius > 0 && aHMR.mSettings.SDFObjectFiles.size() >= 3 )
        {
            // start timer
            tic tTimer;

            // ask mesh for nodes
            Mat<uint> tNodesOnProc = aMesh.get_entities_owned_and_shared_by_current_proc(EntityRank::NODE);
            uint  tNumberOfNodes   = tNodesOnProc.length();

            // create vector of node coordinates
            // get min and max coordinates
            Mat<real> tMinNode( 3, 1,  1e20 );
            Mat<real> tMaxNode( 3, 1, -1e20);
            for (uint k = 0; k < tNumberOfNodes; ++k)
            {
                // get the coordinates of this node
                Mat<real> tNodeCoordinate
                = aMesh.get_selected_nodes_coords(tNodesOnProc.rows(k,k));

                // copy min and max value
                for( uint i=0; i<3; ++i )
                {
                    if( tNodeCoordinate( i ) < tMinNode( i ) )
                    {
                        tMinNode( i ) = tNodeCoordinate( i );
                    }
                    if( tNodeCoordinate( i ) > tMaxNode( i ) )
                    {
                        tMaxNode( i ) = tNodeCoordinate( i );
                    }
                }
            }

            // get length of domain
            Mat<real> tLength( 3, 1 );
            for( uint i=0; i<3; ++i )
            {
                tLength( i ) = tMaxNode( i ) - tMinNode( i );
            }

            // distance between centers
            real tStep = tRelativeDistance*tRadius;

            // calculate number of Holes per direction
            uint tNumberOfHoles = 1;
            Mat< uint > tNumberOfHolesPerDirection ( 3, 1 );
            for( uint i=0; i<3; ++i )
            {
                tNumberOfHolesPerDirection( i ) = floor( ( (real) tLength( i ) )/tStep );

                // calculate total number of Holes
                tNumberOfHoles *= tNumberOfHolesPerDirection( i );
            }

            // create coordinates of Holes
            Mat< real > tHoles( 3, tNumberOfHoles );

            // Hole counter
            uint tCount = 0;

            // create hole centers

            // loop over all k
            real tZ = tMinNode( 2 ) + 0.5*tStep;
            for( uint k=0; k<tNumberOfHolesPerDirection( 2 ); ++k )
            {
                // loop over all j
                real tY = tMinNode( 1 ) + 0.5*tStep;
                for( uint j=0; j<tNumberOfHolesPerDirection( 1 ); ++j )
                {
                    // loop over all i
                    real tX = tMinNode( 0 ) + 0.5*tStep;
                    for( uint i=0; i<tNumberOfHolesPerDirection( 0 ); ++i )
                    {
                        // calculate coordinates of center
                        tHoles( 0, tCount ) = tX;
                        tHoles( 1, tCount ) = tY;
                        tHoles( 2, tCount ) = tZ;

                        // increment counter
                        ++tCount;

                        // increment X
                        tX += tStep;
                    }
                    // increment Y
                    tY += tStep;
                }
                // increment Z
                tZ += tStep;
            }

            // print debug messahe
            std::fprintf(stdout, "Checking criterion for %u holes...\n", ( unsigned int ) tNumberOfHoles );

            // reset output array
            Mat< real > tValues( tNumberOfNodes, 1, -1e20 );
            BoostBitset tFlags( tNumberOfNodes );
            // reset counter
            tCount = 0;

            // dummy value for nodes not regarded
            real tMinVal = 1E20;

            // define epsilon environment for Holes to test
            real tEpsilon = 2.5*tStep;

            for( uint k=0; k<tNumberOfNodes; ++k )
            {
                // get the coordinates of this node
                Mat<real> tNodeCoordinate
                = aMesh.get_selected_nodes_coords(tNodesOnProc.rows(k,k));

                // test if coordinate is inside volume
                if ( tHolesEverywhere || aHMR.mFieldData.ObjectSDF( 0 )( k ) <= 0 )
                {
                    // loop over all Holes
                    for ( uint j=0; j<tNumberOfHoles; ++j )
                    {
                        // calculate distance between point and center
                        real tDeltaX = std::abs( tNodeCoordinate( 0 ) - tHoles( 0, j ) );
                        real tDeltaY = std::abs( tNodeCoordinate( 1 ) - tHoles( 1, j ) );
                        real tDeltaZ = std::abs( tNodeCoordinate( 2 ) - tHoles( 2, j ) );

                        // test if Hole is in the proximity of the Hole
                        // ( we don't want to do this too often. Pow is expensive )
                        if ( tDeltaX < tEpsilon && tDeltaY < tEpsilon && tDeltaZ < tEpsilon )
                        {
                            // calculate value
                            real tValue = tRadius - std::pow(
                                    std::pow( tDeltaX , tExpo )
                            + std::pow( tDeltaY , tExpo )
                            + std::pow( tDeltaZ , tExpo ),
                            1/tExpo );

                            // take value if it is bigger
                            if ( tValue > tValues( k ) )
                            {
                                tValues( k ) = tValue;
                            }

                            // set have data switch to true
                            tFlags.set( k );

                            // remember minumum value
                            if ( tValue < tMinVal )
                            {
                                tMinVal = tValue;
                            }
                        }
                    }

                    // test for other SDFs
                    if (       aHMR.mFieldData.ObjectSDF( 1 )( k ) <= tDelta1
                            || aHMR.mFieldData.ObjectSDF( 2 )( k ) <= tDelta2 )
                    {
                        if ( tValues( k ) > 0 )
                        {
                            // flip value
                            tValues( k ) *= -1;
                        }
                    }
                }
            }

            // fill nodes outside of object 0 with dummu values
            for ( uint k=0; k<tNumberOfNodes; ++k )
            {
                // test if value exists
                if ( ! tFlags.test( k ) )
                {
                    tValues( k ) = tMinVal;
                }
            }

            // copy values into HMR
            // get number of objects
            uint tNumberOfObjects = aHMR.mSettings.SDFObjectFiles.size() ;
            if (tNumberOfObjects == aHMR.mFieldData.ObjectSDF.size() )
            {
                // if field was not created, do a pushback
                aHMR.mFieldData.ObjectSDF.push_back( tValues );
                aHMR.mFieldData.ObjectFlagSDF.push_back( tFlags );
            }
            else
            {
                // overwrite field
                aHMR.mFieldData.ObjectSDF( tNumberOfObjects ) = tValues;
                aHMR.mFieldData.ObjectFlagSDF( tNumberOfObjects ) = tFlags;
            }

            // stop the timer

            moris::real tElapsedTime
            = tTimer.toc<moris::chronos::milliseconds>().wall;

            // print elapsed time
            std::fprintf(stdout, "Time for creating holes           : %5.3f [sec]\n",
                    tElapsedTime/1000 );

        }
    }
}

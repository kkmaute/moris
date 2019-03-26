#include "cl_FEM_Integration_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp"       //FEM/INT/src
#include "op_times.hpp"
#include "fn_trans.hpp"
#include "fn_reshape.hpp"


namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Integrator::Integrator( const Integration_Rule & aIntegrationRule )
        {
            // create space rule
            mSpaceCoeffs = aIntegrationRule.create_space_coeffs();

            // create time rule
            mTimeCoeffs  = aIntegrationRule.create_time_coeffs();
        }

//------------------------------------------------------------------------------

        Integrator::~Integrator()
        {
            // delete space coeffs if they exist
            if ( mSpaceCoeffs != NULL )
            {
                delete mSpaceCoeffs;
            }

            // delete time coeffs if they exist
            if( mTimeCoeffs != NULL )
            {
                delete mTimeCoeffs;
            }
        }

//------------------------------------------------------------------------------

        uint Integrator::get_number_of_points()
        {
            return mSpaceCoeffs->get_number_of_points() * mTimeCoeffs->get_number_of_points();
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Integrator::get_points()
        {
            // matrix with space points
            Matrix< DDRMat > tSpacePoints  = mSpaceCoeffs->get_points();

            // matrix with time points
            Matrix< DDRMat > tTimePoints   = mTimeCoeffs->get_points();

            // get number of points in space
            uint tNumOfSpacePoints = mSpaceCoeffs->get_number_of_points();

            // get number of points in time
            uint tNumOfTimePoints  = mTimeCoeffs->get_number_of_points();

            // get number of dimensions in space
            uint tNumOfSpaceDim = mSpaceCoeffs->get_number_of_dimensions();

            // create output matrix for space time
            Matrix< DDRMat > tPoints( tNumOfSpaceDim + 1,
                                      tNumOfSpacePoints*tNumOfTimePoints );

//            // initialize counter
//            uint tCount = 0;
//            // loop over time
//            for( uint k = 0; k < tNumOfTimePoints; ++k )
//            {
//                // loop over space
//                for( uint j = 0; j < tNumOfSpacePoints; ++j )
//                {
//                    // loop over dimensions in space
//                    for( uint i = 0; i < tNumOfSpaceDim; ++i )
//                    {
//                        tPoints( i, tCount ) = tSpacePoints( i, j );
//                    }
//                    // time
//                    tPoints( tNumOfSpaceDim, tCount++ ) = tTimePoints( k );
//                }
//            }

            Matrix< DDRMat > tOnes( 1, tNumOfSpacePoints, 1.0 );

            // loop over time
            uint startCol, stopCol;
            for( uint k = 0; k < tNumOfTimePoints; ++k )
            {
                // indices for columns
                startCol = k * tNumOfSpacePoints;
                stopCol  = ( k + 1 ) * tNumOfSpacePoints - 1;

                // fill in the space points coordinates
                tPoints( { 0, tNumOfSpaceDim-1 }, { startCol, stopCol } )
                    = tSpacePoints.matrix_data();

                // fill in the time point coordinates
                tPoints( { tNumOfSpaceDim, tNumOfSpaceDim }, { startCol, stopCol } )
                    = tTimePoints( k ) * tOnes;
            }
            return tPoints;
        }

//------------------------------------------------------------------------------

//        uint Integrator::get_number_of_dimensions()
//        {
//            return mSpaceCoeffs->get_number_of_dimensions() + 1;
//        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Integrator::get_weights()
        {

            // matrix with space weights
            Matrix< DDRMat > tSpaceWeights = mSpaceCoeffs->get_weights();

            // matrix with time weights
            Matrix< DDRMat > tTimeWeights  = mTimeCoeffs->get_weights();

            // get number of points in space
            uint tNumOfSpacePoints = mSpaceCoeffs->get_number_of_points();

            // get number of points in time
            uint tNumOfTimePoints  = mTimeCoeffs->get_number_of_points();

//            // create output matrix for weights
//            Matrix< DDRMat > tWeights( 1, tNumOfSpacePoints*tNumOfTimePoints );
//
//            // initialize counter
//            uint tCount = 0;
//
//            // loop over time
//            for( uint k = 0; k < tNumberOfTimePoints; ++k )
//            {
//                // loop over space
//                for( uint j = 0; j < tNumberOfSpacePoints; ++j )
//                {
//                    // weight
//                    tWeights( tCount++ ) = tSpaceWeights( j ) * tTimeWeights( k );
//                }
//            }

            Matrix< DDRMat > tWeights = reshape ( trans( tSpaceWeights ) * tTimeWeights, 1, tNumOfSpacePoints*tNumOfTimePoints );

            return tWeights;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

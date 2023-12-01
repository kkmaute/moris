/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integrator.cpp
 *
 */

#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Integrator.hpp"

//LINALG/src
#include "op_times.hpp"
#include "fn_trans.hpp"
#include "fn_vectorize.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        Integrator::Integrator( const Integration_Rule & aIntegrationRule )
        {
            // create space rule
            mSpaceCoeffs = aIntegrationRule.create_space_coeffs();

            // create time rule
            mTimeCoeffs  = aIntegrationRule.create_time_coeffs();

            // get number of points in space
            mNumOfSpacePoints = mSpaceCoeffs->get_number_of_points();

            // get number of points in time
            mNumOfTimePoints  = mTimeCoeffs->get_number_of_points();

            // matrix with space points
            mSpaceCoeffs->get_points( mSpacePoints );

            // matrix with time points
            mTimeCoeffs->get_points( mTimePoints );

            // matrix with space weights
            mSpaceCoeffs->get_weights( mSpaceWeights );

            // matrix with time weights
            mTimeCoeffs->get_weights( mTimeWeights );
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
            return mNumOfSpacePoints * mNumOfTimePoints;
        }

        //------------------------------------------------------------------------------

        void Integrator::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            // get number of dimensions in space
            uint tNumOfSpaceDim = mSpaceCoeffs->get_number_of_dimensions();

            // set output matrix size for space time
            aIntegrationPoints.set_size( tNumOfSpaceDim + 1,
                    mNumOfSpacePoints * mNumOfTimePoints );

            Matrix< DDRMat > tOnes( 1, mNumOfSpacePoints, 1.0 );

            // loop over time
            uint startCol, stopCol;
            for( uint k = 0; k < mNumOfTimePoints; ++k )
            {
                // indices for columns
                startCol = k * mNumOfSpacePoints;
                stopCol  = ( k + 1 ) * mNumOfSpacePoints - 1;

                // fill in the space points coordinates
                aIntegrationPoints( { 0, tNumOfSpaceDim - 1 }, { startCol, stopCol } ) =
                        mSpacePoints.matrix_data();

                // fill in the time point coordinates
                aIntegrationPoints( { tNumOfSpaceDim, tNumOfSpaceDim }, { startCol, stopCol } ) =
                        mTimePoints( k ) * tOnes;
            }
        }

        //------------------------------------------------------------------------------

        void Integrator::get_weights( Matrix< DDRMat > & aIntegrationWeights )
        {
            // get weights
            aIntegrationWeights = trans(
                vectorize ( trans( mSpaceWeights ) * mTimeWeights ) );
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */


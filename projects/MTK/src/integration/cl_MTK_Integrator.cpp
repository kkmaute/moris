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

// LINALG/src
#include "op_times.hpp"
#include "fn_trans.hpp"
#include "fn_vectorize.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        Integrator::Integrator( const Integration_Rule &aIntegrationRule )
                : mSpaceCoeffs( aIntegrationRule.create_space_coeffs() )
                , mTimeCoeffs( aIntegrationRule.create_time_coeffs() )
                , mNumOfSpacePoints( mSpaceCoeffs->get_number_of_points() )
                , mNumOfTimePoints( mTimeCoeffs->get_number_of_points() )
        {
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

        // Integrator::~Integrator()
        // {
        //     // delete space coeffs if they exist
        //     if ( mSpaceCoeffs != NULL )
        //     {
        //         delete mSpaceCoeffs;
        //     }
        //
        //     // delete time coeffs if they exist
        //     if ( mTimeCoeffs != NULL )
        //     {
        //         delete mTimeCoeffs;
        //     }
        // }

        //------------------------------------------------------------------------------

        uint Integrator::get_number_of_points() const
        {
            return mNumOfSpacePoints * mNumOfTimePoints;
        }

        //------------------------------------------------------------------------------

        void Integrator::get_points( Matrix< DDRMat > &aIntegrationPoints ) const
        {
            // get number of dimensions in space
            uint tNumOfSpaceDim = mSpaceCoeffs->get_number_of_dimensions();

            // set output matrix size for space time
            aIntegrationPoints.set_size( tNumOfSpaceDim + 1,
                    mNumOfSpacePoints * mNumOfTimePoints );

            Matrix< DDRMat > tOnes( 1, mNumOfSpacePoints, 1.0 );

            // loop over time
            uint startCol, stopCol;
            for ( uint k = 0; k < mNumOfTimePoints; ++k )
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

        Matrix< DDRMat > Integrator::get_points() const
        {
            Matrix< DDRMat > tPoints;
            get_points( tPoints );
            return tPoints;
        }


        //------------------------------------------------------------------------------

        void Integrator::get_weights( Matrix< DDRMat > &aIntegrationWeights ) const
        {
            // get weights
            aIntegrationWeights = trans(
                    vectorize( trans( mSpaceWeights ) * mTimeWeights ) );
        }

        Matrix< DDRMat > Integrator::get_weights() const
        {
            Matrix< DDRMat > tWeights;
            get_weights( tWeights );
            return tWeights;
        }


        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

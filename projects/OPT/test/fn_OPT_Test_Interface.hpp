/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_OPT_Test_Interface.hpp
 *
 */

#pragma once

namespace moris::opt
{

    //--------------------------------------------------------------------------------------------------------------

    void initialize_test_1( Vector< real >& aADVs, Vector< real >& aLowerBounds, Vector< real >& aUpperBounds )
    {
        // Initial Guess
        aADVs.resize( 2 );
        aADVs( 0 ) = 1.0;
        aADVs( 1 ) = 2.0;

        // Lower Bounds
        aLowerBounds.resize( 2 );
        aLowerBounds( 0 ) = 1.0;
        aLowerBounds( 1 ) = 2.0;

        // Upper Bounds
        aUpperBounds.resize( 2 );
        aUpperBounds( 0 ) = 1.0;
        aUpperBounds( 1 ) = 2.0;
    }

    //--------------------------------------------------------------------------------------------------------------

    void initialize_test_2( Vector< real >& aADVs, Vector< real >& aLowerBounds, Vector< real >& aUpperBounds )
    {
        // Initial Guess
        aADVs.resize( 2 );
        aADVs( 0 ) = 3.0;
        aADVs( 1 ) = 4.0;

        // Lower Bounds
        aLowerBounds.resize( 2 );
        aLowerBounds( 0 ) = 3.0;
        aLowerBounds( 1 ) = 4.0;

        // Upper Bounds
        aUpperBounds.resize( 2 );
        aUpperBounds( 0 ) = 3.0;
        aUpperBounds( 1 ) = 4.0;
    }

    //--------------------------------------------------------------------------------------------------------------

    void initialize_test_3( Vector< real >& aADVs, Vector< real >& aLowerBounds, Vector< real >& aUpperBounds )
    {
        // Initial Guess
        aADVs.resize( 2 );
        aADVs( 0 ) = 5.0;
        aADVs( 1 ) = 6.0;

        // Lower Bounds
        aLowerBounds.resize( 2 );
        aLowerBounds( 0 ) = 5.0;
        aLowerBounds( 1 ) = 6.0;

        // Upper Bounds
        aUpperBounds.resize( 2 );
        aUpperBounds( 0 ) = 5.0;
        aUpperBounds( 1 ) = 6.0;
    }

    //--------------------------------------------------------------------------------------------------------------

    void initialize_test_4( Vector< real >& aADVs, Vector< real >& aLowerBounds, Vector< real >& aUpperBounds )
    {
        // Initial Guess
        aADVs.resize( 2 );
        aADVs( 0 ) = 7.0;
        aADVs( 1 ) = 8.0;

        // Lower Bounds
        aLowerBounds.resize( 2 );
        aLowerBounds( 0 ) = 7.0;
        aLowerBounds( 1 ) = 8.0;

        // Upper Bounds
        aUpperBounds.resize( 2 );
        aUpperBounds( 0 ) = 7.0;
        aUpperBounds( 1 ) = 8.0;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< real > get_criteria_test( const Vector< real >& aADVs )
    {
        Vector< real > tCriteria( 2 );
        tCriteria( 0 ) = aADVs( 0 );
        tCriteria( 1 ) = aADVs( 1 );

        return tCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > get_dcriteria_dadv_test( const Vector< real >& aADVs )
    {
        Matrix< DDRMat > tDCriteria( 2, 2, 0.0 );
        tDCriteria( 0, 0 ) = aADVs( 0 );
        tDCriteria( 0, 1 ) = 0;
        tDCriteria( 1, 0 ) = 0;
        tDCriteria( 1, 1 ) = aADVs( 1 );

        return tDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::opt

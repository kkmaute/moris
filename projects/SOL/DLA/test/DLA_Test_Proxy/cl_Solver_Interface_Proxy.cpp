/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Solver_Interface_Proxy.cpp
 *
 */

#include "cl_Solver_Interface_Proxy.hpp"
#include "cl_Communication_Tools.hpp"    // COM/src

using namespace moris;

Solver_Interface_Proxy::Solver_Interface_Proxy()
{
    // Determine process rank
    size_t rank = par_rank();
    size_t size = par_size();

    mUseMatrixMarketFiles = false;

    if ( size == 4 )
    {
        // Set input values
        mNumMyDofs         = 0;
        mNumDofsPerElement = 8;    // dofs per element

        mElementMatrixValues.resize( 64, 1 );
        mMyRHSValues.resize( 1 );
        mMyRHSValues( 0 ).resize( 8, 1 );

        // Define input test values
        switch ( rank )
        {
            case 0:
                mNumMyDofs = 8;

                mMyGlobalElements            = { { 0 }, { 1 }, { 8 }, { 9 }, { 16 }, { 17 }, { 14 }, { 15 } };
                mMyGlobalElementsOverlapping = { { 0 }, { 1 }, { 8 }, { 9 }, { 16 }, { 17 }, { 14 }, { 15 } };

                mNumElements      = 1;
                mMyConstraintDofs = { { 0 }, { 1 } };

                mEleDofConectivity = { { 0 }, { 1 }, { 8 }, { 9 }, { 16 }, { 17 }, { 14 }, { 15 } };

                mMyRHSValues( 0 ).fill( 0.0 );
                break;
            case 1:
                mNumMyDofs = 4;

                mMyGlobalElements            = { { 2 }, { 3 }, { 10 }, { 11 } };
                mMyGlobalElementsOverlapping = { { 8 }, { 9 }, { 2 }, { 3 }, { 10 }, { 11 }, { 16 }, { 17 } };

                mNumElements      = 1;
                mMyConstraintDofs = { { 3 } };

                mEleDofConectivity = { { 8 }, { 9 }, { 2 }, { 3 }, { 10 }, { 11 }, { 16 }, { 17 } };

                mMyRHSValues( 0 ).fill( 0.0 );
                break;
            case 2:
                mNumMyDofs = 2;

                mMyGlobalElements            = { { 4 }, { 5 } };
                mMyGlobalElementsOverlapping = { { 16 }, { 17 }, { 10 }, { 11 }, { 4 }, { 5 }, { 12 }, { 13 } };

                mNumElements       = 1;
                mEleDofConectivity = { { 16 }, { 17 }, { 10 }, { 11 }, { 4 }, { 5 }, { 12 }, { 13 } };

                mMyRHSValues( 0 ) = { { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 1.0 }, { 0.0 }, { 0.0 } };
                break;
            case 3:
                mNumMyDofs = 4;

                mMyGlobalElements            = { { 12 }, { 13 }, { 6 }, { 7 } };
                mMyGlobalElementsOverlapping = { { 14 }, { 15 }, { 16 }, { 17 }, { 12 }, { 13 }, { 6 }, { 7 } };

                mNumElements       = 1;
                mEleDofConectivity = { { 14 }, { 15 }, { 16 }, { 17 }, { 12 }, { 13 }, { 6 }, { 7 } };

                mMyRHSValues( 0 ).fill( 0.0 );
                break;
        }

        // Define elemental matrix values
        mElementMatrixValues( 0, 0 )  = 12;
        mElementMatrixValues( 1, 0 )  = 3;
        mElementMatrixValues( 2, 0 )  = -6;
        mElementMatrixValues( 3, 0 )  = -3;
        mElementMatrixValues( 4, 0 )  = -6;
        mElementMatrixValues( 5, 0 )  = -3;
        mElementMatrixValues( 6, 0 )  = 0;
        mElementMatrixValues( 7, 0 )  = 3;
        mElementMatrixValues( 8, 0 )  = 3;
        mElementMatrixValues( 9, 0 )  = 12;
        mElementMatrixValues( 10, 0 ) = 3;
        mElementMatrixValues( 11, 0 ) = 0;
        mElementMatrixValues( 12, 0 ) = -3;
        mElementMatrixValues( 13, 0 ) = -6;
        mElementMatrixValues( 14, 0 ) = -3;
        mElementMatrixValues( 15, 0 ) = -6;
        mElementMatrixValues( 16, 0 ) = -6;
        mElementMatrixValues( 17, 0 ) = 3;
        mElementMatrixValues( 18, 0 ) = 12;
        mElementMatrixValues( 19, 0 ) = -3;
        mElementMatrixValues( 20, 0 ) = 0;
        mElementMatrixValues( 21, 0 ) = -3;
        mElementMatrixValues( 22, 0 ) = -6;
        mElementMatrixValues( 23, 0 ) = 3;
        mElementMatrixValues( 24, 0 ) = -3;
        mElementMatrixValues( 25, 0 ) = 0;
        mElementMatrixValues( 26, 0 ) = -3;
        mElementMatrixValues( 27, 0 ) = 12;
        mElementMatrixValues( 28, 0 ) = 3;
        mElementMatrixValues( 29, 0 ) = -6;
        mElementMatrixValues( 30, 0 ) = 3;
        mElementMatrixValues( 31, 0 ) = -6;
        mElementMatrixValues( 32, 0 ) = -6;
        mElementMatrixValues( 33, 0 ) = -3;
        mElementMatrixValues( 34, 0 ) = 0;
        mElementMatrixValues( 35, 0 ) = 3;
        mElementMatrixValues( 36, 0 ) = 12;
        mElementMatrixValues( 37, 0 ) = 3;
        mElementMatrixValues( 38, 0 ) = -6;
        mElementMatrixValues( 39, 0 ) = -3;
        mElementMatrixValues( 40, 0 ) = -3;
        mElementMatrixValues( 41, 0 ) = -6;
        mElementMatrixValues( 42, 0 ) = -3;
        mElementMatrixValues( 43, 0 ) = -6;
        mElementMatrixValues( 44, 0 ) = 3;
        mElementMatrixValues( 45, 0 ) = 12;
        mElementMatrixValues( 46, 0 ) = 3;
        mElementMatrixValues( 47, 0 ) = 0;
        mElementMatrixValues( 48, 0 ) = 0;
        mElementMatrixValues( 49, 0 ) = -3;
        mElementMatrixValues( 50, 0 ) = -6;
        mElementMatrixValues( 51, 0 ) = 3;
        mElementMatrixValues( 52, 0 ) = -6;
        mElementMatrixValues( 53, 0 ) = 3;
        mElementMatrixValues( 54, 0 ) = 12;
        mElementMatrixValues( 55, 0 ) = -3;
        mElementMatrixValues( 56, 0 ) = 3;
        mElementMatrixValues( 57, 0 ) = -6;
        mElementMatrixValues( 58, 0 ) = 3;
        mElementMatrixValues( 59, 0 ) = -6;
        mElementMatrixValues( 60, 0 ) = -3;
        mElementMatrixValues( 61, 0 ) = 0;
        mElementMatrixValues( 62, 0 ) = -3;
        mElementMatrixValues( 63, 0 ) = 12;
    }
}

Solver_Interface_Proxy::Solver_Interface_Proxy( uint aNumRHS )
{
    // Determine process rank
    size_t size = par_size();

    mUseMatrixMarketFiles = false;

    if ( size == 1 )
    {
        // check that aNumRHS < 3
        MORIS_ERROR( aNumRHS < 3,
                "Solver_Interface_Proxy::Solver_Interface_Proxy - aNumRHS larger needs to be smaller than 3" );

        // Set input values
        mNumDofsPerElement = 8;    // dofs per element
        mNumRHS            = aNumRHS;

        mElementMatrixValues.resize( 64, 1 );
        mMyRHSValues.resize( mNumRHS );

        mNumMyDofs                   = 18;
        mMyGlobalElements            = { { 0 }, { 1 }, { 2 }, { 3 }, { 4 }, { 5 }, { 6 }, { 7 }, { 8 }, { 9 }, { 10 }, { 11 }, { 12 }, { 13 }, { 14 }, { 15 }, { 16 }, { 17 } };
        mMyGlobalElementsOverlapping = mMyGlobalElements;

        mNumElements = 4;

        mMyConstraintDofs = { { 0 }, { 1 }, { 3 } };

        mEleDofConectivity = { { 0, 8, 16, 14 },
            { 1, 9, 17, 15 },
            { 8, 2, 10, 16 },
            { 9, 3, 11, 17 },
            { 16, 10, 4, 12 },
            { 17, 11, 5, 13 },
            { 14, 16, 12, 6 },
            { 15, 17, 13, 7 } };

        mMyRHSValues( 0 ) = { { 0, 0, 0, 0 },
            { 0, 0, 0, 0 },
            { 0, 0, 0, 0 },
            { 0, 0, 0, 0 },
            { 0, 0, 0, 0 },
            { 0, 0, 1, 0 },
            { 0, 0, 0, 0 },
            { 0, 0, 0, 0 } };

        if ( mNumRHS > 1 )
        {
            mMyRHSValues( 1 ) = { { 0, 0, 0, 0 },
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 },
                { 0, 0, 1, 0 },
                { 0, 0, 0, 0 },
                { 0, 0, 0, 0 } };
        }

        // Define elemental matrix values
        mElementMatrixValues( 0, 0 )  = 12;
        mElementMatrixValues( 1, 0 )  = 3;
        mElementMatrixValues( 2, 0 )  = -6;
        mElementMatrixValues( 3, 0 )  = -3;
        mElementMatrixValues( 4, 0 )  = -6;
        mElementMatrixValues( 5, 0 )  = -3;
        mElementMatrixValues( 6, 0 )  = 0;
        mElementMatrixValues( 7, 0 )  = 3;
        mElementMatrixValues( 8, 0 )  = 3;
        mElementMatrixValues( 9, 0 )  = 12;
        mElementMatrixValues( 10, 0 ) = 3;
        mElementMatrixValues( 11, 0 ) = 0;
        mElementMatrixValues( 12, 0 ) = -3;
        mElementMatrixValues( 13, 0 ) = -6;
        mElementMatrixValues( 14, 0 ) = -3;
        mElementMatrixValues( 15, 0 ) = -6;
        mElementMatrixValues( 16, 0 ) = -6;
        mElementMatrixValues( 17, 0 ) = 3;
        mElementMatrixValues( 18, 0 ) = 12;
        mElementMatrixValues( 19, 0 ) = -3;
        mElementMatrixValues( 20, 0 ) = 0;
        mElementMatrixValues( 21, 0 ) = -3;
        mElementMatrixValues( 22, 0 ) = -6;
        mElementMatrixValues( 23, 0 ) = 3;
        mElementMatrixValues( 24, 0 ) = -3;
        mElementMatrixValues( 25, 0 ) = 0;
        mElementMatrixValues( 26, 0 ) = -3;
        mElementMatrixValues( 27, 0 ) = 12;
        mElementMatrixValues( 28, 0 ) = 3;
        mElementMatrixValues( 29, 0 ) = -6;
        mElementMatrixValues( 30, 0 ) = 3;
        mElementMatrixValues( 31, 0 ) = -6;
        mElementMatrixValues( 32, 0 ) = -6;
        mElementMatrixValues( 33, 0 ) = -3;
        mElementMatrixValues( 34, 0 ) = 0;
        mElementMatrixValues( 35, 0 ) = 3;
        mElementMatrixValues( 36, 0 ) = 12;
        mElementMatrixValues( 37, 0 ) = 3;
        mElementMatrixValues( 38, 0 ) = -6;
        mElementMatrixValues( 39, 0 ) = -3;
        mElementMatrixValues( 40, 0 ) = -3;
        mElementMatrixValues( 41, 0 ) = -6;
        mElementMatrixValues( 42, 0 ) = -3;
        mElementMatrixValues( 43, 0 ) = -6;
        mElementMatrixValues( 44, 0 ) = 3;
        mElementMatrixValues( 45, 0 ) = 12;
        mElementMatrixValues( 46, 0 ) = 3;
        mElementMatrixValues( 47, 0 ) = 0;
        mElementMatrixValues( 48, 0 ) = 0;
        mElementMatrixValues( 49, 0 ) = -3;
        mElementMatrixValues( 50, 0 ) = -6;
        mElementMatrixValues( 51, 0 ) = 3;
        mElementMatrixValues( 52, 0 ) = -6;
        mElementMatrixValues( 53, 0 ) = 3;
        mElementMatrixValues( 54, 0 ) = 12;
        mElementMatrixValues( 55, 0 ) = -3;
        mElementMatrixValues( 56, 0 ) = 3;
        mElementMatrixValues( 57, 0 ) = -6;
        mElementMatrixValues( 58, 0 ) = 3;
        mElementMatrixValues( 59, 0 ) = -6;
        mElementMatrixValues( 60, 0 ) = -3;
        mElementMatrixValues( 61, 0 ) = 0;
        mElementMatrixValues( 62, 0 ) = -3;
        mElementMatrixValues( 63, 0 ) = 12;
    }
}

// ---------------------------------------------------------------------------------------------

Solver_Interface_Proxy::Solver_Interface_Proxy( const std::string& aProblem )
{
    std::cout << "Problem is:" << aProblem << '\n';

    // Determine process rank
    size_t size = par_size();

    std::cout << size << '\n';

    mUseMatrixMarketFiles = false;

    if ( size == 1 )
    {
        mElementMatrixValues.resize( 64, 1 );
        mElementMatrixValues.fill( 1.0 );

        mElementMassMatrixValues.resize( 64, 1 );
        mElementMassMatrixValues.fill( 2.0 );

        mNumMyDofs = 8;

        mMyGlobalElements            = { { 0 }, { 1 }, { 8 }, { 9 }, { 16 }, { 17 }, { 14 }, { 15 } };
        mMyGlobalElementsOverlapping = { { 0 }, { 1 }, { 8 }, { 9 }, { 16 }, { 17 }, { 14 }, { 15 } };

        mNumElements      = 1;
        mMyConstraintDofs = { { 0 }, { 1 } };

        mEleDofConectivity = { { 0 }, { 1 }, { 8 }, { 9 }, { 16 }, { 17 }, { 14 }, { 15 } };
    }
}

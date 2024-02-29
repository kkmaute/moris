/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Solver_Interface_Proxy.cpp
 *
 */

#include "cl_MSI_Solver_Interface_Proxy.hpp"
#include "cl_Communication_Tools.hpp"    // COM/src

using namespace moris::MSI;

MSI_Solver_Interface_Proxy::MSI_Solver_Interface_Proxy()
{
    // Determine process rank
    size_t rank = par_rank();
    size_t size = par_size();

    mUseMatrixMarketFiles = false;

    if ( size )
    {
        // Set input values
        mNumMyDofs         = 0;
        mNumDofsPerElement = 8;    // dofs per element

        mElementMatrixValues.resize( 64, 1 );
        mMyRHSValues.resize( 1 );
        mMyRHSValues( 0 ).resize( 8, 1 );

        mCommTable.set_size( 1, 3 );

        // Define input test values
        switch ( rank )
        {
            case 0:
                mNumMyDofs = 8;

                mMyGlobalElements            = { { 0 }, { 1 }, { 8 }, { 9 }, { 16 }, { 17 }, { 14 }, { 15 } };
                mMyGlobalElementsOverlapping = { { 0 }, { 1 }, { 8 }, { 9 }, { 16 }, { 17 }, { 14 }, { 15 } };

                mNumElements      = 2;
                mMyConstraintDofs = { { 0 }, { 1 } };

                mEleDofConectivity.resize( mNumElements );
                mEleDofConectivity( 0 ) = { { 0 }, { 1 }, { 8 }, { 9 } };
                mEleDofConectivity( 1 ) = { { 16 }, { 17 }, { 14 }, { 15 } };

                mMyRHSValues( 0 ).fill( 0.0 );

                mCommTable                         = { { 1 }, { 2 }, { 3 } };
                mMyGlobalElementsOverlappingOwners = { { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 }, { 0 } };
                break;

            case 1:
                mNumMyDofs = 4;

                mMyGlobalElements            = { { 2 }, { 3 }, { 10 }, { 11 } };
                mMyGlobalElementsOverlapping = { { 8 }, { 9 }, { 2 }, { 3 }, { 10 }, { 11 }, { 16 }, { 17 } };

                mNumElements      = 3;
                mMyConstraintDofs = { { 3 } };

                mEleDofConectivity.resize( mNumElements );
                mEleDofConectivity( 0 ) = { { 8 }, { 9 }, { 2 }, { 3 } };
                mEleDofConectivity( 1 ) = { { 10 }, { 11 } };
                mEleDofConectivity( 2 ) = { { 16 }, { 17 } };

                mMyRHSValues( 0 ).fill( 0.0 );
                mCommTable                         = { { 0 }, { 2 }, { 3 } };
                mMyGlobalElementsOverlappingOwners = { { 0 }, { 0 }, { 1 }, { 1 }, { 1 }, { 1 }, { 0 }, { 0 } };
                break;
            case 2:
                mNumMyDofs = 2;

                mMyGlobalElements            = { { 4 }, { 5 } };
                mMyGlobalElementsOverlapping = { { 16 }, { 17 }, { 10 }, { 11 }, { 4 }, { 5 }, { 12 }, { 13 } };

                mNumElements = 4;
                mEleDofConectivity.resize( mNumElements );
                mEleDofConectivity( 0 ) = { { 17 }, { 10 }, { 11 }, { 4 }, { 12 }, { 13 } };
                mEleDofConectivity( 1 ) = { { 16 }, { 13 } };
                mEleDofConectivity( 2 ) = { { 17 }, { 10 }, { 11 }, { 4 }, { 5 }, { 12 }, { 13 } };
                mEleDofConectivity( 3 ) = { { 16 }, { 17 }, { 10 } };

                mMyRHSValues( 0 ) = { { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 }, { 1.0 }, { 0.0 }, { 0.0 } };

                mCommTable                         = { { 0 }, { 1 }, { 3 } };
                mMyGlobalElementsOverlappingOwners = { { 0 }, { 0 }, { 1 }, { 1 }, { 2 }, { 2 }, { 3 }, { 3 } };

                break;
            case 3:
                mNumMyDofs = 4;

                mMyGlobalElements            = { { 12 }, { 13 }, { 6 }, { 7 } };
                mMyGlobalElementsOverlapping = { { 14 }, { 15 }, { 16 }, { 17 }, { 12 }, { 13 }, { 6 }, { 7 } };

                mNumElements = 1;
                mEleDofConectivity.resize( mNumElements );
                mEleDofConectivity( 0 ) = { { 6 }, { 7 } };

                mMyRHSValues( 0 ).fill( 0.0 );

                mCommTable                         = { { 0 }, { 1 }, { 2 } };
                mMyGlobalElementsOverlappingOwners = { { 0 }, { 0 }, { 0 }, { 0 }, { 3 }, { 3 }, { 3 }, { 3 } };

                break;
        }
    }
}

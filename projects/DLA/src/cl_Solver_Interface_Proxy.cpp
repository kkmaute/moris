/*
 * cl_Solver_Input_Test.cpp
 *
 *  Created on: Jun 18, 2018
 *      Author: schmidt
 */
#include "cl_Solver_Interface_Proxy.hpp"
#include "cl_Communication_Tools.hpp" // COM/src

using namespace moris;

Solver_Interface_Proxy::Solver_Interface_Proxy()
{
    // Determine process rank
    size_t rank = par_rank();
    size_t size = par_size();

    mUseMatrixMarketFiles = false;

    if (size == 4)
    {
    // Set input values
    mNumMyDofs = 0;
    mNumDofsPerElement = 8;                              // dofs per element

    mEleDofConectivity.resize( 8, 1 );
    mElementMatrixValues.resize( 64, 1 );
    mMyRHSValues.resize( 8, 1 );

    // Define input test values
    switch( rank )
        {
        case 0:
          mNumMyDofs = 8;
          mMyGlobalElements.resize( mNumMyDofs, 1 );
          mMyConstraintDofs.resize( 2, 1 );
          mMyGlobalElements(0,0) = 0;    mMyGlobalElements(1,0) = 1;  mMyGlobalElements(2,0) = 8;    mMyGlobalElements(3,0) = 9;    mMyGlobalElements(4,0) = 16;    mMyGlobalElements(5,0) = 17;    mMyGlobalElements(6,0) = 14;    mMyGlobalElements(7,0) = 15;
          mNumElements = 1;
          mMyConstraintDofs(0,0) = 0;    mMyConstraintDofs(1,0) = 1;
          mEleDofConectivity( 0, 0) = 0;   mEleDofConectivity( 1, 0) = 1; mEleDofConectivity( 2, 0) = 8; mEleDofConectivity( 3, 0) = 9; mEleDofConectivity( 4, 0) = 16; mEleDofConectivity( 5, 0) = 17; mEleDofConectivity( 6, 0) = 14; mEleDofConectivity( 7, 0) = 15;
          mMyRHSValues.fill( 0.0 );
          break;
        case 1:
          mNumMyDofs = 4;
          mMyGlobalElements.resize( mNumMyDofs, 1 );
          mMyConstraintDofs.resize( 1, 1 );
          mMyGlobalElements(0,0) = 2;   mMyGlobalElements(1,0) = 3;    mMyGlobalElements(2,0) = 10;    mMyGlobalElements(3,0) = 11;
          mNumElements = 1;
          mMyConstraintDofs(0,0) = 3;
          mEleDofConectivity( 0, 0) = 8;        mEleDofConectivity( 1, 0) = 9; mEleDofConectivity( 2, 0) = 2; mEleDofConectivity( 3, 0) = 3; mEleDofConectivity( 4, 0) = 10; mEleDofConectivity( 5, 0) = 11; mEleDofConectivity( 6, 0) = 16; mEleDofConectivity( 7, 0) = 17;
          mMyRHSValues.fill( 0.0 );
          break;
        case 2:
          mNumMyDofs = 2;
          mMyGlobalElements.resize( mNumMyDofs, 1 );
          mMyGlobalElements(0,0) = 4;    mMyGlobalElements(1,0) = 5;
          mNumElements = 1;
          mEleDofConectivity( 0, 0) = 16;        mEleDofConectivity( 1, 0) = 17;    mEleDofConectivity( 2, 0) = 10;    mEleDofConectivity( 3, 0) = 11;    mEleDofConectivity( 4, 0) = 4;    mEleDofConectivity( 5, 0) = 5;    mEleDofConectivity( 6, 0) = 12;    mEleDofConectivity( 7, 0) = 13;
          mMyRHSValues( 0, 0 )= 0.0;      mMyRHSValues( 1, 0 )= 0.0;    mMyRHSValues( 2, 0 )= 0.0;    mMyRHSValues( 3, 0 )= 0.0;   mMyRHSValues( 4, 0 )= 0.0;
          mMyRHSValues( 5, 0 )= 1.0;    mMyRHSValues( 6, 0 )= 0.0;    mMyRHSValues( 7, 0 )= 0.0;
          break;
        case 3:
          mNumMyDofs = 4;
          mMyGlobalElements.resize( mNumMyDofs, 1 );
          mMyGlobalElements(0,0) = 12;    mMyGlobalElements(1,0) = 13;    mMyGlobalElements(2,0) = 6;    mMyGlobalElements(3,0) = 7;
          mNumElements = 1;
          mEleDofConectivity( 0, 0) = 14;    mEleDofConectivity( 1, 0) = 15;    mEleDofConectivity( 2, 0) = 16;    mEleDofConectivity( 3, 0) = 17;    mEleDofConectivity( 4, 0) = 12;    mEleDofConectivity( 5, 0) = 13;    mEleDofConectivity( 6, 0) = 6;        mEleDofConectivity( 7, 0) = 7;
          mMyRHSValues.fill( 0.0 );
          break;
        }

    // Define elemental matrix values
    mElementMatrixValues( 0, 0 ) = 12;    mElementMatrixValues( 1, 0 ) = 3;     mElementMatrixValues( 2, 0 ) = -6;    mElementMatrixValues( 3, 0 ) = -3;    mElementMatrixValues( 4, 0 ) = -6;    mElementMatrixValues( 5, 0 ) = -3;    mElementMatrixValues( 6, 0 ) = 0;     mElementMatrixValues( 7, 0 ) = 3;
    mElementMatrixValues( 8, 0 ) = 3;     mElementMatrixValues( 9, 0 ) = 12;    mElementMatrixValues( 10, 0 ) = 3;    mElementMatrixValues( 11, 0 ) = 0;    mElementMatrixValues( 12, 0 ) = -3;   mElementMatrixValues( 13, 0 ) = -6;   mElementMatrixValues( 14, 0 ) = -3;   mElementMatrixValues( 15, 0 ) = -6;
    mElementMatrixValues( 16, 0 ) = -6;   mElementMatrixValues( 17, 0 ) = 3;    mElementMatrixValues( 18, 0 ) = 12;   mElementMatrixValues( 19, 0 ) = -3;   mElementMatrixValues( 20, 0 ) = 0;    mElementMatrixValues( 21, 0 ) = -3;   mElementMatrixValues( 22, 0 ) = -6;   mElementMatrixValues( 23, 0 ) = 3;
    mElementMatrixValues( 24, 0 ) = -3;   mElementMatrixValues( 25, 0 ) = 0;    mElementMatrixValues( 26, 0 ) = -3;   mElementMatrixValues( 27, 0 ) = 12;   mElementMatrixValues( 28, 0 ) = 3;    mElementMatrixValues( 29, 0 ) = -6;   mElementMatrixValues( 30, 0 ) = 3;    mElementMatrixValues( 31, 0 ) = -6;
    mElementMatrixValues( 32, 0 ) = -6;   mElementMatrixValues( 33, 0 ) = -3;   mElementMatrixValues( 34, 0 ) = 0;    mElementMatrixValues( 35, 0 ) = 3;    mElementMatrixValues( 36, 0 ) = 12;   mElementMatrixValues( 37, 0 ) = 3;    mElementMatrixValues( 38, 0 ) = -6;   mElementMatrixValues( 39, 0 ) = -3;
    mElementMatrixValues( 40, 0 ) = -3;   mElementMatrixValues( 41, 0 ) = -6;   mElementMatrixValues( 42, 0 ) = -3;   mElementMatrixValues( 43, 0 ) = -6;   mElementMatrixValues( 44, 0 ) = 3;    mElementMatrixValues( 45, 0 ) = 12;   mElementMatrixValues( 46, 0 ) = 3;    mElementMatrixValues( 47, 0 ) = 0;
    mElementMatrixValues( 48, 0 ) = 0;    mElementMatrixValues( 49, 0 ) = -3;   mElementMatrixValues( 50, 0 ) = -6;   mElementMatrixValues( 51, 0 ) = 3;    mElementMatrixValues( 52, 0 ) = -6;   mElementMatrixValues( 53, 0 ) = 3;    mElementMatrixValues( 54, 0 ) = 12;   mElementMatrixValues( 55, 0 ) = -3;
    mElementMatrixValues( 56, 0 ) = 3;    mElementMatrixValues( 57, 0 ) = -6;   mElementMatrixValues( 58, 0 ) = 3;    mElementMatrixValues( 59, 0 ) = -6;   mElementMatrixValues( 60, 0 ) = -3;   mElementMatrixValues( 61, 0 ) = 0;    mElementMatrixValues( 62, 0 ) = -3;   mElementMatrixValues( 63, 0 ) = 12;
}
}


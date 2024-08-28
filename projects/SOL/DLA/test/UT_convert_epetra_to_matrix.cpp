/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Map_Test.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"    // ALG/src
#include "moris_typedefs.hpp"       // COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Tools.hpp"          // COM/src/

#include "cl_Map_Epetra.hpp"                   // DLA/src/
#include "cl_SOL_Matrix_Vector_Factory.hpp"    // DLA/src/
#include "Epetra_Map.h"                        // DLA/src/
#include "Epetra_CrsMatrix.h"                  // DLA/src/
#include "Epetra_RowMatrix.h"

#include "Epetra_MultiVector.h"

#include "fn_convert_epetra_operator_to_matrix.hpp"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

#include "fn_check_equal.hpp"

namespace moris
{
    TEST_CASE( "testing stuff", "[Epetraaa Mat],[tested]" )
    {
        // run the test only on 2 processors
        if ( par_size() != 2 )
        {
            return;
        }

        // Construct the ePetra object for communication
        Epetra_MpiComm tComm( MPI_COMM_WORLD );

        // construct the map
        int tMyElements = 0;

        // get the processor id
        int tMyPID = tComm.MyPID();

        switch ( tMyPID )
        {
            case 0:
                tMyElements = 2;
                break;
            case 1:
                tMyElements = 1;
                break;
        }

        // create a map
        Epetra_Map tMap( -1, tMyElements, 0, tComm );

        // Epetra matrix
        Epetra_CrsMatrix A = Epetra_CrsMatrix( Copy, tMap, 3 );

        // declare a validity of inserted values
        int tSuccess = 0;

        if ( tMyPID == 0 )
        {
            tSuccess = A.InsertGlobalValues( 0, 2, Vector< double >{ 1.0, 2.0 }.memptr(), Vector< int >{ 0, 1 }.memptr() );
            tSuccess = A.InsertGlobalValues( 1, 2, Vector< double >{ 2.0, 3.0 }.memptr(), Vector< int >{ 0, 2 }.memptr() );
        }
        else
        {
            tSuccess = A.InsertGlobalValues( 2, 2, Vector< double >{ 4.0, 5.0 }.memptr(), Vector< int >{ 1, 2 }.memptr() );
        }

        REQUIRE( tSuccess == 0 );

        A.FillComplete();

        Epetra_CrsMatrix B = convert_epetra_operator_to_matrix( &A, tMap );

        // check the rows of A and B are the same
        for ( int iRow = 0; iRow < B.NumMyRows(); ++iRow )
        {
            int     tNumEntriesA;
            double* tValuesA;
            int*    tIndicesA;

            int     tNumEntriesB;
            double* tValuesB;
            int*    tIndicesB;

            A.ExtractMyRowView( iRow, tNumEntriesA, tValuesA, tIndicesA );
            B.ExtractMyRowView( iRow, tNumEntriesB, tValuesB, tIndicesB );

            REQUIRE( tNumEntriesA == tNumEntriesB );
           
            for (int iEntry = 0; iEntry < tNumEntriesA; iEntry++)
            {
                CHECK(equal_to(tValuesA[iEntry], tValuesB[iEntry], 1e8));
                CHECK(tIndicesA[iEntry] == tIndicesB[iEntry]);
            }
        }

        moris::Matrix< DDRMat > tMorisMatrixB = convert_epetra_operator_to_arma_sp_mat< moris::Matrix, DDRMat >( B );
        moris::Matrix< DDRMat > tMorisMatrixA = convert_epetra_operator_to_arma_sp_mat< moris::Matrix, DDRMat >( A );

        if ( par_rank() == 1 )
        {
            return;
        }

        moris::Matrix< DDRMat > tExpectedMatrix = {{1.0, 2.0, 0.0},
                                                   {2.0, 0.0, 3.0},
                                                   {0.0, 4.0, 5.0}};

        CHECK_EQUAL( tMorisMatrixB, tExpectedMatrix , );
        CHECK_EQUAL( tMorisMatrixA, tExpectedMatrix ,  );
    }

}    // namespace moris

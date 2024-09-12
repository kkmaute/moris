#pragma once

#include "linalg_typedefs.hpp"

// include epetra headers
#include "Epetra_ConfigDefs.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"

/**
 * @brief converts an implicit epetra operator to a matrix
 *
 * @param aOperator
 * @param aMatrixToFill
 */
template< typename OperatorType >
Epetra_CrsMatrix
convert_epetra_operator_to_matrix(
        OperatorType *const aOperator,
        Epetra_Map const   &aMap )
{
    // create the identity matrix by glueing multiple unit multivectors tother
    Epetra_MultiVector tIdentity( aMap, aMap.NumGlobalPoints() );
    tIdentity.PutScalar( 0.0 );

    // get the number of GIDs on this processor
    int  tNumMyGIDs = aMap.NumMyPoints();
    int *tMyGIDs    = aMap.MyGlobalElements();

    // create the identity matrix
    for ( int iGID = 0; iGID < tNumMyGIDs; iGID++ )
    {
        tIdentity.ReplaceGlobalValue( tMyGIDs[ iGID ], tMyGIDs[ iGID ], 1.0 );
    }

    // Mutiple the identity multivectort with the matrix to get column of the matrix
    Epetra_MultiVector tOutputVector( tIdentity );
    aOperator->Apply( tIdentity, tOutputVector );

    // // create the row map for the matrix
    Epetra_CrsMatrix tReturnMatrix( Copy, aMap, 0 );

    // Now transfer the data from the multivector to the matrix with max size
    moris::Vector< int >    tColumnIndices( tIdentity.NumVectors() );
    moris::Vector< double > tNonZeroVals( tIdentity.NumVectors() );

    // insert the values into the matrix
    for ( int iRow = 0; iRow < tIdentity.MyLength(); iRow++ )
    {
        // counter to count non-zeros
        int tNnZCount = 0;
        for ( int iCol = 0; iCol < tIdentity.NumVectors(); iCol++ )
        {
            // This tolerance value is arbitrary
            if ( std::abs( tOutputVector[ iCol ][ iRow ] ) > 1e-10 )
            {
                tColumnIndices( tNnZCount ) = iCol;    // tFullMap.GID( iCol );
                tNonZeroVals( tNnZCount )   = tOutputVector[ iCol ][ iRow ];
                tNnZCount++;
            }
        }

        // for debugging
        tReturnMatrix.InsertGlobalValues( aMap.GID( iRow ), tNnZCount, tNonZeroVals.memptr(), tColumnIndices.memptr() );
    }

    tReturnMatrix.FillComplete();

    return tReturnMatrix;
}

// forward declare the write the serial matrix function
template< template< typename > class MatType, typename ElementType >
void write_serial_matrix( Epetra_CrsMatrix const &aEpetraMatrix, MatType< ElementType > &aArmaMatrix );

template< template< typename > class MatType, typename ElementType >
MatType< ElementType > convert_epetra_operator_to_arma_sp_mat( Epetra_CrsMatrix const &aEpetraMatrix )
{
    int                    tMatrixSize = aEpetraMatrix.NumGlobalRows();
    MatType< ElementType > tReturnMatrix;

    // different cases for different matrix types
    if constexpr ( std::is_same< MatType< ElementType >, moris::Matrix< ElementType > >::value )
    {
        tReturnMatrix = MatType< ElementType >( tMatrixSize, tMatrixSize, 0.0 );
    }
    else if constexpr ( std::is_same< MatType< ElementType >, arma::SpMat< ElementType > >::value )
    {
        tReturnMatrix = MatType< ElementType >( tMatrixSize, tMatrixSize );
        tReturnMatrix.mem_resize( aEpetraMatrix.NumGlobalNonzeros() );
    }

    // write a function that extract s
    if ( aEpetraMatrix.RowMap().Comm().NumProc() == 1 )
    {
        write_serial_matrix< MatType, ElementType >( aEpetraMatrix, tReturnMatrix );
        return tReturnMatrix;
    }
    else
    {
        // else case ;  break the matrix into chunks and transfer all of them into PE 0 and use the write_serial_matrix function
        const Epetra_Map &tRowMap    = aEpetraMatrix.RowMatrixRowMap();
        int               tNumMyRows = tRowMap.NumMyElements();

        // create a map on each processor containing all the GIDs on that PE
        Epetra_Map tAllGidsMap( -1, tNumMyRows, 0, tRowMap.Comm() );

        // create an epetra multivector with the GIDs,  parallel vector
        Epetra_IntVector tAllGids( tAllGidsMap );
        for ( int iRow = 0; iRow < tNumMyRows; iRow++ )
        {
            tAllGids[ iRow ] = tRowMap.GID( iRow );
        }

        // Now construct a RowMatrix on PE 0 by strip-mining the rows of the input matrix A.
        int numChunks    = tRowMap.Comm().NumProc();
        int stripSize    = tAllGids.GlobalLength() / numChunks;
        int remainder    = tAllGids.GlobalLength() % numChunks;
        int curStart     = 0;
        int curStripSize = 0;

        // an array that will store the GIDs that needs to be imported from other PEs
        Epetra_IntSerialDenseVector tImportedGIDList;

        // set the size of the
        if ( tRowMap.Comm().MyPID() == 0 )
        {
            tImportedGIDList.Size( stripSize + 1 );
        }

        for ( int iChunk = 0; iChunk < numChunks; iChunk++ )
        {
            if ( tRowMap.Comm().MyPID() == 0 )
            {    // Only PE 0 does this part
                curStripSize = stripSize;
                if ( iChunk < remainder ) curStripSize++;    // handle leftovers
                for ( int j = 0; j < curStripSize; j++ ) tImportedGIDList[ j ] = j + curStart;
                curStart += curStripSize;
            }

            if ( tRowMap.Comm().MyPID() > 0 ) assert( curStripSize == 0 );
            // The following import map will be non-trivial only on PE 0.
            Epetra_Map       tImportGIDMap( -1, curStripSize, tImportedGIDList.Values(), 0, tRowMap.Comm() );
            Epetra_Import    tGIDImporter( tImportGIDMap, tAllGidsMap );
            Epetra_IntVector tImportGids( tImportGIDMap );
            if ( tImportGids.Import( tAllGids, tGIDImporter, Insert ) != 0 ) { std::cout << "BAD" << '\n'; }

            // importGids now has a list of GIDs for the current strip of matrix rows.
            // Use these values to build another importer that will get rows of the matrix.

            // The following import map will be non-trivial only on PE 0.
            Epetra_Map       tImportMapMatrixMap( -1, tImportGids.MyLength(), tImportGids.Values(), 0, tRowMap.Comm() );
            Epetra_Import    tImporterMatrix( tImportMapMatrixMap, tRowMap );
            Epetra_CrsMatrix tImportedMatrix( Copy, tImportMapMatrixMap, 0 );
            if ( tImportedMatrix.Import( aEpetraMatrix, tImporterMatrix, Insert ) != 0 ) { std::cout << "BAD" << '\n'; }
            if ( tImportedMatrix.FillComplete( aEpetraMatrix.OperatorDomainMap(), tImportMapMatrixMap ) != 0 ) { std::cout << "BAD" << '\n'; }

            // Finally we are ready to write this strip of the matrix to ostream
            write_serial_matrix< MatType, ElementType >( tImportedMatrix, tReturnMatrix );
        }
    }
    return tReturnMatrix;
}

template< template< typename > class MatType, typename ElementType >
void write_serial_matrix(
        Epetra_CrsMatrix const &aEpetraMatrix,
        MatType< ElementType > &aArmaMatrix )
{

    int               tNumRows = aEpetraMatrix.NumMyRows();
    const Epetra_Map &tRowMap  = aEpetraMatrix.RowMap();
    const Epetra_Map &tColMap  = aEpetraMatrix.ColMap();

    // create arrays to store the values
    moris::real *tValues     = new moris::real[ aEpetraMatrix.MaxNumEntries() ];
    int         *tColIndices = new int[ aEpetraMatrix.MaxNumEntries() ];

    for ( int iRow = 0; iRow < tNumRows; iRow++ )
    {
        int I = tRowMap.GID( iRow );
        int tNonZeros;
        aEpetraMatrix.ExtractMyRowCopy( iRow, aEpetraMatrix.MaxNumEntries(), tNonZeros, tValues, tColIndices );

        // now loop over the columns that are non-zero and insert the i,j value into the matrix
        for ( int iCol = 0; iCol < tNonZeros; iCol++ )
        {
            int J               = tColMap.GID( tColIndices[ iCol ] );
            aArmaMatrix( I, J ) = tValues[ iCol ];
        }
    }

    // delete the pointers
    delete[] tValues;
    delete[] tColIndices;
}
// #endif

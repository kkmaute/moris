/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_load_matrix_from_binary_file.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_LOAD_MATRIX_FROM_BINARY_FILE_HPP_
#define PROJECTS_LINALG_SRC_FN_LOAD_MATRIX_FROM_BINARY_FILE_HPP_

#include <fstream>
#include <iostream>
#include <string>

// MORIS header files.
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
//--------------------------------------------------------------------------------

    template< typename Matrix_Type >
    void
    load_matrix_from_binary_file(
            Matrix< Matrix_Type > & aMatrix,
            const std::string     & aFilePath )
    {
        typedef typename Matrix< Matrix_Type >::Data_Type Type;

        // size of buffer in bit
        const uint tSizeOfBuffer = 512 * 1024 * 8;

        // samples in buffer to write
        const uint tNumberOfSamplesInBuffer = tSizeOfBuffer/sizeof( Type );

        // input file object
        std::ifstream tFile;

        // open input file
        tFile.open( aFilePath, std::ios::binary );

        // throw error if file can not be opened
        if( ! tFile )
        {
            std::string tMessage = "Could not open file " + aFilePath + " ." ;
            throw std::runtime_error( tMessage );
        }

        // number of rows of input matrix
        uint tNumberOfRows = 0;

        // number of cols of input matrix
        uint tNumberOfCols = 0;

        // write number of rows to file
        tFile.read( ( char* ) &tNumberOfRows, sizeof( uint ) );

        // write number of cols to file
        tFile.read( ( char* ) &tNumberOfCols, sizeof( uint ) );

        // initialize output matrix
        aMatrix.set_size( tNumberOfRows, tNumberOfCols );

        // create buffer
        std::vector<Type> tBuffer( tNumberOfSamplesInBuffer );

        // row index
        uint i = 0;

        // column index
        uint j = 0;

        // calculate number of full blocks
        uint tNumberOfBlocks = tNumberOfRows*tNumberOfCols / tNumberOfSamplesInBuffer ;

        // loop over all full blocks
        for ( uint b=0; b< tNumberOfBlocks; ++b )
        {
            tFile.read( ( char *) & tBuffer,  tSizeOfBuffer );

            for ( uint k=0; k< tNumberOfSamplesInBuffer; ++k )
            {
                // copy buffer to output matrix
                aMatrix( i, j ) = tBuffer[ k ];

                // increment row counter
                ++i;

                // go to next column
                if( i == tNumberOfRows )
                {
                    i = 0;
                    ++j;
                }
            }
        }

        // read final incomplete block
        uint tNumberOfSamples = tNumberOfRows*tNumberOfCols % tNumberOfSamplesInBuffer ;
        for ( uint k=0; k< tNumberOfSamples; ++k )
        {
            // read data and save into matrix
            tFile.read( ( char* ) &aMatrix( i, j ), sizeof( Type ) );

            // increment row counter
            ++i;

            // go to next column
            if( i == tNumberOfRows )
            {
                i = 0;
                ++j;
            }
        }

        // close input file
        tFile.close();
    }

//--------------------------------------------------------------------------------
}
#endif /* PROJECTS_LINALG_SRC_FN_LOAD_MATRIX_FROM_BINARY_FILE_HPP_ */


/*
 * fn_load_matrix_to_binary_file.hpp
 *
 *  Created on: Jul 18, 2018
 *      Author: messe
 */

#ifndef SRC_LINALG_FN_LOAD_MATRIX_FROM_BINARY_FILE_HPP_
#define SRC_LINALG_FN_LOAD_MATRIX_FROM_BINARY_FILE_HPP_

#include <iostream>
#include <string>

// MORIS header files.
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Mat.hpp" //LNA/src

namespace moris
{

//------------------------------------------------------------------------------
    /**
     * laods a matrix from a binary file
     * @param[ in ] a FilePath    string to file path
     * @param[ out ] aMatrix      matrix that is to be loaded
     */
    template< typename T >
    void
    load_matrix_from_binary_file(
            moris::Mat< T >    & aMatrix,
            const std::string  & aFilePath )
    {
        // size of buffer in bit
        const uint tSizeOfBuffer = 512 * 1024 * 8;

        // samples in buffer to write
        const uint tNumberOfSamplesInBuffer = tSizeOfBuffer/sizeof( T );

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
        T tBuffer[ tNumberOfSamplesInBuffer ];

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
            tFile.read( ( char* ) &aMatrix( i, j ), sizeof( T ) );

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

//------------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_LINALG_FN_LAOD_MATRIX_FROM_BINARY_FILE_HPP_ */

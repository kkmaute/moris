/*
 * fn_load_vector_to_binary_file.hpp
 *
 *  Created on: Apr 10, 2018
 *      Author: messe
 */

#ifndef SRC_LINALG_FN_LOAD_VECTOR_FROM_BINARY_FILE_HPP_
#define SRC_LINALG_FN_LOAD_VECTOR_FROM_BINARY_FILE_HPP_

#include <iostream>
#include <fstream>
#include <string>

// MORIS header files.

#include "assert.hpp"
#include "typedefs.hpp" // COR/src
#include "cl_Mat.hpp" // LNA/src


namespace moris
{
    template< typename T >
    void
    load_vector_from_binary_file(
        const std::string     & aFilePath,
              moris::Mat< T > & aVector )
    {

        // input file object
        std::ifstream tFile;

        // size of buffer in bit
        const moris::size_t tSizeOfBuffer = 512 * 1024 * 8;

        // samples in buffer to write
        const moris::size_t tNumberOfSamplesInBuffer = tSizeOfBuffer/sizeof( T );

        // open file
        tFile.open( aFilePath, std::ios::binary );

        // throw error if file can not be opened
        if( ! tFile )
        {
            std::string tMessage = "Could not open file " + aFilePath + " ." ;
            throw std::runtime_error(tMessage);
        }

        // calculate size of input file
        tFile.seekg( 0, std::ios::end );


        // get number of samples
        moris::size_t tNumberOfSamples = tFile.tellg()/sizeof( T );

        // jump to beginning of file
        tFile.seekg( 0, std::ios::beg );

        // assign output vector
        aVector.set_size( tNumberOfSamples, 1 );

        // create buffer
        //T tBuffer[ tNumberOfSamplesInBuffer ];
        T* tBuffer = new T[ tNumberOfSamplesInBuffer ];

        // calculate number of full blocks
        moris::size_t tNumberOfBlocks =  tNumberOfSamples / tNumberOfSamplesInBuffer ;

        // position in vector
        moris::size_t i=0;

        // loop over all full blocks
        for ( moris::size_t b=0; b< tNumberOfBlocks; ++b )
        {
            tFile.read( reinterpret_cast<char*> ( tBuffer ) ,  std::streamsize( tSizeOfBuffer ) );

            // loop over all values in buffer
            for ( moris::size_t k=0; k< tNumberOfSamplesInBuffer; ++k )
            {
                aVector( i++ ) = tBuffer[ k ];
            }
        }

        // read final incomplete block
        moris::size_t tNumberOfRemainingSamples = tNumberOfSamples % tNumberOfSamplesInBuffer ;
        for ( moris::size_t k=0; k< tNumberOfRemainingSamples; ++k )
        {
            // read data and save into matrix
            tFile.read( reinterpret_cast<char*> ( &aVector( i++ ) ), std::streamsize( sizeof( T ) ) );
        }

        // close file
        tFile.close();

        // free pointer to buffer
        delete [] tBuffer;
    }
}


#endif /* SRC_LINALG_FN_SAVE_VECTOR_TO_BINARY_FILE_HPP_ */

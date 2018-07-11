/*
 * fn_save_vector_to_binary_file.hpp
 *
 *  Created on: Apr 10, 2018
 *      Author: messe
 */

#ifndef SRC_LINALG_FN_SAVE_VECTOR_TO_BINARY_FILE_HPP_
#define SRC_LINALG_FN_SAVE_VECTOR_TO_BINARY_FILE_HPP_

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
    save_vector_to_binary_file(
        const std::string     & aFilePath,
        const moris::Mat< T > & aVector )
    {
        // size of buffer in bit
        const moris::size_t tSizeOfBuffer = 512 * 1024 * 8;

        // samples in buffer to write
        const moris::size_t tNumberOfSamplesInBuffer = tSizeOfBuffer / sizeof( T );

        // get length of vector
        moris::size_t tNumberOfSamples = aVector.length();

        // output file object
        std::ofstream tFile;

        // open output file
        tFile.open( aFilePath, std::ios::binary );

        // throw error if file can not be opened
        if( ! tFile )
        {
            std::string tMessage = "Could not create file " + aFilePath + " ." ;
            throw std::runtime_error(tMessage);
        }

        // create buffer
        T* tBuffer = new T[ tNumberOfSamplesInBuffer ];

        // counter for buffer data
        moris::size_t tCount = 0;

        for ( moris::size_t k=0; k< tNumberOfSamples; ++k )
        {
            // write data into buffer
            tBuffer[ tCount++ ] = aVector( k );

            // check if buffer is full
            if( tCount == tNumberOfSamplesInBuffer )
            {
                // write buffer to output file
                tFile.write( reinterpret_cast< char* >( tBuffer ), std::streamsize( tSizeOfBuffer ) );

                // reset counter
                tCount = 0;
            }
        }

        // write last incomplete buffer to output file
        if ( tCount != 0 )
        {
            for ( moris::size_t k=0; k<tCount; ++k )
            {
                tFile.write( reinterpret_cast< char* >( &tBuffer[ k ] ), std::streamsize( sizeof( T ) ) );
            }
        }

        // close output file
        tFile.close();

        // free pointer to buffer
        delete [] tBuffer;

    }
}


#endif /* SRC_LINALG_FN_SAVE_VECTOR_TO_BINARY_FILE_HPP_ */

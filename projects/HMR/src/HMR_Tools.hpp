/*
 * HMR_Tooms.hpp
 *
 *  Created on: May 14, 2018
 *      Author: messe
 */

#ifndef SRC_HMR_HMR_TOOLS_HPP_
#define SRC_HMR_HMR_TOOLS_HPP_
#include <string>
#include "cl_Communication_Tools.hpp" //COM/src
#include "typedefs.hpp" //COR/src
#include "cl_Map.hpp" //CON/src
#include "cl_Mat.hpp" //LNA/src
#include "fn_unique.hpp" //LNA/src

namespace moris
{
    namespace hmr
    {

// -----------------------------------------------------------------------------

    /**
     * This function inverts little endian to big endian and vice versa.
     * Needed for VTK debug files.
     */
    template <typename T> T swap_byte_endian(T aValue)
    {
        T aOutValue;
        auto *tPointer = (char*) &aValue;
        auto *tOutPointer = (char*)&aOutValue;
        int size = sizeof(T);
        for(int i=0; i<size; i++)
        {
            tOutPointer[size - 1 - i] = tPointer[i];
        }
        return aOutValue;
    }

// -----------------------------------------------------------------------------
    // print dots for nice output

    std::string
    proc_string()
    {
        std::string tString = "              ";

        if( par_size() > 1 )
        {
            uint tMyRank = par_rank();
            tString = "  proc " + std::to_string( tMyRank );

            if ( tMyRank < 10 )
            {
                tString +=" ... :" ;
            }
            else if ( tMyRank < 100 )
            {
                tString +=" .. :" ;
            }
            else if ( tMyRank < 1000 )
            {
                tString +=" . :" ;
            }
            else if ( tMyRank < 10000 )
            {
                tString +="  :" ;
            }
            else
            {
                tString +=" :" ;
            }
        }

        return tString;
    }

// -----------------------------------------------------------------------------

    template <typename T, typename U >
    void unique_with_matrix( Mat< T > & aVector, Mat< U > & aMatrix )
    {
        // backup input vector
        Mat< T > tVector = aVector;

        // backup input Matrix
        Mat< U > tMatrix = aMatrix;

        // make vector unique;
        aVector = unique( tVector );

        // get number of rows
        luint tRows = tMatrix.n_rows();

        // get number of cols
        luint tCols = tMatrix.n_cols();

        // get length of unique vector
        luint tLength = aVector.length();

        // reset size of output matrix
        aMatrix.set_size( tRows, tLength  );

        // create map
        map< T, T > tMap;
        for( uint j=0; j<tCols; ++j )
        {
            tMap[ tVector(j) ] = j;
        }

        // loop over all new entries
        for( luint k=0; k<tLength; ++k )
        {
            /* // compare with old entry
            for( luint j=0; j<tCols; ++j )
            {
                // test if equal
                if ( tVector( j ) == aVector( k ) )
                {
                    // copy column
                    aMatrix.cols( k, k ) = tMatrix.cols( j, j );

                    // exit loop
                    break;
                }
            } */

            // find position in old vector
            luint j = tMap.find( aVector( k ) );

            // copy column
            aMatrix.cols( k, k ) = tMatrix.cols( j, j );
        }
    }

// -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_HMR_TOOLS_HPP_ */

/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_Json_Object.cpp
 *
 */

#include "cl_Json_Object.hpp"
#include "cl_Matrix.hpp"
#include "moris_typedefs.hpp"

namespace moris
{
    void write_json( std::string const &aFileName, Json const &aTree )
    {
        std::ofstream tFile( aFileName );
        boost::property_tree::write_json( tFile, aTree );
        tFile.close();
    }

    std::string to_string( Json const &aTree )
    {
        std::ostringstream tStream;
        boost::property_tree::write_json( tStream, aTree );
        return tStream.str();
    }

    Json to_json( Json const &aTree )
    {
        return aTree;
    }

    Json to_json( Matrix< DDRMat > const &aMatrix )
    {
        uint tRows = aMatrix.n_rows();
        uint tCols = aMatrix.n_cols();

        bool isColumnVector = tCols == 1;
        bool isRowVector    = tRows == 1;

        Json tMatrix;
        for ( uint i = 0; i < tRows; ++i )
        {
            Json tRow;
            for ( uint j = 0; j < tCols; ++j )
            {
                Json tEntry;
                if ( isColumnVector )
                {
                    // treat the first column like a row
                    tEntry.put( "", aMatrix( j, i ) );
                }
                else
                {
                    tEntry.put( "", aMatrix( i, j ) );
                }
                tRow.push_back( { "", tEntry } );
            }
            if ( isRowVector || isColumnVector )
            {
                // for a vector, no nested array is needed
                return tRow;
            }
            tMatrix.push_back( { "", tRow } );
        }
        return tMatrix;
    }

    Json to_json( std::string const &aString )
    {
        Json tElement;
        tElement.put_value( aString );
        return tElement;
    }
}    // namespace moris

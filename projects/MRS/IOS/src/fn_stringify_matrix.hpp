/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_stringify_matrix.hpp
 *
 */

#ifndef PROJECTS_MRS_IOS_SRC_FN_STRINGIFY_MATRIX_HPP_
#define PROJECTS_MRS_IOS_SRC_FN_STRINGIFY_MATRIX_HPP_

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <limits>

#include "Log_Constants.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    namespace ios
    {
        // ----------------------------------------------------------------------------

        template< typename Matrix_Type >
        inline std::string
        stringify( const Matrix< Matrix_Type >& aMatrix )
        {
            std::ostringstream out;

            if ( aMatrix.numel() > 0 )
            {
                for ( uint iCol = 0; iCol < aMatrix.n_cols(); iCol++ )
                {
                    if ( iCol > 0 )
                    {
                        out << "; ";
                    }

                    for ( uint iRow = 0; iRow < aMatrix.n_rows() - 1; iRow++ )
                    {
                        out << aMatrix( iRow, iCol ) << ", ";
                    }

                    out << aMatrix( aMatrix.n_rows() - 1, iCol );
                }
            }

            return out.str();
        }

        // ----------------------------------------------------------------------------

        template< typename Matrix_Type >
        inline std::string
        stringify_log( const Matrix< Matrix_Type >& aMatrix )
        {
            // check matrix size being printed
            if ( aMatrix.numel() > LOGGER_MAX_NUMEL_MATRIX_PRINT )
            {
                return "[Matrix has too many elements to print.]";
            }

            // initialize string stream
            std::ostringstream out;
            out << "[";

            if ( aMatrix.numel() > 0 )
            {
                for ( uint iRow = 0; iRow < aMatrix.n_rows(); iRow++ )
                {
                    if ( iRow > 0 )
                    {
                        out << "; ";
                    }

                    for ( uint iCol = 0; iCol < aMatrix.n_cols(); iCol++ )
                    {
                        if ( iCol > 0 )
                        {
                            out << ", ";
                        }

                        out << aMatrix( iRow, iCol );
                    }
                }
            }

            // end matrix with square bracket
            out << "]";

            // return string stream as string
            return out.str();
        }

        // ----------------------------------------------------------------------------

        template<>
        inline std::string
        stringify_log( const Matrix< DDRMat >& aMatrix )
        {
            // check matrix size being printed
            if ( aMatrix.numel() > LOGGER_MAX_NUMEL_MATRIX_PRINT )
                return "[Matrix has too many elements to print.]";

            // initialize string stream
            std::ostringstream out;
            out << "[" << std::setprecision( LOGGER_FLOAT_PRECISION ) << std::scientific;

            if ( aMatrix.numel() > 0 )
            {
                for ( uint iRow = 0; iRow < aMatrix.n_rows(); iRow++ )
                {
                    if ( iRow > 0 )
                    {
                        out << " ; ";
                    }

                    for ( uint iCol = 0; iCol < aMatrix.n_cols(); iCol++ )
                    {
                        if ( iCol > 0 )
                        {
                            out << ", ";
                        }

                        out << aMatrix( iRow, iCol );
                    }
                }
            }

            // end matrix with square bracket
            out << "]";

            // return string stream as string
            return out.str();
        }

        // ----------------------------------------------------------------------------

        template<>
        inline std::string
        stringify_log( const Matrix< DDBMat >& aMatrix )
        {
            // check matrix size being printed
            if ( aMatrix.numel() > LOGGER_MAX_NUMEL_MATRIX_PRINT )
                return "[Matrix has too many elements to print.]";

            // initialize string stream
            std::ostringstream out;
            out << "[" << std::boolalpha;

            if ( aMatrix.numel() > 0 )
            {
                for ( uint iRow = 0; iRow < aMatrix.n_rows(); iRow++ )
                {
                    if ( iRow > 0 )
                    {
                        out << " ; ";
                    }

                    for ( uint iCol = 0; iCol < aMatrix.n_cols(); iCol++ )
                    {
                        if ( iCol > 0 )
                        {
                            out << ", ";
                        }

                        out << aMatrix( iRow, iCol );
                    }
                }
            }

            // end matrix with square bracket
            out << "]";

            // return string stream as string
            return out.str();
        }

        // ----------------------------------------------------------------------------

    }    // end namespace ios
}    // end namespace moris

#endif /* PROJECTS_MRS_IOS_SRC_FN_STRINGIFY_MATRIX_HPP_ */

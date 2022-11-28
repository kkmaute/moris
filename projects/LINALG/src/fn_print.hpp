/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_print.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_PRINT_HPP_
#define PROJECTS_LINALG_SRC_FN_PRINT_HPP_

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
#include <iomanip>      // std::setw

namespace moris
{

    //------------------------------------------------------------------------------

    template< typename Matrix_Type >
    inline void
    print( 
            Matrix< Matrix_Type > aMat,
            std::string           aTitle )
    {
        size_t tNumRows    = aMat.n_rows();
        size_t tNumColumns = aMat.n_cols();
        std::cout << "\n-------------------------------------------------\n";
        std::cout << aTitle + ": \n";
        std::cout << "Num Rows: " << tNumRows << " | Num Cols: " << tNumColumns << "\n";
        for ( size_t r = 0; r < tNumRows; r++ )
        {
            for ( size_t c = 0; c < tNumColumns; c++ )
            {
                std::cout << std::setw( 22 ) << aMat( r, c );
            }

            std::cout << "\n";
        }
        std::cout << "\n-------------------------------------------------\n" << std::endl;
    }

    //------------------------------------------------------------------------------

    template< typename Matrix_Type >
    inline void
    print_as_row_vector( 
            Matrix< Matrix_Type > aMat,
            std::string           aTitle )
    {
        size_t tNumRows    = aMat.n_rows();
        size_t tNumColumns = aMat.n_cols();

        MORIS_ERROR( tNumRows == 1 || tNumColumns == 1, 
                "print_as_row_vector() - Function should only be called on matrices with size 1 x n or n x 1." );

        std::cout << aTitle << " = [ " << std::flush;
        for ( moris::uint i = 0; i < aMat.numel(); i++ )
        {
            if ( i == aMat.numel() - 1 )
            {
                std::cout << aMat( i ) << std::flush;
            }
            else
            {
                std::cout << aMat( i ) << ", " << std::flush;
            }
        }

        std::cout << " ]" << std::endl;
    }

    //------------------------------------------------------------------------------

    /*
     * prints the col and row index as well as the data
     */
    //------------------------------------------------------------------------------
        template<typename Matrix_Type>
        inline void
        print_fancy(Matrix< Matrix_Type > aMat,
                    std::string aTitle)
        {
            size_t tNumRows = aMat.n_rows();
            size_t tNumColumns = aMat.n_cols();
            std::cout << "\n-------------------------------------------------\n";
            std::cout << aTitle + ": \n";
            std::cout << "Num Rows: " << tNumRows << " | Num Cols: " << tNumColumns << "\n";

            std::cout<<std::right<<std::setw(5)<<"";
            for(size_t c = 0; c < tNumColumns; c++)
            {
                std::cout<<std::right<<std::setw(12)<<c;
            }

            std::cout<<"\n";
            std::cout<<"--------";

            for(size_t c = 0; c < tNumColumns+1; c++)
            {
                std::cout<<"-------------";
            }

            std::cout<<"\n";

            for(size_t r = 0; r < tNumRows; r++)
            {
                std::cout<<std::right<<std::setw(5)<<r<<" | ";

                for(size_t c = 0; c < tNumColumns; c++)
                {
                    std::cout<< std::setw(12) << aMat(r, c);
                }

                std::cout << "\n";
            }
            std::cout << "\n-------------------------------------------------\n"<<std::endl;
        }

    //------------------------------------------------------------------------------

    inline
    void
    print(Matrix< DDRMat > const & aMat,
          std::string aTitle)
    {
        FILE * outFile = stdout;

        fprintf( outFile,"%s-------------------------------------------------\n\n","%" );

        if ( aTitle.empty() )
            fprintf( outFile, "morisMat = [ ... \n" );
        else
            fprintf( outFile, "%s = [ ... \n", aTitle.c_str() );

        for( uint ir = 0; ir < aMat.n_rows(); ++ir )
        {
            for( uint ic = 0; ic < aMat.n_cols(); ++ic )
            {
                // FIXME: need to type cast the output of mMat(ir,ic)
                fprintf( outFile, "%+.15e", ( real ) aMat(ir,ic) );

                if(ic < aMat.n_cols() -1)
                    fprintf( outFile, ",  " );
                else
                    fprintf( outFile, ";" );
            }

            if (ir < aMat.n_rows() - 1 )
                fprintf( outFile, ".... \n" );
            else
                fprintf( outFile, "];\n" );
        }
        fprintf( outFile,"%s-------------------------------------------------\n\n", "%") ;
    }

    //------------------------------------------------------------------------------

    template<typename Matrix_Type>
    inline void print_std_initializer_list(Matrix< Matrix_Type > const & aMatrix,
                                    std::string aTitle)
    {
        std::cout << "\n" << aTitle + ": " << std::endl;
        size_t tNumRows = aMatrix.n_rows();
        size_t tNumColumns = aMatrix.n_cols();

        std::cout << "{";
        for(size_t r = 0; r < tNumRows; r++)
        {
            std::cout << "{";
            for(size_t c = 0; c < tNumColumns; c++)
            {
    //            std::cout << std::setprecision(std::numeric_limits<Type>::digits10 + 1) << aMatrix(r, c);
                std::cout << aMatrix(r, c);
                if(c != tNumColumns - 1)
                {
                    std::cout << ", ";
                }
            }
            std::cout << "}";
            if(r != tNumRows - 1)
            {
                std::cout << ", \n";
            }

        }

        std::cout << "}\n" << std::endl;
    }

    /*
     * Cell of matrices
     */
    template< typename T >
    inline void
    print(moris::Cell< moris::Matrix<T> > const & aCell,
          std::string aStr = "Cell")
    {
        std::cout<<"Cell Name: "<<aStr<<"\n";
        std::cout<<"Number of entries = "<<aCell.size()<<"\n";
        for(moris::uint  i = 0; i <aCell.size(); i++)
        {
            moris::print(aCell(i),"Cell " + std::to_string(i));
        }

        std::cout<<std::endl;
    }

    /*
     * Cell of matrices
     */
    template< typename T >
    inline void
    print(moris::Cell< moris::Matrix<T> * > const & aCell,
          std::string aStr = "Cell")
    {
        std::cout<<"Cell Name: "<<aStr<<"\n";
        std::cout<<"Number of entries = "<<aCell.size()<<"\n";
        for(moris::uint  i = 0; i <aCell.size(); i++)
        {
            moris::print(*aCell(i),"Cell " + std::to_string(i));
        }

        std::cout<<std::endl;
    }

} /* namespace moris */

#endif /* PROJECTS_LINALG_SRC_FN_PRINT_HPP_ */


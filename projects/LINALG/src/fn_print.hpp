/*
 * fn_print.hpp
 *
 *  Created on: Aug 27, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_FN_PRINT_HPP_
#define PROJECTS_LINALG_SRC_FN_PRINT_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"


namespace moris
{

    template<typename Matrix_Type>
    void
    print(Matrix< Matrix_Type > aMat,
          std::string aTitle)
    {
        size_t tNumRows = aMat.n_rows();
        size_t tNumColumns = aMat.n_cols();

        std::cout << aTitle + ": " << std::endl;
        std::cout << "Num Rows: " << tNumRows << " | Num Cols: " << tNumColumns << std::endl;
        for(size_t r = 0; r < tNumRows; r++)
        {
            for(size_t c = 0; c < tNumColumns; c++)
            {
                std::cout << aMat(r, c) << " ";
            }

            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

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

        for( moris::uint ir = 0; ir < aMat.n_rows(); ++ir )
        {
            for( moris::uint ic = 0; ic < aMat.n_cols(); ++ic )
            {
                // FIXME: need to type cast the output of mMat(ir,ic)
                fprintf( outFile, "%+.15e", (moris::real)aMat(ir,ic) );

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
}


#endif /* PROJECTS_LINALG_SRC_FN_PRINT_HPP_ */

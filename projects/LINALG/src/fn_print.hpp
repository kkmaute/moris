/*
 * fn_print.hpp
 *
 *  Created on: Aug 27, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_FN_PRINT_HPP_
#define PROJECTS_LINALG_SRC_FN_PRINT_HPP_

#include "cl_Matrix.hpp"


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
}


#endif /* PROJECTS_LINALG_SRC_FN_PRINT_HPP_ */

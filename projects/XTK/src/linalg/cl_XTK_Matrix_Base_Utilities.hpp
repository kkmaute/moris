/*
 * cl_XTK_Matrix_Tools.hpp
 *
 *  Created on: Jun 20, 2017
 *      Author: ktdoble
 */

#ifndef SRC_LINALG_CL_XTK_MATRIXBASE_UTILITIES_HPP_
#define SRC_LINALG_CL_XTK_MATRIXBASE_UTILITIES_HPP_

#include <mpi.h>
#include <iostream>

#include "tools/cl_MPI_Tools.hpp"
#include "cl_Matrix.hpp"
#include "tools/fn_bubble_sort.hpp"
#include "tools/fn_approximate.hpp"

#include "assert/fn_xtk_assert.hpp"

namespace xtk
{
template<typename Type, typename Matrix_Type>
void print(moris::Matrix<Type, Matrix_Type> const & aMatrix, std::string aTitle)
{
    int tProcRank = 0;
    MPI_Comm_rank(get_comm(), &tProcRank);

    size_t tNumRows = aMatrix.n_rows();
    size_t tNumColumns = aMatrix.n_cols();

    std::cout << aTitle + ": " << std::endl;
    std::cout << "Num Rows: " << tNumRows << " | Num Cols: " << tNumColumns << " | Proc Rank: " << tProcRank << std::endl;
    for(size_t r = 0; r < tNumRows; r++)
    {
        for(size_t c = 0; c < tNumColumns; c++)
        {
            std::cout << aMatrix(r, c) << " ";
        }

        std::cout << std::endl;
    }
    std::cout << std::endl;
}




template<typename Type, typename Matrix_Type>
void print_std_initializer_list(moris::Matrix<Type,Matrix_Type> const & aMatrix, std::string aTitle)
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
            std::cout << ", ";
        }

    }

    std::cout << "}\n" << std::endl;
}

template<typename Type, typename Matrix_Type>
void print_to_matlab(moris::Matrix<Type,Matrix_Type> const & aMatrix,
                          std::string const & aVarName)
{
    FILE * outFile = stdout;

     fprintf( outFile,"%s-------------------------------------------------\n\n","%" );

     if ( aVarName.empty() )
         fprintf( outFile, "morisMat = [ ... \n" );
     else
         fprintf( outFile, "%s = [ ... \n", aVarName.c_str() );

     for( size_t ir = 0; ir < aMatrix.n_rows(); ++ir )
     {
         for( size_t ic = 0; ic < aMatrix.n_cols(); ++ic )
         {
             // FIXME: need to type cast the output of mMat(ir,ic)
             fprintf( outFile, "%+.15e", (real)aMatrix(ir,ic) );

             if(ic < aMatrix.n_cols() -1)
                 fprintf( outFile, ",  " );
             else
                 fprintf( outFile, ";" );
         }

         if (ir < aMatrix.n_rows() - 1 )
             fprintf( outFile, ".... \n" );
         else
             fprintf( outFile, "]\n" );
     }
     fprintf( outFile,"%s-------------------------------------------------\n\n", "%") ;
}


template<typename Type,typename Matrix_Type>
bool equal_to(moris::Matrix<Type, Matrix_Type> const & aMatrix1,
              moris::Matrix<Type, Matrix_Type> const & aMatrix2,
              bool aCheckRowsandCols = true)
{
    bool tFlag = true;

    size_t tNumRows1 = aMatrix1.n_rows();
    size_t tNumCols1 = aMatrix1.n_cols();

    size_t tNumRows2 = aMatrix2.n_rows();
    size_t tNumCols2 = aMatrix2.n_cols();

    if(aCheckRowsandCols)
    {
        if(tNumRows1 != tNumRows2)
        {
            XTK_ERROR << "There are a different number of rows in the provided matrices";
            tFlag = false;
        }

        if(tNumCols1 != tNumCols2)
        {
            XTK_ERROR << "There are a different number of columns in the provided matrices";
            tFlag = false;
        }
    }

    for(size_t r = 0; r < tNumRows1; r++)
    {
        for(size_t c = 0; c < tNumCols1; c++)
        {
            if(!approximate(aMatrix2(r,c),aMatrix1(r,c)))
            {
                tFlag = false;
                break;
            }
        }
    }
    return tFlag;
}


template<typename Type,typename Matrix_Type>
bool row_equal(size_t const & aRowIndex1,
               moris::Matrix<Type,Matrix_Type> const & aMatrix1,
               size_t const & aRowIndex2,
               moris::Matrix<Type,Matrix_Type> const & aMatrix2)
{
    bool tEqual = true;
    size_t tNumCols1 = aMatrix1.n_cols();

    if(tNumCols1 > 0)
    {
        for(size_t c = 0; c < tNumCols1; c++)
        {
            if(!approximate(aMatrix1(aRowIndex1,c),aMatrix2(aRowIndex2,c)))
            {
                tEqual =  false;
                break;
            }
        }
    }

    return tEqual;
}


/*
*
* True if there are no duplicate faces
* Order in the row does not matter (sorts in ascending order)
*/
template<typename Type,typename Matrix_Type>
bool
check_for_duplicate_rows(moris::Matrix<Type,Matrix_Type> & aMat,
                         bool aOrderMatters = true)
{

    bool tSame = true;
    size_t tNumCols = aMat.n_cols();
    moris::Matrix<Type,Matrix_Type> tRow1(1,tNumCols);
    moris::Matrix<Type,Matrix_Type> tRow2(1,tNumCols);

    // Sort the rows in ascending order if the order matters
    if(!aOrderMatters){  xtk::row_bubble_sort(aMat); }


    for(size_t i = 0; i<aMat.n_rows(); i++)
    {
        tRow1 = aMat.get_row(i);

        for(size_t j = i+1; j<aMat.n_rows(); j++)
        {
            tRow2 = aMat.get_row(j);

            tSame = xtk::equal_to(tRow1,tRow2);

            if(tSame)
            {
                std::cout<< "i = " << i<< " j = "<< j<<std::endl;
                return false;
            }
        }
    }

    return true;
}

template<typename Type,typename Matrix_Type>
bool
check_for_duplicate_columns(moris::Matrix<Type,Matrix_Type> & aMat,
                            bool aOrderMatters = true)
{

    bool tSame = true;
    size_t tNumRows = aMat.n_rows();
    std::shared_ptr<moris::Matrix<Type,Matrix_Type> > tCol1 = aMat.create(tNumRows,1);
    moris::Matrix<Type,Matrix_Type> tCol2 = aMat.create(tNumRows,1);

    // Sort the rows in ascending order if the order matters
    if(!aOrderMatters){  xtk::row_bubble_sort(aMat); }


    for(size_t i = 0; i<aMat.n_cols(); i++)
    {
        aMat.get_column(i,*tCol1);

        for(size_t j = i+1; j<aMat.n_cols(); j++)
        {
            aMat.get_column(j,*tCol2);

            tSame = xtk::equal_to(*tCol1,*tCol2);

            if(tSame)
            {
                return false;
            }
        }
    }

    return true;
}

/*
 * Sets aMatrix.row(aRowIndex1) as aMatrix2.row(aRowIndex2)
 */
template<typename Type,typename Matrix_Type>
void
replace_row(size_t const & aRowIndex1,
            moris::Matrix<Type,Matrix_Type> const & aMatrix1,
            size_t const & aRowIndex2,
            moris::Matrix<Type,Matrix_Type> & aMatrix2,
            bool aCheckRowsandCols = false)
{
    size_t tNumCols1 = aMatrix1.n_cols();

    if(aCheckRowsandCols)
    {

        size_t tNumRows1 = aMatrix1.n_rows();
        size_t tNumRows2 = aMatrix2.n_rows();
        size_t tNumCols2 = aMatrix2.n_cols();

        XTK_ASSERT(tNumRows1 = tNumRows2, "Different number of rows");
        XTK_ASSERT(tNumCols1 = tNumCols2, "Different number of cols");
        XTK_ASSERT(aRowIndex1<tNumRows1,"Row index out of bounds matrix 1");
        XTK_ASSERT(aRowIndex2<tNumRows2,"Row index out of bounds matrix 2");
    }

    for(size_t i = 0; i<tNumCols1; i++)
    {
        aMatrix2(aRowIndex2,i) = aMatrix1(aRowIndex1,i);
    }
}


template<typename Type,typename Matrix_Type>
void
fill_row(Type const & aFillValue,
         size_t const & aRowIndex,
         moris::Matrix<Type,Matrix_Type> & aMatrix1)
{
    size_t tNumCols1 = aMatrix1.n_cols();

    for(size_t i = 0; i<tNumCols1; i++)
    {
        aMatrix1(aRowIndex,i) = aFillValue;
    }
}

/*
 * Fill Matrix with another matrix without resize
 *
 * Matrix 1 goes into Matrix 2
 */
template<typename Type,typename Matrix_Type>
void conservative_copy(moris::Matrix<Type,Matrix_Type> const & aMatrix1,
                       moris::Matrix<Type,Matrix_Type>  & aMatrix2)
{
    size_t tNumRow1 = aMatrix1.n_rows();
    size_t tNumCol1 = aMatrix1.n_cols();

    for(size_t iR = 0; iR<tNumRow1; iR++)
    {
       for( size_t iC = 0; iC<tNumCol1; iC++ )
        {
            aMatrix2(iR,iC) = aMatrix1(iR,iC);
        }
    }
}


/*
 * Similar to MATLAB operation
 * aMatToReindex(aMatRow,:)(aIndexMap)
 *
 * I.e
 * aMatToReindex = [2,3,4,5;
 *                  1,2,3,4;
 *                  4,2,1,1];
 * aMatRow = 0;
 * aIndexMap = [0,1,2,3;
 *              3,2,1,0];
 *
 * would result in
 *
 * [2,3,4,5;
 *  5,4,3,2];
 *
 */


template<typename Type,typename Matrix_Type, typename Integer_Type, typename Integer_Matrix>
moris::Matrix<Type,Matrix_Type>
reindex_matrix(moris::Matrix<Integer_Type,Integer_Matrix> const & aIndexMap,
               size_t                        const & aMatRow,
               moris::Matrix<Type,Matrix_Type> const & aMatToReindex)
{
    size_t tNumRow = aIndexMap.n_rows();
    size_t tNumCol = aIndexMap.n_cols();
    moris::Matrix<Type,Matrix_Type> tReindexedMat(tNumRow,tNumCol);

    for( size_t i = 0; i<tNumRow; i++)
    {
        for(size_t j =0; j<tNumCol; j++)
        {
            tReindexedMat(i,j) = aMatToReindex(aMatRow,aIndexMap(i,j));
        }
    }

    return tReindexedMat;
}
/*
 *
 */





}

#endif /* SRC_LINALG_CL_XTK_MATRIXBASE_UTILITIES_HPP_ */

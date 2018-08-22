/*
 * fn_matrix_to_offset.hpp
 *
 *  Created on: Apr 23, 2018
 *      Author: ktdoble
 */

#ifndef SRC_TOOLS_FN_MATRIX_TO_OFFSET_HPP_
#define SRC_TOOLS_FN_MATRIX_TO_OFFSET_HPP_


namespace xtk
{

/*
 * Converts the matrix to offsets and indicies format. If there is a dummy value in the matrix indicated an empty space
 * this can be represented through the Dummy inputs. This format is more compact then a sparsely populated matrix
 */
template< typename Type, typename Matrix_Type>
void
convert_matrix_to_offsets( Matrix_Base<Type, Matrix_Type> const & aMatrix,
                           Matrix_Base<Type, Matrix_Type> & aIndices,
                           Matrix_Base<Type, Matrix_Type> & aOffsets,
                           bool aHasDummyVal = false,
                           Type aDummyVal = std::numeric_limits<Type>::max())
{
    bool tDummyFlag = false;
    size_t tDummyCountIndices = 0;
    size_t tNumCols = aMatrix.get_num_columns();
    size_t tNumRows = aMatrix.get_num_rows();
    size_t tSize = tNumCols*tNumRows;
    aIndices.resize(1,tSize);
    aOffsets.resize(1,tNumRows+1);


    aOffsets(0,0) = 0;
    size_t tCount = 0;
    for(size_t i = 0; i<tNumRows; i++)
    {
        for(size_t j = 0; j <tNumCols; j++)
        {
            if(aHasDummyVal && aMatrix(i,j) ==aDummyVal)
            {
                tDummyFlag = true;
                tDummyCountIndices++;
            }

            if(!tDummyFlag)
            {
                aIndices(0,tCount) = aMatrix(i,j);
                tCount++;
            }

            tDummyFlag = false;
        }
        aOffsets(0,i+1) = tCount;
    }

    aIndices.resize(1,tSize-tDummyCountIndices);
    aOffsets.resize(1,tNumRows+1);
}

/**
 * Converts from the mesh storage structure to the matrix
 * This assumes a uniform offset
 */
template<typename Type, typename Matrix_Type>
void
convert_offsets_to_matrix(Matrix_Base<Type, Matrix_Type> const & aIndices,
                          Matrix_Base<Type, Matrix_Type> const & aOffsets,
                          Matrix_Base<Type, Matrix_Type> & aMatrix)
{
    size_t tNumCols = aOffsets(0,1) - aOffsets(0,0);
    size_t tNumRows;
    if(tNumCols!=0)
    {
        tNumRows = aIndices.get_num_columns()/tNumCols;
    }
    aMatrix.resize(tNumRows,tNumCols);
    aMatrix.fill(std::numeric_limits<Type>::max());

    Type tCount = 0;
    for(size_t i = 0; i<tNumRows; i++)
    {
        for(size_t j = 0; j <tNumCols; j++)
        {

            aMatrix(i,j) = aIndices(0,tCount);
            tCount++;
        }
  }

}


}


#endif /* SRC_TOOLS_FN_MATRIX_TO_OFFSET_HPP_ */

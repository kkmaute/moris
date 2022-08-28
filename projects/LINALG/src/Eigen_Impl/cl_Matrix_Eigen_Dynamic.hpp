/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Matrix_Eigen_Dynamic.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_DYNAMIC_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_DYNAMIC_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"
#include "fn_iscol.hpp"
#include "fn_isvector.hpp"

namespace moris
{
    template<typename Type>
    class Matrix<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>>
    {
        private:
            Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> mMatrix;

            void fill_with_NANs()
            {
#ifdef MATRIX_FILL
                if (std::numeric_limits<Data_Type>::has_quiet_NaN)
                {
                    mMatrix.fill( std::numeric_limits<Data_Type>::quiet_NaN() );
                }
                else
                {
                    mMatrix.fill( std::numeric_limits<Data_Type>::max() );
                }
#endif
            }

        public:
            typedef Type Data_Type;

            Matrix()
            {

            };

            Matrix(size_t const & aNumRows,
                    size_t const & aNumCols):
                        mMatrix(aNumRows,aNumCols)
            {
                this->fill_with_NANs();
            }

            // template constructor
            Matrix(Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> const & X ):
                mMatrix(X)
            {

            }

            Matrix( size_t const & aNumEl):
                mMatrix(aNumEl,1)
            {
                this->fill_with_NANs();
            }

            // template constructor
            template< typename A >
            Matrix(A const & X ):
            mMatrix(X)
            {

            }

            Matrix( size_t const & aNumRows,
                    size_t const & aNumCols,
                    Type   const & aFillVal):
                        mMatrix( aNumRows, aNumCols )
            {
                mMatrix.fill(aFillVal);
            }

            Matrix( Type*          & aArray,
                    moris::size_t const & aNumRows,
                    moris::size_t const & aNumCols):
                        mMatrix(aNumRows,aNumCols)
            {
                this->mMatrix = Eigen::Map< Eigen::Matrix< Type, Eigen::Dynamic, Eigen::Dynamic > >(aArray, aNumRows, aNumCols);
            }

            Matrix(std::initializer_list<std::initializer_list<Type> > const & aInitList)
            {
                size_t i = 0;
                size_t j = 0;
                size_t tNumRows = aInitList.size();
                size_t tNumColumns = aInitList.begin()->size();

                mMatrix = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>(tNumRows, tNumColumns);

                for(const auto tRow : aInitList) // loop over number of rows
                {
                    MORIS_ASSERT(tRow.size() == aInitList.begin()->size(),
                            "The number of elements in one of the rows does not equal the number of columns.");

                    for(const auto tCol : tRow) // loop over every value in the row
                    {
                        mMatrix(i, j) = tCol;
                        ++j;
                    }
                    j = 0;
                    ++i;
                }
            }

            // Copy operations
            Matrix<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>>
            copy() const
            {
                Matrix<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> tMatCopy(this->n_rows(),this->n_cols());
                tMatCopy.matrix_data() = mMatrix;
                return tMatCopy;
            }

            void
            resize(const size_t & aNumRows,
                    const size_t & aNumCols)
            {
                mMatrix.conservativeResize(aNumRows, aNumCols);
            }

            void
            set_size(const size_t & aNumRows,
                    const size_t & aNumCols)
            {
                mMatrix.resize( aNumRows, aNumCols );

                this->fill_with_NANs();
            }

            void
            set_size(const size_t & aNumRows,
                    const size_t & aNumCols,
                    const Type   & aFillValue )
            {
                mMatrix.resize( aNumRows, aNumCols );
                mMatrix.fill( aFillValue );
            }

            void
            reshape(const size_t & aNumRows,
                    const size_t & aNumCols)
            {
                Eigen::Map<Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> tMatrix(mMatrix.data(),aNumRows,aNumCols);

                mMatrix=tMatrix;
            }

            void
            inplace_trans()
            {
                mMatrix.transposeInPlace();
            }

            void
            fill(const Type & aFillValue)
            {
                mMatrix.fill( aFillValue );
            }

            /**
             * Get the number of columns in a data set, similar to Matlab cols().
             *
             * @return Number of columns.
             */
            size_t
            n_cols() const
            {
                return mMatrix.cols();
            }

            /**
             * Get the number of rows in a data set, similar to Matlab rows().
             *
             * @return Number of rows.
             */
            size_t
            n_rows() const
            {
                return mMatrix.rows();
            }

            size_t
            size(size_t aDim)
            {
                if(aDim == 0)
                {
                    return this->n_rows();
                }
                else if(aDim == 1)
                {
                    return this->n_cols();
                }
                else
                {
                    MORIS_ASSERT(false,"Invalid matrix dimension specified, 0-for n_rows, 1- for n_cols");
                    return 0;
                }
            }

            /**
             * Returns the number of elements in the %matrix.
             *
             * @return Number of elements in the %matrix.
             *
             */
            size_t
            numel() const
            {
                return mMatrix.size();
            }

            void
            set_row(size_t aRowIndex,
                    Matrix< Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> const & aVec)
            {
                MORIS_ASSERT(isvector(aVec), "aVec needs to be a vector");
                MORIS_ASSERT(aRowIndex < this->n_rows(), "Specified row index out of bounds");
                MORIS_ASSERT(aVec.numel() == this->n_cols(),
                        "Dimension mismatch (argument matrix and member matrix do not have same number of columns)");

                size_t tROW_INDEX = 0;

                if(!iscol(aVec))
                {
                    mMatrix.row(aRowIndex) = aVec.matrix_data().row(tROW_INDEX);
                }
                else
                {
                    mMatrix.row(aRowIndex) = aVec.matrix_data().col(tROW_INDEX);
                }
            }

            void set_column(size_t aColumnIndex, const Matrix< Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> & aColumn)
            {

                MORIS_ASSERT(aColumn.n_cols() == 1, "aColumn needs to be a column matrix");
                MORIS_ASSERT(aColumnIndex < this->n_cols(), "Specified column index out of bounds");
                MORIS_ASSERT(aColumn.n_rows() == this->n_rows(),
                        "Dimension mismatch (argument matrix and member matrix do not have same number of rows)");

                size_t tCOLUMN_INDEX = 0;
                mMatrix.col(aColumnIndex) = aColumn.matrix_data().col(tCOLUMN_INDEX);
            }

            void
            get_column(size_t aColumnIndex,
                    Matrix< Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> & aColumn) const
            {
                MORIS_ASSERT(aColumn.n_cols() == 1,"aColumn needs to be a column matrix");
                MORIS_ASSERT(aColumnIndex < this->n_cols(),"Specified column index out of bounds");
                MORIS_ASSERT(aColumn.n_rows() == this->n_rows(),"Dimension mismatch (argument matrix and member matrix do not have same number of rows)");
                const size_t tCOLUMN_INDEX = 0;
                aColumn.matrix_data().col(tCOLUMN_INDEX) = mMatrix.col(aColumnIndex);
            }

            Matrix< Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>>
            get_column(size_t aColumnIndex) const
            {
                MORIS_ASSERT(aColumnIndex < this->n_cols(),"Specified column index out of bounds");
                return mMatrix.col(aColumnIndex);
            }

            void get_row(size_t aRowIndex, Matrix< Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> & aRow) const
            {
                MORIS_ASSERT(aRow.n_rows() == 1,"aRow needs to be a row matrix");
                MORIS_ASSERT(aRowIndex < this->n_rows(),"Specified row index out of bounds");
                MORIS_ASSERT(aRow.n_cols() == this->n_cols(),"Dimension mismatch (argument matrix and member matrix do not have same number of columns)");

                const size_t tROW_INDEX = 0;
                aRow.mMatrix.row(tROW_INDEX) = mMatrix.row(aRowIndex);
            }

            auto
            get_row(size_t aRowIndex) const
            ->decltype(mMatrix.row(aRowIndex) )
            {
                MORIS_ASSERT(aRowIndex < this->n_rows(),"Specified row index out of bounds");
                return mMatrix.row(aRowIndex);
            }

            auto
            get_row(size_t aRowIndex)
            ->decltype(mMatrix.row(aRowIndex) )
            {
                MORIS_ASSERT(aRowIndex < this->n_rows(),"Specified row index out of bounds");
                return mMatrix.row(aRowIndex);
            }

            const Type*
            data() const
            {
                return mMatrix.data();
            }

            Type*
            data()
            {
                return mMatrix.data();
            }

            inline
            Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> &
            matrix_data()
            {
                return mMatrix;
            }

            inline
            Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> const &
            matrix_data() const
            {
                return mMatrix;
            }

            Type
            max() const
            {
                return mMatrix.maxCoeff( );
            }

            Type
            min() const
            {
                return mMatrix.minCoeff( ) ;
            }

            Type
            max(uint & aRowIndex,
                    uint & aColIndex) const
            {
                auto val = this->mMatrix.maxCoeff( &aRowIndex, &aColIndex );

                return val;
            }

            Type
            min(uint & aRowIndex,
                    uint & aColIndex) const
            {
                auto val = this->mMatrix.minCoeff( &aRowIndex, &aColIndex );

                return val;
            }

            /**
             * @brief Overloaded moris::Matrix_Base::operator()
             *
             * @param[in] aRowIndex Row index for which data should be accessed.
             * @param[in] aColIndex Column index for which data should be accessed.
             */
            inline
            Type &
            operator()( size_t const & aRowIndex,
                    size_t const & aColIndex )
            {
                return mMatrix(aRowIndex,aColIndex);
            }

            /**
             * @brief Overloaded moris::Matrix_Base::operator()
             *
             * @param[in] aRowIndex Row index for which data should be accessed.
             * @param[in] aColIndex Column index for which data should be accessed.
             */
            const Type &
            operator()(const size_t & aRowIndex,
                    const size_t & aColIndex) const
            {
                return mMatrix(aRowIndex,aColIndex);
            }

            /**
             * @brief Overloaded moris::Matrix_Base::operator()
             *
             * @param[in] aIndex Index for which data should be accessed.
             */
            inline
            Type &
            operator()( const size_t & aIndex )
            {
                return mMatrix( aIndex );
            }

            /**
             * @brief Overloaded moris::Matrix_Base::operator()
             *
             * @param[in] aIndex Index for which data should be accessed.
             */
            const Type &
            operator()( const size_t & aIndex ) const
            {
                return mMatrix( aIndex );
            }

            /*
             * Block operations
             */
            auto
            operator()(
                    std::pair< moris::size_t, moris::size_t > const & aI,
                    std::pair< moris::size_t, moris::size_t > const & aJ )
            ->decltype(mMatrix.block( aI.first, aJ.first, aI.second-aI.first+1,aJ.second-aJ.first+1 ) )
            {
                return mMatrix.block( aI.first, aJ.first, aI.second-aI.first+1,aJ.second-aJ.first+1 ) ;
            }

            /**
             * Assignment operator.
             *
             * @param[in] X Matrix or column vector or row vector.
             *
             * @return Assignment.
             */
            Matrix< Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > &
            operator=( Matrix< Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > const & X )
            {
                mMatrix = X.matrix_data();
                return *this;
            }

            /**
             * Returns the length of a vector. Thows error neither rows nor cols are equal 1.
             */
            size_t
            length() const
            {
                // get number of rows from matrix implementation
                size_t n_rows = this->n_rows();

                // get number of cols from matrix implementation
                size_t n_cols = this->n_cols();

                // catch special case of zero length
                if( n_rows == 0 || n_cols == 0 )
                {
                    return 0;
                }
                else
                {
                    // assert that this is really a vector
                    MORIS_ASSERT(  n_rows == 1 || n_cols == 1,
                            "Tried to get length of a matrix. Check dimensions." );

                    // return the smaller of both values
                    return ( n_rows < n_cols ) ? n_cols : n_rows;
                }
            }
    };
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_DYNAMIC_HPP_ */


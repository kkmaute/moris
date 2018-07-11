#ifndef MORIS_LINALG_CL_BASE_MAT_HPP_
#define MORIS_LINALG_CL_BASE_MAT_HPP_

// C++ header files.
#include <utility>

// MORIS library header files.
#include "typedefs.hpp"
#include "cl_Tuple.hpp"
#include "ios.hpp"

// Class forward declarations.
namespace moris {
    template< typename T >
    class Base_Mat;

    template< typename T >
    class Mat;

    template< typename T >
    class Sp_Mat;
}

/**
 * @brief The base matrix class.
 *
 * For dense matrices, see moris::Mat; for sparse matrices, see moris::Sp_Mat.
 */
template< typename T >
class moris::Base_Mat
{

protected:
    /**
     * @brief Type of wrapped matrix.
     *
     * When moris::Base_Mat is inherited by moris::Mat, T is either:
     * -# moris::arma_Mat
     * -# moris::eigen_Mat
     *
     * when moris::Base_Mat is inherited by moris::Sp_Mat, T is either:
     * -# moris::arma_Sp_Mat
     * -# moris::eigen_Sp_Mat
     */
    T  mMat;

public:

    typedef moris::Tuple< moris::size_t , moris::size_t > size_tuple_t;
    typedef typename T::type_t type_t;

    /**
     * Base_Mat default constructor.
     */
    inline
    Base_Mat() = default;

    // -------------------------------------------------------------------------

    /**
     * Base_Mat copy constructor for copying moris::Mat.
     */
    inline
    Base_Mat(
            moris::Mat<type_t> const & aMat )
    : mMat( aMat.data() )
    {
    }

    // -------------------------------------------------------------------------

    /**
     * Base_Mat copy constructor for copying moris::Sp_Mat.
     */
    inline
    Base_Mat(
            moris::Sp_Mat<type_t> const & aMat )
    : mMat( aMat.data() )
    {
    }

    // -------------------------------------------------------------------------

    /**
     * Base_Mat initialization constructor.
     *
     * @param[in] aIElems Number of rows in matrix.
     * @param[in] aJElems Number of cols in matrix.
     */
    inline
    Base_Mat(
            moris::size_t const & aIElems,
            moris::size_t const & aJElems )
        : mMat( aIElems, aJElems )
    {
    }

    // -------------------------------------------------------------------------

    /**
     * Base_Mat copy constructor.
     */
    template< typename A >
    Base_Mat(
            A const & X )
            : mMat( X )
      {
      }

    // -------------------------------------------------------------------------

    /**
     * Sp_Mat default constructor.
     */
    inline
    Base_Mat(
            moris::Mat< moris::uint > const & aRowInd,
            moris::Mat< moris::uint > const & aColInd,
            moris::Mat< type_t >      const & aValues )
    : mMat( aRowInd, aColInd, aValues )
      {
      }

    // -------------------------------------------------------------------------

    /**
     * Sp_Mat default constructor.
     */
    inline
    Base_Mat(
            moris::Mat< moris::uint > const & aRowInd,
            moris::Mat< moris::uint > const & aColInd,
            moris::Mat< type_t >      const & aValues,
            moris::size_t             const & aIElems,
            moris::size_t             const & aJElems )
    : mMat( aRowInd, aColInd, aValues, aIElems, aJElems )
      {
      }

    // -------------------------------------------------------------------------

    /**
     * @brief Initialization constructor - Will not work (and doesn't make sense) for empty matrices
     *
     * @note a 2D initialization list is used to construct the matrices.
     * For example - { {1,2,3}, {4,5,6} } will give a 2 by 3 matrix.
     * list.size() corresponds to the number of rows. list.begin()->size()
     * corresponds to the number of colunms. We use 'range based for loops' to ensure
     * that the constructor works well for complex numbers too.
     */
    Base_Mat(
            std::initializer_list< std::initializer_list< type_t > > list)
        : mMat( list )
    {
    }

    // -------------------------------------------------------------------------

    /**
     * Base_Mat destructor.
     */
    inline
    ~Base_Mat() = default;

    // -------------------------------------------------------------------------

    /**
     * @brief constant accessor
     *
     * Returns a reference to the underlying vector.
     *
     * @return Reference to the underlying vector.
     */
    auto
    data() const
    -> decltype( mMat.data() )
    {
        return mMat.data();
    }

    // -------------------------------------------------------------------------

    /**
     * Returns a non-const reference to the underlying vector.
     *
     * @return Reference to the underlying vector.
     */
    auto
    data()
    -> decltype( mMat.data() )
    {
        return mMat.data();
    }

    // -------------------------------------------------------------------------

    /**
     * Assignment operator.
     *
     * @param[in] X Matrix or column vector or row vector.
     *
     * @return Assignment.
     */
    template< typename A >
    const moris::Base_Mat< T > &
    operator=(
            A const & X )
    {
        mMat.operator=( X );

        return *this;
    }

    // -------------------------------------------------------------------------

    /**
     * Assignment operator.
     *
     * @param[in] X Matrix or column vector or row vector.
     *
     * @return Assignment.
     */
    const moris::Base_Mat< T > &
    operator=(
            moris::Base_Mat< T > const & X )
    {
        mMat = X.data();

        return *this;
    }

    // -------------------------------------------------------------------------

    /**
     * Get column vector, similar to Matlab A(:,aK).
     *
     * @param[in] aK Column index.
     *
     * @return Column vector.
     *
     * Example:
     * @include LNA/src/cl_Mat/Mat_col.inc
     */
    auto
    col(
        moris::size_t const & aK )
    -> decltype( mMat.col( aK ) )
    {
        return mMat.col( aK );
    }

    // -------------------------------------------------------------------------

    /**
     * Get column vector, similar to Matlab A(:,aK).
     *
     * @param[in] aK Column index.
     *
     * @return Column vector.
     *
     * Example:
     * @include LNA/src/cl_Mat/Mat_col.inc
     */
    auto
    col(
        moris::size_t const & aK ) const
    -> decltype( mMat.col( aK ) )
    {
        return mMat.col( aK );
    }

    // -------------------------------------------------------------------------

    /**
     * Get column vectors, similar to Matlab A(:,aP:aQ).
     *
     * @param[in] aP column beginning index
     * @param[in] aQ column ending index
     *
     * @return Matrix of specified columns.
     *
     * Example:
     * @include LNA/src/cl_Mat/Mat_cols.inc
     */
    auto
    cols(
        moris::size_t const & aP,
        moris::size_t const & aQ )
    -> decltype( mMat.cols( aP, aQ ) )
    {
        return mMat.cols( aP, aQ );
    }

    // -------------------------------------------------------------------------

    /**
     * Get column vectors, similar to Matlab A(:,aP:aQ).
     *
     * @param[in] aP column beginning index
     * @param[in] aQ column ending index
     *
     * @return Matrix of specified columns.
     *
     * Example:
     * @include LNA/src/cl_Mat/Mat_cols.inc
     */
    auto
    cols(
        moris::size_t const & aP,
        moris::size_t const & aQ ) const
    -> decltype( mMat.cols( aP, aQ ) )
    {
        return mMat.cols( aP, aQ );
    }

    // -------------------------------------------------------------------------

    /**
     * Get the number of columns in a data set, similar to Matlab cols().
     *
     * @return Number of columns.
     */
    moris::size_t
    n_cols() const
    {
        return mMat.n_cols();
    }

    // -------------------------------------------------------------------------

    /**
     * Get row vector, similar to Matlab A(aK,:).
     *
     * @param[in] aK row index.
     *
     * @return row vector.
     *
     * Example:
     * @include LNA/src/cl_Mat/Mat_row.inc
     */
    auto
    row(
            moris::size_t const & aK )
    -> decltype( mMat.row( aK ) )
    {
        return mMat.row( aK );
    }

    // -------------------------------------------------------------------------

    /**
     * Get row vector, similar to Matlab A(aK,:).
     *
     * @param[in] aK row index.
     *
     * @return row vector.
     *
     * Example:
     * @include LNA/src/cl_Mat/Mat_row.inc
     */
    auto
    row(
            moris::size_t const & aK ) const
    -> decltype( mMat.row( aK ) )
    {
        return mMat.row( aK );
    }

    // -------------------------------------------------------------------------

    /**
     * Get row vectors, similar to Matlab A(aP:aQ,:).
     *
     * @param[in] aP row beginning index
     * @param[in] aQ row ending index
     *
     * @return Matrix of specified rows.
     *
     * Example:
     * @include LNA/src/cl_Mat/Mat_rows.inc
     */
    auto
    rows(
        moris::size_t const & aP,
        moris::size_t const & aQ )
    -> decltype( mMat.rows( aP, aQ ) )
    {
        return mMat.rows( aP, aQ );
    }

    // -------------------------------------------------------------------------

    /**
     * Get row vectors, similar to Matlab A(aP:aQ,:).
     *
     * @param[in] aP row beginning index
     * @param[in] aQ row ending index
     *
     * @return Matrix of specified rows.
     *
     * Example:
     * @include LNA/src/cl_Mat/Mat_rows.inc
     */
    auto
    rows(
        moris::size_t const & aP,
        moris::size_t const & aQ ) const
   -> decltype( mMat.rows( aP, aQ ) )
   {
        return mMat.rows( aP, aQ );
   }

    // -------------------------------------------------------------------------

    /**
     * Get the number of rows in a data set, similar to Matlab rows().
     *
     * @return Number of rows.
     */
    moris::size_t
    n_rows() const
    {
        return mMat.n_rows();
    }

    // -------------------------------------------------------------------------

    /**
     * Returns the number of elements in the %matrix.
     *
     * @return Number of elements in the %matrix.
     *
     */
    auto
    numel() const
    -> decltype( mMat.numel() )
    {
        return mMat.numel();
    }

    // -------------------------------------------------------------------------

    /**
     * Returns the number of rows and columns of a %matrix.
     * If a direction is specified (0 = rows, 1 = columns) the size of the
     * matrix in this direction is returned.
     *
     * @return Number of rows and columns of a matrix.
     */
    size_tuple_t
    size()
    {
        return size_tuple_t( mMat.n_rows() , mMat.n_cols() );
    }

    // -------------------------------------------------------------------------

    moris::size_t
    size(
            moris::size_t const & aDim ) const
    {
        return mMat.size( aDim );
    }

    // -------------------------------------------------------------------------

    moris::size_t
    length() const
    {
        moris::size_t n_rows = this->n_rows();
        moris::size_t n_cols = this->n_cols();

        if( n_rows == 0 || n_cols == 0)
        {
            return 0;
        }
        else
        {
            // Check if a mxn or empty matrix was provided
            if ( n_rows != 1 && n_cols != 1 )
            {
                MORIS_LOG_ERROR << "Tried to get length of a matrix. Check dimensions.";
            }
            return ( n_rows < n_cols ) ? n_cols : n_rows;
        }
    }

    // -------------------------------------------------------------------------

    /**
     * Print function
     *
     * @param[in] aVarName as the variable name to be printed
     *
     * @return Prints moris::Mat or moris::Sp_Mat in MATLAB style to the screen
     */
    void
    print( const std::string & aVarName = std::string() )
    {
        FILE * outFile = stdout;

        fprintf( outFile,"%s-------------------------------------------------\n\n","%" );

        if ( aVarName.empty() )
            fprintf( outFile, "morisMat = [ ... \n" );
        else
            fprintf( outFile, "%s = [ ... \n", aVarName.c_str() );

        for( moris::uint ir = 0; ir < mMat.n_rows(); ++ir )
        {
            for( moris::uint ic = 0; ic < mMat.n_cols(); ++ic )
            {
                // FIXME: need to type cast the output of mMat(ir,ic)
                fprintf( outFile, "%+.15e", (moris::real)mMat(ir,ic) );

                if(ic < mMat.n_cols() -1)
                    fprintf( outFile, ",  " );
                else
                    fprintf( outFile, ";" );
            }

            if (ir < mMat.n_rows() - 1 )
                fprintf( outFile, ".... \n" );
            else
                fprintf( outFile, "];\n" );
        }
        fprintf( outFile,"%s-------------------------------------------------\n\n", "%") ;
    }

    /**
     * Print function for constant matrices
     */
    void
    print( const std::string & aVarName = std::string() ) const
    {
        FILE * outFile = stdout;

        fprintf( outFile,"%s-------------------------------------------------\n\n","%" );

        if ( aVarName.empty() )
            fprintf( outFile, "morisMat = [ ... \n" );
        else
            fprintf( outFile, "%s = [ ... \n", aVarName.c_str() );

        for( moris::uint ir = 0; ir < mMat.n_rows(); ++ir )
        {
            for( moris::uint ic = 0; ic < mMat.n_cols(); ++ic )
            {
                // FIXME: need to type cast the output of mMat(ir,ic)
                fprintf( outFile, "%+.15e", (moris::real)mMat(ir,ic) );

                if(ic < mMat.n_cols() -1)
                    fprintf( outFile, ",  " );
                else
                    fprintf( outFile, ";" );
            }

            if (ir < mMat.n_rows() - 1 )
                fprintf( outFile, ".... \n" );
            else
                fprintf( outFile, "]\n" );
        }
        fprintf( outFile,"%s-------------------------------------------------\n\n", "%") ;
    }
};

#endif  /* MORIS_LINALG_CL_BASE_MAT_HPP_ */

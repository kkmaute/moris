#ifndef MORIS_LINALG_CL_BASE_EIGEN_MAT_HPP_
#define MORIS_LINALG_CL_BASE_EIGEN_MAT_HPP_

// C++ header files.
#include <utility>

// MORIS library header files.
#include "typedefs.hpp" // COR/src
#include "assert.hpp"

// Class forward declarations.
namespace moris {

    /**
     * The Eigen base matrix class
     */
    template< typename T >
    class Base_Eigen_Mat;
}

template< typename T >
class moris::Base_Eigen_Mat
{

protected:

    // in child classes, T is either Eigen::Matrix or Eigen::SparseMatrix
    T mMat;

public:

    typedef moris::Tuple< moris::size_t , moris::size_t > size_tuple_t;
    typedef typename T::Scalar type_t;

    // default constructor
    inline
    Base_Eigen_Mat() = default;

    // -------------------------------------------------------------------------

    // size constructor
    inline
    Base_Eigen_Mat(
            moris::size_t const & i_elems,
            moris::size_t const & j_elems )
        : mMat(i_elems, j_elems )
    {
    }

    // -------------------------------------------------------------------------

    // initializer list is allowed for eigen::SpMat. see moris::Mat for documentation
    Base_Eigen_Mat(
            std::initializer_list< std::initializer_list< type_t > > list)
        : mMat( list.size(), list.begin()->size() )
    {
        MORIS_ASSERT( list.size() > 0, "The initializer list is empty." );

      moris::size_t i = 0, j = 0;
      for (const auto row : list) // loop over number of rows
      {
          MORIS_ASSERT( row.size() == list.begin()->size(),
          "The number of elements in one of the rows does not equal the number of columns." );

          for (const auto col : row) // loop over every value in the row
          {
              mMat( i, j ) = col;
              ++j;
          }
          j = 0;
          ++i;
      }
    }

    // -------------------------------------------------------------------------

    // default destructor
    inline
    ~Base_Eigen_Mat() = default;

    // -------------------------------------------------------------------------

    // two index access (different for dense and sparse matrices, so it is virtual here)
    virtual
    typename T::Scalar &
    operator()(
            const moris::size_t & i_index,
            const moris::size_t & j_index ) = 0;

    // -------------------------------------------------------------------------

    // data access
    auto
    data()
    -> decltype( mMat ) &
    {
        return mMat;
    }

    // const data access
    auto
    data() const
    -> decltype( mMat ) const &
    {
        return mMat;
    }

    // -------------------------------------------------------------------------

    // const assignment operator
    template< typename A >
    const moris::Base_Eigen_Mat< T > &
    operator=(
            A const & X )
    {
        mMat = X;

        return *this;
    }

    // -------------------------------------------------------------------------

    // copy constructor
    const moris::Base_Eigen_Mat< T > &
    operator=(
            moris::Base_Eigen_Mat< T > const & X )
    {
        mMat = X.data();

        return *this;
    }

    // -------------------------------------------------------------------------

    // column access
    auto
    col(
        moris::size_t const & aK )
    -> decltype( mMat.col( aK ) )
    {
        return mMat.col( aK );
    }

    // const column access
    auto
    col(
        moris::size_t const & aK ) const
    -> decltype( mMat.col( aK ) )
    {
        return mMat.col( aK );
    }

    // -------------------------------------------------------------------------

    // columns access
    auto
    cols(
        moris::size_t const & aP,
        moris::size_t const & aQ )
    -> decltype( mMat.block( 0, aP, mMat.rows(), aQ-aP+1 ) )
    {
        return mMat.block( 0, aP, mMat.rows(), aQ-aP+1 );
    }

    // const columns access
    auto
    cols(
        moris::size_t const & aP,
        moris::size_t const & aQ ) const
    -> decltype( mMat.block( 0, aP, mMat.rows(), aQ-aP+1 ) )
    {
        return mMat.block( 0, aP, mMat.rows(), aQ-aP+1 );
    }

    // -------------------------------------------------------------------------

    // number of cols
    moris::size_t
    n_cols() const
    {
        return mMat.cols();
    }

    // -------------------------------------------------------------------------

    // row access
    auto
    row(
            moris::size_t const & aK )
    -> decltype( mMat.row( aK ) )
    {
        return mMat.row( aK );
    }

    // const row access
    auto
    row(
            moris::size_t const & aK ) const
    -> decltype( mMat.row( aK ) )
    {
        return mMat.row( aK );
    }

    // -------------------------------------------------------------------------

    // rows access
    auto
    rows(
        moris::size_t const & aP,
        moris::size_t const & aQ )
    -> decltype( mMat.block( aP, 0, aQ-aP+1, mMat.cols() ) )
    {
        return mMat.block( aP, 0, aQ-aP+1, mMat.cols() );
    }

    // const rows access
    auto
    rows(
        moris::size_t const & aP,
        moris::size_t const & aQ ) const
   -> decltype( mMat.block( aP, 0, aQ-aP+1, mMat.cols() ) )
   {
        return mMat.block( aP, 0, aQ-aP+1, mMat.cols() );
   }

    // -------------------------------------------------------------------------

    // number of rows
    moris::size_t
    n_rows() const
    {
        return mMat.rows();
    }

    // -------------------------------------------------------------------------

    // number of elements
    auto
    numel() const
    -> decltype( mMat.size() )
    {
        return mMat.size();
    }

    // -------------------------------------------------------------------------

    // size as a tuple (similar to matlab)
    size_tuple_t
    size()
    {
        return size_tuple_t( mMat.rows() , mMat.cols() );
    }

    // size in a particular direction
    moris::size_t
    size(
        moris::size_t const & aDim ) const
    {
        if (aDim == 0)
        {
            return mMat.rows();
        }
        else if(aDim == 1)
        {
            return mMat.cols();
        }
        else
        {
            throw std::runtime_error ("Matrix size was requested, "
                    "for a dimension that has not been set up");
        }
    }

    // -------------------------------------------------------------------------

    moris::size_t
    length() const
    {
        moris::size_t n_rows = mMat.rows();
        moris::size_t n_cols = mMat.cols();

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

    // fast set_size function
    auto
    set_size(
        moris::size_t const & aNumRows,
        moris::size_t const & aNumCols )
    -> decltype( mMat.resize( aNumRows, aNumCols ) )
    {
        return mMat.resize( aNumRows, aNumCols );
    }
};

#endif  /* MORIS_LINALG_CL_BASE_EIGEN_MAT_HPP_ */

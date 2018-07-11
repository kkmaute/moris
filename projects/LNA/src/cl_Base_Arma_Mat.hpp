#ifndef MORIS_LINALG_CL_BASE_ARMA_MAT_HPP_
#define MORIS_LINALG_CL_BASE_ARMA_MAT_HPP_

// C++ header files.
#include <utility>

// MORIS library header files.
#include "typedefs.hpp"
#include "assert.hpp"

// Class forward declarations.
namespace moris {

    /**
     * Armadillo base matrix class
     */
    template< typename T >
    class Base_Arma_Mat;
}

template< typename T >
class moris::Base_Arma_Mat
{

protected:

    // in child classes, T is either arma::Mat or arma::SpMat
    T mMat;

public:

    typedef moris::Tuple< moris::size_t , moris::size_t > size_tuple_t;
    typedef typename T::elem_type type_t;

    // default constructor
    inline
    Base_Arma_Mat() = default;

    // -------------------------------------------------------------------------

    // size constructor
    inline
    Base_Arma_Mat(
            moris::size_t const & i_elems,
            moris::size_t const & j_elems )
        : mMat(i_elems, j_elems )
    {
    }

    // -------------------------------------------------------------------------

    // initializer list is allowed for arma::SpMat. see moris::Mat for documentation
    Base_Arma_Mat(
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

    //default destructor
    inline
    ~Base_Arma_Mat() = default;

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
    const moris::Base_Arma_Mat< T > &
    operator=(
            A const & X )
    {
        mMat = X;

        return *this;
    }

    // -------------------------------------------------------------------------

    // copy constructor
    const moris::Base_Arma_Mat< T > &
    operator=(
            moris::Base_Arma_Mat< T > const & X )
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
    -> decltype( mMat.cols( aP, aQ ) )
    {
        return mMat.cols( aP, aQ );
    }

    // const columns access
    auto
    cols(
        moris::size_t const & aP,
        moris::size_t const & aQ ) const
    -> decltype( mMat.cols( aP, aQ ) )
    {
        return mMat.cols( aP, aQ );
    }

    // -------------------------------------------------------------------------

    // number of columns
    moris::size_t
    n_cols() const
    {
        return mMat.n_cols;
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
    -> decltype( mMat.rows( aP, aQ ) )
    {
        return mMat.rows( aP, aQ );
    }

    // const rows access
    auto
    rows(
        moris::size_t const & aP,
        moris::size_t const & aQ ) const
   -> decltype( mMat.rows( aP, aQ ) )
   {
        return mMat.rows( aP, aQ );
   }

    // -------------------------------------------------------------------------

    // number of rows
    moris::size_t
    n_rows() const
    {
        return mMat.n_rows;
    }

    // -------------------------------------------------------------------------

    // number of elements
    auto
    numel() const
    -> decltype( mMat.n_elem )
    {
        return mMat.n_elem;
    }

    // -------------------------------------------------------------------------

    // size as a tuple (similar to matlab)
    size_tuple_t
    size()
    {
        return size_tuple_t( mMat.n_rows , mMat.n_cols );
    }

    // size in a particular direction
    moris::size_t
    size(
        moris::size_t const & aDim ) const
    {
        if (aDim == 0)
        {
            return mMat.n_rows;
        }
        else if(aDim == 1)
        {
            return mMat.n_cols;
        }
        else
        {
            throw std::runtime_error ("Matrix size was requested,"
                      "for a dimension that has not been set up");
        }
    }

    // -------------------------------------------------------------------------

    moris::size_t
    length() const
    {
        moris::size_t n_rows = mMat.n_rows;
        moris::size_t n_cols = mMat.n_cols;

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
    -> decltype( mMat.set_size( aNumRows, aNumCols ) )
    {
        return mMat.set_size( aNumRows, aNumCols );
    }
};

#endif  /* MORIS_LINALG_CL_BASE_ARMA_MAT_HPP_ */

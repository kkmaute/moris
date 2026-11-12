/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_DYNAMIC_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_DYNAMIC_HPP_

#include "cl_Matrix.hpp"
#include <Eigen/Dense>
#include "fn_Eigen_Aliasing.hpp"
#include <algorithm>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace moris
{
    // A stable native Map exposes both owned and borrowed storage to Eigen's
    // expression engine. Fixed and dynamic matrices share the same public API.
    template< typename Type, int Rows, int Cols, int Options, int MaxRows, int MaxCols >
    class Matrix< Eigen::Matrix< Type, Rows, Cols, Options, MaxRows, MaxCols > >
    {
        using Native = Eigen::Matrix< Type, Rows, Cols, Options, MaxRows, MaxCols >;
        using View   = Eigen::Map< Native >;
        Native mOwned;
        View   mMatrix;
        Type*  mExternal = nullptr;
        bool   mStrict   = false;

        void remap( Type* aData, size_t aRows, size_t aCols )
        {
            // Reconstruct the Map to rebind storage; assignment would copy coefficients.
            mMatrix.~View();
            new ( &mMatrix ) View( aData, aRows, aCols );
        }

        void dimensions( size_t aRows, size_t aCols ) const
        {
            MORIS_ERROR( ( Rows == Eigen::Dynamic || aRows == static_cast< size_t >( Rows ) ) && ( Cols == Eigen::Dynamic || aCols == static_cast< size_t >( Cols ) ),
                    "Matrix - incompatible fixed matrix dimensions." );
        }

        void fill_with_NANs()
        {
#ifdef MATRIX_FILL
            if constexpr ( std::numeric_limits< Type >::has_quiet_NaN )
            {
                mMatrix.fill( std::numeric_limits< Type >::quiet_NaN() );
            }
            else
            {
                mMatrix.fill( std::numeric_limits< Type >::max() );
            }
#endif
        }

        void resize_storage( size_t aRows, size_t aCols )
        {
            dimensions( aRows, aCols );
            if ( n_rows() == aRows && n_cols() == aCols )
            {
                return;
            }
            if ( mExternal && aRows * aCols == numel() )
            {
                remap( mExternal, aRows, aCols );
                return;
            }
            MORIS_ERROR( !mExternal || !mStrict, "Matrix - strict borrowed storage cannot change its number of elements." );
            mOwned.resize( aRows, aCols );
            mExternal = nullptr;
            mStrict   = false;
            remap( mOwned.data(), aRows, aCols );
        }

        template< bool Maximum >
        Type complex_extreme( uint& aRow, uint& aCol ) const
        {
            MORIS_ASSERT( numel() > 0, "Matrix - extrema require a nonempty matrix." );
            size_t tIndex = 0;
            for ( size_t i = 1; i < numel(); ++i )
            {
                if ( Maximum ? std::abs( data()[ i ] ) > std::abs( data()[ tIndex ] )
                             : std::abs( data()[ i ] ) < std::abs( data()[ tIndex ] ) )
                {
                    tIndex = i;
                }
            }
            aRow = tIndex % n_rows();
            aCol = tIndex / n_rows();
            return data()[ tIndex ];
        }

        template< typename E >
        void assign( const E& aExpression )
        {
            if ( linalg_internal::eigen_assignment_hazard( aExpression, data(), n_rows(), n_cols() )
                    || ( ( aExpression.rows() != mMatrix.rows() || aExpression.cols() != mMatrix.cols() )
                            && linalg_internal::eigen_references( aExpression, data(), numel() ) ) )
            {
                // Evaluate before resize or reordered writes invalidate the RHS.
                typename E::PlainObject tCopy( aExpression );
                resize_storage( tCopy.rows(), tCopy.cols() );
                mMatrix = tCopy;
            }
            else
            {
                resize_storage( aExpression.rows(), aExpression.cols() );
                mMatrix = aExpression;
            }
        }

      public:
        using Data_Type = Type;
        Matrix()
                : mMatrix( mOwned.data(), mOwned.rows(), mOwned.cols() )
        {
        }
        Matrix( size_t aRows, size_t aCols )
                : Matrix()
        {
            resize_storage( aRows, aCols );
            fill_with_NANs();
        }
        explicit Matrix( size_t aSize )
                : Matrix( aSize, 1 )
        {
        }
        Matrix( size_t aRows, size_t aCols, const Type& aValue )
                : Matrix( aRows, aCols )
        {
            fill( aValue );
        }
        Matrix( const Matrix& aOther )
                : Matrix()
        {
            assign( aOther.matrix_data() );
        }
        template< typename Other >
        Matrix( const Matrix< Other >& aOther )
                : Matrix()
        {
            assign( aOther.matrix_data() );
        }
        template< typename E >
        Matrix( const Eigen::EigenBase< E >& aExpression )
                : Matrix()
        {
            assign( aExpression.derived() );
        }
        Matrix( Matrix&& aOther ) noexcept( std::is_nothrow_move_constructible< Native >::value )
                : mOwned( std::move( aOther.mOwned ) )
                , mMatrix( aOther.mExternal ? aOther.mExternal : mOwned.data(), aOther.n_rows(), aOther.n_cols() )
                , mExternal( aOther.mExternal )
                , mStrict( aOther.mStrict )
        {
            aOther.mExternal = nullptr;
            aOther.mStrict   = false;
            aOther.remap( aOther.mOwned.data(), aOther.mOwned.rows(), aOther.mOwned.cols() );
        }
        Matrix( const Type* const & aData, size_t aRows, size_t aCols )
                : Matrix()
        {
            resize_storage( aRows, aCols );
            mMatrix = Eigen::Map< const Native >( aData, aRows, aCols );
        }
        Matrix( Type* aData, size_t aRows, size_t aCols, bool aCopy, bool aStrict = false )
                : Matrix()
        {
            dimensions( aRows, aCols );
            if ( aCopy )
            {
                resize_storage( aRows, aCols );
                mMatrix = Eigen::Map< const Native >( aData, aRows, aCols );
            }
            else
            {
                mExternal = aData;
                mStrict   = aStrict;
                remap( aData, aRows, aCols );
            }
        }
        Matrix( std::initializer_list< std::initializer_list< Type > > aValues )
                : Matrix()
        {
            const size_t tRows = aValues.size();
            const size_t tCols = tRows ? aValues.begin()->size() : 0;
            resize_storage( tRows, tCols );
            size_t i = 0;
            for ( const auto& tRow : aValues )
            {
                MORIS_ERROR( tRow.size() == tCols, "Matrix - inconsistent initializer list rows." );
                size_t j = 0;
                for ( const auto& tValue : tRow )
                {
                    mMatrix( i, j++ ) = tValue;
                }
                ++i;
            }
        }

        Matrix& operator=( const Matrix& aOther )
        {
            if ( this != &aOther )
            {
                assign( aOther.matrix_data() );
            }
            return *this;
        }
        Matrix& operator=( Matrix&& aOther )
        {
            if ( this != &aOther )
            {
                // Retain borrowed destinations and owning storage referenced
                // by the source view; releasing that allocation would dangle it.
                if ( mExternal || linalg_internal::eigen_references( aOther.matrix_data(), data(), numel() ) )
                {
                    assign( aOther.matrix_data() );
                }
                else
                {
                    const size_t tRows = aOther.n_rows(), tCols = aOther.n_cols();
                    mOwned    = std::move( aOther.mOwned );
                    mExternal = aOther.mExternal;
                    mStrict   = aOther.mStrict;
                    remap( mExternal ? mExternal : mOwned.data(), tRows, tCols );
                    aOther.mExternal = nullptr;
                    aOther.mStrict   = false;
                    aOther.remap( aOther.mOwned.data(), aOther.mOwned.rows(), aOther.mOwned.cols() );
                }
            }
            return *this;
        }
        template< typename Other >
        Matrix& operator=( const Matrix< Other >& aOther )
        {
            assign( aOther.matrix_data() );
            return *this;
        }
        template< typename E >
        Matrix& operator=( const Eigen::EigenBase< E >& aExpression )
        {
            assign( aExpression.derived() );
            return *this;
        }
        Matrix& operator=( const Type& aValue )
        {
            resize_storage( 1, 1 );
            mMatrix( 0, 0 ) = aValue;
            return *this;
        }
        Matrix& operator+=( const Matrix& aOther )
        {
            return *this += aOther.matrix_data();
        }
        Matrix& operator-=( const Matrix& aOther )
        {
            return *this -= aOther.matrix_data();
        }
        template< typename E >
        Matrix& operator+=( const Eigen::MatrixBase< E >& aExpression )
        {
            if ( linalg_internal::eigen_assignment_hazard( aExpression.derived(), data(), n_rows(), n_cols() ) )
            {
                mMatrix += aExpression.eval();
            }
            else
            {
                mMatrix += aExpression;
            }
            return *this;
        }
        template< typename E >
        Matrix& operator-=( const Eigen::MatrixBase< E >& aExpression )
        {
            if ( linalg_internal::eigen_assignment_hazard( aExpression.derived(), data(), n_rows(), n_cols() ) )
            {
                mMatrix -= aExpression.eval();
            }
            else
            {
                mMatrix -= aExpression;
            }
            return *this;
        }
        Matrix& operator+=( const Type& aValue )
        {
            mMatrix.array() += aValue;
            return *this;
        }
        Matrix& operator-=( const Type& aValue )
        {
            mMatrix.array() -= aValue;
            return *this;
        }
        size_t n_rows() const
        {
            return mMatrix.rows();
        }
        size_t n_cols() const
        {
            return mMatrix.cols();
        }
        size_t numel() const
        {
            return mMatrix.size();
        }
        size_t capacity() const
        {
            return numel() * sizeof( Type );
        }
        size_t size( size_t aDim ) const
        {
            MORIS_ASSERT( aDim < 2, "Matrix - invalid dimension." );
            return aDim == 0 ? n_rows() : n_cols();
        }
        size_t length() const
        {
            MORIS_ASSERT( !numel() || n_rows() == 1 || n_cols() == 1, "Matrix - length requires a vector." );
            return numel() ? std::max( n_rows(), n_cols() ) : 0;
        }
        Matrix copy() const
        {
            return Matrix( *this );
        }
        void set_size( size_t aRows, size_t aCols )
        {
            resize_storage( aRows, aCols );
            fill_with_NANs();
        }
        void set_size( size_t aRows, size_t aCols, const Type& aValue )
        {
            resize_storage( aRows, aCols );
            fill( aValue );
        }
        void resize( size_t aRows, size_t aCols )
        {
            if ( aRows == n_rows() && aCols == n_cols() )
            {
                return;
            }
            Native tOld( mMatrix );
            resize_storage( aRows, aCols );
            mMatrix.setZero();
            const size_t tRows                  = std::min( aRows, static_cast< size_t >( tOld.rows() ) );
            const size_t tCols                  = std::min( aCols, static_cast< size_t >( tOld.cols() ) );
            mMatrix.block( 0, 0, tRows, tCols ) = tOld.block( 0, 0, tRows, tCols );
        }
        void reshape( size_t aRows, size_t aCols )
        {
            Native tOld( mMatrix );
            resize_storage( aRows, aCols );
            mMatrix.setZero();
            std::copy_n( tOld.data(), std::min( numel(), static_cast< size_t >( tOld.size() ) ), data() );
        }
        void inplace_trans()
        {
            Eigen::Matrix< Type, Eigen::Dynamic, Eigen::Dynamic > tCopy( mMatrix.adjoint() );
            assign( tCopy );
        }
        void fill( const Type& aValue )
        {
            mMatrix.fill( aValue );
        }
        View& matrix_data()
        {
            return mMatrix;
        }
        const View& matrix_data() const
        {
            return mMatrix;
        }
        Type* data()
        {
            return mMatrix.data();
        }
        const Type* data() const
        {
            return mMatrix.data();
        }
        Type* colptr( size_t aCol )
        {
            MORIS_ASSERT( aCol < n_cols(), "Matrix - invalid column." );
            return data() + aCol * n_rows();
        }
        const Type* colptr( size_t aCol ) const
        {
            MORIS_ASSERT( aCol < n_cols(), "Matrix - invalid column." );
            return data() + aCol * n_rows();
        }
        Type* begin()
        {
            return data();
        }
        const Type* begin() const
        {
            return data();
        }
        const Type* cbegin() const
        {
            return data();
        }
        Type* end()
        {
            return numel() ? data() + numel() : data();
        }
        const Type* end() const
        {
            return numel() ? data() + numel() : data();
        }
        const Type* cend() const
        {
            return end();
        }
        Type& operator()( size_t i, size_t j )
        {
            MORIS_ASSERT( i < n_rows() && j < n_cols(), "Matrix - index out of bounds." );
            return mMatrix( i, j );
        }
        const Type& operator()( size_t i, size_t j ) const
        {
            MORIS_ASSERT( i < n_rows() && j < n_cols(), "Matrix - index out of bounds." );
            return mMatrix( i, j );
        }
        Type& operator()( size_t i )
        {
            MORIS_ASSERT( ( n_rows() == 1 || n_cols() == 1 ) && i < numel(), "Matrix - invalid vector index." );
            return data()[ i ];
        }
        const Type& operator()( size_t i ) const
        {
            MORIS_ASSERT( ( n_rows() == 1 || n_cols() == 1 ) && i < numel(), "Matrix - invalid vector index." );
            return data()[ i ];
        }
        auto operator()( const std::pair< size_t, size_t >& i, const std::pair< size_t, size_t >& j )
        {
            return mMatrix.block( i.first, j.first, i.second - i.first + 1, j.second - j.first + 1 );
        }
        auto operator()( const std::pair< size_t, size_t >& i, const std::pair< size_t, size_t >& j ) const
        {
            return mMatrix.block( i.first, j.first, i.second - i.first + 1, j.second - j.first + 1 );
        }
        auto operator()( const std::pair< size_t, size_t >& i, size_t j )
        {
            return mMatrix.col( j ).segment( i.first, i.second - i.first + 1 );
        }
        auto operator()( const std::pair< size_t, size_t >& i, size_t j ) const
        {
            return mMatrix.col( j ).segment( i.first, i.second - i.first + 1 );
        }
        auto operator()( const std::pair< size_t, size_t >& i )
        {
            return ( *this )( i, 0 );
        }
        auto operator()( const std::pair< size_t, size_t >& i ) const
        {
            return ( *this )( i, 0 );
        }
        auto get_row( size_t i )
        {
            MORIS_ASSERT( i < n_rows(), "Matrix - invalid row." );
            return mMatrix.row( i );
        }
        auto get_row( size_t i ) const
        {
            MORIS_ASSERT( i < n_rows(), "Matrix - invalid row." );
            return mMatrix.row( i );
        }
        auto get_column( size_t j )
        {
            MORIS_ASSERT( j < n_cols(), "Matrix - invalid column." );
            return mMatrix.col( j );
        }
        auto get_column( size_t j ) const
        {
            MORIS_ASSERT( j < n_cols(), "Matrix - invalid column." );
            return mMatrix.col( j );
        }
        template< typename Other >
        void get_row( size_t i, Matrix< Other >& aRow ) const
        {
            aRow = get_row( i );
        }
        template< typename Other >
        void get_column( size_t j, Matrix< Other >& aColumn ) const
        {
            aColumn = get_column( j );
        }
        template< typename Other >
        void set_row( size_t i, const Matrix< Other >& aRow )
        {
            set_row( i, aRow.matrix_data() );
        }
        template< typename E >
        void set_row( size_t i, const Eigen::MatrixBase< E >& aRow )
        {
            MORIS_ASSERT( i < n_rows() && aRow.size() == static_cast< Eigen::Index >( n_cols() ) && ( aRow.rows() == 1 || aRow.cols() == 1 ), "Matrix - invalid row." );

            // Preserve native expressions unless overlapping storage requires a snapshot.
            if ( linalg_internal::eigen_references( aRow.derived(), data(), numel() ) )
            {
                const typename E::PlainObject tRow = aRow.eval();
                set_row( i, tRow );
            }
            else if ( aRow.rows() == 1 )
            {
                mMatrix.row( i ) = aRow.block( 0, 0, 1, aRow.cols() );
            }
            else
            {
                mMatrix.row( i ) = aRow.block( 0, 0, aRow.rows(), 1 ).transpose();
            }
        }
        template< typename Other >
        void set_column( size_t j, const Matrix< Other >& aColumn )
        {
            set_column( j, aColumn.matrix_data() );
        }
        template< typename E >
        void set_column( size_t j, const Eigen::MatrixBase< E >& aColumn )
        {
            MORIS_ASSERT( j < n_cols() && aColumn.cols() == 1 && aColumn.rows() == static_cast< Eigen::Index >( n_rows() ), "Matrix - invalid column." );

            if ( linalg_internal::eigen_references( aColumn.derived(), data(), numel() ) )
            {
                mMatrix.col( j ) = aColumn.eval();
            }
            else
            {
                mMatrix.col( j ) = aColumn;
            }
        }
        Type min() const
        {
            uint tRow = 0, tCol = 0;
            return min( tRow, tCol );
        }
        Type max() const
        {
            uint tRow = 0, tCol = 0;
            return max( tRow, tCol );
        }
        Type min( uint& i, uint& j ) const
        {
            if constexpr ( Eigen::NumTraits< Type >::IsComplex )
            {
                return complex_extreme< false >( i, j );
            }
            else
            {
                return mMatrix.minCoeff( &i, &j );
            }
        }
        Type max( uint& i, uint& j ) const
        {
            if constexpr ( Eigen::NumTraits< Type >::IsComplex )
            {
                return complex_extreme< true >( i, j );
            }
            else
            {
                return mMatrix.maxCoeff( &i, &j );
            }
        }
    };
}    // namespace moris
#endif

/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Lumped_Matrix.hpp
 *
 */

#ifndef MORIS_CL_FEM_LUMPED_MATRIX_HPP
#define MORIS_CL_FEM_LUMPED_MATRIX_HPP

#include <utility>
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "assert.hpp"

namespace moris::fem
{
    class Scaled_Lumped_Matrix;

    /**
     * Spatially diagonal blocks with full coupling between time bases.
     * Coefficients are ordered by field, then time basis, then spatial basis.
     * No dense conversion is provided: products and assembly visit stored entries.
     */
    class Lumped_Matrix
    {
        uint             mNumSpaceBases = 0;
        uint             mNumTimeBases  = 0;
        uint             mNumFields     = 0;
        Matrix< DDRMat > mDiagonals;

        template< typename Target >
        void add_entries( Target& aTarget, real aScale ) const
        {
            for ( uint iField = 0; iField < mNumFields; ++iField )
            {
                const uint tOffset = iField * mNumTimeBases * mNumSpaceBases;

                for ( uint iTime = 0; iTime < mNumTimeBases; ++iTime )
                {
                    for ( uint jTime = 0; jTime < mNumTimeBases; ++jTime )
                    {
                        for ( uint iSpace = 0; iSpace < mNumSpaceBases; ++iSpace )
                        {
                            aTarget( tOffset + iTime * mNumSpaceBases + iSpace,
                                    tOffset + jTime * mNumSpaceBases + iSpace ) +=
                                    aScale * diagonal( iField, iTime, jTime, iSpace );
                        }
                    }
                }
            }
        }

      public:
        // These are lightweight views in both supported dense backends.
        using Matrix_Block = decltype( std::declval< Matrix< DDRMat >& >()(
                std::pair< moris::size_t, moris::size_t >{},
                std::pair< moris::size_t, moris::size_t >{} ) );

        //------------------------------------------------------------------------------

        void set_size( uint aNumSpaceBases, uint aNumTimeBases, uint aNumFields )
        {
            MORIS_ERROR( aNumSpaceBases > 0 && aNumTimeBases > 0 && aNumFields > 0,
                    "Lumped_Matrix::set_size - Basis and field counts must be positive." );

            mNumSpaceBases = aNumSpaceBases;
            mNumTimeBases  = aNumTimeBases;
            mNumFields     = aNumFields;
            mDiagonals.set_size( aNumSpaceBases * aNumFields, aNumTimeBases * aNumTimeBases, 0.0 );
        }

        uint n_rows() const { return mNumSpaceBases * mNumTimeBases * mNumFields; }
        uint n_cols() const { return n_rows(); }
        uint num_stored_entries() const { return mDiagonals.numel(); }

        //------------------------------------------------------------------------------

        real& diagonal( uint aField, uint aRowTime, uint aColumnTime, uint aSpace )
        {
            MORIS_ASSERT( aField < mNumFields && aRowTime < mNumTimeBases
                                  && aColumnTime < mNumTimeBases && aSpace < mNumSpaceBases,
                    "Lumped_Matrix::diagonal - Index out of bounds." );

            return mDiagonals( aField * mNumSpaceBases + aSpace, aRowTime * mNumTimeBases + aColumnTime );
        }

        //------------------------------------------------------------------------------

        real diagonal( uint aField, uint aRowTime, uint aColumnTime, uint aSpace ) const
        {
            MORIS_ASSERT( aField < mNumFields && aRowTime < mNumTimeBases
                                  && aColumnTime < mNumTimeBases && aSpace < mNumSpaceBases,
                    "Lumped_Matrix::diagonal - Index out of bounds." );

            return mDiagonals( aField * mNumSpaceBases + aSpace, aRowTime * mNumTimeBases + aColumnTime );
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > multiply( const Matrix< DDRMat >& aValues, real aScale = 1.0 ) const
        {
            MORIS_ERROR( n_rows() > 0 && aValues.n_rows() == n_cols(),
                    "Lumped_Matrix::multiply - Incompatible dimensions." );

            Matrix< DDRMat > tResult( n_rows(), aValues.n_cols(), 0.0 );

            for ( uint iField = 0; iField < mNumFields; ++iField )
            {
                const uint tOffset = iField * mNumTimeBases * mNumSpaceBases;

                for ( uint iTime = 0; iTime < mNumTimeBases; ++iTime )
                {
                    for ( uint jTime = 0; jTime < mNumTimeBases; ++jTime )
                    {
                        for ( uint iSpace = 0; iSpace < mNumSpaceBases; ++iSpace )
                        {
                            const real tValue = aScale * diagonal( iField, iTime, jTime, iSpace );

                            for ( uint iVector = 0; iVector < aValues.n_cols(); ++iVector )
                            {
                                tResult( tOffset + iTime * mNumSpaceBases + iSpace, iVector ) +=
                                        tValue * aValues( tOffset + jTime * mNumSpaceBases + iSpace, iVector );
                            }
                        }
                    }
                }
            }

            return tResult;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > operator*( const Matrix< DDRMat >& aValues ) const
        {
            return multiply( aValues );
        }

        //------------------------------------------------------------------------------

        void add_to( Matrix< DDRMat >& aTarget, real aScale = 1.0 ) const
        {
            MORIS_ERROR( n_rows() > 0 && aTarget.n_rows() == n_rows() && aTarget.n_cols() == n_cols(),
                    "Lumped_Matrix::add_to - Incompatible matrix dimensions." );

            add_entries( aTarget, aScale );
        }

        //------------------------------------------------------------------------------

        void add_to( Matrix_Block aTarget, real aScale = 1.0 ) const
        {
#ifdef MORIS_USE_ARMA
            const bool tMatches = aTarget.n_rows == n_rows() && aTarget.n_cols == n_cols();
#else
            const bool tMatches = aTarget.rows() == n_rows() && aTarget.cols() == n_cols();
#endif
            MORIS_ERROR( n_rows() > 0 && tMatches,
                    "Lumped_Matrix::add_to - Incompatible block dimensions." );

            add_entries( aTarget, aScale );
        }
    };

    // Like backend matrix expressions, scaled operators borrow their operand.
    // Keep the source lumped matrix alive until the expression is evaluated.
    class Scaled_Lumped_Matrix
    {
      public:
        const Lumped_Matrix& mMatrix;
        real                 mScale;

        //------------------------------------------------------------------------------

        Matrix< DDRMat > operator*( const Matrix< DDRMat >& aValues ) const
        {
            return mMatrix.multiply( aValues, mScale );
        }
    };

    //------------------------------------------------------------------------------

    inline Scaled_Lumped_Matrix operator*( real aScale, const Lumped_Matrix& aMatrix )
    {
        return { aMatrix, aScale };
    }

    //------------------------------------------------------------------------------

    inline Scaled_Lumped_Matrix operator*( const Lumped_Matrix& aMatrix, real aScale )
    {
        return { aMatrix, aScale };
    }

    //------------------------------------------------------------------------------

    inline Scaled_Lumped_Matrix operator*( real aScale, const Scaled_Lumped_Matrix& aMatrix )
    {
        return { aMatrix.mMatrix, aScale * aMatrix.mScale };
    }

    //------------------------------------------------------------------------------

    inline Scaled_Lumped_Matrix operator*( const Scaled_Lumped_Matrix& aMatrix, real aScale )
    {
        return { aMatrix.mMatrix, aMatrix.mScale * aScale };
    }

    // Exact view overloads avoid extending generic backend expression machinery.
    //------------------------------------------------------------------------------

    inline Lumped_Matrix::Matrix_Block operator+=( Lumped_Matrix::Matrix_Block aTarget, const Scaled_Lumped_Matrix& aMatrix )
    {
        aMatrix.mMatrix.add_to( aTarget, aMatrix.mScale );

        return aTarget;
    }

    //------------------------------------------------------------------------------

    inline Lumped_Matrix::Matrix_Block operator-=( Lumped_Matrix::Matrix_Block aTarget, const Scaled_Lumped_Matrix& aMatrix )
    {
        aMatrix.mMatrix.add_to( aTarget, -aMatrix.mScale );

        return aTarget;
    }

    //------------------------------------------------------------------------------

    inline Lumped_Matrix::Matrix_Block operator+=( Lumped_Matrix::Matrix_Block aTarget, const Lumped_Matrix& aMatrix )
    {
        aMatrix.add_to( aTarget );

        return aTarget;
    }

    //------------------------------------------------------------------------------

    inline Lumped_Matrix::Matrix_Block operator-=( Lumped_Matrix::Matrix_Block aTarget, const Lumped_Matrix& aMatrix )
    {
        aMatrix.add_to( aTarget, -1.0 );

        return aTarget;
    }

    //------------------------------------------------------------------------------

    inline Matrix< DDRMat >& operator+=( Matrix< DDRMat >& aTarget, const Scaled_Lumped_Matrix& aMatrix )
    {
        aMatrix.mMatrix.add_to( aTarget, aMatrix.mScale );

        return aTarget;
    }

    //------------------------------------------------------------------------------

    inline Matrix< DDRMat >& operator-=( Matrix< DDRMat >& aTarget, const Scaled_Lumped_Matrix& aMatrix )
    {
        aMatrix.mMatrix.add_to( aTarget, -aMatrix.mScale );

        return aTarget;
    }

    //------------------------------------------------------------------------------

    inline Matrix< DDRMat >& operator+=( Matrix< DDRMat >& aTarget, const Lumped_Matrix& aMatrix )
    {
        aMatrix.add_to( aTarget );

        return aTarget;
    }

    //------------------------------------------------------------------------------

    inline Matrix< DDRMat >& operator-=( Matrix< DDRMat >& aTarget, const Lumped_Matrix& aMatrix )
    {
        aMatrix.add_to( aTarget, -1.0 );

        return aTarget;
    }

}    // namespace moris::fem
#endif

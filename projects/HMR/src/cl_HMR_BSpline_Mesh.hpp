/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_BSpline_Mesh.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_BSPLINE_MESH_HPP_
#define SRC_HMR_CL_HMR_BSPLINE_MESH_HPP_

#include "cl_HMR_Background_Element_Base.hpp" //HMR/src
#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_BSpline_Element.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "HMR_Globals.hpp" //HMR/src
#include "typedefs.hpp" //COR/src
#include "cl_Stopwatch.hpp" //CHR/src

namespace moris::hmr
{
    // ----------------------------------------------------------------------------

    /**
     * B-spline element class
     *
     * @tparam P Polynomial degree in x-direction
     * @tparam Q Polynomial degree in y-direction
     * @tparam R Polynomial degree in z-direction
     */
    template< uint P, uint Q, uint R >
    class BSpline_Mesh : public BSpline_Mesh_Base
    {
        //! Number of dimensions
        static constexpr uint N = ( P > 0 ) + ( Q > 0 ) + ( R > 0 );

        //! Container of degrees, used to avoid repeated conditionals
        static constexpr uint PQR[ 3 ] = { P, Q, R };

        //! Lookup table containing offset for node IDs
        luint mBasisLevelOffset[ gMaxNumberOfLevels ];

        //! Lookup table containing number of elements per dimension for each level
        luint mNumberOfElementsPerDimensionIncludingAura[ gMaxNumberOfLevels ][ N ];

        //! Lookup table for node IDs
        luint mMySubdomainOffset[ gMaxNumberOfLevels ][ N ];

        //! Lookup table containing number of basis per dimension for each level
        luint mNumberOfBasisPerDimensionIncludingPadding[ gMaxNumberOfLevels ][ N ];

    public:

        // ----------------------------------------------------------------------------

        /**
         * Constructor for Lagrange Mesh
         *
         * @param[in] aParameters       ref to container of user defined settings
         * @param[in] aBackgroundMesh pointer to background mesh
         *
         */
        BSpline_Mesh(
                const Parameters      * aParameters,
                Background_Mesh_Base  * aBackgroundMesh,
                uint aActivationPattern )
        : BSpline_Mesh_Base(
                aParameters,
                aBackgroundMesh,
                P,
                aActivationPattern )
        {
            // ask background mesh for number of elements per ijk-direction
            this->get_number_of_elements_per_dimension();

            // calculate lookup table mBasisLevelOffset
            this->calculate_lookup_tables();

            // calculate subdomain offset
            this->calculate_subdomain_offset();

            // calculate any value that can change after refinement
            this->update_mesh();
        }

        // ----------------------------------------------------------------------------

        /**
         * Default destructor.
         */
        ~BSpline_Mesh() override
        {
            mActiveBasisOnProc.clear();
            this->delete_pointers();
        }

    protected:

        // ----------------------------------------------------------------------------

        /**
         * creates a basis depending on polynomial order and dimension
         */
        Basis * create_basis(
                const luint * aIJK,
                uint          aLevel,
                uint          aOwner ) override;

        // ----------------------------------------------------------------------------

        /**
         * calculates domain wide unique basisID (1D case)
         * Useful for debugging.
         *
         * @param[in]  aLevel    level of basis
         * @param[in]  aI        proc local i-position of basis
         * @return uint          domain wide unique ID
         */
        luint calculate_basis_id(
                uint  aLevel,
                luint aI ) override
        {
            if( aLevel < gMaxNumberOfLevels && N == 1 )
            {
                return  aI + mMySubdomainOffset[ aLevel ][ 0 ];
            }
            else
            {
                return gNoEntityID;
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * calculates domain wide unique basisID (2D case)
         * Useful for debugging.
         *
         * @param[in]  aLevel    level of basis
         * @param[in]  aI        proc local i-position of basis
         * @param[in]  aJ        proc local j-position of basis
         * @return uint          domain wide unique ID
         */
        luint calculate_basis_id(
                uint  aLevel,
                luint aI,
                luint aJ ) override
        {
            if( aLevel < gMaxNumberOfLevels && N == 2 )
            {
                return  aI + mMySubdomainOffset[ aLevel ][ 0 ] +
                        ( aJ + mMySubdomainOffset[ aLevel ][ 1 ] ) *
                        (  mNumberOfBasisPerDimensionIncludingPadding[ aLevel ][ 0 ] ) +
                        mBasisLevelOffset[ aLevel ];
            }
            else
            {
                return gNoEntityID;
            }
        }

        // ----------------------------------------------------------------------------
        /**
         * calculates domain wide unique basisID (3D case)
         * Useful for debugging.
         *
         * @param[in]  aLevel    level of basis
         * @param[in]  aI        proc local i-position of basis
         * @param[in]  aJ        proc local j-position of basis
         * @param[in]  aK        proc local k-position of basis
         * @return uint          domain wide unique ID
         */
        luint calculate_basis_id(
                uint  aLevel,
                luint aI,
                luint aJ,
                luint aK ) override
        {
            if( aLevel < gMaxNumberOfLevels && N == 3 )
            {
                return aI + mMySubdomainOffset[ aLevel ][ 0 ] +
                        ( mNumberOfBasisPerDimensionIncludingPadding[ aLevel ][ 0 ] ) *
                        ( ( aJ + mMySubdomainOffset[ aLevel ][ 1 ] ) +
                                ( aK + mMySubdomainOffset[ aLevel ][ 2 ] ) *
                                ( mNumberOfBasisPerDimensionIncludingPadding[ aLevel ][ 1 ] ) ) +
                                mBasisLevelOffset[ aLevel ];
            }
            else
            {
                return gNoEntityID;
            }
        }

        // ----------------------------------------------------------------------------
    private:
        // ----------------------------------------------------------------------------

        /**
         *  Private function, creates the mNodeLevelOffset lookup table.
         *
         *  @return void
         */
        void calculate_lookup_tables()
        {
            // Number of basis on level 0
            for ( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                mNumberOfBasisPerDimensionIncludingPadding[ 0 ][ iDimension ] = mNumberOfElementsPerDimensionIncludingAura[ 0 ][ iDimension ] + PQR[ iDimension ];
            }

            // Basis level offset on level 0
            mBasisLevelOffset[ 0 ] = 0;

            for( uint iLevel = 1; iLevel < gMaxNumberOfLevels; iLevel++ )
            {
                // Start counting number of bases
                luint tNumberOfBasis = 1;

                for ( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    // Calculate number of basis on higher levels
                    mNumberOfBasisPerDimensionIncludingPadding[ iLevel ][ iDimension ] =
                            2 * mNumberOfBasisPerDimensionIncludingPadding[ iLevel - 1 ][ iDimension ] + PQR[ iDimension ];

                    // Calculate number of nodes on this level
                    tNumberOfBasis *= mNumberOfBasisPerDimensionIncludingPadding[ iLevel - 1 ][ iDimension ];
                }

                // Add number of nodes to offset table
                mBasisLevelOffset[ iLevel ] = mBasisLevelOffset[ iLevel - 1 ] + tNumberOfBasis;
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * Private function calculates the mMySubdomainOffset lookup table
         *
         * @return void
         */
        void calculate_subdomain_offset()
        {
            Matrix< DDLUMat > tIJK = mBackgroundMesh->get_subdomain_offset_of_proc();

            for( uint iLevel = 0; iLevel < gMaxNumberOfLevels; iLevel++ )
            {
                for( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    mMySubdomainOffset[ iLevel ][ iDimension ] = tIJK(iDimension, iLevel );
                }
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * Internal function. Asks the background mesh for number of elements
         * per direction, stores result in  mAuraNumberOfElementsPerDimension
         *
         * @return void
         */
        void get_number_of_elements_per_dimension()
        {
            // get elements per level from background mesh
            Matrix< DDLUMat > tMat = mBackgroundMesh->get_number_of_elements_per_direction();

            // convert matrix to fixed size array
            for( uint iLevel = 0; iLevel < gMaxNumberOfLevels; iLevel++ )
            {
                for( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    mNumberOfElementsPerDimensionIncludingAura[ iLevel ][ iDimension ] = tMat(iDimension, iLevel );
                }
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * calculates XZY coordinates for each basis
         *
         * @return void
         */
        void calculate_basis_coordinates() override
        {
            // get domain dimensions from settings
            Matrix< DDRMat > tDomainDimensions = mParameters->get_domain_dimensions();

            // get number of elements on coarsest level from settings
            Matrix< DDLUMat > tNumberOfElements = mParameters->get_number_of_elements_per_dimension();

            // calculate step width
            real tDeltaX[ gMaxNumberOfLevels ][ N ];

            // calculate width for first level
            for( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                tDeltaX[ 0 ][ iDimension ] = tDomainDimensions(iDimension ) / ( ( real ) ( tNumberOfElements(iDimension ) ) );
            }

            // loop over all higher levels
            for( uint iLevel = 1; iLevel < gMaxNumberOfLevels; iLevel++ )
            {
                for( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    tDeltaX[ iLevel ][ iDimension ] = 0.5 * tDeltaX[ iLevel-1 ][ iDimension ];
                }
            }

            // get domain offset
            Matrix< DDRMat > tParametersOffset = mParameters->get_domain_offset();

            // domain offset
            real tOffset[ N ];

            // get coordinates from background mesh
            Matrix< DDRMat > tOffsetCoords = mBackgroundMesh->get_domain_offset();

            // unflatten coordinates to a normal array
            for( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                tOffset[ iDimension ] = tOffsetCoords(iDimension );
            }

            // coordinate shift for B-Spline
            // 0.5 - 0.5*P
            real tShift[ gMaxNumberOfLevels ][ N ];
            for ( uint iLevel = 0; iLevel < gMaxNumberOfLevels; iLevel++ )
            {
                for ( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    tShift[ iLevel ][ iDimension ] = 0.5 * ( PQR[ iDimension ] + 1 ) - std::pow( 2, iLevel ) * PQR[ iDimension ];
                }
            }

            // loop over all nodes
            for( auto tBasis : mAllBasisOnProc )
            {
                // get ijk position of node
                const luint* tIJK = tBasis->get_ijk();

                // get level of node
                luint tLevel = tBasis->get_level();

                // array containing coordinate
                real tXYZ[ N ];

                // loop over all dimensions
                for( uint iDimension = 0; iDimension < N; iDimension++ )
                {
                    tXYZ[ iDimension ] = ( tShift[ tLevel ][ iDimension ] + ( real ) ( tIJK[ iDimension ]
                                                                      + mMySubdomainOffset[ tLevel ][ iDimension ] ) )
                                                                      * tDeltaX[ tLevel ][ iDimension ] + tOffset[ iDimension ];
                }

                // write XYZ coordinate into node
                tBasis->set_xyz( tXYZ );
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * Creates a new B-spline element based on the provided background element, N, and P
         *
         * @param aBackgroundElement Background element
         * @return Created B-spline element
         */
        Element * create_element( Background_Element_Base* aBackgroundElement ) override
        {
            // Check for dimension less than or equal to 3 and degree less than or equal to 3
            MORIS_ERROR( P <= 3 and Q <= 3 and R <= 3, "Don't know how to create B-Spline element.");

            // Create element
            Element * aBSplineElement = new BSpline_Element< P, Q, R >( aBackgroundElement, mActivationPattern );

            // Return element
            return aBSplineElement;
        }
    };

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    Basis * BSpline_Mesh< 1, 1, 0 >::create_basis(
            const luint* aIJK,
            uint         aLevel,
            uint         aOwner )
    {
        return new BSpline< 2, 9, 8 >( aIJK, aLevel, aOwner );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    Basis * BSpline_Mesh< 2, 2, 0 >::create_basis(
            const luint* aIJK,
            uint         aLevel,
            uint         aOwner )
    {
        return new BSpline< 2, 16, 8 >( aIJK, aLevel, aOwner );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    Basis * BSpline_Mesh< 3, 3, 0 >::create_basis(
            const luint* aIJK,
            uint         aLevel,
            uint         aOwner )
    {
        return new BSpline< 2, 25, 8 >( aIJK, aLevel, aOwner );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    Basis * BSpline_Mesh< 4, 4, 0 >::create_basis(
            const luint* aIJK,
            uint         aLevel,
            uint         aOwner )
    {
        return new BSpline< 2, 36, 8 >( aIJK, aLevel, aOwner );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    Basis * BSpline_Mesh< 5, 5, 0 >::create_basis(
            const luint* aIJK,
            uint         aLevel,
            uint         aOwner )
    {
        return new BSpline< 2, 49, 8 >( aIJK, aLevel, aOwner );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    Basis * BSpline_Mesh< 1, 1, 1 >::create_basis(
            const luint* aIJK,
            uint         aLevel,
            uint         aOwner )
    {
        return new BSpline< 3, 27, 26 >( aIJK, aLevel, aOwner );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    Basis*
    BSpline_Mesh< 2, 2, 2 >::create_basis(
            const luint* aIJK,
            uint         aLevel,
            uint         aOwner )
    {
        return new BSpline< 3, 64, 26 >( aIJK, aLevel, aOwner );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    Basis * BSpline_Mesh< 3, 3, 3 >::create_basis(
            const luint* aIJK,
            uint         aLevel,
            uint         aOwner )
    {
        return new BSpline< 3, 125, 26 >( aIJK, aLevel, aOwner );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    Basis * BSpline_Mesh< 4, 4, 4 >::create_basis(
            const luint* aIJK,
            uint         aLevel,
            uint         aOwner )
    {
        return new BSpline< 3, 216, 26 >( aIJK, aLevel, aOwner );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    Basis * BSpline_Mesh< 5, 5, 5 >::create_basis(
            const luint* aIJK,
            uint         aLevel,
            uint         aOwner )
    {
        return new BSpline< 3, 343, 26 >( aIJK, aLevel, aOwner );
    }

    // ----------------------------------------------------------------------------

} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_BSPLINE_MESH_HPP_ */


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

namespace moris
{
    namespace hmr
    {
        // ----------------------------------------------------------------------------

        /**
         * \brief the BSpline_Mesh class calculates B-Splines for a given
         *  background mesh.
         *
         */
        template< uint N, uint P >
        class BSpline_Mesh : public BSpline_Mesh_Base
        {
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
                        const uint            & aActivationPattern )
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

                    this->calculate_subdomain_offset();

                    // calculate any value that can change after refinement
                    this->update_mesh();
            }

                // ----------------------------------------------------------------------------

                /**
                 * Default destructor.
                 */

                ~BSpline_Mesh()
                {
                    mActiveBasisOnProc.clear();

                    this->delete_pointers();
                }

                //------------------------------------------------------------------------------

                /**
                 * string needed for gmsh output
                 *
                 * @return std::string
                 *
                 */

                std::string get_gmsh_string();

                //------------------------------------------------------------------------------

                /**
                 * VTK ID needed for VTK output
                 *
                 * @return uint
                 */

                uint get_vtk_type();

                //------------------------------------------------------------------------------

                /**
                 * node IDs needed for VTK output
                 *
                 * @param[out] moris::Matrix< DDLUMat >
                 *
                 * @return void
                 *
                 */

                void get_basis_indices_for_vtk( Matrix< DDLUMat > & aNodes );

            protected:

                // ----------------------------------------------------------------------------

                Element * create_element( Background_Element_Base* aElement );

                // ----------------------------------------------------------------------------

                /**
                 * creates a basis depending on polynomial order and dimension
                 */
                Basis * create_basis(
                        const luint * aIJK,
                        const  uint & aLevel,
                        const  uint & aOwner );

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
                        const uint  & aLevel,
                        const luint & aI )
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
                        const uint  & aLevel,
                        const luint & aI,
                        const luint & aJ )
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
                        const uint  & aLevel,
                        const luint & aI,
                        const luint & aJ,
                        const luint & aK )
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
                    // calculate number of basis on first level
                    for( uint k = 0; k < N; ++k )
                    {
                        mNumberOfBasisPerDimensionIncludingPadding[ 0 ][ k ] =
                                mNumberOfElementsPerDimensionIncludingAura[ 0 ][ k ] + P;
                    }

                    // calculate number of basis on higher levels
                    for( uint l = 1; l < gMaxNumberOfLevels; ++l )
                    {
                        for( uint k = 0; k < N; ++k )
                        {
                            mNumberOfBasisPerDimensionIncludingPadding[ l ][ k ] =
                                    2*mNumberOfBasisPerDimensionIncludingPadding[ l-1 ][ k ] + P;
                        }
                    }

                    // calculate basis level offset
                    mBasisLevelOffset[ 0 ] = 0;

                    for( uint l = 1; l < gMaxNumberOfLevels; ++l )
                    {
                        // calculate number of nodes on this level
                        luint tNumberOfBasis = 1;
                        for( uint k = 0; k < N; ++k )
                        {
                            tNumberOfBasis *= mNumberOfBasisPerDimensionIncludingPadding[ l-1 ][ k ];
                        }

                        // add number of nodes to offset table
                        mBasisLevelOffset[ l ] = mBasisLevelOffset[ l-1 ] + tNumberOfBasis;
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

                    for( uint l = 0; l < gMaxNumberOfLevels; ++l )
                    {
                        for( uint k = 0; k < N; ++k )
                        {
                            mMySubdomainOffset[ l ][ k ] = tIJK( k, l );
                        }
                    }
                }

                // ----------------------------------------------------------------------------

                /**
                 * Internal function. Asks the background mesh for number of elements
                 * per direction, stores result in  mAuraNumberOfElementsPerDimension
                 *
                 * @return void
                 *
                 */

                void get_number_of_elements_per_dimension()
                {
                    // get elements per level from background mesh
                    Matrix< DDLUMat > tMat = mBackgroundMesh->get_number_of_elements_per_direction();

                    // convert matrix to fixed size array
                    for( uint l = 0; l < gMaxNumberOfLevels; ++l )
                    {
                        for( uint k = 0; k < N; ++k )
                        {
                            mNumberOfElementsPerDimensionIncludingAura[ l ][ k ] = tMat( k, l );
                        }
                    }
                }

                // ----------------------------------------------------------------------------

                /**
                 * calculates XZY coordinates for each basis
                 *
                 * @return void
                 */

                void calculate_basis_coordinates()
                {
                    // get domain dimensions from settings
                    Matrix< DDRMat > tDomainDimensions = mParameters->get_domain_dimensions();

                    // get number of elements on coarsest level from settings
                    Matrix< DDLUMat > tNumberOfElements = mParameters->get_number_of_elements_per_dimension();

                    // calculate step width
                    real tDeltaX[ gMaxNumberOfLevels ][ N ];

                    // calculate width for first level
                    for( uint k = 0; k < N; ++k )
                    {
                        tDeltaX[ 0 ][ k ] = tDomainDimensions( k ) / ( ( real ) ( tNumberOfElements( k ) ) );
                    }

                    // loop over all higher levels
                    for( uint l = 1; l < gMaxNumberOfLevels; ++l )
                    {
                        for( uint k = 0; k < N; ++k )
                        {
                            tDeltaX[ l ][ k ] = 0.5*tDeltaX[ l-1 ][ k ];
                        }
                    }

                    // get domain offset
                    Matrix< DDRMat > tParametersOffset = mParameters->get_domain_offset();

                    // domain offset
                    real tOffset[ N ];

                    // get coordinates from background mesh
                    Matrix< DDRMat > tOffsetCoords = mBackgroundMesh->get_domain_offset();

                    // unflatten coordinates to a normal array
                    for( uint k = 0; k < N; ++k )
                    {
                        tOffset[ k ] = tOffsetCoords( k );
                    }

                    // coordinate shift for B-Spline
                    // 0.5 - 0.5*P

                    real tOff = real( P );

                    real tShift[ gMaxNumberOfLevels ];

                    for( uint i = 0; i < gMaxNumberOfLevels; ++i )
                    {
                        tShift[ i ] = 0.5 -tOff + 0.5*(real)  P;
                        tOff *= 2;
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
                        for( uint k = 0; k < N; ++k )
                        {
                            tXYZ[ k ] = ( tShift[ tLevel ] + ( real ) ( tIJK[ k ]
                                                                              + mMySubdomainOffset[ tLevel ][ k ] ) )
                                                                              * tDeltaX[ tLevel ][ k ] + tOffset[ k ];
                        }

                        // write XYZ coordinate into node
                        tBasis->set_xyz( tXYZ );
                    }
                }

                // ----------------------------------------------------------------------------
        };

        // ----------------------------------------------------------------------------

        template < uint N, uint P >
        Element * BSpline_Mesh< N, P >::create_element( Background_Element_Base* aElement )
        {
            MORIS_ERROR( false, "Don't know how to create B-Spline element.");
            return nullptr;
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element * BSpline_Mesh< 2, 1 >::create_element( Background_Element_Base* aElement )
        {
            Element * aBSplineElement = new BSpline_Element< 2, 4 >( aElement, mActivationPattern );

            return aBSplineElement;
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element * BSpline_Mesh< 2, 2 >::create_element( Background_Element_Base* aElement )
        {
            Element * aBSplineElement = new BSpline_Element< 2,9 >( aElement, mActivationPattern );

            return aBSplineElement;
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element * BSpline_Mesh< 2, 3 >::create_element( Background_Element_Base* aElement )
        {
            Element * aBSplineElement = new BSpline_Element< 2, 16 >( aElement, mActivationPattern );

            return aBSplineElement;
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element * BSpline_Mesh< 2, 4 >::create_element( Background_Element_Base* aElement )
        {
            Element * aBSplineElement = new BSpline_Element< 2, 25 >( aElement, mActivationPattern );

            return aBSplineElement;
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element * BSpline_Mesh< 2, 5 >::create_element( Background_Element_Base* aElement )
        {
            Element * aBSplineElement = new BSpline_Element< 2, 36 >( aElement, mActivationPattern );

            return aBSplineElement;
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element * BSpline_Mesh< 3, 1 >::create_element( Background_Element_Base* aElement )
        {
            Element * aBSplineElement = new BSpline_Element< 3, 8 >( aElement, mActivationPattern );

            return aBSplineElement;
        }
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element * BSpline_Mesh< 3, 2 >::create_element( Background_Element_Base* aElement )
        {
            Element * aBSplineElement = new BSpline_Element< 3, 27 >( aElement, mActivationPattern );

            return aBSplineElement;
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element * BSpline_Mesh< 3, 3 >::create_element( Background_Element_Base* aElement )
        {
            Element * aBSplineElement = new BSpline_Element< 3, 64 >( aElement, mActivationPattern );

            return aBSplineElement;
        }
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element * BSpline_Mesh< 3, 4 >::create_element( Background_Element_Base* aElement )
        {
            Element * aBSplineElement = new BSpline_Element< 3, 125 >( aElement, mActivationPattern );

            return aBSplineElement;
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        Element * BSpline_Mesh< 3, 5 >::create_element( Background_Element_Base* aElement )
        {
            Element * aBSplineElement = new BSpline_Element< 3, 216 >( aElement, mActivationPattern );

            return aBSplineElement;
        }

        // ----------------------------------------------------------------------------

        template < uint N, uint P >
        Basis * BSpline_Mesh< N, P >::create_basis( const luint * aIJK,
                const uint  & aLevel,
                const uint  & aOwner )
        {
            MORIS_ERROR( false, "Don't know how to create B-Spline element.");

            return nullptr;
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        Basis * BSpline_Mesh< 2, 1 >::create_basis( const luint * aIJK,
                const uint  & aLevel,
                const uint  & aOwner )
        {
            return new BSpline< 2, 9, 8 >( aIJK, aLevel, aOwner );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        Basis * BSpline_Mesh< 2, 2 >::create_basis(
                const luint * aIJK,
                const  uint & aLevel,
                const  uint & aOwner )
        {
            return new BSpline< 2, 16, 8 >( aIJK, aLevel, aOwner );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        Basis * BSpline_Mesh< 2, 3 >::create_basis(
                const luint * aIJK,
                const  uint & aLevel,
                const  uint & aOwner )
        {
            return new BSpline< 2, 25, 8 >( aIJK, aLevel, aOwner );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        Basis * BSpline_Mesh< 2, 4 >::create_basis(
                const luint * aIJK,
                const  uint & aLevel,
                const  uint & aOwner )
        {
            return new BSpline< 2, 36, 8 >( aIJK, aLevel, aOwner );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        Basis * BSpline_Mesh< 2, 5 >::create_basis(
                const luint * aIJK,
                const  uint & aLevel,
                const  uint & aOwner )
        {
            return new BSpline< 2, 49, 8 >( aIJK, aLevel, aOwner );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        Basis * BSpline_Mesh< 3, 1 >::create_basis(
                const luint * aIJK,
                const  uint & aLevel,
                const  uint & aOwner )
        {
            return new BSpline< 3, 27, 26 >( aIJK, aLevel, aOwner );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        Basis*
        BSpline_Mesh< 3, 2 >::create_basis(
                const luint * aIJK,
                const  uint & aLevel,
                const  uint & aOwner )
        {
            return new BSpline< 3, 64, 26 >( aIJK, aLevel, aOwner );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        Basis * BSpline_Mesh< 3, 3 >::create_basis(
                const luint * aIJK,
                const  uint & aLevel,
                const  uint & aOwner )
        {
            return new BSpline< 3, 125, 26 >( aIJK, aLevel, aOwner );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        Basis * BSpline_Mesh< 3, 4 >::create_basis(
                const luint * aIJK,
                const  uint & aLevel,
                const  uint & aOwner )
        {
            return new BSpline< 3, 216, 26 >( aIJK, aLevel, aOwner );
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        Basis * BSpline_Mesh< 3, 5 >::create_basis(
                const luint * aIJK,
                const  uint & aLevel,
                const  uint & aOwner )
        {
            return new BSpline< 3, 343, 26 >( aIJK, aLevel, aOwner );
        }

        // ----------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_BSPLINE_MESH_HPP_ */


/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_BSpline_Element.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_BSPLINE_ELEMENT_HPP_
#define SRC_HMR_CL_HMR_BSPLINE_ELEMENT_HPP_

#include "cl_HMR_Background_Element_Base.hpp"    //HMR/src
#include "cl_HMR_BSpline.hpp"                    //HMR/src
#include "cl_HMR_Element.hpp"                    //HMR/src
#include "moris_typedefs.hpp"                          //COR/src

namespace moris::hmr
{
    //------------------------------------------------------------------------------

    /**
     * B-spline element class
     *
     * @tparam P Polynomial degree in x-direction
     * @tparam Q Polynomial degree in y-direction
     * @tparam R Polynomial degree in z-direction
     */
    template< uint P, uint Q, uint R >
    class BSpline_Element : public Element
    {
        //! Number of dimensions
        static constexpr uint N = ( P > 0 ) + ( Q > 0 ) + ( R > 0 );

        //! Number of bases
        static constexpr uint B = ( P + 1 ) * ( Q + 1 ) * ( R + 1 );

        //! Container of bases in each direction, used to avoid repeated conditionals
        static constexpr uint BXYZ[ 3 ] = { P + 1, Q + 1, R + 1 };

        //! pointer to nodes
        Basis** mBasis;

      public:
        /**
         * default Lagrange Element constructor
         */
        BSpline_Element(
                Background_Element_Base* aElement,
                uint                     aActivationPattern )
                : Element( aElement, aActivationPattern )
        {
        }

        //------------------------------------------------------------------------------

        /**
         * default destructor
         */
        ~BSpline_Element() override
        {
            this->delete_basis_container();
        }

        //------------------------------------------------------------------------------

        void
        init_basis_container() override
        {
            MORIS_ASSERT( !mHaveBasis, "Basis container of element already initiated" );

            mHaveBasis = true;
            mBasis     = new Basis*[ B ];
            for ( uint iBasisIndex = 0; iBasisIndex < B; iBasisIndex++ )
            {
                mBasis[ iBasisIndex ] = nullptr;
            }
        }

        //------------------------------------------------------------------------------

        void
        delete_basis_container() override
        {
            if ( mHaveBasis )
            {
                mHaveBasis = false;
                delete[] mBasis;
            }
        }

        //------------------------------------------------------------------------------

        /**
         * get pointer to basis
         *
         * @param[in]    aIndex   element local index of basis
         *
         * @return Basis* pointer to Lagrange node
         *
         */
        Basis*
        get_basis( uint aIndex ) override
        {
            MORIS_ASSERT( mHaveBasis, "can't return basis if container is not initialized" );
            return mBasis[ aIndex ];
        }

        const Basis*
        get_basis( uint aIndex ) const override
        {
            MORIS_ASSERT( mHaveBasis, "can't return basis if container is not initialized" );
            return mBasis[ aIndex ];
        }
        //------------------------------------------------------------------------------

        /**
         * set pointer of basis to specified index and object
         *
         * @param[in]    aIndex  element local index of node
         * @param[in]    aBasis  pointer to Lagrange node
         *
         * @return void
         *
         */
        void
        insert_basis(
                uint   aIndex,
                Basis* aBasis ) override
        {
            MORIS_ASSERT( mHaveBasis, "can't insert basis if container is not initialized" );
            mBasis[ aIndex ] = aBasis;
        }

        // -----------------------------------------------------------------------------

        /**
         * Creates all bases on the coarsest level. Called by B-spline mesh.
         *
         * @param aAllElementsOnProc Cell containing all B-spline elements including the aura
         * @return Number of created bases
         */
        luint create_basis_on_level_zero( Vector< Element* >& aAllElementsOnProc ) override;

        //------------------------------------------------------------------------------

        /**
         * Creates nodes for children of refined elements. Called by B-spline mesh.
         *
         * @param aAllElementsOnProc Cell containing all B-spline elements including the aura
         * @return Number of created bases
         */
        luint create_basis_for_children( Vector< Element* >& aAllElementsOnProc ) override;

        // -----------------------------------------------------------------------------

        /**
         * for debugging
         *
         * @return void
         */
        void
        print_connectivity() override
        {
            std::fprintf( stdout,
                    "connectivity of element %4lu ( ID %4lu, parent %4lu ):\n",
                    (long unsigned int)mElement->get_hmr_index( mActivationPattern ),
                    (long unsigned int)mElement->get_hmr_id(),
                    (long unsigned int)mElement->get_parent()->get_hmr_id() );
            for ( uint iBasisIndex = 0; iBasisIndex < B; iBasisIndex++ )
            {
                // get node
                Basis* tNode = this->get_basis( iBasisIndex );
                std::fprintf( stdout,
                        "    %2u :  Basis %lu , ID %lu, MEM %lu \n",
                        (unsigned int)iBasisIndex,
                        (long unsigned int)tNode->get_hmr_index(),
                        (long unsigned int)tNode->get_hmr_id(),
                        (long unsigned int)tNode->get_memory_index() );
            }
            std::fprintf( stdout, "\n" );
        }

        //------------------------------------------------------------------------------

        /**
         * string needed for gmsh output
         *
         * @return std::string
         *
         */
        std::string get_gmsh_string() override;

        //------------------------------------------------------------------------------

        /**
         * VTK ID needed for VTK output
         *
         * @return uint
         */
        uint get_vtk_type() override;

        //------------------------------------------------------------------------------

        /**
         * Gets the basis indices needed for basis output
         *
         * @param aBasis Basis indices
         */
        void get_basis_indices_for_vtk( Matrix< DDLUMat >& aBasisIndices ) override;

        //------------------------------------------------------------------------------

        /**
         * returns the ijk position of a given basis
         *
         * @param[in]  aBasisNumber   element local number of basis
         * @param[out] aIJK           proc local ijk position of this basis
         *
         * @return void
         */
        void
        get_ijk_of_basis(
                uint   aBasisNumber,
                luint* aIJK ) override
        {
            // Get position of element on background mesh
            const luint* tElIJK = mElement->get_ijk();

            // Get element local coordinate with element offset
            for ( uint iDimension = 0; iDimension < N; iDimension++ )
            {
                uint tIJK = aBasisNumber;
                for ( uint iPreviousDimension = 0; iPreviousDimension < iDimension; iPreviousDimension++ )
                {
                    tIJK /= BXYZ[ iPreviousDimension ];
                }

                aIJK[ iDimension ] = ( tIJK % BXYZ[ iDimension ] ) + tElIJK[ iDimension ];
            }
        }

        //------------------------------------------------------------------------------

        /**
         * Links each basis of an element with neighbor basis.
         *
         * @param[inout] aAllElementsOnProc   cell containing all Bspline
         *                                    elements including the aura
         * @return void
         */
        void link_basis_with_neighbors( Vector< Element* >& aAllElementsOnProc ) override;

        //------------------------------------------------------------------------------

        /**
         * Refines this element.
         *
         * @param aAllElementsOnProc Cell containing all B-spline elements including the aura
         * @return Number of created bases
         */
        luint refine( Vector< Element* >& aAllElementsOnProc ) override;

        //------------------------------------------------------------------------------

        /**
         * Returns the geometry type of the element
         */
        mtk::Geometry_Type get_geometry_type() const override;

        //------------------------------------------------------------------------------

        /**
         * Returns the number of basis connected to this element.
         */
        uint
        get_number_of_vertices() const override
        {
            return B;
        }

        /**
         * Flags all bases in this element
         */
        void
        flag_all_bases() override
        {
            if ( mHaveBasis )
            {
                for ( uint iNodeIndex = 0; iNodeIndex < B; iNodeIndex++ )
                {
                    mBasis[ iNodeIndex ]->flag();
                }
            }
        }

        //------------------------------------------------------------------------------

      protected:
        //------------------------------------------------------------------------------

        /**
         * create new node at position
         *
         * @param[in]    aBasisNumber   element local index of new basis
         *
         * @return void
         */
        void create_basis( uint aBasisNumber ) override;

        //------------------------------------------------------------------------------

        /**
         * Refines a basis of this element
         *
         * @param aBasisNumber Index of the basis to refine
         * @return Number of created bases
         */
        luint refine_basis( uint aBasisNumber ) override;

        //------------------------------------------------------------------------------

        /**
         * returns the interpolation order of this element
         */
        mtk::Interpolation_Order get_interpolation_order() const override;

        //------------------------------------------------------------------------------

        /**
         * returns a Mat with the basis coords
         */
        Matrix< DDRMat >
        get_vertex_coords() const override
        {
            MORIS_ERROR( false, "get_vertex_coords(): to make this Bspline_Element function work, turn on calculate_basis_coordinates() in collect_basis()" );

            Matrix< DDRMat > aCoords( B, N );
            for ( uint iBasisIndex = 0; iBasisIndex < B; iBasisIndex++ )
            {
                const real* tXYZ = mBasis[ iBasisIndex ]->get_xyz();

                for ( uint iCoordinateIndex = 0; iCoordinateIndex < N; iCoordinateIndex++ )
                {
                    aCoords( iBasisIndex, iCoordinateIndex ) = tXYZ[ iCoordinateIndex ];
                }
            }
            return aCoords;
        }

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------

    template< uint P, uint Q, uint R >
    inline std::string
    BSpline_Element< P, Q, R >::get_gmsh_string()
    {
        std::string aString = "GMSH not implemented for this element";
        return aString;
    }

    //------------------------------------------------------------------------------

    template< uint P, uint Q, uint R >
    inline uint
    BSpline_Element< P, Q, R >::get_vtk_type()
    {
        // this element has no VTK id
        return 0;
    }

    //------------------------------------------------------------------------------

    /**
     * Gets the node IDs needed for VTK output
     *
     * @tparam P Polynomial degree in x-direction
     * @tparam Q Polynomial degree in y-direction
     * @tparam R Polynomial degree in z-direction
     * @param aBasis
     */
    template< uint P, uint Q, uint R >
    void
    BSpline_Element< P, Q, R >::get_basis_indices_for_vtk( Matrix< DDLUMat >& aBasis )
    {
        // Loop over all bases
        for ( uint iBasisIndex = 0; iBasisIndex < B; iBasisIndex++ )
        {
            aBasis( iBasisIndex ) = mBasis[ iBasisIndex ]->get_memory_index();
        }
    }

    //------------------------------------------------------------------------------
    template< uint P, uint Q, uint R >
    void
    BSpline_Element< P, Q, R >::create_basis( uint aBasisNumber )
    {
        MORIS_ERROR( false, "Don't know how to create B-Spline." );
    }

    //------------------------------------------------------------------------------

    template< uint P, uint Q, uint R >
    luint
    BSpline_Element< P, Q, R >::create_basis_on_level_zero( Vector< Element* >& aAllElementsOnProc )
    {
        MORIS_ERROR( false, "Don't know how to create B-Splines on level zero." );
        return 0;
    }

    //------------------------------------------------------------------------------

    template< uint P, uint Q, uint R >
    luint
    BSpline_Element< P, Q, R >::create_basis_for_children( Vector< Element* >& aAllElementsOnProc )
    {
        MORIS_ERROR( false, "Don't know how to create B-Splines for children." );
        return 0;
    }

    //------------------------------------------------------------------------------

    template< uint P, uint Q, uint R >
    inline void
    BSpline_Element< P, Q, R >::link_basis_with_neighbors(
            Vector< Element* >& aAllElementsOnProc )
    {
        MORIS_ERROR( false,
                "Link basis with neighbors not available for this element." );
    }

    //------------------------------------------------------------------------------

    template< uint P, uint Q, uint R >
    luint
    BSpline_Element< P, Q, R >::refine_basis( uint aBasisNumber )
    {
        MORIS_ERROR( false, "refine_basis() not available for this element." );
        return 0;
    }

    //------------------------------------------------------------------------------

    template< uint P, uint Q, uint R >
    luint
    BSpline_Element< P, Q, R >::refine(
            Vector< Element* >& aAllElementsOnProc )
    {
        MORIS_ERROR( false, "refine() not available for this element." );
        return 0;
    }

    //------------------------------------------------------------------------------

    template< uint P, uint Q, uint R >
    inline mtk::Geometry_Type
    BSpline_Element< P, Q, R >::get_geometry_type() const
    {
        MORIS_ERROR( false, "get_geometry_type() not available for this element." );
        return mtk::Geometry_Type::UNDEFINED;
    }

    //------------------------------------------------------------------------------

    template< uint P, uint Q, uint R >
    inline mtk::Interpolation_Order
    BSpline_Element< P, Q, R >::get_interpolation_order() const
    {
        MORIS_ERROR( false, "get_interpolation_order() not available for this element." );
        return mtk::Interpolation_Order::UNDEFINED;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::hmr

#endif /* SRC_HMR_CL_HMR_BSPLINE_ELEMENT_HPP_ */

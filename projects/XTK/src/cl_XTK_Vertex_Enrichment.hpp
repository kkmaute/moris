/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Vertex_Enrichment.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_VERTEX_ENRICHMENT_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_VERTEX_ENRICHMENT_HPP_

#include "cl_Matrix.hpp"
#include "containers.hpp"
#include "assert.h"
#include <unordered_map>
#include "cl_MTK_Vertex_Interpolation_XTK_Impl.hpp"

using namespace moris;

namespace moris::xtk
{
    class Enriched_Interpolation_Mesh;
}

namespace moris::xtk
{
    class Vertex_Enrichment : public mtk::Vertex_Interpolation
    {
        friend class Enriched_Interpolation_Mesh;

        //------------------------------------------------------------------------------

      protected:
        moris::moris_index         mNodeIndex;
        Matrix< IndexMat >         mBasisIndices;
        Matrix< IndexMat >         mBasisIds;
        Matrix< IndexMat >         mBasisOwners;
        Matrix< DDRMat >           mBasisWeights;
        mtk::Vertex_Interpolation* mBaseVertexInterp;
        IndexMap                   mBasisMap; /*From basis to local index*/

                                              //------------------------------------------------------------------------------

      public:
        /**
         * @brief Default Constructor
         */
        Vertex_Enrichment();

        //------------------------------------------------------------------------------
        // Functions to interface with vertex interpolation in MTK
        //------------------------------------------------------------------------------

        /**
         * returns the IDs of the interpolation coefficients
         */
        Matrix< IdMat >
        get_ids() const;

        //------------------------------------------------------------------------------

        /**
         * returns the indices of the interpolation coefficients
         */
        Matrix< IndexMat >
        get_indices() const;

        //------------------------------------------------------------------------------

        /**
         * returns the proc owners of the IDs of this vertex
         */
        Matrix< IdMat >
        get_owners() const;

        //------------------------------------------------------------------------------
        /**
         * set the interpolation weights
         */
        void
        set_weights( const Matrix< DDRMat >& aWeights )
        {
            MORIS_ERROR( 0, "set_weights not implemented in xtk vertex interpolation" );
        }

        //------------------------------------------------------------------------------

        /**
         * returns the interpolation weights
         */
        const Matrix< DDRMat >*
        get_weights() const;

        //------------------------------------------------------------------------------

        /**
         * set the coefficient objects
         */
        void
        set_coefficients( Vector< mtk::Vertex* >& aCoefficients )
        {
            MORIS_ERROR( 0, "set_coefficients not implemented in xtk vertex interpolation" );
        }

        //------------------------------------------------------------------------------

        /**
         * returns the pointers to the coefficient objects
         */
        Vector< mtk::Vertex* >&
        get_coefficients()
        {
            MORIS_ERROR( 0, "get_coefficients not implemented in xtk vertex interpolation" );
            return mCoefficients;
        }

        //------------------------------------------------------------------------------

        /**
         * returns the pointers to the coefficient objects (const version)
         */
        const Vector< mtk::Vertex* >&
        get_coefficients() const
        {
            MORIS_ERROR( 0, "get_coefficients not implemented in xtk vertex interpolation" );
            return mCoefficients;
        }

        //------------------------------------------------------------------------------

        /**
         * returns the number of coefficients attributed to this basis
         */
        uint
        get_number_of_coefficients() const
        {
            MORIS_ERROR( 0, "get_number_of_coefficients not implemented in xtk vertex interpolation" );
            return 0;
        }

        //------------------------------------------------------------------------------

        void
        set_node_index( moris::moris_index aNodeIndex );

        //------------------------------------------------------------------------------

        mtk::Vertex_Interpolation const *
        get_base_vertex_interpolation() const;

        //------------------------------------------------------------------------------

        /*
         * Add the basis information which includes weights, enrichment level, and basis index.
         * There is no "smartness" in this function. Duplicates should have been removed prior to call
         * An assertion will catch duplicates in debug mode
         */
        void
        add_basis_information(
                Matrix< IndexMat > const & aBasisIndices,
                Matrix< IndexMat > const & aBasisId );

        void
        add_basis_weights(
                Matrix< IndexMat > const & aBasisIndices,
                Matrix< DDRMat > const &   aBasisWeight );

        void
        add_basis_owners(
                Matrix< IndexMat > const & aBasisIndices,
                Matrix< IndexMat > const & aBasisOwners );

        void
        add_base_vertex_interpolation( mtk::Vertex_Interpolation* aBaseVertInterp );

        IndexMap&
        get_basis_map();

        uint
        get_num_bases_in_map() const;

        uint
        local_basis_index( uint aBasisIndex ) const;

        void
        condense_out_basis_with_0_weight();

        Matrix< IndexMat > const &
        get_basis_indices() const;

        Matrix< IndexMat > const &
        get_basis_ids() const;

        Matrix< DDRMat > const &
        get_basis_weights() const;

        Matrix< DDRMat >&
        get_basis_weights();

        bool
        basis_exists_in_enrichment( moris_index aBasisIndex ) const;

        //------------------------------------------------------------------------------

        bool
        has_interpolation() const;

        //------------------------------------------------------------------------------

        size_t
        capacity();

        //------------------------------------------------------------------------------

        void
        print() const;

        //------------------------------------------------------------------------------

    };    // class Vertex_Enrichment

    //------------------------------------------------------------------------------
    // OPERATORS ASSOCIATED WITH THIS CLASS
    //------------------------------------------------------------------------------

    /*
     * They are considered the same if they have the same basis weights and basis coefficients,
     * not necessarily in the same order
     */
    inline bool
    operator==( const Vertex_Enrichment& aA, const Vertex_Enrichment& aB )
    {
        // check basis indices ...............

        // get basis indices of aA,aB
        Matrix< IndexMat > const & tBasisIndicesA = aA.get_basis_indices();
        Matrix< IndexMat > const & tBasisIndicesB = aB.get_basis_indices();

        uint tSizeA = tBasisIndicesA.numel();
        uint tSizeB = tBasisIndicesB.numel();

        // if they do not have the same number of basis they cannot be equal
        if ( tSizeA != tSizeB )
        {
            return false;
        }

        // if their index maps have different length something must be wrong
        MORIS_ERROR(
                ( aA.get_num_bases_in_map() == aB.get_num_bases_in_map() ) && ( aA.get_num_bases_in_map() == tSizeA ),
                "xtk::Vertex_Enrichment::operator== - "
                "Basis maps of Vertex Enrichments compared have same number of indices, but different size of basis index maps." );

        // iterate through basis in aA and ask aB if they exist in their basis
        for ( uint iBasis = 0; iBasis < tSizeA; iBasis++ )
        {
            moris_index tCurrentBasisIndex   = tBasisIndicesA( iBasis );
            bool        tBasisIndexExistsInB = aB.basis_exists_in_enrichment( tCurrentBasisIndex );

            if ( !tBasisIndexExistsInB )
            {
                return false;
            }
        }

        // check weights ...............

        // get basis weights of aA,aB
        Matrix< DDRMat > const & tWeightsA = aA.get_basis_weights();
        Matrix< DDRMat > const & tWeightsB = aB.get_basis_weights();

        // check size
        MORIS_ASSERT( tWeightsA.numel() == tSizeA, "xtk::Vertex_Enrichment::operator== - Size of weights and indices don't match for A." );
        MORIS_ASSERT( tWeightsB.numel() == tSizeB, "xtk::Vertex_Enrichment::operator== - Size of weights and indices don't match for B." );

        // check the weights
        for ( uint iWeight = 0; iWeight < tSizeA; iWeight++ )
        {
            // get the basis index
            moris_index tBasisIndex = tBasisIndicesA( iWeight );

            // find the position of this index in B
            uint tBasisIndexPositionInB = aB.local_basis_index( tBasisIndex );

            // get the two weights
            real tWeightA = tWeightsA( iWeight );
            real tWeightB = tWeightsB( tBasisIndexPositionInB );

            // check if the two weights are equal
            if ( std::abs( tWeightA - tWeightB ) > 10.0 * MORIS_REAL_EPS )
            {
                return false;
            }
        }

        // if all checks have passed the T-matrices are equal
        return true;
    }

    //------------------------------------------------------------------------------

    inline std::ostream&
    operator<<( std::ostream& os, const xtk::Vertex_Enrichment& dt )
    {
        Matrix< IndexMat > const & tBasisIndices = dt.get_basis_indices();
        Matrix< IndexMat > const & tBasisOwner   = dt.get_owners();
        Matrix< DDRMat > const &   tBasisWeights = dt.get_basis_weights();

        // base vertex
        mtk::Vertex_Interpolation const * tBaseVertIp       = dt.get_base_vertex_interpolation();
        Matrix< IndexMat >                tBackBasisIndices = tBaseVertIp->get_indices();
        for ( uint iBasis = 0; iBasis < tBasisIndices.numel(); iBasis++ )
        {
            os << "Basis Index: " << std::setw( 9 ) << tBasisIndices( iBasis );
            os << " | Basis Weight: " << std::setw( 9 ) << tBasisWeights( iBasis );
            os << " | Basis Owner: " << std::setw( 9 ) << tBasisOwner( iBasis );
            os << " | Back Basis Index:" << std::setw( 9 ) << tBackBasisIndices( iBasis ) << std::endl;
            os << std::endl;
        }

        return os;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::xtk

//------------------------------------------------------------------------------

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_VERTEX_ENRICHMENT_HPP_ */

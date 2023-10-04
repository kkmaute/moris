/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Vertex_Enrichment.cpp
 *
 */

#include "cl_XTK_Vertex_Enrichment.hpp"
#include "fn_equal_to.hpp"
#include "fn_unique.hpp"

using namespace moris;

namespace xtk
{
    //------------------------------------------------------------------------------

    Vertex_Enrichment::Vertex_Enrichment()
            : mNodeIndex( MORIS_INDEX_MAX )
            , mBaseVertexInterp( nullptr )
    {
    }

    //------------------------------------------------------------------------------

    Matrix< IdMat >
    Vertex_Enrichment::get_ids() const
    {
        return mBasisIds;
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat >
    Vertex_Enrichment::get_indices() const
    {
        return this->get_basis_indices();
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >*
    Vertex_Enrichment::get_weights() const
    {
        return &this->get_basis_weights();
    }

    //------------------------------------------------------------------------------

    Matrix< IdMat >
    Vertex_Enrichment::get_owners() const
    {
        return mBasisOwners;
    }

    //------------------------------------------------------------------------------

    void
    Vertex_Enrichment::set_node_index( moris::moris_index aNodeIndex )
    {
        MORIS_ASSERT( mNodeIndex == MORIS_INDEX_MAX, "Node index already set for Vertex Enrichment" );
        mNodeIndex = aNodeIndex;
    }

    //------------------------------------------------------------------------------

    void
    Vertex_Enrichment::add_basis_information( moris::Matrix< moris::IndexMat > const & aBasisIndices,
            moris::Matrix< moris::IndexMat > const &                                   aBasisId )
    {
#ifdef MORIS_HAVE_DEBUG
        // since I can't write these functions in one line, need to have ifdef
        moris::Matrix< moris::IndexMat > tUniqueBasis;
        moris::unique( aBasisIndices, tUniqueBasis );

        MORIS_ASSERT( tUniqueBasis.numel() == aBasisIndices.numel(), "duplicate basis indices detected" );
#endif

        // num basis
        uint tNumBasis = aBasisIndices.numel();

        // allocate space
        mBasisIndices.resize( tNumBasis, 1 );
        mBasisIds.resize( tNumBasis, 1 );
        mBasisWeights.resize( tNumBasis, 1 );
        mBasisOwners.resize( tNumBasis, 1 );

        // iterate to store data
        for ( uint i = 0; i < aBasisIndices.numel(); i++ )
        {
            uint tBasisLocInd             = this->local_basis_index( aBasisIndices( i ) );
            mBasisIndices( tBasisLocInd ) = aBasisIndices( i );
            mBasisIds( tBasisLocInd )     = aBasisId( i );
        }
    }

    //------------------------------------------------------------------------------

    void
    Vertex_Enrichment::add_basis_weights( moris::Matrix< moris::IndexMat > const & aBasisIndices,
            moris::Matrix< moris::DDRMat > const &                                 aBasisWeight )
    {
        for ( uint i = 0; i < aBasisIndices.numel(); i++ )
        {
            uint tBasisLocInd             = this->local_basis_index( aBasisIndices( i ) );
            mBasisWeights( tBasisLocInd ) = aBasisWeight( i );
        }
    }

    //------------------------------------------------------------------------------

    void
    Vertex_Enrichment::add_basis_owners( moris::Matrix< moris::IndexMat > const & aBasisIndices,
            moris::Matrix< moris::IndexMat > const &                              aBasisOwners )
    {
        for ( uint i = 0; i < aBasisIndices.numel(); i++ )
        {
            uint tBasisLocInd            = this->local_basis_index( aBasisIndices( i ) );
            mBasisOwners( tBasisLocInd ) = aBasisOwners( i );
        }
    }

    //------------------------------------------------------------------------------

    void
    Vertex_Enrichment::add_base_vertex_interpolation( mtk::Vertex_Interpolation* aBaseVertInterp )
    {
        mBaseVertexInterp = aBaseVertInterp;
    }

    //------------------------------------------------------------------------------

    mtk::Vertex_Interpolation const *
    Vertex_Enrichment::get_base_vertex_interpolation() const
    {
        return mBaseVertexInterp;
    }

    //------------------------------------------------------------------------------

    Mini_Map< moris::moris_index, moris::moris_index >&
    Vertex_Enrichment::get_basis_map()
    {
        return mBasisMap;
    }

    //------------------------------------------------------------------------------

    uint
    Vertex_Enrichment::get_num_bases_in_map() const
    {
        return mBasisMap.size();
    }

    //------------------------------------------------------------------------------

    uint
    Vertex_Enrichment::local_basis_index( uint aBasisIndex ) const
    {
        auto tIter = mBasisMap.find( aBasisIndex );

        MORIS_ASSERT( mBasisMap.size() > 0,
                "Vertex_Enrichment::local_basis_index() - Basis map not constructed yet." );
        MORIS_ASSERT( tIter != mBasisMap.end(),
                "Vertex_Enrichment::local_basis_index() - Provided basis index not found in map." );

        return tIter->second;
    }

    //------------------------------------------------------------------------------

    void
    Vertex_Enrichment::condense_out_basis_with_0_weight()
    {
        uint tCount = 0;

        for ( uint i = 0; i < mBasisIndices.numel(); i++ )
        {
            if ( moris::equal_to( mBasisWeights( i ), 0 ) )
            {
                mBasisMap.erase( mBasisIndices( i ) );
            }

            else
            {
                mBasisIndices( tCount ) = mBasisIndices( i );
                mBasisWeights( tCount ) = mBasisWeights( i );

                // change map index
                mBasisMap[ mBasisIndices( i ) ] = tCount;
                tCount++;
            }
        }

        // remove excess space
        mBasisIndices.resize( tCount, 1 );
        mBasisWeights.resize( tCount, 1 );
    }

    //------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat > const &
    Vertex_Enrichment::get_basis_indices() const
    {
        return mBasisIndices;
    }

    //------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat > const &
    Vertex_Enrichment::get_basis_weights() const
    {
        return mBasisWeights;
    }

    moris::Matrix< moris::IndexMat > const &
    Vertex_Enrichment::get_basis_ids() const
    {
        return mBasisIds;
    }

    //------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >&
    Vertex_Enrichment::get_basis_weights()
    {
        return mBasisWeights;
    }

    //------------------------------------------------------------------------------

    bool
    Vertex_Enrichment::basis_exists_in_enrichment( moris_index aBasisIndex ) const
    {
        return mBasisMap.find( aBasisIndex ) != mBasisMap.end();
    }

    //------------------------------------------------------------------------------

    bool
    Vertex_Enrichment::has_interpolation() const
    {
        if ( mBasisIndices.numel() == 0 )
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    //------------------------------------------------------------------------------

    size_t
    Vertex_Enrichment::capacity()
    {
        size_t tTotal = 0;
        tTotal += sizeof( mNodeIndex );
        tTotal += mBasisIndices.capacity();
        tTotal += mBasisIds.capacity();
        tTotal += mBasisOwners.capacity();
        tTotal += mBasisWeights.capacity();
        tTotal += sizeof( mBaseVertexInterp );
        // FIXME: add mBasisMap
        return tTotal;
    }

    //------------------------------------------------------------------------------

    void
    Vertex_Enrichment::print() const
    {
        // get number of entries
        uint tNumEntries = mBasisIndices.numel();

        // print header
        std::cout << "\nVertex Enrichment #" << mNodeIndex << " ( Enr. Basis Index | Weight ) = " << std::endl;
        for ( uint iEntry = 0; iEntry < tNumEntries; iEntry++ )
        {
            std::cout << mBasisIndices( iEntry ) << " | " << mBasisWeights( iEntry ) << std::endl;
        }

        // std::cout << "\n" << std::endl;
        std::cout << "" << std::endl;
    }

    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------

}    // namespace xtk

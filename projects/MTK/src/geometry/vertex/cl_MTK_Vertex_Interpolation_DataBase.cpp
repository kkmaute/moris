/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_Interpolation_DataBase.cpp
 *
 */

#include "cl_MTK_Vertex_Interpolation_DataBase.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "fn_TOL_Capacities.hpp"

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    Vertex_Interpolation_DataBase::Vertex_Interpolation_DataBase(
            moris_index aVertexIndex,
            moris_index aVertexOrder,
            mtk::Mesh*  aMesh )
            : mVertexIndex( aVertexIndex )
            , mVertexOrder( aVertexOrder )
            , mMesh( aMesh )
    {
        this->set_outward_data();
    }

    //------------------------------------------------------------------------------

    Matrix< IdMat >
    Vertex_Interpolation_DataBase::get_ids() const
    {
        // get the basis pointer and number of elements
        moris_id*   tBasisIdPtr = mMesh->get_basis_ids( mVertexIndex, mVertexOrder );
        moris_index tOffset     = mMesh->get_basis_length( mVertexIndex, mVertexOrder );

        Matrix< IdMat > tBasisIds( tBasisIdPtr, tOffset, 1, false, true );

        return tBasisIds;
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat >
    Vertex_Interpolation_DataBase::get_indices() const
    {
        // get the basis pointer and number of elements
        moris_index* tBasisIndicesPtr = mMesh->get_basis_indicies( mVertexIndex, mVertexOrder );
        moris_index  tOffset          = mMesh->get_basis_length( mVertexIndex, mVertexOrder );

        Matrix< IdMat > tBasisIndices( tBasisIndicesPtr, tOffset, 1, false, true );

        return tBasisIndices;
    }

    //------------------------------------------------------------------------------

    Matrix< IdMat >
    Vertex_Interpolation_DataBase::get_owners() const
    {
        // get the basis pointer and number of elements
        moris_id*   tBasisOwnersPtr = mMesh->get_basis_owners( mVertexIndex, mVertexOrder );
        moris_index tOffset         = mMesh->get_basis_length( mVertexIndex, mVertexOrder );

        Matrix< IdMat > tBasisOwners( tBasisOwnersPtr, tOffset, 1, false, true );

        return tBasisOwners;
    }

    //------------------------------------------------------------------------------

    void
    Vertex_Interpolation_DataBase::set_weights( const moris::Matrix< DDRMat >& aWeights )
    {
        MORIS_ERROR( 0, "set_weights not implemented in fem vertex interpolation" );
    }

    //------------------------------------------------------------------------------

    void
    Vertex_Interpolation_DataBase::set_basis_extension_flag( bool aHasExtendedBasis )
    {
        mHasExtendedBasis = aHasExtendedBasis;
    }

    //------------------------------------------------------------------------------

    bool
    Vertex_Interpolation_DataBase::contains_extended_basis() const
    {
        return mHasExtendedBasis;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >*
    Vertex_Interpolation_DataBase::get_weights() const
    {
        return &mWeights;
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< moris_index > >&
    Vertex_Interpolation_DataBase::get_basis_weight_numerators() const
    {
        return mWeightNumerators;
    }

    //------------------------------------------------------------------------------
    const Vector< Vector< moris_index > >&
    Vertex_Interpolation_DataBase::get_basis_weight_denominators() const
    {
        return mWeightDenominators;
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< real > >&
    Vertex_Interpolation_DataBase::get_basis_weight_per_leader_basis() const
    {
        return mWeightsPerLeader;
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< uint > >&
    Vertex_Interpolation_DataBase::get_basis_weight_numerator_count() const
    {
        return mWeightNumeratorIndices;
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< uint > >&
    Vertex_Interpolation_DataBase::get_basis_weight_denominator_count() const
    {
        return mWeightDenominatorIndices;
    }

    //------------------------------------------------------------------------------

    void
    Vertex_Interpolation_DataBase::set_coefficients( Vector< Vertex* >& aCoefficients )
    {
        MORIS_ERROR( 0, "set_coefficients not implemented in fem vertex interpolation" );
    }

    //------------------------------------------------------------------------------

    Vector< Vertex* >&
    Vertex_Interpolation_DataBase::get_coefficients()
    {
        MORIS_ERROR( 0, "get_coefficients not implemented in fem vertex interpolation" );
        return mCoefficients;
    }

    //------------------------------------------------------------------------------

    const Vector< Vertex* >&
    Vertex_Interpolation_DataBase::get_coefficients() const
    {
        MORIS_ERROR( 0, "get_coefficients not implemented in fem vertex interpolation" );
        return mCoefficients;
    }

    //------------------------------------------------------------------------------

    uint
    Vertex_Interpolation_DataBase::get_number_of_coefficients() const
    {
        MORIS_ERROR( 0, "get_number_of_coefficients not implemented in fem vertex interpolation" );
        return 0;
    }

    //------------------------------------------------------------------------------

    void
    Vertex_Interpolation_DataBase::set_outward_data()
    {
        this->set_basis_extension_flag( mMesh->get_basis_extension_flag( mVertexIndex, mVertexOrder ) );
        moris::real* tBasisWightsPtr = mMesh->get_basis_weights( mVertexIndex, mVertexOrder );
        moris_index  tOffset         = mMesh->get_basis_length( mVertexIndex, mVertexOrder );
        //Vector< moris::moris_index >* tBasisWeightNumeratorsPtr = mMesh->get_basis_weight_numerators( mVertexIndex, mVertexOrder );
        //Vector< moris::moris_index >* tBasisWeightDenominatorsPtr = mMesh->get_basis_weight_denominators( mVertexIndex, mVertexOrder );
        //Vector< moris::real >* tBasisWeightsPerLeaderBasisPtr = mMesh->get_basis_weights_per_leader_basis( mVertexIndex, mVertexOrder );
        
        const auto* tBasisWeightNumeratorsPtr        = mMesh->get_basis_weight_numerators( mVertexIndex, mVertexOrder );
        const auto* tBasisWeightDenominatorsPtr      = mMesh->get_basis_weight_denominators( mVertexIndex, mVertexOrder );

        const auto* tBasisWeightsPerLeaderBasisPtr   = mMesh->get_basis_weights_per_leader_basis( mVertexIndex, mVertexOrder );
        const auto* tBasisWeightsNumeratorIndices    = mMesh->get_basis_weight_numerator_indices( mVertexIndex, mVertexOrder );
        const auto* tBasisWeightsDenominatorIndices  = mMesh->get_basis_weight_denominator_indices( mVertexIndex, mVertexOrder );

        mWeights = Matrix< DDRMat >( tBasisWightsPtr, tOffset, 1, false, true );

        mWeightNumerators.resize( tOffset );
        mWeightDenominators.resize( tOffset );
        mWeightsPerLeader.resize( tOffset );
        mWeightNumeratorIndices.resize( tOffset );
        mWeightDenominatorIndices.resize( tOffset );

        for (uint i = 0; i < tOffset; ++i )
        {
            mWeightNumerators( i )         = tBasisWeightNumeratorsPtr[ i ];
            mWeightDenominators( i )       = tBasisWeightDenominatorsPtr[ i ];
            mWeightsPerLeader( i )         = tBasisWeightsPerLeaderBasisPtr[ i ];
            mWeightNumeratorIndices( i )   = tBasisWeightsNumeratorIndices[ i ];
            mWeightDenominatorIndices( i ) = tBasisWeightsDenominatorIndices[ i ];
        }
        


        // 2. In-place lambda for numerators
        // mWeightNumerators = [ this, tOffset ]() {
        //     auto*                           tPtr = mMesh->get_basis_weight_numerators( mVertexIndex, mVertexOrder );
        //     Vector< Vector< moris_index > > tResult;
        //     tResult.reserve( tOffset );
        //     for ( moris_index i = 0; i < tOffset; ++i )
        //     {
        //         tResult.push_back( tPtr[ i ] );
        //     }
        //     return tResult;
        // }();

        // mWeightNumerators = [ this, tOffset ]() {
        //     // Reference directly to the target inner vector
        //     const auto& tInnerVec = mMesh->get_basis_weight_numerators( mVertexIndex, mVertexOrder );

        //     Vector< Vector< moris_index > > tResult;
        //     tResult.reserve( 1 );
        //     tResult.push_back( tInnerVec );    // Single push or copy

        //     return tResult;
        // }();

        // mWeightNumerators = [ this, tOffset ]() {
        //     Vector< Vector< moris_index > > tResult;
        //     tResult.reserve( tOffset );

        //     // Loop through tOffset vertices starting from mVertexIndex
        //     for ( moris_index i = 0; i < tOffset; ++i )
        //     {
        //         // Fetch the inner vector for each vertex index sequentially
        //         const auto& tInnerVec = mMesh->get_basis_weight_numerators( mVertexIndex + i, mVertexOrder );

        //         tResult.push_back( tInnerVec );
        //     }

        //     return tResult;
        // }();

        // // 3. In-place lambda for denominators
        // mWeightDenominators = [ this, tOffset ]() {
        //     auto*                           tPtr = mMesh->get_basis_weight_denominators( mVertexIndex, mVertexOrder );
        //     Vector< Vector< moris_index > > tResult;
        //     tResult.reserve( tOffset );
        //     for ( moris_index i = 0; i < tOffset; ++i )
        //     {
        //         tResult.push_back( tPtr[ i ] );
        //     }
        //     return tResult;
        // }();

        // mWeightDenominators = [ this, tOffset ]() {
        //     // Reference directly to the target inner vector
        //     const auto& tInnerVec = mMesh->get_basis_weight_denominators( mVertexIndex, mVertexOrder );

        //     Vector< Vector< moris_index > > tResult;
        //     tResult.reserve( 1 );
        //     tResult.push_back( tInnerVec );    // Single push or copy

        //     return tResult;
        // }();

        // mWeightDenominators = [ this, tOffset ]() {
        //     Vector< Vector< moris_index > > tResult;
        //     tResult.reserve( tOffset );

        //     // Loop through tOffset vertices starting from mVertexIndex
        //     for ( moris_index i = 0; i < tOffset; ++i )
        //     {
        //         // Fetch the inner vector for each vertex index sequentially
        //         const auto& tInnerVec = mMesh->get_basis_weight_denominators( mVertexIndex + i, mVertexOrder );

        //         tResult.push_back( tInnerVec );
        //     }

        //     return tResult;
        // }();

        // // 4. In-place lambda for weights per leader basis
        // mWeightsPerLeader = [ this, tOffset ]() {
        //     auto*                           tPtr = mMesh->get_basis_weights_per_leader_basis( mVertexIndex, mVertexOrder );
        //     Vector< Vector< moris::real > > tResult;
        //     tResult.reserve( tOffset );
        //     for ( moris_index i = 0; i < tOffset; ++i )
        //     {
        //         tResult.push_back( tPtr[ i ] );
        //     }
        //     return tResult;
        // }();

        // mWeightsPerLeader = [ this, tOffset ]() {
        //     // Reference directly to the target inner vector
        //     const auto& tInnerVec = mMesh->get_basis_weights_per_leader_basis( mVertexIndex, mVertexOrder );

        //     Vector< Vector< real > > tResult;
        //     tResult.reserve( 1 );
        //     tResult.push_back( tInnerVec );    // Single push or copy

        //     return tResult;
        // }();

        // mWeightsPerLeader = [ this, tOffset ]() {
        //     Vector< Vector< real > > tResult;
        //     tResult.reserve( tOffset );

        //     // Loop through tOffset vertices starting from mVertexIndex
        //     for ( moris_index i = 0; i < tOffset; ++i )
        //     {
        //         // Fetch the inner vector for each vertex index sequentially
        //         const auto& tInnerVec = mMesh->get_basis_weights_per_leader_basis( mVertexIndex + i, mVertexOrder );

        //         tResult.push_back( tInnerVec );
        //     }

        //     return tResult;
        // }();

        
        
    }

    //------------------------------------------------------------------------------

    size_t
    Vertex_Interpolation_DataBase::capacity()
    {
        size_t tCapacity = 0;

        tCapacity += sizeof( mVertexIndex );
        tCapacity += sizeof( mVertexOrder );
        tCapacity += sizeof( mMesh );

        tCapacity += sizeof( mWeights ) + mWeights.capacity();

        return tCapacity;
    }
}    // namespace moris::mtk

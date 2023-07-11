/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_HMR_Helper.cpp
 *
 */

#include "cl_XTK_HMR_Helper.hpp"
#include "cl_XTK_Model.hpp"
#include "fn_XTK_convert_cell_to_map.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_WRK_perform_refinement.hpp"
#include "HMR_Globals.hpp"
#include "cl_HMR_Background_Element_Base.hpp"


using namespace moris;

namespace xtk
{
    HMR_Helper::HMR_Helper( xtk::Model* aXTKModel, moris_index aMeshIndex )
            : mXTKModelPtr( aXTKModel )
            , mBsplineMeshIndex( aMeshIndex )
    {
        // get the hmr lagrange mesh through background mesh
        mHMRLagrangeMesh = dynamic_cast< hmr::Mesh& >( mXTKModelPtr->get_background_mesh() ).get_lagrange_mesh();

        mSpatialDimension = mXTKModelPtr->get_spatial_dim();

        // get the bspline order of the mesh
        mBSplineOrder = mHMRLagrangeMesh->get_bspline_order( mBsplineMeshIndex );

        // get discretization mesh index of the background mesh
        // number of bspline coefficients per direction
        //! FIXME: need to update 0 to the correct discretization mesh index
        mNumberOfNodesPerDimension = mBSplineOrder + 1;

        // calculate number of basis per element
        mNumberOfBasis = std::pow( mNumberOfNodesPerDimension, mSpatialDimension );

        // resize the cell of enriched basis IDs, to prevent reallocation only will be overwrites
        mEnrichedBasisIDs.resize( mNumberOfBasis, gNoID );
        mEnrichedBasisOwners.resize( mNumberOfBasis, gNoID );

        // resize the L2 projection matrix
        mL2ProjectionMatrix.set_size( mNumberOfBasis, mNumberOfBasis, 0.0 );
        
        // initialize the basis index ( local ordering of the basis functions in the element based on the dimension and order)
        this->init_basis_index();

        // initialize 1D matrices to find the projection matrix in 1D
        mMatrices1D = moris::Cell< Matrix< DDRMat > >( mSpatialDimension, Matrix< DDRMat >( mNumberOfNodesPerDimension, mNumberOfNodesPerDimension ) );

        // intermediate varibale to access the enrichement data in the next line
        moris::Cell< Enrichment_Data > const & tEnrichmentData = mXTKModelPtr->get_basis_enrichment().get_enrichment_data();

        //  get the BG basis indices and their level for the corresponding SPG
        mEnrichmentData = &tEnrichmentData(mBsplineMeshIndex); 
    }

    //------------------------------------------------------------------------------------

    HMR_Helper::~HMR_Helper(){}

    //------------------------------------------------------------------------------------

    luint
    HMR_Helper::get_global_domain_id_of_cell( const mtk::Cell* aCell )
    {
        const hmr::Element* tHMRElement = dynamic_cast< const hmr::Element* >( aCell );

        // check with moris assert the result is not a nullptr
        MORIS_ASSERT( tHMRElement != nullptr, "The cell is not a hmr element" );

        // get the global domain id of the cell
        return tHMRElement->get_hmr_id();
    }

    //------------------------------------------------------------------------------------

    moris::Cell< moris_id > const &
    HMR_Helper::get_enriched_basis_id_of_cell( const mtk::Cell* aCell, moris_index aSPGIndex )
    {
        // reassign values in for the return cell
        mEnrichedBasisIDs.assign( mNumberOfBasis, gNoID );
        mEnrichedBasisOwners.assign( mNumberOfBasis, gNoID );

        // get the hmr element
        const hmr::Element* tHMRElement = dynamic_cast< const hmr::Element* >( aCell );
        
        // initailzie n index counter to count the basis that are active on the element
        uint iIndexCounter = 0;

        // fill out the cell data
        for ( uint i = 0; i < mNumberOfBasis; i++ )
        {
            // get the ith basis
            const hmr::Basis* tBasis = tHMRElement->get_basis( i );

            // if it is active add it to the cell
            if ( tBasis->is_active() )
            {
                mEnrichedBasisIDs( iIndexCounter )      = tBasis->get_index();
                mEnrichedBasisOwners( iIndexCounter++ ) = tBasis->get_owner();
            }
        }

        // get the BG basis and level for the corresponding SPG, this data is necessary to idenity the enriched basis
        moris::Cell< moris_index > const &             tBGBasisIndices = mEnrichmentData->mSubphaseGroupBGBasisIndices( aSPGIndex );
        moris::Cell< moris_index > const &             tBGBasisLevels  = mEnrichmentData->mSubphaseGroupBGBasisEnrLev( aSPGIndex );
        std::unordered_map< moris_index, moris_index > tBasisToLocalIndexMapRoot;

        // convert bg basis indices to a map
        convert_cell_to_map( tBGBasisIndices, tBasisToLocalIndexMapRoot );

        // The enriched basis indices based on the background and level
        moris::Cell< Matrix< IndexMat > > const & tBasisEnrichmentIndices = mEnrichmentData->mBasisEnrichmentIndices;

        // loop over the unenriched basis indices and convert them to enriched ones based on the background and level
        for ( uint iBasisOrd = 0; iBasisOrd < mNumberOfBasis; iBasisOrd++ )
        {
            // get the background basis in order
            moris_index tBGBasisIndex = mEnrichedBasisIDs( iBasisOrd );

            // get the basis index
            auto        tIter            = tBasisToLocalIndexMapRoot.find( tBGBasisIndex );
            moris_index tLocalBasisOrd = tIter->second;

            // Now use the local basis index to get the enrichment level
            moris_index tEnrLevel = tBGBasisLevels( tLocalBasisOrd );

            // get the enriched basis index of the root basis
            moris_index tEnrichedBasisIndex = tBasisEnrichmentIndices( tBGBasisIndex )( tEnrLevel );

            // switch the BG basis index with the enriched one
            mEnrichedBasisIDs( iBasisOrd ) = tEnrichedBasisIndex;
        }

        // transform the basis indices to basis ids
        moris::Matrix< IdMat > const & tEnrichedBasisIndexToId = mEnrichmentData->mEnrichedBasisIndexToId;

        // use the transform function to get the enriched basis ids from the enriched basis indices
        std::transform( mEnrichedBasisIDs.begin(), mEnrichedBasisIDs.end(), mEnrichedBasisIDs.begin(),    //
                [ &tEnrichedBasisIndexToId ]( moris_index aBasisIndex ) -> moris_id                       //
                {
                    return tEnrichedBasisIndexToId( aBasisIndex );
                } );

        // return the basis ids
        return mEnrichedBasisIDs;
    }

    //------------------------------------------------------------------------------------

    uint
    HMR_Helper::get_number_of_basis_per_element()
    {
        return mNumberOfBasis;
    }

    //------------------------------------------------------------------------------------

    moris::Cell< moris_id > const &
    HMR_Helper::get_bg_basis_indices_of_cell( const mtk::Cell* aCell )
    {
        // reassign values in for the return cell to be gNOID
        mEnrichedBasisIDs.assign( mNumberOfBasis, gNoID );

        // get the BG basis indices and their level for the corresponding SPG
        // get the hmr element
        const hmr::Element* tHMRElement = dynamic_cast< const hmr::Element* >( aCell );

        // initailzie n index counter to count the basis that are active on the element
        uint iIndexCounter = 0;

        // fill out the cell data
        for ( uint i = 0; i < mNumberOfBasis; i++ )
        {
            // get the basis
            const hmr::Basis* tBasis = tHMRElement->get_basis( i );

            // if it is active add it to the cell
            if ( tBasis->is_active() )
            {
                mEnrichedBasisIDs( iIndexCounter++ ) = tBasis->get_index();
            }
        }

        // return the basis indics
        return mEnrichedBasisIDs;
    }

    //------------------------------------------------------------------------------------

    Matrix< DDRMat > const &
    HMR_Helper::get_l2_projection_matrix( const mtk::Cell* aExtendedCell, moris_id aRootBsplineId )
    {
        // reset the L2 projection matrix and individual 1d matrices
        mL2ProjectionMatrix.set_size( mNumberOfBasis, mNumberOfBasis, 0.0 );

        // get the extended hmr element and the corresponding ijk
        const hmr::Element* tHMRExtenedElement = dynamic_cast< const hmr::Element* >( aExtendedCell );

        // get the ijk of the extended element
        luint tExtendedIJK[ 3 ];
        mHMRLagrangeMesh->get_background_mesh()->calc_ijk_from_global_id( tHMRExtenedElement->get_level(), tHMRExtenedElement->get_hmr_id(), tExtendedIJK );

        // NOTE: it is assumed that the root element has the same level as the extended element
        luint tRootIJK[ 3 ];
        mHMRLagrangeMesh->get_background_mesh()->calc_ijk_from_global_id( tHMRExtenedElement->get_level(), aRootBsplineId, tRootIJK );

        // loop over the dimensions and create the i,j,k 1D matrices
        for ( uint iDim = 0; iDim < mSpatialDimension; iDim++ )
        {
            // obtain the shift, a custom subtraction needs to be defined to account for negative values
            real tShift = tRootIJK[ iDim ] < tExtendedIJK[ iDim ] ? real( -tRootIJK[ iDim ] + tExtendedIJK[ iDim ] ) : -real( -tExtendedIJK[ iDim ] + tRootIJK[ iDim ] );

            // find the shift in each direction
            this->get_extention_matrix_1d( tShift, mMatrices1D( iDim ) );
        }

        // Apply a tensor product to get the final weights
        if ( mSpatialDimension == 2 )
        {
            uint b = 0;
            for ( uint l = 0; l < mNumberOfNodesPerDimension; ++l )
            {
                for ( uint k = 0; k < mNumberOfNodesPerDimension; ++k )
                {
                    uint a = 0;
                    for ( uint j = 0; j < mNumberOfNodesPerDimension; ++j )
                    {
                        for ( uint i = 0; i < mNumberOfNodesPerDimension; ++i )
                        {
                            mL2ProjectionMatrix( mBasisIndex( a ), mBasisIndex( b ) ) = mMatrices1D( 0 )( i, k ) * mMatrices1D( 1 )( j, l );
                            ++a;
                        }
                    }
                    ++b;
                }
            }
        }
        else if ( mSpatialDimension == 3 )
        {
            uint b = 0;
            for ( uint p = 0; p < mNumberOfNodesPerDimension; ++p )
            {
                for ( uint q = 0; q < mNumberOfNodesPerDimension; ++q )
                {
                    for ( uint l = 0; l < mNumberOfNodesPerDimension; ++l )
                    {
                        uint a = 0;
                        for ( uint k = 0; k < mNumberOfNodesPerDimension; ++k )
                        {
                            for ( uint j = 0; j < mNumberOfNodesPerDimension; ++j )
                            {
                                for ( uint i = 0; i < mNumberOfNodesPerDimension; ++i )
                                {
                                    mL2ProjectionMatrix( mBasisIndex( a ), mBasisIndex( b ) ) = mMatrices1D( 0 )( i, l ) * mMatrices1D( 1 )( j, q ) * mMatrices1D( 2 )( k, p );
                                    ++a;
                                }
                            }
                        }
                        ++b;
                    }
                }
            }
        }

        // return the matrix
        return mL2ProjectionMatrix;
    }

    //------------------------------------------------------------------------------------

    hmr::Background_Element_Base*
    HMR_Helper::create_background_element()
    {
        hmr::Background_Element_Base* aBackElement = nullptr;

        // create a prototype for a background element
        switch ( mSpatialDimension )
        {
            case ( 2 ):
            {
                // create a prototype for a background element
                luint tIJK[ 2 ] = { 0 };
                return new hmr::Background_Element< 2 >(
                        nullptr,
                        0,
                        tIJK,
                        0,
                        0,
                        0,
                        gNoProcID );
                break;
            }
            case ( 3 ):
            {
                // create a prototype for a background element
                luint tIJK[ 3 ] = { 0 };
                return new hmr::Background_Element< 3 >(
                        nullptr,
                        0,
                        tIJK,
                        0,
                        0,
                        0,
                        gNoProcID );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "unknown number of dimensions." );
            }
        }

        return aBackElement;
    }

    //------------------------------------------------------------------------------------
    void
    HMR_Helper::init_basis_index()
    {
        hmr::Background_Element_Base* tBackElement = this->create_background_element();

        // create a prototype of a B-Spline Element
        hmr::Element* tElement = mHMRLagrangeMesh->get_bspline_mesh( mBsplineMeshIndex )->create_element( tBackElement );

        // initialize index matrix
        mBasisIndex.set_size( mNumberOfBasis, 1 );

        // loop over all basis
        if ( mSpatialDimension == 2 )
        {
            // container for ijk position of basis
            luint tIJ[ 2 ];
            for ( uint k = 0; k < mNumberOfBasis; ++k )
            {
                // get position from element
                tElement->get_ijk_of_basis( k, tIJ );

                // calculate index in matrix
                uint tIndex = tIJ[ 0 ] + tIJ[ 1 ] * mNumberOfNodesPerDimension;

                mBasisIndex( tIndex ) = k;
            }
        }
        else
        {
            // container for ijk position of basis
            luint tIJK[ 3 ];
            for ( uint k = 0; k < mNumberOfBasis; ++k )
            {
                // get position from element
                tElement->get_ijk_of_basis( k, tIJK );

                // calculate index in matrix
                uint tIndex =
                        tIJK[ 0 ] + mNumberOfNodesPerDimension * ( tIJK[ 1 ] + tIJK[ 2 ] * mNumberOfNodesPerDimension );

                mBasisIndex( tIndex ) = k;
            }
        }

        // tidy up
        delete tElement;
        delete tBackElement;
    }

    //------------------------------------------------------------------------------------

    void
    HMR_Helper::get_extention_matrix_1d( real aShift, Matrix< DDRMat >& aExtentionMatrix )
    {
        switch ( mBSplineOrder )
        {
            case 1:
            {
                // fill out the matrix for a given data based on the shift
                aExtentionMatrix( 0, 0 ) = 1.0 - aShift;
                aExtentionMatrix( 0, 1 ) = aShift;
                aExtentionMatrix( 1, 0 ) = -aShift;
                aExtentionMatrix( 1, 1 ) = 1.0 + aShift;
                break;
            }
            case 2:
            {
                // fill out the matrix for a given data based on the shift
                aExtentionMatrix( 0, 0 ) = 0.5 * ( aShift - 1.0 ) * ( aShift - 2.0 );
                aExtentionMatrix( 0, 1 ) = aShift * ( -aShift + 2.0 );
                aExtentionMatrix( 0, 2 ) = 0.5 * aShift * ( aShift - 1.0 );
                aExtentionMatrix( 1, 0 ) = 0.5 * aShift * ( aShift - 1.0 );
                aExtentionMatrix( 1, 1 ) = -( aShift - 1.0 ) * ( aShift + 1.0 );
                aExtentionMatrix( 1, 2 ) = 0.5 * aShift * ( aShift + 1 );
                aExtentionMatrix( 2, 0 ) = 0.5 * aShift * ( aShift + 1.0 );
                aExtentionMatrix( 2, 1 ) = -aShift * ( aShift + 2.0 );
                aExtentionMatrix( 2, 2 ) = 0.5 * ( aShift + 1 ) * ( aShift + 2 );
                break;
            }

            case 3:
            {
                // fill out the matrix for a given data based on the shift
                aExtentionMatrix( 0, 0 ) = -1.0 / 6.0 * std::pow( aShift, 3.0 ) + 1.0 * std::pow( aShift, 2.0 ) - 11.0 / 6.0 * aShift + 1.0;    //
                aExtentionMatrix( 0, 1 ) = +0.5 * std::pow( aShift, 3.0 ) - 2.5 * std::pow( aShift, 2.0 ) + 3.0 * aShift + 0.0;                 //
                aExtentionMatrix( 0, 2 ) = -0.5 * std::pow( aShift, 3.0 ) + 2.0 * std::pow( aShift, 2.0 ) - 1.5 * aShift + 0.0;                 //
                aExtentionMatrix( 0, 3 ) = +1.0 / 6.0 * std::pow( aShift, 3.0 ) - 0.5 * std::pow( aShift, 2.0 ) + 1.0 / 3.0 * aShift + 0.0;     //

                aExtentionMatrix( 1, 0 ) = -1.0 / 6.0 * std::pow( aShift, 3.0 ) + 0.5 * std::pow( aShift, 2.0 ) - 1.0 / 3.0 * aShift + 0.0;     //
                aExtentionMatrix( 1, 1 ) = +0.5 * std::pow( aShift, 3.0 ) - 1.0 * std::pow( aShift, 2.0 ) - 0.5 * aShift + 1.0;                 //
                aExtentionMatrix( 1, 2 ) = -0.5 * std::pow( aShift, 3.0 ) + 0.5 * std::pow( aShift, 2.0 ) + 1.0 * aShift + 0.0;                 //
                aExtentionMatrix( 1, 3 ) = +1.0 / 6.0 * std::pow( aShift, 3.0 ) - 0.0 * std::pow( aShift, 2.0 ) - 1.0 / 6.0 * aShift + 0.0;

                aExtentionMatrix( 2, 0 ) = -1.0 / 6.0 * std::pow( aShift, 3.0 ) - 0.0 * std::pow( aShift, 2.0 ) + 1.0 / 6.0 * aShift + 0.0;    //
                aExtentionMatrix( 2, 1 ) = +0.5 * std::pow( aShift, 3.0 ) + 0.5 * std::pow( aShift, 2.0 ) - 1.0 * aShift + 0.0;                //
                aExtentionMatrix( 2, 2 ) = -0.5 * std::pow( aShift, 3.0 ) - 1.0 * std::pow( aShift, 2.0 ) + 0.5 * aShift + 1.0;                //
                aExtentionMatrix( 2, 3 ) = +1.0 / 6.0 * std::pow( aShift, 3.0 ) + 0.5 * std::pow( aShift, 2.0 ) + 1.0 / 3.0 * aShift + 0.0;

                aExtentionMatrix( 3, 0 ) = -1.0 / 6.0 * std::pow( aShift, 3.0 ) - 0.5 * std::pow( aShift, 2.0 ) - 1.0 / 3.0 * aShift + 0.0;    //
                aExtentionMatrix( 3, 1 ) = +0.5 * std::pow( aShift, 3.0 ) + 2.0 * std::pow( aShift, 2.0 ) + 1.5 * aShift + 0.0;                //
                aExtentionMatrix( 3, 2 ) = -0.5 * std::pow( aShift, 3.0 ) - 2.5 * std::pow( aShift, 2.0 ) - 3.0 * aShift + 0.0;                //
                aExtentionMatrix( 3, 3 ) = +1.0 / 6.0 * std::pow( aShift, 3.0 ) + 1.0 * std::pow( aShift, 2.0 ) + 11.0 / 6.0 * aShift + 1.0;

                break;
            }

            default:
            {
                MORIS_ASSERT( false, "HMR_Helper::get_extention_matrix_1d, The degree of the bpsline basis %u is not supported", mBSplineOrder );
                break;
            }
        }
    }

    //------------------------------------------------------------------------------------

    Matrix< DDRMat > const &
    HMR_Helper::get_l2_projection_matrix( const mtk::Cell* aExtendedCell,const  mtk::Cell* aRootCell )
    {
        // reset the L2 projection matrix and individual 1d matrices
        mL2ProjectionMatrix.set_size( mNumberOfBasis, mNumberOfBasis, 0.0 );

        // get the extended hmr element and the corresponding ijk
        const hmr::Element* tHMRExtenedElement = dynamic_cast< const hmr::Element* >( aExtendedCell );
        const hmr::Element* tHMRRootElement    = dynamic_cast< const hmr::Element* >( aRootCell );
        
        // obtain the ijk of those elements 
        const luint* tExtendedIJK = tHMRExtenedElement->get_ijk();
        const luint* tRootIJK = tHMRRootElement->get_ijk();

        // loop over the dimensions and create the i,j,k 1D matrices
        for ( uint iDim = 0; iDim < mSpatialDimension; iDim++ )
        {
            real tShift = tRootIJK[ iDim ] < tExtendedIJK[ iDim ] ? real( -tRootIJK[ iDim ] + tExtendedIJK[ iDim ] ) : -real( -tExtendedIJK[ iDim ] + tRootIJK[ iDim ] );
            
            //  find the shift in each direction
            this->get_extention_matrix_1d( tShift, mMatrices1D( iDim ) );
        }

        // Apply a tensor product to get the final weights
        if ( mSpatialDimension == 2 )
        {
            uint b = 0;
            for ( uint l = 0; l < mNumberOfNodesPerDimension; ++l )
            {
                for ( uint k = 0; k < mNumberOfNodesPerDimension; ++k )
                {
                    uint a = 0;
                    for ( uint j = 0; j < mNumberOfNodesPerDimension; ++j )
                    {
                        for ( uint i = 0; i < mNumberOfNodesPerDimension; ++i )
                        {
                            mL2ProjectionMatrix( mBasisIndex( a ), mBasisIndex( b ) ) = mMatrices1D( 0 )( i, k ) * mMatrices1D( 1 )( j, l );
                            ++a;
                        }
                    }
                    ++b;
                }
            }
        }
        else if ( mSpatialDimension == 3 )
        {
            uint b = 0;
            for ( uint p = 0; p < mNumberOfNodesPerDimension; ++p )
            {
                for ( uint q = 0; q < mNumberOfNodesPerDimension; ++q )
                {
                    for ( uint l = 0; l < mNumberOfNodesPerDimension; ++l )
                    {
                        uint a = 0;
                        for ( uint k = 0; k < mNumberOfNodesPerDimension; ++k )
                        {
                            for ( uint j = 0; j < mNumberOfNodesPerDimension; ++j )
                            {
                                for ( uint i = 0; i < mNumberOfNodesPerDimension; ++i )
                                {
                                    mL2ProjectionMatrix( mBasisIndex( a ), mBasisIndex( b ) ) = mMatrices1D( 0 )( i, l ) * mMatrices1D( 1 )( j, q ) * mMatrices1D( 2 )( k, p );
                                    ++a;
                                }
                            }
                        }
                        ++b;
                    }
                }
            }
        }

        // return the matrix
        return mL2ProjectionMatrix;
    }


    //------------------------------------------------------------------------------------
    moris::Cell< moris_id >&
    HMR_Helper::get_enriched_basis_indicies_of_cell( const mtk::Cell* aCell, moris_index aSPGIndex )
    {
        // reassign values in for the return cell
        mEnrichedBasisIDs.assign( mNumberOfBasis, gNoID );

        // get the hmr element
        const hmr::Element* tHMRElement = dynamic_cast< const hmr::Element* >( aCell );

        uint iIndexCounter = 0;
        // fill out the cell data
        for ( uint i = 0; i < mNumberOfBasis; i++ )
        {
            // get the basis
            const hmr::Basis* tBasis = tHMRElement->get_basis( i );

            // if it is active add it to the cell
            if ( tBasis->is_active() )
            {
                mEnrichedBasisIDs( iIndexCounter++ ) = tBasis->get_index();
            }
        }

        // get the BG basis anf their level
        moris::Cell< moris_index > const &             tBGBasisIndices = mEnrichmentData->mSubphaseGroupBGBasisIndices( aSPGIndex );
        moris::Cell< moris_index > const &             tBGBasisLevels  = mEnrichmentData->mSubphaseGroupBGBasisEnrLev( aSPGIndex );
        std::unordered_map< moris_index, moris_index > tBasisToLocalIndexMapRoot;

        // need to convert the first one to map
        convert_cell_to_map( tBGBasisIndices, tBasisToLocalIndexMapRoot );

        // The enriched basis indices based on the background and level
        moris::Cell< Matrix< IndexMat > > const & tBasisEnrichmentIndices = mEnrichmentData->mBasisEnrichmentIndices;

        // loop over the unenriched basis indices and convert them to enriched ones based on the background and level
        for ( uint iBasisOrd = 0; iBasisOrd < mNumberOfBasis; iBasisOrd++ )
        {
            // get the background basis in order
            moris_index tBGBasisIndex = mEnrichedBasisIDs( iBasisOrd );

            // get the basis index
            auto        tIter            = tBasisToLocalIndexMapRoot.find( tBGBasisIndex );
            moris_index tLocalBasisIndex = tIter->second;

            // Now use the local basis index to get the enrichment level
            moris_index tEnrLevel = tBGBasisLevels( tLocalBasisIndex );

            // get the enriched basis index of the root basis
            moris_index tEnrichedBasisIndex = tBasisEnrichmentIndices( tBGBasisIndex )( tEnrLevel );

            // switch the BG basis index with the enriched one
            mEnrichedBasisIDs( iBasisOrd ) = tEnrichedBasisIndex;
        }

        return mEnrichedBasisIDs;
    }

    //------------------------------------------------------------------------------------

    const luint*
    HMR_Helper::get_ijk_bspline_cell( const mtk::Cell* aCell )
    {   
        // cast the element into hmr element
        const hmr::Element* tHMRElement = dynamic_cast< const hmr::Element* >( aCell );

        // get the ijk of the element
        const luint* tIJK = tHMRElement->get_ijk();

        // return the ijk 
        return tIJK;
    }


}    // namespace xtk

/*
 * cl_Lagrange_Filter.cpp
 *
 *  Created on: Mar 7, 2018
 *      Author: gleim
 */

#include "cl_Lagrange_Filter.hpp" // STK/src/Hierarchical
using namespace moris;

void Lagrange_Filter::Update_IDandTMatrix_design(
        uint const &      aDim,
        uint const &      aPolynomialDesign,
        Mat<uint> const & aNumElements,
        uint const &      aLevel,
        real const &      aFilterRadius,
        BoostBitset &     aElementActiveDesign,
        const BoostBitset & aDesignBasisActive,
        const BoostBitset & aLagrangeBasisActive,
        const uint        & aNumberOfActiveDesignLagrangeBasis,
        //const map<uint, uint> & aDesignBSplineListMap,
        Mat<real> const & aDimensions,
        Mat<real> const & aDimensions_Offset,
        Mat<real> const & aDimensionsOriginal,
        Mat<uint> &       aIdFieldDesignNodalField,
        Mat<real> &       aTMatrixDesignNodalField,
        bool const &      aFilterLevelWeight,
        bool const &      aPerformNormalization,
        Mat<real> const & aPointOfOrigin,
        Mat<uint> &       IdFieldDesignNodalList)
{

    uint tMaxIDField = 0;

    Mat<real> tNewTMatrix;
    Mat<uint> tBasisFilterIDs;
    Mat<real> tBasisFilterAlpha;
    Mat<real> tBasisFilterBeta;
    Mat<uint> tNewIdField;
    Mat<uint> tFindUniqueList;

    Mat<uint> tBasisIDs;
    Mat<real> tBasisWeights;

    Mat<real> tHelp;
    Mat<real> tHelpa;
    Mat<real> tHelpb;
    Mat<uint> tHelpc;

    Cell<Mat<real>> tNewTMatrixCell(aIdFieldDesignNodalField.n_rows());
    Cell<Mat<uint>> tNewIdFieldCell(aIdFieldDesignNodalField.n_rows());

    // FIXME: should use same for parallel and serial
    // Serial uses a global numbering and saving for the BasisFilterCell (faster in serial), MPI uses a local numbering

    if( par_size() == 1)
    {
        // create map object
        map < uint, uint > tIdFieldMap;
        for(uint j = 0; j < IdFieldDesignNodalList.length(); j++)
        {
            tIdFieldMap[ IdFieldDesignNodalList( j ) ] = j;
        }

        for(uint i = 0; i< IdFieldDesignNodalList.length(); i++)
        {
            // Compute active design basis in radius
            filter_for_smoothing(
                    IdFieldDesignNodalList( i ),
                    aDim,
                    aPolynomialDesign,
                    aNumElements,
                    aLevel,
                    aFilterRadius,
                    aLagrangeBasisActive,
                    aNumberOfActiveDesignLagrangeBasis,
                    aDimensions,
                    aDimensions_Offset,
                    tBasisFilterIDs,
                    tBasisFilterAlpha,
                    tBasisFilterBeta
                    );

            // Check that there is at least one basis function is in filter list
            MORIS_ASSERT( tBasisFilterIDs.length() > 0, "No Lagrange basis found within filter radius");

            uint tNumFilterBasis = 0;
            for(uint j = 0; j< tBasisFilterIDs.length(); ++j)
            {
                tNumFilterBasis += (uint) aIdFieldDesignNodalField( tIdFieldMap[ tBasisFilterIDs( j ) ] );
            }

            // Normalize filter weights w/o accounting for level of basis functions
            if( aFilterLevelWeight == true)
            {
                tBasisFilterAlpha = tBasisFilterAlpha % tBasisFilterBeta;
            }
            real tWeight = sum( tBasisFilterAlpha );

            // Divide weights by sum of all weights for normalization
            tBasisFilterAlpha.col(0) /= tWeight;

            tBasisIDs.set_size( tNumFilterBasis, 1 );
            tBasisWeights.set_size( tNumFilterBasis, 1 );

            // Get the design basis functions from all Lagrange nodes, which are found within the filter radius
            tNumFilterBasis = 0;

            // Loop over all Lagrange nodes within filter radius
            for ( uint j = 0; j < tBasisFilterIDs.length(); j++ )
            {
                // ID of Lagrange node
                uint tLagrangeNodeID    = tBasisFilterIDs( j );

                // Index of Lagrange node
                uint tLagrangeNodeIndex = tIdFieldMap[ tLagrangeNodeID ];

                // Number of design basis of Lagrange node
                uint tNumDesignBasis    = aIdFieldDesignNodalField( tLagrangeNodeIndex, 0);

                // loop over all design basis of a Lagrange node
                for(uint k = 0; k < tNumDesignBasis; k++)
                {
                    // get ID
                    tBasisIDs ( tNumFilterBasis ) = aIdFieldDesignNodalField( tLagrangeNodeIndex, k + 1 );

                    // calculate weight
                    tBasisWeights( tNumFilterBasis ) = tBasisFilterAlpha( j ) * aTMatrixDesignNodalField( tLagrangeNodeIndex, k + 1 );

                    // increment counter
                    tNumFilterBasis++;
                }
            }

            // adapt size of tBasisIDs for unique finding
            tBasisIDs.resize( tNumFilterBasis, 1);

            // Determine unique design basis IDs
            tFindUniqueList = unique( tBasisIDs );

            // create map object
            map < uint, uint > tUniqueMap;
            for(uint j = 0; j < tFindUniqueList.length(); j++)
            {
                tUniqueMap[ tFindUniqueList( j ) ] = j+1;
            }

            // Set size for compressed Tmatrix and design basis ID vector
            tNewTMatrix.set_size(1,tFindUniqueList.length()+1,0);
            tNewIdField.set_size(1,tFindUniqueList.length()+1,0);

            // Build compressed Tmatrix and design basis ID vector
            for(uint j = 0; j < tNumFilterBasis; j++)
            {
                uint tCurrentBasisID         = tBasisIDs ( j );
                uint tUniqueMapID            = tUniqueMap[ tCurrentBasisID ];
                tNewIdField( tUniqueMapID )  = tCurrentBasisID;
                tNewTMatrix( tUniqueMapID ) += tBasisWeights( j );
            }

            // Normalize filter weights
            if( aPerformNormalization == true)
            {
                tWeight     = sum( tNewTMatrix );
                tNewTMatrix = tNewTMatrix / tWeight;
            }

            // Write number of unique design basis on compressed Tmatrix and design basis ID vectors
            tNewTMatrix(0) = tFindUniqueList.length();
            tNewIdField(0) = tFindUniqueList.length();

            // Save maximum number of unique design basis
            tMaxIDField = std::max( tMaxIDField, (uint)tNewIdField.length() );

            // Save compressed Tmatrix and design basis ID vectors for current Lagrange node
            tNewTMatrixCell(i) = tNewTMatrix;
            tNewIdFieldCell(i) = tNewIdField;
        }

        // std::cout << "Candidates ratio " << (real) gCandidateHitCounter / (real) gCandidateCounter << std::endl;
    }
    else
    {
        MORIS_ASSERT( false , "To be removed - parallel implementation not functional");
    }

    aIdFieldDesignNodalField.set_size(aIdFieldDesignNodalField.n_rows(),tMaxIDField,0);
    aTMatrixDesignNodalField.set_size(aIdFieldDesignNodalField.n_rows(),tMaxIDField,0);

    for(uint i=0; i < aIdFieldDesignNodalField.n_rows(); i++)
    {
        tHelpc =  tNewIdFieldCell(i).row(0);
        tHelpa =  tNewTMatrixCell(i).row(0);

        tHelpc.resize(1,tMaxIDField);
        tHelpa.resize(1,tMaxIDField);

        aIdFieldDesignNodalField.row(i) = tHelpc.row(0);
        aTMatrixDesignNodalField.row(i) = tHelpa.row(0);
    }
}

void
Lagrange_Filter::filter_for_smoothing(
        const uint               & aBasis,
        const uint               & aDim,
        const uint               & aPolynomialDesign,
        const Mat<uint>          & aNumElements,
        const uint               & aLevel,
        const real               & aFilterRadius,
        const BoostBitset        & aLagrangeBasisActive,
        const uint               & aNumberOfActiveDesignLagrangeBasis,
        const Mat<real>          & aDimensions,
        const Mat<real>          & aDimensions_Offset,
        Mat< uint >              & aFoundBasis,
        Mat< real >              & aAlpha,
        Mat< real >              & aBeta)
{
    // calculate ijk-position of aBasis
    Mat<uint> tBasisPosition   = mLagrangeBasis.give_position_of_basis(aBasis,aDim,aPolynomialDesign,aNumElements);

    // calculate coordinate of aBasis
    Mat<real> tBasisCoordinate = mLagrangeBasis.give_coordinate_from_basis(aBasis,aDim,aPolynomialDesign,aNumElements,aDimensions,aDimensions_Offset);

    // calculate level of aBasis
    uint taBasisLevel = mLagrangeBasis.give_basis_level(aBasis,aDim,aPolynomialDesign,aNumElements);

    // ---------------------------------------------------------------
    // Pre-selection algorithm

    // reserve memory for list of basis
    Mat< uint > tListOfBasisCandidates(aNumberOfActiveDesignLagrangeBasis, 1);

    // initialize counter
    uint tNumberOfBasis = 0;

    // loop over all levels
    for ( uint tLevel=0; tLevel<= aLevel; ++tLevel )
    {
        // calculate candidates for level and append them to tListOfBasisCandidates
        give_active_basis_on_level_within_bounding_box(
                aBasis,
                tLevel,
                aFilterRadius,
                aDimensions,
                aLagrangeBasisActive,
                aDim,
                aPolynomialDesign,
                aNumElements,
                tListOfBasisCandidates,
                tNumberOfBasis);
    }

    // adapt size of found candidates
    tListOfBasisCandidates.resize( tNumberOfBasis, 1);

    // add number of basis to global counter ( for debugging )
    gCandidateCounter += tNumberOfBasis;

    // ---------------------------------------------------------------

    aFoundBasis.set_size(tNumberOfBasis,1);
    aAlpha.set_size(tNumberOfBasis,1);
    aBeta.set_size(tNumberOfBasis,1);

    uint tCount = 0;
    real tDistance;
    uint tBasisLevel;
    Mat<real> tBasisCoordinates;

    //Search for active basis design functions and check the radius
    MORIS_ASSERT( aLagrangeBasisActive.size() > tListOfBasisCandidates.max(), "Wrong basis function is found, maybe outside the domain ?");

    // Loop over all active basis
    for( uint i = 0; i < tNumberOfBasis; ++i )
    {
        tBasisLevel       = mLagrangeBasis.give_basis_level( tListOfBasisCandidates( i ), aDim,aPolynomialDesign, aNumElements );
        tBasisCoordinates = mLagrangeBasis.give_coordinate_from_basis( tListOfBasisCandidates( i ) , aDim, aPolynomialDesign, aNumElements, aDimensions, aDimensions_Offset );

        tBasisCoordinates = tBasisCoordinate - tBasisCoordinates;
        tDistance         = tBasisCoordinates.norm();

        if( (aFilterRadius - tDistance) > 0.0 )
        {
            // copy ID of basis
            aFoundBasis( tCount ) = tListOfBasisCandidates( i );

            // calculate factor alpha
            aAlpha( tCount )      = aFilterRadius - tDistance;

            // calculate factor beta
            aBeta( tCount )       = ( (real) taBasisLevel + 1.0 )/( (real) tBasisLevel + 1.0 );

            // increment counter
            ++tCount;

            // add number of basis to global counter ( for debugging )
            ++gCandidateHitCounter;
        }
    }

    // reduce array to used size
    aFoundBasis.resize( tCount, 1 );
    aAlpha.resize( tCount, 1 );
    aBeta.resize( tCount, 1 );
}

void
Lagrange_Filter::give_active_basis_on_level_within_bounding_box(
        const uint        & aBasisID,
        const uint        & aLevelToLook,
        const real        & aRadius,
        const Mat<real>   & aDimensions,
        const BoostBitset & aLagrangeBasisActive,
        const uint        & aModelDim,
        const uint        & aPolynomial,
        const Mat<uint>   & aNumberOfElementsPerDirection,
        Mat< uint >       & aBasisWithinBoundingBox,
        uint              & aBasisCounter
)
{
    // assert proper input
    MORIS_ASSERT ( aModelDim > 1, "give_active_basis_on_level_within_bounding_box only works for dimensions 2 and 3");

    // minumum IJK position of interest
    Mat< uint > tMinIJK( aModelDim, 1);

    // maximum IJK position of interest
    Mat< uint > tMaxIJK( aModelDim, 1);

    calculate_min_and_max_ijk(
            aBasisID,
            aLevelToLook,
            aRadius,
            aDimensions,
            aModelDim,
            aPolynomial,
            aNumberOfElementsPerDirection,
            tMinIJK,
            tMaxIJK);

    // IJK position of basis to be investigated
    Mat< uint > tBasisIJK( aModelDim, 1);

    // scan for candidates
    if ( aModelDim == 2)
    {
        // 2D Case

        // loop along j-direction of bounding box
        for ( uint j=tMinIJK( 1 ); j<= tMaxIJK( 1 ); ++j )
        {
            // save second coordinate of basis
            tBasisIJK(1) = j;

            // loop along i-direction of bounding box
            for (uint i=tMinIJK( 0 ); i<= tMaxIJK( 0 ); ++i )
            {
                // save first coordinate of basis
                tBasisIJK( 0 ) = i;

                // get ID of candidate
                uint tCandidate = mLagrangeBasis.give_basis_of_position(
                        aLevelToLook,
                        aModelDim,
                        aPolynomial,
                        aNumberOfElementsPerDirection,
                        tBasisIJK );

                // add to list if basis is active
                if ( aLagrangeBasisActive.test( tCandidate ) )
                {
                    // add basis to list if active
                    aBasisWithinBoundingBox ( aBasisCounter ) = tCandidate;

                    // increment counter
                    ++aBasisCounter;

                }
            }
        }
    }
    else
    {
        // 3D Case

        // loop along k-direction of bounding box
        for ( uint k=tMinIJK( 2 ); k<= tMaxIJK(2); ++k )
        {
            // save third coordinate of basis
            tBasisIJK( 2 ) = k;

            // loop along j-direction of bounding box
            for ( uint j=tMinIJK( 1 ); j<= tMaxIJK(1); ++j )
            {
                // save second coordinate of basis
                tBasisIJK( 1 ) = j;

                // loop along i-direction of bounding box
                for (uint i=tMinIJK( 0 ); i<= tMaxIJK( 0 ); ++i)
                {
                    // save first coordinate of basis
                    tBasisIJK( 0 ) = i;

                    // get ID of candidate
                    uint tCandidate = mLagrangeBasis.give_basis_of_position(
                            aLevelToLook,
                            aModelDim,
                            aPolynomial,
                            aNumberOfElementsPerDirection,
                            tBasisIJK );

                    if ( aLagrangeBasisActive.test( tCandidate )  )
                    {
                        // add basis to list if active
                        aBasisWithinBoundingBox ( aBasisCounter ) = tCandidate;

                        // increment counter
                        ++aBasisCounter;
                    }
                }
            }
        }
    }
}

void
Lagrange_Filter::calculate_min_and_max_ijk(
        const uint        & aBasisID,
        const uint        & aLevelToLook,
        const real        & aRadius,
        const Mat<real>   & aDimensions,
        const uint        & aModelDim,
        const uint        & aPolynomial,
        const Mat<uint>   & aNumberOfElementsPerDirection,
        Mat < uint >      & aMinIJK,
        Mat < uint >      & aMaxIJK )
{

    // get level of basis
    uint tLevelOfBasis = mLagrangeBasis.give_basis_level(
            aBasisID,
            aModelDim,
            aPolynomial,
            aNumberOfElementsPerDirection );

    // level scale factor
    real tDeltaLevelScale = pow( 2, (real) aLevelToLook - (real)  tLevelOfBasis);

    // get position of basis
    Mat<uint> tBasisIJK   = mLagrangeBasis.give_position_of_basis(
            aBasisID,
            aModelDim,
            aPolynomial,
            aNumberOfElementsPerDirection );

    // translate position to level
    Mat < real > tCenterIJK( aModelDim, 1 );

    // calculate ijk position of center
    for ( uint k=0; k< aModelDim; ++k )
    {
        tCenterIJK( k ) = tDeltaLevelScale * tBasisIJK( k );
    }

    // level scale factor
    real tLevelScale = pow( 2, aLevelToLook );

    // IJK difference
    Mat < real > tDeltaIJK ( aModelDim, 1 );
    // loop over all dimensions
    for ( uint k=0; k< aModelDim; ++k )
    {
        // calculate half edge length of bounding box around aBasisID
        tDeltaIJK( k ) = aRadius * tLevelScale * aNumberOfElementsPerDirection( k ) / aDimensions( k ) ;
    }

    // minimum possible i-position
    uint tLowerI = tLevelScale*aPolynomial;

    // loop over all dimensions and calculate min and max IJK
    for ( uint k=0; k< aModelDim; ++k )
    {
        // minimum i-position of interest
        aMinIJK( k ) = ceil(( ( tCenterIJK( k ) )  < ( tLowerI  + tDeltaIJK( k ) ) ) ? ( tLowerI ) : ( tCenterIJK( k ) - tDeltaIJK( k ) ) );

        // maximum possible i-position
        uint tUpperI = aNumberOfElementsPerDirection( k ) * tLevelScale;

        // maximum i-position of interest
        aMaxIJK ( k ) = floor(( ( tCenterIJK( k ) + tDeltaIJK( k ) ) > tUpperI ) ? ( tUpperI ) : ( tCenterIJK( k ) + tDeltaIJK( k ) ) );
    }
}

/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Bilinear.cpp
 *
 */

#include "cl_GEN_Intersection_Node_Bilinear.hpp"
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Interpolation.hpp"

#include "cl_MTK_Interpolation_Function_Base.hpp"       //MTK/src
#include "cl_MTK_Interpolation_Function_Factory.hpp"    //MTK/src
#include "cl_MTK_Enums.hpp"                             //MTK/src

#include "fn_dot.hpp"
#include "fn_trans.hpp"

namespace moris
{
    namespace ge
    {
        //--------------------------------------------------------------------------------------------------------------

        Intersection_Node_Bilinear::Intersection_Node_Bilinear(
                std::shared_ptr< Intersection_Node > aFirstParentNode,
                std::shared_ptr< Intersection_Node > aSecondParentNode,
                uint                                 aFirstParentNodeIndex,
                uint                                 aSecondParentNodeIndex,
                const Matrix< DDRMat >&              aFirstParentNodeLocalCoordinates,
                const Matrix< DDRMat >&              aSecondParentNodeLocalCoordinates,
                const Matrix< DDUMat >&              aAncestorNodeIndices,
                const Cell< Matrix< DDRMat > >&      aAncestorNodeCoordinates,
                const Element_Intersection_Type      aInterpolationType,
                std::shared_ptr< Geometry >          aInterfaceGeometry,
                real                                 aIsocontourThreshold,
                real                                 aIsocontourTolerance,
                real                                 aIntersectionTolerance )
                : Intersection_Node(
                        get_local_coordinate(
                                aFirstParentNodeLocalCoordinates,
                                aSecondParentNodeLocalCoordinates,
                                aAncestorNodeIndices,
                                aAncestorNodeCoordinates,
                                aInterfaceGeometry,
                                aIsocontourThreshold ),
                        aFirstParentNode,
                        aSecondParentNode,
                        aFirstParentNodeIndex,
                        aSecondParentNodeIndex,
                        aFirstParentNodeLocalCoordinates,
                        aSecondParentNodeLocalCoordinates,
                        aAncestorNodeIndices,
                        aAncestorNodeCoordinates,
                        aInterpolationType,
                        aInterfaceGeometry,
                        aIsocontourThreshold,
                        aIsocontourTolerance,
                        aIntersectionTolerance,
                        false )
        {
            // initialize number of basis used by interpolation
            uint tNumBases;

            // build interpolator
            mtk::Interpolation_Function_Factory tFactory;

            mtk::Interpolation_Function_Base* tInterpolation;

            // create interpolation function based on spatial dimension  of problem
            switch ( aAncestorNodeCoordinates( 0 ).numel() )
            {
                case 2:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::QUAD,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    tNumBases = 4;
                    break;
                }
                case 3:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::HEX,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    tNumBases = 8;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false,
                            "Intersection_Node_Bilinear::Intersection_Node_Bilinear - Interpolation type not implemented." );
                }
            }

            // check that number of node indices and number of nodal coodinates of ancestor nodes are identical
            MORIS_ASSERT( aAncestorNodeCoordinates.size() == aAncestorNodeIndices.numel(),
                    "Intersection_Node_Bilinear::compute_intersection - inconsistent ancestor node information." );

            // check that number of bases to be used less or equal number of ancestor nodes
            MORIS_ASSERT( aAncestorNodeCoordinates.size() >= tNumBases,
                    "Intersection_Node_Bilinear::compute_intersection - number of ancestor nodes insufficient." );

            // allocate matrix for level set values at background cell nodes
            Matrix< DDRMat > tPhiBCNodes( tNumBases, 1 );

            // get level set values of corner nodes
            for ( uint in = 0; in < tNumBases; ++in )
            {
                tPhiBCNodes( in ) = aInterfaceGeometry->get_field_value( aAncestorNodeIndices( in ), aAncestorNodeCoordinates( in ) );
            }

            // check that dimension of ancestor node coordinate equals dimension of parent node coordinates
            MORIS_ASSERT( aFirstParentNodeLocalCoordinates.numel() == aAncestorNodeCoordinates( 0 ).numel(),
                    "Intersection_Node_Bilinear::compute_intersection - inconsistent coordinate dimensions." );

            // compute level set value at parent nodes
            Matrix< DDRMat > tFirstParentBasis;
            Matrix< DDRMat > tSecondParentBasis;

            tInterpolation->eval_N( aFirstParentNodeLocalCoordinates, tFirstParentBasis );
            tInterpolation->eval_N( aSecondParentNodeLocalCoordinates, tSecondParentBasis );

            real tFirstParentPhi  = dot( tFirstParentBasis, tPhiBCNodes );
            real tSecondParentPhi = dot( tSecondParentBasis, tPhiBCNodes );

            real tFirstDiffFromThreshold  = tFirstParentPhi - aIsocontourThreshold;
            real tSecondDiffFromThreshold = tSecondParentPhi - aIsocontourThreshold;

            // Overwrite basis values of intersection nodes
            Matrix< DDRMat > tCellCoordinate = 0.5 * ( 1.0 - mLocalCoordinate ) * aFirstParentNodeLocalCoordinates
                                             + 0.5 * ( 1.0 + mLocalCoordinate ) * aSecondParentNodeLocalCoordinates;

            tInterpolation->eval_N( tCellCoordinate, mBasisValues );

            // Global coordinates of intersection and parents
            mGlobalCoordinates = mBasisValues( 0 ) * aAncestorNodeCoordinates( 0 );

            Matrix< DDRMat > tFirstParentGlobalCoordinates  = tFirstParentBasis( 0 ) * aAncestorNodeCoordinates( 0 );
            Matrix< DDRMat > tSecondParentGlobalCoordinates = tSecondParentBasis( 0 ) * aAncestorNodeCoordinates( 0 );

            for ( uint tBasisIndex = 1; tBasisIndex < mBasisValues.length(); tBasisIndex++ )
            {
                mGlobalCoordinates += mBasisValues( tBasisIndex ) * aAncestorNodeCoordinates( tBasisIndex );
                tFirstParentGlobalCoordinates += tFirstParentBasis( tBasisIndex ) * aAncestorNodeCoordinates( tBasisIndex );
                tSecondParentGlobalCoordinates += tSecondParentBasis( tBasisIndex ) * aAncestorNodeCoordinates( tBasisIndex );
            }

            mParentVector      = trans( tSecondParentGlobalCoordinates - tFirstParentGlobalCoordinates );
            real tParentLength = norm( mParentVector );

            mFirstParentOnInterface = std::abs( tFirstDiffFromThreshold ) < aIntersectionTolerance
                                   or 0.5 * tParentLength * std::abs( 1 + mLocalCoordinate ) < aIntersectionTolerance;

            mSecondParentOnInterface = std::abs( tSecondDiffFromThreshold ) < aIntersectionTolerance
                                    or 0.5 * tParentLength * std::abs( 1 - mLocalCoordinate ) < aIntersectionTolerance;

            // Determine if edge is intersected
            if ( mFirstParentOnInterface or mSecondParentOnInterface )
            {
                mIsIntersected = true;
            }
            // FIXME: This check should be unnecessary as the local edge coordinate should be sufficient
            // to determine whether edge is intersected; it is only "useful" if parent node's level set value
            // is determined by method that is different from intersection nodes; for example level set value child node
            // of child node is computed via analytic field and intersection node via bi-linear interpolation
            else if ( tFirstDiffFromThreshold * tSecondDiffFromThreshold > 0 )
            {
                mIsIntersected = false;

                // check for consistency of parent values and local coordinate
                MORIS_ASSERT( std::abs( mLocalCoordinate ) > 1,
                        "Intersection_Node_Bilinear::Intersection_Node_Bilinear - inconsistent parent level set values versus local coordinate - p1 %e p2 %e loc %e.",
                        tFirstDiffFromThreshold,
                        tSecondDiffFromThreshold,
                        mLocalCoordinate );
            }
            else
            {
                mIsIntersected = ( std::abs( mLocalCoordinate ) <= 1.0 );

                // check for consistency with parent values
                // this check is currently useless but should be performed is inconsistency issue (see comment above) is resolved
                MORIS_ASSERT( mIsIntersected ? tFirstDiffFromThreshold * tSecondDiffFromThreshold < 0 : tFirstDiffFromThreshold * tSecondDiffFromThreshold > 0,
                        "Intersection_Node_Bilinear::Intersection_Node_Bilinear - inconsistent parent level set values - p1 %e p2 %e loc %e.",
                        tFirstDiffFromThreshold,
                        tSecondDiffFromThreshold,
                        mLocalCoordinate );
            }

            // allocate matrix for coordinates of parent nodes nodes on edge in local background cell CS
            mParentLocalCoordinates.set_size( aFirstParentNodeLocalCoordinates.numel(), 2 );

            // get coordinates of parent nodes nodes on edge in local background cell CS
            for ( uint iI = 0; iI < aFirstParentNodeLocalCoordinates.numel(); ++iI )
            {
                mParentLocalCoordinates( iI, 0 ) = aFirstParentNodeLocalCoordinates( iI );
                mParentLocalCoordinates( iI, 1 ) = aSecondParentNodeLocalCoordinates( iI );
            }

            // delete interpolator
            delete tInterpolation;
        }

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Node_Bilinear::~Intersection_Node_Bilinear()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Intersection_Node_Bilinear::get_dxi_dfield_from_ancestor( uint aAncestorIndex )
        {
            return this->compute_intersection_derivative( aAncestorIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Intersection_Node_Bilinear::get_dxi_dcoordinate_first_parent()
        {
            MORIS_ERROR( false, "Intersections on intersections not implemented yet for bilinear case." );
            return { {} };
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Intersection_Node_Bilinear::get_dxi_dcoordinate_second_parent()
        {
            MORIS_ERROR( false, "Intersections on intersections not implemented yet for bilinear case." );
            return { {} };
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Intersection_Node_Bilinear::get_local_coordinate(
                const Matrix< DDRMat >&         aFirstParentNodeLocalCoordinates,
                const Matrix< DDRMat >&         aSecondParentNodeLocalCoordinates,
                const Matrix< DDUMat >&         aAncestorNodeIndices,
                const Cell< Matrix< DDRMat > >& aAncestorNodeCoordinates,
                std::shared_ptr< Geometry >     aInterfaceGeometry,
                real                            aIsocontourThreshold )
        {
            // use Newton to compute intersection
            return this->compute_intersection(
                    aFirstParentNodeLocalCoordinates,
                    aSecondParentNodeLocalCoordinates,
                    aAncestorNodeIndices,
                    aAncestorNodeCoordinates,
                    aInterfaceGeometry,
                    aIsocontourThreshold );
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Intersection_Node_Bilinear::compute_intersection(
                const Matrix< DDRMat >&         aFirstParentNodeLocalCoordinates,
                const Matrix< DDRMat >&         aSecondParentNodeLocalCoordinates,
                const Matrix< DDUMat >&         aAncestorNodeIndices,
                const Cell< Matrix< DDRMat > >& aAncestorNodeCoordinates,
                std::shared_ptr< Geometry >     aInterfaceGeometry,
                real                            aIsocontourThreshold )
        {
            // spatial dimension
            uint tNumDims = aAncestorNodeCoordinates( 0 ).numel();

            // number of nodes to be used for interpolation
            uint tNumBases;

            // build interpolator
            mtk::Interpolation_Function_Factory tFactory;

            mtk::Interpolation_Function_Base* tInterpolation;

            // create interpolation function based on spatial dimension of problem
            switch ( tNumDims )
            {
                case 2:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::QUAD,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    tNumBases = 4;
                    break;
                }
                case 3:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::HEX,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    tNumBases = 8;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false,
                            "Intersection_Node_Bilinear::compute_intersection - Interpolation type not implemented." );
                }
            }

            // check that number of node indices and number of nodal coodinates of ancestor nodes are identical
            MORIS_ASSERT( aAncestorNodeCoordinates.size() == aAncestorNodeIndices.numel(),
                    "Intersection_Node_Bilinear::compute_intersection - inconsistent ancestor node information." );

            // check that number of bases to be used less or equal number of ancestor nodes
            MORIS_ASSERT( aAncestorNodeCoordinates.size() >= tNumBases,
                    "Intersection_Node_Bilinear::compute_intersection - number of ancestor nodes insufficient." );

            // allocate matrix for level set values at background cell nodes
            Matrix< DDRMat > tPhiBCNodes( tNumBases, 1 );

            // get level set values of corner nodes
            for ( uint in = 0; in < tNumBases; ++in )
            {
                tPhiBCNodes( in ) = aInterfaceGeometry->get_field_value(
                        aAncestorNodeIndices( in ),
                        aAncestorNodeCoordinates( in ) );
            }

            // scale element level set field such that norm equals 1.0
            const real tPhiScaling = 1.0 / norm( tPhiBCNodes );

            tPhiBCNodes = tPhiScaling * tPhiBCNodes;

            // scale threshold
            aIsocontourThreshold *= tPhiScaling;

            // check that dimension of ancestor node coordinate equals dimension of parent node coordinates
            MORIS_ASSERT( aFirstParentNodeLocalCoordinates.numel() == aAncestorNodeCoordinates( 0 ).numel(),
                    "Intersection_Node_Bilinear::compute_intersection - inconsistent coordinate dimensions." );

            // compute level set value at parent nodes
            Matrix< DDRMat > tFirstParentBasis;
            Matrix< DDRMat > tSecondParentBasis;

            tInterpolation->eval_N( aFirstParentNodeLocalCoordinates, tFirstParentBasis );
            tInterpolation->eval_N( aSecondParentNodeLocalCoordinates, tSecondParentBasis );

            const real tFirstParentPhi  = dot( tFirstParentBasis, tPhiBCNodes );
            const real tSecondParentPhi = dot( tSecondParentBasis, tPhiBCNodes );

            // check that line is intersected
            if ( ( tFirstParentPhi - aIsocontourThreshold ) * ( tSecondParentPhi - aIsocontourThreshold ) > 0 )
            {
                delete tInterpolation;

                return MORIS_REAL_MAX;
            }

            // check that difference between first and second parent is not smaller than MORIS_REAL_EPS
            MORIS_ERROR( std::abs( tFirstParentPhi - tSecondParentPhi ) > MORIS_REAL_EPS,
                    "Intersection_Node_Bilinear::compute_intersection - level set values of edge end points identical" );

            // set Newton parameters
            const uint tNewMaxIter  = 20;    // maximum number of iterations in Newton
            const uint tCurvMaxIter = 10;    // maximum number of iterations using curvature information

            const real tNewRelax  = 1.0;      // relaxation factor for solution update
            const real tNewRelEps = 1e-8;     // required relative residual drop
            const real tNewAbsEps = 1e-12;    // required absolute residual drop
            const real tCurvMin   = 1e-8;     // threshold for curvature magnitude above which curvature is used

            // allocate matrices used within Newton loop
            Matrix< DDRMat > tCellCoordinate( aFirstParentNodeLocalCoordinates.numel(), 1 );
            Matrix< DDRMat > tBasis;
            Matrix< DDRMat > tDBasisDxi;
            Matrix< DDRMat > tDBasisDxi2;
            Matrix< DDRMat > tD2PhiDxi2;

            // vector from first to second parent in local coordinates
            Matrix< DDRMat > tSecondToFirstParent = aSecondParentNodeLocalCoordinates - aFirstParentNodeLocalCoordinates;

            // initialized reference residual
            real tReferenceResidual = 0.0;
            real tResidual          = 0.0;

            // compute initial guess: location of intersection point along edge in edge CS
            real tEdgeCoordinate = ( 2.0 * aIsocontourThreshold - tFirstParentPhi - tSecondParentPhi )    //
                                 / ( tSecondParentPhi - tFirstParentPhi );

            Matrix< DDRMat > tInitialGuess = { { std::min( 1.0, std::max( tEdgeCoordinate, -1.0 ) ), -1.0, 1.0 } };

            // loop over initial guess trials
            for ( uint iGuess = 0; iGuess < tInitialGuess.numel(); iGuess++ )
            {
                // set initial guess
                tEdgeCoordinate = tInitialGuess( iGuess );

                // perform iterations
                for ( uint iNew = 0; iNew < tNewMaxIter; ++iNew )
                {
                    // compute local coordinate in background cell CS
                    tCellCoordinate = 0.5 * ( 1.0 - tEdgeCoordinate ) * aFirstParentNodeLocalCoordinates
                                    + 0.5 * ( 1.0 + tEdgeCoordinate ) * aSecondParentNodeLocalCoordinates;

                    // compute basis function
                    tInterpolation->eval_N( tCellCoordinate, tBasis );

                    // compute residual
                    tResidual = dot( tBasis, tPhiBCNodes ) - aIsocontourThreshold;

                    // check convergence against absolute residual
                    if ( std::abs( tResidual ) < tNewAbsEps )
                    {
                        delete tInterpolation;

                        return tEdgeCoordinate;
                    }

                    // store reference residual
                    if ( iNew == 0 )
                    {
                        tReferenceResidual = std::abs( tResidual );
                    }

                    // check convergence against relative residual
                    if ( std::abs( tResidual ) < tNewRelEps * tReferenceResidual )
                    {
                        delete tInterpolation;

                        return tEdgeCoordinate;
                    }

                    // compute Jacobian
                    tInterpolation->eval_dNdXi( tCellCoordinate, tDBasisDxi );

                    // compute first order gradient of residual with respect to edge coordinate
                    real tGradRes = 0.5 * dot( tSecondToFirstParent, tDBasisDxi * tPhiBCNodes );

                    // initialize solution increment
                    real tSolIncrement;

                    // initialize curvature variable
                    real tSqrt2   = -1.0;
                    real tCurvRes = 0.0;

                    // compute second order derivatives
                    if ( iNew < tCurvMaxIter )
                    {
                        // compute Hessian
                        tInterpolation->eval_d2NdXi2( tCellCoordinate, tDBasisDxi2 );

                        tD2PhiDxi2 = tDBasisDxi2 * tPhiBCNodes;

                        // compute second order derivative of residual with respect to edge coordinate
                        if ( tNumDims == 2 )
                        {
                            tCurvRes =                                                                                   //
                                    0.125 * tD2PhiDxi2( 0 ) * tSecondToFirstParent( 0 ) * tSecondToFirstParent( 0 ) +    //
                                    0.125 * tD2PhiDxi2( 1 ) * tSecondToFirstParent( 1 ) * tSecondToFirstParent( 1 ) +    //
                                    0.250 * tD2PhiDxi2( 2 ) * tSecondToFirstParent( 0 ) * tSecondToFirstParent( 1 );
                        }
                        else
                        {
                            tCurvRes =                                                                                   //
                                    0.125 * tD2PhiDxi2( 0 ) * tSecondToFirstParent( 0 ) * tSecondToFirstParent( 0 ) +    //
                                    0.125 * tD2PhiDxi2( 1 ) * tSecondToFirstParent( 1 ) * tSecondToFirstParent( 1 ) +    //
                                    0.125 * tD2PhiDxi2( 2 ) * tSecondToFirstParent( 2 ) * tSecondToFirstParent( 2 ) +    //
                                    0.250 * tD2PhiDxi2( 3 ) * tSecondToFirstParent( 1 ) * tSecondToFirstParent( 2 ) +    //
                                    0.250 * tD2PhiDxi2( 4 ) * tSecondToFirstParent( 0 ) * tSecondToFirstParent( 2 ) +    //
                                    0.250 * tD2PhiDxi2( 5 ) * tSecondToFirstParent( 0 ) * tSecondToFirstParent( 1 );
                        }

                        // compute square of sqrt term in root finding formula for quadratic equations
                        tSqrt2 = tGradRes * tGradRes - 4.0 * tCurvRes * tResidual;
                    }

                    // update solution depending on curvature
                    if ( iNew < tCurvMaxIter && std::abs( tCurvRes ) > tCurvMin && tSqrt2 >= 0.0 )
                    {
                        // add roots of quadratic function to current edge coordinate
                        real tSqrt = std::sqrt( tSqrt2 );

                        real tRootMin = ( -1.0 * tGradRes - tSqrt ) / ( 2.0 * tCurvRes );
                        real tRootMax = ( -1.0 * tGradRes + tSqrt ) / ( 2.0 * tCurvRes );

                        // check with root is in [-1,1] interval
                        bool tRootMinIsValid = std::abs( tEdgeCoordinate + tRootMin ) <= 1.0;
                        bool tRootMaxIsValid = std::abs( tEdgeCoordinate + tRootMax ) <= 1.0;

                        // compute increment
                        if ( tRootMinIsValid && tRootMaxIsValid )
                        {
                            tSolIncrement = std::abs( tRootMin ) < std::abs( tRootMax ) ? tRootMin : tRootMax;
                        }
                        else if ( tRootMinIsValid )
                        {
                            tSolIncrement = tRootMin;
                        }
                        else if ( tRootMaxIsValid )
                        {
                            tSolIncrement = tRootMax;
                        }
                        else
                        {
                            tSolIncrement = -1.0 * tResidual / ( tGradRes + MORIS_REAL_EPS );
                        }
                    }
                    else
                    {
                        // compute increment
                        tSolIncrement = -1.0 * tResidual / ( tGradRes + MORIS_REAL_EPS );
                    }

                    // update solution
                    tEdgeCoordinate += tNewRelax * tSolIncrement;

                    // trim solution
                    tEdgeCoordinate = std::min( 1.0, std::max( tEdgeCoordinate, -1.0 ) );
                }
            }

            // print debug information
            for ( uint in = 0; in < tNumBases; ++in )
            {
                std::string tStrg = "Anchestor_Node_" + std::to_string( aAncestorNodeIndices( in ) );
                print( aAncestorNodeCoordinates( in ), tStrg );
            }

            print( tPhiBCNodes, "tPhiBCNodes" );

            fprintf( stderr, "tFirstParentPhi =%e   tSecondParentPhi = %e\n", tFirstParentPhi, tSecondParentPhi );

            print( aFirstParentNodeLocalCoordinates, "aFirstParentNodeLocalCoordinates" );
            print( aSecondParentNodeLocalCoordinates, "aSecondParentNodeLocalCoordinates" );

            MORIS_ERROR( false,
                    "Intersection_Node_Bilinear::compute_intersection - Newton did not convergence: %s %e %s %e %s %e",
                    "Reference residual",
                    tReferenceResidual,
                    "Current residual",
                    std::abs( tResidual ),
                    "Edge coordinate",
                    tEdgeCoordinate );

            delete tInterpolation;

            return tEdgeCoordinate;
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Intersection_Node_Bilinear::compute_intersection_derivative( uint aAncestorIndex )
        {
            // Locked interface geometry
            std::shared_ptr< Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // number of nodes to be used for interpolation
            uint tNumBases;

            // build interpolator
            mtk::Interpolation_Function_Factory tFactory;

            mtk::Interpolation_Function_Base* tInterpolation;

            // create interpolation function based on spatial dimension  of problem
            switch ( mAncestorNodeCoordinates( 0 ).numel() )
            {
                case 2:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::QUAD,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    tNumBases = 4;
                    break;
                }
                case 3:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::HEX,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    tNumBases = 8;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false,
                            "Intersection_Node_Bilinear::compute_intersection_derivative - Interpolation type not implemented." );
                }
            }

            // allocate matrix for level set values at background cell nodes
            Matrix< DDRMat > tPhiBCNodes( tNumBases, 1 );

            // get level set values of corner nodes
            for ( uint in = 0; in < tNumBases; ++in )
            {
                tPhiBCNodes( in ) = tLockedInterfaceGeometry->get_field_value( mAncestorNodeIndices( in ), mAncestorNodeCoordinates( in ) );
            }

            // compute level set value at parent nodes
            Matrix< DDRMat > aFirstParentNodeLocalCoordinates  = mParentLocalCoordinates.get_column( 0 );
            Matrix< DDRMat > aSecondParentNodeLocalCoordinates = mParentLocalCoordinates.get_column( 1 );

            // compute local coordinate in background cell CS
            Matrix< DDRMat > tCellCoordinate = 0.5 * ( 1.0 - mLocalCoordinate ) * aFirstParentNodeLocalCoordinates
                                             + 0.5 * ( 1.0 + mLocalCoordinate ) * aSecondParentNodeLocalCoordinates;

            // compute basis function
            Matrix< DDRMat > tBasis;
            tInterpolation->eval_N( tCellCoordinate, tBasis );

            // compute derivative of residual
            real tDResidualDPhi = tBasis( aAncestorIndex );

            // compute Jacobian
            Matrix< DDRMat > tDBasisDxi;
            tInterpolation->eval_dNdXi( tCellCoordinate, tDBasisDxi );

            Matrix< DDRMat > tJac = 0.5 * trans( aSecondParentNodeLocalCoordinates - aFirstParentNodeLocalCoordinates ) * tDBasisDxi * tPhiBCNodes;

            // delete interpolator
            delete tInterpolation;

            // compute derivative of edge coordinate wrt. ancestor level set value
            return -1.0 * tDResidualDPhi / ( tJac( 0 ) + MORIS_REAL_EPS );
        }

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris


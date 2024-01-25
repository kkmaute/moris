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
#include "cl_GEN_Parent_Node.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Interpolation.hpp"

#include "cl_MTK_Interpolation_Function_Base.hpp"
#include "cl_MTK_Interpolation_Function_Factory.hpp"
#include "cl_MTK_Enums.hpp"

#include "fn_dot.hpp"
#include "fn_norm.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node_Bilinear::Intersection_Node_Bilinear(
            uint                     aNodeIndex,
            const Cell< Node* >&     aBackgroundNodes,
            const Parent_Node&       aFirstParentNode,
            const Parent_Node&       aSecondParentNode,
            mtk::Geometry_Type       aBackgroundGeometryType,
            mtk::Interpolation_Order aBackgroundInterpolationOrder,
            Level_Set_Geometry&      aInterfaceGeometry )
            : Intersection_Node_Level_Set(
                    aNodeIndex,
                    aBackgroundNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    Intersection_Node_Bilinear::compute_local_coordinate( aBackgroundNodes, aFirstParentNode, aSecondParentNode, aInterfaceGeometry ),
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder,
                    aInterfaceGeometry )
            , mParametricParentVector( aSecondParentNode.get_parametric_coordinates() - aFirstParentNode.get_parametric_coordinates() )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const Cell< Basis_Node >& Intersection_Node_Bilinear::get_field_basis_nodes() const
    {
        return this->get_background_nodes();
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Intersection_Node_Bilinear::get_dxi_dfield_from_ancestor( uint aAncestorIndex ) const
    {
        // number of nodes to be used for interpolation
        uint tNumBases;

        // determine number of basis for bi or tri-linear interpolation based on spatial dimensions
        switch ( this->get_global_coordinates().length() )
        {
            case 2:
            {
                tNumBases = 4;
                break;
            }
            case 3:
            {
                tNumBases = 8;
                break;
            }
            default:
            {
                MORIS_ERROR( false,
                        "Intersection_Node_Bilinear::compute_intersection_derivative - Improper spatial dimensions." );
            }
        }

        // check that aAncestorIndex <= number of basis
        // note: here only bi and tri-linear interpolation used irrespective of interpolation of background cell; thus only
        // corner nodes values are used; level set values of other nodes do not influence intersection position
        if ( aAncestorIndex >= tNumBases )
        {
            // return zero as only corner node level set values influence intersection
            return 0.0;
        }

        // build interpolator
        mtk::Interpolation_Function_Factory tFactory;
        mtk::Interpolation_Function_Base* tInterpolation;

        // create interpolation function based on spatial dimension  of problem
        switch ( tNumBases )
        {
            case 4:
            {
                tInterpolation = tFactory.create_interpolation_function(
                        mtk::Geometry_Type::QUAD,
                        mtk::Interpolation_Type::LAGRANGE,
                        mtk::Interpolation_Order::LINEAR );
                break;
            }
            case 8:
            {
                tInterpolation = tFactory.create_interpolation_function(
                        mtk::Geometry_Type::HEX,
                        mtk::Interpolation_Type::LAGRANGE,
                        mtk::Interpolation_Order::LINEAR );
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
        const Cell< Basis_Node >& tBackgroundNodes = this->get_background_nodes();
        for ( uint iBackgroundNodeIndex = 0; iBackgroundNodeIndex < tNumBases; ++iBackgroundNodeIndex )
        {
            tPhiBCNodes( iBackgroundNodeIndex ) = mInterfaceGeometry.get_field_value(
                    tBackgroundNodes( iBackgroundNodeIndex ).get_index(),
                    tBackgroundNodes( iBackgroundNodeIndex ).get_global_coordinates() );
        }

        // compute local coordinate in background cell CS
        const Matrix< DDRMat >& tCellCoordinate = this->get_parametric_coordinates();

        // compute basis function
        Matrix< DDRMat > tBasis;
        tInterpolation->eval_N( tCellCoordinate, tBasis );

        // compute derivative of residual
        real tDResidualDPhi = tBasis( aAncestorIndex );

        // compute Jacobian
        Matrix< DDRMat > tDBasisDxi;
        tInterpolation->eval_dNdXi( tCellCoordinate, tDBasisDxi );
        Matrix< DDRMat > tJac = 0.5 * mParametricParentVector * tDBasisDxi * tPhiBCNodes;

        // delete interpolator
        delete tInterpolation;

        // compute derivative of edge coordinate wrt. ancestor level set value
        return -1.0 * tDResidualDPhi / ( tJac( 0 ) + MORIS_REAL_EPS );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Bilinear::get_dxi_dcoordinate_first_parent() const
    {
        MORIS_ERROR( false, "Intersections on intersections not implemented yet for bilinear case." );
        return { {} };
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Bilinear::get_dxi_dcoordinate_second_parent() const
    {
        MORIS_ERROR( false, "Intersections on intersections not implemented yet for bilinear case." );
        return { {} };
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Intersection_Node_Bilinear::compute_local_coordinate(
            const Cell< Node* >&      aBackgroundNodes,
            const Parent_Node&        aFirstParentNode,
            const Parent_Node&        aSecondParentNode,
            const Level_Set_Geometry& aInterfaceGeometry )
    {
        // get isocontour threshold from geometry
        real tIsocontourThreshold = aInterfaceGeometry.get_isocontour_threshold();

        // spatial dimension
        uint tNumDims = aFirstParentNode.get_global_coordinates().length();

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

        // allocate matrix for level set values at background cell nodes
        Matrix< DDRMat > tPhiBCNodes( tNumBases, 1 );

        // get level set values of corner nodes
        for ( uint iBackgroundNode = 0; iBackgroundNode < tNumBases; ++iBackgroundNode )
        {
            tPhiBCNodes( iBackgroundNode ) = aInterfaceGeometry.get_field_value(
                    aBackgroundNodes( iBackgroundNode )->get_index(),
                    aBackgroundNodes( iBackgroundNode )->get_global_coordinates() );
        }

        // Scale element level set field such that norm equals 1.0
        real tPhiScaling = 1.0 / norm( tPhiBCNodes );
        tPhiBCNodes = tPhiScaling * tPhiBCNodes;
        tIsocontourThreshold *= tPhiScaling;

        // Get scaled parent level set values
        real tFirstParentPhi  = tPhiScaling * aInterfaceGeometry.get_field_value( aFirstParentNode.get_index(), aFirstParentNode.get_global_coordinates() );
        real tSecondParentPhi = tPhiScaling * aInterfaceGeometry.get_field_value( aSecondParentNode.get_index(), aSecondParentNode.get_global_coordinates() );

        // check that line is intersected
        if ( ( tFirstParentPhi - tIsocontourThreshold ) * ( tSecondParentPhi - tIsocontourThreshold ) > 0 )
        {
            delete tInterpolation;

            return MORIS_REAL_MAX;
        }

        // check that difference between first and second parent is not smaller than MORIS_REAL_EPS
        if ( std::abs( tFirstParentPhi - tSecondParentPhi ) < MORIS_REAL_EPS )
        {
            delete tInterpolation;

            // return that intersection node is at center of edge
            return 0.0;
        }

        // set Newton parameters
        const uint tNewMaxIter  = 20;    // maximum number of iterations in Newton
        const uint tCurvMaxIter = 10;    // maximum number of iterations using curvature information

        const real tNewRelax  = 1.0;      // relaxation factor for solution update
        const real tNewRelEps = 1e-8;     // required relative residual drop
        const real tNewAbsEps = 1e-12;    // required absolute residual drop
        const real tCurvMin   = 1e-8;     // threshold for curvature magnitude above which curvature is used

        // allocate matrices used within Newton loop
        Matrix< DDRMat > tCellCoordinate( tNumDims, 1 );
        Matrix< DDRMat > tBasis;
        Matrix< DDRMat > tDBasisDxi;
        Matrix< DDRMat > tDBasisDxi2;
        Matrix< DDRMat > tD2PhiDxi2;

        // vector from first to second parent in parametric coordinates
        Matrix< DDRMat > tSecondToFirstParent = aSecondParentNode.get_parametric_coordinates() - aFirstParentNode.get_parametric_coordinates();

        // initialized reference residual
        real tReferenceResidual = 0.0;
        real tResidual          = 0.0;

        // compute initial guess: location of intersection point along edge in edge CS
        real tEdgeCoordinate = ( 2.0 * tIsocontourThreshold - tFirstParentPhi - tSecondParentPhi )
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
                tCellCoordinate = 0.5 * ( 1.0 - tEdgeCoordinate ) * aFirstParentNode.get_parametric_coordinates()
                                + 0.5 * ( 1.0 + tEdgeCoordinate ) * aSecondParentNode.get_parametric_coordinates();

                // compute basis function
                tInterpolation->eval_N( tCellCoordinate, tBasis );

                // compute residual
                tResidual = dot( tBasis, tPhiBCNodes ) - tIsocontourThreshold;

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
        // for ( uint in = 0; in < tNumBases; ++in )
        // {
        //     std::string tStrg = "Anchestor_Node_" + std::to_string( aAncestorNodeIndices( in ) );
        //     print( aAncestorNodeCoordinates( in ), tStrg );
        // }

        // print( tPhiBCNodes, "tPhiBCNodes" );

        // fprintf( stderr, "tFirstParentPhi =%e   tSecondParentPhi = %e\n", tFirstParentPhi, tSecondParentPhi );

        // print( aFirstParentNodeLocalCoordinates, "aFirstParentNodeLocalCoordinates" );
        // print( aSecondParentNodeLocalCoordinates, "aSecondParentNodeLocalCoordinates" );

        MORIS_ERROR( false,
                "Intersection_Node_Bilinear::compute_intersection - Newton did not converge: %s %e %s %e %s %e",
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

}

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

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node_Bilinear::Intersection_Node_Bilinear(
            uint                     aNodeIndex,
            const Vector< Background_Node* >& aBackgroundNodes,
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
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder,
                    aInterfaceGeometry )
            , mParametricParentVector( aSecondParentNode.get_parametric_coordinates() - aFirstParentNode.get_parametric_coordinates() )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< Basis_Node >& Intersection_Node_Bilinear::get_field_basis_nodes() const
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
        const Vector< Basis_Node >& tBackgroundNodes = this->get_background_nodes();
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

}

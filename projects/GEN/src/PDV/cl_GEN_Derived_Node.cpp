/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Derived_Node.cpp
 *
 */

#include "cl_GEN_Derived_Node.hpp"
#include "cl_GEN_Basis_Node.hpp"
#include "cl_MTK_Interpolation_Function_Factory.hpp"
#include "cl_MTK_Interpolation_Function.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Derived_Node::Derived_Node(
            uint                     aIndex,
            const Cell< Node* >&     aBaseNodes,
            const Matrix< DDRMat >&  aParametricCoordinates,
            mtk::Geometry_Type       aGeometryType,
            mtk::Interpolation_Order aInterpolationOrder )
            : Node( aIndex )
            , mParametricCoordinates( aParametricCoordinates )
    {
        // Check that at least one base node was given
        MORIS_ASSERT( aBaseNodes.size() > 0, "A derived GEN node must have at least one basis node." );

        // Create interpolator
        mtk::Interpolation_Function_Factory tInterpolationFactory;
        mtk::Interpolation_Function_Base* tInterpolation = tInterpolationFactory.create_interpolation_function(
                aGeometryType,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // Perform interpolation using parametric coordinates
        Matrix< DDRMat > tBasis;
        tInterpolation->eval_N( mParametricCoordinates, tBasis );

        // Get number of bases
        uint tNumberOfBases = tInterpolation->get_number_of_bases();

        // Create background nodes
        mBackgroundNodes.reserve( tNumberOfBases );
        for ( uint iBasisIndex = 0; iBasisIndex < tNumberOfBases; iBasisIndex++ )
        {
            mBackgroundNodes.emplace_back( aBaseNodes( iBasisIndex ), tBasis( iBasisIndex ) );
        }

        // Size global coordinates based on first locator
        mGlobalCoordinates = Matrix< DDRMat >( 1, tInterpolation->get_number_of_dimensions(), 0.0 );
        delete tInterpolation;

        // Add contributions from all locators
        for ( auto iBasisNode : mBackgroundNodes )
        {
            mGlobalCoordinates += iBasisNode.get_global_coordinates() * iBasisNode.get_basis();
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Derived_Node::get_global_coordinates()
    {
        return mGlobalCoordinates;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Derived_Node::get_parametric_coordinates()
    {
        return mParametricCoordinates;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Cell< Basis_Node >& Derived_Node::get_background_nodes()
    {
        return mBackgroundNodes;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Cell< Basis_Node >& Derived_Node::get_locator_nodes()
    {
        return mBackgroundNodes;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Derived_Node::is_on_interface( Geometry* aGeometry )
    {
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

}

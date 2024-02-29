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
#include "cl_GEN_Background_Node.hpp"
#include "cl_GEN_Basis_Node.hpp"
#include "cl_MTK_Interpolation_Function_Factory.hpp"
#include "cl_MTK_Interpolation_Function.hpp"
#include "cl_Communication_Tools.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Derived_Node::Derived_Node(
            uint                     aIndex,
            const Vector< Background_Node* >& aBackgroundNodes,
            const Matrix< DDRMat >&  aParametricCoordinates,
            mtk::Geometry_Type       aGeometryType,
            mtk::Interpolation_Order aInterpolationOrder )
            : Node( aIndex )
            , mParametricCoordinates( aParametricCoordinates )
    {
        // Check that at least one background node was given
        MORIS_ASSERT( aBackgroundNodes.size() > 0, "A derived GEN node must have at least one basis node." );

        // Override linear interpolation if desired
        if ( gOverrideLinearInterpolation )
        {
            aInterpolationOrder = mtk::Interpolation_Order::LINEAR;
        }

        // Create interpolator
        mtk::Interpolation_Function_Factory tInterpolationFactory;
        mtk::Interpolation_Function_Base* tInterpolation = tInterpolationFactory.create_interpolation_function(
                aGeometryType,
                mtk::Interpolation_Type::LAGRANGE,
                aInterpolationOrder );

        // Perform interpolation using parametric coordinates
        Matrix< DDRMat > tBasis;
        tInterpolation->eval_N( mParametricCoordinates, tBasis );

        // Get number of bases
        uint tNumberOfBases = tInterpolation->get_number_of_bases();

        // Create background nodes
        mBackgroundNodes.reserve( tNumberOfBases );
        for ( uint iBasisIndex = 0; iBasisIndex < tNumberOfBases; iBasisIndex++ )
        {
            mBackgroundNodes.emplace_back( *aBackgroundNodes( iBasisIndex ), tBasis( iBasisIndex ) );
        }

        // Size global coordinates based on first locator
        mGlobalCoordinates = Matrix< DDRMat >( 1, mBackgroundNodes( 0 ).get_global_coordinates().length(), 0.0 );
        delete tInterpolation;

        // Add contributions from all locators
        for ( auto iBasisNode : mBackgroundNodes )
        {
            mGlobalCoordinates += iBasisNode.get_global_coordinates() * iBasisNode.get_basis();
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Derived_Node::get_global_coordinates() const
    {
        return mGlobalCoordinates;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Derived_Node::get_parametric_coordinates() const
    {
        return mParametricCoordinates;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< Basis_Node >& Derived_Node::get_background_nodes() const
    {
        return mBackgroundNodes;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< Basis_Node >& Derived_Node::get_locator_nodes() const
    {
        return mBackgroundNodes;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Derived_Node::is_on_interface( const Geometry& aGeometry ) const
    {
        return false;
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    real Derived_Node::get_coordinate_value( uint aCoordinateIndex ) const
    {
        return mGlobalCoordinates( aCoordinateIndex );
    }
        
    //--------------------------------------------------------------------------------------------------------------
    
    uint Derived_Node::get_num_pdvs()
    {
        return 0;
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    void Derived_Node::set_starting_pdv_id( moris_id aPDVStartingID )
    {
    }

    //--------------------------------------------------------------------------------------------------------------
    
    moris_id Derived_Node::get_starting_pdv_id()
    {
        return -1;
    }

    //--------------------------------------------------------------------------------------------------------------
    
    void Derived_Node::set_id( moris_id aNodeID )
    {
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    void Derived_Node::set_owner( moris_index aNodeOwner )
    {
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    moris_id Derived_Node::get_id()
    {
        return -1;
    }

    //--------------------------------------------------------------------------------------------------------------
    
    moris_index Derived_Node::get_owner()
    {
        return par_rank();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Derived_Node::set_override_linear_interpolation()
    {
        gOverrideLinearInterpolation = true;
    }

    //--------------------------------------------------------------------------------------------------------------

}

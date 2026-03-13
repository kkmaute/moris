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
#include "cl_GEN_Node_Manager.hpp"
#include "cl_MTK_Interpolation_Function_Factory.hpp"
#include "cl_Communication_Tools.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Derived_Node::Derived_Node(
            uint aIndex,
            // const Vector< Background_Node* >& aBackgroundNodes, brendan delete
            const mtk::Cell&        aBackgroundElement,
            const Node_Manager&     aNodeManager,
            const Matrix< DDRMat >& aParametricCoordinates )
            : Node( aIndex )
            , mBackgroundElement( aBackgroundElement )
            , mNodeManager( aNodeManager )
            , mParametricCoordinates( aParametricCoordinates )
    {
        // Size global coordinates based on first locator
        mGlobalCoordinates = Matrix< DDRMat >( 1, aBackgroundElement.get_cell_info()->get_loc_coord_dim(), 0 );
        // delete tInterpolation;

        // Add contributions from all locators
        for ( auto iBasisNode : this->get_locator_nodes() )
        {
            mGlobalCoordinates += iBasisNode.get_global_coordinates() * iBasisNode.get_basis();
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Derived_Node::~Derived_Node()
    {
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

    Vector< const Background_Node* > Derived_Node::get_background_nodes() const
    {
        // Get the cell indices of the background element
        Matrix< IndexMat > tBackgroundElementNodeIndices = mBackgroundElement.get_vertex_inds();
        uint               tNumNodes                     = tBackgroundElementNodeIndices.length();

        // Initialize background nodes
        Vector< const Background_Node* > tBackgroundNodes( tNumNodes );
        tBackgroundNodes.reserve( tNumNodes );
        for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
        {
            tBackgroundNodes( iNode ) = &mNodeManager.get_background_node( tBackgroundElementNodeIndices( iNode ) );
        }

        return tBackgroundNodes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Basis_Node > Derived_Node::get_locator_nodes() const
    {
        // Create interpolator
        mtk::Interpolation_Function_Factory tInterpolationFactory;
        mtk::Interpolation_Function_Base*   tInterpolator = tInterpolationFactory.create_interpolation_function(
                mBackgroundElement.get_geometry_type(),
                mtk::Interpolation_Type::LAGRANGE,
                get_locator_interpolation_order() );


        // Evaluate the basis functions at the parametric coordinates
        Matrix< DDRMat > tBasis;
        tInterpolator->eval_N( mParametricCoordinates, tBasis );
        uint tNumNodes = tInterpolator->get_number_of_bases();

        // Get the cell indices of the background element
        Matrix< IndexMat > tBackgroundElementNodeIndices = mBackgroundElement.get_vertex_inds();

        // Initialize background nodes
        Vector< Basis_Node > tLocatorNodes;
        tLocatorNodes.reserve( tNumNodes );
        for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
        {
            tLocatorNodes.emplace_back( mNodeManager.get_background_node( tBackgroundElementNodeIndices( iNode ) ), tBasis( iNode ) );
        }

        return tLocatorNodes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Basis_Node > Derived_Node::get_field_basis_nodes() const
    {
        return this->get_locator_nodes();
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

    uint Derived_Node::get_num_pdvs() const
    {
        return 0;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Derived_Node::set_starting_pdv_id( moris_id aPDVStartingID )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id Derived_Node::get_starting_pdv_id() const
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

    moris_id Derived_Node::get_id() const
    {
        return -1;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index Derived_Node::get_owner() const
    {
        return par_rank();
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Derived_Node::get_background_element_nodal_coordinates() const
    {
        return mBackgroundElement.get_vertex_coords();
    }

    //--------------------------------------------------------------------------------------------------------------

    mtk::Geometry_Type Derived_Node::get_background_element_geometry_type() const
    {
        return mBackgroundElement.get_geometry_type();
    }

    //--------------------------------------------------------------------------------------------------------------

    mtk::Interpolation_Order Derived_Node::get_background_element_interpolation_order() const
    {
        return mBackgroundElement.get_interpolation_order();
    }

    //--------------------------------------------------------------------------------------------------------------

    mtk::Interpolation_Order Derived_Node::get_locator_interpolation_order() const
    {
        return gOverrideLinearInterpolation ? mtk::Interpolation_Order::LINEAR : mBackgroundElement.get_interpolation_order();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Derived_Node::override_linear_interpolation( bool aOverride )
    {
        gOverrideLinearInterpolation = aOverride;
    }

}    // namespace moris::gen

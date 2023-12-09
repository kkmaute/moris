/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Level_Set_Geometry.cpp
 *
 */

#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Intersection_Node_Linear.hpp"
#include "cl_GEN_Intersection_Node_Bilinear.hpp"
#include "cl_GEN_Derived_Node.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Level_Set_Parameters::Level_Set_Parameters( const ParameterList& aParameterList )
            : Field_Parameters( aParameterList )
            , mIsocontourThreshold( aParameterList.get< real >( "isocontour_threshold" ) )
            , mIsocontourTolerance( aParameterList.get< real >( "isocontour_tolerance" ) )
            , mIntersectionTolerance( aParameterList.get< real >( "intersection_tolerance" ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Level_Set_Geometry::Level_Set_Geometry(
            std::shared_ptr< Field > aField,
            Level_Set_Parameters     aParameters,
            Node_Manager&            aNodeManager )
            : Design_Field( aField, aParameters, aNodeManager )
            , mParameters( aParameters )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void Level_Set_Geometry::set_node_manager( Node_Manager& aNodeManager )
    {
        mNodeManager = &aNodeManager;
    }

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Mode
    Level_Set_Geometry::get_intersection_mode()
    {
        return Intersection_Mode::LEVEL_SET;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Level_Set_Geometry::get_isocontour_threshold()
    {
        return mParameters.mIsocontourThreshold;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Level_Set_Geometry::get_isocontour_tolerance()
    {
        return mParameters.mIsocontourTolerance;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Level_Set_Geometry::get_intersection_tolerance()
    {
        return mParameters.mIntersectionTolerance;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Level_Set_Geometry::depends_on_advs()
    {
        return mField->has_advs();
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Level_Set_Geometry::get_geometric_region(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aNodeCoordinates )
    {
        // If it's a base node, can get geometric region from field value
        if ( mNodeManager->is_base_node( aNodeIndex ) )
        {
            return this->determine_geometric_region( this->get_field_value( aNodeIndex, aNodeCoordinates ) );
        }
        else
        {
            // Get derived node
            Derived_Node* tDerivedNode = mNodeManager->get_derived_node( aNodeIndex );

            // If derived node knows it is on this interface, can return interface
            if ( tDerivedNode->is_on_interface( this ) )
            {
                return Geometric_Region::INTERFACE;
            }

            // Initialize possible region
            Geometric_Region tPossibleRegion = Geometric_Region::INTERFACE;

            // Test for locators
            for ( auto iLocator : tDerivedNode->get_locator_nodes() )
            {
                // Get locator region
                Geometric_Region tLocatorRegion = this->get_geometric_region( iLocator.get_index(), iLocator.get_global_coordinates() );

                // Update possible region
                if ( tPossibleRegion == Geometric_Region::INTERFACE )
                {
                    // Can be any possible region, so set as locator region
                    tPossibleRegion = tLocatorRegion;
                }
                else if ( tLocatorRegion == Geometric_Region::INTERFACE )
                {
                    // No change needed to possible region
                    continue;
                }
                else if ( tLocatorRegion != tPossibleRegion )
                {
                    // Resort to field value
                    return this->determine_geometric_region( this->get_field_value( aNodeIndex, aNodeCoordinates ) );
                }
            }

            // If nothing returned yet, possible region is definite region
            return tPossibleRegion;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node* Level_Set_Geometry::create_intersection_node(
            uint                 aNodeIndex,
            const Cell< Node* >& aBaseNodes,
            const Parent_Node&   aFirstParentNode,
            const Parent_Node&   aSecondParentNode,
            mtk::Geometry_Type   aBaseGeometryType )
    {
        if ( this->use_multilinear_interpolation() )
        {
            // Create multilinear intersection node
            return new Intersection_Node_Bilinear(
                    aNodeIndex,
                    aBaseNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    aBaseGeometryType,
                    shared_from_this() );
        }
        else
        {
            // Create linear intersection node
            return new Intersection_Node_Linear(
                    aNodeIndex,
                    aBaseNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    aBaseGeometryType,
                    shared_from_this() );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< std::shared_ptr< mtk::Field > > Level_Set_Geometry::get_mtk_fields()
    {
        // Get MTK field
        std::shared_ptr< mtk::Field > tMTKField = this->get_mtk_field();

        // Add to cell if it exists
        if ( tMTKField )
        {
            return { this->get_mtk_field() };
        }
        else
        {
            return {};
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometric_Region Level_Set_Geometry::determine_geometric_region( real aLevelSetValue )
    {
        // Determine if value indicates that this point is on the interface
        if ( std::abs( aLevelSetValue - mParameters.mIsocontourThreshold ) <= mParameters.mIsocontourTolerance )
        {
            return Geometric_Region::INTERFACE;
        }

        // Otherwise, give the region relative to the isocontour threshold
        else if ( aLevelSetValue - mParameters.mIsocontourThreshold < 0 )
        {
            return Geometric_Region::NEGATIVE;
        }
        else
        {
            return Geometric_Region::POSITIVE;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}

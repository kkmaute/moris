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

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Level_Set_Parameters::Level_Set_Parameters( const ParameterList& aParameterList )
            : Field_Parameters( aParameterList )
            , mIntersectionInterpolation( aParameterList.get< bool >( "multilinear_intersections" )
                    ? Int_Interpolation::MULTILINEAR : Int_Interpolation::LINEAR )
            , mIsocontourThreshold( aParameterList.get< real >( "isocontour_threshold" ) )
            , mIsocontourTolerance( aParameterList.get< real >( "isocontour_tolerance" ) )
            , mIntersectionTolerance( aParameterList.get< real >( "intersection_tolerance" ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Level_Set_Geometry::Level_Set_Geometry(
            std::shared_ptr< Field > aField,
            Level_Set_Parameters     aParameters )
            : Design_Field( aField, aParameters )
            , mParameters( aParameters )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Int_Interpolation
    Level_Set_Geometry::get_intersection_interpolation()
    {
        return mParameters.mIntersectionInterpolation;
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

    std::shared_ptr< Intersection_Node > Level_Set_Geometry::create_intersection_node(
            uint                                 aEdgeFirstNodeIndex,
            uint                                 aEdgeSecondNodeIndex,
            std::shared_ptr< Intersection_Node > aEdgeFirstIntersectionNode,
            std::shared_ptr< Intersection_Node > aEdgeSecondIntersectionNode,
            const Matrix< DDRMat >&              aEdgeFirstNodeLocalCoordinates,
            const Matrix< DDRMat >&              aEdgeSecondNodeLocalCoordinates,
            const Matrix< DDRMat >&              aEdgeFirstNodeGlobalCoordinates,
            const Matrix< DDRMat >&              aEdgeSecondNodeGlobalCoordinates,
            const Matrix< DDUMat >&              aBackgroundElementNodeIndices,
            const Cell< Matrix< DDRMat > >&      aBackgroundElementNodeCoordinates )
    {
        if ( mParameters.mIntersectionInterpolation == Int_Interpolation::LINEAR )
        {
            // Create linear intersection node
            return std::make_shared< Intersection_Node_Linear >(
                    aEdgeFirstIntersectionNode,
                    aEdgeSecondIntersectionNode,
                    aEdgeFirstNodeIndex,
                    aEdgeSecondNodeIndex,
                    aEdgeFirstNodeGlobalCoordinates,
                    aEdgeSecondNodeGlobalCoordinates,
                    shared_from_this() );
        }
        else
        {
            // Get element interpolation type TODO move this into the bilinear intersection node, it should know what spatial dimension it is
            Element_Interpolation_Type tInterpolationType =
                    aEdgeFirstNodeLocalCoordinates.length() == 2 ? Element_Interpolation_Type::Linear_2D : Element_Interpolation_Type::Linear_3D;

            // Create multilinear intersection node
            return std::make_shared< Intersection_Node_Bilinear >(
                    aEdgeFirstIntersectionNode,
                    aEdgeSecondIntersectionNode,
                    aEdgeFirstNodeIndex,
                    aEdgeSecondNodeIndex,
                    aEdgeFirstNodeLocalCoordinates,
                    aEdgeSecondNodeLocalCoordinates,
                    aBackgroundElementNodeIndices,
                    aBackgroundElementNodeCoordinates,
                    tInterpolationType,
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

}

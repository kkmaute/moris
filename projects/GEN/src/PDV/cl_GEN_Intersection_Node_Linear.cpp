/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Linear.cpp
 *
 */

#include "cl_GEN_Intersection_Node_Linear.hpp"
#include "cl_GEN_Parent_Node.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Interpolation.hpp"

#include "cl_MTK_Interpolation_Function_Base.hpp"
#include "cl_MTK_Interpolation_Function_Factory.hpp"
#include "cl_MTK_Enums.hpp"

#include "fn_dot.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node_Linear::Intersection_Node_Linear(
            uint                                  aNodeIndex,
            const Cell< Node* >&                  aBaseNodes,
            const Parent_Node&                    aFirstParentNode,
            const Parent_Node&                    aSecondParentNode,
            mtk::Geometry_Type                    aBackgroundGeometryType,
            mtk::Interpolation_Order              aBackgroundInterpolationOrder,
            std::shared_ptr< Level_Set_Geometry > aInterfaceGeometry )
            : Intersection_Node_Level_Set(
                    aNodeIndex,
                    aBaseNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    Intersection_Node_Linear::compute_local_coordinate( aFirstParentNode, aSecondParentNode, aInterfaceGeometry ),
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder,
                    aInterfaceGeometry )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const Cell< Basis_Node >& Intersection_Node_Linear::get_field_basis_nodes() const
    {
        return this->get_locator_nodes();
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Intersection_Node_Linear::get_dxi_dfield_from_ancestor( uint aAncestorIndex ) const
    {
        // Locked interface geometry
        std::shared_ptr< Level_Set_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

        // get isocontour threshold from geometry
        real tIsocontourThreshold = tLockedInterfaceGeometry->get_isocontour_threshold();

        // Get geometry field values
        real tPhi0 = tLockedInterfaceGeometry->get_field_value( this->get_first_parent_node().get_index(), this->get_first_parent_node().get_global_coordinates() );
        real tPhi1 = tLockedInterfaceGeometry->get_field_value( this->get_second_parent_node().get_index(), this->get_second_parent_node().get_global_coordinates() );

        // Compute sensitivity of the local coordinate with respect to the field value
        return 2 * ( ( tPhi0 - tIsocontourThreshold ) * ( aAncestorIndex == 1 ) - ( tPhi1 - tIsocontourThreshold ) * ( aAncestorIndex == 0 ) ) / std::pow( ( tPhi1 - tPhi0 ), 2 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Linear::get_dxi_dcoordinate_first_parent() const
    {
        // Locked interface geometry
        std::shared_ptr< Level_Set_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

        // Compute sensitivity of the local coordinate with respect to the ancestor coordinates
        Matrix< DDRMat > tCoordinateSensitivities( 1, this->get_global_coordinates().n_cols(), 0.0 );
        tLockedInterfaceGeometry->get_dfield_dcoordinates(
                this->get_first_parent_node(),
                tCoordinateSensitivities );
        return this->get_dxi_dfield_from_ancestor( 0 ) * tCoordinateSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Linear::get_dxi_dcoordinate_second_parent() const
    {
        // Locked interface geometry
        std::shared_ptr< Level_Set_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

        // Compute sensitivity of the local coordinate with respect to the ancestor coordinates
        Matrix< DDRMat > tCoordinateSensitivities( 1, this->get_global_coordinates().n_cols() );
        tLockedInterfaceGeometry->get_dfield_dcoordinates(
                this->get_second_parent_node(),
                tCoordinateSensitivities );

        return this->get_dxi_dfield_from_ancestor( 1 ) * tCoordinateSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Intersection_Node_Linear::compute_local_coordinate(
            const Parent_Node&                    aFirstParentNode,
            const Parent_Node&                    aSecondParentNode,
            std::shared_ptr< Level_Set_Geometry > aInterfaceGeometry )
    {
        // Interface geometry values
        Matrix< DDRMat > tInterfaceGeometryValues = { { aInterfaceGeometry->get_field_value( aFirstParentNode.get_index(), aFirstParentNode.get_global_coordinates() ) },
            { aInterfaceGeometry->get_field_value( aSecondParentNode.get_index(), aSecondParentNode.get_global_coordinates() ) } };

        // Get isocontour threshold
        real tIsocontourThreshold = aInterfaceGeometry->get_isocontour_threshold();

        // Interpolate
        Matrix< DDRMat > tLocalCoordinates = Interpolation::linear_interpolation_value( tInterfaceGeometryValues, tIsocontourThreshold );

        return tLocalCoordinates( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

}

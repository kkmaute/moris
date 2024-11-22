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
#include "cl_GEN_Basis_Node.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node_Linear::Intersection_Node_Linear(
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
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< Basis_Node >& Intersection_Node_Linear::get_field_basis_nodes() const
    {
        return this->get_locator_nodes();
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Intersection_Node_Linear::get_dxi_dfield_from_ancestor( uint aAncestorIndex ) const
    {
        // get isocontour threshold from geometry
        real tIsocontourThreshold = mInterfaceGeometry.get_isocontour_threshold();

        // Get geometry field values
        real tPhi0 = mInterfaceGeometry.get_field_value( this->get_first_parent_node().get_index(), this->get_first_parent_node().get_global_coordinates() );
        real tPhi1 = mInterfaceGeometry.get_field_value( this->get_second_parent_node().get_index(), this->get_second_parent_node().get_global_coordinates() );

        // Compute sensitivity of the local coordinate with respect to the field value
        return 2 * ( ( tPhi0 - tIsocontourThreshold ) * ( aAncestorIndex == 1 ) - ( tPhi1 - tIsocontourThreshold ) * ( aAncestorIndex == 0 ) ) / std::pow( ( tPhi1 - tPhi0 ), 2 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Linear::get_dxi_dcoordinate_first_parent() const
    {
        // Compute sensitivity of the local coordinate with respect to the ancestor coordinates
        Matrix< DDRMat > tCoordinateSensitivities( 1, this->get_global_coordinates().length(), 0.0 );
        mInterfaceGeometry.get_dfield_dcoordinates(
                this->get_first_parent_node(),
                tCoordinateSensitivities );
        return this->get_dxi_dfield_from_ancestor( 0 ) * tCoordinateSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Linear::get_dxi_dcoordinate_second_parent() const
    {
        // Compute sensitivity of the local coordinate with respect to the ancestor coordinates
        Matrix< DDRMat > tCoordinateSensitivities( 1, this->get_global_coordinates().length() );
        mInterfaceGeometry.get_dfield_dcoordinates(
                this->get_second_parent_node(),
                tCoordinateSensitivities );

        return this->get_dxi_dfield_from_ancestor( 1 ) * tCoordinateSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

}

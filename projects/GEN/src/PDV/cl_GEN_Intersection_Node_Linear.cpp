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
            mtk::Geometry_Type                    aBaseGeometryType,
            std::shared_ptr< Level_Set_Geometry > aInterfaceGeometry )
            : Intersection_Node_Level_Set(
                    aNodeIndex,
                    aBaseNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    Intersection_Node_Linear::compute_local_coordinate( aFirstParentNode, aSecondParentNode, aInterfaceGeometry ),
                    aBaseGeometryType,
                    aInterfaceGeometry )
    {
        // call required setup function
        this->initialize();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Intersection_Node_Linear::get_dfield_dcoordinates(
            Field*            aField,
            Matrix< DDRMat >& aSensitivities )
    {
        // Geometry values
        real tDeltaPhi = aField->get_field_value( mSecondParentNode.get_index(), mSecondParentNode.get_global_coordinates() )
                       - aField->get_field_value( mFirstParentNode.get_index(), mFirstParentNode.get_global_coordinates() );

        // get number of spatial dimensions;
        uint tNumDim = mParentVector.length();

        // Compute square of length of parent vector
        real tParentLengthSquared = dot( mParentVector, mParentVector );

        // Sensitivities: dPhi/dx_i  = delta(Phi) / L_i where L_i = PaerentVectorLenth^2 / (ParentVector * e_i)
        for ( uint tCoordinateIndex = 0; tCoordinateIndex < tNumDim; tCoordinateIndex++ )
        {
            aSensitivities( tCoordinateIndex ) = tDeltaPhi * mParentVector( tCoordinateIndex ) / ( tParentLengthSquared + MORIS_REAL_EPS );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Intersection_Node_Linear::get_dxi_dfield_from_ancestor( uint aAncestorIndex )
    {
        // Locked interface geometry
        std::shared_ptr< Level_Set_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

        // get isocontour threshold from geometry
        real tIsocontourThreshold = tLockedInterfaceGeometry->get_isocontour_threshold();

        // Get geometry field values
        real tPhi0 = tLockedInterfaceGeometry->get_field_value( mFirstParentNode.get_index(), mFirstParentNode.get_global_coordinates() );
        real tPhi1 = tLockedInterfaceGeometry->get_field_value( mSecondParentNode.get_index(), mSecondParentNode.get_global_coordinates() );

        // Compute sensitivity of the local coordinate with respect to the field value
        return 2 * ( ( tPhi0 - tIsocontourThreshold ) * ( aAncestorIndex == 1 ) - ( tPhi1 - tIsocontourThreshold ) * ( aAncestorIndex == 0 ) ) / std::pow( ( tPhi1 - tPhi0 ), 2 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Linear::get_dxi_dcoordinate_first_parent()
    {
        // Locked interface geometry
        std::shared_ptr< Level_Set_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

        // Compute sensitivity of the local coordinate with respect to the ancestor coordinates
        Matrix< DDRMat > tCoordinateSensitivities( 1, this->get_global_coordinates().n_cols() );
        tLockedInterfaceGeometry->get_dfield_dcoordinates(
                mFirstParentNode.get_index(),
                mSecondParentNode.get_global_coordinates(),
                tCoordinateSensitivities );
        return this->get_dxi_dfield_from_ancestor( 0 ) * tCoordinateSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Linear::get_dxi_dcoordinate_second_parent()
    {
        // Locked interface geometry
        std::shared_ptr< Level_Set_Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

        // Compute sensitivity of the local coordinate with respect to the ancestor coordinates
        Matrix< DDRMat > tCoordinateSensitivities( 1, this->get_global_coordinates().n_cols() );
        tLockedInterfaceGeometry->get_dfield_dcoordinates(
                mSecondParentNode.get_index(),
                mSecondParentNode.get_global_coordinates(),
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

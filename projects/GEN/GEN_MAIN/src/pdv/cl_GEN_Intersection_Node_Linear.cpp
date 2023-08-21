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
#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Interpolation.hpp"
#include "cl_XTK_Linear_Basis_Functions.hpp"
#include "fn_trans.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Node_Linear::Intersection_Node_Linear(
                std::shared_ptr< Intersection_Node > aFirstNode,
                std::shared_ptr< Intersection_Node > aSecondNode,
                uint                                 aFirstNodeIndex,
                uint                                 aSecondNodeIndex,
                const Matrix< DDRMat >&              aFirstNodeCoordinates,
                const Matrix< DDRMat >&              aSecondNodeCoordinates,
                std::shared_ptr< Geometry >          aInterfaceGeometry )
                : Intersection_Node(
                        get_local_coordinate(
                                aFirstNodeIndex,
                                aSecondNodeIndex,
                                aFirstNodeCoordinates,
                                aSecondNodeCoordinates,
                                aInterfaceGeometry ),
                        aFirstNode,
                        aSecondNode,
                        aFirstNodeIndex,
                        aSecondNodeIndex,
                        { { -1 } },
                        { { 1 } },
                        { { aFirstNodeIndex, aSecondNodeIndex } },
                        { aFirstNodeCoordinates, aSecondNodeCoordinates },
                        Element_Intersection_Type::Linear_1D,
                        aInterfaceGeometry )
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Intersection_Node_Linear::get_dfield_dcoordinates(
                Field*            aField,
                Matrix< DDRMat >& aSensitivities )
        {
            // Geometry values
            real tDeltaPhi = aField->get_field_value( mAncestorNodeIndices( 1 ), mAncestorNodeCoordinates( 1 ) )    //
                           - aField->get_field_value( mAncestorNodeIndices( 0 ), mAncestorNodeCoordinates( 0 ) );

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
            std::shared_ptr< Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // get isocontour threshold from geometry
            real tIsocontourThreshold = tLockedInterfaceGeometry->get_isocontour_threshold();

            // Get geometry field values
            real tPhi0 = tLockedInterfaceGeometry->get_field_value( mAncestorNodeIndices( 0 ), mAncestorNodeCoordinates( 0 ) );
            real tPhi1 = tLockedInterfaceGeometry->get_field_value( mAncestorNodeIndices( 1 ), mAncestorNodeCoordinates( 1 ) );

            // Compute sensitivity of the local coordinate with respect to the field value
            return 2 * ( ( tPhi0 - tIsocontourThreshold ) * ( aAncestorIndex == 1 ) - ( tPhi1 - tIsocontourThreshold ) * ( aAncestorIndex == 0 ) ) / std::pow( ( tPhi1 - tPhi0 ), 2 );
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Intersection_Node_Linear::get_dxi_dcoordinate_first_parent()
        {
            // Locked interface geometry
            std::shared_ptr< Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // Compute sensitivity of the local coordinate with respect to the ancestor coordinates
            Matrix< DDRMat > tCoordinateSensitivities( 1, mGlobalCoordinates.n_cols() );
            tLockedInterfaceGeometry->get_dfield_dcoordinates(
                    mAncestorNodeIndices( 0 ),
                    mAncestorNodeCoordinates( 0 ),
                    tCoordinateSensitivities );
            return this->get_dxi_dfield_from_ancestor( 0 ) * tCoordinateSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Intersection_Node_Linear::get_dxi_dcoordinate_second_parent()
        {
            // Locked interface geometry
            std::shared_ptr< Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // Compute sensitivity of the local coordinate with respect to the ancestor coordinates
            Matrix< DDRMat > tCoordinateSensitivities( 1, mGlobalCoordinates.n_cols() );
            tLockedInterfaceGeometry->get_dfield_dcoordinates(
                    mAncestorNodeIndices( 1 ),
                    mAncestorNodeCoordinates( 1 ),
                    tCoordinateSensitivities );

            return this->get_dxi_dfield_from_ancestor( 1 ) * tCoordinateSensitivities;
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Intersection_Node_Linear::get_local_coordinate(
                uint                        aFirstNodeIndex,
                uint                        aSecondNodeIndex,
                const Matrix< DDRMat >&     aFirstNodeCoordinates,
                const Matrix< DDRMat >&     aSecondNodeCoordinates,
                std::shared_ptr< Geometry > aInterfaceGeometry )
        {
            // Interface geometry values
            Matrix< DDRMat > tInterfaceGeometryValues = { { aInterfaceGeometry->get_field_value( aFirstNodeIndex, aFirstNodeCoordinates ) },
                { aInterfaceGeometry->get_field_value( aSecondNodeIndex, aSecondNodeCoordinates ) } };

            // Get isocontour threshold
            real tIsocontourThreshold = aInterfaceGeometry->get_isocontour_threshold();

            // Interpolate
            Matrix< DDRMat > tLocalCoordinates = Interpolation::linear_interpolation_value( tInterfaceGeometryValues, tIsocontourThreshold );

            return tLocalCoordinates( 0 );
        }

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris

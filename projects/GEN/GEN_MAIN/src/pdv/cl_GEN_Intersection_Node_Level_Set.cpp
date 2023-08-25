/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Level_Set.cpp
 *
 */

#include "cl_GEN_Intersection_Node_Level_Set.hpp"
#include "cl_GEN_Geometry.hpp"

#include "fn_eye.hpp"

namespace moris
{
    namespace ge
    {
        Intersection_Node_Level_Set::Intersection_Node_Level_Set(
                real                                 aLocalCoordinate,
                std::shared_ptr< Intersection_Node > aFirstParentNode,
                std::shared_ptr< Intersection_Node > aSecondParentNode,
                moris_index                          aFirstParentNodeIndex,
                moris_index                          aSecondParentNodeIndex,
                const Matrix< DDRMat >&              aFirstParentNodeLocalCoordinates,
                const Matrix< DDRMat >&              aSecondParentNodeLocalCoordinates,
                Matrix< DDUMat >                     aAncestorNodeIndices,
                Cell< Matrix< DDRMat > >             aAncestorNodeCoordinates,
                const Element_Intersection_Type      aAncestorBasisFunction,
                std::shared_ptr< Geometry >          aInterfaceGeometry,
                bool                                 aDetermineIsIntersected )
                : Intersection_Node(
                        aLocalCoordinate,
                        aFirstParentNode,
                        aSecondParentNode,
                        aFirstParentNodeIndex,
                        aSecondParentNodeIndex,
                        aFirstParentNodeLocalCoordinates,
                        aSecondParentNodeLocalCoordinates,
                        aAncestorNodeIndices,
                        aAncestorNodeCoordinates,
                        aAncestorBasisFunction,
                        aInterfaceGeometry,
                        aDetermineIsIntersected )
        {
        }




        void
        Intersection_Node_Level_Set::get_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, Matrix< DDRMat > const & aSensitivityFactor )
        {
            // Locked interface geometry
            std::shared_ptr< Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // Get sensitivity values from other ancestors
            Matrix< DDRMat > tSensitivitiesToAdd;
            for ( uint tAncestorNode = 0; tAncestorNode < mAncestorNodeIndices.length(); tAncestorNode++ )
            {
                // Get geometry field sensitivity with respect to ADVs
                const Matrix< DDRMat >& tFieldSensitivities = tLockedInterfaceGeometry->get_dfield_dadvs(
                        mAncestorNodeIndices( tAncestorNode ),
                        mAncestorNodeCoordinates( tAncestorNode ) );

                // Ancestor sensitivities
                tSensitivitiesToAdd =
                        0.5 * aSensitivityFactor * this->get_dxi_dfield_from_ancestor( tAncestorNode ) * mParentVector * tFieldSensitivities;

                // Resize sensitivities
                uint tJoinedSensitivityLength = aCoordinateSensitivities.n_cols();
                aCoordinateSensitivities.resize( tSensitivitiesToAdd.n_rows(),
                        tJoinedSensitivityLength + tSensitivitiesToAdd.n_cols() );

                // Join sensitivities
                for ( uint tCoordinateIndex = 0; tCoordinateIndex < tSensitivitiesToAdd.n_rows(); tCoordinateIndex++ )
                {
                    for ( uint tAddedSensitivity = 0; tAddedSensitivity < tSensitivitiesToAdd.n_cols(); tAddedSensitivity++ )
                    {
                        aCoordinateSensitivities( tCoordinateIndex, tJoinedSensitivityLength + tAddedSensitivity ) =
                                tSensitivitiesToAdd( tCoordinateIndex, tAddedSensitivity );
                    }
                }
            }

            // Add first parent sensitivities, if needed
            if ( mFirstParentNode )
            {
                Matrix< DDRMat > tLocCoord = ( 1.0 - mLocalCoordinate ) *    //
                                             eye( mParentVector.n_rows(), mParentVector.n_rows() );

                Matrix< DDRMat > tSensitivityFactor = 0.5 * ( tLocCoord + mParentVector * this->get_dxi_dcoordinate_first_parent() );
                mFirstParentNode->get_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );
            }

            // Add second parent sensitivities, if needed
            if ( mSecondParentNode )
            {
                Matrix< DDRMat > tLocCoord = ( 1.0 + mLocalCoordinate ) *    //
                                             eye( mParentVector.n_rows(), mParentVector.n_rows() );

                Matrix< DDRMat > tSensitivityFactor = 0.5 * ( tLocCoord + mParentVector * this->get_dxi_dcoordinate_second_parent() );
                mSecondParentNode->get_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDSMat >
        Intersection_Node_Level_Set::get_coordinate_determining_adv_ids()
        {
            // Set size
            mCoordinateDeterminingADVIDs.set_size( 0, 0 );

            // Locked interface geometry
            std::shared_ptr< Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // Get sensitivity values from other ancestors
            for ( uint tAncestorNode = 0; tAncestorNode < mAncestorNodeIndices.length(); tAncestorNode++ )
            {
                // Get geometry field sensitivity with respect to ADVs
                const Matrix< DDSMat >& tAncestorADVIDs = tLockedInterfaceGeometry->get_determining_adv_ids(
                        mAncestorNodeIndices( tAncestorNode ),
                        mAncestorNodeCoordinates( tAncestorNode ) );

                // Join IDs
                this->join_adv_ids( tAncestorADVIDs );
            }

            // Add first parent IDs, if needed
            if ( mFirstParentNode )
            {
                this->join_adv_ids( mFirstParentNode->get_coordinate_determining_adv_ids() );
            }

            // Add second parent IDs, if needed
            if ( mSecondParentNode )
            {
                this->join_adv_ids( mSecondParentNode->get_coordinate_determining_adv_ids() );
            }

            return mCoordinateDeterminingADVIDs;
        }
        
        //--------------------------------------------------------------------------------------------------------------

    }    // namespace ge

}    // namespace moris
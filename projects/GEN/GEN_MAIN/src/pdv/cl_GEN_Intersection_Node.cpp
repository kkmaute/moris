/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node.cpp
 *
 */

#include "cl_GEN_Intersection_Node.hpp"
#include "cl_GEN_Geometry.hpp"

#include "cl_MTK_Interpolation_Function_Base.hpp"       //MTK/src
#include "cl_MTK_Interpolation_Function_Factory.hpp"    //MTK/src
#include "cl_MTK_Enums.hpp"                             //MTK/src

#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Intersection_Node::Intersection_Node(
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
                real                                 aIsocontourThreshold,
                real                                 aIsocontourTolerance,
                real                                 aIntersectionTolerance,
                bool                                 aDetermineIsIntersected )
                : Child_Node(
                        aAncestorNodeIndices,
                        aAncestorNodeCoordinates,
                        aAncestorBasisFunction,
                        0.5 * ( ( 1 - aLocalCoordinate ) * aFirstParentNodeLocalCoordinates + ( 1 + aLocalCoordinate ) * aSecondParentNodeLocalCoordinates ) )
                , mLocalCoordinate( aLocalCoordinate )
                , mIsIntersected( false )
                , mInterfaceGeometry( aInterfaceGeometry )
                , mIsocontourThreshold( aIsocontourThreshold )
                , mFirstParentNode( aFirstParentNode )
                , mSecondParentNode( aSecondParentNode )
                , mFirstParentNodeIndex( aFirstParentNodeIndex )
                , mSecondParentNodeIndex( aSecondParentNodeIndex )
        {
            // check whether intersection computation should be performed
            if ( !aDetermineIsIntersected )
            {
                return;
            }

            // construct interpolator
            mtk::Interpolation_Function_Factory tFactory;

            mtk::Interpolation_Function_Base* tInterpolation;

            switch ( aAncestorBasisFunction )
            {
                case Element_Intersection_Type::Linear_1D:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::LINE,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );
                    break;
                }
                case Element_Intersection_Type::Linear_2D:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::QUAD,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );
                    break;
                }
                case Element_Intersection_Type::Linear_3D:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::HEX,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );
                    break;
                }
                default:
                    MORIS_ERROR( false,
                            "Intersection_Node::Intersection_Node - Interpolation type not implemented." );
            }

            // Parent basis
            Matrix< DDRMat > tFirstParentBasisValues;
            Matrix< DDRMat > tSecondParentBasisValues;

            tInterpolation->eval_N( aFirstParentNodeLocalCoordinates, tFirstParentBasisValues );
            tInterpolation->eval_N( aSecondParentNodeLocalCoordinates, tSecondParentBasisValues );

            // Global coordinates of intersection and parents
            mGlobalCoordinates = mBasisValues( 0 ) * aAncestorNodeCoordinates( 0 );

            Matrix< DDRMat > tFirstParentGlobalCoordinates  = tFirstParentBasisValues( 0 ) * aAncestorNodeCoordinates( 0 );
            Matrix< DDRMat > tSecondParentGlobalCoordinates = tSecondParentBasisValues( 0 ) * aAncestorNodeCoordinates( 0 );

            for ( uint tBasisIndex = 1; tBasisIndex < mBasisValues.length(); tBasisIndex++ )
            {
                mGlobalCoordinates += mBasisValues( tBasisIndex ) * aAncestorNodeCoordinates( tBasisIndex );
                tFirstParentGlobalCoordinates += tFirstParentBasisValues( tBasisIndex ) * aAncestorNodeCoordinates( tBasisIndex );
                tSecondParentGlobalCoordinates += tSecondParentBasisValues( tBasisIndex ) * aAncestorNodeCoordinates( tBasisIndex );
            }

            mParentVector = trans( tSecondParentGlobalCoordinates - tFirstParentGlobalCoordinates );

            // Parents on interface
            real tFirstParentPhi  = aInterfaceGeometry->get_field_value( aFirstParentNodeIndex, tFirstParentGlobalCoordinates );
            real tSecondParentPhi = aInterfaceGeometry->get_field_value( aSecondParentNodeIndex, tSecondParentGlobalCoordinates );

            real tParentLength = norm( mParentVector );

            real tFirstDiffFromThreshold  = tFirstParentPhi - aIsocontourThreshold;
            real tSecondDiffFromThreshold = tSecondParentPhi - aIsocontourThreshold;

            mFirstParentOnInterface = std::abs( tFirstDiffFromThreshold ) < aIntersectionTolerance
                                   or 0.5 * tParentLength * std::abs( 1 + aLocalCoordinate ) < aIntersectionTolerance;

            mSecondParentOnInterface = std::abs( tSecondDiffFromThreshold ) < aIntersectionTolerance
                                    or 0.5 * tParentLength * std::abs( 1 - aLocalCoordinate ) < aIntersectionTolerance;

            // Determine if edge is intersected
            if ( mFirstParentOnInterface or mSecondParentOnInterface )
            {
                mIsIntersected = true;
            }
            // FIXME: This check should be unnecessary as the local edge coordinate should be sufficient
            // to determine whether edge is intersected; it is only "useful" if parent node's level set value
            // is determined by method that is different from intersection nodes; for example level set value child node
            // of child node is computed via analytic field and intersection node via bi-linear interpolation
            else if ( tFirstDiffFromThreshold * tSecondDiffFromThreshold > 0 )
            {
                mIsIntersected = false;

                // check for consistency of parent values and local coordinate
                MORIS_ASSERT( std::abs( mLocalCoordinate ) > 1,
                        "Intersection_Node::Intersection_Node - inconsistent parent level set values versus local coordinate - p1 %e p2 %e loc %e.",
                        tFirstDiffFromThreshold,
                        tSecondDiffFromThreshold,
                        mLocalCoordinate );
            }
            else
            {
                mIsIntersected = ( std::abs( mLocalCoordinate ) <= 1.0 );

                // check for consistency with parent values
                // this check is currently useless but should be performed is inconsistency issue (see comment above) is resolved
                MORIS_ASSERT( mIsIntersected ? tFirstDiffFromThreshold * tSecondDiffFromThreshold < 0 : tFirstDiffFromThreshold * tSecondDiffFromThreshold > 0,
                        "Intersection_Node::Intersection_Node - inconsistent parent level set values - p1 %e p2 %e loc %e.",
                        tFirstDiffFromThreshold,
                        tSecondDiffFromThreshold,
                        mLocalCoordinate );
            }

            // Delete interpolator
            delete tInterpolation;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Intersection_Node::get_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, Matrix< DDRMat > const & aSensitivityFactor )
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
        Intersection_Node::get_coordinate_determining_adv_ids()
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

        bool
        Intersection_Node::parent_edge_is_intersected()
        {
            return mIsIntersected;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool
        Intersection_Node::first_parent_on_interface()
        {
            return mFirstParentOnInterface;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool
        Intersection_Node::second_parent_on_interface()
        {
            return mSecondParentOnInterface;
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Intersection_Node::get_local_coordinate()
        {
            return mLocalCoordinate;
        }

        //--------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Intersection_Node::get_global_coordinates()
        {
            return mGlobalCoordinates;
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Intersection_Node::get_coordinate_value( uint aCoordinateIndex )
        {
            return mGlobalCoordinates( aCoordinateIndex );
        }

        //--------------------------------------------------------------------------------------------------------------

        uint
        Intersection_Node::get_num_pdvs()
        {
            return mGlobalCoordinates.numel();
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Intersection_Node::set_starting_pdv_id( moris_id aPDVStartingID )
        {
            mPDVStartingID    = aPDVStartingID;
            mPDVStartingIDSet = true;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_id
        Intersection_Node::get_starting_pdv_id()
        {
            MORIS_ASSERT( mPDVStartingIDSet, "PDV Starting ID must be set for an intersection." );
            return mPDVStartingID;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Intersection_Node::set_id( moris_id aNodeID )
        {
            mNodeID = aNodeID;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Intersection_Node::set_owner( moris_index aNodeOwner )
        {
            mNodeOwner = aNodeOwner;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_id
        Intersection_Node::get_id()
        {
            return mNodeID;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_index
        Intersection_Node::get_owner()
        {
            return mNodeOwner;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Intersection_Node::join_adv_ids( const Matrix< DDSMat >& aIDsToAdd )
        {
            // Resize IDs
            uint tJoinedSensitivityLength = mCoordinateDeterminingADVIDs.n_cols();
            mCoordinateDeterminingADVIDs.resize( 1, tJoinedSensitivityLength + aIDsToAdd.length() );

            // Join IDs
            for ( uint tAddedSensitivity = 0; tAddedSensitivity < aIDsToAdd.length(); tAddedSensitivity++ )
            {
                mCoordinateDeterminingADVIDs( tJoinedSensitivityLength + tAddedSensitivity ) = aIDsToAdd( tAddedSensitivity );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris

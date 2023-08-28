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
                bool                                 aDetermineIsIntersected )
                : Child_Node(
                        aAncestorNodeIndices,
                        aAncestorNodeCoordinates,
                        aAncestorBasisFunction,
                        0.5 * ( ( 1 - aLocalCoordinate ) * aFirstParentNodeLocalCoordinates + ( 1 + aLocalCoordinate ) * aSecondParentNodeLocalCoordinates ) )
                , mLocalCoordinate( aLocalCoordinate )
                , mIsIntersected( false )
                , mFirstParentNode( aFirstParentNode )
                , mSecondParentNode( aSecondParentNode )
                , mFirstParentNodeIndex( aFirstParentNodeIndex )
                , mSecondParentNodeIndex( aSecondParentNodeIndex )
        {
        }

        void
        Intersection_Node::initialize(
                const Element_Intersection_Type aAncestorBasisFunction,
                const Matrix< DDRMat >&         aFirstParentLocalCoordinates,
                const Matrix< DDRMat >&         aSecondParentLocalCoordinates )
        {
            mGlobalCoordinates       = compute_global_coordinates();
            mParentVector            = compute_parent_vector( aAncestorBasisFunction, aFirstParentLocalCoordinates, aSecondParentLocalCoordinates );
            mFirstParentOnInterface  = determine_first_parent_on_interface( aAncestorBasisFunction, aFirstParentLocalCoordinates );
            mSecondParentOnInterface = determine_second_parent_on_interface( aAncestorBasisFunction, aSecondParentLocalCoordinates );
            mIsIntersected           = determine_is_intersected();
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Intersection_Node::compute_parent_vector(
                const Element_Intersection_Type aAncestorBasisFunction,
                const Matrix< DDRMat >&         aFirstParentNodeLocalCoordinates,
                const Matrix< DDRMat >&         aSecondParentNodeLocalCoordinates )
        {
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

            Matrix< DDRMat > tFirstParentGlobalCoordinates  = tFirstParentBasisValues( 0 ) * mAncestorNodeCoordinates( 0 );
            Matrix< DDRMat > tSecondParentGlobalCoordinates = tSecondParentBasisValues( 0 ) * mAncestorNodeCoordinates( 0 );

            for ( uint tBasisIndex = 1; tBasisIndex < mBasisValues.length(); tBasisIndex++ )
            {
                tFirstParentGlobalCoordinates += tFirstParentBasisValues( tBasisIndex ) * mAncestorNodeCoordinates( tBasisIndex );
                tSecondParentGlobalCoordinates += tSecondParentBasisValues( tBasisIndex ) * mAncestorNodeCoordinates( tBasisIndex );
            }

            Matrix< DDRMat > tParentVector = trans( tSecondParentGlobalCoordinates - tFirstParentGlobalCoordinates );

            // delete interpolator
            delete tInterpolation;

            return tParentVector;
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

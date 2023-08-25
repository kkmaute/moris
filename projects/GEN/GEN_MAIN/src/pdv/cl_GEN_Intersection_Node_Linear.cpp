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

#include "cl_MTK_Interpolation_Function_Base.hpp"       //MTK/src
#include "cl_MTK_Interpolation_Function_Factory.hpp"    //MTK/src
#include "cl_MTK_Enums.hpp"                             //MTK/src

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
                : Intersection_Node_Level_Set(
                        compute_local_coordinate(
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
            // call required setup function
            this->initialize( Element_Intersection_Type::Linear_1D, { { -1 } }, { { 1 } } );
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Intersection_Node_Linear::compute_global_coordinates()
        {
            // Global coordinates of intersection and parents
            Matrix< DDRMat > tGlobalCoordinates = mBasisValues( 0 ) * mAncestorNodeCoordinates( 0 );

            for ( uint tBasisIndex = 1; tBasisIndex < mBasisValues.length(); tBasisIndex++ )
            {
                tGlobalCoordinates += mBasisValues( tBasisIndex ) * mAncestorNodeCoordinates( tBasisIndex );
            }

            return tGlobalCoordinates;
        }

        //--------------------------------------------------------------------------------------------------------------

        

        //--------------------------------------------------------------------------------------------------------------

        bool
        Intersection_Node_Linear::determine_first_parent_on_interface( 
            const Element_Intersection_Type aAncestorBasisFunction,
            const Matrix< DDRMat >& aFirstParentNodeLocalCoordinates )
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

            // lock the interface geometry
            std::shared_ptr< Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // get the isocontour thresholds from the geometry
            real tIsocontourThreshold   = tLockedInterfaceGeometry->get_isocontour_threshold();
            real tIntersectionTolerance = tLockedInterfaceGeometry->get_intersection_tolerance();

            // Parent basis
            Matrix< DDRMat > tFirstParentBasisValues;

            tInterpolation->eval_N( aFirstParentNodeLocalCoordinates, tFirstParentBasisValues );

            Matrix< DDRMat > tFirstParentGlobalCoordinates = tFirstParentBasisValues( 0 ) * mAncestorNodeCoordinates( 0 );

            for ( uint tBasisIndex = 1; tBasisIndex < mBasisValues.length(); tBasisIndex++ )
            {
                tFirstParentGlobalCoordinates += tFirstParentBasisValues( tBasisIndex ) * mAncestorNodeCoordinates( tBasisIndex );
            }

            // First parent on interface
            real tFirstParentPhi = tLockedInterfaceGeometry->get_field_value( mFirstParentNodeIndex, tFirstParentGlobalCoordinates );

            mFirstDiffFromThreshold = tFirstParentPhi - tIsocontourThreshold;

            bool tFirstParentOnInterface = std::abs( mFirstDiffFromThreshold ) < tIntersectionTolerance
                                        or 0.5 * norm( mParentVector ) * std::abs( 1 + mLocalCoordinate ) < tIntersectionTolerance;

            // delete interpolator
            delete tInterpolation;

            return tFirstParentOnInterface;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool
        Intersection_Node_Linear::determine_second_parent_on_interface( 
            const Element_Intersection_Type aAncestorBasisFunction,
            const Matrix< DDRMat >& aSecondParentNodeLocalCoordinates )
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

            // lock the interface geometry
            std::shared_ptr< Geometry > tLockedInterfaceGeometry = mInterfaceGeometry.lock();

            // get the isocontour thresholds from the geometry
            real tIsocontourThreshold   = tLockedInterfaceGeometry->get_isocontour_threshold();
            real tIntersectionTolerance = tLockedInterfaceGeometry->get_intersection_tolerance();

            // Parent basis
            Matrix< DDRMat > tSecondParentBasisValues;

            tInterpolation->eval_N( aSecondParentNodeLocalCoordinates, tSecondParentBasisValues );

            Matrix< DDRMat > tSecondParentGlobalCoordinates = tSecondParentBasisValues( 0 ) * mAncestorNodeCoordinates( 0 );

            for ( uint tBasisIndex = 1; tBasisIndex < mBasisValues.length(); tBasisIndex++ )
            {
                tSecondParentGlobalCoordinates += tSecondParentBasisValues( tBasisIndex ) * mAncestorNodeCoordinates( tBasisIndex );
            }

            // Second parent on interface
            real tSecondParentPhi = tLockedInterfaceGeometry->get_field_value( mSecondParentNodeIndex, tSecondParentGlobalCoordinates );

            mSecondDiffFromThreshold = tSecondParentPhi - tIsocontourThreshold;

            bool tSecondParentOnInterface = std::abs( mSecondDiffFromThreshold ) < tIntersectionTolerance
                                    or 0.5 * norm( mParentVector ) * std::abs( 1 - mLocalCoordinate ) < tIntersectionTolerance;
            
            // delete interpolator
            delete tInterpolation;

            return tSecondParentOnInterface;
        }

        //--------------------------------------------------------------------------------------------------------------

        bool
        Intersection_Node_Linear::determine_is_intersected()
        {
            bool tIsIntersected;
            
            // Determine if edge is intersected
            if ( mFirstParentOnInterface or mSecondParentOnInterface )
            {
                tIsIntersected = true;
            }
            // FIXME: This check should be unnecessary as the local edge coordinate should be sufficient
            // to determine whether edge is intersected; it is only "useful" if parent node's level set value
            // is determined by method that is different from intersection nodes; for example level set value child node
            // of child node is computed via analytic field and intersection node via bi-linear interpolation
            else if ( mFirstDiffFromThreshold * mSecondDiffFromThreshold > 0 )
            {
                tIsIntersected = false;

                // check for consistency of parent values and local coordinate
                MORIS_ASSERT( std::abs( mLocalCoordinate ) > 1,
                        "Intersection_Node::Intersection_Node - inconsistent parent level set values versus local coordinate - p1 %e p2 %e loc %e.",
                        mFirstDiffFromThreshold,
                        mSecondDiffFromThreshold,
                        mLocalCoordinate );
            }
            else
            {
                tIsIntersected = ( std::abs( mLocalCoordinate ) <= 1.0 );

                // check for consistency with parent values
                // this check is currently useless but should be performed is inconsistency issue (see comment above) is resolved
                MORIS_ASSERT( tIsIntersected ? mFirstDiffFromThreshold * mSecondDiffFromThreshold < 0 : mFirstDiffFromThreshold * mSecondDiffFromThreshold > 0,
                        "Intersection_Node::Intersection_Node - inconsistent parent level set values - p1 %e p2 %e loc %e.",
                        mFirstDiffFromThreshold,
                        mSecondDiffFromThreshold,
                        mLocalCoordinate );
            }

            return tIsIntersected;
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
        Intersection_Node_Linear::compute_local_coordinate(
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

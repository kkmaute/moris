/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Space_Interpolator.cpp
 *
 */

#include "fn_norm.hpp"
#include "fn_cross.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"
#include "fn_inv.hpp"
#include "op_div.hpp"
#include "fn_linsolve.hpp"

#include "cl_MTK_Space_Interpolator.hpp"

namespace moris
{
    namespace mtk
    {
        // smallest acceptable value for DetJ
        // note: should be consistent with Geometry_Interpolator::sDetJLowerLimit
        const real Space_Interpolator::sDetJLowerLimit = -1.0e-6;

        // smallest acceptable value for DetJ used in building inverse of Jacobian
        // note: should be consistent with Geometry_Interpolator::sDetJInvJacLowerLimit
        const real Space_Interpolator::sDetJInvJacLowerLimit = 1.0e-24;

        //------------------------------------------------------------------------------

        Space_Interpolator::Space_Interpolator(
                const Interpolation_Rule& aInterpolationRule,
                const CellShape&          aInterpolationShape,
                const bool                aSpaceSideset )
        {
            // set bool for side interpolation to true
            mSpaceSideset = aSpaceSideset;

            // create member pointer to space interpolation function
            mSpaceInterpolation = aInterpolationRule.create_space_interpolation_function();

            // number of space bases and dimensions
            mNumSpaceBases    = mSpaceInterpolation->get_number_of_bases();
            mNumSpaceDim      = mSpaceInterpolation->get_number_of_dimensions();
            mNumSpaceParamDim = mSpaceInterpolation->get_number_of_param_dimensions();

            // set member geometry type
            mGeometryType = aInterpolationRule.get_geometry_type();

            // set the interpolation shape
            mInterpolationShape = aInterpolationShape;

            // Assuming the interpolation cell geometry is the same as the interpolation rule geometry
            mIPMappingGeometryType     = mGeometryType;
            mIPMappingNumSpaceParamDim = mNumSpaceParamDim;

            // set pointers for second derivative depending on space and time dimensions
            this->set_function_pointers();
        }

        //------------------------------------------------------------------------------

        Space_Interpolator::Space_Interpolator(
                const Interpolation_Rule& aInterpolationRule,
                const Interpolation_Rule& aIPMapInterpolationRule,
                const CellShape&          aInterpolationShape,
                const bool                aSpaceSideset )
        {
            // set bool for side interpolation to true
            mSpaceSideset = aSpaceSideset;

            // create member pointer to space interpolation function
            mSpaceInterpolation = aInterpolationRule.create_space_interpolation_function();

            // number of space bases and dimensions
            mNumSpaceBases    = mSpaceInterpolation->get_number_of_bases();
            mNumSpaceDim      = mSpaceInterpolation->get_number_of_dimensions();
            mNumSpaceParamDim = mSpaceInterpolation->get_number_of_param_dimensions();

            // set member geometry type
            mGeometryType = aInterpolationRule.get_geometry_type();

            // set the interpolation shape
            mInterpolationShape = aInterpolationShape;

            // Interpolation cell geometry type and space param  dim.  This will be used
            // to determine an appropriate mapping size.
            mIPMappingGeometryType = aIPMapInterpolationRule.get_geometry_type();

            // getting param dimensions
            mIPMappingNumSpaceParamDim = aIPMapInterpolationRule.get_number_of_param_dimensions();

            // set pointers for second derivative depending on space and time dimensions
            this->set_function_pointers();
        }

        //------------------------------------------------------------------------------

        Space_Interpolator::~Space_Interpolator()
        {
            // delete interpolation functions
            if ( mSpaceInterpolation != NULL )
            {
                delete mSpaceInterpolation;
            }
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::reset_eval_flags()
        {
            // reset booleans for evaluation
            mValxEval = true;

            mNXiEval     = true;
            mdNdXiEval   = true;
            md2NdXi2Eval = true;
            md3NdXi3Eval = true;

            mSpaceDetJEval      = true;
            mSpaceJacEval       = true;
            mInvSpaceJacEval    = true;
            mSpaceJacDerivEval  = true;
            mSpaceDetJDerivEval = true;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::reset_eval_flags_coordinates()
        {
            mValxEval = true;

            mSpaceDetJEval      = true;
            mSpaceJacEval       = true;
            mInvSpaceJacEval    = true;
            mSpaceJacDerivEval  = true;
            mSpaceDetJDerivEval = true;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::reset_eval_flags_deriv()
        {
            mSpaceJacDerivEval  = true;
            mSpaceDetJDerivEval = true;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::set_space_coeff( const Matrix< DDRMat >& aXHat )
        {
            // check the space coefficients input size
            //  fixme can not check the number of cols for aXHat
            MORIS_ASSERT( aXHat.n_rows() == mNumSpaceBases,
                    " Space_Interpolator::set_space_coeff - Wrong input size (aXHat). " );

            // set the space coefficients
            mXHat = aXHat;

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::set_param_coeff()
        {
            // default implementation
            // set space and time param coords
            mSpaceInterpolation->get_param_coords( mXiHat );
            mXiHat = trans( mXiHat );

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::set_space_param_coeff( const Matrix< DDRMat >& aXiHat )
        {
            // check the space param coefficients input size
            //  fixme can not check the number of cols for aXiHat

            MORIS_ASSERT( aXiHat.n_rows() == mNumSpaceBases,
                    " Space_Interpolator::set_space_param_coeff - Wrong input size (aXiHat). %-5i vs %-5i ",
                    aXiHat.n_rows(),
                    mNumSpaceBases );

            // set the space coefficients
            mXiHat = aXiHat;

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::set_space_time( const Matrix< DDRMat >& aParamPoint )
        {

            // set input values
            mXiLocal = aParamPoint( { 0, mNumSpaceParamDim - 1 }, { 0, 0 } );

            // if no mapping required
            if ( !mMapFlag )
            {
                mMappedPoint = aParamPoint;
            }

            // reset bool for evaluation
            this->reset_eval_flags();
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::set_space( const Matrix< DDRMat >& aSpaceParamPoint )
        {

            // check input size aParamPoint
            MORIS_ASSERT( ( ( aSpaceParamPoint.n_cols() == 1 ) && ( aSpaceParamPoint.n_rows() == mNumSpaceParamDim ) ),
                    "Space_Interpolator::set_space - Wrong input size ( aSpaceParamPoint )." );

            // check input values are between -1 and 1
            // fixme what about TRI and TET
            for ( uint Ik = 0; Ik < mNumSpaceParamDim; Ik++ )
            {
                MORIS_ASSERT( ( ( aSpaceParamPoint( Ik ) <= 1.0 + mEpsilon ) && ( aSpaceParamPoint( Ik ) >= -1.0 - mEpsilon ) ),
                        "Space_Interpolator::set_space - Wrong input value ( aSpaceParamPoint )." );
            }

            // set input values
            mXiLocal = aSpaceParamPoint;

            // if no mapping required
            if ( !mMapFlag )
            {
                mMappedPoint( { 0, mNumSpaceParamDim - 1 }, { 0, 0 } ) =
                        aSpaceParamPoint.matrix_data();
            }

            // reset bool for evaluation
            this->reset_eval_flags();
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Space_Interpolator::NXi()
        {
            // if shape functions need to be evaluated
            if ( mNXiEval )
            {
                // evaluate the shape functions
                this->eval_NXi();

                // set bool for evaluation
                mNXiEval = false;
            }

            // return member value
            return mNXi;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_NXi()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Space_Interpolator::eval_NXi - mXiLocal is not set." );

            // pass data through interpolation function
            mSpaceInterpolation->eval_N( mXiLocal, mNXi );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Space_Interpolator::dNdXi()
        {
            // if shape functions need to be evaluated
            if ( mdNdXiEval )
            {
                // evaluate the shape functions 1st derivative
                this->eval_dNdXi();

                // set bool for evaluation
                mdNdXiEval = false;
            }

            // return member value
            return mdNdXi;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_dNdXi()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Space_Interpolator::eval_dNdXi - mXiLocal is not set." );

            // pass data through interpolation function
            mSpaceInterpolation->eval_dNdXi( mXiLocal, mdNdXi );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Space_Interpolator::d2NdXi2()
        {
            // if shape functions need to be evaluated
            if ( md2NdXi2Eval )
            {
                // evaluate the shape functions 2nd derivative
                this->eval_d2NdXi2();

                // set bool for evaluation
                md2NdXi2Eval = false;
            }

            // return member value
            return md2NdXi2;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_d2NdXi2()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Space_Interpolator::eval_d2NdXi2 - mXiLocal is not set." );

            // pass data through interpolation function
            mSpaceInterpolation->eval_d2NdXi2( mXiLocal, md2NdXi2 );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Space_Interpolator::d3NdXi3()
        {
            // if shape functions need to be evaluated
            if ( md3NdXi3Eval )
            {
                // evaluate the shape functions 3rd derivative
                this->eval_d3NdXi3();

                // set bool for evaluation
                md3NdXi3Eval = false;
            }

            // return member value
            return md3NdXi3;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_d3NdXi3()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Space_Interpolator::eval_d3NdXi3 - mXiLocal is not set." );

            // pass data through interpolation function
            mSpaceInterpolation->eval_d3NdXi3( mXiLocal, md3NdXi3 );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Space_Interpolator::space_jacobian()
        {
            // if space Jacobian needs to be evaluated
            if ( mSpaceJacEval )
            {
                // evaluate the space Jacobian
                this->eval_space_jacobian();

                // set bool for evaluation
                mSpaceJacEval = false;
            }

            // return member value
            return mSpaceJac;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_space_jacobian()
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel() > 0,
                    "Space_Interpolator::space_jacobian - mXHat is not set." );

            // compute the Jacobian
            mSpaceJac = this->dNdXi() * mXHat;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Space_Interpolator::space_jacobian_deriv(
                const uint& aLocalVertexID,
                const uint& aDirection )
        {
            // if space Jacobian needs to be evaluated
            if ( mSpaceJacDerivEval )
            {
                // evaluate the space Jacobian
                this->eval_space_jacobian_deriv( aLocalVertexID, aDirection );

                // set bool for evaluation
                mSpaceJacDerivEval = false;
            }

            // return member value
            return mSpaceJacDeriv;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_space_jacobian_deriv(
                const uint& aLocalVertexID,
                const uint& aDirection )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel() > 0,
                    "Space_Interpolator::eval_space_jacobian_deriv - mXHat is not set." );

            // check inputs wrt to xHat
            MORIS_ASSERT( aDirection < mXHat.n_cols(),
                    "Space_Interpolator::eval_space_jacobian_deriv - invalid direction." );
            MORIS_ASSERT( aLocalVertexID < mXHat.n_rows(),
                    "Space_Interpolator::eval_space_jacobian_deriv - invalid vertex ID." );

            // get derivative of space coefficient. Note that only one element of this matrix will have a value
            Matrix< DDRMat > tXHatDeriv = mXHat;
            Matrix< DDRMat > tZeroVector( mXHat.n_rows(), 1, 0.0 );

            // set all other elements in the derivative direction to 0
            tXHatDeriv( { 0, mXHat.n_rows() - 1 }, { aDirection, aDirection } ) = tZeroVector( { 0, mXHat.n_rows() - 1 }, { 0, 0 } );

            // what dof are we taking a derivative wrt? this will be the only used value in that column
            tXHatDeriv( aLocalVertexID, aDirection ) = 1.0;

            // compute the Jacobian derivative
            mSpaceJacDeriv = this->dNdXi() * tXHatDeriv;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Space_Interpolator::inverse_space_jacobian()
        {
            // if inverse of the space Jacobian needs to be evaluated
            if ( mInvSpaceJacEval )
            {
                // evaluate the inverse of the space Jacobian
                this->eval_inverse_space_jacobian();

                // set bool for evaluation
                mInvSpaceJacEval = false;
            }

            // return member value
            return mInvSpaceJac;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_inverse_space_jacobian()
        {
            // compute the standard inv Jacobian
            ( this->*mInvSpaceJacFunc )();
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_inverse_space_jacobian_1d()
        {
            // get the space Jacobian
            const Matrix< DDRMat >& tSpacJac = this->space_jacobian();

            MORIS_ASSERT( tSpacJac( 0, 0 ) > sDetJInvJacLowerLimit,
                    "Space determinate (1D) close to zero or negative: %e\n",
                    tSpacJac( 0, 0 ) );

            mInvSpaceJac.set_size( 1, 1 );

            mInvSpaceJac( 0, 0 ) = 1.0 / tSpacJac( 0, 0 );

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvSpaceJac - inv( tSpacJac ) ) < 1e-8 * norm( mInvSpaceJac ),
                    "Inconsistent space Jacobian (1D)\n" );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_inverse_space_jacobian_2d()
        {
            // get the space Jacobian
            const Matrix< DDRMat >& tSpacJac = this->space_jacobian();

            MORIS_ASSERT( this->space_det_J() > sDetJInvJacLowerLimit,
                    "Space determinate (2D) close to zero or negative: %e\n",
                    this->space_det_J() );

            // compute inverse of 3x3 matrix
            real tInvDet = 1.0 / ( this->space_det_J() );

            // compute inverse
            mInvSpaceJac.set_size( 2, 2 );

            mInvSpaceJac( 0, 0 ) = tSpacJac( 1, 1 ) * tInvDet;
            mInvSpaceJac( 0, 1 ) = -tSpacJac( 0, 1 ) * tInvDet;
            mInvSpaceJac( 1, 0 ) = -tSpacJac( 1, 0 ) * tInvDet;
            mInvSpaceJac( 1, 1 ) = tSpacJac( 0, 0 ) * tInvDet;

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvSpaceJac - inv( tSpacJac ) ) < 1e-8 * norm( mInvSpaceJac ),
                    "Inconsistent space Jacobian (2D)" );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_inverse_space_jacobian_2d_tri()
        {
            // get the space Jacobian
            const Matrix< DDRMat >& tSpacJac = this->space_jacobian();

            MORIS_ASSERT( this->space_det_J() > sDetJInvJacLowerLimit,
                    "Space determinate (2D) close to zero or negative: %e\n",
                    this->space_det_J() );

            // compute inv det J * 1/2
            real tInvDet = 0.5 / ( this->space_det_J() );

            // compute inverse
            mInvSpaceJac.set_size( 2, 2 );

            mInvSpaceJac( 0, 0 ) = tSpacJac( 1, 1 ) * tInvDet;
            mInvSpaceJac( 0, 1 ) = -tSpacJac( 0, 1 ) * tInvDet;
            mInvSpaceJac( 1, 0 ) = -tSpacJac( 1, 0 ) * tInvDet;
            mInvSpaceJac( 1, 1 ) = tSpacJac( 0, 0 ) * tInvDet;

            // no generic checks available for this calculation;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_inverse_space_jacobian_2d_rect()
        {
            // get the space Jacobian
            const Matrix< DDRMat >& tSpacJac = this->space_jacobian();

            MORIS_ASSERT( this->space_det_J() > sDetJInvJacLowerLimit,
                    "Space determinate (2D) close to zero or negative: %e\n",
                    this->space_det_J() );

            MORIS_ASSERT(
                    std::abs( tSpacJac( 0, 1 ) ) < mEpsilon || std::abs( tSpacJac( 1, 0 ) ) < mEpsilon,
                    "Space_Interpolator::eval_inverse_space_jacobian_2d_rect - Jacobian is not diagonal" );

            // compute inverse
            mInvSpaceJac.set_size( 2, 2, 0.0 );

            // reciprocals
            mInvSpaceJac( 0, 0 ) = 1.0 / tSpacJac( 0, 0 );
            mInvSpaceJac( 1, 1 ) = 1.0 / tSpacJac( 1, 1 );

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvSpaceJac - inv( tSpacJac ) ) < 1e-8 * norm( mInvSpaceJac ),
                    "Inconsistent space Jacobian (2D)" );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_inverse_space_jacobian_3d()
        {
            // get the space Jacobian
            const Matrix< DDRMat >& tSpacJac = this->space_jacobian();

            MORIS_ASSERT( this->space_det_J() > sDetJInvJacLowerLimit,
                    "Space determinate (3D) close to zero or negative: %e\n",
                    this->space_det_J() );

            // compute inverse of 3x3 matrix
            real tInvDet = 1.0 / ( this->space_det_J() );

            // compute inverse
            mInvSpaceJac.set_size( 3, 3 );

            mInvSpaceJac( 0, 0 ) = ( tSpacJac( 1, 1 ) * tSpacJac( 2, 2 ) - tSpacJac( 2, 1 ) * tSpacJac( 1, 2 ) ) * tInvDet;
            mInvSpaceJac( 0, 1 ) = ( tSpacJac( 0, 2 ) * tSpacJac( 2, 1 ) - tSpacJac( 0, 1 ) * tSpacJac( 2, 2 ) ) * tInvDet;
            mInvSpaceJac( 0, 2 ) = ( tSpacJac( 0, 1 ) * tSpacJac( 1, 2 ) - tSpacJac( 0, 2 ) * tSpacJac( 1, 1 ) ) * tInvDet;
            mInvSpaceJac( 1, 0 ) = ( tSpacJac( 1, 2 ) * tSpacJac( 2, 0 ) - tSpacJac( 1, 0 ) * tSpacJac( 2, 2 ) ) * tInvDet;
            mInvSpaceJac( 1, 1 ) = ( tSpacJac( 0, 0 ) * tSpacJac( 2, 2 ) - tSpacJac( 0, 2 ) * tSpacJac( 2, 0 ) ) * tInvDet;
            mInvSpaceJac( 1, 2 ) = ( tSpacJac( 1, 0 ) * tSpacJac( 0, 2 ) - tSpacJac( 0, 0 ) * tSpacJac( 1, 2 ) ) * tInvDet;
            mInvSpaceJac( 2, 0 ) = ( tSpacJac( 1, 0 ) * tSpacJac( 2, 1 ) - tSpacJac( 2, 0 ) * tSpacJac( 1, 1 ) ) * tInvDet;
            mInvSpaceJac( 2, 1 ) = ( tSpacJac( 2, 0 ) * tSpacJac( 0, 1 ) - tSpacJac( 0, 0 ) * tSpacJac( 2, 1 ) ) * tInvDet;
            mInvSpaceJac( 2, 2 ) = ( tSpacJac( 0, 0 ) * tSpacJac( 1, 1 ) - tSpacJac( 1, 0 ) * tSpacJac( 0, 1 ) ) * tInvDet;

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvSpaceJac - inv( tSpacJac ) ) < 1e-8 * norm( mInvSpaceJac ),
                    "Inconsistent space Jacobian (3D)" );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_inverse_space_jacobian_3d_rect()
        {
            // get the space Jacobian
            const Matrix< DDRMat >& tSpacJac = this->space_jacobian();

            MORIS_ASSERT( this->space_det_J() > sDetJInvJacLowerLimit,
                    "Space determinate (3D) close to zero or negative: %e\n",
                    this->space_det_J() );

            MORIS_ASSERT(
                    std::abs( tSpacJac( 0, 1 ) ) < mEpsilon ||            //
                            std::abs( tSpacJac( 0, 2 ) ) < mEpsilon ||    //
                            std::abs( tSpacJac( 1, 0 ) ) < mEpsilon ||    //
                            std::abs( tSpacJac( 1, 2 ) ) < mEpsilon ||    //
                            std::abs( tSpacJac( 2, 0 ) ) < mEpsilon ||    //
                            std::abs( tSpacJac( 2, 1 ) ) < mEpsilon,
                    "Space_Interpolator::eval_inverse_space_jacobian_3d_rect - Jacobian is not diagonal" );

            // compute inverse
            mInvSpaceJac.set_size( 3, 3, 0.0 );

            // reciprocals
            mInvSpaceJac( 0, 0 ) = 1.0 / tSpacJac( 0, 0 );
            mInvSpaceJac( 1, 1 ) = 1.0 / tSpacJac( 1, 1 );
            mInvSpaceJac( 2, 2 ) = 1.0 / tSpacJac( 2, 2 );

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvSpaceJac - inv( tSpacJac ) ) < 1e-8 * norm( mInvSpaceJac ),
                    "Inconsistent space Jacobian (3D)" );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::second_space_jacobian( Matrix< DDRMat >& aJ2bt )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel() > 0,
                    "Space_Interpolator::second_space_jacobian - mXHat is not set." );

            // compute the second order Jacobian
            aJ2bt = this->d2NdXi2() * mXHat;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::third_space_jacobian( Matrix< DDRMat >& aJ3ct )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel() > 0,
                    "Space_Interpolator::third_space_jacobian - mXHat is not set." );

            // compute the third order Jacobian
            aJ3ct = this->d3NdXi3() * mXHat;
        }

        //------------------------------------------------------------------------------

        const real&
        Space_Interpolator::space_det_J()
        {
            // if determinant of space Jacobian needs to be evaluated
            if ( mSpaceDetJEval )
            {
                // get the space Jacobian
                const Matrix< DDRMat >& tSpaceJt = this->space_jacobian();

                // det j function pointer
                mSpaceDetJ = ( this->*mSpaceDetJFunc )( tSpaceJt );

                // set bool for evaluation
                mSpaceDetJEval = false;
            }

            // return member value
            return mSpaceDetJ;
        }

        //------------------------------------------------------------------------------

        const real&
        Space_Interpolator::space_det_J_deriv(
                const uint& aLocalVertexID,
                const uint& aDirection )
        {
            // if determinant of space Jacobian derivative needs to be evaluated
            if ( mSpaceDetJDerivEval )
            {
                // get the space Jacobian
                const Matrix< DDRMat >& tSpaceJtDeriv = this->space_jacobian_deriv( aLocalVertexID, aDirection );

                //  det J deriv function pointer
                mSpaceDetJDeriv = ( this->*mSpaceDetJDerivFunc )( tSpaceJtDeriv );

                // set bool for evaluation
                mSpaceDetJDerivEval = false;
            }

            // return member value
            return mSpaceDetJDeriv;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_side_line(
                const Matrix< DDRMat >& aSpaceJt )
        {
            real tDetJ = norm( aSpaceJt );

            MORIS_ASSERT( tDetJ > sDetJLowerLimit,
                    "Space determinant (side line) close to zero or negative: %e\n",
                    tDetJ );

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_deriv_side_line(
                const Matrix< DDRMat >& aSpaceJDerivt )
        {
            real tDetJDeriv = norm( aSpaceJDerivt );

            return tDetJDeriv;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_side_tri(
                const Matrix< DDRMat >& aSpaceJt )
        {
            MORIS_ASSERT( aSpaceJt.n_rows() == 2 && aSpaceJt.n_cols() == 3,
                    "Space_Interpolator::eval_space_detJ_side_tri - incorrect dimension of space Jacobian" );

            // note: space Jacobian contains vectors x1-x3 and x2-x3
            real tDetJ = norm( cross( aSpaceJt.get_row( 0 ), aSpaceJt.get_row( 1 ) ) ) / 2.0;

            MORIS_ASSERT( tDetJ > sDetJLowerLimit,
                    "Space_Interpolator::eval_space_detJ_side_tri - Space determinant (side tri) close to zero or negative: %e\n",
                    tDetJ );

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_deriv_side_tri(
                const Matrix< DDRMat >& aSpaceJDerivt )
        {
            MORIS_ASSERT( aSpaceJDerivt.n_rows() == 2 && aSpaceJDerivt.n_cols() == 3,
                    "Space_Interpolator::eval_space_detJ_deriv_side_tri - incorrect dimension of space Jacobian" );

            real tDetJDeriv = norm( cross( aSpaceJDerivt.get_row( 0 ), aSpaceJDerivt.get_row( 1 ) ) ) / 2.0;

            return tDetJDeriv;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_side_quad(
                const Matrix< DDRMat >& aSpaceJt )
        {
            real tDetJ = norm( cross( aSpaceJt.get_row( 0 ), aSpaceJt.get_row( 1 ) ) );

            MORIS_ASSERT( tDetJ > sDetJLowerLimit,
                    "Space determinant (side quad) close to zero or negative: %e\n",
                    tDetJ );

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_deriv_side_quad(
                const Matrix< DDRMat >& aSpaceJDerivt )
        {
            real tDetJDeriv = norm( cross( aSpaceJDerivt.get_row( 0 ), aSpaceJDerivt.get_row( 1 ) ) );

            return tDetJDeriv;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_bulk_line(
                const Matrix< DDRMat >& aSpaceJt )
        {
            real tDetJ = aSpaceJt( 0, 0 );

            MORIS_ASSERT( tDetJ > sDetJLowerLimit,
                    "Space determinant (bulk 1D) close to zero or negative: %e\n",
                    tDetJ );

            MORIS_ASSERT( std::abs( det( aSpaceJt ) - tDetJ ) < 1e-8 * tDetJ,
                    "Inconsistent space determinant (bulk 1D): %e vs %e\n",
                    tDetJ,
                    det( aSpaceJt ) );

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_deriv_bulk_line(
                const Matrix< DDRMat >& aSpaceJDerivt )
        {
            real tDetJDeriv = aSpaceJDerivt( 0, 0 );

            return tDetJDeriv;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_bulk_quad(
                const Matrix< DDRMat >& aSpaceJt )
        {
            real tDetJ = aSpaceJt( 0, 0 ) * aSpaceJt( 1, 1 ) - aSpaceJt( 0, 1 ) * aSpaceJt( 1, 0 );

            MORIS_ASSERT( tDetJ > sDetJLowerLimit,
                    "Space determinant (bulk 2D) close to zero or negative: %e\n",
                    tDetJ );

            MORIS_ASSERT( std::abs( det( aSpaceJt ) - tDetJ ) < 1e-8 * tDetJ,
                    "Inconsistent space determinant (bulk 2D): %e vs %e\n",
                    tDetJ,
                    det( aSpaceJt ) );

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_deriv_bulk_quad(
                const Matrix< DDRMat >& aSpaceJDerivt )
        {
            real tDetJDeriv = aSpaceJDerivt( 0, 0 ) * aSpaceJDerivt( 1, 1 ) -    //
                              aSpaceJDerivt( 0, 1 ) * aSpaceJDerivt( 1, 0 );

            return tDetJDeriv;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_bulk_quad_rect(
                const Matrix< DDRMat >& aSpaceJt )
        {
            MORIS_ASSERT(
                    std::abs( aSpaceJt( 0, 1 ) ) < mEpsilon || std::abs( aSpaceJt( 1, 0 ) ) < mEpsilon,
                    "Space_Interpolator::eval_space_detJ_bulk_quad_rect - Jacobian is not diagonal" );

            real tDetJ = aSpaceJt( 0, 0 ) * aSpaceJt( 1, 1 );

            MORIS_ASSERT( tDetJ > sDetJLowerLimit,
                    "Space determinant (bulk 2D) close to zero or negative: %e\n",
                    tDetJ );

            MORIS_ASSERT( std::abs( det( aSpaceJt ) - tDetJ ) < 1e-8 * tDetJ,
                    "Inconsistent space determinant (bulk 2D): %e vs %e\n",
                    tDetJ,
                    det( aSpaceJt ) );

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_deriv_bulk_quad_rect(
                const Matrix< DDRMat >& aSpaceJDerivt )
        {
            MORIS_ASSERT(
                    std::abs( aSpaceJDerivt( 0, 1 ) ) < mEpsilon || std::abs( aSpaceJDerivt( 1, 0 ) ) < mEpsilon,
                    "Space_Interpolator::eval_space_detJ_deriv_bulk_quad_rect - Jacobian is not diagonal" );

            real tDetJDeriv = aSpaceJDerivt( 0, 0 ) * aSpaceJDerivt( 1, 1 );

            return tDetJDeriv;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_bulk_hex(
                const Matrix< DDRMat >& aSpaceJt )
        {
            real tDetJ =
                    +aSpaceJt( 0, 0 ) * ( aSpaceJt( 1, 1 ) * aSpaceJt( 2, 2 ) - aSpaceJt( 2, 1 ) * aSpaceJt( 1, 2 ) )
                    - aSpaceJt( 0, 1 ) * ( aSpaceJt( 1, 0 ) * aSpaceJt( 2, 2 ) - aSpaceJt( 1, 2 ) * aSpaceJt( 2, 0 ) )
                    + aSpaceJt( 0, 2 ) * ( aSpaceJt( 1, 0 ) * aSpaceJt( 2, 1 ) - aSpaceJt( 1, 1 ) * aSpaceJt( 2, 0 ) );

            MORIS_ASSERT( tDetJ > sDetJLowerLimit,
                    "Space determinant (bulk 3D) close to zero or negative: %e\n",
                    tDetJ );

            MORIS_ASSERT( std::abs( det( aSpaceJt ) - tDetJ ) < 1e-8 * tDetJ,
                    "Inconsistent space determinant (bulk 3D): %e vs %e\n",
                    tDetJ,
                    det( aSpaceJt ) );

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_deriv_bulk_hex(
                const Matrix< DDRMat >& aSpaceJDerivt )
        {
            real tDetJDeriv =
                    +aSpaceJDerivt( 0, 0 ) * ( aSpaceJDerivt( 1, 1 ) * aSpaceJDerivt( 2, 2 ) - aSpaceJDerivt( 2, 1 ) * aSpaceJDerivt( 1, 2 ) )
                    - aSpaceJDerivt( 0, 1 ) * ( aSpaceJDerivt( 1, 0 ) * aSpaceJDerivt( 2, 2 ) - aSpaceJDerivt( 1, 2 ) * aSpaceJDerivt( 2, 0 ) )
                    + aSpaceJDerivt( 0, 2 ) * ( aSpaceJDerivt( 1, 0 ) * aSpaceJDerivt( 2, 1 ) - aSpaceJDerivt( 1, 1 ) * aSpaceJDerivt( 2, 0 ) );

            return tDetJDeriv;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_bulk_hex_rect(
                const Matrix< DDRMat >& aSpaceJt )
        {
            MORIS_ASSERT(
                    std::abs( aSpaceJt( 0, 1 ) ) < mEpsilon ||            //
                            std::abs( aSpaceJt( 0, 2 ) ) < mEpsilon ||    //
                            std::abs( aSpaceJt( 1, 0 ) ) < mEpsilon ||    //
                            std::abs( aSpaceJt( 1, 2 ) ) < mEpsilon ||    //
                            std::abs( aSpaceJt( 2, 0 ) ) < mEpsilon ||    //
                            std::abs( aSpaceJt( 2, 1 ) ) < mEpsilon,
                    "Space_Interpolator::eval_space_detJ_bulk_hex_rect - Jacobian is not diagonal" );

            // init tDet J
            real tDetJ;

            // get trace of jacobian
            tDetJ = aSpaceJt( 1, 1 ) * aSpaceJt( 0, 0 ) * aSpaceJt( 2, 2 );

            MORIS_ASSERT( tDetJ > sDetJLowerLimit,
                    "Space determinant (bulk 3D) close to zero or negative: %e\n",
                    tDetJ );

            MORIS_ASSERT( std::abs( det( aSpaceJt ) - tDetJ ) < 1e-8 * tDetJ,
                    "Inconsistent space determinant (bulk 3D): %e vs %e\n",
                    tDetJ,
                    det( aSpaceJt ) );

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_deriv_bulk_hex_rect(
                const Matrix< DDRMat >& aSpaceJDerivt )
        {
            MORIS_ASSERT(
                    std::abs( aSpaceJDerivt( 0, 1 ) ) < mEpsilon ||            //
                            std::abs( aSpaceJDerivt( 0, 2 ) ) < mEpsilon ||    //
                            std::abs( aSpaceJDerivt( 1, 0 ) ) < mEpsilon ||    //
                            std::abs( aSpaceJDerivt( 1, 2 ) ) < mEpsilon ||    //
                            std::abs( aSpaceJDerivt( 2, 0 ) ) < mEpsilon ||    //
                            std::abs( aSpaceJDerivt( 2, 1 ) ) < mEpsilon,
                    "Space_Interpolator::eval_space_detJ_deriv_bulk_hex_rect - Jacobian is not diagonal" );

            real tDetJDeriv = aSpaceJDerivt( 0, 0 ) * aSpaceJDerivt( 1, 1 ) * aSpaceJDerivt( 2, 2 );

            return tDetJDeriv;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_bulk_tri_param_2(
                const Matrix< DDRMat >& aSpaceJt )
        {
            real tDetJ = ( aSpaceJt( 0, 0 ) * aSpaceJt( 1, 1 ) - aSpaceJt( 0, 1 ) * aSpaceJt( 1, 0 ) ) / 2.0;

            MORIS_ASSERT( tDetJ > sDetJLowerLimit,
                    "Space determinant (Tri-P2) close to zero or negative: %e\n",
                    tDetJ );

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_deriv_bulk_tri_param_2(
                const Matrix< DDRMat >& aSpaceJDerivt )
        {
            real tDetJDeriv = ( aSpaceJDerivt( 0, 0 ) * aSpaceJDerivt( 1, 1 ) - aSpaceJDerivt( 0, 1 ) * aSpaceJDerivt( 1, 0 ) ) / 2.0;

            return tDetJDeriv;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_bulk_tri_param_3(
                const Matrix< DDRMat >& aSpaceJt )
        {
            real tDetJ = (    //
                                 +( aSpaceJt( 1, 0 ) * aSpaceJt( 2, 1 ) - aSpaceJt( 1, 1 ) * aSpaceJt( 2, 0 ) )
                                 - ( aSpaceJt( 0, 0 ) * aSpaceJt( 2, 1 ) - aSpaceJt( 0, 1 ) * aSpaceJt( 2, 0 ) )
                                 + ( aSpaceJt( 0, 0 ) * aSpaceJt( 1, 1 ) - aSpaceJt( 0, 1 ) * aSpaceJt( 1, 0 ) ) )
                       / 2.0;

            MORIS_ASSERT( tDetJ > sDetJLowerLimit,
                    "Space determinant (Tri-P3) close to zero or negative: %e\n",
                    tDetJ );

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_deriv_bulk_tri_param_3(
                const Matrix< DDRMat >& aSpaceJDerivt )
        {
            real tDetJDeriv = (    //
                                      +( aSpaceJDerivt( 1, 0 ) * aSpaceJDerivt( 2, 1 ) - aSpaceJDerivt( 1, 1 ) * aSpaceJDerivt( 2, 0 ) )
                                      - ( aSpaceJDerivt( 0, 0 ) * aSpaceJDerivt( 2, 1 ) - aSpaceJDerivt( 0, 1 ) * aSpaceJDerivt( 2, 0 ) )
                                      + ( aSpaceJDerivt( 0, 0 ) * aSpaceJDerivt( 1, 1 ) - aSpaceJDerivt( 0, 1 ) * aSpaceJDerivt( 1, 0 ) ) )
                            / 2.0;

            return tDetJDeriv;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_bulk_tet_param_3(
                const Matrix< DDRMat >& aSpaceJt )
        {
            real tDetJ = (    //
                                 +aSpaceJt( 0, 0 ) * ( aSpaceJt( 1, 1 ) * aSpaceJt( 2, 2 ) - aSpaceJt( 2, 1 ) * aSpaceJt( 1, 2 ) )
                                 - aSpaceJt( 0, 1 ) * ( aSpaceJt( 1, 0 ) * aSpaceJt( 2, 2 ) - aSpaceJt( 1, 2 ) * aSpaceJt( 2, 0 ) )
                                 + aSpaceJt( 0, 2 ) * ( aSpaceJt( 1, 0 ) * aSpaceJt( 2, 1 ) - aSpaceJt( 1, 1 ) * aSpaceJt( 2, 0 ) ) )
                       / 6.0;

            MORIS_ASSERT( tDetJ > sDetJLowerLimit,
                    "Space determinant (Tet-P3) close to zero or negative: %e\n",
                    tDetJ );

            MORIS_ASSERT( std::abs( det( aSpaceJt ) - tDetJ ) < 1e-8 * tDetJ,
                    "Inconsistent space determinant (Tet-P3): %e vs %e\n",
                    tDetJ,
                    det( aSpaceJt ) );

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_deriv_bulk_tet_param_3(
                const Matrix< DDRMat >& aSpaceJDerivt )
        {
            real tDetJDeriv = (    //
                                      +aSpaceJDerivt( 0, 0 ) * ( aSpaceJDerivt( 1, 1 ) * aSpaceJDerivt( 2, 2 ) - aSpaceJDerivt( 2, 1 ) * aSpaceJDerivt( 1, 2 ) )
                                      - aSpaceJDerivt( 0, 1 ) * ( aSpaceJDerivt( 1, 0 ) * aSpaceJDerivt( 2, 2 ) - aSpaceJDerivt( 1, 2 ) * aSpaceJDerivt( 2, 0 ) )
                                      + aSpaceJDerivt( 0, 2 ) * ( aSpaceJDerivt( 1, 0 ) * aSpaceJDerivt( 2, 1 ) - aSpaceJDerivt( 1, 1 ) * aSpaceJDerivt( 2, 0 ) ) )
                            / 6.0;

            return tDetJDeriv;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_bulk_tet_param_4(
                const Matrix< DDRMat >& aSpaceJt )
        {
            real tSubDet1 =
                    +aSpaceJt( 1, 0 ) * ( aSpaceJt( 2, 1 ) * aSpaceJt( 3, 2 ) - aSpaceJt( 2, 2 ) * aSpaceJt( 3, 1 ) )
                    - aSpaceJt( 1, 1 ) * ( aSpaceJt( 2, 0 ) * aSpaceJt( 3, 2 ) - aSpaceJt( 2, 2 ) * aSpaceJt( 3, 0 ) )
                    + aSpaceJt( 1, 2 ) * ( aSpaceJt( 2, 0 ) * aSpaceJt( 3, 1 ) - aSpaceJt( 2, 1 ) * aSpaceJt( 3, 0 ) );

            real tSubDet2 =
                    +aSpaceJt( 0, 0 ) * ( aSpaceJt( 2, 1 ) * aSpaceJt( 3, 2 ) - aSpaceJt( 2, 2 ) * aSpaceJt( 3, 1 ) )
                    - aSpaceJt( 0, 1 ) * ( aSpaceJt( 2, 0 ) * aSpaceJt( 3, 2 ) - aSpaceJt( 2, 2 ) * aSpaceJt( 3, 0 ) )
                    + aSpaceJt( 0, 2 ) * ( aSpaceJt( 2, 0 ) * aSpaceJt( 3, 1 ) - aSpaceJt( 2, 1 ) * aSpaceJt( 3, 0 ) );

            real tSubDet3 =
                    +aSpaceJt( 0, 0 ) * ( aSpaceJt( 1, 1 ) * aSpaceJt( 3, 2 ) - aSpaceJt( 1, 2 ) * aSpaceJt( 3, 1 ) )
                    - aSpaceJt( 0, 1 ) * ( aSpaceJt( 1, 0 ) * aSpaceJt( 3, 2 ) - aSpaceJt( 1, 2 ) * aSpaceJt( 3, 0 ) )
                    + aSpaceJt( 0, 2 ) * ( aSpaceJt( 1, 0 ) * aSpaceJt( 3, 1 ) - aSpaceJt( 1, 1 ) * aSpaceJt( 3, 0 ) );

            real tSubDet4 =
                    +aSpaceJt( 0, 0 ) * ( aSpaceJt( 1, 1 ) * aSpaceJt( 2, 2 ) - aSpaceJt( 1, 2 ) * aSpaceJt( 2, 1 ) )
                    - aSpaceJt( 0, 1 ) * ( aSpaceJt( 1, 0 ) * aSpaceJt( 2, 2 ) - aSpaceJt( 1, 2 ) * aSpaceJt( 2, 0 ) )
                    + aSpaceJt( 0, 2 ) * ( aSpaceJt( 1, 0 ) * aSpaceJt( 2, 1 ) - aSpaceJt( 1, 1 ) * aSpaceJt( 2, 0 ) );

            real tDetJ = ( tSubDet1 - tSubDet2 + tSubDet3 - tSubDet4 ) / 6.0;

            MORIS_ASSERT( tDetJ > sDetJLowerLimit,
                    "Space determinant (Tet-P4) close to zero or negative: %e\n",
                    tDetJ );

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real
        Space_Interpolator::eval_space_detJ_deriv_bulk_tet_param_4(
                const Matrix< DDRMat >& aSpaceJDerivt )
        {
            real tSubDet1 =
                    +aSpaceJDerivt( 1, 0 ) * ( aSpaceJDerivt( 2, 1 ) * aSpaceJDerivt( 3, 2 ) - aSpaceJDerivt( 2, 2 ) * aSpaceJDerivt( 3, 1 ) )
                    - aSpaceJDerivt( 1, 1 ) * ( aSpaceJDerivt( 2, 0 ) * aSpaceJDerivt( 3, 2 ) - aSpaceJDerivt( 2, 2 ) * aSpaceJDerivt( 3, 0 ) )
                    + aSpaceJDerivt( 1, 2 ) * ( aSpaceJDerivt( 2, 0 ) * aSpaceJDerivt( 3, 1 ) - aSpaceJDerivt( 2, 1 ) * aSpaceJDerivt( 3, 0 ) );

            real tSubDet2 =
                    +aSpaceJDerivt( 0, 0 ) * ( aSpaceJDerivt( 2, 1 ) * aSpaceJDerivt( 3, 2 ) - aSpaceJDerivt( 2, 2 ) * aSpaceJDerivt( 3, 1 ) )
                    - aSpaceJDerivt( 0, 1 ) * ( aSpaceJDerivt( 2, 0 ) * aSpaceJDerivt( 3, 2 ) - aSpaceJDerivt( 2, 2 ) * aSpaceJDerivt( 3, 0 ) )
                    + aSpaceJDerivt( 0, 2 ) * ( aSpaceJDerivt( 2, 0 ) * aSpaceJDerivt( 3, 1 ) - aSpaceJDerivt( 2, 1 ) * aSpaceJDerivt( 3, 0 ) );

            real tSubDet3 =
                    +aSpaceJDerivt( 0, 0 ) * ( aSpaceJDerivt( 1, 1 ) * aSpaceJDerivt( 3, 2 ) - aSpaceJDerivt( 1, 2 ) * aSpaceJDerivt( 3, 1 ) )
                    - aSpaceJDerivt( 0, 1 ) * ( aSpaceJDerivt( 1, 0 ) * aSpaceJDerivt( 3, 2 ) - aSpaceJDerivt( 1, 2 ) * aSpaceJDerivt( 3, 0 ) )
                    + aSpaceJDerivt( 0, 2 ) * ( aSpaceJDerivt( 1, 0 ) * aSpaceJDerivt( 3, 1 ) - aSpaceJDerivt( 1, 1 ) * aSpaceJDerivt( 3, 0 ) );

            real tSubDet4 =
                    +aSpaceJDerivt( 0, 0 ) * ( aSpaceJDerivt( 1, 1 ) * aSpaceJDerivt( 2, 2 ) - aSpaceJDerivt( 1, 2 ) * aSpaceJDerivt( 2, 1 ) )
                    - aSpaceJDerivt( 0, 1 ) * ( aSpaceJDerivt( 1, 0 ) * aSpaceJDerivt( 2, 2 ) - aSpaceJDerivt( 1, 2 ) * aSpaceJDerivt( 2, 0 ) )
                    + aSpaceJDerivt( 0, 2 ) * ( aSpaceJDerivt( 1, 0 ) * aSpaceJDerivt( 2, 1 ) - aSpaceJDerivt( 1, 1 ) * aSpaceJDerivt( 2, 0 ) );

            real tDetJDeriv = ( tSubDet1 - tSubDet2 + tSubDet3 - tSubDet4 ) / 6.0;

            return tDetJDeriv;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::get_normal( Matrix< DDRMat >& aNormal )
        {
            // check that there is a side interpolation
            MORIS_ASSERT( mSpaceSideset,
                    "Space_Interpolator::normal - not a side." );

            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Space_Interpolator::normal - mXiLocal is not set." );

            // evaluate side space interpolation shape functions first parametric derivatives at aParamPoint
            Matrix< DDRMat > tdNSpacedXi;
            mSpaceInterpolation->eval_dNdXi( mXiLocal, tdNSpacedXi );

            // evaluation of tangent vectors to the space side in the physical space
            Matrix< DDRMat > tRealTangents = trans( tdNSpacedXi * mXHat );

            // computing the normal from the real tangent vectors
            ( this->*mNormalFunc )( tRealTangents, aNormal );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_normal_side_line(
                Matrix< DDRMat >& aRealTangents,
                Matrix< DDRMat >& aNormal )
        {
            aNormal = {
                { aRealTangents( 1 ) },
                { -aRealTangents( 0 ) }
            };

            aNormal = aNormal / norm( aNormal );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_normal_side_quad(
                Matrix< DDRMat >& aRealTangents,
                Matrix< DDRMat >& aNormal )
        {
            aNormal = cross( aRealTangents.get_column( 0 ), aRealTangents.get_column( 1 ) );

            aNormal = aNormal / norm( aNormal );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_normal_side_tri(
                Matrix< DDRMat >& aRealTangents,
                Matrix< DDRMat >& aNormal )
        {
            aNormal = cross( aRealTangents.get_column( 0 ), aRealTangents.get_column( 1 ) );

            aNormal = aNormal / norm( aNormal );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Space_Interpolator::valx()
        {
            if ( mValxEval )
            {
                // check that mTHat is set
                MORIS_ASSERT( mXHat.numel() > 0,
                        "Space_Interpolator::valx - mXHat is not set." );

                // evaluate the field
                mValx = this->NXi() * mXHat;

                mValxEval = false;
            }

            return mValx;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Space_Interpolator::map_integration_point()
        {
            // if eval mapping
            if ( mMapFlag )
            {
                // check that mXiHat and mTauHat are set
                MORIS_ASSERT( mXiHat.numel() > 0,
                        "Space_Interpolator::eval_mapping - mXiHat is not set." );

                uint tSize = mMappedPoint.numel() - 1;

                // set mapped space coordinates
                mMappedPoint( { 0, tSize - 1 } ) = trans( this->NXi() * mXiHat );
                mMappedPoint( tSize )            = 0;
            }
            return mMappedPoint;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::update_local_coordinates(
                Matrix< DDRMat >& aPhysCoordinates,
                Matrix< DDRMat >& aParamCoordinates )
        {
            // set max iteration
            const uint tMaxIter = 5;

            // set convergence criterion
            const real tConvCrit = 1e-12;

            // Newton loop
            for ( uint iIter = 0; iIter <= tMaxIter; iIter++ )
            {
                // compute NXi
                Matrix< DDRMat > tNXi;
                mSpaceInterpolation->eval_N( aParamCoordinates, tNXi );

                // compute residual
                Matrix< DDRMat > tR = trans( aPhysCoordinates - tNXi * mXHat );

                // check for convergence
                if ( norm( tR ) < tConvCrit )
                {
                    return;
                }

                // compute dNdXI
                Matrix< DDRMat > tdNdXi;
                mSpaceInterpolation->eval_dNdXi( aParamCoordinates, tdNdXi );

                // solve
                aParamCoordinates += trans( inv( trans( tdNdXi * mXHat ) ) * tR );
            }

            // getting here means that iterations exceeded maximum number
            MORIS_ASSERT( true, "Space_Interpolator::update_local_coordinates - No convergence." );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::space_jacobian_and_matrices_for_second_derivatives(
                Matrix< DDRMat >&       aJt,
                Matrix< DDRMat >&       aKt,
                Matrix< DDRMat >&       aLt,
                const Matrix< DDRMat >& adNdXi,
                const Matrix< DDRMat >& ad2NdXi2 )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel() > 0,
                    "Space_Interpolator::space_jacobian_and_matrices_for_second_derivatives - mXHat is not set." );

            // evaluate transposed of geometry Jacobian
            aJt = this->space_jacobian();

            // call calculator for second derivatives
            this->mSecondDerivativeMatricesSpace(
                    aJt,    // contains first geometric derivs
                    aKt,
                    aLt,
                    ad2NdXi2,
                    mXHat );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::space_jacobian_and_matrices_for_third_derivatives(
                Matrix< DDRMat >&       aJt,      // contains first geometric derivs
                Matrix< DDRMat >&       aJ2bt,    // contains second geometric derivs = second help matrix for 2nd field derivs
                Matrix< DDRMat >&       aJ3at,    // first help matrix for 3rd field derivs
                Matrix< DDRMat >&       aJ3bt,    // second help matrix for 3rd field derivs
                Matrix< DDRMat >&       aJ3ct,    // third help matrix for 3rd field derivs
                const Matrix< DDRMat >& adNdXi,
                const Matrix< DDRMat >& ad2NdXi2,
                const Matrix< DDRMat >& ad3NdXi3 )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel() > 0,
                    "Space_Interpolator::space_jacobian_and_matrices_for_third_derivatives - mXHat is not set." );

            // evaluate geometry Jacobians
            aJt = this->space_jacobian();
            this->second_space_jacobian( aJ2bt );

            // call calculator for second derivatives
            this->mThirdDerivativeMatricesSpace(
                    aJt,
                    aJ2bt,
                    aJ3at,
                    aJ3bt,
                    aJ3ct,
                    ad3NdXi3,
                    mXHat );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_matrices_for_second_derivative_1d(
                const Matrix< DDRMat >& aJt,
                Matrix< DDRMat >&       aKt,
                Matrix< DDRMat >&       aLt,
                const Matrix< DDRMat >& ad2NdXi2,
                const Matrix< DDRMat >& aXHat )
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 1, 1 );
            aLt( 0, 0 ) = std::pow( aJt( 0, 0 ), 2 );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_matrices_for_third_derivative_1d(
                const Matrix< DDRMat >& aJt,
                const Matrix< DDRMat >& aJ2bt,
                Matrix< DDRMat >&       aJ3at,
                Matrix< DDRMat >&       aJ3bt,
                Matrix< DDRMat >&       aJ3ct,
                const Matrix< DDRMat >& ad3NdXi3,
                const Matrix< DDRMat >& aXHat )
        {
            // first help matrix
            aJ3at.set_size( 1, 1, std::pow( aJt( 0, 0 ), 3 ) );

            // second help matrix
            aJ3bt.set_size( 1, 1, 3 * aJ2bt( 0, 0 ) * aJt( 0, 0 ) );

            // third help matrix
            aJ3ct = ad3NdXi3 * aXHat;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_matrices_for_second_derivative_2d(
                const Matrix< DDRMat >& aJt,
                Matrix< DDRMat >&       aKt,
                Matrix< DDRMat >&       aLt,
                const Matrix< DDRMat >& ad2NdXi2,
                const Matrix< DDRMat >& aXHat )
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 3, 3 );
            aLt( 0, 0 ) = std::pow( aJt( 0, 0 ), 2 );
            aLt( 1, 0 ) = std::pow( aJt( 1, 0 ), 2 );
            aLt( 2, 0 ) = aJt( 0, 0 ) * aJt( 1, 0 );

            aLt( 0, 1 ) = std::pow( aJt( 0, 1 ), 2 );
            aLt( 1, 1 ) = std::pow( aJt( 1, 1 ), 2 );
            aLt( 2, 1 ) = aJt( 0, 1 ) * aJt( 1, 1 );

            aLt( 0, 2 ) = 2.0 * aJt( 0, 0 ) * aJt( 0, 1 );
            aLt( 1, 2 ) = 2.0 * aJt( 1, 0 ) * aJt( 1, 1 );
            aLt( 2, 2 ) = aJt( 0, 0 ) * aJt( 1, 1 ) + aJt( 0, 1 ) * aJt( 1, 0 );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_matrices_for_third_derivative_2d(
                const Matrix< DDRMat >& aJt,
                const Matrix< DDRMat >& aJ2bt,
                Matrix< DDRMat >&       aJ3at,
                Matrix< DDRMat >&       aJ3bt,
                Matrix< DDRMat >&       aJ3ct,
                const Matrix< DDRMat >& ad3NdXi3,
                const Matrix< DDRMat >& aXHat )
        {
            // first help matrix
            aJ3at.set_size( 4, 4 );

            /* matrix structured into 4 parts
             *  _____________     ________
             *  |(1)* |(2)* |     | ,xxx |
             *  |_*_*_|_*_*_|  *  | ,yyy |
             *  |(3)* |(4)* |     | ,xxy |
             *  |_*_*_|_*_*_|     |_,xyy_|
             */

            // Block (1) ------------------------------------------------
            for ( uint j = 0; j < 2; ++j )
            {
                aJ3at( 0, j ) = std::pow( aJt( 0, j ), 3.0 );
                aJ3at( 1, j ) = std::pow( aJt( 1, j ), 3.0 );
            }

            // Block (2) ------------------------------------------------
            aJ3at( 0, 2 ) = 3.0 * std::pow( aJt( 0, 0 ), 2.0 ) * aJt( 0, 1 );
            aJ3at( 1, 2 ) = 3.0 * std::pow( aJt( 1, 0 ), 2.0 ) * aJt( 1, 1 );

            aJ3at( 0, 3 ) = 3.0 * std::pow( aJt( 0, 1 ), 2.0 ) * aJt( 0, 0 );
            aJ3at( 1, 3 ) = 3.0 * std::pow( aJt( 1, 1 ), 2.0 ) * aJt( 1, 0 );

            // Block (3) ------------------------------------------------
            for ( uint j = 0; j < 2; ++j )
            {
                aJ3at( 2, j ) = std::pow( aJt( 0, j ), 2.0 ) * aJt( 1, j );
                aJ3at( 3, j ) = std::pow( aJt( 1, j ), 2.0 ) * aJt( 0, j );
            }

            // Block (4) ------------------------------------------------
            aJ3at( 2, 2 ) = std::pow( aJt( 0, 0 ), 2.0 ) * aJt( 1, 1 ) + 2 * aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 0, 1 );
            aJ3at( 3, 2 ) = std::pow( aJt( 1, 0 ), 2.0 ) * aJt( 0, 1 ) + 2 * aJt( 1, 0 ) * aJt( 0, 0 ) * aJt( 1, 1 );

            aJ3at( 2, 3 ) = std::pow( aJt( 0, 1 ), 2.0 ) * aJt( 1, 0 ) + 2 * aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 0, 0 );
            aJ3at( 3, 3 ) = std::pow( aJt( 1, 1 ), 2.0 ) * aJt( 0, 0 ) + 2 * aJt( 1, 1 ) * aJt( 0, 1 ) * aJt( 1, 0 );

            // second help matrix
            aJ3bt.set_size( 4, 3 );

            /* matrix structured into 4 parts
             *  ___________     _______
             *  |(1)* |(2)|     | ,xx |
             *  |_*_*_|_*_|  *  | ,yy |
             *  |(3)* |(4)|     |_,xy_|
             *  |_*_*_|_*_|
             */

            // Block (1) ------------------------------------------------
            for ( uint j = 0; j < 2; ++j )
            {
                aJ3bt( 0, j ) = 3.0 * aJ2bt( 0, j ) * aJt( 0, j );
                aJ3bt( 1, j ) = 3.0 * aJ2bt( 1, j ) * aJt( 1, j );
            }

            // Block (2) ------------------------------------------------
            aJ3bt( 0, 2 ) = 3.0 * aJ2bt( 0, 0 ) * aJt( 0, 1 ) + 3 * aJ2bt( 0, 1 ) * aJt( 0, 0 );
            aJ3bt( 1, 2 ) = 3.0 * aJ2bt( 1, 0 ) * aJt( 1, 1 ) + 3 * aJ2bt( 1, 1 ) * aJt( 1, 0 );

            // Block (3) ------------------------------------------------
            for ( uint j = 0; j < 2; ++j )
            {
                aJ3bt( 2, j ) = 2.0 * aJ2bt( 2, j ) * aJt( 0, j ) + aJ2bt( 0, j ) * aJt( 1, j );
                aJ3bt( 3, j ) = 2.0 * aJ2bt( 2, j ) * aJt( 1, j ) + aJ2bt( 1, j ) * aJt( 0, j );
            }

            // Block (4) ------------------------------------------------
            aJ3bt( 2, 2 ) =
                    2.0 * aJ2bt( 2, 0 ) * aJt( 0, 1 ) + 2.0 * aJ2bt( 2, 1 ) * aJt( 0, 0 ) +    //
                    1.0 * aJ2bt( 0, 1 ) * aJt( 1, 0 ) + 1.0 * aJ2bt( 0, 0 ) * aJt( 1, 1 );
            aJ3bt( 3, 2 ) =
                    2.0 * aJ2bt( 2, 0 ) * aJt( 1, 1 ) + 2.0 * aJ2bt( 2, 1 ) * aJt( 1, 0 ) +    //
                    1.0 * aJ2bt( 1, 1 ) * aJt( 0, 0 ) + 1.0 * aJ2bt( 1, 0 ) * aJt( 0, 1 );

            // third help matrix
            aJ3ct = ad3NdXi3 * aXHat;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_matrices_for_second_derivative_3d(
                const Matrix< DDRMat >& aJt,
                Matrix< DDRMat >&       aKt,
                Matrix< DDRMat >&       aLt,
                const Matrix< DDRMat >& ad2NdXi2,
                const Matrix< DDRMat >& aXHat )
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 6, 6 );
            for ( uint j = 0; j < 3; ++j )
            {
                aLt( 0, j ) = std::pow( aJt( 0, j ), 2.0 );
                aLt( 1, j ) = std::pow( aJt( 1, j ), 2.0 );
                aLt( 2, j ) = std::pow( aJt( 2, j ), 2.0 );
                aLt( 3, j ) = aJt( 1, j ) * aJt( 2, j );
                aLt( 4, j ) = aJt( 0, j ) * aJt( 2, j );
                aLt( 5, j ) = aJt( 0, j ) * aJt( 1, j );
            }

            aLt( 0, 3 ) = 2.0 * aJt( 0, 1 ) * aJt( 0, 2 );
            aLt( 1, 3 ) = 2.0 * aJt( 1, 1 ) * aJt( 1, 2 );
            aLt( 2, 3 ) = 2.0 * aJt( 2, 1 ) * aJt( 2, 2 );
            aLt( 3, 3 ) = aJt( 1, 1 ) * aJt( 2, 2 ) + aJt( 2, 1 ) * aJt( 1, 2 );
            aLt( 4, 3 ) = aJt( 0, 1 ) * aJt( 2, 2 ) + aJt( 2, 1 ) * aJt( 0, 2 );
            aLt( 5, 3 ) = aJt( 0, 1 ) * aJt( 1, 2 ) + aJt( 1, 1 ) * aJt( 0, 2 );

            aLt( 0, 4 ) = 2.0 * aJt( 0, 0 ) * aJt( 0, 2 );
            aLt( 1, 4 ) = 2.0 * aJt( 1, 0 ) * aJt( 1, 2 );
            aLt( 2, 4 ) = 2.0 * aJt( 2, 0 ) * aJt( 2, 2 );
            aLt( 3, 4 ) = aJt( 1, 0 ) * aJt( 2, 2 ) + aJt( 2, 0 ) * aJt( 1, 2 );
            aLt( 4, 4 ) = aJt( 0, 0 ) * aJt( 2, 2 ) + aJt( 2, 0 ) * aJt( 0, 2 );
            aLt( 5, 4 ) = aJt( 0, 0 ) * aJt( 1, 2 ) + aJt( 1, 0 ) * aJt( 0, 2 );

            aLt( 0, 5 ) = 2.0 * aJt( 0, 0 ) * aJt( 0, 1 );
            aLt( 1, 5 ) = 2.0 * aJt( 1, 0 ) * aJt( 1, 1 );
            aLt( 2, 5 ) = 2.0 * aJt( 2, 0 ) * aJt( 2, 1 );
            aLt( 3, 5 ) = aJt( 1, 0 ) * aJt( 2, 1 ) + aJt( 2, 0 ) * aJt( 1, 1 );
            aLt( 4, 5 ) = aJt( 0, 0 ) * aJt( 2, 1 ) + aJt( 2, 0 ) * aJt( 0, 1 );
            aLt( 5, 5 ) = aJt( 0, 0 ) * aJt( 1, 1 ) + aJt( 1, 0 ) * aJt( 0, 1 );
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::eval_matrices_for_third_derivative_3d(
                const Matrix< DDRMat >& aJt,
                const Matrix< DDRMat >& aJ2bt,
                Matrix< DDRMat >&       aJ3at,
                Matrix< DDRMat >&       aJ3bt,
                Matrix< DDRMat >&       aJ3ct,
                const Matrix< DDRMat >& ad3NdXi3,
                const Matrix< DDRMat >& aXHat )
        {
            // first help matrix
            aJ3at.set_size( 10, 10 );

            /* matrix structured into 9 parts
             *  ___________________________     ________
             *  | * * * | * * * * * * | * |     | ,xxx |
             *  | *(1)* | * *(2)* * * |(3)|     | ,yyy |
             *  |_*_*_*_|_*_*_*_*_*_*_|_*_|     | ,zzz |
             *  | * * * | * * * * * * | * |     | ,xxy |
             *  | * * * | * * * * * * | * |     | ,xxz |
             *  | *(4)* | * *(5)* * * |(6)|  *  | ,xyy |
             *  | * * * | * * * * * * | * |     | ,yyz |
             *  | * * * | * * * * * * | * |     | ,xzz |
             *  |_*_*_*_|_*_*_*_*_*_*_|_*_|     | ,yzz |
             *  |_*(7)*_|_*_*(8)*_*_*_|(9)|     |_,xyz_|
             */

            // Block (1) ------------------------------------------------
            for ( uint j = 0; j < 3; ++j )
            {
                aJ3at( 0, j ) = std::pow( aJt( 0, j ), 3 );
                aJ3at( 1, j ) = std::pow( aJt( 1, j ), 3 );
                aJ3at( 2, j ) = std::pow( aJt( 2, j ), 3 );
            }

            // Block (2) ------------------------------------------------
            aJ3at( 0, 3 ) = 3.0 * std::pow( aJt( 0, 0 ), 2.0 ) * aJt( 0, 1 );
            aJ3at( 1, 3 ) = 3.0 * std::pow( aJt( 1, 0 ), 2.0 ) * aJt( 1, 1 );
            aJ3at( 2, 3 ) = 3.0 * std::pow( aJt( 2, 0 ), 2.0 ) * aJt( 2, 1 );

            aJ3at( 0, 4 ) = 3.0 * std::pow( aJt( 0, 0 ), 2.0 ) * aJt( 0, 2 );
            aJ3at( 1, 4 ) = 3.0 * std::pow( aJt( 1, 0 ), 2.0 ) * aJt( 1, 2 );
            aJ3at( 2, 4 ) = 3.0 * std::pow( aJt( 2, 0 ), 2.0 ) * aJt( 2, 2 );

            aJ3at( 0, 5 ) = 3.0 * std::pow( aJt( 0, 1 ), 2.0 ) * aJt( 0, 0 );
            aJ3at( 1, 5 ) = 3.0 * std::pow( aJt( 1, 1 ), 2.0 ) * aJt( 1, 0 );
            aJ3at( 2, 5 ) = 3.0 * std::pow( aJt( 2, 1 ), 2.0 ) * aJt( 2, 0 );

            aJ3at( 0, 6 ) = 3.0 * std::pow( aJt( 0, 1 ), 2.0 ) * aJt( 0, 2 );
            aJ3at( 1, 6 ) = 3.0 * std::pow( aJt( 1, 1 ), 2.0 ) * aJt( 1, 2 );
            aJ3at( 2, 6 ) = 3.0 * std::pow( aJt( 2, 1 ), 2.0 ) * aJt( 2, 2 );

            aJ3at( 0, 7 ) = 3.0 * std::pow( aJt( 0, 2 ), 2.0 ) * aJt( 0, 0 );
            aJ3at( 1, 7 ) = 3.0 * std::pow( aJt( 1, 2 ), 2.0 ) * aJt( 1, 0 );
            aJ3at( 2, 7 ) = 3.0 * std::pow( aJt( 2, 2 ), 2.0 ) * aJt( 2, 0 );

            aJ3at( 0, 8 ) = 3.0 * std::pow( aJt( 0, 2 ), 2.0 ) * aJt( 0, 1 );
            aJ3at( 1, 8 ) = 3.0 * std::pow( aJt( 1, 2 ), 2.0 ) * aJt( 1, 1 );
            aJ3at( 2, 8 ) = 3.0 * std::pow( aJt( 2, 2 ), 2.0 ) * aJt( 2, 1 );

            // Block (3) ------------------------------------------------
            aJ3at( 0, 9 ) = 6.0 * aJt( 0, 0 ) * aJt( 0, 1 ) * aJt( 0, 2 );
            aJ3at( 1, 9 ) = 6.0 * aJt( 1, 0 ) * aJt( 1, 1 ) * aJt( 1, 2 );
            aJ3at( 2, 9 ) = 6.0 * aJt( 2, 0 ) * aJt( 2, 1 ) * aJt( 2, 2 );

            // Block (4) ------------------------------------------------
            for ( uint j = 0; j < 3; ++j )
            {
                aJ3at( 3, j ) = std::pow( aJt( 0, j ), 2.0 ) * aJt( 1, j );
                aJ3at( 4, j ) = std::pow( aJt( 0, j ), 2.0 ) * aJt( 2, j );
                aJ3at( 5, j ) = std::pow( aJt( 1, j ), 2.0 ) * aJt( 0, j );
                aJ3at( 6, j ) = std::pow( aJt( 1, j ), 2.0 ) * aJt( 2, j );
                aJ3at( 7, j ) = std::pow( aJt( 2, j ), 2.0 ) * aJt( 0, j );
                aJ3at( 8, j ) = std::pow( aJt( 2, j ), 2.0 ) * aJt( 1, j );
            }

            // Block (5) ------------------------------------------------
            aJ3at( 3, 3 ) = std::pow( aJt( 0, 0 ), 2.0 ) * aJt( 1, 1 ) + 2.0 * aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 0, 1 );
            aJ3at( 4, 3 ) = std::pow( aJt( 0, 0 ), 2.0 ) * aJt( 2, 1 ) + 2.0 * aJt( 0, 0 ) * aJt( 2, 0 ) * aJt( 0, 1 );
            aJ3at( 5, 3 ) = std::pow( aJt( 1, 0 ), 2.0 ) * aJt( 0, 1 ) + 2.0 * aJt( 1, 0 ) * aJt( 0, 0 ) * aJt( 1, 1 );
            aJ3at( 6, 3 ) = std::pow( aJt( 1, 0 ), 2.0 ) * aJt( 2, 1 ) + 2.0 * aJt( 1, 0 ) * aJt( 2, 0 ) * aJt( 1, 1 );
            aJ3at( 7, 3 ) = std::pow( aJt( 2, 0 ), 2.0 ) * aJt( 0, 1 ) + 2.0 * aJt( 2, 0 ) * aJt( 0, 0 ) * aJt( 2, 1 );
            aJ3at( 8, 3 ) = std::pow( aJt( 2, 0 ), 2.0 ) * aJt( 1, 1 ) + 2.0 * aJt( 2, 0 ) * aJt( 1, 0 ) * aJt( 2, 1 );

            aJ3at( 3, 4 ) = std::pow( aJt( 0, 0 ), 2.0 ) * aJt( 1, 2 ) + 2.0 * aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 0, 2 );
            aJ3at( 4, 4 ) = std::pow( aJt( 0, 0 ), 2.0 ) * aJt( 2, 2 ) + 2.0 * aJt( 0, 0 ) * aJt( 2, 0 ) * aJt( 0, 2 );
            aJ3at( 5, 4 ) = std::pow( aJt( 1, 0 ), 2.0 ) * aJt( 0, 2 ) + 2.0 * aJt( 1, 0 ) * aJt( 0, 0 ) * aJt( 1, 2 );
            aJ3at( 6, 4 ) = std::pow( aJt( 1, 0 ), 2.0 ) * aJt( 2, 2 ) + 2.0 * aJt( 1, 0 ) * aJt( 2, 0 ) * aJt( 1, 2 );
            aJ3at( 7, 4 ) = std::pow( aJt( 2, 0 ), 2.0 ) * aJt( 0, 2 ) + 2.0 * aJt( 2, 0 ) * aJt( 0, 0 ) * aJt( 2, 2 );
            aJ3at( 8, 4 ) = std::pow( aJt( 2, 0 ), 2.0 ) * aJt( 1, 2 ) + 2.0 * aJt( 2, 0 ) * aJt( 1, 0 ) * aJt( 2, 2 );

            aJ3at( 3, 5 ) = std::pow( aJt( 0, 1 ), 2.0 ) * aJt( 1, 0 ) + 2.0 * aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 0, 0 );
            aJ3at( 4, 5 ) = std::pow( aJt( 0, 1 ), 2.0 ) * aJt( 2, 0 ) + 2.0 * aJt( 0, 1 ) * aJt( 2, 1 ) * aJt( 0, 0 );
            aJ3at( 5, 5 ) = std::pow( aJt( 1, 1 ), 2.0 ) * aJt( 0, 0 ) + 2.0 * aJt( 1, 1 ) * aJt( 0, 1 ) * aJt( 1, 0 );
            aJ3at( 6, 5 ) = std::pow( aJt( 1, 1 ), 2.0 ) * aJt( 2, 0 ) + 2.0 * aJt( 1, 1 ) * aJt( 2, 1 ) * aJt( 1, 0 );
            aJ3at( 7, 5 ) = std::pow( aJt( 2, 1 ), 2.0 ) * aJt( 0, 0 ) + 2.0 * aJt( 2, 1 ) * aJt( 0, 1 ) * aJt( 2, 0 );
            aJ3at( 8, 5 ) = std::pow( aJt( 2, 1 ), 2.0 ) * aJt( 1, 0 ) + 2.0 * aJt( 2, 1 ) * aJt( 1, 1 ) * aJt( 2, 0 );

            aJ3at( 3, 6 ) = std::pow( aJt( 0, 1 ), 2.0 ) * aJt( 1, 2 ) + 2.0 * aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 0, 2 );
            aJ3at( 4, 6 ) = std::pow( aJt( 0, 1 ), 2.0 ) * aJt( 2, 2 ) + 2.0 * aJt( 0, 1 ) * aJt( 2, 1 ) * aJt( 0, 2 );
            aJ3at( 5, 6 ) = std::pow( aJt( 1, 1 ), 2.0 ) * aJt( 0, 2 ) + 2.0 * aJt( 1, 1 ) * aJt( 0, 1 ) * aJt( 1, 2 );
            aJ3at( 6, 6 ) = std::pow( aJt( 1, 1 ), 2.0 ) * aJt( 2, 2 ) + 2.0 * aJt( 1, 1 ) * aJt( 2, 1 ) * aJt( 1, 2 );
            aJ3at( 7, 6 ) = std::pow( aJt( 2, 1 ), 2.0 ) * aJt( 0, 2 ) + 2.0 * aJt( 2, 1 ) * aJt( 0, 1 ) * aJt( 2, 2 );
            aJ3at( 8, 6 ) = std::pow( aJt( 2, 1 ), 2.0 ) * aJt( 1, 2 ) + 2.0 * aJt( 2, 1 ) * aJt( 1, 1 ) * aJt( 2, 2 );

            aJ3at( 3, 7 ) = std::pow( aJt( 0, 2 ), 2.0 ) * aJt( 1, 0 ) + 2.0 * aJt( 0, 2 ) * aJt( 1, 2 ) * aJt( 0, 0 );
            aJ3at( 4, 7 ) = std::pow( aJt( 0, 2 ), 2.0 ) * aJt( 2, 0 ) + 2.0 * aJt( 0, 2 ) * aJt( 2, 2 ) * aJt( 0, 0 );
            aJ3at( 5, 7 ) = std::pow( aJt( 1, 2 ), 2.0 ) * aJt( 0, 0 ) + 2.0 * aJt( 1, 2 ) * aJt( 0, 2 ) * aJt( 1, 0 );
            aJ3at( 6, 7 ) = std::pow( aJt( 1, 2 ), 2.0 ) * aJt( 2, 0 ) + 2.0 * aJt( 1, 2 ) * aJt( 2, 2 ) * aJt( 1, 0 );
            aJ3at( 7, 7 ) = std::pow( aJt( 2, 2 ), 2.0 ) * aJt( 0, 0 ) + 2.0 * aJt( 2, 2 ) * aJt( 0, 2 ) * aJt( 2, 0 );
            aJ3at( 8, 7 ) = std::pow( aJt( 2, 2 ), 2.0 ) * aJt( 1, 0 ) + 2.0 * aJt( 2, 2 ) * aJt( 1, 2 ) * aJt( 2, 0 );

            aJ3at( 3, 8 ) = std::pow( aJt( 0, 2 ), 2.0 ) * aJt( 1, 1 ) + 2.0 * aJt( 0, 2 ) * aJt( 1, 2 ) * aJt( 0, 1 );
            aJ3at( 4, 8 ) = std::pow( aJt( 0, 2 ), 2.0 ) * aJt( 2, 1 ) + 2.0 * aJt( 0, 2 ) * aJt( 2, 2 ) * aJt( 0, 1 );
            aJ3at( 5, 8 ) = std::pow( aJt( 1, 2 ), 2.0 ) * aJt( 0, 1 ) + 2.0 * aJt( 1, 2 ) * aJt( 0, 2 ) * aJt( 1, 1 );
            aJ3at( 6, 8 ) = std::pow( aJt( 1, 2 ), 2.0 ) * aJt( 2, 1 ) + 2.0 * aJt( 1, 2 ) * aJt( 2, 2 ) * aJt( 1, 1 );
            aJ3at( 7, 8 ) = std::pow( aJt( 2, 2 ), 2.0 ) * aJt( 0, 1 ) + 2.0 * aJt( 2, 2 ) * aJt( 0, 2 ) * aJt( 2, 1 );
            aJ3at( 8, 8 ) = std::pow( aJt( 2, 2 ), 2.0 ) * aJt( 1, 1 ) + 2.0 * aJt( 2, 2 ) * aJt( 1, 2 ) * aJt( 2, 1 );

            // Block (6) ------------------------------------------------
            aJ3at( 3, 9 ) =
                    2.0 * aJt( 0, 0 ) * aJt( 0, 1 ) * aJt( 1, 2 ) +    //
                    2.0 * aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 0, 2 ) +    //
                    2.0 * aJt( 1, 0 ) * aJt( 0, 1 ) * aJt( 0, 2 );
            aJ3at( 4, 9 ) =
                    2.0 * aJt( 0, 0 ) * aJt( 0, 1 ) * aJt( 2, 2 ) +    //
                    2.0 * aJt( 0, 0 ) * aJt( 2, 1 ) * aJt( 0, 2 ) +    //
                    2.0 * aJt( 2, 0 ) * aJt( 0, 1 ) * aJt( 0, 2 );
            aJ3at( 5, 9 ) =
                    2.0 * aJt( 1, 0 ) * aJt( 1, 1 ) * aJt( 0, 2 ) +    //
                    2.0 * aJt( 1, 0 ) * aJt( 0, 1 ) * aJt( 1, 2 ) +    //
                    2.0 * aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 1, 2 );
            aJ3at( 6, 9 ) =
                    2.0 * aJt( 1, 0 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +    //
                    2.0 * aJt( 1, 0 ) * aJt( 2, 1 ) * aJt( 1, 2 ) +    //
                    2.0 * aJt( 2, 0 ) * aJt( 1, 1 ) * aJt( 1, 2 );
            aJ3at( 7, 9 ) =
                    2.0 * aJt( 2, 0 ) * aJt( 2, 1 ) * aJt( 0, 2 ) +    //
                    2.0 * aJt( 2, 0 ) * aJt( 0, 1 ) * aJt( 2, 2 ) +    //
                    2.0 * aJt( 0, 0 ) * aJt( 2, 1 ) * aJt( 2, 2 );
            aJ3at( 8, 9 ) =
                    2.0 * aJt( 2, 0 ) * aJt( 2, 1 ) * aJt( 1, 2 ) +    //
                    2.0 * aJt( 2, 0 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +    //
                    2.0 * aJt( 1, 0 ) * aJt( 2, 1 ) * aJt( 2, 2 );

            // Block (7) ------------------------------------------------
            for ( uint j = 0; j < 3; ++j )
            {
                aJ3at( 9, j ) = aJt( 0, j ) * aJt( 1, j ) * aJt( 2, j );
            }

            // Block (8) ------------------------------------------------
            aJ3at( 9, 3 ) =
                    aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 2, 1 ) +    //
                    aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 2, 0 ) +    //
                    aJt( 0, 1 ) * aJt( 1, 0 ) * aJt( 2, 0 );

            aJ3at( 9, 4 ) =
                    aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 2, 2 ) +    //
                    aJt( 0, 0 ) * aJt( 1, 2 ) * aJt( 2, 0 ) +    //
                    aJt( 0, 2 ) * aJt( 1, 0 ) * aJt( 2, 0 );

            aJ3at( 9, 5 ) =
                    aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 2, 0 ) +    //
                    aJt( 0, 1 ) * aJt( 1, 0 ) * aJt( 2, 1 ) +    //
                    aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 2, 1 );

            aJ3at( 9, 6 ) =
                    aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +    //
                    aJt( 0, 1 ) * aJt( 1, 2 ) * aJt( 2, 1 ) +    //
                    aJt( 0, 2 ) * aJt( 1, 1 ) * aJt( 2, 1 );

            aJ3at( 9, 7 ) =
                    aJt( 0, 2 ) * aJt( 1, 2 ) * aJt( 2, 0 ) +    //
                    aJt( 0, 2 ) * aJt( 1, 0 ) * aJt( 2, 2 ) +    //
                    aJt( 0, 0 ) * aJt( 1, 2 ) * aJt( 2, 2 );

            aJ3at( 9, 8 ) =
                    aJt( 0, 2 ) * aJt( 1, 2 ) * aJt( 2, 1 ) +    //
                    aJt( 0, 2 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +    //
                    aJt( 0, 1 ) * aJt( 1, 2 ) * aJt( 2, 2 );

            // Block (9) ------------------------------------------------
            aJ3at( 9, 9 ) =
                    aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +    //
                    aJt( 0, 2 ) * aJt( 1, 1 ) * aJt( 2, 0 ) +    //
                    aJt( 0, 1 ) * aJt( 1, 2 ) * aJt( 2, 0 ) +    //
                    aJt( 0, 0 ) * aJt( 1, 2 ) * aJt( 2, 1 ) +    //
                    aJt( 0, 2 ) * aJt( 1, 0 ) * aJt( 2, 1 ) +    //
                    aJt( 0, 1 ) * aJt( 1, 0 ) * aJt( 2, 2 );

            // second help matrix
            aJ3bt.set_size( 10, 6 );

            /* matrix structured into 6 parts
             *  _________________
             *  | * * * | * * * |
             *  | *(1)* | *(2)* |    _______
             *  |_*_*_*_|_*_*_*_|    | ,xx |
             *  | * * * | * * * |    | ,yy |
             *  | * * * | * * * |    | ,zz |
             *  | *(3)* | *(4)* |  * | ,yz |
             *  | * * * | * * * |    | ,xz |
             *  | * * * | * * * |    |_,xy_|
             *  |_*_*_*_|_*_*_*_|
             *  |_*(5)*_|_*(6)*_|
             */

            // Block (1) ------------------------------------------------
            for ( uint j = 0; j < 3; ++j )
            {
                aJ3bt( 0, j ) = 3.0 * aJ2bt( 0, j ) * aJt( 0, j );
                aJ3bt( 1, j ) = 3.0 * aJ2bt( 1, j ) * aJt( 1, j );
                aJ3bt( 2, j ) = 3.0 * aJ2bt( 2, j ) * aJt( 2, j );
            }

            // Block (2) ------------------------------------------------
            aJ3bt( 0, 5 ) =
                    3.0 * aJ2bt( 0, 0 ) * aJt( 0, 1 ) +    //
                    3.0 * aJ2bt( 0, 1 ) * aJt( 0, 0 );
            aJ3bt( 1, 5 ) =
                    3.0 * aJ2bt( 1, 0 ) * aJt( 1, 1 ) +    //
                    3.0 * aJ2bt( 1, 1 ) * aJt( 1, 0 );
            aJ3bt( 2, 5 ) =
                    3.0 * aJ2bt( 2, 0 ) * aJt( 2, 1 ) +    //
                    3.0 * aJ2bt( 2, 1 ) * aJt( 2, 0 );

            aJ3bt( 0, 3 ) =
                    3.0 * aJ2bt( 0, 1 ) * aJt( 0, 2 ) +    //
                    3.0 * aJ2bt( 0, 2 ) * aJt( 0, 1 );
            aJ3bt( 1, 3 ) =
                    3.0 * aJ2bt( 1, 1 ) * aJt( 1, 2 ) +    //
                    3.0 * aJ2bt( 1, 2 ) * aJt( 1, 1 );
            aJ3bt( 2, 3 ) =
                    3.0 * aJ2bt( 2, 1 ) * aJt( 2, 2 ) +    //
                    3.0 * aJ2bt( 2, 2 ) * aJt( 2, 1 );

            aJ3bt( 0, 4 ) =
                    3.0 * aJ2bt( 0, 0 ) * aJt( 0, 2 ) +    //
                    3.0 * aJ2bt( 0, 2 ) * aJt( 0, 0 );
            aJ3bt( 1, 4 ) =
                    3.0 * aJ2bt( 1, 0 ) * aJt( 1, 2 ) +    //
                    3.0 * aJ2bt( 1, 2 ) * aJt( 1, 0 );
            aJ3bt( 2, 4 ) =
                    3.0 * aJ2bt( 2, 0 ) * aJt( 2, 2 ) +    //
                    3.0 * aJ2bt( 2, 2 ) * aJt( 2, 0 );

            // Block (3) ------------------------------------------------
            for ( uint j = 0; j < 3; ++j )
            {
                aJ3bt( 3, j ) = 2.0 * aJ2bt( 5, j ) * aJt( 0, j ) + aJ2bt( 0, j ) * aJt( 1, j );
                aJ3bt( 4, j ) = 2.0 * aJ2bt( 4, j ) * aJt( 0, j ) + aJ2bt( 0, j ) * aJt( 2, j );
                aJ3bt( 5, j ) = 2.0 * aJ2bt( 5, j ) * aJt( 1, j ) + aJ2bt( 1, j ) * aJt( 0, j );
                aJ3bt( 6, j ) = 2.0 * aJ2bt( 3, j ) * aJt( 1, j ) + aJ2bt( 1, j ) * aJt( 2, j );
                aJ3bt( 7, j ) = 2.0 * aJ2bt( 4, j ) * aJt( 2, j ) + aJ2bt( 2, j ) * aJt( 0, j );
                aJ3bt( 8, j ) = 2.0 * aJ2bt( 3, j ) * aJt( 2, j ) + aJ2bt( 2, j ) * aJt( 1, j );
            }

            // Block (4) ------------------------------------------------
            aJ3bt( 3, 5 ) = 2.0 * aJ2bt( 5, 0 ) * aJt( 0, 1 ) + 2.0 * aJ2bt( 5, 1 ) * aJt( 0, 0 ) +    //
                            aJ2bt( 0, 1 ) * aJt( 1, 0 ) + aJ2bt( 0, 0 ) * aJt( 1, 1 );
            aJ3bt( 4, 5 ) = 2.0 * aJ2bt( 4, 0 ) * aJt( 0, 1 ) + 2.0 * aJ2bt( 4, 1 ) * aJt( 0, 0 ) +    //
                            aJ2bt( 0, 1 ) * aJt( 2, 0 ) + aJ2bt( 0, 0 ) * aJt( 2, 1 );
            aJ3bt( 5, 5 ) = 2.0 * aJ2bt( 5, 0 ) * aJt( 1, 1 ) + 2.0 * aJ2bt( 5, 1 ) * aJt( 1, 0 ) +    //
                            aJ2bt( 1, 1 ) * aJt( 0, 0 ) + aJ2bt( 1, 0 ) * aJt( 0, 1 );
            aJ3bt( 6, 5 ) = 2.0 * aJ2bt( 3, 0 ) * aJt( 1, 1 ) + 2.0 * aJ2bt( 3, 1 ) * aJt( 1, 0 ) +    //
                            aJ2bt( 1, 1 ) * aJt( 2, 0 ) + aJ2bt( 1, 0 ) * aJt( 2, 1 );
            aJ3bt( 7, 5 ) = 2.0 * aJ2bt( 4, 0 ) * aJt( 2, 1 ) + 2.0 * aJ2bt( 4, 1 ) * aJt( 2, 0 ) +    //
                            aJ2bt( 2, 1 ) * aJt( 0, 0 ) + aJ2bt( 2, 0 ) * aJt( 0, 1 );
            aJ3bt( 8, 5 ) = 2.0 * aJ2bt( 3, 0 ) * aJt( 2, 1 ) + 2.0 * aJ2bt( 3, 1 ) * aJt( 2, 0 ) +    //
                            aJ2bt( 2, 1 ) * aJt( 1, 0 ) + aJ2bt( 2, 0 ) * aJt( 1, 1 );

            aJ3bt( 3, 3 ) = 2.0 * aJ2bt( 5, 1 ) * aJt( 0, 2 ) + 2.0 * aJ2bt( 5, 2 ) * aJt( 0, 1 ) +    //
                            aJ2bt( 0, 2 ) * aJt( 1, 1 ) + aJ2bt( 0, 1 ) * aJt( 1, 2 );
            aJ3bt( 4, 3 ) = 2.0 * aJ2bt( 4, 1 ) * aJt( 0, 2 ) + 2.0 * aJ2bt( 4, 2 ) * aJt( 0, 1 ) +    //
                            aJ2bt( 0, 2 ) * aJt( 2, 1 ) + aJ2bt( 0, 1 ) * aJt( 2, 2 );
            aJ3bt( 5, 3 ) = 2.0 * aJ2bt( 5, 1 ) * aJt( 1, 2 ) + 2.0 * aJ2bt( 5, 2 ) * aJt( 1, 1 ) +    //
                            aJ2bt( 1, 2 ) * aJt( 0, 1 ) + aJ2bt( 1, 1 ) * aJt( 0, 2 );
            aJ3bt( 6, 3 ) = 2.0 * aJ2bt( 3, 1 ) * aJt( 1, 2 ) + 2.0 * aJ2bt( 3, 2 ) * aJt( 1, 1 ) +    //
                            aJ2bt( 1, 2 ) * aJt( 2, 1 ) + aJ2bt( 1, 1 ) * aJt( 2, 2 );
            aJ3bt( 7, 3 ) = 2.0 * aJ2bt( 4, 1 ) * aJt( 2, 2 ) + 2.0 * aJ2bt( 4, 2 ) * aJt( 2, 1 ) +    //
                            aJ2bt( 2, 2 ) * aJt( 0, 1 ) + aJ2bt( 2, 1 ) * aJt( 0, 2 );
            aJ3bt( 8, 3 ) = 2.0 * aJ2bt( 3, 1 ) * aJt( 2, 2 ) + 2.0 * aJ2bt( 3, 2 ) * aJt( 2, 1 ) +    //
                            aJ2bt( 2, 2 ) * aJt( 1, 1 ) + aJ2bt( 2, 1 ) * aJt( 1, 2 );

            aJ3bt( 3, 4 ) = 2.0 * aJ2bt( 5, 0 ) * aJt( 0, 2 ) + 2.0 * aJ2bt( 5, 2 ) * aJt( 0, 0 ) +    //
                            aJ2bt( 0, 2 ) * aJt( 1, 0 ) + aJ2bt( 0, 0 ) * aJt( 1, 2 );
            aJ3bt( 4, 4 ) = 2.0 * aJ2bt( 4, 0 ) * aJt( 0, 2 ) + 2.0 * aJ2bt( 4, 2 ) * aJt( 0, 0 ) +    //
                            aJ2bt( 0, 2 ) * aJt( 2, 0 ) + aJ2bt( 0, 0 ) * aJt( 2, 2 );
            aJ3bt( 5, 4 ) = 2.0 * aJ2bt( 5, 0 ) * aJt( 1, 2 ) + 2.0 * aJ2bt( 5, 2 ) * aJt( 1, 0 ) +    //
                            aJ2bt( 1, 2 ) * aJt( 0, 0 ) + aJ2bt( 1, 0 ) * aJt( 0, 2 );
            aJ3bt( 6, 4 ) = 2.0 * aJ2bt( 3, 0 ) * aJt( 1, 2 ) + 2.0 * aJ2bt( 3, 2 ) * aJt( 1, 0 ) +    //
                            aJ2bt( 1, 2 ) * aJt( 2, 0 ) + aJ2bt( 1, 0 ) * aJt( 2, 2 );
            aJ3bt( 7, 4 ) = 2.0 * aJ2bt( 4, 0 ) * aJt( 2, 2 ) + 2.0 * aJ2bt( 4, 2 ) * aJt( 2, 0 ) +    //
                            aJ2bt( 2, 2 ) * aJt( 0, 0 ) + aJ2bt( 2, 0 ) * aJt( 0, 2 );
            aJ3bt( 8, 4 ) = 2.0 * aJ2bt( 3, 0 ) * aJt( 2, 2 ) + 2.0 * aJ2bt( 3, 2 ) * aJt( 2, 0 ) +    //
                            aJ2bt( 2, 2 ) * aJt( 1, 0 ) + aJ2bt( 2, 0 ) * aJt( 1, 2 );

            // Block (5) ------------------------------------------------
            for ( uint j = 0; j < 3; ++j )
            {
                aJ3bt( 9, j ) =
                        aJ2bt( 4, j ) * aJt( 1, j ) + aJ2bt( 3, j ) * aJt( 0, j ) + aJ2bt( 5, j ) * aJt( 2, j );
            }

            // Block (6) ------------------------------------------------
            aJ3bt( 9, 5 ) =
                    aJ2bt( 5, 0 ) * aJt( 2, 1 ) + aJ2bt( 5, 1 ) * aJt( 2, 0 ) + aJ2bt( 3, 0 ) * aJt( 0, 1 ) +    //
                    aJ2bt( 3, 1 ) * aJt( 0, 0 ) + aJ2bt( 4, 0 ) * aJt( 1, 1 ) + aJ2bt( 4, 1 ) * aJt( 1, 0 );

            aJ3bt( 9, 3 ) =
                    aJ2bt( 5, 1 ) * aJt( 2, 2 ) + aJ2bt( 5, 2 ) * aJt( 2, 1 ) + aJ2bt( 3, 1 ) * aJt( 0, 2 ) +    //
                    aJ2bt( 3, 2 ) * aJt( 0, 1 ) + aJ2bt( 4, 1 ) * aJt( 1, 2 ) + aJ2bt( 4, 2 ) * aJt( 1, 1 );

            aJ3bt( 9, 4 ) =
                    aJ2bt( 5, 0 ) * aJt( 2, 2 ) + aJ2bt( 5, 2 ) * aJt( 2, 0 ) + aJ2bt( 3, 0 ) * aJt( 0, 2 ) +    //
                    aJ2bt( 3, 2 ) * aJt( 0, 0 ) + aJ2bt( 4, 0 ) * aJt( 1, 2 ) + aJ2bt( 4, 2 ) * aJt( 1, 0 );

            // third help matrix
            aJ3ct = ad3NdXi3 * aXHat;
        }

        //------------------------------------------------------------------------------

        void
        Space_Interpolator::set_function_pointers()
        {
            switch ( mInterpolationShape )
            {
                case CellShape::GENERAL:
                case CellShape::PARALLEL:
                case CellShape::STRAIGHT:
                {
                    // get number of dimensions and set pointer to function
                    // for second space derivative
                    switch ( mNumSpaceDim )
                    {
                        case 1:
                        {
                            mInvSpaceJacFunc               = &Space_Interpolator::eval_inverse_space_jacobian_1d;
                            mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_1d;
                            mThirdDerivativeMatricesSpace  = this->eval_matrices_for_third_derivative_1d;
                            break;
                        }
                        case 2:
                        {
                            switch ( mGeometryType )
                            {
                                case Geometry_Type::LINE:
                                case Geometry_Type::QUAD:
                                {
                                    mInvSpaceJacFunc = &Space_Interpolator::eval_inverse_space_jacobian_2d;
                                    break;
                                }
                                // effectively inverse for non-square triangle element jacobian
                                case Geometry_Type::TRI:
                                {
                                    mInvSpaceJacFunc = &Space_Interpolator::eval_inverse_space_jacobian_2d_tri;
                                    break;
                                }
                                default:
                                {
                                    MORIS_ERROR( false,
                                            " Space_Interpolator::set_function_pointers - only line, quad, tri allowed for 2D." );
                                }
                            }
                            mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_2d;
                            mThirdDerivativeMatricesSpace  = this->eval_matrices_for_third_derivative_2d;
                            break;
                        }

                        case 3:
                        {
                            mInvSpaceJacFunc               = &Space_Interpolator::eval_inverse_space_jacobian_3d;
                            mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_3d;
                            mThirdDerivativeMatricesSpace  = this->eval_matrices_for_third_derivative_3d;
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, " Space_Interpolator::set_function_pointers - unknown number of dimensions." );
                        }
                    }

                    // switch on geometry type
                    if ( mSpaceSideset )
                    {
                        switch ( mGeometryType )
                        {
                            case Geometry_Type::LINE:
                            {
                                // switching based on background element
                                switch ( mIPMappingGeometryType )
                                {
                                    // FIXME: Geometry_Type of LINE only used to prevent
                                    // existing tests that would fail with this change. This is only
                                    // a logical inconsistency though. Can't have a line sideset of a
                                    // line IP cell.
                                    case Geometry_Type::LINE:
                                    case Geometry_Type::QUAD:
                                    {
                                        mMappedPoint.set_size( 3, 1, 0.0 );
                                        break;
                                    }
                                    // if the IP cell is a tri, then it could use 2 or 3 parameterization dimensions
                                    case Geometry_Type::TRI:
                                    {
                                        switch ( mIPMappingNumSpaceParamDim )
                                        {
                                            case 2:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 3, 1, 0.0 );
                                                break;
                                            }
                                            case 3:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 4, 1, 0.0 );
                                                break;
                                            }
                                            default:
                                            {
                                                MORIS_ERROR( false,
                                                        " Space_Interpolator::set_function_pointers - Parametric "
                                                        "space dimensions can only be 2 or 3." );
                                            }
                                        }
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - Line sidesets can only be "
                                                "based on quad or tri interpolation cells." );
                                    }
                                }
                                mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_side_line;
                                mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_side_line;
                                mNormalFunc         = &Space_Interpolator::eval_normal_side_line;
                                // set size for storage
                                mMapFlag = true;
                                break;
                            }
                            case Geometry_Type::TRI:
                            {
                                // the following switch structure determines the mapping size
                                // based on the IP cell used.
                                switch ( mIPMappingGeometryType )
                                {
                                    // Tri geometry type only used to prevent failing tests
                                    case Geometry_Type::TRI:
                                    case Geometry_Type::HEX:
                                    {
                                        mMappedPoint.set_size( 4, 1, 0.0 );
                                        break;
                                    }
                                    case Geometry_Type::TET:
                                    {
                                        switch ( mIPMappingNumSpaceParamDim )
                                        {
                                            case 3:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 4, 1, 0.0 );
                                                break;
                                            }
                                            case 4:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 5, 1, 0.0 );
                                                break;
                                            }
                                            default:
                                            {
                                                MORIS_ERROR( false,
                                                        " Space_Interpolator::set_function_pointers - Parametric "
                                                        "space dimensions can only be 2 or 3." );
                                            }
                                        }
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - tri sidesets can only "
                                                "come from hex or tet IP cells." );
                                    }
                                }

                                // function pointers will not change based on the IP switch structure
                                mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_side_tri;
                                mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_side_tri;
                                mNormalFunc         = &Space_Interpolator::eval_normal_side_tri;

                                // set size for storage
                                mMapFlag = true;
                                break;
                            }
                            case Geometry_Type::QUAD:
                            {
                                mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_side_quad;
                                mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_side_quad;
                                mNormalFunc         = &Space_Interpolator::eval_normal_side_quad;

                                // set size for storage
                                mMapFlag = true;
                                mMappedPoint.set_size( 4, 1, 0.0 );
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false, "Space_Interpolator::set_function_pointers - unknown or not implemented side space geometry type " );
                            }
                        }
                    }
                    else
                    {
                        switch ( mGeometryType )
                        {
                            case Geometry_Type::LINE:
                            {
                                mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_line;
                                mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_line;

                                // switching based on interpolation cell geometry
                                switch ( mIPMappingGeometryType )
                                {
                                    case Geometry_Type::LINE:
                                    {
                                        // set size for storage
                                        mMappedPoint.set_size( 2, 1, 0.0 );
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - "
                                                "bulk line geometry can only be used within line interpolation cells" );
                                    }
                                }
                                break;
                            }

                            case Geometry_Type::QUAD:
                            {
                                mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_quad;
                                mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_quad;

                                // switching based on interpolation cell geometry
                                switch ( mIPMappingGeometryType )
                                {
                                    case Geometry_Type::QUAD:
                                    {
                                        // set size for storage
                                        mMappedPoint.set_size( 3, 1, 0.0 );
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - "
                                                "bulk quad geometry can only be used within quad interpolation cells" );
                                    }
                                }
                                break;
                            }

                            case Geometry_Type::HEX:
                            {
                                mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_hex;
                                mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_hex;

                                // switching based on interpolation cell geometry
                                switch ( mIPMappingGeometryType )
                                {
                                    case Geometry_Type::HEX:
                                    {
                                        // set size for storage
                                        mMappedPoint.set_size( 4, 1, 0.0 );
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - "
                                                "bulk hex geometry can only be used within hex interpolation cells" );
                                    }
                                }
                                break;
                            }

                            case Geometry_Type::TRI:
                            {
                                switch ( mNumSpaceParamDim )
                                {
                                    case 2:
                                    {
                                        mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_tri_param_2;
                                        mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_tri_param_2;
                                        break;
                                    }
                                    case 3:
                                    {
                                        mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_tri_param_3;
                                        mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_tri_param_3;
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false, " Space_Interpolator::set_function_pointers - Parametric space dimensions can only be 2 or 3." );
                                    }
                                }

                                // switching based on interpolation cell geometry
                                switch ( mIPMappingGeometryType )
                                {
                                    // interpolating on a tri background cell
                                    case Geometry_Type::TRI:
                                    {
                                        switch ( mIPMappingNumSpaceParamDim )
                                        {
                                            case 2:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 3, 1, 0.0 );
                                                break;
                                            }
                                            case 3:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 4, 1, 0.0 );
                                                break;
                                            }
                                            default:
                                            {
                                                MORIS_ERROR( false, " Space_Interpolator::set_function_pointers - Parametric space dimensions can only be 2 or 3." );
                                            }
                                        }
                                        break;
                                    }

                                    // interpolating on a quad background cell
                                    case Geometry_Type::QUAD:
                                    {
                                        // set size for storage
                                        mMapFlag = true;
                                        mMappedPoint.set_size( 3, 1, 0.0 );
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - "
                                                "bulk tri geometry can only be used within tri or quad interpolation cells" );
                                    }
                                }
                                break;
                            }

                            case Geometry_Type::TET:
                            {
                                switch ( mNumSpaceParamDim )
                                {
                                    case 3:
                                        mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_tet_param_3;
                                        mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_tet_param_3;
                                        break;
                                    case 4:
                                        mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_tet_param_4;
                                        mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_tet_param_4;
                                        break;
                                    default:
                                        MORIS_ERROR( false, " Space_Interpolator::set_function_pointers - Parametric space dimensions can only be 3 or 4." );
                                }

                                // switching based on interpolation cell geometry
                                switch ( mIPMappingGeometryType )
                                {
                                    // interpolating on a tri background cell
                                    case Geometry_Type::TET:
                                    {

                                        switch ( mIPMappingNumSpaceParamDim )
                                        {
                                            case 3:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 4, 1, 0.0 );
                                                break;
                                            }
                                            case 4:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 5, 1, 0.0 );
                                                break;
                                            }
                                            default:
                                            {
                                                MORIS_ERROR( false, " Space_Interpolator::set_function_pointers - Parametric space dimensions can only be 2 or 3." );
                                            }
                                        }
                                        break;
                                    }

                                    // interpolating on a quad background cell
                                    case Geometry_Type::HEX:
                                    {
                                        // set size for storage
                                        mMapFlag = true;
                                        mMappedPoint.set_size( 4, 1, 0.0 );
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - "
                                                "bulk tet geometry can only be used within tet or hex interpolation cells" );
                                    }
                                }
                                break;
                            }

                            default:
                            {
                                MORIS_ERROR( false, "Space_Interpolator::set_function_pointers - unknown or not implemented space geometry type " );
                            }
                        }
                    }
                    break;
                }

                // ----------------------------------------------------------------
                case CellShape::RECTANGULAR:
                {
                    // get number of dimensions and set pointer to function
                    // for second space derivative
                    switch ( mNumSpaceDim )
                    {
                        case 1:
                        {
                            mInvSpaceJacFunc               = &Space_Interpolator::eval_inverse_space_jacobian_1d;
                            mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_1d;
                            mThirdDerivativeMatricesSpace  = this->eval_matrices_for_third_derivative_1d;
                            break;
                        }
                        case 2:
                        {
                            switch ( mGeometryType )
                            {
                                case Geometry_Type::LINE:
                                case Geometry_Type::QUAD:
                                {
                                    mInvSpaceJacFunc = &Space_Interpolator::eval_inverse_space_jacobian_2d_rect;
                                    break;
                                }
                                // effectively inverse for non-square triangle element jacobian
                                case Geometry_Type::TRI:
                                {
                                    mInvSpaceJacFunc = &Space_Interpolator::eval_inverse_space_jacobian_2d_tri;
                                    break;
                                }
                                default:
                                {
                                    MORIS_ERROR( false,
                                            " Space_Interpolator::set_function_pointers - only line, quad, tri allowed for 2D." );
                                }
                            }
                            mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_2d;
                            mThirdDerivativeMatricesSpace  = this->eval_matrices_for_third_derivative_2d;
                            break;
                        }

                        case 3:
                        {
                            mInvSpaceJacFunc               = &Space_Interpolator::eval_inverse_space_jacobian_3d_rect;
                            mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_3d;
                            mThirdDerivativeMatricesSpace  = this->eval_matrices_for_third_derivative_3d;
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, " Space_Interpolator::set_function_pointers - unknown number of dimensions." );
                        }
                    }

                    // switch on geometry type
                    if ( mSpaceSideset )
                    {
                        switch ( mGeometryType )
                        {
                            case Geometry_Type::LINE:
                            {
                                // switching based on background element
                                switch ( mIPMappingGeometryType )
                                {
                                    // FIXME: Geometry_Type of LINE only used to prevent
                                    // existing tests that would fail with this change. This is only
                                    // a logical inconsistency though. Can't have a line sideset of a
                                    // line IP cell.
                                    case Geometry_Type::LINE:
                                    case Geometry_Type::QUAD:
                                    {
                                        mMappedPoint.set_size( 3, 1, 0.0 );
                                        break;
                                    }
                                    // if the IP cell is a tri, then it could use 2 or 3 parameterization dimensions
                                    case Geometry_Type::TRI:
                                    {
                                        switch ( mIPMappingNumSpaceParamDim )
                                        {
                                            case 2:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 3, 1, 0.0 );
                                                break;
                                            }
                                            case 3:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 4, 1, 0.0 );
                                                break;
                                            }
                                            default:
                                            {
                                                MORIS_ERROR( false,
                                                        " Space_Interpolator::set_function_pointers - Parametric "
                                                        "space dimensions can only be 2 or 3." );
                                            }
                                        }
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - Line sidesets can only be "
                                                "based on quad or tri interpolation cells." );
                                    }
                                }
                                mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_side_line;
                                mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_side_line;
                                mNormalFunc         = &Space_Interpolator::eval_normal_side_line;
                                // set size for storage
                                mMapFlag = true;
                                break;
                            }
                            case Geometry_Type::TRI:
                            {
                                // the following switch structure determines the mapping size
                                // based on the IP cell used.
                                switch ( mIPMappingGeometryType )
                                {
                                    // Tri geometry type only used to prevent failing tests
                                    case Geometry_Type::TRI:
                                    case Geometry_Type::HEX:
                                    {
                                        mMappedPoint.set_size( 4, 1, 0.0 );
                                        break;
                                    }
                                    case Geometry_Type::TET:
                                    {
                                        switch ( mIPMappingNumSpaceParamDim )
                                        {
                                            case 3:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 4, 1, 0.0 );
                                                break;
                                            }
                                            case 4:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 5, 1, 0.0 );
                                                break;
                                            }
                                            default:
                                            {
                                                MORIS_ERROR( false,
                                                        " Space_Interpolator::set_function_pointers - Parametric "
                                                        "space dimensions can only be 2 or 3." );
                                            }
                                        }
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - tri sidesets can only "
                                                "come from hex or tet IP cells." );
                                    }
                                }

                                // function pointers will not change based on the IP switch structure
                                mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_side_tri;
                                mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_side_tri;
                                mNormalFunc         = &Space_Interpolator::eval_normal_side_tri;

                                // set size for storage
                                mMapFlag = true;
                                break;
                            }
                            case Geometry_Type::QUAD:
                            {
                                mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_side_quad;
                                mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_side_quad;
                                mNormalFunc         = &Space_Interpolator::eval_normal_side_quad;

                                // set size for storage
                                mMapFlag = true;
                                mMappedPoint.set_size( 4, 1, 0.0 );
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false, "Space_Interpolator::set_function_pointers - unknown or not implemented side space geometry type " );
                            }
                        }
                    }
                    else
                    {
                        switch ( mGeometryType )
                        {
                            case Geometry_Type::LINE:
                            {
                                mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_line;
                                mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_line;

                                // switching based on interpolation cell geometry
                                switch ( mIPMappingGeometryType )
                                {
                                    case Geometry_Type::LINE:
                                    {
                                        // set size for storage
                                        mMappedPoint.set_size( 2, 1, 0.0 );
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - "
                                                "bulk line geometry can only be used within line interpolation cells" );
                                    }
                                }
                                break;
                            }

                            case Geometry_Type::QUAD:
                            {
                                mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_quad_rect;
                                mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_quad_rect;

                                // switching based on interpolation cell geometry
                                switch ( mIPMappingGeometryType )
                                {
                                    case Geometry_Type::QUAD:
                                    {
                                        // set size for storage
                                        mMappedPoint.set_size( 3, 1, 0.0 );
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - "
                                                "bulk quad geometry can only be used within quad interpolation cells" );
                                    }
                                }
                                break;
                            }

                            case Geometry_Type::HEX:
                            {
                                mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_hex_rect;
                                mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_hex_rect;

                                // switching based on interpolation cell geometry
                                switch ( mIPMappingGeometryType )
                                {
                                    case Geometry_Type::HEX:
                                    {
                                        // set size for storage
                                        mMappedPoint.set_size( 4, 1, 0.0 );
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - "
                                                "bulk hex geometry can only be used within hex interpolation cells" );
                                    }
                                }
                                break;
                            }

                            case Geometry_Type::TRI:
                            {
                                switch ( mNumSpaceParamDim )
                                {
                                    case 2:
                                    {
                                        mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_tri_param_2;
                                        mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_tri_param_2;
                                        break;
                                    }
                                    case 3:
                                    {
                                        mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_tri_param_3;
                                        mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_tri_param_3;
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false, " Space_Interpolator::set_function_pointers - Parametric space dimensions can only be 2 or 3." );
                                    }
                                }

                                // switching based on interpolation cell geometry
                                switch ( mIPMappingGeometryType )
                                {
                                    // interpolating on a tri background cell
                                    case Geometry_Type::TRI:
                                    {
                                        switch ( mIPMappingNumSpaceParamDim )
                                        {
                                            case 2:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 3, 1, 0.0 );
                                                break;
                                            }
                                            case 3:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 4, 1, 0.0 );
                                                break;
                                            }
                                            default:
                                            {
                                                MORIS_ERROR( false, " Space_Interpolator::set_function_pointers - Parametric space dimensions can only be 2 or 3." );
                                            }
                                        }
                                        break;
                                    }

                                    // interpolating on a quad background cell
                                    case Geometry_Type::QUAD:
                                    {
                                        // set size for storage
                                        mMapFlag = true;
                                        mMappedPoint.set_size( 3, 1, 0.0 );
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - "
                                                "bulk tri geometry can only be used within tri or quad interpolation cells" );
                                    }
                                }
                                break;
                            }

                            case Geometry_Type::TET:
                            {
                                switch ( mNumSpaceParamDim )
                                {
                                    case 3:
                                        mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_tet_param_3;
                                        mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_tet_param_3;
                                        break;
                                    case 4:
                                        mSpaceDetJFunc      = &Space_Interpolator::eval_space_detJ_bulk_tet_param_4;
                                        mSpaceDetJDerivFunc = &Space_Interpolator::eval_space_detJ_deriv_bulk_tet_param_4;
                                        break;
                                    default:
                                        MORIS_ERROR( false, " Space_Interpolator::set_function_pointers - Parametric space dimensions can only be 3 or 4." );
                                }

                                // switching based on interpolation cell geometry
                                switch ( mIPMappingGeometryType )
                                {
                                    // interpolating on a tri background cell
                                    case Geometry_Type::TET:
                                    {

                                        switch ( mIPMappingNumSpaceParamDim )
                                        {
                                            case 3:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 4, 1, 0.0 );
                                                break;
                                            }
                                            case 4:
                                            {
                                                // set size for storage
                                                mMappedPoint.set_size( 5, 1, 0.0 );
                                                break;
                                            }
                                            default:
                                            {
                                                MORIS_ERROR( false, " Space_Interpolator::set_function_pointers - Parametric space dimensions can only be 2 or 3." );
                                            }
                                        }
                                        break;
                                    }

                                    // interpolating on a quad background cell
                                    case Geometry_Type::HEX:
                                    {
                                        // set size for storage
                                        mMapFlag = true;
                                        mMappedPoint.set_size( 4, 1, 0.0 );
                                        break;
                                    }
                                    default:
                                    {
                                        MORIS_ERROR( false,
                                                "Space_Interpolator::set_function_pointers - "
                                                "bulk tet geometry can only be used within tet or hex interpolation cells" );
                                    }
                                }
                                break;
                            }

                            default:
                            {
                                MORIS_ERROR( false, "Space_Interpolator::set_function_pointers - unknown or not implemented space geometry type " );
                            }
                        }
                    }
                    break;
                }

                default:
                {
                    MORIS_ERROR( false,
                            " Space_Interpolator::set_function_pointers - invalid cellshape used." );
                }
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */


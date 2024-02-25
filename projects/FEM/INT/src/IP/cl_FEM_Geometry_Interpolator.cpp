/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Geometry_Interpolator.cpp
 *
 */

#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "fn_norm.hpp"
#include "fn_cross.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"
#include "fn_inv.hpp"
#include "op_div.hpp"
#include "fn_linsolve.hpp"

#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_FEM_Field_Interpolator.hpp"

namespace moris
{
    namespace fem
    {
        // smallest acceptable value for DetJ
        const real Geometry_Interpolator::sDetJLowerLimit = -1.0e-6;

        // smallest acceptable value for DetJ used in building inverse of Jacobian
        const real Geometry_Interpolator::sDetJInvJacLowerLimit = -1.0e-6;

        //------------------------------------------------------------------------------

        Geometry_Interpolator::Geometry_Interpolator(
                const mtk::Interpolation_Rule& aInterpolationRule,
                const mtk::CellShape&          aInterpolationShape,
                const bool                     aSpaceSideset,
                const bool                     aTimeSideset )
        {
            // set bool for side interpolation to true
            mSpaceSideset = aSpaceSideset;

            // set bool for time side interpolation to true
            mTimeSideset = aTimeSideset;

            // create member pointer to space interpolator
            mSpaceInterpolator = new mtk::Space_Interpolator(
                    aInterpolationRule,
                    aInterpolationShape,
                    aSpaceSideset );

            // getting mapping size and flag from space interpolator
            mMappedPoint = mSpaceInterpolator->get_initialized_mapped_point();
            mMapFlag     = mSpaceInterpolator->get_map_flag();

            // create member pointer to time interpolation function
            mTimeInterpolation = aInterpolationRule.create_time_interpolation_function();

            // number of space bases and dimensions
            mNumSpaceParamDim = mSpaceInterpolator->get_number_of_param_dimensions();

            // number of time bases and dimensions
            mNumTimeBases = mTimeInterpolation->get_number_of_bases();
            mNumTimeDim   = mTimeInterpolation->get_number_of_dimensions();

            // set member geometry type
            mTimeGeometryType = aInterpolationRule.get_time_geometry_type();

            // set pointers for second derivative depending on space and time dimensions
            this->set_function_pointers();
        }

        //------------------------------------------------------------------------------

        Geometry_Interpolator::Geometry_Interpolator(
                const mtk::Interpolation_Rule& aInterpolationRule,
                const mtk::Interpolation_Rule& aIPMapInterpolationRule,
                const mtk::CellShape&          aInterpolationShape,
                const bool                     aSpaceSideset,
                const bool                     aTimeSideset )
        {
            // set bool for side interpolation to true
            mSpaceSideset = aSpaceSideset;

            // set bool for time side interpolation to true
            mTimeSideset = aTimeSideset;

            // create member pointer to space interpolator
            mSpaceInterpolator = new mtk::Space_Interpolator(
                    aInterpolationRule,
                    aIPMapInterpolationRule,
                    aInterpolationShape,
                    aSpaceSideset );

            // getting mapping size and flag from space interpolator
            mMappedPoint = mSpaceInterpolator->get_initialized_mapped_point();
            mMapFlag     = mSpaceInterpolator->get_map_flag();

            // create member pointer to time interpolation function
            mTimeInterpolation = aInterpolationRule.create_time_interpolation_function();

            // number of space bases and dimensions
            mNumSpaceParamDim = mSpaceInterpolator->get_number_of_param_dimensions();

            // number of time bases and dimensions
            mNumTimeBases = mTimeInterpolation->get_number_of_bases();
            mNumTimeDim   = mTimeInterpolation->get_number_of_dimensions();

            // set member geometry type
            mTimeGeometryType = aInterpolationRule.get_time_geometry_type();

            // set pointers for second derivative depending on space and time dimensions
            this->set_function_pointers();
        }

        //------------------------------------------------------------------------------

        Geometry_Interpolator::~Geometry_Interpolator()
        {
            if ( mTimeInterpolation )
            {
                delete mTimeInterpolation;
                mTimeInterpolation = nullptr;
            }

            if ( mSpaceInterpolator )
            {
                delete mSpaceInterpolator;
                mSpaceInterpolator = nullptr;
            }
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::reset_eval_flags()
        {
            // reset booleans for evaluation
            //            mValtEval = true;
            mValt.reset();

            mNTauEval     = true;
            mdNdTauEval   = true;
            md2NdTau2Eval = true;
            md3NdTau3Eval = true;

            mTimeDetJEval   = true;
            mTimeJacEval    = true;
            mInvTimeJacEval = true;
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::reset_eval_flags_coordinates()
        {
            //            mValtEval = true;
            mValt.reset();
            mDeformedNodes.reset();

            mTimeDetJEval   = true;
            mTimeJacEval    = true;
            mInvTimeJacEval = true;
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::set_coeff(
                const Matrix< DDRMat >& aXHat,
                const Matrix< DDRMat >& aTHat )
        {
            // set the space coefficients
            mSpaceInterpolator->set_space_coeff( aXHat );

            // check the time coefficients input size
            MORIS_ASSERT( aTHat.n_rows() == mNumTimeBases,
                    " Geometry_Interpolator::set_coeff - Wrong input size (aTHat). " );

            // set the time coefficients
            mTHat = aTHat;

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::set_space_coeff( const Matrix< DDRMat >& aXHat )
        {
            // set the space coefficients
            mSpaceInterpolator->set_space_coeff( aXHat );

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::set_time_coeff( const Matrix< DDRMat >& aTHat )
        {
            // check the time coefficients input size
            //  fixme can not check the number of cols for aTHat
            MORIS_ASSERT( aTHat.n_rows() == mNumTimeBases,
                    " Geometry_Interpolator::set_time_coeff - Wrong input size (aTHat). " );

            // set the time coefficients
            mTHat = aTHat;

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::set_param_coeff()
        {
            // default implementation
            // set space and time param coords
            mSpaceInterpolator->set_param_coeff();

            mTimeInterpolation->get_param_coords( mTauHat );
            mTauHat = trans( mTauHat );

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::set_param_coeff(
                const Matrix< DDRMat >& aXiHat,
                const Matrix< DDRMat >& aTauHat )
        {
            // set the space coefficients
            mSpaceInterpolator->set_space_param_coeff( aXiHat );

            // check the time param coefficients input size
            MORIS_ASSERT( aTauHat.n_rows() == mNumTimeBases,
                    " Geometry_Interpolator::set_time_coeff - Wrong input size (aTauHat). " );

            // set the time coefficients
            mTauHat = aTauHat;

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::set_space_param_coeff( const Matrix< DDRMat >& aXiHat )
        {
            // set the space coefficients
            mSpaceInterpolator->set_space_param_coeff( aXiHat );

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::set_time_param_coeff( const Matrix< DDRMat >& aTauHat )
        {
            // check the time param coefficients input size
            //  fixme can not check the number of cols for aTauHat
            MORIS_ASSERT( aTauHat.n_rows() == mNumTimeBases,
                    " Geometry_Interpolator::set_time_coeff - Wrong input size (aTauHat). " );

            // set the time coefficients
            mTauHat = aTauHat;

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::set_space_time( const Matrix< DDRMat >& aParamPoint )
        {
            // check input size aParamPoint
            MORIS_ASSERT( ( ( aParamPoint.n_cols() == 1 ) && ( aParamPoint.n_rows() == mNumSpaceParamDim + mNumTimeDim ) ),
                    "Geometry_Interpolator::set_space_time - Wrong input size ( aParamPoint )." );

            // check input values are between -1 and 1
            // fixme what about TRI and TET
            for ( uint Ik = 0; Ik < mNumSpaceParamDim + mNumTimeDim; Ik++ )
            {
                MORIS_ASSERT( ( ( aParamPoint( Ik ) <= 1.0 + mEpsilon ) && ( aParamPoint( Ik ) >= -1.0 - mEpsilon ) ),
                        "Geometry_Interpolator::set_space_time - Wrong input value ( aParamPoint )." );
            }

            // set input values
            mSpaceInterpolator->set_space_time( aParamPoint );
            mTauLocal = aParamPoint( mNumSpaceParamDim );

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
        Geometry_Interpolator::set_space( const Matrix< DDRMat >& aSpaceParamPoint )
        {
            // set input values
            mSpaceInterpolator->set_space( aSpaceParamPoint );

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

        void
        Geometry_Interpolator::set_time( const Matrix< DDRMat >& aTimeParamPoint )
        {
            // check input size aParamPoint
            MORIS_ASSERT( ( ( aTimeParamPoint.n_cols() == 1 ) && ( aTimeParamPoint.n_rows() == mNumTimeDim ) ),
                    "Geometry_Interpolator::set_space - Wrong input size ( aTimeParamPoint )." );

            // check input values are between -1 and 1
            // fixme what about TRI and TET
            for ( uint Ik = 0; Ik < mNumTimeDim; Ik++ )
            {
                MORIS_ASSERT( ( ( aTimeParamPoint( Ik ) <= 1.0 + mEpsilon ) && ( aTimeParamPoint( Ik ) >= -1.0 - mEpsilon ) ),
                        "Geometry_Interpolator::set_time - Wrong input value ( aTimeParamPoint )." );
            }

            // set input values
            mTauLocal = aTimeParamPoint;

            // if no mapping required
            if ( !mMapFlag )
            {
                mMappedPoint( mNumSpaceParamDim ) = aTimeParamPoint( 0 );
            }

            // reset bool for evaluation
            this->reset_eval_flags();
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        real
        Geometry_Interpolator::get_time_step()
        {
            // check that mTHat is set
            MORIS_ASSERT( mTHat.numel() > 0, "Geometry_Interpolator::get_time_step - mTHat is not set." );

            // compute time increment deltat
            return mTHat.max() - mTHat.min();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::NXi()
        {
            // return member value
            return mSpaceInterpolator->NXi();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::NTau()
        {
            // if shape functions need to be evaluated
            if ( mNTauEval )
            {
                // evaluate the shape functions
                this->eval_NTau();

                // set bool for evaluation
                mNTauEval = false;
            }

            // return member value
            return mNTau;
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::eval_NTau()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mTauLocal.numel() > 0,
                    "Geometry_Interpolator::eval_NTau - mTauLocal is not set." );

            // pass data through interpolation function
            mTimeInterpolation->eval_N( mTauLocal, mNTau );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::dNdXi()
        {
            // return member value
            return mSpaceInterpolator->dNdXi();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::dNdTau()
        {
            // if shape functions need to be evaluated
            if ( mdNdTauEval )
            {
                // evaluate the shape functions 1st derivative
                this->eval_dNdTau();

                // set bool for evaluation
                mdNdTauEval = false;
            }

            // return member value
            return mdNdTau;
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::eval_dNdTau()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mTauLocal.numel() > 0,
                    "Geometry_Interpolator::eval_dNdTau - mTauLocal is not set." );

            // pass data through interpolation function
            mTimeInterpolation->eval_dNdXi( mTauLocal, mdNdTau );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::d2NdXi2()
        {
            // return member value
            return mSpaceInterpolator->d2NdXi2();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::d3NdXi3()
        {
            // return member value
            return mSpaceInterpolator->d3NdXi3();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::d2NdTau2()
        {
            // if shape functions need to be evaluated
            if ( md2NdTau2Eval )
            {
                // evaluate the shape functions 2nd derivative
                this->eval_d2NdTau2();

                // set bool for evaluation
                md2NdTau2Eval = false;
            }

            // return member value
            return md2NdTau2;
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::eval_d2NdTau2()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mTauLocal.numel() > 0,
                    "Geometry_Interpolator::eval_d2NdTau2 - mTauLocal is not set." );

            // pass data through interpolation function
            mTimeInterpolation->eval_d2NdXi2( mTauLocal, md2NdTau2 );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::space_jacobian()
        {
            // return member value
            return mSpaceInterpolator->space_jacobian();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::inverse_space_jacobian()
        {
            // return member value
            return mSpaceInterpolator->inverse_space_jacobian();
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::second_space_jacobian( Matrix< DDRMat >& aJ2bt )
        {
            // compute the second order Jacobian
            mSpaceInterpolator->second_space_jacobian( aJ2bt );
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::third_space_jacobian( Matrix< DDRMat >& aJ3ct )
        {
            // compute the third order Jacobian
            mSpaceInterpolator->third_space_jacobian( aJ3ct );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::time_jacobian()
        {
            // if space Jacobian needs to be evaluated
            if ( mTimeJacEval )
            {
                // evaluate the space Jacobian
                this->eval_time_jacobian();

                // set bool for evaluation
                mTimeJacEval = false;
            }

            // return member value
            return mTimeJac;
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::eval_time_jacobian()
        {
            // check that mTHat is set
            MORIS_ASSERT( mTHat.numel() > 0,
                    "Geometry_Interpolator::time_jacobian - mTHat is not set." );

            // compute the Jacobian
            mTimeJac = this->dNdTau() * mTHat;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::inverse_time_jacobian()
        {
            // if inverse of the time Jacobian needs to be evaluated
            if ( mInvTimeJacEval )
            {
                // evaluate the inverse of the time Jacobian
                this->eval_inverse_time_jacobian();

                // set bool for evaluation
                mInvTimeJacEval = false;
            }

            // return member value
            return mInvTimeJac;
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::eval_inverse_time_jacobian()
        {
            // compute inverse the time Jacobian
            ( this->*mInvTimeJacFunc )();
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::eval_inverse_time_jacobian_1d()
        {
            // get the space Jacobian
            const Matrix< DDRMat >& tTimeJac = this->time_jacobian();

            MORIS_ASSERT( tTimeJac( 0, 0 ) > sDetJInvJacLowerLimit,
                    "Time determinate (1D) close to zero or negative: %e\n",
                    tTimeJac( 0, 0 ) );

            mInvTimeJac.set_size( 1, 1 );

            mInvTimeJac = 1.0 / tTimeJac( 0, 0 );

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvTimeJac - inv( tTimeJac ) ) < 1e-8 * norm( mInvTimeJac ),
                    "Inconsistent time Jacobian (1D)" );
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::eval_inverse_time_jacobian_2d()
        {
            // get the space Jacobian
            const Matrix< DDRMat >& tTimeJac = this->time_jacobian();

            MORIS_ASSERT( this->time_det_J() > sDetJInvJacLowerLimit,
                    "Time determinate (2D) close to zero or negative: %e\n",
                    this->time_det_J() );

            // compute inverse of 2x2 matrix
            real tInvDet = 1.0 / ( this->time_det_J() );

            // compute inverse
            mInvTimeJac.set_size( 2, 2 );

            mInvTimeJac( 0, 0 ) = tTimeJac( 1, 1 ) * tInvDet;
            mInvTimeJac( 0, 1 ) = -tTimeJac( 0, 1 ) * tInvDet;
            mInvTimeJac( 1, 0 ) = -tTimeJac( 1, 0 ) * tInvDet;
            mInvTimeJac( 1, 1 ) = tTimeJac( 0, 0 ) * tInvDet;

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvTimeJac - inv( tTimeJac ) ) < 1e-8 * norm( mInvTimeJac ),
                    "Inconsistent time Jacobian (2D)" );
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::eval_inverse_time_jacobian_3d()
        {
            // get the space Jacobian
            const Matrix< DDRMat >& tTimeJac = this->time_jacobian();

            MORIS_ASSERT( this->time_det_J() > sDetJInvJacLowerLimit,
                    "Time determinate (3D) close to zero or negative: %e\n",
                    this->time_det_J() );

            // compute inverse of 3x3 matrix
            real tInvDet = 1.0 / ( this->time_det_J() );

            // compute inverse
            mInvTimeJac.set_size( 3, 3 );

            mInvTimeJac( 0, 0 ) = ( tTimeJac( 1, 1 ) * tTimeJac( 2, 2 ) - tTimeJac( 2, 1 ) * tTimeJac( 1, 2 ) ) * tInvDet;
            mInvTimeJac( 0, 1 ) = ( tTimeJac( 0, 2 ) * tTimeJac( 2, 1 ) - tTimeJac( 0, 1 ) * tTimeJac( 2, 2 ) ) * tInvDet;
            mInvTimeJac( 0, 2 ) = ( tTimeJac( 0, 1 ) * tTimeJac( 1, 2 ) - tTimeJac( 0, 2 ) * tTimeJac( 1, 1 ) ) * tInvDet;
            mInvTimeJac( 1, 0 ) = ( tTimeJac( 1, 2 ) * tTimeJac( 2, 0 ) - tTimeJac( 1, 0 ) * tTimeJac( 2, 2 ) ) * tInvDet;
            mInvTimeJac( 1, 1 ) = ( tTimeJac( 0, 0 ) * tTimeJac( 2, 2 ) - tTimeJac( 0, 2 ) * tTimeJac( 2, 0 ) ) * tInvDet;
            mInvTimeJac( 1, 2 ) = ( tTimeJac( 1, 0 ) * tTimeJac( 0, 2 ) - tTimeJac( 0, 0 ) * tTimeJac( 1, 2 ) ) * tInvDet;
            mInvTimeJac( 2, 0 ) = ( tTimeJac( 1, 0 ) * tTimeJac( 2, 1 ) - tTimeJac( 2, 0 ) * tTimeJac( 1, 1 ) ) * tInvDet;
            mInvTimeJac( 2, 1 ) = ( tTimeJac( 2, 0 ) * tTimeJac( 0, 1 ) - tTimeJac( 0, 0 ) * tTimeJac( 2, 1 ) ) * tInvDet;
            mInvTimeJac( 2, 2 ) = ( tTimeJac( 0, 0 ) * tTimeJac( 1, 1 ) - tTimeJac( 1, 0 ) * tTimeJac( 0, 1 ) ) * tInvDet;

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvTimeJac - inv( tTimeJac ) ) < 1e-8 * norm( mInvTimeJac ),
                    "Inconsistent time Jacobian (3D)" );
        }

        //------------------------------------------------------------------------------

        real
        Geometry_Interpolator::det_J()
        {
            // compute the determinant of the space time Jacobian
            return this->space_det_J() * this->time_det_J();
        }

        //------------------------------------------------------------------------------

        const real&
        Geometry_Interpolator::space_det_J()
        {
            // return member value
            return mSpaceInterpolator->space_det_J();
        }

        //------------------------------------------------------------------------------

        const real&
        Geometry_Interpolator::time_det_J()
        {
            // if determinant of time Jacobian needs to be evaluated
            if ( mTimeDetJEval )
            {
                // get the space Jacobian
                const Matrix< DDRMat >& tTimeJt = this->time_jacobian();

                // compute detJ for space
                mTimeDetJ = ( this->*mTimeDetJFunc )( tTimeJt );

                // set bool for evaluation
                mTimeDetJEval = false;
            }

            // return member value
            return mTimeDetJ;
        }

        //------------------------------------------------------------------------------

        real
        Geometry_Interpolator::eval_time_detJ_side(
                const Matrix< DDRMat >& aTimeJt )
        {
            return 1.0;
        }

        //------------------------------------------------------------------------------

        real
        Geometry_Interpolator::eval_time_detJ_bulk(
                const Matrix< DDRMat >& aTimeJt )
        {
            // FIXME: Needs to be specialized
            real tDetJ = det( aTimeJt );

            MORIS_ASSERT( tDetJ > sDetJLowerLimit,
                    "Time determinant (bulk) close to zero or negative: %e\n",
                    tDetJ );

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::get_normal( Matrix< DDRMat >& aNormal )
        {
            // call space interpolator
            mSpaceInterpolator->get_normal( aNormal );
        }

        Matrix< DDRMat >
        Geometry_Interpolator::get_normal()
        {
            // call space interpolator
            Matrix< DDRMat > tNormal;
            mSpaceInterpolator->get_normal( tNormal );
            return tNormal;
        }

        Matrix< DDRMat >
        Geometry_Interpolator::get_normal_current( Field_Interpolator* aFieldInterpolator )
        {
            Matrix< DDRMat > const tPreviousSpaceCoeffs = mSpaceInterpolator->get_space_coeff();           // store the previous state
            mSpaceInterpolator->set_space_coeff( this->get_space_coeff_current( aFieldInterpolator ) );    // set the deformed coefficients
            Matrix< DDRMat > tNormal;
            mSpaceInterpolator->get_normal( tNormal );
            mSpaceInterpolator->set_space_coeff( tPreviousSpaceCoeffs );    // reset the space interpolator to the previous state
            return tNormal;
        }


        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::valx()
        {
            // call space interpolator
            return mSpaceInterpolator->valx();
        }

        Matrix< DDRMat > Geometry_Interpolator::valx_current( Field_Interpolator* aFieldInterpolator )
        {
            // get the current parametric space coordinate
            Matrix< DDRMat > const tPreviousSpaceCoeffs = mSpaceInterpolator->get_space_coeff();
            mSpaceInterpolator->set_space_coeff( this->get_space_coeff_current( aFieldInterpolator ) );
            Matrix< DDRMat > const tXCurrent = mSpaceInterpolator->valx();
            mSpaceInterpolator->set_space_coeff( tPreviousSpaceCoeffs );
            return tXCurrent;
        }

        const Matrix< DDRMat >&
        Geometry_Interpolator::get_space_coeff_current( Field_Interpolator* aFieldInterpolator )
        {
            // calculate the current coordinates of the element nodes (if not cached)
            if ( !mDeformedNodes.has_value() )
            {
                // Since the internal state of the field interpolator is changed for each evaluation of the element node,
                // we have to store the previous state and restore it after the evaluation
                Matrix< DDRMat > const tPreviousFIPoint = aFieldInterpolator->get_space_time();

                Matrix< DDRMat > const tXhat  = mSpaceInterpolator->get_space_coeff();          // the physical coordinates of the element nodes (n_nodes x n_dim)
                Matrix< DDRMat > const tXiHat = mSpaceInterpolator->get_space_param_coeff();    // the parametric coordinates of the element nodes (n_nodes x n_dim)
                mDeformedNodes                = Matrix< DDRMat >( tXhat.n_rows(), tXhat.n_cols() );
                for ( size_t iNode = 0; iNode < tXhat.n_rows(); ++iNode )
                {
                    // Set the Field Interpolator to the parametric coordinates of the current node to get the value of the displacement field for this node.
                    Matrix< DDRMat > tSpaceTime( tXiHat.n_cols() + 1, 1 );
                    tSpaceTime( { 0, tXiHat.n_cols() - 1 }, { 0, 0 } ) = trans( tXiHat.get_row( iNode ) );
                    aFieldInterpolator->set_space_time( tSpaceTime );
                    mDeformedNodes->set_row( iNode, tXhat.get_row( iNode ) + trans( aFieldInterpolator->val() ) );
                }
                // Reset the Field Interpolator to the previous space coordinates.
                aFieldInterpolator->set_space_time( tPreviousFIPoint );
            }
            return mDeformedNodes.value();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::valt()
        {
            if ( !mValt.has_value() )
            {
                // check that mTHat is set
                MORIS_ASSERT( mTHat.numel() > 0,
                        "Geometry_Interpolator::valt - mTHat is not set." );

                // evaluate the field
                mValt = this->NTau() * mTHat;
            }

            return mValt.value();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::map_integration_point()
        {
            // if eval mapping
            if ( mMapFlag )
            {
                MORIS_ASSERT( mTauHat.numel() > 0,
                        "Geometry_Interpolator::eval_mapping - mTauHat is not set." );

                uint tSize = mMappedPoint.numel() - 1;

                // set mapped space coordinates
                mMappedPoint = mSpaceInterpolator->map_integration_point();

                // set mapped time coordinates
                mMappedPoint( tSize ) = dot( this->NTau(), mTauHat );
            }

            return mMappedPoint;
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::update_parametric_coordinates(
                Matrix< DDRMat > const & aPhysCoordinates,
                Matrix< DDRMat >&        aParamCoordinates ) const
        {
            // call space interpolator
            mSpaceInterpolator->update_parametric_coordinates( aPhysCoordinates, aParamCoordinates );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Geometry_Interpolator::metric_tensor()
        {
            // call space interpolator
            return mSpaceInterpolator->metric_tensor();
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::space_jacobian_and_matrices_for_second_derivatives(
                Matrix< DDRMat >&       aJt,
                Matrix< DDRMat >&       aKt,
                Matrix< DDRMat >&       aLt,
                const Matrix< DDRMat >& adNdXi,
                const Matrix< DDRMat >& ad2NdXi2 )
        {
            // call space interpolator
            mSpaceInterpolator->space_jacobian_and_matrices_for_second_derivatives(
                    aJt,
                    aKt,
                    aLt,
                    adNdXi,
                    ad2NdXi2 );
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::space_jacobian_and_matrices_for_third_derivatives(
                Matrix< DDRMat >&       aJt,      // contains first geometric derivs
                Matrix< DDRMat >&       aJ2bt,    // contains second geometric derivs = second help matrix for 2nd field derivs
                Matrix< DDRMat >&       aJ3at,    // first help matrix for 3rd field derivs
                Matrix< DDRMat >&       aJ3bt,    // second help matrix for 3rd field derivs
                Matrix< DDRMat >&       aJ3ct,    // third help matrix for 3rd field derivs
                const Matrix< DDRMat >& adNdXi,
                const Matrix< DDRMat >& ad2NdXi2,
                const Matrix< DDRMat >& ad3NdXi3 )
        {
            // call space interpolator
            mSpaceInterpolator->space_jacobian_and_matrices_for_third_derivatives(
                    aJt,
                    aJ2bt,
                    aJ3at,
                    aJ3bt,
                    aJ3ct,
                    adNdXi,
                    ad2NdXi2,
                    ad3NdXi3 );
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::time_jacobian_and_matrices_for_second_derivatives(
                Matrix< DDRMat >&       aJt,
                Matrix< DDRMat >&       aKt,
                Matrix< DDRMat >&       aLt,
                const Matrix< DDRMat >& adNdTau,
                const Matrix< DDRMat >& ad2NdTau2 )
        {
            // check that mTHat is set
            MORIS_ASSERT( mTHat.numel() > 0,
                    "Geometry_Interpolator::time_jacobian_and_matrices_for_second_derivatives - mTHat is not set." );

            // evaluate transposed of geometry Jacobian
            aJt = this->time_jacobian();

            // call calculator for second derivatives
            this->mSecondDerivativeMatricesTime(
                    aJt,
                    aKt,
                    aLt,
                    ad2NdTau2,
                    mTHat );
        }

        //------------------------------------------------------------------------------

        void
        Geometry_Interpolator::set_function_pointers()
        {
            // get number of dimensions and set pointer to function
            // for second time derivative
            switch ( mNumTimeDim )
            {
                case 1:
                {
                    mInvTimeJacFunc               = &Geometry_Interpolator::eval_inverse_time_jacobian_1d;
                    mSecondDerivativeMatricesTime = mSpaceInterpolator->eval_matrices_for_second_derivative_1d;
                    break;
                }
                case 2:
                {
                    mInvTimeJacFunc               = &Geometry_Interpolator::eval_inverse_time_jacobian_2d;
                    mSecondDerivativeMatricesTime = mSpaceInterpolator->eval_matrices_for_second_derivative_2d;
                    break;
                }
                case 3:
                {
                    mInvTimeJacFunc               = &Geometry_Interpolator::eval_inverse_time_jacobian_3d;
                    mSecondDerivativeMatricesTime = mSpaceInterpolator->eval_matrices_for_second_derivative_3d;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, " Geometry_Interpolator::set_function_pointers - unknown number of dimensions. " );
                }
            }

            // switch Geometry_Type
            switch ( mTimeGeometryType )
            {
                case mtk::Geometry_Type::POINT:
                {
                    mMapFlag      = true;
                    mTimeDetJFunc = &Geometry_Interpolator::eval_time_detJ_side;
                    break;
                }
                case mtk::Geometry_Type::LINE:
                {
                    mTimeDetJFunc = &Geometry_Interpolator::eval_time_detJ_bulk;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Geometry_Interpolator::set_function_pointers - unknown or not implemented time geometry type " );
                }
            }
        }
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#include "fn_norm.hpp"
#include "fn_cross.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"
#include "fn_inv.hpp"
#include "op_div.hpp"
#include "fn_linsolve.hpp"

#include "cl_FEM_Geometry_Interpolator.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        Geometry_Interpolator::Geometry_Interpolator(
                const Interpolation_Rule & aInterpolationRule,
                const bool                 aSpaceSideset,
                const bool                 aTimeSideset )
        {
            // set bool for side interpolation to true
            mSpaceSideset = aSpaceSideset;

            // set bool for time side interpolation to true
            mTimeSideset = aTimeSideset;

            // create member pointer to space interpolation function
            mSpaceInterpolation = aInterpolationRule.create_space_interpolation_function();

            // create member pointer to time interpolation function
            mTimeInterpolation  = aInterpolationRule.create_time_interpolation_function();

            // number of space bases and dimensions
            mNumSpaceBases    = mSpaceInterpolation->get_number_of_bases();
            mNumSpaceDim      = mSpaceInterpolation->get_number_of_dimensions();
            mNumSpaceParamDim = mSpaceInterpolation->get_number_of_param_dimensions();

            // number of time bases and dimensions
            mNumTimeBases = mTimeInterpolation->get_number_of_bases();
            mNumTimeDim   = mTimeInterpolation->get_number_of_dimensions();

            // set member geometry type
            mGeometryType     = aInterpolationRule.get_geometry_type();
            mTimeGeometryType = aInterpolationRule.get_time_geometry_type();

            // set pointers for second derivative depending on space and time dimensions
            this->set_function_pointers();
        }

        //------------------------------------------------------------------------------

        Geometry_Interpolator::~Geometry_Interpolator()
        {
            // delete interpolation functions
            if( mSpaceInterpolation != NULL )
            {
                delete mSpaceInterpolation;
            }

            if( mTimeInterpolation != NULL )
            {
                delete mTimeInterpolation;
            }
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::reset_eval_flags()
        {
            // reset bool for evaluation
            mNXiEval      = true;
            mdNdXiEval    = true;
            md2NdXi2Eval  = true;
            md3NdXi3Eval  = true;
            mNTauEval     = true;
            mdNdTauEval   = true;
            md2NdTau2Eval = true;
            md3NdTau3Eval = true;

            mSpaceDetJEval   = true;
            mSpaceJacEval    = true;
            mInvSpaceJacEval = true;

            mTimeDetJEval    = true;
            mTimeJacEval     = true;
            mInvTimeJacEval  = true;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::reset_eval_flags_coordinates()
        {
            mSpaceDetJEval   = true;
            mSpaceJacEval    = true;
            mInvSpaceJacEval = true;

            mTimeDetJEval    = true;
            mTimeJacEval     = true;
            mInvTimeJacEval  = true;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_coeff(
                const Matrix< DDRMat > & aXHat,
                const Matrix< DDRMat > & aTHat )
        {
            // check the space coefficients input size
            // fixme can not check the number of cols for aXHat and aTHat
            MORIS_ASSERT( aXHat.n_rows() == mNumSpaceBases ,
                    " Geometry_Interpolator::set_coeff - Wrong input size (aXHat). ");

            // set the space coefficients
            mXHat = aXHat;

            //check the time coefficients input size
            MORIS_ASSERT( aTHat.n_rows() == mNumTimeBases,
                    " Geometry_Interpolator::set_coeff - Wrong input size (aTHat). ");

            // set the time coefficients
            mTHat = aTHat;

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_space_coeff( const Matrix< DDRMat > & aXHat )
        {
            //check the space coefficients input size
            // fixme can not check the number of cols for aXHat
            MORIS_ASSERT( aXHat.n_rows() == mNumSpaceBases,
                    " Geometry_Interpolator::set_space_coeff - Wrong input size (aXHat). ");

            // set the space coefficients
            mXHat = aXHat;

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_time_coeff( const Matrix< DDRMat > & aTHat )
        {
            //check the time coefficients input size
            // fixme can not check the number of cols for aTHat
            MORIS_ASSERT( aTHat.n_rows() == mNumTimeBases,
                    " Geometry_Interpolator::set_time_coeff - Wrong input size (aTHat). ");

            // set the time coefficients
            mTHat = aTHat;

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_param_coeff()
        {
            // default implementation
            // set space and time param coords
            mSpaceInterpolation->get_param_coords( mXiHat );
            mXiHat = trans( mXiHat );

            mTimeInterpolation->get_param_coords( mTauHat );
            mTauHat = trans( mTauHat );

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_param_coeff(
                const Matrix< DDRMat > & aXiHat,
                const Matrix< DDRMat > & aTauHat )
        {
            //check the space param coefficients input size
            MORIS_ASSERT( aXiHat.n_rows() == mNumSpaceBases,
                    " Geometry_Interpolator::set_space_param_coeff - Wrong input size (aXiHat). ");

            // set the space coefficients
            mXiHat = aXiHat;

            //check the time param coefficients input size
            MORIS_ASSERT( aTauHat.n_rows() == mNumTimeBases,
                    " Geometry_Interpolator::set_time_coeff - Wrong input size (aTauHat). ");

            // set the time coefficients
            mTauHat = aTauHat;

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_space_param_coeff( const Matrix< DDRMat > & aXiHat )
        {
            //check the space param coefficients input size
            // fixme can not check the number of cols for aXiHat
            MORIS_ASSERT( aXiHat.n_rows() == mNumSpaceBases,
                    " Geometry_Interpolator::set_space_param_coeff - Wrong input size (aXiHat). ");

            // set the space coefficients
            mXiHat = aXiHat;

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_time_param_coeff( const Matrix< DDRMat > & aTauHat )
        {
            //check the time param coefficients input size
            // fixme can not check the number of cols for aTauHat
            MORIS_ASSERT( aTauHat.n_rows() == mNumTimeBases,
                    " Geometry_Interpolator::set_time_coeff - Wrong input size (aTauHat). ");

            // set the time coefficients
            mTauHat = aTauHat;

            // reset evaluation flags
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_space_time( const Matrix< DDRMat > & aParamPoint )
        {
            // check input size aParamPoint
            MORIS_ASSERT( ( ( aParamPoint.n_cols() == 1 ) && ( aParamPoint.n_rows() == mNumSpaceParamDim + mNumTimeDim )),
                    "Geometry_Interpolator::set_space_time - Wrong input size ( aParamPoint ).");

            // check input values are between -1 and 1
            // fixme what about TRI and TET
            for ( uint Ik = 0; Ik < mNumSpaceParamDim + mNumTimeDim; Ik++ )
            {
                MORIS_ASSERT( ( ( aParamPoint( Ik ) <= 1.0 + mEpsilon ) && ( aParamPoint( Ik ) >= -1.0 - mEpsilon ) ),
                        "Geometry_Interpolator::set_space_time - Wrong input value ( aParamPoint ).");
            }

            // set input values
            mXiLocal  = aParamPoint( { 0, mNumSpaceParamDim-1 }, { 0, 0 } );
            mTauLocal = aParamPoint( mNumSpaceParamDim );

            // if no mapping required
            if( !mMapFlag )
            {
                mMappedPoint = aParamPoint;
            }

            // reset bool for evaluation
            this->reset_eval_flags();
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_space( const Matrix< DDRMat > & aSpaceParamPoint )
        {
            // check input size aParamPoint
            MORIS_ASSERT( ( ( aSpaceParamPoint.n_cols() == 1 ) && ( aSpaceParamPoint.n_rows() == mNumSpaceParamDim )),
                    "Geometry_Interpolator::set_space - Wrong input size ( aSpaceParamPoint ).");

            // check input values are between -1 and 1
            // fixme what about TRI and TET
            for ( uint Ik = 0; Ik < mNumSpaceParamDim; Ik++ )
            {
                MORIS_ASSERT( ( ( aSpaceParamPoint( Ik ) <= 1.0 + mEpsilon ) && ( aSpaceParamPoint( Ik ) >= -1.0 - mEpsilon ) ),
                        "Geometry_Interpolator::set_space - Wrong input value ( aSpaceParamPoint ).");
            }

            // set input values
            mXiLocal  = aSpaceParamPoint;

            // if no mapping required
            if( !mMapFlag )
            {
                mMappedPoint( { 0, mNumSpaceParamDim-1 }, { 0, 0 } ) =
                        aSpaceParamPoint.matrix_data();
            }

            // reset bool for evaluation
            this->reset_eval_flags();
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_time( const Matrix< DDRMat > & aTimeParamPoint )
        {
            // check input size aParamPoint
            MORIS_ASSERT( ( ( aTimeParamPoint.n_cols() == 1 ) && ( aTimeParamPoint.n_rows() == mNumTimeDim )),
                    "Geometry_Interpolator::set_space - Wrong input size ( aTimeParamPoint ).");

            // check input values are between -1 and 1
            // fixme what about TRI and TET
            for ( uint Ik = 0; Ik < mNumTimeDim; Ik++ )
            {
                MORIS_ASSERT( ( ( aTimeParamPoint( Ik ) <= 1.0 + mEpsilon ) && ( aTimeParamPoint( Ik ) >= -1.0 - mEpsilon ) ),
                        "Geometry_Interpolator::set_time - Wrong input value ( aTimeParamPoint ).");
            }

            // set input values
            mTauLocal  = aTimeParamPoint;

            // if no mapping required
            if( !mMapFlag )
            {
                mMappedPoint( mNumSpaceParamDim ) = aTimeParamPoint( 0 );
            }

            // reset bool for evaluation
            this->reset_eval_flags();
            this->reset_eval_flags_coordinates();
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::get_time_step()
        {
            // check that mTHat is set
            MORIS_ASSERT( mTHat.numel()>0, "Geometry_Interpolator::get_time_step - mTHat is not set." );

            // compute time increment deltat
            return mTHat.max() - mTHat.min();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::NXi()
        {
            // if shape functions need to be evaluated
            if( mNXiEval )
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

        void Geometry_Interpolator::eval_NXi()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Geometry_Interpolator::eval_NXi - mXiLocal is not set." );

            // pass data through interpolation function
            mSpaceInterpolation->eval_N( mXiLocal, mNXi );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::NTau()
        {
            // if shape functions need to be evaluated
            if( mNTauEval )
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

        void Geometry_Interpolator::eval_NTau()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mTauLocal.numel() > 0,
                    "Geometry_Interpolator::eval_NTau - mTauLocal is not set." );

            // pass data through interpolation function
            mTimeInterpolation->eval_N( mTauLocal, mNTau );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::dNdXi()
        {
            // if shape functions need to be evaluated
            if( mdNdXiEval )
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

        void Geometry_Interpolator::eval_dNdXi()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Geometry_Interpolator::eval_dNdXi - mXiLocal is not set." );

            // pass data through interpolation function
            mSpaceInterpolation->eval_dNdXi( mXiLocal, mdNdXi );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::dNdTau()
        {
            // if shape functions need to be evaluated
            if( mdNdTauEval )
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

        void Geometry_Interpolator::eval_dNdTau()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mTauLocal.numel() > 0,
                    "Geometry_Interpolator::eval_dNdTau - mTauLocal is not set." );

            // pass data through interpolation function
            mTimeInterpolation->eval_dNdXi( mTauLocal, mdNdTau );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::d2NdXi2()
        {
            // if shape functions need to be evaluated
            if( md2NdXi2Eval )
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

        void Geometry_Interpolator::eval_d2NdXi2()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Geometry_Interpolator::eval_d2NdXi2 - mXiLocal is not set." );

            // pass data through interpolation function
            mSpaceInterpolation->eval_d2NdXi2( mXiLocal, md2NdXi2 );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::d3NdXi3()
        {
            // if shape functions need to be evaluated
            if( md3NdXi3Eval )
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

        void Geometry_Interpolator::eval_d3NdXi3()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Geometry_Interpolator::eval_d3NdXi3 - mXiLocal is not set." );

            // pass data through interpolation function
            mSpaceInterpolation->eval_d3NdXi3( mXiLocal, md3NdXi3 );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::d2NdTau2()
        {
            // if shape functions need to be evaluated
            if( md2NdTau2Eval )
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

        void Geometry_Interpolator::eval_d2NdTau2()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mTauLocal.numel() > 0,
                    "Geometry_Interpolator::eval_d2NdTau2 - mTauLocal is not set." );

            // pass data through interpolation function
            mTimeInterpolation->eval_d2NdXi2( mTauLocal, md2NdTau2 );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::space_jacobian()
        {
            // if space Jacobian needs to be evaluated
            if( mSpaceJacEval )
            {
                // evaluate the space Jacobian
                this->eval_space_jacobian();

                // set bool for evaluation
                mSpaceJacEval = false;
            }

            // return member value
            return mSpaceJac;
        }

        void Geometry_Interpolator::eval_space_jacobian()
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel() > 0,
                    "Geometry_Interpolator::space_jacobian - mXHat is not set." );

            // compute the Jacobian
            mSpaceJac = this->dNdXi() * mXHat;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::inverse_space_jacobian()
        {
            // if inverse of the space Jacobian needs to be evaluated
            if( mInvSpaceJacEval )
            {
                // evaluate the inverse of the space Jacobian
                this->eval_inverse_space_jacobian();

                // set bool for evaluation
                mInvSpaceJacEval = false;
            }

            // return member value
            return mInvSpaceJac;
        }

        void Geometry_Interpolator::eval_inverse_space_jacobian()
        {
            // compute the Jacobian
            (this->*mInvSpaceJacFunc)();
        }

        void Geometry_Interpolator::eval_inverse_space_jacobian_1d()
        {
            // get the space Jacobian
            const Matrix< DDRMat > & tSpacJac = this->space_jacobian();

            MORIS_ASSERT(tSpacJac(0,0) > 1e-18,
                    "Space determinate (1D) close to zero or negative: %e\n", tSpacJac(0,0));

            mInvSpaceJac.set_size(1,1);

            mInvSpaceJac(0,0) = 1.0/tSpacJac(0,0);

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvSpaceJac-inv( tSpacJac ) ) < 1e-8*norm(mInvSpaceJac ),
                    "Inconsistent space Jacobian (1D)\n");
        }

        void Geometry_Interpolator::eval_inverse_space_jacobian_2d()
        {
            // get the space Jacobian
            const Matrix< DDRMat > & tSpacJac = this->space_jacobian();

            // compute inverse of 3x3 matrix
            real tInvDet = 1.0/(this->space_det_J());

            MORIS_ASSERT(std::abs(tInvDet) < 1e+18,
                    "Space determinate (2D) close to zero or negative: %e\n", this->space_det_J());

            // compute inverse
            mInvSpaceJac.set_size(2,2);

            mInvSpaceJac(0, 0) =  tSpacJac(1, 1) * tInvDet;
            mInvSpaceJac(0, 1) = -tSpacJac(0, 1) * tInvDet;
            mInvSpaceJac(1, 0) = -tSpacJac(1, 0) * tInvDet;
            mInvSpaceJac(1, 1) =  tSpacJac(0, 0) * tInvDet;

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvSpaceJac-inv( tSpacJac ) ) < 1e-8*norm(mInvSpaceJac ),
                    "Inconsistent space Jacobian (2D)");
        }

        void Geometry_Interpolator::eval_inverse_space_jacobian_3d()
        {
            // get the space Jacobian
            const Matrix< DDRMat > & tSpacJac = this->space_jacobian();

            // compute inverse of 3x3 matrix
            real tInvDet = 1.0/(this->space_det_J());

            MORIS_ASSERT(std::abs(tInvDet) < 1e+18,
                    "Space determinate (3D) close to zero or negative: %e\n", this->space_det_J());

            // compute inverse
            mInvSpaceJac.set_size(3,3);

            mInvSpaceJac(0, 0) = (tSpacJac(1, 1) * tSpacJac(2, 2) - tSpacJac(2, 1) * tSpacJac(1, 2)) * tInvDet;
            mInvSpaceJac(0, 1) = (tSpacJac(0, 2) * tSpacJac(2, 1) - tSpacJac(0, 1) * tSpacJac(2, 2)) * tInvDet;
            mInvSpaceJac(0, 2) = (tSpacJac(0, 1) * tSpacJac(1, 2) - tSpacJac(0, 2) * tSpacJac(1, 1)) * tInvDet;
            mInvSpaceJac(1, 0) = (tSpacJac(1, 2) * tSpacJac(2, 0) - tSpacJac(1, 0) * tSpacJac(2, 2)) * tInvDet;
            mInvSpaceJac(1, 1) = (tSpacJac(0, 0) * tSpacJac(2, 2) - tSpacJac(0, 2) * tSpacJac(2, 0)) * tInvDet;
            mInvSpaceJac(1, 2) = (tSpacJac(1, 0) * tSpacJac(0, 2) - tSpacJac(0, 0) * tSpacJac(1, 2)) * tInvDet;
            mInvSpaceJac(2, 0) = (tSpacJac(1, 0) * tSpacJac(2, 1) - tSpacJac(2, 0) * tSpacJac(1, 1)) * tInvDet;
            mInvSpaceJac(2, 1) = (tSpacJac(2, 0) * tSpacJac(0, 1) - tSpacJac(0, 0) * tSpacJac(2, 1)) * tInvDet;
            mInvSpaceJac(2, 2) = (tSpacJac(0, 0) * tSpacJac(1, 1) - tSpacJac(1, 0) * tSpacJac(0, 1)) * tInvDet;

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvSpaceJac-inv( tSpacJac ) ) < 1e-8*norm(mInvSpaceJac ),
                    "Inconsistent space Jacobian (3D)");
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::second_space_jacobian( Matrix< DDRMat > & aJ2bt )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel()>0,
                    "Geometry_Interpolator::second_space_jacobian - mXHat is not set." );

            // compute the second order Jacobian
            aJ2bt = this->d2NdXi2() * mXHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::third_space_jacobian( Matrix< DDRMat > & aJ3ct )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel()>0,
                    "Geometry_Interpolator::third_space_jacobian - mXHat is not set." );

            // compute the third order Jacobian
            aJ3ct = this->d3NdXi3() * mXHat;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::time_jacobian()
        {
            // if space Jacobian needs to be evaluated
            if( mTimeJacEval )
            {
                // evaluate the space Jacobian
                this->eval_time_jacobian();

                // set bool for evaluation
                mTimeJacEval = false;
            }

            // return member value
            return mTimeJac;
        }

        void Geometry_Interpolator::eval_time_jacobian()
        {
            // check that mTHat is set
            MORIS_ASSERT( mTHat.numel()>0,
                    "Geometry_Interpolator::time_jacobian - mTHat is not set." );

            // compute the Jacobian
            mTimeJac = this->dNdTau() * mTHat;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::inverse_time_jacobian()
        {
            // if inverse of the time Jacobian needs to be evaluated
            if( mInvTimeJacEval )
            {
                // evaluate the inverse of the time Jacobian
                this->eval_inverse_time_jacobian();

                // set bool for evaluation
                mInvTimeJacEval = false;
            }

            // return member value
            return mInvTimeJac;
        }

        void Geometry_Interpolator::eval_inverse_time_jacobian()
        {
            // compute inverse the time Jacobian
            (this->*mInvTimeJacFunc)();
        }

        void Geometry_Interpolator::eval_inverse_time_jacobian_1d()
        {
            // get the space Jacobian
            const Matrix< DDRMat > & tTimeJac = this->time_jacobian();

            MORIS_ASSERT(tTimeJac(0,0) > 1e-18,
                    "Time determinate (1D) close to zero or negative: %e\n", tTimeJac(0,0));

            mInvTimeJac.set_size(1,1);

            mInvTimeJac = 1.0/tTimeJac(0,0);

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvTimeJac-inv( tTimeJac ) ) < 1e-8*norm(mInvTimeJac ),
                    "Inconsistent time Jacobian (1D)");
        }

        void Geometry_Interpolator::eval_inverse_time_jacobian_2d()
        {
            // get the space Jacobian
            const Matrix< DDRMat > & tTimeJac = this->time_jacobian();

            // compute inverse of 2x2 matrix
            real tInvDet = 1.0/(this->time_det_J());

            MORIS_ASSERT(std::abs(tInvDet) < 1e+18,
                    "Time determinate (2D) close to zero or negative: %e\n", this->time_det_J());

            // compute inverse
            mInvTimeJac.set_size(2,2);

            mInvTimeJac(0, 0) =  tTimeJac(1, 1) * tInvDet;
            mInvTimeJac(0, 1) = -tTimeJac(0, 1) * tInvDet;
            mInvTimeJac(1, 0) = -tTimeJac(1, 0) * tInvDet;
            mInvTimeJac(1, 1) =  tTimeJac(0, 0) * tInvDet;

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvTimeJac-inv( tTimeJac ) ) < 1e-8*norm(mInvTimeJac ),
                    "Inconsistent time Jacobian (2D)");
        }

        void Geometry_Interpolator::eval_inverse_time_jacobian_3d()
        {
            // get the space Jacobian
            const Matrix< DDRMat > & tTimeJac = this->time_jacobian();

            // compute inverse of 3x3 matrix
            real tInvDet = 1.0/(this->time_det_J());

            MORIS_ASSERT(std::abs(tInvDet) < 1e+18,
                    "Time determinate (3D) close to zero or negative: %e\n", this->time_det_J());

            // compute inverse
            mInvTimeJac.set_size(3,3);

            mInvTimeJac(0, 0) = (tTimeJac(1, 1) * tTimeJac(2, 2) - tTimeJac(2, 1) * tTimeJac(1, 2)) * tInvDet;
            mInvTimeJac(0, 1) = (tTimeJac(0, 2) * tTimeJac(2, 1) - tTimeJac(0, 1) * tTimeJac(2, 2)) * tInvDet;
            mInvTimeJac(0, 2) = (tTimeJac(0, 1) * tTimeJac(1, 2) - tTimeJac(0, 2) * tTimeJac(1, 1)) * tInvDet;
            mInvTimeJac(1, 0) = (tTimeJac(1, 2) * tTimeJac(2, 0) - tTimeJac(1, 0) * tTimeJac(2, 2)) * tInvDet;
            mInvTimeJac(1, 1) = (tTimeJac(0, 0) * tTimeJac(2, 2) - tTimeJac(0, 2) * tTimeJac(2, 0)) * tInvDet;
            mInvTimeJac(1, 2) = (tTimeJac(1, 0) * tTimeJac(0, 2) - tTimeJac(0, 0) * tTimeJac(1, 2)) * tInvDet;
            mInvTimeJac(2, 0) = (tTimeJac(1, 0) * tTimeJac(2, 1) - tTimeJac(2, 0) * tTimeJac(1, 1)) * tInvDet;
            mInvTimeJac(2, 1) = (tTimeJac(2, 0) * tTimeJac(0, 1) - tTimeJac(0, 0) * tTimeJac(2, 1)) * tInvDet;
            mInvTimeJac(2, 2) = (tTimeJac(0, 0) * tTimeJac(1, 1) - tTimeJac(1, 0) * tTimeJac(0, 1)) * tInvDet;

            // check results against generic inverse operator
            MORIS_ASSERT( norm( mInvTimeJac-inv( tTimeJac ) ) < 1e-8*norm(mInvTimeJac ),
                    "Inconsistent time Jacobian (3D)");
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::det_J()
        {
            // compute the determinant of the space time Jacobian
            return this->space_det_J() * this->time_det_J();
        }

        //------------------------------------------------------------------------------

        const real & Geometry_Interpolator::space_det_J()
        {
            // if determinant of space Jacobian needs to be evaluated
            if( mSpaceDetJEval )
            {
                // get the space Jacobian
                const Matrix< DDRMat > & tSpaceJt = this->space_jacobian();

                // compute detJ for space
                mSpaceDetJ = ( this->*mSpaceDetJFunc )( tSpaceJt );

                // set bool for evaluation
                mSpaceDetJEval = false;
            }

            // return member value
            return mSpaceDetJ;
        }

        //------------------------------------------------------------------------------

        const real & Geometry_Interpolator::time_det_J()
        {
            // if determinant of time Jacobian needs to be evaluated
            if( mTimeDetJEval )
            {
                // get the space Jacobian
                const Matrix< DDRMat > & tTimeJt = this->time_jacobian();

                // compute detJ for space
                mTimeDetJ = ( this->*mTimeDetJFunc )( tTimeJt );

                // set bool for evaluation
                mTimeDetJEval = false;
            }

            // return member value
            return mTimeDetJ;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_space_detJ_side_line(
                const Matrix< DDRMat > & aSpaceJt )
        {
            real tDetJ = norm( aSpaceJt );

            MORIS_ASSERT( tDetJ > 1e-18,
                    "Space determinant (side line) close to zero or negative: %e\n", tDetJ);

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_space_detJ_side_tri(
                const Matrix< DDRMat > & aSpaceJt )
        {
            real tDetJ =
                    norm( cross(
                            aSpaceJt.get_row( 0 ) - aSpaceJt.get_row( 2 ),
                            aSpaceJt.get_row( 1 ) - aSpaceJt.get_row( 2 ) ) ) / 2.0;

            MORIS_ASSERT( tDetJ > 1e-18,
                    "Space determinant (side tri) close to zero or negative: %e\n", tDetJ);

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_space_detJ_side_quad(
                const Matrix< DDRMat > & aSpaceJt )
        {
            real tDetJ = norm( cross( aSpaceJt.get_row( 0 ), aSpaceJt.get_row( 1 ) ) );

            MORIS_ASSERT( tDetJ > 1e-18,
                    "Space determinant (side quad) close to zero or negative: %e\n", tDetJ);

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_space_detJ_bulk_line(
                const Matrix< DDRMat > & aSpaceJt )
        {
            real tDetJ = aSpaceJt(0,0);

            MORIS_ASSERT( tDetJ > 1e-18,
                    "Space determinant (bulk 1D) close to zero or negative: %e\n", tDetJ);

            MORIS_ASSERT( std::abs( det(aSpaceJt) - tDetJ ) < 1e-8*tDetJ,
                    "Inconsistent space determinant (bulk 1D): %e vs %e\n",
                    tDetJ,det(aSpaceJt));

            return tDetJ;
        }

        real Geometry_Interpolator::eval_space_detJ_bulk_quad(
                const Matrix< DDRMat > & aSpaceJt )
        {
            real tDetJ = aSpaceJt(0,0)*aSpaceJt(1,1)-aSpaceJt(0,1)*aSpaceJt(1,0);

            MORIS_ASSERT( tDetJ > 1e-18,
                    "Space determinant (bulk 2D) close to zero or negative: %e\n", tDetJ);

            MORIS_ASSERT( std::abs( det(aSpaceJt) - tDetJ ) < 1e-8*tDetJ,
                    "Inconsistent space determinant (bulk 2D): %e vs %e\n",
                    tDetJ,det(aSpaceJt));

            return tDetJ;
        }

        real Geometry_Interpolator::eval_space_detJ_bulk_hex(
                const Matrix< DDRMat > & aSpaceJt )
        {
            real tDetJ =
                    +aSpaceJt(0,0)*(aSpaceJt(1,1)*aSpaceJt(2,2)-aSpaceJt(2,1)*aSpaceJt(1,2))
                    -aSpaceJt(0,1)*(aSpaceJt(1,0)*aSpaceJt(2,2)-aSpaceJt(1,2)*aSpaceJt(2,0))
                    +aSpaceJt(0,2)*(aSpaceJt(1,0)*aSpaceJt(2,1)-aSpaceJt(1,1)*aSpaceJt(2,0));

            MORIS_ASSERT( tDetJ > 1e-18,
                    "Space determinant (bulk 3D) close to zero or negative: %e\n", tDetJ);

            MORIS_ASSERT( std::abs( det(aSpaceJt) - tDetJ ) < 1e-8*tDetJ,
                    "Inconsistent space determinant (bulk 3D): %e vs %e\n",
                    tDetJ,det(aSpaceJt));

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_space_detJ_bulk_tri_param_2(
                const Matrix< DDRMat > & aSpaceJt )
        {
            real tDetJ = (aSpaceJt(0,0)*aSpaceJt(1,1)-aSpaceJt(0,1)*aSpaceJt(1,0))/ 2.0;

            MORIS_ASSERT( tDetJ > 1e-18,
                    "Space determinant (Tri-P2) close to zero or negative: %e\n", tDetJ);

            MORIS_ASSERT( std::abs( det(aSpaceJt) - tDetJ ) < 1e-8*tDetJ,
                    "Inconsistent space determinant (Tri-P2): %e vs %e\n",
                    tDetJ,det(aSpaceJt));

            return tDetJ;
        }

        real Geometry_Interpolator::eval_space_detJ_bulk_tri_param_3(
                const Matrix< DDRMat > & aSpaceJt )
        {
            real tDetJ = (
                    +(aSpaceJt(1,0)*aSpaceJt(2,1)-aSpaceJt(1,1)*aSpaceJt(2,0))
                    -(aSpaceJt(0,0)*aSpaceJt(2,1)-aSpaceJt(0,1)*aSpaceJt(2,0))
                    +(aSpaceJt(0,0)*aSpaceJt(1,1)-aSpaceJt(0,1)*aSpaceJt(1,0))) / 2.0;

            MORIS_ASSERT( tDetJ > 1e-18,
                    "Space determinant (Tri-P3) close to zero or negative: %e\n", tDetJ);

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_space_detJ_bulk_tet_param_3(
                const Matrix< DDRMat > & aSpaceJt )
        {
            real tDetJ = (
                    +aSpaceJt(0,0)*(aSpaceJt(1,1)*aSpaceJt(2,2)-aSpaceJt(2,1)*aSpaceJt(1,2))
                    -aSpaceJt(0,1)*(aSpaceJt(1,0)*aSpaceJt(2,2)-aSpaceJt(1,2)*aSpaceJt(2,0))
                    +aSpaceJt(0,2)*(aSpaceJt(1,0)*aSpaceJt(2,1)-aSpaceJt(1,1)*aSpaceJt(2,0))) / 6.0;

            MORIS_ASSERT( tDetJ > 1e-18,
                    "Space determinant (Tet-P3) close to zero or negative: %e\n", tDetJ);

            MORIS_ASSERT( std::abs( det(aSpaceJt) - tDetJ ) < 1e-8*tDetJ,
                    "Inconsistent space determinant (Tet-P3): %e vs %e\n",
                    tDetJ,det(aSpaceJt));

            return tDetJ;
        }

        real Geometry_Interpolator::eval_space_detJ_bulk_tet_param_4(
                const Matrix< DDRMat > & aSpaceJt )
        {
            real tSubDet1 =
                    +aSpaceJt(1,0)*(aSpaceJt(2,1)*aSpaceJt(3,2)-aSpaceJt(2,2)*aSpaceJt(3,1))
                    -aSpaceJt(1,1)*(aSpaceJt(2,0)*aSpaceJt(3,2)-aSpaceJt(2,2)*aSpaceJt(3,0))
                    +aSpaceJt(1,2)*(aSpaceJt(2,0)*aSpaceJt(3,1)-aSpaceJt(2,1)*aSpaceJt(3,0));

            real tSubDet2 =
                    +aSpaceJt(0,0)*(aSpaceJt(2,1)*aSpaceJt(3,2)-aSpaceJt(2,2)*aSpaceJt(3,1))
                    -aSpaceJt(0,1)*(aSpaceJt(2,0)*aSpaceJt(3,2)-aSpaceJt(2,2)*aSpaceJt(3,0))
                    +aSpaceJt(0,2)*(aSpaceJt(2,0)*aSpaceJt(3,1)-aSpaceJt(2,1)*aSpaceJt(3,0));

            real tSubDet3 =
                    +aSpaceJt(0,0)*(aSpaceJt(1,1)*aSpaceJt(3,2)-aSpaceJt(1,2)*aSpaceJt(3,1))
                    -aSpaceJt(0,1)*(aSpaceJt(1,0)*aSpaceJt(3,2)-aSpaceJt(1,2)*aSpaceJt(3,0))
                    +aSpaceJt(0,2)*(aSpaceJt(1,0)*aSpaceJt(3,1)-aSpaceJt(1,1)*aSpaceJt(3,0));

            real tSubDet4 =
                    +aSpaceJt(0,0)*(aSpaceJt(1,1)*aSpaceJt(2,2)-aSpaceJt(1,2)*aSpaceJt(2,1))
                    -aSpaceJt(0,1)*(aSpaceJt(1,0)*aSpaceJt(2,2)-aSpaceJt(1,2)*aSpaceJt(2,0))
                    +aSpaceJt(0,2)*(aSpaceJt(1,0)*aSpaceJt(2,1)-aSpaceJt(1,1)*aSpaceJt(2,0));

            real tDetJ = ( tSubDet1 - tSubDet2 + tSubDet3 - tSubDet4 ) / 6.0;

            MORIS_ASSERT( tDetJ > 1e-18,
                    "Space determinant (Tet-P4) close to zero or negative: %e\n", tDetJ);

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_time_detJ_side(
                const Matrix< DDRMat > & aTimeJt )
        {
            return 1.0;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_time_detJ_bulk(
                const Matrix< DDRMat > & aTimeJt )
        {
            //FIXME: Needs to be specialized
            real tDetJ = det(aTimeJt);

            MORIS_ASSERT( tDetJ > 1e-18,
                    "Time determinant (bulk) close to zero or negative: %e\n", tDetJ);

            return tDetJ;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::get_normal( Matrix< DDRMat > & aNormal )
        {
            // check that there is a side interpolation
            MORIS_ASSERT( mSpaceSideset,
                    "Geometry_Interpolator::normal - not a side." );

            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Geometry_Interpolator::normal - mXiLocal is not set." );

            // evaluate side space interpolation shape functions first parametric derivatives at aParamPoint
            Matrix< DDRMat > tdNSpacedXi;
            mSpaceInterpolation->eval_dNdXi( mXiLocal, tdNSpacedXi );

            // evaluation of tangent vectors to the space side in the physical space
            Matrix< DDRMat > tRealTangents = trans( tdNSpacedXi * mXHat );

            // computing the normal from the real tangent vectors
            ( this->*mNormalFunc )( tRealTangents, aNormal );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_normal_side_line(
                Matrix< DDRMat > & aRealTangents,
                Matrix< DDRMat > & aNormal )
        {
            aNormal = {
                    {   aRealTangents( 1 ) },
                    { - aRealTangents( 0 ) }};

            aNormal = aNormal / norm( aNormal );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_normal_side_quad(
                Matrix< DDRMat > & aRealTangents,
                Matrix< DDRMat > & aNormal )
        {
            aNormal = cross( aRealTangents.get_column( 0 ), aRealTangents.get_column( 1 ) );
            aNormal = aNormal / norm( aNormal );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_normal_side_tri(
                Matrix< DDRMat > & aRealTangents,
                Matrix< DDRMat > & aNormal )
        {
            aNormal = cross(
                    aRealTangents.get_column( 0 ) - aRealTangents.get_column( 2 ),
                    aRealTangents.get_column( 1 ) - aRealTangents.get_column( 2 ) );

            aNormal = aNormal / norm( aNormal );
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::valx()
        {
            // check that mTHat is set
            MORIS_ASSERT( mXHat.numel() > 0,
                    "Geometry_Interpolator::valx - mXHat is not set." );

            //evaluate the field
            return this->NXi() * mXHat ;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::valt()
        {
            // check that mTHat is set
            MORIS_ASSERT( mTHat.numel()>0,
                    "Geometry_Interpolator::valt - mTHat is not set." );

            //evaluate the field
            return this->NTau() * mTHat;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::map_integration_point()
        {
            // if eval mapping
            if( mMapFlag )
            {
                // check that mXiHat and mTauHat are set
                MORIS_ASSERT( mXiHat.numel() > 0,
                        "Geometry_Interpolator::eval_mapping - mXiHat is not set." );
                MORIS_ASSERT( mTauHat.numel() > 0,
                        "Geometry_Interpolator::eval_mapping - mTauHat is not set." );

                uint tSize = mXiHat.n_cols();

                // set mapped space coordinates
                mMappedPoint( { 0, tSize - 1 } ) = trans( this->NXi() * mXiHat );

                // set mapped time coordinates
                mMappedPoint( tSize ) = dot( this->NTau(), mTauHat );
            }

            return mMappedPoint;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::update_local_coordinates(
                Matrix< DDRMat > & aPhysCoordinates,
                Matrix< DDRMat > & aParamCoordinates )
        {
            // set max iteration
            uint tMaxIter = 10;

            // init previous param coords
            Matrix< DDRMat > tPreviousParamCoordinates(
                    aParamCoordinates.n_rows(),
                    aParamCoordinates.n_cols(),
                    10.0 );

            // Newton loop
            for( uint iIter = 0; iIter <= tMaxIter; iIter++ )
            {
                // break update
                real tError = norm( tPreviousParamCoordinates - aParamCoordinates );
                if( tError < 1e-12 )
                {
                    break;
                }

                // compute NXi
                Matrix< DDRMat > tNXi;
                mSpaceInterpolation->eval_N( aParamCoordinates, tNXi );

                // compute residual
                Matrix< DDRMat > tR = trans( aPhysCoordinates - tNXi * mXHat );

                // compute dNdXI
                Matrix< DDRMat > tdNdXi;
                mSpaceInterpolation->eval_dNdXi( aParamCoordinates, tdNdXi );

                // compute Jacobian
                Matrix< DDRMat > tJ = trans( tdNdXi * mXHat );

                // solve
                aParamCoordinates += trans( inv( tJ ) * tR );

                // update previous coords
                tPreviousParamCoordinates = aParamCoordinates;

                if( iIter == tMaxIter )
                {
                    MORIS_LOG_INFO("Geometry_Interpolator::update_local_coordinates - No convergence.");
                }
            }
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::space_jacobian_and_matrices_for_second_derivatives(
                Matrix< DDRMat >       & aJt,
                Matrix< DDRMat >       & aKt,
                Matrix< DDRMat >       & aLt,
                const Matrix< DDRMat > & adNdXi,
                const Matrix< DDRMat > & ad2NdXi2 )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel() > 0,
                    "Geometry_Interpolator::space_jacobian_and_matrices_for_second_derivatives - mXHat is not set." );

            // evaluate transposed of geometry Jacobian
            aJt = this->space_jacobian();

            // call calculator for second derivatives
            this->mSecondDerivativeMatricesSpace(
                    aJt, // contains first geometric derivs
                    aKt,
                    aLt,
                    ad2NdXi2,
                    mXHat );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::space_jacobian_and_matrices_for_third_derivatives(
                Matrix< DDRMat >       & aJt,   // contains first geometric derivs
                Matrix< DDRMat >       & aJ2bt, // contains second geometric derivs = second help matrix for 2nd field derivs
                Matrix< DDRMat >       & aJ3at, // first help matrix for 3rd field derivs
                Matrix< DDRMat >       & aJ3bt, // second help matrix for 3rd field derivs
                Matrix< DDRMat >       & aJ3ct, // third help matrix for 3rd field derivs
                const Matrix< DDRMat > & adNdXi,
                const Matrix< DDRMat > & ad2NdXi2,
                const Matrix< DDRMat > & ad3NdXi3 )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel() > 0,
                    "Geometry_Interpolator::space_jacobian_and_matrices_for_third_derivatives - mXHat is not set." );

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

        void Geometry_Interpolator::time_jacobian_and_matrices_for_second_derivatives(
                Matrix< DDRMat >       & aJt,
                Matrix< DDRMat >       & aKt,
                Matrix< DDRMat >       & aLt,
                const Matrix< DDRMat > & adNdTau,
                const Matrix< DDRMat > & ad2NdTau2 )
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

        void Geometry_Interpolator::eval_matrices_for_second_derivative_1d(
                const Matrix< DDRMat > & aJt,
                Matrix< DDRMat >       & aKt,
                Matrix< DDRMat >       & aLt,
                const Matrix< DDRMat > & ad2NdXi2,
                const Matrix< DDRMat > & aXHat )
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 1, 1 );
            aLt( 0, 0 ) = std::pow( aJt( 0, 0 ), 2 );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_third_derivative_1d(
                const Matrix< DDRMat > & aJt,
                const Matrix< DDRMat > & aJ2bt,
                Matrix< DDRMat >       & aJ3at,
                Matrix< DDRMat >       & aJ3bt,
                Matrix< DDRMat >       & aJ3ct,
                const Matrix< DDRMat > & ad3NdXi3,
                const Matrix< DDRMat > & aXHat )
        {
            // first help matrix
            aJ3at.set_size( 1, 1, std::pow( aJt( 0, 0 ), 3 ) );

            // second help matrix
            aJ3bt.set_size( 1, 1, 3 * aJ2bt( 0, 0 ) * aJt( 0, 0 ) );

            // third help matrix
            aJ3ct = ad3NdXi3 * aXHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_2d(
                const Matrix< DDRMat > & aJt,
                Matrix< DDRMat >       & aKt,
                Matrix< DDRMat >       & aLt,
                const Matrix< DDRMat > & ad2NdXi2,
                const Matrix< DDRMat > & aXHat )
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 3, 3 );
            aLt( 0, 0 ) = std::pow( aJt( 0, 0 ), 2 );
            aLt( 1, 0 ) = std::pow( aJt( 1, 0 ), 2 );
            aLt( 2, 0 ) = aJt( 0 , 0 ) * aJt( 1 , 0 );

            aLt( 0, 1 ) = std::pow( aJt( 0, 1 ), 2 );
            aLt( 1, 1 ) = std::pow( aJt( 1, 1 ), 2 );
            aLt( 2, 1 ) = aJt( 0 , 1 ) * aJt( 1, 1 );

            aLt( 0, 2 ) = 2.0 * aJt( 0, 0 ) * aJt( 0, 1 );
            aLt( 1, 2 ) = 2.0 * aJt( 1, 0 ) * aJt( 1, 1 );
            aLt( 2, 2 ) = aJt( 0, 0 ) * aJt( 1, 1 ) + aJt( 0, 1 ) * aJt( 1, 0 );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_third_derivative_2d(
                const Matrix< DDRMat > & aJt,
                const Matrix< DDRMat > & aJ2bt,
                Matrix< DDRMat >       & aJ3at,
                Matrix< DDRMat >       & aJ3bt,
                Matrix< DDRMat >       & aJ3ct,
                const Matrix< DDRMat > & ad3NdXi3,
                const Matrix< DDRMat > & aXHat )
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
            for( uint j = 0; j < 2; ++j )
            {
                aJ3at( 0, j ) = std::pow( aJt( 0, j ), 3 );
                aJ3at( 1, j ) = std::pow( aJt( 1, j ), 3 );
            }

            // Block (2) ------------------------------------------------
            aJ3at( 0, 2 ) = 3 * std::pow( aJt( 0, 0 ), 2 ) * aJt( 0, 1 );
            aJ3at( 1, 2 ) = 3 * std::pow( aJt( 1, 0 ), 2 ) * aJt( 1, 1 );

            aJ3at( 0, 3 ) = 3 * std::pow( aJt( 0, 1 ), 2 ) * aJt( 0, 0 );
            aJ3at( 1, 3 ) = 3 * std::pow( aJt( 1, 1 ), 2 ) * aJt( 1, 0 );

            // Block (3) ------------------------------------------------
            for( uint j = 0; j < 2; ++j )
            {
                aJ3at( 2, j ) = std::pow( aJt( 0, j ), 2 ) * aJt( 1, j );
                aJ3at( 3, j ) = std::pow( aJt( 1, j ), 2 ) * aJt( 0, j );
            }

            // Block (4) ------------------------------------------------
            aJ3at( 2, 2 ) = std::pow( aJt( 0, 0 ), 2 ) * aJt( 1, 1 )  +  2 * aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 0, 1 );
            aJ3at( 3, 2 ) = std::pow( aJt( 1, 0 ), 2 ) * aJt( 0, 1 )  +  2 * aJt( 1, 0 ) * aJt( 0, 0 ) * aJt( 1, 1 );

            aJ3at( 2, 3 ) = std::pow( aJt( 0, 1 ), 2 ) * aJt( 1, 0 )  +  2 * aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 0, 0 );
            aJ3at( 3, 3 ) = std::pow( aJt( 1, 1 ), 2 ) * aJt( 0, 0 )  +  2 * aJt( 1, 1 ) * aJt( 0, 1 ) * aJt( 1, 0 );

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
            for( uint j=0; j<2; ++j )
            {
                aJ3bt( 0, j ) = 3 * aJ2bt( 0, j ) * aJt( 0, j );
                aJ3bt( 1, j ) = 3 * aJ2bt( 1, j ) * aJt( 1, j );
            }

            // Block (2) ------------------------------------------------
            aJ3bt( 0, 2 ) =   3 * aJ2bt( 0, 0 ) * aJt( 0, 1 ) + 3 * aJ2bt( 0, 1 ) * aJt( 0, 0 );
            aJ3bt( 1, 2 ) =   3 * aJ2bt( 1, 0 ) * aJt( 1, 1 ) + 3 * aJ2bt( 1, 1 ) * aJt( 1, 0 );

            // Block (3) ------------------------------------------------
            for( uint j=0; j<2; ++j )
            {
                aJ3bt( 2, j ) = 2 * aJ2bt( 2, j ) * aJt( 0, j )  +  aJ2bt( 0, j ) * aJt( 1, j );
                aJ3bt( 3, j ) = 2 * aJ2bt( 2, j ) * aJt( 1, j )  +  aJ2bt( 1, j ) * aJt( 0, j );
            }

            // Block (4) ------------------------------------------------
            aJ3bt( 2, 2 ) =
                    2 * aJ2bt( 2, 0 ) * aJt( 0, 1 )  +  2 * aJ2bt( 2, 1 ) * aJt( 0, 0 ) +
                    1 * aJ2bt( 0, 1 ) * aJt( 1, 0 )  +  1 * aJ2bt( 0, 0 ) * aJt( 1, 1 );
            aJ3bt( 3, 2 ) =
                    2 * aJ2bt( 2, 0 ) * aJt( 1, 1 )  +  2 * aJ2bt( 2, 1 ) * aJt( 1, 0 ) +
                    1 * aJ2bt( 1, 1 ) * aJt( 0, 0 )  +  1 * aJ2bt( 1, 0 ) * aJt( 0, 1 );

            // third help matrix
            aJ3ct = ad3NdXi3 * aXHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_3d(
                const Matrix< DDRMat > & aJt,
                Matrix< DDRMat >       & aKt,
                Matrix< DDRMat >       & aLt,
                const Matrix< DDRMat > & ad2NdXi2,
                const Matrix< DDRMat > & aXHat )
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 6, 6 );
            for( uint j = 0; j < 3; ++j )
            {
                aLt( 0, j ) = std::pow( aJt( 0, j ), 2 );
                aLt( 1, j ) = std::pow( aJt( 1, j ), 2 );
                aLt( 2, j ) = std::pow( aJt( 2, j ), 2 );
                aLt( 3, j ) = aJt( 1 , j ) * aJt( 2 , j );
                aLt( 4, j ) = aJt( 0 , j ) * aJt( 2 , j );
                aLt( 5, j ) = aJt( 0 , j ) * aJt( 1 , j );
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

        void Geometry_Interpolator::eval_matrices_for_third_derivative_3d(
                const Matrix< DDRMat > & aJt,
                const Matrix< DDRMat > & aJ2bt,
                Matrix< DDRMat >       & aJ3at,
                Matrix< DDRMat >       & aJ3bt,
                Matrix< DDRMat >       & aJ3ct,
                const Matrix< DDRMat > & ad3NdXi3,
                const Matrix< DDRMat > & aXHat )
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
            for( uint j=0; j<3; ++j )
            {
                aJ3at( 0, j ) = std::pow( aJt( 0, j ), 3 );
                aJ3at( 1, j ) = std::pow( aJt( 1, j ), 3 );
                aJ3at( 2, j ) = std::pow( aJt( 2, j ), 3 );
            }

            // Block (2) ------------------------------------------------
            aJ3at( 0, 3 ) = 3 * std::pow( aJt( 0, 0 ), 2 ) * aJt( 0, 1 );
            aJ3at( 1, 3 ) = 3 * std::pow( aJt( 1, 0 ), 2 ) * aJt( 1, 1 );
            aJ3at( 2, 3 ) = 3 * std::pow( aJt( 2, 0 ), 2 ) * aJt( 2, 1 );

            aJ3at( 0, 4 ) = 3 * std::pow( aJt( 0, 0 ), 2 ) * aJt( 0, 2 );
            aJ3at( 1, 4 ) = 3 * std::pow( aJt( 1, 0 ), 2 ) * aJt( 1, 2 );
            aJ3at( 2, 4 ) = 3 * std::pow( aJt( 2, 0 ), 2 ) * aJt( 2, 2 );

            aJ3at( 0, 5 ) = 3 * std::pow( aJt( 0, 1 ), 2 ) * aJt( 0, 0 );
            aJ3at( 1, 5 ) = 3 * std::pow( aJt( 1, 1 ), 2 ) * aJt( 1, 0 );
            aJ3at( 2, 5 ) = 3 * std::pow( aJt( 2, 1 ), 2 ) * aJt( 2, 0 );

            aJ3at( 0, 6 ) = 3 * std::pow( aJt( 0, 1 ), 2 ) * aJt( 0, 2 );
            aJ3at( 1, 6 ) = 3 * std::pow( aJt( 1, 1 ), 2 ) * aJt( 1, 2 );
            aJ3at( 2, 6 ) = 3 * std::pow( aJt( 2, 1 ), 2 ) * aJt( 2, 2 );

            aJ3at( 0, 7 ) = 3 * std::pow( aJt( 0, 2 ), 2 ) * aJt( 0, 0 );
            aJ3at( 1, 7 ) = 3 * std::pow( aJt( 1, 2 ), 2 ) * aJt( 1, 0 );
            aJ3at( 2, 7 ) = 3 * std::pow( aJt( 2, 2 ), 2 ) * aJt( 2, 0 );

            aJ3at( 0, 8 ) = 3 * std::pow( aJt( 0, 2 ), 2 ) * aJt( 0, 1 );
            aJ3at( 1, 8 ) = 3 * std::pow( aJt( 1, 2 ), 2 ) * aJt( 1, 1 );
            aJ3at( 2, 8 ) = 3 * std::pow( aJt( 2, 2 ), 2 ) * aJt( 2, 1 );


            // Block (3) ------------------------------------------------
            aJ3at( 0, 9 ) = 6 * aJt( 0, 0 ) * aJt( 0, 1 ) * aJt( 0, 2 );
            aJ3at( 1, 9 ) = 6 * aJt( 1, 0 ) * aJt( 1, 1 ) * aJt( 1, 2 );
            aJ3at( 2, 9 ) = 6 * aJt( 2, 0 ) * aJt( 2, 1 ) * aJt( 2, 2 );

            // Block (4) ------------------------------------------------
            for( uint j=0; j<3; ++j )
            {
                aJ3at( 3, j ) = std::pow( aJt( 0, j ) , 2 ) * aJt( 1, j );
                aJ3at( 4, j ) = std::pow( aJt( 0, j ) , 2 ) * aJt( 2, j );
                aJ3at( 5, j ) = std::pow( aJt( 1, j ) , 2 ) * aJt( 0, j );
                aJ3at( 6, j ) = std::pow( aJt( 1, j ) , 2 ) * aJt( 2, j );
                aJ3at( 7, j ) = std::pow( aJt( 2, j ) , 2 ) * aJt( 0, j );
                aJ3at( 8, j ) = std::pow( aJt( 2, j ) , 2 ) * aJt( 1, j );
            }

            // Block (5) ------------------------------------------------
            aJ3at( 3, 3 ) = std::pow( aJt( 0, 0 ), 2 ) * aJt( 1, 1 ) + 2 * aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 0, 1 );
            aJ3at( 4, 3 ) = std::pow( aJt( 0, 0 ), 2 ) * aJt( 2, 1 ) + 2 * aJt( 0, 0 ) * aJt( 2, 0 ) * aJt( 0, 1 );
            aJ3at( 5, 3 ) = std::pow( aJt( 1, 0 ), 2 ) * aJt( 0, 1 ) + 2 * aJt( 1, 0 ) * aJt( 0, 0 ) * aJt( 1, 1 );
            aJ3at( 6, 3 ) = std::pow( aJt( 1, 0 ), 2 ) * aJt( 2, 1 ) + 2 * aJt( 1, 0 ) * aJt( 2, 0 ) * aJt( 1, 1 );
            aJ3at( 7, 3 ) = std::pow( aJt( 2, 0 ), 2 ) * aJt( 0, 1 ) + 2 * aJt( 2, 0 ) * aJt( 0, 0 ) * aJt( 2, 1 );
            aJ3at( 8, 3 ) = std::pow( aJt( 2, 0 ), 2 ) * aJt( 1, 1 ) + 2 * aJt( 2, 0 ) * aJt( 1, 0 ) * aJt( 2, 1 );

            aJ3at( 3, 4 ) = std::pow( aJt( 0, 0 ), 2 ) * aJt( 1, 2 ) + 2 * aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 0, 2 );
            aJ3at( 4, 4 ) = std::pow( aJt( 0, 0 ), 2 ) * aJt( 2, 2 ) + 2 * aJt( 0, 0 ) * aJt( 2, 0 ) * aJt( 0, 2 );
            aJ3at( 5, 4 ) = std::pow( aJt( 1, 0 ), 2 ) * aJt( 0, 2 ) + 2 * aJt( 1, 0 ) * aJt( 0, 0 ) * aJt( 1, 2 );
            aJ3at( 6, 4 ) = std::pow( aJt( 1, 0 ), 2 ) * aJt( 2, 2 ) + 2 * aJt( 1, 0 ) * aJt( 2, 0 ) * aJt( 1, 2 );
            aJ3at( 7, 4 ) = std::pow( aJt( 2, 0 ), 2 ) * aJt( 0, 2 ) + 2 * aJt( 2, 0 ) * aJt( 0, 0 ) * aJt( 2, 2 );
            aJ3at( 8, 4 ) = std::pow( aJt( 2, 0 ), 2 ) * aJt( 1, 2 ) + 2 * aJt( 2, 0 ) * aJt( 1, 0 ) * aJt( 2, 2 );

            aJ3at( 3, 5 ) = std::pow( aJt( 0, 1 ), 2 ) * aJt( 1, 0 ) + 2 * aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 0, 0 );
            aJ3at( 4, 5 ) = std::pow( aJt( 0, 1 ), 2 ) * aJt( 2, 0 ) + 2 * aJt( 0, 1 ) * aJt( 2, 1 ) * aJt( 0, 0 );
            aJ3at( 5, 5 ) = std::pow( aJt( 1, 1 ), 2 ) * aJt( 0, 0 ) + 2 * aJt( 1, 1 ) * aJt( 0, 1 ) * aJt( 1, 0 );
            aJ3at( 6, 5 ) = std::pow( aJt( 1, 1 ), 2 ) * aJt( 2, 0 ) + 2 * aJt( 1, 1 ) * aJt( 2, 1 ) * aJt( 1, 0 );
            aJ3at( 7, 5 ) = std::pow( aJt( 2, 1 ), 2 ) * aJt( 0, 0 ) + 2 * aJt( 2, 1 ) * aJt( 0, 1 ) * aJt( 2, 0 );
            aJ3at( 8, 5 ) = std::pow( aJt( 2, 1 ), 2 ) * aJt( 1, 0 ) + 2 * aJt( 2, 1 ) * aJt( 1, 1 ) * aJt( 2, 0 );

            aJ3at( 3, 6 ) = std::pow( aJt( 0, 1 ), 2 ) * aJt( 1, 2 ) + 2 * aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 0, 2 );
            aJ3at( 4, 6 ) = std::pow( aJt( 0, 1 ), 2 ) * aJt( 2, 2 ) + 2 * aJt( 0, 1 ) * aJt( 2, 1 ) * aJt( 0, 2 );
            aJ3at( 5, 6 ) = std::pow( aJt( 1, 1 ), 2 ) * aJt( 0, 2 ) + 2 * aJt( 1, 1 ) * aJt( 0, 1 ) * aJt( 1, 2 );
            aJ3at( 6, 6 ) = std::pow( aJt( 1, 1 ), 2 ) * aJt( 2, 2 ) + 2 * aJt( 1, 1 ) * aJt( 2, 1 ) * aJt( 1, 2 );
            aJ3at( 7, 6 ) = std::pow( aJt( 2, 1 ), 2 ) * aJt( 0, 2 ) + 2 * aJt( 2, 1 ) * aJt( 0, 1 ) * aJt( 2, 2 );
            aJ3at( 8, 6 ) = std::pow( aJt( 2, 1 ), 2 ) * aJt( 1, 2 ) + 2 * aJt( 2, 1 ) * aJt( 1, 1 ) * aJt( 2, 2 );

            aJ3at( 3, 7 ) = std::pow( aJt( 0, 2 ), 2 ) * aJt( 1, 0 ) + 2 * aJt( 0, 2 ) * aJt( 1, 2 ) * aJt( 0, 0 );
            aJ3at( 4, 7 ) = std::pow( aJt( 0, 2 ), 2 ) * aJt( 2, 0 ) + 2 * aJt( 0, 2 ) * aJt( 2, 2 ) * aJt( 0, 0 );
            aJ3at( 5, 7 ) = std::pow( aJt( 1, 2 ), 2 ) * aJt( 0, 0 ) + 2 * aJt( 1, 2 ) * aJt( 0, 2 ) * aJt( 1, 0 );
            aJ3at( 6, 7 ) = std::pow( aJt( 1, 2 ), 2 ) * aJt( 2, 0 ) + 2 * aJt( 1, 2 ) * aJt( 2, 2 ) * aJt( 1, 0 );
            aJ3at( 7, 7 ) = std::pow( aJt( 2, 2 ), 2 ) * aJt( 0, 0 ) + 2 * aJt( 2, 2 ) * aJt( 0, 2 ) * aJt( 2, 0 );
            aJ3at( 8, 7 ) = std::pow( aJt( 2, 2 ), 2 ) * aJt( 1, 0 ) + 2 * aJt( 2, 2 ) * aJt( 1, 2 ) * aJt( 2, 0 );

            aJ3at( 3, 8 ) = std::pow( aJt( 0, 2 ), 2 ) * aJt( 1, 1 ) + 2 * aJt( 0, 2 ) * aJt( 1, 2 ) * aJt( 0, 1 );
            aJ3at( 4, 8 ) = std::pow( aJt( 0, 2 ), 2 ) * aJt( 2, 1 ) + 2 * aJt( 0, 2 ) * aJt( 2, 2 ) * aJt( 0, 1 );
            aJ3at( 5, 8 ) = std::pow( aJt( 1, 2 ), 2 ) * aJt( 0, 1 ) + 2 * aJt( 1, 2 ) * aJt( 0, 2 ) * aJt( 1, 1 );
            aJ3at( 6, 8 ) = std::pow( aJt( 1, 2 ), 2 ) * aJt( 2, 1 ) + 2 * aJt( 1, 2 ) * aJt( 2, 2 ) * aJt( 1, 1 );
            aJ3at( 7, 8 ) = std::pow( aJt( 2, 2 ), 2 ) * aJt( 0, 1 ) + 2 * aJt( 2, 2 ) * aJt( 0, 2 ) * aJt( 2, 1 );
            aJ3at( 8, 8 ) = std::pow( aJt( 2, 2 ), 2 ) * aJt( 1, 1 ) + 2 * aJt( 2, 2 ) * aJt( 1, 2 ) * aJt( 2, 1 );


            // Block (6) ------------------------------------------------
            aJ3at( 3, 9 ) =
                    2 * aJt( 0, 0 ) * aJt( 0, 1 ) * aJt( 1, 2 ) +
                    2 * aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 0, 2 ) +
                    2 * aJt( 1, 0 ) * aJt( 0, 1 ) * aJt( 0, 2 );
            aJ3at( 4, 9 ) =
                    2 * aJt( 0, 0 ) * aJt( 0, 1 ) * aJt( 2, 2 ) +
                    2 * aJt( 0, 0 ) * aJt( 2, 1 ) * aJt( 0, 2 ) +
                    2 * aJt( 2, 0 ) * aJt( 0, 1 ) * aJt( 0, 2 );
            aJ3at( 5, 9 ) =
                    2 * aJt( 1, 0 ) * aJt( 1, 1 ) * aJt( 0, 2 ) +
                    2 * aJt( 1, 0 ) * aJt( 0, 1 ) * aJt( 1, 2 ) +
                    2 * aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 1, 2 );
            aJ3at( 6, 9 ) =
                    2 * aJt( 1, 0 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +
                    2 * aJt( 1, 0 ) * aJt( 2, 1 ) * aJt( 1, 2 ) +
                    2 * aJt( 2, 0 ) * aJt( 1, 1 ) * aJt( 1, 2 );
            aJ3at( 7, 9 ) =
                    2 * aJt( 2, 0 ) * aJt( 2, 1 ) * aJt( 0, 2 ) +
                    2 * aJt( 2, 0 ) * aJt( 0, 1 ) * aJt( 2, 2 ) +
                    2 * aJt( 0, 0 ) * aJt( 2, 1 ) * aJt( 2, 2 );
            aJ3at( 8, 9 ) =
                    2 * aJt( 2, 0 ) * aJt( 2, 1 ) * aJt( 1, 2 ) +
                    2 * aJt( 2, 0 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +
                    2 * aJt( 1, 0 ) * aJt( 2, 1 ) * aJt( 2, 2 );

            // Block (7) ------------------------------------------------
            for( uint j=0; j<3; ++j )
            {
                aJ3at( 9, j ) = aJt( 0, j ) * aJt( 1, j ) * aJt( 2, j );
            }

            // Block (8) ------------------------------------------------
            aJ3at( 9, 3 ) =
                    aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 2, 1 ) +
                    aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 2, 0 ) +
                    aJt( 0, 1 ) * aJt( 1, 0 ) * aJt( 2, 0 ) ;

            aJ3at( 9, 4 ) =
                    aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 2, 2 ) +
                    aJt( 0, 0 ) * aJt( 1, 2 ) * aJt( 2, 0 ) +
                    aJt( 0, 2 ) * aJt( 1, 0 ) * aJt( 2, 0 );

            aJ3at( 9, 5 ) =
                    aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 2, 0 ) +
                    aJt( 0, 1 ) * aJt( 1, 0 ) * aJt( 2, 1 ) +
                    aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 2, 1 );

            aJ3at( 9, 6 ) =
                    aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +
                    aJt( 0, 1 ) * aJt( 1, 2 ) * aJt( 2, 1 ) +
                    aJt( 0, 2 ) * aJt( 1, 1 ) * aJt( 2, 1 );

            aJ3at( 9, 7 ) =
                    aJt( 0, 2 ) * aJt( 1, 2 ) * aJt( 2, 0 ) +
                    aJt( 0, 2 ) * aJt( 1, 0 ) * aJt( 2, 2 ) +
                    aJt( 0, 0 ) * aJt( 1, 2 ) * aJt( 2, 2 );

            aJ3at( 9, 8 ) =
                    aJt( 0, 2 ) * aJt( 1, 2 ) * aJt( 2, 1 ) +
                    aJt( 0, 2 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +
                    aJt( 0, 1 ) * aJt( 1, 2 ) * aJt( 2, 2 );

            // Block (9) ------------------------------------------------
            aJ3at( 9, 9 ) =
                    aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +
                    aJt( 0, 2 ) * aJt( 1, 1 ) * aJt( 2, 0 ) +
                    aJt( 0, 1 ) * aJt( 1, 2 ) * aJt( 2, 0 ) +
                    aJt( 0, 0 ) * aJt( 1, 2 ) * aJt( 2, 1 ) +
                    aJt( 0, 2 ) * aJt( 1, 0 ) * aJt( 2, 1 ) +
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
            for( uint j=0; j<3; ++j )
            {
                aJ3bt( 0, j ) = 3 * aJ2bt( 0, j ) * aJt( 0, j );
                aJ3bt( 1, j ) = 3 * aJ2bt( 1, j ) * aJt( 1, j );
                aJ3bt( 2, j ) = 3 * aJ2bt( 2, j ) * aJt( 2, j );
            }

            // Block (2) ------------------------------------------------
            aJ3bt( 0, 5 ) =
                    3 * aJ2bt( 0, 0 ) * aJt( 0, 1 ) +
                    3 * aJ2bt( 0, 1 ) * aJt( 0, 0 );
            aJ3bt( 1, 5 ) =
                    3 * aJ2bt( 1, 0 ) * aJt( 1, 1 ) +
                    3 * aJ2bt( 1, 1 ) * aJt( 1, 0 );
            aJ3bt( 2, 5 ) =
                    3 * aJ2bt( 2, 0 ) * aJt( 2, 1 ) +
                    3 * aJ2bt( 2, 1 ) * aJt( 2, 0 );

            aJ3bt( 0, 3 ) =
                    3 * aJ2bt( 0, 1 ) * aJt( 0, 2 ) +
                    3 * aJ2bt( 0, 2 ) * aJt( 0, 1 );
            aJ3bt( 1, 3 ) =
                    3 * aJ2bt( 1, 1 ) * aJt( 1, 2 ) +
                    3 * aJ2bt( 1, 2 ) * aJt( 1, 1 );
            aJ3bt( 2, 3 ) =
                    3 * aJ2bt( 2, 1 ) * aJt( 2, 2 ) +
                    3 * aJ2bt( 2, 2 ) * aJt( 2, 1 );

            aJ3bt( 0, 4 ) =
                    3 * aJ2bt( 0, 0 ) * aJt( 0, 2 ) +
                    3 * aJ2bt( 0, 2 ) * aJt( 0, 0 );
            aJ3bt( 1, 4 ) =
                    3 * aJ2bt( 1, 0 ) * aJt( 1, 2 ) +
                    3 * aJ2bt( 1, 2 ) * aJt( 1, 0 );
            aJ3bt( 2, 4 ) =
                    3 * aJ2bt( 2, 0 ) * aJt( 2, 2 ) +
                    3 * aJ2bt( 2, 2 ) * aJt( 2, 0 );

            // Block (3) ------------------------------------------------
            for( uint j=0; j<3; ++j )
            {
                aJ3bt( 3, j ) = 2 * aJ2bt( 5, j ) * aJt( 0, j )  +  aJ2bt( 0, j ) * aJt( 1, j );
                aJ3bt( 4, j ) = 2 * aJ2bt( 4, j ) * aJt( 0, j )  +  aJ2bt( 0, j ) * aJt( 2, j );
                aJ3bt( 5, j ) = 2 * aJ2bt( 5, j ) * aJt( 1, j )  +  aJ2bt( 1, j ) * aJt( 0, j );
                aJ3bt( 6, j ) = 2 * aJ2bt( 3, j ) * aJt( 1, j )  +  aJ2bt( 1, j ) * aJt( 2, j );
                aJ3bt( 7, j ) = 2 * aJ2bt( 4, j ) * aJt( 2, j )  +  aJ2bt( 2, j ) * aJt( 0, j );
                aJ3bt( 8, j ) = 2 * aJ2bt( 3, j ) * aJt( 2, j )  +  aJ2bt( 2, j ) * aJt( 1, j );
            }

            // Block (4) ------------------------------------------------
            aJ3bt( 3, 5 ) =  2 * aJ2bt( 5, 0 ) * aJt( 0, 1 )  +  2 * aJ2bt( 5, 1 ) * aJt( 0, 0 )  +  aJ2bt( 0, 1 ) * aJt( 1, 0 )  +  aJ2bt( 0, 0 ) * aJt( 1, 1 );
            aJ3bt( 4, 5 ) =  2 * aJ2bt( 4, 0 ) * aJt( 0, 1 )  +  2 * aJ2bt( 4, 1 ) * aJt( 0, 0 )  +  aJ2bt( 0, 1 ) * aJt( 2, 0 )  +  aJ2bt( 0, 0 ) * aJt( 2, 1 );
            aJ3bt( 5, 5 ) =  2 * aJ2bt( 5, 0 ) * aJt( 1, 1 )  +  2 * aJ2bt( 5, 1 ) * aJt( 1, 0 )  +  aJ2bt( 1, 1 ) * aJt( 0, 0 )  +  aJ2bt( 1, 0 ) * aJt( 0, 1 );
            aJ3bt( 6, 5 ) =  2 * aJ2bt( 3, 0 ) * aJt( 1, 1 )  +  2 * aJ2bt( 3, 1 ) * aJt( 1, 0 )  +  aJ2bt( 1, 1 ) * aJt( 2, 0 )  +  aJ2bt( 1, 0 ) * aJt( 2, 1 );
            aJ3bt( 7, 5 ) =  2 * aJ2bt( 4, 0 ) * aJt( 2, 1 )  +  2 * aJ2bt( 4, 1 ) * aJt( 2, 0 )  +  aJ2bt( 2, 1 ) * aJt( 0, 0 )  +  aJ2bt( 2, 0 ) * aJt( 0, 1 );
            aJ3bt( 8, 5 ) =  2 * aJ2bt( 3, 0 ) * aJt( 2, 1 )  +  2 * aJ2bt( 3, 1 ) * aJt( 2, 0 )  +  aJ2bt( 2, 1 ) * aJt( 1, 0 )  +  aJ2bt( 2, 0 ) * aJt( 1, 1 );

            aJ3bt( 3, 3 ) =  2 * aJ2bt( 5, 1 ) * aJt( 0, 2 )  +  2 * aJ2bt( 5, 2 ) * aJt( 0, 1 )  +  aJ2bt( 0, 2 ) * aJt( 1, 1 )  +  aJ2bt( 0, 1 ) * aJt( 1, 2 );
            aJ3bt( 4, 3 ) =  2 * aJ2bt( 4, 1 ) * aJt( 0, 2 )  +  2 * aJ2bt( 4, 2 ) * aJt( 0, 1 )  +  aJ2bt( 0, 2 ) * aJt( 2, 1 )  +  aJ2bt( 0, 1 ) * aJt( 2, 2 );
            aJ3bt( 5, 3 ) =  2 * aJ2bt( 5, 1 ) * aJt( 1, 2 )  +  2 * aJ2bt( 5, 2 ) * aJt( 1, 1 )  +  aJ2bt( 1, 2 ) * aJt( 0, 1 )  +  aJ2bt( 1, 1 ) * aJt( 0, 2 );
            aJ3bt( 6, 3 ) =  2 * aJ2bt( 3, 1 ) * aJt( 1, 2 )  +  2 * aJ2bt( 3, 2 ) * aJt( 1, 1 )  +  aJ2bt( 1, 2 ) * aJt( 2, 1 )  +  aJ2bt( 1, 1 ) * aJt( 2, 2 );
            aJ3bt( 7, 3 ) =  2 * aJ2bt( 4, 1 ) * aJt( 2, 2 )  +  2 * aJ2bt( 4, 2 ) * aJt( 2, 1 )  +  aJ2bt( 2, 2 ) * aJt( 0, 1 )  +  aJ2bt( 2, 1 ) * aJt( 0, 2 );
            aJ3bt( 8, 3 ) =  2 * aJ2bt( 3, 1 ) * aJt( 2, 2 )  +  2 * aJ2bt( 3, 2 ) * aJt( 2, 1 )  +  aJ2bt( 2, 2 ) * aJt( 1, 1 )  +  aJ2bt( 2, 1 ) * aJt( 1, 2 );

            aJ3bt( 3, 4 ) =  2 * aJ2bt( 5, 0 ) * aJt( 0, 2 )  +  2 * aJ2bt( 5, 2 ) * aJt( 0, 0 )  +  aJ2bt( 0, 2 ) * aJt( 1, 0 )  +  aJ2bt( 0, 0 ) * aJt( 1, 2 );
            aJ3bt( 4, 4 ) =  2 * aJ2bt( 4, 0 ) * aJt( 0, 2 )  +  2 * aJ2bt( 4, 2 ) * aJt( 0, 0 )  +  aJ2bt( 0, 2 ) * aJt( 2, 0 )  +  aJ2bt( 0, 0 ) * aJt( 2, 2 );
            aJ3bt( 5, 4 ) =  2 * aJ2bt( 5, 0 ) * aJt( 1, 2 )  +  2 * aJ2bt( 5, 2 ) * aJt( 1, 0 )  +  aJ2bt( 1, 2 ) * aJt( 0, 0 )  +  aJ2bt( 1, 0 ) * aJt( 0, 2 );
            aJ3bt( 6, 4 ) =  2 * aJ2bt( 3, 0 ) * aJt( 1, 2 )  +  2 * aJ2bt( 3, 2 ) * aJt( 1, 0 )  +  aJ2bt( 1, 2 ) * aJt( 2, 0 )  +  aJ2bt( 1, 0 ) * aJt( 2, 2 );
            aJ3bt( 7, 4 ) =  2 * aJ2bt( 4, 0 ) * aJt( 2, 2 )  +  2 * aJ2bt( 4, 2 ) * aJt( 2, 0 )  +  aJ2bt( 2, 2 ) * aJt( 0, 0 )  +  aJ2bt( 2, 0 ) * aJt( 0, 2 );
            aJ3bt( 8, 4 ) =  2 * aJ2bt( 3, 0 ) * aJt( 2, 2 )  +  2 * aJ2bt( 3, 2 ) * aJt( 2, 0 )  +  aJ2bt( 2, 2 ) * aJt( 1, 0 )  +  aJ2bt( 2, 0 ) * aJt( 1, 2 );

            // Block (5) ------------------------------------------------
            for( uint j=0; j<3; ++j )
            {
                aJ3bt( 9, j ) =
                        aJ2bt( 4, j ) * aJt( 1, j ) +
                        aJ2bt( 3, j ) * aJt( 0, j ) +
                        aJ2bt( 5, j ) * aJt( 2, j );
            }

            // Block (6) ------------------------------------------------
            aJ3bt( 9, 5 ) =
                    aJ2bt( 5, 0 ) * aJt( 2, 1 ) + aJ2bt( 5, 1 ) * aJt( 2, 0 ) +
                    aJ2bt( 3, 0 ) * aJt( 0, 1 ) + aJ2bt( 3, 1 ) * aJt( 0, 0 ) +
                    aJ2bt( 4, 0 ) * aJt( 1, 1 ) + aJ2bt( 4, 1 ) * aJt( 1, 0 );

            aJ3bt( 9, 3 ) =
                    aJ2bt( 5, 1 ) * aJt( 2, 2 ) + aJ2bt( 5, 2 ) * aJt( 2, 1 ) +
                    aJ2bt( 3, 1 ) * aJt( 0, 2 ) + aJ2bt( 3, 2 ) * aJt( 0, 1 ) +
                    aJ2bt( 4, 1 ) * aJt( 1, 2 ) + aJ2bt( 4, 2 ) * aJt( 1, 1 );

            aJ3bt( 9, 4 ) =
                    aJ2bt( 5, 0 ) * aJt( 2, 2 ) + aJ2bt( 5, 2 ) * aJt( 2, 0 ) +
                    aJ2bt( 3, 0 ) * aJt( 0, 2 ) + aJ2bt( 3, 2 ) * aJt( 0, 0 ) +
                    aJ2bt( 4, 0 ) * aJt( 1, 2 ) + aJ2bt( 4, 2 ) * aJt( 1, 0 );

            // third help matrix
            aJ3ct = ad3NdXi3 * aXHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_function_pointers()
        {
            // get number of dimensions and set pointer to function
            // for second space derivative
            switch ( mNumSpaceDim )
            {
                case 1 :
                {
                    mInvSpaceJacFunc               = &Geometry_Interpolator::eval_inverse_space_jacobian_1d;
                    mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_1d;
                    mThirdDerivativeMatricesSpace  = this->eval_matrices_for_third_derivative_1d;
                    break;
                }
                case 2 :
                {
                    mInvSpaceJacFunc               = &Geometry_Interpolator::eval_inverse_space_jacobian_2d;
                    mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_2d;
                    mThirdDerivativeMatricesSpace  = this->eval_matrices_for_third_derivative_2d;
                    break;
                }
                case 3 :
                {
                    mInvSpaceJacFunc               = &Geometry_Interpolator::eval_inverse_space_jacobian_3d;
                    mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_3d;
                    mThirdDerivativeMatricesSpace  = this->eval_matrices_for_third_derivative_3d;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, " Geometry_Interpolator::set_function_pointers - unknown number of dimensions. " );
                }
            }

            // get number of dimensions and set pointer to function
            // for second time derivative
            switch ( mNumTimeDim )
            {
                case 1 :
                {
                    mInvTimeJacFunc               = &Geometry_Interpolator::eval_inverse_time_jacobian_1d;
                    mSecondDerivativeMatricesTime = this->eval_matrices_for_second_derivative_1d;
                    break;
                }
                case 2 :
                {
                    mInvTimeJacFunc               = &Geometry_Interpolator::eval_inverse_time_jacobian_2d;
                    mSecondDerivativeMatricesTime = this->eval_matrices_for_second_derivative_2d;
                    break;
                }
                case 3 :
                {
                    mInvTimeJacFunc               = &Geometry_Interpolator::eval_inverse_time_jacobian_3d;
                    mSecondDerivativeMatricesTime = this->eval_matrices_for_second_derivative_3d;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, " Geometry_Interpolator::set_function_pointers - unknown number of dimensions. " );
                }
            }

            // switch on geometry type
            if( mSpaceSideset )
            {
                switch( mGeometryType )
                {
                    case mtk::Geometry_Type::LINE :
                    {
                        mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_side_line;
                        mNormalFunc    = &Geometry_Interpolator::eval_normal_side_line;

                        // set size for storage
                        mMapFlag = true;
                        mMappedPoint.set_size( 3, 1, 0.0 );
                        break;
                    }
                    case mtk::Geometry_Type::TRI :
                    {
                        mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_side_tri;
                        mNormalFunc    = &Geometry_Interpolator::eval_normal_side_tri;

                        // set size for storage
                        mMapFlag = true;
                        mMappedPoint.set_size( 4, 1, 0.0 );
                        break;
                    }
                    case mtk::Geometry_Type::QUAD :
                    {
                        mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_side_quad;
                        mNormalFunc    = &Geometry_Interpolator::eval_normal_side_quad;

                        // set size for storage
                        mMapFlag = true;
                        mMappedPoint.set_size( 4, 1, 0.0 );
                        break;
                    }
                    default :
                    {
                        MORIS_ERROR( false, "Geometry_Interpolator::set_function_pointers - unknown or not implemented side space geometry type ");
                    }
                }
            }
            else
            {
                switch( mGeometryType )
                {
                    case mtk::Geometry_Type::LINE :
                    {
                        mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_bulk_line;

                        // set size for storage
                        mMappedPoint.set_size( 2, 1, 0.0 );
                        break;
                    }

                    case mtk::Geometry_Type::QUAD :
                    {
                        mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_bulk_quad;

                        // set size for storage
                        mMappedPoint.set_size( 3, 1, 0.0 );
                        break;
                    }

                    case mtk::Geometry_Type::HEX :
                    {
                        mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_bulk_hex;

                        // set size for storage
                        mMappedPoint.set_size( 4, 1, 0.0 );
                        break;
                    }

                    case mtk::Geometry_Type::TRI :
                    {
                        switch( mNumSpaceParamDim )
                        {
                            case 2 :
                                mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_bulk_tri_param_2;
                                break;
                            case 3 :
                                mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_bulk_tri_param_3;
                                break;
                            default :
                                MORIS_ERROR( false, " Geometry_Interpolator::set_function_pointers - Parametric space dimensions can only be 2 or 3." );
                        }
                        // set size for storage
                        mMapFlag = true;
                        mMappedPoint.set_size( 3, 1, 0.0 );
                        break;
                    }

                    case mtk::Geometry_Type::TET :
                    {
                        switch( mNumSpaceParamDim )
                        {
                            case 3 :
                                mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_bulk_tet_param_3;
                                break;
                            case 4 :
                                mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_bulk_tet_param_4;
                                break;
                            default :
                                MORIS_ERROR( false, " Geometry_Interpolator::set_function_pointers - Parametric space dimensions can only be 3 or 4." );
                        }
                        // set size for storage
                        mMapFlag = true;
                        mMappedPoint.set_size( 4, 1, 0.0 );
                        break;
                    }

                    default :
                    {
                        MORIS_ERROR( false, "Geometry_Interpolator::set_function_pointers - unknown or not implemented space geometry type ");
                    }
                }
            }

            // switch Geometry_Type
            switch ( mTimeGeometryType )
            {
                case mtk::Geometry_Type::POINT :
                {
                    mMapFlag = true;
                    mTimeDetJFunc = &Geometry_Interpolator::eval_time_detJ_side;
                    break;
                }
                case mtk::Geometry_Type::LINE :
                {
                    mTimeDetJFunc = &Geometry_Interpolator::eval_time_detJ_bulk;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Geometry_Interpolator::set_function_pointers - unknown or not implemented time geometry type ");
                }
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

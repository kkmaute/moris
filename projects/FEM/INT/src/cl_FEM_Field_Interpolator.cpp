
#include "cl_Matrix.hpp"                //LNA/src
#include "linalg_typedefs.hpp"
#include "fn_linsolve.hpp"
#include "fn_det.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_reshape.hpp"
#include "fn_det.hpp"
#include "op_times.hpp"
#include "op_equal_equal.hpp"
#include "op_less_equal.hpp"
#include "op_greater_equal.hpp"

#include "cl_FEM_Field_Interpolator.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        Field_Interpolator::Field_Interpolator( const uint                   & aNumberOfFields,
                                                const Interpolation_Rule     & aFieldInterpolationRule,
                                                      Geometry_Interpolator*   aGeometryInterpolator,
                                                const moris::Cell< MSI::Dof_Type >            aDofType )
                                              : mNumberOfFields( aNumberOfFields ),
                                                mGeometryInterpolator( aGeometryInterpolator ),
                                                mDofType( aDofType )
        {
            // create space and time interpolation function
            mSpaceInterpolation = aFieldInterpolationRule.create_space_interpolation_function();
            mTimeInterpolation  = aFieldInterpolationRule.create_time_interpolation_function();

            // get number of space, time dimensions
            mNSpaceDim = mSpaceInterpolation->get_number_of_dimensions();
            mNTimeDim  = mTimeInterpolation ->get_number_of_dimensions();

            // get number of space parametric dimensions
            mNSpaceParamDim = mSpaceInterpolation->get_number_of_param_dimensions();

            // check dimensions consistency
            MORIS_ERROR( ( mNSpaceDim == mGeometryInterpolator->get_number_of_space_dimensions() ) ,
                         "Field_Interpolator - Space dimension inconsistency." );
            MORIS_ERROR( ( mNTimeDim  == mGeometryInterpolator->get_number_of_time_dimensions() ),
                         "Field_Interpolator - Time dimension inconsistency.");

            // get number of space, time, and space time basis
            mNSpaceBases = mSpaceInterpolation->get_number_of_bases();
            mNTimeBases  = mTimeInterpolation->get_number_of_bases();
            mNFieldBases = mNSpaceBases * mNTimeBases;

            // get number of coefficients
            mNFieldCoeff = mNFieldBases * mNumberOfFields;
        }

//------------------------------------------------------------------------------

        Field_Interpolator::~Field_Interpolator()
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

        void Field_Interpolator::set_space_time( const Matrix< DDRMat > & aParamPoint )
        {
            // check input size aParamPoint
            MORIS_ASSERT( ( ( aParamPoint.n_cols() == 1 ) && ( aParamPoint.n_rows() == mNSpaceParamDim + mNTimeDim )),
                         "Field_Interpolator::set_space_time - Wrong input size ( aParamPoint ).");

            // check input values are between -1 and 1
            // fixme what about TRI and TET
            for ( uint Ik = 0; Ik < mNSpaceParamDim + mNTimeDim; Ik++ )
            {

                MORIS_ASSERT( ( ( aParamPoint( Ik ) <= 1.0 + 1E-12 ) && ( aParamPoint( Ik ) >= -1.0 - 1E-12 ) ),
                             "Field_Interpolator::set_space_time - Wrong input value ( aParamPoint ).");
            }

            // set input values
            mXi.matrix_data()  = aParamPoint( { 0, mNSpaceParamDim-1 }, { 0, 0 } );
            mTau.matrix_data() = aParamPoint( mNSpaceParamDim );

            // set input values for geometry interpolator
            mGeometryInterpolator->set_space_time( aParamPoint );

            // set bool for evaluation
            mNEval      = true;
            mBxEval     = true;
            md2Ndx2Eval = true;
            md3Ndx3Eval = true;
            mBtEval     = true;
            md2Ndt2Eval = true;
        }

//------------------------------------------------------------------------------

        void Field_Interpolator::set_coeff( const Matrix< DDRMat > & aUHat )
        {
            // check the input size
            MORIS_ASSERT( ( ( aUHat.n_cols() == mNumberOfFields ) && ( aUHat.n_rows() == mNFieldBases ) ),
                          "Field_Interpolator::set_coeff - Wrong input size (aUHat).");

            // set the coefficients
            mUHat = aUHat;
        }

//------------------------------------------------------------------------------

        const Matrix< DDRMat > & Field_Interpolator::N()
        {
            // if shape functions need to be evaluated
            if( mNEval )
            {
                // evaluate the shape functions
                this->eval_N();
            }

            // return member value
            return mN;
        }

         void Field_Interpolator::eval_N()
         {
             // check that mXi and mTau are set
             MORIS_ASSERT( mXi.numel()  > 0, "Field_Interpolator::eval_N - mXi  is not set." );
             MORIS_ASSERT( mTau.numel() > 0, "Field_Interpolator::eval_N - mTau is not set." );

             //evaluate space and time SF at Xi, Tau
             Matrix < DDRMat > tNSpace;
             Matrix < DDRMat > tNTime;
             mSpaceInterpolation->eval_N( mXi, tNSpace );
             mTimeInterpolation ->eval_N( mTau, tNTime );

             //evaluate space time SF by multiplying space and time SF
             mN = reshape( trans( tNSpace ) * tNTime, 1, mNFieldBases );

             // set bool for evaluation
             mNEval = false;
         }

//------------------------------------------------------------------------------
         const Matrix< DDRMat > & Field_Interpolator::dnNdxn( const uint & aDerivativeOrder )
         {
             // switch on derivative order
             switch ( aDerivativeOrder )
             {
                 // 1st order spatial derivatives of the shape functions
                 case( 1 ) :
                 {
                     // if dNdx needs to be evaluated
                     if( mBxEval )
                     {
                         // evaluate Bx
                         this->eval_Bx();
                     }
                     // return member data
                     return mBx;
                 }
                 // 2nd order spatial derivatives of the shape functions
                 case ( 2 ) :
                 {
                     // if d2Ndx2 needs to be evaluated
                     if( md2Ndx2Eval )
                     {
                         // evaluate d2Ndx2
                         this->eval_d2Ndx2();
                     }
                     // return member data
                     return md2Ndx2;
                 }
                 case ( 3 ) :
                 {
                     // if d3Ndx3 needs to be evaluated
                     if( md3Ndx3Eval )
                     {
                         // evaluate d3Ndx3
                         this->eval_d3Ndx3();
                     }
                     // return member data
                     return md3Ndx3;
                 }

                 default :
                     MORIS_ERROR( false, " Field_Interpolator::dnNdxn - Derivative order not implemented. " );
                     return mBx;
                     break;
             }
         }

//------------------------------------------------------------------------------

        const Matrix< DDRMat > & Field_Interpolator::Bx()
        {
            if( mBxEval )
            {
                // evaluate Bx
                this->eval_Bx();
            }
            return mBx;
        }

        void Field_Interpolator::eval_Bx()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel()>0,  "Field_Interpolator::eval_Bx - mXi  is not set." );
            MORIS_ASSERT( mTau.numel()>0, "Field_Interpolator::eval_Bx - mTau is not set." );

            // evaluate dNSpacedXi for the field time interpolation and transpose
            Matrix< DDRMat> tdNSpacedXi;
            mSpaceInterpolation->eval_dNdXi( mXi, tdNSpacedXi );
            tdNSpacedXi = trans( tdNSpacedXi );

            // evaluate NTime for the field time interpolation
            Matrix < DDRMat > tNTime;
            mTimeInterpolation->eval_N( mTau, tNTime );

            // set size dNFielddXi for the field
            Matrix< DDRMat> tdNFielddXi ( mNSpaceDim, mNFieldBases );

            // build the space time dNFielddXi row by row
            for ( moris::uint Ik = 0; Ik < mNSpaceDim; Ik++ )
            {
                tdNFielddXi.get_row( Ik ) = reshape( tdNSpacedXi.get_column(Ik) * tNTime , 1, mNFieldBases );
            }

            // evaluate the space Jacobian from the geometry interpolator
            Matrix< DDRMat > tdNGeodXi;
            mGeometryInterpolator->dNdXi( tdNGeodXi );
            Matrix< DDRMat > tJGeot;
            mGeometryInterpolator->space_jacobian( tdNGeodXi, tJGeot );

            // compute first derivative of the SF wrt x
            mBx = solve( tJGeot, tdNFielddXi );

            // set bool for evaluation
            mBxEval = false;
        }

//------------------------------------------------------------------------------

        const Matrix< DDRMat > & Field_Interpolator::d2Ndx2()
        {
            // if d2Ndx2 needs to be evaluated
            if( md2Ndx2Eval )
            {
                // evaluate the shape functions second order derivatives
                this->eval_d2Ndx2();
            }

            // return member data
            return md2Ndx2;
        }

        void Field_Interpolator::eval_d2Ndx2()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel()>0,  "Field_Interpolator::eval_d2Ndx2 - mXi  is not set." );
            MORIS_ASSERT( mTau.numel()>0, "Field_Interpolator::eval_d2Ndx2 - mTau is not set." );

            // get first and second derivatives of space SF wrt xi
            Matrix< DDRMat > tdNGeodxi;
            mGeometryInterpolator->dNdXi( tdNGeodxi );
            Matrix< DDRMat > td2NGeodxi2;
            mGeometryInterpolator->d2NdXi2( td2NGeodxi2 );

            // get matrices for second space derivatives from geometry interpolator
            Matrix< DDRMat > tJGeot, tKGeot, tLGeot;
            mGeometryInterpolator->space_jacobian_and_matrices_for_second_derivatives( tJGeot,
                                                                                       tKGeot,
                                                                                       tLGeot,
                                                                                       tdNGeodxi,
                                                                                       td2NGeodxi2 );

            // get the derivatives of the space time SF wrt x
            Matrix< DDRMat > tdNFielddx = this->Bx();

            // evaluate N for the field time interpolation
            Matrix< DDRMat > tNTime;
            mTimeInterpolation->eval_N( mTau, tNTime );

            // evaluate d2Ndxi2 for the field space interpolation
            Matrix< DDRMat > td2NSpacedxi2;
            mSpaceInterpolation->eval_d2NdXi2( mXi, td2NSpacedxi2 );
            td2NSpacedxi2 = trans( td2NSpacedxi2 );

            // get the number of rows for td2NFielddxi2
            uint tNSecondDerivatives = td2NGeodxi2.n_rows();

            // compute td2NFielddxi2 row by row
            Matrix< DDRMat > td2NFielddxi2( tNSecondDerivatives, mNFieldBases );
            for ( moris::uint Ik = 0; Ik < tNSecondDerivatives; Ik++ )
            {
                td2NFielddxi2.get_row( Ik ) = reshape( td2NSpacedxi2.get_column(Ik) * tNTime , 1, mNFieldBases );
            }

            //build the second derivatives of the space time SF wrt x
            Matrix< DDRMat > td2NFielddXi2 = td2NFielddxi2 - tKGeot * tdNFielddx;
            md2Ndx2 = solve( tLGeot, td2NFielddXi2 );

            // set bool for evaluation
            md2Ndx2Eval = false;
        }

//------------------------------------------------------------------------------

        const Matrix< DDRMat > & Field_Interpolator::d3Ndx3()
        {
            // if d3Ndx3 needs to be evaluated
            if( md3Ndx3Eval )
            {
                // evaluate the shape functions third order derivatives
                this->eval_d3Ndx3();
            }

            // return member date
            return md3Ndx3;
        }

        //FIXME: Function not complete
        void Field_Interpolator::eval_d3Ndx3()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel()  > 0, "Field_Interpolator::eval_d3Ndx3 - mXi is not set." );
            MORIS_ASSERT( mTau.numel() > 0, "Field_Interpolator::eval_d3Ndx3 - mTau is not set." );

            // get first and second derivatives of space SF wrt xi
            Matrix< DDRMat > tdNGeodxi;
            mGeometryInterpolator->dNdXi( tdNGeodxi );
            Matrix< DDRMat > td2NGeodxi2;
            mGeometryInterpolator->d2NdXi2( td2NGeodxi2 );
            Matrix< DDRMat > td3NGeodxi3;
            mGeometryInterpolator->d3NdXi3( td3NGeodxi3 );

            // get matrices for second space derivatives from geometry interpolator
            Matrix< DDRMat > tJGeot, tJ2bGeot, tJ3aGeot, tJ3bGeot, tJ3cGeot;
            mGeometryInterpolator->space_jacobian_and_matrices_for_third_derivatives( tJGeot,
                                                                                      tJ2bGeot,
                                                                                      tJ3aGeot,
                                                                                      tJ3bGeot,
                                                                                      tJ3cGeot,
                                                                                      tdNGeodxi,
                                                                                      td2NGeodxi2,
                                                                                      td3NGeodxi3 );

            // get the derivatives of the space time SF wrt x
            Matrix< DDRMat > tdNFielddx   = this->Bx();
            Matrix< DDRMat > td2NFielddx2 = this->d2Ndx2();

            Matrix< DDRMat > tNTime;

            // evaluate N for the field time interpolation
            mTimeInterpolation->eval_N( mTau, tNTime );

            // evaluate derivatives of the field space interpolation
            Matrix< DDRMat > td2NSpacedxi2;
            mSpaceInterpolation->eval_d2NdXi2( mXi, td2NSpacedxi2 );
            td2NSpacedxi2 = trans( td2NSpacedxi2 );
            Matrix< DDRMat > td3NSpacedxi3;
            mSpaceInterpolation->eval_d3NdXi3( mXi, td3NSpacedxi3 );
            td3NSpacedxi3 = trans( td3NSpacedxi3 );

            // get the number 3rd derivatives
            uint tNThirdDerivatives  = td3NGeodxi3.n_rows();

            // compute td3NFielddxi3 row by row
            Matrix< DDRMat > td3NFielddxi3( tNThirdDerivatives, mNFieldBases );
            for ( moris::uint Ik = 0; Ik < tNThirdDerivatives; Ik++ )
            {
                td3NFielddxi3.get_row( Ik ) = reshape( td3NSpacedxi3.get_column(Ik) * tNTime , 1, mNFieldBases );
            }

            //build the third derivatives of the space time SF wrt x
            Matrix< DDRMat > td3NFielddXi3 = td3NFielddxi3 - tJ3bGeot * td2NFielddx2 - tJ3cGeot * tdNFielddx;
            md3Ndx3 = solve( tJ3aGeot, td3NFielddXi3 );

            // set bool for evaluation
            md3Ndx3Eval = false;
        }

//------------------------------------------------------------------------------

        const Matrix< DDRMat > & Field_Interpolator::Bt()
        {
            // if Bt needs to be evaluated
            if( mBtEval )
            {
                // evaluate the shape functions first order derivatives
                this->eval_Bt();
            }

            // return member date
            return mBt;
        }

        void Field_Interpolator::eval_Bt()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel()>0,  "Field_Interpolator::eval_Bt - mXi  is not set." );
            MORIS_ASSERT( mTau.numel()>0, "Field_Interpolator::eval_Bt - mTau is not set." );

            // evaluate dNdTau for the field time interpolation
            Matrix< DDRMat > tdNTimedTau;
            mTimeInterpolation->eval_dNdXi( mTau, tdNTimedTau );

            Matrix < DDRMat > tNSpace;

            // evaluate N for the field space interpolation
            mSpaceInterpolation->eval_N( mXi, tNSpace );
            tNSpace = trans( tNSpace );

            // set size tdNFielddTau for the field
            Matrix< DDRMat> tdNFielddTau ( mNTimeDim, mNFieldBases );

            // build the space time dNdTau row by row
            for ( moris::uint Ik = 0; Ik < mNTimeDim; Ik++ )
            {
                tdNFielddTau.get_row( Ik ) = reshape( tNSpace * tdNTimedTau.get_row(Ik), 1, mNFieldBases);
            }

            // evaluate the Jacobian from the space geometry interpolator
            Matrix< DDRMat > tdNGeodTau;
            mGeometryInterpolator->dNdTau( tdNGeodTau );
            Matrix< DDRMat > tJGeot;
            mGeometryInterpolator->time_jacobian( tdNGeodTau, tJGeot );

            // transform output matrix to dNdX
            mBt = solve( tJGeot, tdNFielddTau );

            // set bool for evaluation
            mBtEval = false;
        }

//------------------------------------------------------------------------------

        const Matrix< DDRMat > & Field_Interpolator::d2Ndt2()
        {
            // if d2Ndt2 needs to be evaluated
            if( md2Ndt2Eval )
            {
                // evaluate the shape functions second order derivatives
                this->eval_d2Ndt2();
            }

            // return member date
            return md2Ndt2;
        }

        void Field_Interpolator::eval_d2Ndt2()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel()>0,  "Field_Interpolator::eval_d2Ndt2 - mXi  is not set." );
            MORIS_ASSERT( mTau.numel()>0, "Field_Interpolator::eval_d2Ndt2 - mTau is not set." );

            // get first and second derivatives of space SF wrt tau
            Matrix< DDRMat > tdNGeodtau;
            mGeometryInterpolator->dNdTau( tdNGeodtau );
            Matrix< DDRMat > td2NGeodtau2;
            mGeometryInterpolator->d2NdTau2( td2NGeodtau2 );

            // get matrices for second derivatives from space geometry interpolator
            Matrix< DDRMat > tJGeot, tKGeot, tLGeot;
            mGeometryInterpolator->time_jacobian_and_matrices_for_second_derivatives( tJGeot,
                                                                                      tKGeot,
                                                                                      tLGeot,
                                                                                      tdNGeodtau,
                                                                                      td2NGeodtau2 );

            // get the derivatives of the space time SF wrt t
            Matrix< DDRMat > tdNFielddt = this->Bt();

            //get the second derivatives of the space time SF wrt tau
            Matrix< DDRMat > tNSpace;

            // get N for the field space interpolation and transpose
            mSpaceInterpolation->eval_N( mXi, tNSpace );
            tNSpace = trans( tNSpace );

            // get d2Ndtau2 for the field time interpolation
            Matrix< DDRMat > td2NTimedtau2;
            mTimeInterpolation->eval_d2NdXi2( mTau, td2NTimedtau2 );

            // get the number of rows for td2NFielddtau2
            uint tNSecondDerivatives = td2NGeodtau2.n_rows();

            // compute second derivatives of the space time SF wrt tau row by row
            Matrix< DDRMat > td2NFielddtau2( tNSecondDerivatives, mNFieldBases );
            for ( uint Ik = 0; Ik < tNSecondDerivatives; Ik++ )
            {
                td2NFielddtau2.get_row(Ik) = reshape( tNSpace * td2NTimedtau2.get_row(Ik) , 1, mNFieldBases );
            }

            //build the second derivatives of the space time SF wrt t
            Matrix< DDRMat > td2NFielddTau2 = td2NFielddtau2 - tKGeot * tdNFielddt;
            md2Ndt2 = solve( tLGeot, td2NFielddTau2 );

            // set bool for evaluation
            md2Ndt2Eval = false;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Field_Interpolator::val()
        {
            // check that mUHat is set
            MORIS_ASSERT( mUHat.numel() > 0, "Field_Interpolator::val - mUHat is not set." );

            //evaluate the field value
            return this->N() * mUHat ;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Field_Interpolator::gradx( const uint & aDerivativeOrder )
        {
            // check that mUHat is set
            MORIS_ASSERT( mUHat.numel()>0,  "Field_Interpolator::gradx - mUHat  is not set." );

            switch ( aDerivativeOrder )
            {
                case( 1 ) :
                    //evaluate the field first space derivative
                    return this->Bx() * mUHat ;
                    break;

                case ( 2 ) :
                    //evaluate the field second space derivative
                    return this->d2Ndx2() * mUHat ;
                    break;

                case ( 3 ) :
                    //evaluate the field second space derivative
                    return this->d3Ndx3() * mUHat ;
                    break;

                default :
                    MORIS_ERROR( false, " Field_Interpolator::gradx - Derivative order not implemented. " );
                    Matrix< DDRMat > tEmpty;
                    return tEmpty;
                    break;
            }
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Field_Interpolator::gradt( const uint & aDerivativeOrder )
        {
            // check that mUHat is set
            MORIS_ASSERT( mUHat.numel()>0,  "Field_Interpolator::gradt - mUHat  is not set." );

            switch ( aDerivativeOrder )
            {
                case( 1 ) :
                    //evaluate the field first time derivative
                    return this->Bt() * mUHat ;
                    break;

                case ( 2 ) :
                    //evaluate the field second time derivative
                    return this->d2Ndt2() * mUHat ;
                    break;

                default :
                    MORIS_ERROR( false, " Field_Interpolator::gradt - Derivative order not implemented. " );
                    Matrix< DDRMat > tEmpty;
                    return tEmpty;
                    break;
            }
        }

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */


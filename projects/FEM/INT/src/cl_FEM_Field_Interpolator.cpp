
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
#include "cl_FEM_Property.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        Field_Interpolator::Field_Interpolator( const uint                   & aNumberOfFields,
                                                const Interpolation_Rule     & aFieldInterpolationRule,
                                                const Geometry_Interpolator*   aGeometryInterpolator,
                                                const MSI::Dof_Type            aDofType )
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

        Field_Interpolator::Field_Interpolator( const uint                   & aNumberOfFields,
                                                const Interpolation_Rule     & aFieldInterpolationRule,
                                                const Geometry_Interpolator*   aGeometryInterpolator,
                                                const Property *               aProperty,
                                                const fem::Property_Type       aPropertyType )
                                              : mNumberOfFields( aNumberOfFields ),
                                                mGeometryInterpolator( aGeometryInterpolator ),
                                                mProperty( aProperty ),
                                                mPropertyType( aPropertyType )
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
            //print(aParamPoint, "aParamPoint" );

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

         Matrix < DDRMat> Field_Interpolator::N()
         {
             // check that mXi and mTau are set
             MORIS_ASSERT( mXi.numel()  > 0, "Field_Interpolator::N - mXi  is not set." );
             MORIS_ASSERT( mTau.numel() > 0, "Field_Interpolator::N - mTau is not set." );

             Matrix < DDRMat > tNSpace;
             Matrix < DDRMat > tNTime;

             //evaluate space and time SF at Xi, Tau
             mSpaceInterpolation->eval_N( mXi, tNSpace );
             mTimeInterpolation ->eval_N( mTau, tNTime );

             //evaluate space time SF by multiplying space and time SF
             return reshape( trans( tNSpace ) * tNTime, 1, mNFieldBases );
         }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Field_Interpolator::Bx()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel()>0,  "Field_Interpolator::Bx - mXi  is not set." );
            MORIS_ASSERT( mTau.numel()>0, "Field_Interpolator::Bx - mTau is not set." );

            // evaluate dNSpacedXi for the field time interpolation and transpose
            Matrix< DDRMat> tdNSpacedXi = mSpaceInterpolation->eval_dNdXi( mXi );
            tdNSpacedXi = trans( tdNSpacedXi );
            Matrix < DDRMat > tNTime;

            // evaluate NTime for the field time interpolation
            mTimeInterpolation->eval_N( mTau, tNTime );

            // set size dNFielddXi for the field
            Matrix< DDRMat> tdNFielddXi ( mNSpaceDim, mNFieldBases );

            // build the space time dNFielddXi row by row
            for ( moris::uint Ik = 0; Ik < mNSpaceDim; Ik++ )
            {
                tdNFielddXi.get_row(Ik) = reshape( tdNSpacedXi.get_column(Ik) * tNTime , 1, mNFieldBases );
            }

            // evaluate the space Jacobian from the geometry interpolator
            Matrix< DDRMat > tdNGeodXi = mGeometryInterpolator->dNdXi( mXi );
            Matrix< DDRMat > tJGeot    = mGeometryInterpolator->space_jacobian( tdNGeodXi );

            // compute first derivative of the SF wrt x
            return solve( tJGeot, tdNFielddXi );
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Field_Interpolator::eval_d2Ndx2()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel()>0,  "Field_Interpolator::eval_d2Ndx2 - mXi  is not set." );
            MORIS_ASSERT( mTau.numel()>0, "Field_Interpolator::eval_d2Ndx2 - mTau is not set." );

            // get first and second derivatives of space SF wrt xi
            Matrix< DDRMat > tdNGeodxi   = mGeometryInterpolator->dNdXi( mXi );
            Matrix< DDRMat > td2NGeodxi2 = mGeometryInterpolator->d2NdXi2( mXi );

            // get matrices for second space derivatives from geometry interpolator
            Matrix< DDRMat > tJGeot, tKGeot, tLGeot;
            mGeometryInterpolator->space_jacobian_and_matrices_for_second_derivatives( tJGeot,
                                                                                       tKGeot,
                                                                                       tLGeot,
                                                                                       tdNGeodxi,
                                                                                       td2NGeodxi2 );

            // get the derivatives of the space time SF wrt x
            Matrix< DDRMat > tdNFielddx = this->Bx();

            Matrix< DDRMat > tNTime;

            // evaluate N for the field time interpolation
           mTimeInterpolation->eval_N( mTau, tNTime );

            // evaluate d2Ndxi2 for the field space interpolation
            Matrix< DDRMat > td2NSpacedxi2 = mSpaceInterpolation->eval_d2NdXi2( mXi );
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
            return solve( tLGeot, td2NFielddXi2 );
        }

//------------------------------------------------------------------------------
//FIXME: Function not complete

        Matrix< DDRMat > Field_Interpolator::eval_d3Ndx3()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel()  > 0, "Field_Interpolator::eval_d3Ndx3 - mXi is not set." );
            MORIS_ASSERT( mTau.numel() > 0, "Field_Interpolator::eval_d3Ndx3 - mTau is not set." );

            // get first and second derivatives of space SF wrt xi
            Matrix< DDRMat > tdNGeodxi   = mGeometryInterpolator->dNdXi( mXi );
            Matrix< DDRMat > td2NGeodxi2 = mGeometryInterpolator->d2NdXi2( mXi );
            Matrix< DDRMat > td3NGeodxi3 = mGeometryInterpolator->d3NdXi3( mXi );

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
            Matrix< DDRMat > tdNFielddx = this->Bx();
            Matrix< DDRMat > td2NFielddx2 = this->eval_d2Ndx2();

            Matrix< DDRMat > tNTime;

            // evaluate N for the field time interpolation
            mTimeInterpolation->eval_N( mTau, tNTime );

            // evaluate derivatives of the field space interpolation
            Matrix< DDRMat > td2NSpacedxi2 = mSpaceInterpolation->eval_d2NdXi2( mXi );
            td2NSpacedxi2 = trans( td2NSpacedxi2 );
            Matrix< DDRMat > td3NSpacedxi3 = mSpaceInterpolation->eval_d3NdXi3( mXi );
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
            return solve( tJ3aGeot, td3NFielddXi3 );
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Field_Interpolator::Bt()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel()>0,  "Field_Interpolator::Bt - mXi  is not set." );
            MORIS_ASSERT( mTau.numel()>0, "Field_Interpolator::Bt - mTau is not set." );

            // evaluate dNdTau for the field time interpolation
            Matrix< DDRMat > tdNTimedTau = mTimeInterpolation->eval_dNdXi( mTau );

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
            Matrix< DDRMat > tdNGeodTau = mGeometryInterpolator->dNdTau( mTau );
            Matrix< DDRMat > tJGeot = mGeometryInterpolator->time_jacobian( tdNGeodTau );

            // transform output matrix to dNdX
            return solve( tJGeot, tdNFielddTau );
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Field_Interpolator::eval_d2Ndt2()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel()>0,  "Field_Interpolator::eval_d2Ndt2 - mXi  is not set." );
            MORIS_ASSERT( mTau.numel()>0, "Field_Interpolator::eval_d2Ndt2 - mTau is not set." );

            // get first and second derivatives of space SF wrt tau
            Matrix< DDRMat > tdNGeodtau   = mGeometryInterpolator->dNdTau( mTau );
            Matrix< DDRMat > td2NGeodtau2 = mGeometryInterpolator->d2NdTau2( mTau );

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
            Matrix< DDRMat > td2NTimedtau2 = mTimeInterpolation->eval_d2NdXi2( mTau );

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
            return solve( tLGeot, td2NFielddTau2 );
        }
//------------------------------------------------------------------------------

        Matrix< DDRMat > Field_Interpolator::val()
        {
            // check that mUHat is set
            MORIS_ASSERT( mUHat.numel() > 0, "Field_Interpolator::val - mUHat  is not set." );

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
                    return this->eval_d2Ndx2() * mUHat ;
                    break;

                case ( 3 ) :
                    //evaluate the field second space derivative
                    return this->eval_d3Ndx3() * mUHat ;
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
                    return this->eval_d2Ndt2() * mUHat ;
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


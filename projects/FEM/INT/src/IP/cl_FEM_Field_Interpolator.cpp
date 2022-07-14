
#include "cl_Matrix.hpp"    //LNA/src
#include "linalg_typedefs.hpp"
#include "fn_inv.hpp"
#include "fn_det.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_vectorize.hpp"
#include "fn_det.hpp"
#include "fn_sum.hpp"
#include "fn_diag_vec.hpp"
#include "op_times.hpp"
#include "op_equal_equal.hpp"
#include "op_less_equal.hpp"
#include "op_greater_equal.hpp"


#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_MTK_Enums.hpp"                 //MTK/src

#include <iostream>

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        Field_Interpolator::Field_Interpolator(
                const uint&                        aNumberOfFields,
                const mtk::Interpolation_Rule&     aFieldInterpolationRule,
                Geometry_Interpolator*             aGeometryInterpolator,
                const moris::Cell< MSI::Dof_Type > aDofType )
                : mNumberOfFields( aNumberOfFields )
                , mGeometryInterpolator( aGeometryInterpolator )
                , mDofType( aDofType )
        {
            // create space and time interpolation function
            mSpaceInterpolation = aFieldInterpolationRule.create_space_interpolation_function();
            mTimeInterpolation  = aFieldInterpolationRule.create_time_interpolation_function();

            // get number of space, time dimensions
            mNSpaceDim = mSpaceInterpolation->get_number_of_dimensions();
            mNTimeDim  = mTimeInterpolation->get_number_of_dimensions();

            // get number of space parametric dimensions
            mNSpaceParamDim = mSpaceInterpolation->get_number_of_param_dimensions();

            // check dimensions consistency
            MORIS_ERROR( ( mNSpaceDim == mGeometryInterpolator->get_number_of_space_dimensions() ),
                    "Field_Interpolator - Space dimension inconsistency." );
            MORIS_ERROR( ( mNTimeDim == mGeometryInterpolator->get_number_of_time_dimensions() ),
                    "Field_Interpolator - Time dimension inconsistency." );

            // get number of space, time, and space time basis
            mNSpaceBases = mSpaceInterpolation->get_number_of_bases();
            mNTimeBases  = mTimeInterpolation->get_number_of_bases();
            mNFieldBases = mNSpaceBases * mNTimeBases;

            // get number of coefficients
            mNFieldCoeff = mNFieldBases * mNumberOfFields;

            // init storage size
            mN.set_size( mNumberOfFields, mNFieldCoeff, 0.0 );
            mdNdx.set_size( mNSpaceDim, mNFieldBases, 0.0 );
            mdNdt.set_size( mNTimeDim, mNFieldBases, 0.0 );
            md2Ndxt.set_size( mNSpaceDim, mNFieldBases, 0.0 );
        }

        //------------------------------------------------------------------------------

        Field_Interpolator::Field_Interpolator(
                const uint&                    aNumberOfFields,
                const mtk::Interpolation_Rule& aFieldInterpolationRule,
                Geometry_Interpolator*         aGeometryInterpolator,
                const moris::Cell< PDV_Type >  aDvType )
                : mNumberOfFields( aNumberOfFields )
                , mGeometryInterpolator( aGeometryInterpolator )
                , mDvType( aDvType )
        {
            // create space and time interpolation function
            mSpaceInterpolation = aFieldInterpolationRule.create_space_interpolation_function();
            mTimeInterpolation  = aFieldInterpolationRule.create_time_interpolation_function();

            // get number of space, time dimensions
            mNSpaceDim = mSpaceInterpolation->get_number_of_dimensions();
            mNTimeDim  = mTimeInterpolation->get_number_of_dimensions();

            // get number of space parametric dimensions
            mNSpaceParamDim = mSpaceInterpolation->get_number_of_param_dimensions();

            // check dimensions consistency
            MORIS_ERROR( ( mNSpaceDim == mGeometryInterpolator->get_number_of_space_dimensions() ),
                    "Field_Interpolator - Space dimension inconsistency." );
            MORIS_ERROR( ( mNTimeDim == mGeometryInterpolator->get_number_of_time_dimensions() ),
                    "Field_Interpolator - Time dimension inconsistency." );

            // get number of space, time, and space time basis
            mNSpaceBases = mSpaceInterpolation->get_number_of_bases();
            mNTimeBases  = mTimeInterpolation->get_number_of_bases();
            mNFieldBases = mNSpaceBases * mNTimeBases;

            // get number of coefficients
            mNFieldCoeff = mNFieldBases * mNumberOfFields;

            // init storage size
            mN.set_size( mNumberOfFields, mNFieldCoeff, 0.0 );
            mdNdx.set_size( mNSpaceDim, mNFieldBases, 0.0 );
            mdNdt.set_size( mNTimeDim, mNFieldBases, 0.0 );
            md2Ndxt.set_size( mNSpaceDim, mNFieldBases, 0.0 );
        }

        //------------------------------------------------------------------------------

        Field_Interpolator::Field_Interpolator(
                const uint&                          aNumberOfFields,
                const mtk::Interpolation_Rule&       aFieldInterpolationRule,
                Geometry_Interpolator*               aGeometryInterpolator,
                const moris::Cell< mtk::Field_Type > aFieldType )
                : mNumberOfFields( aNumberOfFields )
                , mGeometryInterpolator( aGeometryInterpolator )
                , mFieldType( aFieldType )
        {
            // create space and time interpolation function
            mSpaceInterpolation = aFieldInterpolationRule.create_space_interpolation_function();
            mTimeInterpolation  = aFieldInterpolationRule.create_time_interpolation_function();

            // get number of space, time dimensions
            mNSpaceDim = mSpaceInterpolation->get_number_of_dimensions();
            mNTimeDim  = mTimeInterpolation->get_number_of_dimensions();

            // get number of space parametric dimensions
            mNSpaceParamDim = mSpaceInterpolation->get_number_of_param_dimensions();

            // check dimensions consistency
            MORIS_ERROR( ( mNSpaceDim == mGeometryInterpolator->get_number_of_space_dimensions() ),
                    "Field_Interpolator - Space dimension inconsistency." );
            MORIS_ERROR( ( mNTimeDim == mGeometryInterpolator->get_number_of_time_dimensions() ),
                    "Field_Interpolator - Time dimension inconsistency." );

            // get number of space, time, and space time basis
            mNSpaceBases = mSpaceInterpolation->get_number_of_bases();
            mNTimeBases  = mTimeInterpolation->get_number_of_bases();
            mNFieldBases = mNSpaceBases * mNTimeBases;

            // get number of coefficients
            mNFieldCoeff = mNFieldBases * mNumberOfFields;

            // init storage size
            mN.set_size( mNumberOfFields, mNFieldCoeff, 0.0 );
            mdNdx.set_size( mNSpaceDim, mNFieldBases, 0.0 );
            mdNdt.set_size( mNTimeDim, mNFieldBases, 0.0 );
            md2Ndxt.set_size( mNSpaceDim, mNFieldBases, 0.0 );
        }

        //------------------------------------------------------------------------------

        Field_Interpolator::~Field_Interpolator()
        {
            // delete interpolation functions
            if ( mSpaceInterpolation != NULL )
            {
                delete mSpaceInterpolation;
            }

            if ( mTimeInterpolation != NULL )
            {
                delete mTimeInterpolation;
            }
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::reset_eval_flags()
        {
            // reset bool for evaluation
            mNEval           = true;
            mNTransEval      = true;
            mNBuildEval      = true;
            mdNdxEval        = true;
            md2Ndx2Eval      = true;
            md3Ndx3Eval      = true;
            mdNdtEval        = true;
            md2Ndt2Eval      = true;
            md2NdxtEval      = true;
            mDivOperatorEval = true;

            mValEval      = true;
            mValTransEval = true;

            mGradx1Eval = true;
            mGradx2Eval = true;
            mGradx3Eval = true;

            mGradt1Eval = true;
            mGradt2Eval = true;
            mGradt3Eval = true;

            mGradxtEval = true;
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::reset_eval_flags_coefficients()
        {
            // reset bool for evaluation
            mValEval      = true;
            mValTransEval = true;

            mGradx1Eval = true;
            mGradx2Eval = true;
            mGradx3Eval = true;

            mGradt1Eval = true;
            mGradt2Eval = true;
            mGradt3Eval = true;

            mGradxtEval = true;
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::set_discretization_mesh_index( const moris_index aDiscretizationMeshIndex )
        {
            mDiscretizationMeshIndex = aDiscretizationMeshIndex;
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::set_space_time( const Matrix< DDRMat >& aParamPoint )
        {
            // check input size aParamPoint
            MORIS_ASSERT( ( ( aParamPoint.n_cols() == 1 ) && ( aParamPoint.n_rows() == mNSpaceParamDim + mNTimeDim ) ),
                    "Field_Interpolator::set_space_time - Wrong input size ( aParamPoint )." );

            // check input values are between -1 and 1
            for ( uint Ik = 0; Ik < mNSpaceParamDim; Ik++ )
            {
                switch ( mGeometryInterpolator->get_space_geometry_type() )
                {
                    case mtk::Geometry_Type::LINE:
                    case mtk::Geometry_Type::QUAD:
                    case mtk::Geometry_Type::HEX:
                    {
                        MORIS_ASSERT( ( ( aParamPoint( Ik ) <= 1.0 + mEpsilon ) && ( aParamPoint( Ik ) >= -1.0 - mEpsilon ) ),
                                "Field_Interpolator::set_space_time - Wrong input value space line/quad/hex ( aParamPoint ): %f \n",
                                aParamPoint( Ik ) );
                        break;
                    }

                    case mtk::Geometry_Type::TRI:
                    case mtk::Geometry_Type::TET:
                    {
                        MORIS_ASSERT( ( ( aParamPoint( Ik ) <= 1.0 + mEpsilon ) && ( aParamPoint( Ik ) >= 0.0 - mEpsilon ) ),
                                "Field_Interpolator::set_space_time - Wrong input value space tri/tet ( aParamPoint ): %f \n",
                                aParamPoint( Ik ) );
                        break;
                    }

                    default:
                        MORIS_ERROR( false, "Field_Interpolator::set_space_time - unknown geometry type." );
                }
            }

            MORIS_ASSERT( ( ( aParamPoint( mNSpaceParamDim ) <= 1.0 + mEpsilon ) && ( aParamPoint( mNSpaceParamDim ) >= -1.0 - mEpsilon ) ),
                    "Field_Interpolator::set_space_time - Wrong input value time line ( aParamPoint ) " );

            // set input values
            mXi  = aParamPoint( { 0, mNSpaceParamDim - 1 }, { 0, 0 } );
            mTau = aParamPoint( mNSpaceParamDim );

            // reset bool for evaluation
            this->reset_eval_flags();
            this->reset_eval_flags_coefficients();
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::set_coeff( const Matrix< DDRMat >& aUHat )
        {
            // check the input size
            MORIS_ASSERT( ( ( aUHat.n_cols() == mNumberOfFields ) && ( aUHat.n_rows() == mNFieldBases ) ),
                    "Field_Interpolator::set_coeff - Wrong input size (aUHat)." );

            // set the coefficients
            mUHat = aUHat;

            // reset bool for evaluation
            this->reset_eval_flags_coefficients();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Field_Interpolator::NBuild()
        {
            // if shape functions need to be evaluated
            if ( mNBuildEval )
            {
                // evaluate the shape functions
                this->eval_NBuild();

                // set bool for evaluation
                mNBuildEval = false;
            }

            // return member value
            return mNBuild;
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::eval_NBuild()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel() > 0,
                    "Field_Interpolator::eval_NBuild - mXi  is not set." );

            MORIS_ASSERT( mTau.numel() > 0,
                    "Field_Interpolator::eval_NBuild - mTau is not set." );

            // evaluate space and time SF at Xi, Tau
            Matrix< DDRMat > tNSpace;
            Matrix< DDRMat > tNTime;

            mSpaceInterpolation->eval_N( mXi, tNSpace );
            mTimeInterpolation->eval_N( mTau, tNTime );

            // evaluate space time SF by multiplying space and time SF and create row vector
            mNBuild = trans( vectorize( trans( tNSpace ) * tNTime ) );
        }

        //------------------------------------------------------------------------------

        const Matrix< SDRMat >&
        Field_Interpolator::N()
        {
            // if shape functions need to be evaluated
            if ( mNEval )
            {
                // evaluate the shape functions
                this->eval_N();

                // set bool for evaluation
                mNEval = false;
            }

            // return member value
            return mN;
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::eval_N()
        {
            // loop over the fields
            for ( uint iField = 0; iField < mNumberOfFields; iField++ )
            {
                // fill the matrix for each field
                mN(
                        { iField, iField },
                        { iField * mNFieldBases, ( iField + 1 ) * mNFieldBases - 1 } ) =
                        this->NBuild().matrix_data();
            }
        }
        //------------------------------------------------------------------------------

        const Matrix< SDRMat >&
        Field_Interpolator::N_trans()
        {
            // if shape functions need to be evaluated
            if ( mNTransEval )
            {
                // evaluate the shape functions
                mNTrans = trans( this->N() );

                // set bool for evaluation
                mNTransEval = false;
            }

            // return member value
            return mNTrans;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Field_Interpolator::dnNdxn( const uint& aDerivativeOrder )
        {
            // switch on derivative order
            switch ( aDerivativeOrder )
            {
                // 1st order spatial derivatives of the shape functions
                case 1:
                {
                    // if dNdx needs to be evaluated
                    if ( mdNdxEval )
                    {
                        // evaluate d1Ndx1
                        this->eval_d1Ndx1();

                        // set bool for evaluation
                        mdNdxEval = false;
                    }
                    // return member data
                    return mdNdx;
                }
                // 2nd order spatial derivatives of the shape functions
                case 2:
                {
                    // if d2Ndx2 needs to be evaluated
                    if ( md2Ndx2Eval )
                    {
                        // evaluate d2Ndx2
                        this->eval_d2Ndx2();

                        // set bool for evaluation
                        md2Ndx2Eval = false;
                    }
                    // return member data
                    return md2Ndx2;
                }
                case 3:
                {
                    // if d3Ndx3 needs to be evaluated
                    if ( md3Ndx3Eval )
                    {
                        // evaluate d3Ndx3
                        this->eval_d3Ndx3();

                        // set bool for evaluation
                        md3Ndx3Eval = false;
                    }
                    // return member data
                    return md3Ndx3;
                }

                default:
                    MORIS_ERROR( false, " Field_Interpolator::dnNdxn - Derivative order not implemented. " );
                    return mdNdx;
            }
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::eval_d1Ndx1()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel() > 0,
                    "Field_Interpolator::eval_d1Ndx1 - mXi  is not set." );
            MORIS_ASSERT( mTau.numel() > 0,
                    "Field_Interpolator::eval_d1Ndx1 - mTau is not set." );

            // evaluate dNSpacedXi for the space interpolation
            Matrix< DDRMat > tdNSpacedXi;
            mSpaceInterpolation->eval_dNdXi( mXi, tdNSpacedXi );

            // evaluate the space Jacobian from the geometry interpolator
            const Matrix< DDRMat >& tInvJGeot = mGeometryInterpolator->inverse_space_jacobian();

            // compute first derivative of the space shape function wrt x
            auto tdNSpacedX = tInvJGeot * tdNSpacedXi;

            // evaluate NTime for the time interpolation
            Matrix< DDRMat > tNTime;
            mTimeInterpolation->eval_N( mTau, tNTime );

            // build the space time dNFielddXi row by row
            for ( moris::uint Ik = 0; Ik < mNTimeBases; Ik++ )
            {
                mdNdx(
                        { 0, mNSpaceDim - 1 },
                        { Ik * mNSpaceBases, ( Ik + 1 ) * mNSpaceBases - 1 } ) =
                        tdNSpacedX * tNTime( Ik );
            }
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::eval_d2Ndx2()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel() > 0,
                    "Field_Interpolator::eval_d2Ndx2 - mXi  is not set." );
            MORIS_ASSERT( mTau.numel() > 0,
                    "Field_Interpolator::eval_d2Ndx2 - mTau is not set." );

            // get matrices for second space derivatives from geometry interpolator
            Matrix< DDRMat > tJGeot, tKGeot, tLGeot;
            mGeometryInterpolator->space_jacobian_and_matrices_for_second_derivatives(
                    tJGeot,
                    tKGeot,
                    tLGeot,
                    mGeometryInterpolator->dNdXi(),
                    mGeometryInterpolator->d2NdXi2() );

            // compute first derivative of the field shape function wrt x
            const Matrix< DDRMat >& tdNFielddx = this->dnNdxn( 1 );

            // evaluate d2Ndxi2 for the field space interpolation
            Matrix< DDRMat > td2NSpacedxi2;
            mSpaceInterpolation->eval_d2NdXi2( mXi, td2NSpacedxi2 );

            // evaluate NTime for the time interpolation
            Matrix< DDRMat > tNTime;
            mTimeInterpolation->eval_N( mTau, tNTime );

            // set size d2NFielddxi2 for the field
            uint             tNumRows = td2NSpacedxi2.n_rows();
            Matrix< DDRMat > td2NFielddxi2( tNumRows, mNFieldBases, 0.0 );

            // build the space time d2NFielddxi2 row by row
            for ( moris::uint Ik = 0; Ik < mNTimeBases; Ik++ )
            {
                td2NFielddxi2(
                        { 0, tNumRows - 1 },
                        { Ik * mNSpaceBases, ( Ik + 1 ) * mNSpaceBases - 1 } ) +=
                        td2NSpacedxi2 * tNTime( Ik );
            }

            // build the second derivatives of the space time SF wrt x
            td2NFielddxi2 -= tKGeot * tdNFielddx;
            md2Ndx2 = inv( tLGeot ) * td2NFielddxi2;
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::eval_d3Ndx3()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel() > 0,
                    "Field_Interpolator::eval_d3Ndx3 - mXi is not set." );
            MORIS_ASSERT( mTau.numel() > 0,
                    "Field_Interpolator::eval_d3Ndx3 - mTau is not set." );

            // get matrices for second space derivatives from geometry interpolator
            Matrix< DDRMat > tJGeot, tJ2bGeot, tJ3aGeot, tJ3bGeot, tJ3cGeot;
            mGeometryInterpolator->space_jacobian_and_matrices_for_third_derivatives(
                    tJGeot,
                    tJ2bGeot,
                    tJ3aGeot,
                    tJ3bGeot,
                    tJ3cGeot,
                    mGeometryInterpolator->dNdXi(),
                    mGeometryInterpolator->d2NdXi2(),
                    mGeometryInterpolator->d3NdXi3() );

            // get the derivatives of the space time SF wrt x
            const Matrix< DDRMat >& tdNFielddx   = this->dnNdxn( 1 );
            const Matrix< DDRMat >& td2NFielddx2 = this->dnNdxn( 2 );

            // evaluate N for the field time interpolation
            Matrix< DDRMat > tNTime;
            mTimeInterpolation->eval_N( mTau, tNTime );

            // evaluate derivatives of the field space interpolation
            Matrix< DDRMat > td3NSpacedxi3;
            mSpaceInterpolation->eval_d3NdXi3( mXi, td3NSpacedxi3 );

            // set size for td3NFielddxi3
            uint             tNumRows = td3NSpacedxi3.n_rows();
            Matrix< DDRMat > td3NFielddxi3( tNumRows, mNFieldBases, 0.0 );

            // compute td3NFielddXi3
            for ( moris::uint Ik = 0; Ik < mNTimeBases; Ik++ )
            {
                td3NFielddxi3(
                        { 0, tNumRows - 1 },
                        { Ik * mNSpaceBases, ( Ik + 1 ) * mNSpaceBases - 1 } ) +=
                        td3NSpacedxi3 * tNTime( Ik );
            }

            // build the third derivatives of the space time SF wrt x
            td3NFielddxi3 -= tJ3bGeot * td2NFielddx2 + tJ3cGeot * tdNFielddx;
            md3Ndx3 = inv( tJ3aGeot ) * td3NFielddxi3;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Field_Interpolator::dnNdtn( const uint& aDerivativeOrder )
        {
            // switch on derivative order
            switch ( aDerivativeOrder )
            {
                // 1st order time derivatives of the shape functions
                case 1:
                {
                    // if dNdx needs to be evaluated
                    if ( mdNdtEval )
                    {
                        // evaluate d1Ndt1
                        this->eval_d1Ndt1();

                        // set bool for evaluation
                        mdNdtEval = false;
                    }
                    // return member data
                    return mdNdt;
                }
                // 2nd order time derivatives of the shape functions
                case 2:
                {
                    // if d2Ndt2 needs to be evaluated
                    if ( md2Ndt2Eval )
                    {
                        // evaluate d2Ndt2
                        this->eval_d2Ndt2();

                        // set bool for evaluation
                        md2Ndt2Eval = false;
                    }
                    // return member data
                    return md2Ndt2;
                }
                default:
                    MORIS_ERROR( false, "Field_Interpolator::dnNdtn - Derivative order not implemented." );
                    return mdNdt;
            }
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::eval_d1Ndt1()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel() > 0,
                    "Field_Interpolator::eval_d1Ndt1 - mXi  is not set." );
            MORIS_ASSERT( mTau.numel() > 0,
                    "Field_Interpolator::eval_d1Ndt1 - mTau is not set." );

            // evaluate dNTimedtau for the time interpolation
            Matrix< DDRMat > tdNTimedtau;
            mTimeInterpolation->eval_dNdXi( mTau, tdNTimedtau );

            // evaluate the Jacobian from the time geometry interpolator
            const Matrix< DDRMat >& tInvJGeot = mGeometryInterpolator->inverse_time_jacobian();

            // evaluate dNTimedt
            const Matrix< DDRMat > tdNTimedt = tInvJGeot * tdNTimedtau;

            // evaluate N for the field space interpolation
            Matrix< DDRMat > tNSpace;
            mSpaceInterpolation->eval_N( mXi, tNSpace );

            // build the space time dNdTau row by row
            for ( moris::uint Ik = 0; Ik < mNTimeBases; Ik++ )
            {
                mdNdt(
                        { 0, mNTimeDim - 1 },
                        { Ik * mNSpaceBases, ( Ik + 1 ) * mNSpaceBases - 1 } ) =
                        tNSpace * tdNTimedt( Ik );
            }
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::eval_d2Ndt2()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel() > 0,
                    "Field_Interpolator::eval_d2Ndt2 - mXi  is not set." );
            MORIS_ASSERT( mTau.numel() > 0,
                    "Field_Interpolator::eval_d2Ndt2 - mTau is not set." );

            // get matrices for second derivatives from space geometry interpolator
            Matrix< DDRMat > tJGeot, tKGeot, tLGeot;
            mGeometryInterpolator->time_jacobian_and_matrices_for_second_derivatives(
                    tJGeot,
                    tKGeot,
                    tLGeot,
                    mGeometryInterpolator->dNdTau(),
                    mGeometryInterpolator->d2NdTau2() );

            // get the derivatives of the space time SF wrt t
            Matrix< DDRMat > tdNFielddt = this->dnNdtn( 1 );

            // get space SF from the space interpolation
            Matrix< DDRMat > tNSpace;
            mSpaceInterpolation->eval_N( mXi, tNSpace );

            // get d2Ndtau2 for the time interpolation
            Matrix< DDRMat > td2NTimedtau2;
            mTimeInterpolation->eval_d2NdXi2( mTau, td2NTimedtau2 );

            // get the number of rows for td2NFielddtau2
            uint             tNSecondDerivatives = td2NTimedtau2.n_rows();
            Matrix< DDRMat > td2NFielddtau2( tNSecondDerivatives, mNFieldBases, 0.0 );

            // compute second derivatives of the space time SF wrt tau row by row
            for ( uint Ik = 0; Ik < mNTimeBases; Ik++ )
            {
                td2NFielddtau2(
                        { 0, mNTimeDim - 1 },
                        { Ik * mNSpaceBases, ( Ik + 1 ) * mNSpaceBases - 1 } ) +=
                        tNSpace * td2NTimedtau2( Ik );
            }

            // build the second derivatives of the space time SF wrt t
            td2NFielddtau2 -= tKGeot * tdNFielddt;
            md2Ndt2 = inv( tLGeot ) * td2NFielddtau2;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Field_Interpolator::d2Ndxt()
        {
            // if d2Ndxt needs to be evaluated
            if ( md2NdxtEval )
            {
                // evaluate d1Ndt1
                this->eval_d2Ndxt();

                // set bool for evaluation
                md2NdxtEval = false;
            }
            // return member data
            return md2Ndxt;
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::eval_d2Ndxt()
        {
            // check that mXi and mTau are set
            MORIS_ASSERT( mXi.numel() > 0,
                    "Field_Interpolator::eval_d2Ndxt - mXi  is not set." );
            MORIS_ASSERT( mTau.numel() > 0,
                    "Field_Interpolator::eval_d2Ndxt - mTau is not set." );

            // evaluate dNdTau for the field time interpolation
            Matrix< DDRMat > tdNTimedTau;
            mTimeInterpolation->eval_dNdXi( mTau, tdNTimedTau );

            // evaluate the time Jacobian from the geometry interpolator
            const Matrix< DDRMat >& tJGeoTimet = mGeometryInterpolator->time_jacobian();

            // compute first derivative of the space shape function wrt x
            Matrix< DDRMat > tdNTimedT = tdNTimedTau / tJGeoTimet( 0 );

            // evaluate dNSpacedXi for the field space interpolation
            Matrix< DDRMat > tdNSpacedXi;
            mSpaceInterpolation->eval_dNdXi( mXi, tdNSpacedXi );

            // evaluate the space Jacobian from the geometry interpolator
            const Matrix< DDRMat >& tInvJGeoSpacet = mGeometryInterpolator->inverse_space_jacobian();

            // compute first derivative of the space shape function wrt x
            Matrix< DDRMat > tdNSpacedX = tInvJGeoSpacet * tdNSpacedXi;

            // build the space time d2Ndxt
            for ( moris::uint Ik = 0; Ik < mNTimeBases; Ik++ )
            {
                md2Ndxt(
                        { 0, mNSpaceDim - 1 },
                        { Ik * mNSpaceBases, ( Ik + 1 ) * mNSpaceBases - 1 } ) =
                        tdNSpacedX * tdNTimedT( Ik );
            }
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Field_Interpolator::val()
        {
            // if field value needs to be evaluated
            if ( mValEval )
            {
                // evaluate field value
                this->eval_val();

                // set bool for evaluation
                mValEval = false;
            }

            // return member data
            return mVal;
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::eval_val()
        {
            // check that mUHat is set
            MORIS_ASSERT( mUHat.numel() > 0,
                    "Field_Interpolator::eval_val - mUHat is not set." );

            // evaluate the field value
            mVal = trans( this->NBuild() * mUHat );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Field_Interpolator::val_trans()
        {
            // if field value needs to be evaluated
            if ( mValTransEval )
            {
                // evaluate transpose of field value
                mValTrans = trans( this->val() );

                // set bool for evaluation
                mValTransEval = false;
            }

            // return member data
            return mValTrans;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Field_Interpolator::gradx( const uint& aDerivativeOrder )
        {
            switch ( aDerivativeOrder )
            {
                case 1:
                {
                    // if field value needs to be evaluated
                    if ( mGradx1Eval )
                    {
                        // evaluate gradient
                        this->eval_gradx( 1 );

                        // set bool for evaluation
                        mGradx1Eval = false;
                    }
                    return mGradx1;
                    break;
                }
                case 2:
                {
                    // if field value needs to be evaluated
                    if ( mGradx2Eval )
                    {
                        // evaluate gradient
                        this->eval_gradx( 2 );

                        // set bool for evaluation
                        mGradx2Eval = false;
                    }
                    return mGradx2;
                    break;
                }
                case 3:
                {
                    // if field value needs to be evaluated
                    if ( mGradx3Eval )
                    {
                        // evaluate gradient
                        this->eval_gradx( 3 );

                        // set bool for evaluation
                        mGradx3Eval = false;
                    }
                    return mGradx3;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Field_Interpolator:: gradx - Spatial derivative order not implemented." );
                    return mGradx1;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::eval_gradx( const uint& aDerivativeOrder )
        {
            // check that mUHat is set
            MORIS_ASSERT( mUHat.numel() > 0,
                    "Field_Interpolator::eval_gradx - mUHat is not set." );

            switch ( aDerivativeOrder )
            {
                case 1:
                {
                    mGradx1 = this->dnNdxn( 1 ) * mUHat;
                    return;
                    break;
                }
                case 2:
                {
                    mGradx2 = this->dnNdxn( 2 ) * mUHat;
                    return;
                    break;
                }
                case 3:
                {
                    mGradx3 = this->dnNdxn( 3 ) * mUHat;
                    return;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Field_Interpolator:: eval_gradx - Spatial derivative order not implemented." );
                }
            }
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Field_Interpolator::gradxt()
        {
            if ( mGradxtEval )
            {
                // evaluate gradient
                this->eval_gradxt();

                // set bool for evaluation
                mGradxtEval = false;
            }

            return mGradxt;
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::eval_gradxt()
        {
            // check that mUHat is set
            MORIS_ASSERT( mUHat.numel() > 0,
                    "Field_Interpolator::gradxt - mUHat is not set." );

            // evaluate and return gradient
            mGradxt = this->d2Ndxt() * mUHat;
        }

        //------------------------------------------------------------------------------

        moris::real
        Field_Interpolator::div()
        {
            // evaluate spatial divergence from gradient
            return sum( diag_vec( this->gradx( 1 ) ) );
        }

        //------------------------------------------------------------------------------

        const moris::Matrix< DDRMat >&
        Field_Interpolator::div_operator()
        {
            // if div_operator needs to be evaluated
            if ( mDivOperatorEval )
            {
                // evaluate d1Ndt1
                this->eval_div_operator();

                // set bool for evaluation
                mDivOperatorEval = false;
            }
            // return member data
            return mDivOperator;
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::eval_div_operator()
        {
            // evaluate spatial divergence operator from dNdx and flatten to row vector
            mDivOperator = trans( vectorize( trans( this->dnNdxn( 1 ) ) ) );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Field_Interpolator::gradt( const uint& aDerivativeOrder )
        {
            switch ( aDerivativeOrder )
            {
                case 1:
                {
                    // if field value needs to be evaluated
                    if ( mGradt1Eval )
                    {
                        // evaluate gradient
                        this->eval_gradt( 1 );

                        // set bool for evaluation
                        mGradt1Eval = false;
                    }
                    return mGradt1;
                    break;
                }
                case 2:
                {
                    // if field value needs to be evaluated
                    if ( mGradt2Eval )
                    {
                        // evaluate gradient
                        this->eval_gradt( 2 );

                        // set bool for evaluation
                        mGradt2Eval = false;
                    }
                    return mGradt2;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Field_Interpolator:: gradt - Temporal derivative order not implemented." );
                    return mGradt1;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        Field_Interpolator::eval_gradt( const uint& aDerivativeOrder )
        {
            // check that mUHat is set
            MORIS_ASSERT( mUHat.numel() > 0,
                    "Field_Interpolator::eval_gradt - mUHat is not set." );

            switch ( aDerivativeOrder )
            {
                case 1:
                {
                    mGradt1 = this->dnNdtn( 1 ) * mUHat;
                    return;
                    break;
                }
                case 2:
                {
                    mGradt2 = this->dnNdtn( 2 ) * mUHat;
                    return;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Field_Interpolator:: eval_gradt - Temporal derivative order not implemented." );
                }
            }
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

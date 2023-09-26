/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Geometry_Interpolator.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_
#define SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_

// MRS/COR/src
#include "typedefs.hpp"
// MTK/src
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Interpolation_Rule.hpp"
#include "cl_MTK_Space_Interpolator.hpp"
// FEM/INT/src
#include "cl_FEM_Enums.hpp"
// LINALG/src
#include "linalg_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "op_times.hpp"
#include "op_plus.hpp"
#include "op_minus.hpp"
#include "fn_trans.hpp"
#include "fn_det.hpp"
#include "fn_inv.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        /**
         * \brief a special interpolation class for geometry
         */
        class Geometry_Interpolator
        {
            // tolerance for check
            const real mEpsilon = 1e-12;

            // pointer to space interpolator object
            mtk::Space_Interpolator* mSpaceInterpolator = nullptr;

            // pointer to time interpolation function object
            mtk::Interpolation_Function_Base* mTimeInterpolation = nullptr;

            // number of space bases, number of physical and parametric dimensions
            uint mNumSpaceParamDim;

            // number of time bases and dimensions
            uint mNumTimeBases;
            uint mNumTimeDim;

            // matrix of space coefficients xHat
            // and matrix of time coefficients tHat
            Matrix< DDRMat > mTHat;

            // matrix of space param coefficients xiHat in the interpolation param space
            // and matrix of time param coefficients tauHat in the interpolation param space
            Matrix< DDRMat > mTauHat;

            // matrix of space param coords Xi in the local param space
            // and matrix of time param coords Tau in the local param space
            Matrix< DDRMat > mTauLocal;

            // element geometry type
            mtk::Geometry_Type mTimeGeometryType;

            // boolean true if side interpolation
            bool mSpaceSideset = false;

            // boolean true if time side interpolation
            bool mTimeSideset = false;

            // flag for evaluation
            bool mValtEval = true;

            bool mNTauEval     = true;
            bool mdNdTauEval   = true;
            bool md2NdTau2Eval = true;
            bool md3NdTau3Eval = true;

            bool mTimeDetJEval   = true;
            bool mTimeJacEval    = true;
            bool mInvTimeJacEval = true;

            // storage
            Matrix< DDRMat > mValt;

            Matrix< DDRMat > mNTau;
            Matrix< DDRMat > mdNdTau;
            Matrix< DDRMat > md2NdTau2;
            Matrix< DDRMat > md3NdTau3;

            Matrix< DDRMat > mTimeJac;
            Matrix< DDRMat > mInvTimeJac;
            real             mTimeDetJ;

            Matrix< DDRMat > mMappedPoint;

            // flag for mapping evaluation point
            bool mMapFlag = false;

            // pointer to function for time detJ
            real ( Geometry_Interpolator::*mTimeDetJFunc )(
                    const Matrix< DDRMat >& aTimeJt ) = nullptr;

            // point to function for inverse of time Jacobian
            void ( Geometry_Interpolator::*mInvTimeJacFunc )() = nullptr;

            // pointer to function for time second derivative
            void ( *mSecondDerivativeMatricesTime )(
                    const Matrix< DDRMat >& aJt,
                    Matrix< DDRMat >&       aKt,
                    Matrix< DDRMat >&       aLt,
                    const Matrix< DDRMat >& ad2NdTau2,
                    const Matrix< DDRMat >& aTHat ) = nullptr;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            // smallest acceptable value for DetJ
            // note: needs to allow for small negative numbers as used to identify degenerated elements
            static const real sDetJLowerLimit;

            // smallest acceptable value for DetJ used in building inverse of Jacobian
            // note: should be strictly positive as used for building inverse of matrix
            // note: this value should be used to identify and skip degenerated cells
            static const real sDetJInvJacLowerLimit;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------
            /**
             * trivial constructor
             */
            Geometry_Interpolator(){};

            /**
             * constructor
             * @param[ in ] interpolation rule for geometry
             * @param[ in ] flag true if side interpolation
             * @param[ in ] flag true if time side interpolation
             */
            Geometry_Interpolator(
                    const mtk::Interpolation_Rule& aInterpolationRule,
                    const mtk::CellShape&               aInterpolationShape = mtk::CellShape::GENERAL,
                    const bool                     aSpaceSideset       = false,
                    const bool                     aTimeSideset        = false );

            /**
             * constructor
             * @param[ in ] interpolation rule for geometry
             * @param[ in ] interpolation rule for geometry mapping
             * @param[ in ] flag true if side interpolation
             * @param[ in ] flag true if time side interpolation
             */
            Geometry_Interpolator(
                    const mtk::Interpolation_Rule& aInterpolationRule,
                    const mtk::Interpolation_Rule& aIPMapInterpolationRule,
                    const mtk::CellShape&               aInterpolationShape = mtk::CellShape::GENERAL,
                    const bool                     aSpaceSideset       = false,
                    const bool                     aTimeSideset        = false );

            //------------------------------------------------------------------------------
            /**
             * destructor
             */
            ~Geometry_Interpolator();

            //------------------------------------------------------------------------------
            /**
             * reset evaluation flags when xiLocal, tauLocal coordinates of
             * evaluation point are changed
             */
            void reset_eval_flags();

            //------------------------------------------------------------------------------
            /**
             * reset evaluation flags when global x,t coordinates or
             * local xi, tau coordinates are changed
             */
            void reset_eval_flags_coordinates();

            //------------------------------------------------------------------------------
            /**
             * returns the order of the space interpolation
             */
            mtk::Interpolation_Order
            get_space_interpolation_order() const
            {
                return mSpaceInterpolator->get_space_interpolation_order();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the order of the time interpolation
             */
            mtk::Interpolation_Order
            get_time_interpolation_order() const
            {
                return mTimeInterpolation->get_interpolation_order();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the type of the space interpolation
             */
            mtk::Interpolation_Type
            get_space_interpolation_type() const
            {
                return mSpaceInterpolator->get_space_interpolation_type();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the type of the time interpolation
             */
            mtk::Interpolation_Type
            get_time_interpolation_type() const
            {
                return mTimeInterpolation->get_interpolation_type();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the number of space dimensions
             */
            uint
            get_number_of_space_dimensions() const
            {
                return mSpaceInterpolator->get_number_of_space_dimensions();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the number of time dimensions
             */
            uint
            get_number_of_time_dimensions() const
            {
                return mTimeInterpolation->get_number_of_dimensions();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the number of bases for the space function
             */
            uint
            get_number_of_space_bases() const
            {
                return mSpaceInterpolator->get_number_of_space_bases();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the number of bases for the space function
             */
            uint
            get_number_of_time_bases() const
            {
                return mTimeInterpolation->get_number_of_bases();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the space geometry type
             */
            mtk::Geometry_Type
            get_space_geometry_type()
            {
                return mSpaceInterpolator->get_space_geometry_type();
            }

            //------------------------------------------------------------------------------
            /**
             * set the space and time coefficients of the geometry field xHat, tHat
             * @param[ in ] space coefficients
             * @param[ in ] time coefficients
             */
            void set_coeff(
                    const Matrix< DDRMat >& aXHat,
                    const Matrix< DDRMat >& aTHat );

            //------------------------------------------------------------------------------
            /**
             * set the space coefficients of the geometry field xHat
             * @param[ in ] space coefficients
             */
            void set_space_coeff( const Matrix< DDRMat >& aXHat );

            //------------------------------------------------------------------------------
            /**
             * set the time coefficients of the geometry field tHat
             * @param[ in ] time coefficients
             */
            void set_time_coeff( const Matrix< DDRMat >& aTHat );

            //------------------------------------------------------------------------------
            /**
             * get the space coefficients of the geometry field xHat
             */
            const Matrix< DDRMat >&
            get_space_coeff() const
            {
                return mSpaceInterpolator->get_space_coeff();
            }

            //------------------------------------------------------------------------------
            /**
             * get the time coefficients of the geometry field tHat
             */
            const Matrix< DDRMat >&
            get_time_coeff() const
            {
                // check that mTHat is set
                MORIS_ASSERT( mTHat.numel() > 0,
                        "Geometry_Interpolator::get_time_coeff - mTHat is not set." );

                return mTHat;
            }

            //------------------------------------------------------------------------------
            /**
             * get the time step delta t
             */
            real get_time_step();

            //------------------------------------------------------------------------------
            /**
             * set the space param coefficients of the geometry field xiHat
             * @param[ in ] space coefficients
             */
            // default implementation
            void set_param_coeff();

            // set implementation
            void set_param_coeff(
                    const Matrix< DDRMat >& aXiHat,
                    const Matrix< DDRMat >& aTauHat );

            //------------------------------------------------------------------------------
            /**
             * set the space param coefficients of the geometry field xiHat
             * @param[ in ] space coefficients
             */
            void set_space_param_coeff( const Matrix< DDRMat >& aXiHat );

            //------------------------------------------------------------------------------
            /**
             * set the time param coefficients of the geometry field tHat
             * @param[ in ] time coefficients
             */
            void set_time_param_coeff( const Matrix< DDRMat >& aTauHat );

            //------------------------------------------------------------------------------
            /**
             * get the space parametric coefficients of the geometry field xiHat
             */
            const Matrix< DDRMat >&
            get_space_param_coeff() const
            {
                // call space interpolator
                return mSpaceInterpolator->get_space_param_coeff();
            }

            //------------------------------------------------------------------------------
            /**
             * get the time parametric coefficients of the geometry field tauHat
             */
            const Matrix< DDRMat >&
            get_time_param_coeff() const
            {
                // check that mTauHat is set
                MORIS_ASSERT( mTauHat.numel() > 0,
                        "Geometry_Interpolator::get_time_param_coeff - mTauHat is not set." );

                return mTauHat;
            }

            //------------------------------------------------------------------------------
            /**
             * set the parametric point where geometry is interpolated
             * @param[ in ] aParamPoint evaluation point in space and time
             */
            void set_space_time( const Matrix< DDRMat >& aParamPoint );
            void set_space( const Matrix< DDRMat >& aSpaceParamPoint );
            void set_time( const Matrix< DDRMat >& aTimeParamPoint );

            void
            get_space_time( Matrix< DDRMat >& aParamPoint )
            {
                aParamPoint.set_size( mNumSpaceParamDim + 1, 1 );

                // setting space portion
                mSpaceInterpolator->get_space_time( aParamPoint );

                // setting time portion
                aParamPoint( mNumSpaceParamDim ) = mTauLocal( 0 );
            }

            //------------------------------------------------------------------------------
            /**
             * gets the space shape functions at a given evaluation point
             * @param[ out ] aNXi shape functions ( 1 x <number of nodes> )
             */
            const Matrix< DDRMat >& NXi();

            //------------------------------------------------------------------------------
            /**
             * gets the time shape functions at a given evaluation point
             * @param[ out ] aNTau shape functions ( 1 x <number of nodes> )
             */
            const Matrix< DDRMat >& NTau();

            /**
             * evaluates the time shape functions at a given evaluation point
             */
            void eval_NTau();

            //------------------------------------------------------------------------------
            /**
             * gets the first derivatives of the space shape functions
             * wrt parametric coordinates at a given evaluation point
             * @param[ out ] adNdXi derivatives ( <number of dimensions> x <number of nodes> )
             */
            const Matrix< DDRMat >& dNdXi();

            //------------------------------------------------------------------------------
            /**
             * gets the first derivatives of the time shape functions
             * wrt parametric coordinates at a given evaluation point
             * @param[ out ] adNdTau first order derivatives ( <number of dimensions> x <number of nodes> )
             */
            const Matrix< DDRMat >& dNdTau();

            /**
             * evaluates the first derivatives of the time shape functions
             * wrt parametric coordinates at a given evaluation point
             */
            void eval_dNdTau();

            //------------------------------------------------------------------------------
            /**
             * gets the second derivatives of the space shape functions
             * wrt parametric coordinates at a given evaluation point
             * @param[ out ] ad2NdXi2 second order derivatives ( <1D:1, 2D:3, 3D:6> x <number of nodes> )
             */
            const Matrix< DDRMat >& d2NdXi2();

            //------------------------------------------------------------------------------
            /**
             * gets the third derivatives of the space shape functions
             * wrt parametric coordinates at a given evaluation point
             * @param[ in ] ad3NdXi3 third order derivatives ( <1D:1, 2D:4, 3D: 10> x <number of nodes> )
             */
            const Matrix< DDRMat >& d3NdXi3();

            //------------------------------------------------------------------------------
            /**
             * gets the second derivatives of the time shape functions
             * wrt parametric coordinates at a given evaluation point
             * @param[ in ] ad2NdTau2 second order derivatives ( <number of dimensions> x <number of nodes> )
             */
            const Matrix< DDRMat >& d2NdTau2();

            /**
             * evaluates the second derivatives of the time shape functions
             * wrt parametric coordinates at a given evaluation point
             */
            void eval_d2NdTau2();

            //------------------------------------------------------------------------------
            /**
             * get the geometry Jacobian in space
             * @param[ out ] mSpaceJac transposed of geometry Jacobian in space
             */
            const Matrix< DDRMat >& space_jacobian();

            //------------------------------------------------------------------------------

            /**
             * get the inverse of the geometry Jacobian in space
             * @param[ out ] mInvSpaceJac inverse of the transposed of geometry Jacobian in space
             */
            const Matrix< DDRMat >& inverse_space_jacobian();

            //------------------------------------------------------------------------------
            /**
             * evaluates the 2nd geometry Jacobian in space
             * @param[ in ] aJ2bt     2nd geometry Jacobian in space
             */
            void second_space_jacobian( Matrix< DDRMat >& aJ2bt );

            //------------------------------------------------------------------------------
            /**
             * evaluates the geometry Jacobian in space
             * @param[ in ] aJ3ct     3rd geometry Jacobian in space
             */
            void third_space_jacobian( Matrix< DDRMat >& aJ3ct );

            //------------------------------------------------------------------------------
            /**
             * get the geometry Jacobian in time
             * @param[ out ] mTimeJac transposed of geometry Jacobian in time
             */
            const Matrix< DDRMat >& time_jacobian();

            /**
             * evaluates the geometry Jacobian in time
             */
            void eval_time_jacobian();

            //------------------------------------------------------------------------------
            /**
             * get the inverse of the geometry Jacobian in time
             * @param[ out ] mInvTimeJac inverse of the transposed of geometry Jacobian in time
             */
            const Matrix< DDRMat >& inverse_time_jacobian();

            /**
             * evaluates the inverse of the geometry Jacobian in time
             */
            void eval_inverse_time_jacobian();

            //------------------------------------------------------------------------------
            /**
             * evaluates the determinant of the Jacobian mapping
             * at given space and time evaluation point
             */
            real        det_J();
            const real& space_det_J();
            const real& time_det_J();

            //------------------------------------------------------------------------------
            /**
             * evaluates the normal to the side
             * in the case of a space side interpolation
             * at given space and time evaluation point
             * @param[ in ]  aNormal normal to be filled
             */
            void get_normal( Matrix< DDRMat >& aNormal );

            //------------------------------------------------------------------------------
            /**
             * evaluates the geometry Jacobian and the matrices needed for the second
             * derivatives wrt to space in physical space
             * @param[ in ]  aJt      transposed of geometry Jacobian
             * @param[ out ] aKt      matrix for second derivatives in physical space
             * @param[ out ] aLt      matrix for second derivatives in physical space
             * @param[ in ]  adNdXi   first derivatives of N in parameter space
             * @param[ in ]  ad2NdXi2 second derivatives of N in parameter space
             *
             */
            void space_jacobian_and_matrices_for_second_derivatives(
                    Matrix< DDRMat >&       aJt,
                    Matrix< DDRMat >&       aKt,
                    Matrix< DDRMat >&       aLt,
                    const Matrix< DDRMat >& adNdXi,
                    const Matrix< DDRMat >& ad2NdXi2 );

            //------------------------------------------------------------------------------
            /**
             * evaluates the geometry Jacobian and the matrices needed for the second
             * derivatives wrt to time in physical space
             * @param[ in ]  aJt       transposed of geometry Jacobian
             * @param[ out ] aKt       matrix for second derivatives in physical space
             * @param[ out ] aLt       matrix for second derivatives in physical space
             * @param[ in ]  adNdTau   first derivatives of N in parameter space
             * @param[ in ]  ad2NdTau2 second derivatives of N in parameter space
             *
             */
            void time_jacobian_and_matrices_for_second_derivatives(
                    Matrix< DDRMat >&       aJt,
                    Matrix< DDRMat >&       aKt,
                    Matrix< DDRMat >&       aLt,
                    const Matrix< DDRMat >& adNdTau,
                    const Matrix< DDRMat >& ad2NdTau2 );

            //------------------------------------------------------------------------------
            /**
             * evaluates the geometry Jacobian and the matrices needed for the second
             * derivatives wrt to space in physical space
             * @param[ in ]  aJ1t       transposed of 1st geometry Jacobian
             * @param[ in ]  aJ2bt      2nd geometry Jacobian = 2nd help matrix for 2nd derivs
             *
             * @param[ out ] aJ3at      first help matrix for 3rd field derivs
             * @param[ out ] aJ3bt      second help matrix for 3rd field derivs
             * @param[ out ] aJ3ct      third help matrix for 3rd field derivs
             *
             * @param[ in ]  adNdXi   first derivatives of N in parameter space
             * @param[ in ]  ad2NdXi2 second derivatives of N in parameter space
             * @param[ in ]  ad3NdXi3 third derivatives of N in parameter space
             *
             */
            void space_jacobian_and_matrices_for_third_derivatives(
                    Matrix< DDRMat >&       aJt,
                    Matrix< DDRMat >&       aJ2bt,
                    Matrix< DDRMat >&       aJ3at,
                    Matrix< DDRMat >&       aJ3bt,
                    Matrix< DDRMat >&       aJ3ct,
                    const Matrix< DDRMat >& adNdXi,
                    const Matrix< DDRMat >& ad2NdXi2,
                    const Matrix< DDRMat >& ad3NdXi3 );

            //------------------------------------------------------------------------------
            /**
             * evaluates the space geometry field at a given evaluation point in space
             * @param[ out ] aX   location in space
             */
            const Matrix< DDRMat >& valx();

            //------------------------------------------------------------------------------
            /**
             * evaluates the space geometry field at a given evaluation point in time
             * @param[ out ] aT   location in time
             */
            const Matrix< DDRMat >& valt();

            //------------------------------------------------------------------------------
            /**
             * map an integration point from local param coords to global param coords
             * @param[ in ] aGlobalParamPoint param coords in global parametric space
             */
            const Matrix< DDRMat >& map_integration_point();

            //------------------------------------------------------------------------------
            /**
             * update parametric coordinates (xi, eta, zeta)
             * for given physical coordinates (x, y, z)
             * @param[ in ] aPhysCoordinates  coords in physical space
             * @param[ in ] aParamCoordinates coords in parametric space
             */
            void update_local_coordinates(
                    Matrix< DDRMat >& aPhysCoordinates,
                    Matrix< DDRMat >& aParamCoordinates );

            //------------------------------------------------------------------------------
            /**
             * get metric tensor
             * @param[ out ] aG   metric tensor
             * Gij = sum_d dxi_d/dx_i dxi_d/dx_j, d = 1, ..., nSpaceDim
             */
            const Matrix< DDRMat >& metric_tensor();

            //------------------------------------------------------------------------------

          private:
            //------------------------------------------------------------------------------
            /**
             * sets the function pointers for 2D and 3D. Called during construction.
             */
            void set_function_pointers();

            //------------------------------------------------------------------------------
            /**
             * evaluate time detJ
             */
            real eval_time_detJ_side( const Matrix< DDRMat >& aTimeJt );
            real eval_time_detJ_bulk( const Matrix< DDRMat >& aTimeJt );

            //------------------------------------------------------------------------------
            /**
             * evaluate time Jacobians
             */
            void eval_inverse_time_jacobian_1d();
            void eval_inverse_time_jacobian_2d();
            void eval_inverse_time_jacobian_3d();
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_ */

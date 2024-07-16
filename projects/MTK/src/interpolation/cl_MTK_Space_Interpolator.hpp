/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Space_Interpolator.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_SPACE_INTERPOLATOR_HPP_
#define SRC_MTK_CL_MTK_SPACE_INTERPOLATOR_HPP_

// MRS/COR/src
#include "moris_typedefs.hpp"
// MTK/src
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Interpolation_Rule.hpp"
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
    namespace mtk
    {
        //------------------------------------------------------------------------------
        /**
         * \brief a special interpolation class for Space
         */
        class Space_Interpolator
        {
            // tolerance for check
            const real mEpsilon = 1e-12;

            // pointer to space interpolation function object
            Interpolation_Function_Base* mSpaceInterpolation = nullptr;

            // number of space bases, number of physical and parametric dimensions
            uint mNumSpaceBases;
            uint mNumSpaceDim;
            uint mNumSpaceParamDim;

            // matrix of space coefficients xHat
            Matrix< DDRMat > mXHat;

            // matrix of space param coefficients xiHat in the interpolation param space
            Matrix< DDRMat > mXiHat;

            // matrix of space param coords Xi in the local param space
            Matrix< DDRMat > mXiLocal;

            // element geometry type
            Geometry_Type mGeometryType;

            // element shape
            CellShape mInterpolationShape;

            // interpolation cell geometry type
            Geometry_Type mIPMappingGeometryType;
            uint          mIPMappingNumSpaceParamDim;

            // boolean true if side interpolation
            bool mSpaceSideset = false;

            // flag for evaluation
            bool mValxEval = true;

            bool mNXiEval     = true;
            bool mdNdXiEval   = true;
            bool md2NdXi2Eval = true;
            bool md3NdXi3Eval = true;

            bool mSpaceDetJEval      = true;
            bool mSpaceJacEval       = true;
            bool mInvSpaceJacEval    = true;
            bool mSpaceJacDerivEval  = true;
            bool mSpaceDetJDerivEval = true;

            bool mMetricTensorEval = true;

            // storage
            Matrix< DDRMat > mValx;

            Matrix< DDRMat > mNXi;
            Matrix< DDRMat > mdNdXi;
            Matrix< DDRMat > md2NdXi2;
            Matrix< DDRMat > md3NdXi3;

            Matrix< DDRMat > mSpaceJac;
            Matrix< DDRMat > mInvSpaceJac;
            real             mSpaceDetJ;
            Matrix< DDRMat > mSpaceJacDeriv;
            real             mSpaceDetJDeriv;

            Matrix< DDRMat > mMetricTensor;

            Matrix< DDRMat > mMappedPoint;

            // flag for mapping evaluation point
            bool mMapFlag = false;

            // pointer to function for space detJ
            real ( Space_Interpolator::*mSpaceDetJFunc )(
                    const Matrix< DDRMat >& aSpaceJt ) = nullptr;

            // pointers for derivatives of space detJ wrt a single dof
            real ( Space_Interpolator::*mSpaceDetJDerivFunc )(
                    const Matrix< DDRMat >& aSpaceJt ) = nullptr;

            // point to function for inverse of space Jacobian
            void ( Space_Interpolator::*mInvSpaceJacFunc )() = nullptr;

            // pointer to function for normal
            void ( Space_Interpolator::*mNormalFunc )(
                    Matrix< DDRMat >& aTangent,
                    Matrix< DDRMat >& aNormal ) = nullptr;

            // pointer to function for space second derivative
            void ( *mSecondDerivativeMatricesSpace )(
                    const Matrix< DDRMat >& aJt,
                    Matrix< DDRMat >&       aKt,
                    Matrix< DDRMat >&       aLt,
                    const Matrix< DDRMat >& ad2NdXi2,
                    const Matrix< DDRMat >& aXHat );

            // pointer to function for space third derivatives
            void ( *mThirdDerivativeMatricesSpace )(
                    const Matrix< DDRMat >& aJt,
                    const Matrix< DDRMat >& aJ2bt,
                    Matrix< DDRMat >&       aJ3at,
                    Matrix< DDRMat >&       aJ3bt,
                    Matrix< DDRMat >&       aJ3ct,
                    const Matrix< DDRMat >& ad3NdXi3,
                    const Matrix< DDRMat >& aXHat );

            // pointer to function for metric tensor
            void ( Space_Interpolator::*mMetricTensorFunc )(
                    const Matrix< DDRMat >& aInvSpaceJacobian ) = nullptr;

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
            Space_Interpolator(){};

            /**
             * constructor
             * @param[ in ] interpolation rule for geometry
             * @param[ in ] flag true if side interpolation
             */
            Space_Interpolator(
                    const Interpolation_Rule& aInterpolationRule,
                    const CellShape&          aInterpolationShape = CellShape::GENERAL,
                    const bool                aSpaceSideset       = false );

            /**
             * constructor
             * @param[ in ] interpolation rule for geometry
             * @param[ in ] interpolation rule for geometry mapping
             * @param[ in ] flag true if side interpolation
             */
            Space_Interpolator(
                    const Interpolation_Rule& aInterpolationRule,
                    const Interpolation_Rule& aIPMapInterpolationRule,
                    const CellShape&          aInterpolationShape = CellShape::GENERAL,
                    const bool                aSpaceSideset       = false );

            //------------------------------------------------------------------------------
            /**
             * destructor
             */
            ~Space_Interpolator();

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
             * reset evaluation flags when needing to get derivative of same detJ wrt to a
             * different dof
             */
            void reset_eval_flags_deriv();

            //------------------------------------------------------------------------------
            /**
             * returns the order of the space interpolation
             */
            Interpolation_Order
            get_space_interpolation_order() const
            {
                return mSpaceInterpolation->get_interpolation_order();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the type of the space interpolation
             */
            Interpolation_Type
            get_space_interpolation_type() const
            {
                return mSpaceInterpolation->get_interpolation_type();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the number of space dimensions
             */
            uint
            get_number_of_space_dimensions() const
            {
                return mSpaceInterpolation->get_number_of_dimensions();
            }

            //------------------------------------------------------------------------------

            /**
             * returns the number of parametric dimensions
             */
            uint
            get_number_of_param_dimensions() const
            {
                return mSpaceInterpolation->get_number_of_param_dimensions();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the number of bases for the space function
             */
            uint
            get_number_of_space_bases() const
            {
                return mSpaceInterpolation->get_number_of_bases();
            }

            //------------------------------------------------------------------------------
            /**
             * returns the space geometry type
             */
            Geometry_Type
            get_space_geometry_type()
            {
                return mGeometryType;
            }

            //------------------------------------------------------------------------------
            /**
             * returns the boolean map flag
             */
            bool
            get_map_flag()
            {
                return mMapFlag;
            }

            //------------------------------------------------------------------------------
            /**
             * returns the initialed mapped point
             * will not necessarily be filled with a mapping yet
             */
            Matrix< DDRMat >
            get_initialized_mapped_point()
            {
                return mMappedPoint;
            }

            //------------------------------------------------------------------------------
            /**
             * set the space coefficients of the geometry field xHat
             * @param[ in ] space coefficients
             */
            void set_space_coeff( const Matrix< DDRMat >& aXHat );

            //------------------------------------------------------------------------------
            /**
             * get the space coefficients of the geometry field xHat
             */
            const Matrix< DDRMat >&
            get_space_coeff() const
            {
                // check that mXHat is set
                MORIS_ASSERT( mXHat.numel() > 0,
                        "Space_Interpolator::get_space_coeff - mXHat is not set." );

                return mXHat;
            }

            //------------------------------------------------------------------------------
            /**
             * set the space param coefficients of the geometry field xiHat
             * @param[ in ] space coefficients
             */
            // default implementation
            void set_param_coeff();

            //------------------------------------------------------------------------------
            /**
             * set the space param coefficients of the geometry field xiHat
             * @param[ in ] space coefficients
             */
            void set_space_param_coeff( const Matrix< DDRMat >& aXiHat );

            //------------------------------------------------------------------------------
            /**
             * get the space parametric coefficients of the geometry field xiHat
             */
            const Matrix< DDRMat >&
            get_space_param_coeff() const
            {
                // check that mXiHat is set
                MORIS_ASSERT( mXiHat.numel() > 0,
                        "Space_Interpolator::get_space_param_coeff - mXiHat is not set." );

                return mXiHat;
            }

            //------------------------------------------------------------------------------
            /**
             * set the parametric point where geometry is interpolated
             * @param[ in ] aParamPoint evaluation point in space and time
             */
            void set_space_time( const Matrix< DDRMat >& aParamPoint );
            void set_space( const Matrix< DDRMat >& aSpaceParamPoint );

            void
            get_space_time( Matrix< DDRMat >& aParamPoint )
            {
                aParamPoint.set_size( mNumSpaceParamDim + 1, 1 );

                aParamPoint( { 0, mNumSpaceParamDim - 1 } ) = mXiLocal.matrix_data();
            }

            //------------------------------------------------------------------------------
            /**
             * gets the space shape functions at a given evaluation point
             * @param[ out ] aNXi shape functions ( 1 x <number of nodes> )
             */
            const Matrix< DDRMat >& NXi();

            /**
             * evaluates the space shape functions at a given evaluation point
             */
            void eval_NXi();

            //------------------------------------------------------------------------------
            /**
             * gets the first derivatives of the space shape functions
             * wrt parametric coordinates at a given evaluation point
             * @param[ out ] adNdXi derivatives ( <number of dimensions> x <number of nodes> )
             */
            const Matrix< DDRMat >& dNdXi();

            /**
             * evaluates the first derivatives of the space shape functions
             * wrt parametric coordinates at a given evaluation point
             */
            void eval_dNdXi();

            //------------------------------------------------------------------------------
            /**
             * gets the second derivatives of the space shape functions
             * wrt parametric coordinates at a given evaluation point
             * @param[ out ] ad2NdXi2 second order derivatives ( <1D:1, 2D:3, 3D:6> x <number of nodes> )
             */
            const Matrix< DDRMat >& d2NdXi2();

            /**
             * evaluates the second derivatives of the space shape functions
             * wrt parametric coordinates at a given evaluation point
             */
            void eval_d2NdXi2();

            //------------------------------------------------------------------------------
            /**
             * gets the third derivatives of the space shape functions
             * wrt parametric coordinates at a given evaluation point
             * @param[ in ] ad3NdXi3 third order derivatives ( <1D:1, 2D:4, 3D: 10> x <number of nodes> )
             */
            const Matrix< DDRMat >& d3NdXi3();

            /**
             * evaluates the third derivatives of the space shape functions
             * wrt parametric coordinates at a given evaluation point
             */
            void eval_d3NdXi3();

            //------------------------------------------------------------------------------
            /**
             * get the geometry Jacobian in space
             * @param[ out ] mSpaceJac transposed of geometry Jacobian in space
             */
            const Matrix< DDRMat >& space_jacobian();

            /**
             * evaluates the geometry Jacobian in space
             */
            void eval_space_jacobian();

            //------------------------------------------------------------------------------
            /**
             * get the transformation matrix that is used to determine d(detJ)/dxHat_i
             * this is not a true derivative of the space jacobian
             * @param[ out ] aLocalVertexID local vertex to take derivative wrt
             * @param[ out ] aDirection     direction to take derivative wrt (0, 1, or 2)
             */
            const Matrix< DDRMat >& space_jacobian_deriv( const uint& aLocalVertexID, const uint& aDirection );

            /**
             * evaluates the transformation matrix used to determine d(detJ)/dxHat_i.
             * In a 2D situation, with a derivative wrt to y_0, this will produce a matrix
             * that looks like the following
             * { { d[N]/dXi .* [x],  dN_0/dXi },
             *   { d[N]/dEta.* [x],  dN_0/dEta}}
             * @param[ out ] aLocalVertexID local vertex to take derivative wrt
             * @param[ out ] aDirection     direction to take derivative wrt (0, 1, or 2)
             */
            void eval_space_jacobian_deriv( const uint& aLocalVertexID, const uint& aDirection );

            //------------------------------------------------------------------------------

            /**
             * get the inverse of the geometry Jacobian in space
             * @param[ out ] mInvSpaceJac inverse of the transposed of geometry Jacobian in space
             */
            const Matrix< DDRMat >& inverse_space_jacobian();

            /**
             * evaluates the inverse of the geometry Jacobian in space
             */
            void eval_inverse_space_jacobian();

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
             * evaluates the determinant of the Jacobian mapping
             * at given space and time evaluation point
             */
            const real& space_det_J();

            //------------------------------------------------------------------------------
            /**
             * evaluates the determinant of the Jacobian Derivative wrt a single dof mapping
             * at given space and time evaluation point
             */
            const real& space_det_J_deriv( const uint& aLocalVertexID, const uint& aDirection );

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
             * evaluates matrices that are needed for the second derivative
             * in space, 2D version. It is
             *
             * \f[
             *      \mathbf{L}^T \, \mathbf{\frac{\partial^2 N}{\partial x^2}}
             *      = \mathbf{\frac{\partial^2 N}{\partial \xi^2}}
             *      - K^T \, mathbf{\frac{\partial N}{\partial x}}
             * \f]
             *
             * @param[ in ]  aJt          transposed of geometry Jacobian
             * @param[ out ] aKt          transposed help matrix K
             * @param[ out ] aLt          transposed help matrix L
             * @param[ in ]  adNdXi       first derivatives in parameter space
             * @param[ in ]  ad2NdX2i     second derivatives in parameter space
             *
             */
            static void eval_matrices_for_second_derivative_1d(
                    const Matrix< DDRMat >& aJt,
                    Matrix< DDRMat >&       aKt,
                    Matrix< DDRMat >&       aLt,
                    const Matrix< DDRMat >& ad2NdXi2,
                    const Matrix< DDRMat >& aXHat );

            //------------------------------------------------------------------------------
            /**
             * evaluates matrices that are needed for the second derivative
             * in space, 2D version. It is
             *
             * \f[
             *      \mathbf{L}^T \, \mathbf{\frac{\partian^2 N}{\partial x^2}}
             *      = \mathbf{\frac{\partian^2 N}{\partial \xi^2}}
             *      - K^T \, mathbf{\frac{\partian N}{\partial x}}
             * \f]
             *
             * @param[ in ]  aJt          transposed of geometry Jacobian
             * @param[ out ] aKt          transposed help matrix K
             * @param[ out ] aLt          transposed help matrix L
             * @param[ in ]  adNdXi       first derivatives in parameter space
             * @param[ in ]  ad2NdX2i     second derivatives in parameter space
             *
             */
            static void eval_matrices_for_second_derivative_2d(
                    const Matrix< DDRMat >& aJt,
                    Matrix< DDRMat >&       aKt,
                    Matrix< DDRMat >&       aLt,
                    const Matrix< DDRMat >& ad2NdXi2,
                    const Matrix< DDRMat >& aXHat );

            //------------------------------------------------------------------------------
            /**
             * evaluates matrices that are needed for the second derivative
             * in space, 3D version. It is
             *
             * \f[
             *      \mathbf{L}^T \, \mathbf{\frac{\partian^2 N}{\partial x^2}}
             *      = \mathbf{\frac{\partian^2 N}{\partial \xi^2}}
             *      - K^T \, mathbf{\frac{\partian N}{\partial x}}
             * \f]
             *
             * @param[ in ]  aJt          transposed of geometry Jacobian
             * @param[ out ] aKt          transposed help matrix K
             * @param[ out ] aLt          transposed help matrix L
             * @param[ in ]  adNdXi       first derivatives in parameter space
             * @param[ in ]  ad2NdX2i     second derivatives in parameter space
             *
             */
            static void eval_matrices_for_second_derivative_3d(
                    const Matrix< DDRMat >& aJt,
                    Matrix< DDRMat >&       aKt,
                    Matrix< DDRMat >&       aLt,
                    const Matrix< DDRMat >& ad2NdXi2,
                    const Matrix< DDRMat >& aXHat );

            //------------------------------------------------------------------------------
            /**
             * get the metric tensor at a given evaluation point in space
             * where Gij = sum_d dxi_d/dx_i dxi_d/dx_j, d = 1, ..., nSpaceDim
             * @param[ out ] mMetricTensor   metric tensor
             */
            const Matrix< DDRMat >& metric_tensor();

            /**
             * get the metric tensor at a given evaluation point in space
             * @param[ in ] aInvSpaceJacobian  inverse of space jacobian
             */
            void eval_metric_tensor_1d(
                    const Matrix< DDRMat >& aInvSpaceJacobian );

            void eval_metric_tensor_2d(
                    const Matrix< DDRMat >& aInvSpaceJacobian );

            void eval_metric_tensor_3d(
                    const Matrix< DDRMat >& aInvSpaceJacobian );

            //------------------------------------------------------------------------------

          private:
            //------------------------------------------------------------------------------
            /**
             * sets the function pointers for 2D and 3D. Called during construction.
             */
            void set_function_pointers();

            //------------------------------------------------------------------------------
            /**
             * evaluate space detJ.
             */
            real eval_space_detJ_side_line( const Matrix< DDRMat >& aSpaceJt );
            real eval_space_detJ_side_tri( const Matrix< DDRMat >& aSpaceJt );
            real eval_space_detJ_side_quad( const Matrix< DDRMat >& aSpaceJt );

            real eval_space_detJ_bulk_line( const Matrix< DDRMat >& aSpaceJt );
            real eval_space_detJ_bulk_quad( const Matrix< DDRMat >& aSpaceJt );
            real eval_space_detJ_bulk_quad_rect( const Matrix< DDRMat >& aSpaceJt );
            real eval_space_detJ_bulk_hex( const Matrix< DDRMat >& aSpaceJt );
            real eval_space_detJ_bulk_hex_rect( const Matrix< DDRMat >& aSpaceJt );
            real eval_space_detJ_bulk_tri_param_2( const Matrix< DDRMat >& aSpaceJt );
            real eval_space_detJ_bulk_tri_param_3( const Matrix< DDRMat >& aSpaceJt );
            real eval_space_detJ_bulk_tet_param_3( const Matrix< DDRMat >& aSpaceJt );
            real eval_space_detJ_bulk_tet_param_4( const Matrix< DDRMat >& aSpaceJt );

            //------------------------------------------------------------------------------
            /**
             * evaluate derivative of the space detJ wrt to a single dof.
             * These functions are very similar to the standard detJ calcs with different asserts
             */
            real eval_space_detJ_deriv_side_line( const Matrix< DDRMat >& aSpaceJDerivt );
            real eval_space_detJ_deriv_side_tri( const Matrix< DDRMat >& aSpaceJDerivt );
            real eval_space_detJ_deriv_side_quad( const Matrix< DDRMat >& aSpaceJDerivt );

            real eval_space_detJ_deriv_bulk_line( const Matrix< DDRMat >& aSpaceJDerivt );
            real eval_space_detJ_deriv_bulk_quad( const Matrix< DDRMat >& aSpaceJDerivt );
            real eval_space_detJ_deriv_bulk_quad_rect( const Matrix< DDRMat >& aSpaceJDerivt );
            real eval_space_detJ_deriv_bulk_hex( const Matrix< DDRMat >& aSpaceJDerivt );
            real eval_space_detJ_deriv_bulk_hex_rect( const Matrix< DDRMat >& aSpaceJDerivt );
            real eval_space_detJ_deriv_bulk_tri_param_2( const Matrix< DDRMat >& aSpaceJDerivt );
            real eval_space_detJ_deriv_bulk_tri_param_3( const Matrix< DDRMat >& aSpaceJDerivt );
            real eval_space_detJ_deriv_bulk_tet_param_3( const Matrix< DDRMat >& aSpaceJDerivt );
            real eval_space_detJ_deriv_bulk_tet_param_4( const Matrix< DDRMat >& aSpaceJDerivt );

            //------------------------------------------------------------------------------
            /**
             * evaluate space Jacobians
             */
            void eval_inverse_space_jacobian_1d();
            void eval_inverse_space_jacobian_2d();
            void eval_inverse_space_jacobian_3d();
            void eval_inverse_space_jacobian_2d_rect();
            void eval_inverse_space_jacobian_3d_rect();
            void eval_inverse_space_jacobian_2d_tri();
            void eval_inverse_space_jacobian_3d_tri();

            //------------------------------------------------------------------------------
            /**
             * evaluate normal to side
             */
            void eval_normal_side_line(
                    Matrix< DDRMat >& aRealTangents,
                    Matrix< DDRMat >& aNormal );
            void eval_normal_side_quad(
                    Matrix< DDRMat >& aRealTangents,
                    Matrix< DDRMat >& aNormal );
            void eval_normal_side_tri(
                    Matrix< DDRMat >& aRealTangents,
                    Matrix< DDRMat >& aNormal );

            //------------------------------------------------------------------------------
            /**
             * evaluates the geometry Jacobian and the matrices needed for the second
             * derivatives wrt to space in 1D physical space
             * @param[ in ]  aJ1t       transposed of 1st geometry Jacobian
             * @param[ in ]  aJ2bt      2nd geometry Jacobian = 2nd help matrix for 2nd derivs
             *
             * @param[ out ] aJ3at      first help matrix for 3rd field derivs
             * @param[ out ] aJ3bt      second help matrix for 3rd field derivs
             * @param[ out ] aJ3ct      third help matrix for 3rd field derivs
             *
             * @param[ in ]  ad3NdXi3 third derivatives of N in parameter space
             *
             */
            static void eval_matrices_for_third_derivative_1d(
                    const Matrix< DDRMat >& aJt,
                    const Matrix< DDRMat >& aJ2bt,
                    Matrix< DDRMat >&       aJ3at,
                    Matrix< DDRMat >&       aJ3bt,
                    Matrix< DDRMat >&       aJ3ct,
                    const Matrix< DDRMat >& ad3NdXi3,
                    const Matrix< DDRMat >& aXHat );

            //------------------------------------------------------------------------------
            /**
             * evaluates the geometry Jacobian and the matrices needed for the second
             * derivatives wrt to space in 2D physical space
             * @param[ in ]  aJ1t       transposed of 1st geometry Jacobian
             * @param[ in ]  aJ2bt      2nd geometry Jacobian = 2nd help matrix for 2nd derivs
             *
             * @param[ out ] aJ3at      first help matrix for 3rd field derivs
             * @param[ out ] aJ3bt      second help matrix for 3rd field derivs
             * @param[ out ] aJ3ct      third help matrix for 3rd field derivs
             *
             * @param[ in ]  ad3NdXi3 third derivatives of N in parameter space
             *
             */
            static void eval_matrices_for_third_derivative_2d(
                    const Matrix< DDRMat >& aJt,
                    const Matrix< DDRMat >& aJ2bt,
                    Matrix< DDRMat >&       aJ3at,
                    Matrix< DDRMat >&       aJ3bt,
                    Matrix< DDRMat >&       aJ3ct,
                    const Matrix< DDRMat >& ad3NdXi3,
                    const Matrix< DDRMat >& aXHat );

            //------------------------------------------------------------------------------
            /**
             * evaluates the geometry Jacobian and the matrices needed for the second
             * derivatives wrt to space in 3D physical space
             * @param[ in ]  aJ1t       transposed of 1st geometry Jacobian
             * @param[ in ]  aJ2bt      2nd geometry Jacobian = 2nd help matrix for 2nd derivs
             *
             * @param[ out ] aJ3at      first help matrix for 3rd field derivs
             * @param[ out ] aJ3bt      second help matrix for 3rd field derivs
             * @param[ out ] aJ3ct      third help matrix for 3rd field derivs
             *
             * @param[ in ]  ad3NdXi3 third derivatives of N in parameter space
             *
             */
            static void eval_matrices_for_third_derivative_3d(
                    const Matrix< DDRMat >& aJt,
                    const Matrix< DDRMat >& aJ2bt,
                    Matrix< DDRMat >&       aJ3at,
                    Matrix< DDRMat >&       aJ3bt,
                    Matrix< DDRMat >&       aJ3ct,
                    const Matrix< DDRMat >& ad3NdXi3,
                    const Matrix< DDRMat >& aXHat );
        };

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* SRC_MTK_CL_MTK_SPACE_INTERPOLATOR_HPP_ */

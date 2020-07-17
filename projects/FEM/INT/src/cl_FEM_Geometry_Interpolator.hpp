/*
 * cl_FEM_Geometry_Interpolator.hpp
 *
 *  Created on: Jan 31, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_
#define SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_

//MRS/COR/src
#include "typedefs.hpp"
//MTK/src
#include "cl_MTK_Enums.hpp"
//FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Interpolation_Rule.hpp"
//LINALG/src
#include "linalg_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "op_times.hpp"
#include "op_plus.hpp"
#include "op_minus.hpp"
#include "fn_trans.hpp"
#include "fn_det.hpp"

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
                real mEpsilon = 1e-12;

                // pointer to space interpolation function object
                Interpolation_Function_Base * mSpaceInterpolation = nullptr;

                // pointer to time interpolation function object
                Interpolation_Function_Base * mTimeInterpolation = nullptr;

                // number of space bases, number of physical and parametric dimensions
                uint mNumSpaceBases;
                uint mNumSpaceDim;
                uint mNumSpaceParamDim;

                // number of time bases and dimensions
                uint mNumTimeBases;
                uint mNumTimeDim;

                // matrix of space coefficients xHat
                // and matrix of time coefficients tHat
                Matrix < DDRMat > mXHat;
                Matrix < DDRMat > mTHat;

                // matrix of space param coefficients xiHat in the interpolation param space
                // and matrix of time param coefficients tauHat in the interpolation param space
                Matrix < DDRMat > mXiHat;
                Matrix < DDRMat > mTauHat;

                // matrix of space param coords Xi in the local param space
                // and matrix of time param coords Tau in the local param space
                Matrix< DDRMat > mXiLocal;
                Matrix< DDRMat > mTauLocal;

                // element geometry type
                mtk::Geometry_Type mGeometryType;
                mtk::Geometry_Type mTimeGeometryType;

                // boolean true if side interpolation
                bool mSpaceSideset = false;

                // boolean true if time side interpolation
                bool mTimeSideset = false;

                // flag for evaluation
                bool mNXiEval     = true;
                bool mdNdXiEval   = true;
                bool md2NdXi2Eval = true;
                bool md3NdXi3Eval = true;

                bool mNTauEval     = true;
                bool mdNdTauEval   = true;
                bool md2NdTau2Eval = true;
                bool md3NdTau3Eval = true;

                // storage
                Matrix< DDRMat > mNXi;
                Matrix< DDRMat > mdNdXi;
                Matrix< DDRMat > md2NdXi2;
                Matrix< DDRMat > md3NdXi3;

                Matrix< DDRMat > mNTau;
                Matrix< DDRMat > mdNdTau;
                Matrix< DDRMat > md2NdTau2;
                Matrix< DDRMat > md3NdTau3;

                // pointer to function for space detJ
                real ( Geometry_Interpolator:: * mSpaceDetJFunc )(
                        Matrix< DDRMat > & aSpaceJt ) = nullptr;

                // pointer to function for time detJ
                real (  Geometry_Interpolator:: * mTimeDetJFunc )(
                        Matrix< DDRMat > & aTimeJt ) = nullptr;

                // pointer to function for normal
                void (  Geometry_Interpolator:: * mNormalFunc )(
                        Matrix< DDRMat > & aTangent,
                        Matrix< DDRMat > & aNormal ) = nullptr;

                // pointer to function for space second derivative
                void ( * mSecondDerivativeMatricesSpace )(
                        const Matrix< DDRMat > & aJt,
                        Matrix< DDRMat > & aKt,
                        Matrix< DDRMat > & aLt,
                        const Matrix< DDRMat > & ad2NdXi2,
                        const Matrix< DDRMat > & aXHat );

                // pointer to function for space third derivatives
                void ( * mThirdDerivativeMatricesSpace )(
                        const Matrix< DDRMat > & aJt,
                        const Matrix< DDRMat > & aJ2bt,
                        Matrix< DDRMat > & aJ3at,
                        Matrix< DDRMat > & aJ3bt,
                        Matrix< DDRMat > & aJ3ct,
                        const Matrix< DDRMat > & ad3NdXi3,
                        const Matrix< DDRMat > & aXHat );

                // pointer to function for time second derivative
                void ( * mSecondDerivativeMatricesTime )(
                        const Matrix< DDRMat > & aJt,
                        Matrix< DDRMat > & aKt,
                        Matrix< DDRMat > & aLt,
                        const Matrix< DDRMat > & ad2NdTau2,
                        const Matrix< DDRMat > & aTHat );

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
                        const Interpolation_Rule & aInterpolationRule,
                        const bool                 aSpaceSideset = false,
                        const bool                 aTimeSideset  = false );

                //------------------------------------------------------------------------------
                /**
                 * destructor
                 */
                ~Geometry_Interpolator();

                //------------------------------------------------------------------------------
                /**
                 * reset evaluation flags
                 */
                void reset_eval_flags();

                //------------------------------------------------------------------------------
                /**
                 * returns the order of the space interpolation
                 */
                mtk::Interpolation_Order get_space_interpolation_order() const
                {
                    return mSpaceInterpolation->get_interpolation_order();
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the order of the time interpolation
                 */
                mtk::Interpolation_Order get_time_interpolation_order() const
                {
                    return mTimeInterpolation->get_interpolation_order();
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the type of the space interpolation
                 */
                Interpolation_Type get_space_interpolation_type() const
                {
                    return mSpaceInterpolation->get_interpolation_type();
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the type of the time interpolation
                 */
                Interpolation_Type get_time_interpolation_type() const
                {
                    return mTimeInterpolation->get_interpolation_type();
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the number of space dimensions
                 */
                uint get_number_of_space_dimensions() const
                {
                    return mSpaceInterpolation->get_number_of_dimensions();
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the number of time dimensions
                 */
                uint get_number_of_time_dimensions() const
                {
                    return mTimeInterpolation->get_number_of_dimensions();
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the number of bases for the space function
                 */
                uint get_number_of_space_bases() const
                {
                    return mSpaceInterpolation->get_number_of_bases();
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the number of bases for the space function
                 */
                uint get_number_of_time_bases() const
                {
                    return mTimeInterpolation->get_number_of_bases();
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the space geometry type
                 */
                mtk::Geometry_Type get_space_geometry_type()
                {
                    return mGeometryType;
                }

                //------------------------------------------------------------------------------
                /**
                 * set the space and time coefficients of the geometry field xHat, tHat
                 * @param[ in ] space coefficients
                 * @param[ in ] time coefficients
                 */
                void set_coeff(
                        const Matrix< DDRMat > & aXHat,
                        const Matrix< DDRMat > & aTHat );

                //------------------------------------------------------------------------------
                /**
                 * set the space coefficients of the geometry field xHat
                 * @param[ in ] space coefficients
                 */
                void set_space_coeff( const Matrix< DDRMat > & aXHat );

                //------------------------------------------------------------------------------
                /**
                 * set the time coefficients of the geometry field tHat
                 * @param[ in ] time coefficients
                 */
                void set_time_coeff( const Matrix< DDRMat > & aTHat );

                //------------------------------------------------------------------------------
                /**
                 * get the space coefficients of the geometry field xHat
                 */
                const Matrix< DDRMat > & get_space_coeff() const
                  {
                    // check that mXHat is set
                    MORIS_ASSERT( mXHat.numel()>0, "Geometry_Interpolator::get_space_coeff - mXHat is not set." );

                    return mXHat;
                  }

                //------------------------------------------------------------------------------
                /**
                 * get the time coefficients of the geometry field tHat
                 */
                const Matrix< DDRMat > & get_time_coeff() const
                   {
                    // check that mTHat is set
                    MORIS_ASSERT( mTHat.numel()>0, "Geometry_Interpolator::get_time_coeff - mTHat is not set." );

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
                        const Matrix< DDRMat > & aXiHat,
                        const Matrix< DDRMat > & aTauHat );

                //------------------------------------------------------------------------------
                /**
                 * set the space param coefficients of the geometry field xiHat
                 * @param[ in ] space coefficients
                 */
                void set_space_param_coeff( const Matrix< DDRMat > & aXiHat );

                //------------------------------------------------------------------------------
                /**
                 * set the time param coefficients of the geometry field tHat
                 * @param[ in ] time coefficients
                 */
                void set_time_param_coeff( const Matrix< DDRMat > & aTauHat );

                //------------------------------------------------------------------------------
                /**
                 * get the space parametric coefficients of the geometry field xiHat
                 */
                const Matrix< DDRMat > & get_space_param_coeff() const
                {
                    // check that mXiHat is set
                    MORIS_ASSERT( mXiHat.numel()>0, "Geometry_Interpolator::get_space_param_coeff - mXiHat is not set." );

                    return mXiHat;
                }

                //------------------------------------------------------------------------------
                /**
                 * get the time parametric coefficients of the geometry field tauHat
                 */
                const Matrix< DDRMat > & get_time_param_coeff() const
                {
                    // check that mTauHat is set
                    MORIS_ASSERT( mTauHat.numel()>0, "Geometry_Interpolator::get_time_param_coeff - mTauHat is not set." );

                    return mTauHat;
                }

                //------------------------------------------------------------------------------
                /**
                 * set the parametric point where geometry is interpolated
                 * @param[ in ] aParamPoint evaluation point in space and time
                 */
                void set_space_time( const Matrix< DDRMat > & aParamPoint );
                void set_space( const Matrix< DDRMat > & aSpaceParamPoint );
                void set_time( const Matrix< DDRMat > & aTimeParamPoint );

                void get_space_time( Matrix< DDRMat > & aParamPoint )
                {
                    aParamPoint.set_size( mNumSpaceParamDim + 1, 1 );
                    aParamPoint( { 0, mNumSpaceParamDim-1 }, { 0, 0 } ) = mXiLocal.matrix_data();
                    aParamPoint( mNumSpaceParamDim ) = mTauLocal( 0 );
                }

                //------------------------------------------------------------------------------
                /**
                 * get the vertices ordinals of each face of the parent element
                 */
                moris::Cell< moris::Cell< moris_index > > get_face_vertices_ordinals( uint aNumSpaceBases );

                //------------------------------------------------------------------------------
                /**
                 * get the parametric coordinates of a space side
                 * @param[ in ] a space face ordinal
                 */
                Matrix< DDRMat > extract_space_side_space_param_coeff(
                        moris_index              aSpaceOrdinal,
                        mtk::Interpolation_Order aInterpolationOrder );

                //------------------------------------------------------------------------------
                /**
                 * get the parametric coordinates in space
                 */
                Matrix< DDRMat > extract_space_param_coeff( mtk::Interpolation_Order aInterpolationOrder );

                //------------------------------------------------------------------------------
                /**
                 * gets the space shape functions at a given evaluation point
                 * @param[ out ] aNXi shape functions ( 1 x <number of nodes> )
                 */
                const Matrix< DDRMat > & NXi();

                /**
                 * evaluates the space shape functions at a given evaluation point
                 */
                void eval_NXi();

                //------------------------------------------------------------------------------
                /**
                 * gets the time shape functions at a given evaluation point
                 * @param[ out ] aNTau shape functions ( 1 x <number of nodes> )
                 */
                const Matrix< DDRMat > & NTau();

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
                const Matrix< DDRMat > & dNdXi();

                /**
                 * evaluates the first derivatives of the space shape functions
                 * wrt parametric coordinates at a given evaluation point
                 */
                void eval_dNdXi();

                //------------------------------------------------------------------------------
                /**
                 * gets the first derivatives of the time shape functions
                 * wrt parametric coordinates at a given evaluation point
                 * @param[ out ] adNdTau first order derivatives ( <number of dimensions> x <number of nodes> )
                 */
                const Matrix< DDRMat > & dNdTau();

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
                const Matrix< DDRMat > & d2NdXi2();

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
                const Matrix< DDRMat > & d3NdXi3();

                /**
                 * evaluates the third derivatives of the space shape functions
                 * wrt parametric coordinates at a given evaluation point
                 */
                void eval_d3NdXi3();

                //------------------------------------------------------------------------------
                /**
                 * gets the second derivatives of the time shape functions
                 * wrt parametric coordinates at a given evaluation point
                 * @param[ in ] ad2NdTau2 second order derivatives ( <number of dimensions> x <number of nodes> )
                 */
                const Matrix< DDRMat > & d2NdTau2();

                /**
                 * evaluates the second derivatives of the time shape functions
                 * wrt parametric coordinates at a given evaluation point
                 */
                void eval_d2NdTau2();

                //------------------------------------------------------------------------------
                /**
                 * evaluates the geometry Jacobian in space
                 * @param[ in ] aJt     transposed of geometry Jacobian in space
                 */
                void space_jacobian( Matrix< DDRMat > & aJt );

                //------------------------------------------------------------------------------
                /**
                 * evaluates the 2nd geometry Jacobian in space
                 * @param[ in ] aJ2bt     2nd geometry Jacobian in space
                 */
                void second_space_jacobian( Matrix< DDRMat > & aJ2bt );

                //------------------------------------------------------------------------------
                /**
                 * evaluates the geometry Jacobian in space
                 * @param[ in ] aJ3ct     3rd geometry Jacobian in space
                 */
                void third_space_jacobian( Matrix< DDRMat > & aJ3ct );

                //------------------------------------------------------------------------------
                /**
                 * evaluates the geometry Jacobian in time
                 * @param[ in ] aJt     transposed of geometry Jacobian in time
                 */
                void time_jacobian( Matrix< DDRMat > & aJt );

                //------------------------------------------------------------------------------
                /**
                 * evaluates the determinant of the Jacobian mapping
                 * at given space and time evaluation point
                 */
                real det_J();

                //------------------------------------------------------------------------------
                /**
                 * evaluates the normal to the side
                 * in the case of a space side interpolation
                 * at given space and time evaluation point
                 * @param[ in ]  aNormal normal to be filled
                 */
                void get_normal( Matrix< DDRMat > & aNormal );

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
                        Matrix< DDRMat > & aJt,
                        Matrix< DDRMat > & aKt,
                        Matrix< DDRMat > & aLt,
                        const Matrix< DDRMat > & adNdXi,
                        const Matrix< DDRMat > & ad2NdXi2 );


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
                        Matrix< DDRMat > & aJt,
                        Matrix< DDRMat > & aKt,
                        Matrix< DDRMat > & aLt,
                        const Matrix< DDRMat > & adNdTau,
                        const Matrix< DDRMat > & ad2NdTau2 );

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
                        Matrix< DDRMat > & aJt,
                        Matrix< DDRMat > & aJ2bt,
                        Matrix< DDRMat > & aJ3at,
                        Matrix< DDRMat > & aJ3bt,
                        Matrix< DDRMat > & aJ3ct,
                        const Matrix< DDRMat > & adNdXi,
                        const Matrix< DDRMat > & ad2NdXi2,
                        const Matrix< DDRMat > & ad3NdXi3 );

                //------------------------------------------------------------------------------
                /**
                 * evaluates the space geometry field at a given evaluation point in space
                 * @param[ out ] aX   location in space
                 */
                Matrix< DDRMat > valx();

                //------------------------------------------------------------------------------
                /**
                 * evaluates the space geometry field at a given evaluation point in time
                 * @param[ out ] aT   location in time
                 */
                Matrix< DDRMat > valt();

                //------------------------------------------------------------------------------
                /**
                 * map an integration point from local param coords to global param coords
                 * @param[ in ] aGlobalParamPoint param coords in global parametric space
                 */
                void map_integration_point( Matrix< DDRMat > & aGlobalParamPoint );

                //------------------------------------------------------------------------------
                /**
                 * update parametric coordinates (xi, eta, zeta)
                 * for given physical coordinates (x, y, z)
                 * @param[ in ] aPhysCoordinates  coords in physical space
                 * @param[ in ] aParamCoordinates coords in parametric space
                 */
                void update_local_coordinates(
                        Matrix< DDRMat > & aPhysCoordinates,
                        Matrix< DDRMat > & aParamCoordinates );

                //------------------------------------------------------------------------------
            private:
                //------------------------------------------------------------------------------
                /**
                 * sets the function pointers for 2D and 3D. Called during construction.
                 */
                void set_function_pointers();

                //------------------------------------------------------------------------------
                /**
                 * evaluate space detJ
                 */
                real eval_space_detJ_side_line( Matrix< DDRMat > & aSpaceJt );
                real eval_space_detJ_side_tri( Matrix< DDRMat > & aSpaceJt );
                real eval_space_detJ_side_quad( Matrix< DDRMat > & aSpaceJt );

                real eval_space_detJ_bulk_line_quad_hex( Matrix< DDRMat > & aSpaceJt );
                real eval_space_detJ_bulk_tri( Matrix< DDRMat > & aSpaceJt );
                real eval_space_detJ_bulk_tet( Matrix< DDRMat > & aSpaceJt );

                //------------------------------------------------------------------------------
                /**
                 * evaluate time detJ
                 */
                real eval_time_detJ_side( Matrix< DDRMat > & aTimeJt );
                real eval_time_detJ_bulk( Matrix< DDRMat > & aTimeJt );

                //------------------------------------------------------------------------------
                /**
                 * evaluate normal to side
                 */
                void eval_normal_side_line(
                        Matrix< DDRMat > & aRealTangents,
                        Matrix< DDRMat > & aNormal );
                void eval_normal_side_quad(
                        Matrix< DDRMat > & aRealTangents,
                        Matrix< DDRMat > & aNormal );
                void eval_normal_side_tri(
                        Matrix< DDRMat > & aRealTangents,
                        Matrix< DDRMat > & aNormal );

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
                        const Matrix< DDRMat > & aJt,
                        Matrix< DDRMat > & aKt,
                        Matrix< DDRMat > & aLt,
                        const Matrix< DDRMat > & ad2NdXi2,
                        const Matrix< DDRMat > & aXHat );

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
                        const Matrix< DDRMat > & aJt,
                        Matrix< DDRMat > & aKt,
                        Matrix< DDRMat > & aLt,
                        const Matrix< DDRMat > & ad2NdXi2,
                        const Matrix< DDRMat > & aXHat );

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
                        const Matrix< DDRMat > & aJt,
                        Matrix< DDRMat > & aKt,
                        Matrix< DDRMat > & aLt,
                        const Matrix< DDRMat > & ad2NdXi2,
                        const Matrix< DDRMat > & aXHat );

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
                        const Matrix< DDRMat > & aJt,
                        const Matrix< DDRMat > & aJ2bt,
                        Matrix< DDRMat > & aJ3at,
                        Matrix< DDRMat > & aJ3bt,
                        Matrix< DDRMat > & aJ3ct,
                        const Matrix< DDRMat > & ad3NdXi3,
                        const Matrix< DDRMat > & aXHat );

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
                        const Matrix< DDRMat > & aJt,
                        const Matrix< DDRMat > & aJ2bt,
                        Matrix< DDRMat > & aJ3at,
                        Matrix< DDRMat > & aJ3bt,
                        Matrix< DDRMat > & aJ3ct,
                        const Matrix< DDRMat > & ad3NdXi3,
                        const Matrix< DDRMat > & aXHat );

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
                        const Matrix< DDRMat > & aJt,
                        const Matrix< DDRMat > & aJ2bt,
                        Matrix< DDRMat > & aJ3at,
                        Matrix< DDRMat > & aJ3bt,
                        Matrix< DDRMat > & aJ3ct,
                        const Matrix< DDRMat > & ad3NdXi3,
                        const Matrix< DDRMat > & aXHat );

                //------------------------------------------------------------------------------

        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_ */

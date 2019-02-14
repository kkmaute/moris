/*
 * cl_FEM_Geometry_Interpolator_Bis.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_BIS_HPP_
#define SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_BIS_HPP_


#include "cl_FEM_Interpolation_Matrix.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Rule_Bis.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
    /**
     * \brief a special interpolation class for geometry
     */
    class Geometry_Interpolator_Bis
    {
        //! pointer to interpolation function object
        Interpolation_Function_Base * mInterpolation;


        //! pointer to function for second derivative
        void
        ( * mSecondDerivativeMatrices )(
                        const Matrix< DDRMat > 		& aJt,
                              Matrix< DDRMat > 		& aKt,
                              Matrix< DDRMat > 		& aLt,
                        const Interpolation_Matrix  & adNdXi,
                        const Interpolation_Matrix  & ad2NdXi2,
                        const Matrix< DDRMat > 		& aXhat );

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * default constructor
         *
         * @param[ in ] interpolation rule
         */
        Geometry_Interpolator_Bis(const Interpolation_Rule_Bis  & aInterpolationRule,
        					      const bool mIsSpace);
//------------------------------------------------------------------------------

        /**
         * destructor
         */
        ~Geometry_Interpolator_Bis();

//------------------------------------------------------------------------------

        /**
         * returns the order of the interpolation
         */
        mtk::Interpolation_Order
        get_interpolation_order() const;

//------------------------------------------------------------------------------

        /**
         * returns the order of the interpolation
         */
        Interpolation_Type
        get_interpolation_type() const;

 //------------------------------------------------------------------------------

        /**
         * evaluates the shape function at a given point
         *
         * @param[ out ] aN  shape function as
         *                   ( 1 x <number of nodes> )
         *
         * @param[ in ]  aXi parameter coordinates
         *                   ( <number of dimensions>  x 1 )
         */
        void
        eval_N(       Interpolation_Matrix  & aN,
                const Matrix< DDRMat >      & aXi ) const;

//------------------------------------------------------------------------------
        /**
        * evaluates the shape function at a given point
        *
        * @param[ in ]  aN  shape function as
        *                   ( 1 x <number of nodes> )
        * @param[ in ]  aXi parameter coordinates
        *                   ( <number of dimensions>  x 1 )
        * @param[ in ]  aXHat node coordinates
        *                   ( <number of nodes>  x 1 )
        * @param[ out ] aX interpolated field value
        *                   ( <number of dimensions>  x 1 )
        */
        void
        eval_geometryField(		Matrix< DDRMat > 		& aX,
                    			Interpolation_Matrix    & aN,
						  const Matrix< DDRMat > 		& aXi,
						  const Matrix< DDRMat > 		& aXHat) const;

//------------------------------------------------------------------------------

        /**
         * calculates the first derivative of the shape function
         * in parameter space
         *
         * @param[ out ] adNdXi ( <number of dimensions> x <number of nodes> )
         *
         * @param[ in ] aXi    point where function is evaluated
         *                     ( <number of dimensions>  x 1 )
         *
         */
        void
        eval_dNdXi(
                      Interpolation_Matrix & adNdXi,
                const Matrix< DDRMat >     & aXi ) const;

//------------------------------------------------------------------------------

        /**
         * calculates the second derivative of the shape function
         * in parameter space
         *
         * @param[ out ] ad2NdXi2 ( <number of dimensions> x <number of nodes> )
         *
         * @param[ in ] aXi    point where function is evaluated
         *                     ( <number of dimensions>  x 1 )
         *
         */
        void
        eval_d2NdXi2 (
                      Interpolation_Matrix 	& ad2NdXi2,
                const Matrix< DDRMat >      & aXi ) const;

//------------------------------------------------------------------------------

        /**
         * evaluates the geometry Jacobian
         *
         * @param[ out ] aJt    transposed of geometry Jacobian
         *
         * @param[ in ] adNdXi  derivatives of N in parameter space
         *
         * @param[ in ] aXhat   pointer to object containing node coordinates
         *
         */
        void
        eval_jacobian(
                      Matrix< DDRMat >      & aJt,
                const Interpolation_Matrix  & adNdXi,
                const Matrix< DDRMat >      & aXhat ) const;

//------------------------------------------------------------------------------

        /**
         * evaluates the geometry Jacobian and the matrices needed for the second
         * derivative
         *
         * @param[ out ] aJt    transposed of geometry Jacobian
         *
         * @param[ in ] adNdXi  derivatives of N in parameter space
         *
         * @param[ in ] aXhat   pointer to object containing node coordinates
         *
         */
        void
        eval_jacobian_and_matrices_for_second_derivatives(
                              Matrix< DDRMat > & aJt,
                              Matrix< DDRMat > & aKt,
                              Matrix< DDRMat > & aLt,
                        const Interpolation_Matrix    & adNdXi,
                        const Interpolation_Matrix    & ad2NdXi2,
                        const Matrix< DDRMat >   & aXhat ) const;

//------------------------------------------------------------------------------

        /**
         * creates a pointer to one of the J, K or L matrices
         *
         * @param[ in ]  aDerivativeInSpace, 0, 1 or 2
         * @param[ in ]  aDerivativeInTime   0, 1 or 2
         * @param[ in ]  aCoeffsSwitch :
         *                      0: evaluated value
         *                      1: vector N, N_x or N_x2
         *                      2: Matrix J, L
         *                      3: Matrix K
         */
        Interpolation_Matrix *
        create_matrix_pointer(
                const uint & aDerivativeInSpace,
                const uint & aDerivativeInTime ) const;

//------------------------------------------------------------------------------

        /**
         * returns the number of basis for this function
         */
        uint
        get_number_of_basis() const;

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        /**
         * sets the function pointers for 2D and 3D. Called during construction.
         */
        void
        set_function_pointers();

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
         * @param[ in ]  adNdXi       first derivative in parameter space
         * @param[ in ]  ad2NdX2i     second derivative in parameter space
         *
         */
        static void
        eval_matrices_for_second_derivative_1d(
                const Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Interpolation_Matrix    & adNdXi,
                const Interpolation_Matrix    & ad2NdXi2,
                const Matrix< DDRMat > & aXhat );

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
         * @param[ in ]  adNdXi       first derivative in parameter space
         * @param[ in ]  ad2NdX2i     second derivative in parameter space
         *
         */
        static void
        eval_matrices_for_second_derivative_2d(
                const Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Interpolation_Matrix    & adNdXi,
                const Interpolation_Matrix    & ad2NdXi2,
                const Matrix< DDRMat > & aXhat );

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
         * @param[ in ]  adNdXi       first derivative in parameter space
         * @param[ in ]  ad2NdX2i     second derivative in parameter space
         *
         */
        static void
        eval_matrices_for_second_derivative_3d(
                const Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Interpolation_Matrix    & adNdXi,
                const Interpolation_Matrix    & ad2NdXi2,
                const Matrix< DDRMat > & aXhat );

//------------------------------------------------------------------------------

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_BIS_HPP_ */

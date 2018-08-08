/*
 * cl_FEM_Geometry_Interpolator.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_
#define SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_


#include "cl_FEM_Interpolation_Matrix.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

    // forward declaration of element
    class Element ;

//------------------------------------------------------------------------------
    /**
     * \brief a special interpolation class for geometry
     */
    class Geometry_Interpolator
    {
        //! pointer to interpolation function object
        Interpolation_Function_Base * mInterpolation;

        //! pointer to function for second derivative
        void
        ( * mSecondDerivativeMatrices )(
                        const Mat<real> & aJt,
                              Mat<real> & aKt,
                              Mat<real> & aLt,
                        const Interpolation_Matrix  & adNdXi,
                        const Interpolation_Matrix  & ad2NdXi2,
                        const Mat<real> & aXhat );

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * default constructor
         *
         * @param[ in ] element this function refers to
         * @param[ in ] order of interpolation
         */
        Geometry_Interpolator(
            Element                   * aElement,
            const Interpolation_Rule  & aInterpolationRule );

//------------------------------------------------------------------------------

        /**
         * destructor
         */
        ~Geometry_Interpolator();

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
                const Mat<real>             & aXi ) const;

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
                const Mat<real>            & aXi ) const;

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
                      Interpolation_Matrix & ad2NdXi2,
                const Mat<real>            & aXi ) const;

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
                      Mat<real>             & aJt,
                const Interpolation_Matrix  & adNdXi,
                const Mat<real>             & aXhat ) const;

//------------------------------------------------------------------------------

        /**
         * evaluates the geometry Jacobian and the matrices needed for the second
         * derovative
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
                              Mat< real > & aJt,
                              Mat< real > & aKt,
                              Mat< real > & aLt,
                        const Interpolation_Matrix    & adNdXi,
                        const Interpolation_Matrix    & ad2NdXi2,
                        const Mat<real>   & aXhat ) const;

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
                const uint & aDerivativeInTime, // does nothing right now
                const uint & aCoeffsSwitch ) const;

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
        eval_matrices_for_second_derivative_1d(
                const Mat< real > & aJt,
                      Mat< real > & aKt,
                      Mat< real > & aLt,
                const Interpolation_Matrix    & adNdXi,
                const Interpolation_Matrix    & ad2NdXi2,
                const Mat< real > & aXhat );

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
                const Mat< real > & aJt,
                      Mat< real > & aKt,
                      Mat< real > & aLt,
                const Interpolation_Matrix    & adNdXi,
                const Interpolation_Matrix    & ad2NdXi2,
                const Mat< real > & aXhat );

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
         * @param[ in ] aJt          transposed of geometry Jacobian
         * @param[ out ] aKt          transposed help matrix K
         * @param[ out ] aLt          transposed help matrix L
         * @param[ in ]  adNdXi       first derivative in parameter space
         * @param[ in ]  ad2NdX2i     second derivative in parameter space
         *
         */
        static void
        eval_matrices_for_second_derivative_3d(
                const Mat< real > & aJt,
                      Mat< real > & aKt,
                      Mat< real > & aLt,
                const Interpolation_Matrix    & adNdXi,
                const Interpolation_Matrix    & ad2NdXi2,
                const Mat< real > & aXhat );

//------------------------------------------------------------------------------

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_ */

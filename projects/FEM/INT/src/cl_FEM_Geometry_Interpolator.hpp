/*
 * cl_FEM_Geometry_Interpolator.hpp
 *
 *  Created on: Jan 31, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_
#define SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_

#include "typedefs.hpp" //MRS/COR/src

#include "cl_MTK_Enums.hpp" //MTK/src

#include "cl_FEM_Enums.hpp"              //FEM/INT/src
#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src

#include "linalg_typedefs.hpp"
#include "cl_Matrix.hpp" //LINALG/src
#include "op_times.hpp"  //LINALG/src
#include "op_plus.hpp"   //LINALG/src
#include "op_minus.hpp"  //LINALG/src
#include "fn_trans.hpp"  //LINALG/src
#include "fn_det.hpp"    //LINALG/src
#include "fn_inv.hpp"    //LINALG/src

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
        // pointer to space interpolation function object
        Interpolation_Function_Base * mSpaceInterpolation = nullptr;

        // pointer to time interpolation function object
        Interpolation_Function_Base * mTimeInterpolation = nullptr;

        // matrix of space coefficients xHat
        Matrix < DDRMat > mXHat;

        // matrix of time coefficients tHat
        Matrix < DDRMat > mTHat;

        // pointer to function for second derivative
        void ( * mSecondDerivativeMatricesSpace )( const Matrix< DDRMat > & aJt,
                                                         Matrix< DDRMat > & aKt,
                                                         Matrix< DDRMat > & aLt,
                                                   const Matrix< DDRMat > & ad2NdXi2,
                                                   const Matrix< DDRMat > & aXHat );

        void ( * mSecondDerivativeMatricesTime )( const Matrix< DDRMat > & aJt,
                                                        Matrix< DDRMat > & aKt,
                                                        Matrix< DDRMat > & aLt,
                                                  const Matrix< DDRMat > & ad2NdTau2,
                                                  const Matrix< DDRMat > & aTHat );


//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * default constructor
         *
         * @param[ in ] interpolation rule for geometry
         */
        Geometry_Interpolator( const Interpolation_Rule & aInterpolationRule );

//------------------------------------------------------------------------------
        /**
         * destructor
         */
        ~Geometry_Interpolator();

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
          * returns the number of basis for the space function
          */
          uint get_number_of_space_bases() const
          {
              return mSpaceInterpolation->get_number_of_bases();
          }

//------------------------------------------------------------------------------
          /**
           * returns the number of basis for the space function
           */
           uint get_number_of_time_bases() const
           {
               return mTimeInterpolation->get_number_of_bases();
           }

//------------------------------------------------------------------------------
         /**
          * set the coefficients of the geometry field xHat, tHat
          */
          void set_coeff( const Matrix< DDRMat > & aXHat,
                          const Matrix< DDRMat > & aTHat );

          void set_coeff( const Matrix< DDRMat > & aXHat );

//------------------------------------------------------------------------------
         /**
          * get the space coefficients of the geometry field xHat
          */
          Matrix< DDRMat > get_space_coeff() const
          {
              return mXHat;
          }

//------------------------------------------------------------------------------
          /**
           * get the time coefficients of the geometry field xHat
           */
           Matrix< DDRMat > get_time_coeff() const
           {
               return mTHat;
           }
//------------------------------------------------------------------------------
        /**
         * evaluates the space shape function at a given point
         * @param[ out ] aNXi shape function ( 1 x <number of nodes> )
         * @param[ in ]  aXi  parameter coordinates ( <number of dimensions> x 1 )
         */
        Matrix < DDRMat > NXi( const Matrix< DDRMat > & aXi ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the time shape function at a given point
         * @param[ out ] aNTau shape function as ( 1 x <number of nodes> )
         * @param[ in ]  aTau  parameter coordinates ( <number of dimensions> x 1 )
         */
        Matrix < DDRMat > NTau( const Matrix< DDRMat > & aTau ) const;

//------------------------------------------------------------------------------
        /**
         * calculates the first derivative of the shape function
         * in parameter space
         *
         * @param[ out ] adNdXi ( <number of dimensions> x <number of nodes> )
         *
         * @param[ in ]  aXi    point where function is evaluated
         *                     ( <number of dimensions>  x 1 )
         */
        Matrix< DDRMat > dNdXi( const Matrix< DDRMat > & aXi ) const;

//------------------------------------------------------------------------------
        /**
         * calculates the first derivative of the shape function
         * in parameter space
         *
         * @param[ out ] adNdTau ( <number of dimensions> x <number of nodes> )
         *
         * @param[ in ] aTau    point where function is evaluated
         *                     ( <number of dimensions>  x 1 )
         */
         Matrix< DDRMat > dNdTau( const Matrix< DDRMat > & aTau ) const;

//------------------------------------------------------------------------------
        /**
         * calculates the second derivative of the shape function
         * in parameter space
         *
         * @param[ out ] ad2NdXi2 ( <number of dimensions> x <number of nodes> )
         *
         * @param[ in ] aXi    point where function is evaluated
         *                     ( <number of dimensions>  x 1 )
         */
        Matrix< DDRMat > d2NdXi2 ( const Matrix< DDRMat > & aXi ) const;

//------------------------------------------------------------------------------
        /**
         * calculates the second derivative of the shape function
         * in parameter space
         *
         * @param[ out ] ad2NdTau2 ( <number of dimensions> x <number of nodes> )
         *
         * @param[ in ] aTau    point where function is evaluated
         *                     ( <number of dimensions>  x 1 )
         */
        Matrix< DDRMat > d2NdTau2 ( const Matrix< DDRMat > & aTau ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the geometry Jacobian
         *
         * @param[ out ] tJt    transposed of geometry Jacobian
         *
         * @param[ in ] adNdXi  derivatives of N in parameter space
         *         *
         */
        Matrix< DDRMat > space_jacobian( const Matrix< DDRMat > & adNdXi ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the geometry Jacobian
         *
         * @param[ out ] tJt    transposed of geometry Jacobian
         *
         * @param[ in ] adNdTau  derivatives of N in parameter space
         *         *
         */
        Matrix< DDRMat > time_jacobian( const Matrix< DDRMat > & adNdTau ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the determinant of the Jacobian mapping
         * at given space and time Xi, Tau
         */
        real det_J( const Matrix< DDRMat > & aParamPoint );

//------------------------------------------------------------------------------
        /**
         * evaluates the geometry Jacobian and the matrices needed for the second
         * derivative
         *
         * @param[ out ] aJt    transposed of geometry Jacobian
         *
         * @param[ in ] adNdXi  derivatives of N in parameter space
         *
         */
        void space_jacobian_and_matrices_for_second_derivatives(       Matrix< DDRMat > & aJt,
                                                                       Matrix< DDRMat > & aKt,
                                                                       Matrix< DDRMat > & aLt,
                                                                 const Matrix< DDRMat > & adNdXi,
                                                                 const Matrix< DDRMat > & ad2NdXi2 ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the geometry Jacobian and the matrices needed for the second
         * derivative
         *
         * @param[ out ] aJt    transposed of geometry Jacobian
         * @param[ in ] adNdXi  derivatives of N in parameter space
         *
         */
        void time_jacobian_and_matrices_for_second_derivatives(       Matrix< DDRMat > & aJt,
                                                                      Matrix< DDRMat > & aKt,
                                                                      Matrix< DDRMat > & aLt,
                                                                const Matrix< DDRMat > & adNdTau,
                                                                const Matrix< DDRMat > & ad2NdTau2 ) const;

//------------------------------------------------------------------------------
        /**
         * evaluates the space geometry field at xi
         */
         Matrix< DDRMat > valx( const Matrix< DDRMat > & aXi );

//------------------------------------------------------------------------------
        /**
         * evaluates the time geometry field at tau
         */
        Matrix< DDRMat > valt( const Matrix< DDRMat > & aTau );

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------
        /**
         * sets the function pointers for 2D and 3D. Called during construction.
         */
        void set_function_pointers();

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
        static void eval_matrices_for_second_derivative_1d( const Matrix< DDRMat > & aJt,
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
         * @param[ in ]  adNdXi       first derivative in parameter space
         * @param[ in ]  ad2NdX2i     second derivative in parameter space
         *
         */
        static void eval_matrices_for_second_derivative_2d( const Matrix< DDRMat > & aJt,
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
         * @param[ in ]  adNdXi       first derivative in parameter space
         * @param[ in ]  ad2NdX2i     second derivative in parameter space
         *
         */
        static void eval_matrices_for_second_derivative_3d( const Matrix< DDRMat > & aJt,
                                                                  Matrix< DDRMat > & aKt,
                                                                  Matrix< DDRMat > & aLt,
                                                            const Matrix< DDRMat > & ad2NdXi2,
                                                            const Matrix< DDRMat > & aXHat );

//------------------------------------------------------------------------------

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_GEOMETRY_INTERPOLATOR_HPP_ */

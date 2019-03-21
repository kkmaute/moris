/*
 * cl_FEM_IWG_L2.hpp
 *
 *  Created on: Aug 13, 2018
 *      Author: messe
 */
#ifndef SRC_FEM_CL_FEM_IWG_L2_HPP_
#define SRC_FEM_CL_FEM_IWG_L2_HPP_

#include "typedefs.hpp"                     //MRS/COR/src

#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class IWG_L2 : public IWG
        {

            // Alpha-Parameter, for J = M + alpha*K
            real mAlpha;

            // pointer to interpolator
            Field_Interpolator * mInterpolator = nullptr;

            void
            ( IWG_L2:: * mComputeFunction )( Cell< Matrix< DDRMat > >    & aJacobians,
                                             Matrix< DDRMat >            & aResidual,
                                             Cell< Field_Interpolator* > & aFieldInterpolators);

            void
            ( IWG_L2:: * mComputeJacFunction )( Cell< Matrix< DDRMat > >    & aJacobians,
                                                Cell< Field_Interpolator* > & aFielInterpolators );

            void
            ( IWG_L2:: * mComputeResFunction )( Matrix< DDRMat >            & aResidual,
                                                Cell< Field_Interpolator* > & aFielInterpolators );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            //constructor
             IWG_L2( const real aAlpha = 0.0 );

//------------------------------------------------------------------------------
            // trivial destructor
            ~IWG_L2(){};

//------------------------------------------------------------------------------
            /**
             * set the alpha value for the L2 parameter so that
             * A = M + alpha * K
             */
            void set_alpha( const real aAlpha );

//------------------------------------------------------------------------------
            /**
             * get the alpha value for the L2 parameter
             */
            real get_alpha() const
            {
                return mAlpha;
            }

//------------------------------------------------------------------------------
            /**
             * compute jacobian and residual
             */
            void compute_jacobian_and_residual( Cell< Matrix< DDRMat > >    & aJacobian,
                                                Matrix< DDRMat >            & aResidual,
                                                Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
            /**
             * compute jacobian
             */
            void compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                   Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
            /**
             * compute residual
             */
            void compute_residual( Matrix< DDRMat >            & aResidual,
                                   Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
//            /**
//             * calculates the square of the error at a given point
//             */
//            real
//            compute_integration_error( const Matrix< DDRMat > & aNodalDOF,
//                    real (*aFunction)( const Matrix< DDRMat > & aPoint ) ,
//                                       const uint             & aPointIndex );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------
            /**
             * j = N'*N, r =
             */
            void compute_jacobian_and_residual_without_alpha( Cell< Matrix< DDRMat > >    & aJacobians,
                                                              Matrix< DDRMat >            & aResidual,
                                                              Cell< Field_Interpolator* > & aFieldInterpolators );
//------------------------------------------------------------------------------
            /**
             * j = N'*N + alpha * B'*B
             */
            void compute_jacobian_and_residual_with_alpha( Cell< Matrix< DDRMat > >    & aJacobians,
                                                           Matrix< DDRMat >            & aResidual,
                                                           Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
            /**
             * j = N'*N
             */
            void compute_jacobian_without_alpha( Cell< Matrix< DDRMat > >    & aJacobians,
                                                 Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
            /**
             * j = N'*N + alpha * B'*B
             */
            void compute_jacobian_with_alpha( Cell< Matrix< DDRMat > >    & aJacobians,
                                              Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------

            void compute_residual_without_alpha( Matrix< DDRMat >            & aResidual,
                                                 Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------

            void compute_residual_with_alpha( Matrix< DDRMat >            & aResidual,
                                              Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
//            real
//            interpolate_scalar_at_point(
//                                const Matrix< DDRMat > & aNodalWeakBC,
//                                const uint             & aPointIndex );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_L2_HPP_ */

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

//            // N-Matrix
//            //Interpolation_Matrix * mN = nullptr;
//            Matrix< DDRMat > * mN = nullptr;
//
//            // B-Matrix
//            //Interpolation_Matrix * mB = nullptr;
//            Matrix< DDRMat > * mB = nullptr;

//            void
//            ( IWG_L2:: * mComputeFunction )(
//                    Matrix< DDRMat >       & aJacobian,
//                    Matrix< DDRMat >       & aResidual,
//                    const Matrix< DDRMat > & aNodalDOF,
//                    const Matrix< DDRMat > & aNodalWeakBC,
//                    const uint             & aPointIndex );

            void
            ( IWG_L2:: * mComputeFunction )(       Matrix< DDRMat > & aJacobian,
                                                   Matrix< DDRMat > & aResidual,
                                             const Matrix< DDRMat > & aNodalDOF,
                                             const Matrix< DDRMat > & aNodalWeakBC );

            void
            ( IWG_L2:: * mComputeJacFunction )( Cell< Matrix< DDRMat > >    & aJacobians,
                                                Cell< Field_Interpolator* > & aFielInterpolators );

            void
            ( IWG_L2:: * mComputeResFunction )( Matrix< DDRMat >            & aResidual,
                                                Cell< Field_Interpolator* > & aFielInterpolators );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /*
             *  constructor
             *
             *  J = N'*N + alpha* B'*B
             */
             IWG_L2( const real aAlpha = 0.0 );

//------------------------------------------------------------------------------

            // trivial destructor
            ~IWG_L2(){};

//------------------------------------------------------------------------------

//            /**
//             * returns a cell with the dof types, assuming that all nodes
//             * have the same type
//             */
//            Cell< MSI::Dof_Type > get_dof_types()
//            {
//                Cell< MSI::Dof_Type > aDofTypes( 1, MSI::Dof_Type::L2 );
//                return aDofTypes;
//            }

//------------------------------------------------------------------------------

            /**
             * set the alpha value for the L2 parameter so that
             *
             * A = M + alpha * K
             */
            void set_alpha( const real aAlpha );

//------------------------------------------------------------------------------

            real get_alpha() const
            {
                return mAlpha;
            }

//------------------------------------------------------------------------------

//            void create_matrices( Field_Interpolator * aInterpolator );

//------------------------------------------------------------------------------

//            void delete_matrices();

//------------------------------------------------------------------------------
//
//            void compute_jacobian_and_residual(
//                    Matrix< DDRMat >       & aJacobian,
//                    Matrix< DDRMat >       & aResidual,
//                    const Matrix< DDRMat > & aNodalDOF,
//                    const Matrix< DDRMat > & aNodalWeakBC,
//                    const uint             & aPointIndex );

            void compute_jacobian_and_residual(       Matrix< DDRMat > & aJacobian,
                                                      Matrix< DDRMat > & aResidual,
                                                const Matrix< DDRMat > & aNodalDOF,
                                                const Matrix< DDRMat > & aNodalWeakBC );
//------------------------------------------------------------------------------
//
//            /**
//             * calculates the square of the error at a given point
//             */
//            real
//            compute_integration_error( const Matrix< DDRMat > & aNodalDOF,
//                    real (*aFunction)( const Matrix< DDRMat > & aPoint ) ,
//                                       const uint             & aPointIndex );

//------------------------------------------------------------------------------

            void compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                   Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------

            void compute_residual( Matrix< DDRMat >            & aResidual,
                                   Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------
//
//            /**
//             * J = N'*N
//             */
//            void
//            compute_jacobian_and_residual_without_alpha(
//                    Matrix< DDRMat >       & aJacobian,
//                    Matrix< DDRMat >       & aResidual,
//                    const Matrix< DDRMat > & aNodalDOF,
//                    const Matrix< DDRMat > & aNodalWeakBC,
//                    const uint             & aPointIndex );


            /**
             * J = N'*N
             */
            void compute_jacobian_and_residual_without_alpha(       Matrix< DDRMat > & aJacobian,
                                                                    Matrix< DDRMat > & aResidual,
                                                              const Matrix< DDRMat > & aNodalDOF,
                                                              const Matrix< DDRMat > & aNodalWeakBC );

//------------------------------------------------------------------------------
//
//            /**
//             * J = N'*N + alpha * B'*B
//             */
//            void
//            compute_jacobian_and_residual_with_alpha(
//                    Matrix< DDRMat >       & aJacobian,
//                    Matrix< DDRMat >       & aResidual,
//                    const Matrix< DDRMat > & aNodalDOF,
//                    const Matrix< DDRMat > & aNodalWeakBC,
//                    const uint        & aPointIndex );

            /**
             * J = N'*N + alpha * B'*B
             */
            void compute_jacobian_and_residual_with_alpha(       Matrix< DDRMat > & aJacobian,
                                                                 Matrix< DDRMat > & aResidual,
                                                           const Matrix< DDRMat > & aNodalDOF,
                                                           const Matrix< DDRMat > & aNodalWeakBC );

//------------------------------------------------------------------------------
//
//            real
//            interpolate_scalar_at_point(
//                                const Matrix< DDRMat > & aNodalWeakBC,
//                                const uint             & aPointIndex );

//------------------------------------------------------------------------------
            void compute_jacobian_without_alpha( Cell< Matrix< DDRMat > >    & aJacobians,
                                                 Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------

            void compute_jacobian_with_alpha( Cell< Matrix< DDRMat > >    & aJacobians,
                                              Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------

            void compute_residual_without_alpha( Matrix< DDRMat >            & aResidual,
                                                 Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------

            void compute_residual_with_alpha( Matrix< DDRMat >            & aResidual,
                                              Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_L2_HPP_ */

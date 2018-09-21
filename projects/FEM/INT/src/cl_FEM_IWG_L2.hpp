/*
 * cl_FEM_Interpolator.hpp
 *
 *  Created on: Aug 13, 2018
 *      Author: messe
 */
#ifndef SRC_FEM_CL_FEM_IWG_L2_TEST_HPP_
#define SRC_FEM_CL_FEM_IWG_L2_TEST_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "fn_norm.hpp"
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src
#include "cl_FEM_Interpolation_Matrix.hpp"  //FEM/INT/src
#include "cl_FEM_Interpolator.hpp"          //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src
#include "op_times.hpp"                     //LINALG/src
#include "fn_trans.hpp"                     //LINALG/src
namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class IWG_L2 : public IWG
        {
            // pointer to interpolator
            Interpolator         * mInterpolator = nullptr;

            // N-Matrix
            Interpolation_Matrix * mN = nullptr;

            // B-Matrix
            Interpolation_Matrix * mB = nullptr;

            // Alpha-Parameter, for J = M + alpha*K
            real               mAlpha = 0.0; //10000.0;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            // trivial constructor
            IWG_L2() // : IWG()
            {};

//------------------------------------------------------------------------------

            // trivial destructor
            ~IWG_L2(){};

//------------------------------------------------------------------------------

            /**
             * returns a cell with the dof types, assuming that all nodes
             * have the same type
             */
            Cell< MSI::Dof_Type >
            get_dof_types()
            {
                Cell< MSI::Dof_Type > aDofTypes( 1, MSI::Dof_Type::L2 );

                return aDofTypes;
            }

//------------------------------------------------------------------------------

            void
            create_matrices( Interpolator * aInterpolator )
            {
                // copy pointer to interpolator class
                mInterpolator = aInterpolator;

                // create N-Matrix
                mN = aInterpolator->create_matrix( 0, 0 );

                // create B-Matrix
                mB = aInterpolator->create_matrix( 1, 0 );

            }

//------------------------------------------------------------------------------

            void
            delete_matrices()
            {
                delete mN;
                delete mB;
            }

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_residual(
                    Matrix< DDRMat >       & aJacobian,
                    Matrix< DDRMat >       & aResidual,
                    const Matrix< DDRMat > & aNodalDOF,
                    const Matrix< DDRMat > & aNodalWeakBC,
                    const uint        & aPointIndex )
            {

                // get shape function
                mN->compute( aPointIndex );

                if ( mAlpha > 0.0 )
                {
                    // compute derivative
                    mB->compute( aPointIndex );

                    // calculate Jacobian
                    aJacobian = trans( mN ) * mN
                            + mAlpha * ( trans( mB ) * mB );
                }
                else
                {
                    // calculate Jacobian
                    aJacobian = trans( mN ) * mN;
                }

                aResidual = aJacobian * ( aNodalWeakBC - aNodalDOF );

            }
//------------------------------------------------------------------------------

            /**
             * calculates the square of the error at a given point
             */
            real
            compute_integration_error(
                    const Matrix< DDRMat > & aNodalDOF,
                    real (*aFunction)( const Matrix< DDRMat > & aPoint ) ,
                    const uint        & aPointIndex )
            {
                mN->compute( aPointIndex );

                //Matrix< DDRMat > tPoint = mInterpolator->get_point( aPointIndex );
                //Matrix< DDRMat > tCoords = mInterpolator->eval_geometry_coords( tPoint );
                //Matrix< DDRMat > tCoords =  mN->matrix_data() * mInterpolator->get_node_coords().matrix_data()
                Matrix< DDRMat > tCoords = mN * mInterpolator->get_node_coords();
                // get shape function

                Matrix< DDRMat > tPhiHat = mN->matrix_data() * aNodalDOF.matrix_data();

                return std::pow(
                        tPhiHat( 0 )
                        - aFunction( tCoords ), 2 );

            }
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_L2_TEST_HPP_ */

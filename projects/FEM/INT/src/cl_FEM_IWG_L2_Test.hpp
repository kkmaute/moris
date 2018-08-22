/*
 * cl_FEM_Interpolator.hpp
 *
 *  Created on: Aug 13, 2018
 *      Author: messe
 */
#ifndef SRC_FEM_CL_FEM_IWG_L2_TEST_HPP_
#define SRC_FEM_CL_FEM_IWG_L2_TEST_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Mat.hpp"                       //LNA/src
#include "cl_FEM_Interpolation_Matrix.hpp"  //FEM/INT/src
#include "cl_FEM_Interpolator.hpp"          //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class IWG_L2_Test : public IWG
        {
            // pointer to interpolator
            Interpolator         * mInterpolator = nullptr;

            // N-Matrix
            Interpolation_Matrix * mN = nullptr;
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            // trivial constructor
            IWG_L2_Test() // : IWG()
            {};

//------------------------------------------------------------------------------

            // trivial destructor
            ~IWG_L2_Test(){};

//------------------------------------------------------------------------------

            /**
             * returns a cell with the dof types, assuming that all nodes
             * have the same type
             */
            Cell< MSI::Dof_Type >
            get_dof_types()
            {
                Cell< MSI::Dof_Type > aDofTypes( 1, MSI::Dof_Type::TEMP );

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
            }

//------------------------------------------------------------------------------

            void
            delete_matrices()
            {
                delete mN;
            }

//------------------------------------------------------------------------------
            void
            compute_jacobian_and_residual(
                    Mat< real >       & aJacobian,
                    Mat< real >       & aResidual,
                    const Mat< real > & aNodalDOF,
                    const Mat< real > & aNodalWeakBC,
                    const uint        & aPointIndex )
            {
                // get shape function
                mN->compute( aPointIndex );


                // calculate Jacobian
                aJacobian = trans( mN->data() ) * mN->data();

                // get point
                //auto tPoint = mInterpolator->eval_geometry_coords( aPointIndex );

                // circle function
                //real tCircle = tPoint.norm(); // - 1.0;
                // residual ( sign ? )
                //aR = trans( mN->data() )*tCircle - aJ*aU;
                //aR = trans( mN->data() )*tCircle;
                aResidual = aJacobian * ( aNodalWeakBC - aNodalDOF );

            }
//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_L2_TEST_HPP_ */

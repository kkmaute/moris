/*
 * cl_FEM_Interpolator.hpp
 *
 *  Created on: Aug 13, 2018
 *      Author: messe
 */
#ifndef SRC_FEM_CL_FEM_IWG_HPP_
#define SRC_FEM_CL_FEM_IWG_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Mat.hpp"                       //LNA/src
#include "cl_FEM_Interpolation_Matrix.hpp"  //FEM/INT/src
#include "cl_FEM_Interpolator.hpp"          //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        /**
         * Integrand of Weak Form of Governing Euqations
         */
        class IWG
        {
//------------------------------------------------------------------------------
        public :
//------------------------------------------------------------------------------

            /**
             * constructor that defines which interpolator is used
             */
            IWG(){};

//------------------------------------------------------------------------------

            /**
             * virtual destructor
             */
            virtual ~IWG(){};

//------------------------------------------------------------------------------

            /**
             * creates the required interpolation matrices as member variables
             */
            virtual void
            create_matrices( Interpolator * aInterpolator ) = 0;

//------------------------------------------------------------------------------

            /**
             * deletes the matrices
             */
            virtual void
            delete_matrices() = 0;

//------------------------------------------------------------------------------

            /**
             * evaluates the Jacobian matrix
             */
            virtual void
            compute_jacobian_and_residual(
                    Mat< real > & aJ,
                    Mat< real > & aR,
                    const Mat< real > & aU,
                    const uint        & aPoint )
            {
                MORIS_ERROR( false, "This function does nothing" );
            }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IWG_HPP_ */

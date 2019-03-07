/*
 * cl_FEM_IWG_Test.hpp
 *
 *  Created on: Feb 11, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_TEST_HPP_
#define SRC_FEM_CL_FEM_IWG_TEST_HPP_

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

        class IWG_Test : public IWG
        {
            // pointer to interpolator
            Field_Interpolator * mFieldInterpolator = nullptr;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             *  r = Nt * u
             */
            IWG_Test( Field_Interpolator * aFieldInterpolator);

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Test(){};

//------------------------------------------------------------------------------
            /**
             * compute the residual
             */
             void compute_residual( Matrix< DDRMat > & aResidual,
                                    uint             & aSpaceDerivativeOrder,
                                    uint             & aTimeDerivativeOrder);

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             */
            void compute_jacobian( Matrix< DDRMat > & aJacobian,
                                   uint             & aSpaceDerivativeOrder,
                                   uint             & aTimeDerivativeOrder);

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             */
             void compute_jacobian_and_residual( Matrix< DDRMat > & aJacobian,
                                                 Matrix< DDRMat > & aResidual,
                                                 uint             & aSpaceDerivativeOrder,
                                                 uint             & aTimeDerivativeOrder);

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_TEST_HPP_ */

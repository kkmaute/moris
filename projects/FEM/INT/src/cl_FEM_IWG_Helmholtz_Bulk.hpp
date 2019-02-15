/*
 * cl_FEM_IWG_Helmoltz_Bulk.hpp
 *
 *  Created on: Feb 13, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_HELMHOLTZ_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_HELMHOLTZ_BULK_HPP_

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

        class IWG_Helmholtz_Bulk : public IWG
        {
            // pointer to interpolator
            Field_Interpolator * mFieldInterpolator = nullptr;

            // Helmholtz filter length parameter
            real mFilterParam ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Helmholtz_Bulk(       Field_Interpolator * aFieldInterpolator,
                                const real                 aFilterParam );

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Helmholtz_Bulk(){};

//------------------------------------------------------------------------------
            /**
             * compute the residual
             * r = kappa * Bxt * gradx(v) + Nt * v - Nt * vHat
             *
             * @param[ in ] aResidual residual vector to fill
             * @param[ in ] aVHat     unfiltered value at evaluation point
             */
            void compute_residual( Matrix< DDRMat > & aResidual,
                                   real               aVHat );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * j = kappa * Bxt * Bx + Nt * N
             *
             * @param[ in ] aJacobian jacobian matrix to fill
             */
            void compute_jacobian( Matrix< DDRMat > & aJacobian );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             *
             * @param[ in ] aJacobian jacobian matrix to fill
             * @param[ in ] aResidual residual vector to fill
             * @param[ in ] aVHat     unfiltered value at evaluation point
             *
             */
            void compute_jacobian_and_residual( Matrix< DDRMat > & aJacobian,
                                                Matrix< DDRMat > & aResidual,
                                                real               aVHat );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_HELMHOLTZ_BULK_HPP_ */

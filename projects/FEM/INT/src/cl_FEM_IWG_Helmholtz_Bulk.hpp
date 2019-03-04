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

            // Helmholtz filter length parameter
            real mFilterParam ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
//            IWG_Helmholtz_Bulk( const real aFilterParam );
            IWG_Helmholtz_Bulk();

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
             * @param[ in ] aResidual            residual vector to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void compute_residual( Matrix< DDRMat >            & aResidual,
                                   Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * j = kappa * Bxt * Bx + Nt * N
             *
             * @param[ in ] aJacobians           list of jacobian matrices to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void compute_jacobian( Cell< Matrix< DDRMat > > & aJacobians,
                                   Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             *
             * @param[ in ] aJacobians           list of jacobian matrices to fill
             * @param[ in ] aResidual            residual vector to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void compute_jacobian_and_residual( Cell< Matrix< DDRMat > >    & aJacobians,
                                                Matrix< DDRMat >            & aResidual,
                                                Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_HELMHOLTZ_BULK_HPP_ */

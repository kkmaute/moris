/*
 * cl_FEM_IWG_LSNormal_Bulk.hpp
 *
 *  Created on: Feb 13, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_LSNORMAL_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_LSNORMAL_BULK_HPP_

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

        class IWG_LSNormal_Bulk : public IWG
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
             IWG_LSNormal_Bulk();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_LSNormal_Bulk(){};

//------------------------------------------------------------------------------
            /**
             * compute the residual
             * r = Nt *( N * n_phiHat - grad( phi ) / norm( grad( phi ) ) )
             *
             * @param[ in ] aResidual            residual vector to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void compute_residual( Matrix< DDRMat >                   & aResidual );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * j = Nt * N
             *
             * @param[ in ] aJacobians           list of jacobian matrices to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             *
             */
            void compute_jacobian( moris::Cell< Matrix< DDRMat > >    & aJacobians );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             *
             * @param[ in ] aJacobian            list of jacobian matrices to fill
             * @param[ in ] aResidual            residual vector to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void compute_jacobian_and_residual( moris::Cell< Matrix< DDRMat > >    & aJacobians,
                                                Matrix< DDRMat >                   & aResidual );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_HELMHOLTZ_BULK_HPP_ */

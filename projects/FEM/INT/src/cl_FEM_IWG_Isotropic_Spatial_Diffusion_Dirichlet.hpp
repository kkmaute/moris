/*
 * cl_FEM_IWG_Isotropic_Spatial_Diffusion_Dirichlet.hpp
 *
 *  Created on: Mar 22, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_SPATIAL_DIFFUSION_DIRICHLET_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_SPATIAL_DIFFUSION_DIRICHLET_HPP_

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

        class IWG_Isotropic_Spatial_Diffusion_Dirichlet : public IWG
        {
            // Nitsche's penalty parameter
            real mGamma;

            // diffusion tensor
            Matrix< DDRMat > mKappa;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Isotropic_Spatial_Diffusion_Dirichlet();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Isotropic_Spatial_Diffusion_Dirichlet(){};

//------------------------------------------------------------------------------
            /**
             * compute the residual
             * r =
             *
             * @param[ in ] aResidual            residual vector to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void
            compute_residual( Matrix< DDRMat >                   & aResidual );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * j =
             *
             * @param[ in ] aJacobians           list of jacobian matrices to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void
            compute_jacobian( moris::Cell< Matrix< DDRMat > >     & aJacobians );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             *
             * @param[ in ] aJacobians           list of jacobian matrices to fill
             * @param[ in ] aResidual            residual vector to fill
             * @param[ in ] aFieldInterpolators  list of active field interpolators
             */
            void
            compute_jacobian_and_residual( moris::Cell< Matrix< DDRMat > >    & aJacobians,
                                           Matrix< DDRMat >                   & aResidual );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_SPATIAL_DIFFUSION_DIRICHLET_HPP_ */

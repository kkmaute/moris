/*
 * cl_FEM_IWG_Isotropic_Spatial_Diffusion_Ghost.hpp
 *
 *  Created on: May 2, 2019
 *      Author: wunsch/noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_SPATIAL_DIFFUSION_GHOST_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_SPATIAL_DIFFUSION_GHOST_HPP_

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

        class IWG_Isotropic_Spatial_Diffusion_Ghost : public IWG
        {
            // Ghost penalty parameter
            real mGammaGhost;

            // mesh parameter describing length of elements
            real mMeshParameter;

            // order of Shape functions
            uint mOrder;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Isotropic_Spatial_Diffusion_Ghost();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Isotropic_Spatial_Diffusion_Ghost(){};

//------------------------------------------------------------------------------
            /**
             * compute the residual
             * @param[ in ] aResidual cell of residual vectors to fill
             */
            void compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
             */
            void compute_jacobian( moris::Cell< Cell< Matrix< DDRMat > > > & aJacobians );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
             * @param[ in ] aResidual  cell of residual vectors to fill
             */
            void compute_jacobian_and_residual( moris::Cell< Cell< Matrix< DDRMat > > > & aJacobians,
                                                moris::Cell< Matrix< DDRMat > >         & aResidual );

//------------------------------------------------------------------------------
            /**
             * method to assemble "normal matrix" from normal vector needed for
             * 2nd and 3rd order Ghost formulations
             * @param[ in ] aOrderGhost Order of derivatives and ghost formulation
             */
            Matrix< DDRMat > get_normal_matrix ( uint aOrderGhost );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_SPATIAL_DIFFUSION_GHOST_HPP_ */

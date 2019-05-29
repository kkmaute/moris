/*
 * cl_FEM_IWG_Isotropic_Spatial_Diffusion_Ghost.hpp
 *
 *  Created on: May 2, 2019
 *      Author: wunsch
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

            // diffusion tensor
            Matrix< DDRMat > mKappa;

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
             * r =
             *
             * @param[ in ] aResidual                residual vector to fill
             * @param[ in ] aLeftFieldInterpolators  list of active left field interpolators
             * @param[ in ] aRightFieldInterpolators list of active right field interpolators
             *
             */
            void
            compute_residual( Matrix< DDRMat >                   & aResidual,
                              moris::Cell< Field_Interpolator* > & aLeftFieldInterpolators,
                              moris::Cell< Field_Interpolator* > & aRightFieldInterpolators );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * j =
             *
             * @param[ in ] aJacobians               list of jacobian matrices to fill
             * @param[ in ] aLeftFieldInterpolators  list of active left field interpolators
             * @param[ in ] aRightFieldInterpolators list of active right field interpolators
             */
            void
            compute_jacobian( Cell< Matrix< DDRMat > >           & aJacobians,
                              moris::Cell< Field_Interpolator* > & aLeftFieldInterpolators,
                              moris::Cell< Field_Interpolator* > & aRightFieldInterpolators );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             *
             * @param[ in ] aJacobians              list of jacobian matrices to fill
             * @param[ in ] aResidual               residual vector to fill
             * @param[ in ] aLeftFieldInterpolators  list of active left field interpolators
             * @param[ in ] aRightFieldInterpolators list of active right field interpolators
             */
            void
            compute_jacobian_and_residual( Cell< Matrix< DDRMat > >           & aJacobians,
                                           Matrix< DDRMat >                   & aResidual,
                                           moris::Cell< Field_Interpolator* > & aLeftFieldInterpolators,
                                           moris::Cell< Field_Interpolator* > & aRightFieldInterpolators );


//------------------------------------------------------------------------------
            /**
             * method to assemble "normal matrix" from normal vector needed for
             * 2nd and 3rd order Ghost formulations
             * @param[ in ] aOrderGhost Order of derivatives and ghost formulation
             */

            Matrix<DDRMat> get_normal_matrix ( uint aOrderGhost );


//------------------------------------------------------------------------------
            /**
             * method to feed order of interpolation functions to the Ghost element
             * @param[ in ] aOrder Order of interpolation functions used on Ghost element
             */
            void set_interpolation_order ( uint aOrder );

//------------------------------------------------------------------------------
            /**
             * method to set the penalty factor of the Ghost element
             * @param[ in ] aGamma penalty factor
             */
            void set_penalty_factor ( real aGamma );


//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_SPATIAL_DIFFUSION_GHOST_HPP_ */

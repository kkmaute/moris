/*
 * cl_FEM_IWG_Isotropic_Struc_Linear_Neumann.hpp
 *
 *  Created on: Okt 06, 2019
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_NEUMANN_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_NEUMANN_HPP_

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

        class IWG_Isotropic_Struc_Linear_Neumann : public IWG
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Isotropic_Struc_Linear_Neumann(){};

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Isotropic_Struc_Linear_Neumann(){};

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
            void compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
             * @param[ in ] aResidual  cell of residual vectors to fill
             */
            void compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                moris::Cell< Matrix< DDRMat > >                & aResidual );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_NEUMANN_HPP_ */

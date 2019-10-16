/*
 * cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.hpp
 *
 *  Created on: Okt 06, 2019
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_BULK_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class IWG_Isotropic_Struc_Linear_Bulk : public IWG
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
            IWG_Isotropic_Struc_Linear_Bulk(){};

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Isotropic_Struc_Linear_Bulk(){};

//------------------------------------------------------------------------------
            /**
             * compute the residual
             * r =
             * @param[ in ] aResidual residual vector to fill
             */
            void compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual );

//------------------------------------------------------------------------------
            /**
             * compute the jacobian
             * j =
             * @param[ in ] aJacobians list of jacobian matrices to fill
             */
            void compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians );

//------------------------------------------------------------------------------
            /**
             * compute the residual and the jacobian
             * @param[ in ] aJacobians list of jacobian matrices to fill
             * @param[ in ] aResidual  residual vector to fill
             */
            void compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                moris::Cell< Matrix< DDRMat > >                & aResidual );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_BULK_HPP_ */

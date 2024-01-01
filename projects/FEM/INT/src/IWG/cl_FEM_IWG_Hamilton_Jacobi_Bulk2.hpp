/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Hamilton_Jacobi_Bulk2.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK2_HPP_
#define SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK2_HPP_

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                      //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Hamilton_Jacobi_Bulk2 : public IWG
        {

                //------------------------------------------------------------------------------
            public:
                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Hamilton_Jacobi_Bulk2();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Hamilton_Jacobi_Bulk2(){};

                //------------------------------------------------------------------------------
                /**
                 * compute the residual
                 * @param[ in ] aResidual cell of residual vector to fill
                 */
                void compute_residual( real tWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the jacobian
                 * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
                 */
                void compute_jacobian( real tWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the residual and the jacobian
                 * @param[ in ] aJacobians cell of cell of jacobian matrices to fill
                 * @param[ in ] aResidual  cell of residual vectors to fill
                 */
                void compute_jacobian_and_residual( Vector< Vector< Matrix< DDRMat > > > & aJacobians,
                        Vector< Matrix< DDRMat > >                & aResidual );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK_HPP_ */


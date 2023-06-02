/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Hamilton_Jacobi_Bulk.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Hamilton_Jacobi_Bulk : public IWG
        {

                //------------------------------------------------------------------------------
            public:
                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Hamilton_Jacobi_Bulk();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Hamilton_Jacobi_Bulk(){};

                //------------------------------------------------------------------------------
                /**
                 * compute the residual
                 * @param[ in ] aResidual cell of residual vectors to fill
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
                void compute_jacobian_and_residual(
                        moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                        moris::Cell< Matrix< DDRMat > >                & aResidual );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_HAMILTON_JACOBI_BULK_HPP_ */


/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Helmholtz_Bulk2.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_HELMHOLTZ_BULK2_HPP_
#define SRC_FEM_CL_FEM_IWG_HELMHOLTZ_BULK2_HPP_

#include "moris_typedefs.hpp"                     //MRS/COR/src
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

        class IWG_Helmholtz_Bulk2 : public IWG
        {

            // Helmholtz filter length parameter
            real mFilterParam;

            // sharpness parameter for smoothed Heaviside function
            real mSharpParam;

            // element size
            real mHe;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /*
             *  constructor
             */
//            IWG_Helmholtz_Bulk( const real aFilterParam );
            IWG_Helmholtz_Bulk2();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Helmholtz_Bulk2(){};

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
            void compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                moris::Cell< Matrix< DDRMat > >                & aResidual );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_HELMHOLTZ_BULK_HPP_ */


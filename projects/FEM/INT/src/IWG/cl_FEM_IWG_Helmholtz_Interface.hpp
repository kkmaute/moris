/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Helmholtz_Interface.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_HELMHOLTZ_INTERFACE_HPP_
#define SRC_FEM_CL_FEM_IWG_HELMHOLTZ_INTERFACE_HPP_

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

        class IWG_Helmholtz_Interface : public IWG
        {
                // Helmholtz filter length parameter
                real mFilterParam ;

                //------------------------------------------------------------------------------
            public:
                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                //            IWG_Helmholtz_Interface( const real aFilterParam );
                IWG_Helmholtz_Interface();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Helmholtz_Interface(){};

                //------------------------------------------------------------------------------
                /**
                 * compute the residual
                 * @param[ in ] aResidual        residual vector to fill
                 */
                void compute_residual( real tWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the jacobian
                 * @param[ in ] aJacobians           list of jacobian matrices to fill
                 */
                void compute_jacobian( real tWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the residual and the jacobian
                 * @param[ in ] aJacobians           list of jacobian matrices to fill
                 * @param[ in ] aResidual            residual vector to fill
                 */
                void compute_jacobian_and_residual(
                        moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                        moris::Cell< Matrix< DDRMat > >                & aResidual );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_HELMHOLTZ_INTERFACE_HPP_ */


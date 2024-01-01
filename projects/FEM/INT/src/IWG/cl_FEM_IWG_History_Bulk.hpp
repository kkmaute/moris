/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_History_Bulk.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_HISTORY_BULK_HPP_
#define SRC_FEM_CL_FEM_IWG_HISTORY_BULK_HPP_

#include <map>
#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_History_Bulk : public IWG
        {

            private:
              //------------------------------------------------------------------------------

              enum class IWG_Property_Type
              {
                  WEIGHT,
                  LUMP,
                  MAX_ENUM
              };

              enum class IWG_Constitutive_Type
              {
                  ELASTIC_DAMAGE,
                  MAX_ENUM
              };

              enum class IWG_Stabilization_Type
              {
                  MAX_ENUM
              };

            public:
              //------------------------------------------------------------------------------
              /*
               *  constructor
               */
              IWG_History_Bulk();

              //------------------------------------------------------------------------------
              /**
               * trivial destructor
               */
              ~IWG_History_Bulk(){};

              //------------------------------------------------------------------------------
              /**
               * compute the residual
               * @param[ in ] aWStar weight associated to the evaluation point
               */
              void compute_residual( real tWStar );

              //------------------------------------------------------------------------------
              /**
               * compute the jacobian
               * @param[ in ] aWStar weight associated to the evaluation point
               */
              void compute_jacobian( real tWStar );

              //------------------------------------------------------------------------------
              /**
               * compute the residual and the jacobian
               * @param[ in ] aWStar weight associated to the evaluation point
               */
              void compute_jacobian_and_residual( real aWStar );

              //------------------------------------------------------------------------------
              /**
               * compute the derivative of the residual wrt design variables
               * @param[ in ] aWStar weight associated to the evaluation point
               */
              void compute_dRdp( real aWStar );

              //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_HISTORY_BULK_HPP_ */

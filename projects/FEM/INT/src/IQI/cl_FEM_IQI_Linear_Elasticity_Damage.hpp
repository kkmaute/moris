/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Linear_Elasticity_Damage.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_LINEAR_ELASTICITY_DAMAGE_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_LINEAR_ELASTICITY_DAMAGE_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IQI.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        /*
         * Output the coefficient from the linear elasticity damage constitutive model
         * 0 - equivalent strain
         * 1 - damage
         */
        class IQI_Linear_Elasticity_Damage : public IQI
        {
            private:
              enum class IQI_Constitutive_Type
              {
                  ELASTIC_DAMAGE,
                  MAX_ENUM
              };

              //------------------------------------------------------------------------------

            public:
                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
              IQI_Linear_Elasticity_Damage();

              //------------------------------------------------------------------------------
              /**
               * trivial destructor
               */
              ~IQI_Linear_Elasticity_Damage(){};

            private:

                //------------------------------------------------------------------------------
                /**
                 * compute the quantity of interest,
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_QI( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * Evaluate the quantity of interest and fill aQI with value
                 * @param[ in ] aQI IQI value at evaluation point
                 */
                void compute_QI( Matrix< DDRMat > & aQI );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the quantity of interest wrt dof types
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_dQIdu( real aWStar )
                {
                    MORIS_ERROR( false, "IQI_Linear_Elasticity_Damage::compute_dQIdu - not implemented." );
                }

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the quantity of interest wrt dof types
                 * @param[ in ] aDofType group of dof types wrt which derivatives are evaluated
                 * @param[ in ] adQIdu   derivative of quantity of interest matrix to fill
                 */
                void compute_dQIdu(
                        moris::Cell< MSI::Dof_Type > & aDofType,
                        Matrix< DDRMat >             & adQIdu )
                {
                    MORIS_ERROR( false, "IQI_Linear_Elasticity_Damage::compute_dQIdu() - not implemented for a drag/lift coefficient IQI." );
                }
        };
    }/* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_LINEAR_ELASTICITY_DAMAGE_HPP_ */

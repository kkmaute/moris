/*
 * cl_FEM_IQI_Thermal_Energy.hpp
 *
 *  Created on: Sept 27, 2020
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_THERMAL_ENERGY_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_THERMAL_ENERGY_HPP_

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

        class IQI_Thermal_Energy : public IQI
        {
                //------------------------------------------------------------------------------

                enum class IQI_Constitutive_Type
                {
                        FLUID,
                        MAX_ENUM
                };

            public:
                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                IQI_Thermal_Energy();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IQI_Thermal_Energy(){};

                //------------------------------------------------------------------------------

            private:

                //------------------------------------------------------------------------------
                /**
                 * compute the quantity of interest
                 * @param[ in ] aQI quantity of interest matrix to fill
                 */
                void compute_QI( Matrix< DDRMat > & aQI );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the quantity of interest wrt dof types
                 * @param[ in ] aDofType group of dof types wrt which derivatives are evaluated
                 * @param[ in ] adQIdu   derivative of quantity of interest matrix to fill
                 */
                void compute_dQIdu(
                        moris::Cell< MSI::Dof_Type > & aDofType,
                        Matrix< DDRMat >             & adQIdu );

                //------------------------------------------------------------------------------
        };
    }/* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_THERMAL_ENERGY_HPP_ */

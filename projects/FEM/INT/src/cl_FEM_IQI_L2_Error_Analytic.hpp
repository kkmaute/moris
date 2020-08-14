/*
 * cl_FEM_IQI_L2_Error_Analytic.hpp
 *
 *  Created on: Feb 3, 2020
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_L2_ERROR_ANALYTIC_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_L2_ERROR_ANALYTIC_HPP_

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

        class IQI_L2_Error_Analytic : public IQI
        {
                //------------------------------------------------------------------------------

                enum class IQI_Property_Type
                {
                        L2_CHECK,
                        MAX_ENUM
                };

                // Local string to property enum map
                std::map< std::string, IQI_Property_Type > mPropertyMap;

            public:
                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                IQI_L2_Error_Analytic();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IQI_L2_Error_Analytic(){};

                //------------------------------------------------------------------------------
                /**
                 * set property
                 * @param[ in ] aProperty       a property pointer
                 * @param[ in ] aPropertyString a string defining the property
                 * @param[ in ] aIsMaster       an enum for master or slave
                 */
                void set_property(
                        std::shared_ptr< Property > aProperty,
                        std::string                 aPropertyString,
                        mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * compute the quantity of interest
                 * @param[ in ] aQI quantity of interest matrix to fill
                 */
                void compute_QI( Matrix< DDRMat > & aQI );

                //------------------------------------------------------------------------------
                /**
                 * compute the quantity of interest
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_QI( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the quantity of interest wrt dof types
                 * @param[ in ] adQIdDof derivative of quantity of interest matrix to fill
                 */
                void compute_dQIdu( Matrix< DDRMat > & adQIdDof );

                //------------------------------------------------------------------------------
        };
    }/* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRAIN_ENERGY_HPP_ */

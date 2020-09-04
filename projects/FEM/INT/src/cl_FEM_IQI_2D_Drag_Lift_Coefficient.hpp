/*
 * cl_FEM_IQI_2D_Drag_Lift_Coefficient.hpp
 *
 *  Created on: May 11, 2020
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_DRAG_LIFT_COEFFICIENT_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_DRAG_LIFT_COEFFICIENT_HPP_

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

        class IQI_Drag_Lift_Coefficient : public IQI
        {
                //------------------------------------------------------------------------------
            public:

                // property type for the IQI
                enum class Property_Type
                {
                    DENSITY,    // fluid density
                    VISCOSITY,
                    VELOCITY_MAX,
                    DIAMETER,
                    MAX_ENUM
                };

                // local string to property enum map
                std::map< std::string, Property_Type > mPropertyMap;

                // sign for drag/lift
                sint mBeta;

                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IQI_Drag_Lift_Coefficient( sint aBeta );

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IQI_Drag_Lift_Coefficient(){};

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

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_DRAG_LIFT_HPP_ */

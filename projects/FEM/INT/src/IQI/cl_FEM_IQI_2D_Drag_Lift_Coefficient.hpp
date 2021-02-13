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
                    MORIS_ERROR( false, "IQI_Drag_Lift_Coefficient::compute_dQIdu - not implemented." );
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
                    MORIS_ERROR( false, "IQI_Drag_Lift_Coefficient::compute_dQIdu() - not implemented for a drag/lift coefficient IQI.");
                }

                //------------------------------------------------------------------------------
        };
    }/* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_DRAG_LIFT_HPP_ */
/*
 * cl_FEM_IQI_J_Integral.hpp
 *
 *  Created on: Feb 11, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_J_INTEGRAL_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_J_INTEGRAL_HPP_

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_FEM_IQI.hpp"


namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IQI_J_Integral : public IQI
        {

                enum class IQI_Constitutive_Type
                {
                        ELAST_LIN_ISO,
                        MAX_ENUM
                };

            public :
                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                IQI_J_Integral();

                //------------------------------------------------------------------------------
                /*
                 * destructor
                 */
                ~IQI_J_Integral(){};

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
        //------------------------------------------------------------------------------
    }   // end fem namespace
}       // end moris namespace




#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_J_INTEGRAL_HPP_ */

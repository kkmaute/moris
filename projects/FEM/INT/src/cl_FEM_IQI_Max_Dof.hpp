/*
 * cl_FEM_IQI_Max_Dof.hpp
 *
 *  Created on: Jul 10, 2020
 *      Author: wunsch
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_MAX_DOF_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_MAX_DOF_HPP_

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

        class IQI_Max_Dof : public IQI
        {

                //------------------------------------------------------------------------------

            public:
                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                IQI_Max_Dof();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IQI_Max_Dof(){};

            private:
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

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_MAX_DOF_HPP_ */

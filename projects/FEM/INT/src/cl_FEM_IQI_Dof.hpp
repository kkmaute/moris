/*
 * cl_FEM_IQI_Dof.hpp
 *
 *  Created on: Jan 23, 2020
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_DOF_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_DOF_HPP_

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

        class IQI_Dof : public IQI
        {

        public:
            //------------------------------------------------------------------------------
            /*
             * constructor
             */
            IQI_Dof();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IQI_Dof(){};

            //------------------------------------------------------------------------------

        private:
            /**
             * compute the quantity of interest
             * @param[ in ] aQI quantity of interest matrix to fill
             */
            void compute_QI( Matrix< DDRMat > & aQI );

            //------------------------------------------------------------------------------
            /**
             * compute the derivative of the quantity of interest wrt dof types
             * @param[ in ] adQIdu derivative of quantity of interest matrix to fill
             */
            void compute_dQIdu( MSI::Dof_Type aDofType, Matrix< DDRMat > & adQIdu );

            //------------------------------------------------------------------------------
        };
    }/* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRAIN_ENERGY_HPP_ */

/*
 * cl_FEM_IQI_H1_Semi_Error.hpp
 *
 *  Created on: Feb 3, 2020
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_H1_SEMI_ERROR_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_H1_SEMI_ERROR_HPP_

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

        class IQI_H1_Semi_Error : public IQI
        {
//------------------------------------------------------------------------------

            enum class IQI_Property_Type
            {
                MAX_ENUM
            };

            // Local string to property enum map
            std::map< std::string, IQI_Property_Type > mPropertyMap;

            enum class IQI_Constitutive_Type
            {
                MAX_ENUM
            };

            // Local string to constitutive enum map
            std::map< std::string, IQI_Constitutive_Type > mConstitutiveMap;

            enum class IQI_Stabilization_Type
            {
                MAX_ENUM
            };

            // Local string to constitutive enum map
            std::map< std::string, IQI_Stabilization_Type > mStabilizationMap;

        public:
//------------------------------------------------------------------------------
            /*
             * constructor
             */
            IQI_H1_Semi_Error();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IQI_H1_Semi_Error(){};

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

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_H1_SEMI_ERROR_HPP_ */

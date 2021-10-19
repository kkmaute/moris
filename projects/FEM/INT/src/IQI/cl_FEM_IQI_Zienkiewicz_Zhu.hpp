/*
 * cl_FEM_IQI_Zienkiewicz_Zhu.hpp
 *
 *  Created on: Feb, 2020
 *      Author: Schmidt
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_ZIENKIEWICZ_ZHU_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_ZIENKIEWICZ_ZHU_HPP_

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

        class IQI_Zienkiewicz_Zhu : public IQI
        {

                //------------------------------------------------------------------------------

            private:

                //------------------------------------------------------------------------------

                // stress type to evaluate
                enum Stress_Type mStressType;

                enum class IQI_Constitutive_Type
                {
                        ELAST_LIN_ISO,
                        MAX_ENUM
                };

                //------------------------------------------------------------------------------

            protected:

                //------------------------------------------------------------------------------
                /*
                 * default constructor
                 */
                IQI_Zienkiewicz_Zhu()
                {
                    MORIS_ERROR( false, " IQI_Zienkiewicz_Zhu: Default constructor unavailable. Use enum to construct IQI. ");
                };

                //------------------------------------------------------------------------------

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                IQI_Zienkiewicz_Zhu( enum Stress_Type aStressType );

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IQI_Zienkiewicz_Zhu(){};

                //------------------------------------------------------------------------------

            private:

                //------------------------------------------------------------------------------
                /**
                 * gets the stress vector from the constitutive model and sorts it into a format standardized
                 * (independent of 2D-plane strain, 2D-plain stress, and 3D)
                 * @param[ out ] aStressVector vector with 6 entries containing the normal and shear stress values
                 */
                void get_stress_vector( Matrix< DDRMat > & aStressVector );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the Von-Mises stress using the stress vector provided by the constitutive model
                 * @param[ out ] tStressValue the value of the Von-Mises stress
                 */
                real eval_Von_Mises_stress();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the i-th principal stress using the stress vector provided by the constitutive model
                 * @param[ in ] aPrincipalStressIndex number i of the principal stress to get,
                 *               i.e. 1,2,3 for the first, second and third principal stress, respectively
                 * @param[ out ] tStressValue the value of the requested principal stress
                 */
                real eval_principal_stress( uint aPrincipalStressIndex );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the i-th normal stress using the stress vector provided by the constitutive model
                 * @param[ in ] aStressIndex number of the principal stress to get,
                 *               i.e. 1,2,3 for the x-, y-, and z-normal stresses, respectively
                 * @param[ out ] tStressValue the value of the requested principal stress
                 */
                real eval_normal_stress( uint aStressIndex );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the i-th shear stress using the stress vector provided by the constitutive model
                 * @param[ in ] aStressIndex number of the principal stress to get,
                 *               i.e. 1,2,3 for the yz-, xz-, and xy-shear stresses, respectively
                 *               in 2D the input is ignored, as there's only one shear stress
                 * @param[ out ] tStressValue the value of the requested principal stress
                 */
                real eval_shear_stress( uint aStressIndex );

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
                    MORIS_ERROR( false, "IQI_Zienkiewicz_Zhu::compute_dQIdu - not implemented." );
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
                    MORIS_ERROR( false, "IQI_Zienkiewicz_Zhu::compute_dQIdu() - not implemented for a drag/lift coefficient IQI.");
                }

                //------------------------------------------------------------------------------

        };
    }/* end namespace fem */
} /* end namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_ZIENKIEWICZ_ZHU_HPP_ */

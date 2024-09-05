/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Stress.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRESS_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRESS_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                          //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_IQI.hpp"                   //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class IQI_Stress : public IQI
    {

        //------------------------------------------------------------------------------

      private:
        //------------------------------------------------------------------------------

        // stress type to evaluate
        enum Stress_Type mStressType;

        // flux type to evaluate
        enum CM_Function_Type mFluxType;

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
        IQI_Stress()
        {
            MORIS_ERROR( false, " IQI_Stress: Default constructor unavailable. Use enum to construct IQI. " );
        };

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         * @param[ in ] aStressType stress type to evaluate stress in IQI
         * @param[ in ] aFluxType   flux type to evaluate for the CM to evaluate stress in IQI
         */
        IQI_Stress(
                enum Stress_Type      aStressType,
                enum CM_Function_Type aFluxType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~IQI_Stress() override{};

        //------------------------------------------------------------------------------

      private:
        //------------------------------------------------------------------------------
        /**
         * gets the stress vector from the constitutive model and sorts it into a format standardized
         * (independent of 2D-plane strain, 2D-plain stress, and 3D)
         * @param[ out ] aStressVector vector with 6 entries containing the normal and shear stress values
         */
        void get_stress_vector( Matrix< DDRMat > &aStressVector );

        //------------------------------------------------------------------------------
        /**
         * evaluate the Von-Mises stress using the stress vector provided by the constitutive model
         * @param[ out ] aStressVector the value of the Von-Mises stress
         */
        void eval_Von_Mises_stress( Matrix< DDRMat > &aStressVector );

        //------------------------------------------------------------------------------
        /**
         * evaluate the i-th principal stress using the stress vector provided by the constitutive model
         * @param[ in ] aPrincipalStressIndex number i of the principal stress to get,
         *               i.e. 1,2,3 for the first, second and third principal stress, respectively
         * @param[ out ] tStressValue the value of the requested principal stress
         */
        void eval_principal_stress( uint aPrincipalStressIndex, Matrix< DDRMat > &aStressVector );

        //------------------------------------------------------------------------------
        /**
         * evaluate the i-th normal stress using the stress vector provided by the constitutive model
         * @param[ in ] aStressIndex number of the principal stress to get,
         *               i.e. 1,2,3 for the x-, y-, and z-normal stresses, respectively
         * @param[ out ] tStressValue the value of the requested principal stress
         */
        void eval_normal_stress( uint aStressIndex, Matrix< DDRMat > &aStressVector );

        //------------------------------------------------------------------------------
        /**
         * evaluate the i-th shear stress using the stress vector provided by the constitutive model
         * @param[ in ] aStressIndex number of the principal stress to get,
         *               i.e. 1,2,3 for the yz-, xz-, and xy-shear stresses, respectively
         *               in 2D the input is ignored, as there's only one shear stress
         * @param[ out ] tStressValue the value of the requested principal stress
         */
        void eval_shear_stress( uint aStressIndex, Matrix< DDRMat > &aStressVector );

        //------------------------------------------------------------------------------
        /**
         * compute the quantity of interest
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_QI( real aWStar ) override;

        //------------------------------------------------------------------------------
        /**
         * Evaluate the quantity of interest and fill aQI with value
         * @param[ in ] aQI IQI value at evaluation point
         */
        void compute_QI( Matrix< DDRMat > &aQI ) override;

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        void compute_dQIdu( real aWStar ) override
        {
            MORIS_ERROR( false, "IQI_Stress::compute_dQIdu - not implemented." );
        }

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the quantity of interest wrt dof types
         * @param[ in ] aDofType group of dof types wrt which derivatives are evaluated
         * @param[ in ] adQIdu   derivative of quantity of interest matrix to fill
         */
        void compute_dQIdu(
                Vector< MSI::Dof_Type > &aDofType,
                Matrix< DDRMat >        &adQIdu ) override
        {
            MORIS_ERROR( false, "IQI_Stress::compute_dQIdu() - not implemented for a drag/lift coefficient IQI." );
        }

        //------------------------------------------------------------------------------
        /**
         * compute matrix dimension of the IQI
         * @param[ out ] space dimension of the IQI
         */
        std::pair< uint, uint > get_matrix_dim() override;
    };
}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_IQI_STRESS_HPP_ */


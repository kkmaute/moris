/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Diffusion_Linear_Isotropic.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_HPP_
#define SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_HPP_

#include <iostream>
#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class CM_Diffusion_Linear_Isotropic : public Constitutive_Model
    {
        //------------------------------------------------------------------------------

      protected:
        // default local properties
        std::shared_ptr< Property > mPropConductivity = nullptr;
        std::shared_ptr< Property > mPropHeatCapacity = nullptr;
        std::shared_ptr< Property > mPropDensity      = nullptr;
        std::shared_ptr< Property > mPropEigenStrain  = nullptr;

      private:
        // Default dof type for CM
        MSI::Dof_Type mTempDof  = MSI::Dof_Type::TEMP;
        MSI::Dof_Type mThetaDof = MSI::Dof_Type::UNDEFINED;

        // property type for CM
        enum class CM_Property_Type
        {
            CONDUCTIVITY,
            HEAT_CAPACITY,
            DENSITY,
            EIGEN_STRAIN,
            MAX_ENUM
        };

        //------------------------------------------------------------------------------

      public:
        /*
         * trivial constructor
         */
        CM_Diffusion_Linear_Isotropic();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~CM_Diffusion_Linear_Isotropic() override {};

        //------------------------------------------------------------------------------
        /*
         * @return constitutive_type
         */
        Constitutive_Type
        get_constitutive_type() const override
        {
            return Constitutive_Type::DIFF_LIN_ISO;
        }

        //------------------------------------------------------------------------------
        /**
         * set constitutive model dof types
         * @param[ in ] aDofTypes a list of group of dof types
         * @param[ in ] aDofStrings a list of strings to describe the dof types
         */
        void set_dof_type_list(
                const Vector< Vector< MSI::Dof_Type > > &aDofTypes,
                const Vector< std::string >             &aDofStrings ) override;

        //------------------------------------------------------------------------------
        /**
         * set constitutive model dv types
         * @param[ in ] aDvTypes a list of group of dv types
         * @param[ in ] aDvStrings a list of strings to describe the dv types
         */
        void set_dv_type_list(
                const Vector< gen::PDV_Type > &aDvTypes,
                const Vector< std::string >   &aDvStrings ) override
        {
            Constitutive_Model::set_dv_type_list( aDvTypes );
        }

        //------------------------------------------------------------------------------
        /**
         * set local properties
         */
        void set_local_properties() override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model flux
         */
        void eval_flux() override;

        //------------------------------------------------------------------------------
        /**
         * evaluates the constitutive model enthalpy
         */
        void eval_Energy() override;

        //------------------------------------------------------------------------------
        /**
         * evaluates the constitutive model change rate of enthalpy
         */
        void eval_EnergyDot() override;

        //------------------------------------------------------------------------------
        /**
         * evaluates the constitutive model spatial gradient of enthalpy
         */
        void eval_gradEnergy() override;

        //------------------------------------------------------------------------------
        /**
         * evaluates the constitutive model change rate of spatial gradient of enthalpy (needed for GGLS-stabilization)
         */
        void eval_gradEnergyDot() override;

        //------------------------------------------------------------------------------
        /**
         * evaluates the gradient of the divergence of the flux (needed for GGLS-stabilization)
         */
        void eval_graddivflux() override;
        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model test flux
         * flux ( mSpaceDim x 1)
         */
        void eval_testFlux();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the divergence of the flux
         */
        void eval_divflux() override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model traction
         * @param[ in ] aNormal normal
         * traction ( 1 x 1 )
         */
        void eval_traction( const Matrix< DDRMat > &aNormal ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model test traction
         * @param[ in ] aNormal normal
         * test traction ( numDof x 1 )
         */
        void eval_testTraction(
                const Matrix< DDRMat >        &aNormal,
                const Vector< MSI::Dof_Type > &aTestDofType ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model strain
         * strain ( mSpaceDim x 1 )
         */
        void eval_strain() override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the divergence of the strain
         */
        void eval_divstrain() override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model test strain
         * test strain ( mSpaceDim x numDof  )
         */
        void eval_testStrain() override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model matrix
         * constitutive matrix ( mSpaceDim x mSpaceDim )
         */
        void eval_const() override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model flux derivative wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * dFluxdDOF ( mSpaceDim x numDerDof )
         */
        void eval_dFluxdDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model enthalpy wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * dEnergydDOF ( 1 x numDerDof )
         */
        void eval_dEnergydDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model enthalpy change rate wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * dEnergyDotdDOF ( 1 x numDerDof )
         */
        void eval_dEnergyDotdDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model gradient of enthalpy wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * dGradEnergydDOF ( mSpaceDim x numDerDof )
         */
        void eval_dGradEnergydDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model gradient of enthalpy change rate wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * dgradEnergyDotdDOF ( mSpaceDim x numDerDof )
         */
        void eval_dGradEnergyDotdDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model gradient of divergence of flux wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * dGradDivFluxdDOF ( mSpaceDim x numDerDof )
         */
        void eval_dGradDivFluxdDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the divergence of the flux wrt dof type
         */
        void eval_ddivfluxdu( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model traction derivative wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aNormal   normal
         * dTractiondDOF ( 1 x numDerDof )
         */
        void eval_dTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >        &aNormal ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model test traction derivative wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aNormal   normal
         * dTestTractiondDOF ( numDof x numDerDof )
         */
        void eval_dTestTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >        &aNormal,
                const Vector< MSI::Dof_Type > &aTestDofTypes ) override;

        /**
         * evaluate the constitutive model test traction derivative wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aNormal   normal
         * dTestTractiondDOF ( numDof x numDerDof )
         */
        void eval_dTestTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >        &aNormal,
                const Matrix< DDRMat >        &aJump,
                const Vector< MSI::Dof_Type > &aTestDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model strain derivative wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * dStraindDOF ( mSpaceDim x numDerDof )
         */
        void eval_dStraindDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the divergence of the strain wrt dof type
         */
        void eval_ddivstraindu( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the divergence of the strain wrt dof type
         */
        void eval_dConstdDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model flux derivative wrt to a dv type
         * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
         * dFluxdDV ( mSpaceDim x numDerDv )
         */
        void eval_dFluxdDV( const Vector< gen::PDV_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model strain derivative wrt to a dv type
         * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
         * dStraindDV ( mSpaceDim x numDerDV )
         */
        void eval_dStraindDV( const Vector< gen::PDV_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_HPP_ */

/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_TURBULENCE_HPP_
#define SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_TURBULENCE_HPP_

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class CM_Diffusion_Linear_Isotropic_Turbulence : public CM_Diffusion_Linear_Isotropic
    {
        //------------------------------------------------------------------------------

      protected:
        // property type for CM
        enum class CM_Property_Type
        {
            CONDUCTIVITY,
            HEAT_CAPACITY,
            DENSITY,
            EIGEN_STRAIN,
            TURBULENT_PRANDTL,
            MAX_ENUM
        };

        // default local properties
        std::shared_ptr< Property > mPropPrandtlT = nullptr;

        // flags for turbulence dynamic viscosity related evaluation
        bool                    mTurbDynViscEval = true;
        moris::Matrix< DDBMat > mdTurbDynViscduEval;
        moris::Matrix< DDBMat > mdTurbDynViscdxEval;
        moris::Matrix< DDBMat > mdTurbDynViscdxduEval;
        moris::Matrix< DDBMat > mTestTurbDynViscEval;
        moris::Matrix< DDBMat > mdTestTurbDynViscduEval;

        // storage for turbulence dynamic viscosity related evaluation
        Matrix< DDRMat >                     mTurbDynVisc;
        Vector< Matrix< DDRMat > >           mdTurbDynViscdu;
        Vector< Matrix< DDRMat > >           mdTurbDynViscdx;
        Vector< Vector< Matrix< DDRMat > > > mdTurbDynViscdxdu;
        Vector< Matrix< DDRMat > >           mTestTurbDynVisc;
        Vector< Vector< Matrix< DDRMat > > > mdTestTurbDynViscdu;

        // flags for effective conductivity related evaluation
        bool                    mEffCondEval = true;
        moris::Matrix< DDBMat > mdEffCondduEval;
        moris::Matrix< DDBMat > mdEffConddxEval;
        moris::Matrix< DDBMat > mdEffConddxduEval;
        moris::Matrix< DDBMat > mTestEffCondEval;
        moris::Matrix< DDBMat > mdTestEffCondduEval;

        // storage for effective conductivity related evaluation
        Matrix< DDRMat >                     mEffCond;
        Vector< Matrix< DDRMat > >           mdEffConddu;
        Vector< Matrix< DDRMat > >           mdEffConddx;
        Vector< Vector< Matrix< DDRMat > > > mdEffConddxdu;
        Vector< Matrix< DDRMat > >           mTestEffCond;
        Vector< Vector< Matrix< DDRMat > > > mdTestEffConddu;

        //------------------------------------------------------------------------------

      public:
        /*
         * trivial constructor
         */
        CM_Diffusion_Linear_Isotropic_Turbulence();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~CM_Diffusion_Linear_Isotropic_Turbulence() override {};

        //------------------------------------------------------------------------------
        /*
         * @return constitutive_type
         */
        Constitutive_Type get_constitutive_type() const override
        {
            return Constitutive_Type::DIFF_LIN_ISO_TURBULENCE;
        }

        //------------------------------------------------------------------------------
        /*
         * reset specific evaluation flag
         * (child implementation)
         */
        void reset_specific_eval_flags() override;

        //------------------------------------------------------------------------------
        /*
         * initialize specific storage and evaluation flag
         * (child implementation)
         */
        void initialize_spec_storage_vars_and_eval_flags() override;

        //------------------------------------------------------------------------------
        /**
         * set constitutive model dv types
         * @param[ in ] aDvTypes a list of group of dv types
         * @param[ in ] aDvStrings a list of strings to describe the dv types
         */
        void set_dv_type_list(
                const Vector< Vector< gen::PDV_Type > > &aDvTypes,
                const Vector< std::string >             &aDvStrings ) override
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
         * evaluate the flux derivative wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * dFluxdDOF ( mSpaceDim x numDerDof )
         */
        void eval_dFluxdDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the divergence of the flux
         */
        void eval_divflux() override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the divergence of the flux wrt dof type
         */
        void eval_ddivfluxdu( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model traction
         * @param[ in ] aNormal normal
         * traction ( 1 x 1 )
         */
        void eval_traction( const Matrix< DDRMat > &aNormal ) override;

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
         * evaluate the constitutive model test traction
         * @param[ in ] aNormal normal
         * test traction ( numDof x 1 )
         */
        void eval_testTraction(
                const Matrix< DDRMat >        &aNormal,
                const Vector< MSI::Dof_Type > &aTestDofType ) override;

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
         * evaluate the constitutive model matrix
         * constitutive matrix ( mSpaceDim x mSpaceDim )
         */
        void eval_const() override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the divergence of the strain wrt dof type
         */
        void eval_dConstdDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the the effective conductivity k_eff = k + k_t
         * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
         *               if there are several
         * @param[ out ] mEffCond effective conductivity
         */
        const Matrix< DDRMat > &effective_conductivity(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the effective conductivity k_eff = k + k_t
         */
        void eval_effective_conductivity();

        //------------------------------------------------------------------------------
        /**
         * get the derivative of the effective conductivity wrt space
         * @param[ in ] aOrder order of the derivative
         * @param[ in ] aCMFunctionType enum for specific type of which effective conductivity
         * @param[ out ] mdeffconddx derivative of the effective conductivity wrt space
         */
        const Matrix< DDRMat > &deffconddx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the effective conductivity wrt space
         * @param[ in ] aOrder order of the derivative
         */
        void eval_deffconddx( uint aOrder );

        //------------------------------------------------------------------------------
        /**
         * get the derivative of the effective conductivity wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which effective conductivity
         * @param[ out ] mdeffconddu derivative of the effective conductivity wrt dof types
         */
        const Matrix< DDRMat > &deffconddu(
                const Vector< MSI::Dof_Type > &aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the effective conductivity derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_deffconddu( const Vector< MSI::Dof_Type > &aDofTypes );

        //------------------------------------------------------------------------------
        /**
         * evaluate the effective conductivity derivative wrt to space and dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         */
        void eval_deffconddxdu( const Vector< MSI::Dof_Type > &aDofTypes,
                uint                                           aOrder );

        /**
         * get the derivative of the effective conductivity wrt space and dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         * @param[ in ] aCMFunctionType enum for specific type of which effective conductivity
         * @param[ out ] mdeffconddxdu derivative of the effective conductivity wrt dof types
         */
        const Matrix< DDRMat > &deffconddxdu(
                const Vector< MSI::Dof_Type > &aDofType,
                uint                           aOrder,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get the test of the effective conductivity wrt test dof type
         * @param[ in ] aTestDofTypes a dof type wrt which the test is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which turbulent dynamic viscosity
         * @param[ out ] mdTestTurbDynViscdu test of the effective conductivity wrt dof types
         */
        const Matrix< DDRMat > &
        testeffcond(
                const Vector< MSI::Dof_Type > &aTestDofTypes,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the test effective conductivity wrt to test dof type
         * @param[ in ] aTestDofTypes a dof type wrt which the test is evaluated
         */
        virtual void eval_testeffcond(
                const Vector< MSI::Dof_Type > &aTestDofTypes );
        
        //------------------------------------------------------------------------------
        /**
         * get the derivative of the test effective conductivity wrt dof
         * @param[ in ]  aDofTypes group of dof type for derivative
         * @param[ in ]  aTestDofTypes group of dof type for test
         * @param[ out ] mdTestEffConddu derivative of the test effective conductivity wrt dof
         */
        const Matrix< DDRMat > &
        dtesteffconddu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Vector< MSI::Dof_Type > &aTestDofTypes,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative wrt der dof type of
         * the test effective conductivity wrt to test dof type
         * @param[ in ] aTestDofTypes a dof type wrt which the test is evaluated
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        virtual void eval_dtesteffconddu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Vector< MSI::Dof_Type > &aTestDofTypes );

        //------------------------------------------------------------------------------
        /**
         * get the the turbulent dynamic viscosity mu_t = rho * vtilde * tf1
         * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
         *               if there are several
         * @param[ out ] mTurbDynVisc effective conductivity
         */
        const Matrix< DDRMat > &turbulent_dynamic_viscosity(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the turbulent dynamic viscosity mu_t = rho * vtilde * tf1
         */
        virtual void eval_turbulent_dynamic_viscosity();

        //------------------------------------------------------------------------------
        /**
         * evaluate the turbulent dynamic viscosity derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        virtual void eval_dturbdynviscdu( const Vector< MSI::Dof_Type > &aDofTypes );

        /**
         * get the derivative of the turbulent dynamic viscosity wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which turbulent dynamic viscosity
         * @param[ out ] mdturbdynviscdu derivative of the turbulent dynamic viscosity wrt dof types
         */
        const Matrix< DDRMat > &dturbdynviscdu(
                const Vector< MSI::Dof_Type > &aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the turbulent dynamic viscosity wrt space
         * @param[ in ] aOrder order of the derivative
         */
        virtual void eval_dturbdynviscdx( uint aOrder );

        /**
         * get the derivative of the turbulent dynamic viscosity wrt space
         * @param[ in ] aOrder order of the derivative
         * @param[ in ] aCMFunctionType enum for specific type of which turbulent dynamic viscosity
         * @param[ out ] mdeffconddx derivative of the turbulent dynamic viscosity wrt space
         */
        const Matrix< DDRMat > &dturbdynviscdx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the turbulent dynamic viscosity derivative wrt to space and dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         */
        virtual void eval_dturbdynviscdxdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                uint                           aOrder );

        /**
         * get the derivative of the turbulent dynamic viscosity wrt space and dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         * @param[ in ] aCMFunctionType enum for specific type of which turbulent dynamic viscosity
         * @param[ out ] mdeffconddxdu derivative of the turbulent dynamic viscosity wrt dof types
         */
        const Matrix< DDRMat > &dturbdynviscdxdu(
                const Vector< MSI::Dof_Type > &aDofType,
                uint                           aOrder,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get the test of the turbulent dynamic viscosity wrt test dof type
         * @param[ in ] aTestDofTypes a dof type wrt which the test is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which turbulent dynamic viscosity
         * @param[ out ] mdTestTurbDynViscdu derivative of the turbulent dynamic viscosity wrt dof types
         */
        const Matrix< DDRMat > &
        testturbdynvisc(
                const Vector< MSI::Dof_Type > &aTestDofTypes,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the test turbulent dynamic viscosity wrt to test dof type
         * @param[ in ] aTestDofTypes a dof type wrt which the test is evaluated
         */
        virtual void eval_testturbdynvisc(
                const Vector< MSI::Dof_Type > &aTestDofTypes );

        //------------------------------------------------------------------------------
        /**
         * get the derivative of the test effective conductivity wrt dof
         * @param[ in ]  aDofTypes group of dof type for derivative
         * @param[ in ]  aTestDofTypes group of dof type for test
         * @param[ out ] mdTestTurbDynViscdu derivative of the test turbulent dyn. viscosity wrt dof
         */
        const Matrix< DDRMat > &
        dtestturbdynviscdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Vector< MSI::Dof_Type > &aTestDofTypes,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative wrt der dof type of
         * the test turbulent dynamic viscosity wrt to test dof type
         * @param[ in ] aTestDofTypes a dof type wrt which the test is evaluated
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        virtual void eval_dtestturbdynviscdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Vector< MSI::Dof_Type > &aTestDofTypes );

        //------------------------------------------------------------------------------
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_TURBULENCE_HPP_ */

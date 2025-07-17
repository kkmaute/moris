/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Fluid_Incompressible_Turbulence.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_TURBULENCE_HPP_
#define SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_TURBULENCE_HPP_

// MRS/COR/src
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
// LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
// FEM/INT/src
#include "cl_FEM_CM_Fluid_Incompressible.hpp"

namespace moris::fem
{

    //--------------------------------------------------------------------------------------------------------------
    class CM_Fluid_Incompressible_Turbulence : public CM_Fluid_Incompressible
    {

        //--------------------------------------------------------------------------------------------------------------

      private:
        // property type for CM
        enum class CM_Property_Type
        {
            DENSITY,          // fluid density
            VISCOSITY,        // fluid dynamic viscosity
            MAX_ENUM
        };

      protected:
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
        bool                    mEffDynViscEval = true;
        moris::Matrix< DDBMat > mdEffDynViscduEval;
        moris::Matrix< DDBMat > mdEffDynViscdxEval;
        moris::Matrix< DDBMat > mdEffDynViscdxduEval;
        moris::Matrix< DDBMat > mTestEffDynViscEval;
        moris::Matrix< DDBMat > mdTestEffDynViscduEval;

        // storage for effective conductivity related evaluation
        Matrix< DDRMat >                     mEffDynVisc;
        Vector< Matrix< DDRMat > >           mdEffDynViscdu;
        Vector< Matrix< DDRMat > >           mdEffDynViscdx;
        Vector< Vector< Matrix< DDRMat > > > mdEffDynViscdxdu;
        Vector< Matrix< DDRMat > >           mTestEffDynVisc;
        Vector< Vector< Matrix< DDRMat > > > mdTestEffDynViscdu;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /*
         * constructor
         */
        CM_Fluid_Incompressible_Turbulence();

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~CM_Fluid_Incompressible_Turbulence() override{};

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
        /*
         * @return constitutive_type
         */
        Constitutive_Type get_constitutive_type() const override
        {
            return Constitutive_Type::FLUID_INCOMPRESSIBLE_TURBULENCE;
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
         * evaluate the traction
         * @param[ in ] aNormal normal
         */
        void eval_traction( const Matrix< DDRMat > &aNormal ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the test traction
         * @param[ in ] aNormal   normal
         */
        void eval_testTraction(
                const Matrix< DDRMat >        &aNormal,
                const Vector< MSI::Dof_Type > &aTestDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model flux derivative wrt to a dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         */
        void eval_dFluxdDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model traction derivative wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aNormal   normal
         */
        void eval_dTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >        &aNormal ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model test traction derivative wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aNormal   normal
         */
        void eval_dTestTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >        &aNormal,
                const Matrix< DDRMat >        &aJump,
                const Vector< MSI::Dof_Type > &aTestDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * get the the effective dynamic viscosity mu_eff = mu + mu_t
         * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
         *               if there are several
         * @param[ out ] mEffDynVisc effective conductivity
         */
        const Matrix< DDRMat > &effective_dynamic_viscosity(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the effective dynamic viscosity
         * mu_eff = mu + mu_t
         */
        void eval_effective_dynamic_viscosity();

        //------------------------------------------------------------------------------
        /**
         * get the derivative of the effective dynamic viscosity wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which effective dynamic viscosity
         * @param[ out ] mdeffconddu derivative of the effective dynamic viscosity wrt dof types
         *               dimensions ( mDof x 1 )
         */
        const Matrix< DDRMat > &deffdynviscdu(
                const Vector< MSI::Dof_Type > &aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the effective dynamic viscosity derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_deffdynviscdu( const Vector< MSI::Dof_Type > &aDofTypes );

        //------------------------------------------------------------------------------
        /**
         * get the derivative of the effective dynamic viscosity wrt space
         * @param[ in ] aOrder order of the derivative
         * @param[ in ] aCMFunctionType enum for specific type of which effective dynamic viscosity
         * @param[ out ] mdeffconddx derivative of the effective dynamic viscosity wrt space
         *               dimensions ( mSpaceDim x 1 )
         */
        const Matrix< DDRMat > &deffdynviscdx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the effective dynamic viscosity wrt space
         * @param[ in ] aOrder order of the derivative
         */
        void eval_deffdynviscdx( uint aOrder );

        //------------------------------------------------------------------------------
        /**
         * get the derivative of the effective dynamic viscosity wrt space and dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         * @param[ in ] aCMFunctionType enum for specific type of which effective dynamic viscosity
         * @param[ out ] mdeffconddxdu derivative of the effective dynamic viscosity wrt dof types
         *               dimensions ( mSpaceDim x mDof )
         */
        const Matrix< DDRMat > &deffdynviscdxdu(
                const Vector< MSI::Dof_Type > &aDofType,
                uint                           aOrder,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the effective dynamic viscosity derivative wrt to space and dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         */
        void eval_deffdynviscdxdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                uint                           aOrder );

        //------------------------------------------------------------------------------
        /**
         * get the test of the effective dynamic viscosity wrt test dof type
         * @param[ in ] aTestDofTypes a dof type wrt which the test is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which turbulent dynamic viscosity
         * @param[ out ] mdTestEffDynViscdu test of the effective dynamic viscosity wrt dof types
         */
        const Matrix< DDRMat > &
        testeffdynvisc(
                const Vector< MSI::Dof_Type > &aTestDofTypes,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the test effective dynamic viscosity wrt to test dof type
         * @param[ in ] aTestDofTypes a dof type wrt which the test is evaluated
         */
        virtual void eval_testeffdynvisc(
                const Vector< MSI::Dof_Type > &aTestDofTypes );

        //------------------------------------------------------------------------------
        /**
         * get the derivative of the test effective dynamic viscosity wrt dof
         * @param[ in ]  aDofTypes group of dof type for derivative
         * @param[ in ]  aTestDofTypes group of dof type for test
         * @param[ out ] mdTestEffDynViscdu derivative of the test effective dynamic viscosity  wrt dof
         */
        const Matrix< DDRMat > &
        dtesteffdynviscdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Vector< MSI::Dof_Type > &aTestDofTypes,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative wrt der dof type of
         * the test effective  dynamic viscosity wrt to test dof type
         * @param[ in ] aTestDofTypes a dof type wrt which the test is evaluated
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        virtual void eval_dtesteffdynviscdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Vector< MSI::Dof_Type > &aTestDofTypes );

        //--------------------------------------------------------------------------------------------------------------
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
         * evaluate the turbulent dynamic viscosity derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        virtual void eval_dturbdynviscdu( const Vector< MSI::Dof_Type > &aDofTypes );

        //------------------------------------------------------------------------------
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
         * evaluate the derivative of the turbulent dynamic viscosity wrt space
         * @param[ in ] aOrder order of the derivative
         */
        virtual void eval_dturbdynviscdx( uint aOrder );

        //------------------------------------------------------------------------------
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
         * evaluate the turbulent dynamic viscosity derivative wrt to space and dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         */
        virtual void eval_dturbdynviscdxdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                uint                           aOrder );

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
         * get the derivative of the test turbulent dyn. viscosity wrt dof
         * @param[ in ]  aDofTypes group of dof type for deriavtive
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

        //--------------------------------------------------------------------------------------------------------------
        /**
         * select derivative wrt to a dof type
         * @param[ in ] aCMRequestType  a type for required derivative
         * @param[ in ] aTestDofTypes   a test dof type wrt which the test traction is evaluated
         * @param[ in ] aNormal         a normal
         * @param[ in ] aJump           a jump
         * @param[ in ] aCMFunctionType
         * Rem: child implementation
         */
        const Matrix< DDRMat > &select_derivative_FD(
                enum CM_Request_Type           aCMRequestType,
                const Vector< MSI::Dof_Type > &aTestDofTypes,
                const Matrix< DDRMat >        &aNormal,
                const Matrix< DDRMat >        &aJump,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * select derivative wrt to a dof type
         * @param[ in ] aCMRequestType  a type for required derivative
         * @param[ in ] aDerivativeFD   a derivative value to set to storage
         * @param[ in ] aTestDofTypes   a test dof type wrt which the test traction is evaluated
         * @param[ in ] aCMFunctionType
         * Rem: child implementation
         */
        void set_derivative_FD(
                enum CM_Request_Type           aCMRequestType,
                Matrix< DDRMat >              &aDerivativeFD,
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Vector< MSI::Dof_Type > &aTestDofTypes,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT ) override;
    };

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_TURBULENCE_HPP_ */

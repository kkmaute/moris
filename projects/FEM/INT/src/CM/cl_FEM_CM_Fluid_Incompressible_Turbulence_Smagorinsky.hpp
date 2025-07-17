/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Fluid_Incompressible_Turbulence_Smagorinsky.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_TURBULENCE_SMAGORINSKY_HPP_
#define SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_TURBULENCE_SMAGORINSKY_HPP_

// MRS/COR/src
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
// LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
// FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_CM_Fluid_Incompressible_Turbulence.hpp"

namespace moris::fem
{

    //--------------------------------------------------------------------------------------------------------------
    class CM_Fluid_Incompressible_Turbulence_Smagorinsky : public CM_Fluid_Incompressible_Turbulence
    {

        //--------------------------------------------------------------------------------------------------------------

      private:
        // property type for CM
        enum class CM_Property_Type
        {
            DENSITY,          // fluid density
            VISCOSITY,        // fluid dynamic viscosity
            WALL_DISTANCE,    // distance to closest wall
            MAX_ENUM
        };

        // storage for fluid strain rate evaluation
        Matrix< DDRMat >                     mFluidStrainRate;
        Vector< Matrix< DDRMat > >           mdFluidStrainRatedu;
        Vector< Matrix< DDRMat > >           mdFluidStrainRatedx;
        Vector< Vector< Matrix< DDRMat > > > mdFluidStrainRatedxdu;

        // flag for fluid strain rate related evaluation

        bool                    mFluidStrainRateEval = true;
        moris::Matrix< DDBMat > mdFluidStrainRateduEval;
        moris::Matrix< DDBMat > mdFluidStrainRatedxEval;
        moris::Matrix< DDBMat > mdFluidStrainRatedxduEval;

      protected:
        // default properties
        std::shared_ptr< Property > mPropWallDistance = nullptr;

        // default parameters
        real mKappa       = 0.4;        // von Karman constant
        real mFluidStrainRateTol = 1.0e-12;    // tolerance on the distance to closest wall

        //--------------------------------------------------------------------------------------------------------------

      public:
        //--------------------------------------------------------------------------------------------------------------
        /*
         * constructor
         */
        CM_Fluid_Incompressible_Turbulence_Smagorinsky();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~CM_Fluid_Incompressible_Turbulence_Smagorinsky() override{};

        //------------------------------------------------------------------------------
        /*
         * @return constitutive_type
         */
        Constitutive_Type get_constitutive_type() const override
        {
            return Constitutive_Type::FLUID_INCOMPRESSIBLE_TURBULENCE_SMAGORINSKY;
        }

        //------------------------------------------------------------------------------
        /**
         * set local properties
         */
        void set_local_properties() override;

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
         * evaluate the turbulent dynamic viscosity mu_t = rho * vtilde * tf1
         */
        void eval_turbulent_dynamic_viscosity() override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the turbulent dynamic viscosity derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dturbdynviscdu(
                const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the turbulent dynamic viscosity wrt space
         * @param[ in ] aOrder order of the derivative
         */
        void eval_dturbdynviscdx(
                uint aOrder ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the turbulent dynamic viscosity derivative wrt to space and dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         */
        void eval_dturbdynviscdxdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                uint                           aOrder ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the test turbulent dynamic viscosity wrt to test dof type
         * @param[ in ] aTestDofTypes a dof type wrt which the test is evaluated
         */
        void eval_testturbdynvisc(
                const Vector< MSI::Dof_Type > &aTestDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative wrt der dof type of
         * the test turbulent dynamic viscosity wrt to test dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aTestDofTypes a dof type wrt which the test is evaluated
         */
        void eval_dtestturbdynviscdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Vector< MSI::Dof_Type > &aTestDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the fluid strain rate
         * 2 * strain : strain
         * @param[ out ] mFluidStrainRate fluid strain rate
         */
        const Matrix< DDRMat > &fluid_strain_rate(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the fluid strain rate
         */
        void eval_fluid_strain_rate();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the derivative of the fluid strain rate wrt dof
         * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
         * @param[ out ] mdFluidStraindu derivative of the fluid strain wrt dof
         */
        const Matrix< DDRMat > &dfluidstrainratedu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the fluid strain rate wrt dof
         */
        void eval_dfluidstrainratedu(
                const Vector< MSI::Dof_Type > &aDofTypes );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the derivative of the fluid strain rate wrt space
         * @param[ in ] aOrder order of the derivative
         */
        const Matrix< DDRMat > &dfluidstrainratedx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the fluid strain ratewrt space
         * @param[ in ] aOrder order of the derivative
         */
        void eval_dfluidstrainratedx( uint aOrder );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the derivative of the fluid strain rate wrt space and dof type
         * @param[ in ] aOrder order of the derivative
         * @param[ in ] aDofTypes vector of dof type for derivative
         *
         */
        const Matrix< DDRMat > &dfluidstrainratedxdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                uint                           aOrder,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the strain rate wrt space and dof type
         * @param[ in ] aOrder    order of the derivative
         * @param[ in ] aDofTypes vector of dof type for derivative
         */
        void eval_dfluidstrainratedxdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                uint                           aOrder );

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
        //------------------------------------------------------------------------------
    };
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_TURBULENCE_SMAGORINSKY_HPP_ */

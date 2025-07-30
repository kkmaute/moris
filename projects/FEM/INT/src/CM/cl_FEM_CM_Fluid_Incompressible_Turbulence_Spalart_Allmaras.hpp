/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_TURBULENCE_SPALART_ALLMARAS_HPP_
#define SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_TURBULENCE_SPALART_ALLMARAS_HPP_

#include <map>
// MRS/COR/src
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
// LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
// FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_CM_Fluid_Incompressible.hpp"
#include "cl_FEM_CM_Fluid_Incompressible_Turbulence.hpp"
#include "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp"

namespace moris::fem
{

    //--------------------------------------------------------------------------------------------------------------
    // Spalart-Allmaras turbulence model
    // where turbulent dynamic viscosity mu_t = rho * vtilde * tf1
    class CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras : public CM_Fluid_Incompressible_Turbulence
    {

        //--------------------------------------------------------------------------------------------------------------

      private:
        // default dof type
        MSI::Dof_Type mDofViscosity = MSI::Dof_Type::VISCOSITY;

        // default properties
        std::shared_ptr< Property > mPropKinViscosity = nullptr;

        // property type for CM
        enum class CM_Property_Type
        {
            DENSITY,          // fluid density
            VISCOSITY,        // fluid dynamic viscosity
            KIN_VISCOSITY,    // fluid kinematic viscosity
            MAX_ENUM
        };

        // flags for chi related evaluation
        bool                    mChiEval = true;
        moris::Matrix< DDBMat > mdChiduEval;
        moris::Matrix< DDBMat > mdChidxEval;
        moris::Matrix< DDBMat > mdChidxduEval;

        // storage for chi
        moris::real                          mChi;
        Vector< Matrix< DDRMat > >           mdChidu;
        Vector< Matrix< DDRMat > >           mdChidx;
        Vector< Vector< Matrix< DDRMat > > > mdChidxdu;

        // flags for fv1 related evaluation
        bool                    mFv1Eval = true;
        moris::Matrix< DDBMat > mdFv1duEval;
        moris::Matrix< DDBMat > mdFv1dxEval;
        moris::Matrix< DDBMat > mdFv1dxduEval;

        // storage for fv1
        moris::real                          mFv1;
        Vector< Matrix< DDRMat > >           mdFv1du;
        Vector< Matrix< DDRMat > >           mdFv1dx;
        Vector< Vector< Matrix< DDRMat > > > mdFv1dxdu;

        //--------------------------------------------------------------------------------------------------------------

      public:
        //--------------------------------------------------------------------------------------------------------------
        /*
         * constructor
         */
        CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~CM_Fluid_Incompressible_Turbulence_Spalart_Allmaras() override{};

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
            return Constitutive_Type::FLUID_INCOMPRESSIBLE_TURBULENCE_SPALART_ALLMARAS;
        }

        //------------------------------------------------------------------------------
        /**
         * set constitutive model dof types
         * @param[ in ] aDofTypes   a list of group of dof types
         * @param[ in ] aDofStrings a list of strings to describe the dof types
         */
        void set_dof_type_list(
                const Vector< Vector< MSI::Dof_Type > > &aDofTypes,
                const Vector< std::string >             &aDofStrings ) override;
        
        //------------------------------------------------------------------------------
        /**
         * set local properties
         */
        void set_local_properties() override;

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
        void eval_dturbdynviscdu( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the turbulent dynamic viscosity wrt space
         * @param[ in ] aOrder order of the derivative
         */
        void eval_dturbdynviscdx( uint aOrder ) override;

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
        void
        eval_testturbdynvisc(
                const Vector< MSI::Dof_Type > &aTestDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative wrt der dof type of
         * the test turbulent dynamic viscosity wrt to test dof type
         * @param[ in ] aTestDofTypes a dof type wrt which the test is evaluated
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dtestturbdynviscdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Vector< MSI::Dof_Type > &aTestDofTypes ) override;

      private:
        //------------------------------------------------------------------------------
        /**
         * get chi = nu_tilde/nu
         * @param[ in ]  aCMFunctionType enum indicating which function if several
         * @param[ out ] mChi chi
         */
        real chi(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * get the derivative of chi wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which chi
         * @param[ out ] mdchidu derivative of the chi wrt dof types
         */
        const Matrix< DDRMat > &
        dchidu(
                const Vector< MSI::Dof_Type > &aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * get the derivative of chi wrt space
         * @param[ in ] aOrder order of the derivative
         * @param[ in ] aCMFunctionType enum for specific type of which chi
         * @param[ out ] mdchidx derivative of the chi wrt space
         */
        const Matrix< DDRMat > &
        dchidx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * get the derivative of chi wrt space and dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         * @param[ in ] aCMFunctionType enum for specific type of which chi
         * @param[ out ] mdchidxdu derivative of the chi wrt dof types
         */
        const Matrix< DDRMat > &
        dchidxdu(
                const Vector< MSI::Dof_Type > &aDofType,
                uint                           aOrder,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get fv1 = chi^3 / ( chi^3 + cv1^3 )
         * @param[ in ]  aCMFunctionType enum for specific type of fv2
         * @param[ out ] mFv1 fv1
         */
        real fv1(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * get the derivative of fv1 wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of fv1
         * @param[ out ] mdfv1du derivative of fv1 wrt dof types
         */
        const Matrix< DDRMat > &
        dfv1du(
                const Vector< MSI::Dof_Type > &aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * get the derivative of fv1 wrt space
         * @param[ in ] aOrder order of the derivative
         * @param[ in ] aCMFunctionType enum for specific type of which fv1
         * @param[ out ] mdfv1dx derivative of the fv1 wrt space
         */
        const Matrix< DDRMat > &
        dfv1dx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * get the derivative of fv1 wrt space and dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         * @param[ in ] aCMFunctionType enum for specific type of which fv1
         * @param[ out ] mdfv1dxdu derivative of the fv1 wrt dof types
         */
        const Matrix< DDRMat > &
        dfv1dxdu(
                const Vector< MSI::Dof_Type > &aDofType,
                uint                           aOrder,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

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

#endif /* SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_TURBULENCE_SPALART_ALLMARAS_HPP_ */

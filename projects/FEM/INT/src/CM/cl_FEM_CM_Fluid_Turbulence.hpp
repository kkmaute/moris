/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Fluid_Turbulence.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_FLUID_TURBULENCE_HPP_
#define SRC_FEM_CL_FEM_CM_FLUID_TURBULENCE_HPP_

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
#include "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp"

namespace moris::fem
{

    //--------------------------------------------------------------------------------------------------------------
    class CM_Fluid_Turbulence : public CM_Fluid_Incompressible
    {

        //--------------------------------------------------------------------------------------------------------------

      private:
        // default dof type
        MSI::Dof_Type mDofVelocity  = MSI::Dof_Type::VX;
        MSI::Dof_Type mDofPressure  = MSI::Dof_Type::P;
        MSI::Dof_Type mDofViscosity = MSI::Dof_Type::VISCOSITY;

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

        // flags for turbulence dynamic viscosity related evaluation
        bool                    mTurbDynViscEval = true;
        moris::Matrix< DDBMat > mdTurbDynViscduEval;
        moris::Matrix< DDBMat > mdTurbDynViscdxEval;
        moris::Matrix< DDBMat > mdTurbDynViscdxduEval;

        // storage for turbulence dynamic viscosity related evaluation
        Matrix< DDRMat >                     mTurbDynVisc;
        Vector< Matrix< DDRMat > >           mdTurbDynViscdu;
        Vector< Matrix< DDRMat > >           mdTurbDynViscdx;
        Vector< Vector< Matrix< DDRMat > > > mdTurbDynViscdxdu;

        // flags for effective conductivity related evaluation
        bool                    mEffDynViscEval = true;
        moris::Matrix< DDBMat > mdEffDynViscduEval;
        moris::Matrix< DDBMat > mdEffDynViscdxEval;
        moris::Matrix< DDBMat > mdEffDynViscdxduEval;

        // storage for effective conductivity related evaluation
        Matrix< DDRMat >                     mEffDynVisc;
        Vector< Matrix< DDRMat > >           mdEffDynViscdu;
        Vector< Matrix< DDRMat > >           mdEffDynViscdx;
        Vector< Vector< Matrix< DDRMat > > > mdEffDynViscdxdu;

      protected:
        // default properties
        std::shared_ptr< Property > mPropKinViscosity = nullptr;

        //--------------------------------------------------------------------------------------------------------------

      public:
        //--------------------------------------------------------------------------------------------------------------
        /*
         * constructor
         */
        CM_Fluid_Turbulence();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~CM_Fluid_Turbulence() override {};

        //------------------------------------------------------------------------------
        /*
         * reset evaluation flag
         * (child implementation)
         */
        void reset_eval_flags() override;

        //------------------------------------------------------------------------------
        /*
         * create a global dof type list including constitutive and property dependencies
         * (child implementation)
         */
        void build_global_dof_type_list() override;

        //------------------------------------------------------------------------------
        /*
         * @return constitutive_type
         */
        Constitutive_Type get_constitutive_type() const override
        {
            return Constitutive_Type::FLUID_TURBULENCE;
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
         * set constitutive model dv types
         * @param[ in ] aDvTypes   a list of group of dv types
         * @param[ in ] aDvStrings a list of strings to describe the dv types
         */
        void set_dv_type_list(
                const Vector< Vector< gen::PDV_Type > > & aDvTypes,
                const Vector< std::string >             & aDvStrings ) override
        {
            Constitutive_Model::set_dv_type_list( aDvTypes );
        }

        //------------------------------------------------------------------------------
        /**
         * set local properties
         */
        void set_local_properties() override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model flux
         */
        void eval_flux() override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the divergence of the flux
         */
        void eval_divflux() override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the divergence of the flux wrt dof type
         */
        void eval_ddivfluxdu( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the flux wrt space
         * @param[ in ] aOrder order of the derivative
         */
        void eval_dfluxdx( uint aOrder ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the traction
         * @param[ in ] aNormal normal
         */
        void eval_traction( const Matrix< DDRMat > &aNormal ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the test traction
         * @param[ in ] aNormal   normal
         */
        void eval_testTraction(
                const Matrix< DDRMat >        &aNormal,
                const Vector< MSI::Dof_Type > &aTestDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model flux derivative wrt to a dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         */
        void eval_dFluxdDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model traction derivative wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aNormal   normal
         */
        void eval_dTractiondDOF(
                const Vector< MSI::Dof_Type > &aDofTypes,
                const Matrix< DDRMat >        &aNormal ) override;

        //--------------------------------------------------------------------------------------------------------------
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
        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the the turbulent dynamic viscosity mu_t = rho * vtilde * tf1
         * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
         *               if there are several
         * @param[ out ] mTurbDynVisc effective conductivity
         */
        const Matrix< DDRMat > &turbulent_dynamic_viscosity(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the the effective dynamic viscosity mu_eff = mu + mu_t
         * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
         *               if there are several
         * @param[ out ] mEffDynVisc effective conductivity
         */
        const Matrix< DDRMat > &effective_dynamic_viscosity(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        /**
         * get the derivative of the effective dynamic viscosity wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which effective dynamic viscosity
         * @param[ out ] mdeffconddu derivative of the effective dynamic viscosity wrt dof types
         */
        const Matrix< DDRMat > &deffdynviscdu(
                const Vector< MSI::Dof_Type > &aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //--------------------------------------------------------------------------------------------------------------

      private:
        //------------------------------------------------------------------------------
        /**
         * evaluate the effective dynamic viscosity mu_eff = mu + mu_t
         */
        void eval_effective_dynamic_viscosity();

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the effective dynamic viscosity wrt space
         * @param[ in ] aOrder order of the derivative
         */
        void eval_deffdynviscdx( uint aOrder );

        /**
         * get the derivative of the effective dynamic viscosity wrt space
         * @param[ in ] aOrder order of the derivative
         * @param[ in ] aCMFunctionType enum for specific type of which effective dynamic viscosity
         * @param[ out ] mdeffconddx derivative of the effective dynamic viscosity wrt space
         */
        const Matrix< DDRMat > &deffdynviscdx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the effective dynamic viscosity derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_deffdynviscdu( const Vector< MSI::Dof_Type > &aDofTypes );

        //------------------------------------------------------------------------------
        /**
         * evaluate the effective dynamic viscosity derivative wrt to space and dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         */
        void eval_deffdynviscdxdu(
                const Vector< MSI::Dof_Type > &aDofTypes,
                uint                           aOrder );

        /**
         * get the derivative of the effective dynamic viscosity wrt space and dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         * @param[ in ] aCMFunctionType enum for specific type of which effective dynamic viscosity
         * @param[ out ] mdeffconddxdu derivative of the effective dynamic viscosity wrt dof types
         */
        const Matrix< DDRMat > &deffdynviscdxdu(
                const Vector< MSI::Dof_Type > &aDofType,
                uint                           aOrder,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the turbulent dynamic viscosity mu_t = rho * vtilde * tf1
         */
        void eval_turbulent_dynamic_viscosity();

        //                /**
        //                 * get the the turbulent dynamic viscosity mu_t = rho * vtilde * tf1
        //                 * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
        //                 *               if there are several
        //                 * @param[ out ] mTurbDynVisc effective conductivity
        //                 */
        //                const Matrix< DDRMat > & turbulent_dynamic_viscosity(
        //                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate the turbulent dynamic viscosity derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dturbdynviscdu( const Vector< MSI::Dof_Type > &aDofTypes );

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
        void eval_dturbdynviscdx( uint aOrder );

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
        void eval_dturbdynviscdxdu(
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
        const Matrix< DDRMat > &dchidu(
                const Vector< MSI::Dof_Type > &aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * get the derivative of chi wrt space
         * @param[ in ] aOrder order of the derivative
         * @param[ in ] aCMFunctionType enum for specific type of which chi
         * @param[ out ] mdchidx derivative of the chi wrt space
         */
        const Matrix< DDRMat > &dchidx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * get the derivative of chi wrt space and dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         * @param[ in ] aCMFunctionType enum for specific type of which chi
         * @param[ out ] mdchidxdu derivative of the chi wrt dof types
         */
        const Matrix< DDRMat > &dchidxdu(
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
        const Matrix< DDRMat > &dfv1du(
                const Vector< MSI::Dof_Type > &aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * get the derivative of fv1 wrt space
         * @param[ in ] aOrder order of the derivative
         * @param[ in ] aCMFunctionType enum for specific type of which fv1
         * @param[ out ] mdfv1dx derivative of the fv1 wrt space
         */
        const Matrix< DDRMat > &dfv1dx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * get the derivative of fv1 wrt space and dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         * @param[ in ] aCMFunctionType enum for specific type of which fv1
         * @param[ out ] mdfv1dxdu derivative of the fv1 wrt dof types
         */
        const Matrix< DDRMat > &dfv1dxdu(
                const Vector< MSI::Dof_Type > &aDofType,
                uint                           aOrder,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );
    };

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_CM_FLUID_TURBULENCE_HPP_ */

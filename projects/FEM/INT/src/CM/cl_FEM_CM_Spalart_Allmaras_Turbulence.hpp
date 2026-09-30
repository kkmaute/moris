/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Spalart_Allmaras_Turbulence.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_SPALART_ALLMARAS_TURBULENCE_HPP_
#define SRC_FEM_CL_FEM_CM_SPALART_ALLMARAS_TURBULENCE_HPP_

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp"

namespace moris::fem
{

    //--------------------------------------------------------------------------------------------------------------

    class CM_Spalart_Allmaras_Turbulence : public Constitutive_Model
    {
      protected:
        // default properties
        std::shared_ptr< Property > mPropKinViscosity = nullptr;
        std::shared_ptr< Property > mPropWallDistance = nullptr;

      private:
        // internal threshold for zero
        // needs to be consistent with threshold set in fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.cpp
        const real mEpsilon      = 1e-10;
        const real mEpsilonDeriv = 1e-60;

        // Spalart-Allmaras model constants
        const real mCb1   = 0.1355;
        const real mCb2   = 0.6220;
        const real mSigma = 2.0 / 3.0;
        const real mKappa = 0.41;
        const real mCw1   = mCb1 / std::pow( mKappa, 2.0 ) + ( 1.0 + mCb2 ) / mSigma;
        const real mCw2   = 0.3;
        const real mCw3   = 2.0;
        const real mCt3   = 1.2;
        const real mCt4   = 0.5;
        const real mCv1   = 7.1;
        const real mCv2   = 0.7;
        const real mCv3   = 0.9;
        const real mRLim  = 10.0;
        const real mCn1   = 16.0;

        // correction to destruction term for negative viscosity
        real mAlpha = 1.0;

        // flag to turn on/off ft2
        bool mUseFt2 = true;

        // default dof type
        MSI::Dof_Type mDofVelocity  = MSI::Dof_Type::VX;
        MSI::Dof_Type mDofViscosity = MSI::Dof_Type::VISCOSITY;

        // property type for CM
        enum class CM_Property_Type
        {
            KINVISCOSITY,    // fluid viscosity
            WALL_DISTANCE,
            MAX_ENUM
        };

        // flags for production coefficient related evaluation
        bool                    mProductionCoeffEval = true;
        moris::Matrix< DDBMat > mdProductionCoeffduEval;
        bool                    mProductionTermEval = true;
        moris::Matrix< DDBMat > mdProductionTermduEval;

        bool                    mWallDestructionCoeffEval = true;
        moris::Matrix< DDBMat > mdWallDestructionCoeffduEval;
        bool                    mWallDestructionTermEval = true;
        moris::Matrix< DDBMat > mdWallDestructionTermduEval;

        bool                    mDiffusionCoeffEval = true;
        moris::Matrix< DDBMat > mdDiffusionCoeffduEval;
        moris::Matrix< DDBMat > mdDiffusionCoeffdxEval;
        moris::Matrix< DDBMat > mdDiffusionCoeffdxduEval;

        bool                    mModVelocityEval = true;
        moris::Matrix< DDBMat > mdModVelocityduEval;

        bool                    mModVelocityLinearizedEval = true;
        moris::Matrix< DDBMat > mdModVelocityLinearizedduEval;

        bool                    mChiEval = true;
        moris::Matrix< DDBMat > mdChiduEval;
        moris::Matrix< DDBMat > mdChidxEval;
        moris::Matrix< DDBMat > mdChidxduEval;

        bool                    mFt2Eval = true;
        moris::Matrix< DDBMat > mdFt2duEval;

        bool                    mSEval = true;
        moris::Matrix< DDBMat > mdSduEval;
        bool                    mSBarEval = true;
        moris::Matrix< DDBMat > mdSBarduEval;
        bool                    mSModEval = true;
        moris::Matrix< DDBMat > mdSModduEval;
        bool                    mSTildeEval = true;
        moris::Matrix< DDBMat > mdSTildeduEval;
        bool                    mWEval = true;
        moris::Matrix< DDBMat > mdWduEval;
        bool                    mFv1Eval = true;
        moris::Matrix< DDBMat > mdFv1duEval;
        bool                    mFv2Eval = true;
        moris::Matrix< DDBMat > mdFv2duEval;

        bool                    mFwEval = true;
        moris::Matrix< DDBMat > mdFwduEval;
        bool                    mGEval = true;
        moris::Matrix< DDBMat > mdGduEval;
        bool                    mREval = true;
        moris::Matrix< DDBMat > mdRduEval;

        bool                    mFnEval = true;
        moris::Matrix< DDBMat > mdFnduEval;
        moris::Matrix< DDBMat > mdFndxEval;
        moris::Matrix< DDBMat > mdFndxduEval;

        // storage for production coefficient related evaluation
        Matrix< DDRMat >           mProductionCoeff;
        Vector< Matrix< DDRMat > > mdProductionCoeffdu;
        Matrix< DDRMat >           mProductionTerm;
        Vector< Matrix< DDRMat > > mdProductionTermdu;

        Matrix< DDRMat >           mWallDestructionCoeff;
        Vector< Matrix< DDRMat > > mdWallDestructionCoeffdu;
        Matrix< DDRMat >           mWallDestructionTerm;
        Vector< Matrix< DDRMat > > mdWallDestructionTermdu;

        Matrix< DDRMat >                     mDiffusionCoeff;
        Vector< Matrix< DDRMat > >           mdDiffusionCoeffdu;
        Vector< Matrix< DDRMat > >           mdDiffusionCoeffdx;
        Vector< Vector< Matrix< DDRMat > > > mdDiffusionCoeffdxdu;

        Matrix< DDRMat >           mModVelocity;
        Vector< Matrix< DDRMat > > mdModVelocitydu;

        Matrix< DDRMat >           mModVelocityLinearized;
        Vector< Matrix< DDRMat > > mdModVelocityLinearizeddu;

        // storage for chi
        moris::real                          mChi;
        Vector< Matrix< DDRMat > >           mdChidu;
        Vector< Matrix< DDRMat > >           mdChidx;
        Vector< Vector< Matrix< DDRMat > > > mdChidxdu;

        moris::real                mFt2;
        Vector< Matrix< DDRMat > > mdFt2du;

        Matrix< DDRMat >           mW;
        Vector< Matrix< DDRMat > > mdWdu;
        moris::real                mS;
        Vector< Matrix< DDRMat > > mdSdu;
        moris::real                mSBar;
        Vector< Matrix< DDRMat > > mdSBardu;
        moris::real                mSMod;
        Vector< Matrix< DDRMat > > mdSModdu;
        moris::real                mSTilde;
        Vector< Matrix< DDRMat > > mdSTildedu;
        moris::real                mFv1;
        Vector< Matrix< DDRMat > > mdFv1du;
        moris::real                mFv2;
        Vector< Matrix< DDRMat > > mdFv2du;

        moris::real                mFw;
        Vector< Matrix< DDRMat > > mdFwdu;
        moris::real                mG;
        Vector< Matrix< DDRMat > > mdGdu;
        moris::real                mR;
        Vector< Matrix< DDRMat > > mdRdu;

        moris::real                          mFn;
        Vector< Matrix< DDRMat > >           mdFndu;
        Vector< Matrix< DDRMat > >           mdFndx;
        Vector< Vector< Matrix< DDRMat > > > mdFndxdu;

        friend class IQI_Spalart_Allmaras_Coefficient;

        //--------------------------------------------------------------------------------------------------------------

      public:
        //--------------------------------------------------------------------------------------------------------------
        /*
         * constructor
         */
        CM_Spalart_Allmaras_Turbulence();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~CM_Spalart_Allmaras_Turbulence() override {};

        //------------------------------------------------------------------------------
        /*
         * @return constitutive_type
         */
        Constitutive_Type
        get_constitutive_type() const override
        {
            return Constitutive_Type::SPALART_ALLMARAS_TURBULENCE;
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * set space dim
         */
        void
        set_space_dim( uint aSpaceDim ) override
        {
            mSpaceDim = aSpaceDim;
        }

        //------------------------------------------------------------------------------
        /**
         * set parameters
         * @param[ in ] aParameters a list of parameters
         */
        void set_parameters( const Vector< Matrix< DDRMat > >& aParameters ) override;

        //------------------------------------------------------------------------------
        /**
         * set constitutive model dof types
         * @param[ in ] aDofTypes a list of group of dof types
         * @param[ in ] aDofStrings a list of strings to describe the dof types
         */
        void set_dof_type_list(
                const Vector< Vector< MSI::Dof_Type > >& aDofTypes,
                const Vector< std::string >&             aDofStrings ) override;

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
         * @param[ in ] aDvTypes   a list of group of dv types
         * @param[ in ] aDvStrings a list of strings to describe the dv types
         */
        void
        set_dv_type_list(
                const Vector< Vector< gen::PDV_Type > >& aDvTypes,
                const Vector< std::string >&             aDvStrings ) override
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
         * get the the production coefficient
         * @param[ in ]  aCMFunctionType enum for specific production term if several
         * @param[ out ] mProductionTerm production term
         */
        const Matrix< DDRMat >&
        production_term(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get the derivative of the production term wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which production term
         * @param[ out ] mdProductionTermdu derivative of the production term wrt dof types
         */
        const Matrix< DDRMat >& 
        dproductiontermdu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );
        
        //------------------------------------------------------------------------------
        /**
         * get the the production coefficient
         * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
         *               if there are several
         * @param[ out ] mProductionCoeff production coefficient
         */
        const Matrix< DDRMat >& 
        production_coefficient(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );
        
        //------------------------------------------------------------------------------
        /**
         * get the derivative of the production coefficient wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which production coefficient
         * @param[ out ] mdproductioncoeffdu derivative of the production coefficient wrt dof types
         */
        const Matrix< DDRMat >& 
        dproductioncoeffdu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get the the wall destruction term
         * @param[ in ]  aCMFunctionType enum for specific wall destruction term if several
         * @param[ out ] mWallDestructionTerm wall destruction term
         */
        const Matrix< DDRMat >& 
        wall_destruction_term(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get the derivative of the wall destruction term wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which wall destruction term
         * @param[ out ] mdwalldestructiontermdu derivative of the wall destruction term wrt dof types
         */
        const Matrix< DDRMat >&
        dwalldestructiontermdu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get the the wall destruction coefficient
         * @param[ in ]  aCMFunctionType enum indicating which wall destruction coefficient
         *               if there are several
         * @param[ out ] mWallDestructionCoeff wall destruction coefficient
         */
        const Matrix< DDRMat >& 
        wall_destruction_coefficient(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get the derivative of the wall destruction coefficient wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which wall destruction coefficient
         * @param[ out ] mdwalldestructioncoeffdu derivative of the wall destruction coefficient wrt dof types
         */
        const Matrix< DDRMat >& 
        dwalldestructioncoeffdu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );
        
        //------------------------------------------------------------------------------
        /**
         * get the diffusion coefficient
         * @param[ in ]  aCMFunctionType enum for specific diffusion coefficient if several
         * @param[ out ] mDiffusionCoeff diffusion coefficient
         */
        const Matrix< DDRMat >& 
        diffusion_coefficient(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );
        
        //------------------------------------------------------------------------------
        /**
         * get the derivative of the diffusion coefficient wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which diffusion
         * @param[ out ] mdproductioncoeffdu derivative of the diffusion wrt dof types
         */
        const Matrix< DDRMat >& 
        ddiffusioncoeffdu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get the modified velocity u_tilde = u - ( cb2 / sigma ) * dnu_tilde/dx
         * @param[ in ]  aCMFunctionType enum for modified velocity if several
         * @param[ out ] mModVelocity modified velocity
         */
        const Matrix< DDRMat >& 
        modified_velocity(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get the derivative of the modified velocity wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which diffusion
         * @param[ out ] mdmodvelocitydu derivative of the modified velocity wrt dof types
         */
        const Matrix< DDRMat >& 
        dmodvelocitydu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get the linearized version of modified velocity u_tilde = u - ( cb2 / sigma ) * dnu_tilde/dx
         * @param[ in ]  aCMFunctionType enum for modified velocity if several
         * @param[ out ] mModVelocity modified velocity
         */
        const Matrix< DDRMat >& modified_velocity_linearized(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get the derivative of the linearized version of the modified velocity wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which diffusion
         * @param[ out ] mdmodvelocitydu derivative of the modified velocity wrt dof types
         */
        const Matrix< DDRMat >& dmodvelocitylinearizeddu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get the derivative of diffusion coefficient wrt space
         * @param[ in ] aOrder order of the derivative
         * @param[ in ] aCMFunctionType enum for specific type of diffusion coefficient
         * @param[ out ] mdDiffusionCoeffdx derivative of the diffusion coefficient wrt space
         */
        const Matrix< DDRMat >& ddiffusioncoeffdx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * get the derivative of the diffusion coefficient wrt space and dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         * @param[ in ] aCMFunctionType enum for specific type of diffusion coefficient
         * @param[ out ] mdDiffusionCoeffdxdu derivative of the diffusion coefficient wrt dof types
         */
        const Matrix< DDRMat >& ddiffusioncoeffdxdu(
                const Vector< MSI::Dof_Type >& aDofType,
                uint                           aOrder,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

      private:
        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the traction
         * @param[ in ] aNormal normal
         */
        void eval_traction( const Matrix< DDRMat >& aNormal ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the test traction
         * @param[ in ] aNormal   normal
         */
        void eval_testTraction(
                const Matrix< DDRMat >&        aNormal,
                const Vector< MSI::Dof_Type >& aTestDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model traction derivative wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aNormal   normal
         */
        void eval_dTractiondDOF(
                const Vector< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&        aNormal ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model test traction derivative wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aNormal   normal
         */
        void eval_dTestTractiondDOF(
                const Vector< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&        aNormal,
                const Vector< MSI::Dof_Type >& aTestDofTypes ) override;
        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the divergence of the flux
         */
        void eval_divflux() override;

        /**
         * evaluate the derivative of the divergence of the flux wrt to dof type
         */
        void eval_ddivfluxdu(
                const Vector< MSI::Dof_Type >& aDofTypes ) override;

        //------------------------------------------------------------------------------
        /**
         * evaluate the production coefficient
         */
        void eval_production_coefficient();

        //------------------------------------------------------------------------------
        /**
         * evaluate the production coefficient derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dproductioncoeffdu( const Vector< MSI::Dof_Type >& aDofTypes );

        //------------------------------------------------------------------------------
        /**
         * evaluate the production term
         */
        void eval_production_term();

        //------------------------------------------------------------------------------
        /**
         * evaluate the production term derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dproductiontermdu( const Vector< MSI::Dof_Type >& aDofTypes );

        //------------------------------------------------------------------------------
        /**
         * evaluate the wall destruction coefficient
         */
        void eval_wall_destruction_coefficient();

        //------------------------------------------------------------------------------
        /**
         * evaluate the wall destruction coefficient derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dwalldestructioncoeffdu( const Vector< MSI::Dof_Type >& aDofTypes );

        //------------------------------------------------------------------------------
        /**
         * evaluate the wall destruction term
         */
        void eval_wall_destruction_term();

        //------------------------------------------------------------------------------
        /**
         * evaluate the wall destruction term derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dwalldestructiontermdu( const Vector< MSI::Dof_Type >& aDofTypes );

        //------------------------------------------------------------------------------
        /**
         * evaluate the diffusion coefficient
         */
        void eval_diffusion_coefficient();

        //------------------------------------------------------------------------------
        /**
         * evaluate the diffusion coefficient derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_ddiffusioncoeffdu( const Vector< MSI::Dof_Type >& aDofTypes );

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative of diffusion coeffcient wrt space
         * @param[ in ] aOrder order of the derivative
         */
        void eval_ddiffusioncoeffdx( uint aOrder );

        //------------------------------------------------------------------------------
        /**
         * evaluate the effective diffusion coefficient wrt to space and dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         */
        void eval_ddiffusioncoeffdxdu(
                const Vector< MSI::Dof_Type >& aDofTypes,
                uint                           aOrder );

        //------------------------------------------------------------------------------
        /**
         * evaluate the modified velocity
         */
        void eval_modified_velocity();

        //------------------------------------------------------------------------------
        /**
         * evaluate the modified velocity derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dmodvelocitydu( const Vector< MSI::Dof_Type >& aDofTypes );

        //------------------------------------------------------------------------------
        /**
         * evaluate the modified velocity
         */
        void eval_modified_velocity_linearized();

        //------------------------------------------------------------------------------
        /**
         * evaluate the modified velocity derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dmodvelocitylinearizeddu( const Vector< MSI::Dof_Type >& aDofTypes );

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
        const Matrix< DDRMat >& dchidu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * get the derivative of chi wrt space
         * @param[ in ] aOrder order of the derivative
         * @param[ in ] aCMFunctionType enum for specific type of which chi
         * @param[ out ] mdeffconddx derivative of the chi wrt space
         */
        const Matrix< DDRMat >& dchidx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * get the derivative of chi wrt space and dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         * @param[ in ] aCMFunctionType enum for specific type of which chi
         * @param[ out ] mdeffconddxdu derivative of the chi wrt dof types
         */
        const Matrix< DDRMat >& dchidxdu(
                const Vector< MSI::Dof_Type >& aDofType,
                uint                           aOrder,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate ft2
         */
        void eval_ft2();

        /**
         * get ft2 = ct3 * exp(-ct4*chi^2)
         * @param[ in ]  aCMFunctionType enum indicating which function if several
         * @param[ out ] mFt2 ft2
         */
        real ft2(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * evaluate ft2 derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dft2du( const Vector< MSI::Dof_Type >& aDofTypes );

        /**
         * get the derivative of ft2 wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which ft2
         * @param[ out ] mdchidu derivative of ft2 wrt dof types
         */
        const Matrix< DDRMat >& dft2du(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate w = 1/2 *(dui/dxk - duj/dxi)
         */
        void eval_w();

        /**
         * get w =
         * @param[ in ]  aCMFunctionType  enum indicating which W
         *               if there are several
         * @param[ out ] mW diffusion coefficient
         */
        const Matrix< DDRMat >& w(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * evaluate the w derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dwdu( const Vector< MSI::Dof_Type >& aDofTypes );

        /**
         * get the derivative of w wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which w
         * @param[ out ] mdwdu derivative of w wrt dof types
         */
        const Matrix< DDRMat >& dwdu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate S
         */
        void eval_s();

        /**
         * get S = sqrt( 2 W : W )
         * @param[ in ]  aCMFunctionType  enum indicating which S
         *               if there are several
         * @param[ out ] mS S
         */
        real s(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * evaluate the s derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dsdu( const Vector< MSI::Dof_Type >& aDofTypes );

        /**
         * get the derivative of S wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which S
         * @param[ out ] mdSdu derivative of S wrt dof types
         */
        const Matrix< DDRMat >& dsdu(
                const Vector< MSI::Dof_Type >& aDofType,
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
        const Matrix< DDRMat >& dfv1du(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate fv2
         */
        void eval_fv2();

        /**
         * get fv2 = 1 - chi / ( 1 + chi * fv1 )
         * @param[ in ]  aCMFunctionType enum for specific type of fv2
         * @param[ out ] mFv2 fv2
         */
        real fv2(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * evaluate the fv2 derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dfv2du( const Vector< MSI::Dof_Type >& aDofTypes );

        /**
         * get the derivative of fv2 wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of fv2
         * @param[ out ] mdFv2du derivative of fv2 wrt dof types
         */
        const Matrix< DDRMat >& dfv2du(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate Sbar
         */
        void eval_sbar();

        /**
         * get Sbar = nu_tilde * fv2 / ( kappa^2 * d^2)
         * @param[ in ]  aCMFunctionType  enum indicating which S if several
         * @param[ out ] mSbar sbar
         */
        real sbar(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * evaluate the Sbar derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dsbardu( const Vector< MSI::Dof_Type >& aDofTypes );

        /**
         * get the derivative of Sbar wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of which S
         * @param[ out ] mdSbardu derivative of Sbar wrt dof types
         */
        const Matrix< DDRMat >& dsbardu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate SMod
         */
        void eval_smod();

        /**
         * get SMod = S * (cv2^2 * S + cv3 * SBar) / ((cv3 - 2 * cv2) * S - SBar)
         * @param[ in ]  aCMFunctionType enum for specific type of SMod
         * @param[ out ] mSMod SMod
         */
        real smod(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * evaluate the SMod derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dsmoddu( const Vector< MSI::Dof_Type >& aDofTypes );

        /**
         * get the derivative of SMod wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of SMod
         * @param[ out ] mdSModdu derivative of SMod wrt dof types
         */
        const Matrix< DDRMat >& dsmoddu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate STilde
         */
        void eval_stilde();

        /**
         * get STilde = S + sBar if SBar <= -cv2 * S
         *            = S + SMod if SBar <  -cv2 * S
         * @param[ in ]  aCMFunctionType enum for specific type of SMod
         * @param[ out ] mSTilde STilde
         */
        real stilde(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * evaluate the STilde derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dstildedu( const Vector< MSI::Dof_Type >& aDofTypes );

        /**
         * get the derivative of STilde wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of STilde
         * @param[ out ] mdSTildedu derivative of STilde wrt dof types
         */
        const Matrix< DDRMat >& dstildedu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate r
         */
        void eval_r();

        /**
         * get r = min( nutilde / (stilde * kappa^2 * d^2), rlim)
         * @param[ in ]  aCMFunctionType enum for specific type of SMod
         * @param[ out ] mR R
         */
        real r(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * evaluate the r derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_drdu( const Vector< MSI::Dof_Type >& aDofTypes );

        /**
         * get the derivative of r wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of r
         * @param[ out ] mdRdu derivative of r wrt dof types
         */
        const Matrix< DDRMat >& drdu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate g
         */
        void eval_g();

        /**
         * get g = r + mCw2 * ( std::pow( r, 6.0 ) - r );
         * @param[ in ]  aCMFunctionType enum for specific type of SMod
         * @param[ out ] mG g
         */
        real g(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * evaluate the g derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dgdu( const Vector< MSI::Dof_Type >& aDofTypes );

        /**
         * get the derivative of g wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of r
         * @param[ out ] mdGdu derivative of g wrt dof types
         */
        const Matrix< DDRMat >& dgdu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate fw
         */
        void eval_fw();

        /**
         * get fw = g * ( 1 + cw3^6) / (g^6 + cw3^6)^(1/6)
         * @param[ in ]  aCMFunctionType enum for specific type of SMod
         * @param[ out ] mFw fw
         */
        real fw(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * evaluate the fw derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dfwdu( const Vector< MSI::Dof_Type >& aDofTypes );

        /**
         * get the derivative of fw wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of r
         * @param[ out ] mdFwdu derivative of fw wrt dof types
         */
        const Matrix< DDRMat >& dfwdu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //------------------------------------------------------------------------------
        /**
         * evaluate fn = (cn1 + chi^3)/(cn1 - chi^3)
         */
        void eval_fn();

        /**
         * get fn =
         * @param[ in ]  aCMFunctionType enum for specific type of fn
         * @param[ out ] mFn fn
         */
        real fn(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * evaluate the fn derivative wrt to dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         */
        void eval_dfndu( const Vector< MSI::Dof_Type >& aDofTypes );

        /**
         * get the derivative of fn wrt dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aCMFunctionType enum for specific type of fn
         * @param[ out ] mdFwdu derivative of fn wrt dof types
         */
        const Matrix< DDRMat >& dfndu(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * evaluate the derivative of fn wrt space
         * @param[ in ] aOrder order of the derivative
         */
        void eval_dfndx( uint aOrder );

        /**
         * get the derivative of fn wrt space
         * @param[ in ] aOrder order of the derivative
         * @param[ in ] aCMFunctionType enum for specific type of fn
         * @param[ out ] mdfndx derivative of the fn wrt space
         */
        const Matrix< DDRMat >& dfndx(
                uint                  aOrder,
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        /**
         * evaluate the effective fn wrt to space and dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         */
        void eval_dfndxdu(
                const Vector< MSI::Dof_Type >& aDofTypes,
                uint                           aOrder );

        /**
         * get the derivative of fn wrt space and dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ in ] aOrder order of the space derivative
         * @param[ in ] aCMFunctionType enum for specific type of fn
         * @param[ out ] mdfndxdu derivative of the fn wrt dof types
         */
        const Matrix< DDRMat >& dfndxdu(
                const Vector< MSI::Dof_Type >& aDofType,
                uint                           aOrder,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //--------------------------------------------------------------------------------------------------------------
    };

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_CM_SPALART_ALLMARAS_TURBULENCE_HPP_ */

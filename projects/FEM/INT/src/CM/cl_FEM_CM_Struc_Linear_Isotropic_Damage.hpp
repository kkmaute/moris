/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Struc_Linear_Isotropic_Damage.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_STRUC_LINEAR_ISOTROPIC_DAMAGE_HPP_
#define SRC_FEM_CL_FEM_CM_STRUC_LINEAR_ISOTROPIC_DAMAGE_HPP_

#include <map>

#include "typedefs.hpp"    //MRS/COR/src
#include "cl_Cell.hpp"     //MRS/CON/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_CM_Struc_Linear_Isotropic.hpp"

namespace moris
{
    namespace fem
    {
        //--------------------------------------------------------------------------------------------------------------

        class CM_Struc_Linear_Isotropic_Damage : public CM_Struc_Linear_Isotropic
        {

          protected:
            // default dof
            MSI::Dof_Type mDofNonlocalEqStrain = MSI::Dof_Type::NLEQSTRAIN;
            MSI::Dof_Type mDofHistory          = MSI::Dof_Type::HISTORY;

          private:
            // property type for CM
            enum class CM_Property_Type_Iso
            {
                EMOD,
                NU,
                MAX_ENUM,
            };

            // parameters for local equivalent strain, damage law, smoothing
            Matrix< DDRMat > mLEqStrainParam;
            Matrix< DDRMat > mDamageParam;
            Matrix< DDRMat > mSmoothParam;

            // function pointer for local equivalent strain
            void ( CM_Struc_Linear_Isotropic_Damage::*m_eval_equivalent_strain )(
                    const Matrix< DDRMat >& aLEqStrainParam ) = nullptr;
            void ( CM_Struc_Linear_Isotropic_Damage::*m_eval_dEqStraindu )(
                    const Cell< MSI::Dof_Type >& aDofTypes ) = nullptr;

            // function pointer for damage law
            void ( CM_Struc_Linear_Isotropic_Damage::*m_eval_damage )(
                    const Matrix< DDRMat >& aDamageParam ) = nullptr;
            void ( CM_Struc_Linear_Isotropic_Damage::*m_eval_dDamagedu )(
                    const Cell< MSI::Dof_Type >& aDofTypes ) = nullptr;

            // function pointer for smoothed damage
            void ( CM_Struc_Linear_Isotropic_Damage::*m_eval_smooth_damage )(
                    const Matrix< DDRMat >& mSmoothParam ) = nullptr;
            void ( CM_Struc_Linear_Isotropic_Damage::*m_eval_dSmoothDamagedu )(
                    const Cell< MSI::Dof_Type >& aDofTypes ) = nullptr;

            // parameters for damage law
            real mKappa0 = 1.0;
            real mK      = 1.0;
            real mAlpha  = 1.0;
            real mBeta   = 1.0;

            // tolerance for zero
            real mEpsilon = 1.0e-12;

            // parameters for smoothing
            real mSmoothC = 1.0;

            // storage for equivalent strain value
            Matrix< DDRMat >                mEqStrain;
            moris::Cell< Matrix< DDRMat > > mdEqStraindu;

            // storage for nonlocal equivalent strain history value
            Matrix< DDRMat >                mHistory;
            moris::Cell< Matrix< DDRMat > > mdHistorydu;

            // storage for nonlocal equivalent strain history value at beginning of time slab
            Matrix< DDRMat >                mHistoryRef;
            moris::Cell< Matrix< DDRMat > > mdHistoryRefdu;

            // storage for damage value
            Matrix< DDRMat >                mDamage;
            moris::Cell< Matrix< DDRMat > > mdDamagedu;

            // storage for smoothed damage value
            Matrix< DDRMat >                mSmoothDamage;
            moris::Cell< Matrix< DDRMat > > mdSmoothDamagedu;

            // flags for equivalent strain related evaluation
            bool                    mEqStrainEval = true;
            moris::Matrix< DDBMat > mdEqStrainduEval;

            // flags for nonlocal equivalent strain history related evaluation
            bool                    mHistoryEval = true;
            moris::Matrix< DDBMat > mdHistoryduEval;

            // flags for nonlocal equivalent strain history at beginning of time slab related evaluation
            bool                    mHistoryRefEval = true;
            moris::Matrix< DDBMat > mdHistoryRefduEval;

            // flags for damage related evaluation
            bool                    mDamageEval = true;
            moris::Matrix< DDBMat > mdDamageduEval;

            // flags for smooth damage related evaluation
            bool                    mSmoothDamageEval = true;
            moris::Matrix< DDBMat > mdSmoothDamageduEval;

            //--------------------------------------------------------------------------------------------------------------

          public:
            /*
             * trivial constructor
             */
            CM_Struc_Linear_Isotropic_Damage(){};

            //--------------------------------------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~CM_Struc_Linear_Isotropic_Damage(){};

            //------------------------------------------------------------------------------
            /**
             * set parameters - specific child implementation
             * @param[ in ] aParameters a list of parameters
             */
            void set_parameters( moris::Cell< Matrix< DDRMat > > aParameters );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * set dof types
             * @param[ in ] aDofTypes a list of group of dof types
             * @param[ in ] aDofStrings a list of strings to describe the dof types
             */
            void set_dof_type_list(
                    Cell< Cell< MSI::Dof_Type > > aDofTypes,
                    Cell< std::string >           aDofStrings );

            //--------------------------------------------------------------------------------------------------------------

            Constitutive_Type
            get_constitutive_type() const
            {
                return Constitutive_Type::STRUC_LIN_ISO_DAMAGE;
            }

            //------------------------------------------------------------------------------
            /*
             * reset evaluation flag
             * Rem: child implementation
             */
            void reset_eval_flags();

            //------------------------------------------------------------------------------
            /*
             * build global dof type list
             * Rem: child implementation
             */
            void build_global_dof_type_list();

            //--------------------------------------------------------------------------------------------------------------
            /**
             * get the equivalent strain
             * @param[ in ]  aCMFunctionType  enum indicating which equivalent strain is called,
             *               if there are several
             * @param[ out ] mEqStrain equivalent strain
             */
            const Matrix< DDRMat >& equivalent_strain(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * get the derivative of the equivalent strain wrt dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ in ] aCMFunctionType enum for specific type of which equivalent strain
             * @param[ out ] mdEqStraindu derivative of the equivalent strain wrt dof types
             */
            const Matrix< DDRMat >& dEqStraindu(
                    const moris::Cell< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type               aCMFunctionType = CM_Function_Type::DEFAULT );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * get the damage
             * @param[ in ]  aCMFunctionType  enum indicating which damage is called,
             *               if there are several
             * @param[ out ] mDamage damage
             */
            const Matrix< DDRMat >& damage(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * get the smoothed damage
             * @param[ in ]  aCMFunctionType  enum indicating which smooth damage is called,
             *               if there are several
             * @param[ out ] mSmoothDamage smoothed damage
             */
            const Matrix< DDRMat >& smooth_damage(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * get the derivative of the smoothed damage wrt dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ in ] aCMFunctionType enum for specific type of which smoothed damage
             * @param[ out ] mdSmoothDamagedu derivative of the smoothed damage wrt dof types
             */
            const Matrix< DDRMat >& dSmoothDamagedu(
                    const moris::Cell< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type               aCMFunctionType = CM_Function_Type::DEFAULT );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * get the nonlocal equivalent strain history
             * @param[ in ]  aCMFunctionType  enum indicating which equivalent strain is called,
             *               if there are several
             * @param[ out ] mHistory nonlocal equivalent strain history
             */
            const Matrix< DDRMat >& history(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * get the derivative of the nonlocal equivalent strain history wrt dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ in ] aCMFunctionType enum for specific type of which equivalent strain
             * @param[ out ] mdHistorydu derivative of the nonlocal equivalent strain history wrt dof types
             */
            const Matrix< DDRMat >& dHistorydu(
                    const moris::Cell< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type               aCMFunctionType = CM_Function_Type::DEFAULT );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * get the nonlocal equivalent strain history at beginning of time slab
             * @param[ in ]  aCMFunctionType  enum indicating which equivalent strain is called,
             *               if there are several
             * @param[ out ] mHistoryRef nonlocal equivalent strain history
             */
            const Matrix< DDRMat >& history_ref(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * get the derivative of the nonlocal equivalent strain history
             * at beginning of time slab wrt dof type
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] aCMFunctionType enum for specific type of which equivalent strain
             * @param[ out ] mdHistoryRefdu derivative of the nonlocal equivalent strain history
             *                              at beginning of time slab wrt dof types
             */
            const Matrix< DDRMat >& dHistoryRefdu(
                    const moris::Cell< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type               aCMFunctionType = CM_Function_Type::DEFAULT );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * get the derivative of the damage wrt dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ in ] aCMFunctionType enum for specific type of which damage
             * @param[ out ] mdDamagedu derivative of the damage wrt dof types
             */
            const Matrix< DDRMat >& dDamagedu(
                    const moris::Cell< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type               aCMFunctionType = CM_Function_Type::DEFAULT );

          protected:
            //--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux
             */
            void eval_flux();

            //--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux derivative wrt to a dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             */
            void eval_dFluxdDOF(
                    const Cell< MSI::Dof_Type >& aDofTypes );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model traction derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * @param[ in ] aNormal   normal
             */
            void eval_dTractiondDOF(
                    const Cell< MSI::Dof_Type >& aDofTypes,
                    const Matrix< DDRMat >&      aNormal );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test traction
             * @param[ in ] aNormal   normal
             */
            void eval_testTraction(
                    const Matrix< DDRMat >&      aNormal,
                    const Cell< MSI::Dof_Type >& aTestDofTypes );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test traction derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * @param[ in ] aNormal   normal
             */
             void eval_dTestTractiondDOF(
                    const Cell< MSI::Dof_Type >& aDofTypes,
                    const Matrix< DDRMat >&      aNormal,
                    const Matrix< DDRMat >&      aJump,
                    const Cell< MSI::Dof_Type >& aTestDofTypes );

           private:
             //--------------------------------------------------------------------------------------------------------------
             /**
              * evaluate the constitutive model equivalent strain
              */
             void eval_equivalent_strain();

             /**
              * formulation based on energy release rate
              * Lemaitre and Chaboche (1990)
              */
             void eval_equivalent_strain_LemaitreChaboche(
                     const Matrix< DDRMat >& aLEqStrainParam );

             /**
              * formulation based on difference in tensile and compressive strength
              * de Vree et al (1995)
              */
             void
             eval_equivalent_strain_deVree_2d_plane_stress(
                     const Matrix< DDRMat >& aLEqStrainParam );
             void
             eval_equivalent_strain_deVree_2d_plane_strain(
                     const Matrix< DDRMat >& aLEqStrainParam );
             void
             eval_equivalent_strain_deVree_3d(
                     const Matrix< DDRMat >& aLEqStrainParam );

             //--------------------------------------------------------------------------------------------------------------
             /**
              * evaluate the constitutive model equivalent strain derivative wrt to a dof type
              * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
              */
             void eval_dEqStraindu(
                     const Cell< MSI::Dof_Type >& aDofTypes );

             /**
              * formulation based on energy release rate
              * Lemaitre and Chaboche (1990)
              */
             void eval_dEqStraindu_LemaitreChaboche(
                     const Cell< MSI::Dof_Type >& aDofTypes );
             /**
              * formulation based on difference in tensile and compressive strength
              * de Vree et al. (1995)
              */
             void eval_dEqStraindu_deVree_2d_plane_stress(
                     const Cell< MSI::Dof_Type >& aDofTypes );
             void eval_dEqStraindu_deVree_2d_plane_strain(
                     const Cell< MSI::Dof_Type >& aDofTypes );
             void eval_dEqStraindu_deVree_3d(
                     const Cell< MSI::Dof_Type >& aDofTypes );

             //--------------------------------------------------------------------------------------------------------------
             /**
              * evaluate the damage
              */
             void eval_damage();

             /**
              * formulation based on linear law
              * mDamage = ( mAlpha / tHistory ) * ( tHistory - mKappa0 ) / ( mAlpha - mKappa0 )
              */
             void eval_damage_linear(
                     const Matrix< DDRMat >& aDamageParam );

             /**
              * formulation based on exponential law
              * mDamage = 1.0 - ( mKappa0 / tHistory ) * ( 1.0 - mAlpha + mAlpha * std::exp( mBeta * ( mKappa0 - tHistory ) ) )
              */
             void eval_damage_exponential(
                     const Matrix< DDRMat >& aDamageParam );

             /**
              * formulation based on smooth exponential law
              * mDamage = std::exp( -std::exp( mAlpha * ( mKappa0 - this->history()( 0 ) ) ) )
              */
             void eval_damage_smooth_exponential(
                     const Matrix< DDRMat >& aDamageParam );

             //--------------------------------------------------------------------------------------------------------------
             /**
              * evaluate the damage derivative wrt to a dof type in current slab
              * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
              */
             void eval_dDamagedu(
                     const Cell< MSI::Dof_Type >& aDofTypes );

             /**
              * formulation based on linear law
              * mDamage = ( mAlpha / tHistory ) * ( tHistory - mKappa0 ) / ( mAlpha - mKappa0 )
              */
             void eval_dDamagedu_linear(
                     const Cell< MSI::Dof_Type >& aDofTypes );

             /**
              * formulation based on exponential law
              * mDamage = 1.0 - ( mKappa0 / tHistory ) * ( 1.0 - mAlpha + mAlpha * std::exp( mBeta * ( mKappa0 - tHistory ) ) )
              */
             void eval_dDamagedu_exponential(
                     const Cell< MSI::Dof_Type >& aDofTypes );

             /**
              * formulation based on smooth exponential law
              * mDamage = std::exp( -std::exp( mAlpha * ( mKappa0 - this->history()( 0 ) ) ) )
              */
             void eval_dDamagedu_smooth_exponential(
                     const Cell< MSI::Dof_Type >& aDofTypes );

             //--------------------------------------------------------------------------------------------------------------
             /**
              * evaluate the smooth damage
              */
             void eval_smooth_damage();

             /**
              * no smoothing
              */
             void eval_smooth_damage_noSmoothing(
                     const Matrix< DDRMat >& aSmoothParam );

             /**
              * Kreisselmeier-Steinhauser smoothing
              */
             void eval_smooth_damage_ks(
                     const Matrix< DDRMat >& aSmoothParam );

             /**
              * Kreisselmeier-Steinhauser corrected smoothing
              */
             void eval_smooth_damage_corrected_ks(
                     const Matrix< DDRMat >& aSmoothParam );

             //--------------------------------------------------------------------------------------------------------------
             /**
              * evaluate the smooth damage derivative wrt to a dof type in current slab
              * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
              */
             void eval_dSmoothDamagedu(
                     const Cell< MSI::Dof_Type >& aDofTypes );

             /**
              * no smoothing
              */
             void eval_dSmoothDamagedu_noSmoothing(
                     const Cell< MSI::Dof_Type >& aDofTypes );

             /**
              * Kreisselmeier-Steinhauser smoothing
              */
             void eval_dSmoothDamagedu_ks(
                     const Cell< MSI::Dof_Type >& aDofTypes );

             /**
              * Kreisselmeier-Steinhauser corrected smoothing
              */
             void eval_dSmoothDamagedu_corrected_ks(
                     const Cell< MSI::Dof_Type >& aDofTypes );

             //--------------------------------------------------------------------------------------------------------------
             /**
              * evaluate the nonlocal equivalent strain history
              */
             void eval_history();

             //--------------------------------------------------------------------------------------------------------------
             /**
              * evaluate the nonlocal equivalent strain history derivative wrt to a dof type in current slab
              * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
              */
             void eval_dHistorydu(
                     const Cell< MSI::Dof_Type >& aDofTypes );

             //--------------------------------------------------------------------------------------------------------------
             /**
              * evaluate the nonlocal equivalent strain history at beginning of time slab
              */
             void eval_history_ref();

             //--------------------------------------------------------------------------------------------------------------
             /**
              * evaluate the nonlocal equivalent strain history derivative at beginning of time slab
              * wrt to a dof type in current slab
              * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
              */
             void eval_dHistoryRefdu(
                     const Cell< MSI::Dof_Type >& aDofTypes );

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
             const Matrix< DDRMat >& select_derivative_FD(
                     enum CM_Request_Type                aCMRequestType,
                     const moris::Cell< MSI::Dof_Type >& aTestDofTypes,
                     const Matrix< DDRMat >&             aNormal,
                     const Matrix< DDRMat >&             aJump,
                     enum CM_Function_Type               aCMFunctionType = CM_Function_Type::DEFAULT );

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
                     enum CM_Request_Type                aCMRequestType,
                     Matrix< DDRMat >&                   aDerivativeFD,
                     const moris::Cell< MSI::Dof_Type >& aDofTypes,
                     const moris::Cell< MSI::Dof_Type >& aTestDofTypes,
                     enum CM_Function_Type               aCMFunctionType = CM_Function_Type::DEFAULT );

             //--------------------------------------------------------------------------------------------------------------
        };
        //--------------------------------------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_STRUC_LINEAR_ISOTROPIC_DAMAGE_HPP_ */

/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Struc_Nonlinear_Isotropic.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_STRUC_NONLINEAR_ISOTROPIC_HPP_
#define SRC_FEM_CL_FEM_CM_STRUC_NONLINEAR_ISOTROPIC_HPP_

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src

namespace moris::fem
{
    //--------------------------------------------------------------------------------------------------------------

    class CM_Struc_Nonlinear_Isotropic : public Constitutive_Model
    {

      protected:
        const real mEpsilon = 1e-10;

        // default dof
        MSI::Dof_Type mDofDispl    = MSI::Dof_Type::UX;
        MSI::Dof_Type mDofTemp     = MSI::Dof_Type::UNDEFINED;
        MSI::Dof_Type mDofPressure = MSI::Dof_Type::UNDEFINED;

        // property type for CM
        enum class CM_Property_Type
        {
            EIGEN_STRAIN,
            MAX_ENUM
        };

        // default local properties
        std::shared_ptr< Property > mPropEigenStrain = nullptr;

        // storage for deformation related evaluation
        Matrix< DDRMat > mLGStrain;    // Lagrangian/Green strain tensor
        Matrix< DDRMat > mEAStrain;    // Eulerian/Almansi strain tensor
        Matrix< DDRMat > mDGStrain;    // Deformation gradient
        // FIXME: storage for deformation gradient is set two times: Strain - Vector notation DefGrad: Matrix notation
        Matrix< DDRMat > mDefGrad;         // deformation gradient
        Matrix< DDRMat > mRCGStrain;       // right Cauchy-Green deformation tensor
        Matrix< DDRMat > mInvRCGStrain;    // inverse of right Cauchy-Green deformation tensor
        Matrix< DDRMat > mLCGStrain;       // left Cauchy-Green deformation tensor
        Matrix< DDRMat > mInvLCGStrain;    // inverse of left Cauchy-Green deformation tensor

        // storage for derivative of deformation related evaluation
        Vector< Matrix< DDRMat > > mdLGStraindu;
        Vector< Matrix< DDRMat > > mdEAStraindu;
        Vector< Matrix< DDRMat > > mdDGStraindu;
        // FIXME: storage for derivative of deformation gradient is set two times: Strain - Vector notation DefGrad: Matrix notation
        Vector< Matrix< DDRMat > > mdDefGraddu;

        // storage for test deformation related evaluation
        Matrix< DDRMat > mLGTestStrain;
        Matrix< DDRMat > mEATestStrain;
        Matrix< DDRMat > mDGTestStrain;
        // FIXME: storage for deformation gradient is set two times: Strain - Vector notation DefGrad: Matrix notation
        Matrix< DDRMat > mTestDefGrad;

        Vector< Matrix< DDRMat > > mdLGTestStraindu;
        Vector< Matrix< DDRMat > > mdEATestStraindu;
        Vector< Matrix< DDRMat > > mdDGTestStraindu;

        // storage for volume change jacobian evaluation
        real mVolumeChangeJ;

        // storage for projected stress
        Matrix< DDRMat > mFluxProj;

        // storage for stress related evaluation
        Matrix< DDRMat > m1PKStress;
        Matrix< DDRMat > m2PKStress;
        Matrix< DDRMat > mCauchyStress;

        Vector< Matrix< DDRMat > > md1PKStressdu;
        Vector< Matrix< DDRMat > > md2PKStressdu;
        Vector< Matrix< DDRMat > > mdCauchyStressdu;

        Matrix< DDRMat > m1PKTraction;
        Matrix< DDRMat > m2PKTraction;
        Matrix< DDRMat > mCauchyTraction;

        Vector< Matrix< DDRMat > > md1PKTractiondu;
        Vector< Matrix< DDRMat > > md2PKTractiondu;
        Vector< Matrix< DDRMat > > mdCauchyTractiondu;

        Vector< Matrix< DDRMat > > m1PKTestTraction;
        Vector< Matrix< DDRMat > > m2PKTestTraction;
        Vector< Matrix< DDRMat > > mCauchyTestTraction;

        Vector< Vector< Matrix< DDRMat > > > md1PKTestTractiondu;
        Vector< Vector< Matrix< DDRMat > > > md2PKTestTractiondu;
        Vector< Vector< Matrix< DDRMat > > > mdCauchyTestTractiondu;

        Matrix< DDRMat > mGeometricStiffness;

        Model_Type mPlaneType = Model_Type::PLANE_STRAIN;    // Plane stress or plane strain, only used in 2d

        Model_Type mTensorType = Model_Type::FULL;    // Hydrostatic or deviatoric (default: full tensor)

      private:
        // map for Voigt notation
        moris::Matrix< DDRMat > mVoigtNonSymMap;
        moris::Matrix< DDRMat > mVoigtSymMap;

        // flag for deformation related evaluation
        bool mDefGradEval      = true;    // FIXME full version of deformation tensor evaluation flag
        bool mDGStrainEval     = true;
        bool mRCGStrainEval    = true;
        bool mInvRCGStrainEval = true;
        bool mLCGStrainEval    = true;
        bool mInvLCGStrainEval = true;
        bool mLGStrainEval     = true;
        bool mEAStrainEval     = true;

        bool                    mLGTestStrainEval = true;
        bool                    mEATestStrainEval = true;
        bool                    mDGTestStrainEval = true;
        moris::Matrix< DDBMat > mdLGTestStrainduEval;
        moris::Matrix< DDBMat > mdEATestStrainduEval;
        moris::Matrix< DDBMat > mdDGTestStrainduEval;

        bool                    mTestDefGradEval = true;
        moris::Matrix< DDBMat > mdDefGradduEval;
        moris::Matrix< DDBMat > mdLGStrainduEval;
        moris::Matrix< DDBMat > mdEAStrainduEval;
        moris::Matrix< DDBMat > mdDGStrainduEval;

        // flag for volume change jacobian evaluation
        bool mVolumeChangeJEval = true;

        // flag for stress related evaluation
        bool m1PKStressEval    = true;
        bool m2PKStressEval    = true;
        bool mCauchyStressEval = true;

        bool mFluxProjEval = true;

        moris::Matrix< DDBMat > md1PKStressduEval;
        moris::Matrix< DDBMat > md2PKStressduEval;
        moris::Matrix< DDBMat > mdCauchyStressduEval;

        bool m1PKTractionEval    = true;
        bool m2PKTractionEval    = true;
        bool mCauchyTractionEval = true;

        moris::Matrix< DDBMat > md1PKTractionduEval;
        moris::Matrix< DDBMat > md2PKTractionduEval;
        moris::Matrix< DDBMat > mdCauchyTractionduEval;

        moris::Matrix< DDBMat > m1PKTestTractionEval;
        moris::Matrix< DDBMat > m2PKTestTractionEval;
        moris::Matrix< DDBMat > mCauchyTestTractionEval;

        moris::Matrix< DDBMat > md1PKTestTractionduEval;
        moris::Matrix< DDBMat > md2PKTestTractionduEval;
        moris::Matrix< DDBMat > mdCauchyTestTractionduEval;

        moris::Matrix< DDBMat > mGeometricStiffnessEval;

        // function pointers

        void ( CM_Struc_Nonlinear_Isotropic::* m_eval_test_deformation_gradient )(
                const Vector< MSI::Dof_Type >& aTestDofTypes ) = nullptr;

        void ( CM_Struc_Nonlinear_Isotropic::* m_eval_dLGStraindDOF )(
                const Vector< MSI::Dof_Type >& aDofTypes ) = nullptr;
        void ( CM_Struc_Nonlinear_Isotropic::* m_eval_dEAStraindDOF )(
                const Vector< MSI::Dof_Type >& aDofTypes ) = nullptr;

        void ( CM_Struc_Nonlinear_Isotropic::* m_eval_dLGTestStraindDOF )(
                const Vector< MSI::Dof_Type >& aDofTypes ) = nullptr;
        void ( CM_Struc_Nonlinear_Isotropic::* m_eval_dEATestStraindDOF )(
                const Vector< MSI::Dof_Type >& aDofTypes ) = nullptr;

        void ( CM_Struc_Nonlinear_Isotropic::* m_eval_inv_right_cauchy_green_deformation_tensor )() = nullptr;
        void ( CM_Struc_Nonlinear_Isotropic::* m_eval_inv_left_cauchy_green_deformation_tensor )()  = nullptr;

        void ( CM_Struc_Nonlinear_Isotropic::* m_eval_volume_change_jacobian )() = nullptr;

        void ( CM_Struc_Nonlinear_Isotropic::* m_flatten_normal )(
                const Matrix< DDRMat >& aNormal,
                Matrix< DDRMat >&       aFlatNormal ) = nullptr;

        void ( CM_Struc_Nonlinear_Isotropic::* m_flatten_normal_nonsym )(
                const Matrix< DDRMat >& aNormal,
                Matrix< DDRMat >&       aFlatNormal ) = nullptr;

        void ( CM_Struc_Nonlinear_Isotropic::* m_proj_jump )(
                const Matrix< DDRMat >& aVector,
                Matrix< DDRMat >&       aProjMatrix ) = nullptr;

        void ( CM_Struc_Nonlinear_Isotropic::* m_proj_sym )(
                const Matrix< DDRMat >& aVector,
                Matrix< DDRMat >&       aProjMatrix ) = nullptr;

        void ( CM_Struc_Nonlinear_Isotropic::* m_proj_nsym )(
                const Matrix< DDRMat >& aVector,
                Matrix< DDRMat >&       aProjMatrix ) = nullptr;

        // number of normal stresses and strains in the tensors
        uint mNumNormalStress;
        uint mNumNormalStrain;

        //--------------------------------------------------------------------------------------------------------------

      public:
        /*
         * trivial constructor
         */
        CM_Struc_Nonlinear_Isotropic(
                enum Constitutive_Type aConstType = Constitutive_Type::STRUC_NON_LIN_ISO );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~CM_Struc_Nonlinear_Isotropic() override {};

        //------------------------------------------------------------------------------
        /*
         * @return constitutive_type
         */
        Constitutive_Type
        get_constitutive_type() const override
        {
            return Constitutive_Type::STRUC_NON_LIN_ISO;
        }

        //------------------------------------------------------------------------------
        /**
         * reset local evaluation flags
         */
        void reset_eval_flags() override;

        //------------------------------------------------------------------------------
        /**
         * create a global dof type list including local constitutive dependencies
         */
        void build_global_dof_type_list() override;

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
        /**
         * set constitutive model dv types
         * @param[ in ] aDvTypes a list of group of dv types
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
         * set spatial dimensions. Modified from base to set function pointers to the appropriate eval_const()
         * @param[ in ] aSpaceDim a spatial dimension
         */
        void set_space_dim( uint aSpaceDim ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Sets the option to use a modified elasticity model (e.g. plane stress, plane stress, etc.)
         * @param aModelType Nonlinear isotropic elasticity supports combinations of Model_Type::PLANE_STRESS or
         * Model_Type::PLANE_STRAIN, and Model_Type::HYDROSTATIC or
         * Model_Type::DEVIATORIC
         */
        void set_model_type( Model_Type aModelType ) override;

        //------------------------------------------------------------------------------
        /*
         * @return plane_type
         */
        Model_Type
        get_plane_type() const override
        {
            return mPlaneType;
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the test deformation gradient
         * @param[ in ]  aTestDofTypes a dof type wrt which the test is evaluated
         * @param[ out ] mTestDefGrad  test deformation gradient delta F wrt test dof
         */
        const Matrix< DDRMat >&
        test_deformation_gradient(
                const Vector< MSI::Dof_Type >& aTestDofTypes );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the deformation gradient
         * @param[ out ] mDefGrad        deformation gradient F = dx/dX = I + du/dX
         */
        const Matrix< DDRMat >&
        deformation_gradient();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the volume change Jacobian
         * @param[ in ]  aCMFunctionType enum indicating which deformation gradient is called, if there are several
         * @param[ out ] mVolumeChangeJ        volume change Jacobina J = det( F )
         */
        const real&
        volume_change_jacobian(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the strain
         * @param[ in ]  aCMFunctionType enum indicating which strain is called
         * @param[ out ] mStrain         strain
         */
        const Matrix< DDRMat >& strain(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the derivative of the strain wrt dof
         * @param[ in ] aDofTypes        a dof type wrt which the derivative is evaluated
         * @param[ in ]  aCMFunctionType enum indicating which strain is called
         * @param[ out ] mdStraindDof    derivative of the strain wrt dof
         */
        const Matrix< DDRMat >& dStraindDOF(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the test strain
         * @param[ in ]  aCMFunctionType enum indicating which strain is called
         * @param[ out ] mTestStrain     test strain
         */
        const Matrix< DDRMat >& testStrain(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the derivative of the test strain wrt dof
         * @param[ in ]  aCMFunctionType  enum indicating which flux is called, if there are several
         * @param[ out ] mdTeststraindDof Derivative of the teststrain wrt. DOF
         */
        const Matrix< DDRMat >& dTestStraindDOF(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the flux
         * @param[ in ]  aCMFunctionType  enum indicating which flux is called, if there are several
         * @param[ out ] mFlux            evaluated flux
         */
        const Matrix< DDRMat >& flux(
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        const Matrix< DDRMat >& flux(
                int                   aFlatType,
                enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the derivative of the flux wrt dof
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ out ] mFluxDofDer derivative of the flux wrt dof
         */
        const Matrix< DDRMat >& dFluxdDOF(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the constitutive model traction
         * @param[ in ]  aNormal   normal
         * @param[ out ] mTraction constitutive model traction
         */
        const Matrix< DDRMat >& traction(
                const Matrix< DDRMat >& aNormal,
                enum CM_Function_Type   aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the derivative of the constitutive model traction wrt dof
         * @param[ in ]  aNormal   normal
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         * @param[ out ] mFluxDofDer derivative of the flux wrt dof
         */
        const Matrix< DDRMat >& dTractiondDOF(
                const Vector< MSI::Dof_Type >& aDofType,
                const Matrix< DDRMat >&        aNormal,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the constitutive model test traction
         * @param[ in ]  aNormal       normal
         * @param[ out ] mTestTraction constitutive model test traction
         */
        const Matrix< DDRMat >& testTraction(
                const Matrix< DDRMat >&        aNormal,
                const Vector< MSI::Dof_Type >& aTestDofTypes,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * get the derivative of the test traction wrt dof
         * @param[ in ]  aDofType           group of dof type
         * @param[ in ]  aNormal            normal
         * @param[ out ] mdTestTractiondDof derivative of the traction wrt dof
         */
        const Matrix< DDRMat >& dTestTractiondDOF(
                const Vector< MSI::Dof_Type >& aDofType,
                const Matrix< DDRMat >&        aNormal,
                const Matrix< DDRMat >&        aJump,
                const Vector< MSI::Dof_Type >& aTestDofTypes,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate geometric jacobian
         * @param[ in ] aDofTypes a dof type wrt which the jacobian is evaluated
         * @param[out ] mGeometricStiffness geometric stiffness contribution to jacobian
         */
        const Matrix< DDRMat >& GeometricStiffness(
                const Vector< MSI::Dof_Type >& aDofTypes,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

        //--------------------------------------------------------------------------------------------------------------

      private:
        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the deformation gradient
         */
        void
        eval_deformation_gradient();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the test deformation gradient wrt to test dof type
         * @param[ in ]  aTestDofTypes a dof type wrt which the test is evaluated
         */
        void
        eval_test_deformation_gradient(
                const Vector< MSI::Dof_Type >& aTestDofTypes )
        {
            ( this->*m_eval_test_deformation_gradient )( aTestDofTypes );
        }

        void eval_test_deformation_gradient_2d( const Vector< MSI::Dof_Type >& aTestDofTypes );
        void eval_test_deformation_gradient_3d( const Vector< MSI::Dof_Type >& aTestDofTypes );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the volume change Jacobian
         */
        void eval_volume_change_jacobian();

        void eval_volume_change_jacobian_2d();
        void eval_volume_change_jacobian_3d();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the right Cauchy-Green deformation tensor C
         */
        void
        eval_right_cauchy_green_deformation_tensor();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the inverse of the right Cauchy-Green deformation tensor C
         */
        void
        eval_inv_right_cauchy_green_deformation_tensor()
        {
            ( this->*m_eval_inv_right_cauchy_green_deformation_tensor )();
        }

        void eval_inv_right_cauchy_green_deformation_tensor_2d();
        void eval_inv_right_cauchy_green_deformation_tensor_3d();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the left Cauchy-Green deformation tensor b
         */
        void
        eval_left_cauchy_green_deformation_tensor();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the inverse of the left Cauchy-Green deformation tensor b
         */
        void
        eval_inv_left_cauchy_green_deformation_tensor()
        {
            ( this->*m_eval_inv_left_cauchy_green_deformation_tensor )();
        }

        void eval_inv_left_cauchy_green_deformation_tensor_2d();
        void eval_inv_left_cauchy_green_deformation_tensor_3d();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the Lagrangian/Green strain tensor E
         */
        void
        eval_lagrangian_green_strain_tensor();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the Eulerian/Almansi strain tensor e
         */
        void
        eval_eulerian_almansi_strain_tensor();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the deformation gradient strain tensor F
         */
        void
        eval_deformation_gradient_strain_tensor();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the test strain based on the Lagrangian/Green strain tensor E
         */
        void
        eval_testStrain_lagrangian_green_strain();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the test strain based on the Eulerian/Almansi strain tensor e
         */
        void
        eval_testStrain_eulerian_almansi_strain();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the test strain based on the deformation gradient strain tensor F
         */
        void
        eval_testStrain_deformation_gradient_strain();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of Lagrangian/Green strain tensor wrt dof
         * @param[ in ]  aDofTypes     a dof type wrt which the derivative is evaluated
         */
        void
        eval_dLGStraindDOF( const Vector< MSI::Dof_Type >& aDofTypes )
        {
            ( this->*m_eval_dLGStraindDOF )( aDofTypes );
        }

        void eval_dLGStraindDOF_2d( const Vector< MSI::Dof_Type >& aDofTypes );
        void eval_dLGStraindDOF_3d( const Vector< MSI::Dof_Type >& aDofTypes );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of Eulerian/Almansi strain tensor wrt dof
         * @param[ in ]  aDofTypes     a dof type wrt which the derivative is evaluated
         */
        void
        eval_dEAStraindDOF( const Vector< MSI::Dof_Type >& aDofTypes )
        {
            ( this->*m_eval_dLGStraindDOF )( aDofTypes );
        }

        void eval_dEAStraindDOF_2d( const Vector< MSI::Dof_Type >& aDofTypes );
        void eval_dEAStraindDOF_3d( const Vector< MSI::Dof_Type >& aDofTypes );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the deformation gradient strain tensor wrt dof
         * @param[ in ]  aDofTypes     a dof type wrt which the derivative is evaluated
         */
        void
        eval_dDGStraindDOF( const Vector< MSI::Dof_Type >& aDofTypes );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the first Piola-Kirchhoff stress tensor
         */
        virtual void eval_flux_first_piola_kirchhoff()
        {
            MORIS_ERROR( false, "eval_flux_first_piola_kirchhoff - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the second Piola-Kirchhoff stress tensor
         */
        virtual void eval_flux_second_piola_kirchhoff()
        {
            MORIS_ERROR( false, "eval_flux_second_piola_kirchhoff - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the first Piola-Kirchhoff stress tensor
         */
        virtual void eval_geometric_stiffness_second_piola_kirchhoff( const Vector< MSI::Dof_Type >& aDofTypes )
        {
            MORIS_ERROR( false, "eval_geometric_stiffness_second_piola_kirchhoff - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the second Piola-Kirchhoff stress tensor
         */
        virtual void eval_flux_cauchy()
        {
            MORIS_ERROR( false, "eval_flux_cauchy - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of first Piola-Kirchhoff stress tensor wrt dof
         * @param[ in ]  aDofTypes     a dof type wrt which the derivative is evaluated
         */
        virtual void eval_d1PKStressdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
        {
            MORIS_ERROR( false, "eval_d1PKStressdDOF - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of second Piola-Kirchhoff stress tensor wrt dof
         * @param[ in ]  aDofTypes     a dof type wrt which the derivative is evaluated
         */
        virtual void eval_d2PKStressdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
        {
            MORIS_ERROR( false, "eval_d2PKStressdDOF - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of cauchy stress tensor wrt dof
         * @param[ in ]  aDofTypes     a dof type wrt which the derivative is evaluated
         */
        virtual void eval_dCauchyStressdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
        {
            MORIS_ERROR( false, "eval_dCauchyStressdDOF - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model matrix
         */
        void eval_const() override
        {
            MORIS_ERROR( false, "eval_const - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the traction based on first Piola-Kirchhoff stress Tp = P N
         * @param[ in ] aNormal normal in the reference configuration
         */
        virtual void eval_traction_first_piola_kirchhoff( const Matrix< DDRMat >& aNormal )
        {
            MORIS_ERROR( false, "eval_traction_first_piola_kirchhoff - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the traction based on second Piola-Kirchhoff stress Ts = S N
         * @param[ in ] aNormal normal in the reference configuration
         */
        virtual void eval_traction_second_piola_kirchhoff( const Matrix< DDRMat >& aNormal )
        {
            MORIS_ERROR( false, "eval_traction_second_piola_kirchhoff - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the traction based on cauchy stress t = sigma n
         * @param[ in ] aNormal normal in the reference configuration
         */
        virtual void eval_traction_cauchy( const Matrix< DDRMat >& aNormal )
        {
            MORIS_ERROR( false, "eval_traction_cauchy - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of traction based on first Piola-Kirchhoff stress tensor wrt dof
         * @param[ in ]  aDofTypes     a dof type wrt which the derivative is evaluated
         */
        virtual void eval_dTractiondDOF_first_piola_kirchhoff(
                const Matrix< DDRMat >&        aNormal,
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            MORIS_ERROR( false, "eval_d1PKTractiondDOF - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of traction based on second Piola-Kirchhoff stress tensor wrt dof
         * @param[ in ]  aDofTypes     a dof type wrt which the derivative is evaluated
         */
        virtual void eval_dTractiondDOF_second_piola_kirchhoff(
                const Matrix< DDRMat >&        aNormal,
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            MORIS_ERROR( false, "eval_d2PKTractiondDOF - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of traction based on cauchy stress tensor wrt dof
         * @param[ in ]  aDofTypes     a dof type wrt which the derivative is evaluated
         */
        virtual void eval_dTractiondDOF_cauchy(
                const Matrix< DDRMat >&        aNormal,
                const Vector< MSI::Dof_Type >& aDofTypes )
        {
            MORIS_ERROR( false, "eval_dCauchyTractiondDOF - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the test traction based on the first Piola-Kirchhoff stress
         * @param[ in ] aNormal normal in the reference configuration
         */
        virtual void eval_testTraction_first_piola_kirchhoff(
                const Matrix< DDRMat >&        aNormal,
                const Vector< MSI::Dof_Type >& aTestDofTypes )
        {
            MORIS_ERROR( false, "eval_testTraction_first_piola_kirchhoff - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the test traction based on the second Piola-Kirchhoff stress
         * @param[ in ] aNormal normal in the reference configuration
         */
        virtual void eval_testTraction_second_piola_kirchhoff(
                const Matrix< DDRMat >&        aNormal,
                const Vector< MSI::Dof_Type >& aTestDofTypes )
        {
            MORIS_ERROR( false, "eval_testTraction_second_piola_kirchhoff - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the test traction based on the Cauchy stress
         * @param[ in ] aNormal normal in the reference configuration
         */
        virtual void eval_testTraction_cauchy(
                const Matrix< DDRMat >&        aNormal,
                const Vector< MSI::Dof_Type >& aTestDofTypes )
        {
            MORIS_ERROR( false, "eval_testTraction_cauchy - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the test traction derivative based on the first Piola-Kirchhoff stress
         * wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aNormal   normal in the reference configuration
         */
        virtual void eval_dTestTractiondDOF_first_piola_kirchhoff(
                const Vector< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&        aNormal,
                const Matrix< DDRMat >&        aJump,
                const Vector< MSI::Dof_Type >& aTestDofTypes )
        {
            MORIS_ERROR( false, "eval_dTestTractiondDOF_first_piola_kirchhoff - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the test traction derivative based on the second Piola-Kirchhoff stress
         * wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aNormal   normal in the reference configuration
         */
        virtual void eval_dTestTractiondDOF_second_piola_kirchhoff(
                const Vector< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&        aNormal,
                const Matrix< DDRMat >&        aJump,
                const Vector< MSI::Dof_Type >& aTestDofTypes )
        {
            MORIS_ERROR( false, "eval_dTestTractiondDOF_second_piola_kirchhoff - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the test traction derivative based on the cauchy stress
         * wrt to a dof type
         * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
         * @param[ in ] aNormal   normal in the reference configuration
         */
        virtual void eval_dTestTractiondDOF_cauchy(
                const Vector< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&        aNormal,
                const Matrix< DDRMat >&        aJump,
                const Vector< MSI::Dof_Type >& aTestDofTypes )
        {
            MORIS_ERROR( false, "eval_dTestTractiondDOF_cauchy - Not implemented in parent. " );
        }

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the Lagrangian/Green test strain wrt to a dof type
         */
        void eval_dLGTestStraindDOF(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT )
        {
            ( this->*m_eval_dLGTestStraindDOF )( aDofType );
        }

        void eval_dLGTestStraindDOF_2d( const Vector< MSI::Dof_Type >& aDofType );
        void eval_dLGTestStraindDOF_3d( const Vector< MSI::Dof_Type >& aDofType );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the Eulerian/Almansi test strain wrt to a dof type
         */
        void eval_dEATestStraindDOF(
                const Vector< MSI::Dof_Type >& aDofType,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT )
        {
            ( this->*m_eval_dEATestStraindDOF )( aDofType );
        }

        void eval_dEATestStraindDOF_2d( const Vector< MSI::Dof_Type >& aDofType );
        void eval_dEATestStraindDOF_3d( const Vector< MSI::Dof_Type >& aDofType );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the deformation gradient test strain wrt to a dof type
         */
        void eval_dDGTestStraindDOF(
                const Vector< MSI::Dof_Type >& aDofType );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model matrix derivative wrt to a dof type
         * @param[ in ] aDofTypes   a dof type wrt which the derivative is evaluated
         * @param[ in ] adConstdDOF a matrix to fill with derivative evaluation
         */
        void eval_dConstdDOF( const Vector< MSI::Dof_Type >& aDofTypes ) override
        {
            MORIS_ERROR( false, "eval_dConstdDOF - Not implemented in parent. " );
        }

      protected:
        //--------------------------------------------------------------------------------------------------------------
        /**
         *  Sets the appropriate function pointers
         *  based on the current member data
         *  for spatial dimensions and model types
         */
        void set_function_pointers();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * flatten normal vector
         * @param[ in ] aNormal          a normal vector
         * @param[ in ] aFlattenedNormal a matrix for flattened normal to fill
         */
        void flatten_normal(
                const Matrix< DDRMat >& aNormal,
                Matrix< DDRMat >&       aFlatNormal )
        {
            ( this->*m_flatten_normal )( aNormal, aFlatNormal );
        }

        void flatten_normal_2d(
                const Matrix< DDRMat >& aNormal,
                Matrix< DDRMat >&       aFlatNormal );

        void flatten_normal_3d(
                const Matrix< DDRMat >& aNormal,
                Matrix< DDRMat >&       aFlatNormal );

        void flatten_normal_nonsym(
                const Matrix< DDRMat >& aNormal,
                Matrix< DDRMat >&       aFlatNormal )
        {
            ( this->*m_flatten_normal_nonsym )( aNormal, aFlatNormal );
        }

        void flatten_normal_nonsym_2d(
                const Matrix< DDRMat >& aNormal,
                Matrix< DDRMat >&       aFlatNormal );

        void flatten_normal_nonsym_3d(
                const Matrix< DDRMat >& aNormal,
                Matrix< DDRMat >&       aFlatNormal );

        //--------------------------------------------------------------------------------------------------------------
        // FIXME implementation fixed for 2d/3d
        /**
         * transform Voigt vector into full matrix for non-symmetric case
         * @param[ in ] aVoigtVector a matrix in Voigt notation provided
         * @param[ in ] aFullMatrix  a matrix in full notation to fill
         */
        void voigt_to_full_nonsym(
                const Matrix< DDRMat >& aVoigtVector,
                Matrix< DDRMat >&       aFullMatrix );

        /**
         * transform Voigt vector into full matrix for symmetric case for strain
         * @param[ in ] aVoigtVector a matrix in Voigt notation provided
         * @param[ in ] aFullMatrix  a matrix in full notation to fill
         */
        void voigt_to_full_sym_strain(
                const Matrix< DDRMat >& aVoigtVector,
                Matrix< DDRMat >&       aFullMatrix );

        /**
         * transform Voigt vector into full matrix for symmetric case for stress
         * @param[ in ] aVoigtVector a matrix in Voigt notation provided
         * @param[ in ] aFullMatrix  a matrix in full notation to fill
         */
        void voigt_to_full_sym_stress(
                const Matrix< DDRMat >& aVoigtVector,
                Matrix< DDRMat >&       aFullMatrix );

        /**
         * transform full matrix into Voigt vector for non-symmetric case
         * @param[ in ] aFullMatrix  a matrix in full notation provided
         * @param[ in ] aVoigtVector a matrix in Voigt notation to fill
         */
        void full_to_voigt_nonsym(
                const Matrix< DDRMat >& aFullMatrix,
                Matrix< DDRMat >&       aVoigtVector );

        /**
         * transform full matrix into Voigt vector for symmetric case for strain
         * @param[ in ] aFullMatrix  a matrix in full notation provided
         * @param[ in ] aVoigtVector a matrix in Voigt notation to fill
         */
        void full_to_voigt_sym_strain(
                const Matrix< DDRMat >& aFullMatrix,
                Matrix< DDRMat >&       aVoigtVector );

        /**
         * transform full matrix into Voigt vector for symmetric case for stress
         * @param[ in ] aFullMatrix  a matrix in full notation provided
         * @param[ in ] aVoigtVector a matrix in Voigt notation to fill
         */
        void full_to_voigt_sym_stress(
                const Matrix< DDRMat >& aFullMatrix,
                Matrix< DDRMat >&       aVoigtVector );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * derive the projected Matrix
         * @param[ in ] aVector          a vector
         * @param[ in ] aProjMatrix      a projected matrix
         */
        void proj_sym(
                const Matrix< DDRMat >& aVector,
                Matrix< DDRMat >&       aProjMatrix )
        {
            ( this->*m_proj_sym )( aVector, aProjMatrix );
        }

        void proj_sym_2d(
                const Matrix< DDRMat >& aVector,
                Matrix< DDRMat >&       aProjMatrix );

        void proj_sym_3d(
                const Matrix< DDRMat >& aVector,
                Matrix< DDRMat >&       aProjMatrix );

        void eval_flux_proj_sym(
                enum CM_Function_Type aCMFunctionType );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * derive the projected Matrix
         * @param[ in ] aVector          a vector
         * @param[ in ] aProjMatrix      a projected matrix
         */
        void proj_nsym(
                const Matrix< DDRMat >& aVector,
                Matrix< DDRMat >&       aProjMatrix )
        {
            ( this->*m_proj_nsym )( aVector, aProjMatrix );
        }

        void proj_nsym_2d(
                const Matrix< DDRMat >& aVector,
                Matrix< DDRMat >&       aProjMatrix );

        void proj_nsym_3d(
                const Matrix< DDRMat >& aVector,
                Matrix< DDRMat >&       aProjMatrix );

        void eval_flux_proj_nsym(
                enum CM_Function_Type aCMFunctionType );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Flattening of the jump
         * @param aVector     jump in form of vector
         * @param aProjMatrix jump in form of matrix
         */
        void proj_jump(
                const Matrix< DDRMat >& aVector,
                Matrix< DDRMat >&       aProjMatrix )
        {
            ( this->*m_proj_jump )( aVector, aProjMatrix );
        }

        void proj_jump_2d(
                const Matrix< DDRMat >& aVector,
                Matrix< DDRMat >&       aProjMatrix );

        void proj_jump_3d(
                const Matrix< DDRMat >& aVector,
                Matrix< DDRMat >&       aProjMatrix );

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
                enum CM_Request_Type           aCMRequestType,
                const Vector< MSI::Dof_Type >& aTestDofTypes,
                const Matrix< DDRMat >&        aNormal,
                const Matrix< DDRMat >&        aJump,
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
                Matrix< DDRMat >&              aDerivativeFD,
                const Vector< MSI::Dof_Type >& aDofTypes,
                const Vector< MSI::Dof_Type >& aTestDofTypes,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT ) override;
    };
    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_CM_STRUC_NONLINEAR_ISOTROPIC_HPP_ */

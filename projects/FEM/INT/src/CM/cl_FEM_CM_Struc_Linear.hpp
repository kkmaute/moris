/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Struc_Linear.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_STRUC_LINEAR_HPP_
#define SRC_FEM_CL_FEM_CM_STRUC_LINEAR_HPP_

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

    class CM_Struc_Linear : public Constitutive_Model
    {
      protected:
        std::shared_ptr< Property > mPropCTE         = nullptr;
        std::shared_ptr< Property > mPropTemp        = nullptr;
        std::shared_ptr< Property > mPropTRef        = nullptr;
        std::shared_ptr< Property > mPropRotAxis     = nullptr;
        std::shared_ptr< Property > mPropEigenStrain = nullptr;

        // default dof
        MSI::Dof_Type mDofDispl    = MSI::Dof_Type::UX;
        MSI::Dof_Type mDofTemp     = MSI::Dof_Type::UNDEFINED;
        MSI::Dof_Type mDofPressure = MSI::Dof_Type::UNDEFINED;

        // property type for CM
        enum class CM_Property_Type_Lin
        {
            CTE,
            TEMP_PROP,
            TEMP_REF,
            ROT_AXI,
            EIGEN_STRAIN,
            MAX_ENUM
        };

        // function pointers
        void ( CM_Struc_Linear::* m_eval_strain )() = nullptr;

        void ( CM_Struc_Linear::* m_eval_teststrain )() = nullptr;

        void ( CM_Struc_Linear::* m_flatten_normal )(
                const Matrix< DDRMat >& aNormal,
                Matrix< DDRMat >&       aFlatNormal ) = nullptr;

        void ( CM_Struc_Linear::* mConstFunc )(
                std::initializer_list< const real >&& tParams ) = &CM_Struc_Linear::full_3d;

        void ( CM_Struc_Linear::* m_eval_inv_bulk_modulus )(
                const real& aNu,
                const real& aEMod,
                real&       aInvBulkModulus ) = &CM_Struc_Linear::eval_inv_bulk_modulus_generic;

        Model_Type mPlaneType = Model_Type::PLANE_STRESS;    // Plane stress or plane strain, only used in 2d

        Model_Type mTensorType = Model_Type::FULL;    // Hydrostatic or deviatoric (default: full tensor)

        // number of normal stresses and strains in the tensors
        uint mNumNormalStress = MORIS_UINT_MAX;
        uint mNumNormalStrain = MORIS_UINT_MAX;

        // geometric stiffness
        Matrix< DDRMat > mGeometricStiffness;

        moris::Matrix< DDBMat > mGeometricStiffnessEval;

        //--------------------------------------------------------------------------------------------------------------

      public:
        /*
         * trivial constructor
         */
        CM_Struc_Linear();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~CM_Struc_Linear() override {};

        //------------------------------------------------------------------------------
        /*
         * @return constitutive_type
         */
        Constitutive_Type
        get_constitutive_type() const override
        {
            return Constitutive_Type::STRUC_LIN_ISO;
        }

        //------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate geometric jacobian
         * @param[ in ] aDofTypes a dof type wrt which the jacobian is evaluated
         * @param[out ] mGeometricStiffness geometric stiffness contribution to jacobian
         */
        const Matrix< DDRMat >& GeometricStiffness(
                const Vector< MSI::Dof_Type >& aDofTypes,
                enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

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

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the inverse of the bulk modulus, K
         * @return 1/K
         */
        virtual real eval_inv_bulk_modulus();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the inverse of the bulk modulus, K
         * wrt dof type
         * @param[ in ] aDofTypes            a dof type wrt which the derivative is evaluated
         * @param[ out ] dInvBulkModulusdDOF derivative of K
         */
        virtual Matrix< DDRMat > eval_dInvBulkModulusdDOF( const Vector< MSI::Dof_Type >& aDofTypes );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * set spatial dimensions. Modified from base to set function pointers to the appropriate eval_const()
         * @param[ in ] aSpaceDim a spatial dimension
         */
        void set_space_dim( uint aSpaceDim ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Sets the option to use a modified elasticity model (e.g. plane stress, plane stress, etc.)
         * @param aModelType Linear isotropic elasticity supports combinations of Model_Type::PLANE_STRESS or
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
         * returns the E prime values used in the computation of the Stress Intensity Factor(s)
         */
        //     virtual
        //     real get_e_prime();

        //--------------------------------------------------------------------------------------------------------------

      protected:
        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the inverse of the bulk modulus, K
         * @param[ in ]  aNu             Poisson ratio
         * @param[ in ]  aEMod           Elasticity modulus
         * @param[ out ] aInvBulkModulus 1/K
         */
        virtual void eval_inv_bulk_modulus_generic(
                const real& aNu,
                const real& aEMod,
                real&       aInvBulkModulus );

        virtual void eval_inv_bulk_modulus_plane_stress(
                const real& aNu,
                const real& aEMod,
                real&       aInvBulkModulus );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model flux
         */
        void eval_flux() override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model test flux
         */
        void eval_testFlux();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model traction
         * @param[ in ] aNormal normal
         */
        void eval_traction( const Matrix< DDRMat >& aNormal ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model test traction
         * @param[ in ] aNormal   normal
         */
        void eval_testTraction(
                const Matrix< DDRMat >&        aNormal,
                const Vector< MSI::Dof_Type >& aTestDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the strain
         */
        void
        eval_strain() override
        {
            ( this->*m_eval_strain )();
        }

        void eval_strain_2d();
        void eval_strain_3d();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the test strain
         */
        void
        eval_testStrain() override
        {
            ( this->*m_eval_teststrain )();
        }
        void eval_teststrain_2d();
        void eval_teststrain_3d();

        void eval_geometric_stiffness( const Vector< MSI::Dof_Type >& aDofTypes );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model flux derivative wrt to a dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         */
        void eval_dFluxdDOF( const Vector< MSI::Dof_Type >& aDofTypes ) override;

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
                const Matrix< DDRMat >&        aJump,
                const Vector< MSI::Dof_Type >& aTestDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model strain derivative wrt to a dof type
         * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
         * @param[ in ] adStraindDOF a matrix to fill with derivative evaluation
         */
        void eval_dStraindDOF( const Vector< MSI::Dof_Type >& aDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model matrix derivative wrt to a dof type
         * @param[ in ] aDofTypes   a dof type wrt which the derivative is evaluated
         * @param[ in ] adConstdDOF a matrix to fill with derivative evaluation
         */
        void eval_dConstdDOF( const Vector< MSI::Dof_Type >& aDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         *  Sets the appropriate function pointers
         *  based on the current member data
         *  for spatial dimensions and model types
         */
        virtual void set_function_pointers();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * flatten normal vector
         * @param[ in ] aNormal          a normal vector
         * @param[ in ] aFlattenedNormal a matrix for flattened normal to fill
         */
        void
        flatten_normal(
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

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model matrix
         */
        void eval_const() override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Full plane stress tensor
         * @param[ in ] aEmod Elastic modulus
         * @param[ in ] aNu   Poisson ratio
         */
        virtual void full_plane_stress(
                std::initializer_list< const real >&& tParams );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Deviatoric plane stress tensor
         * @param aEmod Elastic modulus
         * @param aNu   Poisson ratio
         */
        virtual void deviatoric_plane_stress(
                std::initializer_list< const real >&& tParams );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Full plane strain tensor
         * @param aEmod Elastic modulus
         * @param aNu   Poisson ratio
         */
        virtual void full_plane_strain(
                std::initializer_list< const real >&& tParams );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Deviatoric plane strain tensor
         * @param aEmod Elastic modulus
         * @param aNu Poisson ratio
         */
        virtual void deviatoric_plane_strain(
                std::initializer_list< const real >&& tParams );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Full axisymmetric tensor
         * @param[ in ] aEmod Elastic modulus
         * @param[ in ] aNu   Poisson ratio
         */
        virtual void full_axisymmetric(
                std::initializer_list< const real >&& tParams );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Deviatoric axisymmetric tensor
         * @param aEmod Elastic modulus
         * @param aNu   Poisson ratio
         */
        virtual void deviatoric_axisymmetric( std::initializer_list< const real >&& tParams );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Full 3d linear isotropic tensor
         * @param aEmod Elastic modulus
         * @param aNu Poisson ratio
         */
        virtual void full_3d( std::initializer_list< const real >&& tParams );

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Deviatoric 3d tensor
         * @param aEmod Elastic modulus
         * @param aNu Poisson ratio
         */
        virtual void deviatoric_3d( std::initializer_list< const real >&& tParams );

        //--------------------------------------------------------------------------------------------------------------
    };
    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_CM_STRUC_LINEAR_ISOTROPIC_HPP_ */

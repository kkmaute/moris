/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Struc_Linear_MoriTanaka.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_STRUC_LINEAR_MoriTanaka_HPP_
#define SRC_FEM_CL_FEM_CM_STRUC_LINEAR_MoriTanaka_HPP_

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_CM_Struc_Linear.hpp"

namespace moris::fem
{
    //--------------------------------------------------------------------------------------------------------------

    class CM_Struc_Linear_MoriTanaka : public CM_Struc_Linear
    {
      protected:
        // default local properties
        std::shared_ptr< Property > mPropEModMat        = nullptr;
        std::shared_ptr< Property > mPropEModFib        = nullptr;
        std::shared_ptr< Property > mPropPoissonMat     = nullptr;
        std::shared_ptr< Property > mPropPoissonFib     = nullptr;
        std::shared_ptr< Property > mPropVolumeFraction = nullptr;
        std::shared_ptr< Property > mPropThetaIp        = nullptr;
        std::shared_ptr< Property > mPropThetaOp        = nullptr;
        std::shared_ptr< Property > mPropAspectRatio    = nullptr;

        // matrices used in the homogenization process
        Matrix< DDRMat > mEshelbyTensor;
        Matrix< DDRMat > mConstMatrix;
        Matrix< DDRMat > mConstFiber;
        Matrix< DDRMat > mRotation;
        Matrix< DDRMat > mRotationDerInPlane;
        Matrix< DDRMat > mRotationDerOutPlane;
        Matrix< DDRMat > mConstPrime;

        // property type for CM
        enum class CM_Property_Type_MT
        {
            EMOD1,
            NU1,
            EMOD2,
            NU2,
            VF,
            THETA_IP,
            THETA_OP,
            AR,
            MAX_ENUM
        };

        //--------------------------------------------------------------------------------------------------------------

      public:
        /*
         * trivial constructor
         */
        CM_Struc_Linear_MoriTanaka();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~CM_Struc_Linear_MoriTanaka() override{};

        //------------------------------------------------------------------------------
        /*
         * @return constitutive_type
         */
        Constitutive_Type
        get_constitutive_type() const override
        {
            return Constitutive_Type::STRUC_LIN_MT;
        }

        //------------------------------------------------------------------------------
        /**
         * set local properties
         */
        void set_local_properties() override;

      private:
        /**
         * evaluate the constitutive model matrix
         */
        void eval_const() override;

        //--------------------------------------------------------------------------------------------------------------

        /**
         * @brief computes the plane stress constitutive tensor
         *
         * @param tParams
         */

        void full_plane_stress(
                std::initializer_list< const real >&& tParams ) override;

        //--------------------------------------------------------------------------------------------------------------

        /**
         * @brief computes the 3d constitutive tensor
         *
         * @param tParams
         */
        void full_3d(
                std::initializer_list< const real >&& tParams ) override;

        //--------------------------------------------------------------------------------------------------------------

        /**
         * @brief Set the isotopic tensor for matrix and fiber
         *
         * @param aEmod
         * @param aNu
         * @param aIsotropicConst
         */

        void
        set_isotopic_tensor( const real& aEmod,
                const real&              aNu,
                Matrix< DDRMat >&        aIsotropicConst );

        //--------------------------------------------------------------------------------------------------------------

        /**
         * @brief Set the eshelby tensor object
         *
         * @param AR
         * @param v
         */

        void
        set_eshelby_tensor( real const & aAspectRatio, real const & aNu );

        //--------------------------------------------------------------------------------------------------------------
        /**
         *  Sets the appropriate function pointers
         *  based on the current member data
         *  for spatial dimensions and model types
         */
        void set_function_pointers() override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model flux derivative wrt to a dof type
         *
         * @param[ in ] aDofTypes  dof type wrt which the derivative is evaluated
         */
        void eval_dFluxdDOF( const Vector< MSI::Dof_Type >& aDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model test traction derivative wrt to a dof type
         *
         * @param[ in ] aDofTypes      dof type wrt which the derivative is evaluated
         * @param[ in ] aNormal        surface normal
         * @param[ in ] aJump          displacement jump
         * @param[ in ] aTestDofTypes  dof type of test function
         */
        void eval_dTestTractiondDOF(
                const Vector< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&        aNormal,
                const Matrix< DDRMat >&        aJump,
                const Vector< MSI::Dof_Type >& aTestDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the rotation tensor wrt in plane and out of plane angeles
         *
         */
        void eval_inplane_rotation_derivative();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the rotation tensor wrt in plane and out of plane angeles
         *
         */
        void eval_outplane_rotation_derivative();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model derivative wrt to a dof type
         *
         * @param[ in ] aDofTypes  dof type wrt which the derivative is evaluated
         */
        void eval_dConstdDOF( const Vector< MSI::Dof_Type >& aDofTypes ) override;
    };
    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_CM_STRUC_LINEAR_ISOTROPIC_HPP_ */


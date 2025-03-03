/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Struc_Linear_Isotropic.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_STRUC_LINEAR_ISOTROPIC_HPP_
#define SRC_FEM_CL_FEM_CM_STRUC_LINEAR_ISOTROPIC_HPP_

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_CM_Struc_Linear.hpp"

namespace moris::fem
{
    //--------------------------------------------------------------------------------------------------------------

    class CM_Struc_Linear_Isotropic : public CM_Struc_Linear
    {

      protected:
        // default local properties
        std::shared_ptr< Property > mPropEMod    = nullptr;
        std::shared_ptr< Property > mPropPoisson = nullptr;

      private:
        // property type for CM
        enum class CM_Property_Type_Iso
        {
            EMOD,
            NU,
            MAX_ENUM,
        };

        //--------------------------------------------------------------------------------------------------------------

      public:
        /*
         * trivial constructor
         */
        CM_Struc_Linear_Isotropic();

        //--------------------------------------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~CM_Struc_Linear_Isotropic() override{};

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
        real eval_inv_bulk_modulus() override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the inverse of the bulk modulus, K
         * wrt dof type
         * @param[ in ] aDofTypes            a dof type wrt which the derivative is evaluated
         * @param[ out ] dInvBulkModulusdDOF derivative of K
         */
        Matrix< DDRMat > eval_dInvBulkModulusdDOF( const Vector< MSI::Dof_Type >& aDofTypes ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * returns the E prime values used in the computation of the Stress Intensity Factor(s)
         */

        real get_e_prime() override;

        //--------------------------------------------------------------------------------------------------------------

        Constitutive_Type
        get_constitutive_type() const override
        {
            return Constitutive_Type::STRUC_LIN_ISO;
        }

      protected:
        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the constitutive model flux derivative wrt to a dof type
         * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
         */
        void eval_dFluxdDOF( const Vector< MSI::Dof_Type >& aDofTypes ) override;

      private:
        //--------------------------------------------------------------------------------------------------------------
        /**
         * evaluate the inverse of the bulk modulus, K
         * @param[ in ]  aNu             Poisson ratio
         * @param[ in ]  aEMod           Elasticity modulus
         * @param[ out ] aInvBulkModulus 1/K
         */
        void eval_inv_bulk_modulus_generic(
                const real& aNu,
                const real& aEMod,
                real&       aInvBulkModulus ) override;

        void eval_inv_bulk_modulus_plane_stress(
                const real& aNu,
                const real& aEMod,
                real&       aInvBulkModulus ) override;

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
         * evaluate the constitutive model matrix
         */
        void eval_const() override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Full plane stress tensor
         * @param[ in ] aEmod Elastic modulus
         * @param[ in ] aNu   Poisson ratio
         */
        void full_plane_stress(
                std::initializer_list< const real >&& tParams ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Deviatoric plane stress tensor
         * @param aEmod Elastic modulus
         * @param aNu   Poisson ratio
         */
        void deviatoric_plane_stress(
                std::initializer_list< const real >&& tParams ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Full plane strain tensor
         * @param aEmod Elastic modulus
         * @param aNu   Poisson ratio
         */
        void full_plane_strain(
                std::initializer_list< const real >&& tParams ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Deviatoric plane strain tensor
         * @param aEmod Elastic modulus
         * @param aNu Poisson ratio
         */
        void deviatoric_plane_strain(
                std::initializer_list< const real >&& tParams ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Full axisymmetric tensor
         * @param[ in ] aEmod Elastic modulus
         * @param[ in ] aNu   Poisson ratio
         */
        void full_axisymmetric(
                std::initializer_list< const real >&& tParams ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Deviatoric axisymmetric tensor
         * @param aEmod Elastic modulus
         * @param aNu   Poisson ratio
         */
        void deviatoric_axisymmetric(
                std::initializer_list< const real >&& tParams ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Full 3d linear isotropic tensor
         * @param aEmod Elastic modulus
         * @param aNu Poisson ratio
         */
        void full_3d(
                std::initializer_list< const real >&& tParams ) override;

        //--------------------------------------------------------------------------------------------------------------
        /**
         * Deviatoric 3d tensor
         * @param aEmod Elastic modulus
         * @param aNu Poisson ratio
         */
        void deviatoric_3d(
                std::initializer_list< const real >&& tParams ) override;

        //--------------------------------------------------------------------------------------------------------------
    };
    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_CM_STRUC_LINEAR_ISOTROPIC_HPP_ */


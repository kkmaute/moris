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
#include "cl_Cell.hpp"     //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_CM_Struc_Linear.hpp"

namespace moris
{
    namespace fem
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
            ~CM_Struc_Linear_Isotropic(){};

            //------------------------------------------------------------------------------
            /**
             * set local properties
             */
            virtual void set_local_properties();

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
            virtual Matrix< DDRMat > eval_dInvBulkModulusdDOF( const Cell< MSI::Dof_Type >& aDofTypes );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * returns the E prime values used in the computation of the Stress Intensity Factor(s)
             */

            real get_e_prime();

            //--------------------------------------------------------------------------------------------------------------

            Constitutive_Type
            get_constitutive_type() const
            {
                return Constitutive_Type::STRUC_LIN_ISO;
            }

          protected:
            //--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux derivative wrt to a dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             */
            virtual void eval_dFluxdDOF( const Cell< MSI::Dof_Type >& aDofTypes );

          private:
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
             * evaluate the constitutive model test traction derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * @param[ in ] aNormal   normal
             */
            virtual void eval_dTestTractiondDOF(
                    const Cell< MSI::Dof_Type >& aDofTypes,
                    const Matrix< DDRMat >&      aNormal,
                    const Matrix< DDRMat >&      aJump,
                    const Cell< MSI::Dof_Type >& aTestDofTypes );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model matrix
             */
            virtual void eval_const() override;

            //--------------------------------------------------------------------------------------------------------------
            /**
             * Full plane stress tensor
             * @param[ in ] aEmod Elastic modulus
             * @param[ in ] aNu   Poisson ratio
             */
            virtual void full_plane_stress(
                    std::initializer_list< const real >&& tParams ) override;

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
            virtual void deviatoric_axisymmetric(
                    std::initializer_list< const real >&& tParams );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * Full 3d linear isotropic tensor
             * @param aEmod Elastic modulus
             * @param aNu Poisson ratio
             */
            virtual void full_3d(
                    std::initializer_list< const real >&& tParams );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * Deviatoric 3d tensor
             * @param aEmod Elastic modulus
             * @param aNu Poisson ratio
             */
            virtual void deviatoric_3d(
                    std::initializer_list< const real >&& tParams );

            //--------------------------------------------------------------------------------------------------------------
        };
        //--------------------------------------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_STRUC_LINEAR_ISOTROPIC_HPP_ */


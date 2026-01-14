/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Fluid_Incompressible.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_HPP_
#define SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_HPP_

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

    class CM_Fluid_Incompressible : public Constitutive_Model
    {
        protected:
            // default properties
            std::shared_ptr< Property > mPropViscosity = nullptr;
            std::shared_ptr< Property > mPropDensity   = nullptr;

            // default dof type
            MSI::Dof_Type mDofVelocity = MSI::Dof_Type::VX;
            MSI::Dof_Type mDofPressure = MSI::Dof_Type::P;

            // flags for additional strain related evaluation
            moris::Matrix< DDBMat > mdStraindxduEval;

            // storage for additional strain related evaluation
            Vector< Vector< Matrix< DDRMat > > > mdStraindxdu;

        private:
            // property type for CM
            enum class CM_Property_Type
            {
                DENSITY,      // fluid density
                VISCOSITY,    // fluid dynamic viscosity
                MAX_ENUM
            };

            // function pointers
            void ( CM_Fluid_Incompressible:: *m_eval_strain )()                 = nullptr;
            void ( CM_Fluid_Incompressible:: *m_eval_divstrain )()              = nullptr;
            void ( CM_Fluid_Incompressible:: *m_eval_teststrain )()             = nullptr;
            void ( CM_Fluid_Incompressible:: *m_eval_dstraindx )( uint aOrder ) = nullptr;
            void ( CM_Fluid_Incompressible:: *m_eval_ddivstraindu )(
                    const Vector< MSI::Dof_Type > &aDofTypes ) = nullptr;
            void ( CM_Fluid_Incompressible:: *m_flatten_normal )(
                    const Matrix< DDRMat > &aNormal,
                    Matrix< DDRMat >       &aFlatNormal ) = nullptr;
            void ( CM_Fluid_Incompressible::*m_eval_dstraindxdu )(
                    const Vector< MSI::Dof_Type > &aDofTypes,
                    uint                           aOrder,
                    const Matrix< DDRMat >        &aJump ) = nullptr;

            //--------------------------------------------------------------------------------------------------------------

        public:
            //--------------------------------------------------------------------------------------------------------------
            /*
             * constructor
             */
            CM_Fluid_Incompressible();

            //--------------------------------------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~CM_Fluid_Incompressible() override {};

            //------------------------------------------------------------------------------
            /*
             * @return constitutive_type
             */
            Constitutive_Type
            get_constitutive_type() const override
            {
                return Constitutive_Type::FLUID_INCOMPRESSIBLE;
            }

            //--------------------------------------------------------------------------------------------------------------
            /**
             * set space dim
             */
            void set_space_dim( uint aSpaceDim ) override
            {
                mSpaceDim = aSpaceDim;
                this->set_function_pointers();
            }

            //--------------------------------------------------------------------------------------------------------------
            /**
             * set function pointers
             */
            void set_function_pointers();

            //------------------------------------------------------------------------------
            /**
             * set constitutive model dof types
             * @param[ in ] aDofTypes a list of group of dof types
             * @param[ in ] aDofStrings a list of strings to describe the dof types
             */
            void set_dof_type_list(
                    const Vector< Vector< MSI::Dof_Type > > &aDofTypes,
                    const Vector< std::string >             &aDofStrings ) override;

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
             * evaluate the strain template
             */
            void eval_strain() override
            {
                ( this->*m_eval_strain )();
            }
            void eval_strain_2d();
            void eval_strain_3d();

            //--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the divergence of the strain
             */
            void eval_divstrain() override
            {
                ( this->*m_eval_divstrain )();
            }
            void eval_divstrain_2d();
            void eval_divstrain_3d();

            //--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the divergence of the strain wrt dof type
             */
            void eval_ddivstraindu( const Vector< MSI::Dof_Type > &aDofTypes ) override
            {
                ( this->*m_eval_ddivstraindu )( aDofTypes );
            }
            void eval_ddivstraindu_2d( const Vector< MSI::Dof_Type > &aDofTypes );
            void eval_ddivstraindu_3d( const Vector< MSI::Dof_Type > &aDofTypes );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the strain wrt space
             * @param[ in ] aOrder order of the derivative
             */
            void eval_dstraindx( uint aOrder ) override
            {
                ( this->*m_eval_dstraindx )( aOrder );
            }
            void eval_dstraindx_2d( uint aOrder );
            void eval_dstraindx_3d( uint aOrder );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * get the derivative of the strain wrt space and dof type
             * @param[ in ] aOrder order of the derivative
             * @param[ in ] aDofTypes vector of dof type for derivative
             *
             */
            const Matrix< DDRMat > &dstraindxdu(
                    const Vector< MSI::Dof_Type > &aDofTypes,
                    uint                           aOrder,
                    const Matrix< DDRMat >        &aJump,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the strain wrt space and dof type
             * @param[ in ] aOrder    order of the derivative
             * @param[ in ] aDofTypes vector of dof type for derivative
             */
            void eval_dstraindxdu(
                    const Vector< MSI::Dof_Type > &aDofTypes,
                    uint                           aOrder,
                    const Matrix< DDRMat >        &aJump )
            {
                ( this->*m_eval_dstraindxdu )( aDofTypes, aOrder, aJump );
            }
            void eval_dstraindxdu_2d( const Vector< MSI::Dof_Type > &aDofTypes, uint aOrder, const Matrix< DDRMat > &aJump );
            void eval_dstraindxdu_3d( const Vector< MSI::Dof_Type > &aDofTypes, uint aOrder, const Matrix< DDRMat > &aJump );

            //--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the test strain
             */
            void eval_testStrain() override
            {
                ( this->*m_eval_teststrain )();
            }
            void eval_teststrain_2d();
            void eval_teststrain_3d();

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
             * evaluate the constitutive model strain derivative wrt to a dof type
             * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
             * @param[ in ] adStraindDOF a matrix to fill with derivative evaluation
             */
            void eval_dStraindDOF( const Vector< MSI::Dof_Type > &aDofTypes ) override;

            //--------------------------------------------------------------------------------------------------------------
            /**
             * flatten normal vector
             * @param[ in ] aNormal          a normal vector
             * @param[ in ] aFlattenedNormal a matrix for flattened normal to fill
             */
            void flatten_normal( const Matrix< DDRMat > &aNormal,
                    Matrix< DDRMat >                    &aFlatNormal )
            {
                ( this->*m_flatten_normal )( aNormal, aFlatNormal );
            }
            void flatten_normal_2d(
                    const Matrix< DDRMat > &aNormal,
                    Matrix< DDRMat >       &aFlatNormal );
            void flatten_normal_3d(
                    const Matrix< DDRMat > &aNormal,
                    Matrix< DDRMat >       &aFlatNormal );

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

            //--------------------------------------------------------------------------------------------------------------
    };

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_HPP_ */

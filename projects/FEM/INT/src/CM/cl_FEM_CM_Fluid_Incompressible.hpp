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

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                      //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {

        //--------------------------------------------------------------------------------------------------------------

        class CM_Fluid_Incompressible : public Constitutive_Model
        {
            protected:

                // default properties
                std::shared_ptr< Property > mPropViscosity = nullptr;
                std::shared_ptr< Property > mPropDensity   = nullptr;

            private:

                // default dof type
                MSI::Dof_Type mDofVelocity = MSI::Dof_Type::VX;
                MSI::Dof_Type mDofPressure = MSI::Dof_Type::P;

                // property type for CM
                enum class CM_Property_Type
                {
                    DENSITY,   // fluid density
                    VISCOSITY, // fluid viscosity
                    MAX_ENUM
                };

                // function pointers
                void ( CM_Fluid_Incompressible:: * m_eval_strain )() = nullptr;
                void ( CM_Fluid_Incompressible:: * m_eval_divstrain )() = nullptr;
                void ( CM_Fluid_Incompressible:: * m_eval_teststrain )() = nullptr;
                void ( CM_Fluid_Incompressible:: * m_eval_dstraindx )( uint aOrder ) = nullptr;
                void ( CM_Fluid_Incompressible:: * m_eval_ddivstraindu )(
                        const Vector< MSI::Dof_Type > & aDofTypes ) = nullptr;
                void ( CM_Fluid_Incompressible:: * m_flatten_normal )(
                        const Matrix< DDRMat > & aNormal,
                        Matrix< DDRMat >       & aFlatNormal ) = nullptr;

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
                ~CM_Fluid_Incompressible(){};

                //------------------------------------------------------------------------------
                /*
                 * @return constitutive_type
                 */
                Constitutive_Type
                get_constitutive_type() const
                {
                    return Constitutive_Type::FLUID_INCOMPRESSIBLE;
                }

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * set space dim
                 */
                void set_space_dim( uint aSpaceDim )
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
                        Vector< Vector< MSI::Dof_Type > > aDofTypes,
                        Vector< std::string >                  aDofStrings );

                //------------------------------------------------------------------------------
                /**
                 * set constitutive model dv types
                 * @param[ in ] aDvTypes a list of group of dv types
                 * @param[ in ] aDvStrings a list of strings to describe the dv types
                 */
                void set_dv_type_list(
                        Vector< Vector< PDV_Type > > aDvTypes,
                        Vector< std::string >             aDvStrings )
                {
                    Constitutive_Model::set_dv_type_list( aDvTypes );
                }

                //------------------------------------------------------------------------------
                /**
                 * set local properties
                 */
                void set_local_properties();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model flux
                 */
                void eval_flux();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the divergence of the flux
                 */
                void eval_divflux();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the divergence of the flux wrt dof type
                 */
                void eval_ddivfluxdu( const Vector< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the flux wrt space
                 * @param[ in ] aOrder order of the derivative
                 */
                void eval_dfluxdx( uint aOrder );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the traction
                 * @param[ in ] aNormal normal
                 */
                void eval_traction( const Matrix< DDRMat > & aNormal );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the test traction
                 * @param[ in ] aNormal   normal
                 */
                void eval_testTraction(
                        const Matrix< DDRMat >             & aNormal,
                        const Vector< MSI::Dof_Type > & aTestDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the strain template
                 */
                void eval_strain()
                {
                    ( this->*m_eval_strain )();
                }
                void eval_strain_2d();
                void eval_strain_3d();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the divergence of the strain
                 */
                void eval_divstrain()
                {
                    ( this->*m_eval_divstrain )();
                }
                void eval_divstrain_2d();
                void eval_divstrain_3d();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the divergence of the strain wrt dof type
                 */
                void eval_ddivstraindu( const Vector< MSI::Dof_Type > & aDofTypes )
                {
                    ( this->*m_eval_ddivstraindu)( aDofTypes );
                }
                void eval_ddivstraindu_2d( const Vector< MSI::Dof_Type > & aDofTypes );
                void eval_ddivstraindu_3d( const Vector< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the strain wrt space
                 * @param[ in ] aOrder order of the derivative
                 */
                void eval_dstraindx( uint aOrder )
                {
                    ( this->*m_eval_dstraindx )( aOrder );
                }
                void eval_dstraindx_2d( uint aOrder );
                void eval_dstraindx_3d( uint aOrder );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the test strain
                 */
                void eval_testStrain()
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
                void eval_dFluxdDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 */
                void eval_dTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model test traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 */
                void eval_dTestTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal,
                        const Matrix< DDRMat >             & aJump,
                        const Vector< MSI::Dof_Type > & aTestDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model strain derivative wrt to a dof type
                 * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
                 * @param[ in ] adStraindDOF a matrix to fill with derivative evaluation
                 */
                void eval_dStraindDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * flatten normal vector
                 * @param[ in ] aNormal          a normal vector
                 * @param[ in ] aFlattenedNormal a matrix for flattened normal to fill
                 */
                void flatten_normal( const Matrix< DDRMat > & aNormal,
                        Matrix< DDRMat > & aFlatNormal )
                {
                    ( this->*m_flatten_normal )( aNormal, aFlatNormal );
                }
                void flatten_normal_2d(
                        const Matrix< DDRMat > & aNormal,
                        Matrix< DDRMat >       & aFlatNormal );
                void flatten_normal_3d(
                        const Matrix< DDRMat > & aNormal,
                        Matrix< DDRMat >       & aFlatNormal );

                //--------------------------------------------------------------------------------------------------------------
            private:

                //--------------------------------------------------------------------------------------------------------------
        };

        //--------------------------------------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_HPP_ */


/*
 * cl_FEM_CM_Fluid_Compressible_Ideal.hpp
 *
 *  Created on: Jul 24, 2020
 *  Author: wunsch
 */

#ifndef SRC_FEM_CL_FEM_CM_FLUID_COMPRESSIBLE_IDEAL_HPP_
#define SRC_FEM_CL_FEM_CM_FLUID_COMPRESSIBLE_IDEAL_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {

        //--------------------------------------------------------------------------------------

        class CM_Fluid_Compressible_Ideal : public Constitutive_Model
        {

                //------------------------------------------------------------------------------

            private:

                // Pressure
                Matrix< DDRMat > mPressure;
                Matrix< DDRMat > mPressureDof;

                // Thermal Flux
                Matrix< DDRMat > mThermalFlux;
                Matrix< DDRMat > mThermalFluxDof;

                // default dof types
                MSI::Dof_Type mDofDensity     = MSI::Dof_Type::RHO;
                MSI::Dof_Type mDofVelocity    = MSI::Dof_Type::VX;
                MSI::Dof_Type mDofTemperature = MSI::Dof_Type::TEMP;

                // flattened identity matrix
                Matrix< DDRMat > mFlatIdentity;

                // velocity matrix for flattened tensors
                Matrix< DDRMat > mVelocityMatrix;

                // property type for CM
                enum class CM_Property_Type
                {
                        ISOCHORIC_HEAT_CAPACITY,  // heat capacity at constant density
                        SPECIFIC_GAS_CONSTANT,    // specific gas constant for fluid (R = R_0/M)
                        DYNAMIC_VISCOSITY,        // dynamic viscosity
                        THERMAL_CONDUCTIVITY,     // thermal conductivity
                        MAX_ENUM
                };

                // local string to property enum map
                std::map< std::string, CM_Property_Type > mPropertyMap;

                // function pointer for functions depending spatial dimension
                void ( CM_Fluid_Compressible_Ideal:: * m_eval_strain )() = nullptr;
                void ( CM_Fluid_Compressible_Ideal:: * m_eval_teststrain )() = nullptr;
                void ( CM_Fluid_Compressible_Ideal:: * m_eval_velocitymatrix )() = nullptr;
                void ( CM_Fluid_Compressible_Ideal:: * m_unfold_tensor )(
                        const Matrix< DDRMat > & aFlattenedTensor,
                        Matrix< DDRMat > & aExpandedTensor) = nullptr;

                //------------------------------------------------------------------------------

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                CM_Fluid_Compressible_Ideal();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~CM_Fluid_Compressible_Ideal(){};

                //------------------------------------------------------------------------------
                /**
                 * set space dim
                 */
                void set_space_dim( uint aSpaceDim )
                {
                    mSpaceDim = aSpaceDim;
                    this->set_function_pointers();
                }

                //------------------------------------------------------------------------------
                /**
                 * set function pointers for 2D and 3D
                 */
                void set_function_pointers();

                //------------------------------------------------------------------------------
                /**
                 * set constitutive model dof types
                 * @param[ in ] aDofTypes a list of group of dof types
                 * @param[ in ] aDofStrings a list of strings to describe the dof types
                 */
                void set_dof_type_list(
                        moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                        moris::Cell< std::string >                  aDofStrings );

                //------------------------------------------------------------------------------
                /**
                 * set constitutive model dv types
                 * @param[ in ] aDvTypes a list of group of dv types
                 * @param[ in ] aDvStrings a list of strings to describe the dv types
                 */
                void set_dv_type_list(
                        moris::Cell< moris::Cell< PDV_Type > > aDvTypes,
                        moris::Cell< std::string >             aDvStrings )
                {
                    Constitutive_Model::set_dv_type_list( aDvTypes );
                }

                //------------------------------------------------------------------------------
                /**
                 * set a property pointer
                 * @param[ in ] aProperty     a property pointer
                 * @param[ in ] aPropertyType a string defining the property
                 */
                void set_property(
                        std::shared_ptr< fem::Property > aProperty,
                        std::string                      aPropertyString );

                //------------------------------------------------------------------------------
                /**
                 * get a property pointer
                 * @param[ in ]  aPropertyType a string defining the property
                 * @param[ out ] aProperty     a property pointer
                 */
                std::shared_ptr< Property > get_property( std::string aPropertyString );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model flux (in this case: the thermal flux)
                 */
                void eval_flux();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model flux derivative wrt to a dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermal part of the constitutive model flux
                 */
                void eval_thermal_flux();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermal part of the constitutive model flux derivative wrt to a dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_thermal_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluates the constitutive model energy
                 */
                void eval_Energy();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model energy change rate wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dEnergyDotdDOF ( 1 x numDerDof )
                 */
                void eval_dEnergydDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluates the constitutive model change rate of energy
                 */
                void eval_EnergyDot();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model enthalpy change rate wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dEnergyDotdDOF ( 1 x numDerDof )
                 */
                void eval_dEnergyDotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model stress
                 */
                void eval_stress();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model stress derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 */
                void eval_dStressdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the traction
                 * @param[ in ] aNormal normal
                 */
                void eval_traction( const Matrix< DDRMat > & aNormal );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 */
                void eval_dTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal );

                //------------------------------------------------------------------------------
                /**
                 * evaluate and get the constitutive model / equation of state pressure
                 * @param[ out ] constitutive model / equation of state pressure
                 */
                const Matrix< DDRMat > & pressure();

                //------------------------------------------------------------------------------
                /**
                 * get the derivative of the Pressure wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] PressureDofDer derivative of the pressure wrt dof
                 */
                const Matrix< DDRMat > & dPressuredDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluates the test strain (is equal to dStrain/dDoF in this case)
                 */
                void eval_testStrain()
                {
                    ( this->*m_eval_teststrain )();
                }
                void eval_teststrain_2d();
                void eval_teststrain_3d();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluates the strain for 2D and 3D
                 */
                void eval_strain()
                {
                    ( this->*m_eval_strain )();
                }
                void eval_strain_2d();
                void eval_strain_3d();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluates velocity matrix for operation on flattended stress and strain vectors
                 */
                void eval_velocityMatrix()
                {
                    ( this->*m_eval_velocitymatrix )();
                }
                void eval_velocitymatrix_2d();
                void eval_velocitymatrix_3d();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluates the constitutive model strain derivative wrt to a dof type
                 * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
                 * @param[ in ] adStraindDOF a matrix to fill with derivative evaluation
                 */
                void eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * unfolds the flattened 1st order strain and stress tensors into their 2nd order form
                 * @param[ in ] aFlattenedTensor a matrix containing the flattened tensor to be expanded
                 * @param[ in ] aExpandedTensor  a matrix to fill with expanded form of the tensor
                 */
                void unfold(
                        const Matrix< DDRMat > & aFlattenedTensor,
                        Matrix< DDRMat > & aExpandedTensor )
                {
                    ( this->*m_unfold_tensor )( aFlattenedTensor, aExpandedTensor );
                }
                void unfold_2d(
                        const Matrix< DDRMat > & aFlattenedTensor,
                        Matrix< DDRMat > & aExpandedTensor);
                void unfold_3d(
                        const Matrix< DDRMat > & aFlattenedTensor,
                        Matrix< DDRMat > & aExpandedTensor);

                //------------------------------------------------------------------------------
            private:

                //------------------------------------------------------------------------------
        };

        //--------------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_FLUID_COMPRESSIBLE_IDEAL_HPP_ */

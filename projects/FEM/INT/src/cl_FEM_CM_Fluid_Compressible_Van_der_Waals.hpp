/*
 * cl_FEM_CM_Fluid_Compressible_Van_der_Waals.hpp
 *
 *  Created on: Jul 25, 2020
 *  Author: wunsch
 */

#ifndef SRC_FEM_CL_FEM_CM_FLUID_COMPRESSIBLE_VAN_DER_WAALS_HPP_
#define SRC_FEM_CL_FEM_CM_FLUID_COMPRESSIBLE_VAN_DER_WAALS_HPP_

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

        class CM_Fluid_Compressible_Van_der_Waals : public Constitutive_Model
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

                // strain contribution due to density surface effects
                Matrix< DDRMat > mDensityStrain;
                Matrix< DDRMat > mDensityStrainDof;

                // laplacian of density
                Matrix< DDRMat > mLaplaceDensity;
                Matrix< DDRMat > mLaplaceDensityDof;

                // property type for CM
                enum class CM_Property_Type
                {
                        ISOCHORIC_HEAT_CAPACITY,  // heat capacity at constant density
                        SPECIFIC_GAS_CONSTANT,    // specific gas constant for fluid (R = R_0/M)
                        DYNAMIC_VISCOSITY,        // dynamic viscosity
                        THERMAL_CONDUCTIVITY,     // thermal conductivity
                        CAPILLARITY_COEFFICIENT,  // capillarity coefficient
                        FIRST_VDW_CONSTANT,       // first Van-der-Waals constant
                        SECOND_VDW_CONSTANT,      // second Van-der-Waals constant
                        MAX_ENUM
                };

                // local string to property enum map
                std::map< std::string, CM_Property_Type > mPropertyMap;

                // function pointer for functions depending spatial dimension
                void ( CM_Fluid_Compressible_Van_der_Waals:: * m_eval_strain )() = nullptr;
                void ( CM_Fluid_Compressible_Van_der_Waals:: * m_eval_teststrain )() = nullptr;
                void ( CM_Fluid_Compressible_Van_der_Waals:: * m_eval_densitystrain )() = nullptr;
                void ( CM_Fluid_Compressible_Van_der_Waals:: * m_eval_densitystraindof )() = nullptr;
                void ( CM_Fluid_Compressible_Van_der_Waals:: * m_eval_laplacedensity )() = nullptr;
                void ( CM_Fluid_Compressible_Van_der_Waals:: * m_eval_laplacedensitydof )() = nullptr;
                void ( CM_Fluid_Compressible_Van_der_Waals:: * m_eval_velocitymatrix )() = nullptr;
                void ( CM_Fluid_Compressible_Van_der_Waals:: * m_unfold_tensor )(
                        const Matrix< DDRMat > & aFlattenedTensor,
                        Matrix< DDRMat > & aExpandedTensor) = nullptr;

                //------------------------------------------------------------------------------

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                CM_Fluid_Compressible_Van_der_Waals();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~CM_Fluid_Compressible_Van_der_Waals(){};

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
                 * evaluates the strain contribution of the density for the stress tensor for 2D and 3D
                 */
                void eval_densityStrain()
                {
                    ( this->*m_eval_densitystrain )();
                }
                void eval_densitystrain_2d();
                void eval_densitystrain_3d();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluates the dof deriv of the strain contribution of the density for the stress tensor for 2D and 3D
                 */
                void eval_densityStrainDof()
                {
                    ( this->*m_eval_densitystraindof )();
                }
                void eval_densitystraindof_2d();
                void eval_densitystraindof_3d();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluates the laplacian of the density for 2D and 3D
                 */
                void eval_laplaceDensity()
                {
                    ( this->*m_eval_laplacedensity )();
                }
                void eval_laplacedensity_2d();
                void eval_laplacedensity_3d();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluates the dof deriv of the laplacian of the density for 2D and 3D
                 */
                void eval_laplaceDensityDof()
                {
                    ( this->*m_eval_laplacedensitydof )();
                }
                void eval_laplacedensitydof_2d();
                void eval_laplacedensitydof_3d();

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

#endif /* SRC_FEM_CL_FEM_CM_FLUID_COMPRESSIBLE_VAN_DER_WAALS_HPP_ */

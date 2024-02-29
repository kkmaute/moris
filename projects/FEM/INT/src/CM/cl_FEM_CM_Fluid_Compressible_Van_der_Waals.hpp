/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Fluid_Compressible_Van_der_Waals.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_FLUID_COMPRESSIBLE_VAN_DER_WAALS_HPP_
#define SRC_FEM_CL_FEM_CM_FLUID_COMPRESSIBLE_VAN_DER_WAALS_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                          //MRS/CNT/src

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
            protected:
                // default local properties
                std::shared_ptr< Property > mPropIsochoricHeatCapacity  = nullptr;
                std::shared_ptr< Property > mPropSpecificGasConstant    = nullptr;
                std::shared_ptr< Property > mPropDynamicViscosity       = nullptr;
                std::shared_ptr< Property > mPropThermalConductivity    = nullptr;
                std::shared_ptr< Property > mPropCapillarityCoefficient = nullptr;
                std::shared_ptr< Property > mPropFirstVdWconstant       = nullptr;
                std::shared_ptr< Property > mPropSecondVdWconstant      = nullptr;

            private:

                // Pressure
                Matrix< DDRMat > mPressure;
                bool mPressureEval = true;
                Vector< Matrix< DDRMat > > mPressureDof;
                moris::Matrix< DDBMat > mPressureDofEval;

                // Thermal Flux ----------------------------------------
                Matrix< DDRMat > mThermalFlux;
                bool mThermalFluxEval = true;
                Vector< Matrix< DDRMat > > mThermalFluxDof;
                moris::Matrix< DDBMat > mThermalFluxDofEval;

                // Work Flux
                Matrix< DDRMat > mWorkFlux;
                bool mWorkFluxEval = true;
                Vector< Matrix< DDRMat > > mWorkFluxDof;
                moris::Matrix< DDBMat > mWorkFluxDofEval;

                // Energy Flux
                Matrix< DDRMat > mEnergyFlux;
                bool mEnergyFluxEval = true;
                Vector< Matrix< DDRMat > > mEnergyFluxDof;
                moris::Matrix< DDBMat > mEnergyFluxDofEval;

//                // Mechanical Flux
//                Matrix< DDRMat > mStress;
//                bool mStressEval = true;
//                Vector< Matrix< DDRMat > > mStressDof;
//                moris::Matrix< DDBMat > mStressDofEval;

                // Thermal Traction ------------------------------------
                Matrix< DDRMat > mThermalTraction;
                bool mThermalTractionEval = true;
                Vector< Matrix< DDRMat > > mThermalTractionDof;
                moris::Matrix< DDBMat > mThermalTractionDofEval;

                // Work Traction
                Matrix< DDRMat > mWorkTraction;
                bool mWorkTractionEval = true;
                Vector< Matrix< DDRMat > > mWorkTractionDof;
                moris::Matrix< DDBMat > mWorkTractionDofEval;

                // Energy Traction
                Matrix< DDRMat > mEnergyTraction;
                bool mEnergyTractionEval = true;
                Vector< Matrix< DDRMat > > mEnergyTractionDof;
                moris::Matrix< DDBMat > mEnergyTractionDofEval;

                // Mechanical Traction
                Matrix< DDRMat > mMechanicalTraction;
                bool mMechanicalTractionEval = true;
                Vector< Matrix< DDRMat > > mMechanicalTractionDof;
                moris::Matrix< DDBMat > mMechanicalTractionDofEval;

                // Thermal Test Traction ------------------------------
                Vector< Matrix< DDRMat > > mThermalTestTraction;
                moris::Matrix< DDBMat > mThermalTestTractionEval;
                Vector< Vector< Matrix< DDRMat > > > mdThermalTestTractiondDof;
                moris::Matrix< DDBMat > mdThermalTestTractiondDofEval;

                // Mechanical Test Traction
                Vector< Matrix< DDRMat > > mMechanicalTestTraction;
                moris::Matrix< DDBMat > mMechanicalTestTractionEval;
                Vector< Vector< Matrix< DDRMat > > > mdMechanicalTestTractiondDof;
                moris::Matrix< DDBMat > mdMechanicalTestTractiondDofEval;

                // DoF derivative of du/dt ----------------------------
                Matrix< DDRMat > mdNveldt;
                bool mdNveldtEval = true;

                // velocity matrix for flattened tensors
                Matrix< DDRMat > mVelocityMatrix;
                bool mVelocityMatrixEval = true;

                // strain contribution due to density surface effects
                Matrix< DDRMat > mDensityStrain;
                bool mDensityStrainEval = true;
                Matrix< DDRMat > mDensityStrainDof;
                bool mDensityStrainDofEval = true;

                // laplacian of density
                Matrix< DDRMat > mLaplaceDensity;
                bool mLaplaceDensityEval = true;
                Matrix< DDRMat > mLaplaceDensityDof;
                bool mLaplaceDensityDofEval = true;

                // default dof types
                MSI::Dof_Type mDofDensity     = MSI::Dof_Type::RHO;
                MSI::Dof_Type mDofVelocity    = MSI::Dof_Type::VX;
                MSI::Dof_Type mDofTemperature = MSI::Dof_Type::TEMP;

                // flattened identity matrix
                Matrix< DDRMat > mFlatIdentity;

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
                void ( CM_Fluid_Compressible_Van_der_Waals:: * m_flatten_normal )(
                        const Matrix< DDRMat > & aNormal,
                        Matrix< DDRMat > & aFlatNormal ) = nullptr;

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
                /*
                 * @return constitutive_type
                 */
                Constitutive_Type
                get_constitutive_type() const
                {
                    return Constitutive_Type::FLUID_COMPRESSIBLE_VDW;
                }

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
                 * reset evaluation flags specific to this constitutive models
                 */
                void reset_specific_eval_flags();

                //------------------------------------------------------------------------------
                /**
                 * initialize storage variables and evaluation flags specific to this child CM
                 * function is called in the build_global_dof_type_list() in parent class
                 */
                void initialize_spec_storage_vars_and_eval_flags();

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
                        Vector< Vector< gen::PDV_Type > > aDvTypes,
                        Vector< std::string >             aDvStrings )
                {
                    Constitutive_Model::set_dv_type_list( aDvTypes );
                }

                //------------------------------------------------------------------------------
                /**
                 * set local properties
                 */
                void set_local_properties();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model fluxes
                 */
                void eval_flux()
                {
                    MORIS_ERROR( false,
                            "CM_Fluid_Compressible_Ideal::eval_flux - not implemented, use specific flux functions." );
                };

                /**
                 * get the constitutive model fluxes
                 * @param[ in ]  aCMFunctionType  enum indicating which flux is called, if there are several
                 * @param[ out ] mFlux constitutive model fluxes
                 */
                const Matrix< DDRMat > & flux( enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                /**
                 * evaluate the constitutive model fluxes derivatives wrt to a dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_dFluxdDOF(
                        const Vector< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false,
                            "CM_Fluid_Compressible_Ideal::eval_dFluxdDOF - not implemented, use specific flux functions." );
                };

                /**
                 * get the derivative of the fluxes wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] mFluxDofDer derivative of the fluxes wrt dof
                 */
                virtual const Matrix< DDRMat > & dFluxdDOF(
                        const Vector< MSI::Dof_Type > & aDofTypes,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model advective energy flux
                 */
                void eval_energy_flux();

                /**
                 *  get the constitutive model advective energy flux
                 * @param[ out ] constitutive model / equation of state pressure
                 */
                const Matrix< DDRMat > & energy_flux();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model advective energy flux derivative wrt to a dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_energy_dFluxdDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the derivative of the advective energy flux wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] ThermalFluxDofDer derivative of the thermal flux wrt dof
                 */
                const Matrix< DDRMat > & energy_dFluxdDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermal part of the constitutive model flux
                 */
                void eval_thermal_flux();

                /**
                 *  get the constitutive model thermal flux
                 * @param[ out ] constitutive model / equation of state pressure
                 */
                const Matrix< DDRMat > & thermal_flux();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the thermal part of the constitutive model flux derivative wrt to a dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_thermal_dFluxdDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the derivative of the thermal flux wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] ThermalFluxDofDer derivative of the thermal flux wrt dof
                 */
                const Matrix< DDRMat > & thermal_dFluxdDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the mechanical work part of the constitutive model flux
                 */
                void eval_work_flux();

                /**
                 *  get the constitutive model mechanical work flux
                 * @param[ out ] constitutive model / equation of state pressure
                 */
                const Matrix< DDRMat > & work_flux();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the mechanical work part of the constitutive model flux derivative wrt to a dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_work_dFluxdDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the derivative of the mechanical work flux wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] ThermalFluxDofDer derivative of the thermal flux wrt dof
                 */
                const Matrix< DDRMat > & work_dFluxdDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
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
                void eval_dStressdDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the traction
                 * @param[ in ] aNormal normal
                 */
                void eval_traction( const Matrix< DDRMat > & aNormal )
                {
                    MORIS_ERROR( false,
                            "CM_Fluid_Compressible_Ideal::eval_traction - not implemented, use specific traction functions." );
                };

                /**
                 * get the constitutive model traction
                 * @param[ in ]  aNormal   normal
                 * @param[ out ] mTraction constitutive model traction
                 */
                const Matrix< DDRMat > & traction(
                        const Matrix< DDRMat > & aNormal,
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 */
                void eval_dTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal )
                {
                    MORIS_ERROR( false,
                            "CM_Fluid_Compressible_Ideal::eval_dTractiondDOF - not implemented, use specific traction functions." );
                };

                /**
                 * get the derivative of the traction wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ in ]  aNormal         normal
                 * @param[ out ] mTractionDofDer derivative of the traction wrt dof
                 */
                const Matrix< DDRMat > & dTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofType,
                        const Matrix< DDRMat >             & aNormal,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model thermal flux traction
                 * @param[ in ] aNormal normal
                 */
                void eval_thermal_traction( const Matrix< DDRMat > & aNormal );

                /**
                 * get the constitutive model thermal flux traction
                 * @param[ in ]  aNormal   normal
                 * @param[ out ] mTraction constitutive model thermal flux traction
                 */
                const Matrix< DDRMat > & thermal_traction(
                        const Matrix< DDRMat > & aNormal);

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model thermal flux traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 */
                void eval_thermal_dTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal );

                /**
                 * get the derivative of the thermal flux traction wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ in ]  aNormal         normal
                 * @param[ out ] mTractionDofDer derivative of the thermal flux traction wrt dof
                 */
                const Matrix< DDRMat > & thermal_dTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofType,
                        const Matrix< DDRMat >             & aNormal );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model advective energy flux traction
                 * @param[ in ] aNormal normal
                 */
                void eval_energy_traction( const Matrix< DDRMat > & aNormal );

                /**
                 * get the constitutive model advective energy flux traction
                 * @param[ in ]  aNormal   normal
                 * @param[ out ] mTraction constitutive model advective energy flux traction
                 */
                const Matrix< DDRMat > & energy_traction(
                        const Matrix< DDRMat > & aNormal);

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model advective energy flux traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 */
                void eval_energy_dTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal );

                /**
                 * get the derivative of the advective energy flux traction wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ in ]  aNormal         normal
                 * @param[ out ] mTractionDofDer derivative of the traction wrt dof
                 */
                const Matrix< DDRMat > & energy_dTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofType,
                        const Matrix< DDRMat >             & aNormal );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model mechanical work flux traction
                 * @param[ in ] aNormal normal
                 */
                void eval_work_traction( const Matrix< DDRMat > & aNormal );

                /**
                 * get the constitutive model mechanical work flux traction
                 * @param[ in ]  aNormal   normal
                 * @param[ out ] mTraction constitutive model mechanical work flux traction
                 */
                const Matrix< DDRMat > & work_traction(
                        const Matrix< DDRMat > & aNormal);

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model mechanical work flux traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 */
                void eval_work_dTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal );

                /**
                 * get the derivative of the mechanical work flux traction wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ in ]  aNormal         normal
                 * @param[ out ] mTractionDofDer derivative of the mechanical work flux traction wrt dof
                 */
                const Matrix< DDRMat > & work_dTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofType,
                        const Matrix< DDRMat >             & aNormal );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model mechanical traction
                 * @param[ in ] aNormal normal
                 */
                void eval_mechanical_traction( const Matrix< DDRMat > & aNormal );

                /**
                 * get the constitutive model mechanical traction
                 * @param[ in ]  aNormal   normal
                 * @param[ out ] mTraction constitutive model mechanical traction
                 */
                const Matrix< DDRMat > & mechanical_traction(
                        const Matrix< DDRMat > & aNormal);

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model mechanical traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 */
                void eval_mechanical_dTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal );

                /**
                 * get the derivative of the mechanical traction wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ in ]  aNormal         normal
                 * @param[ out ] mTractionDofDer derivative of the mechanical traction wrt dof
                 */
                const Matrix< DDRMat > & mechanical_dTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofType,
                        const Matrix< DDRMat >             & aNormal );

                //------------------------------------------------------------------------------
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
                void eval_dEnergydDOF( const Vector< MSI::Dof_Type > & aDofTypes );

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
                void eval_dEnergyDotdDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * get the constitutive model test traction
                 * @param[ in ]  aNormal       normal
                 * @param[ in ]  aTestDofTypes      Dof type of the test functions
                 * @param[ out ] mTestTraction constitutive model test traction
                 */
                const Matrix< DDRMat > & testTraction(
                        const Matrix< DDRMat >             & aNormal,
                        const Vector< MSI::Dof_Type > & aTestDofTypes,
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT);

                /**
                 * get the derivative of the test traction wrt dof
                 * @param[ in ]  aDofType           group of dof type
                 * @param[ in ]  aNormal            normal
                 * @param[ in ]  aJump              jump in field values on boundary
                 * @param[ in ]  aTestDofTypes      Dof type of the test functions
                 * @param[ out ] mdTestTractiondDof derivative of the traction wrt dof
                 */
                const Matrix< DDRMat > & dTestTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofType,
                        const Matrix< DDRMat >             & aNormal,
                        const Matrix< DDRMat >             & aJump,
                        const Vector< MSI::Dof_Type > & aTestDofTypes,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                  * evaluate the constitutive model thermal test traction
                  * @param[ in ]  aNormal normal
                  * @param[ in ]  aTestDofTypes      Dof type of the test functions
                  */
                 void eval_thermal_testTraction(
                         const Matrix< DDRMat >             & aNormal,
                         const Vector< MSI::Dof_Type > & aTestDofTypes );

                 /**
                  * get the constitutive model thermal test traction
                  * @param[ in ]  aNormal       normal
                  * @param[ in ]  aTestDofTypes      Dof type of the test functions
                  * @param[ out ] mTestTraction constitutive model test traction
                  */
                 const Matrix< DDRMat > & thermal_testTraction(
                         const Matrix< DDRMat >             & aNormal,
                         const Vector< MSI::Dof_Type > & aTestDofTypes );

                 //------------------------------------------------------------------------------
                 /**
                  * evaluate the constitutive model thermal test traction derivative wrt to a dof type
                  * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
                  * @param[ in ] aNormal         normal
                  * @param[ in ] aJump           jump in field values on boundary
                  * @param[ in ] aTestDofTypes   Dof type of the test functions
                  */
                 void eval_thermal_dTestTractiondDOF(
                         const Vector< MSI::Dof_Type > & aDofTypes,
                         const Matrix< DDRMat >             & aNormal,
                         const Matrix< DDRMat >             & aJump,
                         const Vector< MSI::Dof_Type > & aTestDofTypes );

                 /**
                  * get the derivative of the thermal test traction wrt dof
                  * @param[ in ]  aDofType           group of dof type
                  * @param[ in ]  aNormal            normal
                  * @param[ in ]  aJump              jump in field values on boundary
                  * @param[ in ]  aTestDofTypes      Dof type of the test functions
                  * @param[ out ] mdTestTractiondDof derivative of the traction wrt dof
                  */
                 const Matrix< DDRMat > & thermal_dTestTractiondDOF(
                         const Vector< MSI::Dof_Type > & aDofType,
                         const Matrix< DDRMat >             & aNormal,
                         const Matrix< DDRMat >             & aJump,
                         const Vector< MSI::Dof_Type > & aTestDofTypes );

                 //------------------------------------------------------------------------------
                 //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model mechanical test traction
                 * @param[ in ]  aNormal normal
                 * @param[ in ]  aTestDofTypes      Dof type of the test functions
                 */
                void eval_mechanical_testTraction(
                        const Matrix< DDRMat >             & aNormal,
                        const Vector< MSI::Dof_Type > & aTestDofTypes );

                /**
                 * get the constitutive model mechanical test traction
                 * @param[ in ]  aNormal       normal
                 * @param[ in ]  aTestDofTypes      Dof type of the test functions
                 * @param[ out ] mTestTraction constitutive model test traction
                 */
                const Matrix< DDRMat > & mechanical_testTraction(
                        const Matrix< DDRMat >             & aNormal,
                        const Vector< MSI::Dof_Type > & aTestDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model mechanical test traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal         normal
                 * @param[ in ] aJump           jump in field values on boundary
                 * @param[ in ] aTestDofTypes   Dof type of the test functions
                 */
                void eval_mechanical_dTestTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal,
                        const Matrix< DDRMat >             & aJump,
                        const Vector< MSI::Dof_Type > & aTestDofTypes );

                /**
                 * get the derivative of the mechanical test traction wrt dof
                 * @param[ in ]  aDofType           group of dof type
                 * @param[ in ]  aNormal            normal
                 * @param[ in ]  aJump              jump in field values on boundary
                 * @param[ in ]  aTestDofTypes      Dof type of the test functions
                 * @param[ out ] mdTestTractiondDof derivative of the traction wrt dof
                 */
                const Matrix< DDRMat > & mechanical_dTestTractiondDOF(
                        const Vector< MSI::Dof_Type > & aDofType,
                        const Matrix< DDRMat >             & aNormal,
                        const Matrix< DDRMat >             & aJump,
                        const Vector< MSI::Dof_Type > & aTestDofTypes );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model pressure
                 */
                void eval_pressure();

                /**
                 * evaluate and get the constitutive model / equation of state pressure
                 * @param[ out ] constitutive model / equation of state pressure
                 */
                const Matrix< DDRMat > & pressure();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model pressure
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_dPressuredDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the derivative of the Pressure wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] PressureDofDer derivative of the pressure wrt dof
                 */
                const Matrix< DDRMat > & dPressuredDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluates the test strain (is equal to dStrain/dDoF in this case)
                 */
                void eval_testStrain()
                {
                    ( this->*m_eval_teststrain )();
                }
                void eval_teststrain_2d();
                void eval_teststrain_3d();

                /**
                 * evaluates the constitutive model strain derivative wrt to a dof type
                 * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
                 * @param[ in ] adStraindDOF a matrix to fill with derivative evaluation
                 */
                void eval_dStraindDOF( const Vector< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
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
                 * evaluates the DoF derivative formated for usage of the velocity time rate of change du/dt
                 * FIXME: The output format of the dnNdtn() function in the FI should probably be changed such that this is the
                 * standard format of this derivative
                 */
                void eval_dNveldt();

                /**
                 * get DoF derivative of the velocity time rate of change du/dt
                 */
                const Matrix< DDRMat > & dNveldt();

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluates velocity matrix for operation on flattended stress and strain vectors
                 */
                void eval_velocityMatrix()
                {
                    ( this->*m_eval_velocitymatrix )();
                }
                void eval_velocitymatrix_2d();
                void eval_velocitymatrix_3d();

                /**
                 * get the velocity matrix for operation on flattened stress and strain vectors
                 * @param[ out ] mVelocityMatrix constitutive model stress
                 */
                const Matrix< DDRMat > & velocityMatrix();

                //--------------------------------------------------------------------------------------------------------------
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

                /**
                 * get the strain contribution of the density for the stress tensor
                 * @param[ out ] mDensityStrain constitutive model density strain contribution
                 */
                const Matrix< DDRMat > & densityStrain();

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

                /**
                 * get the strain contribution of the density for the stress tensor, derivative wrt. the DoFs
                 * @param[ out ] mDensityStrainDoF constitutive model density strain contribution
                 */
                const Matrix< DDRMat > & densityStrainDof();

                //--------------------------------------------------------------------------------------------------------------
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

                /**
                 * get the laplacian of the density
                 * @param[ out ] mLaplaceDensity laplacian of the density
                 */
                const Matrix< DDRMat > & laplaceDensity();

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

                /**
                 * get the derivative of the laplacian of the density with respect to the DoFs
                 * @param[ out ] mLaplaceDensityDoF derivative of the laplacian of the density with respect to the DoFs
                 */
                const Matrix< DDRMat > & laplaceDensityDof();

                //--------------------------------------------------------------------------------------------------------------
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
                void flatten_normal_2d( const Matrix< DDRMat > & aNormal,
                        Matrix< DDRMat > & aFlatNormal );
                void flatten_normal_3d( const Matrix< DDRMat > & aNormal,
                        Matrix< DDRMat > & aFlatNormal );

                //------------------------------------------------------------------------------
            private:

                //------------------------------------------------------------------------------
        };

        //--------------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_FLUID_COMPRESSIBLE_VAN_DER_WAALS_HPP_ */


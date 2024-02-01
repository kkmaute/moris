/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Compressible_Newtonian_Fluid.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_COMPRESSIBLE_NEWTONIAN_FLUID_HPP_
#define SRC_FEM_CL_FEM_CM_COMPRESSIBLE_NEWTONIAN_FLUID_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {

        //--------------------------------------------------------------------------------------

        class CM_Compressible_Newtonian_Fluid : public Constitutive_Model
        {

                //------------------------------------------------------------------------------
            protected:
                // default local properties
                std::shared_ptr< Property > mPropDynamicViscosity       = nullptr;
                std::shared_ptr< Property > mPropThermalConductivity    = nullptr;

                // default thermodynamic material model
                std::shared_ptr< Material_Model > mMaterialModel  = nullptr;

            private:

                // Thermal Flux ---------------------------------------
                Matrix< DDRMat > mThermalFlux;
                bool mThermalFluxEval = true;
                moris::Cell< Matrix< DDRMat > > mThermalFluxDof;
                moris::Matrix< DDBMat > mThermalFluxDofEval;

                // Work Flux
                Matrix< DDRMat > mWorkFlux;
                bool mWorkFluxEval = true;
                moris::Cell< Matrix< DDRMat > > mWorkFluxDof;
                moris::Matrix< DDBMat > mWorkFluxDofEval;

                // Energy Flux
                Matrix< DDRMat > mEnergyFlux;
                bool mEnergyFluxEval = true;
                moris::Cell< Matrix< DDRMat > > mEnergyFluxDof;
                moris::Matrix< DDBMat > mEnergyFluxDofEval;

//                // Mechanical Flux
//                Matrix< DDRMat > mStress;
//                bool mStressEval = true;
//                moris::Cell< Matrix< DDRMat > > mStressDof;
//                moris::Cell< bool > mStressDofEval;

                // Thermal Div Flux -----------------------------------
                Matrix< DDRMat > mThermalDivFlux;
                bool mThermalDivFluxEval = true;
                moris::Cell< Matrix< DDRMat > > mThermalDivFluxDof;
                moris::Matrix< DDBMat > mThermalDivFluxDofEval;

                // Work Div Flux
                Matrix< DDRMat > mWorkDivFlux;
                bool mWorkDivFluxEval = true;
                moris::Cell< Matrix< DDRMat > > mWorkDivFluxDof;
                moris::Matrix< DDBMat > mWorkDivFluxDofEval;

                // Mechanical Div Flux (aka div-viscous-stress)
                Matrix< DDRMat > mMechanicalDivFlux;
                bool mMechanicalDivFluxEval = true;
                moris::Cell< Matrix< DDRMat > > mMechanicalDivFluxDof;
                moris::Matrix< DDBMat > mMechanicalDivFluxDofEval;

                // Thermal Traction -----------------------------------
                Matrix< DDRMat > mThermalTraction;
                bool mThermalTractionEval = true;
                moris::Cell< Matrix< DDRMat > > mThermalTractionDof;
                moris::Matrix< DDBMat > mThermalTractionDofEval;

                // Work Traction
                Matrix< DDRMat > mWorkTraction;
                bool mWorkTractionEval = true;
                moris::Cell< Matrix< DDRMat > > mWorkTractionDof;
                moris::Matrix< DDBMat > mWorkTractionDofEval;

                // Energy Traction
                Matrix< DDRMat > mEnergyTraction;
                bool mEnergyTractionEval = true;
                moris::Cell< Matrix< DDRMat > > mEnergyTractionDof;
                moris::Matrix< DDBMat > mEnergyTractionDofEval;

                // Mechanical Traction
                Matrix< DDRMat > mMechanicalTraction;
                bool mMechanicalTractionEval = true;
                moris::Cell< Matrix< DDRMat > > mMechanicalTractionDof;
                moris::Matrix< DDBMat > mMechanicalTractionDofEval;

                // Thermal Test Traction ------------------------------
                moris::Cell< Matrix< DDRMat > > mThermalTestTraction;
                moris::Matrix< DDBMat > mThermalTestTractionEval;
                moris::Cell< moris::Cell< Matrix< DDRMat > > > mdThermalTestTractiondDof;
                moris::Matrix< DDBMat > mdThermalTestTractiondDofEval;

                // Mechanical Test Traction
                moris::Cell< Matrix< DDRMat > > mMechanicalTestTraction;
                moris::Matrix< DDBMat > mMechanicalTestTractionEval;
                moris::Cell< moris::Cell< Matrix< DDRMat > > > mdMechanicalTestTractiondDof;
                moris::Matrix< DDBMat > mdMechanicalTestTractiondDofEval;

                // themal div-strain ----------------------------------
                Matrix< DDRMat > mThermalDivStrain;
                bool mThermalDivStrainEval = true;
                moris::Cell< Matrix< DDRMat > > mThermalDivStrainDof;
                moris::Matrix< DDBMat > mThermalDivStrainDofEval;

                // div-strain-rate
                Matrix< DDRMat > mDivStrainRate;
                bool mDivStrainRateEval = true;
                moris::Cell< Matrix< DDRMat > > mDivStrainRateDof;
                moris::Matrix< DDBMat > mDivStrainRateDofEval;

                // DoF derivative of du/dt ----------------------------
                Matrix< DDRMat > mdNveldt;
                bool mdNveldtEval = true;

                // div(div(u)*I)
                Matrix< DDRMat > mDivDivVel;
                bool mDivDivVelEval = true;
                moris::Cell< Matrix< DDRMat > > mDivDivVelDof;
                moris::Matrix< DDBMat > mDivDivVelDofEval;

                // velocity matrix for flattened tensors
                Matrix< DDRMat > mVelocityMatrix;
                bool mVelocityMatrixEval = true;

                // flattened identity matrix
                Matrix< DDRMat > mFlatIdentity;

                // multiplication matrices for condensed tensors
                const Matrix< DDRMat > mMultipMat2D = {
                        { 1.0, 0.0, 0.0 },
                        { 0.0, 1.0, 0.0 },
                        { 0.0, 0.0, 2.0 } };

                const Matrix< DDRMat > mMultipMat3D = {
                        { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
                        { 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 },
                        { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },
                        { 0.0, 0.0, 0.0, 2.0, 0.0, 0.0 },
                        { 0.0, 0.0, 0.0, 0.0, 2.0, 0.0 },
                        { 0.0, 0.0, 0.0, 0.0, 0.0, 2.0 } };

                // default dof types
                MSI::Dof_Type mDofDensity     = MSI::Dof_Type::RHO;
                MSI::Dof_Type mDofPressure    = MSI::Dof_Type::P;
                MSI::Dof_Type mDofVelocity    = MSI::Dof_Type::VX;
                MSI::Dof_Type mDofTemperature = MSI::Dof_Type::TEMP;

                // property type for CM
                enum class CM_Property_Type
                {
                        DYNAMIC_VISCOSITY,        // dynamic viscosity
                        THERMAL_CONDUCTIVITY,     // thermal conductivity
                        MAX_ENUM
                };

                // material model type for CM
                enum class MM_Type
                {
                        THERMODYNAMIC_MATERIAL_MODEL,  // thermodynamic material model for fluid
                        MAX_ENUM
                };

                // function pointer for functions depending spatial dimension
                void ( CM_Compressible_Newtonian_Fluid:: * m_eval_strain )() = nullptr;
                void ( CM_Compressible_Newtonian_Fluid:: * m_eval_teststrain )() = nullptr;
                void ( CM_Compressible_Newtonian_Fluid:: * m_eval_divstrainrate )() = nullptr;
                void ( CM_Compressible_Newtonian_Fluid:: * m_eval_ddivstrainratedu )( const moris::Cell< MSI::Dof_Type > & aDofTypes ) = nullptr;
                void ( CM_Compressible_Newtonian_Fluid:: * m_eval_divDivVel )() = nullptr;
                void ( CM_Compressible_Newtonian_Fluid:: * m_eval_dDivDivVeldu )( const moris::Cell< MSI::Dof_Type > & aDofTypes ) = nullptr;
                void ( CM_Compressible_Newtonian_Fluid:: * m_eval_velocitymatrix )() = nullptr;
                void ( CM_Compressible_Newtonian_Fluid:: * m_unfold_tensor )(
                        const Matrix< DDRMat > & aFlattenedTensor,
                        Matrix< DDRMat > & aExpandedTensor) = nullptr;
                void ( CM_Compressible_Newtonian_Fluid:: * m_flatten_normal )(
                        const Matrix< DDRMat > & aNormal,
                        Matrix< DDRMat > & aFlatNormal ) = nullptr;

                //------------------------------------------------------------------------------

            public:

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                CM_Compressible_Newtonian_Fluid();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~CM_Compressible_Newtonian_Fluid(){};

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
                        moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                        moris::Cell< std::string >                  aDofStrings );

                //------------------------------------------------------------------------------
                /**
                 * set constitutive model dv types
                 * @param[ in ] aDvTypes a list of group of dv types
                 * @param[ in ] aDvStrings a list of strings to describe the dv types
                 */
                void set_dv_type_list(
                        moris::Cell< moris::Cell< gen::PDV_Type > > aDvTypes,
                        moris::Cell< std::string >             aDvStrings )
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
                 * set thermodynamic material model
                 */
                void set_local_material_model();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model fluxes
                 */
                void eval_flux()
                {
                    MORIS_ERROR( false,
                            "CM_Compressible_Newtonian_Fluid::eval_flux - not implemented, use specific flux functions." );
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
                        const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false,
                            "CM_Compressible_Newtonian_Fluid::eval_dFluxdDOF - not implemented, use specific flux functions." );
                };

                /**
                 * get the derivative of the fluxes wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] mFluxDofDer derivative of the fluxes wrt dof
                 */
                virtual const Matrix< DDRMat > & dFluxdDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
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
                void eval_energy_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the derivative of the advective energy flux wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] ThermalFluxDofDer derivative of the thermal flux wrt dof
                 */
                const Matrix< DDRMat > & energy_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

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
                void eval_thermal_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the derivative of the thermal flux wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] ThermalFluxDofDer derivative of the thermal flux wrt dof
                 */
                const Matrix< DDRMat > & thermal_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

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
                void eval_work_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the derivative of the mechanical work flux wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] ThermalFluxDofDer derivative of the thermal flux wrt dof
                 */
                const Matrix< DDRMat > & work_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

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
                void eval_dStressdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model divergence of the fluxes
                 */
                void eval_divflux()
                {
                    MORIS_ERROR( false,
                            "CM_Compressible_Newtonian_Fluid::eval_divflux - not implemented, use specific flux functions." );
                };

                /**
                 * get the constitutive model diveregence of the fluxes
                 * @param[ in ]  aCMFunctionType  enum indicating which flux is called, if there are several
                 * @param[ out ] mDivFlux constitutive model fluxes
                 */
                const Matrix< DDRMat > & divflux( enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                /**
                 * evaluate the constitutive model divergence of the fluxes derivatives wrt to a dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_ddivfluxdu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false,
                            "CM_Compressible_Newtonian_Fluid::eval_ddivfluxdu - not implemented, use specific flux functions." );
                };

                /**
                 * get the derivative of the divergence of the fluxes wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] mDivFluxDofDer derivative of the fluxes wrt dof
                 */
                const Matrix< DDRMat > & ddivfluxdu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model divergence of the mechanical flux
                 */
                void eval_mechanical_divflux();

                /**
                 *  get the constitutive model divergence of the mechanical flux
                 * @param[ out ] mDivMechanicalFlux constitutive model divergence of the mechanical flux
                 */
                const Matrix< DDRMat > & mechanical_divflux();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model divergence of the mechanical flux derivative wrt to a dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_mechanical_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the derivative of the divergence of the mechanical flux wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] mDivMechanicalFluxDofDer derivative of the divergence of the mechanical flux wrt dof
                 */
                const Matrix< DDRMat > & mechanical_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model divergence of the thermal flux
                 */
                void eval_thermal_divflux();

                /**
                 *  get the constitutive model divergence of the thermal flux
                 * @param[ out ] mThermalDivFlux constitutive model divergence of the thermal flux
                 */
                const Matrix< DDRMat > & thermal_divflux();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model divergence of the thermal flux derivative wrt to a dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_thermal_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the derivative of the divergence of the thermal flux wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] ThermalDivFluxDofDer derivative of the divergence of the thermal flux wrt dof
                 */
                const Matrix< DDRMat > & thermal_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model divergence of the work flux
                 */
                void eval_work_divflux();

                /**
                 *  get the constitutive model divergence of the work flux
                 * @param[ out ] mWorkDivFlux constitutive model divergence of the work flux
                 */
                const Matrix< DDRMat > & work_divflux();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model divergence of the work flux derivative wrt to a dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_work_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the derivative of the divergence of the work flux wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] mWorkDivFluxDofDer derivative of the divergence of the work flux wrt dof
                 */
                const Matrix< DDRMat > & work_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

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
                //------------------------------------------------------------------------------
                /**
                 * evaluate the traction
                 * @param[ in ] aNormal normal
                 */
                void eval_traction( const Matrix< DDRMat > & aNormal )
                {
                    MORIS_ERROR( false,
                            "CM_Compressible_Newtonian_Fluid::eval_traction - not implemented, use specific traction functions." );
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
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal )
                {
                    MORIS_ERROR( false,
                            "CM_Compressible_Newtonian_Fluid::eval_dTractiondDOF - not implemented, use specific traction functions." );
                };

                /**
                 * get the derivative of the traction wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ in ]  aNormal         normal
                 * @param[ out ] mTractionDofDer derivative of the traction wrt dof
                 */
                const Matrix< DDRMat > & dTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
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
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal );

                /**
                 * get the derivative of the thermal flux traction wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ in ]  aNormal         normal
                 * @param[ out ] mTractionDofDer derivative of the thermal flux traction wrt dof
                 */
                const Matrix< DDRMat > & thermal_dTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
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
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal );

                /**
                 * get the derivative of the advective energy flux traction wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ in ]  aNormal         normal
                 * @param[ out ] mTractionDofDer derivative of the traction wrt dof
                 */
                const Matrix< DDRMat > & energy_dTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
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
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal );

                /**
                 * get the derivative of the mechanical work flux traction wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ in ]  aNormal         normal
                 * @param[ out ] mTractionDofDer derivative of the mechanical work flux traction wrt dof
                 */
                const Matrix< DDRMat > & work_dTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
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
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal );

                /**
                 * get the derivative of the mechanical traction wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ in ]  aNormal         normal
                 * @param[ out ] mTractionDofDer derivative of the mechanical traction wrt dof
                 */
                const Matrix< DDRMat > & mechanical_dTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        const Matrix< DDRMat >             & aNormal );

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
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
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
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        const Matrix< DDRMat >             & aNormal,
                        const Matrix< DDRMat >             & aJump,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
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
                         const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

                 /**
                  * get the constitutive model thermal test traction
                  * @param[ in ]  aNormal       normal
                  * @param[ in ]  aTestDofTypes      Dof type of the test functions
                  * @param[ out ] mTestTraction constitutive model test traction
                  */
                 const Matrix< DDRMat > & thermal_testTraction(
                         const Matrix< DDRMat >             & aNormal,
                         const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

                 //------------------------------------------------------------------------------
                 /**
                  * evaluate the constitutive model thermal test traction derivative wrt to a dof type
                  * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
                  * @param[ in ] aNormal         normal
                  * @param[ in ] aJump           jump in field values on boundary
                  * @param[ in ] aTestDofTypes   Dof type of the test functions
                  */
                 void eval_thermal_dTestTractiondDOF(
                         const moris::Cell< MSI::Dof_Type > & aDofTypes,
                         const Matrix< DDRMat >             & aNormal,
                         const Matrix< DDRMat >             & aJump,
                         const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

                 /**
                  * get the derivative of the thermal test traction wrt dof
                  * @param[ in ]  aDofType           group of dof type
                  * @param[ in ]  aNormal            normal
                  * @param[ in ]  aJump              jump in field values on boundary
                  * @param[ in ]  aTestDofTypes      Dof type of the test functions
                  * @param[ out ] mdTestTractiondDof derivative of the traction wrt dof
                  */
                 const Matrix< DDRMat > & thermal_dTestTractiondDOF(
                         const moris::Cell< MSI::Dof_Type > & aDofType,
                         const Matrix< DDRMat >             & aNormal,
                         const Matrix< DDRMat >             & aJump,
                         const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

                 //------------------------------------------------------------------------------
                 //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model mechanical test traction
                 * @param[ in ]  aNormal normal
                 * @param[ in ]  aTestDofTypes      Dof type of the test functions
                 */
                void eval_mechanical_testTraction(
                        const Matrix< DDRMat >             & aNormal,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

                /**
                 * get the constitutive model mechanical test traction
                 * @param[ in ]  aNormal       normal
                 * @param[ in ]  aTestDofTypes      Dof type of the test functions
                 * @param[ out ] mTestTraction constitutive model test traction
                 */
                const Matrix< DDRMat > & mechanical_testTraction(
                        const Matrix< DDRMat >             & aNormal,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model mechanical test traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal         normal
                 * @param[ in ] aJump           jump in field values on boundary
                 * @param[ in ] aTestDofTypes   Dof type of the test functions
                 */
                void eval_mechanical_dTestTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal,
                        const Matrix< DDRMat >             & aJump,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

                /**
                 * get the derivative of the mechanical test traction wrt dof
                 * @param[ in ]  aDofType           group of dof type
                 * @param[ in ]  aNormal            normal
                 * @param[ in ]  aJump              jump in field values on boundary
                 * @param[ in ]  aTestDofTypes      Dof type of the test functions
                 * @param[ out ] mdTestTractiondDof derivative of the traction wrt dof
                 */
                const Matrix< DDRMat > & mechanical_dTestTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        const Matrix< DDRMat >             & aNormal,
                        const Matrix< DDRMat >             & aJump,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

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

                //------------------------------------------------------------------------------
                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model divergence of the strains
                 */
                void eval_divstrain()
                {
                    MORIS_ERROR( false,
                            "CM_Compressible_Newtonian_Fluid::eval_divstrain - not implemented, use specific flux functions." );
                };

                /**
                 * get the constitutive model divergence of the strains
                 * @param[ in ]  aCMFunctionType  enum indicating which flux is called, if there are several
                 * @param[ out ] mDivStrain constitutive model divergence of the strains
                 */
                const Matrix< DDRMat > & divstrain( enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                /**
                 * evaluate the constitutive model divergence of the strains derivatives wrt to a dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_ddivstraindu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false,
                            "CM_Compressible_Newtonian_Fluid::eval_ddivstraindu - not implemented, use specific flux functions." );
                };

                /**
                 * get the derivative of the divergence of the strains wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] mDivStrainDofDer derivative of the strains wrt dof
                 */
                const Matrix< DDRMat > & ddivstraindu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //--------------------------------------------------------------------------------------------------------------
                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the divergence of the thermal strain
                 */
                void eval_thermal_divstrain();

                /**
                 * get the divergence of the divergence of the thermal strain
                 * @param[ out ] mThermalDivStrain constitutive model divergence of the thermal strain
                 */
                const Matrix< DDRMat > & thermal_divstrain();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the divergence of the thermal strain wrt dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_thermal_ddivstraindu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the divergence of the thermal strain
                 * @param[ in ]  aDofTypes  DoF type wrt with the derivative are computed
                 * @param[ out ] mThermalDivStrainDof constitutive model div(grad(T)) wrt to the dofs
                 */
                const Matrix< DDRMat > & thermal_ddivstraindu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the divergence of the strain
                 */
                void eval_divstrainrate()
                {
                    ( this->*m_eval_divstrainrate )();
                };
                void eval_divstrainrate_2d();
                void eval_divstrainrate_3d();

                /**
                 * get the  divergence of the strain rate
                 * @param[ out ] mDivStrainRate constitutive model divergence of the strain rate
                 */
                const Matrix< DDRMat > & divstrainrate();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the divergence of the strain wrt dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_ddivstrainratedu( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    ( this->*m_eval_ddivstrainratedu)( aDofTypes );
                };
                void eval_ddivstrainratedu_2d( const moris::Cell< MSI::Dof_Type > & aDofTypes );
                void eval_ddivstrainratedu_3d( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the divergence of the strain rate
                 * @param[ in ]  aDofTypes  DoF type wrt with the derivative are computed
                 * @param[ out ] mDivStrainRateDof constitutive model divergence of the strain rate wrt to the dofs
                 */
                const Matrix< DDRMat > & ddivstrainratedu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the divergence of the divergence of the velocity div(div(u)*I)
                 */
                void eval_divDivVel()
                {
                    ( this->*m_eval_divDivVel )();
                }
                void eval_divDivVel_2d();
                void eval_divDivVel_3d();

                /**
                 * get the divergence of the divergence of the velocity div(div(u)*I)
                 * @param[ out ] mDivDivVel constitutive model div(div(u)*I)
                 */
                const Matrix< DDRMat > & divDivVel();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the divergence of the divergnece of the velocity wrt dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_dDivDivVeldu( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    ( this->*m_eval_dDivDivVeldu)( aDofTypes );
                }
                void eval_dDivDivVeldu_2d( const moris::Cell< MSI::Dof_Type > & aDofTypes );
                void eval_dDivDivVeldu_3d( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the divergence of the divergence of the velocity
                 * @param[ in ]  aDofTypes  DoF type wrt with the derivative are computed
                 * @param[ out ] mdDivDivVeldu constitutive model div(div(u)*I) wrt to the dofs
                 */
                const Matrix< DDRMat > & dDivDivVeldu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
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

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluates velocity matrix for operation on flattened stress and strain vectors
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
                /**
                 * get the multiplication matrix for condensed tensors
                 * @param[ out ] mMultipMat multiplication matrix for condensed tensors
                 */
                const Matrix< DDRMat > & MultipMat();

                //------------------------------------------------------------------------------
        };

        //--------------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_COMPRESSIBLE_NEWTONIAN_FLUID_HPP_ */


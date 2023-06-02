/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_TURBULENCE_HPP_
#define SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_TURBULENCE_HPP_

#include <iostream>
#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp"
#include "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp"//FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class CM_Diffusion_Linear_Isotropic_Turbulence : public CM_Diffusion_Linear_Isotropic
        {
                //------------------------------------------------------------------------------
            protected:
                // default local properties
                std::shared_ptr< Property > mPropKinViscosity = nullptr;
                std::shared_ptr< Property > mPropPrandtlT     = nullptr;

            private:
                // Default dof type for CM
                MSI::Dof_Type mTempDof       = MSI::Dof_Type::TEMP;
                MSI::Dof_Type mDofViscosity  = MSI::Dof_Type::VISCOSITY;
                MSI::Dof_Type mThetaDof      = MSI::Dof_Type::UNDEFINED;

                // property type for CM
                enum class CM_Property_Type
                {
                    CONDUCTIVITY,
                    HEAT_CAPACITY,
                    DENSITY,
                    EIGEN_STRAIN,
                    KIN_VISCOSITY,
                    TURBULENT_PRANDTL,
                    MAX_ENUM
                };

                // flags for chi related evaluation
                bool mChiEval = true;
                moris::Matrix< DDBMat > mdChiduEval;
                moris::Matrix< DDBMat > mdChidxEval;
                moris::Matrix< DDBMat > mdChidxduEval;

                // storage for chi
                moris::real mChi;
                moris::Cell< Matrix< DDRMat > > mdChidu;
                moris::Cell< Matrix< DDRMat > > mdChidx;
                moris::Cell< Cell< Matrix< DDRMat > > > mdChidxdu;

                // flags for fv1 related evaluation
                bool mFv1Eval = true;
                moris::Matrix< DDBMat > mdFv1duEval;
                moris::Matrix< DDBMat > mdFv1dxEval;
                moris::Matrix< DDBMat > mdFv1dxduEval;

                // storage for fv1
                moris::real mFv1;
                moris::Cell< Matrix< DDRMat > > mdFv1du;
                moris::Cell< Matrix< DDRMat > > mdFv1dx;
                moris::Cell< Cell< Matrix< DDRMat > > > mdFv1dxdu;

                // flags for turbulence dynamic viscosity related evaluation
                bool mTurbDynViscEval = true;
                moris::Matrix< DDBMat > mdTurbDynViscduEval;
                moris::Matrix< DDBMat > mdTurbDynViscdxEval;
                moris::Matrix< DDBMat > mdTurbDynViscdxduEval;

                // storage for turbulence dynamic viscosity related evaluation
                Matrix< DDRMat > mTurbDynVisc;
                moris::Cell< Matrix< DDRMat > > mdTurbDynViscdu;
                moris::Cell< Matrix< DDRMat > > mdTurbDynViscdx;
                moris::Cell< Cell< Matrix< DDRMat > > > mdTurbDynViscdxdu;

                // flags for effective conductivity related evaluation
                bool mEffCondEval = true;
                moris::Matrix< DDBMat > mdEffCondduEval;
                moris::Matrix< DDBMat > mdEffConddxEval;
                moris::Matrix< DDBMat > mdEffConddxduEval;

                // storage for effective conductivity related evaluation
                Matrix< DDRMat > mEffCond;
                moris::Cell< Matrix< DDRMat > > mdEffConddu;
                moris::Cell< Matrix< DDRMat > > mdEffConddx;
                moris::Cell< Cell< Matrix< DDRMat > > > mdEffConddxdu;

                //------------------------------------------------------------------------------
            public:
                /*
                 * trivial constructor
                 */
                CM_Diffusion_Linear_Isotropic_Turbulence();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~CM_Diffusion_Linear_Isotropic_Turbulence(){};

                //------------------------------------------------------------------------------
                /*
                 * @return constitutive_type
                 */
                Constitutive_Type get_constitutive_type() const
                {
                    return Constitutive_Type::DIFF_LIN_ISO_TURBULENCE;
                }

                //------------------------------------------------------------------------------
                /*
                 * reset evaluation flag
                 * (child implementation)
                 */
                void reset_eval_flags();

                //------------------------------------------------------------------------------
                /*
                 * create a global dof type list including constitutive and property dependencies
                 * (child implementation)
                 */
                void build_global_dof_type_list();

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
                 * set local properties
                 */
                void set_local_properties();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model flux
                 */
                void eval_flux();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the flux derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dFluxdDOF ( mSpaceDim x numDerDof )
                 */
                void eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the divergence of the flux
                 */
                void eval_divflux();

                //------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the divergence of the flux wrt dof type
                 */
                void eval_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model traction
                 * @param[ in ] aNormal normal
                 * traction ( 1 x 1 )
                 */
                void eval_traction( const Matrix< DDRMat > & aNormal );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 * dTractiondDOF ( 1 x numDerDof )
                 */
                void eval_dTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model test traction
                 * @param[ in ] aNormal normal
                 * test traction ( numDof x 1 )
                 */
                void eval_testTraction(
                        const Matrix< DDRMat >             & aNormal,
                        const moris::Cell< MSI::Dof_Type > & aTestDofType );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model test traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 * dTestTractiondDOF ( numDof x numDerDof )
                 */
                void eval_dTestTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

                /**
                 * evaluate the constitutive model test traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 * dTestTractiondDOF ( numDof x numDerDof )
                 */
                void eval_dTestTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal,
                        const Matrix< DDRMat >             & aJump,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model matrix
                 * constitutive matrix ( mSpaceDim x mSpaceDim )
                 */
                void eval_const();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the divergence of the strain wrt dof type
                 */
                void eval_dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * get the the effective conductivity k_eff = k + k_t
                 * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
                 *               if there are several
                 * @param[ out ] mEffCond effective conductivity
                 */
                const Matrix< DDRMat > & effective_conductivity(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                /**
                 * get the the turbulent dynamic viscosity mu_t = rho * vtilde * tf1
                 * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
                 *               if there are several
                 * @param[ out ] mTurbDynVisc effective conductivity
                 */
                const Matrix< DDRMat > & turbulent_dynamic_viscosity(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
            private:

                //------------------------------------------------------------------------------
                /**
                 * evaluate the effective conductivity k_eff = k + k_t
                 */
                void eval_effective_conductivity();

//                /**
//                 * get the the effective conductivity k_eff = k + k_t
//                 * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
//                 *               if there are several
//                 * @param[ out ] mEffCond effective conductivity
//                 */
//                const Matrix< DDRMat > & effective_conductivity(
//                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the effective conductivity wrt space
                 * @param[ in ] aOrder order of the derivative
                 */
                void eval_deffconddx( uint aOrder );

                /**
                 * get the derivative of the effective conductivity wrt space
                 * @param[ in ] aOrder order of the derivative
                 * @param[ in ] aCMFunctionType enum for specific type of which effective conductivity
                 * @param[ out ] mdeffconddx derivative of the effective conductivity wrt space
                 */
                const Matrix< DDRMat > & deffconddx(
                        uint                  aOrder,
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the effective conductivity derivative wrt to dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 */
                void eval_deffconddu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the derivative of the effective conductivity wrt dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ in ] aCMFunctionType enum for specific type of which effective conductivity
                 * @param[ out ] mdeffconddu derivative of the effective conductivity wrt dof types
                 */
                const Matrix< DDRMat > & deffconddu(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the effective conductivity derivative wrt to space and dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aOrder order of the space derivative
                 */
                void eval_deffconddxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                        uint                                 aOrder);

                /**
                 * get the derivative of the effective conductivity wrt space and dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ in ] aOrder order of the space derivative
                 * @param[ in ] aCMFunctionType enum for specific type of which effective conductivity
                 * @param[ out ] mdeffconddxdu derivative of the effective conductivity wrt dof types
                 */
                const Matrix< DDRMat > & deffconddxdu(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        uint                                 aOrder,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the turbulent dynamic viscosity mu_t = rho * vtilde * tf1
                 */
                void eval_turbulent_dynamic_viscosity();

//                /**
//                 * get the the turbulent dynamic viscosity mu_t = rho * vtilde * tf1
//                 * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
//                 *               if there are several
//                 * @param[ out ] mTurbDynVisc effective conductivity
//                 */
//                const Matrix< DDRMat > & turbulent_dynamic_viscosity(
//                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the turbulent dynamic viscosity derivative wrt to dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 */
                void eval_dturbdynviscdu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                /**
                 * get the derivative of the turbulent dynamic viscosity wrt dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ in ] aCMFunctionType enum for specific type of which turbulent dynamic viscosity
                 * @param[ out ] mdturbdynviscdu derivative of the turbulent dynamic viscosity wrt dof types
                 */
                const Matrix< DDRMat > & dturbdynviscdu(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the turbulent dynamic viscosity wrt space
                 * @param[ in ] aOrder order of the derivative
                 */
                void eval_dturbdynviscdx( uint aOrder );

                /**
                 * get the derivative of the turbulent dynamic viscosity wrt space
                 * @param[ in ] aOrder order of the derivative
                 * @param[ in ] aCMFunctionType enum for specific type of which turbulent dynamic viscosity
                 * @param[ out ] mdeffconddx derivative of the turbulent dynamic viscosity wrt space
                 */
                const Matrix< DDRMat > & dturbdynviscdx(
                        uint                  aOrder,
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the turbulent dynamic viscosity derivative wrt to space and dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aOrder order of the space derivative
                 */
                void eval_dturbdynviscdxdu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        uint                                 aOrder);

                /**
                 * get the derivative of the turbulent dynamic viscosity wrt space and dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ in ] aOrder order of the space derivative
                 * @param[ in ] aCMFunctionType enum for specific type of which turbulent dynamic viscosity
                 * @param[ out ] mdeffconddxdu derivative of the turbulent dynamic viscosity wrt dof types
                 */
                const Matrix< DDRMat > & dturbdynviscdxdu(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        uint                                 aOrder,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * get chi = nu_tilde/nu
                 * @param[ in ]  aCMFunctionType enum indicating which function if several
                 * @param[ out ] mChi chi
                 */
                real chi(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                /**
                 * get the derivative of chi wrt dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ in ] aCMFunctionType enum for specific type of which chi
                 * @param[ out ] mdchidu derivative of the chi wrt dof types
                 */
                const Matrix< DDRMat > & dchidu(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                /**
                 * get the derivative of chi wrt space
                 * @param[ in ] aOrder order of the derivative
                 * @param[ in ] aCMFunctionType enum for specific type of which chi
                 * @param[ out ] mdchidx derivative of the chi wrt space
                 */
                const Matrix< DDRMat > & dchidx(
                        uint                  aOrder,
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                /**
                 * get the derivative of chi wrt space and dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ in ] aOrder order of the space derivative
                 * @param[ in ] aCMFunctionType enum for specific type of which chi
                 * @param[ out ] mdchidxdu derivative of the chi wrt dof types
                 */
                const Matrix< DDRMat > & dchidxdu(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        uint                                 aOrder,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * get fv1 = chi^3 / ( chi^3 + cv1^3 )
                 * @param[ in ]  aCMFunctionType enum for specific type of fv2
                 * @param[ out ] mFv1 fv1
                 */
                real fv1(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                /**
                 * get the derivative of fv1 wrt dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ in ] aCMFunctionType enum for specific type of fv1
                 * @param[ out ] mdfv1du derivative of fv1 wrt dof types
                 */
                const Matrix< DDRMat > & dfv1du(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                /**
                 * get the derivative of fv1 wrt space
                 * @param[ in ] aOrder order of the derivative
                 * @param[ in ] aCMFunctionType enum for specific type of which fv1
                 * @param[ out ] mdfv1dx derivative of the fv1 wrt space
                 */
                const Matrix< DDRMat > & dfv1dx(
                        uint                  aOrder,
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                /**
                 * get the derivative of fv1 wrt space and dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ in ] aOrder order of the space derivative
                 * @param[ in ] aCMFunctionType enum for specific type of which fv1
                 * @param[ out ] mdfv1dxdu derivative of the fv1 wrt dof types
                 */
                const Matrix< DDRMat > & dfv1dxdu(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        uint                                 aOrder,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_TURBULENCE_HPP_ */


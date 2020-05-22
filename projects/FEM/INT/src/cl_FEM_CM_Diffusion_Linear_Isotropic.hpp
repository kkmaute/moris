/*
 * cl_FEM_CM_Diffusion_Linear_Isotropic.hpp
 *
 *  Created on: Sep 17, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_HPP_
#define SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_HPP_

#include <iostream>
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
//------------------------------------------------------------------------------

        class CM_Diffusion_Linear_Isotropic : public Constitutive_Model
        {

//------------------------------------------------------------------------------
        private:

            // property type for CM
            enum class Property_Type
            {
                CONDUCTIVITY,
                HEAT_CAPACITY,
                DENSITY,
                MAX_ENUM
            };

            // Local string to property enum map
            std::map< std::string, CM_Diffusion_Linear_Isotropic::Property_Type > mPropertyMap;

//------------------------------------------------------------------------------
        public:
            /*
             * trivial constructor
             */
            CM_Diffusion_Linear_Isotropic();

//------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~CM_Diffusion_Linear_Isotropic(){};

//------------------------------------------------------------------------------
            /**
             * set a property pointer
             * @param[ in ] aProperty     a property pointer
             * @param[ in ] aPropertyType a char
             */
            void set_property( std::shared_ptr< fem::Property > aProperty,
                               std::string                      aPropertyString )
            {
                // check that aPropertyString makes sense
                MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                             "CM_Diffusion_Linear_Isotropic::set_property - Unknown aPropertyString." );

                // set the property in the property cell
                mProperties( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
            };

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux
             */
            void eval_flux();

//------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model change rate of enthalpy
             */
            virtual void eval_Hdot();

//------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model change rate of spatial gradient of enthalpy (needed for GGLS-stabilization)
             */
            virtual void eval_gradHdot();

//------------------------------------------------------------------------------
            /**
             * evaluates the gradient of the divergence of the flux (needed for GGLS-stabilization)
             */
            void eval_graddivflux();
//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test flux
             * flux ( mSpaceDim x 1)
             */
            void eval_testFlux();

//--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the divergence of the flux
             */
            void eval_divflux();

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model traction
             * @param[ in ] aNormal normal
             * traction ( 1 x 1 )
             */
            void eval_traction( const Matrix< DDRMat > & aNormal );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test traction
             * @param[ in ] aNormal normal
             * test traction ( numDof x 1 )
             */
            void eval_testTraction( const Matrix< DDRMat >             & aNormal,
                                    const moris::Cell< MSI::Dof_Type > & aTestDofType );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain
             * strain ( mSpaceDim x 1 )
             */
            void eval_strain();

//--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the divergence of the strain
             */
            void eval_divstrain();

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test strain
             * test strain ( mSpaceDim x numDof  )
             */
            void eval_testStrain();

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model matrix
             * constitutive matrix ( mSpaceDim x mSpaceDim )
             */
            void eval_const();

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dFluxdDOF ( mSpaceDim x numDerDof )
             */
            void eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model enthalpy change rate wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dHdotdDOF ( 1 x numDerDof )
             */
            virtual void eval_dHdotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model gradient of enthalpy change rate wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dgradHdotdDOF ( mSpaceDim x numDerDof )
             */
            virtual void eval_dGradHdotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model gradient of divergence of flux wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dGradDivFluxdDOF ( mSpaceDim x numDerDof )
             */
            void eval_dGradDivFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the divergence of the flux wrt dof type
             */
            void eval_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model traction derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * @param[ in ] aNormal   normal
             * dTractiondDOF ( 1 x numDerDof )
             */
            void eval_dTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                     const Matrix< DDRMat >             & aNormal );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test traction derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * @param[ in ] aNormal   normal
             * dTestTractiondDOF ( numDof x numDerDof )
             */
            void eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                         const Matrix< DDRMat >             & aNormal,
                                         const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

            /**
             * evaluate the constitutive model test traction derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * @param[ in ] aNormal   normal
             * dTestTractiondDOF ( numDof x numDerDof )
             */
            void eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                         const Matrix< DDRMat >             & aNormal,
                                         const Matrix< DDRMat >             & aJump,
                                         const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dStraindDOF ( mSpaceDim x numDerDof )
             */
            void eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the divergence of the strain wrt dof type
             */
            void eval_ddivstraindu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the divergence of the strain wrt dof type
             */
            void eval_dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux derivative wrt to a dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             * dFluxdDV ( mSpaceDim x numDerDv )
             */
            void eval_dFluxdDV( const moris::Cell< PDV_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain derivative wrt to a dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             * dStraindDV ( mSpaceDim x numDerDV )
             */
            void eval_dStraindDV( const moris::Cell< PDV_Type > & aDofTypes );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_DIFFUSION_LINEAR_ISOTROPIC_HPP_ */

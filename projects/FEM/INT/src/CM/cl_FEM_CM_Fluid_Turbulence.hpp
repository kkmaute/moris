#ifndef SRC_FEM_CL_FEM_CM_FLUID_TURBULENCE_HPP_
#define SRC_FEM_CL_FEM_CM_FLUID_TURBULENCE_HPP_

#include <map>
//MRS/COR/src
#include "typedefs.hpp"
#include "cl_Cell.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_CM_Fluid_Incompressible.hpp"
#include "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp"

namespace moris
{
    namespace fem
    {

        //--------------------------------------------------------------------------------------------------------------
        class CM_Fluid_Turbulence : public CM_Fluid_Incompressible
        {

                //--------------------------------------------------------------------------------------------------------------
            private:

                // default dof type
                MSI::Dof_Type mDofVelocity  = MSI::Dof_Type::VX;
                MSI::Dof_Type mDofPressure  = MSI::Dof_Type::P;
                MSI::Dof_Type mDofViscosity = MSI::Dof_Type::VISCOSITY;

                // property type for CM
                enum class CM_Property_Type
                {
                    DENSITY,   // fluid density
                    VISCOSITY, // fluid viscosity
                    MAX_ENUM
                };

                // Spalart Allmaras model constants
                real mCv1 = 7.1;

                //--------------------------------------------------------------------------------------------------------------
            public:

                //--------------------------------------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                CM_Fluid_Turbulence();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~CM_Fluid_Turbulence(){};

                //------------------------------------------------------------------------------
                /*
                 * @return constitutive_type
                 */
                Constitutive_Type
                get_constitutive_type() const
                {
                    return Constitutive_Type::FLUID_TURBULENCE;
                }

                //------------------------------------------------------------------------------
                /**
                 * set constitutive model dof types
                 * @param[ in ] aDofTypes   a list of group of dof types
                 * @param[ in ] aDofStrings a list of strings to describe the dof types
                 */
                void set_dof_type_list(
                        moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                        moris::Cell< std::string >                  aDofStrings );

                //------------------------------------------------------------------------------
                /**
                 * set constitutive model dv types
                 * @param[ in ] aDvTypes   a list of group of dv types
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
                void eval_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the flux wrt space
                 * @param[ in ] aOrder order of the derivative
                 */
                void eval_dfluxdx( uint aOrder );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the test flux
                 */
                void eval_testFlux();

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
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model flux derivative wrt to a dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                void eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 */
                void eval_dTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model test traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 */
                void eval_dTestTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal,
                        const Matrix< DDRMat >             & aJump,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

                //--------------------------------------------------------------------------------------------------------------
            private:

                //------------------------------------------------------------------------------
                /**
                 * compute turbulence viscosity = rho * vtilde * tf1
                 * @param[ out ] turbulence viscosity
                 */
                real compute_turbulence_viscosity();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of turbulence viscosity
                 *  wrt to a dof type
                 * @param[ in ] aDofTypes      a list of dof type wrt which
                 *                             the derivative is requested
                 * @param[ in ] adviscositytdu a matrix to fill with dviscositytdu
                 */
                void compute_dviscositytdu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adviscositytdu );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of viscosityt wrt to x
                 * @param[ in ] adviscositytdx a matrix to fill with dviscositytdx
                 */
                void compute_dviscositytdx( Matrix< DDRMat > & adviscositytdx );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the gradx of turbulence viscosity
                 *  wrt to a dof type
                 * @param[ in ] aDofTypes        a list of dof type wrt which
                 *                               the derivative is requested
                 * @param[ in ] adviscositytdxdu a matrix to fill with dviscositytdxdu
                 */
                void compute_dviscositytdxdu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adviscositytdxdu );

//                //------------------------------------------------------------------------------
//                /**
//                 * compute chi = viscosityDof / viscosityProp
//                 * @param[ out ] chi
//                 */
//                real compute_chi();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of chi wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adchidu    a matrix to fill with dchidu
//                 */
//                void compute_dchidu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adchidu );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of chi wrt to x
//                 * @param[ in ] adchidx a matrix to fill with dchidx
//                 */
//                void compute_dchidx( Matrix< DDRMat > & adchidx );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of dchidx wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adchidxdu  a matrix to fill with dchidxdu
//                 *                  */
//                void compute_dchidxdu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adchidxdu );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute fv1 = chi³ / ( chi³ + cv1³)
//                 * @param[ out ] fv1
//                 */
//                real compute_fv1();
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of fv1 wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adfv1du    a matrix to fill with dfv1du
//                 */
//                void compute_dfv1du(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adfv1du );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of fv1 wrt to x
//                 * @param[ in ] adfv1dx a matrix to fill with dfv1dx
//                 */
//                void compute_dfv1dx( Matrix< DDRMat > & adfv1dx );
//
//                //------------------------------------------------------------------------------
//                /**
//                 * compute the derivative of dfv1dx wrt to a dof type
//                 * @param[ in ] aDofTypes  a list of dof type wrt which
//                 *                         the derivative is requested
//                 * @param[ in ] adfv1dxdu  a matrix to fill with dfv1dxdu
//                 */
//                void compute_dfv1dxdu(
//                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                        Matrix< DDRMat >                   & adfv1dxdu );

                //--------------------------------------------------------------------------------------------------------------
        };

        //--------------------------------------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_FLUID_TURBULENCE_HPP_ */

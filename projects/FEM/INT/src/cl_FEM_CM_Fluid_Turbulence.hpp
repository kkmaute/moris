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

namespace moris
{
    namespace fem
    {

        //--------------------------------------------------------------------------------------------------------------
        class CM_Fluid_Turbulence : public Constitutive_Model
        {

                //--------------------------------------------------------------------------------------------------------------
            private:

                // default dof type
                MSI::Dof_Type mDofVelocity = MSI::Dof_Type::VX;
                MSI::Dof_Type mDofViscosity = MSI::Dof_Type::VISCOSITY;

                // property type for CM
                enum class CM_Property_Type
                {
                    VISCOSITY, // fluid viscosity
                    DENSITY,   // fluid density
                    MAX_ENUM
                };

                // local string to property enum map
                std::map< std::string, CM_Property_Type > mPropertyMap;

                // function pointers
                void ( CM_Fluid_Turbulence:: * m_eval_strain )() = nullptr;
                void ( CM_Fluid_Turbulence:: * m_eval_divstrain )() = nullptr;
                void ( CM_Fluid_Turbulence:: * m_eval_teststrain )() = nullptr;
                void ( CM_Fluid_Turbulence:: * m_eval_dstraindx )( uint aOrder ) = nullptr;
                void ( CM_Fluid_Turbulence:: * m_eval_ddivstraindu )(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes ) = nullptr;
                void ( CM_Fluid_Turbulence:: * m_flatten_normal )(
                        const Matrix< DDRMat > & aNormal,
                        Matrix< DDRMat > & aFlatNormal ) = nullptr;

                // FIXME temp all the constants
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

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * set a property pointer
                 * @param[ in ] aProperty     a property pointer
                 * @param[ in ] aPropertyType a string defining the property
                 */
                void set_property(
                        std::shared_ptr< fem::Property > aProperty,
                        std::string                      aPropertyString );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * get a property pointer
                 * @param[ in ]  aPropertyType a string defining the property
                 * @param[ out ] aProperty     a property pointer
                 */
                std::shared_ptr< Property > get_property( std::string aPropertyString );

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
                void eval_ddivstraindu( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    ( this->*m_eval_ddivstraindu)( aDofTypes );
                }
                void eval_ddivstraindu_2d( const moris::Cell< MSI::Dof_Type > & aDofTypes );
                void eval_ddivstraindu_3d( const moris::Cell< MSI::Dof_Type > & aDofTypes );

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
                /**
                 * evaluate the constitutive model strain derivative wrt to a dof type
                 * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
                 * @param[ in ] adStraindDOF a matrix to fill with derivative evaluation
                 */
                void eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * flatten normal vector
                 * @param[ in ] aNormal          a normal vector
                 * @param[ in ] aFlattenedNormal a matrix for flattened normal to fill
                 */
                void flatten_normal(
                        const Matrix< DDRMat > & aNormal,
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

                //------------------------------------------------------------------------------
                /**
                 * compute chi = viscosityDof / viscosityPtop
                 * @param[ out ] chi
                 */
                real compute_chi();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of chi wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adchidu    a matrix to fill with dchidu
                 */
                void compute_dchidu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adchidu );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of chi wrt to x
                 * @param[ in ] adchidx a matrix to fill with dchidx
                 */
                void compute_dchidx( Matrix< DDRMat > & adchidx );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of dchidx wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adchidxdu  a matrix to fill with dchidxdu
                 *                  */
                void compute_dchidxdu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adchidxdu );

                //------------------------------------------------------------------------------
                /**
                 * compute fv1 = chi³ / ( chi³ + cv1³)
                 * @param[ out ] fv1
                 */
                real compute_fv1();

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of fv1 wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adfv1du    a matrix to fill with dfv1du
                 */
                void compute_dfv1du(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adfv1du );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of fv1 wrt to x
                 * @param[ in ] adfv1dx a matrix to fill with dfv1dx
                 */
                void compute_dfv1dx( Matrix< DDRMat > & adfv1dx );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of dfv1dx wrt to a dof type
                 * @param[ in ] aDofTypes  a list of dof type wrt which
                 *                         the derivative is requested
                 * @param[ in ] adfv1dxdu  a matrix to fill with dfv1dxdu
                 */
                void compute_dfv1dxdu(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adfv1dxdu );

                //--------------------------------------------------------------------------------------------------------------
        };

        //--------------------------------------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_FLUID_TURBULENCE_HPP_ */

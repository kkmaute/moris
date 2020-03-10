#ifndef SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_HPP_
#define SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_HPP_

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

//--------------------------------------------------------------------------------------------------------------

        class CM_Fluid_Incompressible : public Constitutive_Model
        {

//--------------------------------------------------------------------------------------------------------------
        private:

            // property type for CM
            enum class CM_Property_Type
            {
                DENSITY,   // fluid density
                VISCOSITY, // fluid viscosity
                MAX_ENUM
            };

            // local string to property enum map
            std::map< std::string, CM_Property_Type > mPropertyMap;

            // function pointers
            void ( CM_Fluid_Incompressible:: * m_eval_strain )() = nullptr;
            void ( CM_Fluid_Incompressible:: * m_eval_divstrain )() = nullptr;
            void ( CM_Fluid_Incompressible:: * m_eval_teststrain )() = nullptr;
            void ( CM_Fluid_Incompressible:: * m_eval_dstraindx )( uint aOrder ) = nullptr;
            void ( CM_Fluid_Incompressible:: * m_eval_ddivstraindu )( const moris::Cell< MSI::Dof_Type > & aDofTypes ) = nullptr;
            void ( CM_Fluid_Incompressible:: * m_flatten_normal )( const Matrix< DDRMat > & aNormal,
                                                                         Matrix< DDRMat > & aFlatNormal ) = nullptr;

//--------------------------------------------------------------------------------------------------------------
        public:
            /*
             * constructor
             */
            CM_Fluid_Incompressible()
            {
                // set the property pointer cell size
                mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

                // populate the map
                mPropertyMap[ "Density" ]   = CM_Property_Type::DENSITY;
                mPropertyMap[ "Viscosity" ] = CM_Property_Type::VISCOSITY;

                // populate the dof map (default)
                mDofMap[ "Velocity" ] = MSI::Dof_Type::VX;
                mDofMap[ "Pressure" ] = MSI::Dof_Type::P;
            };

//--------------------------------------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~CM_Fluid_Incompressible(){};

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

//--------------------------------------------------------------------------------------------------------------
            /**
             * set a property pointer
             * @param[ in ] aProperty     a property pointer
             * @param[ in ] aPropertyType a string defining the property
             */
             void set_property( std::shared_ptr< fem::Property > aProperty,
                                std::string                      aPropertyString )
             {
                 // check that aPropertyString makes sense
                 MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                              "CM_Fluid_Incompressible::set_property - Unknown aPropertyString." );

                 // set the property in the property cell
                 mProperties( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
             };

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
             * evaluate the constitutive model test flux
             */
            void eval_testFlux();

//--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model traction
             * @param[ in ] aNormal normal
             */
            void eval_traction( const Matrix< DDRMat > & aNormal );

//--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test traction
             * @param[ in ] aNormal   normal
             */
            void eval_testTraction( const Matrix< DDRMat > & aNormal );

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
            void eval_dTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                     const Matrix< DDRMat >             & aNormal );

//--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test traction derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * @param[ in ] aNormal   normal
             */
            void eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                         const Matrix< DDRMat >             & aNormal,
                                         const Matrix< DDRMat >             & aJump );

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
            void flatten_normal( const Matrix< DDRMat > & aNormal,
                                       Matrix< DDRMat > & aFlatNormal )
            {
                ( this->*m_flatten_normal )( aNormal, aFlatNormal );
            }
            void flatten_normal_2d( const Matrix< DDRMat > & aNormal,
                                          Matrix< DDRMat > & aFlatNormal );
            void flatten_normal_3d( const Matrix< DDRMat > & aNormal,
                                          Matrix< DDRMat > & aFlatNormal );

//--------------------------------------------------------------------------------------------------------------
        private:

//--------------------------------------------------------------------------------------------------------------
        };

//--------------------------------------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_HPP_ */

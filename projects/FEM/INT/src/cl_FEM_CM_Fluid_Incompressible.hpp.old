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
            enum class Property_Type
            {
                RHO, // fluid density
                MU,  // fluid viscosity
                MAX_ENUM
            };

            // local string to property enum map
            std::map< std::string, CM_Fluid_Incompressible::Property_Type > mPropertyMap;

//--------------------------------------------------------------------------------------------------------------
        public:
            /*
             * constructor
             */
            CM_Fluid_Incompressible()
            {
                // set the property pointer cell size
                mProperties.resize( static_cast< uint >( CM_Fluid_Incompressible::Property_Type::MAX_ENUM ), nullptr );

                // populate the map
                mPropertyMap[ "Density" ]   = CM_Fluid_Incompressible::Property_Type::RHO;
                mPropertyMap[ "Viscosity" ] = CM_Fluid_Incompressible::Property_Type::MU;

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
             * evaluate the constitutive model strain
             */
            void eval_strain();

//--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test strain
             */
            void eval_testStrain();

//--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model matrix
             */
            void eval_const();

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
             * evaluate the constitutive model matrix derivative wrt to a dof type
             * @param[ in ] aDofTypes   a dof type wrt which the derivative is evaluated
             * @param[ in ] adConstdDOF a matrix to fill with derivative evaluation
             */
            void eval_dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

//--------------------------------------------------------------------------------------------------------------
            /**
             * flatten normal vector
             * @param[ in ] aNormal          a normal vector
             * @param[ in ] aFlattenedNormal a matrix for flattened normal to fill
             */
            void flatten_normal( const Matrix< DDRMat > & aNormal,
                                       Matrix< DDRMat > & aFlattenedNormal );

//--------------------------------------------------------------------------------------------------------------
            /**
             * set spatial dimensions. Modified from base to set function pointers to the appropriate eval_const()
             * @param[ in ] aSpaceDim a spatial dimension
             */
            void set_space_dim( uint aSpaceDim );

//--------------------------------------------------------------------------------------------------------------
        private:

//--------------------------------------------------------------------------------------------------------------
        };
//--------------------------------------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_FLUID_INCOMPRESSIBLE_HPP_ */

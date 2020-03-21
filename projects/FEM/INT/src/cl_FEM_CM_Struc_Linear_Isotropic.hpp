#ifndef SRC_FEM_CL_FEM_CM_STRUC_LINEAR_ISOTROPIC_HPP_
#define SRC_FEM_CL_FEM_CM_STRUC_LINEAR_ISOTROPIC_HPP_

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

        class CM_Struc_Linear_Isotropic : public Constitutive_Model
        {

//--------------------------------------------------------------------------------------------------------------
        private:

            // property type for CM
            enum class Property_Type
            {
                EMOD,
                NU,
                CTE,
                TEMP_REF,
                MAX_ENUM
            };

            // Local string to property enum map
            std::map< std::string, CM_Struc_Linear_Isotropic::Property_Type > mPropertyMap;

            Model_Type mPlaneType  = Model_Type::PLANE_STRESS; // Plane stress or plane strain, only used in 2d
            Model_Type mTensorType = Model_Type::FULL; // Hydrostatic or deviatoric (default: full tensor)
            void ( moris::fem::CM_Struc_Linear_Isotropic:: *mConstFunc)(moris::real, moris::real) = &CM_Struc_Linear_Isotropic::full_3d;

//--------------------------------------------------------------------------------------------------------------
        public:
            /*
             * trivial constructor
             */
            CM_Struc_Linear_Isotropic()
            {
                // set the property pointer cell size
                mProperties.resize( static_cast< uint >( CM_Struc_Linear_Isotropic::Property_Type::MAX_ENUM ), nullptr );

                // populate the map
                mPropertyMap[ "YoungsModulus" ]         = CM_Struc_Linear_Isotropic::Property_Type::EMOD;
                mPropertyMap[ "PoissonRatio" ]          = CM_Struc_Linear_Isotropic::Property_Type::NU;
                mPropertyMap[ "CTE" ]                   = CM_Struc_Linear_Isotropic::Property_Type::CTE;
                mPropertyMap[ "ReferenceTemperature" ]  = CM_Struc_Linear_Isotropic::Property_Type::TEMP_REF;

                // populate the dof map
                mDofMap[ "Displacement" ] = MSI::Dof_Type::UX;
                mDofMap[ "Temperature" ]  = MSI::Dof_Type::UNDEFINED;
                mDofMap[ "Pressure" ]     = MSI::Dof_Type::UNDEFINED;
            };

//--------------------------------------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~CM_Struc_Linear_Isotropic(){};

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
                              "CM_Struc_Linear_Isotropic::set_property - Unknown aPropertyString." );

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
             * evaluate the inverse of the bulk modulus, K
             * @return 1/K
             */
            moris::real eval_inv_bulk_modulus();

//--------------------------------------------------------------------------------------------------------------
            /**
             * evaluate the B matrix from the field interpolator and make it flattened to a 1 x ?
             * @return "flattened" B (not voigt notation B)
             */
            moris::Matrix<moris::DDRMat> eval_B_flat();

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
        
            void flatten_normal( const Matrix< DDRMat > & aNormal,
                                       Matrix< DDRMat > & aFlattenedNormal );

//--------------------------------------------------------------------------------------------------------------

            void get_isotropic_thermal_expansion_vector( Matrix< DDRMat > & aTheramlExpansionVector );

//--------------------------------------------------------------------------------------------------------------
            /**
             * set spatial dimensions. Modified from base to set function pointers to the appropriate eval_const()
             * @param[ in ] aSpaceDim a spatial dimension
             */
            void set_space_dim( uint aSpaceDim );

//--------------------------------------------------------------------------------------------------------------
            /**
             * Sets the option to use a modified elasticity model (e.g. plane stress, plane stress, etc.)
             *
             * @param aModelType Linear isotropic elasticity supports combinations of Model_Type::PLANE_STRESS or
             * Model_Type::PLANE_STRAIN, and Model_Type::HYDROSTATIC or
             * Model_Type::DEVIATORIC
             */
            void set_model_type( Model_Type aModelType);

//--------------------------------------------------------------------------------------------------------------
        private:

            /**
             *  Sets the appropriate function pointers based on the current member data for spatial dimensions and
             *  model types
             */
            void set_function_pointers();

//--------------------------------------------------------------------------------------------------------------
            /**
             * Full plane stress tensor
             *              
             * @param aEmod Elastic modulus
             * @param aNu Poisson ratio
             */
            void full_plane_stress(moris::real aEmod, moris::real aNu);

//--------------------------------------------------------------------------------------------------------------
            /**
             * Deviatoric plane stress tensor
             *
             * @param aEmod Elastic modulus
             * @param aNu Poisson ratio
             */
            void deviatoric_plane_stress(moris::real aEmod, moris::real aNu);

//--------------------------------------------------------------------------------------------------------------
            /**
             * Full plane strain tensor
             *               
             * @param aEmod Elastic modulus
             * @param aNu Poisson ratio
             */
            void full_plane_strain(moris::real aEmod, moris::real aNu);

//--------------------------------------------------------------------------------------------------------------
            /**
             * Deviatoric plane strain tensor
             *
             * @param aEmod Elastic modulus
             * @param aNu Poisson ratio
             */
            void deviatoric_plane_strain(moris::real aEmod, moris::real aNu);

//--------------------------------------------------------------------------------------------------------------
            /**
             * Full 3d linear isotropic tensor
             * 
             * @param aEmod Elastic modulus
             * @param aNu Poisson ratio
             */
            void full_3d(moris::real aEmod, moris::real aNu);

//--------------------------------------------------------------------------------------------------------------
            /**
             * Deviatoric 3d tensor
             *
             * @param aEmod Elastic modulus
             * @param aNu Poisson ratio
             */
            void deviatoric_3d(moris::real aEmod, moris::real aNu);

//--------------------------------------------------------------------------------------------------------------

        };
//--------------------------------------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_STRUC_LINEAR_ISOTROPIC_HPP_ */

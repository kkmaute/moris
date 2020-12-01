#ifndef SRC_FEM_CL_FEM_CM_STRUC_NONLINEAR_ISOTROPIC_HPP_
#define SRC_FEM_CL_FEM_CM_STRUC_NONLINEAR_ISOTROPIC_HPP_

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

        class CM_Struc_Nonlinear_Isotropic : public Constitutive_Model
        {

            protected:

                // default local properties
                std::shared_ptr< Property > mPropEMod    = nullptr;
                std::shared_ptr< Property > mPropPoisson = nullptr;

            private:

                // default dof
                MSI::Dof_Type mDofDispl    = MSI::Dof_Type::UX;

                // property type for CM
                enum class CM_Property_Type
                {
                        EMOD,
                        NU,
                        MAX_ENUM
                };

                // function pointers
                void ( CM_Struc_Nonlinear_Isotropic:: * m_eval_strain )() = nullptr;

                void ( CM_Struc_Nonlinear_Isotropic:: * m_eval_teststrain )() = nullptr;

                void ( CM_Struc_Nonlinear_Isotropic:: * m_flatten_normal )(
                        const Matrix< DDRMat > & aNormal,
                        Matrix< DDRMat >       & aFlatNormal ) = nullptr;

                void ( CM_Struc_Nonlinear_Isotropic:: * mConstFunc )( moris::real,
                        moris::real ) = &CM_Struc_Nonlinear_Isotropic::full_3d;

                Model_Type mPlaneType  = Model_Type::PLANE_STRESS; // Plane stress or plane strain, only used in 2d
                Model_Type mTensorType = Model_Type::FULL; // Hydrostatic or deviatoric (default: full tensor)


                //--------------------------------------------------------------------------------------------------------------
            public:
                /*
                 * trivial constructor
                 */
                CM_Struc_Nonlinear_Isotropic();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~CM_Struc_Nonlinear_Isotropic(){};

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

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * set spatial dimensions. Modified from base to set function pointers to the appropriate eval_const()
                 * @param[ in ] aSpaceDim a spatial dimension
                 */
                void set_space_dim( uint aSpaceDim );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * Sets the option to use a modified elasticity model (e.g. plane stress, plane stress, etc.)
                 * @param aModelType Linear isotropic elasticity supports combinations of Model_Type::PLANE_STRESS or
                 * Model_Type::PLANE_STRAIN
                 */
                void set_model_type( Model_Type aModelType );

                //--------------------------------------------------------------------------------------------------------------
            private:

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
                void eval_testTraction(
                        const Matrix< DDRMat >             & aNormal,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the strain
                 */
                void eval_strain()
                {
                    ( this->*m_eval_strain )();
                }

                void eval_strain_2d();
                void eval_strain_3d();

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
                 * evaluate the constitutive model matrix derivative wrt to a dof type
                 * @param[ in ] aDofTypes   a dof type wrt which the derivative is evaluated
                 * @param[ in ] adConstdDOF a matrix to fill with derivative evaluation
                 */
                void eval_dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 *  Sets the appropriate function pointers
                 *  based on the current member data
                 *  for spatial dimensions and model types
                 */
                void set_function_pointers();

                //--------------------------------------------------------------------------------------------------------------

                void flatten_normal_2d(
                        const Matrix< DDRMat > & aNormal,
                        Matrix< DDRMat >       & aFlatNormal );

                void flatten_normal_3d(
                        const Matrix< DDRMat > & aNormal,
                        Matrix< DDRMat >       & aFlatNormal );

                 /**
                 * flatten normal vector
                 * @param[ in ] aNormal          a normal vector
                 * @param[ in ] aFlattenedNormal a matrix for flattened normal to fill
                 */
                void flatten_normal(
                        const Matrix< DDRMat > & aNormal,
                        Matrix< DDRMat >       & aFlatNormal )
                {
                    ( this->*m_flatten_normal )( aNormal, aFlatNormal );
                }                        

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model matrix
                 */
                void eval_const();

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * Full plane stress tensor
                 * @param[ in ] aEmod Elastic modulus
                 * @param[ in ] aNu   Poisson ratio
                 */
                void full_plane_stress(
                        moris::real aEmod,
                        moris::real aNu );

                //--------------------------------------------------------------------------------------------------------------
                /**
                 * Full plane strain tensor
                 * @param aEmod Elastic modulus
                 * @param aNu   Poisson ratio
                 */
                void full_plane_strain(
                        moris::real aEmod,
                        moris::real aNu );


                //--------------------------------------------------------------------------------------------------------------
                /**
                 * Full 3d linear isotropic tensor
                 * @param aEmod Elastic modulus
                 * @param aNu Poisson ratio
                 */
                void full_3d(
                        moris::real aEmod,
                        moris::real aNu ){};

                //--------------------------------------------------------------------------------------------------------------

        };
        //--------------------------------------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_CM_STRUC_NONLINEAR_ISOTROPIC_HPP_ */

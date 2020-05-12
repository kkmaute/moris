/*
 * cl_FEM_Constitutive_Model.hpp
 *
 *  Created on: Sep 17, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_
#define SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
//#include "linalg_typedefs.hpp"              //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Enums.hpp"                 //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include <map>

#include "cl_GEN_Pdv_Enums.hpp"

namespace moris
{
    namespace fem
    {
    class Field_Interpolator_Manager;
    class Set;
//------------------------------------------------------------------------------
        /**
         * Constitutive model
         */
        class Constitutive_Model
        {

        protected :

            // field interpolator manager
            Field_Interpolator_Manager * mFIManager = nullptr;

            // fem set pointer
            Set * mSet = nullptr;

            // dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mDofTypes;

            // local string to dof enum map
            std::map< std::string, MSI::Dof_Type > mDofMap;

            // dof type map
            Matrix< DDSMat > mDofTypeMap;

            // global dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mGlobalDofTypes;

            // global dof type map
            Matrix< DDSMat > mGlobalDofTypeMap;

            // dv type list
            moris::Cell< moris::Cell< PDV > > mDvTypes;

            // local string to dv enum map
            std::map< std::string, PDV > mDvMap;

            // global dv type list
            moris::Cell< moris::Cell< PDV > > mGlobalDvTypes;

            // global dv type map
            Matrix< DDSMat > mGlobalDvTypeMap;

            // dv type map
            Matrix< DDSMat > mDvTypeMap;

            // properties
            moris::Cell< std::shared_ptr< Property > > mProperties;

            // spatial dimensions
            uint mSpaceDim;

            // storage for evaluation
            Matrix< DDRMat > mFlux;
            Matrix< DDRMat > mDivFlux;
            moris::Cell< Matrix< DDRMat > > mdFluxdDof;
            moris::Cell< Matrix< DDRMat > > mdFluxdDv;
            moris::Cell< Matrix< DDRMat > > mdFluxdx;
            moris::Cell< Matrix< DDRMat > > mddivfluxdu;

            Matrix< DDRMat > mTraction;
            moris::Cell< Matrix< DDRMat > > mdTractiondDof;
            moris::Cell< Matrix< DDRMat > > mdTractiondDv;

            moris::Cell< Matrix< DDRMat > > mTestTraction;
            moris::Cell< moris::Cell< Matrix< DDRMat > > > mdTestTractiondDof;
            moris::Cell< moris::Cell< Matrix< DDRMat > > > mdTestTractiondDv;

            Matrix< DDRMat > mStrain;
            Matrix< DDRMat > mDivStrain;
            moris::Cell< Matrix< DDRMat > > mdStraindDof;
            moris::Cell< Matrix< DDRMat > > mdStraindDv;
            moris::Cell< Matrix< DDRMat > > mdStraindx;
            moris::Cell< Matrix< DDRMat > > mddivstraindu;

            Matrix< DDRMat > mTestStrain;

            Matrix< DDRMat > mConst;
            moris::Cell< Matrix< DDRMat > > mdConstdDof;
            moris::Cell< Matrix< DDRMat > > mdConstdDv;

            // constitutive model name for input and debug
            std::string mName;

            // ------------------------------------ //
            // FIXME: Remove once unified.
            // enthalpy change rate
            Matrix< DDRMat > mHdot;
            moris::Cell< Matrix< DDRMat > > mHdotDof;

            // enthalpy gradient change rate
            Matrix< DDRMat > mGradHdot;
            moris::Cell< Matrix< DDRMat > > mGradHdotDof;

            // gradient of Divergence of Flux
            Matrix< DDRMat > mGradDivFlux;
            moris::Cell< Matrix< DDRMat > > mGradDivFluxDof;
            // ------------------------------------ //

        private:

            // bool for global dof type list and map build
            bool mGlobalDofBuild = true;
            bool mGlobalDvBuild  = true;
            bool mGlobalDofMapBuild = true;
            bool mGlobalDvMapBuild  = true;

            // flag for evaluation
            bool mFluxEval = true;
            bool mDivFluxEval = true;
            moris::Cell< bool > mdFluxdDofEval;
            moris::Cell< bool > mdFluxdDvEval;
            moris::Cell< bool > mdFluxdxEval;
            moris::Cell< bool > mddivfluxduEval;

            bool mTractionEval = true;
            moris::Cell< bool > mdTractiondDofEval;
            moris::Cell< bool > mdTractiondDvEval;

            moris::Cell< bool > mTestTractionEval;
            moris::Cell< moris::Cell< bool > > mdTestTractiondDofEval;
            moris::Cell< moris::Cell< bool > > mdTestTractiondDvEval;

            bool mStrainEval = true;
            bool mDivStrainEval = true;
            moris::Cell< bool > mdStraindDofEval;
            moris::Cell< bool > mdStraindDvEval;
            moris::Cell< bool > mdStraindxEval;
            moris::Cell< bool > mddivstrainduEval;

            bool mTestStrainEval = true;

            bool mConstEval = true;
            moris::Cell< bool > mdConstdDofEval;
            moris::Cell< bool > mdConstdDvEval;

            // ------------------------------------ //
            // FIXME: Remove once unified.
            // eval flags specific to this CM
            bool mHdotEval = true;
            bool mGradHdotEval = true;
            bool mGradDivFluxEval = true;

            moris::Cell<bool> mHdotDofEval;
            moris::Cell<bool> mGradHdotDofEval;
            moris::Cell<bool> mGradDivFluxDofEval;
            // ------------------------------------ //

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Constitutive_Model()
            {
                // FIXME for now only 1st order allowed
                uint tOrder = 1;

                // set storage for evaluation
                mdFluxdx.resize( tOrder );
                mdStraindx.resize( tOrder );

                // set flag for evaluation
                mdFluxdxEval.assign( tOrder, true );
                mdStraindxEval.assign( tOrder, true );
            };

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            virtual ~Constitutive_Model(){};

//------------------------------------------------------------------------------
            /**
             * set name
             * param[ in ] aName a string for CM name
             */
            void set_name( std::string aName )
            {
                mName = aName;
            }

//------------------------------------------------------------------------------
            /**
             * get name
             * param[ out ] mName a string for CM name
             */
            std::string get_name()
            {
                return mName;
            }

//------------------------------------------------------------------------------
            /**
             * print names
             */
            void print_names()
            {
                std::cout<<"----------"<<std::endl;
                std::cout<<"CM: "<<mName<<std::endl;

                // properties
                for( uint iProp = 0; iProp < mProperties.size(); iProp++ )
                {
                    if( mProperties( iProp ) != nullptr )
                    {
                        std::cout<<"Property: "<<mProperties( iProp )->get_name()<<std::endl;
                    }
                }
                std::cout<<"----------"<<std::endl;
            }

//------------------------------------------------------------------------------
            /*
             * set member set pointer
             * @param[ in ] aSetPointer a FEM set pointer
             */
            void set_set_pointer( Set * aSetPointer )
            {
                mSet = aSetPointer;
            }

//------------------------------------------------------------------------------
            /**
             * set space dimension
             * @param[ in ] aSpaceDim a spatial dimension
             */
            virtual void set_space_dim( uint aSpaceDim )
            {
                // check that space dimension is 1, 2, 3
                MORIS_ERROR( aSpaceDim > 0 && aSpaceDim < 4, "Constitutive_Model::set_space_dim - wrong space dimension." );

                // set space dimension
                mSpaceDim = aSpaceDim;
            }

//------------------------------------------------------------------------------
            /**
             * reset evaluation flags
             */
            void reset_eval_flags();

//------------------------------------------------------------------------------
            /**
             * set constitutive model dof types
             * @param[ in ] aDofTypes a list of group of dof types
             */
            void set_dof_type_list( moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes );

//------------------------------------------------------------------------------
            /**
             * set constitutive model dof types
             * @param[ in ] aDofTypes a list of group of dof types
             * @param[ in ] aDofStrings a list of strings to describe the dof types
             */
            void set_dof_type_list( moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                                    moris::Cell< std::string >                  aDofStrings );

//------------------------------------------------------------------------------
            /**
             * return a cell of dof types
             * @param[ in ] mDofTypes a cell of cell of dof types
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list() const
            {
                return mDofTypes;
            }

//------------------------------------------------------------------------------
            /**
             * build a map for the dof types
             */
            void build_dof_type_map();

//------------------------------------------------------------------------------
            /**
             * get dof type map
             * @param[ out ] mDofTypeMap map for the dof types
             */
            const Matrix< DDSMat > & get_dof_type_map()
            {
                return mDofTypeMap;
            }

//------------------------------------------------------------------------------
            /**
             * set constitutive model dv types
             * @param[ in ] aDvTypes a list of group of dv types
             */
            void set_dv_type_list( moris::Cell< moris::Cell< PDV > > aDvTypes );

//------------------------------------------------------------------------------
            /**
             * set constitutive model dv types
             * @param[ in ] aDvTypes   a list of group of dv types
             * @param[ in ] aDvStrings a list of strings to describe the dv types
             */
            void set_dv_type_list( moris::Cell< moris::Cell< PDV > > aDvTypes,
                                   moris::Cell< std::string >           aDvStrings );

//------------------------------------------------------------------------------
            /**
             * return a cell of dv types
             * @param[ out ] aDvTypes a cell of cell of dv types
             */
            const moris::Cell< moris::Cell< PDV > > & get_dv_type_list() const
            {
                return mDvTypes;
            };

//------------------------------------------------------------------------------
            /**
             * build a map for the dv types
             */
            void build_dv_type_map();

//------------------------------------------------------------------------------
            /**
             * get dv type map
             * @param[ out ] mDvTypeMap map for the dv types
             */
            const Matrix< DDSMat > & get_dv_type_map()
            {
                return mDvTypeMap;
            }

//------------------------------------------------------------------------------
            /**
             * set field interpolator manager
             * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
             */
            void set_field_interpolator_manager( Field_Interpolator_Manager * aFieldInterpolatorManager );

//------------------------------------------------------------------------------
            /**
             * set model type
             * @param[ in ] aModelType an enum for model type
             */
            virtual void set_model_type( fem::Model_Type aModelType )
            {
                MORIS_ERROR( false, "Constitutive_Model::set_model_type - Not implemented for base class." );
            }

//------------------------------------------------------------------------------
            /**
             * set property
             * @param[ in ] aProperty a property pointer
             * @param[ in ] aPropertyString a string describing the property
             */
            virtual void set_property( std::shared_ptr< fem::Property > aProperty,
                                       std::string                      aPropertyType )
            {
                MORIS_ERROR( false, "Constitutive_Model::set_property - Not implemented for base class." );
            }

//------------------------------------------------------------------------------
            /**
             * get properties
             * @param[ out ] mProperties cell of property pointers
             */
            moris::Cell< std::shared_ptr< Property > > & get_properties()
            {
                return mProperties;
            };

//------------------------------------------------------------------------------
            /**
             * create a global dof type list including constitutive and property dependencies
             */
            void build_global_dof_type_list();

//------------------------------------------------------------------------------
            /**
             * get global dof type list
             * @param[ out ] mGlobalDofTypes global list of dof type
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_type_list();

//------------------------------------------------------------------------------
            /**
             * get non unique list of dof type for the constitutive model
             * @param[ in ] aDofTypes a cell of dof type to fill
             */
            void get_non_unique_dof_types( moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * get non unique lists of dof and dv type for the constitutive model
             * @param[ in ] aDofTypes a cell of dof type to fill
             * @param[ in ] aDvTypes  a cell of dv type to fill
             */
            void get_non_unique_dof_and_dv_types( moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                  moris::Cell< PDV >        & aDvTypes );

//------------------------------------------------------------------------------
            /**
             * build global dof type map
             */
            void build_global_dof_type_map();

//------------------------------------------------------------------------------
            /**
             * get global dof type map
             * @param[ out ] mGlobalDofTypeMap global map for the dof types
             */
            const Matrix< DDSMat > & get_global_dof_type_map();

//------------------------------------------------------------------------------
            /**
             * check dependency on a given group of dof types
             * @param[ in ]  aDofType       a group of dof types
             * @param[ out ] tDofDependency a bool true if dependency on dof type
             *
             */
            bool check_dof_dependency( const moris::Cell< MSI::Dof_Type > & aDofType );

//------------------------------------------------------------------------------
            /**
             * create a global dv type list including constitutive and property dependencies
             */
            void build_global_dv_type_list();

//------------------------------------------------------------------------------
            /**
             * get global dv type list
             * @param[ out ] mGlobalDvTypes global list of dv types
             */
            const moris::Cell< moris::Cell< PDV > > & get_global_dv_type_list();

//------------------------------------------------------------------------------
            /**
             * build global dv type list map
             */
            void build_global_dv_type_map();

//------------------------------------------------------------------------------
            /**
             * get global dv type map
             * @param[ out ] mGlobalDvTypeMap global map for the dv types
             */
            const Matrix< DDSMat > & get_global_dv_type_map()
            {
                return mGlobalDvTypeMap;
            }

//------------------------------------------------------------------------------
            /**
             * check dependency on a given group of dv types
             * @param[ in ]  aDvType       a group of dv types
             * @param[ out ] tDvDependency a bool true if dependency on dv type
             *
             */
            bool check_dv_dependency( const moris::Cell< PDV > & aDvType );

//------------------------------------------------------------------------------
            /**
             * get the constitutive model flux
             * @param[ out ] mFlux constitutive model flux
             */
            const Matrix< DDRMat > & flux();

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux
             */
            virtual void eval_flux()
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_flux - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the divergence of the flux
             * @param[ out ] mDivFlux divergence of the flux
             */
            const Matrix< DDRMat > & divflux();

//------------------------------------------------------------------------------
            /**
             * evaluate the divergence of the flux
             */
            virtual void eval_divflux()
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_divflux - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the divergence of the flux wrt to dof type
             * @param[ out ] mddivfluxdu derivative of the divergence of the flux
             *                           wrt to dof type
             */
            const Matrix< DDRMat > & ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofType );

//------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the divergence of the flux wrt to dof type
             */
            virtual void eval_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_ddivfluxdu - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the constitutive model traction
             * @param[ in ]  aNormal   normal
             * @param[ out ] mTraction constitutive model traction
             */
            const Matrix< DDRMat > & traction( const Matrix< DDRMat > & aNormal );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model traction
             * @param[ in ]  aNormal normal
             */
            virtual void eval_traction( const Matrix< DDRMat > & aNormal )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_traction - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the constitutive model test traction
             * @param[ in ]  aNormal       normal
             * @param[ out ] mTestTraction constitutive model test traction
             */
            const Matrix< DDRMat > & testTraction( const Matrix< DDRMat >             & aNormal,
                                                   const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test traction
             * @param[ in ]  aNormal normal
             */
            virtual void eval_testTraction( const Matrix< DDRMat >             & aNormal,
                                            const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_testTraction - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the constitutive model strain
             * @param[ out ] mStrain constitutive model strain
             */
            const Matrix< DDRMat > & strain();

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain
             */
            virtual void eval_strain()
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_strain - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the divergence of the strain
             * @param[ out ] mDivFlux divergence of the strain
             */
            const Matrix< DDRMat > & divstrain();

//------------------------------------------------------------------------------
            /**
             * evaluate the divergence of the strain
             */
            virtual void eval_divstrain()
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_divstrain - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the divergence of the strain wrt to dof type
             * @param[ out ] mddivstraindu derivative of the divergence of the strain
             *                             wrt to dof type
             */
            const Matrix< DDRMat > & ddivstraindu( const moris::Cell< MSI::Dof_Type > & aDofType );

//------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the divergence of the strain wrt to dof type
             */
            virtual void eval_ddivstraindu( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_ddivstraindu - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the constitutive model test strain
             * @param[ out ] mTestStrain constitutive model test strain
             */
            const Matrix< DDRMat > & testStrain();

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test strain
             */
            virtual void eval_testStrain()
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_testStrain - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the constitutive model constitutive matrix
             * @param[ out ] mConst constitutive matrix
             */
            const Matrix< DDRMat > & constitutive();

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model constitutive matrix
             */
            virtual void eval_const()
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_const - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the flux wrt space
             * @param[ in ] aOrder order of the derivative
             */
            const Matrix< DDRMat > & dfluxdx( uint aOrder );

//------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the flux wrt space
             * @param[ in ] aOrder order of the derivative
             */
            virtual void eval_dfluxdx( uint aOrder )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dfluxdx - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the flux wrt dof
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ out ] mFluxDofDer derivative of the flux wrt dof
             */
            const Matrix< DDRMat > & dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofType );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux derivative wrt to a dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             */
            virtual void eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dFluxdDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model change rate of enthalpy
             */
            virtual void eval_Hdot()
            {
                MORIS_ASSERT(false, "eval_Hdot: not implemented in base class.");
            };

            /**
             * get the constitutive model change rate of enthalpy
             * @param[ out ] mHdot change rate of enthalpy
             */
            const Matrix< DDRMat > & Hdot();

//------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model change rate of spatial gradient of enthalpy (needed for GGLS-stabilization)
             */
            virtual void eval_gradHdot()
            {
                MORIS_ASSERT(false, "eval_gradHdot: not implemented in base class.");
            };

            /**
             * get the constitutive model change rate of spatial gradient of enthalpy (needed for GGLS-stabilization)
             * @param[ out ] mGradHdot gradient of change rate of enthalpy
             */
            const Matrix< DDRMat > & gradHdot();

//------------------------------------------------------------------------------
            /**
             * evaluates the gradient of the divergence of the flux (needed for GGLS-stabilization)
             */
            virtual void eval_graddivflux()
            {
                MORIS_ASSERT(false, "eval_graddivflux: not implemented in base class.");
            };

            /**
             * get the gradient of the divergence of the flux (needed for GGLS-stabilization)
             * @param[ out ] mGradDivFlux gradient of divergence of flux
             */
            const Matrix< DDRMat > & graddivflux();

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model enthalpy change rate wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dHdotdDOF ( 1 x numDerDof )
             */
            virtual void eval_dHdotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
            {
                MORIS_ASSERT(false, "eval_dHdotdDOF: not implemented in base class.");
            };

            /**
             * get the enthalpy change rate wrt dof
             * @param[ in ]  aDofType        group of dof type
             * @param[ out ] mHdotDofDer derivative of the traction wrt dof
             */
            const Matrix< DDRMat > & dHdotdDOF( const moris::Cell< MSI::Dof_Type > & aDofType);

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model gradient of enthalpy change rate wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dgradHdotdDOF ( mSpaceDim x numDerDof )
             */
            virtual void eval_dGradHdotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
            {
                MORIS_ASSERT(false, "eval_dGradHdotdDOF: not implemented in base class.");
            };

            /**
             * get the gradient of enthalpy change rate wrt dof
             * @param[ in ]  aDofType        group of dof type
             * @param[ out ] mGradHdotDer derivative of the traction wrt dof
             */
            const Matrix< DDRMat > & dGradHdotdDOF( const moris::Cell< MSI::Dof_Type > & aDofType);

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model gradient of divergence of flux wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dGradDivFluxdDOF ( mSpaceDim x numDerDof )
             */
            virtual void eval_dGradDivFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
            {
                MORIS_ASSERT(false, "eval_dGradDivFluxdDOF: not implemented in base class.");
            };

            /**
             * get the gradient of enthalpy change rate wrt dof
             * @param[ in ]  aDofType        group of dof type
             * @param[ out ] mGradHdotDer derivative of the traction wrt dof
             */
            const Matrix< DDRMat > & dGradDivFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofType);

//------------------------------------------------------------------------------
            /**
             * get the derivative of the traction wrt dof
             * @param[ in ]  aDofType        group of dof type
             * @param[ in ]  aNormal         normal
             * @param[ out ] mTractionDofDer derivative of the traction wrt dof
             */
            const Matrix< DDRMat > & dTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofType,
                                                    const Matrix< DDRMat >             & aNormal );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model traction derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * @param[ in ] aNormal   normal
             */
            virtual void eval_dTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                             const Matrix< DDRMat >             & aNormal)
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dTractiondDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the test traction wrt dof
             * @param[ in ]  aDofType           group of dof type
             * @param[ in ]  aNormal            normal
             * @param[ out ] mdTestTractiondDof derivative of the traction wrt dof
             */
            const Matrix< DDRMat > & dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofType,
                                                        const Matrix< DDRMat >             & aNormal,
                                                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test traction derivative wrt to a dof type
             * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ in ] adTractiondDOF a matrix to fill with derivative evaluation
             */
            virtual void eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                 const Matrix< DDRMat >             & aNormal,
                                                 const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dTestTractiondDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * FIXME this is not a test traction, used for elast lin iso!!!!
             * get the derivative of the test traction wrt dof
             * @param[ in ]  aDofType           group of dof type
             * @param[ in ]  aNormal            normal
             * @param[ out ] mdTestTractiondDof derivative of the traction wrt dof
             */
            const Matrix< DDRMat > & dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofType,
                                                        const Matrix< DDRMat >             & aNormal,
                                                        const Matrix< DDRMat >             & aJump,
                                                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes );

//------------------------------------------------------------------------------
            /**
             * FIXME this is not a test traction, used for elast lin iso!!!!
             * evaluate the constitutive model test traction derivative wrt to a dof type
             * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ in ] adTractiondDOF a matrix to fill with derivative evaluation
             */
            virtual void eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                 const Matrix< DDRMat >             & aNormal,
                                                 const Matrix< DDRMat >             & aJump,
                                                 const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dTestTractiondDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the strain wrt space
             * @param[ in ] aOrder order of the derivative
             */
            const Matrix< DDRMat > & dstraindx( uint aOrder );

//------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the strain wrt space
             * @param[ in ] aOrder order of the derivative
             */
            virtual void eval_dstraindx( uint aOrder )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dstraindx - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the strain wrt dof
             * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ out ] mdStraindDof derivative of the strain wrt dof
             */
            const Matrix< DDRMat > & dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofType );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            virtual void eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dStraindDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the constitutive matrix wrt dof
             * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
             * @param[ out ] mdConstdDof derivative of the constitutive matrix wrt dof
             */
            const Matrix< DDRMat > & dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofType );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model constitutive matrix derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            virtual void eval_dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dConstdDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the flux wrt dv
             * @param[ in ]  aDvTypes  a dv type wrt which the derivative is evaluated
             * @param[ out ] mdFluxdDv derivative of the flux wrt dv
             */
            const Matrix< DDRMat > & dFluxdDV( const moris::Cell< PDV > & aDvType );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux derivative wrt to a dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             */
            virtual void eval_dFluxdDV( const moris::Cell< PDV > & aDvTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dFluxdDV - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the strain wrt dv
             * @param[ in ]  aDvTypes    a dv type wrt which the derivative is evaluated
             * @param[ out ] mdStraindDv derivative of the strain wrt dv
             */
            const Matrix< DDRMat > & dStraindDV( const moris::Cell< PDV > & aDvType );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain derivative wrt to a dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             */
            virtual void eval_dStraindDV( const moris::Cell< PDV > & aDvTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dStraindDV - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the constitutive matrix wrt dv
             * @param[ in ]  aDvTypes   a dv type wrt which the derivative is evaluated
             * @param[ out ] mdConstdDv derivative of the constitutive matrix wrt dv
             */
            const Matrix< DDRMat > & dConstdDV( const moris::Cell< PDV > & aDvType );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model constitutive matrix derivative wrt to a dv type
             * @param[ in ] aDvTypes a dof type wrt which the derivative is evaluated
             */
            virtual void eval_dConstdDV( const moris::Cell< PDV > & aDvTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dConstdDV - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model stress derivative wrt to a dof type
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] adFluxdDOF_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation   real to perturb for FD
             */
            void eval_dFluxdDOF_FD( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                          Matrix< DDRMat >             & adFluxdDOF_FD,
                                          real                           aPerturbation );

//------------------------------------------------------------------------------
            /**
             * evaluate the enthalpy change rate wrt dof using finite differences
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] adFluxdDOF_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation   real to perturb for FD
             */
            void eval_dHdotdDOF_FD( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                          Matrix< DDRMat >             & adHdotdDOF_FD,
                                          real                           aPerturbation );

//------------------------------------------------------------------------------
            /**
             * evaluate the gradient of enthalpy change rate wrt dof using finite differences
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] adFluxdDOF_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation   real to perturb for FD
             */
            void eval_dGradHdotdDOF_FD( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                              Matrix< DDRMat >             & adGradHdotdDOF_FD,
                                              real                           aPerturbation );

//------------------------------------------------------------------------------
            /**
             * evaluate the gradient of enthalpy change rate wrt dof using finite differences
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] adFluxdDOF_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation   real to perturb for FD
             */
            void eval_dGradDivFluxdDOF_FD( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                 Matrix< DDRMat >             & adGradDivFluxdDOF_FD,
                                                 real                           aPerturbation );

//------------------------------------------------------------------------------
            /**
            * evaluate the constitutive model strain derivative wrt to a dof type
            * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
            * @param[ in ] adStraindDOF_FD a matrix to fill with derivative evaluation
            * @param[ in ] aPerturbation   real to perturb for FD
            */
            void eval_dStraindDOF_FD( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                            Matrix< DDRMat >             & adStraindDOF_FD,
                                            real                           aPerturbation );

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model stress derivative wrt to a dv type
             * @param[ in ] aDofTypes     a dof type wrt which the derivative is evaluated
             * @param[ in ] adFluxdDV_FD  a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation real to perturb for FD
             */

            void eval_dFluxdDV_FD( const moris::Cell< PDV > & aDvTypes,
                                         Matrix< DDRMat >      & adFluxdDV_FD,
                                         real                    aPerturbation );
//------------------------------------------------------------------------------
            /**
            * evaluate the constitutive model strain derivative wrt to a dv type
            * @param[ in ] aDvTypes       a dv type wrt which the derivative is evaluated
            * @param[ in ] adStraindDV_FD a matrix to fill with derivative evaluation
            * @param[ in ] aPerturbation  real to perturb for FD
            */

            void eval_dStraindDV_FD( const moris::Cell< PDV > & aDvTypes,
                                     Matrix< DDRMat >            & adStraindDV_FD,
                                     real                          aPerturbation );

//------------------------------------------------------------------------------
            /*
             * evaluates and returns the value of E' which is used in the evaluation of stress intensity factors
             */
            virtual moris::real get_e_prime(  )
            {
                MORIS_ERROR( false, " Constitutive_Model::get_e_prime - This function does nothing. " );
                return 0;
            }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_ */

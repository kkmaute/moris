/*
 * cl_FEM_Constitutive_Model.hpp
 *
 *  Created on: Sep 17, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_
#define SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_

//MRS/COR/src
#include "typedefs.hpp"
//#include "linalg_typedefs.hpp"
//MRS/CON/src
#include "cl_Cell.hpp"
//LNA/src
#include "cl_Matrix.hpp"
//FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Enums.hpp"
#include "fn_FEM_FD_Scheme.hpp"
#include "cl_FEM_Material_Model.hpp"
//FEM/MSI/src
#include "cl_MSI_Dof_Type_Enums.hpp"
#include <map>
//GEN/src
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

                // dof type map
                Matrix< DDSMat > mDofTypeMap;

                // global dof type list
                moris::Cell< moris::Cell< MSI::Dof_Type > > mGlobalDofTypes;

                // global dof type map
                Matrix< DDSMat > mGlobalDofTypeMap;

                // dv type list
                moris::Cell< moris::Cell< PDV_Type > > mDvTypes;

                // local string to dv enum map
                std::map< std::string, PDV_Type > mDvMap;

                // global dv type list
                moris::Cell< moris::Cell< PDV_Type > > mGlobalDvTypes;

                // global dv type map
                Matrix< DDSMat > mGlobalDvTypeMap;

                // dv type map
                Matrix< DDSMat > mDvTypeMap;

                // dof type list
                moris::Cell< moris::Cell< mtk::Field_Type > > mFieldTypes;

                // dof type map
                Matrix< DDSMat > mFieldTypeMap;

                // global dof type list
                moris::Cell< moris::Cell< mtk::Field_Type > > mGlobalFieldTypes;

                // global dof type map
                Matrix< DDSMat > mGlobalFieldTypeMap;

                // properties
                moris::Cell< std::shared_ptr< Property > > mProperties;

                // local string to property enum map
                std::map< std::string, uint > mPropertyMap;

                // Material Model
                moris::Cell< std::shared_ptr< Material_Model > > mMaterialModels;

                // local string to property enum map
                std::map< std::string, uint > mMaterialModelMap;

                // spatial dimensions
                uint mSpaceDim;

                // storage for flux evaluation
                Matrix< DDRMat > mFlux;
                moris::Cell< Matrix< DDRMat > > mdFluxdDof;
                moris::Cell< Matrix< DDRMat > > mdFluxdDv;
                moris::Cell< Matrix< DDRMat > > mdFluxdx;

                // storage for divergence of flux evaluation
                Matrix< DDRMat > mDivFlux;
                moris::Cell< Matrix< DDRMat > > mddivfluxdu;

                // storage for energy evaluation
                Matrix< DDRMat > mEnergy;
                moris::Cell< Matrix< DDRMat > > mEnergyDof;

                // storage for gradient of energy evaluation
                Matrix< DDRMat > mGradEnergy;
                moris::Cell< Matrix< DDRMat > > mGradEnergyDof;

                // storage for energy rate evaluation
                Matrix< DDRMat > mEnergyDot;
                moris::Cell< Matrix< DDRMat > > mEnergyDotDof;

                // storage for gradient of energy rate evaluation
                Matrix< DDRMat > mGradEnergyDot;
                moris::Cell< Matrix< DDRMat > > mGradEnergyDotDof;

                // storage for gradient of div flux rate evaluation
                Matrix< DDRMat > mGradDivFlux;
                moris::Cell< Matrix< DDRMat > > mGradDivFluxDof;

                // storage for traction evaluation
                Matrix< DDRMat > mTraction;
                moris::Cell< Matrix< DDRMat > > mdTractiondDof;
                moris::Cell< Matrix< DDRMat > > mdTractiondDv;

                // storage for test traction evaluation
                moris::Cell< Matrix< DDRMat > > mTestTraction;
                moris::Cell< moris::Cell< Matrix< DDRMat > > > mdTestTractiondDof;
                moris::Cell< moris::Cell< Matrix< DDRMat > > > mdTestTractiondDv;

                // storage for stress evaluation
                Matrix< DDRMat > mStress;
                moris::Cell< Matrix< DDRMat > > mdStressdDof;

                // storage for strain evaluation
                Matrix< DDRMat > mStrain;
                moris::Cell< Matrix< DDRMat > > mdStraindDof;
                moris::Cell< Matrix< DDRMat > > mdStraindDv;
                moris::Cell< Matrix< DDRMat > > mdStraindx;

                // storage for divergence of strain evaluation
                Matrix< DDRMat > mDivStrain;
                moris::Cell< Matrix< DDRMat > > mddivstraindu;

                // storage for test strain evaluation
                Matrix< DDRMat > mTestStrain;
                Matrix< DDRMat > mTestStrainTrans;

                // storage for constitutive matrix evaluation
                Matrix< DDRMat > mConst;
                moris::Cell< Matrix< DDRMat > > mdConstdDof;
                moris::Cell< Matrix< DDRMat > > mdConstdDv;

                // constitutive model name for input and debug
                std::string mName = "Undefined";

            private:

                // bool for global dof type list and map build
                bool mGlobalDofBuild    = true;
                bool mGlobalDvBuild     = true;
                bool mGlobalFieldBuild  = true;
                bool mGlobalDofMapBuild    = true;
                bool mGlobalDvMapBuild     = true;
                bool mGlobalFieldMapBuild  = true;

                // flag for flux related evaluation
                bool mFluxEval = true;
                moris::Matrix< DDBMat > mdFluxdDofEval;
                moris::Matrix< DDBMat > mdFluxdDvEval;
                moris::Matrix< DDBMat > mdFluxdxEval;

                // flag for div flux related evaluation
                bool mDivFluxEval = true;
                moris::Matrix< DDBMat > mddivfluxduEval;

                // flag for grad div flux related evaluation
                bool mGradDivFluxEval = true;
                moris::Matrix< DDBMat > mGradDivFluxDofEval;

                // flag for traction related evaluation
                bool mTractionEval = true;
                moris::Matrix< DDBMat > mdTractiondDofEval;
                moris::Matrix< DDBMat > mdTractiondDvEval;

                // flag for test traction related evaluation
                moris::Matrix< DDBMat > mTestTractionEval;
                moris::Matrix< DDBMat > mdTestTractiondDofEval;
                moris::Matrix< DDBMat > mdTestTractiondDvEval;

                // flag for stress related evaluation
                bool mStressEval = true;
                moris::Matrix< DDBMat > mdStressdDofEval;

                // flag for strain related evaluation
                bool mStrainEval = true;
                moris::Matrix< DDBMat > mdStraindDofEval;
                moris::Matrix< DDBMat > mdStraindDvEval;
                moris::Matrix< DDBMat > mdStraindxEval;

                // flag for div strain related evaluation
                bool mDivStrainEval = true;
                moris::Matrix< DDBMat > mddivstrainduEval;

                // flag for test strain related evaluation
                bool mTestStrainEval      = true;
                bool mTestStrainTransEval = true;

                // flag for constitutive matrix related evaluation
                bool mConstEval = true;
                moris::Matrix< DDBMat > mdConstdDofEval;
                moris::Matrix< DDBMat > mdConstdDvEval;

                // flag for energy related evaluation
                bool mEnergyEval = true;
                moris::Matrix< DDBMat > mEnergyDofEval;

                // flag for energy rate related evaluation
                bool mEnergyDotEval = true;
                moris::Matrix< DDBMat > mEnergyDotDofEval;

                // flag for grad energy related evaluation
                bool mGradEnergyEval = true;
                moris::Matrix< DDBMat > mGradEnergyDofEval;

                // flag for grad energy rate related evaluation
                bool mGradEnergyDotEval = true;
                moris::Matrix< DDBMat > mGradEnergyDotDofEval;

                //------------------------------------------------------------------------------
            public :

                //------------------------------------------------------------------------------
                /**
                 * constructor
                 */
                Constitutive_Model();

                //------------------------------------------------------------------------------
                /**
                 * virtual destructor
                 */
                virtual ~Constitutive_Model(){};

                //------------------------------------------------------------------------------
                /**
                 * @return constitutive type
                 */
                virtual
                Constitutive_Type get_constitutive_type() const
                {
                    // need to define this for every CM
                    return Constitutive_Type::UNDEFINED;
                }

                //------------------------------------------------------------------------------
                /**
                 * set property
                 * @param[ in ] aProperty       property shared pointer to set
                 * @param[ in ] aPropertyString string describing the property to set
                 */
                void set_property(
                        std::shared_ptr< fem::Property > aProperty,
                        std::string                      aPropertyString );

                //------------------------------------------------------------------------------
                /**
                 * get property
                 * @param[ in ]  aPropertyString string describing the property to get
                 * @param[ out ] tProperty       property shared pointer
                 */
                std::shared_ptr< fem::Property > & get_property( std::string aPropertyString );

                //------------------------------------------------------------------------------
                /**
                 * set thermodynamic material model
                 * @param[ in ] aMaterialModel       material model shared pointer to set
                 * @param[ in ] aMaterialModelString string describing the material model to set
                 */
                void set_material_model(
                        std::shared_ptr< fem::Material_Model > aMaterialModel,
                        std::string                            aMaterialModelString );                

                //------------------------------------------------------------------------------
                /**
                 * get thermodynamic material model
                 * @param[ in ]  aMaterialModelString   string describing the Material Model to get
                 * @param[ out ] tMaterialModel         material model shared pointer
                 */
                std::shared_ptr< fem::Material_Model > & get_material_model( std::string aMaterialModelString );                

                //------------------------------------------------------------------------------
                /**
                 * set local properties
                 */
                virtual void set_local_properties(){}

                //------------------------------------------------------------------------------
                /**
                 * set local material model
                 */
                virtual void set_local_material_model(){}

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
                void print_names();

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
                    MORIS_ERROR(
                            aSpaceDim > 0 && aSpaceDim < 4,
                            "Constitutive_Model::set_space_dim - wrong space dimension." );

                    // set space dimension
                    mSpaceDim = aSpaceDim;
                }

                //------------------------------------------------------------------------------
                /**
                 * get number of space dimensions
                 * @param[ in ] NumSpaceDim number of spatial dimensions
                 */
                uint get_num_space_dims()
                {
                    // check that space dimension is 1, 2, 3
                    MORIS_ASSERT( mSpaceDim > 0,
                            "Constitutive_Model::get_num_space_dims() - number of spatial dimensions has not been set." );

                    // return number of space dimensions
                    return mSpaceDim;
                }

                //------------------------------------------------------------------------------
                /**
                 * reset evaluation flags
                 */
                void reset_eval_flags();

                //------------------------------------------------------------------------------
                /**
                 * reset evaluation flags specific to certain constitutive models
                 */
                virtual void reset_specific_eval_flags(){};

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
                virtual void set_dof_type_list(
                        moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                        moris::Cell< std::string >                  aDofStrings )
                {
                    MORIS_ERROR( false, "Constitutive_Model::set_dof_type_list - Not implemented for base class." );
                }

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
                void set_dv_type_list( moris::Cell< moris::Cell< PDV_Type > > aDvTypes );

                //------------------------------------------------------------------------------
                /**
                 * set constitutive model dv types
                 * @param[ in ] aDvTypes   a list of group of dv types
                 * @param[ in ] aDvStrings a list of strings to describe the dv types
                 */
                virtual void set_dv_type_list(
                        moris::Cell< moris::Cell< PDV_Type > > aDvTypes,
                        moris::Cell< std::string >             aDvStrings )
                {
                    MORIS_ERROR( false, "Constitutive_Model::set_model_type - Not implemented for base class." );
                }

                //------------------------------------------------------------------------------
                /**
                 * return a cell of dv types
                 * @param[ out ] aDvTypes a cell of cell of dv types
                 */
                const moris::Cell< moris::Cell< PDV_Type > > & get_dv_type_list() const
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
                 * set constitutive model field types
                 * @param[ in ] aFieldTypes a list of group of field types
                 */
                void set_field_type_list( moris::Cell< moris::Cell< mtk::Field_Type > > aFieldTypes );

                //------------------------------------------------------------------------------
                /**
                 * return a cell of field types
                 * @param[ in ] mFieldTypes a cell of cell of field types
                 */
                const moris::Cell< moris::Cell< mtk::Field_Type > > & get_field_type_list() const
                {
                    return mFieldTypes;
                }

                //------------------------------------------------------------------------------
                /**
                 * build a map for the field types
                 */
                void build_field_type_map();

                //------------------------------------------------------------------------------
                /**
                 * get field type map
                 * @param[ out ] mFieldTypeMap map for the field types
                 */
                const Matrix< DDSMat > & get_field_type_map()
                {
                    return mFieldTypeMap;
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
                 * @return plane type
                 */
                virtual
                Model_Type get_plane_type() const
                {
                    return Model_Type::UNDEFINED;
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
                 * get material model
                 * @param[ out ] mMaterialModel cell of property pointers
                 */
                moris::Cell< std::shared_ptr< Material_Model > > & get_material_models()
                {
                    return mMaterialModels;
                };

                //------------------------------------------------------------------------------
                /**
                 * create a global dof type list including constitutive and property dependencies
                 */
                void build_global_dof_type_list();

                //------------------------------------------------------------------------------
                /**
                 * initialize storage variables and evaluation flags specific to some child CMs
                 * function is called in the build_global_dof_type_list()
                 */
                virtual void initialize_spec_storage_vars_and_eval_flags(){};

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
                 * @param[ in ] aDofTypes   a cell of dof type to fill
                 * @param[ in ] aDvTypes    a cell of dv type to fill
                 * @param[ in ] aFieldTypes a cell of field type to fill
                 */
                void get_non_unique_dof_dv_and_field_types(
                        moris::Cell< MSI::Dof_Type >   & aDofTypes,
                        moris::Cell< PDV_Type >        & aDvTypes,
                        moris::Cell< mtk::Field_Type > & aFieldTypes );

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
                const moris::Cell< moris::Cell< PDV_Type > > & get_global_dv_type_list();

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
                bool check_dv_dependency( const moris::Cell< PDV_Type > & aDvType );

                //------------------------------------------------------------------------------
                /**
                 * create a global field type list including constitutive and property dependencies
                 */
                void build_global_field_type_list();

                //------------------------------------------------------------------------------
                /**
                 * get global field type list
                 * @param[ out ] mGlobalFieldTypes global list of field type
                 */
                const moris::Cell< moris::Cell< mtk::Field_Type > > & get_global_field_type_list();

                //------------------------------------------------------------------------------
                /**
                 * build global field type map
                 */
                void build_global_field_type_map();

                //------------------------------------------------------------------------------
                /**
                 * get global field type map
                 * @param[ out ] mGlobalFieldTypeMap global map for the field types
                 */
                const Matrix< DDSMat > & get_global_field_type_map();

                //------------------------------------------------------------------------------
                /**
                 * check dependency on a given group of field types
                 * @param[ in ]  aFieldType       a group of field types
                 * @param[ out ] tFieldDependency a bool true if dependency on field type
                 *
                 */
                bool check_field_dependency( const moris::Cell< mtk::Field_Type > & aFieldType );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model flux
                 */
                virtual void eval_flux()
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_flux - This function does nothing. " );
                }

                /**
                 * get the constitutive model flux
                 * @param[ in ]  aCMFunctionType  enum indicating which flux is called, if there are several
                 * @param[ out ] mFlux constitutive model flux
                 */
                virtual const Matrix< DDRMat > & flux( enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the divergence of the flux
                 */
                virtual void eval_divflux()
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_divflux - This function does nothing. " );
                }

                /**
                 * get the divergence of the flux
                 * @param[ out ] mDivFlux divergence of the flux
                 */
                virtual const Matrix< DDRMat > & divflux( enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the divergence of the flux wrt to dof type
                 */
                virtual void eval_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofType )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_ddivfluxdu - This function does nothing. " );
                }

                /**
                 * get the derivative of the divergence of the flux wrt to dof type
                 * @param[ out ] mddivfluxdu derivative of the divergence of the flux
                 *                           wrt to dof type
                 */
                virtual const Matrix< DDRMat > & ddivfluxdu(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model traction
                 * @param[ in ]  aNormal normal
                 */
                virtual void eval_traction( const Matrix< DDRMat > & aNormal )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_traction - This function does nothing. " );
                }

                /**
                 * get the constitutive model traction
                 * @param[ in ]  aNormal   normal
                 * @param[ out ] mTraction constitutive model traction
                 */
                virtual const Matrix< DDRMat > & traction(
                        const Matrix< DDRMat > & aNormal,
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model test traction
                 * @param[ in ]  aNormal normal
                 */
                virtual void eval_testTraction(
                        const Matrix< DDRMat >             & aNormal,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_testTraction - This function does nothing. " );
                }

                /**
                 * get the constitutive model test traction
                 * @param[ in ]  aNormal       normal
                 * @param[ out ] mTestTraction constitutive model test traction
                 */
                virtual const Matrix< DDRMat > & testTraction(
                        const Matrix< DDRMat >             & aNormal,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT);

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model stress
                 */
                virtual void eval_stress()
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_stress - This function does nothing. " );
                }

                /**
                 * get the constitutive model stress
                 * @param[ out ] mStress constitutive model stress
                 */
                virtual const Matrix< DDRMat > & stress(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model strain
                 */
                virtual void eval_strain()
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_strain - This function does nothing. " );
                }

                /**
                 * get the constitutive model strain
                 * @param[ out ] mStrain constitutive model strain
                 */
                virtual const Matrix< DDRMat > & strain(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the divergence of the strain
                 */
                virtual void eval_divstrain()
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_divstrain - This function does nothing. " );
                }

                /**
                 * get the divergence of the strain
                 * @param[ out ] mDivFlux divergence of the strain
                 */
                virtual const Matrix< DDRMat > & divstrain(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the divergence of the strain wrt to dof type
                 */
                virtual void eval_ddivstraindu( const moris::Cell< MSI::Dof_Type > & aDofType )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_ddivstraindu - This function does nothing. " );
                }

                /**
                 * get the derivative of the divergence of the strain wrt to dof type
                 * @param[ out ] mddivstraindu derivative of the divergence of the strain
                 *                             wrt to dof type
                 */
                virtual const Matrix< DDRMat > & ddivstraindu(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model test strain
                 */
                virtual void eval_testStrain()
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_testStrain - This function does nothing. " );
                }

                /**
                 * get the constitutive model test strain
                 * @param[ out ] mTestStrain constitutive model test strain
                 */
                virtual const Matrix< DDRMat > & testStrain(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                /**
                 * get the transpose of the constitutive model test strain
                 * @param[ out ] mTestStrain transpose of constitutive model test strain
                 */
                virtual const Matrix< DDRMat > & testStrain_trans(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model constitutive matrix
                 */
                virtual void eval_const()
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_const - This function does nothing. " );
                }

                /**
                 * get the constitutive model constitutive matrix
                 * @param[ out ] mConst constitutive matrix
                 */
                virtual const Matrix< DDRMat > & constitutive(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the flux wrt space
                 * @param[ in ] aOrder order of the derivative
                 */
                virtual void eval_dfluxdx( uint aOrder )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_dfluxdx - This function does nothing. " );
                }

                /**
                 * get the derivative of the flux wrt space
                 * @param[ in ] aOrder order of the derivative
                 */
                virtual const Matrix< DDRMat > & dfluxdx(
                        uint                  aOrder,
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model flux derivative wrt to a dof type
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 */
                virtual void eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_dFluxdDOF - This function does nothing. " );
                }

                /**
                 * get the derivative of the flux wrt dof
                 * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
                 * @param[ out ] mFluxDofDer derivative of the flux wrt dof
                 */
                virtual const Matrix< DDRMat > & dFluxdDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluates the constitutive model energy
                 */
                virtual void eval_Energy()
                {
                    MORIS_ASSERT(false, "eval_Energy: not implemented in base class.");
                };

                /**
                 * get the constitutive model change rate of energy
                 * @param[ out ] mEnergyDot change rate of energy
                 */
                virtual const Matrix< DDRMat > & Energy( enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluates the constitutive model change rate of energy
                 */
                virtual void eval_EnergyDot()
                {
                    MORIS_ASSERT(false, "eval_EnergyDot: not implemented in base class.");
                };

                /**
                 * get the constitutive model change rate of enthalpy
                 * @param[ out ] mEnergyDot change rate of enthalpy
                 */
                virtual const Matrix< DDRMat > & EnergyDot(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluates the constitutive model spatial gradient of enthalpy
                 */
                virtual void eval_gradEnergy()
                {
                    MORIS_ASSERT(false, "eval_gradEnergy: not implemented in base class.");
                };

                /**
                 * get the constitutive model spatial gradient of enthalpy
                 * @param[ out ] mGradEnergy gradient of enthalpy
                 */
                virtual const Matrix< DDRMat > & gradEnergy(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluates the constitutive model change rate of spatial gradient of enthalpy (needed for GGLS-stabilization)
                 */
                virtual void eval_gradEnergyDot()
                {
                    MORIS_ASSERT(false, "eval_gradEnergyDot: not implemented in base class.");
                };

                /**
                 * get the constitutive model change rate of spatial gradient of enthalpy (needed for GGLS-stabilization)
                 * @param[ out ] mGradEnergyDot gradient of change rate of enthalpy
                 */
                virtual const Matrix< DDRMat > & gradEnergyDot(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

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
                virtual const Matrix< DDRMat > & graddivflux(
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model energy change rate wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dEnergyDotdDOF ( 1 x numDerDof )
                 */
                virtual void eval_dEnergydDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ASSERT(false, "eval_dEnergydDOF: not implemented in base class.");
                };

                /**
                 * get the energy change rate wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ out ] mEnergyDotDofDer derivative of the traction wrt dof
                 */
                virtual const Matrix< DDRMat > & dEnergydDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model enthalpy change rate wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dEnergyDotdDOF ( 1 x numDerDof )
                 */
                virtual void eval_dEnergyDotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ASSERT(false, "eval_dEnergyDotdDOF: not implemented in base class.");
                };

                /**
                 * get the enthalpy change rate wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ out ] mEnergyDotDofDer derivative of the traction wrt dof
                 */
                virtual const Matrix< DDRMat > & dEnergyDotdDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model gradient of enthalpy wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dgradEnergydDOF ( mSpaceDim x numDerDof )
                 */
                virtual void eval_dGradEnergydDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ASSERT(false, "eval_dGradEnergydDOF: not implemented in base class.");
                };

                /**
                 * get the gradient of enthalpy wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ out ] mGradEnergyDer derivative of the gradient of enthalpy wrt dof
                 */
                virtual const Matrix< DDRMat > & dGradEnergydDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model gradient of enthalpy change rate wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * dgradEnergyDotdDOF ( mSpaceDim x numDerDof )
                 */
                virtual void eval_dGradEnergyDotdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ASSERT(false, "eval_dGradEnergyDotdDOF: not implemented in base class.");
                };

                /**
                 * get the gradient of enthalpy change rate wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ out ] mGradEnergyDotDer derivative of the traction wrt dof
                 */
                virtual const Matrix< DDRMat > & dGradEnergyDotdDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

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
                 * @param[ out ] mGradEnergyDotDer derivative of the traction wrt dof
                 */
                virtual const Matrix< DDRMat > & dGradDivFluxdDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 * @param[ in ] aNormal   normal
                 */
                virtual void eval_dTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_dTractiondDOF - This function does nothing. " );
                }

                /**
                 * get the derivative of the traction wrt dof
                 * @param[ in ]  aDofType        group of dof type
                 * @param[ in ]  aNormal         normal
                 * @param[ out ] mTractionDofDer derivative of the traction wrt dof
                 */
                virtual const Matrix< DDRMat > & dTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        const Matrix< DDRMat >             & aNormal,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model test traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
                 * @param[ in ] adTractiondDOF a matrix to fill with derivative evaluation
                 */
                virtual void eval_dTestTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_dTestTractiondDOF - This function does nothing. " );
                }

                /**
                 * get the derivative of the test traction wrt dof
                 * @param[ in ]  aDofType           group of dof type
                 * @param[ in ]  aNormal            normal
                 * @param[ out ] mdTestTractiondDof derivative of the traction wrt dof
                 */
                virtual const Matrix< DDRMat > & dTestTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        const Matrix< DDRMat >             & aNormal,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * FIXME this is not a test traction, used for elast lin iso!!!!
                 * evaluate the constitutive model test traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
                 * @param[ in ] adTractiondDOF a matrix to fill with derivative evaluation
                 */
                virtual void eval_dTestTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const Matrix< DDRMat >             & aNormal,
                        const Matrix< DDRMat >             & aJump,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_dTestTractiondDOF - This function does nothing. " );
                }

                /**
                 * FIXME this is not a test traction, used for elast lin iso!!!!
                 * get the derivative of the test traction wrt dof
                 * @param[ in ]  aDofType           group of dof type
                 * @param[ in ]  aNormal            normal
                 * @param[ out ] mdTestTractiondDof derivative of the traction wrt dof
                 */
                virtual const Matrix< DDRMat > & dTestTractiondDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        const Matrix< DDRMat >             & aNormal,
                        const Matrix< DDRMat >             & aJump,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the strain wrt space
                 * @param[ in ] aOrder order of the derivative
                 */
                virtual void eval_dstraindx( uint aOrder )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_dstraindx - This function does nothing. " );
                }

                /**
                 * get the derivative of the strain wrt space
                 * @param[ in ] aOrder order of the derivative
                 */
                virtual const Matrix< DDRMat > & dstraindx(
                        uint                  aOrder,
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model stress derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 */
                virtual void eval_dStressdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_dStressdDOF - This function does nothing. " );
                }

                /**
                 * get the derivative of the stress wrt dof
                 * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
                 * @param[ out ] mdStressdDof derivative of the stress wrt dof
                 */
                virtual const Matrix< DDRMat > & dStressdDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );


                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model strain derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 */
                virtual void eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_dStraindDOF - This function does nothing. " );
                }

                /**
                 * get the derivative of the strain wrt dof
                 * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
                 * @param[ out ] mdStraindDof derivative of the strain wrt dof
                 */
                virtual const Matrix< DDRMat > & dStraindDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model constitutive matrix derivative wrt to a dof type
                 * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
                 */
                virtual void eval_dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_dConstdDOF - This function does nothing. " );
                }

                /**
                 * get the derivative of the constitutive matrix wrt dof
                 * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
                 * @param[ out ] mdConstdDof derivative of the constitutive matrix wrt dof
                 */
                virtual const Matrix< DDRMat > & dConstdDOF(
                        const moris::Cell< MSI::Dof_Type > & aDofType,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model flux derivative wrt to a dv type
                 * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
                 */
                virtual void eval_dFluxdDV( const moris::Cell< PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_dFluxdDV - This function does nothing. " );
                }

                /**
                 * get the derivative of the flux wrt dv
                 * @param[ in ]  aDvTypes  a dv type wrt which the derivative is evaluated
                 * @param[ out ] mdFluxdDv derivative of the flux wrt dv
                 */
                virtual const Matrix< DDRMat > & dFluxdDV(
                        const moris::Cell< PDV_Type > & aDvType,
                        enum CM_Function_Type           aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model strain derivative wrt to a dv type
                 * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
                 */
                virtual void eval_dStraindDV( const moris::Cell< PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_dStraindDV - This function does nothing. " );
                }

                /**
                 * get the derivative of the strain wrt dv
                 * @param[ in ]  aDvTypes    a dv type wrt which the derivative is evaluated
                 * @param[ out ] mdStraindDv derivative of the strain wrt dv
                 */
                virtual const Matrix< DDRMat > & dStraindDV(
                        const moris::Cell< PDV_Type > & aDvType,
                        enum CM_Function_Type           aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model constitutive matrix derivative wrt to a dv type
                 * @param[ in ] aDvTypes a dof type wrt which the derivative is evaluated
                 */
                virtual void eval_dConstdDV( const moris::Cell< PDV_Type > & aDvTypes )
                {
                    MORIS_ERROR( false, " Constitutive_Model::eval_dConstdDV - This function does nothing. " );
                }

                /**
                 * get the derivative of the constitutive matrix wrt dv
                 * @param[ in ]  aDvTypes   a dv type wrt which the derivative is evaluated
                 * @param[ out ] mdConstdDv derivative of the constitutive matrix wrt dv
                 */
                const Matrix< DDRMat > & dConstdDV(
                        const moris::Cell< PDV_Type > & aDvType,
                        enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model stress derivative wrt to a dof type
                 * @param[ in ] aDofTypes       dof type wrt which the derivative is evaluated
                 * @param[ in ] adFluxdDOF_FD   matrix to fill with derivative evaluation
                 * @param[ in ] aPerturbation   real to perturb for FD
                 * @param[ in ] aFDSchemeType   enum for FD scheme
                 */
                virtual void eval_dFluxdDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adFluxdDOF_FD,
                        real                                 aPerturbation,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
                 * @param[ in ] adtractiondu_FD a matrix to fill with derivative evaluation
                 * @param[ in ] aPerturbation   a real to perturb for FD
                 * @param[ in ] aNormal         a normal
                 */
                void eval_dtractiondu_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adtractiondu_FD,
                        real                                 aPerturbation,
                        Matrix< DDRMat >                   & aNormal,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the test traction derivative wrt to a dof type
                 * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
                 * @param[ in ] aTestDofTypes   a test dof type wrt which the test traction is evaluated
                 * @param[ in ] adtestractiondu_FD a matrix to fill with derivative evaluation
                 * @param[ in ] aPerturbation   a real to perturb for FD
                 * @param[ in ] aNormal         a normal
                 * @param[ in ] aJump         a jump
                 */
                void eval_dtesttractiondu_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        const moris::Cell< MSI::Dof_Type > & aTestDofTypes,
                        Matrix< DDRMat >                   & adtesttractiondu_FD,
                        real                                 aPerturbation,
                        const Matrix< DDRMat >             & aNormal,
                        const Matrix< DDRMat >             & aJump,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the div flux derivative wrt to a dof type
                 * @param[ in ] aDofTypes     a dof type wrt which the derivative is evaluated
                 * @param[ in ] ddivfluxdu_FD a matrix to fill with derivative evaluation
                 * @param[ in ] aPerturbation a real to perturb for FD
                 */
                void eval_ddivfluxdu_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & ddivfluxdu_FD,
                        real                                 aPerturbation,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the div strain derivative wrt to a dof type
                 * @param[ in ] aDofTypes        a dof type wrt which the derivative is evaluated
                 * @param[ in ] addivstraindu_FD a matrix to fill with derivative evaluation
                 * @param[ in ] aPerturbation    a real to perturb for FD
                 */
                void eval_ddivstraindu_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & addivstraindu_FD,
                        real                                 aPerturbation,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the internal energy / enthalpy wrt dof using finite differences
                 * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
                 * @param[ in ] adEnergydDOF_FD a matrix to fill with derivative evaluation
                 * @param[ in ] aPerturbation   real to perturb for FD
                 */
                void eval_dEnergydDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adEnergydDOF_FD,
                        real                                 aPerturbation,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the internal energy / enthalpy change rate wrt dof using finite differences
                 * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
                 * @param[ in ] adEnergyDotdDOF_FD a matrix to fill with derivative evaluation
                 * @param[ in ] aPerturbation   real to perturb for FD
                 */
                void eval_dEnergyDotdDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adEnergyDotdDOF_FD,
                        real                                 aPerturbation,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the gradient of enthalpy wrt dof using finite differences
                 * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
                 * @param[ in ] adGradEnergydDOF_FD  a matrix to fill with derivative evaluation
                 * @param[ in ] aPerturbation   real to perturb for FD
                 */
                void eval_dGradEnergydDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adGradEnergydDOF_FD,
                        real                                 aPerturbation,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the gradient of enthalpy change rate wrt dof using finite differences
                 * @param[ in ] aDofTypes          a dof type wrt which the derivative is evaluated
                 * @param[ in ] adGradEnergyDotdDOF_FD  a matrix to fill with derivative evaluation
                 * @param[ in ] aPerturbation      real to perturb for FD
                 */
                void eval_dGradEnergyDotdDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adGradEnergyDotdDOF_FD,
                        real                                 aPerturbation,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the gradient of enthalpy change rate wrt dof using finite differences
                 * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
                 * @param[ in ] adFluxdDOF_FD a matrix to fill with derivative evaluation
                 * @param[ in ] aPerturbation   real to perturb for FD
                 */
                void eval_dGradDivFluxdDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adGradDivFluxdDOF_FD,
                        real                                 aPerturbation,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model strain derivative wrt to a dof type
                 * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
                 * @param[ in ] adStraindDOF_FD a matrix to fill with derivative evaluation
                 * @param[ in ] aPerturbation   real to perturb for FD
                 */
                void eval_dStraindDOF_FD(
                        const moris::Cell< MSI::Dof_Type > & aDofTypes,
                        Matrix< DDRMat >                   & adStraindDOF_FD,
                        real                                 aPerturbation,
                        fem::FDScheme_Type                   aFDSchemeType = fem::FDScheme_Type::POINT_5,
                        enum CM_Function_Type                aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model stress derivative wrt to a dv type
                 * @param[ in ] aDofTypes     a dof type wrt which the derivative is evaluated
                 * @param[ in ] adFluxdDV_FD  a matrix to fill with derivative evaluation
                 * @param[ in ] aPerturbation real to perturb for FD
                 */

                void eval_dFluxdDV_FD(
                        const moris::Cell< PDV_Type > & aDvTypes,
                        Matrix< DDRMat >              & adFluxdDV_FD,
                        real                            aPerturbation,
                        fem::FDScheme_Type              aFDSchemeType = fem::FDScheme_Type::POINT_5,
                        enum CM_Function_Type           aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the constitutive model strain derivative wrt to a dv type
                 * @param[ in ] aDvTypes       a dv type wrt which the derivative is evaluated
                 * @param[ in ] adStraindDV_FD a matrix to fill with derivative evaluation
                 * @param[ in ] aPerturbation  real to perturb for FD
                 */

                void eval_dStraindDV_FD(
                        const moris::Cell< PDV_Type > & aDvTypes,
                        Matrix< DDRMat >              & adStraindDV_FD,
                        real                            aPerturbation,
                        fem::FDScheme_Type              aFDSchemeType = fem::FDScheme_Type::POINT_5,
                        enum CM_Function_Type           aCMFunctionType = CM_Function_Type::DEFAULT );

                //------------------------------------------------------------------------------
                /*
                 * evaluates and returns the value of E' which is used in the evaluation of stress intensity factors
                 */
                virtual moris::real get_e_prime()
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

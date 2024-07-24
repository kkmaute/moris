/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Constitutive_Model.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_
#define SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_

// MRS/COR/src
#include "moris_typedefs.hpp"
// #include "linalg_typedefs.hpp"
//  MRS/CNT/src
#include "cl_Vector.hpp"
// LNA/src
#include "cl_Matrix.hpp"
// FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Enums.hpp"
#include "fn_FEM_FD_Scheme.hpp"
#include "cl_FEM_Material_Model.hpp"
// FEM/MSI/src
#include "cl_MSI_Dof_Type_Enums.hpp"
#include <map>
// GEN/src
#include "GEN_Data_Types.hpp"

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

          protected:
            // field interpolator manager
            Field_Interpolator_Manager* mFIManager = nullptr;

            // fem set pointer
            Set* mSet = nullptr;

            // list of parameters
            Vector< Matrix< DDRMat > > mParameters;

            // dof type list
            Vector< Vector< MSI::Dof_Type > > mDofTypes;

            // dof type map
            Matrix< DDSMat > mDofTypeMap;

            // global dof type list
            Vector< Vector< MSI::Dof_Type > > mGlobalDofTypes;

            // global dof type map
            Matrix< DDSMat > mGlobalDofTypeMap;

            // dv type list
            Vector< Vector< gen::PDV_Type > > mDvTypes;

            // local string to dv enum map
            std::map< std::string, gen::PDV_Type > mDvMap;

            // global dv type list
            Vector< Vector< gen::PDV_Type > > mGlobalDvTypes;

            // global dv type map
            Matrix< DDSMat > mGlobalDvTypeMap;

            // dv type map
            Matrix< DDSMat > mDvTypeMap;

            // dof type list
            Vector< Vector< mtk::Field_Type > > mFieldTypes;

            // dof type map
            Matrix< DDSMat > mFieldTypeMap;

            // global dof type list
            Vector< Vector< mtk::Field_Type > > mGlobalFieldTypes;

            // global dof type map
            Matrix< DDSMat > mGlobalFieldTypeMap;

            // properties
            Vector< std::shared_ptr< Property > > mProperties;

            // local string to property enum map
            std::map< std::string, uint > mPropertyMap;

            // Material Model
            Vector< std::shared_ptr< Material_Model > > mMaterialModels;

            // local string to property enum map
            std::map< std::string, uint > mMaterialModelMap;

            // spatial dimensions
            uint mSpaceDim = 0;

            // storage for flux evaluation
            Matrix< DDRMat >           mFlux;
            Vector< Matrix< DDRMat > > mdFluxdDof;
            Vector< Matrix< DDRMat > > mdFluxdDv;
            Vector< Matrix< DDRMat > > mdFluxdx;

            // storage for divergence of flux evaluation
            Matrix< DDRMat >           mDivFlux;
            Vector< Matrix< DDRMat > > mddivfluxdu;

            // storage for energy evaluation
            Matrix< DDRMat >           mEnergy;
            Vector< Matrix< DDRMat > > mEnergyDof;

            // storage for gradient of energy evaluation
            Matrix< DDRMat >           mGradEnergy;
            Vector< Matrix< DDRMat > > mGradEnergyDof;

            // storage for energy rate evaluation
            Matrix< DDRMat >           mEnergyDot;
            Vector< Matrix< DDRMat > > mEnergyDotDof;

            // storage for gradient of energy rate evaluation
            Matrix< DDRMat >           mGradEnergyDot;
            Vector< Matrix< DDRMat > > mGradEnergyDotDof;

            // storage for gradient of div flux rate evaluation
            Matrix< DDRMat >           mGradDivFlux;
            Vector< Matrix< DDRMat > > mGradDivFluxDof;

            // storage for traction evaluation
            Matrix< DDRMat >           mTraction;
            Vector< Matrix< DDRMat > > mdTractiondDof;
            Vector< Matrix< DDRMat > > mdTractiondDv;

            // storage for test traction evaluation
            Vector< Matrix< DDRMat > >           mTestTraction;
            Vector< Matrix< DDRMat > >           mTestTractionTrans;
            Vector< Vector< Matrix< DDRMat > > > mdTestTractiondDof;
            Vector< Vector< Matrix< DDRMat > > > mdTestTractiondDv;

            // storage for stress evaluation
            Matrix< DDRMat >           mStress;
            Vector< Matrix< DDRMat > > mdStressdDof;

            // storage for strain evaluation
            Matrix< DDRMat >           mStrain;
            Vector< Matrix< DDRMat > > mdStraindDof;
            Vector< Matrix< DDRMat > > mdStraindDv;
            Vector< Matrix< DDRMat > > mdStraindx;

            // storage for divergence of strain evaluation
            Matrix< DDRMat >           mDivStrain;
            Vector< Matrix< DDRMat > > mddivstraindu;

            // storage for test strain evaluation
            Matrix< DDRMat >           mTestStrain;
            Matrix< DDRMat >           mTestStrainTrans;
            Vector< Matrix< DDRMat > > mdTestStraindDof;

            // storage for constitutive matrix evaluation
            Matrix< DDRMat >           mConst;
            Vector< Matrix< DDRMat > > mdConstdDof;
            Vector< Matrix< DDRMat > > mdConstdDv;

            // constitutive model name for input and debug
            std::string mName = "Undefined";

          private:
            // bool for global dof type list and map build
            bool mGlobalDofBuild      = true;
            bool mGlobalDvBuild       = true;
            bool mGlobalFieldBuild    = true;
            bool mGlobalDofMapBuild   = true;
            bool mGlobalDvMapBuild    = true;
            bool mGlobalFieldMapBuild = true;

            // flag for flux related evaluation
            bool                    mFluxEval = true;
            moris::Matrix< DDBMat > mdFluxdDofEval;
            moris::Matrix< DDBMat > mdFluxdDvEval;
            moris::Matrix< DDBMat > mdFluxdxEval;

            // flag for div flux related evaluation
            bool                    mDivFluxEval = true;
            moris::Matrix< DDBMat > mddivfluxduEval;

            // flag for grad div flux related evaluation
            bool                    mGradDivFluxEval = true;
            moris::Matrix< DDBMat > mGradDivFluxDofEval;

            // flag for traction related evaluation
            bool                    mTractionEval = true;
            moris::Matrix< DDBMat > mdTractiondDofEval;
            moris::Matrix< DDBMat > mdTractiondDvEval;

            // flag for test traction related evaluation
            moris::Matrix< DDBMat > mTestTractionEval;
            moris::Matrix< DDBMat > mTestTractionTransEval;
            moris::Matrix< DDBMat > mdTestTractiondDofEval;
            moris::Matrix< DDBMat > mdTestTractiondDvEval;

            // flag for stress related evaluation
            bool                    mStressEval = true;
            moris::Matrix< DDBMat > mdStressdDofEval;

            // flag for strain related evaluation
            bool                    mStrainEval = true;
            moris::Matrix< DDBMat > mdStraindDofEval;
            moris::Matrix< DDBMat > mdStraindDvEval;
            moris::Matrix< DDBMat > mdStraindxEval;

            // flag for div strain related evaluation
            bool                    mDivStrainEval = true;
            moris::Matrix< DDBMat > mddivstrainduEval;

            // flag for test strain related evaluation
            bool                    mTestStrainEval      = true;
            bool                    mTestStrainTransEval = true;
            moris::Matrix< DDBMat > mdTestStraindDofEval;

            // flag for constitutive matrix related evaluation
            bool                    mConstEval = true;
            moris::Matrix< DDBMat > mdConstdDofEval;
            moris::Matrix< DDBMat > mdConstdDvEval;

            // flag for energy related evaluation
            bool                    mEnergyEval = true;
            moris::Matrix< DDBMat > mEnergyDofEval;

            // flag for energy rate related evaluation
            bool                    mEnergyDotEval = true;
            moris::Matrix< DDBMat > mEnergyDotDofEval;

            // flag for grad energy related evaluation
            bool                    mGradEnergyEval = true;
            moris::Matrix< DDBMat > mGradEnergyDofEval;

            // flag for grad energy rate related evaluation
            bool                    mGradEnergyDotEval = true;
            moris::Matrix< DDBMat > mGradEnergyDotDofEval;

            //------------------------------------------------------------------------------

          public:
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
            virtual Constitutive_Type
            get_constitutive_type() const
            {
                // need to define this for every CM
                return Constitutive_Type::UNDEFINED;
            }

            //------------------------------------------------------------------------------
            /**
             * set parameters
             * @param[ in ] aParameters a list of parameters
             */
            virtual void set_parameters( Vector< Matrix< DDRMat > > aParameters );

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
            std::shared_ptr< fem::Property >& get_property( std::string aPropertyString );

            //------------------------------------------------------------------------------
            /**
             * get the entire property map
             * @param[ out ]  mPropertyMap map of string to uint of properties
             */
            std::map< std::string, uint >&
            get_property_map()
            {
                return mPropertyMap;
            }

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
            std::shared_ptr< fem::Material_Model >& get_material_model( std::string aMaterialModelString );

            //------------------------------------------------------------------------------
            /**
             * set local properties
             */
            virtual void
            set_local_properties()
            {
            }

            //------------------------------------------------------------------------------
            /**
             * set local material model
             */
            virtual void
            set_local_material_model()
            {
            }

            //------------------------------------------------------------------------------
            /**
             * set name
             * param[ in ] aName a string for CM name
             */
            void
            set_name( std::string aName )
            {
                mName = aName;
            }

            //------------------------------------------------------------------------------
            /**
             * get name
             * param[ out ] mName a string for CM name
             */
            std::string
            get_name()
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
            void
            set_set_pointer( Set* aSetPointer )
            {
                mSet = aSetPointer;
            }

            //------------------------------------------------------------------------------
            /**
             * set space dimension
             * @param[ in ] aSpaceDim a spatial dimension
             */
            virtual void
            set_space_dim( uint aSpaceDim )
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
            uint
            get_num_space_dims()
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
            virtual void reset_eval_flags();

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
            void set_dof_type_list( Vector< Vector< MSI::Dof_Type > > aDofTypes );

            //------------------------------------------------------------------------------
            /**
             * set constitutive model dof types
             * @param[ in ] aDofTypes a list of group of dof types
             * @param[ in ] aDofStrings a list of strings to describe the dof types
             */
            virtual void
            set_dof_type_list(
                    Vector< Vector< MSI::Dof_Type > > aDofTypes,
                    Vector< std::string >             aDofStrings )
            {
                MORIS_ERROR( false, "Constitutive_Model::set_dof_type_list - Not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * return a cell of dof types
             * @param[ in ] mDofTypes a cell of cell of dof types
             */
            const Vector< Vector< MSI::Dof_Type > >&
            get_dof_type_list() const
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
            const Matrix< DDSMat >&
            get_dof_type_map()
            {
                return mDofTypeMap;
            }

            //------------------------------------------------------------------------------
            /**
             * set constitutive model dv types
             * @param[ in ] aDvTypes a list of group of dv types
             */
            void set_dv_type_list( Vector< Vector< gen::PDV_Type > > aDvTypes );

            //------------------------------------------------------------------------------
            /**
             * set constitutive model dv types
             * @param[ in ] aDvTypes   a list of group of dv types
             * @param[ in ] aDvStrings a list of strings to describe the dv types
             */
            virtual void
            set_dv_type_list(
                    Vector< Vector< gen::PDV_Type > > aDvTypes,
                    Vector< std::string >             aDvStrings )
            {
                MORIS_ERROR( false, "Constitutive_Model::set_model_type - Not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * return a cell of dv types
             * @param[ out ] aDvTypes a cell of cell of dv types
             */
            const Vector< Vector< gen::PDV_Type > >&
            get_dv_type_list() const
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
            const Matrix< DDSMat >&
            get_dv_type_map()
            {
                return mDvTypeMap;
            }

            //------------------------------------------------------------------------------
            /**
             * set constitutive model field types
             * @param[ in ] aFieldTypes a list of group of field types
             */
            void set_field_type_list( Vector< Vector< mtk::Field_Type > > aFieldTypes );

            //------------------------------------------------------------------------------
            /**
             * return a cell of field types
             * @param[ in ] mFieldTypes a cell of cell of field types
             */
            const Vector< Vector< mtk::Field_Type > >&
            get_field_type_list() const
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
            const Matrix< DDSMat >&
            get_field_type_map()
            {
                return mFieldTypeMap;
            }

            //------------------------------------------------------------------------------
            /**
             * set field interpolator manager
             * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
             */
            void set_field_interpolator_manager( Field_Interpolator_Manager* aFieldInterpolatorManager );

            //------------------------------------------------------------------------------
            /**
             * set model type
             * @param[ in ] aModelType an enum for model type
             */
            virtual void
            set_model_type( fem::Model_Type aModelType )
            {
                MORIS_ERROR( false, "Constitutive_Model::set_model_type - Not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * @return plane type
             */
            virtual Model_Type
            get_plane_type() const
            {
                return Model_Type::UNDEFINED;
            }

            //------------------------------------------------------------------------------
            /**
             * get properties
             * @param[ out ] mProperties cell of property pointers
             */
            Vector< std::shared_ptr< Property > >&
            get_properties()
            {
                return mProperties;
            };

            //------------------------------------------------------------------------------
            /**
             * get material model
             * @param[ out ] mMaterialModel cell of property pointers
             */
            Vector< std::shared_ptr< Material_Model > >&
            get_material_models()
            {
                return mMaterialModels;
            };

            //------------------------------------------------------------------------------
            /**
             * create a global dof type list including constitutive and property dependencies
             */
            virtual void build_global_dof_type_list();

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
            const Vector< Vector< MSI::Dof_Type > >& get_global_dof_type_list();

            //------------------------------------------------------------------------------
            /**
             * get non unique list of dof type for the constitutive model
             * @param[ in ] aDofTypes a cell of dof type to fill
             */
            void get_non_unique_dof_types( Vector< MSI::Dof_Type >& aDofTypes );

            //------------------------------------------------------------------------------
            /**
             * get non unique lists of dof and dv type for the constitutive model
             * @param[ in ] aDofTypes   a cell of dof type to fill
             * @param[ in ] aDvTypes    a cell of dv type to fill
             * @param[ in ] aFieldTypes a cell of field type to fill
             */
            void get_non_unique_dof_dv_and_field_types(
                    Vector< MSI::Dof_Type >&   aDofTypes,
                    Vector< gen::PDV_Type >&   aDvTypes,
                    Vector< mtk::Field_Type >& aFieldTypes );

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
            const Matrix< DDSMat >& get_global_dof_type_map();

            //------------------------------------------------------------------------------
            /**
             * check dependency on a given group of dof types
             * @param[ in ]  aDofType       a group of dof types
             * @param[ out ] tDofDependency a bool true if dependency on dof type
             *
             */
            bool check_dof_dependency( const Vector< MSI::Dof_Type >& aDofType );

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
            const Vector< Vector< gen::PDV_Type > >& get_global_dv_type_list();

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
            const Matrix< DDSMat >&
            get_global_dv_type_map()
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
            bool check_dv_dependency( const Vector< gen::PDV_Type >& aDvType );

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
            const Vector< Vector< mtk::Field_Type > >& get_global_field_type_list();

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
            const Matrix< DDSMat >& get_global_field_type_map();

            //------------------------------------------------------------------------------
            /**
             * check dependency on a given group of field types
             * @param[ in ]  aFieldType       a group of field types
             * @param[ out ] tFieldDependency a bool true if dependency on field type
             *
             */
            bool check_field_dependency( const Vector< mtk::Field_Type >& aFieldType );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux
             */
            virtual void
            eval_flux()
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_flux - This function does nothing. " );
            }

            /**
             * get the constitutive model flux
             * @param[ in ]  aCMFunctionType  enum indicating which flux is called, if there are several
             * @param[ out ] mFlux constitutive model flux
             */
            virtual const Matrix< DDRMat >&
            flux( enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            virtual const Matrix< DDRMat >&
            flux( int                     aFlatType,
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            // ---------------------------------------------------------------------------------------------------------------------------------

            //------------------------------------------------------------------------------
            /**
             * evaluate the divergence of the flux
             */
            virtual void
            eval_divflux()
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_divflux - This function does nothing. " );
            }

            /**
             * get the divergence of the flux
             * @param[ out ] mDivFlux divergence of the flux
             */
            virtual const Matrix< DDRMat >&
            divflux( enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the divergence of the flux wrt to dof type
             */
            virtual void
            eval_ddivfluxdu( const Vector< MSI::Dof_Type >& aDofType )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_ddivfluxdu - This function does nothing. " );
            }

            /**
             * get the derivative of the divergence of the flux wrt to dof type
             * @param[ out ] mddivfluxdu derivative of the divergence of the flux
             *                           wrt to dof type
             */
            virtual const Matrix< DDRMat >& ddivfluxdu(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model traction
             * @param[ in ]  aNormal normal
             */
            virtual void
            eval_traction( const Matrix< DDRMat >& aNormal )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_traction - This function does nothing. " );
            }

            /**
             * get the constitutive model traction
             * @param[ in ]  aNormal   normal
             * @param[ out ] mTraction constitutive model traction
             */
            virtual const Matrix< DDRMat >& traction(
                    const Matrix< DDRMat >& aNormal,
                    enum CM_Function_Type   aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test traction
             * @param[ in ]  aNormal normal
             */
            virtual void
            eval_testTraction(
                    const Matrix< DDRMat >&        aNormal,
                    const Vector< MSI::Dof_Type >& aTestDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_testTraction - This function does nothing. " );
            }

            /**
             * get the constitutive model test traction
             * @param[ in ]  aNormal       normal
             * @param[ out ] mTestTraction test traction
             *                             size nSpaceDim x nTestDof
             */
            virtual const Matrix< DDRMat >& testTraction(
                    const Matrix< DDRMat >&        aNormal,
                    const Vector< MSI::Dof_Type >& aTestDofTypes,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            /**
             * get the transpose of the constitutive model test traction
             * @param[ in ]  aNormal            normal
             * @param[ out ] mTestTractionTrans transpose of test traction
             *                                  size nSpaceDim x nTestDof
             */
            virtual const Matrix< DDRMat >& testTraction_trans(
                    const Matrix< DDRMat >&        aNormal,
                    const Vector< MSI::Dof_Type >& aTestDofTypes,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model stress
             */
            virtual void
            eval_stress()
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_stress - This function does nothing. " );
            }

            /**
             * get the constitutive model stress
             * @param[ out ] mStress constitutive model stress
             */
            virtual const Matrix< DDRMat >& stress(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain
             */
            virtual void
            eval_strain()
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_strain - This function does nothing. " );
            }

            /**
             * get the constitutive model strain
             * @param[ out ] mStrain constitutive model strain
             */
            virtual const Matrix< DDRMat >& strain(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the divergence of the strain
             */
            virtual void
            eval_divstrain()
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_divstrain - This function does nothing. " );
            }

            /**
             * get the divergence of the strain
             * @param[ out ] mDivFlux divergence of the strain
             */
            virtual const Matrix< DDRMat >& divstrain(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the divergence of the strain wrt to dof type
             */
            virtual void
            eval_ddivstraindu( const Vector< MSI::Dof_Type >& aDofType )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_ddivstraindu - This function does nothing. " );
            }

            /**
             * get the derivative of the divergence of the strain wrt to dof type
             * @param[ out ] mddivstraindu derivative of the divergence of the strain
             *                             wrt to dof type
             */
            virtual const Matrix< DDRMat >& ddivstraindu(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test strain
             */
            virtual void
            eval_testStrain()
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_testStrain - This function does nothing. " );
            }

            /**
             * get the constitutive model test strain
             * @param[ out ] mTestStrain constitutive model test strain
             */
            virtual const Matrix< DDRMat >& testStrain(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            /**
             * get the transpose of the constitutive model test strain
             * @param[ out ] mTestStrain transpose of constitutive model test strain
             */
            virtual const Matrix< DDRMat >& testStrain_trans(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            /**
             * evaluate the derivative of the constitutive model test strain wrt dof
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            virtual void
            eval_dTestStraindDOF( const Vector< MSI::Dof_Type >& aDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dTestStraindDOF - This function does nothing. " );
            }

            /**
             * get the derivative of the constitutive model test strain wrt dof
             * @param[ in ]  aDofTypes a dof type wrt which the derivative is evaluated
             * @param[ out ] mdTestStraindDOF derivative of constitutive model test strain
             */
            virtual const Matrix< DDRMat >& dTestStraindDOF(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model constitutive matrix
             */
            virtual void
            eval_const()
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_const - This function does nothing. " );
            }

            /**
             * get the constitutive model constitutive matrix
             * @param[ out ] mConst constitutive matrix
             */
            virtual const Matrix< DDRMat >& constitutive(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the flux wrt space
             * @param[ in ] aOrder order of the derivative
             */
            virtual void
            eval_dfluxdx( uint aOrder )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dfluxdx - This function does nothing. " );
            }

            /**
             * get the derivative of the flux wrt space
             * @param[ in ] aOrder order of the derivative
             */
            virtual const Matrix< DDRMat >& dfluxdx(
                    uint                  aOrder,
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux derivative wrt to a dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             */
            virtual void
            eval_dFluxdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dFluxdDOF - This function does nothing. " );
            }

            /**
             * get the derivative of the flux wrt dof
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ out ] mFluxDofDer derivative of the flux wrt dof
             */
            virtual const Matrix< DDRMat >& dFluxdDOF(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model energy
             */
            virtual void
            eval_Energy()
            {
                MORIS_ASSERT( false, "eval_Energy: not implemented in base class." );
            };

            /**
             * get the constitutive model change rate of energy
             * @param[ out ] mEnergyDot change rate of energy
             */
            virtual const Matrix< DDRMat >& Energy( enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model change rate of energy
             */
            virtual void
            eval_EnergyDot()
            {
                MORIS_ASSERT( false, "eval_EnergyDot: not implemented in base class." );
            };

            /**
             * get the constitutive model change rate of enthalpy
             * @param[ out ] mEnergyDot change rate of enthalpy
             */
            virtual const Matrix< DDRMat >& EnergyDot(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model spatial gradient of enthalpy
             */
            virtual void
            eval_gradEnergy()
            {
                MORIS_ASSERT( false, "eval_gradEnergy: not implemented in base class." );
            };

            /**
             * get the constitutive model spatial gradient of enthalpy
             * @param[ out ] mGradEnergy gradient of enthalpy
             */
            virtual const Matrix< DDRMat >& gradEnergy(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluates the constitutive model change rate of spatial gradient of enthalpy (needed for GGLS-stabilization)
             */
            virtual void
            eval_gradEnergyDot()
            {
                MORIS_ASSERT( false, "eval_gradEnergyDot: not implemented in base class." );
            };

            /**
             * get the constitutive model change rate of spatial gradient of enthalpy (needed for GGLS-stabilization)
             * @param[ out ] mGradEnergyDot gradient of change rate of enthalpy
             */
            virtual const Matrix< DDRMat >& gradEnergyDot(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluates the gradient of the divergence of the flux (needed for GGLS-stabilization)
             */
            virtual void
            eval_graddivflux()
            {
                MORIS_ASSERT( false, "eval_graddivflux: not implemented in base class." );
            };

            /**
             * get the gradient of the divergence of the flux (needed for GGLS-stabilization)
             * @param[ out ] mGradDivFlux gradient of divergence of flux
             */
            virtual const Matrix< DDRMat >& graddivflux(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model energy change rate wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dEnergyDotdDOF ( 1 x numDerDof )
             */
            virtual void
            eval_dEnergydDOF( const Vector< MSI::Dof_Type >& aDofTypes )
            {
                MORIS_ASSERT( false, "eval_dEnergydDOF: not implemented in base class." );
            };

            /**
             * get the energy change rate wrt dof
             * @param[ in ]  aDofType        group of dof type
             * @param[ out ] mEnergyDotDofDer derivative of the traction wrt dof
             */
            virtual const Matrix< DDRMat >& dEnergydDOF(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model enthalpy change rate wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dEnergyDotdDOF ( 1 x numDerDof )
             */
            virtual void
            eval_dEnergyDotdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
            {
                MORIS_ASSERT( false, "eval_dEnergyDotdDOF: not implemented in base class." );
            };

            /**
             * get the enthalpy change rate wrt dof
             * @param[ in ]  aDofType        group of dof type
             * @param[ out ] mEnergyDotDofDer derivative of the traction wrt dof
             */
            virtual const Matrix< DDRMat >& dEnergyDotdDOF(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model gradient of enthalpy wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dgradEnergydDOF ( mSpaceDim x numDerDof )
             */
            virtual void
            eval_dGradEnergydDOF( const Vector< MSI::Dof_Type >& aDofTypes )
            {
                MORIS_ASSERT( false, "eval_dGradEnergydDOF: not implemented in base class." );
            };

            /**
             * get the gradient of enthalpy wrt dof
             * @param[ in ]  aDofType        group of dof type
             * @param[ out ] mGradEnergyDer derivative of the gradient of enthalpy wrt dof
             */
            virtual const Matrix< DDRMat >& dGradEnergydDOF(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model gradient of enthalpy change rate wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dgradEnergyDotdDOF ( mSpaceDim x numDerDof )
             */
            virtual void
            eval_dGradEnergyDotdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
            {
                MORIS_ASSERT( false, "eval_dGradEnergyDotdDOF: not implemented in base class." );
            };

            /**
             * get the gradient of enthalpy change rate wrt dof
             * @param[ in ]  aDofType        group of dof type
             * @param[ out ] mGradEnergyDotDer derivative of the traction wrt dof
             */
            virtual const Matrix< DDRMat >& dGradEnergyDotdDOF(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model gradient of divergence of flux wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * dGradDivFluxdDOF ( mSpaceDim x numDerDof )
             */
            virtual void
            eval_dGradDivFluxdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
            {
                MORIS_ASSERT( false, "eval_dGradDivFluxdDOF: not implemented in base class." );
            };

            /**
             * get the gradient of enthalpy change rate wrt dof
             * @param[ in ]  aDofType        group of dof type
             * @param[ out ] mGradEnergyDotDer derivative of the traction wrt dof
             */
            virtual const Matrix< DDRMat >& dGradDivFluxdDOF(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model traction derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             * @param[ in ] aNormal   normal
             */
            virtual void
            eval_dTractiondDOF(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    const Matrix< DDRMat >&        aNormal )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dTractiondDOF - This function does nothing. " );
            }

            /**
             * get the derivative of the traction wrt dof
             * @param[ in ]  aDofType        group of dof type
             * @param[ in ]  aNormal         normal
             * @param[ out ] mTractionDofDer derivative of the traction wrt dof
             */
            virtual const Matrix< DDRMat >& dTractiondDOF(
                    const Vector< MSI::Dof_Type >& aDofType,
                    const Matrix< DDRMat >&        aNormal,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test traction derivative wrt to a dof type
             * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ in ] adTractiondDOF a matrix to fill with derivative evaluation
             */
            virtual void
            eval_dTestTractiondDOF(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    const Matrix< DDRMat >&        aNormal,
                    const Vector< MSI::Dof_Type >& aTestDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dTestTractiondDOF - This function does nothing. " );
            }

            /**
             * get the derivative of the test traction wrt dof
             * @param[ in ]  aDofType           group of dof type
             * @param[ in ]  aNormal            normal
             * @param[ out ] mdTestTractiondDof derivative of the traction wrt dof
             */
            virtual const Matrix< DDRMat >& dTestTractiondDOF(
                    const Vector< MSI::Dof_Type >& aDofType,
                    const Matrix< DDRMat >&        aNormal,
                    const Vector< MSI::Dof_Type >& aTestDofTypes,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * FIXME this is not a test traction, used for elast lin iso!!!!
             * evaluate the constitutive model test traction derivative wrt to a dof type
             * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ in ] adTractiondDOF a matrix to fill with derivative evaluation
             */
            virtual void
            eval_dTestTractiondDOF(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    const Matrix< DDRMat >&        aNormal,
                    const Matrix< DDRMat >&        aJump,
                    const Vector< MSI::Dof_Type >& aTestDofTypes )
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
            virtual const Matrix< DDRMat >& dTestTractiondDOF(
                    const Vector< MSI::Dof_Type >& aDofType,
                    const Matrix< DDRMat >&        aNormal,
                    const Matrix< DDRMat >&        aJump,
                    const Vector< MSI::Dof_Type >& aTestDofTypes,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the strain wrt space
             * @param[ in ] aOrder order of the derivative
             */
            virtual void
            eval_dstraindx( uint aOrder )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dstraindx - This function does nothing. " );
            }

            /**
             * get the derivative of the strain wrt space
             * @param[ in ] aOrder order of the derivative
             */
            virtual const Matrix< DDRMat >& dstraindx(
                    uint                  aOrder,
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model stress derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            virtual void
            eval_dStressdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dStressdDOF - This function does nothing. " );
            }

            /**
             * get the derivative of the stress wrt dof
             * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ out ] mdStressdDof derivative of the stress wrt dof
             */
            virtual const Matrix< DDRMat >& dStressdDOF(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            virtual void
            eval_dStraindDOF( const Vector< MSI::Dof_Type >& aDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dStraindDOF - This function does nothing. " );
            }

            /**
             * get the derivative of the strain wrt dof
             * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ out ] mdStraindDof derivative of the strain wrt dof
             */
            virtual const Matrix< DDRMat >& dStraindDOF(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model constitutive matrix derivative wrt to a dof type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            virtual void
            eval_dConstdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dConstdDOF - This function does nothing. " );
            }

            /**
             * get the derivative of the constitutive matrix wrt dof
             * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
             * @param[ out ] mdConstdDof derivative of the constitutive matrix wrt dof
             */
            virtual const Matrix< DDRMat >& dConstdDOF(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux derivative wrt to a dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             */
            virtual void
            eval_dFluxdDV( const Vector< gen::PDV_Type >& aDvTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dFluxdDV - This function does nothing. " );
            }

            /**
             * get the derivative of the flux wrt dv
             * @param[ in ]  aDvTypes  a dv type wrt which the derivative is evaluated
             * @param[ out ] mdFluxdDv derivative of the flux wrt dv
             */
            virtual const Matrix< DDRMat >& dFluxdDV(
                    const Vector< gen::PDV_Type >& aDvType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain derivative wrt to a dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             */
            virtual void
            eval_dStraindDV( const Vector< gen::PDV_Type >& aDvTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dStraindDV - This function does nothing. " );
            }

            /**
             * get the derivative of the strain wrt dv
             * @param[ in ]  aDvTypes    a dv type wrt which the derivative is evaluated
             * @param[ out ] mdStraindDv derivative of the strain wrt dv
             */
            virtual const Matrix< DDRMat >& dStraindDV(
                    const Vector< gen::PDV_Type >& aDvType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model constitutive matrix derivative wrt to a dv type
             * @param[ in ] aDvTypes a dof type wrt which the derivative is evaluated
             */
            virtual void
            eval_dConstdDV( const Vector< gen::PDV_Type >& aDvTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dConstdDV - This function does nothing. " );
            }

            /**
             * get the derivative of the constitutive matrix wrt dv
             * @param[ in ]  aDvTypes   a dv type wrt which the derivative is evaluated
             * @param[ out ] mdConstdDv derivative of the constitutive matrix wrt dv
             */
            const Matrix< DDRMat >& dConstdDV(
                    const Vector< gen::PDV_Type >& aDvType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate derivative wrt to a dof type
             * @param[ in ] aDerivativeFD   a matrix to fill with derivative evaluation
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] aTestDofTypes   a test dof type wrt which the test traction is evaluated
             * @param[ in ] aNormal         a normal
             * @param[ in ] aJump           a jump
             * @param[ in ] aPerturbation   a real to perturb for FD
             * @param[ in ] aFDSchemeType   enum for FD scheme
             * @param[ in ] aCMFunctionType
             */
            void eval_derivative_FD(
                    enum CM_Request_Type           aCMRequestType,
                    Matrix< DDRMat >&              aDerivativeFD,
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    real                           aPerturbation,
                    const Vector< MSI::Dof_Type >& aTestDofTypes,
                    const Matrix< DDRMat >&        aNormal,
                    const Matrix< DDRMat >&        aJump,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            /**
             * select derivative wrt to a dof type
             * @param[ in ] aCMRequestType  a type for required derivative
             * @param[ in ] aTestDofTypes   a test dof type wrt which the test traction is evaluated
             * @param[ in ] aNormal         a normal
             * @param[ in ] aJump           a jump
             * @param[ in ] aCMFunctionType
             * Rem: implement specific version for child CM
             */
            virtual const Matrix< DDRMat >& select_derivative_FD(
                    enum CM_Request_Type           aCMRequestType,
                    const Vector< MSI::Dof_Type >& aTestDofTypes,
                    const Matrix< DDRMat >&        aNormal,
                    const Matrix< DDRMat >&        aJump,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            /**
             * select derivative wrt to a dof type
             * @param[ in ] aCMRequestType  a type for required derivative
             * @param[ in ] aDerivativeFD   a derivative value to set to storage
             * @param[ in ] aTestDofTypes   a test dof type wrt which the test traction is evaluated
             * @param[ in ] aCMFunctionType
             * Rem: implement specific version for child CM
             */
            virtual void set_derivative_FD(
                    enum CM_Request_Type           aCMRequestType,
                    Matrix< DDRMat >&              aDerivativeFD,
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    const Vector< MSI::Dof_Type >& aTestDofTypes,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model stress derivative wrt to a dof type
             * @param[ in ] aDofTypes       dof type wrt which the derivative is evaluated
             * @param[ in ] adFluxdDOF_FD   matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation   real to perturb for FD
             * @param[ in ] aFDSchemeType   enum for FD scheme
             */
            virtual void eval_dFluxdDOF_FD(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    Matrix< DDRMat >&              adFluxdDOF_FD,
                    real                           aPerturbation,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the traction derivative wrt to a dof type
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] adtractiondu_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation   a real to perturb for FD
             * @param[ in ] aNormal         a normal
             */
            void eval_dtractiondu_FD(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    Matrix< DDRMat >&              adtractiondu_FD,
                    real                           aPerturbation,
                    Matrix< DDRMat >&              aNormal,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

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
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    const Vector< MSI::Dof_Type >& aTestDofTypes,
                    Matrix< DDRMat >&              adtesttractiondu_FD,
                    real                           aPerturbation,
                    const Matrix< DDRMat >&        aNormal,
                    const Matrix< DDRMat >&        aJump,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the test traction derivative wrt to a dof type
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] aTestDofTypes   a test dof type wrt which the test traction is evaluated
             * @param[ in ] adtestractiondu_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation   a real to perturb for FD
             * @param[ in ] aNormal         a normal
             */
            void eval_dtesttractiondu_FD(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    const Vector< MSI::Dof_Type >& aTestDofTypes,
                    Matrix< DDRMat >&              adtesttractiondu_FD,
                    real                           aPerturbation,
                    const Matrix< DDRMat >&        aNormal,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the div flux derivative wrt to a dof type
             * @param[ in ] aDofTypes     a dof type wrt which the derivative is evaluated
             * @param[ in ] ddivfluxdu_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation a real to perturb for FD
             */
            void eval_ddivfluxdu_FD(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    Matrix< DDRMat >&              ddivfluxdu_FD,
                    real                           aPerturbation,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the div strain derivative wrt to a dof type
             * @param[ in ] aDofTypes        a dof type wrt which the derivative is evaluated
             * @param[ in ] addivstraindu_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation    a real to perturb for FD
             */
            void eval_ddivstraindu_FD(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    Matrix< DDRMat >&              addivstraindu_FD,
                    real                           aPerturbation,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the internal energy / enthalpy wrt dof using finite differences
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] adEnergydDOF_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation   real to perturb for FD
             */
            void eval_dEnergydDOF_FD(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    Matrix< DDRMat >&              adEnergydDOF_FD,
                    real                           aPerturbation,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the internal energy / enthalpy change rate wrt dof using finite differences
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] adEnergyDotdDOF_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation   real to perturb for FD
             */
            void eval_dEnergyDotdDOF_FD(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    Matrix< DDRMat >&              adEnergyDotdDOF_FD,
                    real                           aPerturbation,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the gradient of enthalpy wrt dof using finite differences
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] adGradEnergydDOF_FD  a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation   real to perturb for FD
             */
            void eval_dGradEnergydDOF_FD(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    Matrix< DDRMat >&              adGradEnergydDOF_FD,
                    real                           aPerturbation,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the gradient of enthalpy change rate wrt dof using finite differences
             * @param[ in ] aDofTypes          a dof type wrt which the derivative is evaluated
             * @param[ in ] adGradEnergyDotdDOF_FD  a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation      real to perturb for FD
             */
            void eval_dGradEnergyDotdDOF_FD(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    Matrix< DDRMat >&              adGradEnergyDotdDOF_FD,
                    real                           aPerturbation,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the gradient of enthalpy change rate wrt dof using finite differences
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] adFluxdDOF_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation   real to perturb for FD
             */
            void eval_dGradDivFluxdDOF_FD(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    Matrix< DDRMat >&              adGradDivFluxdDOF_FD,
                    real                           aPerturbation,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain derivative wrt to a dof type
             * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
             * @param[ in ] adStraindDOF_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation   real to perturb for FD
             */
            void eval_dStraindDOF_FD(
                    const Vector< MSI::Dof_Type >& aDofTypes,
                    Matrix< DDRMat >&              adStraindDOF_FD,
                    real                           aPerturbation,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model stress derivative wrt to a dv type
             * @param[ in ] aDofTypes     a dof type wrt which the derivative is evaluated
             * @param[ in ] adFluxdDV_FD  a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation real to perturb for FD
             */

            void eval_dFluxdDV_FD(
                    const Vector< gen::PDV_Type >& aDvTypes,
                    Matrix< DDRMat >&              adFluxdDV_FD,
                    real                           aPerturbation,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain derivative wrt to a dv type
             * @param[ in ] aDvTypes       a dv type wrt which the derivative is evaluated
             * @param[ in ] adStraindDV_FD a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation  real to perturb for FD
             */

            void eval_dStraindDV_FD(
                    const Vector< gen::PDV_Type >& aDvTypes,
                    Matrix< DDRMat >&              adStraindDV_FD,
                    real                           aPerturbation,
                    fem::FDScheme_Type             aFDSchemeType   = fem::FDScheme_Type::POINT_5,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT );

            //------------------------------------------------------------------------------
            /*
             * evaluates and returns the value of E' which is used in the evaluation of stress intensity factors
             */
            virtual moris::real
            get_e_prime()
            {
                MORIS_ERROR( false, " Constitutive_Model::get_e_prime - This function does nothing. " );
                return 0;
            }

            //--------------------------------------------------------------------------------------------------------------
            // FIXME to be removed
            /**
             * get the turbulent dynamic viscosity mu_t = rho * vtilde * tf1
             * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
             *               if there are several
             * @param[ out ] mTurbDynVisc effective conductivity
             */
            virtual const Matrix< DDRMat >&
            turbulent_dynamic_viscosity(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::turbulent_dynamic_viscosity - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the effective dynamic viscosity mu_eff = mu + mu_t
             * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
             *               if there are several
             * @param[ out ] mEffDynVisc effective conductivity
             */
            virtual const Matrix< DDRMat >&
            effective_dynamic_viscosity(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::effective_dynamic_viscosity - This function does nothing. " );
                return mFlux;
            }

            //--------------------------------------------------------------------------------------------------------------
            /**
             * get the effective conductivity k_eff = k + k_t
             * @param[ in ]  aCMFunctionType  enum indicating which effective conductivity is called,
             *               if there are several
             * @param[ out ] mEffCond effective conductivity
             */
            virtual const Matrix< DDRMat >&
            effective_conductivity(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::effective_conductivity - This function does nothing. " );
                return mFlux;
            }

            //--------------------------------------------------------------------------------------------------------------
            // FIXME to be removed
            /**
             * get the the production coefficient
             * @param[ in ]  aCMFunctionType enum for specific production term if several
             * @param[ out ] mProductionTerm production term
             */
            virtual const Matrix< DDRMat >&
            production_term(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::production_term - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the derivative of the production term wrt dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ in ] aCMFunctionType enum for specific type of which production term
             * @param[ out ] mdProductionTermdu derivative of the production term wrt dof types
             */
            virtual const Matrix< DDRMat >&
            dproductiontermdu(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::dproductiontermdu - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the the production coefficient
             * @param[ in ]  aCMFunctionType enum for specific production coefficient if several
             * @param[ out ] mProductionTerm production term
             */
            virtual const Matrix< DDRMat >&
            production_coefficient(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::production_coefficient - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the derivative of the production coefficient wrt dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ in ] aCMFunctionType enum for specific type of which production coefficient
             * @param[ out ] mdProductionCoeffdu derivative of the production coefficient wrt dof types
             */
            virtual const Matrix< DDRMat >&
            dproductioncoeffdu(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::dproductioncoeffdu - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the the wall destruction term
             * @param[ in ]  aCMFunctionType enum for specific wall destruction term if several
             * @param[ out ] mWallDestructionTerm wall destruction term
             */
            virtual const Matrix< DDRMat >&
            wall_destruction_term(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::wall_destruction_term - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the derivative of the wall destruction term wrt dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ in ] aCMFunctionType enum for specific type of which wall destruction term
             * @param[ out ] mdwalldestructiontermdu derivative of the wall destruction term wrt dof types
             */
            virtual const Matrix< DDRMat >&
            dwalldestructiontermdu(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::dwalldestructiontermdu - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the the wall destruction coefficient
             * @param[ in ]  aCMFunctionType enum for specific wall destruction coefficient if several
             * @param[ out ] mWallDestructionCoeff wall destruction coefficient
             */
            virtual const Matrix< DDRMat >&
            wall_destruction_coefficient(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::wall_destruction_coefficient - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the derivative of the wall destruction coefficient wrt dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ in ] aCMFunctionType enum for specific type of which wall destruction coefficient
             * @param[ out ] mdwalldestructioncoeffdu derivative of the wall destruction coefficient wrt dof types
             */
            virtual const Matrix< DDRMat >&
            dwalldestructioncoeffdu(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::dwalldestructioncoeffdu - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the diffusion coefficient
             * @param[ in ]  aCMFunctionType enum for specific diffusion coefficient if several
             * @param[ out ] mDiffusionCoeff diffusion coefficient
             */
            virtual const Matrix< DDRMat >&
            diffusion_coefficient(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::diffusion_coefficient - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the derivative of the diffusion coefficient wrt dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ in ] aCMFunctionType enum for specific type of which diffusion
             * @param[ out ] mdproductioncoeffdu derivative of the diffusion wrt dof types
             */
            virtual const Matrix< DDRMat >&
            ddiffusioncoeffdu(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::ddiffusioncoeffdu - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the modified velocity u_tilde = u - ( cb2 / sigma ) * dnu_tilde/dx
             * @param[ in ]  aCMFunctionType enum for modified velocity if several
             * @param[ out ] mModVelocity modified velocity
             */
            virtual const Matrix< DDRMat >&
            modified_velocity(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::modified_velocity - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the derivative of the modified velocity wrt dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ in ] aCMFunctionType enum for specific type of which diffusion
             * @param[ out ] mdmodvelocitydu derivative of the modified velocity wrt dof types
             */
            virtual const Matrix< DDRMat >&
            dmodvelocitydu(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::dmodvelocitydu - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the linearized version of the modified velocity u_tilde = u - ( cb2 / sigma ) * dnu_tilde/dx
             * @param[ in ]  aCMFunctionType enum for modified velocity if several
             * @param[ out ] mModVelocity modified velocity
             */
            virtual const Matrix< DDRMat >&
            modified_velocity_linearized(
                    enum CM_Function_Type aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::modified_velocity - This function does nothing. " );
                return mFlux;
            }

            /**
             * get the derivative of the linearized version of the modified velocity wrt dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ in ] aCMFunctionType enum for specific type of which diffusion
             * @param[ out ] mdmodvelocitydu derivative of the modified velocity wrt dof types
             */
            virtual const Matrix< DDRMat >&
            dmodvelocitylinearizeddu(
                    const Vector< MSI::Dof_Type >& aDofType,
                    enum CM_Function_Type          aCMFunctionType = CM_Function_Type::DEFAULT )
            {
                MORIS_ERROR( false, " Constitutive_Model::dmodvelocitydu - This function does nothing. " );
                return mFlux;
            }

            //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_ */

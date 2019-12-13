/*
 * cl_FEM_Constitutive_Model.hpp
 *
 *  Created on: Sep 17, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_
#define SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "linalg_typedefs.hpp"              //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Enums.hpp"                 //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src

namespace moris
{
    namespace fem
    {
    class Field_Interpolator_Manager;
//------------------------------------------------------------------------------
        /**
         * Constitutive model
         */
        class Constitutive_Model
        {

        private:

            bool mFluxEval = true;

        protected :

            // constitutive model type
            fem::Constitutive_Type mConstitutiveType;

            Field_Interpolator_Manager * mFieldInterpolatorManager = nullptr;

            Set * mSet = nullptr;

            // dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mDofTypes;

            // dof type map
            Matrix< DDSMat > mDofTypeMap;

            // bool for global dof type list and map build
            bool mGlobalDofBuild = true;
            bool mGlobalDvBuild  = true;

            // global dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mGlobalDofTypes;

            // global dof type map
            Matrix< DDSMat > mGlobalDofTypeMap;

            // dof field interpolators
            moris::Cell< Field_Interpolator* > mDofFI;

            // dv type list
            moris::Cell< moris::Cell< MSI::Dv_Type > > mDvTypes;

            // global dv type list
            moris::Cell< moris::Cell< MSI::Dv_Type > > mGlobalDvTypes;

            // global dv type map
            Matrix< DDSMat > mGlobalDvTypeMap;

            // dv type map
            Matrix< DDSMat > mDvTypeMap;

            // dv field interpolators
            moris::Cell< Field_Interpolator* > mDvFI;

            // properties
            moris::Cell< std::shared_ptr< Property > > mProperties;

            // spatial dimensions
            uint mSpaceDim;

            // flag for evaluation
            moris::Cell< bool > mdFluxdDofEval;
            moris::Cell< bool > mdFluxdDvEval;

            bool mTractionEval   = true;
            moris::Cell< bool > mdTractiondDofEval;
            moris::Cell< bool > mdTractiondDvEval;

            bool mTestTractionEval   = true;
            moris::Cell< bool > mdTestTractiondDofEval;
            moris::Cell< bool > mdTestTractiondDvEval;

            bool mStrainEval     = true;
            moris::Cell< bool > mdStraindDofEval;
            moris::Cell< bool > mdStraindDvEval;

            bool mTestStrainEval = true;

            bool mConstEval      = true;
            moris::Cell< bool > mdConstdDofEval;
            moris::Cell< bool > mdConstdDvEval;

            // storage for evaluation
            Matrix< DDRMat > mFlux;
            moris::Cell< Matrix< DDRMat > > mdFluxdDof;
            moris::Cell< Matrix< DDRMat > > mdFluxdDv;

            Matrix< DDRMat > mTraction;
            moris::Cell< Matrix< DDRMat > > mdTractiondDof;
            moris::Cell< Matrix< DDRMat > > mdTractiondDv;

            Matrix< DDRMat > mTestTraction;
            moris::Cell< Matrix< DDRMat > > mdTestTractiondDof;
            moris::Cell< Matrix< DDRMat > > mdTestTractiondDv;

            Matrix< DDRMat > mStrain;
            moris::Cell< Matrix< DDRMat > > mdStraindDof;
            moris::Cell< Matrix< DDRMat > > mdStraindDv;

            Matrix< DDRMat > mTestStrain;

            Matrix< DDRMat > mConst;
            moris::Cell< Matrix< DDRMat > > mdConstdDof;
            moris::Cell< Matrix< DDRMat > > mdConstdDv;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Constitutive_Model(){};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            virtual ~Constitutive_Model(){};

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
                MORIS_ERROR( aSpaceDim > 0 && aSpaceDim < 4, "Constitutive_Model::set_space_dim - wrong space dimension.");

                // set space dimension
                mSpaceDim = aSpaceDim;
            }

//------------------------------------------------------------------------------
            /**
             * reset evaluation flags
             */
            void reset_eval_flags()
            {
                // reset the value flag
                mFluxEval         = true;
                mTractionEval     = true;
                mTestTractionEval = true;
                mStrainEval       = true;
                mTestStrainEval   = true;
                mConstEval        = true;

                // reset the dof derivative flag
                uint tNumDofTypes = mGlobalDofTypes.size();
                mdFluxdDofEval.assign( tNumDofTypes, true );
                mdTractiondDofEval.assign( tNumDofTypes, true );
                mdTestTractiondDofEval.assign( tNumDofTypes, true );
                mdStraindDofEval.assign( tNumDofTypes, true );
                mdConstdDofEval.assign( tNumDofTypes, true );

                // reset the dv derivative flag
                uint tNumDvTypes = mGlobalDvTypes.size();
                mdFluxdDvEval.assign( tNumDvTypes, true );
                mdTractiondDvEval.assign( tNumDvTypes, true );
                mdTestTractiondDvEval.assign( tNumDvTypes, true );
                mdStraindDvEval.assign( tNumDvTypes, true );
                mdConstdDvEval.assign( tNumDvTypes, true );

                // reset underlying properties
                for( std::shared_ptr< Property > tProp : mProperties )
                {
                    if ( tProp != nullptr )
                    {
                        tProp->reset_eval_flags();
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * set constitutive model dof types
             * @param[ in ] aDofTypes a cell of cell of dof types
             */
            void set_dof_type_list( moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes )
            {
                // set the dof types
                mDofTypes = aDofTypes;

                // build a map for the dof types
                this->build_dof_type_map();
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
            void build_dof_type_map()
            {
                // get number of dof types
                uint tNumDofTypes = mDofTypes.size();

                // determine the max Dof_Type enum
                sint tMaxEnum = 0;
                for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
                {
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mDofTypes( iDOF )( 0 ) ) );
                }
                tMaxEnum++;

                // set map size
                mDofTypeMap.set_size( tMaxEnum, 1, -1 );

                // loop over the dof types
                for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
                {
                    // fill the dof type map
                    mDofTypeMap( static_cast< int >( mDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
                }
            }

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
             * @param[ in ] aDvTypes a cell of cell of dv types
             */
            void set_dv_type_list( moris::Cell< moris::Cell< MSI::Dv_Type > > aDvTypes )
            {
                // set the dv types
                mDvTypes = aDvTypes;

                // build a map for the dv types
                this->build_dv_type_map();

            }

//------------------------------------------------------------------------------
            /**
             * return a cell of dv types
             * @param[ out ] aDvTypes a cell of cell of dv types
             */
            const moris::Cell< moris::Cell< MSI::Dv_Type > > & get_dv_type_list() const
            {
                return mDvTypes;
            };

//------------------------------------------------------------------------------
            /**
             * build a map for the dv types
             */
            void build_dv_type_map()
            {
                // get number of dv types
                uint tNumDvTypes = mDvTypes.size();

                // determine the max Dv_Type enum
                sint tMaxEnum = 0;
                for( uint iDV = 0; iDV < tNumDvTypes; iDV++ )
                {
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mDvTypes( iDV )( 0 ) ) );
                }
                tMaxEnum++;

                // set the Dv_Type map size
                mDvTypeMap.set_size( tMaxEnum, 1, -1 );

                // loop over the dv types
                for( uint iDV = 0; iDV < tNumDvTypes; iDV++ )
                {
                    // fill the dv type map
                    mDvTypeMap( static_cast< int >( mDvTypes( iDV )( 0 ) ), 0 ) = iDV;
                }
            }

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
             * set dof field interpolators
             * @param[ in ] aFieldInterpolators cell of dof field interpolator pointers
             */
            void set_dof_field_interpolators( moris::Cell< Field_Interpolator* > aFieldInterpolators )
            {
                // check input size
                MORIS_ASSERT( aFieldInterpolators.size() == this->get_global_dof_type_list().size(),
                              "Constitutive_Model::set_dof_field_interpolators - wrong input size. " );

                // check field interpolator type
                bool tCheckFI = true;
                for( uint iFI = 0; iFI < aFieldInterpolators.size(); iFI++ )
                {
                    tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dof_type()( 0 ) == get_global_dof_type_list()( iFI )( 0 ) );
                }
                MORIS_ASSERT( tCheckFI, "Constitutive_Model::set_dof_field_interpolators - wrong field interpolator dof type. ");

                // set field interpolators
                mDofFI = aFieldInterpolators;

                // set field interpolators for properties
                for( std::shared_ptr< Property > tProp : this->get_properties() )
                {
                    if (tProp != nullptr )
                    {
                        // get the list of dof types for the property
                        moris::Cell< moris::Cell< MSI::Dof_Type > > tPropDofTypes = tProp->get_dof_type_list();

                        // get the number of dof type for the property
                        uint tNumDofTypes = tPropDofTypes.size();

                        // set the size of the field interpolators list for the property
                        moris::Cell< Field_Interpolator* > tPropFIs( tNumDofTypes, nullptr );

                        // loop over the dof types
                        for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
                        {
                            // get the dof type index in the CM
                            uint tDofIndexInCM = this->get_global_dof_type_map()( static_cast< uint >( tPropDofTypes( iDof )( 0 ) ) );

                            // fill the field interpolators list for the property
                            tPropFIs( iDof ) = this->get_dof_field_interpolators()( tDofIndexInCM );
                        }

                        // set the field interpolators for the property
                        tProp->set_dof_field_interpolators( tPropFIs );
                    }
                }
            }

            void set_field_interpolator_manager( Field_Interpolator_Manager * aFieldInterpolatorManager )
            {
                mFieldInterpolatorManager = aFieldInterpolatorManager;
            }

//------------------------------------------------------------------------------
            /**
             * set dof field interpolators
             * @param[ in ] aFieldInterpolators cell of dof field interpolator pointers
             */
            void set_geometry_interpolator( Geometry_Interpolator* aGeometryInterpolator )
            {
                // set geometry interpolator for properties
                for( std::shared_ptr< Property > tProp : this->get_properties() )
                {
                    if( tProp != nullptr )
                    {
                        tProp->set_geometry_interpolator( aGeometryInterpolator );
                    }
                }
            }
//------------------------------------------------------------------------------
            /**
             * get dof field interpolators
             * @param[ out ] mDofFI cell of dof field interpolator pointers
             */
            const moris::Cell< Field_Interpolator* > & get_dof_field_interpolators()
            {
                return mDofFI;
            };

//------------------------------------------------------------------------------
            /**
             * check that dof field interpolators were assigned
             */
            void check_dof_field_interpolators()
            {
                // check field interpolators cell size
                MORIS_ASSERT( mDofFI.size() == this->get_global_dof_type_list().size(),
                              "Constitutive_Model::check_dof_field_interpolators - wrong FI size. " );

               // loop over the field interpolator pointers
               for( uint iFI = 0; iFI < this->get_global_dof_type_list().size(); iFI++ )
               {
                   // check that the field interpolator was set
                   MORIS_ASSERT( mDofFI( iFI ) != nullptr,
                                 "Constitutive_Model::check_dof_field_interpolators - FI missing. " );
               }
            }

//------------------------------------------------------------------------------
            /**
             * set dv field interpolators
             * @param[ in ] aFieldInterpolators cell of dv field interpolator pointers
             */
            void set_dv_field_interpolators( moris::Cell< Field_Interpolator* > aFieldInterpolators )
            {
                // get input size
                uint tNumInputFI = aFieldInterpolators.size();

                // check input size
                MORIS_ASSERT( tNumInputFI == this->get_global_dv_type_list().size(),
                              "Constitutive_Model::set_dv_field_interpolators - wrong input size. " );

                // check field interpolator type
                bool tCheckFI = true;
                for( uint iFI = 0; iFI < tNumInputFI; iFI++ )
                {
                    tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dv_type()( 0 ) == this->get_global_dv_type_list()( iFI )( 0 ) );
                }
                MORIS_ASSERT( tCheckFI, "Constitutive_Model::set_dv_field_interpolators - wrong field interpolator dv type. ");

                // set field interpolators
                mDvFI = aFieldInterpolators;
            }

//------------------------------------------------------------------------------
            /**
             * get dv field interpolators
             * @param[ out ] mDvFI cell of dv field interpolator pointers
             */
            const moris::Cell< Field_Interpolator* > & get_dv_field_interpolators()
            {
                return mDvFI;
            };

//------------------------------------------------------------------------------
            /**
             * check that dv field interpolators were assigned
             */
            void check_dv_field_interpolators()
            {
                // get num of dv types
                uint tNumDvTypes = this->get_global_dv_type_list().size();

                // check field interpolators cell size
                MORIS_ASSERT( mDvFI.size() == tNumDvTypes,
                              "Constitutive_Model::check_dv_field_interpolators - wrong FI size. " );

               // loop over the field interpolator pointers
               for( uint iFI = 0; iFI < tNumDvTypes; iFI++ )
               {
                   // check that the field interpolator was set
                   MORIS_ASSERT( mDvFI( iFI ) != nullptr,
                                 "Constitutive_Model::check_dv_field_interpolators - FI missing. " );
               }
            }

//------------------------------------------------------------------------------
            /**
             * set property
             * @param[ in ] aProperty a property pointer
             * @param[ in ] aPropertyString a string describing the property
             */
            virtual void set_property( std::shared_ptr< fem::Property > aProperty,
                                       std::string                      aPropertyType ) = 0;
//            {
//                MORIS_ERROR( false, "Constitutive_Model::set_property - This function does nothing." );
//            }

            virtual void set_model_type(fem::Model_Type aModelType)
            {
                MORIS_ERROR( false, "Constitutive_Model::set_model_type - This function does nothing." );
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
            void build_global_dof_type_list()
            {
                // get number of dof types
                uint tNumDofTypes = mDofTypes.size();

                // set the size of the dof type list
                uint tCounterMax = tNumDofTypes;

                for ( std::shared_ptr< Property > tProperty : mProperties )
                {
                    if( tProperty != nullptr )
                    {
                        tCounterMax += tProperty->get_dof_type_list().size();
                    }
                }
                mGlobalDofTypes.resize( tCounterMax );
                moris::Cell< sint > tCheckList( tCounterMax, -1 );

                // init total dof counter
                uint tCounter = 0;

                // get active dof type for constitutive model
                for ( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
                {
                    tCheckList( tCounter ) = static_cast< uint >( mDofTypes( iDOF )( 0 ) );
                    mGlobalDofTypes( tCounter ) = mDofTypes( iDOF );
                    tCounter++;
                }

                for ( std::shared_ptr< Property > tProperty : mProperties )
                {
                    if( tProperty != nullptr )
                    {
                        // get active dof types
                        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType
                        = tProperty->get_dof_type_list();

                        for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                        {
                            // check enum is not already in the list
                            bool tCheck = false;
                            for( uint i = 0; i < tCounter; i++ )
                            {
                                tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDofType( iDOF )( 0 ) ) );
                            }

                            // if dof enum not in the list
                            if ( !tCheck )
                            {
                                tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );
                                mGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );
                                tCounter++;
                            }
                        }
                    }
                }

                // get the number of unique dof type groups, i.e. the number of interpolators
                mGlobalDofTypes.resize( tCounter );

                // number of dof types
                uint tNumGlobalDofTypes = mGlobalDofTypes.size();

                // set flag for evaluation
                mdFluxdDofEval.resize( tNumGlobalDofTypes, true );
                mdTractiondDofEval.resize( tNumGlobalDofTypes, true );
                mdTestTractiondDofEval.resize( tNumGlobalDofTypes, true );
                mdStraindDofEval.resize( tNumGlobalDofTypes, true );
                mdConstdDofEval.resize( tNumGlobalDofTypes, true );

                // set storage for evaluation
                mdFluxdDof.resize( tNumGlobalDofTypes );
                mdTractiondDof.resize( tNumGlobalDofTypes );
                mdTestTractiondDof.resize( tNumGlobalDofTypes );
                mdStraindDof.resize( tNumGlobalDofTypes );
                mdConstdDof.resize( tNumGlobalDofTypes );
            };

//------------------------------------------------------------------------------
            /**
             * get global dof type list
             * @param[ out ] mGlobalDofTypes global list of dof type
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_type_list()
            {
                if( mGlobalDofBuild )
                {
                    // build the stabilization parameter global dof type list
                    this->build_global_dof_type_list();

                    // build the stabilization parameter global dof type map
                    this->build_global_dof_type_map();

                    // update build flag
                    mGlobalDofBuild = false;
                }

                return mGlobalDofTypes;
            };

//------------------------------------------------------------------------------
            /**
             * get non unique list of dof type for the constitutive model
             * @param[ in ] aDofTypes a cell of dof type to fill
             */
            void get_non_unique_dof_types( moris::Cell< MSI::Dof_Type > & aDofTypes )
            {
                // init dof counter
                uint tCounter = 0;

                // loop over direct dof dependencies
                for ( uint iDOF = 0; iDOF < mDofTypes.size(); iDOF++ )
                {
                    // update counter
                    tCounter += mDofTypes( iDOF ).size();
                }

                // loop over properties
                for ( std::shared_ptr< Property > tProperty : mProperties )
                {
                    if ( tProperty != nullptr )
                    {
                        // get property dof type list
                        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

                        // loop over property dof types
                        for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                        {
                            // update counter
                            tCounter += tActiveDofType( iDOF ).size();
                        }
                    }
                }

                // reserve memory for the non unique dof type list
                aDofTypes.reserve( tCounter );

                // loop over direct dof dependencies
                for ( uint iDOF = 0; iDOF < mDofTypes.size(); iDOF++ )
                {
                    // populate the dof type list
                    aDofTypes.append( mDofTypes( iDOF ) );
                }

                // loop over the properties
                for ( std::shared_ptr< Property > tProperty : mProperties )
                {
                    if ( tProperty != nullptr )
                    {
                        // get property dof type list
                        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

                        // loop over property dof types
                        for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                        {
                            // populate te dof type list
                            aDofTypes.append( tActiveDofType( iDOF ) );
                        }
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * build global dof type map
             */
            void build_global_dof_type_map()
            {
                // get number of global dof types
                uint tNumDofTypes = mGlobalDofTypes.size();

                // determine the max Dof_Type enum
                sint tMaxEnum = 0;
                for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
                {
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mGlobalDofTypes( iDOF )( 0 ) ) );
                }
                tMaxEnum++;

                // set the Dof_Type map size
                mGlobalDofTypeMap.set_size( tMaxEnum, 1, -1 );

                // fill the Dof_Type map
                for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
                {
                    // fill the property map
                    mGlobalDofTypeMap( static_cast< int >( mGlobalDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
                }
            }

//------------------------------------------------------------------------------
            /**
             * get global dof type map
             * @param[ out ] mGlobalDofTypeMap global map for the dof types
             */
            const Matrix< DDSMat > & get_global_dof_type_map()
            {
                return mGlobalDofTypeMap;
            }

//------------------------------------------------------------------------------
            /**
             * check dependency on a given group of dof types
             * @param[ in ]  aDofType       a group of dof types
             * @param[ out ] tDofDependency a bool true if dependency on dof type
             *
             */
            bool check_dof_dependency( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
                // set bool for dependency
                bool tDofDependency = false;

                // get dof type index
                uint tDofIndex = static_cast< uint >( aDofType( 0 ) );

                // if aDofType is an active dv type for the constitutive model
                if( tDofIndex < mGlobalDofTypeMap.numel() && mGlobalDofTypeMap( tDofIndex ) != -1 )
                {
                    // bool is set to true
                    tDofDependency = true;
                }
                // return bool for dependency
                return tDofDependency;
            }

//------------------------------------------------------------------------------
            /**
             * create a global dv type list including constitutive and property dependencies
             */
            void build_global_dv_type_list()
            {
                // get number of dv types
                uint tNumDvTypes = mDvTypes.size();

                // set the size of the dv type list
                uint tCounterMax = tNumDvTypes;

                for ( std::shared_ptr< Property > tProperty : mProperties )
                {
                    if( tProperty != nullptr )
                    {
                        tCounterMax += tProperty->get_dv_type_list().size();
                    }
                }
                mGlobalDvTypes.resize( tCounterMax );
                moris::Cell< sint > tCheckList( tCounterMax, -1 );

                // init total dv counter
                uint tCounter = 0;

                // get active dv type for constitutive model
                for ( uint iDV = 0; iDV < tNumDvTypes; iDV++ )
                {
                    tCheckList( tCounter ) = static_cast< uint >( mDvTypes( iDV )( 0 ) );
                    mGlobalDvTypes( tCounter ) = mDvTypes( iDV );
                    tCounter++;
                }

                for ( std::shared_ptr< Property > tProperty : mProperties )
                {
                    if( tProperty != nullptr )
                    {
                        // get active dv types
                        moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvType = tProperty->get_dv_type_list();

                        for ( uint iDV = 0; iDV < tActiveDvType.size(); iDV++ )
                        {
                            // check enum is not already in the list
                            bool tCheck = false;
                            for( uint i = 0; i < tCounter; i++ )
                            {
                                tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDV )( 0 ) ) );
                            }

                            // if dof enum not in the list
                            if ( !tCheck )
                            {
                                tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDV )( 0 ) );
                                mGlobalDvTypes( tCounter ) = tActiveDvType( iDV );
                                tCounter++;
                            }
                        }
                    }
                }

                // get the number of unique dv type groups, i.e. the number of interpolators
                mGlobalDvTypes.resize( tCounter );

                // build global dv type map
                this->build_global_dv_type_map();

                // number of dof types
                uint tNumGlobalDvTypes = mGlobalDvTypes.size();

                // set flag for evaluation
                mdFluxdDvEval.assign( tNumGlobalDvTypes, true );
                mdTractiondDvEval.assign( tNumGlobalDvTypes, true );
                mdTestTractiondDvEval.assign( tNumGlobalDvTypes, true );
                mdStraindDvEval.assign( tNumGlobalDvTypes, true );
                mdConstdDvEval.assign( tNumGlobalDvTypes, true );

                // set storage for evaluation
                mdFluxdDv.resize( tNumGlobalDvTypes );
                mdTractiondDv.resize( tNumGlobalDvTypes );
                mdTestTractiondDv.resize( tNumGlobalDvTypes );
                mdStraindDv.resize( tNumGlobalDvTypes );
                mdConstdDv.resize( tNumGlobalDvTypes );
            };

//------------------------------------------------------------------------------
            /**
             * get global dv type list
             * @param[ out ] mGlobalDvTypes global list of dv types
             */
            const moris::Cell< moris::Cell< MSI::Dv_Type > > & get_global_dv_type_list()
            {
                if( mGlobalDvBuild )
                {
                    // build the stabilization parameter global dv type list
                    this->build_global_dv_type_list();

                    // build the stabilization parameter global dv type map
                    this->build_global_dv_type_map();

                    // update build flag
                    mGlobalDvBuild = false;
                }

                return mGlobalDvTypes;
            };

//------------------------------------------------------------------------------
            /**
             * build global dv type list map
             */
            void build_global_dv_type_map()
            {
                // get number of global dv types
                uint tNumDvTypes = mGlobalDvTypes.size();

                // determine the max Dv_Type enum
                sint tMaxEnum = 0;
                for( uint iDOF = 0; iDOF < tNumDvTypes; iDOF++ )
                {
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mGlobalDvTypes( iDOF )( 0 ) ) );
                }
                tMaxEnum++;

                // set the Dv_Type map size
                mGlobalDvTypeMap.set_size( tMaxEnum, 1, -1 );

                // loop over the dv types
                for( uint iDV = 0; iDV < tNumDvTypes; iDV++ )
                {
                    // fill the dv type map
                    mGlobalDvTypeMap( static_cast< int >( mGlobalDvTypes( iDV )( 0 ) ), 0 ) = iDV;
                }
            }

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
            bool check_dv_dependency( const moris::Cell< MSI::Dv_Type > & aDvType )
            {
                // set bool for dependency
                bool tDvDependency = false;

                // get dv type index
                uint tDvIndex = static_cast< uint >( aDvType( 0 ) );

                // if aDvType is an active dv type for the constitutive model
                if( tDvIndex < mGlobalDvTypeMap.numel() && mGlobalDvTypeMap( tDvIndex ) != -1 )
                {
                    // bool is set to true
                    tDvDependency = true;
                }
                // return bool for dependency
                return tDvDependency;
            }

//------------------------------------------------------------------------------
            /**
             * get the constitutive model flux
             * @param[ out ] mFlux constitutive model flux
             */
            const Matrix< DDRMat > & flux()
            {
                // if the flux was not evaluated
                if( mFluxEval )
                {
                    // evaluate the flux
                    this->eval_flux();

                    // set bool for evaluation
                    mFluxEval = false;
                }
                // return the flux value
                return mFlux;
            }

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
             * get the constitutive model traction
             * @param[ in ]  aNormal   normal
             * @param[ out ] mTraction constitutive model traction
             */
            const Matrix< DDRMat > & traction( const Matrix< DDRMat > & aNormal )
            {
                // if the traction was not evaluated
                if( mTractionEval )
                {
                    // evaluate the traction
                    this->eval_traction( aNormal );

                    // set bool for evaluation
                    mTractionEval = false;
                }
                // return the traction value
                return mTraction;
            }

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
            const Matrix< DDRMat > & testTraction( const Matrix< DDRMat > & aNormal )
            {
                // if the test traction was not evaluated
                if( mTestTractionEval )
                {
                    // evaluate the test traction
                    this->eval_testTraction( aNormal );

                    // set bool for evaluation
                    mTestTractionEval = false;
                }
                // return the test traction value
                return mTestTraction;
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test traction
             * @param[ in ]  aNormal normal
             */
            virtual void eval_testTraction( const Matrix< DDRMat > & aNormal )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_testTraction - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the constitutive model strain
             * @param[ out ] mStrain constitutive model strain
             */
            const Matrix< DDRMat > & strain()
            {
                // if the strain was not evaluated
                if( mStrainEval )
                {
                    // evaluate the strain
                    this->eval_strain();

                    // set bool for evaluation
                    mStrainEval = false;
                }
                // return the strain value
                return mStrain;
            }

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
             * get the constitutive model test strain
             * @param[ out ] mTestStrain constitutive model test strain
             */
            const Matrix< DDRMat > & testStrain()
            {
                // if the test strain was not evaluated
                if( mTestStrainEval )
                {
                    // evaluate the test strain
                    this->eval_testStrain();

                    // set bool for evaluation
                    mTestStrainEval = false;
                }
                // return the test strain value
                return mTestStrain;
            }

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
            const Matrix< DDRMat > & constitutive()
            {
                // if the constitutive matrix was not evaluated
                if( mConstEval )
                {
                    // evaluate the constitutive matrix
                    this->eval_const();

                    // set bool for evaluation
                    mConstEval = false;
                }
                // return the constitutive matrix value
                return mConst;
            }

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
             * get the derivative of the flux wrt dof
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ out ] mFluxDofDer derivative of the flux wrt dof
             */
            const Matrix< DDRMat > & dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
               // if aDofType is not an active dof type for the property
               MORIS_ERROR( this->check_dof_dependency( aDofType ), "Constitutive_Model::dFluxdDOF - no dependency in this dof type." );

               // get the dof index
               uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mdFluxdDofEval( tDofIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dFluxdDOF( aDofType );

                   // set bool for evaluation
                   mdFluxdDofEval( tDofIndex ) = false;
               }

               // return the derivative
               return mdFluxdDof( tDofIndex );
            }

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
             * get the derivative of the traction wrt dof
             * @param[ in ]  aDofType        group of dof type
             * @param[ in ]  aNormal         normal
             * @param[ out ] mTractionDofDer derivative of the traction wrt dof
             */
            const Matrix< DDRMat > & dTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofType,
                                                    const Matrix< DDRMat >             & aNormal )
            {
               // if aDofType is not an active dof type for the property
               MORIS_ERROR( this->check_dof_dependency( aDofType ), "Constitutive_Model::dTractiondDOF - no dependency in this dof type." );

               // get the dof index
               uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mdTractiondDofEval( tDofIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dTractiondDOF( aDofType, aNormal );

                   // set bool for evaluation
                   mdTractiondDofEval( tDofIndex ) = false;
               }

               // return the derivative
               return mdTractiondDof( tDofIndex );
            }

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
                                                        const Matrix< DDRMat >             & aNormal )
            {
                // if aDofType is not an active dof type for the property
                MORIS_ERROR( this->check_dof_dependency( aDofType ), "Constitutive_Model::dTestTractiondDOF - no dependency in this dof type." );

                // get the dof index
                uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

                // if the derivative has not been evaluated yet
                if( mdTestTractiondDofEval( tDofIndex ) )
                {
                    // evaluate the derivative
                    this->eval_dTestTractiondDOF( aDofType, aNormal );

                    // set bool for evaluation
                    mdTestTractiondDofEval( tDofIndex ) = false;
                }

                // return the derivative
                return mdTestTractiondDof( tDofIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model test traction derivative wrt to a dof type
             * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ in ] adTractiondDOF a matrix to fill with derivative evaluation
             */
            virtual void eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                 const Matrix< DDRMat >             & aNormal )
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
                                                        const Matrix< DDRMat >             & aJump )
            {
               // if aDofType is not an active dof type for the property
               MORIS_ERROR( this->check_dof_dependency( aDofType ), "Constitutive_Model::dTestTractiondDOF - no dependency in this dof type." );

               // get the dof index
               uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mdTestTractiondDofEval( tDofIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dTestTractiondDOF( aDofType, aNormal, aJump );

                   // set bool for evaluation
                   mdTestTractiondDofEval( tDofIndex ) = false;
               }

               // return the derivative
               return mdTestTractiondDof( tDofIndex );
            }

//------------------------------------------------------------------------------
            /**
             * FIXME this is not a test traction, used for elast lin iso!!!!
             * evaluate the constitutive model test traction derivative wrt to a dof type
             * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ in ] adTractiondDOF a matrix to fill with derivative evaluation
             */
            virtual void eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                 const Matrix< DDRMat >             & aNormal,
                                                 const Matrix< DDRMat >             & aJump )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dTestTractiondDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the strain wrt dof
             * @param[ in ] aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ out ] mdStraindDof derivative of the strain wrt dof
             */
            const Matrix< DDRMat > & dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
               // if aDofType is not an active dof type for the property
               MORIS_ERROR( this->check_dof_dependency( aDofType ), "Constitutive_Model::dStraindDOF - no dependency in this dof type." );

               // get the dof index
               uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mdStraindDofEval( tDofIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dStraindDOF( aDofType );

                   // set bool for evaluation
                   mdStraindDofEval( tDofIndex ) = false;
               }

               // return the derivative
               return mdStraindDof( tDofIndex );
            }

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
            const Matrix< DDRMat > & dConstdDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
               // if aDofType is not an active dof type for the property
               MORIS_ERROR( this->check_dof_dependency( aDofType ), "Constitutive_Model::dConstdDOF - no dependency in this dof type." );

               // get the dof index
               uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mdConstdDofEval( tDofIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dConstdDOF( aDofType );

                   // set bool for evaluation
                   mdConstdDofEval( tDofIndex ) = false;
               }

               // return the derivative
               return mdConstdDof( tDofIndex );
            }

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
            const Matrix< DDRMat > & dFluxdDV( const moris::Cell< MSI::Dv_Type > & aDvType )
            {
               // if aDvType is not an active dv type
               MORIS_ERROR( this->check_dv_dependency( aDvType ), "Constitutive_Model::dFluxdDV - no dependency in this dv type." );

               // get the dv index
               uint tDvIndex = mGlobalDvTypeMap( static_cast< uint >( aDvType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mdFluxdDvEval( tDvIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dFluxdDV( aDvType );

                   // set bool for evaluation
                   mdFluxdDvEval( tDvIndex ) = false;
               }

               // return the derivative
               return mdFluxdDv( tDvIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux derivative wrt to a dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             */
            virtual void eval_dFluxdDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dFluxdDV - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the strain wrt dv
             * @param[ in ]  aDvTypes    a dv type wrt which the derivative is evaluated
             * @param[ out ] mdStraindDv derivative of the strain wrt dv
             */
            const Matrix< DDRMat > & dStraindDV( const moris::Cell< MSI::Dv_Type > & aDvType )
            {
               // if aDvType is not an active dv type for the property
               MORIS_ERROR( this->check_dv_dependency( aDvType ), "Constitutive_Model::dStraindDV - no dependency in this dv type." );

               // get the dv index
               uint tDvIndex = mGlobalDvTypeMap( static_cast< uint >( aDvType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mdStraindDvEval( tDvIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dStraindDV( aDvType );

                   // set bool for evaluation
                   mdStraindDvEval( tDvIndex ) = false;
               }

               // return the derivative
               return mdStraindDv( tDvIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain derivative wrt to a dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             */
            virtual void eval_dStraindDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dStraindDV - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the constitutive matrix wrt dv
             * @param[ in ]  aDvTypes   a dv type wrt which the derivative is evaluated
             * @param[ out ] mdConstdDv derivative of the constitutive matrix wrt dv
             */
            const Matrix< DDRMat > & dConstdDV( const moris::Cell< MSI::Dv_Type > & aDvType )
            {
               // if aDvType is not an active dv type for the property
               MORIS_ERROR( this->check_dv_dependency( aDvType ), "Constitutive_Model::dConstdDV - no dependency in this dv type." );

               // get the dv index
               uint tDvIndex = mGlobalDvTypeMap( static_cast< uint >( aDvType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mdConstdDvEval( tDvIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dConstdDV( aDvType );

                   // set bool for evaluation
                   mdConstdDvEval( tDvIndex ) = false;
               }

               // return the derivative
               return mdConstdDv( tDvIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model constitutive matrix derivative wrt to a dv type
             * @param[ in ] aDvTypes a dof type wrt which the derivative is evaluated
             */
            virtual void eval_dConstdDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
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
                                    Matrix< DDRMat >                   & adFluxdDOF_FD,
                                    real                                 aPerturbation )
            {
                // get the index for the considered dof type
                uint iFI = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ), 0 );

                // get number of coefficients, fields and bases for the considered FI
                uint tDerNumDof    = mDofFI( iFI )->get_number_of_space_time_coefficients();
                uint tDerNumBases  = mDofFI( iFI )->get_number_of_space_time_bases();
                uint tDerNumFields = mDofFI( iFI )->get_number_of_fields();

                // FIXME works only for diffusion
                // set size for derivative
                adFluxdDOF_FD.set_size( mSpaceDim, tDerNumDof, 0.0 );

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = mDofFI( iFI )->get_coeff();

                // init dof counter
                uint tDofCounter = 0;

                // loop over coefficients columns
                for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over coefficients rows
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // perturbation of the coefficent
                        Matrix< DDRMat > tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        mDofFI( iFI )->set_coeff( tCoeffPert );

                        // reset properties
                        uint tNumProps = mProperties.size();
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tFlux_Plus = this->flux();

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        mDofFI( iFI )->set_coeff( tCoeffPert );

                        // reset properties
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tFlux_Minus = this->flux();

                        // evaluate Jacobian
                        adFluxdDOF_FD.get_column( tDofCounter ) = ( tFlux_Plus - tFlux_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                mDofFI( iFI )->set_coeff( tCoeff );
            }

//------------------------------------------------------------------------------
            /**
            * evaluate the constitutive model strain derivative wrt to a dof type
            * @param[ in ] aDofTypes       a dof type wrt which the derivative is evaluated
            * @param[ in ] adStraindDOF_FD a matrix to fill with derivative evaluation
            * @param[ in ] aPerturbation   real to perturb for FD
            */
            void eval_dStraindDOF_FD( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                      Matrix< DDRMat >                   & adStraindDOF_FD,
                                      real                                 aPerturbation )
            {
                // get the index for the considered dof type
                uint iFI = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ), 0 );

                // get number of master dofs wrt which derivative is computed
                uint tDerNumDof    = mDofFI( iFI )->get_number_of_space_time_coefficients();
                uint tDerNumBases  = mDofFI( iFI )->get_number_of_space_time_bases();
                uint tDerNumFields = mDofFI( iFI )->get_number_of_fields();

                // FIXME works only for diffusion
                // set size for derivative
                adStraindDOF_FD.set_size( mSpaceDim, tDerNumDof, 0.0 );

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = mDofFI( iFI )->get_coeff();

                // init dof counter
                uint tDofCounter = 0;

                // loop over coefficients columns
                for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over coefficients rows
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // perturbation of the coefficent
                        Matrix< DDRMat > tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        mDofFI( iFI )->set_coeff( tCoeffPert );

                        // reset properties
                        uint tNumProps = mProperties.size();
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tStrain_Plus = this->strain();

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        mDofFI( iFI )->set_coeff( tCoeffPert );

                        // reset properties
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tStrain_Minus = this->strain();

                        // evaluate Jacobian
                        adStraindDOF_FD.get_column( tDofCounter ) = ( tStrain_Plus - tStrain_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                mDofFI( iFI )->set_coeff( tCoeff );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model stress derivative wrt to a dv type
             * @param[ in ] aDofTypes     a dof type wrt which the derivative is evaluated
             * @param[ in ] adFluxdDV_FD  a matrix to fill with derivative evaluation
             * @param[ in ] aPerturbation real to perturb for FD
             */
            void eval_dFluxdDV_FD( const moris::Cell< MSI::Dv_Type > & aDvTypes,
                                    Matrix< DDRMat >                  & adFluxdDV_FD,
                                    real                                aPerturbation )
            {
                // get the index for the considered dv type
                uint iFI = mGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ), 0 );

                // get number of coefficients, fields and bases for the considered FI
                uint tDerNumDv     = mDvFI( iFI )->get_number_of_space_time_coefficients();
                uint tDerNumBases  = mDvFI( iFI )->get_number_of_space_time_bases();
                uint tDerNumFields = mDvFI( iFI )->get_number_of_fields();

                // FIXME works only for diffusion
                // set size for derivative
                adFluxdDV_FD.set_size( mSpaceDim, tDerNumDv, 0.0 );

                // coefficients for dv type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = mDvFI( iFI )->get_coeff();

                // init dv counter
                uint tDvCounter = 0;

                // loop over coefficients columns
                for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over coefficients rows
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // perturbation of the coefficent
                        Matrix< DDRMat > tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        mDvFI( iFI )->set_coeff( tCoeffPert );

                        // reset properties
                        uint tNumProps = mProperties.size();
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tFlux_Plus = this->flux();

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        mDvFI( iFI )->set_coeff( tCoeffPert );

                        // reset properties
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tFlux_Minus = this->flux();

                        // evaluate Jacobian
                        adFluxdDV_FD.get_column( tDvCounter ) = ( tFlux_Plus - tFlux_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dv counter
                        tDvCounter++;
                    }
                }
                // reset the coefficients values
                mDvFI( iFI )->set_coeff( tCoeff );
            }

//------------------------------------------------------------------------------
            /**
            * evaluate the constitutive model strain derivative wrt to a dv type
            * @param[ in ] aDvTypes       a dv type wrt which the derivative is evaluated
            * @param[ in ] adStraindDV_FD a matrix to fill with derivative evaluation
            * @param[ in ] aPerturbation  real to perturb for FD
            */
            void eval_dStraindDV_FD( const moris::Cell< MSI::Dv_Type > & aDvTypes,
                                     Matrix< DDRMat >                  & adStraindDV_FD,
                                     real                                aPerturbation )
            {
                // get the index for the considered dof type
                uint iFI = mGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ), 0 );

                // get number of coefficients, fields and bases for the considered FI
                uint tDerNumDv     = mDvFI( iFI )->get_number_of_space_time_coefficients();
                uint tDerNumBases  = mDvFI( iFI )->get_number_of_space_time_bases();
                uint tDerNumFields = mDvFI( iFI )->get_number_of_fields();

                // FIXME works only for diffusion
                // set size for derivative
                adStraindDV_FD.set_size( mSpaceDim, tDerNumDv, 0.0 );

                // coefficients for dv type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = mDvFI( iFI )->get_coeff();

                // init dv counter
                uint tDvCounter = 0;

                // loop over coefficients columns
                for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over coefficients rows
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // perturbation of the coefficent
                        Matrix< DDRMat > tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        mDvFI( iFI )->set_coeff( tCoeffPert );

                        // reset properties
                        uint tNumProps = mProperties.size();
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tStrain_Plus = this->strain();

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        mDvFI( iFI )->set_coeff( tCoeffPert );

                        // reset properties
                        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                        {
                            mProperties( iProp )->reset_eval_flags();
                        }

                        // reset constitutive model
                        this->reset_eval_flags();

                        // evaluate the residual
                        Matrix< DDRMat > tStrain_Minus = this->strain();

                        // evaluate Jacobian
                        adStraindDV_FD.get_column( tDvCounter ) = ( tStrain_Plus - tStrain_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dv counter
                        tDvCounter++;
                    }
                }
                // reset the coefficients values
                mDvFI( iFI )->set_coeff( tCoeff );
            }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_CONSTITUTIVE_MODEL_HPP_ */

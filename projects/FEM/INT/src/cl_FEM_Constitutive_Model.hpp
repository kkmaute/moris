/*
 * cl_FEM_IWG.hpp
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
//------------------------------------------------------------------------------
        /**
         * Constitutive model
         */
        class Constitutive_Model
        {

        protected :

            // constitutive model type
            fem::Constitutive_Type mConstitutiveType;

            // dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mDofTypes;

            // dof type map
            Matrix< DDSMat > mDofTypeMap;

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

            // property type list
            moris::Cell< fem::Property_Type > mPropTypes;

            // properties
            moris::Cell< Property* > mProperties;

            // spatial dimensions
            uint mSpaceDim;

            // flag for evaluation
            bool mFluxEval = true;
            bool mStrainEval = true;
            bool mConstEval  = true;
            moris::Cell< bool > mFluxDofDerEval;
            moris::Cell< bool > mFluxDvDerEval;
            moris::Cell< bool > mStrainDofDerEval;
            moris::Cell< bool > mStrainDvDerEval;
            moris::Cell< bool > mConstDofDerEval;
            moris::Cell< bool > mConstDvDerEval;

            // storage for evaluation
            Matrix< DDRMat > mFlux;
            Matrix< DDRMat > mStrain;
            Matrix< DDRMat > mConst;
            moris::Cell< Matrix< DDRMat > > mFluxDofDer;
            moris::Cell< Matrix< DDRMat > > mFluxDvDer;
            moris::Cell< Matrix< DDRMat > > mStrainDofDer;
            moris::Cell< Matrix< DDRMat > > mStrainDvDer;
            moris::Cell< Matrix< DDRMat > > mConstDofDer;
            moris::Cell< Matrix< DDRMat > > mConstDvDer;

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
            /**
             * set a constitutive type for the constitutive model
             * @param[ in ] aConstitutiveType a constitutive type
             */
            void set_constitutive_type( const fem::Constitutive_Type aConstitutiveType )
            {
                mConstitutiveType = aConstitutiveType;
            }

//------------------------------------------------------------------------------
            /**
             * return a constitutive type for the constitutive model
             * @param[ out ] mConstitutiveType a constitutive type
             */
            const fem::Constitutive_Type & get_constitutive_type() const
            {
                return mConstitutiveType;
            }

//------------------------------------------------------------------------------
            /**
             * set space dimension
             * @param[ in ] aSpaceDim a spatial dimension
             */
            void set_space_dim( uint aSpaceDim )
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
                mFluxEval = true;
                mStrainEval = true;
                mConstEval  = true;

                // reset the dof derivative flag
                uint tNumDofTypes = mGlobalDofTypes.size();
                mFluxDofDerEval.resize( tNumDofTypes, true );
                mStrainDofDerEval.resize( tNumDofTypes, true );
                mConstDofDerEval.resize( tNumDofTypes, true );

                // reset the dv derivative flag
                uint tNumDvTypes = mGlobalDvTypes.size();
                mFluxDvDerEval.resize( tNumDvTypes, true );
                mStrainDvDerEval.resize( tNumDvTypes, true );
                mConstDvDerEval.resize( tNumDvTypes, true );
            }

//------------------------------------------------------------------------------
            /**
             * set constitutive model dof types
             * @param[ in ] aDofTypes a cell of cell of dof types
             */
            void set_dof_type_list( const moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes )
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
            void set_dv_type_list( const moris::Cell< moris::Cell< MSI::Dv_Type > > & aDvTypes )
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
             * set constitutive model property types
             * @param[ in ] aPropertyTypes a cell of property types
             */
            void set_property_type_list( const moris::Cell< fem::Property_Type  > & aPropertyTypes )
            {
                // set the property types
                mPropTypes = aPropertyTypes;
            }

//------------------------------------------------------------------------------
            /**
             * return a cell of property type for the constitutive model
             * @param[ out ] mPropTypes a cell of property types
             */
            const moris::Cell< fem::Property_Type > & get_property_type_list() const
            {
                return mPropTypes;
            }

//------------------------------------------------------------------------------
            /**
             * set dof field interpolators
             * @param[ in ] aFieldInterpolators cell of dof field interpolator pointers
             */
            void set_dof_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators )
            {
                // check input size
                MORIS_ASSERT( aFieldInterpolators.size() == mGlobalDofTypes.size(),
                              "Constitutive_Model::set_field_interpolators - wrong input size. " );

                // check field interpolator type
                bool tCheckFI = true;
                for( uint iFI = 0; iFI < aFieldInterpolators.size(); iFI++ )
                {
                    tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dof_type()( 0 ) == mGlobalDofTypes( iFI )( 0 ) );
                }
                MORIS_ASSERT( tCheckFI, "Constitutive_Model::set_field_interpolators - wrong field interpolator dof type. ");

                // set field interpolators
                mDofFI = aFieldInterpolators;
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
                MORIS_ASSERT( mDofFI.size() == mGlobalDofTypes.size(),
                              "Constitutive_Model::check_dof_field_interpolators - wrong FI size. " );

               // loop over the field interpolator pointers
               for( uint iFI = 0; iFI < mGlobalDofTypes.size(); iFI++ )
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
            void set_dv_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators )
            {
                // get input size
                uint tNumInputFI = aFieldInterpolators.size();

                // check input size
                MORIS_ASSERT( tNumInputFI == mGlobalDvTypes.size(),
                              "Constitutive_Model::set_dv_field_interpolators - wrong input size. " );

                // check field interpolator type
                bool tCheckFI = true;
                for( uint iFI = 0; iFI < tNumInputFI; iFI++ )
                {
                    tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dv_type()( 0 ) == mGlobalDvTypes( iFI )( 0 ) );
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
                uint tNumDvTypes = mGlobalDvTypes.size();

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
             * set properties
             * @param[ in ] aProperties cell of property pointers
             */
            void set_properties( moris::Cell< Property* > & aProperties )
            {
                // get input size
                uint tNumInputProp = aProperties.size();

                // check input size
                MORIS_ASSERT( tNumInputProp == mPropTypes.size(),
                              "Constitutive_Model::set_properties - master, wrong input size. " );

                // check property type
                bool tCheckProp = true;
                for( uint iProp = 0; iProp < tNumInputProp; iProp++ )
                {
                    tCheckProp = tCheckProp && ( aProperties( iProp )->get_property_type() == mPropTypes( iProp ) );
                }
                MORIS_ASSERT( tCheckProp, "Constitutive_Model::set_properties - wrong property type. ");

                // set properties
                mProperties = aProperties;

                // create a global dof type list
                this->build_global_dof_type_list();

                // create a global dv type list
                this->build_global_dv_type_list();
            }

//------------------------------------------------------------------------------
            /**
             * get properties
             * @param[ out ] mProperties cell of property pointers
             */
            const moris::Cell< Property* > & get_properties()
            {
                return mProperties;
            };

//------------------------------------------------------------------------------
            /**
             * check that properties were assigned
             */
            void check_properties()
            {
                // get number of property types
                uint tNumPropTypes = mPropTypes.size();

                // check master properties cell size
                MORIS_ASSERT( mProperties.size() == tNumPropTypes,
                              "Constitutive_Model::check_properties - wrong property size. " );

                // loop over all master properties and check that they are assigned
                for( uint iProp = 0; iProp < tNumPropTypes; iProp++ )
                {
                    MORIS_ASSERT( mProperties( iProp ) != nullptr,
                                  "Constitutive_Model::check_properties - property missing. " );
                }
            }

//------------------------------------------------------------------------------
            /**
             * create a global dof type list including constitutive and property dependencies
             */
            void build_global_dof_type_list()
            {
                // get number of property types
                uint tNumPropTypes = mDofTypes.size();

                // set the size of the dof type list
                uint tCounterMax = tNumPropTypes;

                for ( Property* tProperty : mProperties )
                {
                    tCounterMax += tProperty->get_dof_type_list().size();
                }
                mGlobalDofTypes.resize( tCounterMax );
                moris::Cell< sint > tCheckList( tCounterMax, -1 );

                // init total dof counter
                uint tCounter = 0;

                // get active dof type for constitutive model
                for ( uint iDOF = 0; iDOF < tNumPropTypes; iDOF++ )
                {
                    tCheckList( tCounter ) = static_cast< uint >( mDofTypes( iDOF )( 0 ) );
                    mGlobalDofTypes( tCounter ) = mDofTypes( iDOF );
                    tCounter++;
                }

                for ( Property* tProperty : mProperties )
                {
                    // get active dof types
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

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

                // get the number of unique dof type groups, i.e. the number of interpolators
                mGlobalDofTypes.resize( tCounter );

                // build global dof type map
                this->build_global_dof_type_map();

                // number of dof types
                uint tNumGlobalDofTypes = mGlobalDofTypes.size();

                // set flag for evaluation
                mFluxDofDerEval.resize( tNumGlobalDofTypes, true );
                mStrainDofDerEval.resize( tNumGlobalDofTypes, true );
                mConstDofDerEval.resize( tNumGlobalDofTypes, true );

                // set storage for evaluation
                mFluxDofDer.resize( tNumGlobalDofTypes );
                mStrainDofDer.resize( tNumGlobalDofTypes );
                mConstDofDer.resize( tNumGlobalDofTypes );

            };

//------------------------------------------------------------------------------
            /**
             * get global dof type list
             * @param[ out ] mGlobalDofTypes global list of dof type
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_type_list()
            {
                return mGlobalDofTypes;
            };

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

                for ( Property* tProperty : mProperties )
                {
                    tCounterMax += tProperty->get_dv_type_list().size();
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

                for ( Property* tProperty : mProperties )
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

                // get the number of unique dv type groups, i.e. the number of interpolators
                mGlobalDvTypes.resize( tCounter );

                // build global dv type map
                this->build_global_dv_type_map();

                // number of dof types
                uint tNumGlobalDvTypes = mGlobalDvTypes.size();

                // set flag for evaluation
                mFluxDvDerEval.resize( tNumGlobalDvTypes, true );
                mStrainDvDerEval.resize( tNumGlobalDvTypes, true );
                mConstDvDerEval.resize( tNumGlobalDvTypes, true );

                // set storage for evaluation
                mFluxDvDer.resize( tNumGlobalDvTypes );
                mStrainDvDer.resize( tNumGlobalDvTypes );
                mConstDvDer.resize( tNumGlobalDvTypes );
            };

//------------------------------------------------------------------------------
            /**
             * get global dv type list
             * @param[ out ] mGlobalDvTypes global list of dv types
             */
            const moris::Cell< moris::Cell< MSI::Dv_Type > > & get_global_dv_type_list()
            {
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
             * @param[ out ] mFlux constituitve model flux
             */
            const Matrix< DDRMat > & flux()
            {
                // if the flux was not evaluated
                if( mFluxEval )
                {
                    // evaluate the flux
                    this->eval_flux();
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
             * @param[ out ] mFluxDofDer derivative of the flux wrt dof
             */
            const Matrix< DDRMat > & dFluxdDOF( const moris::Cell< MSI::Dof_Type > aDofType )
            {
               // if aDofType is not an active dof type for the property
               MORIS_ERROR( this->check_dof_dependency( aDofType ), "Constitutive_Model::dFluxdDOF - no dependency in this dof type." );

               // get the dof index
               uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mFluxDofDerEval( tDofIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dFluxdDOF( aDofType );
               }

               // return the derivative
               return mFluxDofDer( tDofIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux derivative wrt to a dof type
             * @param[ in ] aDofTypes  a dof type wrt which the derivative is evaluated
             * @param[ in ] adFluxdDOF a matrix to fill with derivative evaluation
             */
            virtual void eval_dFluxdDOF( moris::Cell< MSI::Dof_Type > aDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dFluxdDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the strain wrt dof
             * @param[ out ] mStrainDofDer derivative of the strain wrt dof
             */
            const Matrix< DDRMat > & dStraindDOF( const moris::Cell< MSI::Dof_Type > aDofType )
            {
               // if aDofType is not an active dof type for the property
               MORIS_ERROR( this->check_dof_dependency( aDofType ), "Constitutive_Model::dStraindDOF - no dependency in this dof type." );

               // get the dof index
               uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mStrainDofDerEval( tDofIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dStraindDOF( aDofType );
               }

               // return the derivative
               return mStrainDofDer( tDofIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain derivative wrt to a dof type
             * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
             */
            virtual void eval_dStraindDOF( moris::Cell< MSI::Dof_Type > aDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dStraindDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the constituitve matrix wrt dof
             * @param[ out ] mConstDofDer derivative of the constitutive matrix wrt dof
             */
            const Matrix< DDRMat > & dConstdDOF( const moris::Cell< MSI::Dof_Type > aDofType )
            {
               // if aDofType is not an active dof type for the property
               MORIS_ERROR( this->check_dof_dependency( aDofType ), "Constitutive_Model::dConstdDOF - no dependency in this dof type." );

               // get the dof index
               uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mConstDofDerEval( tDofIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dConstdDOF( aDofType );
               }

               // return the derivative
               return mConstDofDer( tDofIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model constitutive matrix derivative wrt to a dof type
             * @param[ in ] aDofTypes    a dof type wrt which the derivative is evaluated
             */
            virtual void eval_dConstdDOF( moris::Cell< MSI::Dof_Type > aDofTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dConstdDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the flux wrt dv
             * @param[ out ] mFluxDvDer derivative of the flux wrt dv
             */
            const Matrix< DDRMat > & dFluxdDV( const moris::Cell< MSI::Dv_Type > aDvType )
            {
               // if aDvType is not an active dv type
               MORIS_ERROR( this->check_dv_dependency( aDvType ), "Constitutive_Model::dFluxdDV - no dependency in this dv type." );

               // get the dv index
               uint tDvIndex = mGlobalDvTypeMap( static_cast< uint >( aDvType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mFluxDvDerEval( tDvIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dFluxdDV( aDvType );
               }

               // return the derivative
               return mFluxDvDer( tDvIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model flux derivative wrt to a dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             */
            virtual void eval_dFluxdDV( moris::Cell< MSI::Dv_Type > aDvTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dFluxdDV - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the strain wrt dv
             * @param[ out ] mStrainDvDer derivative of the strain wrt dv
             */
            const Matrix< DDRMat > & dStraindDV( const moris::Cell< MSI::Dv_Type > aDvType )
            {
               // if aDvType is not an active dv type for the property
               MORIS_ERROR( this->check_dv_dependency( aDvType ), "Constitutive_Model::dStraindDV - no dependency in this dv type." );

               // get the dv index
               uint tDvIndex = mGlobalDvTypeMap( static_cast< uint >( aDvType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mStrainDvDerEval( tDvIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dStraindDV( aDvType );
               }

               // return the derivative
               return mStrainDvDer( tDvIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model strain derivative wrt to a dv type
             * @param[ in ] aDvTypes a dv type wrt which the derivative is evaluated
             */
            virtual void eval_dStraindDV( moris::Cell< MSI::Dv_Type > aDvTypes )
            {
                MORIS_ERROR( false, " Constitutive_Model::eval_dStraindDV - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the derivative of the constitutive matrix wrt dv
             * @param[ out ] mConstDvDer derivative of the constitutive matrix wrt dv
             */
            const Matrix< DDRMat > & dConstdDV( const moris::Cell< MSI::Dv_Type > aDvType )
            {
               // if aDvType is not an active dv type for the property
               MORIS_ERROR( this->check_dv_dependency( aDvType ), "Constitutive_Model::dConstdDV - no dependency in this dv type." );

               // get the dv index
               uint tDvIndex = mGlobalDvTypeMap( static_cast< uint >( aDvType( 0 ) ) );

               // if the derivative has not been evaluated yet
               if( mConstDvDerEval( tDvIndex ) )
               {
                   // evaluate the derivative
                   this->eval_dConstdDV( aDvType );
               }

               // return the derivative
               return mConstDvDer( tDvIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the constitutive model constitutive matrix derivative wrt to a dv type
             * @param[ in ] aDofTypes a dof type wrt which the derivative is evaluated
             */
            virtual void eval_dConstdDV( moris::Cell< MSI::Dv_Type >   aDvTypes )
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
            void eval_dFluxdDOF_FD( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                      Matrix< DDRMat >           & adFluxdDOF_FD,
                                      real                         aPerturbation )
            {
                // get the index for the considered dof type
                uint iFI = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ), 0 );

                // get number of master dofs wrt which derivative is computed
                uint tDerNumDof = mDofFI( iFI )->get_number_of_space_time_coefficients();

                // set size for derivative
                adFluxdDOF_FD.set_size( mSpaceDim, tDerNumDof, 0.0 );

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = mDofFI( iFI )->get_coeff();

                for( uint iCoeff = 0; iCoeff < tDerNumDof; iCoeff++ )
                {
                    // perturbation of the coefficent
                    Matrix< DDRMat > tCoeffPert = tCoeff;
                    tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) + aPerturbation * tCoeffPert( iCoeff );

                    // setting the perturbed coefficients
                    mDofFI( iFI )->set_coeff( tCoeffPert );

                    // reset properties
                    uint tNumProps = mPropTypes.size();
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
                    tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) - aPerturbation * tCoeffPert( iCoeff );

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
                    adFluxdDOF_FD.get_column( iCoeff ) = ( tFlux_Plus - tFlux_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeff ) );
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
            void eval_dStraindDOF_FD( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                      Matrix< DDRMat >             & adStraindDOF_FD,
                                      real                           aPerturbation )
            {
                // get the index for the considered dof type
                uint iFI = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ), 0 );

                // get number of master dofs wrt which derivative is computed
                uint tDerNumDof = mDofFI( iFI )->get_number_of_space_time_coefficients();

                // set size for derivative
                adStraindDOF_FD.set_size( mSpaceDim, tDerNumDof, 0.0 );

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = mDofFI( iFI )->get_coeff();

                for( uint iCoeff = 0; iCoeff < tDerNumDof; iCoeff++ )
                {
                    // perturbation of the coefficent
                    Matrix< DDRMat > tCoeffPert = tCoeff;
                    tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) + aPerturbation * tCoeffPert( iCoeff );

                    // setting the perturbed coefficients
                    mDofFI( iFI )->set_coeff( tCoeffPert );

                    // reset properties
                    uint tNumProps = mPropTypes.size();
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
                    tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) - aPerturbation * tCoeffPert( iCoeff );

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
                    adStraindDOF_FD.get_column( iCoeff ) = ( tStrain_Plus - tStrain_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeff ) );
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
            void eval_dFluxdDOF_FD( moris::Cell< MSI::Dv_Type >   aDvTypes,
                                      Matrix< DDRMat >          & adFluxdDV_FD,
                                      real                        aPerturbation )
            {
                // get the index for the considered dv type
                uint iFI = mGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ), 0 );

                // get number of master dvs wrt which derivative is computed
                uint tDerNumDv = mDvFI( iFI )->get_number_of_space_time_coefficients();

                // set size for derivative
                adFluxdDV_FD.set_size( mSpaceDim, tDerNumDv, 0.0 );

                // coefficients for dv type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = mDvFI( iFI )->get_coeff();

                for( uint iCoeff = 0; iCoeff < tDerNumDv; iCoeff++ )
                {
                    // perturbation of the coefficent
                    Matrix< DDRMat > tCoeffPert = tCoeff;
                    tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) + aPerturbation * tCoeffPert( iCoeff );

                    // setting the perturbed coefficients
                    mDvFI( iFI )->set_coeff( tCoeffPert );

                    // reset properties
                    uint tNumProps = mPropTypes.size();
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
                    tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) - aPerturbation * tCoeffPert( iCoeff );

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
                    adFluxdDV_FD.get_column( iCoeff ) = ( tFlux_Plus - tFlux_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeff ) );
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
            void eval_dStraindDV_FD( moris::Cell< MSI::Dv_Type >   aDvTypes,
                                     Matrix< DDRMat >             & adStraindDV_FD,
                                     real                           aPerturbation )
            {
                // get the index for the considered dof type
                uint iFI = mGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ), 0 );

                // get number of master dvs wrt which derivative is computed
                uint tDerNumDv = mDvFI( iFI )->get_number_of_space_time_coefficients();

                // set size for derivative
                adStraindDV_FD.set_size( mSpaceDim, tDerNumDv, 0.0 );

                // coefficients for dv type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = mDvFI( iFI )->get_coeff();

                for( uint iCoeff = 0; iCoeff < tDerNumDv; iCoeff++ )
                {
                    // perturbation of the coefficent
                    Matrix< DDRMat > tCoeffPert = tCoeff;
                    tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) + aPerturbation * tCoeffPert( iCoeff );

                    // setting the perturbed coefficients
                    mDvFI( iFI )->set_coeff( tCoeffPert );

                    // reset properties
                    uint tNumProps = mPropTypes.size();
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
                    tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) - aPerturbation * tCoeffPert( iCoeff );

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
                    adStraindDV_FD.get_column( iCoeff ) = ( tStrain_Plus - tStrain_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeff ) );
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

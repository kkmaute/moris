/*
 * cl_FEM_IWG.hpp
 *
 *  Created on: Oct 08, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_OPT_IWG_HPP_
#define SRC_FEM_CL_FEM_OPT_IWG_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "linalg_typedefs.hpp"              //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/MSI/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        /**
         * Integrand of Weak Form of objective and constraint functions
         */
        class Opt_IWG
        {
        protected :

            // master and slave dof type lists
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveDofTypes;

            // master and slave global dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterGlobalDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveGlobalDofTypes;

            // master and slave field interpolators
            moris::Cell< Field_Interpolator* > mMasterDofFI;
            moris::Cell< Field_Interpolator* > mSlaveDofFI;

            // master and slave dv type lists
            moris::Cell< moris::Cell< MSI::Dv_Type > > mMasterDvTypes;
            moris::Cell< moris::Cell< MSI::Dv_Type > > mSlaveDvTypes;

            // master and slave global dv type list
            moris::Cell< moris::Cell< MSI::Dv_Type > > mMasterGlobalDvTypes;
            moris::Cell< moris::Cell< MSI::Dv_Type > > mSlaveGlobalDvTypes;

            // master and slave field interpolators
            moris::Cell< Field_Interpolator* > mMasterDvFI;
            moris::Cell< Field_Interpolator* > mSlaveDvFI;

            // master and slave property type lists
            moris::Cell< fem::Property_Type > mMasterPropTypes;
            moris::Cell< fem::Property_Type > mSlavePropTypes;

            // master and slave global property type lists
            moris::Cell< fem::Property_Type > mMasterGlobalPropTypes;
            moris::Cell< fem::Property_Type > mSlaveGlobalPropTypes;

            // master and slave properties
            moris::Cell< Property* > mMasterProp;
            moris::Cell< Property* > mSlaveProp;

            // master and slave constitutive type lists
            moris::Cell< fem::Constitutive_Type > mMasterConstitutiveTypes;
            moris::Cell< fem::Constitutive_Type > mSlaveConstitutiveTypes;

            // master and slave constitutive models
            moris::Cell< Constitutive_Model* > mMasterCM;
            moris::Cell< Constitutive_Model* > mSlaveCM;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Opt_IWG(){};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            virtual ~Opt_IWG(){};

//------------------------------------------------------------------------------
            /**
             * sets dof types
             * @param[ in ] aDofTypes a cell of cell of dof types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_dof_type_list( const moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                                    mtk::Master_Slave                                   aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case( mtk::Master_Slave::MASTER ) :
                    {
                        mMasterDofTypes = aDofTypes;
                        break;
                    }
                    case( mtk::Master_Slave::SLAVE ) :
                    {
                        mSlaveDofTypes = aDofTypes;
                        break;
                    }
                    default :
                    {
                        MORIS_ERROR( false, "Opt_IWG::set_dof_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * returns a cell of dof types
             * @param[ in ] aIsMaster enum master or slave
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dof type list
                        return mMasterDofTypes;
                        break;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dof type list
                        return mSlaveDofTypes;
                        break;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Opt_IWG::get_dof_type_list - can only be master or slave." );
                        return mMasterDofTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * sets dv types
             * @param[ in ] aDvTypes  a cell of cell of dv types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_dv_type_list( const moris::Cell< moris::Cell< MSI::Dv_Type > > & aDvTypes,
                                   mtk::Master_Slave                                  aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case( mtk::Master_Slave::MASTER ) :
                    {
                        mMasterDvTypes = aDvTypes;
                        break;
                    }
                    case( mtk::Master_Slave::SLAVE ) :
                    {
                        mSlaveDvTypes = aDvTypes;
                        break;
                    }
                    default :
                    {
                        MORIS_ERROR( false, "Opt_IWG::set_dv_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * returns a cell of dv types
             * @param[ in ] aIsMaster enum master or slave
             */
            const moris::Cell< moris::Cell< MSI::Dv_Type > > & get_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dv type list
                        return mMasterDvTypes;
                        break;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dv type list
                        return mSlaveDvTypes;
                        break;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Opt_IWG::get_dv_type_list - can only be master or slave." );
                        return mMasterDvTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * set property types
             * @param[ in ] aPropertyTypes a cell of property types
             * @param[ in ] aIsMaster enum for master or slave
             *
             */
            void set_property_type_list( const moris::Cell< fem::Property_Type  > & aPropertyTypes,
                                         mtk::Master_Slave                          aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case( mtk::Master_Slave::MASTER ) :
                    {
                        mMasterPropTypes = aPropertyTypes;
                        break;
                    }
                    case( mtk::Master_Slave::SLAVE ) :
                    {
                        mSlavePropTypes = aPropertyTypes;
                        break;
                    }
                    default :
                    {
                        MORIS_ERROR( false, "Opt_IWG::set_property_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * returns a cell of property type
             * @param[ in ] aIsMaster enum master or slave
             */
            const moris::Cell< fem::Property_Type > & get_property_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global property type list
                        return mMasterPropTypes;
                        break;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global property type list
                        return mSlavePropTypes;
                        break;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Opt_IWG::get_property_type_list - can only be master or slave." );
                        return mMasterPropTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * sets constitutive types
             * @param[ in ] aConstitutiveTypes a cell of constitutive types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_constitutive_type_list( const moris::Cell< fem::Constitutive_Type  > & aConstitutiveTypes,
                                             mtk::Master_Slave                              aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case( mtk::Master_Slave::MASTER ) :
                    {
                        mMasterConstitutiveTypes = aConstitutiveTypes;
                        break;
                    }
                    case( mtk::Master_Slave::SLAVE ) :
                    {
                        mSlaveConstitutiveTypes = aConstitutiveTypes;
                        break;
                    }
                    default :
                    {
                        MORIS_ERROR( false, "Opt_IWG::set_constitutive_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * returns a cell of constitutive types
             * @param[ in ]  aIsMaster          enum master or slave
             * @param[ out ] aConstitutiveTypes a cell of constitutive types
             */
            const moris::Cell< fem::Constitutive_Type > & get_constitutive_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master constitutive type list
                        return mMasterConstitutiveTypes;
                        break;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave constitutive type list
                        return mSlaveConstitutiveTypes;
                        break;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Opt_IWG::get_constitutive_type_list - can only be master or slave." );
                        return mMasterConstitutiveTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * sets field interpolators
             * @param[ in ] aFieldInterpolators cell of field interpolator pointers
             * @param[ in ] aIsMaster           enum for master or slave
             */
            void set_dof_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators,
                                              mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER )
            {
                // check input size
                MORIS_ASSERT( aFieldInterpolators.size() == this->get_global_dof_type_list( aIsMaster ).size(),
                              "IWG::set_field_interpolators - wrong input size. " );

                // check field interpolator type
                bool tCheckFI = true;
                for( uint iFI = 0; iFI < aFieldInterpolators.size(); iFI++ )
                {
                    tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dof_type()( 0 ) == this->get_global_dof_type_list( aIsMaster )( iFI )( 0 ) );
                }
                MORIS_ASSERT( tCheckFI, "Opt_IWG::set_field_interpolators - wrong field interpolator dof type. ");

                // set field interpolators
                this->get_dof_field_interpolators( aIsMaster ) = aFieldInterpolators;
            }

//------------------------------------------------------------------------------
            /**
             * get field interpolators
             * @param[ in ]  aIsMaster           enum master or slave
             * @param[ out ] aFieldInterpolators cell of field interpolator pointers
             */
            moris::Cell< Field_Interpolator* > & get_dof_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master field interpolator pointers
                        return mMasterDofFI;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave field interpolator pointers
                        return mSlaveDofFI;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Opt_IWG::set_field_interpolators - can only be master or slave." );
                        return mMasterDofFI;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * sets dv field interpolators
             * @param[ in ] aFieldInterpolators cell of field interpolator pointers
             * @param[ in ] aIsMaster           enum for master or slave
             */
            void set_dv_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators,
                                             mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER )
            {
                // check input size
                MORIS_ASSERT( aFieldInterpolators.size() == this->get_global_dv_type_list( aIsMaster ).size(),
                              "IWG::set_dv_field_interpolators - wrong input size. " );

                // check field interpolator type
                bool tCheckFI = true;
                for( uint iFI = 0; iFI < aFieldInterpolators.size(); iFI++ )
                {
                    tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dv_type()( 0 ) == this->get_global_dv_type_list( aIsMaster )( iFI )( 0 ) );
                }
                MORIS_ASSERT( tCheckFI, "Opt_IWG::set_dv_field_interpolators - wrong field interpolator dof type. ");

                // set field interpolators
                this->get_dv_field_interpolators( aIsMaster ) = aFieldInterpolators;
            }

//------------------------------------------------------------------------------
            /**
             * get field interpolators
             * @param[ in ]  aIsMaster           enum master or slave
             * @param[ out ] aFieldInterpolators cell of field interpolator pointers
             */
            moris::Cell< Field_Interpolator* > & get_dv_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master field interpolator pointers
                        return mMasterDvFI;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave field interpolator pointers
                        return mSlaveDvFI;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Opt_IWG::set_dv_field_interpolators - can only be master or slave." );
                        return mMasterDvFI;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * checks that field interpolators were assigned
             * @param[ in ]  aIsMaster enum master or slave
             */
             void check_dv_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // check field interpolators cell size
                 MORIS_ASSERT( this->get_dv_field_interpolators( aIsMaster ).size() == this->get_dv_type_list( aIsMaster ).size(),
                               "Opt_IWG::check_dv_field_interpolators - wrong FI size. " );

                // loop over the field interpolator pointers
                for( uint iFI = 0; iFI < this->get_dv_type_list( aIsMaster ).size(); iFI++ )
                {
                    // check that the field interpolator was set
                    MORIS_ASSERT( this->get_dv_field_interpolators( aIsMaster )( iFI ) != nullptr,
                                  "Opt_IWG::check_dv_field_interpolators - FI missing. " );
                }
             }

//------------------------------------------------------------------------------
             /**
              * sets properties
              * @param[ in ] aProperties cell of property pointers
              * @param[ in ] aIsMaster   enum master or slave
              */
             void set_properties( moris::Cell< Property* > & aProperties,
                                  mtk::Master_Slave          aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // check input size
                 MORIS_ASSERT( aProperties.size() == this->get_global_property_type_list( aIsMaster ).size(),
                               "Opt_IWG::set_properties - master, wrong input size. " );

                 // check property type
                 bool tCheckProp = true;
                 for( uint iProp = 0; iProp < aProperties.size(); iProp++ )
                 {
                     tCheckProp = tCheckProp && ( aProperties( iProp )->get_property_type() == this->get_global_property_type_list( aIsMaster)( iProp ) );
                 }
                 MORIS_ASSERT( tCheckProp, "Opt_IWG::set_properties - wrong property type. ");

                 // set properties
                 this->get_properties( aIsMaster ) = aProperties;

                 // create a global dof type list
                 this->build_global_dof_type_list( aIsMaster );

                 // create a global dv type list
                 this->build_global_dv_type_list( aIsMaster );
             }

//------------------------------------------------------------------------------
             /**
              * gets properties
              * @param[ in ]  aIsMaster   enum master or slave
              * @param[ out ] aProperties cell of property pointers
              */
             moris::Cell< Property* > & get_properties( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // switch on master/slave
                 switch( aIsMaster )
                 {
                     // if master
                     case( mtk::Master_Slave::MASTER ):
                     {
                         // return master property pointers
                         return mMasterProp;
                         break;
                     }
                     // if slave
                     case( mtk::Master_Slave::SLAVE ):
                     {
                         // return slave property pointers
                         return mSlaveProp;
                         break;
                     }
                     // if none
                     default:
                     {
                         MORIS_ASSERT( false, "Opt_IWG::get_properties - can only be master or slave." );
                         return mMasterProp;
                         break;
                     }
                 }
             }

//------------------------------------------------------------------------------
             /**
              * checks that properties are assigned
              * @param[ in ] aIsMaster enum master or slave
              */
              void check_properties( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // check property cell size
                  MORIS_ASSERT( this->get_properties( aIsMaster ).size() == this->get_global_property_type_list( aIsMaster ).size(),
                                "Opt_IWG::check_properties - wrong property size. " );

                  // loop over all properties and check that they are assigned
                  for( uint iProp = 0; iProp < this->get_global_property_type_list( aIsMaster ).size(); iProp++ )
                  {
                      MORIS_ASSERT( this->get_properties( aIsMaster )( iProp ) != nullptr,
                                    "Opt_IWG::check_properties - property missing. " );
                  }
              }

//------------------------------------------------------------------------------
             /**
              * sets constitutive models
              * @param[ in ] aConstitutiveModels cell of constitutive model pointers
              * @param[ in ] aIsMaster           enum master or slave
              */
             void set_constitutive_models( moris::Cell< Constitutive_Model* > & aConstitutiveModels,
                                           mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // check input size
                 MORIS_ASSERT( aConstitutiveModels.size() == this->get_constitutive_type_list( aIsMaster ).size(),
                               "Opt_IWG::set_constitutive_models - wrong input size. " );

                 // check constitutive type
                 bool tCheckConstitutive = true;
                 for( uint iCM = 0; iCM < aConstitutiveModels.size(); iCM++ )
                 {
                     tCheckConstitutive = tCheckConstitutive
                                        && ( aConstitutiveModels( iCM )->get_constitutive_type() == this->get_constitutive_type_list( aIsMaster)( iCM ) );
                 }
                 MORIS_ASSERT( tCheckConstitutive, "IWG::set_constitutive_models - wrong constitutive type. ");

                 // set constitutive models
                 this->get_constitutive_models( aIsMaster ) = aConstitutiveModels;

                 // create a global property type list
                 this->build_global_property_type_list( aIsMaster );
             }

//------------------------------------------------------------------------------
             /**
              * gets constitutive models
              * @param[ in ]  aIsMaster           enum master or slave
              * @param[ out ] aConstitutiveModels cell of constitutive model pointers
              */
             moris::Cell< Constitutive_Model* > & get_constitutive_models( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // switch on master/slave
                 switch( aIsMaster )
                 {
                     // if master
                     case( mtk::Master_Slave::MASTER ):
                     {
                         // return master property pointers
                         return mMasterCM;
                         break;
                     }
                     // if slave
                     case( mtk::Master_Slave::SLAVE ):
                     {
                         // return slave property pointers
                         return mSlaveCM;
                         break;
                     }
                     // if none
                     default:
                     {
                         MORIS_ASSERT( false, "Opt_IWG::get_constitutive_models - can only be master or slave." );
                         return mMasterCM;
                         break;
                     }
                 }
             }

//------------------------------------------------------------------------------
             /**
              * checks that constitutive models are assigned
              * @param[ in ] aIsMaster enum master or slave
              */
              void check_constitutive_models( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // check constitutive model cell size
                  MORIS_ASSERT( this->get_constitutive_models( aIsMaster ).size() == this->get_constitutive_type_list( aIsMaster ).size(),
                                "Opt_IWG::check_constitutive_models - wrong constitutive model size. " );

                  // loop over all constitutive models and check that they are assigned
                  for( uint iCM = 0; iCM < this->get_constitutive_type_list( aIsMaster ).size(); iCM++ )
                  {
                      MORIS_ASSERT( this->get_constitutive_models( aIsMaster )( iCM ) != nullptr,
                                    "Opt_IWG::check_constitutive_models - constitutive model missing. " );
                  }
              }

//------------------------------------------------------------------------------
              /**
               * create a global dof type list including IWG and CM and properties dependencies
               * @param[ in ] aIsMaster enum master or slave
               */
              void build_global_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // set the size of the dof type list
                  uint tCounterMax = this->get_dof_type_list( aIsMaster ).size();

                  for ( Property* tProperty : this->get_properties( aIsMaster ) )
                  {
                      tCounterMax += tProperty->get_dof_type_list().size();
                  }

                  for( Constitutive_Model* tCM : this->get_constitutive_models( aIsMaster ) )
                  {
                      tCounterMax += tCM->get_global_dof_type_list().size();
                  }
                  this->get_global_dof_type_list( aIsMaster ).resize( tCounterMax );
                  moris::Cell< sint > tCheckList( tCounterMax, -1 );

                  // init total dof counter
                  uint tCounter = 0;

                  // get active dof types from IWG
                  for ( uint iDOF = 0; iDOF < this->get_dof_type_list( aIsMaster ).size(); iDOF++ )
                  {
                      tCheckList( tCounter ) = static_cast< uint >( this->get_dof_type_list( aIsMaster )( iDOF )( 0 ) );
                      this->get_global_dof_type_list( aIsMaster )( tCounter ) = this->get_dof_type_list( aIsMaster )( iDOF );
                      tCounter++;
                  }

                  // get active dof types from IWG properties
                  for ( Property* tProperty : this->get_properties( aIsMaster ) )
                  {
                      // get active dof type
                      moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

                      for ( uint iDOF = 0; iDOF < tProperty->get_dof_type_list().size(); iDOF++ )
                      {
                          // check enum is not already in the list
                          bool tCheck = false;
                          for( uint i = 0; i < tCounter; i++ )
                          {
                              tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tProperty->get_dof_type_list()( iDOF )( 0 ) ) );
                          }

                          // if dof enum not in the list
                          if ( !tCheck )
                          {
                              tCheckList( tCounter ) = static_cast< uint >( tProperty->get_dof_type_list()( iDOF )( 0 ) );
                              this->get_global_dof_type_list( aIsMaster )( tCounter ) = tProperty->get_dof_type_list()( iDOF );
                              tCounter++;
                          }
                      }
                  }

                  // get active dof types from IWG constitutive models
                  for( Constitutive_Model* tCM : this->get_constitutive_models( aIsMaster ) )
                  {
                      // get active dof type
                      moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tCM->get_global_dof_type_list();

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
                              this->get_global_dof_type_list( aIsMaster )( tCounter ) = tActiveDofType( iDOF );
                              tCounter++;
                          }
                      }
                  }

                  // get the number of unique dof type groups, i.e. the number of interpolators
                  this->get_global_dof_type_list( aIsMaster ).resize( tCounter );
              };

//------------------------------------------------------------------------------
              /**virtual
               * get global dof type list
               * @param[ in ] aIsMaster enum master or slave
               */
              moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // switch on master/slave
                  switch( aIsMaster )
                  {
                      // if master
                      case( mtk::Master_Slave::MASTER ):
                      {
                          // return master global dof type list
                          return mMasterGlobalDofTypes;
                          break;
                      }
                      // if slave
                      case( mtk::Master_Slave::SLAVE ):
                      {
                          // return slave global dof type list
                          return mSlaveGlobalDofTypes;
                          break;
                      }
                      // if none
                      default:
                      {
                          MORIS_ASSERT( false, "Opt_IWG::get_global_dof_type_list - can only be master or slave." );
                          return mMasterGlobalDofTypes;
                          break;
                      }
                  }
              };

//------------------------------------------------------------------------------
              /**
               * create a global dv type list including IWG and CM and properties dependencies
               * @param[ in ] aIsMaster enum master or slave
               */
              void build_global_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // set the size of the dof type list
                  uint tCounterMax = this->get_dv_type_list( aIsMaster ).size();

                  for ( Property* tProperty : this->get_properties( aIsMaster ) )
                  {
                      tCounterMax += tProperty->get_dv_type_list().size();
                  }

                  for( Constitutive_Model* tCM : this->get_constitutive_models( aIsMaster ) )
                  {
                      tCounterMax += tCM->get_global_dv_type_list().size();
                  }
                  this->get_global_dv_type_list( aIsMaster ).resize( tCounterMax );
                  moris::Cell< sint > tCheckList( tCounterMax, -1 );

                  // init total dof counter
                  uint tCounter = 0;

                  // get active dof types from IWG
                  for ( uint iDV = 0; iDV < this->get_dv_type_list( aIsMaster ).size(); iDV++ )
                  {
                      tCheckList( tCounter ) = static_cast< uint >( this->get_dv_type_list( aIsMaster )( iDV )( 0 ) );
                      this->get_global_dv_type_list( aIsMaster )( tCounter ) = this->get_dv_type_list( aIsMaster )( iDV );
                      tCounter++;
                  }

                  // get active dv types from IWG properties
                  for ( Property* tProperty : this->get_properties( aIsMaster ) )
                  {
                      // get active dv type
                      moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvType = tProperty->get_dv_type_list();

                      for ( uint iDV = 0; iDV < tProperty->get_dv_type_list().size(); iDV++ )
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
                              this->get_global_dv_type_list( aIsMaster )( tCounter ) = tActiveDvType( iDV );
                              tCounter++;
                          }
                      }
                  }

                  // get active dof types from IWG constitutive models
                  for( Constitutive_Model* tCM : this->get_constitutive_models( aIsMaster ) )
                  {
                      // get active dof type
                      moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvType = tCM->get_global_dv_type_list();

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
                              this->get_global_dv_type_list( aIsMaster )( tCounter ) = tActiveDvType( iDV );
                              tCounter++;
                          }
                      }
                  }

                  // get the number of unique dv type groups, i.e. the number of interpolators
                  this->get_global_dv_type_list( aIsMaster ).resize( tCounter );
              };

//------------------------------------------------------------------------------
              /**
               * gets global dv type list
               * @param[ in ] aIsMaster enum master or slave
               */
              moris::Cell< moris::Cell< MSI::Dv_Type > > & get_global_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // switch on master/slave
                  switch( aIsMaster )
                  {
                      // if master
                      case( mtk::Master_Slave::MASTER ):
                      {
                          // return master global dv type list
                          return mMasterGlobalDvTypes;
                          break;
                      }
                      // if slave
                      case( mtk::Master_Slave::SLAVE ):
                      {
                          // return slave global dv type list
                          return mSlaveGlobalDvTypes;
                          break;
                      }
                      // if none
                      default:
                      {
                          MORIS_ASSERT( false, "Opt_IWG::get_global_dv_type_list - can only be master or slave." );
                          return mMasterGlobalDvTypes;
                          break;
                      }
                  }
              };

//------------------------------------------------------------------------------
              /**
               * create a global property type list including IWG and constitutive model dependencies
               * @param[ in ] aIsMaster enum master or slave
               */
              void build_global_property_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // set the size of the dof type list
                  uint tCounterMax = this->get_property_type_list( aIsMaster ).size();

                  for ( Constitutive_Model* tCM : this->get_constitutive_models( aIsMaster ) )
                  {
                      tCounterMax += tCM->get_property_type_list().size();
                  }
                  this->get_global_property_type_list( aIsMaster ).resize( tCounterMax );
                  moris::Cell< sint > tCheckList( tCounterMax, -1 );

                  // init total property counter
                  uint tCounter = 0;

                  // get active property types from IWG
                  for ( uint iProp = 0; iProp < this->get_property_type_list( aIsMaster ).size(); iProp++ )
                  {
                      tCheckList( tCounter ) = static_cast< uint >( this->get_property_type_list( aIsMaster )( iProp ) );
                      this->get_global_property_type_list( aIsMaster )( tCounter ) = this->get_property_type_list( aIsMaster )( iProp );
                      tCounter++;
                  }

                  // get active property types from IWG constitutive models
                  for( Constitutive_Model* tCM : this->get_constitutive_models( aIsMaster ) )
                  {
                      // get active property type
                      moris::Cell< fem::Property_Type > tPropertyType = tCM->get_property_type_list();

                      for ( uint iProp = 0; iProp < tPropertyType.size(); iProp++ )
                      {
                          // check enum is not already in the list
                          bool tCheck = false;
                          for( uint i = 0; i < tCounter; i++ )
                          {
                              tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tPropertyType( iProp ) ) );
                          }

                          // if property enum not in the list
                          if ( !tCheck )
                          {
                              tCheckList( tCounter ) = static_cast< uint >( tPropertyType( iProp ) );
                              this->get_global_property_type_list( aIsMaster )( tCounter ) = tPropertyType( iProp );
                              tCounter++;
                          }
                      }
                  }

                  // resize to number of unique property types
                  this->get_global_property_type_list( aIsMaster ).resize( tCounter );
              };

//------------------------------------------------------------------------------
              /**
               * get global property type list
               * @param[ in ] aIsMaster enum master or slave
               */
              moris::Cell< fem::Property_Type > & get_global_property_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // switch on master/slave
                  switch( aIsMaster )
                  {
                      // if master
                      case( mtk::Master_Slave::MASTER ):
                      {
                          // return master global dof type list
                          return mMasterGlobalPropTypes;
                          break;
                      }
                      // if slave
                      case( mtk::Master_Slave::SLAVE ):
                      {
                          // return slave global dof type list
                          return mSlaveGlobalPropTypes;
                          break;
                      }
                      // if none
                      default:
                      {
                          MORIS_ASSERT( false, "Opt_IWG::get_global_property_type_list - can only be master or slave." );
                          return mMasterGlobalPropTypes;
                          break;
                      }
                  }
              };

//------------------------------------------------------------------------------
            /**
             * evaluates the objective or constraint function
             * @param[ in ] aOptMeasure matrix to fill with residual
             */
            virtual void compute_opt( Matrix< DDRMat > & aOptMeasure )
            {
                MORIS_ERROR( false, "Opt_IWG::compute_opt - This function does nothing. " );
            }

        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_OPT_IWG_HPP_ */

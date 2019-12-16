/*
 * cl_FEM_IWG.hpp
 *
 *  Created on: Aug 13, 2018
 *      Author: messe/noel
 */
#ifndef SRC_FEM_CL_FEM_IWG_HPP_
#define SRC_FEM_CL_FEM_IWG_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "linalg_typedefs.hpp"              //MRS/COR/src           // note: linalg_typedefs.hpp must be included AFTER the cl_Matrix.hpp
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"    //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/MSI/src

#include "fn_reshape.hpp"

namespace moris
{
    namespace fem
    {
        class Set;
        class Field_Interpolator_Manager;
//------------------------------------------------------------------------------
        /**
         * Integrand of Weak Form of Governing Equations
         */
        class IWG
        {
        protected :

            // FEM set pointer
            fem::Set * mSet = nullptr;

            // nodal weak BCs
            Matrix< DDRMat > mNodalWeakBCs;

            // normal
            Matrix< DDRMat > mNormal;

            // residual dof type
            moris::Cell< MSI::Dof_Type > mResidualDofType;

            // bool true if residual dof type requested
            bool mResidualDofTypeRequested = false;

            // master and slave dof type lists
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveDofTypes;

            // bool for building global dof type list and map
            bool mGlobalDofBuild = true;
            bool mGlobalDvBuild = true;

            // master and slave global dof type lists
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterGlobalDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveGlobalDofTypes;

            // master and slave global dof type lists
            moris::Cell< moris::Cell< MSI::Dof_Type > > mRequestedMasterGlobalDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mRequestedSlaveGlobalDofTypes;

            // master and slave dof field interpolators
            moris::Cell< Field_Interpolator* > mMasterFI;
            moris::Cell< Field_Interpolator* > mSlaveFI;

            // master and slave field interpolator managers
            Field_Interpolator_Manager * mMasterFIManager = nullptr;
            Field_Interpolator_Manager * mSlaveFIManager  = nullptr;

            // master and slave dv type lists
            moris::Cell< moris::Cell< MSI::Dv_Type > > mMasterDvTypes;
            moris::Cell< moris::Cell< MSI::Dv_Type > > mSlaveDvTypes;

            // master and slave global dv type list
            moris::Cell< moris::Cell< MSI::Dv_Type > > mMasterGlobalDvTypes;
            moris::Cell< moris::Cell< MSI::Dv_Type > > mSlaveGlobalDvTypes;

            // master and slave global dv type maps
            Matrix< DDSMat > mMasterGlobalDvTypeMap;
            Matrix< DDSMat > mSlaveGlobalDvTypeMap;

            // master and slave dv field interpolators
            moris::Cell< Field_Interpolator* > mMasterDvFI;
            moris::Cell< Field_Interpolator* > mSlaveDvFI;

            // master and slave properties
            moris::Cell< std::shared_ptr< Property > > mMasterProp;
            moris::Cell< std::shared_ptr< Property > > mSlaveProp;

            // master and slave constitutive models
            moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mMasterCM;
            moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mSlaveCM;

            // stabilization parameters
            moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > mStabilizationParam;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            IWG(){};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            virtual ~IWG(){};

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
            /*
             * set field interpolator manager
             * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
             * @param[ in ] aIsMaster                 an enum for master or slave
             */
            void set_field_interpolator_manager( Field_Interpolator_Manager * aFieldInterpolatorManager,
                                                 mtk::Master_Slave            aIsMaster = mtk::Master_Slave::MASTER );

//------------------------------------------------------------------------------
            /*
             * get field interpolator manager
             * @param[ out ] aFieldInterpolatorManager a field interpolator manager pointer
             * @param[ in ]  aIsMaster                 an enum for master or slave
             */
            Field_Interpolator_Manager * get_field_interpolator_manager( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case ( mtk::Master_Slave::MASTER ) :
                        return mMasterFIManager;

                    case ( mtk::Master_Slave::SLAVE ) :
                        return mSlaveFIManager;

                    default :
                    {
                        MORIS_ERROR( false, "IWG::get_field_inetrpolator_manager - can only be master or slave." );
                        return mMasterFIManager;
                    }
                }
            }

//------------------------------------------------------------------------------
            /*
             * free memory
             */
            void free_memory()
            {
                mResidualDofTypeRequested = false;
            }

//------------------------------------------------------------------------------
            /**
             * set nodal weak BCs
             * @param[ in ] aNodalWeakBCs matrix with nodal values
             */
            void set_nodal_weak_bcs( Matrix< DDRMat > & aNodalWeakBCs )
            {
                mNodalWeakBCs = aNodalWeakBCs;
            }

//------------------------------------------------------------------------------
            /**
             * set normal
             * @param[ in ] aNormal normal vector
             */
            void set_normal( Matrix< DDRMat > & aNormal )
            {
                mNormal = aNormal;
            }

//------------------------------------------------------------------------------
            /**
             * set residual dof type
             * @param[ in ] aResidualdofType a cell of residual dof types
             */
            void set_residual_dof_type( const moris::Cell< MSI::Dof_Type > & aResidualDofType )
            {
                mResidualDofType = aResidualDofType;
            }

//------------------------------------------------------------------------------
            /**
             * return a dof type for the residual
             * @param[ out ] aResidualdofType a cell of residual dof types
             */
            const moris::Cell< MSI::Dof_Type > & get_residual_dof_type() const
            {
                return mResidualDofType;
            };

//------------------------------------------------------------------------------
            /**
             * set IWG active dof types
             * @param[ in ] aDofTypes a list of group of dof types
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
                        MORIS_ERROR( false, "IWG::set_dof_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * return a cell of dof types active for the IWG
             * @param[ in ] aIsMaster enum master or slave
             * @param[ out ] aDofTypes a list of group of dof types
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
                        MORIS_ASSERT( false, "IWG::get_dof_type_list - can only be master or slave." );
                        return mMasterDofTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * set IWG active dv types
             * @param[ in ] aDvTypes a list of group of dv types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_dv_type_list( const moris::Cell< moris::Cell< MSI::Dv_Type > > & aDvTypes,
                                    mtk::Master_Slave                                 aIsMaster = mtk::Master_Slave::MASTER )
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
                        MORIS_ERROR( false, "IWG::set_dv_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * return a cell of dv types active for the IWG
             * @param[ in ]  aIsMaster enum master or slave
             * @param[ out ] aDvTypes a list of group of dv types
             */
            const moris::Cell< moris::Cell< MSI::Dv_Type > > & get_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dof type list
                        return mMasterDvTypes;
                        break;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dof type list
                        return mSlaveDvTypes;
                        break;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "IWG::get_dv_type_list - can only be master or slave." );
                        return mMasterDvTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * set geometry interpolator
             * @param[ in ] aGeometryInterpolator geometry interpolator pointers
             * @param[ in ] aIsMaster             enum for master or slave
             */
            void set_geometry_interpolator( Geometry_Interpolator* aGeometryInterpolator,
                                            mtk::Master_Slave      aIsMaster = mtk::Master_Slave::MASTER )
            {
                // set geometry interpolator for the SP
                for( std::shared_ptr< Stabilization_Parameter > tSP : this->get_stabilization_parameters() )
                {
                    if( tSP != nullptr )
                    {
                        tSP->set_geometry_interpolator( aGeometryInterpolator, aIsMaster );
                    }
                }

                // set geometry interpolator for constitutive models
                for( std::shared_ptr< Constitutive_Model > tCM : this->get_constitutive_models( aIsMaster ) )
                {
                    if( tCM != nullptr )
                    {
                        tCM->set_geometry_interpolator( aGeometryInterpolator );
                    }
                }

                // set geometry interpolator for properties
                for( std::shared_ptr< Property > tProp : this->get_properties( aIsMaster ) )
                {
                    if( tProp != nullptr )
                    {
                        tProp->set_geometry_interpolator( aGeometryInterpolator );
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * check that dof field interpolators were assigned
             * @param[ in ]  aIsMaster enum master or slave
             */
             void check_dof_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

//------------------------------------------------------------------------------
            /**
             * check that dv field interpolators were assigned
             * @param[ in ]  aIsMaster enum master or slave
             */
             void check_dv_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
//                 // check field interpolators cell size
//                 MORIS_ASSERT( this->get_dv_field_interpolators( aIsMaster ).size() == this->get_global_dv_type_list( aIsMaster ).size(),
//                               "IWG::check_dv_field_interpolators - wrong FI size. " );
//
//                // loop over the field interpolator pointers
//                for( uint iFI = 0; iFI < this->get_global_dv_type_list( aIsMaster ).size(); iFI++ )
//                {
//                    // check that the field interpolator was set
//                    MORIS_ASSERT( this->get_dv_field_interpolators( aIsMaster )( iFI ) != nullptr,
//                                  "IWG::check_dv_field_interpolators - FI missing. " );
//                }
             }

//------------------------------------------------------------------------------
             /**
              * set property
              * @param[ in ] aProperty       a property pointer
              * @param[ in ] aPropertyString a string describing the property
              * @param[ in ] aIsMaster       enum master or slave
              */
              virtual void set_property( std::shared_ptr< Property > aProperty,
                                         std::string                 aPropertyString,
                                         mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER )
              {
                  MORIS_ASSERT( false, "IWG::set_property - This function does nothing.");
              }

//------------------------------------------------------------------------------
             /**
              * get properties
              * @param[ in ]  aIsMaster   enum master or slave
              * @param[ out ] aProperties cell of property pointers
              */
             moris::Cell< std::shared_ptr< Property > > & get_properties( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // switch on master/slave
                 switch( aIsMaster )
                 {
                     // if master
                     case( mtk::Master_Slave::MASTER ):
                     {
                         // return master property pointers
                         return mMasterProp;
                     }
                     // if slave
                     case( mtk::Master_Slave::SLAVE ):
                     {
                         // return slave property pointers
                         return mSlaveProp;
                     }
                     // if none
                     default:
                     {
                         MORIS_ASSERT( false, "IWG::get_properties - can only be master or slave." );
                         return mMasterProp;
                     }
                 }
             }

//------------------------------------------------------------------------------
             /**
              * set constitutive model
              * @param[ in ] aConstitutiveModel  a constitutive model pointer
              * @param[ in ] aConstitutiveString a string defining the constitutive model
              * @param[ in ] aIsMaster           an enum for master or slave
              */
             virtual void set_constitutive_model( std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                                                  std::string                           aConstitutiveString,
                                                  mtk::Master_Slave                     aIsMaster = mtk::Master_Slave::MASTER )
             {
                 MORIS_ERROR( false, "IWG::set_constitutive_model - This function does nothing." );
             }

//------------------------------------------------------------------------------
             /**
              * get constitutive models
              * @param[ in ]  aIsMaster           enum master or slave
              * @param[ out ] aConstitutiveModels cell of constitutive model pointers
              */
             moris::Cell< std::shared_ptr< Constitutive_Model > > & get_constitutive_models( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // switch on master/slave
                 switch( aIsMaster )
                 {
                     // if master
                     case( mtk::Master_Slave::MASTER ):
                     {
                         // return master property pointers
                         return mMasterCM;
                     }
                     // if slave
                     case( mtk::Master_Slave::SLAVE ):
                     {
                         // return slave property pointers
                         return mSlaveCM;
                     }
                     // if none
                     default:
                     {
                         MORIS_ASSERT( false, "IWG::get_constitutive_models - can only be master or slave." );
                         return mMasterCM;
                     }
                 }
             }

//------------------------------------------------------------------------------
            /**
             * set stabilization parameter
             * @param[ in ] aStabilizationParameter a stabilization parameter pointer
             * @param[ in ] aStabilizationString    a string defining the stabilization parameter
             */
            virtual void set_stabilization_parameter( std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                                                      std::string                                aStabilizationString )
            {
                MORIS_ERROR( false, "IWG::set_stabilization_parameter - This function does nothing." );
            }

//------------------------------------------------------------------------------
             /**
              * get stabilization parameters
              * @param[ out ] mStabilizationParam cell of stabilization parameter pointers
              */
             moris::Cell< std::shared_ptr< Stabilization_Parameter > > & get_stabilization_parameters()
             {
                 // return penalty parameter pointers
                 return mStabilizationParam;
             }

//------------------------------------------------------------------------------
             /**
              * create a global dof type list including
              * IWG, property, constitutive and stabilization dependencies
              */
             void build_global_dof_type_list();

//------------------------------------------------------------------------------
             /**
              * get a non unique list of dof type including
              * IWG, property, constitutive and stabilization dependencies
              * for both master and slave
              */
             void get_non_unique_dof_types( moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
              /**
               * create a global dv type list including IWG, property and constitutive dependencies
               * @param[ in ] aIsMaster enum master or slave
               */
              void build_global_dv_type_list()
              {
                  // MASTER-------------------------------------------------------
                  // get the size of the dv type list
                  uint tCounterMax = 0;

                  // get number of dv types from IWG
                  tCounterMax += mMasterDvTypes.size();

                  // get number of dv types from properties
                  for ( std::shared_ptr< Property > tProperty : mMasterProp )
                  {
                      if( tProperty != nullptr )
                      {
                          tCounterMax += tProperty->get_dv_type_list().size();
                      }
                  }

                  // get number of dv types from constitutive models
                  for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                  {
                      if( tCM != nullptr )
                      {
                          tCounterMax += tCM->get_global_dv_type_list().size();
                      }
                  }

                  // get number of dv types from stabilization parameters
                  for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
                  {
                      if( tSP != nullptr )
                      {
                          tCounterMax += tSP->get_global_dv_type_list( mtk::Master_Slave::MASTER ).size();
                      }
                  }

                  // set size for the global dv type list
                  mMasterGlobalDvTypes.resize( tCounterMax );

                  // set a size for the checkList (used to avoid repeating a dv type)
                  moris::Cell< sint > tCheckList( tCounterMax, -1 );

                  // init total dv counter
                  uint tCounter = 0;

                  // get dv type from penalty parameter
                  for ( uint iDv = 0; iDv < mMasterDvTypes.size(); iDv++ )
                  {
                      // put the dv type in the checklist
                      tCheckList( tCounter ) = static_cast< uint >( mMasterDvTypes( iDv )( 0 ) );

                      // put the dv type in the global type list
                      mMasterGlobalDvTypes( tCounter ) = mMasterDvTypes( iDv );

                      // update the dv counter
                      tCounter++;
                  }

                  // get dv type from properties
                  for ( std::shared_ptr< Property > tProperty : mMasterProp )
                  {
                      if( tProperty != nullptr )
                      {
                          // get dv types for property
                          moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvType = tProperty->get_dv_type_list();

                          // loop on property dv type
                          for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                          {
                              // check enum is not already in the list
                              bool tCheck = false;
                              for( uint i = 0; i < tCounter; i++ )
                              {
                                  tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                              }

                              // if dv enum not in the list
                              if ( !tCheck )
                              {
                                  // put the dv type in the checklist
                                  tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                                  // put the dv type in the global type list
                                  mMasterGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                                  // update dv counter
                                  tCounter++;
                              }
                          }
                      }
                  }

                  // get dv type from constitutive models
                  for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                  {
                      if( tCM != nullptr )
                      {
                          // get dv types for constitutive model
                          moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvType = tCM->get_global_dv_type_list();

                          // loop on property dv type
                          for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                          {
                              // check enum is not already in the list
                              bool tCheck = false;
                              for( uint i = 0; i < tCounter; i++ )
                              {
                                  tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                              }

                              // if dv enum not in the list
                              if ( !tCheck )
                              {
                                  // put the dv type in the checklist
                                  tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                                  // put the dv type in the global type list
                                  mMasterGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                                  // update dv counter
                                  tCounter++;
                              }
                          }
                      }
                  }

                  // get dv type from stabilization parameters
                  for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
                  {
                      if( tSP != nullptr )
                      {
                          // get dv types for constitutive model
                          moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvType = tSP->get_global_dv_type_list( mtk::Master_Slave::MASTER );

                          // loop on property dv type
                          for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                          {
                              // check enum is not already in the list
                              bool tCheck = false;
                              for( uint i = 0; i < tCounter; i++ )
                              {
                                  tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                              }

                              // if dv enum not in the list
                              if ( !tCheck )
                              {
                                  // put the dv type in the checklist
                                  tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                                  // put the dv type in the global type list
                                  mMasterGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                                  // update dv counter
                                  tCounter++;
                              }
                          }
                      }
                  }

                  // get the number of unique dv type groups for the penalty parameter
                  mMasterGlobalDvTypes.resize( tCounter );

                  // SLAVE--------------------------------------------------------
                  // get the size of the dv type list
                  tCounterMax = 0;

                  // get number of dv types from penalty parameter
                  tCounterMax += mSlaveDvTypes.size();

                  // get number of dv types from properties
                  for ( std::shared_ptr< Property > tProperty : mSlaveProp )
                  {
                      if( tProperty != nullptr )
                      {
                          tCounterMax += tProperty->get_dv_type_list().size();
                      }
                  }

                  // get number of dv types from constitutive models
                  for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
                  {
                      if( tCM != nullptr )
                      {
                          tCounterMax += tCM->get_global_dv_type_list().size();
                      }
                  }

                  // get number of dv types from stabilization parameters
                  for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
                  {
                      if( tSP != nullptr )
                      {
                          tCounterMax += tSP->get_global_dv_type_list( mtk::Master_Slave::SLAVE ).size();
                      }
                  }

                  // set size for the global dv type list
                  mSlaveGlobalDvTypes.resize( tCounterMax );

                  // set a size for the checkList (used to avoid repeating a dv type)
                  tCheckList.resize( tCounterMax, -1 );

                  // init total dv counter
                  tCounter = 0;

                  // get dv type from penalty parameter
                  for ( uint iDv = 0; iDv < mSlaveDvTypes.size(); iDv++ )
                  {
                      // put the dv type in the checklist
                      tCheckList( tCounter ) = static_cast< uint >( mSlaveDvTypes( iDv )( 0 ) );

                      // put the dv type in the global type list
                      mSlaveGlobalDvTypes( tCounter ) = mSlaveDvTypes( iDv );

                      // update the dv counter
                      tCounter++;
                  }

                  // get dv type from properties
                  for ( std::shared_ptr< Property > tProperty : mSlaveProp )
                  {
                      if( tProperty != nullptr )
                      {
                          // get dv types for property
                          moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvType = tProperty->get_dv_type_list();

                          // loop on property dv type
                          for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                          {
                              // check enum is not already in the list
                              bool tCheck = false;
                              for( uint i = 0; i < tCounter; i++ )
                              {
                                  tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                              }

                              // if dv enum not in the list
                              if ( !tCheck )
                              {
                                  // put the dv type in the checklist
                                  tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                                  // put the dv type in the global type list
                                  mSlaveGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                                  // update dv counter
                                  tCounter++;
                              }
                          }
                      }
                  }

                  // get dv type from constitutive models
                  for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
                  {
                      if( tCM != nullptr )
                      {
                          // get dv types for constitutive model
                          moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvType = tCM->get_global_dv_type_list();

                          // loop on property dv type
                          for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                          {
                              // check enum is not already in the list
                              bool tCheck = false;
                              for( uint i = 0; i < tCounter; i++ )
                              {
                                  tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                              }

                              // if dv enum not in the list
                              if ( !tCheck )
                              {
                                  // put the dv type in the checklist
                                  tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                                  // put the dv type in the global type list
                                  mSlaveGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                                  // update dv counter
                                  tCounter++;
                              }
                          }
                      }
                  }

                  // get dv type from stabilization parameters
                  for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
                  {
                      if( tSP != nullptr )
                      {
                          // get dv types for constitutive model
                          moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvType = tSP->get_global_dv_type_list( mtk::Master_Slave::SLAVE );

                          // loop on property dv type
                          for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                          {
                              // check enum is not already in the list
                              bool tCheck = false;
                              for( uint i = 0; i < tCounter; i++ )
                              {
                                  tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                              }

                              // if dv enum not in the list
                              if ( !tCheck )
                              {
                                  // put the dv type in the checklist
                                  tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                                  // put the dv type in the global type list
                                  mSlaveGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                                  // update dv counter
                                  tCounter++;
                              }
                          }
                      }
                  }

                  // get the number of unique dv type groups for the penalty parameter
                  mSlaveGlobalDvTypes.resize( tCounter );
              };

//------------------------------------------------------------------------------
            /**
             * build global dof type map
             */
            void build_global_dv_type_map()
            {
                // MASTER-------------------------------------------------------
                // get number of global dv types
                uint tNumDvTypes = mMasterGlobalDvTypes.size();

                // determine the max Dv_Type enum
                sint tMaxEnum = 0;
                for( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
                {
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mMasterGlobalDvTypes( iDv )( 0 ) ) );
                }
                tMaxEnum++;

                // set the Dv_Type map size
                mMasterGlobalDvTypeMap.set_size( tMaxEnum, 1, -1 );

                // fill the Dv_Type map
                for( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
                {
                    // fill the property map
                    mMasterGlobalDvTypeMap( static_cast< int >( mMasterGlobalDvTypes( iDv )( 0 ) ), 0 ) = iDv;
                }

                // SLAVE-------------------------------------------------------
                // get number v global dv types
                tNumDvTypes = mSlaveGlobalDvTypes.size();

                // determine the max Dv_Type enum
                tMaxEnum = 0;
                for( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
                {
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mSlaveGlobalDvTypes( iDv )( 0 ) ) );
                }
                tMaxEnum++;

                // set the dv type map size
                mSlaveGlobalDvTypeMap.set_size( tMaxEnum, 1, -1 );

                // fill the dv type map
                for( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
                {
                    // fill the property map
                    mSlaveGlobalDvTypeMap( static_cast< int >( mSlaveGlobalDvTypes( iDv )( 0 ) ), 0 ) = iDv;
                }
            }

//------------------------------------------------------------------------------
            /**
             * get global dv type map
             */
            const Matrix< DDSMat > & get_global_dv_type_map( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dv type map
                        return mMasterGlobalDvTypeMap;
                        break;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dv type map
                        return mSlaveGlobalDvTypeMap;
                        break;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Penalty_Parameter::get_global_dv_type_map - can only be master or slave." );
                        return mMasterGlobalDvTypeMap;
                        break;
                    }
                }
            }


//------------------------------------------------------------------------------
              /**
               * get global dof type list
               * @param[ in ]  aIsMaster       enum master or slave
               * @param[ out ] mGlobalDofTypes global list of group of dof types
               */
              moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // if the global list was not yet built
                  if( mGlobalDofBuild )
                  {
                      // build the stabilization parameter global dof type list
                      this->build_global_dof_type_list();

                      // update build flag
                      mGlobalDofBuild = false;
                  }

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
                          MORIS_ASSERT( false, "IWG::get_global_dof_type_list - can only be master or slave." );
                          return mMasterGlobalDofTypes;
                          break;
                      }
                  }
              };

//------------------------------------------------------------------------------
              /**
               * get global dv type list
               * @param[ in ]  aIsMaster       enum master or slave
               * @param[ out ] mGlobalDvTypes global list of group of dv types
               */
              moris::Cell< moris::Cell< MSI::Dv_Type > > & get_global_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // if the global list was not yet built
                  if( mGlobalDvBuild )
                  {
                      // build the stabilization parameter global dof type list
                      this->build_global_dv_type_list();

                      // build the stabilization parameter global dof type map
                      this->build_global_dv_type_map();

                      // update build flag
                      mGlobalDvBuild = false;
                  }

                  // switch on master/slave
                  switch( aIsMaster )
                  {
                      // if master
                      case( mtk::Master_Slave::MASTER ):
                      {
                          // return master global dof type list
                          return mMasterGlobalDvTypes;
                          break;
                      }
                      // if slave
                      case( mtk::Master_Slave::SLAVE ):
                      {
                          // return slave global dof type list
                          return mSlaveGlobalDvTypes;
                          break;
                      }
                      // if none
                      default:
                      {
                          MORIS_ASSERT( false, "IWG::get_global_dv_type_list - can only be master or slave." );
                          return mMasterGlobalDvTypes;
                          break;
                      }
                  }
              };

//------------------------------------------------------------------------------
            /**
             * evaluate the residual
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            virtual void compute_residual( real aWStar ) = 0;

//------------------------------------------------------------------------------
            /**
             * evaluate the Jacobian
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            virtual void compute_jacobian( real aWStar ) = 0;

//------------------------------------------------------------------------------
            /**
             * evaluate the residual and the Jacobian
             * @param[ in ] aResidual matrix to fill with residual
             * @param[ in ] aJacobians cell of matrices to fill with Jacobians
             */
            virtual void compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                        moris::Cell< Matrix< DDRMat > >                & aResidual )
            {
                MORIS_ERROR( false, " IWG::compute_jacobian_and_residual - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * check the Jacobian with FD
             * @param[ in ] aPerturbation real to perturb for FD
             * @param[ in ] aEpsilon      real for check
             * @param[ in ] aWStar        real weight associated to evaluation point
             * @param[ in ] aJacobians    cell of cell of matrices to fill with Jacobians
             * @param[ in ] aJacobians_FD cell of cell of matrices to fill with Jacobians by FD
             */
            bool check_jacobian( real                                             aPerturbation,
                                 real                                             aEpsilon,
                                 real                                             aWStar,
                                 moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                 moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFD );

//------------------------------------------------------------------------------
            /**
             * check the Jacobian with FD double
             * @param[ in ] aPerturbation real to perturb for FD
             * @param[ in ] aEpsilon      real for check
             * @param[ in ] aWStar        real weight associated to evaluation point
             * @param[ in ] aJacobians    cell of cell of matrices to fill with Jacobians
             * @param[ in ] aJacobians_FD cell of cell of matrices to fill with Jacobians by FD
             */
            bool check_jacobian_double( real                                             aPerturbation,
                                        real                                             aEpsilon,
                                        real                                             aWStar,
                                        moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                        moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFDs );

////------------------------------------------------------------------------------
//            /**
//             * evaluate the derivative of the residual wrt the design variables
//             * @param[ in ] aWStar weight associated to the evaluation point
//             */
//            virtual void compute_drdpdv( real aWStar ) = 0;

//------------------------------------------------------------------------------
            /**
             * reset evaluation flags
             */
            void reset_eval_flags()
            {
                // reset properties
                for ( std::shared_ptr< Property > tProp : mMasterProp )
                {
                    if ( tProp != nullptr )
                    {
                        tProp->reset_eval_flags();
                    }
                }
                for ( std::shared_ptr< Property > tProp : mSlaveProp )
                {
                    if( tProp != nullptr )
                    {
                        tProp->reset_eval_flags();
                    }
                }

                // reset constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                {
                    if( tCM != nullptr )
                    {
                        tCM->reset_eval_flags();
                    }
                }
                for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
                {
                    if( tCM != nullptr )
                    {
                        tCM->reset_eval_flags();
                    }
                }

                // reset stabilization parameters
                for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
                {
                    if( tSP != nullptr )
                    {
                        tSP->reset_eval_flags();
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the Jacobian by finite difference
             * @param[ in ] aPerturbation real to perturb for FD
             * @param[ in ] aWStar        weight associated to evaluation point
             * @param[ in ] aJacobiansFD  cell of cell of matrices to fill with Jacobians evaluated by FD
             */
            void compute_jacobian_FD( real                                             aWStar,
                                      real                                             aPerturbation,
                                      moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFD );

//------------------------------------------------------------------------------
            /**
            * evaluate the Jacobian by finite difference
            * @param[ in ] aWStar        weight associated with evaluation point
            * @param[ in ] aPerturbation real to perturb for FD
            * @param[ in ] aJacobiansFD  cell of cell of matrices to fill with Jacobians evaluated by FD
            */
            void compute_jacobian_FD_double( real                                             aWStar,
                                             real                                             aPerturbation,
                                             moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFD );

//------------------------------------------------------------------------------
            /**
             * build a list of dof types requested by the solver and owned by the IWG
             * @param[ in ] aItResidual bool true if ???
             */
            void build_requested_dof_type_list( const bool aItResidual );

//------------------------------------------------------------------------------
//
//            virtual real compute_integration_error( const Matrix< DDRMat > & aNodalDOF,
//                                                    real (*aFunction)( const Matrix< DDRMat > & aPoint ) ,
//                                                    const uint        & aPointIndex )
//            {
//                MORIS_ERROR( false, "This function does nothing" );
//                return 0.0;
//            }

//------------------------------------------------------------------------------
//
//            real interpolate_scalar_at_point( const Matrix< DDRMat > & aNodalWeakBC,
//                                              const uint             & aPointIndex )
//            {
//                MORIS_ERROR( false, "This function does nothing" );
//                return 0.0;
//            }


        };
//------------------------------------------------------------------------------


    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IWG_HPP_ */

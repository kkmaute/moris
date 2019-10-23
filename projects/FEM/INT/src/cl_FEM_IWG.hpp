/*
 * cl_FEM_IWG.hpp
 *
 *  Created on: Aug 13, 2018
 *      Author: messe/noel
 */
#ifndef SRC_FEM_CL_FEM_IWG_HPP_
#define SRC_FEM_CL_FEM_IWG_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "linalg_typedefs.hpp"              //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/MSI/src

#include "fn_reshape.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        /**
         * Integrand of Weak Form of Governing Equations
         */
        class IWG
        {
        protected :

            // nodal weak BCs
            Matrix< DDRMat > mNodalWeakBCs;

            // normal
            Matrix< DDRMat > mNormal;

            // residual dof type
            moris::Cell< MSI::Dof_Type > mResidualDofType;

            // master and slave dof type lists
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveDofTypes;

            // master and slave global dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterGlobalDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveGlobalDofTypes;

            // master and slave dof field interpolators
            moris::Cell< Field_Interpolator* > mMasterFI;
            moris::Cell< Field_Interpolator* > mSlaveFI;

            // master and slave dv type lists
            moris::Cell< moris::Cell< MSI::Dv_Type > > mMasterDvTypes;
            moris::Cell< moris::Cell< MSI::Dv_Type > > mSlaveDvTypes;

            // master and slave global dv type list
            moris::Cell< moris::Cell< MSI::Dv_Type > > mMasterGlobalDvTypes;
            moris::Cell< moris::Cell< MSI::Dv_Type > > mSlaveGlobalDvTypes;

            // master and slave dv field interpolators
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
            IWG(){};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            virtual ~IWG(){};

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
             * set IWG active property types
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
                        MORIS_ERROR( false, "IWG::set_property_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * return a cell of property type active for the IWG
             * @param[ in ]  aIsMaster  enum master or slave
             * @param[ out ] mPropTypes a cell of property types
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
                        MORIS_ASSERT( false, "IWG::get_property_type_list - can only be master or slave." );
                        return mMasterPropTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * set IWG constitutive types
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
                        MORIS_ERROR( false, "IWG::set_constitutive_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * return a cell of constitutive types
             * @param[ in ]  aIsMaster enum master or slave
             * @param[ out ] mConstitutiveTypes a cell of constitutive types
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
                        MORIS_ASSERT( false, "IWG::get_constitutive_type_list - can only be master or slave." );
                        return mMasterConstitutiveTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * set dof field interpolators
             * @param[ in ] aFieldInterpolators cell of dof field interpolator pointers
             * @param[ in ] aIsMaster           enum for master or slave
             */
            void set_dof_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators,
                                              mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER )
            {
                // get input size
                uint tInputNumFI = aFieldInterpolators.size();

                // check input size
                MORIS_ASSERT( tInputNumFI == this->get_global_dof_type_list( aIsMaster ).size(),
                              "IWG::set_dof_field_interpolators - wrong input size. " );

                // check dof field interpolator type
                bool tCheckFI = true;
                for( uint iFI = 0; iFI < tInputNumFI; iFI++ )
                {
                    tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dof_type()( 0 ) == this->get_global_dof_type_list( aIsMaster )( iFI )( 0 ) );
                }
                MORIS_ASSERT( tCheckFI, "IWG::set_dof_field_interpolators - wrong field interpolator dof type. ");

                // set field interpolators
                this->get_dof_field_interpolators( aIsMaster ) = aFieldInterpolators;
            }

//------------------------------------------------------------------------------
            /**
             * get dof field interpolators
             * @param[ in ]  aIsMaster           enum master or slave
             * @param[ out ] aFieldInterpolators cell of dof field interpolator pointers
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
                        return mMasterFI;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave field interpolator pointers
                        return mSlaveFI;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "IWG::set_dof_field_interpolators - can only be master or slave." );
                        return mMasterFI;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * check that dof field interpolators were assigned
             * @param[ in ]  aIsMaster enum master or slave
             */
             void check_dof_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // check field interpolators cell size
                 MORIS_ASSERT( this->get_dof_field_interpolators( aIsMaster ).size() == this->get_global_dof_type_list( aIsMaster ).size(),
                               "IWG::check_dof_field_interpolators - wrong FI size. " );

                // loop over the field interpolator pointers
                for( uint iFI = 0; iFI < this->get_global_dof_type_list( aIsMaster ).size(); iFI++ )
                {
                    // check that the field interpolator was set
                    MORIS_ASSERT( this->get_dof_field_interpolators( aIsMaster )( iFI ) != nullptr,
                                  "IWG::check_dof_field_interpolators - FI missing. " );
                }
             }

//------------------------------------------------------------------------------
            /**
             * set dv field interpolators
             * @param[ in ] aFieldInterpolators cell of dv field interpolator pointers
             * @param[ in ] aIsMaster           enum for master or slave
             */
            void set_dv_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators,
                                              mtk::Master_Slave                   aIsMaster = mtk::Master_Slave::MASTER )
            {
                // get input size
                uint tInputNumFI = aFieldInterpolators.size();

                // check input size
                MORIS_ASSERT( tInputNumFI == this->get_global_dv_type_list( aIsMaster ).size(),
                              "IWG::set_dv_field_interpolators - wrong input size. " );

                // check dv field interpolator type
                bool tCheckFI = true;
                for( uint iFI = 0; iFI < tInputNumFI; iFI++ )
                {
                    tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dv_type()( 0 ) == this->get_global_dv_type_list( aIsMaster )( iFI )( 0 ) );
                }
                MORIS_ASSERT( tCheckFI, "IWG::set_dv_field_interpolators - wrong field interpolator dv type. ");

                // set field interpolators
                this->get_dv_field_interpolators( aIsMaster ) = aFieldInterpolators;
            }

//------------------------------------------------------------------------------
            /**
             * get dv field interpolators
             * @param[ in ]  aIsMaster           enum master or slave
             * @param[ out ] aFieldInterpolators cell of dv field interpolator pointers
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
                        MORIS_ASSERT( false, "IWG::set_dv_field_interpolators - can only be master or slave." );
                        return mMasterDvFI;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * check that dv field interpolators were assigned
             * @param[ in ]  aIsMaster enum master or slave
             */
             void check_dv_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // check field interpolators cell size
                 MORIS_ASSERT( this->get_dv_field_interpolators( aIsMaster ).size() == this->get_global_dv_type_list( aIsMaster ).size(),
                               "IWG::check_dof_field_interpolators - wrong FI size. " );

                // loop over the field interpolator pointers
                for( uint iFI = 0; iFI < this->get_global_dv_type_list( aIsMaster ).size(); iFI++ )
                {
                    // check that the field interpolator was set
                    MORIS_ASSERT( this->get_dv_field_interpolators( aIsMaster )( iFI ) != nullptr,
                                  "IWG::check_dv_field_interpolators - FI missing. " );
                }
             }

//------------------------------------------------------------------------------
             /**
              * set properties
              * @param[ in ] aProperties cell of property pointers
              * @param[ in ]  aIsMaster   enum master or slave
              */
             void set_properties( moris::Cell< Property* > & aProperties,
                                  mtk::Master_Slave          aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // check input size
                 MORIS_ASSERT( aProperties.size() == this->get_global_property_type_list( aIsMaster ).size(),
                               "IWG::set_properties - master, wrong input size. " );

                 // check property type
                 bool tCheckProp = true;
                 for( uint iProp = 0; iProp < aProperties.size(); iProp++ )
                 {
                     tCheckProp = tCheckProp && ( aProperties( iProp )->get_property_type() == this->get_global_property_type_list( aIsMaster)( iProp ) );
                 }
                 MORIS_ASSERT( tCheckProp, "IWG::set_properties - wrong property type. ");

                 // set properties
                 this->get_properties( aIsMaster ) = aProperties;

                 // create a global dof type list
                 this->build_global_dof_type_list( aIsMaster );

                 // create a global dv type list
                 this->build_global_dv_type_list( aIsMaster );
             }

//------------------------------------------------------------------------------
             /**
              * get properties
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
                         MORIS_ASSERT( false, "IWG::get_properties - can only be master or slave." );
                         return mMasterProp;
                         break;
                     }
                 }
             }

//------------------------------------------------------------------------------
             /**
              * check that properties are assigned
              * @param[ in ] aIsMaster enum master or slave
              */
              void check_properties( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // check property cell size
                  MORIS_ASSERT( this->get_properties( aIsMaster ).size() == this->get_global_property_type_list( aIsMaster ).size(),
                                "IWG::check_properties - wrong property size. " );

                  // loop over all properties and check that they are assigned
                  for( uint iProp = 0; iProp < this->get_global_property_type_list( aIsMaster ).size(); iProp++ )
                  {
                      MORIS_ASSERT( this->get_properties( aIsMaster )( iProp ) != nullptr,
                                    "IWG::check_properties - property missing. " );
                  }
              }

//------------------------------------------------------------------------------
             /**
              * set constitutive models
              * @param[ in ] aConstitutiveModels cell of constitutive model pointers
              * @param[ in ] aIsMaster           enum master or slave
              */
             void set_constitutive_models( moris::Cell< Constitutive_Model* > & aConstitutiveModels,
                                           mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // check input size
                 MORIS_ASSERT( aConstitutiveModels.size() == this->get_constitutive_type_list( aIsMaster ).size(),
                               "IWG::set_constitutive_models - wrong input size. " );

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
              * get constitutive models
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
                         MORIS_ASSERT( false, "IWG::get_constitutive_models - can only be master or slave." );
                         return mMasterCM;
                         break;
                     }
                 }
             }

//------------------------------------------------------------------------------
             /**
              * check that constitutive models are assigned
              * @param[ in ] aIsMaster enum master or slave
              */
              void check_constitutive_models( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // check constitutive model cell size
                  MORIS_ASSERT( this->get_constitutive_models( aIsMaster ).size() == this->get_constitutive_type_list( aIsMaster ).size(),
                                "IWG::check_constitutive_models - wrong constitutive model size. " );

                  // loop over all constitutive models and check that they are assigned
                  for( uint iCM = 0; iCM < this->get_constitutive_type_list( aIsMaster ).size(); iCM++ )
                  {
                      MORIS_ASSERT( this->get_constitutive_models( aIsMaster )( iCM ) != nullptr,
                                    "IWG::check_constitutive_models - constitutive model missing. " );
                  }
              }

//------------------------------------------------------------------------------
              /**
               * create a global dof type list including IWG, property and constitutive dependencies
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
              /**
               * get global dof type list
               * @param[ in ]  aIsMaster       enum master or slave
               * @param[ out ] mGlobalDofTypes global list of group of dof types
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
                          MORIS_ASSERT( false, "IWG::get_global_dof_type_list - can only be master or slave." );
                          return mMasterGlobalDofTypes;
                          break;
                      }
                  }
              };

//------------------------------------------------------------------------------
              /**
               * create a global dv type list including IWG, property and constitutive dependencies
               * @param[ in ] aIsMaster enum master or slave
               */
              void build_global_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // set the size of the dv type list
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

                  // get active dof types from IWG properties
                  for ( Property* tProperty : this->get_properties( aIsMaster ) )
                  {
                      // get active dof type
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
               * get global dv type list
               * @param[ in ]  aIsMaster       enum master or slave
               * @param[ out ] mGlobalDvTypes global list of group of dv types
               */
              moris::Cell< moris::Cell< MSI::Dv_Type > > & get_global_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
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
               * create a global property type list including IWG, constitutive model dependencies
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
               * @param[ in ] aIsMaster         enum master or slave
               * @param[ out ] mGlobalPropTypes global list of property type
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
                          MORIS_ASSERT( false, "IWG::get_global_property_type_list - can only be master or slave." );
                          return mMasterGlobalPropTypes;
                          break;
                      }
                  }
              };

//------------------------------------------------------------------------------
            /**
             * evaluate the residual
             * @param[ in ] aResidual matrix to fill with residual
             */
            virtual void compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
            {
                MORIS_ERROR( false, "IWG::compute_residual - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the Jacobian
             * @param[ in ] aJacobians cell of matrices to fill with Jacobians
             */
            virtual void compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
            {
                MORIS_ERROR( false, "IWG::compute_jacobian - This function does nothing. " );
            }

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
             * @param[ in ] aJacobians    cell of cell of matrices to fill with Jacobians
             * @param[ in ] aJacobians_FD cell of cell of matrices to fill with Jacobians by FD
             */
            bool check_jacobian( real                                             aPerturbation,
                                 real                                             aEpsilon,
                                 moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                 moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians_FD )
            {
                // set jacobian size
                this->set_jacobian( aJacobians );

                // compute jacobian with IWG
                this->compute_jacobian( aJacobians );

                // set jacobian size
                this->set_jacobian( aJacobians_FD );

                // compute jacobian by FD
                this->compute_jacobian_FD( aJacobians_FD, aPerturbation );

                //define a boolean for check
                bool tCheckJacobian = true;

                // check each components
                for ( uint iJac = 0; iJac < aJacobians.size(); iJac++ )
                {
                    for( uint jJac = 0; jJac < aJacobians( iJac ).size(); jJac++ )
                    {
                        for( uint iiJac = 0; iiJac < aJacobians( iJac )( jJac ).n_rows(); iiJac++ )
                        {
                            for( uint jjJac = 0; jjJac < aJacobians( iJac )( jJac ).n_cols(); jjJac++ )
                            {
                                tCheckJacobian = tCheckJacobian && ( aJacobians( iJac )( jJac )( iiJac, jjJac ) - aJacobians_FD( iJac )( jJac )( iiJac, jjJac ) < aEpsilon );
                            }
                        }
                    }
                }

                // return bool
                return tCheckJacobian;
            }

//------------------------------------------------------------------------------
            /**
             * check the Jacobian with FD double
             * @param[ in ] aPerturbation real to perturb for FD
             * @param[ in ] aEpsilon      real for check
             * @param[ in ] aJacobians    cell of cell of matrices to fill with Jacobians
             * @param[ in ] aJacobians_FD cell of cell of matrices to fill with Jacobians by FD
             */
            bool check_jacobian_double( real aPerturbation,
                                        real aEpsilon,
                                        moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                        moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians_FD )
            {
                // set jacobian size
                this->set_jacobian_double( aJacobians );

                // compute jacobian with IWG
                this->compute_jacobian( aJacobians );

                // set jacobian size
                this->set_jacobian_double( aJacobians_FD );

                // compute jacobian by FD
                this->compute_jacobian_FD_double( aJacobians_FD, aPerturbation );

                //define a boolean for check
                bool tCheckJacobian = true;

                // check each components
                for ( uint iJac = 0; iJac < aJacobians.size(); iJac++ )
                {
                    for( uint jJac = 0; jJac < aJacobians( iJac ).size(); jJac++ )
                    {
                        for( uint iiJac = 0; iiJac < aJacobians( iJac )( jJac ).n_rows(); iiJac++ )
                        {
                            for( uint jjJac = 0; jjJac < aJacobians( iJac )( jJac ).n_cols(); jjJac++ )
                            {
                                tCheckJacobian = tCheckJacobian && ( aJacobians( iJac )( jJac )( iiJac, jjJac ) - aJacobians_FD( iJac )( jJac )( iiJac, jjJac ) < aEpsilon );
                            }
                        }
                    }
                }

                // return bool
                return tCheckJacobian;
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the Jacobian by finite difference
             * @param[ in ] aJacobiansFD  cell of cell of matrices to fill with Jacobians evaluated by FD
             * @param[ in ] aPerturbation real to perturb for FD
             */
            virtual void compute_jacobian_FD( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFD,
                                              real                                             aPerturbation )
            {
                // get master number of dof types
                uint tNumDofType = mMasterGlobalDofTypes.size();

                // get master number of property types
                uint tNumProps = mMasterGlobalPropTypes.size();

                // get master number of constitutive types
                uint tNumCMs = mMasterConstitutiveTypes.size();

                // set the jacobian size
                this->set_jacobian( aJacobiansFD );

                // loop over the IWG dof types
                for( uint iFI = 0; iFI < tNumDofType; iFI++ )
                {
                    // get number of master FI bases and fields
                    uint tDerNumBases  = mMasterFI( iFI )->get_number_of_space_time_bases();
                    uint tDerNumFields = mMasterFI( iFI )->get_number_of_fields();

                    // coefficients for dof type wrt which derivative is computed
                    Matrix< DDRMat > tCoeff = mMasterFI( iFI )->get_coeff();

                    // init dof counter
                    uint tDofCounter = 0;

                    // loop over the coefficient column
                    for( uint iCoeffCol = 0; iCoeffCol< tDerNumFields; iCoeffCol++ )
                    {
                        // loop over the coefficient row
                        for( uint iCoeffRow = 0; iCoeffRow< tDerNumBases; iCoeffRow++  )
                        {
                            // perturbation of the coefficent
                            Matrix< DDRMat > tCoeffPert = tCoeff;
                            tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                            // setting the perturbed coefficients
                            mMasterFI( iFI )->set_coeff( tCoeffPert );

                            // reset properties
                            for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                            {
                                mMasterProp( iProp )->reset_eval_flags();
                            }

                            // reset constitutive model
                            for ( uint iCM = 0; iCM < tNumCMs; iCM++ )
                            {
                                mMasterCM( iCM )->reset_eval_flags();
                            }

                            // evaluate the residual
                            moris::Cell< Matrix< DDRMat > > tResidual_Plus;
                            this->set_residual( tResidual_Plus );
                            this->compute_residual( tResidual_Plus );

                            // perturbation of the coefficent
                            tCoeffPert = tCoeff;
                            tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                            // setting the perturbed coefficients
                            mMasterFI( iFI )->set_coeff( tCoeffPert );

                            // reset properties
                            for ( uint iProp = 0; iProp < tNumProps; iProp++ )
                            {
                                mMasterProp( iProp )->reset_eval_flags();
                            }
                            // reset constitutive model
                            for ( uint iCM = 0; iCM < tNumCMs; iCM++ )
                            {
                                mMasterCM( iCM )->reset_eval_flags();
                            }

                            // evaluate the residual
                            moris::Cell< Matrix< DDRMat > > tResidual_Minus;
                            this->set_residual( tResidual_Minus );
                            this->compute_residual( tResidual_Minus );

                            // evaluate Jacobian
                            aJacobiansFD( 0 )( iFI ).get_column( tDofCounter ) = ( tResidual_Plus( 0 ) - tResidual_Minus( 0 ) ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                            // update dof counter
                            tDofCounter++;
                        }
                    }
                    // reset the coefficients values
                    mMasterFI( iFI )->set_coeff( tCoeff );
                }
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the Jacobian by finite difference double
             * @param[ in ] aJacobiansFD  cell of cell of matrices to fill with Jacobians evaluated by FD
             * @param[ in ] aPerturbation real to perturb for FD
             * @param[ in ] aIsMaster     enum master or slave
             */
            virtual void compute_jacobian_FD_double( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFD,
                                                     real                                             aPerturbation )
            {
                // get master and slave number of dof types
                uint tMasterNumDofType = mMasterGlobalDofTypes.size();
                uint tSlaveNumDofType  = mSlaveGlobalDofTypes.size();

                // get master and slave number of property type
                uint tMasterNumProps = mMasterGlobalPropTypes.size();
                uint tSlaveNumProps  = mSlaveGlobalPropTypes.size();

                // get master and slave number of constitutive types
                uint tMasterNumCMs = mMasterConstitutiveTypes.size();
                uint tSlaveNumCMs = mSlaveConstitutiveTypes.size();

                // set the jacobian size
                this->set_jacobian_double( aJacobiansFD );

                // loop over the master dof types
                for( uint iFI = 0; iFI < tMasterNumDofType; iFI++ )
                {
                    // get number of master FI bases and fields
                    uint tDerNumBases  = mMasterFI( iFI )->get_number_of_space_time_bases();
                    uint tDerNumFields = mMasterFI( iFI )->get_number_of_fields();

                    // coefficients for dof type wrt which derivative is computed
                    Matrix< DDRMat > tCoeff = mMasterFI( iFI )->get_coeff();

                    // init dof counter
                    uint tDofCounter = 0;

                    // loop over the coefficients column
                    for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                    {
                        // loop over the coefficients row
                        for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                        {
                            // perturbation of the coefficent
                            Matrix< DDRMat > tCoeffPert = tCoeff;
                            tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                            // setting the perturbed coefficients
                            mMasterFI( iFI )->set_coeff( tCoeffPert );

                            // reset properties
                            for ( uint iProp = 0; iProp < tMasterNumProps; iProp++ )
                            {
                                mMasterProp( iProp )->reset_eval_flags();
                            }
                            for ( uint iProp = 0; iProp < tSlaveNumProps; iProp++ )
                            {
                                mSlaveProp( iProp )->reset_eval_flags();
                            }

                            // reset constitutive model
                            for ( uint iCM = 0; iCM < tMasterNumCMs; iCM++ )
                            {
                                mMasterCM( iCM )->reset_eval_flags();
                            }
                            for ( uint iCM = 0; iCM < tSlaveNumCMs; iCM++ )
                            {
                                mSlaveCM( iCM )->reset_eval_flags();
                            }

                            // evaluate the residual
                            moris::Cell< Matrix< DDRMat > > tResidual_Plus;
                            this->set_residual_double( tResidual_Plus );
                            this->compute_residual( tResidual_Plus );

                            // perturbation of the coefficent
                            tCoeffPert = tCoeff;

                            tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                            // setting the perturbed coefficients
                            mMasterFI( iFI )->set_coeff( tCoeffPert );

                            // reset properties
                            for ( uint iProp = 0; iProp < tMasterNumProps; iProp++ )
                            {
                                mMasterProp( iProp )->reset_eval_flags();
                            }
                            for ( uint iProp = 0; iProp < tSlaveNumProps; iProp++ )
                            {
                                mSlaveProp( iProp )->reset_eval_flags();
                            }

                            // reset constitutive model
                            for ( uint iCM = 0; iCM < tMasterNumCMs; iCM++ )
                            {
                                mMasterCM( iCM )->reset_eval_flags();
                            }
                            for ( uint iCM = 0; iCM < tSlaveNumCMs; iCM++ )
                            {
                                mSlaveCM( iCM )->reset_eval_flags();
                            }

                            // evaluate the residual
                            moris::Cell< Matrix< DDRMat > > tResidual_Minus;
                            this->set_residual_double( tResidual_Minus );
                            this->compute_residual( tResidual_Minus );

                            // evaluate Jacobian
                            aJacobiansFD( 0 )( iFI ).get_column( tDofCounter ) = ( tResidual_Plus( 0 ) - tResidual_Minus( 0 ) )/ ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );
                            aJacobiansFD( 1 )( iFI ).get_column( tDofCounter ) = ( tResidual_Plus( 1 ) - tResidual_Minus( 1 ) )/ ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                            // update dof counter
                            tDofCounter++;
                        }
                    }
                    // reset the coefficients values
                    mMasterFI( iFI )->set_coeff( tCoeff );
                }

                // loop over the slave dof types
                for( uint iFI = 0; iFI < tSlaveNumDofType; iFI++ )
                {
                    // get number of master FI bases and fields
                    uint tDerNumBases  = mSlaveFI( iFI )->get_number_of_space_time_bases();
                    uint tDerNumFields = mSlaveFI( iFI )->get_number_of_fields();

                    // coefficients for dof type wrt which derivative is computed
                    Matrix< DDRMat > tCoeff = mSlaveFI( iFI )->get_coeff();

                    // init dof counter
                    uint tDofCounter = 0;

                    // loop over the coefficients columns
                    for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                    {
                        // loop over the coefficients rows
                        for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                        {
                            // perturbation of the coefficent
                            Matrix< DDRMat > tCoeffPert = tCoeff;
                            tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                            // setting the perturbed coefficients
                            mSlaveFI( iFI )->set_coeff( tCoeffPert );

                            // reset properties
                            for ( uint iProp = 0; iProp < tMasterNumProps; iProp++ )
                            {
                                mMasterProp( iProp )->reset_eval_flags();
                            }
                            for ( uint iProp = 0; iProp < tSlaveNumProps; iProp++ )
                            {
                                mSlaveProp( iProp )->reset_eval_flags();
                            }

                            // evaluate the residual
                            moris::Cell< Matrix< DDRMat > > tResidual_Plus;
                            this->set_residual_double( tResidual_Plus );
                            this->compute_residual( tResidual_Plus );

                            // perturbation of the coefficent
                            tCoeffPert = tCoeff;
                            tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                            // setting the perturbed coefficients
                            mSlaveFI( iFI )->set_coeff( tCoeffPert );

                            // reset properties
                            for ( uint iProp = 0; iProp < tMasterNumProps; iProp++ )
                            {
                                mMasterProp( iProp )->reset_eval_flags();
                            }
                            for ( uint iProp = 0; iProp < tSlaveNumProps; iProp++ )
                            {
                                mSlaveProp( iProp )->reset_eval_flags();
                            }

                            // evaluate the residual
                            moris::Cell< Matrix< DDRMat > > tResidual_Minus;
                            this->set_residual_double( tResidual_Minus );
                            this->compute_residual( tResidual_Minus );

                            // evaluate Jacobian
                            aJacobiansFD( 0 )( tMasterNumDofType + iFI ).get_column( tDofCounter ) = ( tResidual_Plus( 0 ) - tResidual_Minus( 0 ) ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );
                            aJacobiansFD( 1 )( tMasterNumDofType + iFI ).get_column( tDofCounter ) = ( tResidual_Plus( 1 ) - tResidual_Minus( 1 ) ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                            // update dof counter
                            tDofCounter++;
                        }
                    }
                    // reset the coefficients values
                    mSlaveFI( iFI )->set_coeff( tCoeff );
                }
            }

//------------------------------------------------------------------------------
            /**
             * set the residual size
             * @param[ in ] aResidual matrix to fill with residual
             */
            void set_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
            {
                // set the size of the residual cell
                aResidual.resize( 1 );

                // set size of each residual matrix
                aResidual( 0 ).set_size( mMasterFI( 0 )->get_number_of_space_time_coefficients(), 1, 0.0 );
            }

//------------------------------------------------------------------------------
            /**
             * set the residual size double
             * @param[ in ] aResidual matrix to fill with residual
             */
            void set_residual_double( moris::Cell< Matrix< DDRMat > > & aResidual )
            {
                // set the size of the residual cell
                aResidual.resize( 2 );

                // set size of each residual matrix
                aResidual( 0 ).set_size( mMasterFI( 0 )->get_number_of_space_time_coefficients(), 1, 0.0 );
                aResidual( 1 ).set_size( mSlaveFI( 0 )->get_number_of_space_time_coefficients(), 1, 0.0 );
            }

//------------------------------------------------------------------------------
            /**
             * set the Jacobian size
             * @param[ in ] aJacobians cell of matrices to fill with Jacobians
             */
            void set_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
            {
                // get number of dof types for the IWG
                uint tNumDofType = this->get_global_dof_type_list().size();

                // set the size of the jacobian cell
                aJacobians.resize( 1 );
                aJacobians( 0 ).resize( tNumDofType );

                // get residual dof type number of dofs
                uint tNumResDof = mMasterFI( 0 )->get_number_of_space_time_coefficients();

                // loop over the master dof dependencies
                for( uint iDOF = 0; iDOF < tNumDofType; iDOF++ )
                {
                    // set size for each jacobian matrix
                    aJacobians( 0 )( iDOF ).set_size( tNumResDof, mMasterFI( iDOF )->get_number_of_space_time_coefficients(), 0.0 );
                }
            }

//------------------------------------------------------------------------------
            /**
              * set the Jacobian size double
              * @param[ in ] aJacobians cell of matrices to fill with Jacobians
              */
             void set_jacobian_double( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
             {
                 // get number of dof types for the IWG
                 uint tMasterNumDofType = this->get_global_dof_type_list().size();
                 uint tSlaveNumDofType  = this->get_global_dof_type_list( mtk::Master_Slave::SLAVE ).size();
                 uint tNumDofType = tMasterNumDofType + tSlaveNumDofType;

                 // get master residual number of dofs
                 uint tMasterNumResDof = mMasterFI( 0 )->get_number_of_space_time_coefficients();

                 // get slave residual number of dofs
                 uint tSlaveNumResDof = mSlaveFI( 0 )->get_number_of_space_time_coefficients();

                 // set the size of the jacobian cell
                 aJacobians.resize( 2 );
                 aJacobians( 0 ).resize( tNumDofType );
                 aJacobians( 1 ).resize( tNumDofType );

                 // loop over the master dof dependencies
                 for( uint iMasterDOF = 0; iMasterDOF < tMasterNumDofType; iMasterDOF++ )
                 {
                     // set size for each residual matrix
                     aJacobians( 0 )( iMasterDOF ).set_size( tMasterNumResDof, mMasterFI( iMasterDOF )->get_number_of_space_time_coefficients(), 0.0 );
                     aJacobians( 1 )( iMasterDOF ).set_size( tSlaveNumResDof,  mMasterFI( iMasterDOF )->get_number_of_space_time_coefficients(), 0.0 );
                 }

                 // loop over the slave dof dependencies
                 for( uint iSlaveDOF = 0; iSlaveDOF < tSlaveNumDofType; iSlaveDOF++ )
                 {
                     // set size for each residual matrix
                     aJacobians( 0 )( tMasterNumDofType + iSlaveDOF ).set_size( tMasterNumResDof, mSlaveFI( iSlaveDOF )->get_number_of_space_time_coefficients(), 0.0 );
                     aJacobians( 1 )( tMasterNumDofType + iSlaveDOF ).set_size( tSlaveNumResDof,  mSlaveFI( iSlaveDOF )->get_number_of_space_time_coefficients(), 0.0 );
                 }
             }

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

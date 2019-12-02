/*
 * cl_FEM_IQI.hpp
 *
 *  Created on: Oct 08, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_IQI_HPP_
#define SRC_FEM_CL_FEM_IQI_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "linalg_typedefs.hpp"              //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LNA/src

#include "cl_FEM_Field_Interpolator.hpp"      //FEM/INT/src
#include "cl_FEM_Property.hpp"                //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"      //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp" //FEM/INT/src

#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/MSI/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        /**
         * Integrand of a quantity of interest
         */
        class IQI
        {
        protected :

            // master and slave dof type lists
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveDofTypes;

            // master and slave global dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterGlobalDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveGlobalDofTypes;

            // master and slave field interpolators
            moris::Cell< fem::Field_Interpolator* > mMasterDofFI;
            moris::Cell< fem::Field_Interpolator* > mSlaveDofFI;

            // master and slave dv type lists
            moris::Cell< moris::Cell< MSI::Dv_Type > > mMasterDvTypes;
            moris::Cell< moris::Cell< MSI::Dv_Type > > mSlaveDvTypes;

            // master and slave global dv type list
            moris::Cell< moris::Cell< MSI::Dv_Type > > mMasterGlobalDvTypes;
            moris::Cell< moris::Cell< MSI::Dv_Type > > mSlaveGlobalDvTypes;

            // master and slave field interpolators
            moris::Cell< fem::Field_Interpolator* > mMasterDvFI;
            moris::Cell< fem::Field_Interpolator* > mSlaveDvFI;

            // master and slave properties
            moris::Cell< std::shared_ptr< fem::Property > > mMasterProp;
            moris::Cell< std::shared_ptr< fem::Property > > mSlaveProp;

            // master and slave constitutive models
            moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mMasterCM;
            moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mSlaveCM;

            // stabilization parameters
            moris::Cell< std::shared_ptr< Stabilization_Parameter > > mStabilizationParam;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            IQI(){};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            ~IQI(){};

//------------------------------------------------------------------------------
            /**
             * set dof types
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
                        MORIS_ERROR( false, "IQI::set_dof_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * return a cell of dof types
             * @param[ in ] aIsMaster enum master or slave
             */
            moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dof type list
                        return mMasterDofTypes;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dof type list
                        return mSlaveDofTypes;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "IQI::get_dof_type_list - can only be master or slave." );
                        return mMasterDofTypes;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * set dv types
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
                        MORIS_ERROR( false, "IQI::set_dv_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * return a cell of dv types
             * @param[ in ] aIsMaster enum master or slave
             */
            moris::Cell< moris::Cell< MSI::Dv_Type > > & get_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dv type list
                        return mMasterDvTypes;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dv type list
                        return mSlaveDvTypes;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Opt_IWG::get_dv_type_list - can only be master or slave." );
                        return mMasterDvTypes;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * set dof field interpolators
             * @param[ in ] aFieldInterpolators cell of field interpolator pointers
             * @param[ in ] aIsMaster           enum for master or slave
             */
            void set_dof_field_interpolators( moris::Cell< fem::Field_Interpolator* > & aFieldInterpolators,
                                              mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER )
            {
                // check input size
                MORIS_ASSERT( aFieldInterpolators.size() == this->get_global_dof_type_list( aIsMaster ).size(),
                              "IQI::set_dof_field_interpolators - wrong input size. " );

                // check field interpolator type
                bool tCheckFI = true;
                for( uint iFI = 0; iFI < aFieldInterpolators.size(); iFI++ )
                {
                    tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dof_type()( 0 ) == this->get_global_dof_type_list( aIsMaster )( iFI )( 0 ) );
                }
                MORIS_ASSERT( tCheckFI, "IQI::set_dof_field_interpolators - wrong field interpolator dof type. ");

                // set field interpolators
                this->get_dof_field_interpolators( aIsMaster ) = aFieldInterpolators;
            }

//------------------------------------------------------------------------------
            /**
             * get dof field interpolators
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
                        MORIS_ASSERT( false, "IQI::get_dof_field_interpolators - can only be master or slave." );
                        return mMasterDofFI;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * set dv field interpolators
             * @param[ in ] aFieldInterpolators cell of field interpolator pointers
             * @param[ in ] aIsMaster           enum for master or slave
             */
            void set_dv_field_interpolators( moris::Cell< fem::Field_Interpolator* > & aFieldInterpolators,
                                             mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER )
            {
                // check input size
                MORIS_ASSERT( aFieldInterpolators.size() == this->get_global_dv_type_list( aIsMaster ).size(),
                              "IQI::set_dv_field_interpolators - wrong input size. " );

                // check field interpolator type
                bool tCheckFI = true;
                for( uint iFI = 0; iFI < aFieldInterpolators.size(); iFI++ )
                {
                    tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dv_type()( 0 ) == this->get_global_dv_type_list( aIsMaster )( iFI )( 0 ) );
                }
                MORIS_ASSERT( tCheckFI, "IQI::set_dv_field_interpolators - wrong field interpolator dof type. ");

                // set field interpolators
                this->get_dv_field_interpolators( aIsMaster ) = aFieldInterpolators;
            }

//------------------------------------------------------------------------------
            /**
             * get dv field interpolators
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
                        // return master dv field interpolator pointers
                        return mMasterDvFI;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave dv lator pointers
                        return mSlaveDvFI;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "IQI::set_dv_field_interpolators - can only be master or slave." );
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
                 MORIS_ASSERT( this->get_dv_field_interpolators( aIsMaster ).size() == this->get_dv_type_list( aIsMaster ).size(),
                               "IQI::check_dv_field_interpolators - wrong FI size. " );

                // loop over the field interpolator pointers
                for( uint iFI = 0; iFI < this->get_dv_type_list( aIsMaster ).size(); iFI++ )
                {
                    // check that the field interpolator was set
                    MORIS_ASSERT( this->get_dv_field_interpolators( aIsMaster )( iFI ) != nullptr,
                                  "IQI::check_dv_field_interpolators - FI missing. " );
                }
             }

//------------------------------------------------------------------------------
             /**
              * get master or slave properties
              * @param[ in ]  aIsMaster   enum master or slave
              * @param[ out ] aProperties cell of property pointers
              */
             moris::Cell< std::shared_ptr< fem::Property > > & get_properties( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
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
                         MORIS_ASSERT( false, "IQI::get_properties - can only be master or slave." );
                         return mMasterProp;
                     }
                 }
             }

//------------------------------------------------------------------------------
             /**
              * get master or slave constitutive models
              * @param[ in ]  aIsMaster           enum master or slave
              * @param[ out ] aConstitutiveModels cell of constitutive model pointers
              */
             moris::Cell< std::shared_ptr< fem::Constitutive_Model > > & get_constitutive_models( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
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
                         MORIS_ASSERT( false, "IQI::get_constitutive_models - can only be master or slave." );
                         return mMasterCM;
                     }
                 }
             }

//------------------------------------------------------------------------------
            /**
             * get stabilization parameters
             * @param[ out ] mStabilizationParam cell of stabilization parameter pointers
             */
            moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > & get_stabilization_parameters()
            {
                // return stabilization parameter pointers
                return mStabilizationParam;
            }

//------------------------------------------------------------------------------
            /**
             * create a global dof type list including IWG, property and constitutive dependencies
             */
            void build_global_dof_type_list()
            {
                // MASTER-------------------------------------------------------
                // get the size of the dof type list
                uint tCounterMax = 0;

                // get number of dof types from IWG
                tCounterMax += mMasterDofTypes.size();

                // get number of dof types from properties
                for ( std::shared_ptr< Property > tProperty : mMasterProp )
                {
                    if( tProperty != nullptr )
                    {
                        tCounterMax += tProperty->get_dof_type_list().size();
                    }
                }

                // get number of dof types from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                {
                    if( tCM != nullptr )
                    {
                        tCounterMax += tCM->get_global_dof_type_list().size();
                    }
                }

                // get number of dof types from stabilization parameters
                for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
                {
                    if( tSP != nullptr )
                    {
                        tCounterMax += tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER ).size();
                    }
                }

                // set size for the global dof type list
                mMasterGlobalDofTypes.resize( tCounterMax );

                // set a size for the checkList (used to avoid repeating a dof type)
                moris::Cell< sint > tCheckList( tCounterMax, -1 );

                // init total dof counter
                uint tCounter = 0;

                // get dof type from penalty parameter
                for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
                {
                    // put the dof type in the checklist
                    tCheckList( tCounter ) = static_cast< uint >( mMasterDofTypes( iDOF )( 0 ) );

                    // put the dof type in the global type list
                    mMasterGlobalDofTypes( tCounter ) = mMasterDofTypes( iDOF );

                    // update the dof counter
                    tCounter++;
                }

                // get dof type from properties
                for ( std::shared_ptr< Property > tProperty : mMasterProp )
                {
                    if( tProperty != nullptr )
                    {
                        // get dof types for property
                        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

                        // loop on property dof type
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
                                // put the dof type in the checklist
                                tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                                // put the dof type in the global type list
                                mMasterGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                                // update dof counter
                                tCounter++;
                            }
                        }
                    }
                }

                // get dof type from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                {
                    if( tCM != nullptr )
                    {
                        // get dof types for constitutive model
                        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tCM->get_global_dof_type_list();

                        // loop on property dof type
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
                                // put the dof type in the checklist
                                tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                                // put the dof type in the global type list
                                mMasterGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                                // update dof counter
                                tCounter++;
                            }
                        }
                    }
                }

                // get dof type from stabilization parameters
                for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
                {
                    if ( tSP != nullptr )
                    {
                        // get dof types for constitutive model
                        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

                        // loop on property dof type
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
                                // put the dof type in the checklist
                                tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                                // put the dof type in the global type list
                                mMasterGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                                // update dof counter
                                tCounter++;
                            }
                        }
                    }
                }

                // get the number of unique dof type groups for the penalty parameter
                mMasterGlobalDofTypes.resize( tCounter );

                // SLAVE--------------------------------------------------------
                // get the size of the dof type list
                tCounterMax = 0;

                // get number of dof types from penalty parameter
                tCounterMax += mSlaveDofTypes.size();

                // get number of dof types from properties
                for ( std::shared_ptr< Property > tProperty : mSlaveProp )
                {
                    if( tProperty != nullptr )
                    {
                        tCounterMax += tProperty->get_dof_type_list().size();
                    }
                }

                // get number of dof types from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
                {
                    if ( tCM != nullptr )
                    {
                        tCounterMax += tCM->get_global_dof_type_list().size();
                    }
                }

                // get number of dof types from stabilization parameters
                for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
                {
                    if( tSP != nullptr )
                    {
                        tCounterMax += tSP->get_global_dof_type_list( mtk::Master_Slave::SLAVE ).size();
                    }
                }

                // set size for the global dof type list
                mSlaveGlobalDofTypes.resize( tCounterMax );

                // set a size for the checkList (used to avoid repeating a dof type)
                tCheckList.resize( tCounterMax, -1 );

                // init total dof counter
                tCounter = 0;

                // get dof type from penalty parameter
                for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
                {
                    // put the dof type in the checklist
                    tCheckList( tCounter ) = static_cast< uint >( mSlaveDofTypes( iDOF )( 0 ) );

                    // put the dof type in the global type list
                    mSlaveGlobalDofTypes( tCounter ) = mSlaveDofTypes( iDOF );

                    // update the dof counter
                    tCounter++;
                }

                // get dof type from properties
                for ( std::shared_ptr< Property > tProperty : mSlaveProp )
                {
                    if ( tProperty != nullptr )
                    {
                        // get dof types for property
                        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

                        // loop on property dof type
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
                                // put the dof type in the checklist
                                tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                                // put the dof type in the global type list
                                mSlaveGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                                // update dof counter
                                tCounter++;
                            }
                        }
                    }
                }

                // get dof type from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
                {
                    if( tCM != nullptr )
                    {
                        // get dof types for constitutive model
                        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tCM->get_global_dof_type_list();

                        // loop on property dof type
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
                                // put the dof type in the checklist
                                tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                                // put the dof type in the global type list
                                mSlaveGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                                // update dof counter
                                tCounter++;
                            }
                        }
                    }
                }

                // get dof type from stabilization parameters
                for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
                {
                    if ( tSP != nullptr )
                    {
                        // get dof types for constitutive model
                        moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

                        // loop on property dof type
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
                                // put the dof type in the checklist
                                tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                                // put the dof type in the global type list
                                mSlaveGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                                // update dof counter
                                tCounter++;
                            }
                        }
                    }
                }

                // get the number of unique dof type groups for the penalty parameter
                mSlaveGlobalDofTypes.resize( tCounter );
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
                      }
                      // if slave
                      case( mtk::Master_Slave::SLAVE ):
                      {
                          // return slave global dof type list
                          return mSlaveGlobalDofTypes;
                      }
                      // if none
                      default:
                      {
                          MORIS_ASSERT( false, "Opt_IWG::get_global_dof_type_list - can only be master or slave." );
                          return mMasterGlobalDofTypes;
                      }
                  }
              };

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
                                    if ( tProperty != nullptr )
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
                                    if (tCM != nullptr )
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
                      }
                      // if slave
                      case( mtk::Master_Slave::SLAVE ):
                      {
                          // return slave global dv type list
                          return mSlaveGlobalDvTypes;
                      }
                      // if none
                      default:
                      {
                          MORIS_ASSERT( false, "Opt_IWG::get_global_dv_type_list - can only be master or slave." );
                          return mMasterGlobalDvTypes;
                      }
                  }
              };

//------------------------------------------------------------------------------
            /**
             * evaluate the quantity of interest
             * @param[ in ] aQIVal quantity of interest matrix to fill
             */
            virtual void compute_QI( Matrix< DDRMat > & aQIVal )
            {
                MORIS_ERROR( false, "IQI::compute_QI - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the quantity of interest wrt to dof types
             * @param[ in ] adQIdDof matrix to fill with derivative of the QoI wrt dof types
             */
            virtual void compute_dQIdDof( Matrix< DDRMat > & adQIdDof )
            {
                MORIS_ERROR( false, "IQI::compute_dIQIdDof - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the quantity of interest wrt to dv types
             * @param[ in ] adQIdDv matrix to fill with derivative of the QI wrt dv types
             */
            virtual void compute_dQIdDv( Matrix< DDRMat > & adQIdDv )
            {
                MORIS_ERROR( false, "IQI::compute_dIQIdDv - This function does nothing. " );
            }
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IQI_HPP_ */

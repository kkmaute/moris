/*
 * cl_FEM_IQI.cpp
 *
 *  Created on: Dec 4, 2019
 *      Author: noel
 */

#include "cl_FEM_IQI.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        /**
         * reset evaluation flags
         */
        void IQI::reset_eval_flags()
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
        void IQI::set_field_interpolator_manager( Field_Interpolator_Manager * aFieldInterpolatorManager,
                                                  mtk::Master_Slave            aIsMaster )
        {
            // FIXME why does this not work?
            //this->get_field_interpolator_manager( aIsMaster ) = aFieldInterpolatorManager;

            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ) :
                {
                    mMasterFIManager = aFieldInterpolatorManager;
                    break;
                }

                case ( mtk::Master_Slave::SLAVE ) :
                {
                    mSlaveFIManager = aFieldInterpolatorManager;
                    break;
                }

                default :
                {
                    MORIS_ERROR( false, "IWG::set_field_interpolator_manager - can only be master or slave");
                    break;
                }
            }
            //END FIXME

            // loop over the the SP
            for( std::shared_ptr< Stabilization_Parameter > tSP : this->get_stabilization_parameters() )
            {
                if ( tSP != nullptr )
                {
                    // set the field interpolator manager for the SP
                    tSP->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsMaster ), aIsMaster );

                    // set th efem set pointer for the SP
                    tSP->set_set_pointer( mSet );
                }
            }

            // loop over the constitutive models
            for( std::shared_ptr< Constitutive_Model > tCM : this->get_constitutive_models( aIsMaster ) )
            {
                if ( tCM != nullptr )
                {
                    // set the field interpolator manager for the CM
                    tCM->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsMaster ) );

                    // set the fem set pointe for the CM
                    tCM->set_set_pointer( mSet );
                }
            }

            // loop over the properties
            for( std::shared_ptr< Property > tProp : this->get_properties( aIsMaster ) )
            {
                if ( tProp != nullptr )
                {
                    // set the field interpolator manager for the property
                    tProp->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsMaster ) );

                    // set the fem set pointer for the property
                    tProp->set_set_pointer( mSet );
                }
            }
        }

//------------------------------------------------------------------------------
        void IQI::get_non_unique_global_dof_type_list( moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // init counter for dof types
            uint tCounter = 0;

            // loop over direct master dof dependencies
            for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
            {
                // update counter
                tCounter += mMasterDofTypes( iDOF ).size();
            }

            // loop over direct slave dof dependencies
            for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
            {
                // update counter
                tCounter += mSlaveDofTypes( iDOF ).size();
            }

            // loop over the master properties
            for ( std::shared_ptr< Property > tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tProperty->get_non_unique_dof_types( tActiveDofType );

                    //update counter
                    tCounter += tActiveDofType.size();
                }
            }

            // loop over slave properties
            for ( std::shared_ptr< Property > tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property nn unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tProperty->get_non_unique_dof_types( tActiveDofType );

                    // update counter
                    tCounter += tActiveDofType.size();
                }
            }

            // loop over master constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tCM->get_non_unique_dof_types( tActiveDofType );

                    // update counter
                    tCounter += tActiveDofType.size();
                }
            }

            // loop over slave constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    // get CM non unique dof type list
                    moris::Cell< MSI::Dof_Type >  tActiveDofType;
                    tCM->get_non_unique_dof_types( tActiveDofType );

                    // update counter
                    tCounter += tActiveDofType.size();
                }
            }

            // loop over master stabilization parameters
            for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // FIXME get SP non unique dof type list
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tMasterActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tSlaveActiveDofType  = tSP->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

                    // update counter
                    tCounter += tMasterActiveDofType.size() + tSlaveActiveDofType.size();
                }
            }

            // reserve memory for dof type list
            aDofTypes.reserve( tCounter );

            // loop over master dof direct dependencies
            for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
            {
                // populate the dof list
                aDofTypes.append( mMasterDofTypes( iDOF ) );
            }

            // loop over slave dof direct dependencies
            for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
            {
                //populate the dof list
                aDofTypes.append( mSlaveDofTypes( iDOF )  );
            }

            // loop over master properties
            for ( std::shared_ptr< Property > tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tProperty->get_non_unique_dof_types( tActiveDofType );

                    // populate the dof list
                    aDofTypes.append( tActiveDofType );
                }
            }

            // loop over slave properties
            for ( std::shared_ptr< Property > tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tProperty->get_non_unique_dof_types( tActiveDofType );

                    // populate the dof list
                    aDofTypes.append( tActiveDofType );
                }
            }

            // loop over the master constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tCM->get_non_unique_dof_types( tActiveDofType );

                    // populate the dof list
                    aDofTypes.append( tActiveDofType );
                }
            }

            // loop over the slave constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    // get CM non unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tCM->get_non_unique_dof_types( tActiveDofType );

                    // populate the dof list
                    aDofTypes.append( tActiveDofType );
                }
            }

            // loop over the stabilization parameters
            for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // FIXME get SP non unique master dof type list
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tMasterActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

                    // loop over dof type groups
                    for ( uint iDOF = 0; iDOF < tMasterActiveDofType.size(); iDOF++ )
                    {
                        // populate the dof list
                        aDofTypes.append( tMasterActiveDofType( iDOF ) );
                    }

                    // FIXME get SP non unique dof type list
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tSlaveActiveDofType = tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

                    // loop over dof type groups
                    for ( uint iDOF = 0; iDOF < tSlaveActiveDofType.size(); iDOF++ )
                    {
                        // populate the dof list
                        aDofTypes.append( tSlaveActiveDofType( iDOF ) );
                    }
                }
            }
        }

//------------------------------------------------------------------------------
        void IQI::build_global_dof_type_list()
        {
            // MASTER-------------------------------------------------------
            // set size for the global dof type list
            uint tNumDofTypes = mSet->get_num_dof_types();

            // set size for the global dof type list
            mMasterGlobalDofTypes.reserve( tNumDofTypes );

            // set a size for the checkList ( used to avoid repeating a dof type)
            Matrix< DDSMat > tCheckList( tNumDofTypes, 1, -1 );

            // get dof type from direct dependencies
            for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
            {
                // get set index for dof type
                sint tDofTypeIndex = mSet->get_dof_index_for_type_1( mMasterDofTypes( iDOF )( 0 ) );  //FIXME'

                // put the dof type in the checklist
                tCheckList( tDofTypeIndex ) = 1;

                // put the dof type in the global type list
                mMasterGlobalDofTypes.push_back( mMasterDofTypes( iDOF ) );
            }

            // get dof type from master properties
            for ( std::shared_ptr< Property > tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get dof types for property
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType
                    = tProperty->get_dof_type_list();

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // get set index for dof type
                        sint tDofTypeindex = mSet->get_dof_index_for_type_1( tActiveDofType( iDOF )( 0 ) );

                        // if dof enum not in the list
                        if ( tCheckList( tDofTypeindex) != 1 )
                        {
                            // put the dof type in the checklist
                            tCheckList( tDofTypeindex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                        }
                    }
                }
            }

            // get dof type from master constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get dof types for constitutive model
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType
                    = tCM->get_global_dof_type_list();

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // get set index for dof type
                        sint tDofTypeindex = mSet->get_dof_index_for_type_1( tActiveDofType( iDOF )( 0 ) );

                        // if dof enum not in the list
                        if ( tCheckList( tDofTypeindex) != 1 )
                        {
                            // put the dof type in the checklist
                            tCheckList( tDofTypeindex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                        }
                    }
                }
            }

            // get dof type from master stabilization parameters
            for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get dof types for constitutive model
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType
                    = tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // get set index for dof type
                        sint tDofTypeindex = mSet->get_dof_index_for_type_1( tActiveDofType( iDOF )( 0 ) );

                        // if dof enum not in the list
                        if ( tCheckList( tDofTypeindex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tCheckList( tDofTypeindex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                        }
                    }
                }
            }

            // reduce size of dof list to fit unique list
            mMasterGlobalDofTypes.shrink_to_fit();

            // SLAVE--------------------------------------------------------

            // set size for the global dof type list
            mSlaveGlobalDofTypes.reserve( tNumDofTypes );

            // set a size for the checkList ( used to avoid repeating a dof type)
            tCheckList.fill( -1 );

            // get dof type from slave direct dependencies
            for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
            {
                // get set index for dof type
                sint tDofTypeindex = mSet->get_dof_index_for_type_1( mSlaveDofTypes( iDOF )( 0 ) );

                // put the dof type in the checklist
                tCheckList( tDofTypeindex ) = 1;

                // put the dof type in the global type list
                mSlaveGlobalDofTypes.push_back( mSlaveDofTypes( iDOF ) );
            }

            // get dof type from master properties
            for ( std::shared_ptr< Property > tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get dof types for property
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType
                    = tProperty->get_dof_type_list();

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // get set index for dof type
                        sint tDofTypeindex = mSet->get_dof_index_for_type_1( tActiveDofType( iDOF )( 0 ) );

                        // if dof enum not in the list
                        if ( tCheckList( tDofTypeindex) != 1 )
                        {
                            // put the dof type in the checklist
                            tCheckList( tDofTypeindex ) = 1;

                            // put the dof type in the global type list
                            mSlaveGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                        }
                    }
                }
            }

            // get dof type from slave constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
            {
                if ( tCM != nullptr )
                {
                    // get dof types for constitutive model
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType
                    = tCM->get_global_dof_type_list();

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_dof_index_for_type_1( tActiveDofType( iDOF )( 0 ) );

                        // if dof enum not in the list
                        if ( tCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mSlaveGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
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
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType
                    = tSP->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // get set index for dof type
                        sint tDofTypeindex = mSet->get_dof_index_for_type_1( tActiveDofType( iDOF )( 0 ) );

                        // if dof enum not in the list
                        if ( tCheckList( tDofTypeindex) != 1 )
                        {
                            // put the dof type in the checklist
                            tCheckList( tDofTypeindex ) = 1;

                            // put the dof type in the global type list
                            mSlaveGlobalDofTypes.push_back( tActiveDofType( iDOF ) );
                        }
                    }
                }
            }

            // reduce size of dof list to fit unique list
            mSlaveGlobalDofTypes.shrink_to_fit();
        }

////------------------------------------------------------------------------------
//        // FIXME the mFIManager needs to be use by
//        // properties, CM and SP
//        // FIXME why not master and slave together?
//        void IQI::set_dof_field_interpolators( mtk::Master_Slave aIsMaster )
//        {
//            // set field interpolators for the SP
//            for( std::shared_ptr< Stabilization_Parameter > tSP : this->get_stabilization_parameters() )
//            {
//                if ( tSP != nullptr )
//                {
//                    // get the list of dof types for the SP
//                    moris::Cell< moris::Cell< MSI::Dof_Type > > tSPDofTypes
//                    = tSP->get_global_dof_type_list( aIsMaster );
//
//                    // get the number of dof type for the SP
//                    uint tNumDofTypes = tSPDofTypes.size();
//
//                    // set the size of the field interpolators list for the SP
//                    moris::Cell< Field_Interpolator* > tSPFIs( tNumDofTypes, nullptr );
//
//                    // loop over the dof types
//                    for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
//                    {
//                        // grab the field interpolator for the dof type
//                        tSPFIs( iDof ) = this->get_field_interpolator_manager( aIsMaster )
//                                             ->get_field_interpolators_for_type( tSPDofTypes( iDof )( 0 ) );
//                    }
//
//                    // set the field interpolators for the SP
//                    tSP->set_dof_field_interpolators( tSPFIs, aIsMaster );
//
//                    // set the field interpolator manager for the SP
//                    tSP->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsMaster ), aIsMaster );
//
//                    // set th efem set pointer for the SP
//                    tSP->set_set_pointer( mSet );
//                }
//            }
//
//            // set field interpolators for constitutive models
//            for( std::shared_ptr< Constitutive_Model > tCM : this->get_constitutive_models( aIsMaster ) )
//            {
//                if ( tCM != nullptr )
//                {
//                    // get the list of dof types for the CM
//                    moris::Cell< moris::Cell< MSI::Dof_Type > > tCMDofTypes
//                    = tCM->get_global_dof_type_list();
//
//                    // get the number of dof type for the CM
//                    uint tNumDofTypes = tCMDofTypes.size();
//
//                    // set the size of the field interpolators list for the CM
//                    moris::Cell< Field_Interpolator* > tCMFIs( tNumDofTypes, nullptr );
//
//                    // loop over the dof types
//                    for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
//                    {
//                        // fill the field interpolators list for the CM
//                        tCMFIs( iDof ) = this->get_field_interpolator_manager( aIsMaster )
//                                             ->get_field_interpolators_for_type( tCMDofTypes( iDof )( 0 ) );
//                    }
//
//                    // set the field interpolators for the CM
//                    tCM->set_dof_field_interpolators( tCMFIs );
//
//                    // set the field interpolator manager for the CM
//                    tCM->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsMaster ) );
//
//                    // set the fem set pointe for the CM
//                    tCM->set_set_pointer( mSet );
//                }
//            }
//
//            // set field interpolators for properties
//            for( std::shared_ptr< Property > tProp : this->get_properties( aIsMaster ) )
//            {
//                if ( tProp != nullptr )
//                {
//                    // get the list of dof types for the property
//                    moris::Cell< moris::Cell< MSI::Dof_Type > > tPropDofTypes
//                    = tProp->get_dof_type_list();
//
//                    // get the number of dof type for the property
//                    uint tNumDofTypes = tPropDofTypes.size();
//
//                    // set the size of the field interpolators list for the property
//                    moris::Cell< Field_Interpolator* > tPropFIs( tNumDofTypes, nullptr );
//
//                    // loop over the dof types
//                    for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
//                    {
//                        tPropFIs( iDof ) = this->get_field_interpolator_manager( aIsMaster )
//                                               ->get_field_interpolators_for_type( tPropDofTypes( iDof )( 0 ) );
//                    }
//
//                    // set the field interpolators for the property
//                    tProp->set_dof_field_interpolators( tPropFIs );
//
//                    // set the field interpolator manager for the property
//                    tProp->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsMaster ) );
//
//                    // set the fem set pointer for the property
//                    tProp->set_set_pointer( mSet );
//                }
//            }
//        }

// //------------------------------------------------------------------------------
//        void IQI::get_non_unique_global_dv_type_list( moris::Cell< MSI::Dv_Type > & aDvTypes )
//        {
//            // init counter for dv types
//            uint tCounter = 0;
//
//            // loop over direct master dv dependencies
//            for ( uint iDv = 0; iDv < mMasterDvTypes.size(); iDv++ )
//            {
//                // update counter
//                tCounter += mMasterDvTypes( iDv ).size();
//            }
//
//            // loop over direct slave dv dependencies
//            for ( uint iDv = 0; iDv < mSlaveDvTypes.size(); iDv++ )
//            {
//                // update counter
//                tCounter += mSlaveDvTypes( iDv ).size();
//            }
//
//            // loop over the master properties
//            for ( std::shared_ptr< Property > tProperty : mMasterProp )
//            {
//                if ( tProperty != nullptr )
//                {
//                    // get property non unique dv type list
//                    moris::Cell< MSI::Dv_Type > tActiveDvType;
//                    tProperty->get_non_unique_dv_types( tActiveDvType );
//
//                    //update counter
//                    tCounter += tActiveDvType.size();
//                }
//            }
//
//            // loop over slave properties
//            for ( std::shared_ptr< Property > tProperty : mSlaveProp )
//            {
//                if ( tProperty != nullptr )
//                {
//                    // get property nn unique dv type list
//                    moris::Cell< MSI::Dv_Type > tActiveDvType;
//                    tProperty->get_non_unique_dv_types( tActiveDvType );
//
//                    // update counter
//                    tCounter += tActiveDvType.size();
//                }
//            }
//
//            // loop over master constitutive models
//            for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
//            {
//                if ( tCM != nullptr )
//                {
//                    // get CM non unique dv type list
//                    moris::Cell< MSI::Dv_Type > tActiveDvType;
//                    tCM->get_non_unique_dv_types( tActiveDvType );
//
//                    // update counter
//                    tCounter += tActiveDvType.size();
//                }
//            }
//
//            // loop over slave constitutive models
//            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
//            {
//                if( tCM != nullptr )
//                {
//                    // get CM non unique dv type list
//                    moris::Cell< MSI::Dv_Type >  tActiveDvType;
//                    tCM->get_non_unique_dv_types( tActiveDvType );
//
//                    // update counter
//                    tCounter += tActiveDvType.size();
//                }
//            }
//
//            // loop over master stabilization parameters
//            for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
//            {
//                if ( tSP != nullptr )
//                {
//                    // FIXME get SP non unique dv type list
//                    moris::Cell< moris::Cell< MSI::Dv_Type > > tMasterActiveDvType = tSP->get_global_dv_type_list( mtk::Master_Slave::MASTER );
//                    moris::Cell< moris::Cell< MSI::Dv_Type > > tSlaveActiveDvType  = tSP->get_global_dv_type_list( mtk::Master_Slave::SLAVE );
//
//                    // update counter
//                    tCounter += tMasterActiveDvType.size() + tSlaveActiveDvType.size();
//                }
//            }
//
//            // reserve memory for dv type list
//            aDvTypes.reserve( tCounter );
//
//            // loop over master dv direct dependencies
//            for ( uint iDv = 0; iDv < mMasterDvTypes.size(); iDv++ )
//            {
//                // populate the dv list
//                aDvTypes.append( mMasterDvTypes( iDv ) );
//            }
//
//            // loop over slave dv direct dependencies
//            for ( uint iDv = 0; iDv < mSlaveDvTypes.size(); iDv++ )
//            {
//                //populate the dv list
//                aDvTypes.append( mSlaveDvTypes( iDv )  );
//            }
//
//            // loop over master properties
//            for ( std::shared_ptr< Property > tProperty : mMasterProp )
//            {
//                if ( tProperty != nullptr )
//                {
//                    // get property non unique dv type list
//                    moris::Cell< MSI::Dv_Type > tActiveDvType;
//                    tProperty->get_non_unique_dv_types( tActiveDvType );
//
//                    // populate the dv list
//                    aDvTypes.append( tActiveDvType );
//                }
//            }
//
//            // loop over slave properties
//            for ( std::shared_ptr< Property > tProperty : mSlaveProp )
//            {
//                if ( tProperty != nullptr )
//                {
//                    // get property non unique dv type list
//                    moris::Cell< MSI::Dv_Type > tActiveDvType;
//                    tProperty->get_non_unique_dv_types( tActiveDvType );
//
//                    // populate the dv list
//                    aDvTypes.append( tActiveDvType );
//                }
//            }
//
//            // loop over the master constitutive models
//            for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
//            {
//                if ( tCM != nullptr )
//                {
//                    // get CM non unique dv type list
//                    moris::Cell< MSI::Dv_Type > tActiveDvType;
//                    tCM->get_non_unique_dv_types( tActiveDvType );
//
//                    // populate the dv list
//                    aDvTypes.append( tActiveDvType );
//                }
//            }
//
//            // loop over the slave constitutive models
//            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
//            {
//                if( tCM != nullptr )
//                {
//                    // get CM non unique dv type list
//                    moris::Cell< MSI::Dv_Type > tActiveDvType;
//                    tCM->get_non_unique_dv_types( tActiveDvType );
//
//                    // populate the dv list
//                    aDvTypes.append( tActiveDvType );
//                }
//            }
//
//            // loop over the stabilization parameters
//            for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
//            {
//                if ( tSP != nullptr )
//                {
//                    // FIXME get SP non unique master dv type list
//                    moris::Cell< moris::Cell< MSI::Dv_Type > > tMasterActiveDvType = tSP->get_global_dv_type_list( mtk::Master_Slave::MASTER );
//
//                    // loop over dv type groups
//                    for ( uint iDv = 0; iDv < tMasterActiveDvType.size(); iDv++ )
//                    {
//                        // populate the dv list
//                        aDvTypes.append( tMasterActiveDvType( iDv ) );
//                    }
//
//                    // FIXME get SP non unique dv type list
//                    moris::Cell< moris::Cell< MSI::Dv_Type > > tSlaveActiveDvType = tSP->get_global_dv_type_list( mtk::Master_Slave::MASTER );
//
//                    // loop over dv type groups
//                    for ( uint iDv = 0; iDv < tSlaveActiveDvType.size(); iDv++ )
//                    {
//                        // populate the dv list
//                        aDvTypes.append( tSlaveActiveDvType( iDv ) );
//                    }
//                }
//            }
//        }

//------------------------------------------------------------------------------
        void IQI::build_global_dv_type_list()
        {
            MORIS_ERROR( false, "This function does nothing" );
        }

//------------------------------------------------------------------------------
        void IQI::set_dv_field_interpolators( mtk::Master_Slave aIsMaster )
        {
            MORIS_ERROR( false, "This function does nothing" );
        }

//------------------------------------------------------------------------------
        void IQI::set_geometry_interpolator( Geometry_Interpolator* aGeometryInterpolator,
                                             mtk::Master_Slave      aIsMaster )
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
    }/* end_namespace_fem */
}/* end_namespace_moris */


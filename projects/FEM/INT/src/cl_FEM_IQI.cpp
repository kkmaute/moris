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
        void IQI::get_non_unique_dof_types( moris::Cell< MSI::Dof_Type > & aDofTypes )
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
                    // get SP non unique dof types
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tSP->get_non_unique_dof_types( tActiveDofType );

                    // update counter
                    tCounter += tActiveDofType.size();
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
                    // get SP non unique dof types
                    moris::Cell< MSI::Dof_Type > tActiveDofType;
                    tSP->get_non_unique_dof_types( tActiveDofType );

                    // populate the dof list
                    aDofTypes.append( tActiveDofType );
                }
            }
        }

//------------------------------------------------------------------------------
        void IQI::get_non_unique_dof_and_dv_types( moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                   moris::Cell< MSI::Dv_Type >  & aDvTypes )
        {
            // init counters for dof and dv types
            uint tDofCounter = 0;
            uint tDvCounter  = 0;

            // get number of direct master dof dependencies
            for ( uint iDof = 0; iDof < mMasterDofTypes.size(); iDof++ )
            {
                tDofCounter += mMasterDofTypes( iDof ).size();
            }

            // get number of direct master dv dependencies
            for ( uint iDv = 0; iDv < mMasterDvTypes.size(); iDv++ )
            {
                tDvCounter += mMasterDvTypes( iDv ).size();
            }

            // get number of direct slave dof dependencies
            for ( uint iDof = 0; iDof < mSlaveDofTypes.size(); iDof++ )
            {
                tDofCounter += mSlaveDofTypes( iDof ).size();
            }

            // get number of direct slave dv dependencies
            for ( uint iDv = 0; iDv < mSlaveDvTypes.size(); iDv++ )
            {
                tDvCounter += mSlaveDvTypes( iDv ).size();
            }

            // loop over the master properties
            for ( std::shared_ptr< Property > tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type > tActiveDvTypes;
                    tProperty->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                                tActiveDvTypes );

                    //update dof and dv counters
                    tDofCounter += tActiveDofTypes.size();
                    tDvCounter  += tActiveDvTypes.size();
                }
            }

            // loop over slave properties
            for ( std::shared_ptr< Property > tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type lists
                    // get property non unique dof and dv type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type > tActiveDvTypes;
                    tProperty->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                                tActiveDvTypes );

                    // update dof and dv counter
                    tDofCounter += tActiveDofTypes.size();
                    tDvCounter  += tActiveDvTypes.size();
                }
            }

            // loop over master constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tCM->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                          tActiveDvTypes );

                    // update dof and dv counters
                    tDofCounter += tActiveDofTypes.size();
                    tDvCounter  += tActiveDvTypes.size();

                }
            }

            // loop over slave constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tCM->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                          tActiveDvTypes );

                    // update dof and dv counters
                    tDofCounter += tActiveDofTypes.size();
                    tDvCounter  += tActiveDvTypes.size();
                }
            }

            // loop over master stabilization parameters
            for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get SP non unique dof type list
                    moris::Cell< MSI::Dof_Type >  tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tSP->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                          tActiveDvTypes );

                    // update dof and dv counters
                    tDofCounter += tActiveDofTypes.size();
                    tDvCounter  += tActiveDvTypes.size();
                }
            }

            // reserve memory for dof and dv type lists
            aDofTypes.reserve( tDofCounter );
            aDvTypes.reserve( tDvCounter );

            // loop over master dof direct dependencies
            for ( uint iDof = 0; iDof < mMasterDofTypes.size(); iDof++ )
            {
                // populate the dof list
                aDofTypes.append( mMasterDofTypes( iDof ) );
            }

            // loop over master dv direct dependencies
            for ( uint iDv = 0; iDv < mMasterDvTypes.size(); iDv++ )
            {
                // populate the dv list
                aDvTypes.append( mMasterDvTypes( iDv ) );
            }

            // loop over slave dof direct dependencies
            for ( uint iDof = 0; iDof < mSlaveDofTypes.size(); iDof++ )
            {
                //populate the dof list
                aDofTypes.append( mSlaveDofTypes( iDof )  );
            }

            // loop over slave dv direct dependencies
            for ( uint iDv = 0; iDv < mSlaveDvTypes.size(); iDv++ )
            {
                //populate the dv list
                aDvTypes.append( mSlaveDvTypes( iDv )  );
            }

            // loop over master properties
            for ( std::shared_ptr< Property > tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tProperty->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                                tActiveDvTypes );

                    // populate the dof and dv lists
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }

            // loop over slave properties
            for ( std::shared_ptr< Property > tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tProperty->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                                tActiveDvTypes );

                    // populate the dof and dv lists
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }

            // loop over the master constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tCM->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                          tActiveDvTypes );

                    // populate the dof and dv lists
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }

            // loop over the slave constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tCM->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                          tActiveDvTypes );

                    // populate the dof and dv lists
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }

            // loop over the stabilization parameters
            for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get SP non unique master dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tSP->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                          tActiveDvTypes );

                    // populate the dof and dv lists
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
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
                    // get dof types for stabilization parameters
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

////------------------------------------------------------------------------------
//        void IQI::build_global_dv_type_list()
//        {
//            MORIS_ERROR( false, "This function does nothing" );
//        }

//------------------------------------------------------------------------------
        // FIXME
        moris::Cell < enum MSI::Dof_Type > IQI::get_requested_dof_types()
        {
            this->get_global_dof_type_list();

            moris::Cell < enum MSI::Dof_Type > tRequestedDofTypes( mMasterGlobalDofTypes.size() );

            for ( uint iDofGroup = 0; iDofGroup < mMasterGlobalDofTypes.size(); iDofGroup++ )
            {
                tRequestedDofTypes( iDofGroup ) = mMasterGlobalDofTypes( iDofGroup )( 0 );
            }
            return tRequestedDofTypes;
        }

//------------------------------------------------------------------------------
        void IQI::compute_dQIdDof_FD( Matrix< DDRMat >  & adQIdDofFD,
                                      real                aPerturbation )
        {
            // get number of dof coefficients
            uint tNumCoeff = 0;

            // get the requested dof types
            moris::Cell < enum MSI::Dof_Type > tRequestedDofTypes = this->get_requested_dof_types();

            // loop over the requested dof types
            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                // get the set index for dof type
                sint tDofIndex = mSet->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER );

                // if the dof type was set
                if( tDofIndex != -1 )
                {
                    // update number of coefficients
                    tNumCoeff += mMasterFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )
                                                 ->get_number_of_space_time_coefficients();
                }
            }

            // set size for dQIdDof
            adQIdDofFD.set_size( tNumCoeff, 1, 0.0 );

            // get master number of dof types
            uint tNumDofType = tRequestedDofTypes.size();

            // loop over the dof types
            for( uint iFI = 0; iFI < tNumDofType; iFI++ )
            {
                // get set index for dof type
                uint tDofIndex = mSet->get_dof_index_for_type( tRequestedDofTypes( iFI ), mtk::Master_Slave::MASTER );

                // get position in return vector
                uint tStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );

                // get field interpolator for dof type
                Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( tRequestedDofTypes( iFI ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

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
                        tCoeffPert( iCoeffRow, iCoeffCol ) += aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the QI
                        Matrix< DDRMat > tQIValPlus;
                        this->compute_QI( tQIValPlus );

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) += - aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the QI
                        Matrix< DDRMat > tQIValMinus;
                        this->compute_QI( tQIValMinus );

                        // evaluate dQIdDof
                        adQIdDofFD( tStartRow + tDofCounter )
                        = ( tQIValPlus( 0 ) - tQIValMinus( 0 ) ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }
        }

//------------------------------------------------------------------------------
        bool IQI::check_dQIdDof_FD( real               aPerturbation,
                                    real               aEpsilon,
                                    Matrix< DDRMat > & adQIdDof,
                                    Matrix< DDRMat > & adQIdDofFD )
        {
            // compute dQIdDof with IQI
            this->compute_dQIdDof( adQIdDof );

            // compute dQIdDof by FD
            this->compute_dQIdDof_FD( adQIdDofFD, aPerturbation );

            //define a boolean for check
            bool tCheckdQIdDof = true;

            // check if adQIdDof and adQIdDofFD have the same size
            tCheckdQIdDof = tCheckdQIdDof && ( adQIdDof.n_rows() == adQIdDofFD.n_rows());
            tCheckdQIdDof = tCheckdQIdDof && ( adQIdDof.n_cols() == adQIdDofFD.n_cols());

            // loop over the rows
            for ( uint iRow = 0; iRow < adQIdDof.n_rows(); iRow++ )
            {
                // loop over the columns
                for( uint jCol = 0; jCol < adQIdDof.n_cols(); jCol++ )
                {
                    // check each components
                    tCheckdQIdDof = tCheckdQIdDof && ( adQIdDof( iRow, jCol ) - adQIdDofFD( iRow, jCol ) < aEpsilon );
                }
            }

            // return bool
            return tCheckdQIdDof;
        }

//------------------------------------------------------------------------------
            void IQI::compute_dQIdDv_FD( Matrix< DDRMat > & adQIdpMatFD,
                                         Matrix< DDRMat > & adQIdpGeoFD,
                                         real               aPerturbation )
            {
                // MATERIAL PDV

                // get master number of dv types
                uint tNumDvType = mMasterGlobalDvTypes.size();

                // set size for adrdpdvFD
                uint tNumCols = 0;
                for( uint iFI = 0; iFI < tNumDvType; iFI++ )
                {
                    // get the dv index in the set
                    uint tDvIndex = mSet->get_dv_index_for_type( mMasterGlobalDvTypes( iFI )( 0 ),
                                                                 mtk::Master_Slave::MASTER );

                    // get the number of cols
                    tNumCols += mSet->get_dv_assembly_map()( tDvIndex )( 0, 1 )
                              - mSet->get_dv_assembly_map()( tDvIndex )( 0, 0 ) + 1;
                }

                adQIdpMatFD.set_size( 1, tNumCols, 0.0 );

                // loop over the dv types associated with a FI
                for( uint iFI = 0; iFI < tNumDvType; iFI++ )
                {
                    // get the dv index in the set
                    uint tDvIndex = mSet->get_dv_index_for_type( mMasterGlobalDvTypes( iFI )( 0 ),
                                                                 mtk::Master_Slave::MASTER );

                    // get the FI for the dv type
                    Field_Interpolator * tFI
                    = mMasterFIManager->get_field_interpolators_for_type( mMasterGlobalDvTypes( iFI )( 0 ) );

                    // get number of master FI bases and fields
                    uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                    uint tDerNumFields = tFI->get_number_of_fields();

                    // coefficients for dof type wrt which derivative is computed
                    Matrix< DDRMat > tCoeff = tFI->get_coeff();

                    // init dv coeff counter
                    uint tCounter = 0;

                    // loop over the coefficient column
                    for( uint iCoeffCol = 0; iCoeffCol< tDerNumFields; iCoeffCol++ )
                    {
                        // loop over the coefficient row
                        for( uint iCoeffRow = 0; iCoeffRow< tDerNumBases; iCoeffRow++  )
                        {
                            // perturbation of the coefficent
                            Matrix< DDRMat > tCoeffPert = tCoeff;
                            tCoeffPert( iCoeffRow, iCoeffCol ) += aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                            // setting the perturbed coefficients
                            tFI->set_coeff( tCoeffPert );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // evaluate the QI
                            Matrix< DDRMat > tQIValPlus;
                            this->compute_QI( tQIValPlus );

                            // perturbation of the coefficent
                            tCoeffPert = tCoeff;
                            tCoeffPert( iCoeffRow, iCoeffCol ) += - aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                            // setting the perturbed coefficients
                            tFI->set_coeff( tCoeffPert );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // evaluate the QI
                            Matrix< DDRMat > tQIValMinus;
                            this->compute_QI( tQIValMinus );

                            // evaluate Jacobian
                            uint tDvAssemblyStart = mSet->get_dv_assembly_map()( tDvIndex )( 0, 0 );

                            adQIdpMatFD( tDvAssemblyStart + tCounter )
                            = ( tQIValPlus( 0 ) - tQIValMinus( 0 ) ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                            // update counter
                            tCounter++;
                        }
                    }
                    // reset the coefficients values
                    tFI->set_coeff( tCoeff );
                }

                // GEOMETRY PDV

                // get the GI for the IG element considered
                Geometry_Interpolator * tGI = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

                // get number of master GI bases and space dimensions
                uint tDerNumBases      = tGI->get_number_of_space_bases();
                uint tDerNumDimensions = tGI->get_number_of_space_dimensions();

                // set size for adrdpdvGeoFD
                adQIdpGeoFD.set_size( 1, tDerNumBases * tDerNumDimensions, 0.0 );

                // coefficients for dv type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tGI->get_space_coeff();

                // init dv counter
                uint tDvCounter = 0;

                // loop over the spatial directions
                for( uint iCoeffCol = 0; iCoeffCol< tDerNumDimensions; iCoeffCol++ )
                {
                    // loop over the IG nodes
                    for( uint iCoeffRow = 0; iCoeffRow< tDerNumBases; iCoeffRow++  )
                    {
                        // perturbation of the coefficent
                        Matrix< DDRMat > tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) += aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tGI->set_space_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the QI
                        Matrix< DDRMat > tQIValPlus;
                        this->compute_QI( tQIValPlus );

                        // perturbation of the coefficient
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) += - aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tGI->set_space_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the QI
                        Matrix< DDRMat > tQIValMinus;
                        this->compute_QI( tQIValMinus );

                        // evaluate drdpdvGeo
                        adQIdpGeoFD( tDvCounter )
                        = ( tQIValPlus( 0 ) - tQIValMinus( 0 ) ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dv counter
                        tDvCounter++;
                    }
                }
                // reset the coefficients values
                tGI->set_space_coeff( tCoeff );
            }

//------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */


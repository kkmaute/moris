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

            // loop over the the SP
            for( std::shared_ptr< Stabilization_Parameter > tSP : this->get_stabilization_parameters() )
            {
                if ( tSP != nullptr )
                {
                    // set the field interpolator manager for the SP
                    tSP->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsMaster ), aIsMaster );

                    // set the fem set pointer for the SP
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
        void IQI::get_non_unique_dof_and_dv_types( moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                   moris::Cell< PDV_Type >        & aDvTypes )
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
                    moris::Cell< PDV_Type >        tActiveDvTypes;
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
                    moris::Cell< PDV_Type >        tActiveDvTypes;
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
                    moris::Cell< PDV_Type >        tActiveDvTypes;
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
                    moris::Cell< PDV_Type >        tActiveDvTypes;
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
                    moris::Cell< PDV_Type >         tActiveDvTypes;
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
                    moris::Cell< PDV_Type >        tActiveDvTypes;
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
                    moris::Cell< PDV_Type >        tActiveDvTypes;
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
                    moris::Cell< PDV_Type >        tActiveDvTypes;
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
                    moris::Cell< PDV_Type >        tActiveDvTypes;
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
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    tSP->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                          tActiveDvTypes );

                    // populate the dof and dv lists
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }
        }

//------------------------------------------------------------------------------
        void IQI::build_global_dof_and_dv_type_lists()
        {
            // MASTER-------------------------------------------------------
            // get number of dof and dv types on set
            uint tNumDofTypes = mSet->get_num_unique_dof_types();
            uint tNumDvTypes  = mSet->get_num_unique_dv_types();

            // set size for the global dof and dv type lists
            mMasterGlobalDofTypes.reserve( tNumDofTypes );
            mMasterGlobalDvTypes.reserve( tNumDvTypes );

            // set a size for the dof and dv checkLists
            //( used to avoid repeating a dof or a dv type)
            Matrix< DDSMat > tDofCheckList( tNumDofTypes, 1, -1 );
            Matrix< DDSMat > tDvCheckList( tNumDvTypes, 1, -1 );

            // get dof type from direct dependencies
            for ( uint iDof = 0; iDof < mMasterDofTypes.size(); iDof++ )
            {
                // get set index for dof type
                sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( mMasterDofTypes( iDof )( 0 ) );  //FIXME'

                // put the dof type in the checklist
                tDofCheckList( tDofTypeIndex ) = 1;

                // put the dof type in the global type list
                mMasterGlobalDofTypes.push_back( mMasterDofTypes( iDof ) );
            }

            // get dv type from direct dependencies
            for ( uint iDv = 0; iDv < mMasterDvTypes.size(); iDv++ )
            {
                // get set index for dv type
                sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( mMasterDvTypes( iDv )( 0 ) );  //FIXME'

                // put the dv type in the checklist
                tDvCheckList( tDvTypeIndex ) = 1;

                // put the dv type in the global type list
                mMasterGlobalDvTypes.push_back( mMasterDvTypes( iDv ) );
            }

            // get dof type from master properties
            for ( std::shared_ptr< Property > tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get dof types for property
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofTypes
                    = tProperty->get_dof_type_list();

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for property
                    moris::Cell< moris::Cell< PDV_Type > > tActiveDvTypes
                    = tProperty->get_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                        // if dof enum not in the list
                        if ( tDvCheckList( tDvTypeIndex) != 1 )
                        {
                            // put the dof type in the checklist
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
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
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofTypes
                    = tCM->get_global_dof_type_list();

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex) != 1 )
                        {
                            // put the dof type in the checklist
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for constitutive model
                    moris::Cell< moris::Cell< PDV_Type > > tActiveDvTypes
                    = tCM->get_global_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                        // if dv enum not in the list
                        if ( tDvCheckList( tDvTypeIndex ) != 1 )
                        {
                            // put the dv type in the checklist
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dv type in the global type list
                            mMasterGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
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
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofTypes
                    = tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex) != 1 )
                        {
                            // put the dof type in the checklist
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for constitutive model
                    moris::Cell< moris::Cell< PDV_Type > > tActiveDvTypes
                    = tSP->get_global_dv_type_list( mtk::Master_Slave::MASTER );

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                        // if dv enum not in the list
                        if ( tDvCheckList( tDvTypeIndex ) != 1 )
                        {
                            // put the dv type in the checklist
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dv type in the global type list
                            mMasterGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                        }
                    }
                }
            }

            // reduce size of dof and dv lists to fit unique list
            mMasterGlobalDofTypes.shrink_to_fit();
            mMasterGlobalDvTypes.shrink_to_fit();

            // SLAVE--------------------------------------------------------

            // set size for the global dof type list
            mSlaveGlobalDofTypes.reserve( tNumDofTypes );

            // set a size for the checkList ( used to avoid repeating a dof type)
            tDofCheckList.fill( -1 );
            tDvCheckList.fill( -1 );

            // get dof type from slave direct dependencies
            for ( uint iDof = 0; iDof < mSlaveDofTypes.size(); iDof++ )
            {
                // get set index for dof type
                sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( mSlaveDofTypes( iDof )( 0 ) );

                // put the dof type in the checklist
                tDofCheckList( tDofTypeIndex ) = 1;

                // put the dof type in the global type list
                mSlaveGlobalDofTypes.push_back( mSlaveDofTypes( iDof ) );
            }

            // get dv type from slave direct dependencies
            for ( uint iDv = 0; iDv < mSlaveDvTypes.size(); iDv++ )
            {
                // get set index for dv type
                sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( mSlaveDvTypes( iDv )( 0 ) );

                // put the dv type in the checklist
                tDvCheckList( tDvTypeIndex ) = 1;

                // put the dv type in the global type list
                mSlaveGlobalDvTypes.push_back( mSlaveDvTypes( iDv ) );
            }

            // get dof type from master properties
            for ( std::shared_ptr< Property > tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get dof types for property
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofTypes
                    = tProperty->get_dof_type_list();

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mSlaveGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for property
                    moris::Cell< moris::Cell< PDV_Type > > tActiveDvTypes
                    = tProperty->get_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                        // if dv enum not in the list
                        if ( tDvCheckList( tDvTypeIndex ) != 1 )
                        {
                            // put the dv type in the checklist
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dv type in the global type list
                            mSlaveGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
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
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofTypes
                    = tCM->get_global_dof_type_list();

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mSlaveGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for constitutive model
                    moris::Cell< moris::Cell< PDV_Type > > tActiveDvTypes
                    = tCM->get_global_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                        // if dv enum not in the list
                        if ( tDvCheckList( tDvTypeIndex ) != 1 )
                        {
                            // put the dv type in the checklist
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dv type in the global type list
                            mSlaveGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
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
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofTypes
                    = tSP->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mSlaveGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for stabilization parameter
                     moris::Cell< moris::Cell< PDV_Type > > tActiveDvTypes
                     = tSP->get_global_dv_type_list( mtk::Master_Slave::SLAVE );

                     // loop on property dv type
                     for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                     {
                         // get set index for dv type
                         sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                         // if dv enum not in the list
                         if ( tDvCheckList( tDvTypeIndex ) != 1 )
                         {
                             // put the dv type in the checklist
                             tDvCheckList( tDvTypeIndex ) = 1;

                             // put the dv type in the global type list
                             mSlaveGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                         }
                     }
                }
            }

            // reduce size of dof list to fit unique list
            mSlaveGlobalDofTypes.shrink_to_fit();
            mSlaveGlobalDvTypes.shrink_to_fit();
        }

//------------------------------------------------------------------------------
        void IQI::build_requested_dof_type_lists()
        {
            // clear the dof lists
            mRequestedMasterGlobalDofTypes.clear();
            mRequestedSlaveGlobalDofTypes .clear();

                // get the requested dof types
                Cell < enum MSI::Dof_Type > tRequestedDofTypes = mSet->get_requested_dof_types();

                // reserve possible max size for requested dof lists
                mRequestedMasterGlobalDofTypes.reserve( tRequestedDofTypes.size() );
                mRequestedSlaveGlobalDofTypes .reserve( tRequestedDofTypes.size() );

                // loop over the requested dof types
                for( auto tDofTypes : tRequestedDofTypes )
                {
                    // loop over the IWG master dof types groups
                    for ( uint Ik = 0; Ik < mMasterGlobalDofTypes.size(); Ik++ )
                    {
                        // if requested dof type matches IWG master dof type
                        if( mMasterGlobalDofTypes( Ik )( 0 ) == tDofTypes )
                        {
                            // add the IWG master dof type to the requested dof list
                            mRequestedMasterGlobalDofTypes.push_back( mMasterGlobalDofTypes( Ik ) );
                            break;
                        }
                    }

                    // loop over the IWG slave dof types groups
                    for ( uint Ik = 0; Ik < mSlaveGlobalDofTypes.size(); Ik++ )
                    {
                        // if requested dof type matches IWG slave dof type
                        if( mSlaveGlobalDofTypes( Ik )( 0 ) == tDofTypes )
                        {
                            // add the IWG slave dof type to the requested dof list
                            mRequestedSlaveGlobalDofTypes.push_back( mSlaveGlobalDofTypes( Ik ) );
                            break;
                        }
                    }
                }

            // reduce size for requested dof lists
            mRequestedMasterGlobalDofTypes.shrink_to_fit();
            mRequestedSlaveGlobalDofTypes.shrink_to_fit();
        }

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
        void IQI::compute_dQIdu_FD( Matrix< DDRMat > & adQIduFD,
                                    real               aWStar,
                                    real               aPerturbation )
        {
            // get master number of dof types
            uint tNumMasterDofType = mRequestedMasterGlobalDofTypes.size();
            uint tNumSlaveDofType  = mRequestedSlaveGlobalDofTypes.size();

            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get num of rows for dQIdu
            uint tNumRows = 0;

            // loop over the IQI master dof types
            for( uint iFI = 0; iFI < tNumMasterDofType; iFI++ )
            {
                // get master dof type index in set
                uint tDepIndex = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ),
                                                               mtk::Master_Slave::MASTER );

                // get num of column for the jacobian
                tNumRows += mSet->get_jac_dof_assembly_map()( 0 )( tDepIndex, 1 )
                          - mSet->get_jac_dof_assembly_map()( 0 )( tDepIndex, 0 ) + 1;
            }

            // loop over the IQI slave dof types
            for( uint iFI = 0; iFI < tNumSlaveDofType; iFI++ )
            {
                // get master dof type index in set
                uint tDepIndex = mSet->get_dof_index_for_type( mRequestedSlaveGlobalDofTypes( iFI )( 0 ),
                                                               mtk::Master_Slave::SLAVE );

                // get num of column for the jacobian
                tNumRows += mSet->get_jac_dof_assembly_map()( 0 )( tDepIndex, 1 )
                          - mSet->get_jac_dof_assembly_map()( 0 )( tDepIndex, 0 ) + 1;
            }

            // set size of dQIdu matrix
            adQIduFD.set_size( tNumRows, 1, 0.0 );

            // init dof counter
            uint tDofCounter = 0;

            // loop over the IWG dof types
            for( uint iFI = 0; iFI < tNumMasterDofType; iFI++ )
            {
                // get field interpolator for dependency dof type
                Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // loop over the coefficient column
                for( uint iCoeffCol = 0; iCoeffCol< tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficient row
                    for( uint iCoeffRow = 0; iCoeffRow< tDerNumBases; iCoeffRow++  )
                    {
                        // perturbation of the coefficent
                        Matrix< DDRMat > tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) += aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the residual
                        mSet->get_QI()( tQIIndex ).fill( 0.0 );
                        this->compute_QI( aWStar );
                        Matrix< DDRMat > tQI_Plus = mSet->get_QI()( tQIIndex );

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) -= aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the residual
                        mSet->get_QI()( tQIIndex ).fill( 0.0 );
                        this->compute_QI( aWStar );
                        Matrix< DDRMat > tQI_Minus = mSet->get_QI()( tQIIndex );

                        // evaluate Jacobian
                        adQIduFD( tDofCounter ) =
                                ( tQI_Plus( 0 ) - tQI_Minus( 0 ) ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // loop over the slave dof types
            for( uint iFI = 0; iFI < tNumSlaveDofType; iFI++ )
            {
                // get slave dependency field interpolator
                Field_Interpolator * tFI = mSlaveFIManager->get_field_interpolators_for_type( mRequestedSlaveGlobalDofTypes( iFI )( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // loop over the coefficients columns
                for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficients rows
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // perturbation of the coefficent
                        Matrix< DDRMat > tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) += aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the residual
                        mSet->get_residual()( tQIIndex ).fill( 0.0 );
                        this->compute_QI( aWStar );
                        Matrix< DDRMat > tQI_Plus = mSet->get_QI()( tQIIndex );

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) -= aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the residual
                        mSet->get_QI()( tQIIndex ).fill( 0.0 );
                        this->compute_QI( aWStar );

                        Matrix< DDRMat > tQI_Minus = mSet->get_QI()( tQIIndex );

                        // evaluate Jacobian
                        adQIduFD( tDofCounter ) = ( tQI_Plus( 0 ) - tQI_Minus( 0 ) )/ ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }
        }

        //------------------------------------------------------------------------------
        bool IQI::check_dQIdu_FD(
                real              aWStar,
                real              aPerturbation,
                real              aEpsilon,
                Matrix< DDRMat > & adQIdu,
                Matrix< DDRMat > & adQIduFD )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // compute dQIdDof with IQI
            this->compute_dQIdu( aWStar );
            adQIdu = mSet->get_residual()( tQIIndex );

            // compute dQIdDof by FD
            this->compute_dQIdu_FD( adQIduFD, aWStar, aPerturbation );

            //define a boolean for check
            bool tCheckdQIdDof = true;

            // check if adQIdu and adQIduFD have the same size
            tCheckdQIdDof = tCheckdQIdDof && ( adQIdu.n_rows() == adQIduFD.n_rows());
            tCheckdQIdDof = tCheckdQIdDof && ( adQIdu.n_cols() == adQIduFD.n_cols());

            // loop over the rows
            for ( uint iRow = 0; iRow < adQIdu.n_rows(); iRow++ )
            {
                // loop over the columns
                for( uint jCol = 0; jCol < adQIdu.n_cols(); jCol++ )
                {
                    // check each components
                    tCheckdQIdDof = tCheckdQIdDof && ( adQIdu( iRow, jCol ) - adQIduFD( iRow, jCol ) < aEpsilon );
                }
            }

            // return bool
            return tCheckdQIdDof;
        }

        //------------------------------------------------------------------------------
        void IQI::compute_dQIdp_FD_material(
                moris::real aWStar,
                moris::real aPerturbation )
        {
            // get the requested ip pdv types
            moris::Cell< moris::Cell< PDV_Type > > tRequestedPdvTypes;
            mSet->get_ip_dv_types_for_set( tRequestedPdvTypes );

            // get number of requested dv types
            uint tNumPdvType = tRequestedPdvTypes.size();

            // get the IQI index
            uint tIQIAssemblyIndex = 0;

            // loop over the pdv types
            for( uint iFI = 0; iFI < tNumPdvType; iFI++ )
            {
                // get the FI for the dv type
                Field_Interpolator * tFI
                = mMasterFIManager->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // init pdv coeff counter
                uint tPdvCoeffCounter = 0;

                // loop over the pdv coefficient column
                for( uint iCoeffCol = 0; iCoeffCol< tDerNumFields; iCoeffCol++ )
                {
                    // loop over the pdv coefficient row
                    for( uint iCoeffRow = 0; iCoeffRow< tDerNumBases; iCoeffRow++  )
                    {
                        // perturbation of the pdv coefficent
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
                        tCoeffPert( iCoeffRow, iCoeffCol ) -= aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the QI
                        Matrix< DDRMat > tQIValMinus;
                        this->compute_QI( tQIValMinus );

                        // get the pdv index for assembly
                        uint tPdvAssemblyIndex = mSet->get_mat_pdv_assembly_map()( iFI )( 0, 0 ) + tPdvCoeffCounter;

                        // evaluate dQIqp
                        mSet->get_dqidpmat()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                ( tQIValPlus( 0 ) - tQIValMinus( 0 ) ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update the pdv coefficient counter
                        tPdvCoeffCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }
        }

        //------------------------------------------------------------------------------
        void IQI::compute_dQIdp_FD_geometry(
                moris::real                       aWStar,
                moris::real                       aPerturbation,
                moris::Cell< Matrix< DDSMat > > & aIsActive,
                Matrix< IndexMat >              & aVertexIndices )
        {
            // get requested geometry pdv types
            moris::Cell< PDV_Type > tRequestedGeoPdvType;
            mSet->get_ig_unique_dv_types_for_set( tRequestedGeoPdvType );

            // get the IQI index
            uint tIQIAssemblyIndex = 0;

            // get the GI for the IP and IG element considered
            Geometry_Interpolator * tIPGI = mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();
            Geometry_Interpolator * tIGGI = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // get number of master GI bases and space dimensions
            uint tDerNumBases      = tIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tIGGI->get_number_of_space_dimensions();

            // coefficients for dv type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tIGGI->get_space_coeff();
            Matrix< DDRMat > tParamCoeff = tIGGI->get_space_param_coeff();
            Matrix< DDRMat > tEvaluationPoint;
            tIGGI->get_space_time( tEvaluationPoint );

            // loop over the spatial directions
            for( uint iCoeffCol = 0; iCoeffCol< tDerNumDimensions; iCoeffCol++ )
            {
                // loop over the IG nodes
                for( uint iCoeffRow = 0; iCoeffRow< tDerNumBases; iCoeffRow++ )
                {
                    if ( aIsActive( iCoeffCol )( iCoeffRow ) == 1 )
                    {
                        // perturbation of the coefficent
                        Matrix< DDRMat > tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) += aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tIGGI->set_space_coeff( tCoeffPert );

                        // update local coordinates
                        Matrix< DDRMat > tXCoords  = tCoeffPert.get_row( iCoeffRow );
                        Matrix< DDRMat > tXiCoords = tParamCoeff.get_row( iCoeffRow );
                        tIPGI->update_local_coordinates( tXCoords, tXiCoords );
                        Matrix< DDRMat > tParamCoeffPert = tParamCoeff;
                        tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();
                        tIGGI->set_space_param_coeff( tParamCoeffPert );

                        // set evaluation point for interpolators (FIs and GIs)
                        mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the QI
                        Matrix< DDRMat > tQIValPlus;
                        this->compute_QI( tQIValPlus );

                        // perturbation of the coefficient
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) -= aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tIGGI->set_space_coeff( tCoeffPert );

                        // update local coordinates
                        tXCoords  = tCoeffPert.get_row( iCoeffRow );
                        tXiCoords = tParamCoeff.get_row( iCoeffRow );
                        tIPGI->update_local_coordinates( tXCoords, tXiCoords );
                        tParamCoeffPert = tParamCoeff;
                        tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();
                        tIGGI->set_space_param_coeff( tParamCoeffPert );

                        // set evaluation point for interpolators (FIs and GIs)
                        mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the QI
                        Matrix< DDRMat > tQIValMinus;
                        this->compute_QI( tQIValMinus );

                        // get the geometry pdv assembly index
                        std::pair< moris_index, PDV_Type > tKeyPair =
                                std::make_pair( aVertexIndices( iCoeffRow ), tRequestedGeoPdvType( iCoeffCol ) );
                        uint tPdvAssemblyIndex = mSet->get_geo_pdv_assembly_map()[ tKeyPair ];

                        // evaluate dqidpdv
                        mSet->get_dqidpgeo()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                ( tQIValPlus( 0 ) - tQIValMinus( 0 ) ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );
                    }
                }
                // reset the coefficients values
                tIGGI->set_space_coeff( tCoeff );
                tIGGI->set_space_param_coeff( tParamCoeff );
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */


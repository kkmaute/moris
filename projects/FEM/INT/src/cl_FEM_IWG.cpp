/*
 * cl_FEM_IWG.cpp
 *
 *  Created on: Nov 12, 2019
 *      Author: sonne
 */

#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        void IWG::print_names()
        {
            std::cout<<"----------"<<std::endl;
            std::cout<<"IWG: "<<mName<<std::endl;

            // properties
            for( uint iProp = 0; iProp < mMasterProp.size(); iProp++ )
            {
                if( mMasterProp( iProp ) != nullptr )
                {
                    std::cout<<"Master property: "<<mMasterProp( iProp )->get_name()<<std::endl;
                }
            }
            for( uint iProp = 0; iProp < mSlaveProp.size(); iProp++ )
            {
                if( mSlaveProp( iProp ) != nullptr )
                {
                    std::cout<<"Slave property:  "<<mSlaveProp( iProp )->get_name()<<std::endl;
                }
            }

            // CM
            for( uint iCM = 0; iCM < mMasterCM.size(); iCM++ )
            {
                if( mMasterCM( iCM ) != nullptr )
                {
                    std::cout<<"Master CM:       "<<mMasterCM( iCM )->get_name()<<std::endl;
                }
            }
            for( uint iCM = 0; iCM < mSlaveCM.size(); iCM++ )
            {
                if( mSlaveCM( iCM ) != nullptr )
                {
                    std::cout<<"Slave CM:        "<<mSlaveCM( iCM )->get_name()<<std::endl;
                }
            }

            // SP
            for( uint iSP = 0; iSP < mStabilizationParam.size(); iSP++ )
            {
                if( mStabilizationParam( iSP ) != nullptr )
                {
                    std::cout<<"SP:              "<<mStabilizationParam( iSP )->get_name()<<std::endl;
                }
            }
            std::cout<<"----------"<<std::endl;
        }

//------------------------------------------------------------------------------
        void IWG::reset_eval_flags()
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
        void IWG::set_field_interpolator_manager( Field_Interpolator_Manager * aFieldInterpolatorManager,
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
        void IWG::set_field_interpolator_manager_previous_time( Field_Interpolator_Manager * aFieldInterpolatorManager,
                                                                mtk::Master_Slave            aIsMaster )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ) :
                {
                    mMasterPreviousFIManager = aFieldInterpolatorManager;
                    break;
                }

                default :
                {
                    MORIS_ERROR( false, "IWG::set_field_interpolator_manager - can only be master");
                    break;
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG::get_non_unique_dof_and_dv_types( moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                   moris::Cell< GEN_DV >        & aDvTypes )
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
                    moris::Cell< GEN_DV >        tActiveDvTypes;
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
                    moris::Cell< GEN_DV >        tActiveDvTypes;
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
                    moris::Cell< GEN_DV >        tActiveDvTypes;
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
                    moris::Cell< GEN_DV >        tActiveDvTypes;
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
                    moris::Cell< GEN_DV >         tActiveDvTypes;
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
                    moris::Cell< GEN_DV >        tActiveDvTypes;
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
                    moris::Cell< GEN_DV >        tActiveDvTypes;
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
                    moris::Cell< GEN_DV >        tActiveDvTypes;
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
                    moris::Cell< GEN_DV >        tActiveDvTypes;
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
                    moris::Cell< GEN_DV >        tActiveDvTypes;
                    tSP->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                          tActiveDvTypes );

                    // populate the dof and dv lists
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG::build_global_dof_and_dv_type_list()
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
                    moris::Cell< moris::Cell< GEN_DV > > tActiveDvTypes
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
                    moris::Cell< moris::Cell< GEN_DV > > tActiveDvTypes
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
                    moris::Cell< moris::Cell< GEN_DV > > tActiveDvTypes
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
                    moris::Cell< moris::Cell< GEN_DV > > tActiveDvTypes
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
                    moris::Cell< moris::Cell< GEN_DV > > tActiveDvTypes
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
                     moris::Cell< moris::Cell< GEN_DV > > tActiveDvTypes
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

void IWG::build_requested_dof_type_list( const bool aIsResidual )
{
    // clear the dof lists
    mRequestedMasterGlobalDofTypes.clear();
    mRequestedSlaveGlobalDofTypes .clear();

    // if residual evaluation
    if ( aIsResidual )
    {
        // get the requested dof types
        moris::Cell< moris::Cell< enum MSI::Dof_Type > > tRequestedDofTypes
        = mSet->get_secundary_dof_types();

        // reserve possible max size for requested dof lists
        mRequestedMasterGlobalDofTypes.reserve( tRequestedDofTypes.size() );
        mRequestedSlaveGlobalDofTypes .reserve( tRequestedDofTypes.size() );

        // loop over the requested dof type groups
        for( auto tDofTypes : tRequestedDofTypes )
        {
            // loop over the IWG master dof types groups
            for ( uint Ik = 0; Ik < mMasterGlobalDofTypes.size(); Ik++ )
            {
                // if requested dof type group matches IWG master dof types group
                if( mMasterGlobalDofTypes( Ik )( 0 ) == tDofTypes( 0 ) )
                {
                    // add the IWG master dof types group to the requested dof list
                    mRequestedMasterGlobalDofTypes.push_back( mMasterGlobalDofTypes( Ik ) );
                    break;
                }
            }

            // loop over the IWG slave dof types groups
            for ( uint Ik = 0; Ik < mSlaveGlobalDofTypes.size(); Ik++ )
            {
                // if requested dof type group matches IWG slave dof types group
                if( mSlaveGlobalDofTypes( Ik )( 0 ) == tDofTypes( 0 ) )
                {
                    // add the IWG slave dof types group to the requested dof list
                    mRequestedSlaveGlobalDofTypes.push_back( mSlaveGlobalDofTypes( Ik ) );
                    break;
                }
            }
        }
    }
    // if jacobian evaluation
    else
    {
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
    }

    // reduce size for requested dof lists
    mRequestedMasterGlobalDofTypes.shrink_to_fit();
    mRequestedSlaveGlobalDofTypes.shrink_to_fit();
}

//------------------------------------------------------------------------------
        void IWG::check_field_interpolators( mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case ( mtk::Master_Slave::MASTER ) :
                {
                    // loop over the dof field interpolator pointers
                    for( uint iDofFI = 0; iDofFI < mRequestedMasterGlobalDofTypes.size(); iDofFI++ )
                    {
                        // check that the field interpolator was set
                        MORIS_ASSERT( this->get_field_interpolator_manager( aIsMaster )
                                          ->get_field_interpolators_for_type( mRequestedMasterGlobalDofTypes( iDofFI )( 0 ) ) != nullptr,
                                      "IWG::check_field_interpolators - Master dof FI missing. " );
                    }

                    // loop over the dv field interpolator pointers
                    for( uint iDvFI = 0; iDvFI < mMasterGlobalDvTypes.size(); iDvFI++ )
                    {
                        // check that the field interpolator was set
                        MORIS_ASSERT( this->get_field_interpolator_manager( aIsMaster )
                                          ->get_field_interpolators_for_type( mMasterGlobalDvTypes( iDvFI )( 0 ) ) != nullptr,
                                      "IWG::check_field_interpolators - Master dv FI missing. " );
                    }
                    break;
                }
                case ( mtk::Master_Slave::SLAVE ) :
                {
                    // loop over the dof field interpolator pointers
                    for( uint iDofFI = 0; iDofFI < mRequestedSlaveGlobalDofTypes.size(); iDofFI++ )
                    {
                        // check that the field interpolator was set
                        MORIS_ASSERT( this->get_field_interpolator_manager( aIsMaster )
                                          ->get_field_interpolators_for_type( mRequestedSlaveGlobalDofTypes( iDofFI )( 0 ) ) != nullptr,
                                      "IWG::check_dof_field_interpolators - Slave dof FI missing. " );
                    }

                    // loop over the dv field interpolator pointers
                    for( uint iDvFI = 0; iDvFI < mSlaveGlobalDvTypes.size(); iDvFI++ )
                    {
                        // check that the field interpolator was set
                        MORIS_ASSERT( this->get_field_interpolator_manager( aIsMaster )
                                          ->get_field_interpolators_for_type( mSlaveGlobalDvTypes( iDvFI )( 0 ) ) != nullptr,
                                      "IWG::check_field_interpolators - Slave dv FI missing. " );
                    }
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "IWG::check_field_interpolators - can only be master or slave." );
                    break;
                }
            }
        }

//------------------------------------------------------------------------------
        moris::Cell< moris::Cell< MSI::Dof_Type > > &
        IWG::get_global_dof_type_list( mtk::Master_Slave aIsMaster )
        {
            // if the global list was not yet built
            if( mGlobalDofBuild )
            {
                // build the stabilization parameter global dof type list
                this->build_global_dof_and_dv_type_list();

                // update build flag
                mGlobalDofBuild = false;
                mGlobalDvBuild  = false;
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
        }

//------------------------------------------------------------------------------
        moris::Cell< moris::Cell< GEN_DV > > &
        IWG::get_global_dv_type_list( mtk::Master_Slave aIsMaster )
        {
            // if the global list was not yet built
            if( mGlobalDvBuild )
            {
                // build the stabilization parameter global dof type list
                this->build_global_dof_and_dv_type_list();

                // update build flag
                mGlobalDofBuild = false;
                mGlobalDvBuild  = false;
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
        }

//------------------------------------------------------------------------------
        void IWG::compute_jacobian_FD( real               aWStar,
                                       real               aPerturbation,
                                       Matrix< DDRMat > & aJacobiansFD )
        {
            // get master number of dof types
            uint tNumDofType = mRequestedMasterGlobalDofTypes.size();

            // get residual dof type index in set
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get start and end indices for residual dof type
            uint tResStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get num of rows for jacobian
            uint tNumRows = tResEndRow - tResStartRow + 1;

            // get num of cols for jacobian
            uint tNumCols = 0;

            // loop over the IWG dof types
            for( uint iFI = 0; iFI < tNumDofType; iFI++ )
            {
                // get dependency dof type index in set
                uint tDepIndex = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::MASTER );

                // get num of column for the jacobian
                tNumCols += mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepIndex, 1 )
                          - mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepIndex, 0 ) + 1;
            }

            // set size of jacobian matrix
            aJacobiansFD.set_size( tNumRows, tNumCols, 0.0 );

            // init dof counter
            uint tDofCounter = 0;

            // loop over the IWG dof types
            for( uint iFI = 0; iFI < tNumDofType; iFI++ )
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
                        mSet->get_residual()( 0 ).fill( 0.0 );
                        this->compute_residual( aWStar );

                        Matrix< DDRMat > tResidual_Plus
                        =  mSet->get_residual()( 0 )( { tResStartRow, tResEndRow }, { 0, 0 } );

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) += - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the residual
                        mSet->get_residual()( 0 ).fill( 0.0 );
                        this->compute_residual( aWStar );

                        Matrix< DDRMat > tResidual_Minus
                        =  mSet->get_residual()( 0 )( { tResStartRow, tResEndRow }, { 0, 0 } );

                        // evaluate Jacobian
                        aJacobiansFD.get_column( tDofCounter )
                        = ( tResidual_Plus - tResidual_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }
        }

//------------------------------------------------------------------------------
        void IWG::compute_jacobian_FD_double( real               aWStar,
                                              real               aPerturbation,
                                              Matrix< DDRMat > & aJacobiansFD )
        {
            // get master and slave number of dof types
            uint tMasterNumDofType = mRequestedMasterGlobalDofTypes.size();
            uint tSlaveNumDofType  = mRequestedSlaveGlobalDofTypes.size();

            // get residual dof type index in set
            uint tMasterResDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tSlaveResDofIndex  = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

            // get start and end indices for residual dof type
            uint tMasterResStartRow = mSet->get_res_dof_assembly_map()( tMasterResDofIndex )( 0, 0 );
            uint tMasterResEndRow   = mSet->get_res_dof_assembly_map()( tMasterResDofIndex )( 0, 1 );
            uint tSlaveResStartRow  = mSet->get_res_dof_assembly_map()( tSlaveResDofIndex )( 0, 0 );
            uint tSlaveResEndRow    = mSet->get_res_dof_assembly_map()( tSlaveResDofIndex )( 0, 1 );

            // get num of rows for jacobian
            uint tMasterNumRows = ( tMasterResEndRow - tMasterResStartRow ) + 1;
            uint tSlaveNumRows  = ( tSlaveResEndRow  - tSlaveResStartRow  ) + 1;

            // get num of cols for jacobian
            uint tNumCols = 0;

            // loop over the IWG master dof types
            for( uint iFI = 0; iFI < tMasterNumDofType; iFI++ )
            {
                // get dependency dof type index in set
                uint tDepIndex = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::MASTER );

                // get num of column for the jacobian
                tNumCols += mSet->get_jac_dof_assembly_map()( tMasterResDofIndex )( tDepIndex, 1 )
                          - mSet->get_jac_dof_assembly_map()( tMasterResDofIndex )( tDepIndex, 0 ) + 1;
            }
            // loop over the IWG slave dof types
            for( uint iFI = 0; iFI < tSlaveNumDofType; iFI++ )
            {
                // get dependency dof type index in set
                uint tDepIndex = mSet->get_dof_index_for_type( mRequestedSlaveGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::SLAVE );

                // get num of column for the jacobian
                tNumCols += mSet->get_jac_dof_assembly_map()( tSlaveResDofIndex )( tDepIndex, 1 )
                          - mSet->get_jac_dof_assembly_map()( tSlaveResDofIndex )( tDepIndex, 0 ) + 1;
            }

            // set size of jacobian matrix
            aJacobiansFD.set_size( tMasterNumRows + tSlaveNumRows, tNumCols, 0.0 );

            // init dof counter
            uint tDofCounter = 0;

            // loop over the master dof types
            for( uint iFI = 0; iFI < tMasterNumDofType; iFI++ )
            {
                // get master dependency field interpolator
                Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // loop over the coefficients column
                for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficients row
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
                        mSet->get_residual()( 0 ).fill( 0.0 );
                        this->compute_residual( aWStar );

                        Matrix< DDRMat > tResidual_Plus_Master
                        =  mSet->get_residual()( 0 )( { tMasterResStartRow, tMasterResEndRow }, { 0, 0 } );
                        Matrix< DDRMat > tResidual_Plus_Slave
                        =  mSet->get_residual()( 0 )( { tSlaveResStartRow, tSlaveResEndRow }, { 0, 0 } );

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) -= aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the residual
                        mSet->get_residual()( 0 ).fill( 0.0 );
                        this->compute_residual( aWStar );

                        Matrix< DDRMat > tResidual_Minus_Master
                        =  mSet->get_residual()( 0 )( { tMasterResStartRow, tMasterResEndRow }, { 0, 0 } );
                        Matrix< DDRMat > tResidual_Minus_Slave
                        =  mSet->get_residual()( 0 )( { tSlaveResStartRow, tSlaveResEndRow }, { 0, 0 } );

                        // evaluate Jacobian
                        aJacobiansFD( { 0, tMasterNumRows -1 }, { tDofCounter, tDofCounter } )
                        = ( tResidual_Plus_Master - tResidual_Minus_Master )/ ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );
                        aJacobiansFD( { tMasterNumRows, tMasterNumRows + tSlaveNumRows - 1 }, { tDofCounter, tDofCounter } )
                        = ( tResidual_Plus_Slave  - tResidual_Minus_Slave  )/ ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // loop over the slave dof types
            for( uint iFI = 0; iFI < tSlaveNumDofType; iFI++ )
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
                        mSet->get_residual()( 0 ).fill( 0.0 );
                        this->compute_residual( aWStar );

                        Matrix< DDRMat > tResidual_Plus_Master
                        = mSet->get_residual()( 0 )( { tMasterResStartRow, tMasterResEndRow }, { 0, 0 } );
                        Matrix< DDRMat > tResidual_Plus_Slave
                        = mSet->get_residual()( 0 )( { tSlaveResStartRow, tSlaveResEndRow }, { 0, 0 } );

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) -= aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the residual
                        mSet->get_residual()( 0 ).fill( 0.0 );
                        this->compute_residual( aWStar );

                        Matrix< DDRMat > tResidual_Minus_Master
                        = mSet->get_residual()( 0 )( { tMasterResStartRow, tMasterResEndRow }, { 0, 0 } );
                        Matrix< DDRMat > tResidual_Minus_Slave
                        = mSet->get_residual()( 0 )( { tSlaveResStartRow, tSlaveResEndRow }, { 0, 0 } );

                        // evaluate Jacobian
                        aJacobiansFD( { 0, tMasterNumRows -1 }, { tDofCounter, tDofCounter } )
                        = ( tResidual_Plus_Master - tResidual_Minus_Master )/ ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );
                        aJacobiansFD( { tMasterNumRows, tMasterNumRows + tSlaveNumRows - 1 }, { tDofCounter, tDofCounter } )
                        = ( tResidual_Plus_Slave  - tResidual_Minus_Slave  )/ ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }
        }

//------------------------------------------------------------------------------
        bool IWG::check_jacobian( real               aPerturbation,
                                  real               aEpsilon,
                                  real               aWStar,
                                  Matrix< DDRMat > & aJacobians,
                                  Matrix< DDRMat > & aJacobiansFD )
        {
            // compute jacobian with IWG
            this->compute_jacobian( aWStar );

            // get residual dof type index in set
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get start and end indices for residual dof type
            uint tResStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get number of cols for jacobian
            uint tNumCols = mSet->get_jacobian().n_cols();

            // get the computed jacobian
            aJacobians = mSet->get_jacobian()( { tResStartRow, tResEndRow }, { 0, tNumCols - 1 });

            // compute jacobian by FD
            this->compute_jacobian_FD( aWStar, aPerturbation, aJacobiansFD );

            //define a boolean for check
            bool tCheckJacobian = true;

            // define a real for absolute difference
            real tAbsolute = 0.0;

            // define a real for relative difference
            real tRelative = 0.0;

            for( uint iiJac = 0; iiJac < aJacobians.n_rows(); iiJac++ )
            {
                for( uint jjJac = 0; jjJac < aJacobians.n_cols(); jjJac++ )
                {
                    // get absolute difference
                    tAbsolute = std::abs( aJacobians( iiJac, jjJac ) - aJacobiansFD( iiJac, jjJac ) );

                    // get relative difference
                    tRelative = std::abs( ( aJacobiansFD( iiJac, jjJac ) - aJacobians( iiJac, jjJac ) ) / aJacobiansFD( iiJac, jjJac ) );

                    // update check value
                    tCheckJacobian = tCheckJacobian && ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) );

//                    // for debug
//                    if( !( tAbsolute < aEpsilon ) && ( tRelative < aEpsilon ) )
//                    {
//                        std::cout<<"iiJac "<<iiJac<<" - jjJac "<<jjJac<<std::endl;
//                        std::cout<<"aJacobians( iiJac, jjJac ) "<<aJacobians( iiJac, jjJac )<<std::endl;
//                        std::cout<<"aJacobiansFD( iiJac, jjJac ) "<<aJacobiansFD( iiJac, jjJac )<<std::endl;
//                        std::cout<<"Absolute difference "<<tAbsolute<<std::endl;
//                        std::cout<<"Relative difference "<<tRelative<<std::endl;
//                    }
                }
            }

            // return bool
            return tCheckJacobian;
        }

//------------------------------------------------------------------------------
        bool IWG::check_jacobian_double( real               aPerturbation,
                                         real               aEpsilon,
                                         real               aWStar,
                                         Matrix< DDRMat > & aJacobians,
                                         Matrix< DDRMat > & aJacobiansFD )
        {
            // compute jacobian with IWG
            this->compute_jacobian( aWStar );

            // get residual dof type index in set
            uint tMasterDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tSlaveDofIndex  = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

            // get start and end indices for residual dof type
            uint tMasterResStartRow = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResEndRow   = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );
            uint tSlaveResStartRow  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResEndRow    = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get number of master and slave rows
            uint tMasterNumRows = tMasterResEndRow - tMasterResStartRow + 1;
            uint tSlaveNumRows  = tSlaveResEndRow  - tSlaveResStartRow  + 1;

            // get number of cols for jacobian
            uint tNumCols = mSet->get_jacobian().n_cols();

            // get the computed jacobian
            aJacobians.set_size( tMasterNumRows + tSlaveNumRows, tNumCols, 0.0 );

            // get the computed jacobian
            aJacobians( { 0, tMasterNumRows - 1 }, { 0, tNumCols - 1 } )
            = mSet->get_jacobian()( { tMasterResStartRow, tMasterResEndRow }, { 0, tNumCols - 1 });
            aJacobians( { tMasterNumRows, tMasterNumRows + tSlaveNumRows - 1 }, { 0, tNumCols - 1 } )
            = mSet->get_jacobian()( { tSlaveResStartRow, tSlaveResEndRow }, { 0, tNumCols - 1 });

            // compute jacobian by FD
            this->compute_jacobian_FD_double( aWStar, aPerturbation, aJacobiansFD );

            //define a boolean for check
            bool tCheckJacobian = true;

            // loop over the rows of jacobian
            for( uint iiJac = 0; iiJac < aJacobians.n_rows(); iiJac++ )
            {
                // loop over the columns of jacobian
                for( uint jjJac = 0; jjJac < aJacobians.n_cols(); jjJac++ )
                {
                    // get absolute difference
                    real tAbsolute = std::abs( aJacobians( iiJac, jjJac ) - aJacobiansFD( iiJac, jjJac ) );

                    // get relative difference
                    real tRelative = tAbsolute / aJacobiansFD( iiJac, jjJac );

                    // update check value
                    tCheckJacobian = tCheckJacobian && ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) );

//                    // for debug
//                    if( !( ( tAbsolute < aEpsilon ) && ( tRelative < aEpsilon ) ) )
//                    {
//                        std::cout<<"iiJac "<<iiJac<<" - jjJac "<<jjJac<<std::endl;
//                        std::cout<<"aJacobians( iiJac, jjJac ) "<<aJacobians( iiJac, jjJac )<<std::endl;
//                        std::cout<<"aJacobiansFD( iiJac, jjJac ) "<<aJacobiansFD( iiJac, jjJac )<<std::endl;
//                        std::cout<<"Absolute difference "<<tAbsolute<<std::endl;
//                        std::cout<<"Relative difference "<<tRelative<<std::endl;
//                    }
                }
            }

            // return bool
            return tCheckJacobian;
        }

//------------------------------------------------------------------------------
        void IWG::compute_dRdp_FD_geometry
        ( moris::real                       aWStar,
          moris::real                       aPerturbation,
          moris::Cell< Matrix< DDSMat > > & aIsActive,
          Matrix< IndexMat >              & aVertexIndices )
        {
            // get the GI for the IG element considered
            Geometry_Interpolator * tGI = mSet->get_field_interpolator_manager()
                                              ->get_IG_geometry_interpolator();

            // get the residual dof type index in the set
            uint tResDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ),
                                                              mtk::Master_Slave::MASTER );
            uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
            uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

            // get number of master GI bases and space dimensions
            uint tDerNumBases      = tGI->get_number_of_space_bases();
            uint tDerNumDimensions = tGI->get_number_of_space_dimensions();

            // coefficients for dv type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tGI->get_space_coeff();

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
                        tGI->set_space_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the residual plus
                        mSet->get_residual()( 0 ).fill( 0.0 );
                        this->compute_residual( aWStar );
                        Matrix< DDRMat > tResidual_Plus
                        =  mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } );

                        // perturbation of the coefficient
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) -= aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tGI->set_space_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the residual minus
                        mSet->get_residual()( 0 ).fill( 0.0 );
                        this->compute_residual( aWStar );
                        Matrix< DDRMat > tResidual_Minus
                        = mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } );

                        // FIXME get geo pdv index
                        moris::Cell< GEN_DV > tGeoPdvType = { GEN_DV::XCOORD, GEN_DV::YCOORD, GEN_DV::ZCOORD };

                        std::pair< moris_index, GEN_DV > tKeyPair
                        = std::make_pair( aVertexIndices( iCoeffRow ), tGeoPdvType( iCoeffCol ) );
                        uint tPdvIndex = mSet->get_geo_pdv_assembly_map()[ tKeyPair ];

                        // evaluate dRdpGeo
                        Matrix< DDRMat > tI( tResidual_Plus.numel(), 1, 1.0 );
                        mSet->get_drdpgeo()( { tResDofAssemblyStart, tResDofAssemblyStop }, { tPdvIndex, tPdvIndex } )
                        += ( tResidual_Plus - tResidual_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );
                    }
                }
                // reset the coefficients values
                tGI->set_space_coeff( tCoeff );
            }
        }

//------------------------------------------------------------------------------
        void IWG::compute_dRdp_FD_material( moris::real aWStar,
                                            moris::real aPerturbation )
        {
            // FIXME get master number of dv types
            uint tNumDvType = mMasterGlobalDvTypes.size();

            // get the residual dof type index in the set
            uint tResDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ),
                                                              mtk::Master_Slave::MASTER );
            uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
            uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

            // loop over the dv types associated with a FI
            for( uint iFI = 0; iFI < tNumDvType; iFI++ )
            {
                // get the FI for the dv type
                Field_Interpolator * tFI
                = mMasterFIManager->get_field_interpolators_for_type( mMasterGlobalDvTypes( iFI )( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // init coeff counter
                uint tCoeffCounter = 0;

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

                        // evaluate the residual plus
                        mSet->get_residual()( 0 ).fill( 0.0 );
                        this->compute_residual( aWStar );
                        Matrix< DDRMat > tResidual_Plus
                        =  mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } );

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) -= aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the residual minus
                        mSet->get_residual()( 0 ).fill( 0.0 );
                        this->compute_residual( aWStar );
                        Matrix< DDRMat > tResidual_Minus
                        =  mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } );

                        // get mat pdv index
                        uint tPdvIndex = mSet->get_mat_pdv_assembly_map()( iFI )( 0, 0 ) + tCoeffCounter;

                        // evaluate dRdpMat
//                        std::cout<<"tResDofAssemblyStart "<<tResDofAssemblyStart<<std::endl;
//                        std::cout<<"tResDofAssemblyStop "<<tResDofAssemblyStop<<std::endl;
//                        std::cout<<"tPdvIndex "<<tPdvIndex<<std::endl;
//                        print(mSet->get_drdpmat()( 0 ),"mSet->get_drdpmat()( 0 )");

                        mSet->get_drdpmat()( { tResDofAssemblyStart, tResDofAssemblyStop }, { tPdvIndex, tPdvIndex } )
                        += ( tResidual_Plus - tResidual_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                        // update coefficient counter
                        tCoeffCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }
        }

//------------------------------------------------------------------------------

}   // end fem namespace
}   // end moris namespace



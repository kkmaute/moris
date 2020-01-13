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
            /*
             * set field interpolator manager
             * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
             * @param[ in ] aIsMaster                 an enum for master or slave
             */
            void IWG::set_field_interpolator_manager( Field_Interpolator_Manager * aFieldInterpolatorManager,
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
        void IWG::get_non_unique_dof_types( moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // set the size of the dof type list for the set
            uint tCounter = 0;

            // get number of direct master dof dependencies
            for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
            {
                tCounter += mMasterDofTypes( iDOF ).size();
            }

            // get number of direct slave dof dependencies
            for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
            {
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
                    // get SP non unique dof type list
                    moris::Cell< MSI::Dof_Type >  tActiveDofType;
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
                    // get SP non unique master dof type list
                    moris::Cell< MSI::Dof_Type >  tActiveDofType;
                    tSP->get_non_unique_dof_types( tActiveDofType );

                    // populate the dof list
                    aDofTypes.append( tActiveDofType );
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG::get_non_unique_dof_and_dv_types( moris::Cell< MSI::Dof_Type > & aDofTypes,
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
        void IWG::build_global_dof_type_list()
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

//------------------------------------------------------------------------------
        void IWG::build_global_dof_and_dv_type_list()
        {
            // MASTER-------------------------------------------------------
            // get number of dof and dv types on set
            uint tNumDofTypes = mSet->get_num_dof_types();
            uint tNumDvTypes  = mSet->get_num_dv_types();

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
                sint tDofTypeIndex = mSet->get_dof_index_for_type_1( mMasterDofTypes( iDof )( 0 ) );  //FIXME'

                // put the dof type in the checklist
                tDofCheckList( tDofTypeIndex ) = 1;

                // put the dof type in the global type list
                mMasterGlobalDofTypes.push_back( mMasterDofTypes( iDof ) );
            }

            // get dv type from direct dependencies
            for ( uint iDv = 0; iDv < mMasterDvTypes.size(); iDv++ )
            {
                // get set index for dv type
                sint tDvTypeIndex = mSet->get_dv_index_for_type_1( mMasterDvTypes( iDv )( 0 ) );  //FIXME'

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
                        sint tDofTypeIndex = mSet->get_dof_index_for_type_1( tActiveDofTypes( iDof )( 0 ) );

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
                    moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvTypes
                    = tProperty->get_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_dv_index_for_type_1( tActiveDvTypes( iDv )( 0 ) );

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
                        sint tDofTypeIndex = mSet->get_dof_index_for_type_1( tActiveDofTypes( iDof )( 0 ) );

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
                    moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvTypes
                    = tCM->get_global_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_dv_index_for_type_1( tActiveDvTypes( iDv )( 0 ) );

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
                        sint tDofTypeIndex = mSet->get_dof_index_for_type_1( tActiveDofTypes( iDof )( 0 ) );

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
                    moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvTypes
                    = tSP->get_global_dv_type_list( mtk::Master_Slave::MASTER );

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_dv_index_for_type_1( tActiveDvTypes( iDv )( 0 ) );

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
                sint tDofTypeIndex = mSet->get_dof_index_for_type_1( mSlaveDofTypes( iDof )( 0 ) );

                // put the dof type in the checklist
                tDofCheckList( tDofTypeIndex ) = 1;

                // put the dof type in the global type list
                mSlaveGlobalDofTypes.push_back( mSlaveDofTypes( iDof ) );
            }

            // get dv type from slave direct dependencies
            for ( uint iDv = 0; iDv < mSlaveDvTypes.size(); iDv++ )
            {
                // get set index for dv type
                sint tDvTypeIndex = mSet->get_dv_index_for_type_1( mSlaveDvTypes( iDv )( 0 ) );

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
                        sint tDofTypeIndex = mSet->get_dof_index_for_type_1( tActiveDofTypes( iDof )( 0 ) );

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
                    moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvTypes
                    = tProperty->get_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_dv_index_for_type_1( tActiveDvTypes( iDv )( 0 ) );

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
                        sint tDofTypeIndex = mSet->get_dof_index_for_type_1( tActiveDofTypes( iDof )( 0 ) );

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
                    moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvTypes
                    = tCM->get_global_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_dv_index_for_type_1( tActiveDvTypes( iDv )( 0 ) );

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
                        sint tDofTypeIndex = mSet->get_dof_index_for_type_1( tActiveDofTypes( iDof )( 0 ) );

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
                     moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvTypes
                     = tSP->get_global_dv_type_list( mtk::Master_Slave::SLAVE );

                     // loop on property dv type
                     for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                     {
                         // get set index for dv type
                         sint tDvTypeIndex = mSet->get_dv_index_for_type_1( tActiveDvTypes( iDv )( 0 ) );

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

void IWG::build_requested_dof_type_list( const bool aItResidual )
{
    mRequestedMasterGlobalDofTypes.clear();
    mRequestedSlaveGlobalDofTypes .clear();

    if ( aItResidual )
    {
        moris::Cell< moris::Cell< enum MSI::Dof_Type > > tRequestedDofTypes
        = mSet->get_secundary_dof_types();

        mRequestedMasterGlobalDofTypes.reserve( tRequestedDofTypes.size() );
        mRequestedSlaveGlobalDofTypes .reserve( tRequestedDofTypes.size() );

        for( auto tDofTypes : tRequestedDofTypes )
        {
            for ( uint Ik = 0; Ik < mMasterGlobalDofTypes.size(); Ik++ )
            {
                if( mMasterGlobalDofTypes( Ik )( 0 ) == tDofTypes( 0 ) )
                {
                    mRequestedMasterGlobalDofTypes.push_back( mMasterGlobalDofTypes( Ik ) );

                    break;
                }
            }

            for ( uint Ik = 0; Ik < mSlaveGlobalDofTypes.size(); Ik++ )
            {
                if( mSlaveGlobalDofTypes( Ik )( 0 ) == tDofTypes( 0 ) )
                {
                    mRequestedSlaveGlobalDofTypes.push_back( mSlaveGlobalDofTypes( Ik ) );

                    break;
                }
            }
        }
    }
    else
    {
        Cell < enum MSI::Dof_Type > tRequestedDofTypes = mSet->get_requested_dof_types();

        mRequestedMasterGlobalDofTypes.reserve( tRequestedDofTypes.size() );
        mRequestedSlaveGlobalDofTypes .reserve( tRequestedDofTypes.size() );

        for( auto tDofTypes : tRequestedDofTypes )
        {
            for ( uint Ik = 0; Ik < mMasterGlobalDofTypes.size(); Ik++ )
            {
                if( mMasterGlobalDofTypes( Ik )( 0 ) == tDofTypes )
                {
                    mRequestedMasterGlobalDofTypes.push_back( mMasterGlobalDofTypes( Ik ) );
                    break;
                }
            }

            for ( uint Ik = 0; Ik < mSlaveGlobalDofTypes.size(); Ik++ )
            {
                if( mSlaveGlobalDofTypes( Ik )( 0 ) == tDofTypes )
                {
                    mRequestedSlaveGlobalDofTypes.push_back( mSlaveGlobalDofTypes( Ik ) );

                    break;
                }
            }

            if( mResidualDofType( 0 ) == tDofTypes )
            {
                mResidualDofTypeRequested = true;
            }
        }
    }

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

void IWG::compute_jacobian_FD( real                                             aWStar,
                               real                                             aPerturbation,
                               moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFD )
{
    // get master number of dof types
    uint tNumDofType = mRequestedMasterGlobalDofTypes.size();

    aJacobiansFD.resize( 1 );
    aJacobiansFD( 0 ).resize( tNumDofType );

    // loop over the IWG dof types
    for( uint iFI = 0; iFI < tNumDofType; iFI++ )
    {
        uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
        uint tDepIndex = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::MASTER );

        uint tNumRows = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) - mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ) + 1;
        uint tNumCols = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepIndex, 1 ) - mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepIndex, 0 ) + 1;

        aJacobiansFD( 0 )( iFI ).set_size( tNumRows, tNumCols, 0.0 );

        Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ) );

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
                tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) + aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset properties, CM and SP for IWG
                this->reset_eval_flags();

                // evaluate the residual
                mSet->get_residual().fill( 0.0 );
                this->compute_residual( aWStar );

                Matrix< DDRMat > tResidual_Plus
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } );

                // perturbation of the coefficent
                tCoeffPert = tCoeff;
                tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset properties, CM and SP for IWG
                this->reset_eval_flags();

                // evaluate the residual
                mSet->get_residual().fill( 0.0 );
                this->compute_residual( aWStar );

                Matrix< DDRMat > tResidual_Minus
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } );

                // evaluate Jacobian
                aJacobiansFD( 0 )( iFI ).get_column( tDofCounter )
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

 void IWG::compute_jacobian_FD_double( real                                             aWStar,
                                       real                                             aPerturbation,
                                       moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFD )
{
    // get master and slave number of dof types
    uint tMasterNumDofType = mRequestedMasterGlobalDofTypes.size();
    uint tSlaveNumDofType  = mRequestedSlaveGlobalDofTypes.size();

    aJacobiansFD.resize( 2 );
    aJacobiansFD( 0 ).resize( tMasterNumDofType + tSlaveNumDofType );
    aJacobiansFD( 1 ).resize( tMasterNumDofType + tSlaveNumDofType );

    // loop over the master dof types
    for( uint iFI = 0; iFI < tMasterNumDofType; iFI++ )
    {
        uint tDofIndexMaster = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::MASTER );
        uint tDofIndexSlave  = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::SLAVE  );

        uint tNumRowsMaster = mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) - mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ) + 1;
        uint tNumColsMaster = mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 1 ) - mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 0 ) + 1;
        uint tNumRowsSlave  = mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0 , 1 ) - mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0 , 0 ) + 1;
        uint tNumColsSlave  = mSet->get_jac_dof_assembly_map()( tDofIndexSlave  )( tDofIndexSlave , 1 ) - mSet->get_jac_dof_assembly_map()( tDofIndexSlave  )( tDofIndexSlave , 0 ) + 1;

        aJacobiansFD( 0 )( iFI ).set_size( tNumRowsMaster, tNumColsMaster, 0.0 );
        aJacobiansFD( 1 )( iFI ).set_size( tNumRowsSlave, tNumColsSlave, 0.0 );

        Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ) );

        // get number of master FI bases and fields
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

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
                tFI->set_coeff( tCoeffPert );

                // reset properties, CM and SP for IWG
                this->reset_eval_flags();

                // evaluate the residual
                mSet->get_residual().fill( 0.0 );
                this->compute_residual( aWStar );

                Matrix< DDRMat > tResidual_Plus_Master
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } );
                Matrix< DDRMat > tResidual_Plus_Slave
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 1 ) }, { 0, 0 } );

                // perturbation of the coefficent
                tCoeffPert = tCoeff;
                tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset properties, CM and SP for IWG
                this->reset_eval_flags();

                // evaluate the residual
                mSet->get_residual().fill( 0.0 );
                this->compute_residual( aWStar );

                Matrix< DDRMat > tResidual_Minus_Master
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } );
                Matrix< DDRMat > tResidual_Minus_Slave
                =  mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 1 ) }, { 0, 0 } );

                // evaluate Jacobian
                aJacobiansFD( 0 )( iFI ).get_column( tDofCounter )
                = ( tResidual_Plus_Master - tResidual_Minus_Master )/ ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );
                aJacobiansFD( 1 )( iFI ).get_column( tDofCounter )
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
        uint tDofIndexMaster = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::MASTER );
        uint tDofIndexSlave  = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ), mtk::Master_Slave::SLAVE  );

        uint tNumRowsMaster = mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) - mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ) + 1;
        uint tNumColsMaster = mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 1 ) - mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 0 ) + 1;
        uint tNumRowsSlave  = mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0 , 1 ) - mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0 , 0 ) + 1;
        uint tNumColsSlave  = mSet->get_jac_dof_assembly_map()( tDofIndexSlave  )( tDofIndexSlave , 1 ) - mSet->get_jac_dof_assembly_map()( tDofIndexSlave  )( tDofIndexSlave , 0 ) + 1;

        aJacobiansFD( 0 )( tMasterNumDofType + iFI ).set_size( tNumRowsMaster, tNumColsMaster, 0.0 );
        aJacobiansFD( 1 )( tMasterNumDofType + iFI ).set_size( tNumRowsSlave, tNumColsSlave, 0.0 );

        Field_Interpolator * tFI = mSlaveFIManager->get_field_interpolators_for_type( mRequestedMasterGlobalDofTypes( iFI )( 0 ) );

        // get number of master FI bases and fields
        uint tDerNumBases  = tFI->get_number_of_space_time_bases();
        uint tDerNumFields = tFI->get_number_of_fields();

        // coefficients for dof type wrt which derivative is computed
        Matrix< DDRMat > tCoeff = tFI->get_coeff();

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
                tFI->set_coeff( tCoeffPert );

                // reset properties, CM and SP for IWG
                this->reset_eval_flags();

                // evaluate the residual
                mSet->get_residual().fill( 0.0 );
                this->compute_residual( aWStar );

                Matrix< DDRMat > tResidual_Plus_Master
                = mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } );
                Matrix< DDRMat > tResidual_Plus_Slave
                = mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 1 ) }, { 0, 0 } );

                // perturbation of the coefficent
                tCoeffPert = tCoeff;
                tCoeffPert( iCoeffRow, iCoeffCol ) = tCoeffPert( iCoeffRow, iCoeffCol ) - aPerturbation * tCoeffPert( iCoeffRow, iCoeffCol );

                // setting the perturbed coefficients
                tFI->set_coeff( tCoeffPert );

                // reset properties, CM and SP for IWG
                this->reset_eval_flags();

                // evaluate the residual
                mSet->get_residual().fill( 0.0 );
                this->compute_residual( aWStar );

                Matrix< DDRMat > tResidual_Minus_Master
                = mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } );
                Matrix< DDRMat > tResidual_Minus_Slave
                = mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave  )( 0, 1 ) }, { 0, 0 } );

                // evaluate Jacobian
                aJacobiansFD( 0 )( tMasterNumDofType + iFI ).get_column( tDofCounter )
                = ( tResidual_Plus_Master - tResidual_Minus_Master )/ ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );
                aJacobiansFD( 1 )( tMasterNumDofType + iFI ).get_column( tDofCounter )
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

        bool IWG::check_jacobian( real                                             aPerturbation,
                                  real                                             aEpsilon,
                                  real                                             aWStar,
                                  moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                  moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFDs )
        {
            // compute jacobian with IWG
            this->compute_jacobian( aWStar );

            // compute jacobian by FD
            this->compute_jacobian_FD( aWStar, aPerturbation, aJacobiansFDs );

            // set size for comparison
            aJacobians.resize( aJacobiansFDs.size() );

            //define a boolean for check
            bool tCheckJacobian = true;

            // check each components
            for ( uint iJac = 0; iJac < aJacobiansFDs.size(); iJac++ )
            {
                // set size for comparison
                aJacobians( iJac ).resize( aJacobiansFDs( iJac ).size() );

                for( uint jJac = 0; jJac < aJacobiansFDs( iJac ).size(); jJac++ )
                {
                    uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
                    uint tIndexDep = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( jJac )( 0 ), mtk::Master_Slave::MASTER );

                    aJacobians( iJac )( jJac )
                    = mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                            { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } );

                    for( uint iiJac = 0; iiJac < aJacobiansFDs( iJac )( jJac ).n_rows(); iiJac++ )
                    {
                        for( uint jjJac = 0; jjJac < aJacobiansFDs( iJac )( jJac ).n_cols(); jjJac++ )
                        {
                            tCheckJacobian = tCheckJacobian && ( aJacobians( iJac )( jJac )( iiJac, jjJac ) - aJacobiansFDs( iJac )( jJac )( iiJac, jjJac ) < aEpsilon );
                        }
                    }
                }
            }

            // return bool
            return tCheckJacobian;
        }

//------------------------------------------------------------------------------
        bool IWG::check_jacobian_double( real                                             aPerturbation,
                                         real                                             aEpsilon,
                                         real                                             aWStar,
                                         moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                         moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFDs )
        {
            // compute jacobian with IWG
            this->compute_jacobian( aWStar );

            // compute jacobian by FD
            this->compute_jacobian_FD_double( aWStar, aPerturbation, aJacobiansFDs );

            // set jacobian size for comparison
            aJacobians.resize( 2 );

            //define a boolean for check
            bool tCheckJacobian = true;

            // check each components
            for ( uint iJac = 0; iJac < aJacobiansFDs.size(); iJac++ )
            {
                //set size for comparison
                aJacobians( iJac ).resize( aJacobiansFDs( iJac ).size() );

                uint tDofIndex;

                if( iJac == 0 )
                {
                    tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
                }
                else if( iJac == 1 )
                {
                    tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
                }
                else { MORIS_ERROR( false, "only case 0 and 1 implemented" ); }

                // init dof counter
                uint tCounter = 0;

                // loop over master dof type
                for( uint jJac = 0; jJac < mRequestedMasterGlobalDofTypes.size(); jJac++ )
                {
                    // get set index for master dof type
                    uint tIndexDep = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( jJac )( 0 ), mtk::Master_Slave::MASTER );

                    // fill jacobian matrix for comparison
                    aJacobians( iJac )( jJac )
                    = mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                            { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } );

                    // loop over the rows of jacobian
                    for( uint iiJac = 0; iiJac < aJacobiansFDs( iJac )( tCounter ).n_rows(); iiJac++ )
                    {
                        // loop over the columns of jacobian
                        for( uint jjJac = 0; jjJac < aJacobiansFDs( iJac )( tCounter ).n_cols(); jjJac++ )
                        {
                            // check component
                            tCheckJacobian = tCheckJacobian && ( aJacobians( iJac )( tCounter )( iiJac, jjJac ) - aJacobiansFDs( iJac )( tCounter )( iiJac, jjJac ) < aEpsilon );
                        }
                    }
                    // update dof counter
                    tCounter++;
                }

                // loop over slave dof types
                for( uint jJac = 0; jJac < mRequestedSlaveGlobalDofTypes.size(); jJac++ )
                {
                    // get set index for the slave dof type
                    uint tIndexDep = mSet->get_dof_index_for_type( mRequestedSlaveGlobalDofTypes( jJac )( 0 ), mtk::Master_Slave::SLAVE );

                    // fill jacobian matrix for comparison
                    aJacobians( iJac )( tCounter )
                    = mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                            { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } );

                    // loop over the rows of jacobian
                    for( uint iiJac = 0; iiJac < aJacobiansFDs( iJac )( tCounter ).n_rows(); iiJac++ )
                    {
                        // loop over the columns of jacobian
                        for( uint jjJac = 0; jjJac < aJacobiansFDs( iJac )( tCounter ).n_cols(); jjJac++ )
                        {
                            // check component
                            tCheckJacobian = tCheckJacobian && ( aJacobians( iJac )( tCounter )( iiJac, jjJac ) - aJacobiansFDs( iJac )( tCounter )( iiJac, jjJac ) < aEpsilon );
                        }
                    }
                    // update dof counter
                    tCounter++;
                }
            }

            // return bool
            return tCheckJacobian;
        }

//------------------------------------------------------------------------------
        void IWG::compute_drdpdv_FD( real aWStar,
                                     real aPerturbation,
                                     moris::Cell< Matrix< DDRMat > > & adrdpdvMatFD,
                                     moris::Cell< Matrix< DDRMat > > & adrdpdvGeoFD )
        {
            // MATERIAL PDV

            // get master number of dv types
            uint tNumDvType = mMasterGlobalDvTypes.size();

            // get the residaul dof type index in the set
            uint tDofIndex = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( 0 )( 0 ),
                                                           mtk::Master_Slave::MASTER );
            uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get number of rows
            uint tNumRows = tResDofAssemblyStop - tResDofAssemblyStart + 1;

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

            adrdpdvMatFD.resize( 1 );
            adrdpdvMatFD( 0 ).set_size( tNumRows, tNumCols, 0.0 );

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

                        // evaluate the residual
                        mSet->get_residual().fill( 0.0 );
                        this->compute_residual( aWStar );

                        Matrix< DDRMat > tResidual_Plus
                        =  mSet->get_residual()( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } );

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeffRow, iCoeffCol ) += - aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // evaluate the residual
                        mSet->get_residual().fill( 0.0 );
                        this->compute_residual( aWStar );

                        Matrix< DDRMat > tResidual_Minus
                        =  mSet->get_residual()( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } );

                        // evaluate Jacobian
                        uint tDvAssemblyStart = mSet->get_dv_assembly_map()( tDvIndex )( 0, 0 );

                        adrdpdvMatFD( 0 ).get_column( tDvAssemblyStart + tCounter )
                        = ( tResidual_Plus - tResidual_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

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
            adrdpdvGeoFD.resize( 1 );
            adrdpdvGeoFD( 0 ).set_size( tNumRows, tDerNumBases * tDerNumDimensions, 0.0 );

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

                    // evaluate the residual
                    mSet->get_residual().fill( 0.0 );
                    this->compute_residual( aWStar );

                    Matrix< DDRMat > tResidual_Plus
                    =  mSet->get_residual()( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } );

                    // perturbation of the coefficient
                    tCoeffPert = tCoeff;
                    tCoeffPert( iCoeffRow, iCoeffCol ) += - aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // setting the perturbed coefficients
                    tGI->set_space_coeff( tCoeffPert );

                    // reset properties, CM and SP for IWG
                    this->reset_eval_flags();

                    // evaluate the residual
                    mSet->get_residual().fill( 0.0 );
                    this->compute_residual( aWStar );

                    Matrix< DDRMat > tResidual_Minus
                    =  mSet->get_residual()( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } );

                    // evaluate drdpdvGeo
                    adrdpdvGeoFD( 0 ).get_column( tDvCounter )
                    = ( tResidual_Plus - tResidual_Minus ) / ( 2.0 * aPerturbation * tCoeff( iCoeffRow, iCoeffCol ) );

                    // update dv counter
                    tDvCounter++;
                }
            }
            // reset the coefficients values
            tGI->set_space_coeff( tCoeff );

            //print( adrdpdvMatFD( 0 ), "tdrdpdvMatFD" );
            //print( adrdpdvGeoFD( 0 ), "tdrdpdvGeoFD" );
        }

//------------------------------------------------------------------------------

}   // end fem namespace
}   // end moris namespace



/*
 * cl_FEM_Stabilization_Parameter.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: noel
 */

#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void Stabilization_Parameter::set_field_interpolator_manager( Field_Interpolator_Manager * aFieldInterpolatorManager,
                                                                      mtk::Master_Slave            aIsMaster )
        {
    //        // FIXME why does this not work?
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
                    MORIS_ERROR( false, "Stabilization_Parameter::set_field_interpolator_manager - can only be master or slave");
                    break;
                }
            }

            // FIXME
            // get the list of dof types for the SP
            moris::Cell< moris::Cell< MSI::Dof_Type > > tSPDofTypes
            = this->get_global_dof_type_list( aIsMaster );

            // get the number of dof type for the SP
            uint tNumDofTypes = tSPDofTypes.size();

            // set the size of the field interpolators list for the SP
            this->get_dof_field_interpolators( aIsMaster ).resize( tNumDofTypes, nullptr );

            // loop over the dof types
            for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
            {
                // grab the field interpolator for the dof type
                this->get_dof_field_interpolators( aIsMaster )( iDof )
                = this->get_field_interpolator_manager( aIsMaster )
                      ->get_field_interpolators_for_type( tSPDofTypes( iDof )( 0 ) );
            }
            // END FIXME

            // loop over the underlying constitutive models
            for( std::shared_ptr< Constitutive_Model > tCM : this->get_constitutive_models( aIsMaster ) )
            {
                if ( tCM != nullptr )
                {
                    // set the field interpolator manager for the property
                    tCM->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsMaster ) );
                }
            }

            // loop over the underlying properties
            for( std::shared_ptr< Property > tProp : this->get_properties( aIsMaster ) )
            {
                if ( tProp != nullptr )
                {
                    // set the field interpolator manager for the property
                    tProp->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsMaster ) );
                }
            }
        }

//------------------------------------------------------------------------------
        void Stabilization_Parameter::get_non_unique_dof_types
        ( moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // set the size of the dof type list for the set
            uint tCounter = 0;

            //loop over direct master dof dependencies
            for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
            {
                // update counter
                tCounter += mMasterDofTypes( iDOF ).size();
            }

            // get number of direct slave dof dependencies
            for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
            {
                //update counter
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
        }

//------------------------------------------------------------------------------
        void Stabilization_Parameter::get_non_unique_dof_and_dv_types
        ( moris::Cell< MSI::Dof_Type > & aDofTypes,
          moris::Cell< MSI::Dv_Type >  & aDvTypes )
        {
            // init dof and dv counters
            uint tDofCounter = 0;
            uint tDvCounter  = 0;

            //loop over direct master dof dependencies
            for ( uint iDof = 0; iDof < mMasterDofTypes.size(); iDof++ )
            {
                // update counter
                tDofCounter += mMasterDofTypes( iDof ).size();
            }

            //loop over direct master dv dependencies
            for ( uint iDv = 0; iDv < mMasterDvTypes.size(); iDv++ )
            {
                // update counter
                tDvCounter += mMasterDvTypes( iDv ).size();
            }

            // get number of direct slave dof dependencies
            for ( uint iDof = 0; iDof < mSlaveDofTypes.size(); iDof++ )
            {
                //update counter
                tDofCounter += mSlaveDofTypes( iDof ).size();
            }

            // get number of direct slave dv dependencies
            for ( uint iDv = 0; iDv < mSlaveDvTypes.size(); iDv++ )
            {
                //update counter
                tDvCounter += mSlaveDvTypes( iDv ).size();
            }

            // loop over the master properties
            for ( std::shared_ptr< Property > tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv types
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tProperty->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                                tActiveDvTypes );

                    //update counters
                    tDofCounter += tActiveDofTypes.size();
                    tDvCounter  += tActiveDvTypes.size();

                }
            }

            // loop over slave properties
            for ( std::shared_ptr< Property > tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv types
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tProperty->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                                tActiveDvTypes );

                    //update counters
                    tDofCounter += tActiveDofTypes.size();
                    tDvCounter  += tActiveDvTypes.size();
                }
            }

            // loop over master constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv types
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tCM->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                          tActiveDvTypes );

                    //update counters
                    tDofCounter += tActiveDofTypes.size();
                    tDvCounter  += tActiveDvTypes.size();
                }
            }

            // loop over slave constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    // get CM non unique dof and dv types
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tCM->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                          tActiveDvTypes );

                    //update counters
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
                    // get property non unique dof and dv types
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tProperty->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                                tActiveDvTypes );

                    // populate the dof list
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }

            // loop over slave properties
            for ( std::shared_ptr< Property > tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv types
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tProperty->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                                tActiveDvTypes );

                    // populate the dof list
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }

            // loop over the master constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv types
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tCM->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                          tActiveDvTypes );

                    // populate the dof list
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }

            // loop over the slave constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    // get CM non unique dof and dv types
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tCM->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                          tActiveDvTypes );

                    // populate the dof list
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }
        }


//------------------------------------------------------------------------------
        void Stabilization_Parameter::build_global_dof_type_list()
        {
            // MASTER-------------------------------------------------------
            // get the size of the dof type list
            uint tCounterMax = 0;

            // get number of dof types from penalty parameter
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
                if( tCM != nullptr )
                {
                    tCounterMax += tCM->get_global_dof_type_list().size();
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
                            mSlaveGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

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
                            mSlaveGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                            // update dof counter
                            tCounter++;
                        }
                    }
                }
            }

            // get the number of unique dof type groups for the penalty parameter
            mSlaveGlobalDofTypes.resize( tCounter );

            // number of global master and slave dof types
            uint tNumMasterGlobalDofTypes = mMasterGlobalDofTypes.size();
            uint tNumSlaveGlobalDofTypes  = mSlaveGlobalDofTypes.size();

            // set flag for evaluation
            mdPPdMasterDofEval.assign( tNumMasterGlobalDofTypes, true );
            mdPPdSlaveDofEval.assign( tNumSlaveGlobalDofTypes, true );

            // set storage for evaluation
            mdPPdMasterDof.resize( tNumMasterGlobalDofTypes );
            mdPPdSlaveDof.resize( tNumSlaveGlobalDofTypes );
        }

//------------------------------------------------------------------------------
        void Stabilization_Parameter::reset_eval_flags()
        {
            // reset the value flag
            mPPEval = true;

            // reset the master dof derivative flags
            uint tNumMasterDofTypes = mMasterGlobalDofTypes.size();
            mdPPdMasterDofEval.assign( tNumMasterDofTypes, true );

            // reset the slave dof derivative flags
            uint tNumSlaveDofTypes = mSlaveGlobalDofTypes.size();
            mdPPdSlaveDofEval.assign( tNumSlaveDofTypes, true );

            // reset the master dv derivative flags
            uint tNumMasterDvTypes = mMasterGlobalDvTypes.size();
            mdPPdMasterDvEval.assign( tNumMasterDvTypes, true );

            // reset the slave dv derivative flags
            uint tNumSlaveDvTypes = mSlaveGlobalDvTypes.size();
            mdPPdSlaveDvEval.assign( tNumSlaveDvTypes, true );

            // reset underlying master constitutive models
            for( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
            {
                if( tCM != nullptr )
                {
                    tCM->reset_eval_flags();
                }
            }

            // reset underlying slave constitutive models
            for( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    tCM->reset_eval_flags();
                }
            }

            // reset underlying master properties
            for( std::shared_ptr< Property > tProp : mMasterProp )
            {
                if( tProp != nullptr )
                {
                    tProp->reset_eval_flags();
                }
            }

            // reset underlying slave properties
            for( std::shared_ptr< Property > tProp : mSlaveProp )
            {
                if( tProp != nullptr )
                {
                    tProp->reset_eval_flags();
                }
            }
        }

//------------------------------------------------------------------------------
    }/*end_fem_namespace */
}/*end_moris_namespace */



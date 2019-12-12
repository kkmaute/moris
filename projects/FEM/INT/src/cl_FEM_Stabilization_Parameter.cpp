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
        void Stabilization_Parameter::get_non_unique_global_dof_type_list
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



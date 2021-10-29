/*
 * cl_FEM_Stabilization_Parameter.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: noel
 */

#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Cluster_Measure.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::print_names()
        {
            std::cout<<"----------"<<std::endl;
            std::cout<<"SP: "<<mName<<std::endl;

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
                    std::cout<<"Master CM: "<<mMasterCM( iCM )->get_name()<<std::endl;
                }
            }
            for( uint iCM = 0; iCM < mSlaveCM.size(); iCM++ )
            {
                if( mSlaveCM( iCM ) != nullptr )
                {
                    std::cout<<"Slave CM:  "<<mSlaveCM( iCM )->get_name()<<std::endl;
                }
            }
            std::cout<<"----------"<<std::endl;
        }

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::set_cluster( fem::Cluster * aCluster )
        {
            // set a cluster
            mCluster = aCluster;
        }

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::reset_eval_flags()
        {
            // reset the value flag
            mPPEval = true;

            // reset the master dof derivative flags
            mdPPdMasterDofEval.fill( true );

            // reset the slave dof derivative flags
            mdPPdSlaveDofEval.fill( true );

            // reset the master dv derivative flags
            mdPPdMasterDvEval.fill( true );

            // reset the slave dv derivative flags
            mdPPdSlaveDvEval.fill( true );

            // reset underlying master constitutive models
            for( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
            {
                if( tCM != nullptr )
                {
                    tCM->reset_eval_flags();
                }
            }

            // reset underlying slave constitutive models
            for( const std::shared_ptr< Constitutive_Model > & tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    tCM->reset_eval_flags();
                }
            }

            // reset underlying master properties
            for( const std::shared_ptr< Property > & tProp : mMasterProp )
            {
                if( tProp != nullptr )
                {
                    tProp->reset_eval_flags();
                }
            }

            // reset underlying slave properties
            for( const std::shared_ptr< Property > & tProp : mSlaveProp )
            {
                if( tProp != nullptr )
                {
                    tProp->reset_eval_flags();
                }
            }
        }

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "Stabilization_Parameter::set_property - SP %s - Unknown aPropertyString : %s \n",
                    mName.c_str(),
                    aPropertyString.c_str() );

            // set the property in the property cell
            this->get_properties( aIsMaster )( mPropertyMap[ aPropertyString ] ) = aProperty;
        }

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster )
            {
                // check that aConstitutiveString makes sense
                MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(),
                        "Stabilization_Parameter::set_constitutive_model - SP %s - Unknown aConstitutiveString : %s \n",
                        mName.c_str(),
                        aConstitutiveString.c_str() );

                // set the constitutive model in the CM cell
                this->get_constitutive_models( aIsMaster )( mConstitutiveMap[ aConstitutiveString ] ) = aConstitutiveModel;
            }

        //------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< Property > > & Stabilization_Parameter::get_properties(
                mtk::Master_Slave aIsMaster )
        {
            // switch on master/slave
            switch( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER :
                {
                    // return master property pointers
                    return mMasterProp;
                }

                // if slave
                case mtk::Master_Slave::SLAVE :
                {
                    // return slave property pointers
                    return mSlaveProp;
                }

                // if none
                default:
                {
                    MORIS_ASSERT( false, "Stabilization_Parameter::get_properties - can only be master or slave." );
                    return mMasterProp;
                }
            }
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< Constitutive_Model > > & Stabilization_Parameter::get_constitutive_models(
                mtk::Master_Slave aIsMaster )
        {
            // switch on master/slave
            switch( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER :
                {
                    // return master CM pointers
                    return mMasterCM;
                }

                // if slave
                case mtk::Master_Slave::SLAVE :
                {
                    // return slave CM pointers
                    return mSlaveCM;
                }

                // if none
                default:
                {
                    MORIS_ASSERT( false, "Stabilization_Parameter::get_constitutive_models - can only be master or slave." );
                    return mMasterCM;
                }
            }
        }

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::set_dof_type_list(
                const moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                mtk::Master_Slave                                   aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER  :
                {
                    mMasterDofTypes = aDofTypes;
                    break;
                }
                case mtk::Master_Slave::SLAVE :
                {
                    mSlaveDofTypes = aDofTypes;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "Stabilization_Parameter::set_dof_type_list - can only be MASTER or SLAVE.");
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< MSI::Dof_Type > > & Stabilization_Parameter::get_dof_type_list(
                mtk::Master_Slave aIsMaster ) const
        {
            // switch on master/slave
            switch( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER :
                {
                    // return master global dof type list
                    return mMasterDofTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE :
                {
                    // return slave global dof type list
                    return mSlaveDofTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "Stabilization_Parameter::get_dof_type_list - can only be master or slave." );
                    return mMasterDofTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::set_dv_type_list(
                moris::Cell< moris::Cell< PDV_Type > > & aDvTypes,
                mtk::Master_Slave                        aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    mMasterDvTypes = aDvTypes;
                    break;
                }
                case mtk::Master_Slave::SLAVE :
                {
                    mSlaveDvTypes = aDvTypes;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "Stabilization_Parameter::set_dv_type_list - can only be MASTER or SLAVE.");
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< PDV_Type > > & Stabilization_Parameter::get_dv_type_list(
                mtk::Master_Slave aIsMaster ) const
        {
            // switch on master/slave
            switch( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER :
                {
                    // return master global dof type list
                    return mMasterDvTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE :
                {
                    // return slave global dof type list
                    return mSlaveDvTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "Stabilization_Parameter::get_dv_type_list - can only be master or slave." );
                    return mMasterDvTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::set_field_interpolator_manager(
                Field_Interpolator_Manager * aFieldInterpolatorManager,
                mtk::Master_Slave            aIsMaster )
        {
            // switch on master or slave
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    mMasterFIManager = aFieldInterpolatorManager;
                    break;
                }

                case mtk::Master_Slave::SLAVE :
                {
                    mSlaveFIManager = aFieldInterpolatorManager;
                    break;
                }

                default :
                {
                    MORIS_ERROR( false, "Stabilization_Parameter::set_field_interpolator_manager - can only be master or slave");
                }
            }

            // loop over the underlying constitutive models
            for( const std::shared_ptr< Constitutive_Model > & tCM : this->get_constitutive_models( aIsMaster ) )
            {
                if ( tCM != nullptr )
                {
                    // set the field interpolator manager for the property
                    tCM->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsMaster ) );
                }
            }

            // loop over the underlying properties
            for( const std::shared_ptr< Property > & tProp : this->get_properties( aIsMaster ) )
            {
                if ( tProp != nullptr )
                {
                    // set the field interpolator manager for the property
                    tProp->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsMaster ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        Field_Interpolator_Manager * Stabilization_Parameter::get_field_interpolator_manager(
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    return mMasterFIManager;
                }

                case mtk::Master_Slave::SLAVE :
                {
                    return mSlaveFIManager;
                }

                default :
                {
                    MORIS_ERROR( false, "Stabilization_Parameter::get_field_interpolator_manager - can only be master or slave");
                    return mMasterFIManager;
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< MSI::Dof_Type > > & Stabilization_Parameter::get_global_dof_type_list(
                mtk::Master_Slave aIsMaster )
        {
            if( mGlobalDofBuild )
            {
                // build the stabilization parameter global dof type list
                this->build_global_dof_type_list();

                // build the stabilization parameter global dof type map
                this->build_global_dof_type_map();

                // update build flag
                mGlobalDofBuild = false;
            }

            // switch on master/slave
            switch( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER :
                {
                    // return master global dof type list
                    return mMasterGlobalDofTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE :
                {
                    // return slave global dof type list
                    return mSlaveGlobalDofTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "Stabilization_Parameter::get_global_dof_type_list - can only be master or slave." );
                    return mMasterGlobalDofTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::get_non_unique_dof_types(
                moris::Cell< MSI::Dof_Type > & aDofTypes )
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
            for ( const std::shared_ptr< Property > & tProperty : mMasterProp )
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
            for ( const std::shared_ptr< Property > & tProperty : mSlaveProp )
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
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
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
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mSlaveCM )
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
            for ( const std::shared_ptr< Property > & tProperty : mMasterProp )
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
            for ( const std::shared_ptr< Property > & tProperty : mSlaveProp )
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
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
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
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mSlaveCM )
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

        void Stabilization_Parameter::get_non_unique_dof_and_dv_types(
                moris::Cell< MSI::Dof_Type > & aDofTypes,
                moris::Cell< PDV_Type >      & aDvTypes )
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
            for ( const std::shared_ptr< Property > tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv types
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tProperty->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes);

                    //update counters
                    tDofCounter += tActiveDofTypes.size();
                    tDvCounter  += tActiveDvTypes.size();

                }
            }

            // loop over slave properties
            for ( const std::shared_ptr< Property > & tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv types
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< PDV_Type >      tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tProperty->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes);

                    //update counters
                    tDofCounter += tActiveDofTypes.size();
                    tDvCounter  += tActiveDvTypes.size();
                }
            }

            // loop over master constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv types
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tCM->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    //update counters
                    tDofCounter += tActiveDofTypes.size();
                    tDvCounter  += tActiveDvTypes.size();
                }
            }

            // loop over slave constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    // get CM non unique dof and dv types
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tCM->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

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
            for ( const std::shared_ptr< Property > & tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv types
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< PDV_Type >      tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tProperty->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // populate the dof list
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }

            // loop over slave properties
            for ( const std::shared_ptr< Property > & tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv types
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tProperty->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes);

                    // populate the dof list
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }

            // loop over the master constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv types
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tCM->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // populate the dof list
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }

            // loop over the slave constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    // get CM non unique dof and dv types
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tCM->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

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
            for ( const std::shared_ptr< Property > & tProperty : mMasterProp )
            {
                if( tProperty != nullptr )
                {
                    tCounterMax += tProperty->get_dof_type_list().size();
                }
            }

            // get number of dof types from constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
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
            for ( const std::shared_ptr< Property > & tProperty : mMasterProp )
            {
                if( tProperty != nullptr )
                {
                    // get dof types for property
                    const moris::Cell< moris::Cell< MSI::Dof_Type > > & tActiveDofType =
                            tProperty->get_dof_type_list();

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
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
            {
                if( tCM != nullptr )
                {
                    // get dof types for constitutive model
                    const moris::Cell< moris::Cell< MSI::Dof_Type > > & tActiveDofType =
                            tCM->get_global_dof_type_list();

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
            for ( const std::shared_ptr< Property > & tProperty : mSlaveProp )
            {
                if( tProperty != nullptr )
                {
                    tCounterMax += tProperty->get_dof_type_list().size();
                }
            }

            // get number of dof types from constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mSlaveCM )
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
            for ( const std::shared_ptr< Property > & tProperty : mSlaveProp )
            {
                if( tProperty != nullptr )
                {
                    // get dof types for property
                    const moris::Cell< moris::Cell< MSI::Dof_Type > > & tActiveDofType =
                            tProperty->get_dof_type_list();

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
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    // get dof types for constitutive model
                    const moris::Cell< moris::Cell< MSI::Dof_Type > > & tActiveDofType =
                            tCM->get_global_dof_type_list();

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
            mdPPdMasterDofEval.set_size( tNumMasterGlobalDofTypes, 1, true );
            mdPPdSlaveDofEval.set_size( tNumSlaveGlobalDofTypes, 1, true );

            // set storage for evaluation
            mdPPdMasterDof.resize( tNumMasterGlobalDofTypes );
            mdPPdSlaveDof.resize( tNumSlaveGlobalDofTypes );
        }

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::build_global_dof_type_map()
        {
            // MASTER-------------------------------------------------------
            // get number of global dof types
            uint tNumDofTypes = mMasterGlobalDofTypes.size();

            // determine the max Dof_Type enum
            sint tMaxEnum = 0;
            for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
            {
                tMaxEnum = std::max( tMaxEnum, static_cast< int >( mMasterGlobalDofTypes( iDOF )( 0 ) ) );
            }
            tMaxEnum++;

            // set the Dof_Type map size
            mMasterGlobalDofTypeMap.set_size( tMaxEnum, 1, -1 );

            // fill the Dof_Type map
            for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
            {
                // fill the property map
                mMasterGlobalDofTypeMap( static_cast< int >( mMasterGlobalDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
            }

            // SLAVE-------------------------------------------------------
            // get number of global dof types
            tNumDofTypes = mSlaveGlobalDofTypes.size();

            // determine the max Dof_Type enum
            tMaxEnum = 0;
            for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
            {
                tMaxEnum = std::max( tMaxEnum, static_cast< int >( mSlaveGlobalDofTypes( iDOF )( 0 ) ) );
            }
            tMaxEnum++;

            // set the dof type map size
            mSlaveGlobalDofTypeMap.set_size( tMaxEnum, 1, -1 );

            // fill the dof type map
            for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
            {
                // fill the property map
                mSlaveGlobalDofTypeMap( static_cast< int >( mSlaveGlobalDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
            }
        }

        //------------------------------------------------------------------------------

        const Matrix< DDSMat > & Stabilization_Parameter::get_global_dof_type_map(
                mtk::Master_Slave aIsMaster )
        {
            // switch on master/slave
            switch( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER :
                {
                    // return master global dof type map
                    return mMasterGlobalDofTypeMap;
                }
                // if slave
                case mtk::Master_Slave::SLAVE :
                {
                    // return slave global dof type map
                    return mSlaveGlobalDofTypeMap;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "Stabilization_Parameter::get_global_dof_type_map - can only be master or slave." );
                    return mMasterGlobalDofTypeMap;
                }
            }
        }

        //------------------------------------------------------------------------------

        bool Stabilization_Parameter::check_dof_dependency(
                const moris::Cell< MSI::Dof_Type > & aDofType,
                mtk::Master_Slave                    aIsMaster )
        {
            // set bool for dependency
            bool tDofDependency = false;

            // get dof type index
            uint tDofIndex = static_cast< uint >( aDofType( 0 ) );

            // if aDofType is an active dof type for the stabilization parameter
            if( tDofIndex < this->get_global_dof_type_map( aIsMaster ).numel()
                    && this->get_global_dof_type_map( aIsMaster )( tDofIndex ) != -1 )
            {
                // bool is set to true
                tDofDependency = true;
            }
            // return bool for dependency
            return tDofDependency;
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< PDV_Type > > & Stabilization_Parameter::get_global_dv_type_list(
                mtk::Master_Slave aIsMaster )
        {
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
                case mtk::Master_Slave::MASTER :
                {
                    // return master global dv type list
                    return mMasterGlobalDvTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE :
                {
                    // return slave global dv type list
                    return mSlaveGlobalDvTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "Stabilization_Parameter::get_global_dv_type_list - can only be master or slave." );
                    return mMasterGlobalDvTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::build_global_dv_type_list()
        {
            // MASTER-------------------------------------------------------
            // get the size of the dv type list
            uint tCounterMax = 0;

            // get number of dv types from penalty parameter
            tCounterMax += mMasterDvTypes.size();

            // get number of dv types from properties
            for ( const std::shared_ptr< Property > & tProperty : mMasterProp )
            {
                if( tProperty != nullptr )
                {
                    tCounterMax += tProperty->get_dv_type_list().size();
                }
            }

            // get number of dof types from constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
            {
                if( tCM != nullptr )
                {
                    tCounterMax += tCM->get_global_dv_type_list().size();
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
            for ( const std::shared_ptr< Property > & tProperty : mMasterProp )
            {
                if( tProperty != nullptr )
                {
                    // get dv types for property
                    const moris::Cell< moris::Cell< PDV_Type > > & tActiveDvType =
                            tProperty->get_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                    {
                        // check enum is not already in the list
                        bool tCheck = false;
                        for( uint i = 0; i < tCounter; i++ )
                        {
                            tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                        }

                        // if dof enum not in the list
                        if ( !tCheck )
                        {
                            // put the dv type in the checklist
                            tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                            // put the dv type in the global type list
                            mMasterGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                            // update dof counter
                            tCounter++;
                        }
                    }
                }
            }

            // get dof type from constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
            {
                if( tCM != nullptr )
                {
                    // get dof types for constitutive model
                   const moris::Cell< moris::Cell< PDV_Type > > & tActiveDvType =
                           tCM->get_global_dv_type_list();

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
            for ( const std::shared_ptr< Property > & tProperty : mSlaveProp )
            {
                if( tProperty != nullptr )
                {
                    tCounterMax += tProperty->get_dv_type_list().size();
                }
            }

            // get number of dv types from constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    tCounterMax += tCM->get_global_dv_type_list().size();
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
            for ( const std::shared_ptr< Property > & tProperty : mSlaveProp )
            {
                if( tProperty != nullptr )
                {
                    // get dv types for property
                    const moris::Cell< moris::Cell< PDV_Type > > & tActiveDvType =
                            tProperty->get_dv_type_list();

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
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
            {
                if( tCM != nullptr )
                {
                    // get dv types for constitutive model
                    const moris::Cell< moris::Cell< PDV_Type > > & tActiveDvType =
                            tCM->get_global_dv_type_list();

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

            // build global dv type map
            this->build_global_dv_type_map();

            // number of global master and slave dv types
            uint tNumMasterGlobalDvTypes = mMasterGlobalDvTypes.size();
            uint tNumSlaveGlobalDvTypes  = mSlaveGlobalDvTypes.size();

            // set flag for evaluation
            mdPPdMasterDvEval.set_size( tNumMasterGlobalDvTypes, 1, true );
            mdPPdSlaveDvEval.set_size( tNumSlaveGlobalDvTypes, 1, true );

            // set storage for evaluation
            mdPPdMasterDv.resize( tNumMasterGlobalDvTypes );
            mdPPdSlaveDv.resize( tNumSlaveGlobalDvTypes );
        }

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::build_global_dv_type_map()
        {
            // MASTER-------------------------------------------------------
            // get number of global dof types
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
            // get number of global dv types
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

        bool Stabilization_Parameter::check_master_dv_dependency(
                const moris::Cell< PDV_Type > & aDvType )
        {
            // set bool for dependency
            bool tDvDependency = false;

            // get dv type index
            uint tDvIndex = static_cast< uint >( aDvType( 0 ) );

            // if aDvType is an active dv type for the constitutive model
            if( tDvIndex < mMasterGlobalDvTypeMap.numel() && mMasterGlobalDvTypeMap( tDvIndex ) != -1 )
            {
                // bool is set to true
                tDvDependency = true;
            }
            // return bool for dependency
            return tDvDependency;
        }

        //------------------------------------------------------------------------------

        bool Stabilization_Parameter::check_slave_dv_dependency(
                const moris::Cell< PDV_Type > & aDvType )
        {
            // set bool for dependency
            bool tDvDependency = false;

            // get dv type index
            uint tDvIndex = static_cast< uint >( aDvType( 0 ) );

            // if aDvType is an active dv type for the constitutive model
            if( tDvIndex < mSlaveGlobalDvTypeMap.numel() && mSlaveGlobalDvTypeMap( tDvIndex ) != -1 )
            {
                // bool is set to true
                tDvDependency = true;
            }
            // return bool for dependency
            return tDvDependency;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Stabilization_Parameter::val()
        {
            // if the penalty parameter was not evaluated
            if( mPPEval )
            {
                // evaluate the penalty parameter
                this->eval_SP();

                // set bool for evaluation
                mPPEval = false;
            }
            // return the penalty parameter value
            return mPPVal;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Stabilization_Parameter::dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofType )
        {
            // if aDofType is not an active dof type for the property
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType, mtk::Master_Slave::MASTER ),
                    "Stabilization_Parameter::dSPdMasterDOF - no dependency on this dof type." );

            // get the dof index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdPPdMasterDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dSPdMasterDOF( aDofType );

                // set bool for evaluation
                mdPPdMasterDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdPPdMasterDof( tDofIndex );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Stabilization_Parameter::dSPdSlaveDOF(
                const moris::Cell< MSI::Dof_Type > & aDofType )
        {
            // if aDofType is not an active dof type for the property
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType, mtk::Master_Slave::SLAVE ),
                    "Stabilization_Parameter::dSPdSlaveDOF - no dependency on this dof type." );

            // get the dof index
            uint tDofIndex = mSlaveGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdPPdSlaveDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dSPdSlaveDOF( aDofType );

                // set bool for evaluation
                mdPPdSlaveDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdPPdSlaveDof( tDofIndex );
        }

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::eval_dSPdMasterDOF_FD(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adSPdDOF_FD,
                real                                 aPerturbation,
                fem::FDScheme_Type                   aFDSchemeType )
        {
            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get the dof index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the field interpolator for type
            Field_Interpolator* tFI =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get number of master dofs wrt which derivative is computed
            uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
            uint tDerNumBases  = tFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFI->get_number_of_fields();

            // get unperturbed SP value
            Matrix< DDRMat > tUnperturbedSPValue = this->val();

            // set size for derivative
            adSPdDOF_FD.set_size( tUnperturbedSPValue.n_rows(), tDerNumDof, 0.0 );

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFI->get_coeff();

            // initialize dof counter
            uint tDofCounter = 0;

            // loop over coefficients columns
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over coefficients rows
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // check that perturbation is not zero
                    if( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed traction contribution to dtractiondu
                        adSPdDOF_FD.get_column( tDofCounter ) +=
                                tFDScheme( 1 )( 0 ) * tUnperturbedSPValue /
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // skip first point in FD
                        tStartPoint = 1;
                    }

                    // loop over the points for FD
                    for( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                    {
                        // reset the perturbed coefficients
                        Matrix< DDRMat > tCoeffPert = tCoeff;

                        // perturb the coefficient
                        tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                        // set the perturbed coefficients to FI
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        this->reset_eval_flags();

                        // assemble derivatiev of SP wrt master dof type
                        adSPdDOF_FD.get_column( tDofCounter ) +=
                                tFDScheme( 1 )( iPoint ) * this->val() /
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );

            // set value to storage
            mdPPdMasterDof( tDofIndex ) = adSPdDOF_FD;
        }

        //------------------------------------------------------------------------------

        void Stabilization_Parameter::eval_dSPdSlaveDOF_FD(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adSPdDOF_FD,
                real                                 aPerturbation,
                fem::FDScheme_Type                   aFDSchemeType )
        {
            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get the dof index
            uint tDofIndex = mSlaveGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the field interpolator for type
            Field_Interpolator* tFI =
                    mSlaveFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get number of master dofs wrt which derivative is computed
            uint tDerNumDof    = tFI->get_number_of_space_time_coefficients();
            uint tDerNumBases  = tFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFI->get_number_of_fields();

            // get unperturbed SP value
            Matrix< DDRMat > tUnperturbedSPValue = this->val();

            // set size for derivative
            adSPdDOF_FD.set_size( tUnperturbedSPValue.n_rows(), tDerNumDof, 0.0 );

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFI->get_coeff();

            // initialize dof counter
            uint tDofCounter = 0;

            // loop over coefficients columns
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over coefficients rows
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // check that perturbation is not zero
                    if( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed traction contribution to dtractiondu
                        adSPdDOF_FD.get_column( tDofCounter ) +=
                                tFDScheme( 1 )( 0 ) * tUnperturbedSPValue /
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // skip first point in FD
                        tStartPoint = 1;
                    }

                    // loop over the points for FD
                    for( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                    {
                        // reset the perturbed coefficients
                        Matrix< DDRMat > tCoeffPert = tCoeff;

                        // perturb the coefficient
                        tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                        // set the perturbed coefficients to FI
                        tFI->set_coeff( tCoeffPert );

                        // reset properties
                        this->reset_eval_flags();

                        // assemble the jacobian
                        mdPPdSlaveDof( tDofIndex ).get_column( tDofCounter ) +=
                                tFDScheme( 1 )( iPoint ) * this->val() /
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );

            // set value to storage
            mdPPdSlaveDof( tDofIndex ) = adSPdDOF_FD;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Stabilization_Parameter::dSPdMasterDV(
                const moris::Cell< PDV_Type > & aDvTypes )
        {
            // if aDofType is not an active dv type for the property
            MORIS_ERROR(
                    this->check_master_dv_dependency( aDvTypes ),
                    "Penalty_Parameter::dSPdMasterDV - no dependency on this dv type." );

            // get the dv index
            uint tDvIndex = mMasterGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdPPdMasterDofEval( tDvIndex ) )
            {
                // evaluate the derivative
                this->eval_dSPdMasterDV( aDvTypes );

                // set bool for evaluation
                mdPPdMasterDofEval( tDvIndex ) = false;
            }

            // return the derivative
            return mdPPdMasterDof( tDvIndex );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Stabilization_Parameter::dSPdSlaveDV(
                const moris::Cell< PDV_Type > & aDvTypes )
        {
            // if aDofType is not an active dv type for the property
            MORIS_ERROR(
                    this->check_slave_dv_dependency( aDvTypes ),
                    "Stabilization_Parameter::dSPdSlaveDV - no dependency on this dv type." );

            // get the dv index
            uint tDvIndex = mSlaveGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdPPdSlaveDvEval( tDvIndex ) )
            {
                // evaluate the derivative
                this->eval_dSPdSlaveDV( aDvTypes );

                // set bool for evaluation
                mdPPdSlaveDvEval( tDvIndex ) = false;
            }

            // return the derivative
            return mdPPdSlaveDv( tDvIndex );
        }

        //------------------------------------------------------------------------------
    }/*end_fem_namespace */
}/*end_moris_namespace */



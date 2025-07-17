/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Stabilization_Parameter.cpp
 *
 */

#include "cl_FEM_Stabilization_Parameter.hpp"

#include <utility>
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Cluster_Measure.hpp"
#include "fn_equal_to.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::print_names()
    {
        std::cout << "----------" << '\n';
        std::cout << "SP: " << mName << '\n';

        // properties
        for ( uint iProp = 0; iProp < mLeaderProp.size(); iProp++ )
        {
            if ( mLeaderProp( iProp ) != nullptr )
            {
                std::cout << "Leader property: " << mLeaderProp( iProp )->get_name() << '\n';
            }
        }
        for ( uint iProp = 0; iProp < mFollowerProp.size(); iProp++ )
        {
            if ( mFollowerProp( iProp ) != nullptr )
            {
                std::cout << "Follower property:  " << mFollowerProp( iProp )->get_name() << '\n';
            }
        }

        // CM
        for ( uint iCM = 0; iCM < mLeaderCM.size(); iCM++ )
        {
            if ( mLeaderCM( iCM ) != nullptr )
            {
                std::cout << "Leader CM: " << mLeaderCM( iCM )->get_name() << '\n';
            }
        }
        for ( uint iCM = 0; iCM < mFollowerCM.size(); iCM++ )
        {
            if ( mFollowerCM( iCM ) != nullptr )
            {
                std::cout << "Follower CM:  " << mFollowerCM( iCM )->get_name() << '\n';
            }
        }
        std::cout << "----------" << '\n';
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
    {
        // set a cluster
        mParameters = aParameters;
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::set_cluster( fem::Cluster* aCluster )
    {
        // set a cluster
        mCluster = aCluster;
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::reset_eval_flags()
    {
        // reset the value flag
        mPPEval = true;

        // reset the leader dof derivative flags
        mdPPdLeaderDofEval.fill( true );

        // reset the follower dof derivative flags
        mdPPdFollowerDofEval.fill( true );

        // reset the leader dv derivative flags
        mdPPdLeaderDvEval.fill( true );

        // reset the follower dv derivative flags
        mdPPdFollowerDvEval.fill( true );

        // reset underlying leader constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
        {
            if ( tCM != nullptr )
            {
                tCM->reset_eval_flags();
            }
        }

        // reset underlying follower constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mFollowerCM )
        {
            if ( tCM != nullptr )
            {
                tCM->reset_eval_flags();
            }
        }

        // reset underlying leader properties
        for ( const std::shared_ptr< Property >& tProp : mLeaderProp )
        {
            if ( tProp != nullptr )
            {
                tProp->reset_eval_flags();
            }
        }

        // reset underlying follower properties
        for ( const std::shared_ptr< Property >& tProp : mFollowerProp )
        {
            if ( tProp != nullptr )
            {
                tProp->reset_eval_flags();
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::set_property(
            std::shared_ptr< Property > aProperty,
            const std::string&          aPropertyString,
            mtk::Leader_Follower        aIsLeader )
    {
        // check that aPropertyString makes sense
        MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                "Stabilization_Parameter::set_property - SP %s - Unknown aPropertyString : %s \n",
                mName.c_str(),
                aPropertyString.c_str() );

        // set the property in the property cell
        this->get_properties( aIsLeader )( mPropertyMap[ aPropertyString ] ) = std::move( aProperty );
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::set_constitutive_model(
            std::shared_ptr< Constitutive_Model > aConstitutiveModel,
            const std::string&                    aConstitutiveString,
            mtk::Leader_Follower                  aIsLeader )
    {
        // check that aConstitutiveString makes sense
        MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(),
                "Stabilization_Parameter::set_constitutive_model - SP %s - Unknown aConstitutiveString : %s \n",
                mName.c_str(),
                aConstitutiveString.c_str() );

        // set the constitutive model in the CM cell
        this->get_constitutive_models( aIsLeader )( mConstitutiveMap[ aConstitutiveString ] ) = std::move( aConstitutiveModel );
    }

    //------------------------------------------------------------------------------

    Vector< std::shared_ptr< Property > >&
    Stabilization_Parameter::get_properties(
            mtk::Leader_Follower aIsLeader )
    {
        // switch on leader/follower
        switch ( aIsLeader )
        {
            // if leader
            case mtk::Leader_Follower::LEADER:
            {
                // return leader property pointers
                return mLeaderProp;
            }

            // if follower
            case mtk::Leader_Follower::FOLLOWER:
            {
                // return follower property pointers
                return mFollowerProp;
            }

            // if none
            default:
            {
                MORIS_ASSERT( false, "Stabilization_Parameter::get_properties - can only be leader or follower." );
                return mLeaderProp;
            }
        }
    }

    //------------------------------------------------------------------------------

    Vector< std::shared_ptr< Constitutive_Model > >&
    Stabilization_Parameter::get_constitutive_models(
            mtk::Leader_Follower aIsLeader )
    {
        // switch on leader/follower
        switch ( aIsLeader )
        {
            // if leader
            case mtk::Leader_Follower::LEADER:
            {
                // return leader CM pointers
                return mLeaderCM;
            }

            // if follower
            case mtk::Leader_Follower::FOLLOWER:
            {
                // return follower CM pointers
                return mFollowerCM;
            }

            // if none
            default:
            {
                MORIS_ASSERT( false, "Stabilization_Parameter::get_constitutive_models - can only be leader or follower." );
                return mLeaderCM;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::set_dof_type_list(
            const Vector< Vector< MSI::Dof_Type > >& aDofTypes,
            mtk::Leader_Follower                     aIsLeader )
    {
        switch ( aIsLeader )
        {
            case mtk::Leader_Follower::LEADER:
            {
                mLeaderDofTypes = aDofTypes;
                break;
            }
            case mtk::Leader_Follower::FOLLOWER:
            {
                mFollowerDofTypes = aDofTypes;
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Stabilization_Parameter::set_dof_type_list - can only be LEADER or FOLLOWER." );
            }
        }
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< MSI::Dof_Type > >&
    Stabilization_Parameter::get_dof_type_list(
            mtk::Leader_Follower aIsLeader ) const
    {
        // switch on leader/follower
        switch ( aIsLeader )
        {
            // if leader
            case mtk::Leader_Follower::LEADER:
            {
                // return leader global dof type list
                return mLeaderDofTypes;
            }
            // if follower
            case mtk::Leader_Follower::FOLLOWER:
            {
                // return follower global dof type list
                return mFollowerDofTypes;
            }
            // if none
            default:
            {
                MORIS_ASSERT( false, "Stabilization_Parameter::get_dof_type_list - can only be leader or follower." );
                return mLeaderDofTypes;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::set_dv_type_list(
            Vector< Vector< gen::PDV_Type > >& aDvTypes,
            mtk::Leader_Follower               aIsLeader )
    {
        switch ( aIsLeader )
        {
            case mtk::Leader_Follower::LEADER:
            {
                mLeaderDvTypes = aDvTypes;
                break;
            }
            case mtk::Leader_Follower::FOLLOWER:
            {
                mFollowerDvTypes = aDvTypes;
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Stabilization_Parameter::set_dv_type_list - can only be LEADER or FOLLOWER." );
            }
        }
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< gen::PDV_Type > >&
    Stabilization_Parameter::get_dv_type_list(
            mtk::Leader_Follower aIsLeader ) const
    {
        // switch on leader/follower
        switch ( aIsLeader )
        {
            // if leader
            case mtk::Leader_Follower::LEADER:
            {
                // return leader global dof type list
                return mLeaderDvTypes;
            }
            // if follower
            case mtk::Leader_Follower::FOLLOWER:
            {
                // return follower global dof type list
                return mFollowerDvTypes;
            }
            // if none
            default:
            {
                MORIS_ASSERT( false, "Stabilization_Parameter::get_dv_type_list - can only be leader or follower." );
                return mLeaderDvTypes;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::set_field_interpolator_manager(
            Field_Interpolator_Manager* aFieldInterpolatorManager,
            mtk::Leader_Follower        aIsLeader )
    {
        // switch on leader or follower
        switch ( aIsLeader )
        {
            case mtk::Leader_Follower::LEADER:
            {
                mLeaderFIManager = aFieldInterpolatorManager;
                break;
            }

            case mtk::Leader_Follower::FOLLOWER:
            {
                mFollowerFIManager = aFieldInterpolatorManager;
                break;
            }

            default:
            {
                MORIS_ERROR( false, "Stabilization_Parameter::set_field_interpolator_manager - can only be leader or follower" );
            }
        }

        // loop over the underlying constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : this->get_constitutive_models( aIsLeader ) )
        {
            if ( tCM != nullptr )
            {
                // set the field interpolator manager for the property
                tCM->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsLeader ) );
            }
        }

        // loop over the underlying properties
        for ( const std::shared_ptr< Property >& tProp : this->get_properties( aIsLeader ) )
        {
            if ( tProp != nullptr )
            {
                // set the field interpolator manager for the property
                tProp->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsLeader ) );
            }
        }
    }

    //------------------------------------------------------------------------------

    Field_Interpolator_Manager*
    Stabilization_Parameter::get_field_interpolator_manager(
            mtk::Leader_Follower aIsLeader )
    {
        switch ( aIsLeader )
        {
            case mtk::Leader_Follower::LEADER:
            {
                return mLeaderFIManager;
            }

            case mtk::Leader_Follower::FOLLOWER:
            {
                return mFollowerFIManager;
            }

            default:
            {
                MORIS_ERROR( false, "Stabilization_Parameter::get_field_interpolator_manager - can only be leader or follower" );
                return mLeaderFIManager;
            }
        }
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< MSI::Dof_Type > >&
    Stabilization_Parameter::get_global_dof_type_list(
            mtk::Leader_Follower aIsLeader )
    {
        if ( mGlobalDofBuild )
        {
            // build the stabilization parameter global dof type list
            this->build_global_dof_type_list();

            // build the stabilization parameter global dof type map
            this->build_global_dof_type_map();

            // update build flag
            mGlobalDofBuild = false;
        }

        // switch on leader/follower
        switch ( aIsLeader )
        {
            // if leader
            case mtk::Leader_Follower::LEADER:
            {
                // return leader global dof type list
                return mLeaderGlobalDofTypes;
            }
            // if follower
            case mtk::Leader_Follower::FOLLOWER:
            {
                // return follower global dof type list
                return mFollowerGlobalDofTypes;
            }
            // if none
            default:
            {
                MORIS_ASSERT( false, "Stabilization_Parameter::get_global_dof_type_list - can only be leader or follower." );
                return mLeaderGlobalDofTypes;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::get_non_unique_dof_types(
            Vector< MSI::Dof_Type >& aDofTypes )
    {
        // set the size of the dof type list for the set
        uint tCounter = 0;

        // loop over direct leader dof dependencies
        for ( uint iDOF = 0; iDOF < mLeaderDofTypes.size(); iDOF++ )
        {
            // update counter
            tCounter += mLeaderDofTypes( iDOF ).size();
        }

        // get number of direct follower dof dependencies
        for ( uint iDOF = 0; iDOF < mFollowerDofTypes.size(); iDOF++ )
        {
            // update counter
            tCounter += mFollowerDofTypes( iDOF ).size();
        }

        // loop over the leader properties
        for ( const std::shared_ptr< Property >& tProperty : mLeaderProp )
        {
            if ( tProperty != nullptr )
            {
                // get property non unique dof type list
                Vector< MSI::Dof_Type > tActiveDofType;
                tProperty->get_non_unique_dof_types( tActiveDofType );

                // update counter
                tCounter += tActiveDofType.size();
            }
        }

        // loop over follower properties
        for ( const std::shared_ptr< Property >& tProperty : mFollowerProp )
        {
            if ( tProperty != nullptr )
            {
                // get property nn unique dof type list
                Vector< MSI::Dof_Type > tActiveDofType;
                tProperty->get_non_unique_dof_types( tActiveDofType );

                // update counter
                tCounter += tActiveDofType.size();
            }
        }

        // loop over leader constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
        {
            if ( tCM != nullptr )
            {
                // get CM non unique dof type list
                Vector< MSI::Dof_Type > tActiveDofType;
                tCM->get_non_unique_dof_types( tActiveDofType );

                // update counter
                tCounter += tActiveDofType.size();
            }
        }

        // loop over follower constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mFollowerCM )
        {
            if ( tCM != nullptr )
            {
                // get CM non unique dof type list
                Vector< MSI::Dof_Type > tActiveDofType;
                tCM->get_non_unique_dof_types( tActiveDofType );

                // update counter
                tCounter += tActiveDofType.size();
            }
        }

        // reserve memory for dof type list
        aDofTypes.reserve( tCounter );

        // loop over leader dof direct dependencies
        for ( uint iDOF = 0; iDOF < mLeaderDofTypes.size(); iDOF++ )
        {
            // populate the dof list
            aDofTypes.append( mLeaderDofTypes( iDOF ) );
        }

        // loop over follower dof direct dependencies
        for ( uint iDOF = 0; iDOF < mFollowerDofTypes.size(); iDOF++ )
        {
            // populate the dof list
            aDofTypes.append( mFollowerDofTypes( iDOF ) );
        }

        // loop over leader properties
        for ( const std::shared_ptr< Property >& tProperty : mLeaderProp )
        {
            if ( tProperty != nullptr )
            {
                // get property non unique dof type list
                Vector< MSI::Dof_Type > tActiveDofType;
                tProperty->get_non_unique_dof_types( tActiveDofType );

                // populate the dof list
                aDofTypes.append( tActiveDofType );
            }
        }

        // loop over follower properties
        for ( const std::shared_ptr< Property >& tProperty : mFollowerProp )
        {
            if ( tProperty != nullptr )
            {
                // get property non unique dof type list
                Vector< MSI::Dof_Type > tActiveDofType;
                tProperty->get_non_unique_dof_types( tActiveDofType );

                // populate the dof list
                aDofTypes.append( tActiveDofType );
            }
        }

        // loop over the leader constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
        {
            if ( tCM != nullptr )
            {
                // get CM non unique dof type list
                Vector< MSI::Dof_Type > tActiveDofType;
                tCM->get_non_unique_dof_types( tActiveDofType );

                // populate the dof list
                aDofTypes.append( tActiveDofType );
            }
        }

        // loop over the follower constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mFollowerCM )
        {
            if ( tCM != nullptr )
            {
                // get CM non unique dof type list
                Vector< MSI::Dof_Type > tActiveDofType;
                tCM->get_non_unique_dof_types( tActiveDofType );

                // populate the dof list
                aDofTypes.append( tActiveDofType );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::get_non_unique_dof_and_dv_types(
            Vector< MSI::Dof_Type >& aDofTypes,
            Vector< gen::PDV_Type >& aDvTypes )
    {
        // init dof and dv counters
        uint tDofCounter = 0;
        uint tDvCounter  = 0;

        // loop over direct leader dof dependencies
        for ( uint iDof = 0; iDof < mLeaderDofTypes.size(); iDof++ )
        {
            // update counter
            tDofCounter += mLeaderDofTypes( iDof ).size();
        }

        // loop over direct leader dv dependencies
        for ( uint iDv = 0; iDv < mLeaderDvTypes.size(); iDv++ )
        {
            // update counter
            tDvCounter += mLeaderDvTypes( iDv ).size();
        }

        // get number of direct follower dof dependencies
        for ( uint iDof = 0; iDof < mFollowerDofTypes.size(); iDof++ )
        {
            // update counter
            tDofCounter += mFollowerDofTypes( iDof ).size();
        }

        // get number of direct follower dv dependencies
        for ( uint iDv = 0; iDv < mFollowerDvTypes.size(); iDv++ )
        {
            // update counter
            tDvCounter += mFollowerDvTypes( iDv ).size();
        }

        // loop over the leader properties
        for ( const std::shared_ptr< Property >& tProperty : mLeaderProp )
        {
            if ( tProperty != nullptr )
            {
                // get property non unique dof and dv types
                Vector< MSI::Dof_Type >   tActiveDofTypes;
                Vector< gen::PDV_Type >   tActiveDvTypes;
                Vector< mtk::Field_Type > tActiveFieldTypes;

                tProperty->get_non_unique_dof_dv_and_field_types(
                        tActiveDofTypes,
                        tActiveDvTypes,
                        tActiveFieldTypes );

                // update counters
                tDofCounter += tActiveDofTypes.size();
                tDvCounter += tActiveDvTypes.size();
            }
        }

        // loop over follower properties
        for ( const std::shared_ptr< Property >& tProperty : mFollowerProp )
        {
            if ( tProperty != nullptr )
            {
                // get property non unique dof and dv types
                Vector< MSI::Dof_Type >   tActiveDofTypes;
                Vector< gen::PDV_Type >   tActiveDvTypes;
                Vector< mtk::Field_Type > tActiveFieldTypes;

                tProperty->get_non_unique_dof_dv_and_field_types(
                        tActiveDofTypes,
                        tActiveDvTypes,
                        tActiveFieldTypes );

                // update counters
                tDofCounter += tActiveDofTypes.size();
                tDvCounter += tActiveDvTypes.size();
            }
        }

        // loop over leader constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
        {
            if ( tCM != nullptr )
            {
                // get CM non unique dof and dv types
                Vector< MSI::Dof_Type >   tActiveDofTypes;
                Vector< gen::PDV_Type >   tActiveDvTypes;
                Vector< mtk::Field_Type > tActiveFieldTypes;

                tCM->get_non_unique_dof_dv_and_field_types(
                        tActiveDofTypes,
                        tActiveDvTypes,
                        tActiveFieldTypes );

                // update counters
                tDofCounter += tActiveDofTypes.size();
                tDvCounter += tActiveDvTypes.size();
            }
        }

        // loop over follower constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mFollowerCM )
        {
            if ( tCM != nullptr )
            {
                // get CM non unique dof and dv types
                Vector< MSI::Dof_Type >   tActiveDofTypes;
                Vector< gen::PDV_Type >   tActiveDvTypes;
                Vector< mtk::Field_Type > tActiveFieldTypes;

                tCM->get_non_unique_dof_dv_and_field_types(
                        tActiveDofTypes,
                        tActiveDvTypes,
                        tActiveFieldTypes );

                // update counters
                tDofCounter += tActiveDofTypes.size();
                tDvCounter += tActiveDvTypes.size();
            }
        }

        // reserve memory for dof and dv type lists
        aDofTypes.reserve( tDofCounter );
        aDvTypes.reserve( tDvCounter );

        // loop over leader dof direct dependencies
        for ( uint iDof = 0; iDof < mLeaderDofTypes.size(); iDof++ )
        {
            // populate the dof list
            aDofTypes.append( mLeaderDofTypes( iDof ) );
        }

        // loop over leader dv direct dependencies
        for ( uint iDv = 0; iDv < mLeaderDvTypes.size(); iDv++ )
        {
            // populate the dv list
            aDvTypes.append( mLeaderDvTypes( iDv ) );
        }

        // loop over follower dof direct dependencies
        for ( uint iDof = 0; iDof < mFollowerDofTypes.size(); iDof++ )
        {
            // populate the dof list
            aDofTypes.append( mFollowerDofTypes( iDof ) );
        }

        // loop over follower dv direct dependencies
        for ( uint iDv = 0; iDv < mFollowerDvTypes.size(); iDv++ )
        {
            // populate the dv list
            aDvTypes.append( mFollowerDvTypes( iDv ) );
        }

        // loop over leader properties
        for ( const std::shared_ptr< Property >& tProperty : mLeaderProp )
        {
            if ( tProperty != nullptr )
            {
                // get property non unique dof and dv types
                Vector< MSI::Dof_Type >   tActiveDofTypes;
                Vector< gen::PDV_Type >   tActiveDvTypes;
                Vector< mtk::Field_Type > tActiveFieldTypes;

                tProperty->get_non_unique_dof_dv_and_field_types(
                        tActiveDofTypes,
                        tActiveDvTypes,
                        tActiveFieldTypes );

                // populate the dof list
                aDofTypes.append( tActiveDofTypes );
                aDvTypes.append( tActiveDvTypes );
            }
        }

        // loop over follower properties
        for ( const std::shared_ptr< Property >& tProperty : mFollowerProp )
        {
            if ( tProperty != nullptr )
            {
                // get property non unique dof and dv types
                Vector< MSI::Dof_Type >   tActiveDofTypes;
                Vector< gen::PDV_Type >   tActiveDvTypes;
                Vector< mtk::Field_Type > tActiveFieldTypes;

                tProperty->get_non_unique_dof_dv_and_field_types(
                        tActiveDofTypes,
                        tActiveDvTypes,
                        tActiveFieldTypes );

                // populate the dof list
                aDofTypes.append( tActiveDofTypes );
                aDvTypes.append( tActiveDvTypes );
            }
        }

        // loop over the leader constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
        {
            if ( tCM != nullptr )
            {
                // get CM non unique dof and dv types
                Vector< MSI::Dof_Type >   tActiveDofTypes;
                Vector< gen::PDV_Type >   tActiveDvTypes;
                Vector< mtk::Field_Type > tActiveFieldTypes;

                tCM->get_non_unique_dof_dv_and_field_types(
                        tActiveDofTypes,
                        tActiveDvTypes,
                        tActiveFieldTypes );

                // populate the dof list
                aDofTypes.append( tActiveDofTypes );
                aDvTypes.append( tActiveDvTypes );
            }
        }

        // loop over the follower constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mFollowerCM )
        {
            if ( tCM != nullptr )
            {
                // get CM non unique dof and dv types
                Vector< MSI::Dof_Type >   tActiveDofTypes;
                Vector< gen::PDV_Type >   tActiveDvTypes;
                Vector< mtk::Field_Type > tActiveFieldTypes;

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

    void
    Stabilization_Parameter::build_global_dof_type_list()
    {
        // LEADER-------------------------------------------------------
        // get the size of the dof type list
        uint tCounterMax = 0;

        // get number of dof types from penalty parameter
        tCounterMax += mLeaderDofTypes.size();

        // get number of dof types from properties
        for ( const std::shared_ptr< Property >& tProperty : mLeaderProp )
        {
            if ( tProperty != nullptr )
            {
                tCounterMax += tProperty->get_dof_type_list().size();
            }
        }

        // get number of dof types from constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
        {
            if ( tCM != nullptr )
            {
                tCounterMax += tCM->get_global_dof_type_list().size();
            }
        }

        // set size for the global dof type list
        mLeaderGlobalDofTypes.resize( tCounterMax );

        // set a size for the checkList (used to avoid repeating a dof type)
        Vector< sint > tCheckList( tCounterMax, -1 );

        // init total dof counter
        uint tCounter = 0;

        // get dof type from penalty parameter
        for ( uint iDOF = 0; iDOF < mLeaderDofTypes.size(); iDOF++ )
        {
            // put the dof type in the checklist
            tCheckList( tCounter ) = static_cast< uint >( mLeaderDofTypes( iDOF )( 0 ) );

            // put the dof type in the global type list
            mLeaderGlobalDofTypes( tCounter ) = mLeaderDofTypes( iDOF );

            // update the dof counter
            tCounter++;
        }

        // get dof type from properties
        for ( const std::shared_ptr< Property >& tProperty : mLeaderProp )
        {
            if ( tProperty != nullptr )
            {
                // get dof types for property
                const Vector< Vector< MSI::Dof_Type > >& tActiveDofType =
                        tProperty->get_dof_type_list();

                // loop on property dof type
                for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                {
                    // check enum is not already in the list
                    bool tCheck = false;
                    for ( uint i = 0; i < tCounter; i++ )
                    {
                        tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDofType( iDOF )( 0 ) ) );
                    }

                    // if dof enum not in the list
                    if ( !tCheck )
                    {
                        // put the dof type in the checklist
                        tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                        // put the dof type in the global type list
                        mLeaderGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                        // update dof counter
                        tCounter++;
                    }
                }
            }
        }

        // get dof type from constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
        {
            if ( tCM != nullptr )
            {
                // get dof types for constitutive model
                const Vector< Vector< MSI::Dof_Type > >& tActiveDofType =
                        tCM->get_global_dof_type_list();

                // loop on property dof type
                for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                {
                    // check enum is not already in the list
                    bool tCheck = false;
                    for ( uint i = 0; i < tCounter; i++ )
                    {
                        tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDofType( iDOF )( 0 ) ) );
                    }

                    // if dof enum not in the list
                    if ( !tCheck )
                    {
                        // put the dof type in the checklist
                        tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                        // put the dof type in the global type list
                        mLeaderGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                        // update dof counter
                        tCounter++;
                    }
                }
            }
        }

        // get the number of unique dof type groups for the penalty parameter
        mLeaderGlobalDofTypes.resize( tCounter );

        // FOLLOWER--------------------------------------------------------
        // get the size of the dof type list
        tCounterMax = 0;

        // get number of dof types from penalty parameter
        tCounterMax += mFollowerDofTypes.size();

        // get number of dof types from properties
        for ( const std::shared_ptr< Property >& tProperty : mFollowerProp )
        {
            if ( tProperty != nullptr )
            {
                tCounterMax += tProperty->get_dof_type_list().size();
            }
        }

        // get number of dof types from constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mFollowerCM )
        {
            if ( tCM != nullptr )
            {
                tCounterMax += tCM->get_global_dof_type_list().size();
            }
        }

        // set size for the global dof type list
        mFollowerGlobalDofTypes.resize( tCounterMax );

        // set a size for the checkList (used to avoid repeating a dof type)
        tCheckList.resize( tCounterMax, -1 );

        // init total dof counter
        tCounter = 0;

        // get dof type from penalty parameter
        for ( uint iDOF = 0; iDOF < mFollowerDofTypes.size(); iDOF++ )
        {
            // put the dof type in the checklist
            tCheckList( tCounter ) = static_cast< uint >( mFollowerDofTypes( iDOF )( 0 ) );

            // put the dof type in the global type list
            mFollowerGlobalDofTypes( tCounter ) = mFollowerDofTypes( iDOF );

            // update the dof counter
            tCounter++;
        }

        // get dof type from properties
        for ( const std::shared_ptr< Property >& tProperty : mFollowerProp )
        {
            if ( tProperty != nullptr )
            {
                // get dof types for property
                const Vector< Vector< MSI::Dof_Type > >& tActiveDofType =
                        tProperty->get_dof_type_list();

                // loop on property dof type
                for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                {
                    // check enum is not already in the list
                    bool tCheck = false;
                    for ( uint i = 0; i < tCounter; i++ )
                    {
                        tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDofType( iDOF )( 0 ) ) );
                    }

                    // if dof enum not in the list
                    if ( !tCheck )
                    {
                        // put the dof type in the checklist
                        tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                        // put the dof type in the global type list
                        mFollowerGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                        // update dof counter
                        tCounter++;
                    }
                }
            }
        }

        // get dof type from constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mFollowerCM )
        {
            if ( tCM != nullptr )
            {
                // get dof types for constitutive model
                const Vector< Vector< MSI::Dof_Type > >& tActiveDofType =
                        tCM->get_global_dof_type_list();

                // loop on property dof type
                for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                {
                    // check enum is not already in the list
                    bool tCheck = false;
                    for ( uint i = 0; i < tCounter; i++ )
                    {
                        tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDofType( iDOF )( 0 ) ) );
                    }

                    // if dof enum not in the list
                    if ( !tCheck )
                    {
                        // put the dof type in the checklist
                        tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                        // put the dof type in the global type list
                        mFollowerGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                        // update dof counter
                        tCounter++;
                    }
                }
            }
        }

        // get the number of unique dof type groups for the penalty parameter
        mFollowerGlobalDofTypes.resize( tCounter );

        // number of global leader and follower dof types
        uint tNumLeaderGlobalDofTypes   = mLeaderGlobalDofTypes.size();
        uint tNumFollowerGlobalDofTypes = mFollowerGlobalDofTypes.size();

        // set flag for evaluation
        mdPPdLeaderDofEval.set_size( tNumLeaderGlobalDofTypes, 1, true );
        mdPPdFollowerDofEval.set_size( tNumFollowerGlobalDofTypes, 1, true );

        // set storage for evaluation
        mdPPdLeaderDof.resize( tNumLeaderGlobalDofTypes );
        mdPPdFollowerDof.resize( tNumFollowerGlobalDofTypes );
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::build_global_dof_type_map()
    {
        // LEADER-------------------------------------------------------
        // get number of global dof types
        uint tNumDofTypes = mLeaderGlobalDofTypes.size();

        // determine the max Dof_Type enum
        sint tMaxEnum = 0;
        for ( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
        {
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mLeaderGlobalDofTypes( iDOF )( 0 ) ) );
        }
        tMaxEnum++;

        // set the Dof_Type map size
        mLeaderGlobalDofTypeMap.set_size( tMaxEnum, 1, -1 );

        // fill the Dof_Type map
        for ( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
        {
            // fill the property map
            mLeaderGlobalDofTypeMap( static_cast< int >( mLeaderGlobalDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
        }

        // FOLLOWER-------------------------------------------------------
        // get number of global dof types
        tNumDofTypes = mFollowerGlobalDofTypes.size();

        // determine the max Dof_Type enum
        tMaxEnum = 0;
        for ( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
        {
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mFollowerGlobalDofTypes( iDOF )( 0 ) ) );
        }
        tMaxEnum++;

        // set the dof type map size
        mFollowerGlobalDofTypeMap.set_size( tMaxEnum, 1, -1 );

        // fill the dof type map
        for ( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
        {
            // fill the property map
            mFollowerGlobalDofTypeMap( static_cast< int >( mFollowerGlobalDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
        }
    }

    //------------------------------------------------------------------------------

    const Matrix< DDSMat >&
    Stabilization_Parameter::get_global_dof_type_map(
            mtk::Leader_Follower aIsLeader )
    {
        // switch on leader/follower
        switch ( aIsLeader )
        {
            // if leader
            case mtk::Leader_Follower::LEADER:
            {
                // return leader global dof type map
                return mLeaderGlobalDofTypeMap;
            }
            // if follower
            case mtk::Leader_Follower::FOLLOWER:
            {
                // return follower global dof type map
                return mFollowerGlobalDofTypeMap;
            }
            // if none
            default:
            {
                MORIS_ASSERT( false, "Stabilization_Parameter::get_global_dof_type_map - can only be leader or follower." );
                return mLeaderGlobalDofTypeMap;
            }
        }
    }

    //------------------------------------------------------------------------------

    bool
    Stabilization_Parameter::check_dof_dependency(
            const Vector< MSI::Dof_Type >& aDofType,
            mtk::Leader_Follower           aIsLeader )
    {
        // set bool for dependency
        bool tDofDependency = false;

        // get dof type index
        uint tDofIndex = static_cast< uint >( aDofType( 0 ) );

        // if aDofType is an active dof type for the stabilization parameter
        if ( tDofIndex < this->get_global_dof_type_map( aIsLeader ).numel()
                && this->get_global_dof_type_map( aIsLeader )( tDofIndex ) != -1 )
        {
            // bool is set to true
            tDofDependency = true;
        }
        // return bool for dependency
        return tDofDependency;
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< gen::PDV_Type > >&
    Stabilization_Parameter::get_global_dv_type_list(
            mtk::Leader_Follower aIsLeader )
    {
        if ( mGlobalDvBuild )
        {
            // build the stabilization parameter global dof type list
            this->build_global_dv_type_list();

            // build the stabilization parameter global dof type map
            this->build_global_dv_type_map();

            // update build flag
            mGlobalDvBuild = false;
        }

        // switch on leader/follower
        switch ( aIsLeader )
        {
            // if leader
            case mtk::Leader_Follower::LEADER:
            {
                // return leader global dv type list
                return mLeaderGlobalDvTypes;
            }
            // if follower
            case mtk::Leader_Follower::FOLLOWER:
            {
                // return follower global dv type list
                return mFollowerGlobalDvTypes;
            }
            // if none
            default:
            {
                MORIS_ASSERT( false, "Stabilization_Parameter::get_global_dv_type_list - can only be leader or follower." );
                return mLeaderGlobalDvTypes;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::build_global_dv_type_list()
    {
        // LEADER-------------------------------------------------------
        // get the size of the dv type list
        uint tCounterMax = 0;

        // get number of dv types from penalty parameter
        tCounterMax += mLeaderDvTypes.size();

        // get number of dv types from properties
        for ( const std::shared_ptr< Property >& tProperty : mLeaderProp )
        {
            if ( tProperty != nullptr )
            {
                tCounterMax += tProperty->get_dv_type_list().size();
            }
        }

        // get number of dof types from constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
        {
            if ( tCM != nullptr )
            {
                tCounterMax += tCM->get_global_dv_type_list().size();
            }
        }

        // set size for the global dv type list
        mLeaderGlobalDvTypes.resize( tCounterMax );

        // set a size for the checkList (used to avoid repeating a dv type)
        Vector< sint > tCheckList( tCounterMax, -1 );

        // init total dv counter
        uint tCounter = 0;

        // get dv type from penalty parameter
        for ( uint iDv = 0; iDv < mLeaderDvTypes.size(); iDv++ )
        {
            // put the dv type in the checklist
            tCheckList( tCounter ) = static_cast< uint >( mLeaderDvTypes( iDv )( 0 ) );

            // put the dv type in the global type list
            mLeaderGlobalDvTypes( tCounter ) = mLeaderDvTypes( iDv );

            // update the dv counter
            tCounter++;
        }

        // get dv type from properties
        for ( const std::shared_ptr< Property >& tProperty : mLeaderProp )
        {
            if ( tProperty != nullptr )
            {
                // get dv types for property
                const Vector< Vector< gen::PDV_Type > >& tActiveDvType =
                        tProperty->get_dv_type_list();

                // loop on property dv type
                for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                {
                    // check enum is not already in the list
                    bool tCheck = false;
                    for ( uint i = 0; i < tCounter; i++ )
                    {
                        tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                    }

                    // if dof enum not in the list
                    if ( !tCheck )
                    {
                        // put the dv type in the checklist
                        tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                        // put the dv type in the global type list
                        mLeaderGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                        // update dof counter
                        tCounter++;
                    }
                }
            }
        }

        // get dof type from constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
        {
            if ( tCM != nullptr )
            {
                // get dof types for constitutive model
                const Vector< Vector< gen::PDV_Type > >& tActiveDvType =
                        tCM->get_global_dv_type_list();

                // loop on property dv type
                for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                {
                    // check enum is not already in the list
                    bool tCheck = false;
                    for ( uint i = 0; i < tCounter; i++ )
                    {
                        tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                    }

                    // if dv enum not in the list
                    if ( !tCheck )
                    {
                        // put the dv type in the checklist
                        tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                        // put the dv type in the global type list
                        mLeaderGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                        // update dv counter
                        tCounter++;
                    }
                }
            }
        }

        // get the number of unique dv type groups for the penalty parameter
        mLeaderGlobalDvTypes.resize( tCounter );

        // FOLLOWER--------------------------------------------------------
        // get the size of the dv type list
        tCounterMax = 0;

        // get number of dv types from penalty parameter
        tCounterMax += mFollowerDvTypes.size();

        // get number of dv types from properties
        for ( const std::shared_ptr< Property >& tProperty : mFollowerProp )
        {
            if ( tProperty != nullptr )
            {
                tCounterMax += tProperty->get_dv_type_list().size();
            }
        }

        // get number of dv types from constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mFollowerCM )
        {
            if ( tCM != nullptr )
            {
                tCounterMax += tCM->get_global_dv_type_list().size();
            }
        }

        // set size for the global dv type list
        mFollowerGlobalDvTypes.resize( tCounterMax );

        // set a size for the checkList (used to avoid repeating a dv type)
        tCheckList.resize( tCounterMax, -1 );

        // init total dv counter
        tCounter = 0;

        // get dv type from penalty parameter
        for ( uint iDv = 0; iDv < mFollowerDvTypes.size(); iDv++ )
        {
            // put the dv type in the checklist
            tCheckList( tCounter ) = static_cast< uint >( mFollowerDvTypes( iDv )( 0 ) );

            // put the dv type in the global type list
            mFollowerGlobalDvTypes( tCounter ) = mFollowerDvTypes( iDv );

            // update the dv counter
            tCounter++;
        }

        // get dv type from properties
        for ( const std::shared_ptr< Property >& tProperty : mFollowerProp )
        {
            if ( tProperty != nullptr )
            {
                // get dv types for property
                const Vector< Vector< gen::PDV_Type > >& tActiveDvType =
                        tProperty->get_dv_type_list();

                // loop on property dv type
                for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                {
                    // check enum is not already in the list
                    bool tCheck = false;
                    for ( uint i = 0; i < tCounter; i++ )
                    {
                        tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                    }

                    // if dv enum not in the list
                    if ( !tCheck )
                    {
                        // put the dv type in the checklist
                        tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                        // put the dv type in the global type list
                        mFollowerGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                        // update dv counter
                        tCounter++;
                    }
                }
            }
        }

        // get dv type from constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
        {
            if ( tCM != nullptr )
            {
                // get dv types for constitutive model
                const Vector< Vector< gen::PDV_Type > >& tActiveDvType =
                        tCM->get_global_dv_type_list();

                // loop on property dv type
                for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                {
                    // check enum is not already in the list
                    bool tCheck = false;
                    for ( uint i = 0; i < tCounter; i++ )
                    {
                        tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                    }

                    // if dv enum not in the list
                    if ( !tCheck )
                    {
                        // put the dv type in the checklist
                        tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                        // put the dv type in the global type list
                        mFollowerGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                        // update dv counter
                        tCounter++;
                    }
                }
            }
        }

        // get the number of unique dv type groups for the penalty parameter
        mFollowerGlobalDvTypes.resize( tCounter );

        // build global dv type map
        this->build_global_dv_type_map();

        // number of global leader and follower dv types
        uint tNumLeaderGlobalDvTypes   = mLeaderGlobalDvTypes.size();
        uint tNumFollowerGlobalDvTypes = mFollowerGlobalDvTypes.size();

        // set flag for evaluation
        mdPPdLeaderDvEval.set_size( tNumLeaderGlobalDvTypes, 1, true );
        mdPPdFollowerDvEval.set_size( tNumFollowerGlobalDvTypes, 1, true );

        // set storage for evaluation
        mdPPdLeaderDv.resize( tNumLeaderGlobalDvTypes );
        mdPPdFollowerDv.resize( tNumFollowerGlobalDvTypes );
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::build_global_dv_type_map()
    {
        // LEADER-------------------------------------------------------
        // get number of global dof types
        uint tNumDvTypes = mLeaderGlobalDvTypes.size();

        // determine the max Dv_Type enum
        sint tMaxEnum = 0;
        for ( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
        {
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mLeaderGlobalDvTypes( iDv )( 0 ) ) );
        }
        tMaxEnum++;

        // set the Dv_Type map size
        mLeaderGlobalDvTypeMap.set_size( tMaxEnum, 1, -1 );

        // fill the Dv_Type map
        for ( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
        {
            // fill the property map
            mLeaderGlobalDvTypeMap( static_cast< int >( mLeaderGlobalDvTypes( iDv )( 0 ) ), 0 ) = iDv;
        }

        // FOLLOWER-------------------------------------------------------
        // get number of global dv types
        tNumDvTypes = mFollowerGlobalDvTypes.size();

        // determine the max Dv_Type enum
        tMaxEnum = 0;
        for ( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
        {
            tMaxEnum = std::max( tMaxEnum, static_cast< int >( mFollowerGlobalDvTypes( iDv )( 0 ) ) );
        }
        tMaxEnum++;

        // set the dv type map size
        mFollowerGlobalDvTypeMap.set_size( tMaxEnum, 1, -1 );

        // fill the dv type map
        for ( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
        {
            // fill the property map
            mFollowerGlobalDvTypeMap( static_cast< int >( mFollowerGlobalDvTypes( iDv )( 0 ) ), 0 ) = iDv;
        }
    }

    //------------------------------------------------------------------------------

    bool
    Stabilization_Parameter::check_leader_dv_dependency(
            const Vector< gen::PDV_Type >& aDvType )
    {
        // set bool for dependency
        bool tDvDependency = false;

        // get dv type index
        uint tDvIndex = static_cast< uint >( aDvType( 0 ) );

        // if aDvType is an active dv type for the constitutive model
        if ( tDvIndex < mLeaderGlobalDvTypeMap.numel() && mLeaderGlobalDvTypeMap( tDvIndex ) != -1 )
        {
            // bool is set to true
            tDvDependency = true;
        }
        // return bool for dependency
        return tDvDependency;
    }

    //------------------------------------------------------------------------------

    bool
    Stabilization_Parameter::check_follower_dv_dependency(
            const Vector< gen::PDV_Type >& aDvType )
    {
        // set bool for dependency
        bool tDvDependency = false;

        // get dv type index
        uint tDvIndex = static_cast< uint >( aDvType( 0 ) );

        // if aDvType is an active dv type for the constitutive model
        if ( tDvIndex < mFollowerGlobalDvTypeMap.numel() && mFollowerGlobalDvTypeMap( tDvIndex ) != -1 )
        {
            // bool is set to true
            tDvDependency = true;
        }
        // return bool for dependency
        return tDvDependency;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Stabilization_Parameter::val()
    {
        // if the penalty parameter was not evaluated
        if ( mPPEval )
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

    const Matrix< DDRMat >&
    Stabilization_Parameter::dSPdLeaderDOF(
            const Vector< MSI::Dof_Type >& aDofType )
    {
        // if aDofType is not an active dof type for the property
        MORIS_ERROR(
                this->check_dof_dependency( aDofType, mtk::Leader_Follower::LEADER ),
                "Stabilization_Parameter::dSPdLeaderDOF - no dependency on this dof type." );

        // get the dof index
        uint tDofIndex = mLeaderGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdPPdLeaderDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dSPdLeaderDOF( aDofType );

            // set bool for evaluation
            mdPPdLeaderDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mdPPdLeaderDof( tDofIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Stabilization_Parameter::dSPdFollowerDOF(
            const Vector< MSI::Dof_Type >& aDofType )
    {
        // if aDofType is not an active dof type for the property
        MORIS_ERROR(
                this->check_dof_dependency( aDofType, mtk::Leader_Follower::FOLLOWER ),
                "Stabilization_Parameter::dSPdFollowerDOF - no dependency on this dof type." );

        // get the dof index
        uint tDofIndex = mFollowerGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdPPdFollowerDofEval( tDofIndex ) )
        {
            // evaluate the derivative
            this->eval_dSPdFollowerDOF( aDofType );

            // set bool for evaluation
            mdPPdFollowerDofEval( tDofIndex ) = false;
        }

        // return the derivative
        return mdPPdFollowerDof( tDofIndex );
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::eval_dSPdLeaderDOF_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            Matrix< DDRMat >&              adSPdDOF_FD,
            real                           aPerturbation,
            fem::FDScheme_Type             aFDSchemeType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the dof index
        uint tDofIndex = mLeaderGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFI =
                mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of leader dofs wrt which derivative is computed
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
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed traction contribution to dtractiondu
                    adSPdDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tUnperturbedSPValue / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
                {
                    // reset the perturbed coefficients
                    Matrix< DDRMat > tCoeffPert = tCoeff;

                    // perturb the coefficient
                    tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                    // set the perturbed coefficients to FI
                    tFI->set_coeff( tCoeffPert );

                    // reset properties
                    this->reset_eval_flags();

                    // assemble derivatiev of SP wrt leader dof type
                    adSPdDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->val() / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }
                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );

        // set value to storage
        mdPPdLeaderDof( tDofIndex ) = adSPdDOF_FD;
    }

    //------------------------------------------------------------------------------

    void
    Stabilization_Parameter::eval_dSPdFollowerDOF_FD(
            const Vector< MSI::Dof_Type >& aDofTypes,
            Matrix< DDRMat >&              adSPdDOF_FD,
            real                           aPerturbation,
            fem::FDScheme_Type             aFDSchemeType )
    {
        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumPoints = tFDScheme( 0 ).size();

        // get the dof index
        uint tDofIndex = mFollowerGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // get the field interpolator for type
        Field_Interpolator* tFI =
                mFollowerFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // get number of leader dofs wrt which derivative is computed
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
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
        {
            // loop over coefficients rows
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // compute the perturbation absolute value
                real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                // check that perturbation is not zero
                if ( std::abs( tDeltaH ) < 1e-12 )
                {
                    tDeltaH = aPerturbation;
                }

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed traction contribution to dtractiondu
                    adSPdDOF_FD.get_column( tDofCounter ) +=
                            tFDScheme( 1 )( 0 ) * tUnperturbedSPValue / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over the points for FD
                for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
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
                    mdPPdFollowerDof( tDofIndex ).get_column( tDofCounter ) +=
                            tFDScheme( 1 )( iPoint ) * this->val() / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                }
                // update dof counter
                tDofCounter++;
            }
        }
        // reset the coefficients values
        tFI->set_coeff( tCoeff );

        // set value to storage
        mdPPdFollowerDof( tDofIndex ) = adSPdDOF_FD;
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Stabilization_Parameter::dSPdLeaderDV(
            const Vector< gen::PDV_Type >& aDvTypes )
    {
        // if aDofType is not an active dv type for the property
        MORIS_ERROR(
                this->check_leader_dv_dependency( aDvTypes ),
                "Penalty_Parameter::dSPdLeaderDV - no dependency on this dv type." );

        // get the dv index
        uint tDvIndex = mLeaderGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdPPdLeaderDofEval( tDvIndex ) )
        {
            // evaluate the derivative
            this->eval_dSPdLeaderDV( aDvTypes );

            // set bool for evaluation
            mdPPdLeaderDofEval( tDvIndex ) = false;
        }

        // return the derivative
        return mdPPdLeaderDof( tDvIndex );
    }

    //------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Stabilization_Parameter::dSPdFollowerDV(
            const Vector< gen::PDV_Type >& aDvTypes )
    {
        // if aDofType is not an active dv type for the property
        MORIS_ERROR(
                this->check_follower_dv_dependency( aDvTypes ),
                "Stabilization_Parameter::dSPdFollowerDV - no dependency on this dv type." );

        // get the dv index
        uint tDvIndex = mFollowerGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ) );

        // if the derivative has not been evaluated yet
        if ( mdPPdFollowerDvEval( tDvIndex ) )
        {
            // evaluate the derivative
            this->eval_dSPdFollowerDV( aDvTypes );

            // set bool for evaluation
            mdPPdFollowerDvEval( tDvIndex ) = false;
        }

        // return the derivative
        return mdPPdFollowerDv( tDvIndex );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem

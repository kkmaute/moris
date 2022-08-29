/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI.cpp
 *
 */

#include "cl_FEM_IQI.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_norm.hpp"
#include "fn_min.hpp"
#include "fn_max.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        void
        IQI::print_names()
        {
            std::cout << "----------" << std::endl;
            std::cout << "IQI: " << mName << std::endl;

            // properties
            for ( uint iProp = 0; iProp < mMasterProp.size(); iProp++ )
            {
                if ( mMasterProp( iProp ) != nullptr )
                {
                    std::cout << "Master property: " << mMasterProp( iProp )->get_name() << std::endl;
                }
            }
            for ( uint iProp = 0; iProp < mSlaveProp.size(); iProp++ )
            {
                if ( mSlaveProp( iProp ) != nullptr )
                {
                    std::cout << "Slave property:  " << mSlaveProp( iProp )->get_name() << std::endl;
                }
            }

            // CM
            for ( uint iCM = 0; iCM < mMasterCM.size(); iCM++ )
            {
                if ( mMasterCM( iCM ) != nullptr )
                {
                    std::cout << "Master CM: " << mMasterCM( iCM )->get_name() << std::endl;
                }
            }
            for ( uint iCM = 0; iCM < mSlaveCM.size(); iCM++ )
            {
                if ( mSlaveCM( iCM ) != nullptr )
                {
                    std::cout << "Slave CM:  " << mSlaveCM( iCM )->get_name() << std::endl;
                }
            }

            // SP
            for ( uint iSP = 0; iSP < mStabilizationParam.size(); iSP++ )
            {
                if ( mStabilizationParam( iSP ) != nullptr )
                {
                    std::cout << "SP: " << mStabilizationParam( iSP )->get_name() << std::endl;
                }
            }
            std::cout << "----------" << std::endl;
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_function_pointers()
        {
            // switch on element type
            switch ( mSet->get_element_type() )
            {
                case fem::Element_Type::BULK:
                {
                    m_compute_dQIdu_FD          = &IQI::select_dQIdu_FD;
                    m_compute_dQIdp_FD_material = &IQI::select_dQIdp_FD_material;
                    m_compute_dQIdp_FD_geometry = &IQI::select_dQIdp_FD_geometry_bulk;
                    break;
                }
                case fem::Element_Type::SIDESET:
                {
                    m_compute_dQIdu_FD          = &IQI::select_dQIdu_FD;
                    m_compute_dQIdp_FD_material = &IQI::select_dQIdp_FD_material;
                    m_compute_dQIdp_FD_geometry = &IQI::select_dQIdp_FD_geometry_sideset;
                    break;
                }
                case fem::Element_Type::TIME_SIDESET:
                case fem::Element_Type::TIME_BOUNDARY:
                {
                    m_compute_dQIdu_FD          = &IQI::select_dQIdu_FD;
                    m_compute_dQIdp_FD_material = &IQI::select_dQIdp_FD_material;
                    m_compute_dQIdp_FD_geometry = &IQI::select_dQIdp_FD_geometry_bulk;
                    break;
                }
                case fem::Element_Type::DOUBLE_SIDESET:
                {
                    m_compute_dQIdu_FD          = &IQI::select_dQIdu_FD;
                    m_compute_dQIdp_FD_material = &IQI::select_dQIdp_FD_material_double;
                    m_compute_dQIdp_FD_geometry = &IQI::select_dQIdp_FD_geometry_double;
                    break;
                }
                default:
                    MORIS_ERROR( false, "IWG::set_function_pointers - unknown element type." );
            }

            // switch on perturbation strategy
            switch ( mSet->get_perturbation_strategy() )
            {
                case fem::Perturbation_Type::RELATIVE:
                {
                    m_build_perturbation_size = &IQI::build_perturbation_size_relative;
                    break;
                }
                case fem::Perturbation_Type::ABSOLUTE:
                {
                    m_build_perturbation_size = &IQI::build_perturbation_size_absolute;
                    break;
                }
                default:
                    MORIS_ERROR( false, "IQI::set_function_pointers - unknown perturbation type." );
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::reset_eval_flags()
        {
            // reset properties
            for ( const std::shared_ptr< Property >& tProp : mMasterProp )
            {
                if ( tProp != nullptr )
                {
                    tProp->reset_eval_flags();
                }
            }

            for ( const std::shared_ptr< Property >& tProp : mSlaveProp )
            {
                if ( tProp != nullptr )
                {
                    tProp->reset_eval_flags();
                }
            }

            // reset constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    tCM->reset_eval_flags();
                }
            }

            for ( const std::shared_ptr< Constitutive_Model >& tCM : mSlaveCM )
            {
                if ( tCM != nullptr )
                {
                    tCM->reset_eval_flags();
                }
            }

            // reset stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter >& tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    tSP->reset_eval_flags();
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_phase_name(
                std::string       aPhaseName,
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    mMasterPhaseName = aPhaseName;
                    break;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    mSlavePhaseName = aPhaseName;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IWG::set_phase_name - aIsMaster can only be master or slave." );
                }
            }
        }

        //------------------------------------------------------------------------------

        std::string
        IQI::get_phase_name( mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    return mMasterPhaseName;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    return mSlavePhaseName;
                }
                default:
                {
                    MORIS_ERROR( false, "IWG::get_phase_name - aIsMaster can only be master or slave." );
                    return mMasterPhaseName;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_reference_value( real aReferenceValue )
        {
            if ( not mNormalized )
            {
                mReferenceValue = aReferenceValue;
                mNormalized     = true;
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_field_interpolator_manager(
                Field_Interpolator_Manager* aFieldInterpolatorManager,
                mtk::Master_Slave           aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    mMasterFIManager = aFieldInterpolatorManager;
                    break;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    mSlaveFIManager = aFieldInterpolatorManager;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IWG::set_field_interpolator_manager - can only be master or slave" );
                }
            }

            // loop over the the SP
            for ( const std::shared_ptr< Stabilization_Parameter >& tSP : this->get_stabilization_parameters() )
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
            for ( const std::shared_ptr< Constitutive_Model >& tCM : this->get_constitutive_models( aIsMaster ) )
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
            for ( const std::shared_ptr< Property >& tProp : this->get_properties( aIsMaster ) )
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

        Field_Interpolator_Manager*
        IQI::get_field_interpolator_manager(
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                    return mMasterFIManager;

                case mtk::Master_Slave::SLAVE:
                    return mSlaveFIManager;

                default:
                    MORIS_ERROR( false, "IQI::get_field_inetrpolator_manager - can only be master or slave." );
                    return mMasterFIManager;
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_dof_type_list(
                const moris::Cell< moris::Cell< MSI::Dof_Type > >& aDofTypes,
                mtk::Master_Slave                                  aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    mMasterDofTypes = aDofTypes;
                    break;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    mSlaveDofTypes = aDofTypes;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IQI::set_dof_type_list - can only be MASTER or SLAVE." );
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< MSI::Dof_Type > >&
        IQI::get_dof_type_list(
                mtk::Master_Slave aIsMaster )
        {
            // switch on master/slave
            switch ( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER:
                {
                    // return master global dof type list
                    return mMasterDofTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE:
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
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< MSI::Dof_Type > >&
        IQI::get_global_dof_type_list(
                mtk::Master_Slave aIsMaster )
        {
            // if the global list was not yet built
            if ( mGlobalDofBuild )
            {
                // build global dof type list
                this->build_global_dof_dv_and_field_type_list();

                // update build flag
                mGlobalDofBuild   = false;
                mGlobalDvBuild    = false;
                mGlobalFieldBuild = false;
            }

            // switch on master/slave
            switch ( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER:
                {
                    // return master global dof type list
                    return mMasterGlobalDofTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE:
                {
                    // return slave global dof type list
                    return mSlaveGlobalDofTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IQI::get_global_dof_type_list - can only be master or slave." );
                    return mMasterGlobalDofTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_dv_type_list(
                const moris::Cell< moris::Cell< PDV_Type > >& aDvTypes,
                mtk::Master_Slave                             aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    mMasterDvTypes = aDvTypes;
                    break;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    mSlaveDvTypes = aDvTypes;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IQI::set_dv_type_list - can only be MASTER or SLAVE." );
                }
            }
        }

        //------------------------------------------------------------------------------

        moris::Cell< moris::Cell< PDV_Type > >&
        IQI::get_dv_type_list(
                mtk::Master_Slave aIsMaster )
        {
            // switch on master/slave
            switch ( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER:
                {
                    // return master global dv type list
                    return mMasterDvTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE:
                {
                    // return slave global dv type list
                    return mSlaveDvTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IQI::get_dv_type_list - can only be master or slave." );
                    return mMasterDvTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< PDV_Type > >&
        IQI::get_global_dv_type_list(
                mtk::Master_Slave aIsMaster )
        {
            // if the global list was not yet built
            if ( mGlobalDvBuild )
            {
                // build global dv type list
                this->build_global_dof_dv_and_field_type_list();

                // update build flag
                mGlobalDvBuild    = false;
                mGlobalDofBuild   = false;
                mGlobalFieldBuild = false;
            }

            // switch on master/slave
            switch ( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER:
                {
                    // return master global dv type list
                    return mMasterGlobalDvTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE:
                {
                    // return slave global dv type list
                    return mSlaveGlobalDvTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IQI::get_global_dv_type_list - can only be master or slave." );
                    return mMasterGlobalDvTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< mtk::Field_Type > >&
        IQI::get_global_field_type_list(
                mtk::Master_Slave aIsMaster )
        {
            // if the global list was not yet built
            if ( mGlobalFieldBuild )
            {
                // build global dv type list
                this->build_global_dof_dv_and_field_type_list();

                // update build flag
                mGlobalDvBuild    = false;
                mGlobalDofBuild   = false;
                mGlobalFieldBuild = false;
            }

            // switch on master/slave
            switch ( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER:
                {
                    // return master global dof type list
                    return mMasterGlobalFieldTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE:
                {
                    // return slave global dof type list
                    return mSlaveGlobalFieldTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IWG::get_global_field_type_list - can only be master or slave." );
                    return mMasterGlobalFieldTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_field_type_list(
                const moris::Cell< moris::Cell< mtk::Field_Type > >& aDofTypes,
                mtk::Master_Slave                                    aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    mMasterFieldTypes = aDofTypes;
                    break;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    mSlaveFieldTypes = aDofTypes;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IQI::set_dof_type_list - can only be MASTER or SLAVE." );
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< mtk::Field_Type > >&
        IQI::get_field_type_list(
                mtk::Master_Slave aIsMaster ) const
        {
            // switch on master/slave
            switch ( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER:
                {
                    // return master global dof type list
                    return mMasterFieldTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE:
                {
                    // return slave global dof type list
                    return mSlaveFieldTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IQI::get_dof_type_list - can only be master or slave." );
                    return mMasterFieldTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "IQI::set_property - IQI %s - Unknown aPropertyString: %s ",
                    mName.c_str(),
                    aPropertyString.c_str() );

            // set the property in the property pointer cell
            this->get_properties( aIsMaster )( mPropertyMap[ aPropertyString ] ) = aProperty;
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< fem::Property > >&
        IQI::get_properties(
                mtk::Master_Slave aIsMaster )
        {
            // switch on master/slave
            switch ( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER:
                {
                    // return master property pointers
                    return mMasterProp;
                }
                // if slave
                case mtk::Master_Slave::SLAVE:
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

        void
        IQI::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster )
        {
            // check that aConstitutiveString makes sense
            MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(),
                    "IQI::set_constitutive_model - IQI %s - Unknown aConstitutiveString: %s ",
                    mName.c_str(),
                    aConstitutiveString.c_str() );

            // set the CM in the CM pointer cell
            this->get_constitutive_models( aIsMaster )( mConstitutiveMap[ aConstitutiveString ] ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< fem::Constitutive_Model > >&
        IQI::get_constitutive_models(
                mtk::Master_Slave aIsMaster )
        {
            // switch on master/slave
            switch ( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER:
                {
                    // return master property pointers
                    return mMasterCM;
                }
                // if slave
                case mtk::Master_Slave::SLAVE:
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

        void
        IQI::set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString )
        {
            // check that aStabilizationString makes sense
            MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(),
                    "IQI::set_stabilization_parameter - IQI %s - Unknown aStabilizationString: %s ",
                    mName.c_str(),
                    aStabilizationString.c_str() );

            // set the stabilization parameter in the stabilization parameter cell
            this->get_stabilization_parameters()( mStabilizationMap[ aStabilizationString ] ) = aStabilizationParameter;

            // set active cluster measure on IQI flag on/off
            mActiveCMEAFlag = mActiveCMEAFlag || ( aStabilizationParameter->get_cluster_measure_tuple_list().size() > 0 );
        }

        //------------------------------------------------------------------------------

        void
        IQI::get_non_unique_dof_dv_and_field_types(
                moris::Cell< moris::Cell< MSI::Dof_Type > >&   aDofTypes,
                moris::Cell< moris::Cell< PDV_Type > >&        aDvTypes,
                moris::Cell< moris::Cell< mtk::Field_Type > >& aFieldTypes )

        {
            // init counters for dof and dv types
            uint tMasterDofCounter   = 0;
            uint tSlaveDofCounter    = 0;
            uint tMasterDvCounter    = 0;
            uint tSlaveDvCounter     = 0;
            uint tMasterFieldCounter = 0;
            uint tSlaveFieldCounter  = 0;

            // get number of direct master dof dependencies
            for ( uint iDof = 0; iDof < mMasterDofTypes.size(); iDof++ )
            {
                tMasterDofCounter += mMasterDofTypes( iDof ).size();
            }

            // get number of direct master dv dependencies
            for ( uint iDv = 0; iDv < mMasterDvTypes.size(); iDv++ )
            {
                tMasterDvCounter += mMasterDvTypes( iDv ).size();
            }

            // get number of direct master field dependencies
            for ( uint iFi = 0; iFi < mMasterFieldTypes.size(); iFi++ )
            {
                tMasterFieldCounter += mMasterFieldTypes( iFi ).size();
            }

            // get number of direct slave dof dependencies
            for ( uint iDof = 0; iDof < mSlaveDofTypes.size(); iDof++ )
            {
                tSlaveDofCounter += mSlaveDofTypes( iDof ).size();
            }

            // get number of direct slave dv dependencies
            for ( uint iDv = 0; iDv < mSlaveDvTypes.size(); iDv++ )
            {
                tSlaveDvCounter += mSlaveDvTypes( iDv ).size();
            }

            // get number of direct slave field dependencies
            for ( uint iFi = 0; iFi < mSlaveFieldTypes.size(); iFi++ )
            {
                tSlaveFieldCounter += mSlaveFieldTypes( iFi ).size();
            }

            // loop over the master properties
            for ( const std::shared_ptr< Property >& tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type list
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tProperty->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // update dof and dv counters
                    tMasterDofCounter += tActiveDofTypes.size();
                    tMasterDvCounter += tActiveDvTypes.size();
                    tMasterFieldCounter += tActiveFieldTypes.size();
                }
            }

            // loop over slave properties
            for ( const std::shared_ptr< Property >& tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tProperty->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // update dof and dv counter
                    tSlaveDofCounter += tActiveDofTypes.size();
                    tSlaveDvCounter += tActiveDvTypes.size();
                    tSlaveFieldCounter += tActiveFieldTypes.size();
                }
            }

            // loop over master constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;
                    tCM->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // update dof and dv counters
                    tMasterDofCounter += tActiveDofTypes.size();
                    tMasterDvCounter += tActiveDvTypes.size();
                    tMasterFieldCounter += tActiveFieldTypes.size();
                }
            }

            // loop over slave constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mSlaveCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;
                    tCM->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // update dof and dv counters
                    tSlaveDofCounter += tActiveDofTypes.size();
                    tSlaveDvCounter += tActiveDvTypes.size();
                    tSlaveFieldCounter += tActiveFieldTypes.size();
                }
            }

            // loop over master stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter >& tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get SP non unique dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< PDV_Type >      tActiveDvTypes;
                    tSP->get_non_unique_dof_and_dv_types(
                            tActiveDofTypes,
                            tActiveDvTypes );

                    // update dof and dv counters
                    tMasterDofCounter += tActiveDofTypes.size();
                    tMasterDvCounter += tActiveDvTypes.size();
                    tSlaveDofCounter += tActiveDofTypes.size();
                    tSlaveDvCounter += tActiveDvTypes.size();
                }
            }

            // reserve memory for dof and dv type lists
            aDofTypes.resize( 2 );
            aDvTypes.resize( 2 );
            aFieldTypes.resize( 2 );

            aDofTypes( 0 ).reserve( tMasterDofCounter );
            aDvTypes( 0 ).reserve( tMasterDvCounter );
            aFieldTypes( 0 ).reserve( tMasterFieldCounter );
            aDofTypes( 1 ).reserve( tSlaveDofCounter );
            aDvTypes( 1 ).reserve( tSlaveDvCounter );
            aFieldTypes( 1 ).reserve( tSlaveFieldCounter );

            // loop over master dof direct dependencies
            for ( uint iDof = 0; iDof < mMasterDofTypes.size(); iDof++ )
            {
                // populate the dof list
                aDofTypes( 0 ).append( mMasterDofTypes( iDof ) );
            }

            // loop over master dv direct dependencies
            for ( uint iDv = 0; iDv < mMasterDvTypes.size(); iDv++ )
            {
                // populate the dv list
                aDvTypes( 0 ).append( mMasterDvTypes( iDv ) );
            }

            // loop over master field direct dependencies
            for ( uint iFi = 0; iFi < mMasterFieldTypes.size(); iFi++ )
            {
                // populate the dv list
                aFieldTypes( 0 ).append( mMasterFieldTypes( iFi ) );
            }

            // loop over slave dof direct dependencies
            for ( uint iDof = 0; iDof < mSlaveDofTypes.size(); iDof++ )
            {
                // populate the dof list
                aDofTypes( 1 ).append( mSlaveDofTypes( iDof ) );
            }

            // loop over slave dv direct dependencies
            for ( uint iDv = 0; iDv < mSlaveDvTypes.size(); iDv++ )
            {
                // populate the dv list
                aDvTypes( 1 ).append( mSlaveDvTypes( iDv ) );
            }

            // loop over slave dv direct dependencies
            for ( uint iFi = 0; iFi < mSlaveFieldTypes.size(); iFi++ )
            {
                aFieldTypes( 1 ).append( mSlaveFieldTypes( iFi ) );
            }

            // loop over master properties
            for ( const std::shared_ptr< Property >& tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type list
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tProperty->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // populate the dof and dv lists
                    aDofTypes( 0 ).append( tActiveDofTypes );
                    aDvTypes( 0 ).append( tActiveDvTypes );
                    aFieldTypes( 0 ).append( tActiveFieldTypes );
                }
            }

            // loop over slave properties
            for ( const std::shared_ptr< Property >& tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type list
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tProperty->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // populate the dof and dv lists
                    aDofTypes( 1 ).append( tActiveDofTypes );
                    aDvTypes( 1 ).append( tActiveDvTypes );
                    aFieldTypes( 1 ).append( tActiveFieldTypes );
                }
            }

            // loop over the master constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tCM->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // populate the dof and dv lists
                    aDofTypes( 0 ).append( tActiveDofTypes );
                    aDvTypes( 0 ).append( tActiveDvTypes );
                    aFieldTypes( 0 ).append( tActiveFieldTypes );
                }
            }

            // loop over the slave constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mSlaveCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type >   tActiveDofTypes;
                    moris::Cell< PDV_Type >        tActiveDvTypes;
                    moris::Cell< mtk::Field_Type > tActiveFieldTypes;

                    tCM->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // populate the dof and dv lists
                    aDofTypes( 1 ).append( tActiveDofTypes );
                    aDvTypes( 1 ).append( tActiveDvTypes );
                    aFieldTypes( 1 ).append( tActiveFieldTypes );
                }
            }

            // FIXME this is potentially problematic since it will add slave dependencies even for bulk elements
            // FIXME Ask lise about it. We could ask the set for the element type. should work for DOUBLE_SIDED.
            // FIXME Whats with time boundary
            // loop over the stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter >& tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get SP non unique master dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< PDV_Type >      tActiveDvTypes;

                    tSP->get_non_unique_dof_and_dv_types(
                            tActiveDofTypes,
                            tActiveDvTypes );

                    // populate the dof and dv lists
                    aDofTypes( 0 ).append( tActiveDofTypes );
                    aDvTypes( 0 ).append( tActiveDvTypes );
                    aDofTypes( 1 ).append( tActiveDofTypes );
                    aDvTypes( 1 ).append( tActiveDvTypes );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::build_global_dof_dv_and_field_type_list()
        {
            // MASTER-------------------------------------------------------
            // get number of dof and dv types on set
            uint tNumDofTypes   = mSet->get_num_unique_dof_types();
            uint tNumDvTypes    = mSet->get_num_unique_dv_types();
            uint tNumFieldTypes = mSet->get_num_unique_field_types();

            // set size for the global dof and dv type lists
            mMasterGlobalDofTypes.reserve( tNumDofTypes );
            mMasterGlobalDvTypes.reserve( tNumDvTypes );
            mMasterGlobalFieldTypes.reserve( tNumFieldTypes );

            // set a size for the dof and dv checkLists
            //( used to avoid repeating a dof or a dv type)
            Matrix< DDSMat > tDofCheckList( tNumDofTypes, 1, -1 );
            Matrix< DDSMat > tDvCheckList( tNumDvTypes, 1, -1 );
            Matrix< DDSMat > tFieldCheckList( tNumFieldTypes, 1, -1 );

            // get dof type from direct dependencies
            for ( uint iDof = 0; iDof < mMasterDofTypes.size(); iDof++ )
            {
                // get set index for dof type
                sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( mMasterDofTypes( iDof )( 0 ) );    // FIXME'

                // put the dof type in the checklist
                tDofCheckList( tDofTypeIndex ) = 1;

                // put the dof type in the global type list
                mMasterGlobalDofTypes.push_back( mMasterDofTypes( iDof ) );
            }

            // get dv type from direct dependencies
            for ( uint iDv = 0; iDv < mMasterDvTypes.size(); iDv++ )
            {
                // get set index for dv type
                sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( mMasterDvTypes( iDv )( 0 ) );    // FIXME'

                // put the dv type in the checklist
                tDvCheckList( tDvTypeIndex ) = 1;

                // put the dv type in the global type list
                mMasterGlobalDvTypes.push_back( mMasterDvTypes( iDv ) );
            }

            // get field type from direct dependencies
            for ( uint iFi = 0; iFi < mMasterFieldTypes.size(); iFi++ )
            {
                // get set index for field type
                sint tFieldTypeIndex = mSet->get_index_from_unique_field_type_map( mMasterFieldTypes( iFi )( 0 ) );    // FIXME'

                // put the field type in the checklist
                tFieldCheckList( tFieldTypeIndex ) = 1;

                // put the field type in the global type list
                mMasterGlobalFieldTypes.push_back( mMasterFieldTypes( iFi ) );
            }

            // get dof type from master properties
            for ( const std::shared_ptr< Property >& tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get dof types for property
                    const moris::Cell< moris::Cell< MSI::Dof_Type > >& tActiveDofTypes =
                            tProperty->get_dof_type_list();

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
                    const moris::Cell< moris::Cell< PDV_Type > >& tActiveDvTypes =
                            tProperty->get_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                        // if dof enum not in the list
                        if ( tDvCheckList( tDvTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                        }
                    }

                    // get field types for property
                    const moris::Cell< moris::Cell< mtk::Field_Type > >& tActiveFieldTypes =
                            tProperty->get_field_type_list();

                    // loop on property field type
                    for ( uint iFi = 0; iFi < tActiveFieldTypes.size(); iFi++ )
                    {
                        // get set index for field type
                        sint tFieldTypeIndex = mSet->get_index_from_unique_field_type_map( tActiveFieldTypes( iFi )( 0 ) );

                        // if field enum not in the list
                        if ( tFieldCheckList( tFieldTypeIndex ) != 1 )
                        {
                            // put the field type in the check list
                            tFieldCheckList( tFieldTypeIndex ) = 1;

                            // put the field type in the global type list
                            mMasterGlobalFieldTypes.push_back( tActiveFieldTypes( iFi ) );
                        }
                    }
                }
            }

            // get dof type from master constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get dof types for constitutive model
                    const moris::Cell< moris::Cell< MSI::Dof_Type > >& tActiveDofTypes =
                            tCM->get_global_dof_type_list();

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

                    // get dv types for constitutive model
                    const moris::Cell< moris::Cell< PDV_Type > >& tActiveDvTypes =
                            tCM->get_global_dv_type_list();

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

                    // get field types for property
                    const moris::Cell< moris::Cell< mtk::Field_Type > >& tActiveFieldTypes =
                            tCM->get_global_field_type_list();

                    // loop on property field type
                    for ( uint iFi = 0; iFi < tActiveFieldTypes.size(); iFi++ )
                    {
                        // get set index for field type
                        sint tFieldTypeIndex = mSet->get_index_from_unique_field_type_map( tActiveFieldTypes( iFi )( 0 ) );

                        // if field enum not in the list
                        if ( tFieldCheckList( tFieldTypeIndex ) != 1 )
                        {
                            // put the field type in the check list
                            tFieldCheckList( tFieldTypeIndex ) = 1;

                            // put the field type in the global type list
                            mMasterGlobalFieldTypes.push_back( tActiveFieldTypes( iFi ) );
                        }
                    }
                }
            }

            // get dof type from master stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter >& tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get dof types for constitutive model
                    const moris::Cell< moris::Cell< MSI::Dof_Type > >& tActiveDofTypes =
                            tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

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

                    // get dv types for constitutive model
                    const moris::Cell< moris::Cell< PDV_Type > >& tActiveDvTypes =
                            tSP->get_global_dv_type_list( mtk::Master_Slave::MASTER );

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
            mMasterGlobalFieldTypes.shrink_to_fit();

            // SLAVE--------------------------------------------------------

            // set size for the global dof type list
            mSlaveGlobalDofTypes.reserve( tNumDofTypes );
            mSlaveGlobalDvTypes.reserve( tNumDvTypes );
            mSlaveGlobalFieldTypes.reserve( tNumFieldTypes );

            // set a size for the checkList ( used to avoid repeating a dof type)
            tDofCheckList.fill( -1 );
            tDvCheckList.fill( -1 );
            tFieldCheckList.fill( -1 );

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

            // get field type from slave direct dependencies
            for ( uint iFi = 0; iFi < mSlaveFieldTypes.size(); iFi++ )
            {
                // get set index for field type
                sint tFieldTypeIndex = mSet->get_index_from_unique_field_type_map( mSlaveFieldTypes( iFi )( 0 ) );

                // put the field type in the check list
                tFieldCheckList( tFieldTypeIndex ) = 1;

                // put the field type in the global type list
                mSlaveGlobalFieldTypes.push_back( mSlaveFieldTypes( iFi ) );
            }

            // get dof type from master properties
            for ( const std::shared_ptr< Property >& tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get dof types for property
                    const moris::Cell< moris::Cell< MSI::Dof_Type > >& tActiveDofTypes =
                            tProperty->get_dof_type_list();

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
                    const moris::Cell< moris::Cell< PDV_Type > >& tActiveDvTypes =
                            tProperty->get_dv_type_list();

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

                    // get field types for property
                    const moris::Cell< moris::Cell< mtk::Field_Type > >& tActiveFieldTypes =
                            tProperty->get_field_type_list();

                    // loop on property field type
                    for ( uint iFi = 0; iFi < tActiveFieldTypes.size(); iFi++ )
                    {
                        // get set index for field type
                        sint tFieldTypeIndex = mSet->get_index_from_unique_field_type_map( tActiveFieldTypes( iFi )( 0 ) );

                        // if field enum not in the list
                        if ( tFieldCheckList( tFieldTypeIndex ) != 1 )
                        {
                            // put the field type in the check list
                            tFieldCheckList( tFieldTypeIndex ) = 1;

                            // put the field type in the global type list
                            mSlaveGlobalFieldTypes.push_back( tActiveFieldTypes( iFi ) );
                        }
                    }
                }
            }

            // get dof type from slave constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mSlaveCM )
            {
                if ( tCM != nullptr )
                {
                    // get dof types for constitutive model
                    const moris::Cell< moris::Cell< MSI::Dof_Type > >& tActiveDofTypes =
                            tCM->get_global_dof_type_list();

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
                    const moris::Cell< moris::Cell< PDV_Type > >& tActiveDvTypes =
                            tCM->get_global_dv_type_list();

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

                    // get field types for property
                    const moris::Cell< moris::Cell< mtk::Field_Type > >& tActiveFieldTypes =
                            tCM->get_global_field_type_list();

                    // loop on property field type
                    for ( uint iFi = 0; iFi < tActiveFieldTypes.size(); iFi++ )
                    {
                        // get set index for field type
                        sint tFieldTypeIndex = mSet->get_index_from_unique_field_type_map( tActiveFieldTypes( iFi )( 0 ) );

                        // if field enum not in the list
                        if ( tFieldCheckList( tFieldTypeIndex ) != 1 )
                        {
                            // put the field type in the check list
                            tFieldCheckList( tFieldTypeIndex ) = 1;

                            // put the field type in the global type list
                            mSlaveGlobalFieldTypes.push_back( tActiveFieldTypes( iFi ) );
                        }
                    }
                }
            }

            // get dof type from stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter >& tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get dof types for constitutive model
                    const moris::Cell< moris::Cell< MSI::Dof_Type > >& tActiveDofTypes =
                            tSP->get_global_dof_type_list( mtk::Master_Slave::SLAVE );

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
                    const moris::Cell< moris::Cell< PDV_Type > >& tActiveDvTypes =
                            tSP->get_global_dv_type_list( mtk::Master_Slave::SLAVE );

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
            mSlaveGlobalFieldTypes.shrink_to_fit();
        }

        //------------------------------------------------------------------------------

        void
        IQI::build_requested_dof_type_list()
        {
            // clear the dof lists
            mRequestedMasterGlobalDofTypes.clear();
            mRequestedSlaveGlobalDofTypes.clear();

            // get the requested dof types
            Cell< enum MSI::Dof_Type > tRequestedDofTypes = mSet->get_requested_dof_types();

            // reserve possible max size for requested dof lists
            mRequestedMasterGlobalDofTypes.reserve( tRequestedDofTypes.size() );
            mRequestedSlaveGlobalDofTypes.reserve( tRequestedDofTypes.size() );

            // loop over the requested dof types
            for ( auto tDofTypes : tRequestedDofTypes )
            {
                // loop over the IWG master dof types groups
                for ( uint Ik = 0; Ik < mMasterGlobalDofTypes.size(); Ik++ )
                {
                    // if requested dof type matches IWG master dof type
                    if ( mMasterGlobalDofTypes( Ik )( 0 ) == tDofTypes )
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
                    if ( mSlaveGlobalDofTypes( Ik )( 0 ) == tDofTypes )
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

        void
        IQI::select_dQIdu_FD(
                real               aWStar,
                real               aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tQIIndex );

            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumFDPoints = tFDScheme( 0 ).size();

            // get master number of dof types
            uint tNumMasterDofType = mRequestedMasterGlobalDofTypes.size();

            // reset the QI
            mSet->get_QI()( tQIIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tQIIndex );

            // loop over the IWG dof types
            for ( uint iFI = 0; iFI < tNumMasterDofType; iFI++ )
            {
                // init dof counter
                uint tDofCounter = 0;

                // get the dof type
                Cell< MSI::Dof_Type >& tDofType = mRequestedMasterGlobalDofTypes( iFI );

                // get master index for residual dof type, indices for assembly
                uint tMasterDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );

                // get field interpolator for dependency dof type
                Field_Interpolator* tFI =
                        mMasterFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // loop over the coefficient column
                for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficient row
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

                        // if backward or forward fd
                        if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdu
                            mSet->get_residual()( tQIIndex )(
                                    { tMasterDepStartIndex + tDofCounter, tMasterDepStartIndex + tDofCounter },
                                    { 0, 0 } ) +=
                                    tFDScheme( 1 )( 0 ) * tQI /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over the points for FD
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // set the perturbed coefficients to FI
                            tFI->set_coeff( tCoeffPert );
                            tFI->reset_eval_flags();    // not useful

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the QI
                            mSet->get_QI()( tQIIndex ).fill( 0.0 );

                            // compute the QI
                            this->compute_QI( aWStar );

                            // assemble the dQIdu
                            mSet->get_residual()( tQIIndex )(
                                    { tMasterDepStartIndex + tDofCounter, tMasterDepStartIndex + tDofCounter },
                                    { 0, 0 } ) +=
                                    tFDScheme( 1 )( iPoint ) *      //
                                    mSet->get_QI()( tQIIndex ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // get slave number of dof types
            uint tNumSlaveDofType = mRequestedSlaveGlobalDofTypes.size();

            // loop over the slave dof types
            for ( uint iFI = 0; iFI < tNumSlaveDofType; iFI++ )
            {
                // init dof counter
                uint tDofCounter = 0;

                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iFI );

                // get the index for the dof type
                sint tSlaveDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDepDofIndex )( tSlaveDepDofIndex, 0 );

                // get slave dependency field interpolator
                Field_Interpolator* tFI =
                        mSlaveFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // loop over the coefficients columns
                for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficients rows
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

                        // if backward or forward fd
                        if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdu
                            mSet->get_residual()( tQIIndex )(
                                    { tSlaveDepStartIndex + tDofCounter, tSlaveDepStartIndex + tDofCounter },
                                    { 0, 0 } ) +=
                                    tFDScheme( 1 )( 0 ) * tQI /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over the points for FD
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // set the perturbed coefficients to FI
                            tFI->set_coeff( tCoeffPert );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();
                            tFI->reset_eval_flags();    // not useful

                            // reset the QI
                            mSet->get_QI()( tQIIndex ).fill( 0.0 );

                            // compute the QI
                            this->compute_QI( aWStar );

                            // assemble the dQIdu
                            mSet->get_residual()( tQIIndex )(
                                    { tSlaveDepStartIndex + tDofCounter, tSlaveDepStartIndex + tDofCounter },
                                    { 0, 0 } ) +=
                                    tFDScheme( 1 )( iPoint ) *      //
                                    mSet->get_QI()( tQIIndex ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // reset QI value
            mSet->get_QI()( tQIIndex ) = tQIStore;
        }

        //------------------------------------------------------------------------------

        bool
        IQI::check_dQIdu_FD(
                real               aWStar,
                real               aPerturbation,
                real               aEpsilon,
                Matrix< DDRMat >&  adQIdu,
                Matrix< DDRMat >&  adQIduFD,
                bool               aErrorPrint,
                fem::FDScheme_Type aFDSchemeType )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // compute dQIdu with IQI
            this->compute_dQIdu( aWStar );
            adQIdu = mSet->get_residual()( tQIIndex );

            // reset dQIdu
            mSet->get_residual()( tQIIndex ).fill( 0.0 );

            // compute dQIdu by FD
            this->compute_dQIdu_FD( aWStar, aPerturbation, aFDSchemeType );
            adQIduFD = mSet->get_residual()( tQIIndex );

            // define a boolean for check
            bool tCheckdQIdu = true;

            // check if adQIdu and adQIduFD have the same size
            tCheckdQIdu = tCheckdQIdu && ( adQIdu.n_rows() == adQIduFD.n_rows() );
            tCheckdQIdu = tCheckdQIdu && ( adQIdu.n_cols() == adQIduFD.n_cols() );

            // check that matrices to compare have same size
            MORIS_ERROR(
                    ( adQIdu.n_rows() == adQIduFD.n_rows() ) &&    //
                            ( adQIdu.n_cols() == adQIduFD.n_cols() ),
                    "IWG::check_dQIdu - matrices to check do not share same dimensions." );

            // define a real for absolute difference
            real tAbsolute = 0.0;

            // define a real for relative difference
            real tRelative = 0.0;

            // loop over the rows
            for ( uint iRow = 0; iRow < adQIdu.n_rows(); iRow++ )
            {
                // loop over the columns
                for ( uint jCol = 0; jCol < adQIdu.n_cols(); jCol++ )
                {
                    // get absolute difference
                    tAbsolute = std::abs( adQIduFD( iRow, jCol ) - adQIdu( iRow, jCol ) );

                    // get relative difference
                    tRelative = std::abs( ( adQIduFD( iRow, jCol ) - adQIdu( iRow, jCol ) ) / adQIduFD( iRow, jCol ) );

                    // update check value
                    tCheckdQIdu = tCheckdQIdu && ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) );

                    // debug print
                    if ( ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) ) == false )
                    {
                        if ( aErrorPrint )
                        {
                            std::cout << "iRow " << iRow << " - jCol " << jCol << "\n"
                                      << std::flush;
                            std::cout << "adQIdu( iRow, jCol )    " << std::setprecision( 12 ) << adQIdu( iRow, jCol ) << "\n"
                                      << std::flush;
                            std::cout << "adQIduFD( iRow, jCol )  " << std::setprecision( 12 ) << adQIduFD( iRow, jCol ) << "\n"
                                      << std::flush;
                            std::cout << "Absolute difference " << tAbsolute << "\n"
                                      << std::flush;
                            std::cout << "Relative difference " << tRelative << "\n"
                                      << std::flush;
                        }
                    }
                }
            }
            // return bool
            return tCheckdQIdu;
        }

        //------------------------------------------------------------------------------

        void
        IQI::select_dQIdp_FD_material(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // get the column index to assemble in residual
            sint tIQIAssemblyIndex = mSet->get_QI_assembly_index( mName );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tIQIAssemblyIndex );

            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumFDPoints = tFDScheme( 0 ).size();

            // get the requested ip pdv types
            moris::Cell< moris::Cell< PDV_Type > > tRequestedPdvTypes;
            mSet->get_ip_dv_types_for_set( tRequestedPdvTypes );

            // get number of requested dv types
            uint tNumPdvType = tRequestedPdvTypes.size();

            // reset the QI
            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tIQIAssemblyIndex );

            // loop over the pdv types
            for ( uint iFI = 0; iFI < tNumPdvType; iFI++ )
            {
                // get the FI for the dv type
                Field_Interpolator* tFI =
                        mMasterFIManager->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // init pdv coeff counter
                uint tPdvCoeffCounter = 0;

                // loop over the pdv coefficient column
                for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the pdv coefficient row
                    for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // compute the perturbation absolute value
                        real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if ( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // get the pdv index for assembly
                        uint tPdvAssemblyIndex = mSet->get_mat_pdv_assembly_map()( iFI )( 0, 0 ) + tPdvCoeffCounter;

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward fd
                        if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdp
                            mSet->get_dqidpmat()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( 0 ) * tQI( 0 ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over the points for FD
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // set the perturbed coefficients to FI
                            tFI->set_coeff( tCoeffPert );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the QI
                            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

                            // compute the QI
                            this->compute_QI( aWStar );

                            // assemble the jacobian
                            mSet->get_dqidpmat()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( iPoint ) *                    //
                                    mSet->get_QI()( tIQIAssemblyIndex )( 0 ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update the pdv coefficient counter
                        tPdvCoeffCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // reset QI value
            mSet->get_QI()( tIQIAssemblyIndex ) = tQIStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_dqidpmat()( tIQIAssemblyIndex ) ),
                    "IQI::compute_dQIdp_FD_material - dQIdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        real
        IQI::build_perturbation_size(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real  aTolerance )
        {
            return ( this->*m_build_perturbation_size )(
                    aPerturbation,
                    aCoefficientToPerturb,
                    aMaxPerturbation,
                    aTolerance );
        }

        //------------------------------------------------------------------------------

        real
        IQI::build_perturbation_size_relative(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance )
        {
            // check that maximum perturbation size is larger than tolerance
            MORIS_ASSERT( aMaxPerturbation >= aTolerance,
                    "IQI::build_perturbation_size_relative - maximum perturbation size is smaller than tolerance: max = %e  tol = %e\n",
                    aMaxPerturbation,
                    aTolerance );

            // determine actual tolerance (only useful when above assert inactive)
            const real tActualTol = std::min( aMaxPerturbation, aTolerance );

            // compute the perturbation value using fraction of maximum allowable perturbation
            real tDeltaH = std::abs( aPerturbation * aMaxPerturbation );

            // compute perturbation such that it is not smaller than tolerance
            // and not larger than maximum value
            tDeltaH = std::max( std::min( tDeltaH, aMaxPerturbation ), tActualTol );

            return tDeltaH;
        }

        //------------------------------------------------------------------------------

        real
        IQI::build_perturbation_size_absolute(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance )
        {
            // check that maximum perturbation size is larger than tolerance
            MORIS_ASSERT( aMaxPerturbation >= aTolerance,
                    "IQI::build_perturbation_size_absolute - maximum perturbation size is smaller than tolerance: max = %e  tol = %e\n",
                    aMaxPerturbation,
                    aTolerance );

            // determine actual tolerance (only useful when above assert inactive)
            const real tActualTol = std::min( aMaxPerturbation, aTolerance );

            // check that absolute value of perturbation is not smaller than tolerance
            // and not larger than maximum value
            return std::max( std::min( std::abs( aPerturbation ), aMaxPerturbation ), tActualTol );
        }

        //------------------------------------------------------------------------------

        real
        IQI::check_ig_coordinates_inside_ip_element(
                const real&         aPerturbation,
                const real&         aCoefficientToPerturb,
                const uint&         aSpatialDirection,
                fem::FDScheme_Type& aUsedFDScheme )
        {
            // FIXME: only works for rectangular IP elements
            // FIXME: only works for forward, backward, central, not for higher as 5-point FD

            // get the IP element geometry interpolator
            Geometry_Interpolator* tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();

            // IP element max/min
            real tMaxIP = max( tIPGI->get_space_coeff().get_column( aSpatialDirection ) );    // get maximum values of coordinates of IP nodes
            real tMinIP = min( tIPGI->get_space_coeff().get_column( aSpatialDirection ) );    // get minimum values of coordinates of IP nodes

            // get maximum possible perturbation
            real tMaxPerturb = ( tMaxIP - tMinIP ) / 3.0;

            // compute the perturbation value
            real tDeltaH = build_perturbation_size( aPerturbation, aCoefficientToPerturb, tMaxPerturb );

            // check that IG node coordinate is consistent with minimum and maximum IP coordinates
            MORIS_ASSERT(
                    tMaxIP >= aCoefficientToPerturb - tDeltaH &&    //
                            tMinIP <= aCoefficientToPerturb + tDeltaH,
                    "ERROR: IG coordinates are outside IP element: dim: %d  minIP: %e  maxIP: %e  cordIG: %e  \n",
                    aSpatialDirection,
                    tMinIP,
                    tMaxIP,
                    aCoefficientToPerturb );

            // check point location
            if ( aCoefficientToPerturb + tDeltaH >= tMaxIP )
            {
                aUsedFDScheme = fem::FDScheme_Type::POINT_1_BACKWARD;

                // check for correctness of perturbation size for backward FD
                MORIS_ASSERT( tDeltaH < aCoefficientToPerturb - tMinIP,
                        "ERROR: backward perturbation size exceed limits of interpolation element:\n",
                        "dim: %d  minIP: %e  maxIP: %e  cordIG: %e  maxPert: %e  delta: %e  precPert: %e\n.",
                        aSpatialDirection,
                        tMinIP,
                        tMaxIP,
                        aCoefficientToPerturb,
                        tMaxPerturb,
                        tDeltaH,
                        aPerturbation );
            }
            else
            {
                if ( aCoefficientToPerturb - tDeltaH <= tMinIP )
                {
                    aUsedFDScheme = fem::FDScheme_Type::POINT_1_FORWARD;

                    // check for correctness of perturbation size for backward FD
                    MORIS_ASSERT( tDeltaH < tMaxIP - aCoefficientToPerturb,
                            "ERROR: forward perturbation size exceeds limits of interpolation element: dim: %d  minIP: %e  maxIP: %e  cordIG: %e  maxPert: %e  delta: %e  precPert: %e\n.",
                            aSpatialDirection,
                            tMinIP,
                            tMaxIP,
                            aCoefficientToPerturb,
                            tMaxPerturb,
                            tDeltaH,
                            aPerturbation );
                }
                else
                {
                    // check for correctness of perturbation size for central FD
                    MORIS_ASSERT(
                            tDeltaH < tMaxIP - aCoefficientToPerturb &&    //
                                    tDeltaH < aCoefficientToPerturb - tMinIP,
                            "ERROR: central perturbation size exceed limits of interpolation element: dim: %d  minIP: %e  maxIP: %e  cordIG: %e  maxPert: %e  delta: %e  precPert: %e\n.",
                            aSpatialDirection,
                            tMinIP,
                            tMaxIP,
                            aCoefficientToPerturb,
                            tMaxPerturb,
                            tDeltaH,
                            aPerturbation );
                }
            }

            return tDeltaH;
        }

        //------------------------------------------------------------------------------

        void
        IQI::select_dQIdp_FD_geometry_bulk(
                moris::real                        aWStar,
                moris::real                        aPerturbation,
                fem::FDScheme_Type                 aFDSchemeType,
                Matrix< DDSMat >&                  aGeoLocalAssembly,
                moris::Cell< Matrix< IndexMat > >& aVertexIndices )
        {
            // get the IQI index
            uint tIQIAssemblyIndex = mSet->get_QI_assembly_index( mName );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tIQIAssemblyIndex );

            // get the GI for the IP and IG element considered
            Geometry_Interpolator* tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // get number of master GI bases and space dimensions
            uint tDerNumBases      = tIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

            // coefficients for dv type wrt which derivative is computed
            Matrix< DDRMat > tCoeff      = tIGGI->get_space_coeff();
            Matrix< DDRMat > tParamCoeff = tIGGI->get_space_param_coeff();
            Matrix< DDRMat > tEvaluationPoint;
            tIGGI->get_space_time( tEvaluationPoint );
            real tGPWeight = aWStar / tIGGI->det_J();

            // reset the QI
            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tIQIAssemblyIndex );

            // init FD scheme
            moris::Cell< moris::Cell< real > > tFDScheme;

            // init perturbation
            real tDeltaH = 0.0;

            // loop over the spatial directions/loop on pdv type
            for ( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
            {
                // loop over the IG nodes/loop one nodes
                for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    // if pdv is active
                    if ( tPdvAssemblyIndex != -1 )
                    {
                        // check point location and define perturbation size and FD scheme accordingly
                        fem::FDScheme_Type tUsedFDScheme = aFDSchemeType;
                        tDeltaH                          = this->check_ig_coordinates_inside_ip_element(
                                aPerturbation,
                                tCoeff( iCoeffRow, iCoeffCol ),
                                iCoeffCol,
                                tUsedFDScheme );

                        // finalize FD scheme
                        fd_scheme( tUsedFDScheme, tFDScheme );
                        uint tNumFDPoints = tFDScheme( 0 ).size();

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward fd
                        if ( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdp
                            mSet->get_dqidpgeo()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( 0 ) * tQI( 0 ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over point of FD scheme
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // setting the perturbed coefficients
                            tIGGI->set_space_coeff( tCoeffPert );

                            // update local coordinates
                            Matrix< DDRMat > tXCoords  = tCoeffPert.get_row( iCoeffRow );
                            Matrix< DDRMat > tXiCoords = tParamCoeff.get_row( iCoeffRow );
                            tIPGI->update_local_coordinates( tXCoords, tXiCoords );
                            Matrix< DDRMat > tParamCoeffPert     = tParamCoeff;
                            tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();
                            tIGGI->set_space_param_coeff( tParamCoeffPert );

                            // set evaluation point for interpolators (FIs and GIs)
                            mSet->get_field_interpolator_manager()->    //
                                    set_space_time_from_local_IG_point( tEvaluationPoint );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the QI
                            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

                            // compute the QI
                            real tWStarPert = tGPWeight * tIGGI->det_J();
                            this->compute_QI( tWStarPert );

                            // evaluate dQIdpGeo
                            mSet->get_dqidpgeo()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( iPoint ) *                    //
                                    mSet->get_QI()( tIQIAssemblyIndex )( 0 ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                }
                // reset the coefficients values
                tIGGI->set_space_coeff( tCoeff );
                tIGGI->set_space_param_coeff( tParamCoeff );
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );
            }

            // reset QI value
            mSet->get_QI()( tIQIAssemblyIndex ) = tQIStore;

            // if active cluster measure on IQI
            if ( mActiveCMEAFlag )
            {
                // add their contribution to dQIdp
                this->add_cluster_measure_dQIdp_FD_geometry(
                        aWStar,
                        aPerturbation,
                        aFDSchemeType );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_dqidpgeo()( tIQIAssemblyIndex ) ),
                    "IQI::compute_dQIdp_FD_geometry - dQIdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IQI::select_dQIdp_FD_geometry_sideset(
                moris::real                        aWStar,
                moris::real                        aPerturbation,
                fem::FDScheme_Type                 aFDSchemeType,
                Matrix< DDSMat >&                  aGeoLocalAssembly,
                moris::Cell< Matrix< IndexMat > >& aVertexIndices )
        {
            // get the IQI index
            uint tIQIAssemblyIndex = mSet->get_QI_assembly_index( mName );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tIQIAssemblyIndex );

            // get the GI for the IP and IG element considered
            Geometry_Interpolator* tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // store unperturbed xyz coordinates
            Matrix< DDRMat > tCoeff = tIGGI->get_space_coeff();

            // store unperturbed local coordinates
            Matrix< DDRMat > tParamCoeff = tIGGI->get_space_param_coeff();

            // store unperturbed evaluation point
            Matrix< DDRMat > tEvaluationPoint;
            tIGGI->get_space_time( tEvaluationPoint );

            // store unperturbed evaluation point weight
            real tGPWeight = aWStar / tIGGI->det_J();

            // store unperturbed normal
            Matrix< DDRMat > tNormal;
            tIGGI->get_normal( tNormal );

            // reset the QI
            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tIQIAssemblyIndex );

            // init FD scheme
            moris::Cell< moris::Cell< real > > tFDScheme;

            // init perturbation
            real tDeltaH = 0.0;

            // get number of master GI bases and space dimensions
            uint tDerNumBases      = tIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

            // loop over the spatial directions/loop on pdv type
            for ( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
            {
                // loop over the IG nodes/loop one nodes
                for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    // if pdv is active
                    if ( tPdvAssemblyIndex != -1 )
                    {
                        // check point location and define perturbation size and FD scheme accordingly
                        fem::FDScheme_Type tUsedFDScheme = aFDSchemeType;
                        tDeltaH                          = this->check_ig_coordinates_inside_ip_element(
                                aPerturbation,
                                tCoeff( iCoeffRow, iCoeffCol ),
                                iCoeffCol,
                                tUsedFDScheme );

                        // finalize FD scheme
                        fd_scheme( tUsedFDScheme, tFDScheme );
                        uint tNumFDPoints = tFDScheme( 0 ).size();

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward fd
                        if ( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdp
                            mSet->get_dqidpgeo()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( 0 ) * tQI( 0 ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over point of FD scheme
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // setting the perturbed coefficients
                            tIGGI->set_space_coeff( tCoeffPert );

                            // update local coordinates
                            Matrix< DDRMat > tXCoords  = tCoeffPert.get_row( iCoeffRow );
                            Matrix< DDRMat > tXiCoords = tParamCoeff.get_row( iCoeffRow );
                            tIPGI->update_local_coordinates( tXCoords, tXiCoords );
                            Matrix< DDRMat > tParamCoeffPert     = tParamCoeff;
                            tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();
                            tIGGI->set_space_param_coeff( tParamCoeffPert );

                            // set evaluation point for interpolators (FIs and GIs)
                            mSet->get_field_interpolator_manager()->    //
                                    set_space_time_from_local_IG_point( tEvaluationPoint );

                            // reset the normal
                            Matrix< DDRMat > tNormalPert;
                            tIGGI->get_normal( tNormalPert );
                            this->set_normal( tNormalPert );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the QI
                            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

                            // compute the QI
                            real tWStarPert = tGPWeight * tIGGI->det_J();
                            this->compute_QI( tWStarPert );

                            // evaluate dQIdpGeo
                            mSet->get_dqidpgeo()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_QI()( tIQIAssemblyIndex )( 0 ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                }
                // reset xyz values
                tIGGI->set_space_coeff( tCoeff );

                // reset local coordinates
                tIGGI->set_space_param_coeff( tParamCoeff );

                // reset evaluation point
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );

                // reset normal
                this->set_normal( tNormal );
            }

            // reset QI value
            mSet->get_QI()( tIQIAssemblyIndex ) = tQIStore;

            // if active cluster measure on IQI
            if ( mActiveCMEAFlag )
            {
                // add their contribution to dQIdp
                this->add_cluster_measure_dQIdp_FD_geometry(
                        aWStar,
                        aPerturbation,
                        aFDSchemeType );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_dqidpgeo()( tIQIAssemblyIndex ) ),
                    "IQI::compute_dQIdp_FD_geometry - dQIdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IQI::add_cluster_measure_dQIdp_FD_geometry(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // get the IQI index
            uint tIQIAssemblyIndex = mSet->get_QI_assembly_index( mName );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tIQIAssemblyIndex );

            // reset the QI
            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tIQIAssemblyIndex );

            // init FD scheme
            moris::Cell< moris::Cell< real > > tFDScheme;

            // FIXME init perturbation
            real tDeltaH = aPerturbation;

            // loop over the cluster measures
            for ( uint iCMEA = 0; iCMEA < mCluster->get_cluster_measures().size(); iCMEA++ )
            {
                // get treated cluster measure
                std::shared_ptr< Cluster_Measure >& tClusterMeasure =
                        mCluster->get_cluster_measures()( iCMEA );

                // create FD scheme
                fd_scheme( aFDSchemeType, tFDScheme );
                uint tNumFDPoints = tFDScheme( 0 ).size();

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward fd
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                        ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed QI contribution to dQIdp
                    mSet->get_dqidpgeo()( tIQIAssemblyIndex ) +=
                            tFDScheme( 1 )( 0 ) * tQI( 0 ) * tClusterMeasure->dMEAdPDV() /    //
                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over point of FD scheme
                for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                {
                    // perturb the cluster measure
                    tClusterMeasure->perturb_cluster_measure( tFDScheme( 0 )( iPoint ) * tDeltaH );

                    // reset properties, CM and SP for IWG
                    this->reset_eval_flags();

                    // reset the QI
                    mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

                    // compute the QI
                    this->compute_QI( aWStar );

                    // evaluate dQIdpGeo
                    mSet->get_dqidpgeo()( tIQIAssemblyIndex ) +=
                            tFDScheme( 1 )( iPoint ) *                                                  //
                            mSet->get_QI()( tIQIAssemblyIndex )( 0 ) * tClusterMeasure->dMEAdPDV() /    //
                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // reset cluster measures
                    mCluster->reset_cluster_measure();
                }
            }

            // reset QI value
            mSet->get_QI()( tIQIAssemblyIndex ) = tQIStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_dqidpgeo()( tIQIAssemblyIndex ) ),
                    "IQI::compute_dQIdp_FD_geometry - dQIdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------
    }    // namespace fem
}    // namespace moris

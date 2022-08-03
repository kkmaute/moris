/*
 * cl_FEM_IWG.cpp
 *
 *  Created on: Nov 12, 2019
 *      Author: sonne
 */

#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Cluster_Measure.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Model.hpp"

#include "fn_max.hpp"
#include "fn_min.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        void
        IWG::print_names()
        {
            std::cout << "----------" << std::endl;
            std::cout << "IWG: " << mName << std::endl;

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
                    std::cout << "Master CM:       " << mMasterCM( iCM )->get_name() << std::endl;
                }
            }
            for ( uint iCM = 0; iCM < mSlaveCM.size(); iCM++ )
            {
                if ( mSlaveCM( iCM ) != nullptr )
                {
                    std::cout << "Slave CM:        " << mSlaveCM( iCM )->get_name() << std::endl;
                }
            }

            // SP
            for ( uint iSP = 0; iSP < mStabilizationParam.size(); iSP++ )
            {
                if ( mStabilizationParam( iSP ) != nullptr )
                {
                    std::cout << "SP:              " << mStabilizationParam( iSP )->get_name() << std::endl;
                }
            }
            std::cout << "----------" << std::endl;
        }

        //------------------------------------------------------------------------------

        void
        IWG::reset_eval_flags()
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

            // reset material models
            for ( const std::shared_ptr< Material_Model >& tMM : mMasterMM )
            {
                if ( tMM != nullptr )
                {
                    tMM->reset_eval_flags();
                }
            }
            for ( const std::shared_ptr< Material_Model >& tMM : mSlaveMM )
            {
                if ( tMM != nullptr )
                {
                    tMM->reset_eval_flags();
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

            // reset evaluation flags specific to child implementations
            this->reset_spec_eval_flags();
        }

        //------------------------------------------------------------------------------

        void
        IWG::set_function_pointers()
        {
            // switch on element type
            switch ( mSet->get_element_type() )
            {
                case fem::Element_Type::BULK:
                {
                    m_compute_jacobian_FD      = &IWG::select_jacobian_FD;
                    m_compute_dRdp_FD_material = &IWG::select_dRdp_FD_material;
                    m_compute_dRdp_FD_geometry = &IWG::select_dRdp_FD_geometry_bulk;
                    break;
                }
                case fem::Element_Type::SIDESET:
                {
                    m_compute_jacobian_FD      = &IWG::select_jacobian_FD;
                    m_compute_dRdp_FD_material = &IWG::select_dRdp_FD_material;
                    m_compute_dRdp_FD_geometry = &IWG::select_dRdp_FD_geometry_sideset;
                    break;
                }
                case fem::Element_Type::TIME_SIDESET:
                case fem::Element_Type::TIME_BOUNDARY:
                {
                    m_compute_jacobian_FD      = &IWG::select_jacobian_FD;
                    m_compute_dRdp_FD_material = &IWG::select_dRdp_FD_material;
                    m_compute_dRdp_FD_geometry = &IWG::select_dRdp_FD_geometry_time_sideset;
                    break;
                }
                case fem::Element_Type::DOUBLE_SIDESET:
                {
                    m_compute_jacobian_FD      = &IWG::select_jacobian_FD_double;
                    m_compute_dRdp_FD_material = &IWG::select_dRdp_FD_material_double;
                    m_compute_dRdp_FD_geometry = &IWG::select_dRdp_FD_geometry_double;
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
                    m_build_perturbation_size = &IWG::build_perturbation_size_relative;
                    break;
                }
                case fem::Perturbation_Type::ABSOLUTE:
                {
                    m_build_perturbation_size = &IWG::build_perturbation_size_absolute;
                    break;
                }
                default:
                    MORIS_ERROR( false, "IWG::set_function_pointers - unknown perturbation type." );
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG::set_phase_name(
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
        IWG::get_phase_name( mtk::Master_Slave aIsMaster )
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
        IWG::set_field_interpolator_manager(
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

            // loop over the material models
            for ( const std::shared_ptr< Material_Model >& tMM : this->get_material_models( aIsMaster ) )
            {
                if ( tMM != nullptr )
                {
                    // set the field interpolator manager for the CM
                    tMM->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsMaster ) );

                    // set the fem set pointe for the CM
                    tMM->set_set_pointer( mSet );
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
        IWG::get_field_interpolator_manager(
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                    return mMasterFIManager;

                case mtk::Master_Slave::SLAVE:
                    return mSlaveFIManager;

                default:
                    MORIS_ERROR( false, "IWG::get_field_interpolator_manager - can only be master or slave." );
                    return mMasterFIManager;
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG::set_field_interpolator_manager_previous_time(
                Field_Interpolator_Manager* aFieldInterpolatorManager,
                mtk::Master_Slave           aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    mMasterPreviousFIManager = aFieldInterpolatorManager;
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "IWG::set_field_interpolator_manager - can only be master" );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG::set_normal( Matrix< DDRMat >& aNormal )
        {
            mNormal = aNormal;

            // set normal for SP
            for ( const std::shared_ptr< Stabilization_Parameter >& tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    tSP->set_normal( mNormal );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG::set_interpolation_order()
        {
            // if order is already set
            if ( ( mOrder != MORIS_UINT_MAX ) ) return;

            // get field interpolator manager
            Field_Interpolator_Manager* tFIManager = mSet->get_field_interpolator_manager();

            // define axuxilliary variable for storing interpolation order
            uint tPrevOrder = MORIS_UINT_MAX;

            // loop overall residual dof types
            for ( uint iType = 0; iType < mResidualDofType.size(); ++iType )
            {
                // get field interpolator for current residual dof type;
                Field_Interpolator* tFI = tFIManager->get_field_interpolators_for_type( mResidualDofType( iType )( 0 ) );

                // get residual dof type interpolation order
                mtk::Interpolation_Order tInterpOrder = tFI->get_space_interpolation_order();

                // set the interpolation order for IWG
                switch ( tInterpOrder )
                {
                    case mtk::Interpolation_Order::LINEAR:
                    {
                        mOrder = 1;
                        break;
                    }
                    case mtk::Interpolation_Order::QUADRATIC:
                    {
                        mOrder = 2;
                        break;
                    }
                    case mtk::Interpolation_Order::CUBIC:
                    {
                        mOrder = 3;
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( false, "IWG::set_interpolation_order - order not supported" );
                    }
                }

                // check order is the same for all residual types
                MORIS_ERROR( iType > 0 ? tPrevOrder == mOrder : true,
                        "IWG::set_interpolation_order - Interpolation order of residual fields need to be identical among all residual fields.\n" );

                tPrevOrder = mOrder;
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG::set_dof_type_list(
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
                    MORIS_ERROR( false, "IWG::set_dof_type_list - can only be MASTER or SLAVE." );
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< MSI::Dof_Type > >&
        IWG::get_dof_type_list(
                mtk::Master_Slave aIsMaster ) const
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
                    MORIS_ASSERT( false, "IWG::get_dof_type_list - can only be master or slave." );
                    return mMasterDofTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG::set_dv_type_list(
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
                    MORIS_ERROR( false, "IWG::set_dv_type_list - can only be MASTER or SLAVE." );
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< PDV_Type > >&
        IWG::get_dv_type_list(
                mtk::Master_Slave aIsMaster ) const
        {
            // switch on master/slave
            switch ( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER:
                {
                    // return master global dof type list
                    return mMasterDvTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE:
                {
                    // return slave global dof type list
                    return mSlaveDvTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IWG::get_dv_type_list - can only be master or slave." );
                    return mMasterDvTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG::set_field_type_list(
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
        IWG::get_field_type_list(
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
        IWG::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "IWG::set_property - IWG %s - Unknown aPropertyString: %s ",
                    mName.c_str(),
                    aPropertyString.c_str() );

            // set the property in the property pointer cell
            this->get_properties( aIsMaster )( mPropertyMap[ aPropertyString ] ) = aProperty;
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< Property > >&
        IWG::get_properties(
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
                    MORIS_ASSERT( false, "IWG::get_properties - can only be master or slave." );
                    return mMasterProp;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG::set_material_model(
                std::shared_ptr< Material_Model > aMaterialModel,
                std::string                       aMaterialModelString,
                mtk::Master_Slave                 aIsMaster )
        {
            // check that aConstitutiveString makes sense
            MORIS_ERROR( mMaterialMap.find( aMaterialModelString ) != mMaterialMap.end(),
                    "IWG::set_material_model - IWG %s - Unknown aMaterialModelString: %s ",
                    mName.c_str(),
                    aMaterialModelString.c_str() );

            // set the CM in the CM pointer cell
            this->get_material_models( aIsMaster )( mMaterialMap[ aMaterialModelString ] ) = aMaterialModel;
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< Material_Model > >&
        IWG::get_material_models(
                mtk::Master_Slave aIsMaster )
        {
            // switch on master/slave
            switch ( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER:
                {
                    // return master property pointers
                    return mMasterMM;
                }
                // if slave
                case mtk::Master_Slave::SLAVE:
                {
                    // return slave property pointers
                    return mSlaveMM;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IWG::get_material_models - can only be master or slave." );
                    return mMasterMM;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster )
        {
            // check that aConstitutiveString makes sense
            MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(),
                    "IWG::set_constitutive_model - IWG %s - Unknown aConstitutiveString: %s ",
                    mName.c_str(),
                    aConstitutiveString.c_str() );

            // set the CM in the CM pointer cell
            this->get_constitutive_models( aIsMaster )( mConstitutiveMap[ aConstitutiveString ] ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::shared_ptr< Constitutive_Model > >&
        IWG::get_constitutive_models(
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
                    MORIS_ASSERT( false, "IWG::get_constitutive_models - can only be master or slave." );
                    return mMasterCM;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG::set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString )
        {
            // check that aStabilizationString makes sense
            MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(),
                    "IWG::set_stabilization_parameter - IWG %s - Unknown aStabilizationString: %s ",
                    mName.c_str(),
                    aStabilizationString.c_str() );

            // set the stabilization parameter in the stabilization parameter cell
            this->get_stabilization_parameters()( mStabilizationMap[ aStabilizationString ] ) = aStabilizationParameter;

            // set active cluster measure on IWG flag on/off
            mActiveCMEAFlag = mActiveCMEAFlag || ( aStabilizationParameter->get_cluster_measure_tuple_list().size() > 0 );
        }

        //------------------------------------------------------------------------------

        void
        IWG::get_non_unique_dof_dv_and_field_types(
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
                    tMasterFieldCounter += tActiveDvTypes.size();
                }
            }

            // loop over slave properties
            for ( const std::shared_ptr< Property >& tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type lists
                    // get property non unique dof and dv type list
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

            // loop over master material models
            for ( const std::shared_ptr< Material_Model >& tMM : mMasterMM )
            {
                if ( tMM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< PDV_Type >      tActiveDvTypes;

                    tMM->get_non_unique_dof_types( tActiveDofTypes );

                    // update dof counters (DVs not part of MM yet)
                    tMasterDofCounter += tActiveDofTypes.size();
                }
            }

            // loop over slave material models
            for ( const std::shared_ptr< Material_Model >& tMM : mSlaveMM )
            {
                if ( tMM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;

                    tMM->get_non_unique_dof_types( tActiveDofTypes );

                    // update dof and dv counters
                    tSlaveDofCounter += tActiveDofTypes.size();
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
                // populate the field list
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

            // loop over the master material models
            for ( const std::shared_ptr< Material_Model >& tMM : mMasterMM )
            {
                if ( tMM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;

                    tMM->get_non_unique_dof_types( tActiveDofTypes );

                    // update dof counters (DVs not part of MM yet)
                    tMasterDofCounter += tActiveDofTypes.size();
                }
            }

            // loop over the slave material models
            for ( const std::shared_ptr< Material_Model >& tMM : mSlaveMM )
            {
                if ( tMM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;

                    tMM->get_non_unique_dof_types( tActiveDofTypes );

                    // update dof counters (DVs not part of MM yet)
                    tMasterDofCounter += tActiveDofTypes.size();
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
                    aFieldTypes( 0 ).append( tActiveFieldTypes );
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
        IWG::build_global_dof_dv_and_field_type_list()
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
                            // put the dof type in the check list
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
                            // put the dof type in the check list
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

            // get dof type from master material models
            for ( const std::shared_ptr< Material_Model >& tMM : mMasterMM )
            {
                if ( tMM != nullptr )
                {
                    // get dof types for material modIWG::build_global_dof_dv_and_field_listel
                    const moris::Cell< moris::Cell< MSI::Dof_Type > >& tActiveDofTypes =
                            tMM->get_global_dof_type_list();

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the check list
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }
                    // skip loop on material model dv type - not implemented
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
                            // put the dof type in the check list
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
                            // put the dv type in the check list
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dv type in the global type list
                            mMasterGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                        }
                    }

                    // get field types for constitutive model
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
                            // put the dof type in the check list
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
                            // put the dv type in the check list
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

                // put the dof type in the check list
                tDofCheckList( tDofTypeIndex ) = 1;

                // put the dof type in the global type list
                mSlaveGlobalDofTypes.push_back( mSlaveDofTypes( iDof ) );
            }

            // get dv type from slave direct dependencies
            for ( uint iDv = 0; iDv < mSlaveDvTypes.size(); iDv++ )
            {
                // get set index for dv type
                sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( mSlaveDvTypes( iDv )( 0 ) );

                // put the dv type in the check list
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

            // get dof type from slave properties
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
                            // put the dof type in the check list
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
                            // put the dv type in the check list
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

            // get dof type from slave material models
            for ( std::shared_ptr< Material_Model > tMM : mSlaveMM )
            {
                if ( tMM != nullptr )
                {
                    // get dof types for constitutive model
                    const moris::Cell< moris::Cell< MSI::Dof_Type > >& tActiveDofTypes =
                            tMM->get_global_dof_type_list();

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the check list
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mSlaveGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }
                    // skip loop on material dv type - not implemented
                }
            }

            // get dof type from slave constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
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
                            // put the dof type in the check list
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
                            // put the dv type in the check list
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dv type in the global type list
                            mSlaveGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                        }
                    }

                    // get field types for constitutive model
                    const moris::Cell< moris::Cell< mtk::Field_Type > >& tActiveFieldTypes =
                            tCM->get_field_type_list();

                    // loop on constitutive model field type
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
            for ( std::shared_ptr< Stabilization_Parameter > tSP : mStabilizationParam )
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
                            // put the dof type in the check list
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
                            // put the dv type in the check list
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
        IWG::build_requested_dof_type_list( const bool aIsStaggered )
        {
            // clear the dof lists
            mRequestedMasterGlobalDofTypes.clear();
            mRequestedSlaveGlobalDofTypes.clear();

            moris::Cell< enum MSI::Dof_Type > tRequestedDofTypes;

            // if residual evaluation
            if ( aIsStaggered )
            {
                // get the requested dof types
                tRequestedDofTypes = mSet->get_secondary_dof_types();
            }
            // if Jacobian evaluation
            else
            {
                // get the requested dof types
                tRequestedDofTypes = mSet->get_requested_dof_types();
            }

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
        IWG::check_field_interpolators( mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    // loop over the dof field interpolator pointers
                    for ( uint iDofFI = 0; iDofFI < mRequestedMasterGlobalDofTypes.size(); iDofFI++ )
                    {
                        // check that the field interpolator was set
                        MORIS_ASSERT(
                                this->get_field_interpolator_manager( aIsMaster )->get_field_interpolators_for_type( mRequestedMasterGlobalDofTypes( iDofFI )( 0 ) ) != nullptr,
                                "IWG::check_field_interpolators - Master dof FI missing. " );
                    }

                    // loop over the dv field interpolator pointers
                    for ( uint iDvFI = 0; iDvFI < mMasterGlobalDvTypes.size(); iDvFI++ )
                    {
                        // check that the field interpolator was set
                        MORIS_ASSERT(
                                this->get_field_interpolator_manager( aIsMaster )->get_field_interpolators_for_type( mMasterGlobalDvTypes( iDvFI )( 0 ) ) != nullptr,
                                "IWG::check_field_interpolators - Master dv FI missing. " );
                    }
                    break;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    // loop over the dof field interpolator pointers
                    for ( uint iDofFI = 0; iDofFI < mRequestedSlaveGlobalDofTypes.size(); iDofFI++ )
                    {
                        // check that the field interpolator was set
                        MORIS_ASSERT(
                                this->get_field_interpolator_manager( aIsMaster )->get_field_interpolators_for_type( mRequestedSlaveGlobalDofTypes( iDofFI )( 0 ) ) != nullptr,
                                "IWG::check_dof_field_interpolators - Slave dof FI missing. " );
                    }

                    // loop over the dv field interpolator pointers
                    for ( uint iDvFI = 0; iDvFI < mSlaveGlobalDvTypes.size(); iDvFI++ )
                    {
                        // check that the field interpolator was set
                        MORIS_ASSERT(
                                this->get_field_interpolator_manager( aIsMaster )->get_field_interpolators_for_type( mSlaveGlobalDvTypes( iDvFI )( 0 ) ) != nullptr,
                                "IWG::check_field_interpolators - Slave dv FI missing. " );
                    }
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IWG::check_field_interpolators - can only be master or slave." );
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< MSI::Dof_Type > >&
        IWG::get_global_dof_type_list(
                mtk::Master_Slave aIsMaster )
        {
            // if the global list was not yet built
            if ( mGlobalDofBuild )
            {
                // build the stabilization parameter global dof type list
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
                    break;
                }
                // if slave
                case mtk::Master_Slave::SLAVE:
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
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< PDV_Type > >&
        IWG::get_global_dv_type_list(
                mtk::Master_Slave aIsMaster )
        {
            // if the global list was not yet built
            if ( mGlobalDvBuild )
            {
                // build the stabilization parameter global dof type list
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
                    return mMasterGlobalDvTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE:
                {
                    // return slave global dof type list
                    return mSlaveGlobalDvTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IWG::get_global_dv_type_list - can only be master or slave." );
                    return mMasterGlobalDvTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< mtk::Field_Type > >&
        IWG::get_global_field_type_list(
                mtk::Master_Slave aIsMaster )
        {
            // if the global list was not yet built
            if ( mGlobalFieldBuild )
            {
                // build the stabilization parameter global dof type list
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
                    // return master global field type list
                    return mMasterGlobalFieldTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE:
                {
                    // return slave global field type list
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
        IWG::select_jacobian_FD(
                real               aWStar,
                real               aPerturbation,
                fem::FDScheme_Type aFDSchemeType,
                bool               aUseAbsolutePerturbations )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumFDPoints = tFDScheme( 0 ).size();

            // get master index for residual dof type, indices for assembly
            sint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get master number of dof types
            uint tMasterNumDofTypes = mRequestedMasterGlobalDofTypes.size();

            // reset and evaluate the residual plus
            mSet->get_residual()( 0 ).fill( 0.0 );
            this->compute_residual( aWStar );
            Matrix< DDRMat > tResidual =
                    mSet->get_residual()( 0 )(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { 0, 0 } );

            // loop over the IWG dof types
            for ( uint iFI = 0; iFI < tMasterNumDofTypes; iFI++ )
            {
                // init dof counter
                uint tDofCounter = 0;

                // get the dof type
                Cell< MSI::Dof_Type >& tDofType = mRequestedMasterGlobalDofTypes( iFI );

                // get the index for the dof type
                sint tMasterDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tMasterDepDofIndex, 0 );

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
                        real tDeltaH = aPerturbation;
                        if ( !aUseAbsolutePerturbations )
                        {
                            tDeltaH = tDeltaH * tCoeff( iCoeffRow, iCoeffCol );
                        }

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
                            // add unperturbed residual contribution to dRdp
                            mSet->get_jacobian()(
                                    { tMasterResStartIndex, tMasterResStopIndex },
                                    { tMasterDepStartIndex + tDofCounter, tMasterDepStartIndex + tDofCounter } ) +=
                                    tFDScheme( 1 )( 0 ) * tResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

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

                            // reset the residual
                            mSet->get_residual()( 0 ).fill( 0.0 );

                            // compute the residual
                            this->compute_residual( aWStar );

                            // assemble the Jacobian
                            mSet->get_jacobian()(
                                    { tMasterResStartIndex, tMasterResStopIndex },
                                    { tMasterDepStartIndex + tDofCounter, tMasterDepStartIndex + tDofCounter } ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;
        }

        //------------------------------------------------------------------------------

        void
        IWG::select_jacobian_FD_double(
                real               aWStar,
                real               aPerturbation,
                fem::FDScheme_Type aFDSchemeType,
                bool               aUseAbsolutePerturbations )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumFDPoints = tFDScheme( 0 ).size();

            // get master index for residual dof type, indices for assembly
            sint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            sint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // reset and evaluate the residual plus
            mSet->get_residual()( 0 ).fill( 0.0 );
            this->compute_residual( aWStar );
            Matrix< DDRMat > tMasterResidual =
                    mSet->get_residual()( 0 )(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { 0, 0 } );
            Matrix< DDRMat > tSlaveResidual =
                    mSet->get_residual()( 0 )(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { 0, 0 } );

            // get master number of dof types
            uint tMasterNumDofTypes = mRequestedMasterGlobalDofTypes.size();

            // loop over the IWG dof types
            for ( uint iFI = 0; iFI < tMasterNumDofTypes; iFI++ )
            {
                // init dof counter
                uint tDofCounter = 0;

                // get the dof type
                Cell< MSI::Dof_Type >& tDofType = mRequestedMasterGlobalDofTypes( iFI );

                // get the index for the dof type
                sint tMasterDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tMasterDepDofIndex, 0 );

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
                        real tDeltaH = aPerturbation;
                        if ( !aUseAbsolutePerturbations )
                        {
                            tDeltaH = tDeltaH * tCoeff( iCoeffRow, iCoeffCol );
                        }

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
                            // add unperturbed master residual contribution to dRdp
                            mSet->get_jacobian()(
                                    { tMasterResStartIndex, tMasterResStopIndex },
                                    { tMasterDepStartIndex + tDofCounter, tMasterDepStartIndex + tDofCounter } ) +=
                                    tFDScheme( 1 )( 0 ) * tMasterResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // add unperturbed slave residual contribution to dRdp
                            mSet->get_jacobian()(
                                    { tSlaveResStartIndex, tSlaveResStopIndex },
                                    { tMasterDepStartIndex + tDofCounter, tMasterDepStartIndex + tDofCounter } ) +=
                                    tFDScheme( 1 )( 0 ) * tSlaveResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

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

                            // reset the residual
                            mSet->get_residual()( 0 ).fill( 0.0 );

                            // compute the residual
                            this->compute_residual( aWStar );

                            // assemble master part of the jacobian
                            mSet->get_jacobian()(
                                    { tMasterResStartIndex, tMasterResStopIndex },
                                    { tMasterDepStartIndex + tDofCounter, tMasterDepStartIndex + tDofCounter } ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // assemble slave part of the jacobian
                            mSet->get_jacobian()(
                                    { tSlaveResStartIndex, tSlaveResStopIndex },
                                    { tMasterDepStartIndex + tDofCounter, tMasterDepStartIndex + tDofCounter } ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // get slave number of dof types
            uint tSlaveNumDofTypes = mRequestedSlaveGlobalDofTypes.size();

            // loop over the IWG dof types
            for ( uint iFI = 0; iFI < tSlaveNumDofTypes; iFI++ )
            {
                // init dof counter
                uint tDofCounter = 0;

                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iFI );

                // get the index for the dof type
                sint tSlaveDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tSlaveDepDofIndex, 0 );

                // get field interpolator for dependency dof type
                Field_Interpolator* tFI =
                        mSlaveFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

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
                        real tDeltaH = aPerturbation;
                        if ( !aUseAbsolutePerturbations )
                        {
                            tDeltaH = tDeltaH * tCoeff( iCoeffRow, iCoeffCol );
                        }

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
                            // add unperturbed master residual contribution to dRdp
                            mSet->get_jacobian()(
                                    { tMasterResStartIndex, tMasterResStopIndex },
                                    { tSlaveDepStartIndex + tDofCounter, tSlaveDepStartIndex + tDofCounter } ) +=
                                    tFDScheme( 1 )( 0 ) * tMasterResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // add unperturbed slave residual contribution to dRdp
                            mSet->get_jacobian()(
                                    { tSlaveResStartIndex, tSlaveResStopIndex },
                                    { tSlaveDepStartIndex + tDofCounter, tSlaveDepStartIndex + tDofCounter } ) +=
                                    tFDScheme( 1 )( 0 ) * tSlaveResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

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

                            // reset the residual
                            mSet->get_residual()( 0 ).fill( 0.0 );

                            // compute the residual
                            this->compute_residual( aWStar );

                            // assemble the jacobian
                            mSet->get_jacobian()(
                                    { tMasterResStartIndex, tMasterResStopIndex },
                                    { tSlaveDepStartIndex + tDofCounter, tSlaveDepStartIndex + tDofCounter } ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // assemble the jacobian
                            mSet->get_jacobian()(
                                    { tSlaveResStartIndex, tSlaveResStopIndex },
                                    { tSlaveDepStartIndex + tDofCounter, tSlaveDepStartIndex + tDofCounter } ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;
        }

        //------------------------------------------------------------------------------

        bool
        IWG::check_jacobian(
                real              aPerturbation,
                real              aEpsilon,
                real              aWStar,
                Matrix< DDRMat >& aJacobian,
                Matrix< DDRMat >& aJacobianFD,
                bool              aErrorPrint,
                bool              aUseAbsolutePerturbations )
        {
            // get residual dof type index in set, start and end indices for residual dof type
            uint tMasterDofIndex    = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartRow = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResEndRow   = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            uint tSlaveDofIndex;
            uint tSlaveResStartRow;
            uint tSlaveResEndRow;
            uint tSlaveNumRows = 0;
            if ( mSlaveGlobalDofTypes.size() > 0 )
            {
                tSlaveDofIndex    = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
                tSlaveResStartRow = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
                tSlaveResEndRow   = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );
                tSlaveNumRows     = tSlaveResEndRow - tSlaveResStartRow + 1;
            }

            // get number of master and slave rows
            uint tMasterNumRows = tMasterResEndRow - tMasterResStartRow + 1;

            // get number of cols for jacobian
            uint tNumCols = mSet->get_jacobian().n_cols();

            // set size for analytical and FD jacobians
            aJacobian.set_size( tMasterNumRows + tSlaveNumRows, tNumCols, 0.0 );
            aJacobianFD.set_size( tMasterNumRows + tSlaveNumRows, tNumCols, 0.0 );

            // compute jacobian with IWG
            this->compute_jacobian( aWStar );

            // get the computed jacobian
            aJacobian( { 0, tMasterNumRows - 1 }, { 0, tNumCols - 1 } ) =
                    mSet->get_jacobian()( { tMasterResStartRow, tMasterResEndRow }, { 0, tNumCols - 1 } );
            if ( tSlaveNumRows > 0 )
            {
                aJacobian( { tMasterNumRows, tMasterNumRows + tSlaveNumRows - 1 }, { 0, tNumCols - 1 } ) =
                        mSet->get_jacobian()( { tSlaveResStartRow, tSlaveResEndRow }, { 0, tNumCols - 1 } );
            }

            // reset the jacobian
            mSet->get_jacobian().fill( 0.0 );

            // compute jacobian by FD
            this->compute_jacobian_FD( aWStar, aPerturbation, fem::FDScheme_Type::POINT_5, aUseAbsolutePerturbations );

            // get the computed jacobian
            aJacobianFD( { 0, tMasterNumRows - 1 }, { 0, tNumCols - 1 } ) =
                    mSet->get_jacobian()( { tMasterResStartRow, tMasterResEndRow }, { 0, tNumCols - 1 } );
            if ( tSlaveNumRows > 0 )
            {
                aJacobianFD( { tMasterNumRows, tMasterNumRows + tSlaveNumRows - 1 }, { 0, tNumCols - 1 } ) =
                        mSet->get_jacobian()( { tSlaveResStartRow, tSlaveResEndRow }, { 0, tNumCols - 1 } );
            }

            // check that matrices to compare have same size
            MORIS_ERROR(
                    ( aJacobian.n_rows() == aJacobianFD.n_rows() ) && ( aJacobian.n_cols() == aJacobianFD.n_cols() ),
                    "IWG::check_jacobian - matrices to check do not share same dimensions." );

            // define a boolean for check
            bool tCheckJacobian = true;

            // define a real for absolute difference
            real tAbsolute = 0.0;

            // define a real for relative difference
            real tRelative = 0.0;

            for ( uint iiJac = 0; iiJac < aJacobian.n_rows(); iiJac++ )
            {
                for ( uint jjJac = 0; jjJac < aJacobian.n_cols(); jjJac++ )
                {
                    // get absolute difference
                    tAbsolute = std::abs( aJacobian( iiJac, jjJac ) - aJacobianFD( iiJac, jjJac ) );

                    // get relative difference
                    tRelative = std::abs( ( aJacobianFD( iiJac, jjJac ) - aJacobian( iiJac, jjJac ) ) / aJacobianFD( iiJac, jjJac ) );

                    // update check value
                    tCheckJacobian = tCheckJacobian && ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) );

                    // debug print
                    if ( ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) ) == false )
                    {
                        if ( aErrorPrint )
                        {
                            std::cout << "iiJac " << iiJac << " - jjJac " << jjJac << "\n"
                                      << std::flush;
                            std::cout << "aJacobian( iiJac, jjJac )   " << std::setprecision( 12 ) << aJacobian( iiJac, jjJac ) << "\n"
                                      << std::flush;
                            std::cout << "aJacobianFD( iiJac, jjJac ) " << std::setprecision( 12 ) << aJacobianFD( iiJac, jjJac ) << "\n"
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
            return tCheckJacobian;
        }

        //------------------------------------------------------------------------------

        // FIXME: This function needs to go, functionality will be integrated into the usual check jacobian function
        bool
        IWG::check_jacobian_multi_residual(
                real              aPerturbation,
                real              aEpsilon,
                real              aWStar,
                Matrix< DDRMat >& aJacobian,
                Matrix< DDRMat >& aJacobianFD,
                bool              aErrorPrint,
                bool              aMaxErrorPrint,
                moris::real       aFDtolerance,
                bool              aUseAbsolutePerturbations )
        {
            // check if FD tolerance is defined
            if ( aFDtolerance < 0.0 )
            {
                aFDtolerance = aEpsilon;
            }

            // get residual dof type index in set, start and end indices for residual dof type
            uint tNumResidualDofs     = mResidualDofType.size();
            uint tMasterFirstDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterLastDofIndex  = mSet->get_dof_index_for_type( mResidualDofType( tNumResidualDofs - 1 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartRow   = mSet->get_res_dof_assembly_map()( tMasterFirstDofIndex )( 0, 0 );
            uint tMasterResEndRow     = mSet->get_res_dof_assembly_map()( tMasterLastDofIndex )( 0, 1 );

            // get number of master and slave rows
            uint tNumRows = tMasterResEndRow - tMasterResStartRow + 1;

            // get number of cols for jacobian
            uint tNumCols = mSet->get_jacobian().n_cols();

            // set size for analytical and FD jacobians
            aJacobian.set_size( tNumRows, tNumCols, 0.0 );
            aJacobianFD.set_size( tNumRows, tNumCols, 0.0 );

            // compute jacobian with IWG
            this->compute_jacobian( aWStar );

            // get the computed jacobian
            aJacobian = mSet->get_jacobian();

            // reset the jacobian
            mSet->get_jacobian().fill( 0.0 );

            // compute jacobian by FD
            // this->compute_jacobian_FD( aWStar, aPerturbation );
            {
                // storage residual value
                Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

                // get the FD scheme info
                moris::Cell< moris::Cell< real > > tFDScheme;
                fd_scheme( fem::FDScheme_Type::POINT_5, tFDScheme );
                uint tNumFDPoints = tFDScheme( 0 ).size();

                // get master number of dof types
                uint tMasterNumDofTypes = mRequestedMasterGlobalDofTypes.size();

                // reset and evaluate the residual plus
                mSet->get_residual()( 0 ).fill( 0.0 );
                this->compute_residual( aWStar );
                Matrix< DDRMat > tResidual = mSet->get_residual()( 0 );

                // loop over the IWG dof types
                for ( uint iFI = 0; iFI < tMasterNumDofTypes; iFI++ )
                {
                    // init dof counter
                    uint tDofCounter = 0;

                    // get the dof type
                    Cell< MSI::Dof_Type >& tDofType = mRequestedMasterGlobalDofTypes( iFI );

                    // get the index for the dof type
                    sint tMasterDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                    uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterFirstDofIndex )( tMasterDepDofIndex, 0 );

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
                            if ( aUseAbsolutePerturbations )
                            {
                                tDeltaH = aPerturbation;
                            }

                            // check that perturbation is not zero
                            if ( std::abs( tDeltaH ) < 1e-12 )
                            {
                                tDeltaH = aPerturbation;
                            }

                            // set starting point for FD
                            uint tStartPoint = 0;

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

                                // reset the residual
                                mSet->get_residual()( 0 ).fill( 0.0 );

                                // compute the residual
                                this->compute_residual( aWStar );

                                // assemble the Jacobian
                                mSet->get_jacobian()(
                                        { tMasterResStartRow, tMasterResEndRow },
                                        { tMasterDepStartIndex + tDofCounter, tMasterDepStartIndex + tDofCounter } ) +=
                                        tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                            }
                            // update dof counter
                            tDofCounter++;
                        }
                    }
                    // reset the coefficients values
                    tFI->set_coeff( tCoeff );
                }

                // reset the value of the residual
                mSet->get_residual()( 0 ) = tResidualStore;
            }

            // get the computed jacobian
            aJacobianFD = mSet->get_jacobian();

            // check that matrices to compare have same size
            MORIS_ERROR(
                    ( aJacobian.n_rows() == aJacobianFD.n_rows() ) && ( aJacobian.n_cols() == aJacobianFD.n_cols() ),
                    "IWG::check_jacobian - matrices to check do not share same dimensions." );

            // define a boolean for check
            bool tCheckJacobian = true;

            // define a real for absolute difference
            real tAbsolute    = 0.0;
            real tMaxAbsolute = 0.0;

            // define a real for relative difference
            real tRelative    = 0.0;
            real tMaxRelative = 0.0;

            for ( uint iiJac = 0; iiJac < aJacobian.n_rows(); iiJac++ )
            {
                for ( uint jjJac = 0; jjJac < aJacobian.n_cols(); jjJac++ )
                {
                    // get absolute difference
                    tAbsolute    = std::abs( aJacobian( iiJac, jjJac ) - aJacobianFD( iiJac, jjJac ) );
                    tMaxAbsolute = std::max( tAbsolute, tMaxAbsolute );

                    // get relative difference
                    tRelative = std::abs( ( aJacobianFD( iiJac, jjJac ) - aJacobian( iiJac, jjJac ) ) / aJacobianFD( iiJac, jjJac ) );
                    if ( tAbsolute > aFDtolerance )
                    {
                        tMaxRelative = std::max( tRelative, tMaxRelative );
                    }

                    // update check value
                    tCheckJacobian = tCheckJacobian && ( ( tAbsolute < aFDtolerance ) || ( tRelative < aEpsilon ) );

                    // debug print
                    if ( ( ( tAbsolute < aFDtolerance ) || ( tRelative < aEpsilon ) ) == false )
                    {
                        if ( aErrorPrint )
                        {
                            std::cout << "iiJac " << iiJac << " - jjJac " << jjJac << "\n"
                                      << std::flush;
                            std::cout << "aJacobian( iiJac, jjJac )   " << std::setprecision( 12 ) << aJacobian( iiJac, jjJac ) << "\n"
                                      << std::flush;
                            std::cout << "aJacobianFD( iiJac, jjJac ) " << std::setprecision( 12 ) << aJacobianFD( iiJac, jjJac ) << "\n"
                                      << std::flush;
                            std::cout << "Absolute difference " << tAbsolute << "\n"
                                      << std::flush;
                            std::cout << "Relative difference " << tRelative << "\n"
                                      << std::flush;
                        }
                    }
                }
            }

            // print maximum difference
            if ( aMaxErrorPrint )
            {
                std::cout << "Maximum absolute difference " << tMaxAbsolute << "\n"
                          << std::flush;
                std::cout << "Maximum relative difference " << tMaxRelative << "\n"
                          << std::flush;
            }

            // return bool
            return tCheckJacobian;
        }

        //------------------------------------------------------------------------------

        real
        IWG::build_perturbation_size(
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
        IWG::build_perturbation_size_relative(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance )
        {
            // compute the perturbation value using fraction of maximum allowable perturbation
            real tDeltaH = std::abs( aPerturbation * aMaxPerturbation );

            // check that maximum perturbation size is larger than tolerance
            MORIS_ASSERT( aMaxPerturbation >= aTolerance,
                    "Error: maximum perturbation size is smaller than tolerance.\n" );

            // determine actual tolerance (only useful when above assert inactive)
            const real tActualTol = std::min( aMaxPerturbation, aTolerance );

            // compute perturbation such that it is not smaller than tolerance
            // and not larger than maximum value
            tDeltaH = std::max( std::min( tDeltaH, aMaxPerturbation ), tActualTol );

            return tDeltaH;
        }

        //------------------------------------------------------------------------------

        real
        IWG::build_perturbation_size_absolute(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance )
        {
            // check that maximum perturbation size is larger than tolerance
            MORIS_ASSERT( aMaxPerturbation >= aTolerance,
                    "Error: maximum perturbation size is smaller than tolerance.\n" );

            // determine actual tolerance (only useful when above assert inactive)
            const real tActualTol = std::min( aMaxPerturbation, aTolerance );

            // check that absolute value of perturbation is not smaller than tolerance
            // and not larger than maximum value
            return std::max( std::min( std::abs( aPerturbation ), aMaxPerturbation ), tActualTol );
        }

        //------------------------------------------------------------------------------

        real
        IWG::check_ig_coordinates_inside_ip_element(
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
                    tMaxIP >= aCoefficientToPerturb - tDeltaH && tMinIP <= aCoefficientToPerturb + tDeltaH,
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

                    // check for correctness of perturbation size for forward FD
                    MORIS_ASSERT( tDeltaH < tMaxIP - aCoefficientToPerturb,
                            "ERROR: forward perturbation size exceeds limits of interpolation element:\n",
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
                    // check for correctness of perturbation size for central FD
                    MORIS_ASSERT(
                            tDeltaH < tMaxIP - aCoefficientToPerturb && tDeltaH < aCoefficientToPerturb - tMinIP,
                            "ERROR: central perturbation size exceed limits of interpolation element:\n"
                            "dim: %d  minIP: %e  maxIP: %e  cordIG: %e  maxPert: %e  delta: %e  precPert: %e\n.",
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
        IWG::select_dRdp_FD_geometry_bulk(
                moris::real                        aWStar,
                moris::real                        aPerturbation,
                fem::FDScheme_Type                 aFDSchemeType,
                Matrix< DDSMat >&                  aGeoLocalAssembly,
                moris::Cell< Matrix< IndexMat > >& aVertexIndices )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the GI for the IG element considered
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();
            Geometry_Interpolator* tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();

            // get the residual dof type index in the set
            uint tResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
            uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

            // reset, evaluate and store the residual for unperturbed case
            mSet->get_residual()( 0 ).fill( 0.0 );
            this->compute_residual( aWStar );
            Matrix< DDRMat > tResidual =
                    mSet->get_residual()( 0 )(
                            { tResDofAssemblyStart, tResDofAssemblyStop },
                            { 0, 0 } );

            // get number of master GI bases and space dimensions
            uint tDerNumBases      = tIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

            // coefficients for dv type wrt which derivative is computed
            Matrix< DDRMat > tCoeff      = tIGGI->get_space_coeff();          // get nodal coordinates of integration element
            Matrix< DDRMat > tParamCoeff = tIGGI->get_space_param_coeff();    // get IP natural coordinate of integration nodes

            Matrix< DDRMat > tEvaluationPoint;    // get IG natural coordinates of quadrature points
            tIGGI->get_space_time( tEvaluationPoint );

            real tGPWeight = aWStar / tIGGI->det_J();    // reconstruct integration weight

            // init perturbation
            real tDeltaH = 0.0;

            // init FD scheme
            moris::Cell< moris::Cell< real > > tFDScheme;

            // loop over the spatial directions
            for ( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
            {
                // loop over the IG nodes
                for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    // if pdv is active
                    if ( tPdvAssemblyIndex != -1 )
                    {
                        // provide adapted perturbation and FD scheme considering ip element boundaries
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

                        // if backward or forward add unperturbed contribution
                        if ( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed residual contribution to dRdp
                            mSet->get_drdpgeo()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                    tFDScheme( 1 )( 0 ) * tResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over point of FD scheme
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients, i.e. the nodal coordinate of integration element
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient, i.e. the nodal coordinate of integration element
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // setting the perturbed coefficients, i.e. the nodal coordinate of integration element
                            tIGGI->set_space_coeff( tCoeffPert );

                            // update natural coordinates of IG nodes in IP element
                            Matrix< DDRMat > tXCoords  = tCoeffPert.get_row( iCoeffRow );
                            Matrix< DDRMat > tXiCoords = tParamCoeff.get_row( iCoeffRow );

                            tIPGI->update_local_coordinates( tXCoords, tXiCoords );

                            Matrix< DDRMat > tParamCoeffPert     = tParamCoeff;
                            tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();
                            tIGGI->set_space_param_coeff( tParamCoeffPert );

                            // set evaluation point for interpolators (FIs and GIs)
                            mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset and evaluate the residual plus
                            mSet->get_residual()( 0 ).fill( 0.0 );
                            real tWStarPert = tGPWeight * tIGGI->det_J();
                            this->compute_residual( tWStarPert );

                            // evaluate dRdpGeo
                            mSet->get_drdpgeo()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                }
            }

            // reset the coefficients values
            tIGGI->set_space_coeff( tCoeff );
            tIGGI->set_space_param_coeff( tParamCoeff );
            mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;

            // add contribution of cluster measure to dRdp
            if ( mActiveCMEAFlag )
            {
                // add their contribution to dQIdp
                this->add_cluster_measure_dRdp_FD_geometry(
                        aWStar,
                        aPerturbation,
                        aFDSchemeType );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ),
                    "IWG::compute_dRdp_FD_geometry - dRdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG::select_dRdp_FD_geometry_sideset(
                moris::real                        aWStar,
                moris::real                        aPerturbation,
                fem::FDScheme_Type                 aFDSchemeType,
                Matrix< DDSMat >&                  aGeoLocalAssembly,
                moris::Cell< Matrix< IndexMat > >& aVertexIndices )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the GI for the IG element considered
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();
            Geometry_Interpolator* tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();

            // get the residual dof type index in the set
            uint tResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
            uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

            // reset, evaluate and store the residual for unperturbed case
            mSet->get_residual()( 0 ).fill( 0.0 );
            this->compute_residual( aWStar );
            Matrix< DDRMat > tResidual = mSet->get_residual()( 0 )(
                    { tResDofAssemblyStart, tResDofAssemblyStop },
                    { 0, 0 } );

            // store unperturbed xyz
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

            // init perturbation
            real tDeltaH = 0.0;

            // init FD scheme
            moris::Cell< moris::Cell< real > > tFDScheme;

            // get number of master GI bases and space dimensions
            uint tDerNumBases      = tIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

            // loop over the spatial directions
            for ( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
            {
                // loop over the IG nodes
                for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    // if pdv is active
                    if ( tPdvAssemblyIndex != -1 )
                    {
                        // provide adapted perturbation and FD scheme considering ip element boundaries
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

                        // if backward or forward add unperturbed contribution
                        if ( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed residual contribution to dRdp
                            mSet->get_drdpgeo()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                    tFDScheme( 1 )( 0 ) * tResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

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
                            mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );

                            // reset the normal
                            Matrix< DDRMat > tNormalPert;
                            tIGGI->get_normal( tNormalPert );
                            this->set_normal( tNormalPert );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset and evaluate the residual
                            mSet->get_residual()( 0 ).fill( 0.0 );
                            real tWStarPert = tGPWeight * tIGGI->det_J();
                            this->compute_residual( tWStarPert );

                            // evaluate dRdpGeo
                            mSet->get_drdpgeo()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                }
            }

            // reset xyz values
            tIGGI->set_space_coeff( tCoeff );

            // reset local coordinates values
            tIGGI->set_space_param_coeff( tParamCoeff );

            // reset evaluation point
            mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );

            // reset normal
            this->set_normal( tNormal );

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;

            // add contribution of cluster measure to dRdp
            if ( mActiveCMEAFlag )
            {
                // add their contribution to dQIdp
                this->add_cluster_measure_dRdp_FD_geometry(
                        aWStar,
                        aPerturbation,
                        aFDSchemeType );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ),
                    "IWG::compute_dRdp_FD_geometry - dRdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG::select_dRdp_FD_geometry_time_sideset(
                moris::real                        aWStar,
                moris::real                        aPerturbation,
                fem::FDScheme_Type                 aFDSchemeType,
                Matrix< DDSMat >&                  aGeoLocalAssembly,
                moris::Cell< Matrix< IndexMat > >& aVertexIndices )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the GI for the IG element considered
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();
            Geometry_Interpolator* tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();
            Geometry_Interpolator* tIGGIPrevious =
                    mSet->get_field_interpolator_manager_previous_time()->get_IG_geometry_interpolator();

            // get the residual dof type index in the set
            uint tResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
            uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

            // reset, evaluate and store the residual for unperturbed case
            mSet->get_residual()( 0 ).fill( 0.0 );
            this->compute_residual( aWStar );
            Matrix< DDRMat > tResidual = mSet->get_residual()( 0 )(
                    { tResDofAssemblyStart, tResDofAssemblyStop },
                    { 0, 0 } );

            // store unperturbed xyz
            Matrix< DDRMat > tCoeff = tIGGI->get_space_coeff();

            // store unperturbed local coordinates
            Matrix< DDRMat > tParamCoeff = tIGGI->get_space_param_coeff();

            // store unperturbed evaluation point
            Matrix< DDRMat > tEvaluationPoint;
            tIGGI->get_space_time( tEvaluationPoint );

            // store unperturbed evaluation point weight
            real tGPWeight = aWStar / tIGGI->det_J();

            // init perturbation
            real tDeltaH = 0.0;

            // init FD scheme
            moris::Cell< moris::Cell< real > > tFDScheme;

            // get number of master GI bases and space dimensions
            uint tDerNumBases      = tIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

            // loop over the spatial directions
            for ( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
            {
                // loop over the IG nodes
                for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    // if pdv is active
                    if ( tPdvAssemblyIndex != -1 )
                    {
                        // provide adapted perturbation and FD scheme considering ip element boundaries
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

                        // if backward or forward add unperturbed contribution
                        if ( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed residual contribution to dRdp
                            mSet->get_drdpgeo()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                    tFDScheme( 1 )( 0 ) * tResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

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
                            tIGGIPrevious->set_space_coeff( tCoeffPert );

                            // update local coordinates
                            Matrix< DDRMat > tXCoords  = tCoeffPert.get_row( iCoeffRow );
                            Matrix< DDRMat > tXiCoords = tParamCoeff.get_row( iCoeffRow );
                            tIPGI->update_local_coordinates( tXCoords, tXiCoords );
                            Matrix< DDRMat > tParamCoeffPert     = tParamCoeff;
                            tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();
                            tIGGI->set_space_param_coeff( tParamCoeffPert );
                            tIGGIPrevious->set_space_param_coeff( tParamCoeffPert );

                            // set evaluation point for interpolators (FIs and GIs)
                            mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );
                            mSet->get_field_interpolator_manager_previous_time()->set_space_time_from_local_IG_point( tEvaluationPoint );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset and evaluate the residual plus
                            mSet->get_residual()( 0 ).fill( 0.0 );
                            real tWStarPert = tGPWeight * tIGGI->det_J();
                            this->compute_residual( tWStarPert );

                            // evaluate dRdpGeo
                            mSet->get_drdpgeo()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                }
            }

            // reset xyz values
            tIGGI->set_space_coeff( tCoeff );
            tIGGIPrevious->set_space_coeff( tCoeff );

            // reset local coordinates values
            tIGGI->set_space_param_coeff( tParamCoeff );
            tIGGIPrevious->set_space_param_coeff( tParamCoeff );

            // reset evaluation point
            mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );
            mSet->get_field_interpolator_manager_previous_time()->set_space_time_from_local_IG_point( tEvaluationPoint );

            // reset the coefficients values
            tIGGI->set_space_coeff( tCoeff );
            tIGGIPrevious->set_space_coeff( tCoeff );

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;

            // add contribution of cluster measure to dRdp
            if ( mActiveCMEAFlag )
            {
                // add their contribution to dQIdp
                this->add_cluster_measure_dRdp_FD_geometry(
                        aWStar,
                        aPerturbation,
                        aFDSchemeType );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ),
                    "IWG::compute_dRdp_FD_geometry - dRdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG::select_dRdp_FD_geometry_double(
                moris::real                        aWStar,
                moris::real                        aPerturbation,
                fem::FDScheme_Type                 aFDSchemeType,
                Matrix< DDSMat >&                  aGeoLocalAssembly,
                moris::Cell< Matrix< IndexMat > >& aVertexIndices )
        {
            // unpack vertex indices
            Matrix< IndexMat >& aMasterVertexIndices = aVertexIndices( 0 );
            Matrix< IndexMat >& aSlaveVertexIndices  = aVertexIndices( 1 );

            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get requested geometry pdv types
            moris::Cell< PDV_Type > tRequestedGeoPdvType;
            mSet->get_ig_unique_dv_types_for_set( tRequestedGeoPdvType );

            // get the pdv active flags from the FEM IG nodes
            Matrix< DDSMat > tAssembly;
            mSet->get_equation_model()->get_integration_xyz_pdv_assembly_indices(
                    aMasterVertexIndices,
                    tRequestedGeoPdvType,
                    tAssembly );

            // get the master GI for the IG and IP element considered
            Geometry_Interpolator* tMasterIGGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->get_IG_geometry_interpolator();
            Geometry_Interpolator* tMasterIPGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->get_IP_geometry_interpolator();

            // get the slave GI for the IG and IP element considered
            Geometry_Interpolator* tSlaveIGGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->get_IG_geometry_interpolator();

            // get the master residual dof type index in the set
            uint tMasterResDofIndex = mSet->get_dof_index_for_type(
                    mResidualDofType( 0 )( 0 ),
                    mtk::Master_Slave::MASTER );
            uint tMasterResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tMasterResDofIndex )( 0, 0 );
            uint tMasterResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tMasterResDofIndex )( 0, 1 );

            // get the slave residual dof type index in the set (if exists)
            sint tSlaveResDofIndex = mSet->get_dof_index_for_type(
                    mResidualDofType( 0 )( 0 ),
                    mtk::Master_Slave::SLAVE );

            uint tSlaveResDofAssemblyStart = MORIS_UINT_MAX;
            uint tSlaveResDofAssemblyStop  = MORIS_UINT_MAX;

            if ( tSlaveResDofIndex != -1 )
            {
                tSlaveResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tSlaveResDofIndex )( 0, 0 );
                tSlaveResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tSlaveResDofIndex )( 0, 1 );
            }

            // init perturbation
            real tDeltaH = 0.0;

            // get GP weight
            real tGPWeight = aWStar / tMasterIGGI->det_J();

            // get master coeff
            Matrix< DDRMat > tMasterCoeff      = tMasterIGGI->get_space_coeff();
            Matrix< DDRMat > tMasterParamCoeff = tMasterIGGI->get_space_param_coeff();
            Matrix< DDRMat > tMasterEvaluationPoint;
            tMasterIGGI->get_space_time( tMasterEvaluationPoint );
            Matrix< DDRMat > tMasterNormal;
            tMasterIGGI->get_normal( tMasterNormal );

            // coefficients for dv type wrt which derivative is computed
            Matrix< DDRMat > tSlaveCoeff      = tSlaveIGGI->get_space_coeff();
            Matrix< DDRMat > tSlaveParamCoeff = tSlaveIGGI->get_space_param_coeff();
            Matrix< DDRMat > tSlaveEvaluationPoint;
            tSlaveIGGI->get_space_time( tSlaveEvaluationPoint );

            // reset, evaluate and store the residual for unperturbed case
            mSet->get_residual()( 0 ).fill( 0.0 );
            this->compute_residual( aWStar );
            Matrix< DDRMat > tMasterResidual =
                    mSet->get_residual()( 0 )( { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop }, { 0, 0 } );

            Matrix< DDRMat > tSlaveResidual;
            if ( tSlaveResDofIndex != -1 )
            {
                tSlaveResidual =
                        mSet->get_residual()( 0 )( { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop }, { 0, 0 } );
            }

            // get number of master GI bases and space dimensions
            uint tDerNumBases      = tMasterIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tMasterIPGI->get_number_of_space_dimensions();

            // init FD scheme
            moris::Cell< moris::Cell< real > > tFDScheme;

            // loop over the IG nodes
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // find the node on the slave side
                sint tSlaveNodeLocalIndex = -1;
                for ( uint iNode = 0; iNode < tDerNumBases; iNode++ )
                {
                    if ( aMasterVertexIndices( iCoeffRow ) == aSlaveVertexIndices( iNode ) )
                    {
                        tSlaveNodeLocalIndex = iNode;
                        break;
                    }
                }
                MORIS_ERROR( tSlaveNodeLocalIndex != -1, "IWG::compute_dRdp_FD_geometry_double - slave index not found." );

                // loop over the spatial directions
                for ( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    if ( tPdvAssemblyIndex != -1 )
                    {
                        // provide adapted perturbation and FD scheme considering ip element boundaries
                        fem::FDScheme_Type tUsedFDScheme = aFDSchemeType;
                        tDeltaH                          = this->check_ig_coordinates_inside_ip_element(
                                aPerturbation,
                                tMasterCoeff( iCoeffRow, iCoeffCol ),
                                iCoeffCol,
                                tUsedFDScheme );

                        // finalize FD scheme
                        fd_scheme( tUsedFDScheme, tFDScheme );
                        uint tNumFDPoints = tFDScheme( 0 ).size();

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward add unperturbed contribution
                        if ( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed master residual contribution to dRdp
                            mSet->get_drdpgeo()(
                                    { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                                    { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                    tFDScheme( 1 )( 0 ) * tMasterResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // add unperturbed slave residual contribution to dRdp
                            if ( tSlaveResDofIndex != -1 )
                            {
                                mSet->get_drdpgeo()(
                                        { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                        { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                        tFDScheme( 1 )( 0 ) * tSlaveResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                            }

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over point of FD scheme
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tMasterCoeffPert = tMasterCoeff;
                            Matrix< DDRMat > tSlaveCoeffPert  = tSlaveCoeff;

                            // perturb the coefficient
                            tMasterCoeffPert( iCoeffRow, iCoeffCol ) +=
                                    tFDScheme( 0 )( iPoint ) * tDeltaH;
                            tSlaveCoeffPert( tSlaveNodeLocalIndex, iCoeffCol ) +=
                                    tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // setting the perturbed coefficients
                            tMasterIGGI->set_space_coeff( tMasterCoeffPert );
                            tSlaveIGGI->set_space_coeff( tSlaveCoeffPert );

                            // update local coordinates
                            Matrix< DDRMat > tXCoords  = tMasterCoeffPert.get_row( iCoeffRow );
                            Matrix< DDRMat > tXiCoords = tMasterParamCoeff.get_row( iCoeffRow );
                            tMasterIPGI->update_local_coordinates( tXCoords, tXiCoords );
                            Matrix< DDRMat > tMasterParamCoeffPert               = tMasterParamCoeff;
                            tMasterParamCoeffPert.get_row( iCoeffRow )           = tXiCoords.matrix_data();
                            Matrix< DDRMat > tSlaveParamCoeffPert                = tSlaveParamCoeff;
                            tSlaveParamCoeffPert.get_row( tSlaveNodeLocalIndex ) = tXiCoords.matrix_data();

                            tMasterIGGI->set_space_param_coeff( tMasterParamCoeffPert );
                            tSlaveIGGI->set_space_param_coeff( tSlaveParamCoeffPert );

                            // set evaluation point for interpolators (FIs and GIs)
                            mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->set_space_time_from_local_IG_point( tMasterEvaluationPoint );
                            mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->set_space_time_from_local_IG_point( tSlaveEvaluationPoint );

                            // reset the normal
                            Matrix< DDRMat > tNormalPert;
                            tMasterIGGI->get_normal( tNormalPert );
                            this->set_normal( tNormalPert );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset and evaluate the residual plus
                            mSet->get_residual()( 0 ).fill( 0.0 );
                            real tWStarPert = tGPWeight * tMasterIGGI->det_J();
                            this->compute_residual( tWStarPert );

                            // evaluate dMasterRdpGeo
                            mSet->get_drdpgeo()(
                                    { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                                    { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // evaluate dSlaveRdpGeo
                            if ( tSlaveResDofIndex != -1 )
                            {
                                mSet->get_drdpgeo()(
                                        { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                        { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                        tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                            }
                        }
                    }
                }
            }
            // reset the coefficients values
            tMasterIGGI->set_space_coeff( tMasterCoeff );
            tMasterIGGI->set_space_param_coeff( tMasterParamCoeff );
            tSlaveIGGI->set_space_coeff( tSlaveCoeff );
            tSlaveIGGI->set_space_param_coeff( tSlaveParamCoeff );
            mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->set_space_time_from_local_IG_point( tMasterEvaluationPoint );
            mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->set_space_time_from_local_IG_point( tSlaveEvaluationPoint );
            this->set_normal( tMasterNormal );

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;

            // add contribution of cluster measure to dRdp
            if ( mActiveCMEAFlag )
            {
                // add their contribution to dQIdp
                this->add_cluster_measure_dRdp_FD_geometry_double(
                        aWStar,
                        aPerturbation,
                        aFDSchemeType );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ),
                    "IWG::compute_dRdp_FD_geometry_double - dRdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG::add_cluster_measure_dRdp_FD_geometry(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the residual dof type index in the set
            uint tResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
            uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

            // reset, evaluate and store the residual for unperturbed case
            mSet->get_residual()( 0 ).fill( 0.0 );
            this->compute_residual( aWStar );
            Matrix< DDRMat > tResidual =
                    mSet->get_residual()( 0 )(
                            { tResDofAssemblyStart, tResDofAssemblyStop },
                            { 0, 0 } );

            // init perturbation
            real tDeltaH = 0.0;

            // init FD scheme
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumFDPoints = tFDScheme( 0 ).size();

            // loop over the cluster measures
            for ( uint iCMEA = 0; iCMEA < mCluster->get_cluster_measures().size(); iCMEA++ )
            {
                // get treated cluster measure
                std::shared_ptr< Cluster_Measure >& tClusterMeasure =
                        mCluster->get_cluster_measures()( iCMEA );

                // evaluate the perturbation
                tDeltaH = this->build_perturbation_size(
                        aPerturbation,
                        tClusterMeasure->val()( 0 ),
                        tClusterMeasure->val()( 0 ) );

                // number of pdv to assemble
                uint tNumPdvToAssemble = tClusterMeasure->dMEAdPDV().numel() - 1;

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed residual contribution to dRdp
                    mSet->get_drdpgeo()(
                            { tResDofAssemblyStart, tResDofAssemblyStop },
                            { 0, tNumPdvToAssemble } ) +=
                            tFDScheme( 1 )( 0 ) * tResidual * tClusterMeasure->dMEAdPDV() / ( tFDScheme( 2 )( 0 ) * tDeltaH );

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

                    // reset and evaluate the residual plus
                    mSet->get_residual()( 0 ).fill( 0.0 );
                    this->compute_residual( aWStar );

                    // evaluate dRdpGeo
                    mSet->get_drdpgeo()(
                            { tResDofAssemblyStart, tResDofAssemblyStop },
                            { 0, tNumPdvToAssemble } ) +=
                                    tFDScheme( 1 )( iPoint ) * 
                                    mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) * 
                                    tClusterMeasure->dMEAdPDV() / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // reset cluster measures
                    mCluster->reset_cluster_measure();
                }
            }

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ),
                    "IWG::compute_dRdp_FD_geometry - dRdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG::add_cluster_measure_dRdp_FD_geometry_double(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the residual dof type index in the set
            uint tMasterDofIndex            = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the slave residual dof type index in the set
            uint tSlaveResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tSlaveResDofIndex )( 0, 0 );
            uint tSlaveResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tSlaveResDofIndex )( 0, 1 );

            // reset, evaluate and store the residual for unperturbed case
            mSet->get_residual()( 0 ).fill( 0.0 );
            this->compute_residual( aWStar );
            Matrix< DDRMat > tMasterResidual =
                    mSet->get_residual()( 0 )( { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop }, { 0, 0 } );
            Matrix< DDRMat > tSlaveResidual =
                    mSet->get_residual()( 0 )( { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop }, { 0, 0 } );

            // init perturbation
            real tDeltaH = 0.0;

            // init FD scheme
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumFDPoints = tFDScheme( 0 ).size();

            // loop over the cluster measures
            for ( uint iCMEA = 0; iCMEA < mCluster->get_cluster_measures().size(); iCMEA++ )
            {
                // get treated cluster measure
                std::shared_ptr< Cluster_Measure >& tClusterMeasure =
                        mCluster->get_cluster_measures()( iCMEA );

                // evaluate the perturbation
                tDeltaH = this->build_perturbation_size(
                        aPerturbation,
                        tClusterMeasure->val()( 0 ),
                        tClusterMeasure->val()( 0 ) );

                // get end pdv index
                uint tEndPdvIndex = tClusterMeasure->dMEAdPDV().numel() - 1;

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward add unperturbed contribution
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed master residual contribution to dRdp
                    mSet->get_drdpgeo()(
                            { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                            { 0, tEndPdvIndex } ) +=
                            tFDScheme( 1 )( 0 ) * tMasterResidual * tClusterMeasure->dMEAdPDV() / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // add unperturbed slave residual contribution to dRdp
                    mSet->get_drdpgeo()(
                            { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                            { 0, tEndPdvIndex } ) +=
                            tFDScheme( 1 )( 0 ) * tSlaveResidual * tClusterMeasure->dMEAdPDV() / ( tFDScheme( 2 )( 0 ) * tDeltaH );

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

                    // reset and evaluate the residual plus
                    mSet->get_residual()( 0 ).fill( 0.0 );
                    this->compute_residual( aWStar );

                    // evaluate dMasterRdpGeo
                    mSet->get_drdpgeo()(
                            { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                            { 0, tEndPdvIndex } ) +=
                            tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop }, { 0, 0 } ) * tClusterMeasure->dMEAdPDV() / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // evaluate dSlaveRdpGeo
                    mSet->get_drdpgeo()(
                            { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                            { 0, tEndPdvIndex } ) +=
                            tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop }, { 0, 0 } ) * tClusterMeasure->dMEAdPDV() / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                    // reset cluster measures
                    mCluster->reset_cluster_measure();
                }
            }

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ),
                    "IWG::compute_dRdp_FD_geometry - dRdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG::select_dRdp_FD_material(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumFDPoints = tFDScheme( 0 ).size();

            // get the requested ip pdv types
            moris::Cell< moris::Cell< PDV_Type > > tRequestedPdvTypes;
            mSet->get_ip_dv_types_for_set( tRequestedPdvTypes );

            // get number of requested dv types
            uint tNumDvType = tRequestedPdvTypes.size();

            // get the residual dof type index in the set
            uint tResDofIndex = mSet->get_dof_index_for_type(
                    mResidualDofType( 0 )( 0 ),
                    mtk::Master_Slave::MASTER );
            uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
            uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

            // reset, evaluate and store the residual for unperturbed case
            mSet->get_residual()( 0 ).fill( 0.0 );
            this->compute_residual( aWStar );
            Matrix< DDRMat > tResidual =
                    mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } );

            // init perturbation
            real tDeltaH = 0.0;

            // loop over the dv types associated with a FI
            for ( uint iFI = 0; iFI < tNumDvType; iFI++ )
            {
                // get dv index
                sint tDvDepIndex = mSet->get_dv_index_for_type(
                        tRequestedPdvTypes( iFI )( 0 ),
                        mtk::Master_Slave::MASTER );

                // get the FI for the dv type
                Field_Interpolator* tFI =
                        mMasterFIManager->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // init coeff counter
                uint tCoeffCounter = 0;

                // loop over the coefficient column
                for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficient row
                    for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // compute the perturbation absolute value
                        tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if ( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // get mat pdv index
                        uint tPdvIndex = mSet->get_mat_pdv_assembly_map()( tDvDepIndex )( 0, 0 ) + tCoeffCounter;

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward FD, add unperturbed residual contribution
                        if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed master residual contribution to dRdp
                            mSet->get_drdpmat()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvIndex, tPdvIndex } ) +=
                                    tFDScheme( 1 )( 0 ) * tResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

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

                            // reset the residual
                            mSet->get_residual()( 0 ).fill( 0.0 );

                            // compute the residual
                            this->compute_residual( aWStar );

                            // evaluate dRdpMat
                            mSet->get_drdpmat()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvIndex, tPdvIndex } ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update coefficient counter
                        tCoeffCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_drdpmat() ),
                    "IWG::compute_dRdp_FD_material - dRdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG::select_dRdp_FD_material_double(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumFDPoints = tFDScheme( 0 ).size();

            // get the requested ip pdv types
            moris::Cell< moris::Cell< PDV_Type > > tRequestedPdvTypes;
            mSet->get_ip_dv_types_for_set( tRequestedPdvTypes );

            // get number of requested dv types
            uint tNumDvType = tRequestedPdvTypes.size();

            // get the master residual dof type index in the set
            uint tMasterResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tMasterResDofIndex )( 0, 0 );
            uint tMasterResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tMasterResDofIndex )( 0, 1 );

            // get the slave residual dof type index in the set
            sint tSlaveResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResDofAssemblyStart = MORIS_UINT_MAX;
            uint tSlaveResDofAssemblyStop  = MORIS_UINT_MAX;

            // check that slave side has residual dof type
            if ( tSlaveResDofIndex != -1 )
            {
                tSlaveResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tSlaveResDofIndex )( 0, 0 );
                tSlaveResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tSlaveResDofIndex )( 0, 1 );
            }

            // reset, evaluate and store the residual for unperturbed case
            mSet->get_residual()( 0 ).fill( 0.0 );

            this->compute_residual( aWStar );

            Matrix< DDRMat > tMasterResidual =
                    mSet->get_residual()( 0 )(
                            { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                            { 0, 0 } );

            Matrix< DDRMat > tSlaveResidual;
            if ( tSlaveResDofIndex != -1 )
            {
                tSlaveResidual = mSet->get_residual()( 0 )(
                        { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                        { 0, 0 } );
            }

            // init perturbation
            real tDeltaH = 0.0;

            // loop over the master dv types associated with a FI
            for ( uint iFI = 0; iFI < tNumDvType; iFI++ )
            {
                // get dv index
                sint tDvDepIndex = mSet->get_dv_index_for_type(
                        tRequestedPdvTypes( iFI )( 0 ),
                        mtk::Master_Slave::MASTER );

                // get the FI for the dv type
                Field_Interpolator* tFI =
                        mMasterFIManager->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // init coeff counter
                uint tCoeffCounter = 0;

                // loop over the coefficient column
                for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficient row
                    for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // compute the perturbation absolute value
                        tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if ( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // get mat pdv index
                        uint tPdvIndex = mSet->get_mat_pdv_assembly_map()( tDvDepIndex )( 0, 0 ) + tCoeffCounter;

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward FD, add unperturbed residual contribution
                        if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed master residual contribution to dRdp
                            mSet->get_drdpmat()(
                                    { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                                    { tPdvIndex, tPdvIndex } ) +=
                                    tFDScheme( 1 )( 0 ) * tMasterResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // add unperturbed master residual contribution to dRdp
                            if ( tSlaveResDofIndex != -1 )
                            {
                                mSet->get_drdpmat()(
                                        { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                        { tPdvIndex, tPdvIndex } ) +=
                                        tFDScheme( 1 )( 0 ) * tSlaveResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                            }

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

                            // reset the residual
                            mSet->get_residual()( 0 ).fill( 0.0 );

                            // compute the residual
                            this->compute_residual( aWStar );

                            // assemble dRMasterdpMat
                            mSet->get_drdpmat()(
                                    { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                                    { tPdvIndex, tPdvIndex } ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // assemble dRSlavedpMat
                            if ( tSlaveResDofIndex != -1 )
                            {
                                mSet->get_drdpmat()(
                                        { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                        { tPdvIndex, tPdvIndex } ) +=
                                        tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                            }
                        }
                        // update coefficient counter
                        tCoeffCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // loop over the slave dv types associated with a FI
            for ( uint iFI = 0; iFI < tNumDvType; iFI++ )
            {
                // get dv index
                sint tDvDepIndex = mSet->get_dv_index_for_type(
                        tRequestedPdvTypes( iFI )( 0 ),
                        mtk::Master_Slave::SLAVE );

                // get the FI for the dv type
                Field_Interpolator* tFI =
                        mSlaveFIManager->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // reset properties, CM and SP for IWG
                this->reset_eval_flags();

                // reset and evaluate the residual
                mSet->get_residual()( 0 ).fill( 0.0 );

                this->compute_residual( aWStar );

                Matrix< DDRMat > tMasterResidual =
                        mSet->get_residual()( 0 )(
                                { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                                { 0, 0 } );

                Matrix< DDRMat > tSlaveResidual;
                if ( tSlaveResDofIndex != -1 )
                {
                    tSlaveResidual =
                            mSet->get_residual()( 0 )(
                                    { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                    { 0, 0 } );
                }

                // init coeff counter
                uint tCoeffCounter = 0;

                // loop over the coefficient column
                for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficient row
                    for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // compute the perturbation absolute value
                        tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if ( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // get mat pdv index
                        uint tPdvIndex = mSet->get_mat_pdv_assembly_map()( tDvDepIndex )( 0, 0 ) + tCoeffCounter;

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward FD, add unperturbed residual contribution
                        if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) || ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed master residual contribution to dRdp
                            mSet->get_drdpmat()(
                                    { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                                    { tPdvIndex, tPdvIndex } ) +=
                                    tFDScheme( 1 )( 0 ) * tMasterResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // add unperturbed master residual contribution to dRdp
                            if ( tSlaveResDofIndex != -1 )
                            {
                                mSet->get_drdpmat()(
                                        { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                        { tPdvIndex, tPdvIndex } ) +=
                                        tFDScheme( 1 )( 0 ) * tSlaveResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                            }

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

                            // reset the residual
                            mSet->get_residual()( 0 ).fill( 0.0 );

                            // compute the residual
                            this->compute_residual( aWStar );

                            // assemble dRMasterdpMat
                            mSet->get_drdpmat()(
                                    { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                                    { tPdvIndex, tPdvIndex } ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // assemble dRSlavedpMat
                            if ( tSlaveResDofIndex != -1 )
                            {
                                mSet->get_drdpmat()(
                                        { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                        { tPdvIndex, tPdvIndex } ) +=
                                        tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 )( { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop }, { 0, 0 } ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                            }
                        }
                        // update coefficient counter
                        tCoeffCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_drdpmat() ),
                    "IWG::compute_dRdp_FD_material - dRdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

    }    // namespace fem
}    // namespace moris

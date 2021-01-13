/*
 * cl_FEM_IWG.cpp
 *
 *  Created on: Nov 12, 2019
 *      Author: sonne
 */

#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Model.hpp"

#include "fn_max.hpp"
#include "fn_min.hpp"

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
            for ( const std::shared_ptr< Property > & tProp : mMasterProp )
            {
                if ( tProp != nullptr )
                {
                    tProp->reset_eval_flags();
                }
            }
            
            for ( const std::shared_ptr< Property > & tProp : mSlaveProp )
            {
                if( tProp != nullptr )
                {
                    tProp->reset_eval_flags();
                }
            }

            // reset constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
            {
                if( tCM != nullptr )
                {
                    tCM->reset_eval_flags();
                }
            }
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    tCM->reset_eval_flags();
                }
            }

            // reset stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter > & tSP : mStabilizationParam )
            {
                if( tSP != nullptr )
                {
                    tSP->reset_eval_flags();
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG::set_function_pointers()
        {
            // switch on element type
            switch( mSet->get_element_type() )
            {
                case fem::Element_Type::BULK :
                {
                    m_compute_jacobian_FD      = &IWG::select_jacobian_FD;
                    m_compute_dRdp_FD_material = &IWG::select_dRdp_FD_material;
                    m_compute_dRdp_FD_geometry = &IWG::select_dRdp_FD_geometry_bulk;
                    break;
                }
                case fem::Element_Type::SIDESET :
                {
                    m_compute_jacobian_FD      = &IWG::select_jacobian_FD;
                    m_compute_dRdp_FD_material = &IWG::select_dRdp_FD_material;
                    m_compute_dRdp_FD_geometry = &IWG::select_dRdp_FD_geometry_sideset;
                    break;
                }
                case fem::Element_Type::TIME_SIDESET :
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
                default :
                    MORIS_ERROR( false, "IWG::set_function_pointers - unknown element type.");
            }
        }

        //------------------------------------------------------------------------------

        void IWG::set_phase_name(
                std::string aPhaseName,
                mtk::Master_Slave aIsMaster )
        {
            switch( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    mMasterPhaseName = aPhaseName;
                    break;
                }
                case mtk::Master_Slave::SLAVE :
                {
                    mSlavePhaseName = aPhaseName;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "IWG::set_phase_name - aIsMaster can only be master or slave.");
                }
            }
        }

        //------------------------------------------------------------------------------

        std::string IWG::get_phase_name( mtk::Master_Slave aIsMaster )
        {
            switch( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    return mMasterPhaseName;
                }
                case mtk::Master_Slave::SLAVE :
                {
                    return mSlavePhaseName;
                }
                default :
                {
                    MORIS_ERROR( false, "IWG::get_phase_name - aIsMaster can only be master or slave.");
                    return mMasterPhaseName;
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG::set_field_interpolator_manager(
                Field_Interpolator_Manager * aFieldInterpolatorManager,
                mtk::Master_Slave            aIsMaster )
        {
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
                    MORIS_ERROR( false, "IWG::set_field_interpolator_manager - can only be master or slave");
                }
            }

            // loop over the the SP
            for( const std::shared_ptr< Stabilization_Parameter > & tSP : this->get_stabilization_parameters() )
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
            for( const std::shared_ptr< Constitutive_Model > & tCM : this->get_constitutive_models( aIsMaster ) )
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
            for( const std::shared_ptr< Property > & tProp : this->get_properties( aIsMaster ) )
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

        Field_Interpolator_Manager * IWG::get_field_interpolator_manager(
                mtk::Master_Slave aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                    return mMasterFIManager;

                case mtk::Master_Slave::SLAVE :
                    return mSlaveFIManager;

                default :
                    MORIS_ERROR( false, "IWG::get_field_interpolator_manager - can only be master or slave." );
                    return mMasterFIManager;
            }
        }

        //------------------------------------------------------------------------------

        void IWG::set_field_interpolator_manager_previous_time(
                Field_Interpolator_Manager * aFieldInterpolatorManager,
                mtk::Master_Slave            aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    mMasterPreviousFIManager = aFieldInterpolatorManager;
                    break;
                }

                default :
                {
                    MORIS_ERROR( false, "IWG::set_field_interpolator_manager - can only be master");
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG::set_normal( Matrix< DDRMat > & aNormal )
        {
            mNormal = aNormal;

            // set normal for SP
            for ( const std::shared_ptr< Stabilization_Parameter > & tSP : mStabilizationParam )
            {
                if( tSP != nullptr )
                {
                    tSP->set_normal( mNormal );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG::set_interpolation_order()
        {
            // if order is already set
            if( ( mOrder != MORIS_UINT_MAX ) ) return;

            // get residual dof type interpolation order
            mtk::Interpolation_Order tInterpOrder =
                    mSet->get_field_interpolator_manager()->
                    get_field_interpolators_for_type( mResidualDofType( 0 ) )->
                    get_space_interpolation_order();

            // set the interpolation order for IWG
            switch ( tInterpOrder )
            {
                case mtk::Interpolation_Order::LINEAR :
                {
                    mOrder = 1;
                    break;
                }
                case mtk::Interpolation_Order::QUADRATIC :
                {
                    mOrder = 2;
                    break;
                }
                case mtk::Interpolation_Order::CUBIC :
                {
                    mOrder = 3;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IWG::set_interpolation_order - order not supported" );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG::set_dof_type_list(
                const moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                mtk::Master_Slave                                   aIsMaster )
        {
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
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
                    MORIS_ERROR( false, "IWG::set_dof_type_list - can only be MASTER or SLAVE.");
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< MSI::Dof_Type > > & IWG::get_dof_type_list(
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
                    MORIS_ASSERT( false, "IWG::get_dof_type_list - can only be master or slave." );
                    return mMasterDofTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG::set_dv_type_list(
                const moris::Cell< moris::Cell< PDV_Type > > & aDvTypes,
                mtk::Master_Slave                              aIsMaster )
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
                    MORIS_ERROR( false, "IWG::set_dv_type_list - can only be MASTER or SLAVE.");
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< PDV_Type > > & IWG::get_dv_type_list(
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
                    MORIS_ASSERT( false, "IWG::get_dv_type_list - can only be master or slave." );
                    return mMasterDvTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster  )
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

        moris::Cell< std::shared_ptr< Property > > & IWG::get_properties(
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
                    MORIS_ASSERT( false, "IWG::get_properties - can only be master or slave." );
                    return mMasterProp;
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                 aConstitutiveString,
                mtk::Master_Slave           aIsMaster  )
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

        moris::Cell< std::shared_ptr< Constitutive_Model > > & IWG::get_constitutive_models(
                mtk::Master_Slave aIsMaster )
        {
            // switch on master/slave
            switch( aIsMaster )
            {
                // if master
                case mtk::Master_Slave::MASTER :
                {
                    // return master property pointers
                    return mMasterCM;
                }
                // if slave
                case mtk::Master_Slave::SLAVE :
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

          void IWG::set_stabilization_parameter(
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
          }

        //------------------------------------------------------------------------------

        void IWG::get_non_unique_dof_and_dv_types(
                moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                moris::Cell< moris::Cell< PDV_Type > >      & aDvTypes )
        {
            // init counters for dof and dv types
            uint tMasterDofCounter = 0;
            uint tSlaveDofCounter  = 0;
            uint tMasterDvCounter  = 0;
            uint tSlaveDvCounter   = 0;

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

            // loop over the master properties
            for ( const std::shared_ptr< Property > & tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< PDV_Type >      tActiveDvTypes;

                    tProperty->get_non_unique_dof_and_dv_types(
                            tActiveDofTypes,
                            tActiveDvTypes );

                    //update dof and dv counters
                    tMasterDofCounter += tActiveDofTypes.size();
                    tMasterDvCounter  += tActiveDvTypes.size();
                }
            }

            // loop over slave properties
            for ( const std::shared_ptr< Property > & tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type lists
                    // get property non unique dof and dv type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< PDV_Type >      tActiveDvTypes;

                    tProperty->get_non_unique_dof_and_dv_types(
                            tActiveDofTypes,
                            tActiveDvTypes );

                    // update dof and dv counter
                    tSlaveDofCounter += tActiveDofTypes.size();
                    tSlaveDvCounter  += tActiveDvTypes.size();
                }
            }

            // loop over master constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< PDV_Type >      tActiveDvTypes;

                    tCM->get_non_unique_dof_and_dv_types(
                            tActiveDofTypes,
                            tActiveDvTypes );

                    // update dof and dv counters
                    tMasterDofCounter += tActiveDofTypes.size();
                    tMasterDvCounter  += tActiveDvTypes.size();

                }
            }

            // loop over slave constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< PDV_Type >      tActiveDvTypes;

                    tCM->get_non_unique_dof_and_dv_types(
                            tActiveDofTypes,
                            tActiveDvTypes );

                    // update dof and dv counters
                    tSlaveDofCounter += tActiveDofTypes.size();
                    tSlaveDvCounter  += tActiveDvTypes.size();
                }
            }

            // loop over master stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter > & tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get SP non unique dof type list
                    moris::Cell< MSI::Dof_Type >  tActiveDofTypes;
                    moris::Cell< PDV_Type >       tActiveDvTypes;

                    tSP->get_non_unique_dof_and_dv_types(
                            tActiveDofTypes,
                            tActiveDvTypes );

                    // update dof and dv counters
                    tMasterDofCounter += tActiveDofTypes.size();
                    tMasterDvCounter  += tActiveDvTypes.size();
                    tSlaveDofCounter  += tActiveDofTypes.size();
                    tSlaveDvCounter   += tActiveDvTypes.size();
                }
            }

            // reserve memory for dof and dv type lists
            aDofTypes.resize( 2 );
            aDvTypes.resize( 2 );

            aDofTypes( 0 ).reserve( tMasterDofCounter );
            aDvTypes ( 0 ).reserve( tMasterDvCounter );
            aDofTypes( 1 ).reserve( tSlaveDofCounter );
            aDvTypes ( 1 ).reserve( tSlaveDvCounter );

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

            // loop over slave dof direct dependencies
            for ( uint iDof = 0; iDof < mSlaveDofTypes.size(); iDof++ )
            {
                //populate the dof list
                aDofTypes( 1 ).append( mSlaveDofTypes( iDof )  );
            }

            // loop over slave dv direct dependencies
            for ( uint iDv = 0; iDv < mSlaveDvTypes.size(); iDv++ )
            {
                //populate the dv list
                aDvTypes( 1 ).append( mSlaveDvTypes( iDv )  );
            }

            // loop over master properties
            for ( const std::shared_ptr< Property > & tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< PDV_Type >      tActiveDvTypes;

                    tProperty->get_non_unique_dof_and_dv_types(
                            tActiveDofTypes,
                            tActiveDvTypes );

                    // populate the dof and dv lists
                    aDofTypes( 0 ).append( tActiveDofTypes );
                    aDvTypes( 0 ).append( tActiveDvTypes );
                }
            }

            // loop over slave properties
            for ( const std::shared_ptr< Property > & tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< PDV_Type >      tActiveDvTypes;

                    tProperty->get_non_unique_dof_and_dv_types(
                            tActiveDofTypes,
                            tActiveDvTypes );

                    // populate the dof and dv lists
                    aDofTypes( 1 ).append( tActiveDofTypes );
                    aDvTypes( 1 ).append( tActiveDvTypes );
                }
            }

            // loop over the master constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< PDV_Type >      tActiveDvTypes;

                    tCM->get_non_unique_dof_and_dv_types(
                            tActiveDofTypes,
                            tActiveDvTypes );

                    // populate the dof and dv lists
                    aDofTypes( 0 ).append( tActiveDofTypes );
                    aDvTypes( 0 ).append( tActiveDvTypes );
                }
            }

            // loop over the slave constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mSlaveCM )
            {
                if( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< PDV_Type >      tActiveDvTypes;

                    tCM->get_non_unique_dof_and_dv_types(
                            tActiveDofTypes,
                            tActiveDvTypes );

                    // populate the dof and dv lists
                    aDofTypes( 1 ).append( tActiveDofTypes );
                    aDvTypes( 1 ).append( tActiveDvTypes );
                }
            }

            // FIXME this is potentially problematic since it will add slave dependencies even for bulk elements
            // FIXME Ask lise about it. We could ask the set for the element type. should work for DOUBLE_SIDED.
            // FIXME Whats with time boundary
            // loop over the stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter > & tSP : mStabilizationParam )
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
                    aDvTypes ( 0 ).append( tActiveDvTypes );
                    aDofTypes( 1 ).append( tActiveDofTypes );
                    aDvTypes ( 1 ).append( tActiveDvTypes );
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
            for ( const std::shared_ptr< Property > & tProperty : mMasterProp )
            {
                if ( tProperty != nullptr )
                {
                    // get dof types for property
                    const moris::Cell< moris::Cell< MSI::Dof_Type > > & tActiveDofTypes =
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
                    moris::Cell< moris::Cell< PDV_Type > > tActiveDvTypes =
                            tProperty->get_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                        // if dof enum not in the list
                        if ( tDvCheckList( tDvTypeIndex) != 1 )
                        {
                            // put the dof type in the check list
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                        }
                    }
                }
            }

            // get dof type from master constitutive models
            for ( const std::shared_ptr< Constitutive_Model > & tCM : mMasterCM )
            {
                if ( tCM != nullptr )
                {
                    // get dof types for constitutive model
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofTypes =
                            tCM->get_global_dof_type_list();

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex) != 1 )
                        {
                            // put the dof type in the check list
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for constitutive model
                    moris::Cell< moris::Cell< PDV_Type > > tActiveDvTypes =
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
                }
            }

            // get dof type from master stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter > & tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get dof types for constitutive model
                    const moris::Cell< moris::Cell< MSI::Dof_Type > > & tActiveDofTypes =
                            tSP->get_global_dof_type_list( mtk::Master_Slave::MASTER );

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex) != 1 )
                        {
                            // put the dof type in the check list
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for constitutive model
                    const moris::Cell< moris::Cell< PDV_Type > > & tActiveDvTypes =
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

            // get dof type from master properties
            for ( const std::shared_ptr< Property > & tProperty : mSlaveProp )
            {
                if ( tProperty != nullptr )
                {
                    // get dof types for property
                    const moris::Cell< moris::Cell< MSI::Dof_Type > > & tActiveDofTypes =
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
                    moris::Cell< moris::Cell< PDV_Type > > tActiveDvTypes =
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
                }
            }

            // get dof type from slave constitutive models
            for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
            {
                if ( tCM != nullptr )
                {
                    // get dof types for constitutive model
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofTypes =
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
                    moris::Cell< moris::Cell< PDV_Type > > tActiveDvTypes =
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
                            // put the dof type in the check list
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mSlaveGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for stabilization parameter
                    moris::Cell< moris::Cell< PDV_Type > > tActiveDvTypes =
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
                moris::Cell< moris::Cell< enum MSI::Dof_Type > > tRequestedDofTypes =
                        mSet->get_secondary_dof_types();

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
            // if Jacobian evaluation
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
                case mtk::Master_Slave::MASTER :
                {
                    // loop over the dof field interpolator pointers
                    for( uint iDofFI = 0; iDofFI < mRequestedMasterGlobalDofTypes.size(); iDofFI++ )
                    {
                        // check that the field interpolator was set
                        MORIS_ASSERT(
                                this->get_field_interpolator_manager( aIsMaster )->
                                get_field_interpolators_for_type( mRequestedMasterGlobalDofTypes( iDofFI )( 0 ) ) != nullptr,
                                "IWG::check_field_interpolators - Master dof FI missing. " );
                    }

                    // loop over the dv field interpolator pointers
                    for( uint iDvFI = 0; iDvFI < mMasterGlobalDvTypes.size(); iDvFI++ )
                    {
                        // check that the field interpolator was set
                        MORIS_ASSERT(
                                this->get_field_interpolator_manager( aIsMaster )->
                                get_field_interpolators_for_type( mMasterGlobalDvTypes( iDvFI )( 0 ) ) != nullptr,
                                "IWG::check_field_interpolators - Master dv FI missing. " );
                    }
                    break;
                }
                case mtk::Master_Slave::SLAVE :
                {
                    // loop over the dof field interpolator pointers
                    for( uint iDofFI = 0; iDofFI < mRequestedSlaveGlobalDofTypes.size(); iDofFI++ )
                    {
                        // check that the field interpolator was set
                        MORIS_ASSERT(
                                this->get_field_interpolator_manager( aIsMaster )->
                                get_field_interpolators_for_type( mRequestedSlaveGlobalDofTypes( iDofFI )( 0 ) ) != nullptr,
                                "IWG::check_dof_field_interpolators - Slave dof FI missing. " );
                    }

                    // loop over the dv field interpolator pointers
                    for( uint iDvFI = 0; iDvFI < mSlaveGlobalDvTypes.size(); iDvFI++ )
                    {
                        // check that the field interpolator was set
                        MORIS_ASSERT(
                                this->get_field_interpolator_manager( aIsMaster )->
                                get_field_interpolators_for_type( mSlaveGlobalDvTypes( iDvFI )( 0 ) ) != nullptr,
                                "IWG::check_field_interpolators - Slave dv FI missing. " );
                    }
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "IWG::check_field_interpolators - can only be master or slave." );
                }
            }
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Cell< MSI::Dof_Type > > & IWG::get_global_dof_type_list(
                mtk::Master_Slave aIsMaster )
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
                case mtk::Master_Slave::MASTER :
                {
                    // return master global dof type list
                    return mMasterGlobalDofTypes;
                    break;
                }
                // if slave
                case mtk::Master_Slave::SLAVE :
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

        const moris::Cell< moris::Cell< PDV_Type > > & IWG::get_global_dv_type_list(
                mtk::Master_Slave aIsMaster )
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
                case mtk::Master_Slave::MASTER :
                {
                    // return master global dof type list
                    return mMasterGlobalDvTypes;
                }
                // if slave
                case mtk::Master_Slave::SLAVE :
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

        void IWG::select_jacobian_FD(
                real                aWStar,
                real                aPerturbation,
                fem::FDScheme_Type  aFDSchemeType )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get master index for residual dof type, indices for assembly
            sint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
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
            for( uint iFI = 0; iFI < tMasterNumDofTypes; iFI++ )
            {
                // init dof counter
                uint tDofCounter = 0;

                // get the dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iFI );

                // get the index for the dof type
                sint tMasterDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tMasterDepDofIndex, 0 );

                // get field interpolator for dependency dof type
                Field_Interpolator * tFI =
                        mMasterFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // loop over the coefficient column
                for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficient row
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++  )
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
                            // add unperturbed residual contribution to dRdp
                            mSet->get_jacobian()(
                                    { tMasterResStartIndex, tMasterResStopIndex },
                                    { tMasterDepStartIndex + tDofCounter, tMasterDepStartIndex + tDofCounter } ) +=
                                            tFDScheme( 1 )( 0 ) * tResidual /
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
                            tFI->reset_eval_flags(); // not useful

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
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );
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

        void IWG::select_jacobian_FD_double(
                real                aWStar,
                real                aPerturbation,
                fem::FDScheme_Type  aFDSchemeType )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get master index for residual dof type, indices for assembly
            sint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            sint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
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
            for( uint iFI = 0; iFI < tMasterNumDofTypes; iFI++ )
            {
                // init dof counter
                uint tDofCounter = 0;

                // get the dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iFI );

                // get the index for the dof type
                sint tMasterDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tMasterDepDofIndex, 0 );

                // get field interpolator for dependency dof type
                Field_Interpolator * tFI =
                        mMasterFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // loop over the coefficient column
                for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficient row
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++  )
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
                            // add unperturbed master residual contribution to dRdp
                            mSet->get_jacobian()(
                                    { tMasterResStartIndex, tMasterResStopIndex },
                                    { tMasterDepStartIndex + tDofCounter, tMasterDepStartIndex + tDofCounter } ) +=
                                            tFDScheme( 1 )( 0 ) * tMasterResidual /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // add unperturbed slave residual contribution to dRdp
                            mSet->get_jacobian()(
                                    { tSlaveResStartIndex, tSlaveResStopIndex },
                                    { tMasterDepStartIndex + tDofCounter, tMasterDepStartIndex + tDofCounter } ) +=
                                            tFDScheme( 1 )( 0 ) * tSlaveResidual /
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
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // assemble slave part of the jacobian
                            mSet->get_jacobian()(
                                    { tSlaveResStartIndex, tSlaveResStopIndex },
                                    { tMasterDepStartIndex + tDofCounter, tMasterDepStartIndex + tDofCounter } ) +=
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } ) /
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
            uint tSlaveNumDofTypes = mRequestedSlaveGlobalDofTypes.size();

            // loop over the IWG dof types
            for( uint iFI = 0; iFI < tSlaveNumDofTypes; iFI++ )
            {
                // init dof counter
                uint tDofCounter = 0;

                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iFI );

                // get the index for the dof type
                sint tSlaveDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tSlaveDepDofIndex, 0 );

                // get field interpolator for dependency dof type
                Field_Interpolator * tFI =
                        mSlaveFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

                // get number of master FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // loop over the coefficient column
                for( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficient row
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++  )
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
                            // add unperturbed master residual contribution to dRdp
                            mSet->get_jacobian()(
                                    { tMasterResStartIndex, tMasterResStopIndex },
                                    { tSlaveDepStartIndex + tDofCounter, tSlaveDepStartIndex + tDofCounter } ) +=
                                            tFDScheme( 1 )( 0 ) * tMasterResidual /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // add unperturbed slave residual contribution to dRdp
                            mSet->get_jacobian()(
                                    { tSlaveResStartIndex, tSlaveResStopIndex },
                                    { tSlaveDepStartIndex + tDofCounter, tSlaveDepStartIndex + tDofCounter } ) +=
                                            tFDScheme( 1 )( 0 ) * tSlaveResidual /
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
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // assemble the jacobian
                            mSet->get_jacobian()(
                                    { tSlaveResStartIndex, tSlaveResStopIndex },
                                    { tSlaveDepStartIndex + tDofCounter, tSlaveDepStartIndex + tDofCounter } ) +=
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );
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

        bool IWG::check_jacobian(
                real               aPerturbation,
                real               aEpsilon,
                real               aWStar,
                Matrix< DDRMat > & aJacobian,
                Matrix< DDRMat > & aJacobianFD,
                bool               aErrorPrint )
        {
            // get residual dof type index in set, start and end indices for residual dof type
            uint tMasterDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartRow = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResEndRow   = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            uint tSlaveDofIndex;
            uint tSlaveResStartRow;
            uint tSlaveResEndRow;
            uint tSlaveNumRows = 0;
            if( mSlaveGlobalDofTypes.size() > 0 )
            {
                tSlaveDofIndex     = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
                tSlaveResStartRow  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
                tSlaveResEndRow    = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );
                tSlaveNumRows      = tSlaveResEndRow  - tSlaveResStartRow  + 1;
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
                    mSet->get_jacobian()( { tMasterResStartRow, tMasterResEndRow }, { 0, tNumCols - 1 });
            if( tSlaveNumRows > 0 )
            {
                aJacobian( { tMasterNumRows, tMasterNumRows + tSlaveNumRows - 1 }, { 0, tNumCols - 1 } ) =
                        mSet->get_jacobian()( { tSlaveResStartRow, tSlaveResEndRow }, { 0, tNumCols - 1 });
            }

            // reset the jacobian
            mSet->get_jacobian().fill( 0.0 );

            // compute jacobian by FD
              this->compute_jacobian_FD( aWStar, aPerturbation );


            // get the computed jacobian
            aJacobianFD( { 0, tMasterNumRows - 1 }, { 0, tNumCols - 1 } ) =
                    mSet->get_jacobian()( { tMasterResStartRow, tMasterResEndRow }, { 0, tNumCols - 1 });
            if( tSlaveNumRows > 0 )
            {
                aJacobianFD( { tMasterNumRows, tMasterNumRows + tSlaveNumRows - 1 }, { 0, tNumCols - 1 } ) =
                        mSet->get_jacobian()( { tSlaveResStartRow, tSlaveResEndRow }, { 0, tNumCols - 1 });
            }

            // check that matrices to compare have same size
            MORIS_ERROR(
                    ( aJacobian.n_rows() == aJacobianFD.n_rows() ) &&
                    ( aJacobian.n_cols() == aJacobianFD.n_cols() ),
                    "IWG::check_jacobian - matrices to check do not share same dimensions." );

            //define a boolean for check
            bool tCheckJacobian = true;

            // define a real for absolute difference
            real tAbsolute = 0.0;

            // define a real for relative difference
            real tRelative = 0.0;

            for( uint iiJac = 0; iiJac < aJacobian.n_rows(); iiJac++ )
            {
                for( uint jjJac = 0; jjJac < aJacobian.n_cols(); jjJac++ )
                {
                    // get absolute difference
                    tAbsolute = std::abs( aJacobian( iiJac, jjJac ) - aJacobianFD( iiJac, jjJac ) );

                    // get relative difference
                    tRelative = std::abs( ( aJacobianFD( iiJac, jjJac ) - aJacobian( iiJac, jjJac ) ) / aJacobianFD( iiJac, jjJac ) );

                    // update check value
                    tCheckJacobian = tCheckJacobian && ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) );

                    // debug print
                    if( ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) ) == false )
                    {
                        if( aErrorPrint )
                        {
                            std::cout<<"iiJac "<<iiJac<<" - jjJac "<<jjJac<<"\n"<<std::flush;
                            std::cout<<"aJacobian( iiJac, jjJac )   "<<std::setprecision( 12 )<<aJacobian( iiJac, jjJac )<<"\n"<<std::flush;
                            std::cout<<"aJacobianFD( iiJac, jjJac ) "<<std::setprecision( 12 )<<aJacobianFD( iiJac, jjJac )<<"\n"<<std::flush;
                            std::cout<<"Absolute difference "<<tAbsolute<<"\n"<<std::flush;
                            std::cout<<"Relative difference "<<tRelative<<"\n"<<std::flush;
                        }
                    }
                }
            }

            // return bool
            return tCheckJacobian;
        }

        //------------------------------------------------------------------------------

        void IWG::select_dRdp_FD_geometry_bulk(
                moris::real                         aWStar,
                moris::real                         aPerturbation,
                fem::FDScheme_Type                  aFDSchemeType,
                Matrix< DDSMat >                  & aGeoLocalAssembly,
                moris::Cell< Matrix< IndexMat > > & aVertexIndices )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the GI for the IG element considered
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();
            Geometry_Interpolator * tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();

            // get the residual dof type index in the set
            uint tResDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
            uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

            // reset, evaluate and store the residual for unperturbed case
            mSet->get_residual()( 0 ).fill( 0.0 );
            this->compute_residual( aWStar );
            Matrix< DDRMat > tResidual = mSet->get_residual()( 0 )(
                    { tResDofAssemblyStart, tResDofAssemblyStop },
                    { 0, 0 } );

            // get number of master GI bases and space dimensions
            uint tDerNumBases      = tIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

            // coefficients for dv type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tIGGI->get_space_coeff();
            Matrix< DDRMat > tParamCoeff = tIGGI->get_space_param_coeff();
            Matrix< DDRMat > tEvaluationPoint;
            tIGGI->get_space_time( tEvaluationPoint );
            real tGPWeight = aWStar / tIGGI->det_J();

            // IP element max/min
            Matrix< DDRMat > tMaxIP = max( tIPGI->get_space_coeff() );
            Matrix< DDRMat > tMinIP = min( tIPGI->get_space_coeff() );

            // init perturbation
            real tDeltaH = 0.0;

            // init FD scheme
            moris::Cell< moris::Cell< real > > tFDScheme;

            // loop over the spatial directions
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
            {
                // loop over the IG nodes
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    // if pdv is active
                    if( tPdvAssemblyIndex != -1 )
                    {
                        // compute the perturbation value
                        tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // check point location
                        fem::FDScheme_Type tUsedFDScheme = aFDSchemeType;
                        if( tCoeff( iCoeffRow, iCoeffCol ) + tDeltaH > tMaxIP( iCoeffCol ) )
                        {
                            tUsedFDScheme = fem::FDScheme_Type::POINT_1_BACKWARD;
                        }
                        else if( tCoeff( iCoeffRow, iCoeffCol ) - tDeltaH < tMinIP( iCoeffCol ) )
                        {
                            tUsedFDScheme = fem::FDScheme_Type::POINT_1_FORWARD;
                        }
                        fd_scheme( tUsedFDScheme, tFDScheme );
                        uint tNumPoints = tFDScheme( 0 ).size();

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward add unperturbed contribution
                        if( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                                ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed residual contribution to dRdp
                            mSet->get_drdpgeo()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvAssemblyIndex,    tPdvAssemblyIndex } ) +=
                                            tFDScheme( 1 )( 0 ) * tResidual /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over point of FD scheme
                        for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
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
                            Matrix< DDRMat > tParamCoeffPert = tParamCoeff;
                            tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();
                            tIGGI->set_space_param_coeff( tParamCoeffPert );

                            // set evaluation point for interpolators (FIs and GIs)
                            mSet->get_field_interpolator_manager()->
                                    set_space_time_from_local_IG_point( tEvaluationPoint );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset and evaluate the residual plus
                            mSet->get_residual()( 0 ).fill( 0.0 );
                            real tWStarPert = tGPWeight * tIGGI->det_J();
                            this->compute_residual( tWStarPert );

                            // evaluate dRdpGeo
                            mSet->get_drdpgeo()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvAssemblyIndex,    tPdvAssemblyIndex } ) +=
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                }
            }
            // reset the coefficients values
            tIGGI->set_space_coeff( tCoeff );
            tIGGI->set_space_param_coeff( tParamCoeff );
            mSet->get_field_interpolator_manager()->
                    set_space_time_from_local_IG_point( tEvaluationPoint );

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ) ,
                    "IWG::compute_dRdp_FD_geometry - dRdp contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG::select_dRdp_FD_geometry_sideset(
                moris::real                         aWStar,
                moris::real                         aPerturbation,
                fem::FDScheme_Type                  aFDSchemeType,
                Matrix< DDSMat >                  & aGeoLocalAssembly,
                moris::Cell< Matrix< IndexMat > > & aVertexIndices )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the GI for the IG element considered
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();
            Geometry_Interpolator * tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();

            // get the residual dof type index in the set
            uint tResDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
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

            // IP element max/min
            Matrix< DDRMat > tMaxIP = max( tIPGI->get_space_coeff() );
            Matrix< DDRMat > tMinIP = min( tIPGI->get_space_coeff() );

            // init perturbation
            real tDeltaH = 0.0;

            // init FD scheme
            moris::Cell< moris::Cell< real > > tFDScheme;

            // get number of master GI bases and space dimensions
            uint tDerNumBases      = tIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

            // loop over the spatial directions
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
            {
                // loop over the IG nodes
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    // if pdv is active
                    if( tPdvAssemblyIndex != -1 )
                    {
                        // compute the perturbation value
                        tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // check point location
                        fem::FDScheme_Type tUsedFDScheme = aFDSchemeType;
                        if( tCoeff( iCoeffRow, iCoeffCol ) + tDeltaH > tMaxIP( iCoeffCol ) )
                        {
                            tUsedFDScheme = fem::FDScheme_Type::POINT_1_BACKWARD;
                        }
                        else if( tCoeff( iCoeffRow, iCoeffCol ) - tDeltaH < tMinIP( iCoeffCol ) )
                        {
                            tUsedFDScheme = fem::FDScheme_Type::POINT_1_FORWARD;
                        }
                        fd_scheme( tUsedFDScheme, tFDScheme );
                        uint tNumPoints = tFDScheme( 0 ).size();

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward add unperturbed contribution
                        if( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                                ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed residual contribution to dRdp
                            mSet->get_drdpgeo()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvAssemblyIndex,    tPdvAssemblyIndex } ) +=
                                            tFDScheme( 1 )( 0 ) * tResidual /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over point of FD scheme
                        for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
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
                            Matrix< DDRMat > tParamCoeffPert = tParamCoeff;
                            tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();
                            tIGGI->set_space_param_coeff( tParamCoeffPert );

                            // set evaluation point for interpolators (FIs and GIs)
                            mSet->get_field_interpolator_manager()->
                                    set_space_time_from_local_IG_point( tEvaluationPoint );

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
                                    { tPdvAssemblyIndex,    tPdvAssemblyIndex } ) +=
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                }
            }

            // reset xyz values
            tIGGI->set_space_coeff( tCoeff );

            // reset local coordinates values
            tIGGI->set_space_param_coeff( tParamCoeff );

            // reset evaluation point
            mSet->get_field_interpolator_manager()->
                    set_space_time_from_local_IG_point( tEvaluationPoint );

            // reset normal
            this->set_normal( tNormal );

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ) ,
                    "IWG::compute_dRdp_FD_geometry - dRdp contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG::select_dRdp_FD_geometry_time_sideset(
                moris::real                         aWStar,
                moris::real                         aPerturbation,
                fem::FDScheme_Type                  aFDSchemeType,
                Matrix< DDSMat >                  & aGeoLocalAssembly,
                moris::Cell< Matrix< IndexMat > > & aVertexIndices )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the GI for the IG element considered
            Geometry_Interpolator * tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();
            Geometry_Interpolator * tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();
            Geometry_Interpolator * tIGGIPrevious =
                    mSet->get_field_interpolator_manager_previous_time()->get_IG_geometry_interpolator();

            // get the residual dof type index in the set
            uint tResDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
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

            // IP element max/min
            Matrix< DDRMat > tMaxIP = max( tIPGI->get_space_coeff() );
            Matrix< DDRMat > tMinIP = min( tIPGI->get_space_coeff() );

            // init perturbation
            real tDeltaH = 0.0;

            // init FD scheme
            moris::Cell< moris::Cell< real > > tFDScheme;

            // get number of master GI bases and space dimensions
            uint tDerNumBases      = tIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

            // loop over the spatial directions
            for( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
            {
                // loop over the IG nodes
                for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    // if pdv is active
                    if( tPdvAssemblyIndex != -1 )
                    {
                        // compute the perturbation value
                        tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // check point location
                        fem::FDScheme_Type tUsedFDScheme = aFDSchemeType;
                        if( tCoeff( iCoeffRow, iCoeffCol ) + tDeltaH > tMaxIP( iCoeffCol ) )
                        {
                            tUsedFDScheme = fem::FDScheme_Type::POINT_1_BACKWARD;
                        }
                        else if( tCoeff( iCoeffRow, iCoeffCol ) - tDeltaH < tMinIP( iCoeffCol ) )
                        {
                            tUsedFDScheme = fem::FDScheme_Type::POINT_1_FORWARD;
                        }
                        fd_scheme( tUsedFDScheme, tFDScheme );
                        uint tNumPoints = tFDScheme( 0 ).size();

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward add unperturbed contribution
                        if( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                                ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed residual contribution to dRdp
                            mSet->get_drdpgeo()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvAssemblyIndex,    tPdvAssemblyIndex } ) +=
                                            tFDScheme( 1 )( 0 ) * tResidual /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over point of FD scheme
                        for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
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
                            Matrix< DDRMat > tParamCoeffPert = tParamCoeff;
                            tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();
                            tIGGI->set_space_param_coeff( tParamCoeffPert );
                            tIGGIPrevious->set_space_param_coeff( tParamCoeffPert );

                            // set evaluation point for interpolators (FIs and GIs)
                            mSet->get_field_interpolator_manager()->
                                    set_space_time_from_local_IG_point( tEvaluationPoint );
                            mSet->get_field_interpolator_manager_previous_time()->
                                    set_space_time_from_local_IG_point( tEvaluationPoint );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset and evaluate the residual plus
                            mSet->get_residual()( 0 ).fill( 0.0 );
                            real tWStarPert = tGPWeight * tIGGI->det_J();
                            this->compute_residual( tWStarPert );

                            // evaluate dRdpGeo
                            mSet->get_drdpgeo()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvAssemblyIndex,    tPdvAssemblyIndex } ) +=
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );
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
            mSet->get_field_interpolator_manager()->
                    set_space_time_from_local_IG_point( tEvaluationPoint );
            mSet->get_field_interpolator_manager_previous_time()->
                    set_space_time_from_local_IG_point( tEvaluationPoint );

            // reset the coefficients values
            tIGGI->set_space_coeff( tCoeff );
            tIGGIPrevious->set_space_coeff( tCoeff );

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ) ,
                    "IWG::compute_dRdp_FD_geometry - dRdp contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG::select_dRdp_FD_geometry_double(
                moris::real                         aWStar,
                moris::real                         aPerturbation,
                fem::FDScheme_Type                  aFDSchemeType,
                Matrix< DDSMat >                  & aGeoLocalAssembly,
                moris::Cell< Matrix< IndexMat > > & aVertexIndices )
        {
            // unpack vertex indices
            Matrix< IndexMat > & aMasterVertexIndices = aVertexIndices( 0 );
            Matrix< IndexMat > & aSlaveVertexIndices  = aVertexIndices( 1 );

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
            Geometry_Interpolator * tMasterIGGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->get_IG_geometry_interpolator();
            Geometry_Interpolator * tMasterIPGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->get_IP_geometry_interpolator();

            // get the slave GI for the IG and IP element considered
            Geometry_Interpolator * tSlaveIGGI =
                    mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->get_IG_geometry_interpolator();

            // IP element max/min
            Matrix< DDRMat > tMasterMaxIP = max( tMasterIPGI->get_space_coeff() );
            Matrix< DDRMat > tMasterMinIP = min( tMasterIPGI->get_space_coeff() );

            // get the master residual dof type index in the set
            uint tMasterResDofIndex = mSet->get_dof_index_for_type(
                    mResidualDofType( 0 ),
                    mtk::Master_Slave::MASTER );
            uint tMasterResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tMasterResDofIndex )( 0, 0 );
            uint tMasterResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tMasterResDofIndex )( 0, 1 );

            // get the slave residual dof type index in the set
            uint tSlaveResDofIndex = mSet->get_dof_index_for_type(
                    mResidualDofType( 0 ),
                    mtk::Master_Slave::SLAVE );
            uint tSlaveResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tSlaveResDofIndex )( 0, 0 );
            uint tSlaveResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tSlaveResDofIndex )( 0, 1 );

            // init perturbation
            real tDeltaH = 0.0;

            // get GP weight
            real tGPWeight = aWStar / tMasterIGGI->det_J();

            // get master coeff
            Matrix< DDRMat > tMasterCoeff = tMasterIGGI->get_space_coeff();
            Matrix< DDRMat > tMasterParamCoeff = tMasterIGGI->get_space_param_coeff();
            Matrix< DDRMat > tMasterEvaluationPoint;
            tMasterIGGI->get_space_time( tMasterEvaluationPoint );
            Matrix< DDRMat > tMasterNormal;
            tMasterIGGI->get_normal( tMasterNormal );

            // coefficients for dv type wrt which derivative is computed
            Matrix< DDRMat > tSlaveCoeff = tSlaveIGGI->get_space_coeff();
            Matrix< DDRMat > tSlaveParamCoeff = tSlaveIGGI->get_space_param_coeff();
            Matrix< DDRMat > tSlaveEvaluationPoint;
            tSlaveIGGI->get_space_time( tSlaveEvaluationPoint );

            // reset, evaluate and store the residual for unperturbed case
            mSet->get_residual()( 0 ).fill( 0.0 );
            this->compute_residual( aWStar );
            Matrix< DDRMat > tMasterResidual =
                    mSet->get_residual()( 0 )( { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop }, { 0, 0 } );
            Matrix< DDRMat > tSlaveResidual =
                    mSet->get_residual()( 0 )( { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop }, { 0, 0 } );

            // get number of master GI bases and space dimensions
            uint tDerNumBases      = tMasterIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tMasterIPGI->get_number_of_space_dimensions();

            // init FD scheme
            moris::Cell< moris::Cell< real > > tFDScheme;

            // loop over the IG nodes
            for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // find the node on the slave side
                sint tSlaveNodeLocalIndex = -1;
                for ( uint iNode = 0; iNode < tDerNumBases; iNode++ )
                {
                    if( aMasterVertexIndices( iCoeffRow ) == aSlaveVertexIndices( iNode ) )
                    {
                        tSlaveNodeLocalIndex = iNode;
                        break;
                    }
                }
                MORIS_ERROR( tSlaveNodeLocalIndex != -1, "IWG::compute_dRdp_FD_geometry_double - slave index not found.");

                // loop over the spatial directions
                for( uint iCoeffCol = 0; iCoeffCol< tDerNumDimensions; iCoeffCol++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    if ( tPdvAssemblyIndex != -1 )
                    {
                        // compute the perturbation value
                        tDeltaH = aPerturbation * tMasterCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // check point location
                        fem::FDScheme_Type tUsedFDScheme = aFDSchemeType;
                        if( tMasterCoeff( iCoeffRow, iCoeffCol ) + tDeltaH > tMasterMaxIP( iCoeffCol ) )
                        {
                            tUsedFDScheme = fem::FDScheme_Type::POINT_1_BACKWARD;
                        }
                        else if( tMasterCoeff( iCoeffRow, iCoeffCol ) - tDeltaH < tMasterMinIP( iCoeffCol ) )
                        {
                            tUsedFDScheme = fem::FDScheme_Type::POINT_1_FORWARD;
                        }
                        fd_scheme( tUsedFDScheme, tFDScheme );
                        uint tNumPoints = tFDScheme( 0 ).size();

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward add unperturbed contribution
                        if( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                                ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed master residual contribution to dRdp
                            mSet->get_drdpgeo()(
                                    { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                                    { tPdvAssemblyIndex,          tPdvAssemblyIndex } ) +=
                                            tFDScheme( 1 )( 0 ) * tMasterResidual /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // add unperturbed slave residual contribution to dRdp
                            mSet->get_drdpgeo()(
                                    { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                    { tPdvAssemblyIndex,          tPdvAssemblyIndex } ) +=
                                            tFDScheme( 1 )( 0 ) * tSlaveResidual /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over point of FD scheme
                        for ( uint iPoint = tStartPoint; iPoint < tNumPoints; iPoint++ )
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
                            Matrix< DDRMat > tMasterParamCoeffPert = tMasterParamCoeff;
                            tMasterParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();
                            Matrix< DDRMat > tSlaveParamCoeffPert = tSlaveParamCoeff;
                            tSlaveParamCoeffPert.get_row( tSlaveNodeLocalIndex ) = tXiCoords.matrix_data();

                            tMasterIGGI->set_space_param_coeff( tMasterParamCoeffPert );
                            tSlaveIGGI->set_space_param_coeff( tSlaveParamCoeffPert );

                            // set evaluation point for interpolators (FIs and GIs)
                            mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->
                                    set_space_time_from_local_IG_point( tMasterEvaluationPoint );
                            mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->
                                    set_space_time_from_local_IG_point( tSlaveEvaluationPoint );

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
                                    { tPdvAssemblyIndex,          tPdvAssemblyIndex } ) +=
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // evaluate dSlaveRdpGeo
                            mSet->get_drdpgeo()(
                                    { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                    { tPdvAssemblyIndex,         tPdvAssemblyIndex } ) +=
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                }
            }
            // reset the coefficients values
            tMasterIGGI->set_space_coeff( tMasterCoeff );
            tMasterIGGI->set_space_param_coeff( tMasterParamCoeff );
            tSlaveIGGI->set_space_coeff( tSlaveCoeff );
            tSlaveIGGI->set_space_param_coeff( tSlaveParamCoeff );
            mSet->get_field_interpolator_manager( mtk::Master_Slave::MASTER )->
                    set_space_time_from_local_IG_point( tMasterEvaluationPoint );
            mSet->get_field_interpolator_manager( mtk::Master_Slave::SLAVE )->
                    set_space_time_from_local_IG_point( tSlaveEvaluationPoint );
            this->set_normal( tMasterNormal );

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ),
                    "IWG::compute_dRdp_FD_geometry_double - dRdp contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG::select_dRdp_FD_material(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get the requested ip pdv types
            moris::Cell< moris::Cell< PDV_Type > > tRequestedPdvTypes;
            mSet->get_ip_dv_types_for_set( tRequestedPdvTypes );

            // get number of requested dv types
            uint tNumDvType = tRequestedPdvTypes.size();

            // get the residual dof type index in the set
            uint tResDofIndex = mSet->get_dof_index_for_type(
                    mResidualDofType( 0 ),
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
            for( uint iFI = 0; iFI < tNumDvType; iFI++ )
            {
                // get dv index
                sint tDvDepIndex = mSet->get_dv_index_for_type(
                        tRequestedPdvTypes( iFI )( 0 ),
                        mtk::Master_Slave::MASTER );

                // get the FI for the dv type
                Field_Interpolator * tFI =
                        mMasterFIManager->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) );

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
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++  )
                    {
                        // compute the perturbation absolute value
                        tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // get mat pdv index
                        uint tPdvIndex = mSet->get_mat_pdv_assembly_map()( tDvDepIndex )( 0, 0 ) + tCoeffCounter;

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward FD, add unperturbed residual contribution
                        if( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                                ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed master residual contribution to dRdp
                            mSet->get_drdpmat()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvIndex,                  tPdvIndex } ) +=
                                            tFDScheme( 1 )( 0 ) * tResidual /
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

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the residual
                            mSet->get_residual()( 0 ).fill( 0.0 );

                            // compute the residual
                            this->compute_residual( aWStar );

                            // evaluate dRdpMat
                            mSet->get_drdpmat()(
                                    { tResDofAssemblyStart, tResDofAssemblyStop },
                                    { tPdvIndex,            tPdvIndex } ) +=
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );
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
            MORIS_ASSERT( isfinite( mSet->get_drdpmat() ) ,
                    "IWG::compute_dRdp_FD_material - dRdp contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG::select_dRdp_FD_material_double(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the FD scheme info
            moris::Cell< moris::Cell< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumPoints = tFDScheme( 0 ).size();

            // get the requested ip pdv types
            moris::Cell< moris::Cell< PDV_Type > > tRequestedPdvTypes;
            mSet->get_ip_dv_types_for_set( tRequestedPdvTypes );

            // get number of requested dv types
            uint tNumDvType = tRequestedPdvTypes.size();

            // get the master residual dof type index in the set
            uint tMasterResDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tMasterResDofIndex )( 0, 0 );
            uint tMasterResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tMasterResDofIndex )( 0, 1 );

            // get the slave residual dof type index in the set
            uint tSlaveResDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tSlaveResDofIndex )( 0, 0 );
            uint tSlaveResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tSlaveResDofIndex )( 0, 1 );

            // reset, evaluate and store the residual for unperturbed case
            mSet->get_residual()( 0 ).fill( 0.0 );
            this->compute_residual( aWStar );
            Matrix< DDRMat > tMasterResidual =
                    mSet->get_residual()( 0 )(
                            { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                            { 0, 0 } );
            Matrix< DDRMat > tSlaveResidual =
                    mSet->get_residual()( 0 )(
                            { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                            { 0, 0 } );

            // init perturbation
            real tDeltaH = 0.0;

            // loop over the master dv types associated with a FI
            for( uint iFI = 0; iFI < tNumDvType; iFI++ )
            {
                // get dv index
                sint tDvDepIndex = mSet->get_dv_index_for_type(
                        tRequestedPdvTypes( iFI )( 0 ),
                        mtk::Master_Slave::MASTER );

                // get the FI for the dv type
                Field_Interpolator * tFI =
                        mMasterFIManager->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) );

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
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++  )
                    {
                        // compute the perturbation absolute value
                        tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // get mat pdv index
                        uint tPdvIndex = mSet->get_mat_pdv_assembly_map()( tDvDepIndex )( 0, 0 ) + tCoeffCounter;

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward FD, add unperturbed residual contribution
                        if( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                                ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed master residual contribution to dRdp
                            mSet->get_drdpmat()(
                                    { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                                    { tPdvIndex,                  tPdvIndex } ) +=
                                            tFDScheme( 1 )( 0 ) * tMasterResidual /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // add unperturbed master residual contribution to dRdp
                            mSet->get_drdpmat()(
                                    { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                    { tPdvIndex,                  tPdvIndex } ) +=
                                            tFDScheme( 1 )( 0 ) * tSlaveResidual /
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

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the residual
                            mSet->get_residual()( 0 ).fill( 0.0 );

                            // compute the residual
                            this->compute_residual( aWStar );

                            // assemble dRMasterdpMat
                            mSet->get_drdpmat()(
                                    { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                                    { tPdvIndex,                  tPdvIndex } ) +=
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // assemble dRSlavedpMat
                            mSet->get_drdpmat()(
                                    { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                    { tPdvIndex,                 tPdvIndex } ) +=
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update coefficient counter
                        tCoeffCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // loop over the slave dv types associated with a FI
            for( uint iFI = 0; iFI < tNumDvType; iFI++ )
            {
                // get dv index
                sint tDvDepIndex = mSet->get_dv_index_for_type(
                        tRequestedPdvTypes( iFI )( 0 ),
                        mtk::Master_Slave::SLAVE );

                // get the FI for the dv type
                Field_Interpolator * tFI =
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
                Matrix< DDRMat > tSlaveResidual =
                        mSet->get_residual()( 0 )(
                                { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                { 0, 0 } );

                // init coeff counter
                uint tCoeffCounter = 0;

                // loop over the coefficient column
                for( uint iCoeffCol = 0; iCoeffCol< tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficient row
                    for( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++  )
                    {
                        // compute the perturbation absolute value
                        tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // get mat pdv index
                        uint tPdvIndex = mSet->get_mat_pdv_assembly_map()( tDvDepIndex )( 0, 0 ) + tCoeffCounter;

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward FD, add unperturbed residual contribution
                        if( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||
                                ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed master residual contribution to dRdp
                            mSet->get_drdpmat()(
                                    { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                                    { tPdvIndex,                  tPdvIndex } ) +=
                                            tFDScheme( 1 )( 0 ) * tMasterResidual /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // add unperturbed master residual contribution to dRdp
                            mSet->get_drdpmat()(
                                    { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                    { tPdvIndex,                  tPdvIndex } ) +=
                                            tFDScheme( 1 )( 0 ) * tSlaveResidual /
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

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the residual
                            mSet->get_residual()( 0 ).fill( 0.0 );

                            // compute the residual
                            this->compute_residual( aWStar );

                            // assemble dRMasterdpMat
                            mSet->get_drdpmat()(
                                    { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop },
                                    { tPdvIndex,                  tPdvIndex } ) +=
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tMasterResDofAssemblyStart, tMasterResDofAssemblyStop }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // assemble dRSlavedpMat
                            mSet->get_drdpmat()(
                                    { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop },
                                    { tPdvIndex,                 tPdvIndex } ) +=
                                            tFDScheme( 1 )( iPoint ) *
                                            mSet->get_residual()( 0 )( { tSlaveResDofAssemblyStart, tSlaveResDofAssemblyStop }, { 0, 0 } ) /
                                            ( tFDScheme( 2 )( 0 ) * tDeltaH );
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
            MORIS_ASSERT( isfinite( mSet->get_drdpmat() ) ,
                    "IWG::compute_dRdp_FD_material - dRdp contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

    }   // end fem namespace
}   // end moris namespace



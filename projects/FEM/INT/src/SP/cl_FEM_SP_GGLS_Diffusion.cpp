/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_GGLS_Diffusion.cpp
 *
 */

#include "cl_FEM_SP_GGLS_Diffusion.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_FEM_CM_Phase_State_Functions.hpp"
//LINALG/src
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        SP_GGLS_Diffusion::SP_GGLS_Diffusion()
        {
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Conductivity" ]       = static_cast< uint >( Property_Type::CONDUCTIVITY );
            mPropertyMap[ "Density" ]            = static_cast< uint >( Property_Type::DENSITY );
            mPropertyMap[ "HeatCapacity" ]       = static_cast< uint >( Property_Type::HEAT_CAPACITY );
            mPropertyMap[ "LatentHeat" ]         = static_cast< uint >( Property_Type::LATENT_HEAT );
            mPropertyMap[ "PCTemp" ]             = static_cast< uint >( Property_Type::PC_TEMP );
            mPropertyMap[ "PhaseStateFunction" ] = static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION );
            mPropertyMap[ "PhaseChangeConst" ]   = static_cast< uint >( Property_Type::PHASE_CHANGE_CONST );
        }

        //------------------------------------------------------------------------------

        void SP_GGLS_Diffusion::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                moris::Cell< std::string >                  & aDofStrings,
                mtk::Master_Slave                             aIsMaster )
        {

            // switch on master slave
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    // set dof type list
                    mMasterDofTypes = aDofTypes;

                    // loop on dof type
                    for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        std::string tDofString = aDofStrings( iDof );

                        // get dof type
                        MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if velocity
                        if( tDofString == "Temperature" )
                        {
                            mMasterDofTemp = tDofType;
                        }
                        else
                        {
                            // error unknown dof string
                            MORIS_ERROR( false ,
                                    "SP_GGLS_Diffusion::set_dof_type_list - Unknown aDofString : %s \n",
                                    tDofString.c_str() );
                        }
                    }
                    break;
                }

                case mtk::Master_Slave::SLAVE :
                {
                    // set dof type list
                    mSlaveDofTypes = aDofTypes;
                    break;
                }

                default:
                    MORIS_ERROR( false, "SP_GGLS_Diffusion::set_dof_type_list - unknown master slave type." );
            }
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::tuple<
        fem::Measure_Type,
        mtk::Primary_Void,
        mtk::Master_Slave > > SP_GGLS_Diffusion::get_cluster_measure_tuple_list()
        {
            return { mElementSizeTuple };
        }

        //------------------------------------------------------------------------------

        void SP_GGLS_Diffusion::eval_SP()
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // get the property values
            real tConductivity = mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->val()( 0 );
            real tDensity      = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );
            real tHeatCapacity = mMasterProp( static_cast< uint >( Property_Type::HEAT_CAPACITY ) )->val()( 0 );

            std::shared_ptr< Property > & tPropLatentHeat = mMasterProp( static_cast< uint >( Property_Type::LATENT_HEAT ) );

            // time step size
            real tDeltat = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

            // threshold product of conductivity and time step
            real tProdCondDeltat = std::max(tConductivity * tDeltat, mEpsilon);

            // get alpha
            real tAlpha = 0.0;

            // if there is a phase change
            if (tPropLatentHeat != nullptr)
            {
                // get values of properties
                real tLatentHeat = tPropLatentHeat->val()( 0 );
                real tMeltTemp   = mMasterProp( static_cast< uint >( Property_Type::PC_TEMP ) )->val()( 0 );
                real tPCconst    = mMasterProp( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) )->val()( 0 );
                real tPSfunc     = mMasterProp( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) )->val()( 0 );

                // get the dof type FI
                Field_Interpolator * tFITemp =
                        mMasterFIManager->get_field_interpolators_for_type( mMasterDofTemp );

                // get phase state function
                moris::real tdfdT = eval_dFdTemp( tMeltTemp, tPCconst, tPSfunc, tFITemp );

                tAlpha = ( tDensity * (tHeatCapacity + tLatentHeat * tdfdT) * std::pow(tElementSize, 2.0) ) / ( 6.0 * tProdCondDeltat );
            }
            // if there is no phase change
            else
            {
                tAlpha = ( tDensity * tHeatCapacity * std::pow(tElementSize, 2.0) ) / ( 6.0 * tProdCondDeltat );
            }

            // compute xi-bar value (for tAlpha -> 0: tXibar = 0.5)
            real tXibar = 0.5;

            // check if tAlpha is close to zero
            if ( tAlpha > mEpsilon )
            {
                tXibar = ( std::cosh( std::sqrt(6.0*tAlpha) ) + 2.0 ) / ( std::cosh( std::sqrt(6.0*tAlpha) ) - 1.0 )  -  (1.0/tAlpha);
            }

            // compute stabilization parameter value
            mPPVal = {{  std::pow(tElementSize, 2.0) / 6.0  * tXibar }};

            /* Note:
             * the artificial GGLS conductivity is returned as the stabilization parameter,
             * which is equal to k*tau_GGLS as defined in Alberto Pizzolato's Thesis,
             * or equal to just tau as defined in Franca's 1988 paper on GGLS
             */
        }

        //------------------------------------------------------------------------------

        void SP_GGLS_Diffusion::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // get the properties
            const std::shared_ptr< Property > & tPropConductivity = mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );
            const std::shared_ptr< Property > & tPropDensity      = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );
            const std::shared_ptr< Property > & tPropHeatCapacity = mMasterProp( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );
            const std::shared_ptr< Property > & tPropLatentHeat   = mMasterProp( static_cast< uint >( Property_Type::LATENT_HEAT ) );
            const std::shared_ptr< Property > & tPropMeltTemp     = mMasterProp( static_cast< uint >( Property_Type::PC_TEMP ) );
            const std::shared_ptr< Property > & tPropPCconst      = mMasterProp( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) );
            const std::shared_ptr< Property > & tPropPSfunc       = mMasterProp( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) );

            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFIDer = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // time step size
            real tDeltat = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

            // initialize values
            real tConductivity = tPropConductivity->val()( 0 );
            real tDensity      = tPropDensity->val()( 0 );
            real tHeatCapacity = tPropHeatCapacity->val()( 0 );
            real tLatentHeat   = 0.0;
            real tMeltTemp     = 0.0;
            real tPCconst      = 0.0;
            real tPSfunc       = 0.0;
            real tAlpha        = 0.0;
            real tdfdT         = 0.0;

            // threshold product of conductivity and time step
            real tProdCondDeltat = std::max(tConductivity * tDeltat, mEpsilon);

            // if there is a phase change
            if (tPropLatentHeat != nullptr)
            {
                // get values of properties
                tLatentHeat = tPropLatentHeat->val()( 0 );
                tMeltTemp   = tPropMeltTemp->val()( 0 );
                tPCconst    = tPropPCconst->val()( 0 );
                tPSfunc     = tPropPSfunc->val()( 0 );

                // get phase state function
                tdfdT = eval_dFdTemp( tMeltTemp, tPCconst, tPSfunc, tFIDer );

                tAlpha = ( tDensity * ( tHeatCapacity + tLatentHeat * tdfdT ) * std::pow(tElementSize, 2.0) ) / ( 6.0 * tProdCondDeltat );
            }

            // if there is no phase change
            else
            {
                tAlpha = ( tDensity * tHeatCapacity * std::pow(tElementSize, 2.0) ) / ( 6.0 * tProdCondDeltat );
            }

            // compute derivatives if tAlpha is not thresholded
            if ( tAlpha > mEpsilon )
            {
                // get derivative of SP wrt to alpha d(k*tau)/d(alpha)
                real dXibardAlpha =
                        (4.0 -
                                8.0 * std::cosh( std::sqrt(6.0*tAlpha) ) +
                                4.0 * std::pow( std::cosh( std::sqrt(6.0*tAlpha) ) , 2.0 ) -
                                std::pow( 6.0*tAlpha , 1.5 ) * std::sinh( std::sqrt(6.0*tAlpha) ) )
                                / ( 4.0 * std::pow( tAlpha , 2.0 ) * std::pow( ( std::cosh( std::sqrt(6.0*tAlpha) ) - 1.0 ) , 2.0 ) );

                real tdSPdAlpha = std::pow(tElementSize, 2.0) / 6.0 * dXibardAlpha;

                // indirect contributions for both with and without phase change ---------------------------

                // if indirect dependency on conductivity
                if ( tPropConductivity->check_dof_dependency( aDofTypes ) )
                {
                    if ( tProdCondDeltat > mEpsilon )
                    {
                        mdPPdMasterDof( tDofIndex ) -=
                                tdSPdAlpha *
                                ( tDensity * tHeatCapacity * std::pow(tElementSize, 2.0) /
                                        ( 6.0 * std::pow(tConductivity, 2.0) * tDeltat ) ) *
                                        tPropConductivity->dPropdDOF( aDofTypes );
                    }
                }

                // if indirect dependency on density
                if ( tPropDensity->check_dof_dependency( aDofTypes ) )
                {
                    mdPPdMasterDof( tDofIndex ) +=
                            tdSPdAlpha *
                            ( tHeatCapacity * std::pow(tElementSize, 2.0) /
                                    ( 6.0 * tProdCondDeltat ) ) *
                                    tPropDensity->dPropdDOF( aDofTypes );
                }

                // if indirect dependency on heat capacity
                if ( tPropHeatCapacity->check_dof_dependency( aDofTypes ) )
                {
                    mdPPdMasterDof( tDofIndex ) +=
                            tdSPdAlpha *
                            ( tDensity * std::pow(tElementSize, 2.0) /
                                    ( 6.0 * tProdCondDeltat ) ) *
                                    tPropHeatCapacity->dPropdDOF( aDofTypes );
                }

                // if there is a phase change --------------------------------------------------------------
                if (tPropLatentHeat != nullptr)
                {
                    // if dof type is temperature
                    if( aDofTypes( 0 ) == mMasterDofTemp )
                    {
                        // get Dof-deriv of phase state function
                        moris::Matrix<DDRMat> td2fdTdDOF =
                                eval_dFdTempdDOF( tMeltTemp, tPCconst, tPSfunc, tFIDer);

                        // derivative of tau wrt temperature DOFs
                        mdPPdMasterDof( tDofIndex ) +=
                                tdSPdAlpha *
                                ( tDensity * tLatentHeat * std::pow(tElementSize, 2.0) /
                                        ( 6.0 * tProdCondDeltat ) ) *
                                        td2fdTdDOF;
                    }

                    // if indirect dependency on conductivity
                    if ( tPropConductivity->check_dof_dependency( aDofTypes ) )
                    {
                        if ( tProdCondDeltat > mEpsilon )
                        {
                            mdPPdMasterDof( tDofIndex ) -=
                                    tdSPdAlpha *
                                    ( tDensity * tLatentHeat * tdfdT * std::pow(tElementSize, 2.0) /
                                            ( 6.0 * std::pow(tConductivity, 2.0) * tDeltat ) ) *
                                            tPropConductivity->dPropdDOF( aDofTypes );
                        }
                    }

                    // if indirect dependency on density
                    if ( tPropDensity->check_dof_dependency( aDofTypes ) )
                    {
                        mdPPdMasterDof( tDofIndex ) +=
                                tdSPdAlpha *
                                ( tLatentHeat * tdfdT * std::pow(tElementSize, 2.0) /
                                        ( 6.0 * tProdCondDeltat ) ) *
                                        tPropDensity->dPropdDOF( aDofTypes );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

// BACKUP -----------------------------------------------
// -1.0 * (4.0 * std::cosh(std::sqrt(6.0*tAlpha)) -
//         2.0 * std::pow(std::cosh(std::sqrt(6.0*tAlpha)), 2.0) +
//         3.0 * std::sqrt(6.0 * std::pow(tAlpha,3.0) ) * std::sinh(std::sqrt(6.0*tAlpha)) - 2.0 ) /
//         ( 2.0 * std::pow(tAlpha,2.0) * std::pow((std::cosh(std::sqrt(6.0*tAlpha)) - 1.0), 2.0) );


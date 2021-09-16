
//FEM/INT/src
#include "cl_FEM_SP_SUPG_Advection.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
#include "fn_norm.hpp"
#include "fn_dot.hpp"

#include "fn_FEM_CM_Phase_State_Functions.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        SP_SUPG_Advection::SP_SUPG_Advection()
        {
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Conductivity" ]        = static_cast< uint >( Property_Type::CONDUCTIVITY );
            mPropertyMap[ "Density" ]             = static_cast< uint >( Property_Type::DENSITY );
            mPropertyMap[ "HeatCapacity" ]        = static_cast< uint >( Property_Type::HEAT_CAPACITY );
            mPropertyMap[ "LatentHeat" ]          = static_cast< uint >( Property_Type::LATENT_HEAT );
            mPropertyMap[ "PCTemp" ]              = static_cast< uint >( Property_Type::PC_TEMP );
            mPropertyMap[ "PhaseStateFunction" ]  = static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION );
            mPropertyMap[ "PhaseChangeConst" ]    = static_cast< uint >( Property_Type::PHASE_CHANGE_CONST );
            mPropertyMap[ "Source" ]              = static_cast< uint >( Property_Type::SOURCE );
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Advection::set_dof_type_list(
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
                        if( tDofString == "Velocity" )
                        {
                            mMasterDofVelocity = tDofType;
                        }
                        else if( tDofString == "ScalarField" )
                        {
                            mMasterDofScalarField = tDofType;
                        }
                        else
                        {
                            // error unknown dof string
                            MORIS_ERROR( false ,
                                    "SP_SUPG_Advection::set_dof_type_list - Unknown aDofString : %s \n",
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
                    MORIS_ERROR( false, "SP_SUPG_Advection::set_dof_type_list - unknown master slave type." );
            }
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Advection::compute_effective_conductivity()
        {
            // get the conductivity property
            const std::shared_ptr< Property > & tPropConductivity =
                    mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // get the density property
            const std::shared_ptr< Property > & tPropDensity =
                    mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get the capacity property
            const std::shared_ptr< Property > & tPropCapacity =
                    mMasterProp( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );

            // get the latent heat property
            const std::shared_ptr< Property > & tPropLatentHeat =
                    mMasterProp( static_cast< uint >( Property_Type::LATENT_HEAT ) );

            // check that conductivity is set
            MORIS_ASSERT( tPropConductivity != nullptr,
                    "SP_SUPG_Advection::compute_effective_conductivity - conductivity not defined\n");

            // get contribution of density to effective conductivity
            real tDensity = 1.0;
            if ( tPropDensity != nullptr )
            {
                tDensity = tPropDensity->val()( 0 );
            }

            // get contribution of capacity to effective conductivity
            real tCapacity = 1.0;
            if ( tPropCapacity != nullptr )
            {
                tCapacity = tPropCapacity->val()( 0 );
            }

            // get contribution of latent to effective conductivity
            real tLatentHeatContrib = 0.0;
            if ( tPropLatentHeat != nullptr )
            {
                // get the phase change properties
                const std::shared_ptr< Property > & tPropPCTemp =
                        mMasterProp( static_cast< uint >( Property_Type::PC_TEMP ) );

                const std::shared_ptr< Property > & tPropPhaseChangeFunction =
                        mMasterProp( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) );

                const std::shared_ptr< Property > & tPropPhaseChangeConstant =
                        mMasterProp( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) );

                // check that all phase change properties are set
                MORIS_ASSERT( tPropPCTemp != nullptr and tPropPhaseChangeFunction != nullptr and tPropPhaseChangeConstant!= nullptr ,
                        "SP_SUPG_Advection::compute_effective_conductivity - some or all change properties are not defined\n");

                // get the scalar field FI
                Field_Interpolator * tFIScalarField =
                        mMasterFIManager->get_field_interpolators_for_type( mMasterDofScalarField );

                 // compute derivative of Phase State Function
                 real tdfdT = eval_dFdTemp(
                         tPropPCTemp->val()( 0 ),
                         tPropPhaseChangeConstant->val()( 0 ),
                         tPropPhaseChangeFunction->val()( 0 ),
                         tFIScalarField );

                // compute contribution of latent heat
                 tLatentHeatContrib = tPropLatentHeat->val()( 0 ) * tdfdT;
            }

            // compute effective conductivity
            return tPropConductivity->val()( 0 ) / (tDensity * tCapacity + tLatentHeatContrib);
        }

        //------------------------------------------------------------------------------

        bool SP_SUPG_Advection::compute_derivative_of_effective_conductivity(
                Matrix< DDRMat >                   & aEffectiveConductivitydu,
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the conductivity property
            const std::shared_ptr< Property > & tPropConductivity =
                    mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // get the density property
            const std::shared_ptr< Property > & tPropDensity =
                    mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );

            // get the capacity property
            const std::shared_ptr< Property > & tPropCapacity =
                    mMasterProp( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );

            // get the latent heat property
            const std::shared_ptr< Property > & tPropLatentHeat =
                    mMasterProp( static_cast< uint >( Property_Type::LATENT_HEAT ) );

            // get conductivity
            real tConductivity = tPropConductivity->val()( 0 );

            // get contribution of density to effective conductivity
            real tDensity = 1.0;
            if ( tPropDensity != nullptr )
            {
                tDensity = tPropDensity->val()( 0 );
            }

            // get contribution of capacity to effective conductivity
            real tCapacity = 1.0;
            if ( tPropCapacity != nullptr )
            {
                tCapacity = tPropCapacity->val()( 0 );
            }

            // get contribution of latent to effective conductivity
            real tLatentHeatContrib = 0.0;
            if ( tPropLatentHeat != nullptr )
            {
                // get the phase change properties
                const std::shared_ptr< Property > & tPropPCTemp =
                        mMasterProp( static_cast< uint >( Property_Type::PC_TEMP ) );

                const std::shared_ptr< Property > & tPropPhaseChangeFunction =
                        mMasterProp( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) );

                const std::shared_ptr< Property > & tPropPhaseChangeConstant =
                        mMasterProp( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) );

                // get the temperature FI
                Field_Interpolator * tFIScalarField =
                        mMasterFIManager->get_field_interpolators_for_type( mMasterDofScalarField );

                // compute derivative of Phase State Function
                real tdfdT = eval_dFdTemp(
                        tPropPCTemp->val()( 0 ),
                        tPropPhaseChangeConstant->val()( 0 ),
                        tPropPhaseChangeFunction->val()( 0 ),
                        tFIScalarField );

                // add contribution of latent heat
                tLatentHeatContrib = tPropLatentHeat->val()( 0 ) * tdfdT;
            }

            // set flag for dependency to false
            bool tIsDependent = false;

            // consider dependency of conductivity on dof types
            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                aEffectiveConductivitydu = 1.0/(tDensity * tCapacity + tLatentHeatContrib) * tPropConductivity->dPropdDOF(aDofTypes);
                tIsDependent = true;
            }
            else
            {
                aEffectiveConductivitydu.fill(0.0);
            }

            // consider dependency of density on dof types
            if ( tPropDensity != nullptr )
            {
                if( tPropDensity->check_dof_dependency( aDofTypes ) )
                {
                    const real tFactor =  tConductivity * tCapacity / std::pow(tDensity * tDensity * tCapacity + tLatentHeatContrib, 2.0);

                    aEffectiveConductivitydu -= tFactor * tPropDensity->dPropdDOF(aDofTypes);

                    tIsDependent = true;
                }
            }

            // consider dependency of density on dof types
            if ( tPropCapacity != nullptr )
            {
                if( tPropCapacity->check_dof_dependency( aDofTypes ) )
                {
                    const real tFactor =  tConductivity * tDensity / std::pow(tDensity * tDensity * tCapacity + tLatentHeatContrib, 2.0);

                    aEffectiveConductivitydu -= tFactor * tPropCapacity->dPropdDOF(aDofTypes);

                    tIsDependent = true;
                }
            }

            // consider dependency of latent heat contribution on dof types
            if ( tPropLatentHeat != nullptr )
            {
                // get the phase change properties
                const std::shared_ptr< Property > & tPropPCTemp =
                        mMasterProp( static_cast< uint >( Property_Type::PC_TEMP ) );

                const std::shared_ptr< Property > & tPropPhaseChangeFunction =
                        mMasterProp( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) );

                const std::shared_ptr< Property > & tPropPhaseChangeConstant =
                        mMasterProp( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) );

                // consider dependency of phase state function on scalar field
                if( aDofTypes( 0 ) == mMasterDofScalarField )
                {
                    // get the scalar field FI
                    Field_Interpolator * tFIScalarField =
                            mMasterFIManager->get_field_interpolators_for_type( mMasterDofScalarField );

                    const moris::Matrix< DDRMat > dfdDof = eval_dFdTempdDOF(
                            tPropPCTemp->val()( 0 ),
                            tPropPhaseChangeConstant->val()( 0 ),
                            tPropPhaseChangeFunction->val()( 0 ),
                            tFIScalarField );

                    const real tFactor =  tConductivity * tPropLatentHeat->val()( 0 ) / std::pow(tDensity * tDensity * tCapacity + tLatentHeatContrib, 2.0);

                    aEffectiveConductivitydu -= tFactor * dfdDof;

                    tIsDependent = true;
                }

                // if latent heat depends on the dof type
                if ( tPropLatentHeat->check_dof_dependency( aDofTypes ) )
                {
                    MORIS_ERROR(false, "SP_SUPG_Advection::compute_derivative_of_effective_conductivity - %s\n",
                            "Dof dependence of Latent heat not implemented.\n");
                }
            }

            // return flag whether dof dependence exists
            return tIsDependent;
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Advection::eval_SP()
        {
            // get the velocity FI
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get the mass source
            const std::shared_ptr< Property > & tSourceProp =
                    mMasterProp( static_cast< uint >( Property_Type::SOURCE ) );

            // compute effective conductivity
            const real tEffectiveConductivity = this->compute_effective_conductivity();

            // compute and threshold the velocity norm (thresholding for consistency with derivatives)
            const real tNorm = std::max( norm( tVelocityFI->val() ), mEpsilon);

            // get the abs term
            const uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
            real tAbs = 0.0;

            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                Matrix< DDRMat > tAdd =
                        trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 ).get_column( iNode );
                tAbs += std::abs( tAdd( 0, 0 ) );
            }

            // threshold tAbs
            tAbs = std::max(tAbs, mEpsilon);

            // compute and threshold hugn
            const real tHugn = std::max( 2.0 * tNorm / tAbs, mEpsilon);

            // compute tau1
            const real tTau1 = 2.0 * tNorm / tHugn;

            // compute tau2
            const real tTau2 = 4.0 * tEffectiveConductivity / std::pow( tHugn, 2.0 );

            // compute time increment tDeltaT
            const real tDeltaT = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

            // compute tau3
            const real tTau3 = 2.0 / tDeltaT;

            // compute sum of square terms
            real tSum = std::pow( tTau1, 2.0 ) + std::pow( tTau2, 2.0 ) + std::pow( tTau3, 2.0 );

            // add contribution from source term
            if ( tSourceProp != nullptr )
            {
                tSum += std::pow( tSourceProp->val()(0), 2.0 );
            }

            // threshold sum of square terms
            tSum = std::max(tSum, mEpsilon);

            // compute stabilization parameter value
            mPPVal = {{ std::pow( tSum, -0.5 ) }};
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Advection::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            const uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof, dTau1dDof, dTau3dDof
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            Matrix< DDRMat > tdTau1dDof( 1, tFIDer->get_number_of_space_time_coefficients() );
            Matrix< DDRMat > tdTau2dDof( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the velocity FI
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get the source property
            const std::shared_ptr< Property > & tSourceProp =
                    mMasterProp( static_cast< uint >( Property_Type::SOURCE ) );

            // compute effective conductivity
            const real tEffectiveConductivity = this->compute_effective_conductivity();

            // compute and threshold the velocity norm (thresholding for consistency with derivatives)
            const real tNorm = std::max( norm( tVelocityFI->val() ), mEpsilon);

            // compute the abs term
            const uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
            real tAbs = 0.0;

            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                Matrix< DDRMat > tAdd = trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 ).get_column( iNode );
                tAbs += std::abs( tAdd( 0, 0 ) );
            }

            // threshold tAbs
            tAbs = std::max(tAbs, mEpsilon);

            // compute and threshold hugn
            const real tHugn = std::max( 2.0 * tNorm / tAbs, mEpsilon);

            // compute tau1
            const real tTau1 = 2.0 * tNorm / tHugn;

            // compute tau2
            const real tTau2 = 4.0 * tEffectiveConductivity / std::pow( tHugn, 2.0 );

            // compute time increment deltat
            const real tDeltaT = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

            // compute tau3
            const real tTau3 = 2.0 / tDeltaT;

            // if dof type is velocity
            if( aDofTypes( 0 ) == mMasterDofVelocity )
            {
                // compute derivative of hugn wrt velocity dof
                Matrix< DDRMat > tdHugndu( 1, tVelocityFI->get_number_of_space_time_coefficients() );
                Matrix< DDRMat > tdNormdu( 1, tVelocityFI->get_number_of_space_time_coefficients() );
                Matrix< DDRMat > tdAbsdu(  1, tVelocityFI->get_number_of_space_time_coefficients(), 0.0 );

                // compute derivative of the velocity norm (compute only derivative if not thresholded)
                if ( tNorm > mEpsilon )
                {
                    tdNormdu = trans( tVelocityFI->val() ) * tVelocityFI->N() / tNorm;
                }
                else
                {
                    tdNormdu.fill( 0.0 );
                }

                // compute derivative of the abs term (compute only derivative if not thresholded)
                if ( tAbs > mEpsilon )
                {
                    uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();

                    for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                    {
                        const Matrix< DDRMat > tAdd =
                                trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 ).get_column( iNode );

                        // handle case that tAdd( 0, 0 ) is smaller than threshold
                        if ( std::abs( tAdd( 0, 0 ) ) > mEpsilon )
                        {
                            tdAbsdu +=
                                    tAdd * trans( tVelocityFI->dnNdxn( 1 ).get_column( iNode ) ) *
                                    tVelocityFI->N() / std::abs( tAdd( 0, 0 ) );
                        }
                    }
                }

                // compute derivative of hugn (compute only derivative if not thresholded)
                if ( tHugn > mEpsilon )
                {
                    tdHugndu = 2.0 * ( tdNormdu * tAbs - tdAbsdu * tNorm ) / std::pow( tAbs, 2.0 );
                }
                else
                {
                    tdHugndu.fill( 0.0 );
                }

                // compute dtau1du
                tdTau1dDof = 2.0 * ( tHugn * tdNormdu - tdHugndu * tNorm ) / std::pow( tHugn, 2.0 );

                // compute dtau2du
                tdTau2dDof = - 8.0 * tEffectiveConductivity * tdHugndu / std::pow( tHugn, 3.0 );
            }
            else
            {
                tdTau1dDof.fill( 0.0 );
                tdTau2dDof.fill( 0.0 );
            }

            // consider derivative of effective conductivity on dof types
            Matrix< DDRMat > tEffectiveConductivitydu( 1, tFIDer->get_number_of_space_time_coefficients() );

            bool tIsDependent = this->compute_derivative_of_effective_conductivity(
                    tEffectiveConductivitydu,
                    aDofTypes );

            // compute dtau2du
            if ( tIsDependent )
            {
                tdTau2dDof += 4.0 * tEffectiveConductivitydu / std::pow( tHugn, 2.0 );
            }

            // compute sum of square terms
            real tSum = std::pow( tTau1, 2.0 ) + std::pow( tTau2, 2.0 ) + std::pow( tTau3, 2.0 );

            // add contribution from source term
            if ( tSourceProp != nullptr )
            {
                tSum += std::pow( tSourceProp->val()(0), 2.0 );
            }

            // compute dSPdu
            if ( tSum > mEpsilon )
            {
                const real tPrefactor = - std::pow( tSum, -1.5 );

                mdPPdMasterDof( tDofIndex ) = tPrefactor * ( tTau1 * tdTau1dDof + tTau2 * tdTau2dDof);

                if ( tSourceProp != nullptr )
                {
                    if( tSourceProp->check_dof_dependency( aDofTypes ) )
                    {
                        // compute dtau3du
                        mdPPdMasterDof( tDofIndex ) +=  tPrefactor * tSourceProp->val()(0) * tSourceProp->dPropdDOF( aDofTypes );
                    }
                }
            }
            else
            {
                mdPPdMasterDof( tDofIndex ).fill( 0.0 );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mdPPdMasterDof( tDofIndex ) ),
                    "SP_SUPG_Advection::eval_dSPdMasterDOF - mdPPdMasterDof contains NAN or INF, exiting for tDofIndex = %d !\n",
                    tDofIndex);
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

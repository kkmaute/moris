//FEM/INT/src
#include "cl_FEM_SP_SUPG_Spalart_Allmaras_Turbulence.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "fn_inv.hpp"
#include "op_div.hpp"
#include "fn_diag_vec.hpp"
#include "fn_isfinite.hpp"

namespace moris
{
    namespace fem
    {
        SP_SUPG_Spalart_Allmaras_Turbulence::SP_SUPG_Spalart_Allmaras_Turbulence()
        {
            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "SpalartAllmarasTurbulence" ] = static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::set_dof_type_list(
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
                        else if( tDofString == "Viscosity" )
                        {
                            mMasterDofViscosity = tDofType;
                        }
                        else
                        {
                            // error unknown dof string
                            MORIS_ERROR( false ,
                                    "SP_SUPG_Spalart_Allmaras_Turbulence::set_dof_type_list - Unknown aDofString : %s \n",
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
                    MORIS_ERROR( false, "SP_SUPG_Spalart_Allmaras_Turbulence::set_dof_type_list - unknown master slave type." );
            }
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::tuple<
        fem::Measure_Type,
        mtk::Primary_Void,
        mtk::Master_Slave > > SP_SUPG_Spalart_Allmaras_Turbulence::get_cluster_measure_tuple_list()
        {
            return { mElementSizeTuple };
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::set_function_pointers(){}

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::eval_SP()
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // set size for SP values
            mPPVal.set_size( 1, 1 );

            // get the velocity FI
            Field_Interpolator * tViscosityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the velocity FI
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model > & tCMSATurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // compute uTilde = u - cb2/sigma gradv
            Matrix< DDRMat > tModVelocity =
                    tVelocityFI->val() - mCb2 * tViscosityFI->gradx( 1 ) / mSigma;

            // compute norm( uTilde )
            real tNormA = std::sqrt( dot( tModVelocity, tModVelocity ) );

            // threshold tNormA (for consistency with derivative computation)
            tNormA = std::max(tNormA,mEpsilon);

            // tau A
            real tTauA = 2.0 * tNormA / tElementSize;

            // tau K
            real tTauK = 4.0 * tCMSATurbulence->diffusion_coefficient()( 0 ) / std::pow( tElementSize, 2.0 );

            // tau S
            real tTauS =
                    tCMSATurbulence->production_coefficient()( 0 ) +
                    tCMSATurbulence->wall_destruction_coefficient()( 0 );

            // evaluate tau
            real tTau = std::pow( tTauA, 2 ) + std::pow( tTauK, 2 ) + std::pow( tTauS, 2 );

            // threshold tau
            tTau = std::max(tTau, mEpsilon);

            // set tau
            mPPVal = {{ std::pow( tTau, -0.5 ) }};
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients() );

            // get the velocity FI
            Field_Interpolator * tViscosityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the velocity FI
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model > & tCMSATurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // evaluate uTilde = u - 1/sigma cb2 gradv
            Matrix< DDRMat > tModVelocity =
                    tVelocityFI->val() - mCb2 * tViscosityFI->gradx( 1 ) / mSigma;

            // evaluate norm( uTilde )
            real tNormA = std::sqrt( dot( tModVelocity, tModVelocity ) );

            // threshold tNormA (for consistency with derivative computation)
            tNormA = std::max(tNormA,mEpsilon);

            // tau A
            real tTauA = 2.0 * tNormA / tElementSize;

            // tau K
            real tTauK = 4.0 * tCMSATurbulence->diffusion_coefficient()( 0 ) / std::pow( tElementSize, 2.0 );

            // tau S
            real tTauS =
                    tCMSATurbulence->production_coefficient()( 0 ) +
                    tCMSATurbulence->wall_destruction_coefficient()( 0 );

            // evaluate tau
            real tTau = std::pow( tTauA, 2 ) + std::pow( tTauK, 2 ) + std::pow( tTauS, 2 );

            // compute dSPdu
            if ( tTau > mEpsilon )
            {
                // if normA greater than threshold
                if( tNormA > mEpsilon )
                {
                    // if dof type is velocity
                    if( aDofTypes( 0 ) == mMasterDofVelocity )
                    {
                        // add contribution to mdPPdMasterDof
                        mdPPdMasterDof( tDofIndex ) = tTauA * 2.0 *
                                trans( tModVelocity ) * tVelocityFI->N() /
                                ( tElementSize * tNormA );
                    }
                    // if dof type is viscosity
                    else if( aDofTypes( 0 ) == mMasterDofViscosity )
                    {
                        // add contribution to mdPPdMasterDof
                        mdPPdMasterDof( tDofIndex ) = - tTauA * 2.0 *
                                mCb2 * trans( tModVelocity ) * tViscosityFI->dnNdxn( 1 ) / mSigma /
                                ( tElementSize * tNormA );
                    }
                    else
                    {
                        mdPPdMasterDof( tDofIndex ).fill( 0.0 );
                    }
                }
                else
                {
                    mdPPdMasterDof( tDofIndex ).fill( 0.0 );
                }

                // if turbulence CM depends on dof type
                if( tCMSATurbulence->check_dof_dependency( aDofTypes ) )
                {
                    // compute tdtauKdu
                    mdPPdMasterDof( tDofIndex ) +=
                            tTauK * 4.0 * tCMSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                            std::pow( tElementSize, 2.0 );

                    // compute tdtauSdu
                    mdPPdMasterDof( tDofIndex ) += tTauS * (
                            tCMSATurbulence->dproductioncoeffdu( aDofTypes ) +
                            tCMSATurbulence->dwalldestructioncoeffdu( aDofTypes ) );
                }
            }
            else
            {
                mdPPdMasterDof( tDofIndex ).fill( 0.0 );
            }

            // scale
            mdPPdMasterDof( tDofIndex ) = - std::pow( tTau, -1.5 ) * mdPPdMasterDof( tDofIndex );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


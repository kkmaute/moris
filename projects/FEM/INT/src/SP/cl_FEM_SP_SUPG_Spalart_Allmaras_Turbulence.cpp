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
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Viscosity" ]    = static_cast< uint >( SP_Property_Type::VISCOSITY );
            mPropertyMap[ "WallDistance" ] = static_cast< uint >( SP_Property_Type::WALL_DISTANCE );
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

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );

            // get the wall distance property
            const std::shared_ptr< Property > & tPropWallDistance =
                    mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );

            // compute uTilde = u - cb2/sigma gradv
            Matrix< DDRMat > tModVelocity =
                    tVelocityFI->val() - mCb2 * tViscosityFI->gradx( 1 ) / mSigma;

            // compute norm( uTilde )
            real tNormA = std::sqrt( dot( tModVelocity, tModVelocity ) );

            // threshold tNormA (for consistency with derivative computation)
            tNormA = std::max(tNormA,mEpsilon);

            // tau A
            real tTauA = 2.0 * tNormA / tElementSize;

            // compute diffusion coefficient
            real tDiffusionCoeff = compute_diffusion_coefficient(
                    { mMasterDofViscosity },
                    mMasterFIManager,
                    tPropKinViscosity );

            // tau K
            real tTauK = 4.0 * tDiffusionCoeff / std::pow( tElementSize, 2.0 );

            // compute production coefficient
            real tProductionTerm = compute_production_coefficient(
                    { mMasterDofViscosity },
                    { mMasterDofVelocity },
                    mMasterFIManager,
                    tPropKinViscosity,
                    tPropWallDistance );

            // compute wall destruction coefficient
            real tWallDestructionTerm = compute_wall_destruction_coefficient(
                    { mMasterDofViscosity },
                    { mMasterDofVelocity },
                    mMasterFIManager,
                    tPropKinViscosity,
                    tPropWallDistance );

            // tau S
            real tTauS = tProductionTerm + tWallDestructionTerm;

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

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );

            // get the wall distance property
            const std::shared_ptr< Property > & tPropWallDistance =
                    mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );

            // evaluate uTilde = u - 1/sigma cb2 gradv
            Matrix< DDRMat > tModVelocity =
                    tVelocityFI->val() - mCb2 * tViscosityFI->gradx( 1 ) / mSigma;

            // evaluate norm( uTilde )
            real tNormA = std::sqrt( dot( tModVelocity, tModVelocity ) );

            // threshold tNormA (for consistency with derivative computation)
            tNormA = std::max(tNormA,mEpsilon);

            // tau A
            real tTauA = 2.0 * tNormA / tElementSize;

            // compute diffusion coefficient
            real tDiffusionCoeff = compute_diffusion_coefficient(
                    { mMasterDofViscosity },
                    mMasterFIManager,
                    tPropKinViscosity );

            // tau K
            real tTauK = 4.0 * tDiffusionCoeff / std::pow( tElementSize, 2.0 );

            // compute production term
            real tProductionTerm = compute_production_coefficient(
                    { mMasterDofViscosity },
                    { mMasterDofVelocity },
                    mMasterFIManager,
                    tPropKinViscosity,
                    tPropWallDistance );

            // compute wall destruction term
            real tWallDestructionTerm = compute_wall_destruction_coefficient(
                    { mMasterDofViscosity },
                    { mMasterDofVelocity },
                    mMasterFIManager,
                    tPropKinViscosity,
                    tPropWallDistance );

            // tau S
            real tTauS = tProductionTerm + tWallDestructionTerm;

            // init derivative of tauA, tauK, tauS
            Matrix< DDRMat > tdtauAdu( 1, tFI->get_number_of_space_time_coefficients() );
            Matrix< DDRMat > tdtauKdu( 1, tFI->get_number_of_space_time_coefficients() );
            Matrix< DDRMat > tdtauSdu( 1, tFI->get_number_of_space_time_coefficients() );

            // if normA greater than threshold
            if( tNormA > mEpsilon )
            {
                // if dof type is velocity
                if( aDofTypes( 0 ) == mMasterDofVelocity )
                {
                    // add contribution to dtauAdu
                    tdtauAdu = trans( tModVelocity ) * tVelocityFI->N();
                }
                // if dof type is viscosity
                else if( aDofTypes( 0 ) == mMasterDofViscosity )
                {
                    // add contribution to dtauAdu
                    tdtauAdu = - mCb2 * trans( tModVelocity ) * tViscosityFI->dnNdxn( 1 ) / mSigma;
                }
                else
                {
                    tdtauAdu.fill( 0.0 );
                }

                // multiply by common scaling factor
                tdtauAdu = 2.0 * tdtauAdu / ( tElementSize * tNormA );
            }

            // compute derivative of diffusion coefficient
            Matrix< DDRMat > tddiffusiondu;
            compute_ddiffusiondu(
                    { mMasterDofViscosity },
                    mMasterFIManager,
                    tPropKinViscosity,
                    aDofTypes,
                    tddiffusiondu );

            // compute tdtauKdu
            tdtauKdu = 4.0 * tddiffusiondu / std::pow( tElementSize, 2.0 );

            // compute derivative of production term
            Matrix< DDRMat > tdproductiondu;
            compute_dproductioncoefficientdu(
                    { mMasterDofViscosity },
                    { mMasterDofVelocity },
                    mMasterFIManager,
                    tPropKinViscosity,
                    tPropWallDistance,
                    aDofTypes,
                    tdproductiondu );

            // compute derivative of wall destruction term
            Matrix< DDRMat > tdwalldestructiondu;
            compute_dwalldestructioncoefficientdu(
                    { mMasterDofViscosity },
                    { mMasterDofVelocity },
                    mMasterFIManager,
                    tPropKinViscosity,
                    tPropWallDistance,
                    aDofTypes,
                    tdwalldestructiondu );

            // compute tdtauSdu
            tdtauSdu = tdproductiondu + tdwalldestructiondu;

            // evaluate tau
            real tTau = std::pow( tTauA, 2 ) + std::pow( tTauK, 2 ) + std::pow( tTauS, 2 );

            // compute dSPdu
            if ( tTau > mEpsilon )
            {
                mdPPdMasterDof( tDofIndex ) = - std::pow( tTau, -1.5 ) *
                    ( tTauA * tdtauAdu + tTauK * tdtauKdu + tTauS * tdtauSdu );
            }
            else
            {
                mdPPdMasterDof( tDofIndex ).fill( 0.0 );
            }
        }

//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_production_coefficient()
//        {
//            // init production term
//            real tProductionTerm = 0.0;
//
//            // get the viscosity FI
//            Field_Interpolator * tFIModViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//            // get the viscosity value
//            real tModViscosity = tFIModViscosity->val()( 0 );
//
//            // if viscosity is positive or zero
//            if( tModViscosity >= 0.0 )
//            {
//                // compute ft2
//                real tFt2 = this->compute_ft2();
//
//                // compute Stilde
//                real tSTilde = this->compute_stilde();
//
//                // compute production term
//                tProductionTerm = mCb1 * ( 1.0 - tFt2 ) * tSTilde;
//            }
//            // if viscosity is negative
//            else
//            {
//                // compute S
//                real tS = this->compute_s();
//
//                // compute production term
//                tProductionTerm = mCb1 * ( 1.0 - mCt3 ) * tS;
//            }
//
//            return tProductionTerm;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dproductiondu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adproductiondu )
//        {
//            // get derivative dof type
//            MSI::Dof_Type tDerDofType = aDofTypes( 0 );
//
//            // get the derivative dof FI
//            Field_Interpolator * tFIDer =
//                    mMasterFIManager->get_field_interpolators_for_type( tDerDofType );
//
//            // init adproductiondu
//            adproductiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//
//            // get the viscosity FI
//            Field_Interpolator * tFIModViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//            // get the viscosity value
//            real tModViscosity = tFIModViscosity->val()( 0 );
//
//            // if viscosity is positive or zero
//            if( tModViscosity >= 0.0 )
//            {
//                // compute ft2
//                real tFt2 = this->compute_ft2();
//
//                // compute dft2du
//                Matrix< DDRMat > tdft2du;
//                this->compute_dft2du( aDofTypes, tdft2du );
//
//                // compute STilde
//                real tSTilde = this->compute_stilde();
//
//                // compute dSTildedu
//                Matrix< DDRMat > tdstildedu;
//                this->compute_dstildedu( aDofTypes, tdstildedu );
//
//                // add contribution to dproductiondu
//                adproductiondu += - mCb1 * tSTilde * tdft2du +
//                        mCb1 * ( 1 - tFt2 ) * tdstildedu;
//            }
//            // if viscosity is negative
//            else
//            {
//                // compute dsdu
//                Matrix< DDRMat > tdsdu;
//                this->compute_dsdu( aDofTypes, tdsdu );
//
//                // add contribution to dproductiondu
//                adproductiondu += mCb1 * ( 1.0 - mCt3 ) * tdsdu;
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_wall_destruction_coefficient()
//        {
//            // init wall destruction term
//            real tWallDestructionTerm = 0.0;
//
//            // get the viscosity FI
//            Field_Interpolator * tFIModViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//            // get the wall distance property
//            std::shared_ptr< Property > & tPropWallDistance =
//                    mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );
//
//            // get the viscosity value
//            real tModViscosity = tFIModViscosity->val()( 0 );
//
//            // get the wall distance value
//            real tWallDistance = tPropWallDistance->val()( 0 );
//
////            // check negative/zero wall distance
////            MORIS_ERROR( tWallDistance > 0.0,
////                    "SP_SUPG_Spalart_Allmaras_Turbulence::compute_wall_destruction_coefficient - %s",
////                    "Negative or zero wall distance.\n");
//
//            // threshold wall distance
//            tWallDistance = std::max(tWallDistance,mEpsilon);
//
//            // if viscosity is positive or zero
//            if( tModViscosity >= 0.0 )
//            {
//                // compute fw
//                real tFw = this->compute_fw();
//
//                // compute ft2
//                real tFt2 = this->compute_ft2();
//
//                // compute wall destruction term
//                tWallDestructionTerm =
//                        ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
//                        tModViscosity / std::pow( tWallDistance, 2.0 );
//            }
//            // if viscosity is negative
//            else
//            {
//                // compute wall destruction term
//                tWallDestructionTerm = - mCw1 * tModViscosity / std::pow( tWallDistance, 2.0 );
//            }
//
//            return tWallDestructionTerm;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dwalldestructiondu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adwalldestructiondu )
//        {
//            // get derivative dof type
//            MSI::Dof_Type tDerDofType = aDofTypes( 0 );
//
//            // get the derivative dof FI
//            Field_Interpolator * tFIDer =
//                    mMasterFIManager->get_field_interpolators_for_type( tDerDofType );
//
//            // init adwalldestructiondu
//            adwalldestructiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//
//            // get the viscosity FI
//            Field_Interpolator * tFIModViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//            // get the wall distance property
//            std::shared_ptr< Property > & tPropWallDistance =
//                    mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );
//
//            // get the viscosity value
//            real tModViscosity = tFIModViscosity->val()( 0 );
//
//            // get the wall distance value
//            real tWallDistance = tPropWallDistance->val()( 0 );
//
////            // check negative/zero wall distance
////            MORIS_ERROR( tWallDistance > 0.0,
////                    "SP_SUPG_Spalart_Allmaras_Turbulence::compute_dwalldestructiondu - %s",
////                    "Negative or zero wall distance.\n");
//
//            // threshold wall distance
//            tWallDistance = std::max(tWallDistance,mEpsilon);
//
//            // if viscosity is positive or zero
//            if( tModViscosity >= 0.0 )
//            {
//                // compute fw
//                real tFw = this->compute_fw();
//
//                // compute dfwfu
//                Matrix< DDRMat > tdfwdu;
//                this->compute_dfwdu( aDofTypes, tdfwdu );
//
//                // compute ft2
//                real tFt2 = this->compute_ft2();
//
//                // compute dft2fu
//                Matrix< DDRMat > tdft2du;
//                this->compute_dft2du( aDofTypes, tdft2du );
//
//                // if derivative dof type is viscosity
//                if( tDerDofType == mMasterDofViscosity )
//                {
//                    // add contribution to dwalldestructiondu
//                    adwalldestructiondu +=
//                            ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
//                            tFIModViscosity->N() /
//                            std::pow( tWallDistance, 2.0 );
//                }
//
//                // add contribution to dwalldestructiondu
//                adwalldestructiondu +=
//                        ( mCw1 * tdfwdu - mCb1 * tdft2du / std::pow( mKappa, 2.0 ) ) *
//                        tModViscosity / std::pow( tWallDistance, 2.0 );
//
//                // if wall distance depends on derivative dof type
//                if( ( tPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
//                        ( tWallDistance > mEpsilon ) )
//                {
//                    //MORIS_LOG( "compute_dwalldestructiondu - Wall distance depends on dof type - not tested yet" );
//
//                    // add contribution to dwalldestructiondu
//                    adwalldestructiondu -=
//                            2.0 * ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
//                            std::pow( tModViscosity, 2.0 ) *
//                            tPropWallDistance->dPropdDOF( aDofTypes ) /
//                            std::pow( tWallDistance, 3.0 );
//                }
//            }
//            // if viscosity is negative
//            else
//            {
//                // if derivative dof type is viscosity
//                if( tDerDofType == mMasterDofViscosity )
//                {
//                    // add contribution to dwalldestructiondu
//                    adwalldestructiondu -=
//                            mCw1 * tFIModViscosity->N() /
//                            std::pow( tWallDistance, 2.0 );
//                }
//
//                // if wall distance depends on derivative dof type
//                if( ( tPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
//                        ( tWallDistance > mEpsilon ) )
//                {
//                    //MORIS_LOG( "compute_dwalldestructiondu - Wall distance depends on dof type - not tested yet" );
//
//                    // add contribution to dwalldestructiondu
//                    adwalldestructiondu +=
//                            2.0 * mCw1 * std::pow( tModViscosity, 2.0 ) *
//                            tPropWallDistance->dPropdDOF( aDofTypes ) /
//                            std::pow( tWallDistance, 3.0 );
//                }
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_diffusion_coefficient()
//        {
//            // init diffusion coeff
//            real tDiffusionTerm = 0.0;
//
//            // get the viscosity FI
//            Field_Interpolator * tFIModViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//            // get the viscosity property
//            std::shared_ptr< Property > & tPropKinViscosity =
//                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );
//
//            // get the viscosity value
//            real tModViscosity = tFIModViscosity->val()( 0 );
//
//            // get the fluid kinematic viscosity value
//            real tKinViscosity = tPropKinViscosity->val()( 0 );
//
//            // if viscosity is positive or zero
//            if( tModViscosity >= 0.0 )
//            {
//                // compute diffusion term
//                tDiffusionTerm = ( tKinViscosity + tModViscosity ) / mSigma;
//            }
//            // if viscosity is negative
//            else
//            {
//                // compute fn
//                real tFn = this->compute_fn();
//
//                // compute diffusion term
//                tDiffusionTerm = ( tKinViscosity + tModViscosity * tFn ) / mSigma;
//            }
//
//            return tDiffusionTerm;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_ddiffusiondu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & addiffusiondu )
//        {
//            // get derivative dof type
//            MSI::Dof_Type tDerDofType = aDofTypes( 0 );
//
//            // get the derivative dof FI
//            Field_Interpolator * tFIDer =
//                    mMasterFIManager->get_field_interpolators_for_type( tDerDofType );
//
//            // init ddiffusiondu
//            addiffusiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//
//            // get the viscosity FI
//            Field_Interpolator * tFIModViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//            // get the fluid  kinematic viscosity property
//            std::shared_ptr< Property > & tPropKinViscosity =
//                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );
//
//            // get the viscosity value
//            real tModViscosity = tFIModViscosity->val()( 0 );
//
//            // if viscosity is positive or zero
//            if( tModViscosity >= 0.0 )
//            {
//                // if derivative dof type is viscosity
//                if( tDerDofType == mMasterDofViscosity )
//                {
//                    // add contribution to ddiffusiondu
//                    addiffusiondu += tFIModViscosity->N() / mSigma;
//                }
//
//                // if kinematic viscosity depends on derivative dof type
//                if( tPropKinViscosity->check_dof_dependency( aDofTypes ) )
//                {
//                    // add contribution to ddiffusiondu
//                    addiffusiondu += tPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
//                }
//            }
//            // if viscosity is negative
//            else
//            {
//                // compute fn
//                real tFn = this->compute_fn();
//
//                // compute dfndu
//                Matrix< DDRMat > tdfndu;
//                this->compute_dfndu( aDofTypes, tdfndu );
//
//                // if derivative dof type is viscosity
//                if( tDerDofType == mMasterDofViscosity )
//                {
//                    // add contribution to ddiffusiondu
//                    addiffusiondu += tFn * tFIModViscosity->N() / mSigma;
//                }
//
//                // if kinematic viscosity depends on derivative dof type
//                if( tPropKinViscosity->check_dof_dependency( aDofTypes ) )
//                {
//                    // add contribution to ddiffusiondu
//                    addiffusiondu += tPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
//                }
//
//                // add contribution from fn to ddiffusiondu
//                addiffusiondu += tModViscosity * tdfndu / mSigma;
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_wij( Matrix< DDRMat > & aWij )
//        {
//            // get the velocity FI
//            Field_Interpolator * tFIVelocity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );
//
//            // get gradient of velocity
//            Matrix< DDRMat > tGradVelocity = tFIVelocity->gradx( 1 );
//
//            // switch on space dim
//            switch ( mSpaceDim )
//            {
//                case 2 :
//                {
//                    // init aWij = [ w11 w12 w21 w22]
//                    aWij.set_size( 4, 1, 0.0 );
//
//                    // compute Wij
//                    aWij( 1 ) = 0.5 * ( tGradVelocity( 1, 0 ) - tGradVelocity( 0, 1 ) );
//                    aWij( 2 ) = 0.5 * ( tGradVelocity( 0, 1 ) - tGradVelocity( 1, 0 ) );
//                    break;
//                }
//                case 3 :
//                {
//                    // init aWij = [ w11 w12 w13 w21 w22 w23 w31 w32 w33 ]
//                    aWij.set_size( 9, 1, 0.0 );
//
//                    // compute Wij
//                    aWij( 1 ) = 0.5 * ( tGradVelocity( 1, 0 ) - tGradVelocity( 0, 1 ) );
//                    aWij( 2 ) = 0.5 * ( tGradVelocity( 2, 0 ) - tGradVelocity( 0, 2 ) );
//                    aWij( 3 ) = 0.5 * ( tGradVelocity( 0, 1 ) - tGradVelocity( 1, 0 ) );
//                    aWij( 5 ) = 0.5 * ( tGradVelocity( 2, 1 ) - tGradVelocity( 1, 2 ) );
//                    aWij( 6 ) = 0.5 * ( tGradVelocity( 0, 2 ) - tGradVelocity( 2, 0 ) );
//                    aWij( 7 ) = 0.5 * ( tGradVelocity( 1, 2 ) - tGradVelocity( 2, 1 ) );
//                    break;
//                }
//                default:
//                    MORIS_ERROR( false, "SP_SUPG_Spalart_Allmaras_Turbulence::compute_wij - space dim can only be 2 or 3" );
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dwijdu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adwijdu )
//        {
//            // get the der FI
//            Field_Interpolator * tFIDer =
//                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
//
//            // get the velocity FI
//            Field_Interpolator * tFIVelocity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );
//
//            // switch on space dim
//            switch ( mSpaceDim )
//            {
//                case 2 :
//                {
//                    // init aWij = [ w11 w12 w21 w22]
//                    adwijdu.set_size( 4, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//
//                    if( aDofTypes( 0 ) == MSI::Dof_Type::VX )
//                    {
//                        // get gradient of velocity
//                        Matrix< DDRMat > tdNdxVelocity = tFIVelocity->dnNdxn( 1 );
//
//                        // get number of bases for displacement
//                        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();
//
//                        // compute adwijdu
//                        adwijdu( { 1, 1 }, { 0, tNumBases - 1 } )             =   0.5 * tdNdxVelocity.get_row( 1 );
//                        adwijdu( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 0 );
//                        adwijdu( { 2, 2 }, { 0, tNumBases - 1 } )             = - 0.5 * tdNdxVelocity.get_row( 1 );
//                        adwijdu( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) =   0.5 * tdNdxVelocity.get_row( 0 );
//                    }
//                    break;
//                }
//                case 3 :
//                {
//                    // init aWij = [ w11 w12 w13 w21 w22 w23 w31 w32 w33 ]
//                    adwijdu.set_size( 9, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//
//                    if( aDofTypes( 0 ) == MSI::Dof_Type::VX )
//                    {
//                        // get gradient of velocity
//                        Matrix< DDRMat > tdNdxVelocity = tFIVelocity->dnNdxn( 1 );
//
//                        // get number of bases for displacement
//                        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();
//
//                        // compute adwijdu
//                        adwijdu( { 1, 1 }, { 0, tNumBases - 1 } )                 =   0.5 * tdNdxVelocity.get_row( 1 );
//                        adwijdu( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } )     = - 0.5 * tdNdxVelocity.get_row( 0 );
//
//                        adwijdu( { 2, 2 }, { 0, tNumBases - 1 } )                 =   0.5 * tdNdxVelocity.get_row( 2 );
//                        adwijdu( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 0 );
//
//                        adwijdu( { 3, 3 }, { 0, tNumBases - 1 } )                 = - 0.5 * tdNdxVelocity.get_row( 1 );
//                        adwijdu( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } )     =   0.5 * tdNdxVelocity.get_row( 0 );
//
//                        adwijdu( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } )     =   0.5 * tdNdxVelocity.get_row( 2 );
//                        adwijdu( { 5, 5 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 1 );
//
//                        adwijdu( { 6, 6 }, { 0, tNumBases - 1 } )                 = - 0.5 * tdNdxVelocity.get_row( 2 );
//                        adwijdu( { 6, 6 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =   0.5 * tdNdxVelocity.get_row( 0 );
//
//                        adwijdu( { 7, 7 }, { tNumBases, 2 * tNumBases - 1 } )     = - 0.5 * tdNdxVelocity.get_row( 2 );
//                        adwijdu( { 7, 7 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =   0.5 * tdNdxVelocity.get_row( 1 );
//                    }
//                    break;
//                }
//                default:
//                    MORIS_ERROR( false, "SP_SUPG_Spalart_Allmaras_Turbulence::compute_dwijdu - space dim can only be 2 or 3" );
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_s()
//        {
//            // compute wij
//            Matrix< DDRMat > tWij;
//            this->compute_wij( tWij );
//
//            // compute WijWij
//            Matrix< DDRMat > tWijWij = trans( tWij ) * tWij;
//
//            // compute s
//            real tS = std::sqrt( 2.0 * tWijWij( 0 ) );
//
//            // threshold s for consistency with derivative
//            return std::max(tS,mEpsilon);
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dsdu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adsdu )
//        {
//            // get the derivative dof FIs
//            Field_Interpolator * tDerFI =
//                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
//
//            // initialize dsdu
//            adsdu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );
//
//            // compute wij
//            Matrix< DDRMat > tWij;
//            this->compute_wij( tWij );
//
//            // compute WijWij
//            Matrix< DDRMat > tWijWij = trans( tWij ) * tWij;
//
//            // compute s
//            real tS = std::sqrt( 2.0 * tWijWij( 0 ) );
//
//            // compute derivative
//            if (tS > mEpsilon)
//            {
//                // compute dwijdu
//                Matrix< DDRMat > tdWijdu;
//                this->compute_dwijdu( aDofTypes, tdWijdu );
//
//                // compute dsdu
//                adsdu += 2.0 * trans( tWij ) * tdWijdu / tS;
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_chi()
//        {
//            // get the residual dof FI (here viscosity)
//            Field_Interpolator * tFIViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//            // get the density and gravity properties
//            std::shared_ptr< Property > & tPropViscosity =
//                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );
//
//            // compute chi
//            return tFIViscosity->val()( 0 ) / tPropViscosity->val()( 0 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dchidu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adchidu )
//        {
//            // get the derivative dof FIs
//            Field_Interpolator * tDerFI =
//                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
//
//            // init adchidu
//            adchidu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );
//
//            // get the residual dof FI (here viscosity)
//            Field_Interpolator * tFIViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//            // get the density and gravity properties
//            std::shared_ptr< Property > & tPropViscosity =
//                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );
//
//            // compute chi
//            real tChi = tFIViscosity->val()( 0 ) / tPropViscosity->val()( 0 );
//
//            // if dof type is viscosity
//            if( aDofTypes( 0 ) == mMasterDofViscosity )
//            {
//                adchidu += tDerFI->N() / tPropViscosity->val()( 0 );
//            }
//
//            // if viscosity property depends on dof type
//            if( tPropViscosity->check_dof_dependency( aDofTypes ) )
//            {
//                adchidu -= tChi * tPropViscosity->dPropdDOF( aDofTypes ) / tPropViscosity->val()( 0 );
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_fv1()
//        {
//            // compute chi
//            real tChi = this->compute_chi();
//
//            // compute fv1
//            real tFv1 = std::pow( tChi, 3.0 ) / ( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ) );
//
//            return tFv1;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dfv1du(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adfv1du )
//        {
//            // compute chi
//            real tChi = this->compute_chi();
//
//            // compute dchidu
//            Matrix< DDRMat > tdchidu;
//            this->compute_dchidu( aDofTypes, tdchidu );
//
//            // compute adfv1du
//            adfv1du = 3.0 * std::pow( mCv1, 3.0 ) * std::pow( tChi, 2.0 ) * tdchidu /
//                    std::pow( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ), 2.0 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_fv2()
//        {
//            // compute chi
//            real tChi = this->compute_chi();
//
//            // compute fv1
//            real tFv1 = this->compute_fv1();
//
//            // compute fv2
//            return 1.0 - tChi / ( 1.0 + tChi * tFv1 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dfv2du(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adfv2du )
//        {
//            // compute chi
//            real tChi = this->compute_chi();
//
//            // compute dchidu
//            Matrix< DDRMat > tdchidu;
//            this->compute_dchidu( aDofTypes, tdchidu );
//
//            // compute fv1
//            real tFv1 = this->compute_fv1();
//
//            // compute dfv1du
//            Matrix< DDRMat > tdfv1du;
//            this->compute_dfv1du( aDofTypes, tdfv1du );
//
//            // compute adfv2du
//            adfv2du = ( std::pow( tChi, 2.0 ) * tdfv1du - tdchidu ) / ( std::pow( 1.0 + tChi * tFv1, 2.0 ) );
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_fn()
//        {
//            // compute chi, chi³
//            real tChi  = this->compute_chi();
//            real tChi3 = std::pow( tChi, 3 );
//
//            // compute fn
//            return ( mCn1 + tChi3 ) / ( mCn1 - tChi3 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dfndu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adfndu )
//        {
//            // compute chi, chi², chi³
//            real tChi = this->compute_chi();
//            real tChi2 = std::pow( tChi, 2 );
//            real tChi3 = std::pow( tChi, 3 );
//
//            // compute dchidu
//            Matrix< DDRMat > tdchidu;
//            this->compute_dchidu( aDofTypes, tdchidu );
//
//            // compute adfndu
//            adfndu = 6.0 * mCn1 * tChi2 * tdchidu / std::pow( mCn1 - tChi3, 2 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_sbar()
//        {
//            // get the viscosity FI
//            Field_Interpolator * tFIViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//            // get the wall distance property
//            std::shared_ptr< Property > & tPropWallDistance =
//                    mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );
//
//            // get the wall distance value
//            real tWallDistance = tPropWallDistance->val()( 0 );
//
////            // check negative/zero wall distance
////            MORIS_ERROR( tWallDistance > 0.0,
////                    "SP_SUPG_Spalart_Allmaras_Turbulence::compute_sbar - %s",
////                    "Negative or zero wall distance.\n");
//
//            // threshold wall distance
//            tWallDistance = std::max(tWallDistance,mEpsilon);
//
//            // compute fv2
//            real tFv2 = this->compute_fv2();
//
//            // compute s
//            return tFv2 * tFIViscosity->val()( 0 ) / std::pow( mKappa * tWallDistance, 2.0 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dsbardu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >             & adsbardu )
//        {
//            // get the derivative dof FIs
//            Field_Interpolator * tFIDer =
//                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
//
//            // init dsbardu
//            adsbardu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//
//            // get the viscosity FI
//            Field_Interpolator * tFIViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//            // get the wall distance properties
//            std::shared_ptr< Property > & tPropWallDistance =
//                    mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );
//
//            // get the wall distance value
//            real tWallDistance = tPropWallDistance->val()( 0 );
//
////            // check negative/zero wall distance
////            MORIS_ERROR( tWallDistance > 0.0,
////                    "SP_SUPG_Spalart_Allmaras_Turbulence::compute_dsbardu - Negative or zero wall distance, exiting!");
//
//            // threshold wall distance
//            tWallDistance = std::max(tWallDistance,mEpsilon);
//
//            // compute fv2
//            real tFv2 = this->compute_fv2();
//
//            // compute dfv2du
//            Matrix< DDRMat > tdfv2du;
//            this->compute_dfv2du( aDofTypes, tdfv2du );
//
//            // compute dsbardu
//            adsbardu +=
//                    tFIViscosity->val() * tdfv2du /
//                    std::pow( mKappa * tWallDistance, 2.0 );
//
//            // if dof type is residual dof type
//            if( aDofTypes( 0 ) == mMasterDofViscosity )
//            {
//                // add contribution
//                adsbardu +=
//                        tFv2 * tFIViscosity->N() /
//                        std::pow( mKappa * tWallDistance, 2.0 );
//            }
//
//            // if wall distance depends on derivative dof type
//            if( ( tPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
//                    ( tWallDistance > mEpsilon ) )
//            {
//                //MORIS_LOG( "compute_dsbardu - Wall distance depends on dof type - not tested!" );
//
//                // add contribution to dsbardu
//                adsbardu -=
//                        2.0 * tFv2 * tFIViscosity->val()( 0 ) * tPropWallDistance->dPropdDOF( aDofTypes ) /
//                        ( std::pow( mKappa, 2.0 ) * std::pow( tWallDistance, 3.0 ) );
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_smod()
//        {
//            // compute Sbar
//            real tS = this->compute_s();
//
//            // compute S
//            real tSBar = this->compute_sbar();
//
//            // compute s
//            real tSMod = tS * ( std::pow( mCv2, 2 ) * tS + mCv3 * tSBar ) /
//                    ( ( mCv3 - 2.0 * mCv2 ) * tS - tSBar );
//
//            return tSMod;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dsmoddu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adsmoddu )
//        {
//            // compute S
//            real tS = this->compute_s();
//
//            // compute dsbardu
//            Matrix< DDRMat > tdsdu;
//            this->compute_dsdu( aDofTypes, tdsdu );
//
//            // compute SBar
//            real tSBar = this->compute_sbar();
//
//            // compute dsdu
//            Matrix< DDRMat > tdsbardu;
//            this->compute_dsbardu( aDofTypes, tdsbardu );
//
//            // compute smod num
//            real tSModNum = tS * ( std::pow( mCv2, 2 ) * tS + mCv3 * tSBar );
//
//            // compute smod deno
//            real tSModDeno = ( mCv3 - 2.0 * mCv2 ) * tS - tSBar;
//
//            // compute dsmoddu
//            adsmoddu = ( ( tdsdu * ( std::pow( mCv2, 2 ) * tS + mCv3 * tSBar ) +
//                    tS * ( std::pow( mCv2, 2 ) * tdsdu + mCv3 * tdsbardu ) ) * tSModDeno -
//                    tSModNum * ( ( mCv3 - 2.0 * mCv2 ) * tdsdu - tdsbardu ) ) /
//                            std::pow( tSModDeno, 2 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_stilde()
//        {
//            // compute S
//            real tS = this->compute_s();
//
//            // compute SBar
//            real tSBar = this->compute_sbar();
//
//            // init stilde
//            real tSTilde = tS;
//
//            // compute STilde
//            if( tSBar >= - mCv2 * tS )
//            {
//                tSTilde += tSBar;
//            }
//            else
//            {
//                // compute sMod
//                tSTilde += this->compute_smod();
//            }
//
//            // return STilde
//            return tSTilde;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dstildedu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adstildedu )
//        {
//            // compute SBar
//            real tS = this->compute_s();
//
//            // compute dSBardu
//            Matrix< DDRMat > tdSdu;
//            this->compute_dsdu( aDofTypes, tdSdu );
//
//            // compute SBar
//            real tSBar = this->compute_sbar();
//
//            // init dstildedu
//            adstildedu = tdSdu;
//
//            // compute dstildedu
//            if( tSBar >= - mCv2 * tS )
//            {
//                // compute dSdu
//                Matrix< DDRMat > tdSBardu;
//                this->compute_dsbardu( aDofTypes, tdSBardu );
//
//                // add dsdu
//                adstildedu += tdSBardu;
//            }
//            else
//            {
//                // compute dSModdu
//                Matrix< DDRMat > tdSModdu;
//                this->compute_dsmoddu( aDofTypes, tdSModdu );
//
//                // compute sMod
//                adstildedu += tdSModdu;
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_ft2()
//        {
//            // compute chi
//            real tChi = this->compute_chi();
//
//            // compute ft2
//            return mCt3 * std::exp( - mCt4 * std::pow( tChi, 2.0 ) );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dft2du(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adft2du )
//        {
//            // compute ft2
//            real tFt2 = this->compute_ft2();
//
//            // compute chi
//            real tChi = this->compute_chi();
//
//            // compute dChidu
//            Matrix< DDRMat > tdChidu;
//            this->compute_dchidu( aDofTypes, tdChidu );
//
//            // compute dft2du
//            adft2du = - mCt4 * tFt2 * 2.0 * tChi * tdChidu;
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_r()
//        {
//            // compute stilde
//            real tSTilde = this->compute_stilde();
//
//            // get the residual dof FI (here viscosity)
//            Field_Interpolator * tFIViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//            // get the wall distance property
//            std::shared_ptr< Property > & tPropWallDistance =
//                    mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );
//
//            // get the wall distance value
//            real tWallDistance = tPropWallDistance->val()( 0 );
//
////            // check negative/zero wall distance
////            MORIS_ERROR( tWallDistance > 0.0,
////                    "SP_SUPG_Spalart_Allmaras_Turbulence::compute_r - Negative or zero wall distance, exiting!");
//
//            // threshold wall distance
//            tWallDistance = std::max(tWallDistance,mEpsilon);
//
//            // compute viscosity / ( stilde * kappa² * d² )
//            real tR = tFIViscosity->val()( 0 ) / ( tSTilde * std::pow( mKappa * tWallDistance, 2.0 ) );
//
//            // check that r is finite and greater then zero or set it to mRLim
//            Matrix<DDRMat> tRMatrix( 1, 1, tR );
//            Matrix<DDRMat> tInvRMatrix( 1, 1, 1/tR );
//            if( !isfinite( tRMatrix ) || !isfinite(tInvRMatrix) )
//            {
//                tR = mRLim;
//            }
//
//            return std::min( tR, mRLim );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_drdu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adrdu )
//        {
//            // get the derivative dof FIs
//            Field_Interpolator * tDerFI =
//                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
//
//            // init adrdu
//            adrdu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );
//
//            // compute r
//            real tR = this->compute_r();
//
//            // if r < 10
//            if( tR < mRLim )
//            {
//                // get the residual dof FI (here viscosity)
//                Field_Interpolator * tFIViscosity =
//                        mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//                // get the wall distance property
//                std::shared_ptr< Property > & tPropWallDistance =
//                        mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );
//
//                // get the wall distance value
//                real tWallDistance = tPropWallDistance->val()( 0 );
//
////                // check negative/zero wall distance
////                MORIS_ERROR( tWallDistance > 0.0,
////                        "SP_SUPG_Spalart_Allmaras_Turbulence::compute_drdu - Negative or zero wall distance, exiting!");
//
//                // threshold wall distance
//                tWallDistance = std::max(tWallDistance,mEpsilon);
//
//                // compute stilde
//                real tSTilde = this->compute_stilde();
//
//                // compute dStildedu
//                Matrix< DDRMat > tdSTildedu;
//                this->compute_dstildedu( aDofTypes, tdSTildedu );
//
//                // add contribution from dStildedu
//                adrdu -= tFIViscosity->val() * tdSTildedu / std::pow( tSTilde * mKappa * tWallDistance, 2.0 );
//
//                // if dof type is viscosity
//                if( aDofTypes( 0 ) == mMasterDofViscosity )
//                {
//                    // add contribution from viscosity
//                    adrdu += tSTilde * tDerFI->N() / std::pow( tSTilde * mKappa * tWallDistance, 2.0 );
//                }
//
//                // if wall distance depends on derivative dof type
//                if( ( tPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
//                        ( tWallDistance > mEpsilon ) )
//                {
//                    //MORIS_LOG( "compute_drdu - Wall distance depends on dof type - not tested!" );
//
//                    // add contribution from wall distance
//                    adrdu -= 2.0 * tFIViscosity->val()( 0 ) * tPropWallDistance->dPropdDOF( aDofTypes ) /
//                            ( tSTilde * std::pow( mKappa, 2.0 ) * std::pow( tWallDistance, 3.0 ) );
//                }
//            }
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_g()
//        {
//            // compute r
//            real tR = this->compute_r();
//
//            // compute g
//            return tR + mCw2 * ( std::pow( tR, 6.0 ) - tR );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dgdu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adgdu )
//        {
//            // compute r
//            real tR = this->compute_r();
//
//            // compute drdu
//            Matrix< DDRMat > tdrdu;
//            this->compute_drdu( aDofTypes, tdrdu );
//
//            // compute adgdu
//            adgdu = ( 1.0 + mCw2 * ( 6.0 * std::pow( tR, 5.0 ) - 1.0 ) ) * tdrdu;
//        }
//
//        //------------------------------------------------------------------------------
//
//        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_fw()
//        {
//            // compute g
//            real tG = this->compute_g();
//
//            // compute fw
//            real tFw = ( 1.0 + std::pow( mCw3, 6.0 ) ) / ( std::pow( tG, 6.0 ) + std::pow( mCw3, 6.0 ) );
//            tFw = tG * std::pow( tFw, 1.0 / 6.0 );
//
//            return tFw;
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dfwdu(
//                const moris::Cell< MSI::Dof_Type > & aDofTypes,
//                Matrix< DDRMat >                   & adfwdu )
//        {
//            // compute g
//            real tG = this->compute_g();
//
//            // compute dgdu
//            Matrix< DDRMat > tdgdu;
//            this->compute_dgdu( aDofTypes, tdgdu );
//
//            // compute fw
//            real tFw = this->compute_fw();
//
//            // init adfwdu
//            adfwdu = ( tFw * std::pow( mCw3, 6.0 ) * tdgdu ) / ( tG * ( std::pow( tG, 6.0 ) + std::pow( mCw3, 6.0 ) ) );
//        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


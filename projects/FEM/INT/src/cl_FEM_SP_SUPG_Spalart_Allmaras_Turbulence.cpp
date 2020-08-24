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

        //------------------------------------------------------------------------------

        SP_SUPG_Spalart_Allmaras_Turbulence::SP_SUPG_Spalart_Allmaras_Turbulence()
        {
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Viscosity" ]    = SP_Property_Type::VISCOSITY;
            mPropertyMap[ "WallDistance" ] = SP_Property_Type::WALL_DISTANCE;
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
                            // create error message
                            std::string tErrMsg =
                                    std::string( "SP_SUPG_Spalart_Allmaras_Turbulence::set_dof_type_list - Unknown aDofString : ") +
                                    tDofString;
                            MORIS_ERROR( false , tErrMsg.c_str() );
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

        void SP_SUPG_Spalart_Allmaras_Turbulence::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "SP_SUPG_Spalart_Allmaras_Turbulence::set_property - Unknown aPropertyString: " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(), tErrMsg.c_str() );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::set_function_pointers()
        {
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::reset_cluster_measures()
        {
            // evaluate element size from the cluster
            mElementSize = mCluster->compute_cluster_cell_length_measure(
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::MASTER );
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::eval_SP()
        {
            // set size for SP values
            mPPVal.set_size( 1, 1, 0.0 );

            // get the velocity FI
            Field_Interpolator * tViscosityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the velocity FI
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // compute uTilde = u - cb2/sigma gradv
            Matrix< DDRMat > tModVelocity =
                    tVelocityFI->val() - mCb2 * tViscosityFI->gradx( 1 ) / mSigma;

            // compute norm( uTilde )
            real tNormA = std::sqrt( dot( tModVelocity, tModVelocity ) );

            // compute the diffusion coefficient
            real tK = this->compute_diffusion_coefficient();

            // compute the source
            real tSource = this->compute_production_coefficient() +
                    this->compute_wall_destruction_coefficient();

            // evaluate tau
            real tTau =
                    std::pow( 2.0 * tNormA / mElementSize, 2 ) +
                    std::pow( 4.0 * tK / std::pow( mElementSize, 2.0 ), 2 ) +
                    std::pow( tSource, 2 );

            // set tau
            mPPVal = {{ std::pow( tTau, -0.5 ) }};
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the velocity FI
            Field_Interpolator * tViscosityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the velocity FI
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get viscosity property
            std::shared_ptr< Property > tPropViscosity =
                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );

            // evaluate uTilde = u - 1/sigma cb2 gradv
            Matrix< DDRMat > tModVelocity =
                    tVelocityFI->val() - mCb2 * tViscosityFI->gradx( 1 ) / mSigma;

            // evaluate norm( uTilde )
            real tNormA = std::sqrt( dot( tModVelocity, tModVelocity ) );

            // tau A
            real tTauA = 2.0 * tNormA / mElementSize;

            // tau K
            real tTauK = 4.0 * this->compute_diffusion_coefficient() / std::pow( mElementSize, 2.0 );

            // tau S
            real tTauS = this->compute_production_coefficient() +
                    this->compute_wall_destruction_coefficient();

            // init derivative of tauA, tauK, tauS
            Matrix< DDRMat > tdtauAdu( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );
            Matrix< DDRMat > tdtauKdu( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );
            Matrix< DDRMat > tdtauSdu( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // if dof type is velocity
            if( aDofTypes( 0 ) == mMasterDofVelocity )
            {
                // add contribution to dSPdu
                tdtauAdu.matrix_data() +=
                        2.0 * trans( tModVelocity ) * tVelocityFI->N() / ( mElementSize * tNormA );
            }

            // if dof type is velocity
            if( aDofTypes( 0 ) == mMasterDofViscosity )
            {
                // evaluate dadu
                tdtauAdu.matrix_data() -=
                        2.0 * trans( tModVelocity ) * ( mCb2 * tViscosityFI->dnNdxn( 1 ) / mSigma ) /
                        ( mElementSize * tNormA );
            }

            // compute tdtauKdu
            Matrix< DDRMat > tddiffusiondu;
            this->compute_ddiffusiondu( aDofTypes, tddiffusiondu );
            tdtauKdu = 4.0 * tddiffusiondu / std::pow( mElementSize, 2.0 );

            // compute tdtauSdu
            Matrix< DDRMat > tdproductiondu;
            this->compute_dproductiondu( aDofTypes, tdproductiondu );
            Matrix< DDRMat > tdwalldestructiondu;
            this->compute_dwalldestructiondu( aDofTypes, tdwalldestructiondu );
            tdtauSdu = tdproductiondu + tdwalldestructiondu;

            // evaluate tau
            Matrix< DDRMat > tTau = this->val();

            // scale dSPdu
            mdPPdMasterDof( tDofIndex ) = - 0.5 * std::pow( tTau( 0 ), 3 ) *
                    ( 2.0 * tTauA * tdtauAdu + 2.0 * tTauK * tdtauKdu + 2.0 * tTauS * tdtauSdu );
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_production_coefficient()
        {
            // init production term
            real tProductionTerm = 0.0;

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity >= 0.0 )
            {
                // compute ft2
                real tFt2 = this->compute_ft2();

                // compute Stilde
                real tSTilde = this->compute_stilde();

                // compute production term
                tProductionTerm = mCb1 * ( 1.0 - tFt2 ) * tSTilde;
            }
            // if viscosity is negative
            else
            {
                // compute S
                real tS = this->compute_s();

                // compute production term
                tProductionTerm = mCb1 * ( 1.0 - mCt3 ) * tS;
            }

            return tProductionTerm;
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dproductiondu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adproductiondu )
        {
            // get derivative dof type
            MSI::Dof_Type tDerDofType = aDofTypes( 0 );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( tDerDofType );

            // init adproductiondu
            adproductiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity >= 0.0 )
            {
                // compute ft2
                real tFt2 = this->compute_ft2();

                // compute dft2du
                Matrix< DDRMat > tdft2du;
                this->compute_dft2du( aDofTypes, tdft2du );

                // compute STilde
                real tSTilde = this->compute_stilde();

                // compute dSTildedu
                Matrix< DDRMat > tdstildedu;
                this->compute_dstildedu( aDofTypes, tdstildedu );

                // add contribution to dproductiondu
                adproductiondu.matrix_data() += - mCb1 * tSTilde * tdft2du +
                        mCb1 * ( 1 - tFt2 ) * tdstildedu;
            }
            // if viscosity is negative
            else
            {
                // compute dsdu
                Matrix< DDRMat > tdsdu;
                this->compute_dsdu( aDofTypes, tdsdu );

                // add contribution to dproductiondu
                adproductiondu.matrix_data() +=
                        mCb1 * ( 1.0 - mCt3 ) * tdsdu;
            }
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_wall_destruction_coefficient()
        {
            // init wall destruction term
            real tWallDestructionTerm = 0.0;

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the wall distance property
            std::shared_ptr< Property > tPropWallDistance =
                    mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // get the wall distance value
            real tWallDistance = tPropWallDistance->val()( 0 );

            // check negative/zero wall distance
            MORIS_ERROR( tWallDistance > 0.0,
                    "SP_SUPG_Spalart_Allmaras_Turbulence::compute_wall_destruction_coefficient - Negative or zero wall distance, exiting!");

            // if viscosity is positive or zero
            if( tModViscosity >= 0.0 )
            {
                // compute fw
                real tFw = this->compute_fw();

                // compute ft2
                real tFt2 = this->compute_ft2();

                // compute wall destruction term
                tWallDestructionTerm =
                        ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
                        tModViscosity / std::pow( tWallDistance, 2.0 );
            }
            // if viscosity is negative
            else
            {
                // compute wall destruction term
                tWallDestructionTerm = - mCw1 * tModViscosity / std::pow( tWallDistance, 2.0 );
            }

            return tWallDestructionTerm;
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dwalldestructiondu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adwalldestructiondu )
        {
            // get derivative dof type
            MSI::Dof_Type tDerDofType = aDofTypes( 0 );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( tDerDofType );

            // init adwalldestructiondu
            adwalldestructiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the wall distance property
            std::shared_ptr< Property > tPropWallDistance =
                    mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // get the wall distance value
            real tWallDistance = tPropWallDistance->val()( 0 );

            // check negative/zero wall distance
            MORIS_ERROR( tWallDistance > 0.0,
                    "SP_SUPG_Spalart_Allmaras_Turbulence::compute_dwalldestructiondu - Negative or zero wall distance, exiting!");

            // if viscosity is positive or zero
            if( tModViscosity >= 0.0 )
            {
                // compute fw
                real tFw = this->compute_fw();

                // compute dfwfu
                Matrix< DDRMat > tdfwdu;
                this->compute_dfwdu( aDofTypes, tdfwdu );

                // compute ft2
                real tFt2 = this->compute_ft2();

                // compute dft2fu
                Matrix< DDRMat > tdft2du;
                this->compute_dft2du( aDofTypes, tdft2du );

                // if derivative dof type is viscosity
                if( tDerDofType == mMasterDofViscosity )
                {
                    // add contribution to dwalldestructiondu
                    adwalldestructiondu.matrix_data() +=
                            ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
                            tFIModViscosity->N() /
                            std::pow( tWallDistance, 2.0 );
                }

                // add contribution to dwalldestructiondu
                adwalldestructiondu.matrix_data() +=
                        ( mCw1 * tdfwdu - mCb1 * tdft2du / std::pow( mKappa, 2.0 ) ) *
                        tModViscosity / std::pow( tWallDistance, 2.0 );
            }
            // if viscosity is negative
            else
            {
                // if derivative dof type is viscosity
                if( tDerDofType == mMasterDofViscosity )
                {
                    // add contribution to dwalldestructiondu
                    adwalldestructiondu.matrix_data() -=
                            mCw1 * tFIModViscosity->N() /
                            std::pow( tWallDistance, 2.0 );
                }
            }
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_diffusion_coefficient()
        {
            // init diffusion coeff
            real tDiffusionTerm = 0.0;

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the wall distance property
            std::shared_ptr< Property > tPropKinViscosity =
                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // get the fluid kinematic viscosity value
            real tKinViscosity = tPropKinViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity >= 0.0 )
            {
                // compute diffusion term
                tDiffusionTerm = ( tKinViscosity + tModViscosity ) / mSigma;
            }
            // if viscosity is negative
            else
            {
                // compute fn
                real tFn = this->compute_fn();

                // compute diffusion term
                tDiffusionTerm = ( tKinViscosity + tModViscosity * tFn ) / mSigma;
            }

            return tDiffusionTerm;
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_ddiffusiondu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & addiffusiondu )
        {
            // get derivative dof type
            MSI::Dof_Type tDerDofType = aDofTypes( 0 );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( tDerDofType );

            // init ddiffusiondu
            addiffusiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the fluid  kinematic viscosity property
            std::shared_ptr< Property > tPropKinViscosity =
                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity >= 0.0 )
            {
                // if derivative dof type is viscosity
                if( tDerDofType == mMasterDofViscosity )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu.matrix_data() += tFIModViscosity->N() / mSigma;
                }

                // if kinematic viscosity depends on derivative dof type
                if( tPropKinViscosity->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu.matrix_data() += tPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
                }
            }
            // if viscosity is negative
            else
            {
                // compute fn
                real tFn = this->compute_fn();

                // compute dfndu
                Matrix< DDRMat > tdfndu;
                this->compute_dfndu( aDofTypes, tdfndu );

                // if derivative dof type is viscosity
                if( tDerDofType == mMasterDofViscosity )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu.matrix_data() += tFn * tFIModViscosity->N() / mSigma;
                }

                // if kinematic viscosity depends on derivative dof type
                if( tPropKinViscosity->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu.matrix_data() += tPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
                }

                // add contribution from fn to ddiffusiondu
                addiffusiondu.matrix_data() += tModViscosity * tdfndu / mSigma;
            }
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_wij( Matrix< DDRMat > & aWij )
        {
            // get the velocity FI
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get gradient of velocity
            Matrix< DDRMat > tGradVelocity = tFIVelocity->gradx( 1 );

            // switch on space dim
            switch ( mSpaceDim )
            {
                case 2 :
                {
                    // init aWij = [ w11 w12 w21 w22]
                    aWij.set_size( 4, 1, 0.0 );

                    // compute Wij
                    aWij( 1 ) = 0.5 * ( tGradVelocity( 1, 0 ) - tGradVelocity( 0, 1 ) );
                    aWij( 2 ) = 0.5 * ( tGradVelocity( 0, 1 ) - tGradVelocity( 1, 0 ) );
                    break;
                }
                case 3 :
                {
                    // init aWij = [ w11 w12 w13 w21 w22 w23 w31 w32 w33 ]
                    aWij.set_size( 9, 1, 0.0 );

                    // compute Wij
                    aWij( 1 ) = 0.5 * ( tGradVelocity( 1, 0 ) - tGradVelocity( 0, 1 ) );
                    aWij( 2 ) = 0.5 * ( tGradVelocity( 2, 0 ) - tGradVelocity( 0, 2 ) );
                    aWij( 3 ) = 0.5 * ( tGradVelocity( 0, 1 ) - tGradVelocity( 1, 0 ) );
                    aWij( 5 ) = 0.5 * ( tGradVelocity( 2, 1 ) - tGradVelocity( 1, 2 ) );
                    aWij( 6 ) = 0.5 * ( tGradVelocity( 0, 2 ) - tGradVelocity( 2, 0 ) );
                    aWij( 7 ) = 0.5 * ( tGradVelocity( 1, 2 ) - tGradVelocity( 2, 1 ) );
                    break;
                }
                default:
                    MORIS_ERROR( false, "SP_SUPG_Spalart_Allmaras_Turbulence::compute_wij - space dim can only be 2 or 3" );
            }
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dwijdu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adwijdu )
        {
            // get the der FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // switch on space dim
            switch ( mSpaceDim )
            {
                case 2 :
                {
                    // init aWij = [ w11 w12 w21 w22]
                    adwijdu.set_size( 4, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

                    if( aDofTypes( 0 ) == MSI::Dof_Type::VX )
                    {
                        // get gradient of velocity
                        Matrix< DDRMat > tdNdxVelocity = tFIVelocity->dnNdxn( 1 );

                        // get number of bases for displacement
                        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

                        // compute adwijdu
                        adwijdu( { 1, 1 }, { 0, tNumBases - 1 } )             =   0.5 * tdNdxVelocity.get_row( 1 );
                        adwijdu( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 0 );
                        adwijdu( { 2, 2 }, { 0, tNumBases - 1 } )             = - 0.5 * tdNdxVelocity.get_row( 1 );
                        adwijdu( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) =   0.5 * tdNdxVelocity.get_row( 0 );
                    }
                    break;
                }
                case 3 :
                {
                    // init aWij = [ w11 w12 w13 w21 w22 w23 w31 w32 w33 ]
                    adwijdu.set_size( 9, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

                    if( aDofTypes( 0 ) == MSI::Dof_Type::VX )
                    {
                        // get gradient of velocity
                        Matrix< DDRMat > tdNdxVelocity = tFIVelocity->dnNdxn( 1 );

                        // get number of bases for displacement
                        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

                        // compute adwijdu
                        adwijdu( { 1, 1 }, { 0, tNumBases - 1 } )                 =   0.5 * tdNdxVelocity.get_row( 1 );
                        adwijdu( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } )     = - 0.5 * tdNdxVelocity.get_row( 0 );

                        adwijdu( { 2, 2 }, { 0, tNumBases - 1 } )                 =   0.5 * tdNdxVelocity.get_row( 2 );
                        adwijdu( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 0 );

                        adwijdu( { 3, 3 }, { 0, tNumBases - 1 } )                 = - 0.5 * tdNdxVelocity.get_row( 1 );
                        adwijdu( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } )     =   0.5 * tdNdxVelocity.get_row( 0 );

                        adwijdu( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } )     =   0.5 * tdNdxVelocity.get_row( 2 );
                        adwijdu( { 5, 5 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 1 );

                        adwijdu( { 6, 6 }, { 0, tNumBases - 1 } )                 = - 0.5 * tdNdxVelocity.get_row( 2 );
                        adwijdu( { 6, 6 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =   0.5 * tdNdxVelocity.get_row( 0 );

                        adwijdu( { 7, 7 }, { tNumBases, 2 * tNumBases - 1 } )     = - 0.5 * tdNdxVelocity.get_row( 2 );
                        adwijdu( { 7, 7 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =   0.5 * tdNdxVelocity.get_row( 1 );
                    }
                    break;
                }
                default:
                    MORIS_ERROR( false, "SP_SUPG_Spalart_Allmaras_Turbulence::compute_dwijdu - space dim can only be 2 or 3" );
            }
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_s()
        {
            // compute wij
            Matrix< DDRMat > tWij;
            this->compute_wij( tWij );

            // compute WijWij
            Matrix< DDRMat > tWijWij = trans( tWij ) * tWij;

            // compute s
            return std::sqrt( 2.0 * tWijWij( 0 ) );
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dsdu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adsdu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init dsdu
            adsdu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // compute sbar
            real tS = this->compute_s();

            // if s is greater than zero
            if( tS > 0.0 )
            {
                // compute wij
                Matrix< DDRMat > tWij;
                this->compute_wij( tWij );

                // compute dwijdu
                Matrix< DDRMat > tdWijdu;
                this->compute_dwijdu( aDofTypes, tdWijdu );

                // compute dsdu
                adsdu.matrix_data() += 2.0 * trans( tWij ) * tdWijdu / tS;
            }
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_chi()
        {
            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropViscosity =
                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );

            // compute chi
            return tFIViscosity->val()( 0 ) / tPropViscosity->val()( 0 );
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dchidu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adchidu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adchidu
            adchidu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropViscosity =
                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );

            // compute chi
            real tChi = tFIViscosity->val()( 0 ) / tPropViscosity->val()( 0 );

            // if dof type is viscosity
            if( aDofTypes( 0 ) == mMasterDofViscosity )
            {
                adchidu.matrix_data() += tDerFI->N() / tPropViscosity->val()( 0 );
            }

            // if viscosity property depends on dof type
            if( tPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                adchidu.matrix_data() -= tChi * tPropViscosity->dPropdDOF( aDofTypes ) / tPropViscosity->val()( 0 );
            }
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_fv1()
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute fv1
            real tFv1 = std::pow( tChi, 3.0 ) / ( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ) );

            return tFv1;
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dfv1du(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfv1du )
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            this->compute_dchidu( aDofTypes, tdchidu );

            // compute adfv1du
            adfv1du = 3.0 * std::pow( mCv1, 3.0 ) * std::pow( tChi, 2.0 ) * tdchidu /
                    std::pow( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ), 2.0 );
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_fv2()
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute fv1
            real tFv1 = this->compute_fv1();

            // compute fv2
            return 1.0 - tChi / ( 1.0 + tChi * tFv1 );
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dfv2du(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfv2du )
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            this->compute_dchidu( aDofTypes, tdchidu );

            // compute fv1
            real tFv1 = this->compute_fv1();

            // compute dfv1du
            Matrix< DDRMat > tdfv1du;
            this->compute_dfv1du( aDofTypes, tdfv1du );

            // compute adfv2du
            adfv2du = ( std::pow( tChi, 2.0 ) * tdfv1du - tdchidu ) / ( std::pow( 1.0 + tChi * tFv1, 2.0 ) );
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_fn()
        {
            // compute chi, chi³
            real tChi  = this->compute_chi();
            real tChi3 = std::pow( tChi, 3 );

            // compute fn
            return ( mCn1 + tChi3 ) / ( mCn1 - tChi3 );
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dfndu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfndu )
        {
            // compute chi, chi², chi³
            real tChi = this->compute_chi();
            real tChi2 = std::pow( tChi, 2 );
            real tChi3 = std::pow( tChi, 3 );

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            this->compute_dchidu( aDofTypes, tdchidu );

            // compute adfndu
            adfndu = 6.0 * mCn1 * tChi2 * tdchidu / std::pow( mCn1 - tChi3, 2 );
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_sbar()
        {
            // get the viscosity FI
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the wall distance property
            std::shared_ptr< Property > tPropWallDistance =
                    mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );

            // get the wall distance value
            real tWallDistance = tPropWallDistance->val()( 0 );

            // check negative/zero wall distance
            MORIS_ERROR( tWallDistance > 0.0,
                    "SP_SUPG_Spalart_Allmaras_Turbulence::compute_sbar - Negative or zero wall distance, exiting!");

            // compute fv2
            real tFv2 = this->compute_fv2();

            // compute s
            return tFv2 * tFIViscosity->val()( 0 ) / std::pow( mKappa * tWallDistance, 2.0 );
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dsbardu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >             & adsbardu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init dsbardu
            adsbardu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity FI
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the wall distance properties
            std::shared_ptr< Property > tPropWallDistance =
                    mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );

            // get the wall distance value
            real tWallDistance = tPropWallDistance->val()( 0 );

            // check negative/zero wall distance
            MORIS_ERROR( tWallDistance > 0.0,
                    "SP_SUPG_Spalart_Allmaras_Turbulence::compute_dsbardu - Negative or zero wall distance, exiting!");

            // compute fv2
            real tFv2 = this->compute_fv2();

            // compute dfv2du
            Matrix< DDRMat > tdfv2du;
            this->compute_dfv2du( aDofTypes, tdfv2du );

            // compute dsbardu
            adsbardu.matrix_data() +=
                    tFIViscosity->val() * tdfv2du /
                    std::pow( mKappa * tWallDistance, 2.0 );

            // if dof type is residual dof type
            if( aDofTypes( 0 ) == mMasterDofViscosity )
            {
                // add contribution
                adsbardu.matrix_data() +=
                        tFv2 * tFIViscosity->N() /
                        std::pow( mKappa * tWallDistance, 2.0 );
            }
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_smod()
        {
            // compute Sbar
            real tS = this->compute_s();

            // compute S
            real tSBar = this->compute_sbar();

            // compute s
            real tSMod = tS * ( std::pow( mCv2, 2 ) * tS + mCv3 * tSBar ) /
                    ( ( mCv3 - 2.0 * mCv2 ) * tS - tSBar );

            return tSMod;
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dsmoddu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adsmoddu )
        {
            // compute S
            real tS = this->compute_s();

            // compute dsbardu
            Matrix< DDRMat > tdsdu;
            this->compute_dsdu( aDofTypes, tdsdu );

            // compute SBar
            real tSBar = this->compute_sbar();

            // compute dsdu
            Matrix< DDRMat > tdsbardu;
            this->compute_dsbardu( aDofTypes, tdsbardu );

            // compute smod num
            real tSModNum = tS * ( std::pow( mCv2, 2 ) * tS + mCv3 * tSBar );

            // compute smod deno
            real tSModDeno = ( mCv3 - 2.0 * mCv2 ) * tS - tSBar;

            // compute dsmoddu
            adsmoddu = ( ( tdsdu * ( std::pow( mCv2, 2 ) * tS + mCv3 * tSBar ) +
                    tS * ( std::pow( mCv2, 2 ) * tdsdu + mCv3 * tdsbardu ) ) * tSModDeno -
                    tSModNum * ( ( mCv3 - 2.0 * mCv2 ) * tdsdu - tdsbardu ) ) /
                            std::pow( tSModDeno, 2 );
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_stilde()
        {
            // compute S
            real tS = this->compute_s();

            // compute SBar
            real tSBar = this->compute_sbar();

            // init stilde
            real tSTilde = tS;

            // compute STilde
            if( tSBar >= - mCv2 * tS )
            {
                tSTilde += tSBar;
            }
            else
            {
                // compute sMod
                tSTilde += this->compute_smod();
            }

            // return STilde
            return tSTilde;
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dstildedu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adstildedu )
        {
            // compute SBar
            real tS = this->compute_s();

            // compute dSBardu
            Matrix< DDRMat > tdSdu;
            this->compute_dsdu( aDofTypes, tdSdu );

            // compute SBar
            real tSBar = this->compute_sbar();

            // init dstildedu
            adstildedu = tdSdu;

            // compute dstildedu
            if( tSBar >= - mCv2 * tS )
            {
                // compute dSdu
                Matrix< DDRMat > tdSBardu;
                this->compute_dsbardu( aDofTypes, tdSBardu );

                // add dsdu
                adstildedu.matrix_data() += tdSBardu.matrix_data();
            }
            else
            {
                // compute dSModdu
                Matrix< DDRMat > tdSModdu;
                this->compute_dsmoddu( aDofTypes, tdSModdu );

                // compute sMod
                adstildedu.matrix_data() += tdSModdu.matrix_data();
            }
        }

        //------------------------------------------------------------------------------
        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_ft2()
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute ft2
            return mCt3 * std::exp( - mCt4 * std::pow( tChi, 2.0 ) );
        }

        //------------------------------------------------------------------------------
        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dft2du(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adft2du )
        {
            // compute ft2
            real tFt2 = this->compute_ft2();

            // compute chi
            real tChi = this->compute_chi();

            // compute dChidu
            Matrix< DDRMat > tdChidu;
            this->compute_dchidu( aDofTypes, tdChidu );

            // compute dft2du
            adft2du = - mCt4 * tFt2 * 2.0 * tChi * tdChidu;
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_r()
        {
            // compute stilde
            real tSTilde = this->compute_stilde();

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the wall distance property
            std::shared_ptr< Property > tPropWallDistance =
                    mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );

            // get the wall distance value
            real tWallDistance = tPropWallDistance->val()( 0 );

            // check negative/zero wall distance
            MORIS_ERROR( tWallDistance > 0.0,
                    "SP_SUPG_Spalart_Allmaras_Turbulence::compute_r - Negative or zero wall distance, exiting!");

            // compute viscosity / ( stilde * kappa² * d² )
            real tR = tFIViscosity->val()( 0 ) / ( tSTilde * std::pow( mKappa * tWallDistance, 2.0 ) );

            // check that r is finite or set it to mRLim
            Matrix<DDRMat> tRMatrix( 1, 1, tR );
            if( !isfinite( tRMatrix ) )
            {
                tR = mRLim;
            }

            return std::min( tR, mRLim );
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_drdu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adrdu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adrdu
            adrdu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // compute r
            real tR = this->compute_r();

            // if r < 10
            if( tR < mRLim )
            {
                // get the residual dof FI (here viscosity)
                Field_Interpolator * tFIViscosity =
                        mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

                // get the wall distance property
                std::shared_ptr< Property > tPropWallDistance =
                        mMasterProp( static_cast< uint >( SP_Property_Type::WALL_DISTANCE ) );

                // get the wall distance value
                real tWallDistance = tPropWallDistance->val()( 0 );

                // check negative/zero wall distance
                MORIS_ERROR( tWallDistance > 0.0,
                        "SP_SUPG_Spalart_Allmaras_Turbulence::compute_drdu - Negative or zero wall distance, exiting!");

                // compute stilde
                real tSTilde = this->compute_stilde();

                // compute dStildedu
                Matrix< DDRMat > tdSTildedu;
                this->compute_dstildedu( aDofTypes, tdSTildedu );

                // add contribution from dStildedu
                adrdu.matrix_data() -= tFIViscosity->val() * tdSTildedu;

                // if dof type is viscosity
                if( aDofTypes( 0 ) == mMasterDofViscosity )
                {
                    // add contribution from viscosity
                    adrdu.matrix_data() += tSTilde * tDerFI->N().matrix_data();
                }

                adrdu = adrdu / std::pow( tSTilde * mKappa * tWallDistance, 2.0 );
            }
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_g()
        {
            // compute r
            real tR = this->compute_r();

            // compute g
            return tR + mCw2 * ( std::pow( tR, 6.0 ) - tR );
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dgdu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adgdu )
        {
            // compute r
            real tR = this->compute_r();

            // compute drdu
            Matrix< DDRMat > tdrdu;
            this->compute_drdu( aDofTypes, tdrdu );

            // compute adgdu
            adgdu = ( 1.0 + mCw2 * ( 6.0 * std::pow( tR, 5.0 ) - 1.0 ) ) * tdrdu;
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::compute_fw()
        {
            // compute g
            real tG = this->compute_g();

            // compute fw
            real tFw = ( 1.0 + std::pow( mCw3, 6.0 ) ) / ( std::pow( tG, 6.0 ) + std::pow( mCw3, 6.0 ) );
            tFw = tG * std::pow( tFw, 1.0 / 6.0 );

            return tFw;
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::compute_dfwdu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfwdu )
        {
            // compute g
            real tG = this->compute_g();

            // compute dgdu
            Matrix< DDRMat > tdgdu;
            this->compute_dgdu( aDofTypes, tdgdu );

            // compute fw
            real tFw = this->compute_fw();

            // init adfwdu
            adfwdu = ( tFw * std::pow( mCw3, 6.0 ) * tdgdu ) / ( tG * ( std::pow( tG, 6.0 ) + std::pow( mCw3, 6.0 ) ) );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */



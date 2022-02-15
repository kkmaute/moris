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

        void SP_SUPG_Spalart_Allmaras_Turbulence::build_global_dof_type_list()
        {
            // call parent implementation
            Stabilization_Parameter::build_global_dof_type_list();

            // get number of dof types
            uint tNumMasterGlobalDofTypes = mMasterGlobalDofTypes.size();

            // init child specific eval flags
            mdLengthScaledMasterDofEval.set_size( tNumMasterGlobalDofTypes, 1, true );

            // init child specific storage
            mdLengthScaledMasterDof.resize( tNumMasterGlobalDofTypes );
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::set_function_pointers(){}

        //------------------------------------------------------------------------------
        /**
         * reset evaluation flags
         */
        void SP_SUPG_Spalart_Allmaras_Turbulence::reset_eval_flags()
        {
            // call parent implementation
            Stabilization_Parameter::reset_eval_flags();

            // reset child specific eval flags for chi
            mLengthScaleEval = true;
            mdLengthScaledMasterDofEval.fill( true );
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Spalart_Allmaras_Turbulence::eval_SP()
        {
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

            // compute the length scale
            real tHugn = this->length_scale();

            // compute uTilde = u - cb2/sigma gradv
            Matrix< DDRMat > tModVelocity =
                    tVelocityFI->val() - mCb2 * tViscosityFI->gradx( 1 ) / mSigma;

            // compute norm( uTilde ) and threshold tNormA (for consistency with derivative computation)
            real tNormA = std::max( norm( tModVelocity ), mEpsilon);

            // tau A
            real tTauA = 2.0 * tNormA / tHugn;

            // tau K
            real tTauK = 4.0 * tCMSATurbulence->diffusion_coefficient()( 0 ) / std::pow( tHugn, 2.0 );

            // tau S
            real tTauS =
                    tCMSATurbulence->production_coefficient()( 0 ) -
                    tCMSATurbulence->wall_destruction_coefficient()( 0 );

            // evaluate tau
            real tTau = std::pow( tTauA, 2 ) + std::pow( tTauK, 2 ) + std::pow( tTauS, 2 );

            // threshold tau
            tTau = std::max( tTau, mEpsilon );

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

            // compute the length scale
            real tHugn = this->length_scale();

            // evaluate uTilde = u - 1/sigma cb2 gradv
            Matrix< DDRMat > tModVelocity =
                    tVelocityFI->val() - mCb2 * tViscosityFI->gradx( 1 ) / mSigma;

            // compute norm( uTilde ) and threshold tNormA (for consistency with derivative computation)
            real tNormA = std::max( norm( tModVelocity ), mEpsilon);

            // tau A
            real tTauA = 2.0 * tNormA / tHugn;

            // tau K
            real tTauK = 4.0 * tCMSATurbulence->diffusion_coefficient()( 0 ) / std::pow( tHugn, 2.0 );

            // tau S
            real tTauS =
                    tCMSATurbulence->production_coefficient()( 0 ) -
                    tCMSATurbulence->wall_destruction_coefficient()( 0 );

            // evaluate tau
            real tTau = std::max( std::pow( tTauA, 2 ) + std::pow( tTauK, 2 ) + std::pow( tTauS, 2 ), mEpsilon );

            // compute dSPdu
            if ( tTau > mEpsilon )
            {
                // add contribution from the length scale derivative
                mdPPdMasterDof( tDofIndex ) = - this->dlengthscaledmasteru( aDofTypes ) *
                        ( tTauA * 2.0 * tNormA / std::pow( tHugn, 2.0 ) +
                          tTauK * 8.0 * tCMSATurbulence->diffusion_coefficient()( 0 )  / std::pow( tHugn, 3.0 ) );

                // if normA greater than threshold
                if( tNormA > mEpsilon )
                {
                    // if dof type is velocity
                    if( aDofTypes( 0 ) == mMasterDofVelocity )
                    {
                        // add contribution to mdPPdMasterDof
                        mdPPdMasterDof( tDofIndex ) += tTauA * 2.0 *
                                trans( tModVelocity ) * tVelocityFI->N() /
                                ( tHugn * tNormA );
                    }
                    // if dof type is viscosity
                    else if( aDofTypes( 0 ) == mMasterDofViscosity )
                    {
                        // add contribution to mdPPdMasterDof
                        mdPPdMasterDof( tDofIndex ) -= tTauA * 2.0 *
                                mCb2 * trans( tModVelocity ) * tViscosityFI->dnNdxn( 1 ) / mSigma /
                                ( tHugn * tNormA );
                    }
                }

                // if turbulence CM depends on dof type
                if( tCMSATurbulence->check_dof_dependency( aDofTypes ) )
                {
                    // compute tdtauKdu
                    mdPPdMasterDof( tDofIndex ) +=
                            tTauK * 4.0 * tCMSATurbulence->ddiffusioncoeffdu( aDofTypes ) /
                            std::pow( tHugn, 2.0 );

                    // compute tdtauSdu
                    mdPPdMasterDof( tDofIndex ) += tTauS * (
                            tCMSATurbulence->dproductioncoeffdu( aDofTypes ) -
                            tCMSATurbulence->dwalldestructioncoeffdu( aDofTypes ) );
                }

                // scale
                mdPPdMasterDof( tDofIndex ) = - std::pow( tTau, -1.5 ) * mdPPdMasterDof( tDofIndex );
            }
            else
            {
                mdPPdMasterDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------

        real SP_SUPG_Spalart_Allmaras_Turbulence::length_scale()
        {
            // if the length scale parameter was not evaluated
            if( mLengthScaleEval )
            {
                // evaluate the length scale parameter
                this->eval_length_scale();

                // set bool for evaluation
                mLengthScaleEval = false;
            }
            // return the length scale parameter value
            return mLengthScale;
        }

        void SP_SUPG_Spalart_Allmaras_Turbulence::eval_length_scale()
        {
            // get the velocity FI
            Field_Interpolator * tViscosityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

            // get the velocity FI
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // compute uTilde = u - cb2/sigma gradv
            Matrix< DDRMat > tModVelocity =
                    tVelocityFI->val() - mCb2 * tViscosityFI->gradx( 1 ) / mSigma;

            // compute and threshold the velocity norm (thresholding for consistency with derivatives)
            const real tNorm = std::max( norm( tModVelocity ), mEpsilon );

            // get the abs term
            const uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
            real tAbs = 0.0;
            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                tAbs += std::abs( dot( tModVelocity,tVelocityFI->dnNdxn( 1 ).get_column( iNode ) ) );
            }

            // threshold tAbs
            tAbs = std::max( tAbs, mEpsilon );

            // compute and threshold hugn
            mLengthScale =  std::max( 2.0 * tNorm / tAbs, mEpsilon );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & SP_SUPG_Spalart_Allmaras_Turbulence::dlengthscaledmasteru(
                const moris::Cell< MSI::Dof_Type > & aDofType )
        {
            // if aDofType is not an active dof type for the property
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType, mtk::Master_Slave::MASTER ),
                    "SP_SUPG_Spalart_Allmaras_Turbulence::dlengthscaledmasteru - no dependency on this dof type." );

            // get the dof index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mdLengthScaledMasterDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dlengthscaledmasteru( aDofType );

                // set bool for evaluation
                mdLengthScaledMasterDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdLengthScaledMasterDof( tDofIndex );
        }

        void SP_SUPG_Spalart_Allmaras_Turbulence::eval_dlengthscaledmasteru(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdLengthScaledMasterDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the length scale value
            real tHugn = this->length_scale();

            // compute derivative of hugn (compute only derivative if not thresholded)
            if ( tHugn > mEpsilon )
            {
                // get the velocity FI
                Field_Interpolator * tViscosityFI =
                        mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );

                // get the velocity FI
                Field_Interpolator * tVelocityFI =
                        mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

                // compute uTilde = u - cb2/sigma gradv
                Matrix< DDRMat > tModVelocity =
                        tVelocityFI->val() - mCb2 * tViscosityFI->gradx( 1 ) / mSigma;

                // compute and threshold the velocity norm (thresholding for consistency with derivatives)
                const real tNorm = std::max( norm( tModVelocity ), mEpsilon );

                // compute derivative of the modified velocity norm (compute only derivative if not thresholded)
                Matrix< DDRMat > tdNormdu( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
                if ( tNorm > mEpsilon )
                {
                    if( aDofTypes( 0 ) == mMasterDofVelocity )
                    {
                        tdNormdu += trans( tModVelocity ) * tVelocityFI->N() / tNorm;
                    }
                    else if( aDofTypes( 0 ) == mMasterDofViscosity )
                    {
                        tdNormdu -= mCb2 * trans( tModVelocity ) * tViscosityFI->dnNdxn( 1 ) / tNorm / mSigma;
                    }
                }

                // get the abs term
                const uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
                real tAbs = 0.0;
                for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                {
                    tAbs += std::abs( dot( tModVelocity,tVelocityFI->dnNdxn( 1 ).get_column( iNode ) ) );
                }

                // threshold tAbs
                tAbs = std::max( tAbs, mEpsilon );

                // compute derivative of the abs term (compute only derivative if not thresholded)
                Matrix< DDRMat > tdAbsdu( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
                if ( tAbs > mEpsilon )
                {
                    if( aDofTypes( 0 ) == mMasterDofVelocity )
                    {
                        uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
                        for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                        {
                            real tAdd = dot( tModVelocity,tVelocityFI->dnNdxn( 1 ).get_column( iNode ) );

                            // handle case that tAdd( 0, 0 ) is smaller than threshold
                            if ( std::abs( tAdd ) > mEpsilon )
                            {
                                tdAbsdu +=
                                        tAdd * trans( tVelocityFI->dnNdxn( 1 ).get_column( iNode ) ) *
                                        tVelocityFI->N() / std::abs( tAdd );
                            }
                        }
                    }
                    else if( aDofTypes( 0 ) == mMasterDofViscosity )
                    {
                        uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
                        for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                        {
                            real tAdd = dot( tModVelocity,tVelocityFI->dnNdxn( 1 ).get_column( iNode ) );

                            // handle case that tAdd is smaller than threshold
                            if ( std::abs( tAdd ) > mEpsilon )
                            {
                                tdAbsdu -=
                                        tAdd * trans( tVelocityFI->dnNdxn( 1 ).get_column( iNode ) ) *
                                        mCb2 * tViscosityFI->dnNdxn( 1 ) / std::abs( tAdd ) / mSigma;
                            }
                        }
                    }
                }

                // compute the derivative of the length scale
                mdLengthScaledMasterDof( tDofIndex ) =
                        2.0 * ( tdNormdu * tAbs - tdAbsdu * tNorm ) / std::pow( tAbs, 2.0 );
            }
            else
            {
                mdLengthScaledMasterDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


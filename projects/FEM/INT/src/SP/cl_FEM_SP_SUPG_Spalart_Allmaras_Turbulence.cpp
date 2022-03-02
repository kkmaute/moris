// FEM/INT/src
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
            mConstitutiveMap[ "SpalartAllmarasTurbulence" ] =
                    static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
        }

        //------------------------------------------------------------------------------

        void
        SP_SUPG_Spalart_Allmaras_Turbulence::set_parameters(
                moris::Cell< Matrix< DDRMat > > aParameters )
        {
            // FIXME not necessary
            // set mParameters
            mParameters = aParameters;

            // get number of parameters
            uint tParamSize = aParameters.size();

            // check for proper size of constant function parameters
            MORIS_ERROR( tParamSize <= 2,
                    "SP_SUPG_Spalart_Allmaras_Turbulence::set_parameters - max 2 constant parameters can to be set." );

            // if exponent specification
            if ( tParamSize > 0 )
            {
                // set exponent
                mExponent = aParameters( 0 )( 0 );
            }

            // if including reaction term specification
            if ( tParamSize > 1 )
            {
                // error unknown dof string
                MORIS_ERROR( aParameters( 1 )( 0 ) == 0.0 || aParameters( 1 )( 0 ) == 1.0,
                        "SP_SUPG_Spalart_Allmaras_Turbulence::set_parameters - invalid mHasReaction parameter \n" );

                // set exponent
                mHasReaction = aParameters( 1 )( 0 ) > 0 ? true : false;
            }
        }

        //------------------------------------------------------------------------------

        void
        SP_SUPG_Spalart_Allmaras_Turbulence::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > >& aDofTypes,
                moris::Cell< std::string >&                  aDofStrings,
                mtk::Master_Slave                            aIsMaster )
        {
            // switch on master slave
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER:
                {
                    // set dof type list
                    mMasterDofTypes = aDofTypes;

                    // loop on dof type
                    for ( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        std::string tDofString = aDofStrings( iDof );

                        // get dof type
                        MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if velocity
                        if ( tDofString == "Velocity" )
                        {
                            mMasterDofVelocity = tDofType;
                        }
                        else if ( tDofString == "Viscosity" )
                        {
                            mMasterDofViscosity = tDofType;
                        }
                        else
                        {
                            // error unknown dof string
                            MORIS_ERROR( false,
                                    "SP_SUPG_Spalart_Allmaras_Turbulence::set_dof_type_list - Unknown aDofString : %s \n",
                                    tDofString.c_str() );
                        }
                    }
                    break;
                }
                case mtk::Master_Slave::SLAVE:
                {
                    // set dof type list
                    mSlaveDofTypes = aDofTypes;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false,
                            "SP_SUPG_Spalart_Allmaras_Turbulence::set_dof_type_list - unknown master slave type." );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        SP_SUPG_Spalart_Allmaras_Turbulence::build_global_dof_type_list()
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

        void
        SP_SUPG_Spalart_Allmaras_Turbulence::set_function_pointers()
        {
        }

        //------------------------------------------------------------------------------
        /**
         * reset evaluation flags
         */
        void
        SP_SUPG_Spalart_Allmaras_Turbulence::reset_eval_flags()
        {
            // call parent implementation
            Stabilization_Parameter::reset_eval_flags();

            // reset child specific eval flags for chi
            mLengthScaleEval = true;
            mdLengthScaledMasterDofEval.fill( true );
        }

        //------------------------------------------------------------------------------

        void
        SP_SUPG_Spalart_Allmaras_Turbulence::eval_SP()
        {
            // set size for SP values
            mPPVal.set_size( 1, 1 );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // compute the length scale
            real tHugn = this->length_scale();

            // compute linearized version of uTilde = u - cb2/sigma gradv
            Matrix< DDRMat > tModVelocity = tCMSATurbulence->modified_velocity_linearized();

            // compute norm( uTilde ) and threshold tNormA (for consistency with derivative computation)
            real tNormA = std::max( norm( tModVelocity ), mEpsilon );

            // tau A
            real tTauA = 2.0 * tNormA / tHugn;

            // tau K
            real tTauK = 4.0 * tCMSATurbulence->diffusion_coefficient()( 0 ) / std::pow( tHugn, 2.0 );

            // evaluate tau
            real tTau = std::pow( tTauA, mExponent ) + std::pow( tTauK, mExponent );

            // if use reaction term
            if ( mHasReaction )
            {
                // compute tau S with linearized version of destruction term
                real tTauS =
                        2.0 * tCMSATurbulence->wall_destruction_coefficient()( 0 )
                        - tCMSATurbulence->production_coefficient()( 0 );

                // compute the regularized absolute of source term
                tTauS = std::sqrt( tTauS * tTauS + mEpsilon ) - std::sqrt( mEpsilon );

                // add contribution of reaction term
                tTau += std::pow( tTauS, mExponent );
            }

            // threshold tau
            tTau = std::max( tTau, mEpsilon );

            // set tau
            mPPVal = { { std::pow( tTau, -1.0 / mExponent ) } };
        }

        //------------------------------------------------------------------------------

        void
        SP_SUPG_Spalart_Allmaras_Turbulence::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator* tFI =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients() );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // compute the length scale
            real tHugn = this->length_scale();

            // evaluate linearized version of uTilde = u - 2/sigma cb2 gradv
            Matrix< DDRMat > tModVelocity = tCMSATurbulence->modified_velocity_linearized();

            // compute norm( uTilde ) and threshold tNormA (for consistency with derivative computation)
            real tNormA = std::max( norm( tModVelocity ), mEpsilon );

            // tau A
            real tTauA = 2.0 * tNormA / tHugn;

            // tau K
            real tTauK = 4.0 * tCMSATurbulence->diffusion_coefficient()( 0 ) / std::pow( tHugn, 2.0 );

            // evaluate tau
            real tTau = std::pow( tTauA, mExponent ) + std::pow( tTauK, mExponent );

            // if use reaction term
            real tTauS = 0;
            if ( mHasReaction )
            {
                // compute tau S with linearized version of destruction term
                tTauS =
                        2.0 * tCMSATurbulence->wall_destruction_coefficient()( 0 )
                        - tCMSATurbulence->production_coefficient()( 0 );

                // compute the regularized absolute of source term
                tTauS = std::sqrt( tTauS * tTauS + mEpsilon ) - std::sqrt( mEpsilon );

                // add contribution of reaction term
                tTau += std::pow( std::abs( tTauS ), mExponent );
            }

            // threshold tau
            tTau = std::max( tTau, mEpsilon );

            // compute dSPdu
            if ( tTau > mEpsilon )
            {
                // add contribution from the length scale derivative
                mdPPdMasterDof( tDofIndex ) = -this->dlengthscaledmasteru( aDofTypes )
                                            * ( std::pow( tTauA, mExponent - 1.0 ) * 2.0 * tNormA / std::pow( tHugn, 2.0 )
                                                    + std::pow( tTauK, mExponent - 1.0 ) * 8.0 * tCMSATurbulence->diffusion_coefficient()( 0 ) / std::pow( tHugn, 3.0 ) );

                // if turbulence CM depends on dof type
                if ( tCMSATurbulence->check_dof_dependency( aDofTypes ) )
                {
                    // if normA greater than threshold
                    if ( tNormA > mEpsilon )
                    {
                        // add contribution to mdPPdMasterDof from linearized version of modified velocity
                        mdPPdMasterDof( tDofIndex ) +=
                                std::pow( tTauA, mExponent - 1.0 ) * 2.0
                                * trans( tModVelocity ) * tCMSATurbulence->dmodvelocitylinearizeddu( aDofTypes ) / ( tHugn * tNormA );
                    }

                    // compute tdtauKdu
                    mdPPdMasterDof( tDofIndex ) +=
                            std::pow( tTauK, mExponent - 1.0 ) * 4.0 * tCMSATurbulence->ddiffusioncoeffdu( aDofTypes ) / std::pow( tHugn, 2.0 );

                    // if use reaction term
                    if ( mHasReaction )
                    {
                        // compute tdtauSdu
                        mdPPdMasterDof( tDofIndex ) +=
                                std::pow( tTauS, mExponent - 1.0 )
                                * tTauS / std::sqrt( tTauS * tTauS + mEpsilon )
                                * ( 2.0 * tCMSATurbulence->dwalldestructioncoeffdu( aDofTypes ) - tCMSATurbulence->dproductioncoeffdu( aDofTypes ) );
                    }
                }

                // scale
                mdPPdMasterDof( tDofIndex ) = -std::pow( tTau, -( mExponent + 1.0 ) / mExponent ) * mdPPdMasterDof( tDofIndex );
            }
            else
            {
                mdPPdMasterDof( tDofIndex ).fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------

        real
        SP_SUPG_Spalart_Allmaras_Turbulence::length_scale()
        {
            // if the length scale parameter was not evaluated
            if ( mLengthScaleEval )
            {
                // evaluate the length scale parameter
                this->eval_length_scale();

                // set bool for evaluation
                mLengthScaleEval = false;
            }
            // return the length scale parameter value
            return mLengthScale;
        }

        void
        SP_SUPG_Spalart_Allmaras_Turbulence::eval_length_scale()
        {
            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // compute and threshold the velocity norm (thresholding for consistency with derivatives)
            const real tNorm = std::max( norm( tCMSATurbulence->modified_velocity_linearized() ), mEpsilon );

            // get the velocity type FI
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get the abs term
            const uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
            real       tAbs      = 0.0;

            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                tAbs += std::abs( dot( tCMSATurbulence->modified_velocity_linearized(), tVelocityFI->dnNdxn( 1 ).get_column( iNode ) ) );
            }

            // threshold tAbs
            tAbs = std::max( tAbs, mEpsilon );

            // compute and threshold hugn
            mLengthScale = std::max( 2.0 * tNorm / tAbs, mEpsilon );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        SP_SUPG_Spalart_Allmaras_Turbulence::dlengthscaledmasteru(
                const moris::Cell< MSI::Dof_Type >& aDofType )
        {
            // if aDofType is not an active dof type for the property
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType, mtk::Master_Slave::MASTER ),
                    "SP_SUPG_Spalart_Allmaras_Turbulence::dlengthscaledmasteru - no dependency on this dof type." );

            // get the dof index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if ( mdLengthScaledMasterDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_dlengthscaledmasteru( aDofType );

                // set bool for evaluation
                mdLengthScaledMasterDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mdLengthScaledMasterDof( tDofIndex );
        }

        void
        SP_SUPG_Spalart_Allmaras_Turbulence::eval_dlengthscaledmasteru(
                const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator* tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdLengthScaledMasterDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the length scale value
            real tHugn = this->length_scale();

            // compute derivative of hugn (compute only derivative if not thresholded)
            if ( tHugn > mEpsilon )
            {
                // get the SA turbulence CM
                const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                        mMasterCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

                // get the velocity FI
                Field_Interpolator* tVelocityFI =
                        mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

                // compute and threshold the velocity norm (thresholding for consistency with derivatives)
                const real tNorm = std::max( norm( tCMSATurbulence->modified_velocity_linearized() ), mEpsilon );

                // compute derivative of the modified velocity norm (compute only derivative if not thresholded)
                Matrix< DDRMat > tdNormdu( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
                if ( tNorm > mEpsilon && tCMSATurbulence->check_dof_dependency( aDofTypes ) )
                {
                    tdNormdu += trans( tCMSATurbulence->modified_velocity_linearized() ) * tCMSATurbulence->dmodvelocitylinearizeddu( aDofTypes ) / tNorm;
                }

                // get the abs term
                const uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
                real       tAbs      = 0.0;
                for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                {
                    tAbs += std::abs( dot( tCMSATurbulence->modified_velocity_linearized(), tVelocityFI->dnNdxn( 1 ).get_column( iNode ) ) );
                }

                // threshold tAbs
                tAbs = std::max( tAbs, mEpsilon );

                // compute derivative of the abs term (compute only derivative if not thresholded)
                Matrix< DDRMat > tdAbsdu( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
                if ( tAbs > mEpsilon && tCMSATurbulence->check_dof_dependency( aDofTypes ) )
                {
                    uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
                    for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                    {
                        real tAdd = dot( tCMSATurbulence->modified_velocity_linearized(), tVelocityFI->dnNdxn( 1 ).get_column( iNode ) );

                        // handle case that tAdd is smaller than threshold
                        if ( std::abs( tAdd ) > mEpsilon )
                        {
                            tdAbsdu += tAdd * trans( tVelocityFI->dnNdxn( 1 ).get_column( iNode ) )
                                     * tCMSATurbulence->dmodvelocitylinearizeddu( aDofTypes ) / std::abs( tAdd );
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

//FEM/INT/src
#include "cl_FEM_SP_Turbulence_Dirichlet_Nitsche.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        SP_Turbulence_Dirichlet_Nitsche::SP_Turbulence_Dirichlet_Nitsche()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Viscosity" ] = static_cast< uint >( SP_Property_Type::VISCOSITY );
        }

        //------------------------------------------------------------------------------

        void SP_Turbulence_Dirichlet_Nitsche::set_dof_type_list(
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
                        if( tDofString == "Viscosity" )
                        {
                            mMasterDofViscosity = tDofType;
                        }
                        else
                        {
                            // create error message
                            std::string tErrMsg =
                                    std::string( "SP_Turbulence_Dirichlet_Nitsche::set_dof_type_list - Unknown aDofString : ") +
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
                    MORIS_ERROR( false, "SP_Turbulence_Dirichlet_Nitsche::set_dof_type_list - unknown master slave type." );
            }
        }

        //------------------------------------------------------------------------------

        moris::Cell< std::tuple<
        fem::Measure_Type,
        mtk::Primary_Void,
        mtk::Master_Slave > > SP_Turbulence_Dirichlet_Nitsche::get_cluster_measure_tuple_list()
        {
            return { mElementSizeTuple };
        }

        //------------------------------------------------------------------------------

        void SP_Turbulence_Dirichlet_Nitsche::eval_SP()
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );

            // compute the diffusion coefficient
            real tDiffusionCoeff = compute_diffusion_coefficient(
                    { mMasterDofViscosity },
                    mMasterFIManager,
                    tPropKinViscosity );

            // compute stabilization parameter value
            mPPVal = mParameters( 0 ) * tDiffusionCoeff / tElementSize;
        }

        //------------------------------------------------------------------------------

        void SP_Turbulence_Dirichlet_Nitsche::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get element size cluster measure value
            real tElementSize = mCluster->get_cluster_measure(
                    std::get<0>( mElementSizeTuple ),
                    std::get<1>( mElementSizeTuple ),
                    std::get<2>( mElementSizeTuple ) )->val()( 0 );

            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( tDofType );

            // get FI for derivative dof type
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdMasterDof( tDofIndex ).set_size(
                    1,
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the kinematic viscosity property
            const std::shared_ptr< Property > & tPropKinViscosity =
                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );

            // compute the derivative of the diffusion coefficient
            Matrix< DDRMat > tdDiffusiondu;
            compute_ddiffusiondu(
                    { mMasterDofViscosity },
                    mMasterFIManager,
                    tPropKinViscosity,
                    aDofTypes,
                    tdDiffusiondu );

            // add contribution from diffusion coefficient
            mdPPdMasterDof( tDofIndex ) +=
                    mParameters( 0 ) * tdDiffusiondu / tElementSize;
        }

//        //------------------------------------------------------------------------------
//
//        real SP_Turbulence_Dirichlet_Nitsche::compute_diffusion_coefficient()
//        {
//            // init diffusion coeff
//            real tDiffusionTerm = 0.0;
//
//            // get the viscosity FI
//            Field_Interpolator * tFIModViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//            // get the wall distance property
//            const std::shared_ptr< Property > & tPropKinViscosity =
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
//        void SP_Turbulence_Dirichlet_Nitsche::compute_ddiffusiondu(
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
//            const std::shared_ptr< Property > & tPropKinViscosity =
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
//        real SP_Turbulence_Dirichlet_Nitsche::compute_fn()
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
//        void SP_Turbulence_Dirichlet_Nitsche::compute_dfndu(
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
//        real SP_Turbulence_Dirichlet_Nitsche::compute_chi()
//        {
//            // get the residual dof FI (here viscosity)
//            Field_Interpolator * tFIViscosity =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofViscosity );
//
//            // get the density and gravity properties
//            const std::shared_ptr< Property > & tPropViscosity =
//                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );
//
//            // compute chi
//            return tFIViscosity->val()( 0 ) / tPropViscosity->val()( 0 );
//        }
//
//        //------------------------------------------------------------------------------
//
//        void SP_Turbulence_Dirichlet_Nitsche::compute_dchidu(
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
//            const std::shared_ptr< Property > & tPropViscosity =
//                    mMasterProp( static_cast< uint >( SP_Property_Type::VISCOSITY ) );
//
//            // compute chi
//            real tChi = tFIViscosity->val()( 0 ) / tPropViscosity->val()( 0 );
//
//            // if residual dof type (here viscosity)
//            if( aDofTypes( 0 ) == mMasterDofViscosity )
//            {
//                adchidu += tDerFI->N() / tPropViscosity->val()( 0 );
//            }
//
//            if( tPropViscosity->check_dof_dependency( aDofTypes ) )
//            {
//                adchidu -= tChi * tPropViscosity->dPropdDOF( aDofTypes ) / tPropViscosity->val()( 0 );
//            }
//        }


        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

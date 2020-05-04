
#include "cl_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        IWG_Spalart_Allmaras_Turbulence_Bulk::IWG_Spalart_Allmaras_Turbulence_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "WallDistance" ] = IWG_Property_Type::WALL_DISTANCE;
            mPropertyMap[ "Viscosity" ]    = IWG_Property_Type::VISCOSITY;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
//            mConstitutiveMap[ "IncompressibleFluid" ] = IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "SUPG" ] = IWG_Stabilization_Type::SUPG;

        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity
            = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the velocity FI
            // FISME protect dof type
            Field_Interpolator * tFIVelocity
            = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropWallDistance
            = mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );

            // compute ft2
            real tFt2 = this->compute_ft2();

            // compute STilde
            //real tSTilde = this->compute_stilde();

            // compute the residual weak form
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
            += aWStar * (   trans( tFIViscosity->N() ) * tFIViscosity->gradt( 1 )
                          + trans( tFIViscosity->N() ) * trans( tFIVelocity->val() ) * tFIViscosity->gradx( 1 )
                          - trans( tFIViscosity->N() ) * mCb1 * ( 1 - tFt2 ) * tFIViscosity->val() );
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity
            = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the velocity FI
            // FISME protect dof type
            Field_Interpolator * tFIVelocity
            = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropWallDistance
            = mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );

            // compute ft2
            real tFt2 = this->compute_ft2();

            // compute STilde
            //real tSTilde = this->compute_stilde();

            // get number of dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over the dof dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if residual dof type (here viscosity)
                if( tDofType( 0 ) == mResidualDofType( 0 ) )
                {
                    // compute the jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += aWStar * (   trans( tFIViscosity->N() ) * tFIViscosity->dnNdtn( 1 )
                                  + trans( tFIViscosity->N() ) * trans( tFIVelocity->val() ) * tFIViscosity->dnNdxn( 1 )
                                  - trans( tFIViscosity->N() ) * mCb1 * ( 1 - tFt2 ) * tFIViscosity->N() );
                }

                // if velocity dof type
                // FIXME protect dof type
                if( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    // compute the jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += aWStar * ( trans( tFIViscosity->N() ) * trans( tFIViscosity->gradx( 1 ) ) * tFIVelocity->N() );
                }

                // compute dft2du
                Matrix< DDRMat > tdft2du;
                this->compute_dft2du( tDofType, tdft2du );

//                // compute dSTildedu
//                Matrix< DDRMat > tdstildedu;
//                this->compute_dstildedu( tDofType, tdstildedu );

                // compute the jacobian
                mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                      { tMasterDepStartIndex, tMasterDepStopIndex } )
                += aWStar * (   trans( tFIViscosity->N() ) * mCb1 * tFIViscosity->val() * tdft2du );
                              //- trans( tFIViscosity->N() ) * mCb1 * ( 1 - tFt2 )* tFIViscosity->val() );
            }

        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_dRdp - Not implemented." );
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual_strong_form( Matrix< DDRMat > & aR )
        {
            MORIS_ASSERT( false, "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual_strong_form - not implemented.");
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian_strong_form
        ( moris::Cell< MSI::Dof_Type >   aDofTypes,
          Matrix< DDRMat >             & aJ )
        {
            MORIS_ASSERT( false, "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian_strong_form - not implemented.");
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_wij( Matrix< DDRMat > & aWij )
        {
            // get the velocity FI
            // FIXME protect dof
            Field_Interpolator * tFIVelocity = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get gradient of velocity
            Matrix< DDRMat > tGradVelocity = tFIVelocity->gradx( 1 );

            // switch on space dim
            switch ( tFIVelocity->get_number_of_fields() )
            {
                case ( 2 ):
                {
                    // init aWij = [ w11 w12 w21 w22]
                    aWij.set_size( 4, 1, 0.0 );

                    // compute Wij
                    aWij( 1 ) = 0.5 * ( tGradVelocity( 1, 0 ) - tGradVelocity( 0, 1 ) );
                    aWij( 2 ) = 0.5 * ( tGradVelocity( 0, 1 ) - tGradVelocity( 1, 0 ) );
                    break;
                }
                case ( 3 ):
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
                }
                default:
                    MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_wij - space dim can only be 2 or 3" );
                    break;
            }
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dwijdu( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                                   Matrix< DDRMat >             & adwijdu )
        {
            // get the der FI
            Field_Interpolator * tFIDer = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
            Field_Interpolator * tFIVelocity = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // switch on space dim
            switch ( tFIVelocity->get_number_of_fields() )
            {
                case ( 2 ):
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
                case ( 3 ):
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
                        adwijdu( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } )  =   0.5 * tdNdxVelocity.get_row( 0 );

                        adwijdu( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } )     =   0.5 * tdNdxVelocity.get_row( 2 );
                        adwijdu( { 5, 5 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 1 );

                        adwijdu( { 6, 6 }, { 0, tNumBases - 1 } )                 = - 0.5 * tdNdxVelocity.get_row( 2 );
                        adwijdu( { 6, 6 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =   0.5 * tdNdxVelocity.get_row( 0 );

                        adwijdu( { 7, 7 }, { tNumBases, 2 * tNumBases - 1 } )     =   0.5 * tdNdxVelocity.get_row( 1 );
                        adwijdu( { 7, 7 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 2 );
                    }
                }
                default:
                    MORIS_ERROR( false, "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dwijdu - space dim can only be 2 or 3" );
                    break;
            }
        }

//------------------------------------------------------------------------------
        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_sbar()
        {
            // compute wij
            Matrix< DDRMat > tWij;
            this->compute_wij( tWij );

            // compute WijWij
            Matrix< DDRMat > tWijWij = trans( tWij ) * tWij;

            // compute sbar
            real tSBar = std::sqrt( 2.0 * tWijWij( 0 ) );

            return tSBar;
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dsbardu( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                                    Matrix< DDRMat >             & adsbardu )
        {
            // compute sbar
            real tSBar = this->compute_sbar();

            // compute wij
            Matrix< DDRMat > tWij;
            this->compute_wij( tWij );

            // compute dwijdu
            Matrix< DDRMat > tdWijdu;
            this->compute_dwijdu( aDofTypes, tdWijdu );

            // compute dsbardu
            adsbardu = 2.0 * trans( tWij ) * tdWijdu / tSBar;
        }

//------------------------------------------------------------------------------
        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_chi()
        {
            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity
            = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropViscosity
            = mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // compute chi
            real tChi = tFIViscosity->val()( 0 ) / tPropViscosity->val()( 0 );

            return tChi;
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dchidu( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                                   Matrix< DDRMat >             & adchidu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adchidu
            adchidu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity
            = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropViscosity
            = mMasterProp( static_cast< uint >( IWG_Property_Type::VISCOSITY ) );

            // compute chi
            real tChi = tFIViscosity->val()( 0 ) / tPropViscosity->val()( 0 );

            // if residual dof type (here viscosity)
            if( aDofTypes( 0 ) == mResidualDofType( 0 ) )
            {
                adchidu.matrix_data() += tDerFI->N() / tPropViscosity->val()( 0 );
            }

            if( tPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                adchidu.matrix_data() += tChi * tPropViscosity->dPropdDOF( aDofTypes ) / tPropViscosity->val()( 0 );
            }
        }

//------------------------------------------------------------------------------
        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_fv1()
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute fv1
            real tFv1 = std::pow( tChi, 3.0 ) / ( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ) );

            return tFv1;
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dfv1du( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                                   Matrix< DDRMat >             & adfv1du )
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            this->compute_dchidu( aDofTypes, tdchidu );

            // compute adfv1du
            adfv1du = 3.0 * std::pow( mCv1, 3.0 ) * std::pow( tChi, 2.0 ) * tdchidu / std::pow( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ), 2.0 );
        }

//------------------------------------------------------------------------------
        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_fv2()
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute fv1
            real tFv1 = this->compute_fv1();

            // compute fv2
            real tFv2 = 1.0 - tChi / ( 1 + tChi * tFv1 );

            return tFv2;
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dfv2du( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                                   Matrix< DDRMat >             & adfv2du )
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
            adfv2du = -1.0 - tChi * tFv1 + tChi * tFv1 * tdchidu + std::pow( tChi, 2.0 ) * tdfv1du / ( std::pow( 1.0 + tChi * tFv1, 2.0 ) );
        }

//------------------------------------------------------------------------------
        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_stilde()
        {
            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity
            = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the wall distance property
            std::shared_ptr< Property > tPropWallDistance
            = mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );

            // compute SBar
            real tSBar = this->compute_sbar();

            // compute fv2
            real tFv2 = this->compute_fv2();

            // compute stilde
            real tStilde = tSBar + tFv2 * tFIViscosity->val()( 0 ) / ( std::pow( mKappa, 2.0 ) * std::pow( tPropWallDistance->val()( 0 ), 2.0 ) );

            return tStilde;
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dstildedu( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                                      Matrix< DDRMat >             & adstildedu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adSTildedu
            adstildedu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity
            = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the wall distance properties
            std::shared_ptr< Property > tPropWallDistance
            = mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );

            // compute dSBardu
            Matrix< DDRMat > tdSBardu;
            this->compute_dsbardu( aDofTypes, tdSBardu );

            // compute fv2
            real tFv2 = this->compute_fv2();

            // compute dfv2du
            Matrix< DDRMat > tdfv2du;
            this->compute_dfv2du( aDofTypes, tdfv2du );

            // compute dstildedu
            adstildedu.matrix_data() += tdSBardu + tFIViscosity->val() * tdfv2du / ( std::pow( mKappa, 2.0 ) * std::pow( tPropWallDistance->val()( 0 ), 2.0 ) );

            // if dof type is residual dof type
            if( aDofTypes( 0 ) == mResidualDofType( 0 ) )
            {
                // add contribution
                adstildedu.matrix_data() += tFv2 * tDerFI->N() / ( std::pow( mKappa, 2.0 ) * std::pow( tPropWallDistance->val()( 0 ), 2.0 ) );
            }
        }

//------------------------------------------------------------------------------
        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_ft2()
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute ft2
            real tFt2 = mCt3 * std::exp( - mCt4 * std::pow( tChi, 2.0 ) );

            return tFt2;
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dft2du( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                                   Matrix< DDRMat >             & adft2du )
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
        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_r()
        {
            // compute stilde
            real tSTilde = this->compute_stilde();

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity
            = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the wall distance property
            std::shared_ptr< Property > tPropWallDistance
            = mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );

            // compute viscosity / ( stilde * kappa² * d² )
            real tR = tFIViscosity->val()( 0 ) / ( tSTilde * std::pow( mKappa, 2.0 ) * std::pow( tPropWallDistance->val()( 0 ), 2.0 ) );

            // compute r
            if ( tR > 10.0 )
            {
                tR = 10.0;
            }

            return tR;
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_drdu( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                                  Matrix< DDRMat >             & adrdu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adrdu
            adrdu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // compute r
            real tR = this->compute_r();

            // if r > 10
            if( tR < 10 )
            {
                // get the wall distance property
                std::shared_ptr< Property > tPropWallDistance
                = mMasterProp( static_cast< uint >( IWG_Property_Type::WALL_DISTANCE ) );

                // compute stilde
                real tSTilde = this->compute_stilde();

                // compute dStildedu
                Matrix< DDRMat > tdSTildedu;
                this->compute_dstildedu( aDofTypes, tdSTildedu );

                // add contribution from dStildedu
                adrdu.matrix_data() -= tDerFI->val() * std::pow( mKappa, 2.0 ) * std::pow( tPropWallDistance->val()( 0 ), 2.0 ) * tdSTildedu;

                // if dof type is viscosity
                if( aDofTypes( 0 ) == mResidualDofType( 0 ) )
                {
                    // add contribution from viscosity
                    adrdu.matrix_data() += tDerFI->N().matrix_data();
                }

                adrdu = adrdu / ( std::pow( tSTilde * std::pow( mKappa, 2.0 ) * std::pow( tPropWallDistance->val()( 0 ), 2.0 ), 2.0 ) );
            }
        }

//------------------------------------------------------------------------------
        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_g()
        {
            // compute r
            real tR = this->compute_r();

            // compute g
            real tG = tR + mCw2 * ( std::pow( tR, 6.0 ) - tR );

            return tG;
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dgdu( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                                 Matrix< DDRMat >             & adgdu )
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
        real IWG_Spalart_Allmaras_Turbulence_Bulk::compute_fw()
        {
            // compute g
            real tG = this->compute_g();

            // compute fw
            real tFw = ( 1.0 + std::pow( mCw3, 6.0 ) ) / ( std::pow( tG, 6.0 ) + std::pow( mCw3, 6.0 ) );
            tFw = tG * std::pow( tFw, 1/6 );

            return tFw;
        }

//------------------------------------------------------------------------------
        void IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dfwdu( moris::Cell< MSI::Dof_Type >   aDofTypes,
                                                                  Matrix< DDRMat >             & adfwdu )
        {
            // compute g
            real tG = this->compute_g();

            // compute dgdu
            Matrix< DDRMat > tdgdu;
            this->compute_dgdu( aDofTypes, tdgdu );

            // compute fw
            real tFw = this->compute_fw();

            // init adfwdu
            adfwdu = ( 1.0 - 5.0 * std::pow( tG, 6.0 ) / ( 6.0 * ( std::pow( tG, 6.0 ) + std::pow( mCw3, 6.0 ) ) ) ) * ( tFw / tG ) * tdgdu;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

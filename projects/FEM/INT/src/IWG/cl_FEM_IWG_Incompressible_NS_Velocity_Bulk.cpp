
#include "cl_FEM_IWG_Incompressible_NS_Velocity_Bulk.hpp"
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

        IWG_Incompressible_NS_Velocity_Bulk::IWG_Incompressible_NS_Velocity_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Gravity" ]          = static_cast< uint >( IWG_Property_Type::GRAVITY );
            mPropertyMap[ "ThermalExpansion" ] = static_cast< uint >( IWG_Property_Type::THERMAL_EXPANSION );
            mPropertyMap[ "ReferenceTemp" ]    = static_cast< uint >( IWG_Property_Type::REF_TEMP );
            mPropertyMap[ "InvPermeability" ]  = static_cast< uint >( IWG_Property_Type::INV_PERMEABILITY );
            mPropertyMap[ "MassSource" ]       = static_cast< uint >( IWG_Property_Type::MASS_SOURCE );
            mPropertyMap[ "Load" ]             = static_cast< uint >( IWG_Property_Type::BODY_LOAD );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "IncompressibleFluid" ] =
                    static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "IncompressibleFlow" ] =
                    static_cast< uint >( IWG_Stabilization_Type::INCOMPRESSIBLE_FLOW );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Velocity_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the velocity FI
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the gravity property
            const std::shared_ptr< Property >& tGravityProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::GRAVITY ) );

            // get the thermal expansion property
            const std::shared_ptr< Property >& tThermalExpProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THERMAL_EXPANSION ) );

            // get the reference temperature property
            const std::shared_ptr< Property >& tRefTempProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::REF_TEMP ) );

            // get the inverted permeability (porosity) property
            const std::shared_ptr< Property >& tInvPermeabProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::INV_PERMEABILITY ) );

            // get the body load property
            const std::shared_ptr< Property >& tLoadProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            // get the incompressible fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tIncFluidCM =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the density property
            const std::shared_ptr< Property >& tDensityProp = tIncFluidCM->get_property( "Density" );

            // get the incompressible flow stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPSUPSPSPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::INCOMPRESSIBLE_FLOW ) );

            // compute the residual strong form
            Matrix< DDRMat > tRM;
            real             tRC;
            this->compute_residual_strong_form( tRM, tRC );

            // get the density value
            const real tDensity = tDensityProp->val()( 0 );

            // compute uj vij
            Matrix< DDRMat > tujvij;
            this->compute_ujvij( tujvij );

            // build multiplication matrix for sigma_ij epsilon_ij
            Matrix< DDRMat > tPre;
            if ( tVelocityFI->get_number_of_fields() == 2 )
            {
                tPre.set_size( 3, 3, 0.0 );
                tPre( 0, 0 ) = 1.0;
                tPre( 1, 1 ) = 1.0;
                tPre( 2, 2 ) = 2.0;
            }
            else
            {
                tPre.set_size( 6, 6, 0.0 );
                tPre( 0, 0 ) = 1.0;
                tPre( 1, 1 ) = 1.0;
                tPre( 2, 2 ) = 1.0;
                tPre( 3, 3 ) = 2.0;
                tPre( 4, 4 ) = 2.0;
                tPre( 5, 5 ) = 2.0;
            }

            // get sub-matrix of residual
            auto tRes = mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex } );

            // compute the residual
            tRes += aWStar
                  * ( tDensity * tVelocityFI->N_trans()
                                  * ( trans( tVelocityFI->gradt( 1 ) ) + trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->val() )
                          + trans( tIncFluidCM->testStrain() ) * tPre * tIncFluidCM->flux()
                          + trans( tujvij ) * tDensity * tSPSUPSPSPG->val()( 0 ) * tRM
                          + trans( tVelocityFI->div_operator() ) * tSPSUPSPSPG->val()( 1 ) * tRC );

            // if gravity
            if ( tGravityProp != nullptr )
            {
                // add gravity to residual weak form
                tRes += aWStar * ( tVelocityFI->N_trans() * tDensity * tGravityProp->val() );
            }

            // if thermal expansion and reference temperature
            if ( tGravityProp != nullptr && tThermalExpProp != nullptr && tRefTempProp != nullptr )
            {
                // get the temperature field interpolator
                // FIXME protect FI
                Field_Interpolator* tTempFI = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

                // add contribution to residual
                tRes -= aWStar
                      * ( tVelocityFI->N_trans() * tDensity * tGravityProp->val() * tThermalExpProp->val()
                              * ( tTempFI->val() - tRefTempProp->val() ) );
            }

            // if permeability
            if ( tInvPermeabProp != nullptr )
            {
                // add Brinkman term to residual weak form
                tRes += aWStar * ( tVelocityFI->N_trans() * tInvPermeabProp->val()( 0 ) * tVelocityFI->val() );
            }

            // if body load
            if ( tLoadProp != nullptr ) { tRes -= aWStar * ( tVelocityFI->N_trans() * tLoadProp->val() ); }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Incompressible_NS_Velocity_Bulk::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Velocity_Bulk::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            const uint tMasterDofIndex =
                    mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );

            const uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            const uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get velocity FI
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the gravity property
            const std::shared_ptr< Property >& tGravityProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::GRAVITY ) );

            // get the thermal expansion property
            const std::shared_ptr< Property >& tThermalExpProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THERMAL_EXPANSION ) );

            // get the reference temperature property
            const std::shared_ptr< Property >& tRefTempProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::REF_TEMP ) );

            // get the inverted permeability property
            const std::shared_ptr< Property >& tInvPermeabProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::INV_PERMEABILITY ) );

            // get the body load property
            const std::shared_ptr< Property >& tLoadProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            // get the incompressible fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tIncFluidCM =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the density property from CM
            const std::shared_ptr< Property >& tDensityProp = tIncFluidCM->get_property( "Density" );

            // evaluate the density
            const real tDensity = tDensityProp->val()( 0 );

            // get the incompressible flow stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPSUPSPSPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::INCOMPRESSIBLE_FLOW ) );

            // build multiplication matrix for sigma_ij epsilon_ij
            Matrix< DDRMat > tPre;
            if ( tVelocityFI->get_number_of_fields() == 2 )
            {
                tPre.set_size( 3, 3, 0.0 );
                tPre( 0, 0 ) = 1.0;
                tPre( 1, 1 ) = 1.0;
                tPre( 2, 2 ) = 2.0;
            }
            else
            {
                tPre.set_size( 6, 6, 0.0 );
                tPre( 0, 0 ) = 1.0;
                tPre( 1, 1 ) = 1.0;
                tPre( 2, 2 ) = 1.0;
                tPre( 3, 3 ) = 2.0;
                tPre( 4, 4 ) = 2.0;
                tPre( 5, 5 ) = 2.0;
            }

            // compute the residual strong form
            Matrix< DDRMat > tRM;
            real             tRC;
            this->compute_residual_strong_form( tRM, tRC );

            // compute the Jacobian for dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Cell< MSI::Dof_Type >& tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                const sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                const uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                const uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // compute uj vij
                Matrix< DDRMat > tujvij;
                this->compute_ujvij( tujvij );

                // extract sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex }, { tMasterDepStartIndex, tMasterDepStopIndex } );

                // if residual dof type (velocity)
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute dnNdtn
                    Matrix< DDRMat > tdnNdtn;
                    this->compute_dnNdtn( tdnNdtn );

                    // compute uj vij rM
                    Matrix< DDRMat > tujvijrM;
                    this->compute_ujvijrm( tujvijrM, tRM );

                    // compute the Jacobian
                    tJac += aWStar
                          * ( tDensity * tVelocityFI->N_trans()
                                          * ( tdnNdtn + trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->N() + tujvij )
                                  + tDensity * tujvijrM * tSPSUPSPSPG->val()( 0 ) );

                    // if permeability
                    if ( tInvPermeabProp != nullptr )
                    {
                        // add Brinkman term to Jacobian of weak form
                        tJac += aWStar * ( tVelocityFI->N_trans() * tInvPermeabProp->val()( 0 ) * tVelocityFI->N() );
                    }
                }

                // compute the Jacobian strong form
                Matrix< DDRMat > tJM;
                Matrix< DDRMat > tJC;
                compute_jacobian_strong_form( tDofType, tJM, tJC );

                // compute the Jacobian contribution from strong form
                tJac += aWStar
                      * ( trans( tujvij ) * tDensity * tSPSUPSPSPG->val()( 0 ) * tJM
                              + trans( tVelocityFI->div_operator() ) * tSPSUPSPSPG->val()( 1 ) * tJC );

                // if property has dependency on the dof type
                if ( tDensityProp->check_dof_dependency( tDofType ) )
                {
                    // compute the Jacobian
                    tJac += aWStar
                          * ( tVelocityFI->N_trans() * trans( tVelocityFI->gradt( 1 ) )
                                  + tVelocityFI->N_trans() * trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->val()
                                  + trans( tujvij ) * tSPSUPSPSPG->val()( 0 ) * tRM )
                          * tDensityProp->dPropdDOF( tDofType );
                }

                // if fluid CM depends on dof type
                if ( tIncFluidCM->check_dof_dependency( tDofType ) )
                {
                    // compute the Jacobian
                    tJac += aWStar * ( trans( tIncFluidCM->testStrain() ) * tPre * tIncFluidCM->dFluxdDOF( tDofType ) );
                    // FIXME add dteststrainddof
                }

                // if stabilization parameter has dependency on the dof type
                if ( tSPSUPSPSPG->check_dof_dependency( tDofType ) )
                {
                    // compute the Jacobian
                    tJac += aWStar
                          * ( trans( tujvij ) * tDensity * tRM * tSPSUPSPSPG->dSPdMasterDOF( tDofType ).get_row( 0 )
                                  + trans( tVelocityFI->div_operator() ) * tRC
                                            * tSPSUPSPSPG->dSPdMasterDOF( tDofType ).get_row( 1 ) );
                }

                // if permeability depends on DoF type
                if ( tInvPermeabProp != nullptr )
                {
                    if ( tInvPermeabProp->check_dof_dependency( tDofType ) )
                    {
                        tJac += aWStar
                              * ( tVelocityFI->N_trans() * tVelocityFI->val()
                                      * tInvPermeabProp->dPropdDOF( tDofType ) );
                    }
                }

                // if body load term depends on DoF type
                if ( tLoadProp != nullptr )
                {
                    if ( tLoadProp->check_dof_dependency( tDofType ) )
                    {
                        tJac -= aWStar * ( tVelocityFI->N_trans() * tLoadProp->dPropdDOF( tDofType ) );
                    }
                }

                // if gravity
                if ( tGravityProp != nullptr )
                {
                    // if density property depends on dof type
                    if ( tDensityProp->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJac += aWStar
                              * ( tVelocityFI->N_trans() * tGravityProp->val() * tDensityProp->dPropdDOF( tDofType ) );
                    }

                    // if gravity property depends on dof type
                    if ( tGravityProp->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJac += aWStar * ( tVelocityFI->N_trans() * tDensity * tGravityProp->dPropdDOF( tDofType ) );
                    }
                }

                // if thermal expansion and reference temperature
                if ( tGravityProp != nullptr && tThermalExpProp != nullptr && tRefTempProp != nullptr )
                {
                    // get the temperature field interpolator
                    // FIXME protect FI
                    Field_Interpolator* tTempFI =
                            mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

                    // if dof is temperature
                    if ( tDofType( 0 ) == MSI::Dof_Type::TEMP )
                    {
                        // add contribution to Jacobian
                        tJac -= aWStar
                              * ( tVelocityFI->N_trans() * tDensity * tGravityProp->val() * tThermalExpProp->val()
                                      * tTempFI->N() );
                    }

                    // if density property depends on dof type
                    if ( tDensityProp->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJac -= aWStar
                              * ( tVelocityFI->N_trans() * tGravityProp->val() * tThermalExpProp->val()
                                      * ( tTempFI->val() - tRefTempProp->val() )
                                      * tDensityProp->dPropdDOF( tDofType ) );
                    }

                    // if thermal expansion property depends on dof type
                    if ( tThermalExpProp->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJac -= aWStar
                              * ( tVelocityFI->N_trans() * tDensity * tGravityProp->val()
                                      * ( tTempFI->val() - tRefTempProp->val() )
                                      * tThermalExpProp->dPropdDOF( tDofType ) );
                    }

                    // if reference temperature property depends on dof type
                    if ( tRefTempProp->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJac += aWStar
                              * ( tVelocityFI->N_trans() * tDensity * tGravityProp->val() * tThermalExpProp->val()
                                      * tRefTempProp->dPropdDOF( tDofType ) );
                    }

                    // if gravity property depends on dof type
                    if ( tGravityProp->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJac -= aWStar
                              * ( tVelocityFI->N_trans() * tDensity * tThermalExpProp->val()( 0 )
                                      * ( tTempFI->val()( 0 ) - tRefTempProp->val()( 0 ) )
                                      * tGravityProp->dPropdDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Incompressible_NS_Velocity_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Velocity_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR(
                    false, "IWG_Incompressible_NS_Velocity_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Velocity_Bulk::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Velocity_Bulk::compute_residual_strong_form( Matrix< DDRMat >& aRM, real& aRC )
        {
            // get the velocity and pressure FIs
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the properties
            const std::shared_ptr< Property >& tGravityProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::GRAVITY ) );

            const std::shared_ptr< Property >& tThermalExpProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THERMAL_EXPANSION ) );

            const std::shared_ptr< Property >& tRefTempProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::REF_TEMP ) );

            const std::shared_ptr< Property >& tInvPermeabProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::INV_PERMEABILITY ) );

            const std::shared_ptr< Property >& tMassSourceProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::MASS_SOURCE ) );

            // get the body load property
            const std::shared_ptr< Property >& tLoadProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            // get the incompressible fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tIncFluidCM =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the density property from CM
            const std::shared_ptr< Property >& tDensityProp = tIncFluidCM->get_property( "Density" );

            // get the density value
            real tDensity = tDensityProp->val()( 0 );

            // compute the residual strong form of momentum equation
            aRM = tDensity * trans( tVelocityFI->gradt( 1 ) )
                + tDensity * trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->val()
                - tIncFluidCM->divflux();

            // compute the residual strong form of continuity equation
            aRC = tVelocityFI->div();

            // if body load
            if ( tLoadProp != nullptr )
            {
                // add contribution of body load term to momentum residual
                aRM -= tLoadProp->val();
            }

            // if mass source
            if ( tMassSourceProp != nullptr )
            {
                // add mass source to continuity residual
                aRC -= tMassSourceProp->val()( 0 ) / tDensity;
            }

            // if permeability
            if ( tInvPermeabProp != nullptr )
            {
                // add Brinkman term to residual strong form
                aRM += tInvPermeabProp->val()( 0 ) * tVelocityFI->val();
            }

            // if gravity
            if ( tGravityProp != nullptr )
            {
                // add gravity to residual strong form
                aRM += tDensity * tGravityProp->val();

                // if thermal expansion and reference temperature
                if ( tThermalExpProp != nullptr && tRefTempProp != nullptr )
                {
                    // get the temperature field interpolator
                    // FIXME protect FI
                    Field_Interpolator* tTempFI =
                            mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

                    // add contribution to residual
                    aRM -= tDensity * tGravityProp->val() * tThermalExpProp->val()
                         * ( tTempFI->val() - tRefTempProp->val() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Velocity_Bulk::compute_jacobian_strong_form(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                Matrix< DDRMat >&                   aJM,
                Matrix< DDRMat >&                   aJC )
        {
            // get the res dof and the derivative dof FIs
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );
            Field_Interpolator* tDerFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init aJM and aJC
            aJM.set_size( tVelocityFI->get_number_of_fields(), tDerFI->get_number_of_space_time_coefficients() );

            aJC.set_size( 1, tDerFI->get_number_of_space_time_coefficients() );

            // get the properties
            const std::shared_ptr< Property >& tGravityProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::GRAVITY ) );

            const std::shared_ptr< Property >& tThermalExpProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THERMAL_EXPANSION ) );

            const std::shared_ptr< Property >& tRefTempProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::REF_TEMP ) );

            const std::shared_ptr< Property >& tInvPermeabProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::INV_PERMEABILITY ) );

            const std::shared_ptr< Property >& tMassSourceProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::MASS_SOURCE ) );

            // get the body load property
            const std::shared_ptr< Property >& tLoadProp =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            // get the incompressible fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tIncFluidCM =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the density property from CM
            const std::shared_ptr< Property >& tDensityProp = tIncFluidCM->get_property( "Density" );

            // get the density value
            real tDensity = tDensityProp->val()( 0 );

            // if derivative wrt to residual dof type (here velocity)
            if ( aDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) )
            {
                // compute the term uj vij
                Matrix< DDRMat > tujvij;
                this->compute_ujvij( tujvij );

                // compute the term dnNdtn
                Matrix< DDRMat > tdnNdtn;
                this->compute_dnNdtn( tdnNdtn );

                // compute the Jacobian strong form of momentum equation
                aJM = tDensity * tdnNdtn + tDensity * trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->N()
                    + tDensity * tujvij;

                // compute the Jacobian strong form of continuity equation
                aJC = tVelocityFI->div_operator();
            }
            else
            {
                aJM.fill( 0.0 );
                aJC.fill( 0.0 );
            }

            // if density depends on dof type
            if ( tDensityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution to Jacobian strong form
                aJM += trans( tVelocityFI->gradt( 1 ) ) * tDensityProp->dPropdDOF( aDofTypes )
                     + trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->val() * tDensityProp->dPropdDOF( aDofTypes );
            }

            // if permeability
            if ( tInvPermeabProp != nullptr )
            {
                // if derivative dof type is residual dof type
                if ( aDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // add brinkman term to Jacobian strong form
                    aJM += tInvPermeabProp->val()( 0 ) * tVelocityFI->N();
                }

                // if permeability depends on dof type
                if ( tInvPermeabProp->check_dof_dependency( aDofTypes ) )
                {
                    // add brinkman term to Jacobian strong form
                    aJM += tVelocityFI->val() * tInvPermeabProp->dPropdDOF( aDofTypes );
                }
            }

            // if body load
            if ( tLoadProp != nullptr )
            {
                // if DoF dependency of body load
                if ( tLoadProp->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to momentum Jacobian
                    aJM -= tLoadProp->dPropdDOF( aDofTypes );
                }
            }

            // if mass source
            if ( tMassSourceProp != nullptr )
            {
                // if DoF dependency of source term
                if ( tMassSourceProp->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to continuity Jacobian
                    aJC -= tMassSourceProp->dPropdDOF( aDofTypes ) / tDensity;
                }

                // if density depends on dof type
                if ( tDensityProp->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to Jacobian
                    aJC += tMassSourceProp->val() * tDensityProp->dPropdDOF( aDofTypes ) / std::pow( tDensity, 2 );
                }
            }

            // if CM depends on dof type
            if ( tIncFluidCM->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution to Jacobian strong form
                aJM -= tIncFluidCM->ddivfluxdu( aDofTypes );
            }

            // if gravity
            if ( tGravityProp != nullptr )
            {
                // if gravity depends on dof type
                if ( tGravityProp->check_dof_dependency( aDofTypes ) )
                {
                    // add gravity to Jacobian strong form
                    aJM += tDensity * tGravityProp->dPropdDOF( aDofTypes );
                }

                // if density depends on dof type
                if ( tDensityProp->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to Jacobian
                    aJM += tGravityProp->val() * tDensityProp->dPropdDOF( aDofTypes );
                }

                // if thermal expansion and reference temperature
                if ( tThermalExpProp != nullptr && tRefTempProp != nullptr )
                {
                    // get the temperature field interpolator
                    // FIXME protect FI
                    Field_Interpolator* tTempFI =
                            mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

                    if ( aDofTypes( 0 ) == MSI::Dof_Type::TEMP )
                    {
                        // add contribution to Jacobian
                        aJM -= tDensity * tGravityProp->val() * tThermalExpProp->val() * tTempFI->N();
                    }

                    if ( tThermalExpProp->check_dof_dependency( aDofTypes ) )
                    {
                        // add contribution to Jacobian
                        aJM -= tDensity * tGravityProp->val() * ( tTempFI->val() - tRefTempProp->val() )
                             * tThermalExpProp->dPropdDOF( aDofTypes );
                    }

                    if ( tRefTempProp->check_dof_dependency( aDofTypes ) )
                    {
                        // add contribution to Jacobian
                        aJM += tDensity * tGravityProp->val() * tThermalExpProp->val()
                             * tRefTempProp->dPropdDOF( aDofTypes );
                    }

                    // if gravity property has dependency on the dof type
                    if ( tGravityProp->check_dof_dependency( aDofTypes ) )
                    {
                        // compute the Jacobian
                        aJM -= tDensity * tThermalExpProp->val()( 0 )
                             * ( tTempFI->val()( 0 ) - tRefTempProp->val()( 0 ) )
                             * tGravityProp->dPropdDOF( aDofTypes );
                    }

                    // if density depends on dof type
                    if ( tDensityProp->check_dof_dependency( aDofTypes ) )
                    {
                        // add density contribution to residual strong form
                        aJM -= tGravityProp->val() * tThermalExpProp->val() * ( tTempFI->val() - tRefTempProp->val() )
                             * tDensityProp->dPropdDOF( aDofTypes );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Velocity_Bulk::compute_ujvij( Matrix< DDRMat >& aujvij )
        {
            // get the residual dof type FI (here velocity)
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // init size for uj vij
            uint tNumRow = tVelocityFI->dnNdxn( 1 ).n_rows();
            uint tNumCol = tVelocityFI->dnNdxn( 1 ).n_cols();

            aujvij.set_size( tNumRow, tNumRow * tNumCol, 0.0 );

            // compute: u_j * dv_i/dx_j
            Matrix< DDRMat > tUdVdx = trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 );

            // build test function operator by looping over the number of spatial dimensions
            for ( uint i = 0; i < tNumRow; i++ )
            {
                aujvij( { i, i }, { i * tNumCol, ( i + 1 ) * tNumCol - 1 } ) = tUdVdx.matrix_data();
            }
        }

        //------------------------------------------------------------------------------

        void
        IWG_Incompressible_NS_Velocity_Bulk::compute_ujvijrm( Matrix< DDRMat >& aujvijrm, Matrix< DDRMat >& arm )
        {
            // get the residual dof type FI (here velocity)
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // set size for uj vij rM
            uint tNumField = tVelocityFI->get_number_of_fields();
            uint tNumBases = tVelocityFI->get_number_of_space_time_bases();

            aujvijrm.set_size( tNumField * tNumBases, tNumField * tNumBases, 0.0 );

            // loop over the number of fields
            for ( uint iField = 0; iField < tNumField; iField++ )
            {
                // loop over the number of fields
                for ( uint iField2 = 0; iField2 < tNumField; iField2++ )
                {
                    // compute uj vij rm
                    aujvijrm( { iField * tNumBases, ( iField + 1 ) * tNumBases - 1 },
                            { iField2 * tNumBases, ( iField2 + 1 ) * tNumBases - 1 } ) =
                            trans( tVelocityFI->dnNdxn( 1 ).get_row( iField2 ) ) * tVelocityFI->NBuild()
                            * arm( iField );
                }
            }
        }

        //------------------------------------------------------------------------------

        // FIXME provided directly by the field interpolator?
        void
        IWG_Incompressible_NS_Velocity_Bulk::compute_dnNdtn( Matrix< DDRMat >& adnNdtn )
        {
            // get the residual dof type FI (here velocity)
            Field_Interpolator* tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // init size for dnNdtn
            uint tNumRowt = tVelocityFI->get_number_of_fields();
            uint tNumColt = tVelocityFI->dnNdtn( 1 ).n_cols();

            adnNdtn.set_size( tNumRowt, tNumRowt * tNumColt, 0.0 );

            // loop over the fields
            for ( uint iField = 0; iField < tNumRowt; iField++ )
            {
                // fill the matrix for each dimension
                adnNdtn( { iField, iField }, { iField * tNumColt, ( iField + 1 ) * tNumColt - 1 } ) =
                        tVelocityFI->dnNdtn( 1 ).matrix_data();
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

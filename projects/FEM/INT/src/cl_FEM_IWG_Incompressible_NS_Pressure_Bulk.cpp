
#include "cl_FEM_IWG_Incompressible_NS_Pressure_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------
        IWG_Incompressible_NS_Pressure_Bulk::IWG_Incompressible_NS_Pressure_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Gravity" ]          = IWG_Property_Type::GRAVITY;
            mPropertyMap[ "ThermalExpansion" ] = IWG_Property_Type::THERMAL_EXPANSION;
            mPropertyMap[ "ReferenceTemp" ]    = IWG_Property_Type::REF_TEMP;
            mPropertyMap[ "MassSource" ]       = IWG_Property_Type::MASS_SOURCE;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "IncompressibleFluid" ] = IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID;
            mConstitutiveMap[ "TurbulenceFluid" ]     = IWG_Constitutive_Type::TURBULENCE_FLUID;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "IncompressibleFlow" ] = IWG_Stabilization_Type::INCOMPRESSIBLE_FLOW;
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Incompressible_NS_Pressure_Bulk::set_property - Unknown aPropertyString: " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(), tErrMsg.c_str() );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Incompressible_NS_Pressure_Bulk::set_property - No slave allowed." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster )
        {
            // check that aConstitutiveString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Incompressible_NS_Pressure_Bulk::set_constitutive_model - Unknown aConstitutiveString: " ) +
                    aConstitutiveString;
            MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(), tErrMsg.c_str() );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Incompressible_NS_Pressure_Bulk::set_constitutive_model - No slave allowed." );

            // set the constitutive model in the constitutive model cell
            this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString )
        {
            // check that aStabilizationString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Incompressible_NS_Pressure_Bulk::set_stabilization_parameter - Unknown aStabilizationString: " ) +
                    aStabilizationString;
            MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(), tErrMsg.c_str() );

            // set the stabilization parameter in the stabilization parameter cell
            this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here pressure), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the pressure and velocity FIs
            Field_Interpolator * tPressureFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
            Field_Interpolator * tVelocityFI = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX ); //FIXME this need to be protected

            // get the source property
            std::shared_ptr< Property > tPropSource =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::MASS_SOURCE ) );

            // get the incompressible fluid constitutive model
            std::shared_ptr< Constitutive_Model > tIncFluidCM =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the density
            real tDensity = tIncFluidCM->get_property( "Density" )->val()(0);

            // get the incompressible flow stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tIncFlowSP =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::INCOMPRESSIBLE_FLOW ) );

            // compute residual strong form
            Matrix< DDRMat > tRM;
            this->compute_residual_strong_form( tRM );

            // compute the residual weak form
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) += aWStar * (
                    trans( tPressureFI->N() ) * tVelocityFI->div()
                    + trans( tPressureFI->dnNdxn( 1 ) ) * tIncFlowSP->val()( 0 ) * tRM );

            // if source term is defined
            if (tPropSource != nullptr)
            {
                // add contribution of source term
                mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) += aWStar * (
                        trans( tPressureFI->N() ) * tPropSource->val()( 0 ) / tDensity );
            }

        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here pressure), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the pressure FIs
            Field_Interpolator * tPressureFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the source property
            std::shared_ptr< Property > tPropSource =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::MASS_SOURCE ) );

            // get the incompressible fluid constitutive model
            std::shared_ptr< Constitutive_Model > tIncFluidCM =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the density
            real tDensity = tIncFluidCM->get_property( "Density" )->val()( 0 );

            // get the incompressible flow stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tIncFlowSP =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::INCOMPRESSIBLE_FLOW ) );

            // compute the residual strong form
            Matrix< DDRMat > tRM;
            compute_residual_strong_form( tRM );

            // compute the jacobian for dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if dof type is velocity
                // FIXME protect velocity dof type
                if( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    Field_Interpolator * tVelocityFI = mMasterFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

                    // compute the jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tPressureFI->N() ) * tVelocityFI->div_operator() );
                }

                // compute the jacobian strong form
                Matrix< DDRMat > tJM;
                compute_jacobian_strong_form( tDofType, tJM );

                // compute the jacobian contribution from strong form
                mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                trans( tPressureFI->dnNdxn( 1 ) ) * tIncFlowSP->val()( 0 ) * tJM );

                // if source term has dependency on the dof type
                if ( tPropSource != nullptr )
                {
                    // if source property depends on dof type
                    if ( tPropSource->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                        trans( tPressureFI->N() ) * 1/tDensity * tPropSource->dPropdDOF( tDofType ) );
                    }
                }

                // if stabilization parameter has dependency on the dof type
                if ( tIncFlowSP->check_dof_dependency( tDofType ) )
                {
                    // compute the jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tPressureFI->dnNdxn( 1 ) ) * tRM * tIncFlowSP->dSPdMasterDOF( tDofType ).get_row( 0 ) );
                }
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Pressure_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::compute_residual_strong_form( Matrix< DDRMat > & aRM )
        {
            // get the velocity and pressure FIs
            Field_Interpolator * tVelocityFI = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX ); //FIXME

            // get the density and gravity properties
            std::shared_ptr< Property > tGravityProp = mMasterProp( static_cast< uint >( IWG_Property_Type::GRAVITY ) );
            std::shared_ptr< Property > tThermalExpProp = mMasterProp( static_cast< uint >( IWG_Property_Type::THERMAL_EXPANSION ) );
            std::shared_ptr< Property > tRefTempProp    = mMasterProp( static_cast< uint >( IWG_Property_Type::REF_TEMP ) );

            // get the incompressible fluid constitutive model
            std::shared_ptr< Constitutive_Model > tIncFluidCM =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the turbulence fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMTurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::TURBULENCE_FLUID ) );

            // get the incompressible flow stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tIncFlowSP =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::INCOMPRESSIBLE_FLOW ) );

            // get the density property from CM
            std::shared_ptr< Property > tDensityProp = tIncFluidCM->get_property( "Density" );

            // get the density value
            real tDensity = tDensityProp->val()( 0 );

            // compute the residual strong form
            aRM =
                    tDensity * trans( tVelocityFI->gradt( 1 ) ) +
                    tDensity * trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->val() -
                    tIncFluidCM->divflux();

            // if gravity
            if ( tGravityProp != nullptr )
            {
                // add gravity to residual strong form
                aRM.matrix_data() += tDensity * tGravityProp->val();

                // if thermal expansion and reference temperature
                if( tThermalExpProp != nullptr && tRefTempProp != nullptr )
                {
                    // get the temperature field interpolator
                    // FIXME protect FI
                    Field_Interpolator * tTempFI =
                            mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

                    // add contribution to residual
                    aRM.matrix_data() -=
                            tDensity * tGravityProp->val() * tThermalExpProp->val() *
                            ( tTempFI->val() - tRefTempProp->val() );
                }
            }

            if( tCMTurbulence != nullptr )
            {
                // add contribution to residual
                aRM.matrix_data() -= tCMTurbulence->divflux().matrix_data();
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::compute_jacobian_strong_form(
                moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >             & aJM )
        {
            // get the velocity dof and the derivative dof FIs
            Field_Interpolator * tVelocityFI = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX ); //FIXME
            Field_Interpolator * tDerFI      = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init aJM and aJC
            aJM.set_size(
                    tVelocityFI->get_number_of_fields(),
                    tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the gravity properties
            std::shared_ptr< Property > tGravityProp = mMasterProp( static_cast< uint >( IWG_Property_Type::GRAVITY ) );
            std::shared_ptr< Property > tThermalExpProp = mMasterProp( static_cast< uint >( IWG_Property_Type::THERMAL_EXPANSION ) );
            std::shared_ptr< Property > tRefTempProp    = mMasterProp( static_cast< uint >( IWG_Property_Type::REF_TEMP ) );

            // get the incompressible fluid constitutive model
            std::shared_ptr< Constitutive_Model > tIncFluidCM =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::INCOMPRESSIBLE_FLUID ) );

            // get the turbulence fluid constitutive model
            std::shared_ptr< Constitutive_Model > tCMTurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::TURBULENCE_FLUID ) );

            // get the incompressible flow stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tIncFlowSP =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::INCOMPRESSIBLE_FLOW ) );

            // get the density property from CM
            std::shared_ptr< Property > tDensityProp = tIncFluidCM->get_property( "Density" );

            // get the density value
            real tDensity = tDensityProp->val()( 0 );

            // if derivative wrt to velocity dof type
            // FIXME protect the velocity dof type
            if( aDofTypes( 0 ) == MSI::Dof_Type::VX )
            {
                // compute the term uj vij
                Matrix< DDRMat > tujvij;
                this->compute_ujvij( tujvij );

                // compute the term dnNdtn
                Matrix< DDRMat > tdnNdtn;
                this->compute_dnNdtn( tdnNdtn );

                // compute the jacobian strong form
                aJM.matrix_data() +=
                        tDensity * tdnNdtn
                        + tDensity * trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->N()
                        + tDensity * tujvij ;
            }

            // if density depends on dof type
            if( tDensityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution to jacobian strong form
                aJM.matrix_data() +=
                        trans( tVelocityFI->gradt( 1 ) ) * tDensityProp->dPropdDOF( aDofTypes )
                        + trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->val() * tDensityProp->dPropdDOF( aDofTypes );
            }

            // if CM depends on dof type
            if( tIncFluidCM->check_dof_dependency( aDofTypes ) )
            {
                // add contribution to jacobian
                aJM.matrix_data() -= tIncFluidCM->ddivfluxdu( aDofTypes ).matrix_data();
            }

            // if gravity
            if ( tGravityProp != nullptr )
            {
                // if gravity depends on dof type
                if( tGravityProp->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to jacobian
                    aJM.matrix_data() += tDensity * tGravityProp->dPropdDOF( aDofTypes ).matrix_data();
                }

                // if density depends on dof type
                if( tDensityProp->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to jacobian
                    aJM.matrix_data() += tGravityProp->val() * tDensityProp->dPropdDOF( aDofTypes );
                }

                // if thermal expansion and reference temperature
                if( tThermalExpProp != nullptr && tRefTempProp != nullptr )
                {
                    // get the temperature field interpolator
                    // FIXME protect FI
                    Field_Interpolator * tTempFI
                    = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

                    // if dof type is temperature
                    if( aDofTypes( 0 ) == MSI::Dof_Type::TEMP )
                    {
                        // add contribution to jacobian
                        aJM.matrix_data() -= tDensity * tGravityProp->val() * tThermalExpProp->val() * tTempFI->N();
                    }

                    // if thermal expansion property depends on dof type
                    if( tThermalExpProp->check_dof_dependency( aDofTypes ) )
                    {
                        // add contribution to jacobian
                        aJM.matrix_data() -=
                                tDensity * tGravityProp->val() *
                                ( tTempFI->val() - tRefTempProp->val() ) * tThermalExpProp->dPropdDOF( aDofTypes );
                    }

                    // if reference temperature property depends on dof type
                    if( tRefTempProp->check_dof_dependency( aDofTypes ) )
                    {
                        // add contribution to jacobian
                        aJM.matrix_data() +=
                                tDensity * tGravityProp->val() *
                                tThermalExpProp->val() * tRefTempProp->dPropdDOF( aDofTypes );
                    }

                    // if gravity property has dependency on the dof type
                    if ( tGravityProp->check_dof_dependency( aDofTypes ) )
                    {
                        // compute the jacobian
                        aJM.matrix_data() -=
                                tDensity *
                                tThermalExpProp->val()( 0 ) * ( tTempFI->val()( 0 ) - tRefTempProp->val()( 0 ) ) *
                                tGravityProp->dPropdDOF( aDofTypes );
                    }

                    // if density depends on dof type
                    if( tDensityProp->check_dof_dependency( aDofTypes ) )
                    {
                        // add density contribution to residual strong form
                        aJM.matrix_data() -=
                                tGravityProp->val() *
                                tThermalExpProp->val() * ( tTempFI->val() - tRefTempProp->val() ) *
                                tDensityProp->dPropdDOF( aDofTypes );
                    }
                }

                // if turbulence
                if( tCMTurbulence != nullptr )
                {
                    // if turbulence CM depends on dof type
                    if( tCMTurbulence->check_dof_dependency( aDofTypes ) )
                    {
                        // compute contribution to jacobian strong form
                        aJM.matrix_data() -= tCMTurbulence->ddivfluxdu( aDofTypes ).matrix_data();
                    }
                }
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::compute_ujvij( Matrix< DDRMat > & aujvij )
        {
            // get the velocity dof type FI
            // FIXME protect velocity dof type
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // init size for uj vij
            uint tNumRow = tVelocityFI->dnNdxn( 1 ).n_rows();
            uint tNumCol = tVelocityFI->dnNdxn( 1 ).n_cols();
            aujvij.set_size( tNumRow, tNumRow * tNumCol, 0.0 );

            // loop over the number of rows
            for( uint i = 0; i < tNumRow; i++ )
            {
                // compute uj vij
                aujvij( { i, i }, { i * tNumCol, ( i + 1 ) * tNumCol -1 } ) =
                        trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 );
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Incompressible_NS_Pressure_Bulk::compute_ujvijrm(
                Matrix< DDRMat > & aujvijrm,
                Matrix< DDRMat > & arm )
        {
            // get the velocity dof type FI
            // FIXME protect velocity dof type
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // set size for uj vij rM
            uint tNumField = tVelocityFI->get_number_of_fields();
            uint tNumBases = tVelocityFI->get_number_of_space_time_bases();
            aujvijrm.set_size( tNumField * tNumBases, tNumField * tNumBases , 0.0 );

            // loop over the number of fields
            for( uint iField = 0; iField < tNumField; iField++ )
            {
                // loop over the number of fields
                for( uint iField2 = 0; iField2 < tNumField; iField2++ )
                {
                    // compute uj vij rm
                    aujvijrm(
                            { iField  * tNumBases, ( iField + 1 )  * tNumBases - 1 },
                            { iField2 * tNumBases, ( iField2 + 1 ) * tNumBases - 1 } ) =
                                    trans( tVelocityFI->dnNdxn( 1 ).get_row( iField2 ) ) * tVelocityFI->NBuild() * arm( iField );
                }
            }
        }

        //------------------------------------------------------------------------------
        // FIXME provided directly by the field interpolator?
        void IWG_Incompressible_NS_Pressure_Bulk::compute_dnNdtn( Matrix< DDRMat > & adnNdtn )
        {
            // get the velocity dof type FI
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // init size for dnNdtn
            uint tNumRowt = tVelocityFI->get_number_of_fields();
            uint tNumColt = tVelocityFI->dnNdtn( 1 ).n_cols();
            adnNdtn.set_size( tNumRowt, tNumRowt * tNumColt , 0.0 );

            // loop over the fields
            for( uint iField = 0; iField < tNumRowt; iField++ )
            {
                // fill the matrix for each dimension
                adnNdtn( { iField, iField }, { iField * tNumColt, ( iField + 1 ) * tNumColt - 1 } ) =
                        tVelocityFI->dnNdtn( 1 ).matrix_data();
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */






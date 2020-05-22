//FEM/INT/src
#include "cl_FEM_IWG_Diffusion_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------
        IWG_Diffusion_Bulk::IWG_Diffusion_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Load" ]         = IWG_Property_Type::BODY_LOAD;
            mPropertyMap[ "Density" ]      = IWG_Property_Type::DENSITY;
            mPropertyMap[ "HeatCapacity" ] = IWG_Property_Type::HEAT_CAPACITY;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Diffusion" ] = IWG_Constitutive_Type::DIFFUSION;
        }

        //------------------------------------------------------------------------------
        void IWG_Diffusion_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get residual dof type field interpolator (here temperature)
            Field_Interpolator * tFITemp = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get body load property
            std::shared_ptr< Property > tPropLoad =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            // get density property
            std::shared_ptr< Property > tPropDensity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DENSITY ) );

            // get heat capacity property
            std::shared_ptr< Property > tPropHeatCapacity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::HEAT_CAPACITY ) );

            // get the elasticity CM
            std::shared_ptr< Constitutive_Model > tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // compute the residual
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
                    += aWStar * ( trans( tCMDiffusion->testStrain() ) * tCMDiffusion->flux() );

            // if body load
            if ( tPropLoad != nullptr )
            {
                // compute contribution of body load to residual
                mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
                        -= aWStar * ( trans( tFITemp->N() ) * tPropLoad->val()( 0 ) );
            }

            // if density and heat capacity
            if( tPropDensity != nullptr && tPropHeatCapacity != nullptr )
            {
                // add contribution to residual
                mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
                        += aWStar * ( tPropDensity->val()( 0 ) * tPropHeatCapacity->val()( 0 ) * trans( tFITemp->N() ) * tFITemp->gradt( 1 ) );
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Diffusion_Bulk::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolator for a given dof type
            Field_Interpolator * tFITemp = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get body load property
            std::shared_ptr< Property > tPropLoad =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            // get density property
            std::shared_ptr< Property > tPropDensity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DENSITY ) );

            // get heat capacity property
            std::shared_ptr< Property > tPropHeatCapacity =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::HEAT_CAPACITY ) );

            // get the elasticity CM
            std::shared_ptr< Constitutive_Model > tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get the number of master dof type dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over master dof type dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if body load
                if ( tPropLoad != nullptr )
                {
                    // if body load property has dependency on the dof type
                    if ( tPropLoad->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } )
                                -= aWStar * ( trans( tFITemp->N() ) * tPropLoad->dPropdDOF( tDofType ) );
                    }
                }

                if( tPropDensity != nullptr && tPropHeatCapacity != nullptr )
                {
                    // if residual dof type (here temperature)
                    if( tDofType( 0 ) == mResidualDofType( 0 ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } )
                                += aWStar * ( tPropDensity->val()( 0 ) * tPropHeatCapacity->val()( 0 ) * trans( tFITemp->N() ) * tFITemp->dnNdtn( 1 ) );
                    }

                    // if density property has dependency on the dof type
                    if ( tPropDensity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } )
                                += aWStar * ( tPropHeatCapacity->val()( 0 ) * trans( tFITemp->N() ) * tFITemp->gradt( 1 ) * tPropDensity->dPropdDOF( tDofType ) );
                    }

                    // if heat capacity property has dependency on the dof type
                    if ( tPropHeatCapacity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } )
                                += aWStar * ( tPropDensity->val()( 0 ) * trans( tFITemp->N() ) * tFITemp->gradt( 1 ) * tPropHeatCapacity->dPropdDOF( tDofType ) );
                    }
                }

                // if constitutive model has dependency on the dof type
                if ( tCMDiffusion->check_dof_dependency( tDofType ) )
                {
                    // compute the jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } )
                            += aWStar * ( trans( tCMDiffusion->testStrain() ) * tCMDiffusion->dFluxdDOF( tDofType ) );
                    // FIXME add derivative of the test strain
                }
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Diffusion_Bulk::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Diffusion_Bulk::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Bulk::compute_dRdp - Not implemented." );

            //            // get master index for residual dof type, indices for assembly
            //            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            //            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            //            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );
            //
            //            // get residual dof type field interpolator (here temperature)
            //            Field_Interpolator * tFITemp = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
            //
            //            // get body load property
            //            std::shared_ptr< Property > tPropLoad
            //            = mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );
            //
            //            // get density property
            //            std::shared_ptr< Property > tPropDensity
            //            = mMasterProp( static_cast< uint >( IWG_Property_Type::DENSITY ) );
            //
            //            // get heat capacity property
            //            std::shared_ptr< Property > tPropHeatCapacity
            //            = mMasterProp( static_cast< uint >( IWG_Property_Type::HEAT_CAPACITY ) );
            //
            //            // get the elasticity CM
            //            std::shared_ptr< Constitutive_Model > tCMDiffusion
            //            = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );
            //
            //            // get number of master dv dependencies
            //            uint tNumDvDependencies = mMasterGlobalDvTypes.size();
            //            // loop over master dv dependencies
            //            for( uint iDv = 0; iDv < tNumDvDependencies; iDv++ )
            //            {
            //                // get the treated dv type
            //                Cell< PDV_Type > tDvType = mMasterGlobalDvTypes( iDv );
            //
            //                // get the index for dof type, indices for assembly
            //                sint tDvDepIndex          = ;
            //                uint tMasterDepStartIndex = ;
            //                uint tMasterDepStopIndex  = ;
            //
            //                // get index for the treated dv type
            //                uint tIndexDep = mSet->get_dv_index_for_type( tDvType( 0 ), mtk::Master_Slave::MASTER );
            //
            //                // if body load
            //                if ( tPropLoad != nullptr )
            //                {
            //                    // if load property has dependency on the dv type
            //                    if ( tPropLoad->check_dv_dependency( tDvType ) )
            //                    {
            //                        // compute drdpdv
            //                        mSet->get_drdpdv()( { tMasterResStartIndex, tMasterResStopIndex },
            //                                            { tMasterDepStartIndex, tMasterDepStopIndex } )
            //                        -= aWStar * ( trans( tFI->N() ) * tPropLoad->dPropdDV( tDvType ) );
            //                    }
            //                }
            //
            //                // if diffusion constitutive model has dependency on the dv type
            //                if ( tCMDiffusion->check_dv_dependency( tDvType ) )
            //                {
            //                    // compute the jacobian
            //                    // compute the jacobian
            //                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
            //                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
            //                    += aWStar * ( trans( tCMDiffusion->testStrain() ) * tCMDiffusion->dFluxdDV( tDvType ) );
            //                }
            //            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

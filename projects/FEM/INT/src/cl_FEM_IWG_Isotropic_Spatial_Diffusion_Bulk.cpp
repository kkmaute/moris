
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Bulk.hpp"
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
        IWG_Isotropic_Spatial_Diffusion_Bulk::IWG_Isotropic_Spatial_Diffusion_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Load" ] = IWG_Property_Type::LOAD;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "DiffLinIso" ] = IWG_Constitutive_Type::DIFF_LIN_ISO;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_residual( real tWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get indices for SP, CM and property
            uint tDiffLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );
            uint tLoadIndex       = static_cast< uint >( IWG_Property_Type::LOAD );

            // compute the residual
            mSet->get_residual()( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
            += trans( mMasterCM( tDiffLinIsoIndex )->testStrain() ) * mMasterCM( tDiffLinIsoIndex )->flux() * tWStar;

            // if body load
            if ( mMasterProp( tLoadIndex ) != nullptr )
            {
                // get field interpolator for a given dof type
                Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

                // compute contribution of body load to residual
                mSet->get_residual()( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
                += - trans( tFI->N() ) * mMasterProp( tLoadIndex )->val()( 0 ) * tWStar;
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian( real tWStar )
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
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get indices for SP, CM and property
            uint tDiffLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );
            uint tLoadIndex       = static_cast< uint >( IWG_Property_Type::LOAD );

            // compute the jacobian for direct dof dependencies
            // Here no direct dependencies

            // compute the jacobian for indirect dof dependencies through properties and constitutive model
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if body load
                if ( mMasterProp( tLoadIndex ) != nullptr )
                {
                    // if property has dependency on the dof type
                    if ( mMasterProp( tLoadIndex )->check_dof_dependency( tDofType ) )
                    {
                        // compute the jacobian
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                              { tMasterDepStartIndex, tMasterDepStopIndex } )
                        -= trans( tFI->N() ) * mMasterProp( tLoadIndex )->dPropdDOF( tDofType ) * tWStar;
                    }
                }

                // if constitutive model has dependency on the dof type
                if ( mMasterCM( tDiffLinIsoIndex )->check_dof_dependency( tDofType ) )
                {
                    // compute the jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += ( trans( mMasterCM( tDiffLinIsoIndex )->testStrain() ) * mMasterCM( tDiffLinIsoIndex )->dFluxdDOF( tDofType ) ) * tWStar;
                    // fixme add derivative of the test strain
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_drdpdv( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

//            // get index for a given dof type
//            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
//
//            // get field interpolator for a given dof type
//            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
//
//            // get indices for SP, CM and property
//            uint tDiffLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );
//            uint tLoadIndex       = static_cast< uint >( IWG_Property_Type::LOAD );
//
//            // compute the drdpdv for indirect dv dependencies
//            // through properties, constitutive models and stabilization parameters
//            uint tNumDvDependencies = mMasterGlobalDvTypes.size();
//            for( uint iDv = 0; iDv < tNumDvDependencies; iDv++ )
//            {
//                // get the treated dv type
//                Cell< MSI::Dv_Type > tDvType = mMasterGlobalDvTypes( iDv );
//
//                // get index for the treated dv type
//                uint tIndexDep = mSet->get_dv_index_for_type( tDvType( 0 ), mtk::Master_Slave::MASTER );
//
//                // if body load
//                if ( mMasterProp( tLoadIndex ) != nullptr )
//                {
//                    // if property has dependency on the dv type
//                    if ( mMasterProp( tLoadIndex )->check_dv_dependency( tDvType ) )
//                    {
//                        // compute drdpdv
//                        mSet->get_drdpdv()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
//                                            { mSet->get_dv_assembly_map()( tDvIndex )( tIndexDep, 0 ), mSet->get_dv_assembly_map()( tDvIndex )( tIndexDep, 1 ) } )
//                        += - trans( tFI->N() ) * mMasterProp( tLoadIndex )->dPropdDV( tDvType ) * tWStar;
//                    }
//                }
//
//                // if constitutive model has dependency on the dv type
//                if ( mMasterCM( tDiffLinIsoIndex )->check_dv_dependency( tDvType ) )
//                {
//                    // compute the jacobian
//                    // compute the jacobian
//                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
//                                          { mSet->get_dv_assembly_map()( tDvIndex )( tIndexDep, 0 ), mSet->get_dv_assembly_map()( tDvIndex )( tIndexDep, 1 ) } )
//                    += ( trans( mMasterCM( tDiffLinIsoIndex )->testStrain() ) * mMasterCM( tDiffLinIsoIndex )->dFluxdDV( tDvType ) ) * tWStar;
//                }
//            }

            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Bulk::compute_drdpdv - Not implemented." );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


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
            // check master field interpolators, properties and constitutive models
#ifdef DEBUG
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
#endif

            // get index for a given dof type
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get start and end indices for residual assembly
            uint tStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get indices for SP, CM and property
            uint tDiffLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );
            uint tLoadIndex       = static_cast< uint >( IWG_Property_Type::LOAD );

            // compute the residual
            mSet->get_residual()( { tStartRow, tEndRow }, { 0, 0 } )
            += trans( mMasterCM( tDiffLinIsoIndex )->testStrain() ) * mMasterCM( tDiffLinIsoIndex )->flux() * tWStar;

            // if body load
            if ( mMasterProp( tLoadIndex ) != nullptr )
            {
                // get field interpolator for a given dof type
                Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

                // compute contribution of body load to residual
                mSet->get_residual()( { tStartRow, tEndRow }, { 0, 0 } )
                += - trans( tFI->N() ) * mMasterProp( tLoadIndex )->val()( 0 ) * tWStar;
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian( real tWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
#endif

            // get index for a given dof type
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get field interpolator for a given dof type
            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

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

                // get index for the treated dof type
                uint tIndexDep = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iDOF )( 0 ), mtk::Master_Slave::MASTER );

                // if body load
                if ( mMasterProp( tLoadIndex ) != nullptr )
                {
                    // if property has dependency on the dof type
                    if ( mMasterProp( tLoadIndex )->check_dof_dependency( tDofType ) )
                    {
                        // compute the jacobian
                        mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                              { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                        += - trans( tFI->N() ) * mMasterProp( tLoadIndex )->dPropdDOF( tDofType ) * tWStar;
                    }
                }

                // if constitutive model has dependency on the dof type
                if ( mMasterCM( tDiffLinIsoIndex )->check_dof_dependency( tDofType ) )
                {
                    // compute the jacobian
                    // compute the jacobian
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
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
    } /* namespace fem */
} /* namespace moris */

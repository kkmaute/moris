
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_residual( real tWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
//            this->check_properties();
//            this->check_constitutive_models();
#endif

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            uint tStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // compute the residual
            mSet->get_residual()( { tStartRow, tEndRow }, { 0, 0 } )
                    += trans( mMasterCM( 0 )->testStrain() ) * mMasterCM( 0 )->flux() * tWStar;

            // if body load
            if ( mMasterProp.size() > 0 )
            //if ( mMasterGlobalPropTypes( 0 ) == fem::Property_Type::TEMP_LOAD )
            {
                mSet->get_residual()( { tStartRow, tEndRow }, { 0, 0 } )
                        += - trans( tFI->N() ) * mMasterProp( 0 )->val()( 0 ) * tWStar;
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian( real tWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
//            this->check_properties();
//            this->check_constitutive_models();
#endif

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // compute the jacobian for direct dof dependencies
            // Here no direct dependencies

            // compute the jacobian for indirect dof dependencies through properties and constitutive model
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                uint tIndexDep = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iDOF )( 0 ), mtk::Master_Slave::MASTER );

                // if we have a body load FIXME
                if ( mMasterProp.size() > 0 )
                //if ( mMasterGlobalPropTypes( 0 ) == fem::Property_Type::TEMP_LOAD )
                {
                    // if property has dependency on the dof type
                    if ( mMasterProp( 0 )->check_dof_dependency( tDofType ) )
                    {
                        // compute the jacobian
                        mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                              { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                                += - trans( tFI->N() ) * mMasterProp( 0 )->dPropdDOF( tDofType ) * tWStar;
                    }
                }

                // if constitutive model has dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( tDofType ) )
                {
                    // compute the jacobian
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                            += trans( mMasterCM( 0 )->testStrain() ) * mMasterCM( 0 )->dFluxdDOF( tDofType ) * tWStar;
                    // fixme add derivative of the test strain
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                             moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Bulk.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_CM_Struc_Linear_Isotropic.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
    //------------------------------------------------------------------------------
        IWG_Isotropic_Struc_Linear_Pressure_Bulk::IWG_Isotropic_Struc_Linear_Pressure_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = IWG_Constitutive_Type::ELAST_LIN_ISO;
            mConstitutiveMap[ "ElastLinIsoPressure" ] = IWG_Constitutive_Type::ELAST_LIN_ISO_PRESSURE;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Pressure_Bulk::compute_residual( real aWStar )
        {

#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
#endif

            // get index for a given dof type
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get field interpolator for dof type
            Field_Interpolator* tDisplacementFI = mFieldInterpolatorManager->get_field_interpolators_for_type( MSI::Dof_Type::UX, mtk::Master_Slave::MASTER );
            Field_Interpolator* tPressureFI = mFieldInterpolatorManager->get_field_interpolators_for_type( MSI::Dof_Type::P, mtk::Master_Slave::MASTER );

            // get start and end index for residual assembly
            uint tStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get property, CM, SP indices
            uint tElastLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
            moris::fem::CM_Struc_Linear_Isotropic* tLinearIso = static_cast <moris::fem::CM_Struc_Linear_Isotropic*> (mMasterCM(tElastLinIsoIndex).get());

            // compute the residual
            mSet->get_residual()( { tStartRow, tEndRow }, { 0, 0 } )
            +=  trans(tPressureFI->N()) * (tDisplacementFI->div() + tLinearIso->eval_inv_bulk_modulus() * tPressureFI->val()(0)) * aWStar;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Pressure_Bulk::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
#endif

            // get index for given dof type
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get field interpolator for given dof type
            Field_Interpolator* tPressureFI = mFieldInterpolatorManager->get_field_interpolators_for_type( MSI::Dof_Type::P, mtk::Master_Slave::MASTER );

            // get property, CM, SP indices
            uint tElastLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
            uint tElastLinIsoPressureIndex = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO_PRESSURE );
            moris::fem::CM_Struc_Linear_Isotropic* tLinearIso = static_cast <moris::fem::CM_Struc_Linear_Isotropic*> (mMasterCM(tElastLinIsoIndex).get());

            // compute the jacobian for direct dof dependencies
            // Here no direct dependencies

            // compute the jacobian for indirect dof dependencies through properties and constitutive model
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );
                uint tIndexDep = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iDOF )( 0 ), mtk::Master_Slave::MASTER );

                // if displacement constitutive model has dependency on the dof type (if displacement)
                if ( mMasterCM( tElastLinIsoIndex )->check_dof_dependency( tDofType ) )
                {
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                    += trans(tPressureFI->N()) * tLinearIso->eval_B_flat() * aWStar;
                }

                // if pressure constitutive model has dependency on the dof type (if pressure)
                if ( mMasterCM( tElastLinIsoPressureIndex )->check_dof_dependency( tDofType ) )
                {
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                            += trans(tPressureFI->N()) * tLinearIso->eval_inv_bulk_modulus() * tPressureFI->N() * aWStar;
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Pressure_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                             moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Pressure_Bulk::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

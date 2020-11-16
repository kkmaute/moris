
#include "cl_FEM_IWG_Heat_Method_Boundary.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        void IWG_Heat_Method_Boundary::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get field interpolators
            Field_Interpolator * tFIPhiD = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
            Field_Interpolator * tFITheta = mMasterFIManager->get_field_interpolators_for_type( mDofTheta );

            // compute normalized gradient of Theta
            Matrix< DDRMat > tGradTheta = tFITheta->gradx( 1 );
            Matrix< DDRMat > tNormGradTheta = tGradTheta / norm( tGradTheta );

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tResStartIndex, tResStopIndex },
                    { 0, 0 } ) -= aWStar * (
                            trans( tFIPhiD->N() ) * dot( tNormGradTheta, mNormal ) );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Heat_Method_Boundary::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Heat_Method_Boundary::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get field interpolators
            Field_Interpolator * tFIPhiD = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
            Field_Interpolator * tFITheta = mMasterFIManager->get_field_interpolators_for_type( mDofTheta );

            // compute the jacobian for dof dependencies
            for( uint iDOF = 0; iDOF < mRequestedMasterGlobalDofTypes.size(); iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the dof type indices for assembly
                uint tDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepDofIndex, 0 );
                uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDepDofIndex, 1 );

                // Dof Deriv wrt. to Theta
                if( tDofType( 0 ) == mDofTheta )
                {
                    // compute normalized gradient of Theta
                    Matrix< DDRMat > tBTheta = tFITheta->dnNdxn( 1 );
                    Matrix< DDRMat > tGradTheta = tFITheta->gradx( 1 );
                    real tNorm = norm( tGradTheta );
                    Matrix< DDRMat > tNormGradTheta = tGradTheta / tNorm;
                    Matrix< DDRMat > tNormBTheta = tBTheta / tNorm;

                    // compute dof derivative
                    Matrix< DDRMat > tdHeatFluxdDOF = tNormBTheta - tNormGradTheta * trans( tNormGradTheta ) * tNormBTheta;

                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tResStartIndex, tResStopIndex },
                            { tDepStartIndex, tDepStopIndex } ) -= aWStar * (
                                    trans( tFIPhiD->N() ) * trans( mNormal ) * tdHeatFluxdDOF );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Heat_Method_Boundary::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Heat_Method_Boundary::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, " IWG_Heat_Method_Boundary::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Heat_Method_Boundary::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Heat_Method_Boundary::compute_dRdp - Not implemented.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

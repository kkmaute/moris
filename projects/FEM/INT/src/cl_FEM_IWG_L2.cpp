#include "cl_FEM_IWG_L2.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "op_times.hpp" //LINALG/src
#include "fn_norm.hpp"  //LINALG/src
#include "fn_trans.hpp" //LINALG/src
#include "fn_dot.hpp"   //LINALG/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
    IWG_L2::IWG_L2( const real aAlpha )
    {
//        // set the residual dof type
//        mResidualDofType = { MSI::Dof_Type::L2 };
//
//        // set the active dof types
//        mMasterDofTypes = { { MSI::Dof_Type::L2 } };

        // set alpha
        this->set_alpha( aAlpha );
    }

//------------------------------------------------------------------------------

    void IWG_L2::set_alpha( const real aAlpha )
    {
        mAlpha = aAlpha;

        if(  aAlpha == 0.0 )
        {
            mComputeFunction
                = & IWG_L2::compute_jacobian_and_residual_without_alpha;

            mComputeJacFunction
                = & IWG_L2::compute_jacobian_without_alpha;

            mComputeResFunction
                = & IWG_L2::compute_residual_without_alpha;
        }
        else
        {
            mComputeFunction
                = & IWG_L2::compute_jacobian_and_residual_with_alpha;

            mComputeJacFunction
                = & IWG_L2::compute_jacobian_with_alpha;

            mComputeResFunction
                = & IWG_L2::compute_residual_with_alpha;
        }
    }

//------------------------------------------------------------------------------

    void IWG_L2::compute_jacobian_and_residual( real aWStar )
    {
        // call the right mComputeFunction for residual and jacobian evaluations
        ( this->*mComputeFunction )( aWStar );
    }

//------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_and_residual_without_alpha( real aWStar )
        {
            MORIS_ERROR( false, "will not work because of weights");
//            // check master field interpolators
//            this->check_dof_field_interpolators();
//            this->check_dv_field_interpolators();
//
//            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
//
//            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
//
//            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
//
//            // compute Jacobian
//            aJacobians( 0 ).resize( 1 );
//            aJacobians( 0 )( 0 ) = trans( tFI->N() ) * tFI->N();
//
//            // set the jacobian size
//            aResidual.resize( 1 );
//
//            // compute residual
//            //FIXME mNodalWeakBCs
//            mSet->get_residual()( { mSet->get_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } )
//                                += aJacobians( 0 )( 0 ) * ( tFI->get_coeff() - mNodalWeakBCs );

        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_and_residual_with_alpha( real aWStar )
        {
            MORIS_ERROR( false, "will not work because of weights");
            // check master field interpolators
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // compute Jacobian
            mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) },
                                  { mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 2 ), mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 3 ) } )
                    += (   trans( tFI->N() ) * tFI->N()
                         + mAlpha * ( trans( tFI->dnNdxn( 1 ) ) * tFI->dnNdxn( 1 ) ) ) * aWStar;

            // compute residual
            //FIXME: mNodalWeakBCs
            mSet->get_residual()( { mSet->get_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } )
                    += mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) },
                                             { mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 2 ), mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 3 ) } )
                            * ( tFI->get_coeff() - mNodalWeakBCs ) * aWStar;
        }

//------------------------------------------------------------------------------
//
//        real
//        IWG_L2::interpolate_scalar_at_point(
//                                    const Matrix< DDRMat > & aNodalWeakBC,
//                                    const uint             & aPointIndex )
//        {
//            // get shape function
//            mN->compute( aPointIndex );
//
//            // return interpolation
//            return dot( mN->matrix() , aNodalWeakBC );
//        }

//------------------------------------------------------------------------------
//
//        real
//        IWG_L2::compute_integration_error(
//                const Matrix< DDRMat >                    & aNodalDOF,
//                real (*aFunction)( const Matrix< DDRMat > & aPoint ) ,
//                const uint                                & aPointIndex )
//        {
//            mN->compute( aPointIndex );
//
//            Matrix< DDRMat > tCoords = mN->matrix_data() * mInterpolator->get_node_coords();
//
//            // get shape function
//            Matrix< DDRMat > tPhiHat = mN->matrix_data() * aNodalDOF.matrix_data();
//
//            return std::pow( tPhiHat( 0 ) - aFunction( tCoords ), 2 );
//        }
//------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian( real aWStar )
        {
            // call the right mComputeJacFunction for jacobian evaluation
            ( this->*mComputeJacFunction )( aWStar );
        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_without_alpha( real aWStar )
        {
            // check master field interpolators
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // compute Jacobian
            mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) },
                                  { mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 2 ), mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 3 ) } )
                    += trans( tFI->N() ) * tFI->N() * aWStar;
        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_with_alpha( real aWStar )
        {
            // check master field interpolators
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // compute Jacobian
            mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) },
                                  { mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 2 ), mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 3 ) } )
                    += (   trans( tFI->N() ) * tFI->N()
                         + mAlpha * ( trans( tFI->dnNdxn( 1 ) ) * tFI->dnNdxn( 1 ) ) ) * aWStar;
        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_residual( real aWStar )
        {
            // call the right mComputeResFunction for residual evaluation
            ( this->*mComputeResFunction )( aWStar );
        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_residual_without_alpha( real aWStar )
        {
            // check master field interpolators
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // compute residual
            //FIXME: mNodalWeakBCs
            mSet->get_residual()( { mSet->get_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } )
                                += trans( tFI->N() ) * ( tFI->val() - tFI->N() * mNodalWeakBCs ) * aWStar;

        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_residual_with_alpha( real aWStar )
        {
            // check master field interpolators
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // compute residual
            //FIXME mNodalWeakBCs
            mSet->get_residual()( { mSet->get_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } )
                          += ( trans( tFI->N() ) * ( tFI->val() - tFI->N() * mNodalWeakBCs )
                           + mAlpha * trans( tFI->dnNdxn( 1 ) ) * ( tFI->gradx( 1 ) - tFI->dnNdxn( 1 ) * mNodalWeakBCs ) )        * aWStar;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk_Test.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        IWG_Hamilton_Jacobi_Bulk_Test::IWG_Hamilton_Jacobi_Bulk_Test()
        {
//            // set the residual dof type
//            mResidualDofType = { MSI::Dof_Type::LS1 };
//
//            // set the active dof type
//            mMasterDofTypes = {{ MSI::Dof_Type::LS1 }};
        }

//------------------------------------------------------------------------------
        void IWG_Hamilton_Jacobi_Bulk_Test::compute_residual( real tWStar )
        {
            // velocity field value
            Matrix< DDRMat > aVN( 1, 3, 1.0 );

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

           //compute the residual
            mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } )
                    += trans( mMasterFI( 0 )->N() ) * ( mMasterFI( 0 )->gradt( 1 ) + aVN * mMasterFI( 0 )->gradx( 1 ) ) * tWStar;
        }

//------------------------------------------------------------------------------
        void IWG_Hamilton_Jacobi_Bulk_Test::compute_jacobian( real tWStar )
        {
            // velocity field value
            Matrix< DDRMat > aVN( 1, 3, 1.0 );

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // compute the jacobian
            mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                  { mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) } )
                    += trans( mMasterFI( 0 )->N() ) * ( mMasterFI( 0 )->dnNdtn( 1 ) + aVN * mMasterFI( 0 )->dnNdxn( 1 ) ) * tWStar;

        }

//------------------------------------------------------------------------------
        void IWG_Hamilton_Jacobi_Bulk_Test::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                           moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
//            // velocity field value
//            Matrix< DDRMat > aVN( 1, 3, 1.0 );
//
//            // set the residual size
//            this->set_residual( aResidual );
//
//            //compute the residual
//            aResidual( 0 ) = trans( mMasterFI( 0 )->N() ) * ( mMasterFI( 0 )->gradt( 1 ) + aVN * mMasterFI( 0 )->gradx( 1 ) );
//
//            // set the jacobian size
//            this->set_jacobian( aJacobians );
//
//            // compute the jacobian
//            aJacobians( 0 )( 0 ) = trans( mMasterFI( 0 )->N() ) * ( mMasterFI( 0 )->dnNdtn( 1 ) + aVN * mMasterFI( 0 )->dnNdxn( 1 ) );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

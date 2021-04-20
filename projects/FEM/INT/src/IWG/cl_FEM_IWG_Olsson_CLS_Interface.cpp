
#include "cl_FEM_IWG_Olsson_CLS_Interface.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        IWG_Olsson_CLS_Interface::IWG_Olsson_CLS_Interface()
        {
            //FIXME set field upper and lower bound
            mPhiUB = 1.0;
            mPhiLB = 0.0;

            //FIXME set Olsson CLS epsilon parameter
            mEpsilon = 1.0;

            // set the residual dof type
            //FIXME: level set scalar field not UX
            mResidualDofType = { MSI::Dof_Type::LS2 };

            // set the active dof type
            //FIXME: level set scalar field not UX
            //       level set normal field not UY
            mMasterDofTypes = {{ MSI::Dof_Type::LS2 },
                               { MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY, MSI::Dof_Type::NLSZ } };

        }

//------------------------------------------------------------------------------
        void IWG_Olsson_CLS_Interface::compute_residual( real tWStar )
        {
            // check master field interpolators
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            // set field interpolators
            Field_Interpolator* phi  = mMasterFI( 0 );
            Field_Interpolator* nPhi = mMasterFI( 1 );

            //FIXME set the interface normal
            Matrix< DDRMat > aInterfaceNormal( phi->gradx( 1 ).n_cols() , 1, 1.0 );

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );

            //compute the residual
            mSet->get_residual()( 0 )( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } )
                    += trans( phi->N() )
                           * ( ( phi->val()( 0 ) - mPhiLB ) * ( mPhiUB - phi->val()( 0 ) ) - mEpsilon * dot( phi->gradx( 1 ), nPhi->val() ) )
                           * trans( nPhi->val() ) * aInterfaceNormal * tWStar;
        }

//------------------------------------------------------------------------------

        void IWG_Olsson_CLS_Interface::compute_jacobian( real tWStar )
        {
            // check master field interpolators
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            // set field interpolators
            Field_Interpolator* phi  = mMasterFI( 0 );
            Field_Interpolator* nPhi = mMasterFI( 1 );

            //FIXME set the interface normal
            Matrix< DDRMat > aInterfaceNormal( phi->gradx( 1 ).n_cols() , 1, 1.0 );

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );

            // compute the jacobian
            mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                  { mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ) } )
                    += ( trans( phi->N() )
                                 * ( ( mPhiLB + mPhiUB - 2 * phi->val()( 0 ) ) * phi->N() - mEpsilon * trans( nPhi->val() ) * phi->dnNdxn( 1 ) )
                                 * dot( nPhi->val(), aInterfaceNormal ) ) * tWStar;

//            aJacobians( 0 )( 1 ) = ( trans( phi->N() ) * ( ( ( phi->val()( 0 ) - mPhiLB ) * ( mPhiUB - phi->val()( 0 ) )
//                                 - 2 * mEpsilon * dot( phi->gradx( 1 ), nPhi->val() ) ) * nPhi->N() ) * aInterfaceNormal ) * tWStar;
        }

//------------------------------------------------------------------------------

        void IWG_Olsson_CLS_Interface::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                      moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Olsson_CLS_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

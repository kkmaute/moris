
#include "cl_FEM_IWG_Helmholtz_Interface.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Helmholtz_Interface::IWG_Helmholtz_Interface()
        {
            //FIXME set the Helmholtz filter parameter
            mFilterParam = 1.0;

            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::VX };

            // set the active dof type
            mMasterDofTypes = {{ MSI::Dof_Type::VX }};
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Interface::compute_residual( real tWStar )
        {
            // set the field interpolator
            Field_Interpolator* vN = mMasterFI( 0 );

            //FIXME set the interface normal
            Matrix< DDRMat > aInterfaceNormal( vN->gradx( 1 ).n_cols() , 1, 1.0 );


            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // compute the residual
            mSet->get_residual()( { mSet->get_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_dof_assembly_map()( tDofIndex )( 0, 1 ) }, { 0, 0 } )
                    += - mFilterParam * trans( vN->N() ) * trans( vN->gradx( 1 ) ) * aInterfaceNormal * tWStar;
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Interface::compute_jacobian( real tWStar )
        {
            // set the field interpolator
            Field_Interpolator* vN = mMasterFI( 0 );

            //FIXME set the interface normal
            Matrix< DDRMat > aInterfaceNormal( vN->gradx( 1 ).n_cols() , 1, 1.0 );

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // compute the jacobian
            mSet->get_jacobian()( { mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) },
                                  { mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 2 ), mSet->get_dof_assembly_map()( tDofIndex )( tDofIndex, 3 ) } )
                    += - mFilterParam * trans( vN->N() ) * trans( aInterfaceNormal ) * vN->dnNdxn( 1 ) * tWStar;
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Interface::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                     moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
//            // set the field interpolator
//            Field_Interpolator* vN = mMasterFI( 0 );
//
//            //FIXME set the interface normal
//            Matrix< DDRMat > aInterfaceNormal( vN->gradx( 1 ).n_cols() , 1, 1.0 );
//
//            // set resiaul size
//            this->set_residual( aResidual );
//
//            // compute the residual
//            aResidual( 0 ) = - mFilterParam * trans( vN->N() ) * trans( vN->gradx( 1 ) ) * aInterfaceNormal;
//
//            // set the jacobian size
//            this->set_jacobian( aJacobians );
//
//            // compute the residual
//            aJacobians( 0 )( 0 ) = - mFilterParam * trans( vN->N() ) * trans( aInterfaceNormal ) * vN->dnNdxn( 1 );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#include "cl_FEM_IWG_L2.hpp"

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
        // set the residual dof type
        mResidualDofType = { MSI::Dof_Type::L2 };

        // set the active dof types
        mMasterDofTypes = { { MSI::Dof_Type::L2 } };

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

    void IWG_L2::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                moris::Cell< Matrix< DDRMat > >                & aResidual )
    {
        // call the right mComputeFunction for residual and jacobian evaluations
        ( this->*mComputeFunction )( aJacobians,
                                     aResidual );
    }

//------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_and_residual_without_alpha( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute Jacobian
            aJacobians( 0 ).resize( 1 );
            aJacobians( 0 )( 0 ) = trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();

            // set the jacobian size
            aResidual.resize( 1 );

            // compute residual
            //FIXME mNodalWeakBCs
            aResidual( 0 ) = aJacobians( 0 )( 0 ) * ( mMasterFI( 0 )->get_coeff() - mNodalWeakBCs );

        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_and_residual_with_alpha( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                               moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute Jacobian
            aJacobians( 0 ).resize( 1 );
            aJacobians( 0 )( 0 ) = trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N()
                            + mAlpha * ( trans( mMasterFI( 0 )->Bx() ) * mMasterFI( 0 )->Bx() );

            // set the residual size
            aResidual.resize( 1 );

            // compute residual
            //FIXME: mNodalWeakBCs
            aResidual( 0 ) = aJacobians( 0 )( 0 ) * ( mMasterFI( 0 )->get_coeff() - mNodalWeakBCs );
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

        void IWG_L2::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // call the right mComputeJacFunction for jacobian evaluation
            ( this->*mComputeJacFunction )( aJacobians );
        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_without_alpha( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute Jacobian
            aJacobians( 0 ).resize( 1 );
            aJacobians( 0 )( 0 ) = trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();
        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_with_alpha( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute Jacobian
            aJacobians( 0 ).resize( 1 );
            aJacobians( 0 )( 0 ) = trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N()
                                 + mAlpha * ( trans( mMasterFI( 0 )->Bx() ) * mMasterFI( 0 )->Bx() );
        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // call the right mComputeResFunction for residual evaluation
            ( this->*mComputeResFunction )( aResidual );
        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_residual_without_alpha( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );

            // set size
            aResidual.resize( 1 );

            // compute residual
            //FIXME: mNodalWeakBCs
            aResidual( 0 ) = trans( mMasterFI( 0 )->N() ) * ( mMasterFI( 0 )->val() - mMasterFI( 0 )->N() * mNodalWeakBCs );

        }

//------------------------------------------------------------------------------

        void IWG_L2::compute_residual_with_alpha( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );

            // set size
            aResidual.resize( 1 );

            // compute residual
            //FIXME mNodalWeakBCs
            aResidual( 0 ) = trans( mMasterFI( 0 )->N() ) * ( mMasterFI( 0 )->val() - mMasterFI( 0 )->N() * mNodalWeakBCs )
                           + mAlpha * trans( mMasterFI( 0 )->Bx() ) * ( mMasterFI( 0 )->gradx( 1 ) - mMasterFI( 0 )->Bx() * mNodalWeakBCs );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

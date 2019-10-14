
#include "cl_FEM_IWG_Test.hpp"
#include "fn_trans.hpp"        //LINALG/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Test::IWG_Test( Field_Interpolator * aFieldInterpolator )
        {
//            MORIS_ERROR( !( aFieldInterpolator->space_only() ),
//                         "IWG_Field - not implemented for space only.");

            mFieldInterpolator = aFieldInterpolator;
        }

//------------------------------------------------------------------------------

        void IWG_Test::compute_residual( Matrix< DDRMat > & aResidual,
                                         uint             & aSpaceDerivativeOrder,
                                         uint             & aTimeDerivativeOrder )
        {
            // r = delta v * u
            Matrix< DDRMat > tN = mFieldInterpolator->N();
            Matrix< DDRMat > tu;
            if      ( aSpaceDerivativeOrder == 0 && aTimeDerivativeOrder == 0 )
            {
                tu = mFieldInterpolator->val();
            }
            else if ( aSpaceDerivativeOrder == 1 && aTimeDerivativeOrder == 0 )
            {
                tu = mFieldInterpolator->gradx( aSpaceDerivativeOrder );
            }
            else if ( aSpaceDerivativeOrder == 2 && aTimeDerivativeOrder == 0 )
            {
                tu = mFieldInterpolator->gradx( aSpaceDerivativeOrder );
            }
            else if ( aSpaceDerivativeOrder == 0 && aTimeDerivativeOrder == 1 )
            {
                tu = mFieldInterpolator->gradt( aTimeDerivativeOrder );
            }
            else if ( aSpaceDerivativeOrder == 0 && aTimeDerivativeOrder == 2 )
            {
                tu = mFieldInterpolator->gradt( aTimeDerivativeOrder );
            }
            else
            {
                MORIS_ERROR( false, "IWG_Field::compute_residual - derivative order not implemented.");
            }

            aResidual = trans( tN ) * tu;
        }

//------------------------------------------------------------------------------

        void IWG_Test::compute_jacobian( Matrix< DDRMat > & aJacobian,
                                         uint             & aSpaceDerivativeOrder,
                                         uint             & aTimeDerivativeOrder )
        {
            // j = N * Nt
            Matrix< DDRMat > tN = mFieldInterpolator->N();

            Matrix< DDRMat > tduduHat;
            if      ( aSpaceDerivativeOrder == 0 && aTimeDerivativeOrder == 0 )
            {
                tduduHat = mFieldInterpolator->N();
            }
            else if ( aSpaceDerivativeOrder == 1 && aTimeDerivativeOrder == 0 )
            {
                tduduHat = mFieldInterpolator->dnNdxn( 1 );
            }
            else if ( aSpaceDerivativeOrder == 2 && aTimeDerivativeOrder == 0 )
            {
                tduduHat = mFieldInterpolator->dnNdxn( 2 );
            }
            else if ( aSpaceDerivativeOrder == 0 && aTimeDerivativeOrder == 1 )
            {
                tduduHat = mFieldInterpolator->dnNdtn( 1 );
            }
            else if ( aSpaceDerivativeOrder == 0 && aTimeDerivativeOrder == 2 )
            {
                tduduHat = mFieldInterpolator->dnNdtn( 2 );
            }
            else
            {
                MORIS_ERROR( false, "IWG_Field::compute_residual - derivative order not implemented.");
            }

            aJacobian = trans( tN ) * tduduHat;
        }

//------------------------------------------------------------------------------

        void IWG_Test::compute_jacobian_and_residual( Matrix< DDRMat > & aJacobian,
                                                      Matrix< DDRMat > & aResidual,
                                                      uint             & aSpaceDerivativeOrder,
                                                      uint             & aTimeDerivativeOrder )
        {
            // r = delta v * u & j = N * Nt
            Matrix< DDRMat > tN = mFieldInterpolator->N();

            Matrix< DDRMat > tu, tduduHat;
            if      ( aSpaceDerivativeOrder == 0 && aTimeDerivativeOrder == 0 )
            {
                tu       = mFieldInterpolator->val();
                tduduHat = mFieldInterpolator->N();
            }
            else if ( aSpaceDerivativeOrder == 1 && aTimeDerivativeOrder == 0 )
            {
                tu       = mFieldInterpolator->gradx( aSpaceDerivativeOrder );
                tduduHat = mFieldInterpolator->dnNdxn( 1 );
            }
            else if ( aSpaceDerivativeOrder == 2 && aTimeDerivativeOrder == 0 )
            {
                tu       = mFieldInterpolator->gradx( aSpaceDerivativeOrder );
                tduduHat = mFieldInterpolator->dnNdxn( 2 );
            }
            else if ( aSpaceDerivativeOrder == 0 && aTimeDerivativeOrder == 1 )
            {
                tu       = mFieldInterpolator->gradt( aTimeDerivativeOrder );
                tduduHat = mFieldInterpolator->dnNdtn( 1 );
            }
            else if ( aSpaceDerivativeOrder == 0 && aTimeDerivativeOrder == 2 )
            {
                tu       = mFieldInterpolator->gradt( aTimeDerivativeOrder );
                tduduHat = mFieldInterpolator->dnNdtn( 2 );
            }
            else
            {
                MORIS_ERROR( false, "IWG_Field::compute_residual - derivative order not implemented.");
            }

            aResidual = trans( tN ) * tu;
            aJacobian = trans( tN ) * tduduHat;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

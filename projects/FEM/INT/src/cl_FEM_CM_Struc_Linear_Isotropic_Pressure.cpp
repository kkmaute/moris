
#include "cl_FEM_CM_Struc_Linear_Isotropic_Pressure.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic_Pressure::eval_strain()
        {
            // compute displacement gradient
            Matrix<DDRMat> tN = mDofFI(0)->N();

            switch ( mSpaceDim )
            {
                case ( 2 ):
                {
                    mStrain.set_size(3, tN.numel(), 0.0);
                    for (uint tCol = 0; tCol < tN.numel(); tCol++)
                    {
                        mStrain(0, tCol) = tN(tCol);
                        mStrain(1, tCol) = tN(tCol);
                    }
                    break;
                }
                case( 3 ):
                {
                    mStrain.set_size(6, tN.numel(), 0.0);
                    for (uint tCol = 0; tCol < tN.numel(); tCol++)
                    {
                        mStrain(0, tCol) = tN(tCol);
                        mStrain(1, tCol) = tN(tCol);
                        mStrain(2, tCol) = tN(tCol);
                    }
                    break;
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Linear_Isotropic_Pressure::eval_strain - Flattening of strain tensor only implemented in 2 and 3D");
                }
            }
        }

        void CM_Struc_Linear_Isotropic_Pressure::eval_flux()
        {
            // compute displacement gradient
            moris::real tP = mDofFI(0)->val()(0);

            switch ( mSpaceDim )
            {
                case ( 2 ):
                {
                    mFlux.set_size(3, 1, 0.0);
                    mFlux(0) = tP;
                    mFlux(1) = tP;
                    break;
                }
                case( 3 ):
                {
                    mFlux.set_size(6, 1, 0.0);
                    mFlux(0) = tP;
                    mFlux(1) = tP;
                    mFlux(2) = tP;
                    break;
                }
                default:
                {
                    MORIS_ERROR(false, "CM_Struc_Linear_Isotropic_Pressure::eval_strain - Flattening of strain tensor only implemented in 2 and 3D");
                }
            }
        }
        
        //--------------------------------------------------------------------------------------------------------------

//        void CM_Struc_Linear_Isotropic_Pressure::eval_const()
//        {
//            uint tNuIndex = static_cast< uint >( Property_Type::NU );
//            moris::real tNu = mProperties(tNuIndex)->val()(0);
//            uint tEmodIndex = static_cast< uint >( Property_Type::EMOD );
//            moris::real tEmod = mProperties(tEmodIndex)->val()(0);
//            (this->*mConstFunc)(tEmod, tNu);
//        }

        //--------------------------------------------------------------------------------------------------------------

        void CM_Struc_Linear_Isotropic_Pressure::eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
//            // get the dof type as a uint
//            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );
//
//            // get the dof type index
//            uint tDofIndex = mGlobalDofTypeMap( tDofType );
//
//            // if direct dependency on the dof type
//            if( tDofType < mDofTypeMap.numel() && mDofTypeMap( tDofType ) != -1 && aDofTypes(0) == MSI::Dof_Type::P )
//            {
//                // compute derivative with direct dependency
//                mdStraindDof( tDofIndex ) = mDofFI(0)->N();
//            }
//            else
//            {
//                // reset the matrix
//                mdStraindDof( tDofIndex ).set_size( 1, mDofFI( tDofIndex )->get_number_of_space_time_coefficients(), 0.0 );
//            }
        }

        //--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

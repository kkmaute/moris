#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Spalart_Allmaras_Coefficient.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Spalart_Allmaras_Coefficient::IQI_Spalart_Allmaras_Coefficient()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::SPALART_ALLMARAS_COEFFICIENT;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "SpalartAllmarasTurbulence" ] =
                    static_cast< uint >( IQI_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
        }
        
        //------------------------------------------------------------------------------

        void IQI_Spalart_Allmaras_Coefficient::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMSATurbulence =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // switch on type index
            switch ( mIQITypeIndex )
            {
                // production coefficient
                case 0:
                {
                    // compute production coefficient
                    aQI = tCMSATurbulence->production_coefficient();

                    break;
                }
                // wall destruction coefficient
                case 1:
                {
                    aQI = tCMSATurbulence->wall_destruction_coefficient();
                    break;
                }
                // diffusion coefficient
                case 2:
                {
                    aQI = tCMSATurbulence->diffusion_coefficient();
                    break;
                }
                // if none of the above
                default:
                {
                    MORIS_ERROR(false, "IQI_Spalart_Allmaras_Coefficient::compute_QI - wrong mIQITypeIndex type");
                }
            }
        }
        
        //------------------------------------------------------------------------------

        void IQI_Spalart_Allmaras_Coefficient::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // init QI value matrix
            Matrix<DDRMat> tQI(1,1);
            this->compute_QI( tQI );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * tQI(0);
        }
        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */




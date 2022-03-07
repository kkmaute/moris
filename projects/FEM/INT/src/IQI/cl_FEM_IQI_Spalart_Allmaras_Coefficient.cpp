#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Spalart_Allmaras_Coefficient.hpp"

#include "cl_FEM_CM_Spalart_Allmaras_Turbulence.hpp"

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

        void
        IQI_Spalart_Allmaras_Coefficient::compute_QI( Matrix< DDRMat >& aQI )
        {
            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // cast constitutive model base class pointer to SA constitutive model
            CM_Spalart_Allmaras_Turbulence* tCMSATurbulencePtr =
                    dynamic_cast< CM_Spalart_Allmaras_Turbulence* >( tCMSATurbulence.get() );

            // switch on type index
            switch ( mIQITypeIndex )
            {
                // production coefficient
                case 0:
                {
                    // compute production coefficient
                    aQI = tCMSATurbulencePtr->production_coefficient();

                    break;
                }
                // wall destruction coefficient
                case 1:
                {
                    aQI = tCMSATurbulencePtr->wall_destruction_coefficient();
                    break;
                }
                // diffusion coefficient
                case 2:
                {
                    aQI = tCMSATurbulencePtr->diffusion_coefficient();
                    break;
                }
                // production term
                case 3:
                {
                    aQI = tCMSATurbulencePtr->production_term();
                    break;
                }
                // wall destruction term
                case 4:
                {
                    aQI = tCMSATurbulencePtr->wall_destruction_term();
                    break;
                }
                // s
                case 5:
                {
                    aQI = tCMSATurbulencePtr->s();
                    break;
                }
                // stilde
                case 6:
                {
                    aQI = tCMSATurbulencePtr->stilde();
                    break;
                }
                // smod
                case 7:
                {
                    aQI = tCMSATurbulencePtr->smod();
                    break;
                }
                // chi
                case 8:
                {
                    aQI = tCMSATurbulencePtr->chi();
                    break;
                }
                    // ft2
                case 9:
                {
                    aQI = tCMSATurbulencePtr->ft2();
                    break;
                }
                // if none of the above
                default:
                {
                    MORIS_ERROR( false, "IQI_Spalart_Allmaras_Coefficient::compute_QI - wrong mIQITypeIndex type" );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI_Spalart_Allmaras_Coefficient::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // init QI value matrix
            Matrix< DDRMat > tQI( 1, 1 );
            this->compute_QI( tQI );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * tQI( 0 );
        }
        //------------------------------------------------------------------------------
    }    // namespace fem
}    // namespace moris

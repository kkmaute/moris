
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Stabilization.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Stabilization::IQI_Stabilization()
        {
            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IQI_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "Stabilization" ] = static_cast< uint >( IQI_Stabilization_Type::STABILIZATION );
        }

        //------------------------------------------------------------------------------

        void IQI_Stabilization::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP =
                    mStabilizationParam( static_cast< uint >( IQI_Stabilization_Type::STABILIZATION ) );

            // check if index was set
            if( mQuantityDofType.size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Stabilization::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // evaluate the QI
            aQI = {{ tSP->val()( mIQITypeIndex ) }};
        }

        //------------------------------------------------------------------------------

        void IQI_Stabilization::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP =
                    mStabilizationParam( static_cast< uint >( IQI_Stabilization_Type::STABILIZATION ) );

            // check if index was set
            if( mQuantityDofType.size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Stabilization::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * ( tSP->val()( mIQITypeIndex ) );
        }

        //------------------------------------------------------------------------------

        void IQI_Stabilization::compute_dQIdu( real aWStar )
        {
            MORIS_ERROR(false, "IQI_Stabilization::compute_dQIdu - not implemented\n.");
        }

        //------------------------------------------------------------------------------

        void IQI_Stabilization::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            MORIS_ERROR(false, "IQI_Stabilization::compute_dQIdu - not implemented\n.");
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */




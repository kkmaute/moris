#include "cl_FEM_Cluster.hpp"                       //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"    //FEM/INT/src
#include "cl_FEM_SP_Reciprocal_Total_Volume.hpp"    //FEM/INT/src

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "op_div.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        moris::Cell< std::tuple<
                fem::Measure_Type,
                mtk::Primary_Void,
                mtk::Master_Slave > >
        SP_Reciprocal_Total_Volume::get_cluster_measure_tuple_list()
        {
            return { mMasterVolumeTuple, mSlaveVolumeTuple };
        }

        //------------------------------------------------------------------------------

        void
        SP_Reciprocal_Total_Volume::eval_SP()
        {
            // get master volume cluster measure value
            real tMasterVolume = mCluster->get_cluster_measure(
                                                 std::get< 0 >( mMasterVolumeTuple ),
                                                 std::get< 1 >( mMasterVolumeTuple ),
                                                 std::get< 2 >( mMasterVolumeTuple ) )
                                         ->val()( 0 );

            // get slave volume cluster measure value
            real tSlaveVolume = mCluster->get_cluster_measure(
                                                std::get< 0 >( mSlaveVolumeTuple ),
                                                std::get< 1 >( mSlaveVolumeTuple ),
                                                std::get< 2 >( mSlaveVolumeTuple ) )
                                        ->val()( 0 );

            // compute stabilization parameter value
            mPPVal = { { 1.0 / ( tMasterVolume + tSlaveVolume ) } };
        }

        //------------------------------------------------------------------------------

        void
        SP_Reciprocal_Total_Volume::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            MORIS_ERROR( false, "SP_Volume_Fraction::eval_dSPdMasterDOF(), not implemented for this SP" );
        }

        //------------------------------------------------------------------------------

        void
        SP_Reciprocal_Total_Volume::eval_dSPdSlaveDOF( const moris::Cell< MSI::Dof_Type >& aDofTypes )
        {
            MORIS_ERROR( false, "SP_Volume_Fraction::eval_dSPdSlaveDOF(), not implemented for this SP" );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

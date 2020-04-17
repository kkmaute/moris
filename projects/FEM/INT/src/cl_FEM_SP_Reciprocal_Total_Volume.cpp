#include "cl_FEM_Cluster.hpp"              //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"              //FEM/INT/src
#include "cl_FEM_SP_Reciprocal_Total_Volume.hpp" //FEM/INT/src

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
        SP_Reciprocal_Total_Volume::SP_Reciprocal_Total_Volume()
        {
            // set size for the property pointer cells
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );
            mSlaveProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );
        }

//------------------------------------------------------------------------------
        void SP_Reciprocal_Total_Volume::reset_cluster_measures()
        {
            // evaluate cluster measures from the cluster
            mMasterVolume     = mCluster->compute_cluster_cell_measure( mtk::Primary_Void::PRIMARY, mtk::Master_Slave::MASTER );
            mSlaveVolume      = mCluster->compute_cluster_cell_measure( mtk::Primary_Void::VOID, mtk::Master_Slave::SLAVE );
        }

//------------------------------------------------------------------------------
        void SP_Reciprocal_Total_Volume::eval_SP()
        {
            mPPVal.set_size( 1, 1 );

            // compute stabilization parameter value
            mPPVal( 0 ) = 1 / ( mMasterVolume + mSlaveVolume );
        }

//------------------------------------------------------------------------------
        void SP_Reciprocal_Total_Volume::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            MORIS_ERROR( false, "SP_Volume_Fraction::eval_dSPdMasterDOF(), not implemented for this SP" );
        }

//------------------------------------------------------------------------------
        void SP_Reciprocal_Total_Volume::eval_dSPdSlaveDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            MORIS_ERROR( false, "SP_Volume_Fraction::eval_dSPdSlaveDOF(), not implemented for this SP" );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


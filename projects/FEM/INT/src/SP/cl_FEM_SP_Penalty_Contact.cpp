/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Penalty_Contact.cpp
 *
 */

#include "cl_FEM_SP_Penalty_Contact.hpp"   //FEM/INT/src
#include "cl_FEM_Cluster.hpp"              //FEM/INT/src

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

        SP_Penalty_Contact::SP_Penalty_Contact()
        {
            // set size for the property pointer cells
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );
            mSlaveProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Material" ] = static_cast< uint >( SP_Property_Type::MATERIAL );
        }

        //------------------------------------------------------------------------------

        void SP_Penalty_Contact::reset_cluster_measures()
        {
            // evaluate cluster measures from the cluster
            mMasterVolume = mCluster->compute_cluster_cell_measure(
                    mtk::Primary_Void::INTERP,
                    mtk::Master_Slave::MASTER );
            mSlaveVolume  = mCluster->compute_cluster_cell_measure(
                    mtk::Primary_Void::INTERP,
                    mtk::Master_Slave::SLAVE );
        }

        //------------------------------------------------------------------------------

        void SP_Penalty_Contact::eval_SP()
        {
            moris::real tEMaster = mMasterProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 );
            moris::real tESlave = mSlaveProp( static_cast< uint >( SP_Property_Type::MATERIAL ) )->val()( 0 );

            mPPVal = mParameters( 0 ) * (mMasterVolume /( mMasterVolume / tEMaster + mSlaveVolume / tESlave)
                    + mSlaveVolume /( mMasterVolume / tEMaster + mSlaveVolume / tESlave)) / mElementSize;
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


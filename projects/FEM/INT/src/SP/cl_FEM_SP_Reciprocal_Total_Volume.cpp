/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Reciprocal_Total_Volume.cpp
 *
 */

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

        moris::Vector< std::tuple<
                fem::Measure_Type,
                mtk::Primary_Void,
                mtk::Leader_Follower > >
        SP_Reciprocal_Total_Volume::get_cluster_measure_tuple_list()
        {
            return { mLeaderVolumeTuple, mFollowerVolumeTuple };
        }

        //------------------------------------------------------------------------------

        void
        SP_Reciprocal_Total_Volume::eval_SP()
        {
            // get leader volume cluster measure value
            real tLeaderVolume = mCluster->get_cluster_measure(
                                                 std::get< 0 >( mLeaderVolumeTuple ),
                                                 std::get< 1 >( mLeaderVolumeTuple ),
                                                 std::get< 2 >( mLeaderVolumeTuple ) )
                                         ->val()( 0 );

            // get follower volume cluster measure value
            real tFollowerVolume = mCluster->get_cluster_measure(
                                                std::get< 0 >( mFollowerVolumeTuple ),
                                                std::get< 1 >( mFollowerVolumeTuple ),
                                                std::get< 2 >( mFollowerVolumeTuple ) )
                                        ->val()( 0 );

            // compute stabilization parameter value
            mPPVal = { { 1.0 / ( tLeaderVolume + tFollowerVolume ) } };
        }

        //------------------------------------------------------------------------------

        void
        SP_Reciprocal_Total_Volume::eval_dSPdLeaderDOF( const moris::Vector< MSI::Dof_Type >& aDofTypes )
        {
            MORIS_ERROR( false, "SP_Volume_Fraction::eval_dSPdLeaderDOF(), not implemented for this SP" );
        }

        //------------------------------------------------------------------------------

        void
        SP_Reciprocal_Total_Volume::eval_dSPdFollowerDOF( const moris::Vector< MSI::Dof_Type >& aDofTypes )
        {
            MORIS_ERROR( false, "SP_Volume_Fraction::eval_dSPdFollowerDOF(), not implemented for this SP" );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


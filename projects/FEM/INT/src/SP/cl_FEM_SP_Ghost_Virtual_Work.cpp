/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Ghost_Virtual_Work.cpp
 *
 */

#include "cl_FEM_SP_Ghost_Virtual_Work.hpp" //FEM/INT/src
#include "cl_FEM_Cluster.hpp"              //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"              //FEM/INT/src

namespace moris::fem
{
    //------------------------------------------------------------------------------

    Vector< std::tuple<
            fem::Measure_Type,
            mtk::Primary_Void,
            mtk::Leader_Follower > >
    SP_Ghost_Virtual_Work::get_cluster_measure_tuple_list()
    {
        return { mElementSizeTuple };
    }

    //------------------------------------------------------------------------------

    void SP_Ghost_Virtual_Work::eval_SP()
    {
        // get element size cluster measure value
        real tElementSize = mCluster->get_cluster_measure(
                                            std::get< 0 >( mElementSizeTuple ),
                                            std::get< 1 >( mElementSizeTuple ),
                                            std::get< 2 >( mElementSizeTuple ) )
                                    ->val()( 0 );

        // compute stabilization parameter value
        mPPVal = mParameters( 0 ) * std::pow( tElementSize, 2 * ( mOrder - 1 ) + 1 );
    }

    //------------------------------------------------------------------------------

    void SP_Ghost_Virtual_Work::eval_dSPdLeaderDOF(
            const Vector< MSI::Dof_Type > &aDofTypes )
    {
        // get the dof type as a uint
        uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        uint tDofIndex = mLeaderGlobalDofTypeMap( tDofType );

        // get FI for derivative dof type
        Field_Interpolator *tFIDer =
                mLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // reset the matrix
        mdPPdLeaderDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem

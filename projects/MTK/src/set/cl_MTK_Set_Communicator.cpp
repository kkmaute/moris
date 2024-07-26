/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Set_Communicator.cpp
 *
 */

#include "cl_MTK_Set_Communicator.hpp"

#include "cl_MTK_Set.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        Set_Communicator::Set_Communicator( Vector< mtk::Set* >& aSetsToCommunicate )
        {
            // get the number of sets to be communicated
            uint tNumSetsToComm = aSetsToCommunicate.size();

            /* ---------------------------------------- */
            // 1.) make sure that each proc has the same sets (same quantity and same names)
            MORIS_ASSERT(
                    this->check_sets( aSetsToCommunicate ),
                    "mtk::Set_Communicator::Set_Communicator() - "
                    "List of sets to communicate not consistent across all procs." );

            /* ---------------------------------------- */
            // 2.) collect the important set information on the current proc
            Matrix< DDUMat > tIpCellGeometries( tNumSetsToComm, 1, (uint)mtk::Geometry_Type::UNDEFINED );
            Matrix< DDUMat > tIgCellGeometries( tNumSetsToComm, 1, (uint)mtk::Geometry_Type::UNDEFINED );
            Matrix< DDUMat > tIpInterpolationOrders( tNumSetsToComm, 1, (uint)mtk::Interpolation_Order::UNDEFINED );
            Matrix< DDUMat > tIgInterpolationOrders( tNumSetsToComm, 1, (uint)mtk::Interpolation_Order::UNDEFINED );

            // collect the information from each set and store it in the matrix
            for ( uint iSet = 0; iSet < tNumSetsToComm; iSet++ )
            {
                // get access to the current set
                mtk::Set* tSet = aSetsToCommunicate( iSet );

                // skip and communicate the default if the Set is a nullptr
                if ( !tSet )
                {
                    continue;
                }

                // copy into matrices
                tIpCellGeometries( iSet )      = (uint)tSet->get_interpolation_cell_geometry_type();
                tIgCellGeometries( iSet )      = (uint)tSet->get_integration_cell_geometry_type();
                tIpInterpolationOrders( iSet ) = (uint)tSet->get_interpolation_cell_interpolation_order();
                tIgInterpolationOrders( iSet ) = (uint)tSet->get_integration_cell_interpolation_order();
            }

            /* ---------------------------------------- */
            // 3.) communicate all the information
            min_all_matrix( tIpCellGeometries );
            min_all_matrix( tIgCellGeometries );
            min_all_matrix( tIpInterpolationOrders );
            min_all_matrix( tIgInterpolationOrders );

            /* ---------------------------------------- */
            // 4.) store communicated information on current proc

            // collect the information from each set and store it in the matrix
            for ( uint iSet = 0; iSet < tNumSetsToComm; iSet++ )
            {
                // get access to the current set
                mtk::Set* tSet = aSetsToCommunicate( iSet );

                // skip if the Set is a nullptr
                if ( !tSet )
                {
                    continue;
                }

                // translate the information obtained from communication
                mtk::Geometry_Type       tIpCellGeomType    = static_cast< mtk::Geometry_Type >( tIpCellGeometries( iSet ) );
                mtk::Geometry_Type       tIgCellGeomType    = static_cast< mtk::Geometry_Type >( tIgCellGeometries( iSet ) );
                mtk::Interpolation_Order tIpCellInterpOrder = static_cast< mtk::Interpolation_Order >( tIpInterpolationOrders( iSet ) );
                mtk::Interpolation_Order tIgCellInterpOrder = static_cast< mtk::Interpolation_Order >( tIgInterpolationOrders( iSet ) );

                // store away information
                tSet->set_interpolation_cell_geometry_type( tIpCellGeomType );
                tSet->set_integration_cell_geometry_type( tIgCellGeomType );
                tSet->set_interpolation_cell_interpolation_order( tIpCellInterpOrder );
                tSet->set_integration_cell_interpolation_order( tIgCellInterpOrder );
            }

        }    // end function: Set_Communicator::Set_Communicator()

        //------------------------------------------------------------------------------

        bool
        Set_Communicator::check_sets( Vector< mtk::Set* >& aSetsToCommunicate )
        {
            // get the number of sets to be communicated
            uint tNumSetsToComm = aSetsToCommunicate.size();

            // check that all procs have the same number of sets that they are trying to communicate
            uint tMinNumberOfSetsOnProcs = min_all( tNumSetsToComm );
            MORIS_ERROR(
                    tMinNumberOfSetsOnProcs == tNumSetsToComm,
                    "mtk::Set_Communicator::check_sets() - "
                    "The number of sets in the list for set communication is different across processors." );

            // TODO: check that all sets have the same name

            // return a pass if no assertion has been thrown up to here
            return true;
        }

        //------------------------------------------------------------------------------

    }    // namespace mtk

}    // namespace moris

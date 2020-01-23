/*
 * cl_MTK_Mesh_Checker.cpp
 *
 *  Created on: Jan 15, 2020
 *      Author: doble
 */

#include "cl_MTK_Mesh_Checker.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "fn_norm.hpp"
#include "op_minus.hpp"
namespace moris
{

namespace mtk
{
//--------------------------------------------------------------------------------
Mesh_Checker::Mesh_Checker(){}
//--------------------------------------------------------------------------------
Mesh_Checker::~Mesh_Checker(){}
//--------------------------------------------------------------------------------
bool
Mesh_Checker::verify_double_side_sets(Integration_Mesh const * aIgMesh)
{
    moris::uint tNumDoubleSideSets = aIgMesh->get_num_double_side_set();

    // iterate through the double side sets
    for(moris::uint iSet = 0 ; iSet < tNumDoubleSideSets; iSet++)
    {
        moris::Cell<Cluster const*> tClustersInSet = aIgMesh->get_double_side_set_cluster((moris_index)iSet);

        // iterate through clusters in set
        for(auto iCl:tClustersInSet)
        {
            // verify the master
            bool tMasterValid = this->verify_side_cluster(iCl, Master_Slave::MASTER);

            if(!tMasterValid)
            {
                MORIS_LOG_ERROR("\nInvalid master cluster");
            }

            // verify the slave
            bool tSlaveValid = this->verify_side_cluster(iCl, Master_Slave::SLAVE);

            if(!tSlaveValid)
            {
                MORIS_LOG_ERROR("\nInvalid master cluster");
            }

            if(!tMasterValid || !tSlaveValid)
            {
                return false;
            }
        }
    }

    return true;
}
//--------------------------------------------------------------------------------
bool
Mesh_Checker::verify_side_cluster(Cluster const* aCluster,
                                  enum Master_Slave aMasterSlave)
{
    // check if trivial
    bool tTrivial = aCluster->is_trivial(aMasterSlave);

    if(tTrivial)
    {
        // get the interpolation cell
        mtk::Cell const & tIpCell = aCluster->get_interpolation_cell(aMasterSlave);

        // get the integration primary cells
        moris::Cell<moris::mtk::Cell const *> const & tIgCells = aCluster->get_primary_cells_in_cluster( aMasterSlave );

        if(tIgCells.size() != 1)
        {
            MORIS_LOG_ERROR("\nTrivial cluster needs to have 1 primary integration cell.");
            return false;
        }

        // get the side ordinals
        moris::Matrix<moris::IndexMat> tSideOrd = aCluster->get_cell_side_ordinals(aMasterSlave);

        if(tSideOrd.numel() != 1)
        {
            MORIS_LOG_ERROR("\nTrivial cluster needs to have exactly 1 side ordinal.");
            return false;
        }

        moris::Cell<moris::mtk::Vertex const *> tIpVertsOnSide = tIpCell.get_vertices_on_side_ordinal(tSideOrd(0));
        moris::Cell<moris::mtk::Vertex const *> tIgVertsOnSide = tIgCells(0)->get_vertices_on_side_ordinal(tSideOrd(0));

        if(tIpVertsOnSide.size() != tIgVertsOnSide.size())
        {
            MORIS_LOG_ERROR("\nTrivial cluster ip and ig vertex on side mismatch.");
            return false;
        }


        // check the coordinates
        for(moris::uint iV = 0; iV < tIpVertsOnSide.size(); iV++)
        {
            if(moris::norm(tIpVertsOnSide(iV)->get_coords() - tIgVertsOnSide(iV)->get_coords()) > 1e-8)
            {
                return false;
            }
        }

    }
    else
    {

        // get the interpolation cell
        mtk::Cell const & tIpCell = aCluster->get_interpolation_cell(aMasterSlave);

        // create cell info
        Cell_Info_Factory tCellInfoFactory;
        moris::mtk::Cell_Info* tCellInfo = tCellInfoFactory.create_cell_info(tIpCell.get_geometry_type(),tIpCell.get_interpolation_order());

        // get the vertices in the cluster
        moris::Cell<moris::mtk::Vertex const *> tVertsInCluster = aCluster->get_vertices_in_cluster(aMasterSlave);

        // ip verts
        Matrix<DDRMat> tIpCoords = tIpCell.get_vertex_coords();

        for(size_t i= 0; i<tVertsInCluster.size(); i++)
        {
            // local coordinates
            moris::Matrix<moris::DDRMat> tLocalCoords = aCluster->get_vertex_local_coordinate_wrt_interp_cell( tVertsInCluster(i), aMasterSlave);

            // Evalute the basis function at the point
            moris::Matrix<moris::DDRMat> tN;
            tCellInfo->eval_N(tLocalCoords,tN);

            // Evaluate the nodes global coordinate from the basis weights
            moris::Matrix<moris::DDRMat> tInterpedNodeCoord = tN*tIpCoords;

            // Verify the interpolated coordinate is equal to the node coordinate row
            if(norm(tInterpedNodeCoord - tVertsInCluster(i)->get_coords()) > 1e-8)
            {
                return false;
            }
        }

        delete tCellInfo;
    }

    return true;
}
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

}
}

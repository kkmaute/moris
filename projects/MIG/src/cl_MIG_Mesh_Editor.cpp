/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MIG_Mesh_Editor.cpp
 *
 */

#include "cl_Matrix.hpp"
#include <numeric>

#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_MIG_Mesh_Editor.hpp"
#include "cl_MIG_Periodic_2D.hpp"
#include "cl_MIG_Periodic_3D.hpp"
#include "cl_MTK_Mesh_DataBase_IG.hpp"
#include "cl_Tracer.hpp"

namespace moris::mig
{
    //------------------------------------------------------------------------------------------------------------
    Periodic_Mesh_Editor::Periodic_Mesh_Editor( mtk::Integration_Mesh_DataBase_IG* aIGMesh, mig::Periodic_2D* aPeriodicData2D ) :
        mPeriodicData2D( aPeriodicData2D )
    {
        mOutputMesh = aIGMesh;
    }

    Periodic_Mesh_Editor::Periodic_Mesh_Editor( mtk::Integration_Mesh_DataBase_IG* aIGMesh, mig::Periodic_3D* aPeriodicData3D ) :
        mPeriodicData3D( aPeriodicData3D )
    {
        mOutputMesh = aIGMesh;
    }
    //------------------------------------------------------------------------------------------------------------

    Periodic_Mesh_Editor::~Periodic_Mesh_Editor()
    {
    }

    //------------------------------------------------------------------------------------------------------------

    void
    Periodic_Mesh_Editor::perform()
    {
        Tracer tTracer("MIG", "Editor", "Merging Meshes");

        // link the newly created nodes to the geomtry engine
        this->link_nodes_to_geomtry_engine();

        // depending on the dimension add data to the data base
        if ( mPeriodicData3D == nullptr )
        {
            this->construct_periodic_data_base(
                mPeriodicData2D->mSideClusterToVertexIndices,
                mPeriodicData2D->mVerticesCoords,
                mPeriodicData2D->mSideClusterToCells,
                mPeriodicData2D->mCellToVertexIndices,
                mPeriodicData2D->mSideClusterToIPCell,
                mPeriodicData2D->mVertexParametricCoords,
                mPeriodicData2D->mDoubleSidedClustersIndex,
                mPeriodicData2D->mNumDblSideCluster,
                mGeometryEngine->get_num_phases() );
        }

        else
        {
            this->construct_periodic_data_base(
                mPeriodicData3D->mSideClusterToVertexIndices,
                mPeriodicData3D->mVerticesCoords,
                mPeriodicData3D->mSideClusterToCells,
                mPeriodicData3D->mCellToVertexIndices,
                mPeriodicData3D->mSideClusterToIPCell,
                mPeriodicData3D->mVertexParametricCoords,
                mPeriodicData3D->mDoubleSidedClustersIndex,
                mPeriodicData3D->mNumDblSideCluster,
                mGeometryEngine->get_num_phases() );
        }
    }

    //------------------------------------------------------------------------------------------------------------
    void
    Periodic_Mesh_Editor::set_geometry_engine( moris::gen::Geometry_Engine* aGeometryEngine )
    {
        mGeometryEngine = aGeometryEngine;
    }

    //------------------------------------------------------------------------------------------------------------

    void
    Periodic_Mesh_Editor::link_nodes_to_geomtry_engine()
    {
        // get number of current nodes
        uint tNumPreviousVertices = mOutputMesh->get_num_nodes();

        // update underlying ids and owners of interpolation nodes in GE
        if ( mOutputMesh->get_spatial_dim() == 2 )
        {
            for ( uint iVertex = 0; iVertex < mPeriodicData2D->mNumVertices; iVertex++ )
            {
                moris_index tNodeIndex = iVertex + tNumPreviousVertices;

                // FIXME - it is for serial
                moris_id    tNodeId    = tNodeIndex + 1;
                moris_index tNodeOwner = 0;

                // add nodes to the ge
                mGeometryEngine->update_intersection_node( tNodeIndex, tNodeId, tNodeOwner );
            }
        }
        else
        {
            for ( uint iVertex = 0; iVertex < mPeriodicData3D->mNumVertices; iVertex++ )
            {
                moris_index tNodeIndex = iVertex + tNumPreviousVertices;

                // FIXME - it is for serial
                moris_id    tNodeId    = tNodeIndex + 1;
                moris_index tNodeOwner = 0;

                // add nodes to the ge
                mGeometryEngine->update_intersection_node( tNodeIndex, tNodeId, tNodeOwner );
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------

    void
    Periodic_Mesh_Editor::construct_periodic_data_base(
        moris::Cell< moris::Cell< moris_index > >& aSideClusterToVertexIndices,
        Matrix< DDRMat >                           aVerticesCoords,
        moris::Cell< moris::Cell< moris_index > >& aSideClusterToCells,
        moris::Cell< moris::Cell< moris_index > >& aCellToVertexIndices,
        moris::Cell< moris_index >&                aSideClusterToIPCell,
        Matrix< DDRMat >&                          aVertexParametricCoords,
        moris::Cell< moris_index >&                aDoubleSidedClustersIndex,
        uint                                       mNumDblSideCluster,
        uint                                       aNumGeometry )
    {
        // call the parent class function since the data is in MTK
        Integration_Mesh_Editor::construct_periodic_data_base(
            aSideClusterToVertexIndices,
            aVerticesCoords,
            aSideClusterToCells,
            aCellToVertexIndices,
            aSideClusterToIPCell,
            aVertexParametricCoords,
            aDoubleSidedClustersIndex,
            mNumDblSideCluster,
            aNumGeometry );
    }
}// namespace moris::mig


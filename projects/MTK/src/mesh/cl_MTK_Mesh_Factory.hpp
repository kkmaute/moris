/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_Factory.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MESH_FACTORY_HPP_
#define PROJECTS_MTK_SRC_CL_MESH_FACTORY_HPP_

#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Cell_Cluster_Input.hpp"
#include "cl_MTK_Side_Cluster_Input.hpp"

namespace moris
{
    namespace mtk
    {
        /**
         * Create an interpolation from an Exodus file.
         *
         * @param aMeshType Mesh type to create
         * @param aFileName Exodus file name
         * @param aSuppMeshData Supplementary mesh data to add to the mesh
         * @param aCreateFacesAndEdges Whether or not to create faces and edges
         * @return Created interpolation mesh
         */
        Interpolation_Mesh* create_interpolation_mesh(
                enum MeshType aMeshType,
                std::string   aFileName,
                MtkMeshData*  aSuppMeshData        = nullptr,
                const bool    aCreateFacesAndEdges = true );

        /**
         * Create an interpolation mesh from mesh data.
         *
         * @param aMeshType Mesh type to create
         * @param aMeshData Mesh data for creating the mesh
         * @return Created interpolation mesh
         */
        Interpolation_Mesh* create_interpolation_mesh(
                enum MeshType aMeshType,
                MtkMeshData   aMeshData );

        /**
         * Create an integration mesh from an Exodus file.
         *
         * @param aMeshType Mesh type to create
         * @param aFileName Exodus file name
         * @param aSuppMeshData Supplementary mesh data to add to the mesh
         * @param aCreateFacesAndEdges Whether or not to create faces and edges
         * @return Created integration mesh
         */
        Integration_Mesh* create_integration_mesh(
                enum MeshType aMeshType,
                std::string   aFileName,
                MtkMeshData*  aSuppMeshData        = nullptr,
                const bool    aCreateFacesAndEdges = true );

        /**
         * Create an integration mesh from mesh data.
         *
         * @param aMeshType Mesh type to create
         * @param aMeshData Mesh data for creating the mesh
         * @return Created integration mesh
         */
        Integration_Mesh* create_integration_mesh(
                enum MeshType aMeshType,
                MtkMeshData   aMeshData );

        /**
         * Create an integration mesh from mesh data with a link to an interpolation mesh.
         *
         * @param aMeshType Mesh type to create
         * @param aMeshData Mesh data for creating the mesh
         * @param aInterpMesh Interpolation mesh to link
         * @return Created integration mesh
         */
        Integration_Mesh* create_integration_mesh(
                enum MeshType       aMeshType,
                MtkMeshData         aMeshData,
                Interpolation_Mesh* aInterpMesh );

        /**
         * Create an integration mesh from an existing interpolation mesh.
         *
         * @param aMeshType Mesh type to create
         * @param aInterpMesh Interpolation mesh to get mesh data from
         * @param aCellClusterData Additional cell cluster data
         * @return Created integration mesh
         */
        Integration_Mesh* create_integration_mesh_from_interpolation_mesh(
                enum MeshType       aMeshType,
                Interpolation_Mesh* aInterpMesh,
                uint                aMeshIndex       = 0,
                Cell_Cluster_Input* aCellClusterData = nullptr );

    }    // namespace mtk
}    // namespace moris

#endif /* PROJECTS_MTK_SRC_CL_MESH_FACTORY_HPP_ */

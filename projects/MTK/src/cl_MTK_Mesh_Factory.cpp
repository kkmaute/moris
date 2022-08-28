/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_Factory.cpp
 *
 */

#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"

namespace moris
{
    namespace mtk
    {

        //--------------------------------------------------------------------------------------------------------------

        Interpolation_Mesh* create_interpolation_mesh(
                enum MeshType aMeshType,
                std::string aFileName,
                MtkMeshData *aSuppMeshData,
                const bool aCreateFacesAndEdges)
        {
            // Mark mesh data as being supplementary to an input file
            if(aSuppMeshData != nullptr)
            {
                aSuppMeshData->SupplementaryToFile = true;
            }

            // Create interpolation mesh
            switch (aMeshType)
            {
                case MeshType::STK:
                {
                    return new Interpolation_Mesh_STK( aFileName, aSuppMeshData, aCreateFacesAndEdges );
                }
                default:
                {
                    MORIS_ERROR(false, "Cannot create specified integration mesh type from file.");
                    return nullptr;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Interpolation_Mesh* create_interpolation_mesh(
                enum MeshType aMeshType,
                MtkMeshData aMeshData)
        {
            switch (aMeshType)
            {
                case MeshType::STK:
                {
                    return new Interpolation_Mesh_STK( aMeshData );
                }
                default:
                {
                    MORIS_ERROR(false, "Cannot create specified interpolation mesh type with mesh data.");
                    return nullptr;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Integration_Mesh* create_integration_mesh(
                enum MeshType aMeshType,
                std::string aFileName,
                MtkMeshData *aSuppMeshData,
                const bool aCreateFacesAndEdges)
        {
            // Mark mesh data as being supplementary to an input file
            if (aSuppMeshData != nullptr)
            {
                aSuppMeshData->SupplementaryToFile = true;
            }

            // Create integration mesh
            switch (aMeshType)
            {
                case MeshType::STK:
                {
                    return new Integration_Mesh_STK( aFileName, aSuppMeshData, aCreateFacesAndEdges );
                }
                default:
                {
                    MORIS_ERROR(false, "Cannot create specified integration mesh type from file.");
                    return nullptr;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Integration_Mesh* create_integration_mesh(
                enum MeshType aMeshType,
                MtkMeshData aMeshData)
        {
            switch (aMeshType)
            {
                case MeshType::STK:
                {
                    return new Integration_Mesh_STK( aMeshData );
                }
                default:
                {
                    MORIS_ERROR(false, "Cannot create specified integration mesh type with mesh data.");
                    return nullptr;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Integration_Mesh* create_integration_mesh(
                enum MeshType aMeshType,
                MtkMeshData aMeshData,
                Interpolation_Mesh *aInterpMesh)
        {
            switch (aMeshType)
            {
                case MeshType::STK:
                {
                    return new Integration_Mesh_STK( aMeshData, aInterpMesh );
                }
                default:
                {
                    MORIS_ERROR(false, "Cannot create specified integration mesh type with mesh data plus an interpolation mesh." );
                    return nullptr;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Integration_Mesh* create_integration_mesh_from_interpolation_mesh(
                enum MeshType aMeshType,
                Interpolation_Mesh* aInterpMesh,
                uint                aMeshIndex,
                Cell_Cluster_Input* aCellClusterData)
        {
            switch (aMeshType)
            {
                case MeshType::STK:
                {
                    return new Integration_Mesh_STK(*aInterpMesh, aCellClusterData);
                }
                case MeshType::HMR:
                {
                    return new hmr::Integration_Mesh_HMR(aMeshIndex, static_cast<hmr::Interpolation_Mesh_HMR*>(aInterpMesh));
                }
                default:
                {
                    MORIS_ERROR(false, "Cannot create specified integration mesh type from an interpolation mesh.");
                    return nullptr;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}


/*
 * cl_Mesh_Factory.hpp
 *
 *  Created on: Sep 18, 2018
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MESH_FACTORY_HPP_
#define PROJECTS_MTK_SRC_CL_MESH_FACTORY_HPP_

#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
// implementations
#include "stk_impl/cl_MTK_Mesh_STK.hpp"
#include "stk_impl/cl_MTK_Interpolation_Mesh_STK.hpp"
#include "stk_impl/cl_MTK_Integration_Mesh_STK.hpp"

namespace moris
{
namespace mtk
{

/**
 * Mesh constructor (mesh generated internally or obtained from an Exodus file )
 *
 * @param[in] aFileName  .................    String with mesh file name.
 */

inline Mesh*
create_mesh(enum MeshType  aMeshType,
            std::string    aFileName,
            MtkMeshData*   aSuppMeshData = nullptr,
            const bool     aCreateFacesAndEdges = true )
{
    Mesh* tMeshBase = nullptr;

    // Mark mesh data as being supplementary to an input file
    if(aSuppMeshData != nullptr)
    {
        aSuppMeshData->SupplementaryToFile = true;
    }

    switch (aMeshType)
    {
        case(MeshType::STK):
        {
            tMeshBase = new Mesh_STK( aFileName, aSuppMeshData, aCreateFacesAndEdges );
            break;
        }
        default:
        {
            MORIS_ASSERT( 0, "Specified mesh type not supported by MORIS or this construction method not implemented" );
        }
    }
    return tMeshBase;
}

inline Mesh*
create_mesh(enum MeshType aMeshType,
            MtkMeshData   aMeshData )
{
    Mesh* tMeshBase = nullptr;
    switch (aMeshType)
    {
        case(MeshType::STK):
        {
            tMeshBase = new Mesh_STK( aMeshData );
            break;
        }
        default:
        {
            MORIS_ASSERT( 0, "Specified mesh type not supported by MORIS or this construction method not implemented" );
        }
    }
    return tMeshBase;
}

Interpolation_Mesh*
create_interpolation_mesh(enum MeshType  aMeshType,
                 std::string    aFileName,
                 MtkMeshData*   aSuppMeshData = nullptr,
                 const bool     aCreateFacesAndEdges = true)
{
    // Mark mesh data as being supplementary to an input file
    if(aSuppMeshData != nullptr)
    {
        aSuppMeshData->SupplementaryToFile = true;
    }
    Interpolation_Mesh* tMesh = new Interpolation_Mesh_STK( aFileName, aSuppMeshData, aCreateFacesAndEdges );

    return tMesh;
}


Integration_Mesh*
create_integration_mesh(enum MeshType  aMeshType,
                        std::string    aFileName,
                        MtkMeshData*   aSuppMeshData = nullptr,
                        const bool     aCreateFacesAndEdges = true)
{
    // Mark mesh data as being supplementary to an input file
    if(aSuppMeshData != nullptr)
    {
        aSuppMeshData->SupplementaryToFile = true;
    }
    Integration_Mesh* tMesh = new Integration_Mesh_STK( aFileName, aSuppMeshData, aCreateFacesAndEdges );

    return tMesh;
}

inline Integration_Mesh*
create_integration_mesh(enum MeshType aMeshType,
                        MtkMeshData   aMeshData )
{
    Integration_Mesh* tMeshBase = nullptr;
    switch (aMeshType)
    {
        case(MeshType::STK):
        {
            tMeshBase = new Integration_Mesh_STK( aMeshData );
            break;
        }
        default:
        {
            MORIS_ASSERT( 0, "Specified mesh type not supported by MORIS or this construction method not implemented" );
        }
    }
    return tMeshBase;
}

Integration_Mesh*
create_integration_mesh_from_interpolation_mesh(enum MeshType       aMeshType,
                                                Interpolation_Mesh* aInterpMesh)
{
    MORIS_ERROR(aMeshType == MeshType::STK,"create_integration_mesh_from_interpolation_mesh only currently set up for STK meshes");
    return new Integration_Mesh_STK(*aInterpMesh);
}

}
}



#endif /* PROJECTS_MTK_SRC_CL_MESH_FACTORY_HPP_ */

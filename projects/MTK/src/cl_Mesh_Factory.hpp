/*
 * cl_Mesh_Factory.hpp
 *
 *  Created on: Sep 18, 2018
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MESH_FACTORY_HPP_
#define PROJECTS_MTK_SRC_CL_MESH_FACTORY_HPP_

#include "cl_MTK_Mesh.hpp"
// implementations
#include "stk_impl/cl_MTK_Mesh_STK.hpp"

namespace moris
{
namespace mtk
{

/**
 * Mesh constructor (mesh generated internally or obtained from an Exodus file )
 *
 * @param[in] aFileName  .................    String with mesh file name.
 */

Mesh*
create_mesh(enum MeshType aMeshType,
            std::string    aFileName,
            MtkSetsInfo*   aSetsInfo = nullptr,
            const bool     aCreateFacesAndEdges = true )
{
    Mesh* tMeshBase = nullptr;
    switch (aMeshType)
    {
        case(MeshType::STK):
        {
            tMeshBase = new Mesh_STK( aFileName, aSetsInfo, aCreateFacesAndEdges );
            break;
        }
        default:
        {
            MORIS_ASSERT( 0, "Specified mesh type not supported by MORIS or this construction method not implemented" );
        }
    }
    return tMeshBase;
}

Mesh*
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

}
}



#endif /* PROJECTS_MTK_SRC_CL_MESH_FACTORY_HPP_ */

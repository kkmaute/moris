/*
 * cl_MTK_Mesh_Data_Input.hpp
 *
 *  Created on: Sep 17, 2018
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_MESH_DATA_INPUT_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_MESH_DATA_INPUT_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Sets_Info.hpp"
#include "cl_MTK_Fields_Info.hpp"

namespace moris
{
namespace mtk
{
    //////////////////////////
    // STRUC FOR BLOCK SET  //
    //////////////////////////
    struct MtkMeshData
    {
        uint*        SpatialDim ;
        Matrix < IdMat >* ElemConn;
        Matrix < IdMat >* EntProcOwner;
        Matrix < DDRMat >* NodeCoords;
        Matrix < IdMat >* LocaltoGlobalElemMap;
        Matrix < IdMat >* LocaltoGlobalNodeMap;
        bool         CreateAllEdgesAndFaces;
        MtkFieldsInfo* FieldsInfo;
        MtkSetsInfo* SetsInfo;
        real         TimeStamp = 0.0;

        MtkMeshData():
            SpatialDim(),
            ElemConn(),
            EntProcOwner(),
            NodeCoords(),
            LocaltoGlobalElemMap(),
            LocaltoGlobalNodeMap(),
            CreateAllEdgesAndFaces(false),
            FieldsInfo(),
            SetsInfo(){}
    };
}
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_MESH_DATA_INPUT_HPP_ */

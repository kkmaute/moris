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
        uint*                           SpatialDim ;
        moris::Cell<Matrix < IdMat >*>  ElemConn;
        Matrix < IdMat >*               EntProcOwner;
        Matrix < DDRMat >*              NodeCoords;
        moris::Cell<Matrix < IdMat >*>  LocaltoGlobalElemMap;
        Matrix < IdMat >*               LocaltoGlobalNodeMap;
        bool                            CreateAllEdgesAndFaces;
        MtkFieldsInfo*                  FieldsInfo;
        MtkSetsInfo*                    SetsInfo;
        real                            TimeStamp = 0.0;
        bool                            AutoAuraOptionInSTK = true;

        MtkMeshData(uint aNumElementTypes):
            SpatialDim(),
            ElemConn(aNumElementTypes),
            EntProcOwner(),
            NodeCoords(),
            LocaltoGlobalElemMap(aNumElementTypes),
            LocaltoGlobalNodeMap(),
            CreateAllEdgesAndFaces(true),
            FieldsInfo(),
            SetsInfo(nullptr)
        {

        }

        MtkMeshData():
            SpatialDim(),
            ElemConn(1),
            EntProcOwner(),
            NodeCoords(),
            LocaltoGlobalElemMap(1),
            LocaltoGlobalNodeMap(),
            CreateAllEdgesAndFaces(true),
            FieldsInfo(),
            SetsInfo(nullptr)
        {

        }

        bool
        has_mesh_sets()
        {
            bool tAnswer = true;
            if ( SetsInfo == NULL )
            {
                return false;
            }

            return tAnswer;
        }

        /*
         * count the number of elements provided
         */
        uint
        get_num_elements()
        {
            uint tNumElements = 0;
            for(uint i = 0; i<ElemConn.size(); i++)
            {
                tNumElements += ElemConn(i)->n_rows();
            }
            return tNumElements;
        }

        /*
         * Number of elements in the element map
         */
        uint
        size_local_to_global_elem_map()
        {
            uint tSizeLocalToGlobal = 0;
            for(uint i = 0; i <LocaltoGlobalElemMap.size(); i++)
            {
                tSizeLocalToGlobal += LocaltoGlobalElemMap(i)->numel();
            }
            return tSizeLocalToGlobal;
        }

        /*
         * Collapse the element cell to a continuous vector
         */
        Matrix< IdMat >
        collapse_element_map()
        {
            Matrix< IdMat > tCollapsedMap(this->get_num_elements(),1);

            uint tCount = 0;
            for(uint i = 0; i <LocaltoGlobalElemMap.size(); i++)
            {
                for(uint j = 0; j<LocaltoGlobalElemMap(i)->n_rows(); j++)
                {
                    tCollapsedMap(tCount) = (*LocaltoGlobalElemMap(i))(j);
                    tCount++;
                }
            }

            return tCollapsedMap;
        }

    };
}
}

#endif /* PROJECTS_MTK_SRC_CL_MTK_MESH_DATA_INPUT_HPP_ */

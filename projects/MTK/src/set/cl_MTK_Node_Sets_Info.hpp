/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Node_Sets_Info.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_NODE_SETS_INFO_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_NODE_SETS_INFO_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"
namespace moris
{
    namespace mtk
    {
        //////////////////////////
        // STRUC FOR NODE SET  //
        //////////////////////////
        /*
         * A node set requires the following information
         *
         * Node ids - the ids of the nodes in the node set
         * Node set name - name of the node set
         */
        struct MtkNodeSetInfo
        {
            Matrix< IdMat >* mNodeIds;
            std::string      mNodeSetName;

            MtkNodeSetInfo()
                    : mNodeIds()
                    , mNodeSetName()
            {
            }

            bool
            nodeset_has_name()
            {
                return !mNodeSetName.empty();
            }
        };

        inline std::ostream&
        operator<<( std::ostream& os, mtk::MtkNodeSetInfo const * const & dt )
        {
            os << "Vertex Set Name:" << dt->mNodeSetName << " | Number of Vertices: " << dt->mNodeIds->numel();

            return os;
        }
    }    // namespace mtk
}    // namespace moris

#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_NODE_SETS_INFO_HPP_ */

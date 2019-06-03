/*
 * cl_MTK_Side_Set_STK.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: schmidt
 */

#ifndef SRC_MESH_CL_MTK_SIDE_SET_STK_HPP_
#define SRC_MESH_CL_MTK_SIDE_SET_STK_HPP_

#include <string>

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Map.hpp"
#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Cell.hpp" //MTK/src

#include "cl_MTK_Cell_Cluster.hpp" //MTK/src
#include "cl_MTK_Side_Set.hpp" //MTK/src

namespace moris
{
    namespace mtk
    {

//------------------------------------------------------------------------------
        class Side_Set_STK : public Side_Set
        {
//        private :
//            moris::Matrix< DDUMat > mMyBlockSetClusterInds;
//            moris::Cell<Cell_Cluster const *> mBlockSetClusters;
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Side_Set_STK( moris::Cell<Side_Cluster const *>  aBlockSetClusters ) : Side_Set(aBlockSetClusters)
            {
            };

//------------------------------------------------------------------------------

            /**
             * virtual destructor
             */
            ~Side_Set_STK(){};

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------
#endif /* SRC_MESH_CL_MTK_SIDE_SET_HPP_STK_ */

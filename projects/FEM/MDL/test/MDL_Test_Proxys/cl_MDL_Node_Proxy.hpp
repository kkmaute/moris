/*
 * cl_MDL_Node_Proxy.hpp
 *
 *  Created on: Mai 22, 2019
 *      Author: Schmidt
 */
#ifndef PROJECTS_GEN_SRC_CL_MDL_NODE_PROXY_HPP_
#define PROJECTS_GEN_SRC_CL_MDL_NODE_PROXY_HPP_

#include "cl_Matrix.hpp"
#include "cl_MTK_Vertex.hpp"

namespace moris
{
    namespace mdl
    {
        class Node : public mtk::Vertex
        {
//------------------------------------------------------------------------------
        private:
//            moris::Matrix< DDRMat > coord;
            moris::sint mIndex = -1;
            moris::sint mID = -1;

//------------------------------------------------------------------------------

        public:
//------------------------------------------------------------------------------

            Node( const sint aIndex, const sint aID ) : mIndex( aIndex ),
                                                        mID( aID )
            {};

//------------------------------------------------------------------------------

            ~Node(){};

//------------------------------------------------------------------------------
//            // gives the coordinates of node
//            Matrix< DDRMat > get_coords() const
//            {
//                return coord;
//            };
//
//            //------------------------------------------------------------------------------
//            // gives domain wide ID of node
//            moris_id
//            get_id() const
//            {
//                MORIS_ERROR( false, "get_id() not implemented " );
//                return gNoID;
//            };
//
//            //------------------------------------------------------------------------------
//            // gives index of node
//            moris_index
//            get_index() const
//            {
//                MORIS_ERROR( false, "get_index() not implemented " );
//                return 0;
//            };
//
//            //------------------------------------------------------------------------------
//            // gives owner of node
//            moris_index
//            get_owner() const
//            {
//                MORIS_ERROR( false, "get_index() not implemented " );
//                return 0;
//            };
//
//            //------------------------------------------------------------------------------
//            mtk::Vertex_Interpolation*
//            get_interpolation( const uint aOrder )
//            {
//                MORIS_ERROR( false, "get_interpolation() not implemented " );
//                return nullptr;
//            };
//
//            //------------------------------------------------------------------------------
//            const mtk::Vertex_Interpolation*
//            get_interpolation( const uint aOrder ) const
//            {
//                MORIS_ERROR( false, "get_interpolation() not implemented " );
//                return nullptr;
//            };
//
//            //------------------------------------------------------------------------------
//            uint
//            get_level() const
//            {
//                return( 0 );
//            };
//
//            //------------------------------------------------------------------------------
//            void
//            flag()
//            {
//                MORIS_ERROR( false, "flag() not implemented " );
//            };
//
//            //------------------------------------------------------------------------------
//            void
//            unflag() const
//            {
//                MORIS_ERROR( false, "unflag() not implemented " );
//            };
//
//            //------------------------------------------------------------------------------
//            bool
//            is_flagged() const
//            {
//                MORIS_ERROR( false, "is_flagged() not implemented " );
//                return false;
//            };


        };

    } /* namespace mdl */

} /* namespace moris */




#endif /* PROJECTS_GEN_SRC_CL_MDL_NODE_PROXY_HPP_ */

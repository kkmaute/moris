/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_NodeProxy.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_NODEPROXY_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_NODEPROXY_HPP_

#include "cl_Matrix.hpp"
#include "cl_MTK_Vertex.hpp"

namespace moris::fem
{

    class NodeProxy : public mtk::Vertex
    {
    private:
        moris::Matrix< DDRMat > coord;
        sint mindex = -1;

    public:
        NodeProxy( double xval, double yval , uint index)
        {
            coord.set_size(1,2,1.0);
            coord( 0, 0 ) = xval;
            coord( 0, 1 ) = yval;
            mindex = index;
        }

    //------------------------------------------------------------------------------
        ~NodeProxy() override{};

    //------------------------------------------------------------------------------
        Matrix< DDRMat > get_coords() const override
        {
            return coord;
        };

    //------------------------------------------------------------------------------
        moris_id get_id() const override
        {
            MORIS_ERROR( false," NodeProxy - get_id - Function not implemented. ");
            return gNoID;
        };

    //------------------------------------------------------------------------------
        moris_index get_index() const override
        {
            return mindex;
        }

    //------------------------------------------------------------------------------
        moris_index get_owner() const override
        {
            MORIS_ERROR( false," NodeProxy - get_owner - Function not implemented. " );
            return 0;
        };

    //------------------------------------------------------------------------------
        mtk::Vertex_Interpolation * get_interpolation( const uint aBSplineMeshIndex ) override
        {
            MORIS_ERROR( false," NodeProxy - get_interpolation - Function not implemented. " );
            return nullptr;
        }

    //------------------------------------------------------------------------------
        mtk::Vertex_Interpolation * get_interpolation( const uint aBSplineMeshIndex ) const override
        {
            MORIS_ERROR( false," NodeProxy - get_interpolation - Function not implemented. " );
            return nullptr;
        }

    //------------------------------------------------------------------------------
        uint get_level() const override
        {
            return 0;
        };

    //------------------------------------------------------------------------------
        void flag() override
        {
            MORIS_ERROR( false," NodeProxy - flag - Function not implemented. " );
        };

    //------------------------------------------------------------------------------
        void unflag() override
        {
            MORIS_ERROR( false," NodeProxy - unflag - Function not implemented. " );
        };

    //------------------------------------------------------------------------------
        bool is_flagged() const
        {
            MORIS_ERROR( false," NodeProxy - is_flagged - Function not implemented. " );
            return false;
        };

    //------------------------------------------------------------------------------
    };

}    // namespace moris::fem

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_NODEPROXY_HPP_ */


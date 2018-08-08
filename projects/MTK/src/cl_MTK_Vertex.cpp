#include "assert.hpp"
#include "cl_MTK_Vertex.hpp" //MTK/src

//------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        /**
         * returns a moris::Mat with node coordinates
         */
        Mat< real >
        Vertex::get_coords() const
        {
            MORIS_ERROR( false, "get_coords() not implemented for this vertex.");
            Mat< real > aEmpty;
            return aEmpty;
        }

//-----------------------------------------------------------------------------

        /**
         * returns the domain wide id of this vertex
         */
        luint
        Vertex::get_id() const
        {
            MORIS_ERROR( false, "get_id() not implemented for this vertex.");
            return MORIS_LUINT_MAX;
        }

//------------------------------------------------------------------------------

        /**
         * returns the B-Spline IDs of this vertex
         */
        Mat< luint >
        Vertex::get_bspline_ids() const
        {
            MORIS_ERROR( false, "get_bspline_ids() not implemented for this vertex.");
            Mat< luint > aEmpty;
            return aEmpty;
        }

//------------------------------------------------------------------------------

        /**
         * returns the T-Matrix of this vertex
         */
        Mat< real >
        Vertex::get_t_matrix() const
        {
            MORIS_ERROR( false, "get_t_matrix() not implemented for this vertex.");
            Mat< real > aEmpty;
            return aEmpty;
        }

//------------------------------------------------------------------------------

        /**
         * returns the id of the proc that owns this vertex
         */
        uint
        Vertex::get_owner() const
        {
            MORIS_ERROR( false, "get_owner() not implemented for this vertex.");
            return MORIS_UINT_MAX;
        }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

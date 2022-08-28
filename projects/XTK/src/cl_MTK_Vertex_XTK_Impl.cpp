/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_XTK_Impl.cpp
 *
 */

#include "cl_MTK_Vertex_XTK_Impl.hpp"
#include "cl_XTK_Background_Mesh.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        Vertex_XTK::Vertex_XTK(){}

        Vertex_XTK::Vertex_XTK(
            moris::moris_id           aVertexId,
            moris::moris_index        aVertexIndex,
            moris::moris_index        aOwner,
            std::shared_ptr<moris::Matrix<moris::DDRMat>> aCoordinates)
        :
        mVertexId(aVertexId),
        mVertexIndex(aVertexIndex),
        mVertexOwner(aOwner),
        mCoordinates(aCoordinates)
        {
        }

        //------------------------------------------------------------------------------

        Vertex_XTK::~Vertex_XTK(){}

        //------------------------------------------------------------------------------

        void
        Vertex_XTK::set_vertex_interpolation(Vertex_Interpolation_XTK * aVertexInterpolation)
        {
            MORIS_ASSERT(aVertexInterpolation != nullptr, "aVertexInterpolation provided is a null ptr");
            mVertexInterpolation = aVertexInterpolation;
        }
        //------------------------------------------------------------------------------

        Matrix<DDRMat>
        Vertex_XTK::get_coords() const
        {
            return *mCoordinates;
        }
        //------------------------------------------------------------------------------
        moris_id
        Vertex_XTK::get_id() const
        {
            return mVertexId;
        }

        //------------------------------------------------------------------------------

        moris_index
        Vertex_XTK::get_index() const
        {
            return mVertexIndex;
        }

        //------------------------------------------------------------------------------
        moris_index
        Vertex_XTK::get_owner() const
        {
            return mVertexOwner;
        }
        //------------------------------------------------------------------------------
        Vertex_Interpolation *
        Vertex_XTK::get_interpolation(const uint aOrder)
        {
            MORIS_ASSERT(mVertexInterpolation != nullptr, "mInterpolation is a null ptr");
            return mVertexInterpolation;
        }

        //------------------------------------------------------------------------------
        const Vertex_Interpolation *
        Vertex_XTK::get_interpolation(const uint aOrder) const
        {
            MORIS_ASSERT(mVertexInterpolation != nullptr, "mInterpolation is a null ptr");
            return mVertexInterpolation;
        }
        //------------------------------------------------------------------------------
        uint
        Vertex_XTK::get_level() const
        {
            return 0;
        }

        //------------------------------------------------------------------------------
        size_t
        Vertex_XTK::capacity()
        {
            size_t tTotal = 0;
            tTotal += sizeof(mVertexId);
            tTotal += sizeof(mVertexIndex);
            tTotal += sizeof(mVertexInterpolation);
            tTotal += sizeof(mCoordinates);
            return tTotal;
        }
        //------------------------------------------------------------------------------

    } // namespace mtk
} // namespace moris


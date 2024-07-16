/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_ISC_Impl.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_VERTEX_ISC_IMPL_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_VERTEX_ISC_IMPL_HPP_

#include "cl_MTK_Vertex.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//#include "cl_MTK_Vertex_Interpolation.hpp"

//------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------
        class Vertex_ISC: public Vertex
        {
            private:
                moris::moris_id            mVertexId;
                moris::moris_index         mVertexIndex;
                Matrix<DDRMat>             mVertexCoords;
//                Vertex_Interpolation*      mVertexInterpolation = nullptr;

            public:

                //------------------------------------------------------------------------------
                /**
                 * trivial constructor
                 */
                Vertex_ISC(){};

                //------------------------------------------------------------------------------

                Vertex_ISC(moris::moris_id        aVertexId,
                        moris::moris_index     aVertexIndex,
                        moris::Matrix<moris::DDRMat>         aVertexCoords):
                            mVertexId(aVertexId),
                            mVertexIndex(aVertexIndex),
                            mVertexCoords(aVertexCoords)
                {};

                //------------------------------------------------------------------------------

                /**
                 * Destructor
                 */
                ~Vertex_ISC(){};

                //------------------------------------------------------------------------------
                /**
                 * returns a moris::Matrix with node coordinates
                 */
                Matrix< DDRMat >
                get_coords() const
                {
                    return mVertexCoords;
                }

                //------------------------------------------------------------------------------
                /**
                 * returns a moris::Matrix with node coordinates
                 */
                Matrix< DDRMat > const &
                get_coords_by_ref() const
                {
                    return mVertexCoords;
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the domain wide id of this vertex
                 */
                moris_id
                get_id() const
                {
                    return mVertexId;
                }

                //------------------------------------------------------------------------------
                /**
                 * returns the processor unique index of this vertex
                 */
                moris_index
                get_index() const
                {
                    return mVertexIndex;
                }

                //------------------------------------------------------------------------------
                //
                //                // fixme: change this into moris_id
                //                // FIXME: add owner function in background mesh
                //                moris_index
                //                get_owner() const;

                //------------------------------------------------------------------------------
                //
                //                Vertex_Interpolation *
                //                get_interpolation( const uint aOrder ) const
                //                {
                //                    return mVertexInterpolation;;
                //                }
                //
                //                //------------------------------------------------------------------------------
                //
                //                const Vertex_Interpolation *
                //                get_interpolation( const uint aOrder ) const
                //                {
                //                    return mVertexInterpolation; ;
                //                }

                //------------------------------------------------------------------------------

                //                uint
                //                get_level() const;
                //
                //                //------------------------------------------------------------------------------
                //
                //                size_t
                ////                capacity();
                //
        };
    } /* namespace mtk */
} /* namespace moris */

#endif /* PROJECTS_MTK_SRC_CL_MTK_VERTEX_ISC_IMPL_HPP_ */


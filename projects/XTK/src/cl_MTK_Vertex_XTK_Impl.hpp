/*
 * cl_MTK_Vertex_XTK_Impl.hpp
 *
 *  Created on: Mar 8, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_MTK_VERTEX_XTK_IMPL_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_MTK_VERTEX_XTK_IMPL_HPP_

#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Vertex_Interpolation_XTK_Impl.hpp"

namespace xtk
{
 class Background_Mesh;
}

//------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        class Vertex_XTK: public Vertex
        {
        protected :
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Vertex_XTK();


            Vertex_XTK(moris::moris_id        aVertexId,
                       moris::moris_index     aVertexIndex,
                       xtk::Background_Mesh * aBackgroundMeshPtr);

            Vertex_XTK(moris::moris_id        aVertexId,
                       moris::moris_index     aVertexIndex,
                       std::shared_ptr<moris::Matrix<moris::DDRMat>> aCoordinates);                       
            /*
             * Constructor for a background mesh vertex
             */
            Vertex_XTK(mtk::Vertex* aBackgroundMeshVertex);
//------------------------------------------------------------------------------

            /**
             * Destructor
             */
            ~Vertex_XTK();


//------------------------------------------------------------------------------
            void
            set_vertex_interpolation(Vertex_Interpolation_XTK* aVertexInterpolation);

//------------------------------------------------------------------------------

            /**
             * returns a moris::Matrix with node coordinates
             */
            Matrix< DDRMat >
            get_coords() const;

//------------------------------------------------------------------------------

            /**
             * returns the domain wide id of this vertex
             */
            moris_id
            get_id() const;

//------------------------------------------------------------------------------

            /**
             * returns the processor unique index of this vertex
             */
            moris_index
            get_index() const;

//------------------------------------------------------------------------------

            // fixme: change this into moris_id
            // FIXME: add owner function in background mesh
            moris_index
            get_owner() const;

//------------------------------------------------------------------------------

            Vertex_Interpolation *
            get_interpolation( const uint aOrder );

//------------------------------------------------------------------------------

            const Vertex_Interpolation *
            get_interpolation( const uint aOrder ) const;

//------------------------------------------------------------------------------

            uint
            get_level() const;

//------------------------------------------------------------------------------

            size_t
            capacity();

        private:
            moris::moris_id              mVertexId;
            moris::moris_index           mVertexIndex;
            Vertex_Interpolation_XTK *   mVertexInterpolation = nullptr; // (basis weights and basis identity)

            // If this vertex is craeted by XTK we need a pointer to the background mesh
            xtk::Background_Mesh * mBackgroundMeshPtr = nullptr; // To access coords... (set to null for assertion purposes)

            // If this was a vertex in the background mesh (simply store the pointer)
            Vertex * mBackgroundMeshVertex = nullptr;
            
            std::shared_ptr<moris::Matrix<moris::DDRMat>> mCoordinates;
        };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------



#endif /* PROJECTS_XTK_SRC_XTK_CL_MTK_VERTEX_XTK_IMPL_HPP_ */

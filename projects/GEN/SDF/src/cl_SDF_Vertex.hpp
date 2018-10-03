/*
 * cl_SDF_Vertex.hpp
 *
 *  Created on: Oct 1, 2018
 *      Author: messe
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_VERTEX_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_VERTEX_HPP_

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Vertex.hpp"

namespace moris
{
    namespace sdf
    {
// -----------------------------------------------------------------------------

        /**
         * The sdf vertex is a wrapper around an MTK vertex.
         * It contains a pointer to the MTK vertex and
         * has the ability to flag nodes
         */
        class Vertex
        {
            //! pointer to underlying MTK vertex
            const mtk::Vertex * mVertex;

            //! flag telling if vertex is inside
            bool                mIsInside = false;

            //! flag telling if an SDF has been calculated for this vertex
            bool                mHasSDF = false;

            bool                mIsCandidate = false;

// -----------------------------------------------------------------------------
        public:
// -----------------------------------------------------------------------------

            /**
             * constructor
             */
            Vertex( const mtk::Vertex * aVertex ) : mVertex( aVertex )
            {

            }

// -----------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~Vertex(){};

// -----------------------------------------------------------------------------

            Matrix< F31RMat >
            get_coords() const
            {
                // convert moris DDRMat to F31Mat
                auto tCoords = mVertex->get_coords();

                Matrix< F31RMat > aCoords;
                aCoords( 0 ) = tCoords( 0 );
                aCoords( 1 ) = tCoords( 1 );
                aCoords( 2 ) = tCoords( 2 );
            }

// -----------------------------------------------------------------------------

            void
            set_inside_flag()
            {
                mIsInside = true;
            }

// -----------------------------------------------------------------------------

            void
            unset_inside_flag()
            {
                mIsInside = false;
            }

// -----------------------------------------------------------------------------

            bool
            is_inside() const
            {
                return mIsInside;
            }

// -----------------------------------------------------------------------------

            void
            set_candidate_flag()
            {
                mIsCandidate = true;
            }

// -----------------------------------------------------------------------------

            void
            unset_candidate_flag()
            {
                mIsCandidate = false;
            }

// -----------------------------------------------------------------------------

            bool
            is_candidate() const
            {
                return mIsCandidate;
            }

// -----------------------------------------------------------------------------

            void
            set_sdf_flag()
            {
                mHasSDF = true;
            }

 // -----------------------------------------------------------------------------

            void
            unset_sdf_flag()
            {
                mHasSDF = false;
            }

// -----------------------------------------------------------------------------

            bool
            has_sdf() const
            {
                return mHasSDF;
            }

// -----------------------------------------------------------------------------

            moris_index
            get_index() const
            {
                return mVertex->get_index();
            }

// -----------------------------------------------------------------------------
        };

// -----------------------------------------------------------------------------
    }
}



#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_VERTEX_HPP_ */

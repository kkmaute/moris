/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_Interpolation_XTK_Impl.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_MTK_VERTEX_INTERPOLATION_XTK_IMPL_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_MTK_VERTEX_INTERPOLATION_XTK_IMPL_HPP_

#include "cl_Matrix.hpp"

#include "cl_MTK_Vertex.hpp"

namespace xtk
{
    class Vertex_Enrichment;
}

namespace moris
{
    namespace mtk
    {

        class Vertex_Interpolation_XTK : public Vertex_Interpolation
        {

          public:
            //------------------------------------------------------------------------------

            Vertex_Interpolation_XTK(){};

            //------------------------------------------------------------------------------

            ~Vertex_Interpolation_XTK(){};

            //------------------------------------------------------------------------------

            void
            set_vertex_enrichment( xtk::Vertex_Enrichment* aVertexEnrichment );
            /**
             * returns the IDs of the interpolation coefficients
             */
            Matrix< IdMat >
            get_ids() const
            {
                MORIS_ERROR( 0, "get_ids not implemented in xtk vertex interpolation" );
                return moris::Matrix< IdMat >( 0, 0 );
            }

            //------------------------------------------------------------------------------

            /**
             * returns the indices of the interpolation coefficients
             */
            Matrix< IndexMat >
            get_indices() const;

            //------------------------------------------------------------------------------

            /**
             * returns the proc owners of the IDs of this vertex
             */
            Matrix< IdMat >
            get_owners() const
            {
                MORIS_ERROR( 0, "get_owners not implemented in xtk vertex interpolation" );
                return moris::Matrix< IdMat >( 0, 0 );
            }

            //------------------------------------------------------------------------------
            /**
             * set the interpolation weights
             */
            void
            set_weights( const moris::Matrix< DDRMat >& aWeights )
            {
                MORIS_ERROR( 0, "set_weights not implemented in xtk vertex interpolation" );
            }

            //------------------------------------------------------------------------------

            /**
             * returns the interpolation weights
             */
            const Matrix< DDRMat >*
            get_weights() const;

            //------------------------------------------------------------------------------
            /**
             * set the coefficient objects
             */
            void
            set_coefficients( Vector< Vertex* >& aCoefficients )
            {
                MORIS_ERROR( 0, "set_coefficients not implemented in xtk vertex interpolation" );
            }

            //------------------------------------------------------------------------------

            /**
             * returns the pointers to the coefficient objects
             */
            Vector< Vertex* >&
            get_coefficients()
            {
                MORIS_ERROR( 0, "get_coefficients not implemented in xtk vertex interpolation" );
                return mCoefficients;
            }

            //------------------------------------------------------------------------------

            /**
             * returns the pointers to the coefficient objects (const version)
             */
            const Vector< Vertex* >&
            get_coefficients() const
            {
                MORIS_ERROR( 0, "get_coefficients not implemented in xtk vertex interpolation" );
                return mCoefficients;
            }

            //------------------------------------------------------------------------------

            /**
             * returns the number of coefficients attributed to this basis
             */
            uint
            get_number_of_coefficients() const
            {
                MORIS_ERROR( 0, "get_number_of_coefficients not implemented in xtk vertex interpolation" );
                return 0;
            }

            //------------------------------------------------------------------------------

          private:
            xtk::Vertex_Enrichment* mVertexEnrichment = nullptr;
        };
    }    // namespace mtk
}    // namespace moris

#endif /* PROJECTS_XTK_SRC_XTK_CL_MTK_VERTEX_INTERPOLATION_XTK_IMPL_HPP_ */

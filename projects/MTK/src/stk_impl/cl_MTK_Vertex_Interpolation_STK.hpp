/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_Interpolation_STK.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_MTK_VERTEX_INTERPOLATION_STK_HPP_
#define PROJECTS_HMR_SRC_CL_MTK_VERTEX_INTERPOLATION_STK_HPP_

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"

namespace moris
{
    namespace mtk
    {
        class Vertex_Interpolation_STK : public mtk::Vertex_Interpolation
        {
            mtk::Vertex* mVertex = nullptr;

            // ----------------------------------------------------------------------------
        public:
            // ----------------------------------------------------------------------------
            Vertex_Interpolation_STK(){};

            Vertex_Interpolation_STK( mtk::Vertex * aVertex ) : mVertex(aVertex)
            {
                mWeights.set_size( 1, 1, 1.0 );
            };

            // ----------------------------------------------------------------------------

            ~Vertex_Interpolation_STK(){};

            // ----------------------------------------------------------------------------
            /**
             * sets the values of the T-Matrix
             */
            void
            set_weights( const Matrix< DDRMat > & aWeights );

            // ----------------------------------------------------------------------------

            /**
             * return the interpolation weights
             */
            const Matrix< DDRMat > *
            get_weights() const;

            // ----------------------------------------------------------------------------

            /**
             * sets the coefficients of this basis
             */
            void
            set_coefficients( moris::Cell< mtk::Vertex* > & aCoefficients );

            // ----------------------------------------------------------------------------

            /**
             * returns the coefficients of this basis
             */
            moris::Cell< mtk::Vertex* > &
            get_coefficients();

            // ----------------------------------------------------------------------------

            /**
             * returns the coefficients of this basis ( const version )
             */
            const moris::Cell< mtk::Vertex* > &
            get_coefficients() const;

            // ----------------------------------------------------------------------------

            /**
             * returns the number of coefficients attributed to this basis
             */
            uint
            get_number_of_coefficients() const;

            // ----------------------------------------------------------------------------

            /**
             * returns the IDs of the interpolation coefficients
             */
            Matrix< IdMat >
            get_ids() const;

            // ----------------------------------------------------------------------------

            /**
             * returns the Indices of the interpolation coefficients
             */
            Matrix< IndexMat >
            get_indices() const;
            //            {
            //                MORIS_ERROR( false, "Node_Interpolation_STK::get_indices - not implemented");
            //                moris::Matrix< IndexMat > tEmptyMatrix;
            //                return tEmptyMatrix;
            //             }

            // ----------------------------------------------------------------------------

            /**
             * returns the owners of the interpolation coefficients
             */
            Matrix< IdMat >
            get_owners() const;

            // ----------------------------------------------------------------------------
        };
    }
}

#endif /* PROJECTS_HMR_SRC_CL_MTK_NODE_INTERPOLATION_STK_HPP_ */


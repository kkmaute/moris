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

#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"

namespace moris::mtk
{
    class Vertex_Interpolation_STK : public mtk::Vertex_Interpolation
    {
        mtk::Vertex* mVertex = nullptr;

        // ----------------------------------------------------------------------------

      public:
        // ----------------------------------------------------------------------------
        Vertex_Interpolation_STK(){};

        Vertex_Interpolation_STK( mtk::Vertex* aVertex )
                : mVertex( aVertex )
        {
            mWeights.set_size( 1, 1, 1.0 );
        };

        // ----------------------------------------------------------------------------

        ~Vertex_Interpolation_STK() override{};

        // ----------------------------------------------------------------------------
        /**
         * sets the values of the T-Matrix
         */
        void
        set_weights( const Matrix< DDRMat >& aWeights ) override;

        // ----------------------------------------------------------------------------

        /**
         * return the interpolation weights
         */
        const Matrix< DDRMat >*
        get_weights() const override;

        // ----------------------------------------------------------------------------

        /**
         * sets the coefficients of this basis
         */
        void
        set_coefficients( Vector< mtk::Vertex* >& aCoefficients ) override;

        // ----------------------------------------------------------------------------

        /**
         * returns the coefficients of this basis
         */
        Vector< mtk::Vertex* >&
        get_coefficients() override;

        // ----------------------------------------------------------------------------

        /**
         * returns the coefficients of this basis ( const version )
         */
        const Vector< mtk::Vertex* >&
        get_coefficients() const override;

        // ----------------------------------------------------------------------------

        /**
         * returns the number of coefficients attributed to this basis
         */
        uint
        get_number_of_coefficients() const override;

        // ----------------------------------------------------------------------------

        /**
         * returns the IDs of the interpolation coefficients
         */
        Matrix< IdMat >
        get_ids() const override;

        // ----------------------------------------------------------------------------

        /**
         * returns the Indices of the interpolation coefficients
         */
        Matrix< IndexMat >
        get_indices() const override;
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
        get_owners() const override;

        // ----------------------------------------------------------------------------
    };
    }

#endif /* PROJECTS_HMR_SRC_CL_MTK_NODE_INTERPOLATION_STK_HPP_ */


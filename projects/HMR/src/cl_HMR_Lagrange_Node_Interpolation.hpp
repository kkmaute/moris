/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Node_Interpolation.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_LAGRANGE_NODE_INTERPOLATION_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_LAGRANGE_NODE_INTERPOLATION_HPP_

#include "cl_HMR_Basis.hpp"
#include "moris_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"

namespace moris::hmr
{
    class Lagrange_Node_Interpolation : public mtk::Vertex_Interpolation
    {

// ----------------------------------------------------------------------------
public:
// ----------------------------------------------------------------------------

        Lagrange_Node_Interpolation(){};

// ----------------------------------------------------------------------------

        ~Lagrange_Node_Interpolation(){};

// ----------------------------------------------------------------------------
        /**
         * sets the values of the T-Matrix
         */
        void set_weights( const Matrix< DDRMat > & aWeights );

// ----------------------------------------------------------------------------

        /**
         * return the interpolation weights
         */
        const Matrix< DDRMat > * get_weights() const;

// ----------------------------------------------------------------------------

        /**
         * sets the coefficients of this basis
         */
        void
        set_coefficients( Cell< mtk::Vertex* > & aCoefficients );

// ----------------------------------------------------------------------------

        /**
         * returns the coefficients of this basis
         */
        moris::Cell< mtk::Vertex* > & get_coefficients();

// ----------------------------------------------------------------------------

        /**
         * returns the coefficients of this basis ( const version )
         */
        const moris::Cell< mtk::Vertex* > & get_coefficients() const;

// ----------------------------------------------------------------------------

        /**
         * returns the number of coefficients attributed to this basis
         */
        uint get_number_of_coefficients() const;

// ----------------------------------------------------------------------------

        /**
         * returns the IDs of the interpolation coefficients
         */
        Matrix< IdMat > get_ids() const;

// ----------------------------------------------------------------------------

        /**
         * returns the Indices of the interpolation coefficients
         */
        Matrix< IndexMat > get_indices() const;

// ----------------------------------------------------------------------------

        /**
         * returns the owners of the interpolation coefficients
         */
        Matrix< IdMat > get_owners() const;

        /**
         * returns the owners of the interpolation coefficients
         * these ids are not consecutive and iriginally created as luint.
         * consider using luint here for large problems
         */
        Matrix< IdMat > get_ijkl_id() const;

// ----------------------------------------------------------------------------
    };
}

#endif /* PROJECTS_HMR_SRC_CL_HMR_LAGRANGE_NODE_INTERPOLATION_HPP_ */


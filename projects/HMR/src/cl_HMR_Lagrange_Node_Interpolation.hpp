/*
 * cl_HMR_Lagrange_Node_Interpolation.hpp
 *
 *  Created on: Sep 17, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_LAGRANGE_NODE_INTERPOLATION_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_LAGRANGE_NODE_INTERPOLATION_HPP_

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"
#include "cl_HMR_Basis.hpp"

namespace moris
{
    namespace hmr
    {
        class Lagrange_Node_Interpolation : public mtk::Vertex_Interpolation
        {
            Cell< mtk::Vertex* > mCoefficients;
            Matrix< DDRMat >          mWeights;

// ----------------------------------------------------------------------------
    public:
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
            set_coefficients( Cell< mtk::Vertex* > & aCoefficients );

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
            Mat< moris_id >
            get_ids() const;

// ----------------------------------------------------------------------------

            /**
             * returns the Indices of the interpolation coefficients
             */
            Mat< moris_index >
            get_indices() const;

// ----------------------------------------------------------------------------

            /**
             * returns the owners of the interpolation coefficients
             */
            Matrix< DDUMat >
            get_owners() const;

// ----------------------------------------------------------------------------
        };
    }
}



#endif /* PROJECTS_HMR_SRC_CL_HMR_LAGRANGE_NODE_INTERPOLATION_HPP_ */

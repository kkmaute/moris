/*
 * cl_FEM_Element_Bulk.hpp
 *
 *  Created on: Apr 22, 2018
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_BULK_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_BULK_HPP_

#include "assert.h"
#include "cl_FEM_Element.hpp"               //FEM/INT/src

namespace moris
{
    namespace fem
    {
    class Element_Block;
//------------------------------------------------------------------------------
    /**
     * \brief element class that communicates with the mesh interface
     */
    class Element_Bulk : public Element
    {

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * constructor
         *
         * @param[ in ]     pointer to mesh interface object
         * @param[ in ]     cell of pointers to integrand of weak form of governing eqs.
         * @param[ in ]     cell of pointer to fem nodes
         * @param[ in ]     Pointer to element block
         */
        Element_Bulk( mtk::Cell    const *        aCell,
                      Element_Block      * aElementBlock,
                      Cluster            * aCluster );
//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Element_Bulk();

//------------------------------------------------------------------------------

        void compute_jacobian();

//------------------------------------------------------------------------------

        void compute_residual();

//------------------------------------------------------------------------------

        void compute_jacobian_and_residual();

//------------------------------------------------------------------------------

//        real compute_integration_error( real (*aFunction)( const Matrix< DDRMat > & aPoint ) );

//------------------------------------------------------------------------------

//        real compute_element_average_of_scalar_field();

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ELEMENT_BULK_HPP_ */

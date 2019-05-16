/*
 * cl_FEM_Element_Double_Sideset.hpp
 *
 *  Created on: May 13, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_DOUBLE_SIDESET_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_DOUBLE_SIDESET_HPP_

#include "assert.h"
#include "cl_FEM_Element.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
    class Element_Block;
//------------------------------------------------------------------------------
    /**
     * \brief Element_Sideset class
     */
    class Element_Double_Sideset : public Element
    {

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * constructor
         * @param[ in ]     pointer to mesh cell
         * @param[ in ]     pointer to element block
         * @param[ in ]     pointer to cluster
         *
         */
        Element_Double_Sideset( mtk::Cell     const * aMasterIGCell,
                                mtk::Cell     const * aSlaveIGCell,
                                Element_Block       * aElementBlock,
                                Cluster             * aCluster);

//------------------------------------------------------------------------------
        /**
         * destructor
         */
        ~Element_Double_Sideset();

//------------------------------------------------------------------------------
        /**
         * compute jacobian over the element
         */
        void compute_jacobian();

//------------------------------------------------------------------------------
        /**
         * compute residual over the element
         */
        void compute_residual();

//------------------------------------------------------------------------------
        /**
         * compute jacobian and residual over the element
         */
        void compute_jacobian_and_residual();

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_ELEMENT_DOUBLE_SIDESET_HPP_ */

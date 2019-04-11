/*
 * cl_FEM_Element_Sideset.hpp
 *
 *  Created on: Mar 07, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_SIDESET_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_SIDESET_HPP_

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
    class Element_Sideset : public Element
    {
//------------------------------------------------------------------------------
    protected:
//------------------------------------------------------------------------------
//    // a member list of side ordinals
//    moris::Cell< moris_index > mListOfSideOrdinals;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * constructor
         *
         * @param[ in ]     pointer to mesh interface object
         * @param[ in ]     cell of pointers to integrand of weak form of governing eqs.
         * @param[ in ]     cell of pointer to fem nodes
         */
        Element_Sideset( mtk::Cell                  * aCell,
                         moris::Cell< IWG* >        & aIWGs,
                         moris::Cell< Node_Base* >  & aNodes );

        Element_Sideset( mtk::Cell                  * aCell,
                         moris::Cell< IWG* >        & aIWGs,
                         moris::Cell< Node_Base* >  & aNodes,
                         Element_Block      * aElementBlock);

//------------------------------------------------------------------------------
        /**
         * destructor
         */
        ~Element_Sideset();

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
    protected:
//------------------------------------------------------------------------------
//        /**
//         * create the field interpolators for the element
//         */
//        moris::Cell< Field_Interpolator* > create_field_interpolators
//            ( Geometry_Interpolator* aGeometryInterpolator );

//------------------------------------------------------------------------------
//        /**
//         * set the field interpolators coefficients
//         */
//        void set_field_interpolators_coefficients
//            ( Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_ELEMENT_SIDESET_HPP_ */

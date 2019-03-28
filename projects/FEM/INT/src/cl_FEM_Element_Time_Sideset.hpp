/*
 * cl_FEM_Element_Time_Sideset.hpp
 *
 *  Created on: Mar 19, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_TIME_SIDESET_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_TIME_SIDESET_HPP_

#include "assert.h"
#include "cl_FEM_Element.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
    /**
     * \brief Element_Sideset class
     */
    class Element_Time_Sideset : public Element
    {

//------------------------------------------------------------------------------
    protected:
//------------------------------------------------------------------------------
//        // a member list of time side ordinals
//        moris::Cell< moris_index > mListOfTimeOrdinals;

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
        Element_Time_Sideset( mtk::Cell                  * aCell,
                              moris::Cell< IWG* >        & aIWGs,
                              moris::Cell< Node_Base* >  & aNodes );
//------------------------------------------------------------------------------
        /**
         * destructor
         */
        ~Element_Time_Sideset();

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
//        Cell< Field_Interpolator* > create_field_interpolators
//            ( Geometry_Interpolator* aGeometryInterpolator );

//------------------------------------------------------------------------------
//        /**
//         * set the field interpolators coefficients
//         */
//        void
//        set_field_interpolators_coefficients
//        ( Cell< Field_Interpolator* > & aFieldInterpolators );
//

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_ELEMENT_TIME_SIDESET_HPP_ */

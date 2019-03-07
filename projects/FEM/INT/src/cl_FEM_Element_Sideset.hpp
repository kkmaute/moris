/*
 * cl_FEM_Element_Sideset.hpp
 *
 *  Created on: Mar 07, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_SIDESET_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_SIDESET_HPP_

#include "assert.h"
#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MSI_Equation_Object.hpp"       //FEM/MSI/src

#include "cl_FEM_Enums.hpp"                 //FEM/INT/src
#include "cl_FEM_Node.hpp"                  //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Element.hpp"               //FEM/INT/src

namespace moris
{

    // forward declaration of mtk classes
    namespace mtk
    {
        class Cell;
        enum class Geometry_Type;
        enum class Interpolation_Order;
    }

//------------------------------------------------------------------------------

    namespace fem
    {
//------------------------------------------------------------------------------

    /**
     * \brief element class that communicates with the mesh interface
     */
    class Element_Sideset : public Element
    {

    protected:
        //moris::Cell< mtk::Cell > mSidesetCell;

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
        Element_Sideset( mtk::Cell                 * aCell,
                         moris::Cell< IWG* >       & aIWGs,
                         moris::Cell< Node_Base* > & aNodes );
//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Element_Sideset();

//------------------------------------------------------------------------------

        void compute_jacobian();

//------------------------------------------------------------------------------

        void compute_residual();

//------------------------------------------------------------------------------

        void compute_jacobian_and_residual();

//------------------------------------------------------------------------------
    protected:
//------------------------------------------------------------------------------
        /**
         * create the field interpolators for the element
         */
        Cell< Field_Interpolator* >
        create_element_field_interpolators
        ( Geometry_Interpolator* aGeometryInterpolator );

//------------------------------------------------------------------------------
        /**
         * set the field interpolators coefficients
         */
        void
        set_element_field_interpolators_coefficients
        ( Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
        /**
         * get the field interpolators for an IWG
         */
        Cell< Field_Interpolator* >
        get_IWG_field_interpolators( IWG*                        & aIWG,
                                     Cell< Field_Interpolator* > & aFieldInterpolators );
//------------------------------------------------------------------------------
        /**
         * set the initial sizes and values for mJacobianElement and mResidualElement
         */
        void
        initialize_mJacobianElement_and_mResidualElement
        ( Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_ELEMENT_SIDESET_HPP_ */

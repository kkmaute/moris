/*
 * cl_FEM_Element_Bulk.hpp
 *
 *  Created on: Jul 16, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_BULK_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_BULK_HPP_

#include "assert.h"
#include "cl_FEM_Element.hpp"               //FEM/INT/src

namespace moris
{
    namespace fem
    {

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
         */
        Element_Bulk( mtk::Cell                 * aCell,
                      moris::Cell< IWG* >       & aIWGs,
                      moris::Cell< Node_Base* > & aNodes );
//------------------------------------------------------------------------------


        /**
         * trivial destructor
         */
        ~Element_Bulk();
//------------------------------------------------------------------------------

//        /**
//         * returns the element geometry type
//         */
//        mtk::Geometry_Type get_geometry_type() const ;

//------------------------------------------------------------------------------

//        /**
//         * returns the node coordinates of the element
//         */
//        Matrix< DDRMat > get_node_coords() const ;

//------------------------------------------------------------------------------

//        /**
//         * returns the interpolation order of the underlying cell
//         */
//        mtk::Interpolation_Order get_interpolation_order() const;

//------------------------------------------------------------------------------

//        /**
//         * returns a moris::Mat with ids of vertices that are connected to this element
//         */
//        Matrix< IdMat > get_vertex_ids() const;

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
    protected:
//------------------------------------------------------------------------------
        /**
         * create the field interpolators for the element
         */
        moris::Cell< Field_Interpolator* >
        create_element_field_interpolators
        ( Geometry_Interpolator* aGeometryInterpolator );

//------------------------------------------------------------------------------
        /**
         * set the field interpolators coefficients
         */
        void
        set_element_field_interpolators_coefficients
        ( moris::Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
        /**
         * get the field interpolators for an IWG
         */
        moris::Cell< Field_Interpolator* >
        get_IWG_field_interpolators( IWG*                               & aIWG,
                                     moris::Cell< Field_Interpolator* > & aFieldInterpolators );
//------------------------------------------------------------------------------
        /**
         * set the initial sizes and values for mJacobianElement and mResidualElement
         */
        void
        initialize_mJacobianElement_and_mResidualElement
        ( moris::Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ELEMENT_BULK_HPP_ */

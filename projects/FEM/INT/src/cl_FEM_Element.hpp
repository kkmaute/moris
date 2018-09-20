/*
 * cl_FEM_Element.hpp
 *
 *  Created on: Jul 16, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_HPP_

#include "assert.h"

#include "typedefs.hpp"           //MRS/COR/src

#include "cl_FEM_Interpolation_Matrix.hpp"
#include "cl_FEM_Enums.hpp"       //FEM/INT/src
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#include "cl_FEM_Node.hpp"         //FEM/INT/src
#include "cl_MSI_Equation_Object.hpp"

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
    class Element : public MSI::Equation_Object
    {
        //! pointer to cell on mesh
        mtk::Cell * mCell;

        //! pointer to IWG object
        IWG       * mIWG;

        //! weak BCs of element
        moris::Mat< moris::real > mNodalWeakBCs;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * constructor
         *
         * @param[ in ]     pointer to mesh interface object
         * @param[ in ]     pointer to integrand of weak form of governing eqs.
         */
        Element(
                mtk::Cell * aCell,
                IWG * aIWG,
                Cell< Node_Base* > & aNodes,
                const Matrix< DDRMat >  & aNodalWeakBCs );

//------------------------------------------------------------------------------

        /**
         * trivial destructor
         */
        ~Element(){};
//------------------------------------------------------------------------------

        /**
         * returns the elememt geometry type
         */
        mtk::Geometry_Type
        get_geometry_type() const ;

//------------------------------------------------------------------------------

        /**
         * returns the node coordinates of the element
         */
        Matrix<DDRMat>
        get_node_coords() const ;

//------------------------------------------------------------------------------

        /**
         * returns the interpolation order of the underlying cell
         */
        mtk::Interpolation_Order
        get_interpolation_order() const;

//------------------------------------------------------------------------------

        /**
         * returns a moris::Mat with ids of vertices that are connected to this element
         */
        Matrix< IdMat >
        get_vertex_ids() const;

//------------------------------------------------------------------------------

        void
        compute_jacobian_and_residual();

//------------------------------------------------------------------------------

        real
        compute_integration_error(
                real (*aFunction)( const Matrix< DDRMat > & aPoint ) );

//------------------------------------------------------------------------------

        /**
         * how many nodes are connected to this element
         */
        uint
        get_number_of_nodes() const;

//------------------------------------------------------------------------------

        //Mat< moris_index >
        //get_adof_indices();

//------------------------------------------------------------------------------
    protected:
//------------------------------------------------------------------------------

        /**
         * auto detect full interpolation scheme
         */
        Integration_Order
        get_auto_integration_order();

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */




#endif /* SRC_FEM_CL_FEM_ELEMENT_HPP_ */

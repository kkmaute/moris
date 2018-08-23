#include <iostream>
#include "cl_FEM_Element.hpp" //FEM/INT/src

//#include "cl_FEM_Integrator_old.hpp"

#include "op_times.hpp" //LNA/src
#include "op_plus.hpp" //LNA/src
#include "fn_det.hpp" //LNA/src
#include "cl_MTK_Vertex.hpp"
#include "cl_FEM_Integration_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Matrix.hpp" //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp" //FEM/INT/src
#include "cl_FEM_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Element::Element(
                mtk::Cell * aCell,
                IWG * aIWG,
                Cell< msi::Node* > & aNodes,
                const Mat< real >  & aNodalWeakBCs ) :
                Equation_Object(),
                mCell( aCell ),
                mIWG( aIWG )
        {
            // get vertices from cell
            Cell< mtk::Vertex* > tVertices = aCell->get_vertex_pointers();

            // get number of nodes from Cell
            uint tNumberOfNodes = tVertices.size();

            // assign node object
            mNodeObj.resize( tNumberOfNodes, nullptr );

            // reset counter
            uint tCount = 0;

            // set size of Weak BCs
            mNodalWeakBCs.set_size( tNumberOfNodes, 1, 0.0 );

            // fill node objects
            for( auto tVertex : tVertices )
            {
                // get index from vertex
                auto tIndex = tVertex->get_index();
                mNodeObj( tCount ) = aNodes( tIndex );
                mNodalWeakBCs( tCount++) = aNodalWeakBCs( tIndex );
            }

            // FIXME: Mathias, please comment
            mTimeSteps.set_size( 1, 1, 0 );

            // set the dof types of this element
            mEqnObjDofTypeList = mIWG->get_dof_types();

            // FIXME: Mathias, please comment
            mPdofValues.set_size( tNumberOfNodes, 1, 0.0 );



            //
            this->compute_jacobian_and_residual();
        }

//------------------------------------------------------------------------------

        mtk::Geometry_Type
        Element::get_geometry_type() const
        {
            return mCell->get_geometry_type();
        }
//------------------------------------------------------------------------------

        mtk::Interpolation_Order
        Element::get_interpolation_order() const
        {
            return mCell->get_interpolation_order();
        }

//------------------------------------------------------------------------------

        Mat< real >
        Element::get_node_coords() const
        {
            return mCell->get_vertex_coords();
        }

//------------------------------------------------------------------------------


        Integration_Order
        Element::get_auto_integration_order()
        {
            switch( this->get_geometry_type() )
            {
                case( mtk::Geometry_Type::QUAD ) :
                {
                     return Integration_Order::QUAD_3x3;
                     break;
                }
                case( mtk::Geometry_Type::HEX ) :
                {
                    return Integration_Order::HEX_3x3x3;
                }
                default :
                {
                    MORIS_ERROR( false, "get_integration_order() not defined for this geometry type");
                    return Integration_Order::UNDEFINED;
                }
            }
        }

//------------------------------------------------------------------------------
        void
        Element::compute_jacobian_and_residual()
        {
            // create field interpolation rule
            Interpolation_Rule tFieldInterpolationRule(
                    this->get_geometry_type(),
                    Interpolation_Type::LAGRANGE,
                    this->get_interpolation_order() ); // <- add second type in order
                                                      //    to interpolate in space
                                                      //    and time

            // create geometry interpolation rule
            Interpolation_Rule tGeometryInterpolationRule(
                    this->get_geometry_type(),
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // create integration rule
            Integration_Rule tIntegration_Rule(
                    this->get_geometry_type(),
                    Integration_Type::GAUSS,
                    this->get_auto_integration_order()
                    );

            // set number of fields
            uint tNumberOfFields = 1;

            // create interpolator
            Interpolator tInterpolator(
                    this,
                    tNumberOfFields,
                    tFieldInterpolationRule,
                    tGeometryInterpolationRule,
                    tIntegration_Rule );

            // get number of points
            auto tNumberOfIntegrationPoints
                = tInterpolator.get_number_of_integration_points();

            // get number of nodes
            auto tNumberOfNodes = tInterpolator.get_number_of_dofs();

            // mass matrix
            mJacobian.set_size( tNumberOfNodes, tNumberOfNodes, 0.0 );
            mResidual.set_size( tNumberOfNodes, 1, 0.0 );

            Mat< real > tJacobian( tNumberOfNodes, tNumberOfNodes );
            Mat< real > tResidual( tNumberOfNodes, 1 );

            mIWG->create_matrices( &tInterpolator );

            for( uint k=0; k<tNumberOfIntegrationPoints; ++k )
            {
                // evaluate shape function at given integration point
                mIWG->compute_jacobian_and_residual(
                        tJacobian,
                        tResidual,
                        mPdofValues,
                        mNodalWeakBCs,
                        k );


                mJacobian = mJacobian + tJacobian*tInterpolator.get_det_J( k )
                          *tInterpolator.get_integration_weight( k );

                mResidual = mResidual + tResidual*tInterpolator.get_det_J( k )
                                      * tInterpolator.get_integration_weight( k );

            }

            //mJacobian.print("J");
            //mResidual.print("R");
            // close IWG object
            mIWG->delete_matrices();
        }

//------------------------------------------------------------------------------

        moris::Cell< mtk::Vertex* >
        Element::get_vertex_pointers()
        {
            return mCell->get_vertex_pointers();
        }

//------------------------------------------------------------------------------

        Mat< luint >
        Element::get_vertex_ids() const
        {
            return mCell->get_vertex_ids();
        }

//------------------------------------------------------------------------------

        uint
        Element::get_number_of_nodes() const
        {
            return mCell->get_number_of_vertices();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

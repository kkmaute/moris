#include <iostream>
#include "cl_FEM_Element.hpp" //FEM/INT/src

//#include "cl_FEM_Integrator_old.hpp"

#include "op_times.hpp" //LNA/src
#include "op_plus.hpp" //LNA/src
#include "fn_det.hpp" //LNA/src
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

        Element::Element( mtk::Cell * aCell ) : mCell( aCell )
        {
            // if the FEM element is constructed, we assume that we are
            // interested in the T-Matrix of this cell
            aCell->set_t_matrix_flag();
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

        void
        Element::eval_mass( Mat< real > & aM )
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
                    Integration_Order::QUAD_3x3
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

            // create N-Matrix
            auto tN = tInterpolator.create_matrix( 0, 0, 1 );

            // get number of points
            auto tNumberOfIntegrationPoints
                = tInterpolator.get_number_of_integration_points();

            // get number of nodes
            auto tNumberOfNodes = tInterpolator.get_number_of_dofs();

            // mass matrix
            aM.set_size( tNumberOfNodes, tNumberOfNodes, 0.0 );

            for( uint k=0; k<tNumberOfIntegrationPoints; ++k )
            {
                tInterpolator.evaluate_matrix( tN, k );

                aM = aM + trans( tN )
                          *tN
                          *tInterpolator.get_det_J( k )
                          *tInterpolator.get_integration_weight( k );
            }
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

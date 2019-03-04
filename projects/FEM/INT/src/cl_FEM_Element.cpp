#include <iostream>
#include "cl_FEM_Element.hpp" //FEM/INT/src

//#include "cl_FEM_Integrator_old.hpp"
#include "op_times.hpp"
#include "op_plus.hpp"
#include "fn_det.hpp"
#include "fn_sort.hpp"
#include "fn_eye.hpp"

#include "cl_MTK_Vertex.hpp"
#include "cl_FEM_Integration_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Matrix.hpp" //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp" //FEM/INT/src
#include "cl_FEM_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src

#include "cl_FEM_Node.hpp"         //FEM/INT/src

#include "cl_MTK_Cell.hpp" //MTK/src
#include "cl_Vector.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Element::Element(
                mtk::Cell * aCell,
                IWG * aIWG,
                Cell< Node_Base* > & aNodes ) :
                //Equation_Object(),
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
            mNodalWeakBCs.set_size( tNumberOfNodes, 1 );

            // fill node objects
            for( auto tVertex : tVertices )
            {
                mNodeObj( tCount++ ) = aNodes( tVertex->get_index() );
            }

            // FIXME: Mathias, please comment
            mTimeSteps.set_size( 1, 1, 0 );

            // set the dof types of this element
            mEqnObjDofTypeList = mIWG->get_dof_types();

            // FIXME: Mathias, please comment
            mPdofValues.set_size( tNumberOfNodes, 1, 0.0 );

            //
            //this->compute_jacobian_and_residual();

            eye( tNumberOfNodes,tNumberOfNodes, mJacobian);

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

        Matrix< DDRMat >
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
        Element::compute_jacobian()
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
            Interpolator tInterpolator( this->get_node_coords(),
                                        tNumberOfFields,
                                        tFieldInterpolationRule,
                                        tGeometryInterpolationRule,
                                        tIntegration_Rule );

            // get number of points
            auto tNumberOfIntegrationPoints = tInterpolator.get_number_of_integration_points();

            // get number of nodes
            auto tNumberOfNodes = tInterpolator.get_number_of_dofs();

            // mass matrix
            mJacobian.set_size( tNumberOfNodes, tNumberOfNodes, 0.0 );
            mResidual.set_size( tNumberOfNodes, 1, 0.0 );

            Matrix< DDRMat > tJacobian( tNumberOfNodes, tNumberOfNodes );
            Matrix< DDRMat > tResidual( tNumberOfNodes, 1 );

            mIWG->create_matrices( &tInterpolator );

            this->get_my_pdof_values();
            // end update values

            moris::Cell< enum MSI::Dof_Type > tDofTypes( 1, MSI::Dof_Type::L2 );
            moris::Matrix< DDRMat > tDofTypePdofVal;
            this->get_my_pdof_values( tDofTypes, tDofTypePdofVal );

            for( uint k=0; k<tNumberOfIntegrationPoints; ++k )
            {
                // evaluate shape function at given integration point
                mIWG->compute_jacobian_and_residual(
                        tJacobian,
                        tResidual,
                        tDofTypePdofVal,
                        mNodalWeakBCs,
                        k );


                mJacobian = mJacobian + tJacobian.matrix_data()*tInterpolator.get_det_J( k )
                          *tInterpolator.get_integration_weight( k );

                //mResidual = mResidual + tResidual.matrix_data()*tInterpolator.get_det_J( k )
                //                      * tInterpolator.get_integration_weight( k );
            }

            //print(mJacobian,"J");
            //mResidual.print("R");
            // close IWG object
            mIWG->delete_matrices();
      }

//------------------------------------------------------------------------------

        void
        Element::compute_residual()
        {
            this->get_my_pdof_values();

            moris::Cell< enum MSI::Dof_Type > tDofTypes( 1, MSI::Dof_Type::L2 );
            moris::Matrix< DDRMat > tDofTypePdofVal;
            this->get_my_pdof_values( tDofTypes, tDofTypePdofVal );

            mResidual = mJacobian*( tDofTypePdofVal - mNodalWeakBCs );
        }
//------------------------------------------------------------------------------

        real
        Element::compute_integration_error(
                real (*aFunction)( const Matrix< DDRMat > & aPoint ) )
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
                    this->get_interpolation_order() );
                    //mtk::Interpolation_Order::LINEAR );

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
                    this->get_node_coords(),
                    tNumberOfFields,
                    tFieldInterpolationRule,
                    tGeometryInterpolationRule,
                    tIntegration_Rule );

            // get number of points
            auto tNumberOfIntegrationPoints
                = tInterpolator.get_number_of_integration_points();

            real aError = 0.0;

            mIWG->create_matrices( &tInterpolator );

            moris::Cell< enum MSI::Dof_Type > tDofTypes( 1, MSI::Dof_Type::L2 );
            moris::Matrix< DDRMat > tDofTypePdofVal;
            this->get_my_pdof_values( tDofTypes, tDofTypePdofVal );

            for( uint k=0; k<tNumberOfIntegrationPoints; ++k )
            {
                // evaluate shape function at given integration point
                aError += mIWG->compute_integration_error(
                            tDofTypePdofVal,
                            aFunction,
                            k ) * tInterpolator.get_det_J( k )
                            * tInterpolator.get_integration_weight( k );
            }

            //std::cout << "Element error " << aError << std::endl;
            mIWG->delete_matrices();

            return aError;
        }

//------------------------------------------------------------------------------

        Matrix< IdMat >
        Element::get_vertex_ids() const
        {
            return mCell->get_vertex_ids();
        }

//------------------------------------------------------------------------------

        real
        Element::compute_element_average_of_scalar_field()
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
                    this->get_node_coords(),
                    tNumberOfFields,
                    tFieldInterpolationRule,
                    tGeometryInterpolationRule,
                    tIntegration_Rule );

            // get number of points
            auto tNumberOfIntegrationPoints
                = tInterpolator.get_number_of_integration_points();

            mIWG->create_matrices( &tInterpolator );

            real aValue = 0.0;
            real tWeight = 0.0;

            for( uint k=0; k<tNumberOfIntegrationPoints; ++k )
            {
                real tScale = tInterpolator.get_integration_weight( k )
                              * tInterpolator.get_det_J( k );

                aValue +=
                        mIWG->interpolate_scalar_at_point( mNodalWeakBCs, k )
                        * tScale;

                tWeight += tScale;

            }

            // close IWG object
            mIWG->delete_matrices();

            return aValue / tWeight;

        }

        /*Mat< moris_index >
        Element::get_adof_indices()
        {
            return sort( mCell->get_adof_indices() );
        }*/

//------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

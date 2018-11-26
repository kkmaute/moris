/*
 * fn_HMR_Exec_perform_mapping.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_FN_HMR_EXEC_PERFORM_MAPPING_HPP_
#define PROJECTS_HMR_SRC_FN_HMR_EXEC_PERFORM_MAPPING_HPP_

#include <memory>

#include "typedefs.hpp"
#include "cl_Cell.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_unique.hpp"

#include "cl_MTK_Mapper.hpp"
#include "HMR_Globals.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

    void
    check_for_forbidden_fields( Cell< ParameterList > & aFieldParams )
    {
        for( auto tField : aFieldParams )
        {
            bool tFieldIsForbidden = false;

            std::string tLabel = tField.get< std::string >( "label" );

            /**
             * the following fields can not be mapped, because their data
             * can not be inquired using get_entity_field_value_real_scalar()
             */
            tFieldIsForbidden = tFieldIsForbidden || tLabel == "Element_Level";
            tFieldIsForbidden = tFieldIsForbidden || tLabel == "Element_Owner";
            tFieldIsForbidden = tFieldIsForbidden || tLabel == "Node_IDs";

            std::string tError = "Mapping of field " + tLabel + " is forbidden.";

            MORIS_ERROR( ! ( tFieldIsForbidden && tField.get< sint >( "perform_mapping" ) == 1 ),
                    tError.c_str() );
        }
    }

// -----------------------------------------------------------------------------

        void
        perform_mapping(
                HMR * aHMR,
                Cell< ParameterList >            & aFieldParams,
                Cell< std::shared_ptr< Field > > & aInputFields,
                Cell< std::shared_ptr< Field > > & aOutputFields )
        {
            // make sure that we only map allowed fields

            check_for_forbidden_fields( aFieldParams );

            // - - - - - - - - - - - - - - - - - - - - - -
            // step 1: find out which orders are needed
            // - - - - - - - - - - - - - - - - - - - - - -

            // number of input fields
            uint tNumberOfFields = aInputFields.size();

            // counter
            uint tCount = 0;

            // container for orders of fields
            Matrix< DDUMat > tInputFieldOrders( 2*tNumberOfFields, 1 );

            // loop over all fields
            for( uint f=0; f<tNumberOfFields; ++f )
            {
                // test if we want to map this field
                if( aFieldParams( f ).get< sint >( "perform_mapping" ) == 1 )
                {
                    tInputFieldOrders( tCount++ ) = aInputFields( f )->get_bspline_order();
                    tInputFieldOrders( tCount++ ) = aInputFields( f )->get_lagrange_order();
                }
            }

            // chop container
            tInputFieldOrders.resize( tCount, 1 );

            // make orders unique
            Matrix< DDUMat > tMeshOrders;
            unique( tInputFieldOrders, tMeshOrders );

            uint tNumberOfMappers = tMeshOrders.length();

            // create map for mappers
            Matrix< DDUMat > tMapperIndex( gMaxBSplineOrder+1, 1, MORIS_UINT_MAX );

            for( uint k = 0; k<tNumberOfMappers; ++k )
            {
                tMapperIndex( tMeshOrders( k ) ) = k;
            }

            // - - - - - - - - - - - - - - - - - - - - - -
            // step 2: create union meshes and mappers
            // - - - - - - - - - - - - - - - - - - - - - -

            Cell< std::shared_ptr< Mesh > > tUnionMeshes;
            Cell< std::shared_ptr< Mesh > > tInputMeshes;
            Cell< mapper::Mapper * > tMappers( tNumberOfMappers, nullptr );

            for( uint m=0; m<tNumberOfMappers; ++m )
            {
                // get pointer to input mesh
                tInputMeshes.push_back( aHMR->create_mesh(
                        tMeshOrders( m ),
                        aHMR->get_parameters()->get_lagrange_input_pattern() ) );

                // create union mesh from HMR object
                tUnionMeshes.push_back( aHMR->create_mesh(
                        tMeshOrders( m ),
                        aHMR->get_parameters()->get_union_pattern() ) );

                // create mapper
                tMappers( m ) = new mapper::Mapper( tUnionMeshes( m ) );
            }


            // - - - - - - - - - - - - - - - - - - - - - -
            // step 3: map and project fields
            // - - - - - - - - - - - - - - - - - - - - - -

            for( uint f=0; f<tNumberOfFields; ++f )
            {
                // test if we want to map this field
                if( aFieldParams( f ).get< sint >( "perform_mapping" ) == 1 )
                {

                    // get pointer to input field
                    std::shared_ptr< Field > tInputField = aInputFields( f );

                    // get order
                    uint tBSplineOrder = tInputField->get_bspline_order();

                    // get index of mapper
                    uint m = tMapperIndex( tBSplineOrder );

                    // get pointer to field on union mesh
                    std::shared_ptr< Field > tUnionField =  tUnionMeshes( m )->create_field(
                            tInputField->get_label(),
                            tBSplineOrder );


                    if( tInputField->get_lagrange_order() >= tBSplineOrder )
                    {
                        // interpolate field onto union mesh
                        aHMR->get_database()->interpolate_field(
                                aHMR->get_parameters()->get_lagrange_input_pattern(),
                                tInputField,
                                aHMR->get_parameters()->get_union_pattern(),
                                tUnionField );
                    }
                    else
                    {
                        // first, project field on mesh with correct order
                        std::shared_ptr< Field > tTemporaryField =
                                tInputMeshes( m )->create_field(
                                        tInputField->get_label(),
                                        tBSplineOrder );

                        aHMR->get_database()->change_field_order(
                                tInputField, tTemporaryField );

                        // now, interpolate this field onto the inion
                        aHMR->get_database()->interpolate_field(
                                aHMR->get_parameters()->get_lagrange_input_pattern(),
                                tTemporaryField,
                                aHMR->get_parameters()->get_union_pattern(),
                                tUnionField );
                    }

                    // perform mapping
                    tMappers( m )->perform_mapping(
                            tInputField->get_label(),
                            EntityRank::NODE,
                            tInputField->get_label(),
                            tUnionField->get_bspline_rank() );

                    // a small sanity test
                    MORIS_ASSERT(  tUnionField->get_coefficients().length()
                            == tUnionMeshes( m )->get_num_entities(
                                    mtk::order_to_entity_rank( tBSplineOrder ) ),
                                    "Number of B-Splines does not match" );

                    // get pointer to output mesh
                    std::shared_ptr< Mesh >  tOutputMesh = aHMR->create_mesh(
                            tInputField->get_lagrange_order(),
                            aHMR->get_parameters()->get_lagrange_output_pattern() );

                    // create output field
                    std::shared_ptr< Field >  tOutputField =
                            tOutputMesh->create_field(
                                    tInputField->get_label(),
                                    tBSplineOrder );

                    // move coefficients to output field
                    // fixme: to be tested with Eigen also
                    tOutputField->get_coefficients() = std::move( tUnionField->get_coefficients() );

                    // allocate nodes for output
                    tOutputField->get_node_values().set_size( tOutputMesh->get_num_nodes(), 1 );

                    // evaluate nodes
                    tOutputField->evaluate_node_values();

                    // append field to output cell
                    aOutputFields.push_back( tOutputField );

                }
            }

            // delete mappers
            for( mapper::Mapper * tMapper : tMappers )
            {
                delete tMapper;
            }
        }

// -----------------------------------------------------------------------------

    } /* namespace hmr */
} /* namespace moris */



#endif /* PROJECTS_HMR_SRC_FN_HMR_EXEC_PERFORM_MAPPING_HPP_ */

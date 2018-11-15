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

#include "HMR_Globals.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

        void
        perform_mapping(
                HMR * aHMR,
                Cell< ParameterList >            & aFieldParams,
                Cell< std::shared_ptr< Field > > & aInputFields,
                Cell< std::shared_ptr< Field > > & aOutputFields )
        {

            // reset output cell
            aOutputFields.clear();

            // initialize counter
            uint tCount = 0;

            for( std::shared_ptr< Field >  tField : aInputFields )
            {
                // test if we want to map against this field
                if( aFieldParams( tCount++ ).get< sint >( "perform_mapping" ) == 1 )
                {
                    // get pointer to union mesh
                    std::shared_ptr< Mesh > tUnionMesh = aHMR->create_mesh(
                            tField->get_lagrange_order(),
                            aHMR->get_parameters()->get_union_pattern() );

                    // create pointer to field on union mesh
                    std::shared_ptr< Field > tUnionField =  tUnionMesh->create_field(
                            tField->get_label(),
                            tField->get_bspline_order() );

                    // interpolate field onto union mesh
                    aHMR->get_database()->interpolate_field(
                            aHMR->get_parameters()->get_input_pattern(),
                            tField,
                            aHMR->get_parameters()->get_union_pattern(),
                            tUnionField );

                    // create mapper
                    // fixme: mapper can be changed so that it has to be created only once
                    mapper::Mapper tUnionMapper( tUnionMesh );

                    // perform mapping
                    tUnionMapper.perform_mapping(
                            tField->get_label(),
                            EntityRank::NODE,
                            tField->get_label(),
                            tUnionField->get_bspline_rank() );

                    // get pointer to output mesh
                    std::shared_ptr< Mesh >  tOutputMesh = aHMR->create_mesh(
                            tField->get_lagrange_order(),
                            aHMR->get_parameters()->get_output_pattern() );

                    // create output field
                    std::shared_ptr< Field >  tOutputField =
                            tOutputMesh->create_field(
                                    tField->get_label(),
                                    tField->get_bspline_order() );

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
        }

// -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */



#endif /* PROJECTS_HMR_SRC_FN_HMR_EXEC_PERFORM_MAPPING_HPP_ */

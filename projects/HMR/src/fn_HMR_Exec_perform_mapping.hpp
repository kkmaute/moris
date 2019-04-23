/*
 * fn_HMR_Exec_perform_mapping.hpp
 *
 *  Created on: Nov 14, 2018
 *      Author: messe
 */

#ifndef PROJECTS_HMR_SRC_FN_HMR_EXEC_PERFORM_MAPPING_HPP_
#define PROJECTS_HMR_SRC_FN_HMR_EXEC_PERFORM_MAPPING_HPP_

#include <memory>

#include "cl_HMR.hpp"
#include "cl_HMR_Arguments.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Paramfile.hpp"
#include "HMR_Globals.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_unique.hpp"

#include "cl_MTK_Mapper.hpp"
namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

    void
    check_for_forbidden_fields( Cell< std::shared_ptr< Field > > & aInputFields )
    {
        for( auto tField : aInputFields )
        {
            bool tFieldIsForbidden = false;

            const std::string & tLabel = tField->get_label();

            /**
             * the following fields can not be mapped, because their data
             * can not be inquired using get_entity_field_value_real_scalar()
             */
            tFieldIsForbidden = tFieldIsForbidden || tLabel == "Element_Level";
            tFieldIsForbidden = tFieldIsForbidden || tLabel == "Element_Owner";
            tFieldIsForbidden = tFieldIsForbidden || tLabel == "Node_IDs";

            std::string tError = "Mapping of field " + tLabel + " is forbidden.";

            MORIS_ERROR( ! tFieldIsForbidden, tError.c_str() );
        }
    }

// -----------------------------------------------------------------------------

        void
        perform_mapping(
                const Arguments                  & aArguments,
                const Paramfile                  & aParamfile,
                HMR                              * aHMR,
                Cell< std::shared_ptr< Field > > & aInputFields,
                Cell< std::shared_ptr< Field > > & aOutputFields )
        {
            // make sure that we only map allowed fields

            check_for_forbidden_fields( aInputFields );

            // - - - - - - - - - - - - - - - - - - - - - -
            // step 1: find out which orders are needed
            // - - - - - - - - - - - - - - - - - - - - - -

            // number of input fields
            uint tNumberOfFields = aInputFields.size();

            // counter
            uint tCount = 0;

            // container for orders of fields
            Matrix< DDUMat > tInputFieldOrders( 3*tNumberOfFields
                    + aParamfile.get_number_of_meshes(), 1 );

            // loop over all fields
            for( uint f=0; f<tNumberOfFields; ++f )
            {
                tInputFieldOrders( tCount++ ) = aInputFields( f )->get_bspline_order();
                tInputFieldOrders( tCount++ ) = aInputFields( f )->get_bspline_output_order();
                tInputFieldOrders( tCount++ ) = aInputFields( f )->get_lagrange_order();
            }

            // loop over all defined meshes
            for( uint m=0; m< aParamfile.get_number_of_meshes(); ++m )
            {
                tInputFieldOrders( tCount++ ) =  aParamfile.get_mesh_order( m );
            }

            // make orders unique
            Matrix< DDUMat > tMeshOrders;
            unique( tInputFieldOrders, tMeshOrders );

            uint tNumberOfMappers = tMeshOrders.length();

            // create map for mappers
            Matrix< DDUMat > tMapperIndex( gMaxBSplineOrder+1, 1, MORIS_UINT_MAX );

            // flag telling if we have a linear and quadratic mesh
            bool tHaveLinearMesh = false;
            bool tHaveQuadraticMesh = false;

            for( uint k = 0; k<tNumberOfMappers; ++k )
            {
                tMapperIndex( tMeshOrders( k ) ) = k;

                if( tMeshOrders( k ) == 1 )
                {
                    tHaveLinearMesh = true;
                }
                else if( tMeshOrders( k ) == 2 )
                {
                    tHaveQuadraticMesh = true;
                }
            }

            // - - - - - - - - - - - - - - - - - - - - - -
            // step 2: create union meshes and mappers
            // - - - - - - - - - - - - - - - - - - - - - -
            mtk::Mesh_Manager tMeshManager;
            Cell< std::shared_ptr< Interpolation_Mesh_HMR > > tUnionInterpMeshes;
            Cell< std::shared_ptr< Integration_Mesh_HMR > >   tUnionIntegMeshes;
            Cell< std::shared_ptr< Interpolation_Mesh_HMR > > tInputInterpMeshes;
            Cell< std::shared_ptr< Integration_Mesh_HMR > >   tInputIntegMeshes;
            Cell< mapper::Mapper * > tMappers( tNumberOfMappers, nullptr );

            for( uint m=0; m<tNumberOfMappers; ++m )
            {
                // get pointer to input interpolation mesh
                tInputInterpMeshes.push_back( aHMR->create_interpolation_mesh(
                                            tMeshOrders( m ),
                        aHMR->get_parameters()->get_lagrange_input_pattern() )  );

                // get pointer to input integration mesh
                tInputIntegMeshes.push_back( aHMR->create_integration_mesh(
                        tMeshOrders( m ),
                        aHMR->get_parameters()->get_lagrange_input_pattern()));

                // create union mesh from HMR object
                tUnionInterpMeshes.push_back( aHMR->create_interpolation_mesh(
                        tMeshOrders( m ),
                        aHMR->get_parameters()->get_union_pattern() ) );


                tUnionIntegMeshes.push_back( aHMR->create_integration_mesh(
                        tMeshOrders( m ),
                        aHMR->get_parameters()->get_union_pattern() ) );

                // add pairs to mesh manager
                moris::uint tMeshPairIndex = tMeshManager.register_mesh_pair(tUnionInterpMeshes(m).get(),tUnionIntegMeshes(m).get());

                // create mapper
                tMappers( m ) = new mapper::Mapper( &tMeshManager,tMeshPairIndex );
            }


            // - - - - - - - - - - - - - - - - - - - - - -
            // step 3: map and project fields
            // - - - - - - - - - - - - - - - - - - - - - -

            for( uint f=0; f<tNumberOfFields; ++f )
            {
                // get pointer to input field
                std::shared_ptr< Field > tInputField = aInputFields( f );

                // get order
                uint tBSplineOrder = tInputField->get_bspline_output_order();
                uint tLagrangeOrder = tInputField->get_lagrange_order();

                // get index of mapper, pick mesh with same order as output
                uint m = tMapperIndex( tBSplineOrder );

                // get pointer to field on union mesh
                std::shared_ptr< Field > tUnionField =  tUnionInterpMeshes( m )->create_field(
                        tInputField->get_label(),
                        tBSplineOrder );

                if( tLagrangeOrder == tBSplineOrder )
                {
                    // interpolate field onto union mesh
                    aHMR->get_database()->interpolate_field(
                            aHMR->get_parameters()->get_lagrange_input_pattern(),
                            tInputField,
                            aHMR->get_parameters()->get_union_pattern(),
                            tUnionField );

                    // copy field id
                    tUnionField->set_id( tInputField->get_id() );
                }
                else
                {
                    // first, project field on mesh with correct order
                    std::shared_ptr< Field > tTemporaryField =
                            tInputInterpMeshes( m )->create_field(
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

                    // copy field id
                    tUnionField->set_id( tInputField->get_id() );
                }

                // set alpha parameter of mapper
                tMappers( m )->set_l2_alpha(
                        aParamfile.get_field_params( f ).mL2alpha );

                // perform mapping
                tMappers( m )->perform_mapping(
                        tInputField->get_label(),
                        EntityRank::NODE,
                        tInputField->get_label(),
                        tUnionField->get_bspline_rank() );

                // a small sanity test
                MORIS_ASSERT(  tUnionField->get_coefficients().length()
                        == tUnionInterpMeshes( m )->get_num_entities(
                                mtk::order_to_entity_rank( tBSplineOrder ) ),
                                "Number of B-Splines does not match" );

                // get pointer to output mesh
                std::shared_ptr< Mesh >  tOutputMesh = aHMR->create_mesh(
                        tLagrangeOrder,
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

                // copy field id
                tOutputField->set_id( tInputField->get_id() );

                // append field to output cell
                aOutputFields.push_back( tOutputField );

                //if the field is higher order, also map it to linear
                if( tLagrangeOrder != 1 && tHaveLinearMesh )
                {
                    // get pointer to output mesh
                    std::shared_ptr< Mesh > tLinearMesh = aHMR->create_mesh( 1 );

                    // create a linear field
                    std::shared_ptr< Field > tLinearField = tLinearMesh->create_field( tOutputField->get_label(), tBSplineOrder );

                    // evaluate node values for linear field
                    tLinearField->evaluate_node_values( tOutputField->get_coefficients() );
                }

                // if the field is not quadratic, also map it to quadratic mesh
                if( tLagrangeOrder != 2 && tHaveQuadraticMesh )
                {
                    // get pointer to output mesh
                    std::shared_ptr< Mesh > tQuadraticMesh = aHMR->create_mesh( 2 );

                    // create a linear field
                    std::shared_ptr< Field > tQuadraticField = tQuadraticMesh->create_field( tOutputField->get_label(), tBSplineOrder );

                    // evaluate node values for linear field
                    tQuadraticField->evaluate_node_values( tOutputField->get_coefficients() );
                }
            }

            // dump union meshes

            // test if union mesh is to be dumped
            if ( aParamfile.get_union_mesh_path().size() > 0 )
            {
                // grab basepath
                std::string tBasepath = aParamfile.get_union_mesh_path().substr(
                        0,aParamfile.get_union_mesh_path().find_last_of(".") );

                // loop over all meshes
                for( uint m=0; m<tNumberOfMappers; ++m )
                {
                    // note: can't dump cubic meshes to exodus
                    if( tMeshOrders( m ) < 3 )
                    {
                        // path to mesh
                        std::string tMeshPath = tBasepath + "_" + std::to_string( tMeshOrders( m ) ) + ".exo";

                        // get index of mesh
                        uint tMeshIndex = aHMR->get_mesh_index( tMeshOrders( m ) ,
                                aHMR->get_parameters()->get_union_pattern() );

                        // dump mesh (assume timestep to be zero )
                        aHMR->save_to_exodus( tMeshIndex, tMeshPath, aArguments.get_timestep() );
                    }
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

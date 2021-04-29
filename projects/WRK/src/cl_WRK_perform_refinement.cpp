#include "cl_WRK_Performer.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_WRK_perform_refinement.hpp"
#include "HMR_Globals.hpp"

#include "cl_MTK_Field.hpp"

#include "cl_Param_List.hpp"

namespace moris
{
    namespace wrk
    {

        Refinement_Mini_Performer::Refinement_Mini_Performer(
                ParameterList                 & aParameterlist,
                std::shared_ptr< Library_IO >   aLibrary )
        : mLibrary( aLibrary )
        {
            // set field names
            moris::Cell< std::string > tFieldNames;
            string_to_cell(
                    aParameterlist.get< std::string >( "field_names" ),
                    tFieldNames );

            mParameters.mFieldNames = tFieldNames;

            // set refinement level
            Cell< Matrix< DDSMat > > tRefinementLevel;
            string_to_cell_mat(
                    aParameterlist.get< std::string >( "levels_of_refinement" ),
                    tRefinementLevel );

            mParameters.mRefinementLevel = tRefinementLevel;

            // set refinementpattern
            Cell< Matrix< DDSMat > >  tRefinementPattern;
            string_to_cell_mat(
                    aParameterlist.get< std::string >( "refinement_pattern" ),
                    tRefinementPattern );

            mParameters.mRefinementPattern = tRefinementPattern;

        }

        //--------------------------------------------------------------------------------------------------------------

        void Refinement_Mini_Performer::perform_refinement(
                Cell< std::shared_ptr< mtk::Field > >  & aFields,
                std::shared_ptr<hmr::HMR>                aHMR )
        {
            // create field name to index map
            moris::map< std::string, moris_index >   tFieldNameToIndexMap;

            for( uint Ik = 0; Ik < aFields.size(); Ik++ )
            {
                tFieldNameToIndexMap[ aFields( Ik )->get_label() ] = Ik;
            }

            Cell< moris_index >                       tPattern;
            moris::Cell< moris::Cell< std::string > > tFieldNames;
            moris::Cell< moris::Cell< uint > >        tRefinements;
            moris::Cell< sint >                       tMaxRefinementPerLevel;

            this->prepare_input_for_refinement(
                    tPattern,
                    tFieldNames,
                    tRefinements,
                    tMaxRefinementPerLevel);

            for( uint Ik = 0; Ik< tPattern.size(); Ik++ )
            {
                // get pattern
                uint tActivationPattern = tPattern( Ik );

                for( uint Ii = 0; Ii< tMaxRefinementPerLevel.size(); Ii++ )
                {
                    // Mesh has changed after first refinement and therefore has to be rebuilt
                    // This will only work with analytic fields. When remeshing the L2 id done in the remesh mini performer
                    if( Ii > 0 )
                    {
                        // Interpolation_Mesh_HMR * tInterpolationMesh = aHMR->create_interpolation_mesh( tOrder, tPattern, tPattern );
                        ////   mtk::Mesh_Pair

                        for( uint Ia = 0; Ia< tFieldNames( Ik ).size(); Ia++ )
                        {
                            //  mtk::Field * tField = aFields( aFieldNameToIndexMap.find( tFieldNames( Ia ) ) );
                            //  tField->
                        }
                    }

                    for( uint Ia = 0; Ia< tFieldNames( Ik ).size(); Ia++ )
                    {
                        std::shared_ptr< mtk::Field > tField = aFields( tFieldNameToIndexMap.find( tFieldNames( Ik )( Ia ) ) );

                        uint tOrder = tField->get_lagrange_order();

                        const Matrix< DDRMat > & tFieldValues = tField->get_nodal_values();

                        // Put elements on queue and set flag for refinement
                        aHMR->based_on_field_put_elements_on_queue(
                                tFieldValues,
                                tActivationPattern,
                                tOrder,
                                -1);    //FieldValues, Patter,Order, function pointer index
                    }
                }

                aHMR->perform_refinement( tActivationPattern );
                aHMR->update_refinement_pattern( tActivationPattern );

                //                while ( true )
                //                {
                //                    uint tNumElements = 0;
                //
                //                    for( uint Ia = 0; Ia< tFieldNames.size(); Ia++ )
                //                    {
                //                        mtk::Field * tField = aFields( tFieldNameToIndexMap.find( tFieldNames( Ia ) ) );
                //
                //                        const Matrix< DDRMat > & tFieldValues = tField->get_nodal_values();
                //
                //                        tNumElements += aHMR->based_on_field_put_low_level_elements_on_queue(
                //                                tFieldValues,
                //                                tActivationPattern,
                //                                tOrder,
                //                                0);
                //                    }
                //
                //                    if( tNumElements == 0 )
                //                    {
                //                        break;
                //                    }
                //                    else
                //                    {
                //                        aHMR->perform_refinement( tActivationPattern );
                //                        aHMR->update_refinement_pattern( tActivationPattern );
                //                    }
                //
                //                    break;
                //                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Refinement_Mini_Performer::prepare_input_for_refinement(
                Cell< moris_index >                       & aPatternForRefinement,
                moris::Cell< moris::Cell< std::string > > & aFieldsForRefinement,
                moris::Cell< moris::Cell< uint > >        & aRefinements,
                moris::Cell< sint >                       & aMaxRefinementPerPattern )
        {

            // produce unique list of pattern which will be refined
            for( uint Ik = 0; Ik< mParameters.mRefinementPattern.size(); Ik++ )
            {
                for( uint Ii = 0; Ii< mParameters.mRefinementPattern( Ik ).numel(); Ii++ )
                {
                    aPatternForRefinement.push_back( mParameters.mRefinementPattern( Ik )( Ii ) );
                }
            }

            // Sort this created list
            std::sort( ( aPatternForRefinement.data() ).data(), ( aPatternForRefinement.data() ).data() + aPatternForRefinement.size() );

            // use std::unique and std::distance to create list containing all used dof types. This list is unique
            auto last = std::unique( ( aPatternForRefinement.data() ).data(), ( aPatternForRefinement.data() ).data() + aPatternForRefinement.size() );
            auto pos  = std::distance( ( aPatternForRefinement.data() ).data(), last );

            aPatternForRefinement.resize( pos );

            uint tNumberOfRefinementPattern = aPatternForRefinement.size();

            // resize
            aFieldsForRefinement    .resize( tNumberOfRefinementPattern );
            aRefinements            .resize( tNumberOfRefinementPattern );
            aMaxRefinementPerPattern.resize( tNumberOfRefinementPattern, 0 );

            // create list with field pointers and refinements per pattern
            for( uint Ik = 0; Ik< tNumberOfRefinementPattern; Ik++ )
            {
                moris_index tPattern = aPatternForRefinement( Ik );

                // loop over all fields and corresponding patterns. Find the pattern which corresponds to tPattern and put it in list.
                // This is kind of a brute force algorithm. however there will be only a few fields
                for( uint Ii = 0; Ii< mParameters.mRefinementPattern.size(); Ii++ )
                {
                    for( uint Ia = 0; Ia< mParameters.mRefinementPattern( Ii ).numel(); Ia++ )
                    {
                        if( tPattern == mParameters.mRefinementPattern( Ii )( Ia ) )
                        {
                            aFieldsForRefinement( Ik ).push_back( mParameters.mFieldNames( Ii ) );
                            aRefinements( Ik )        .push_back( mParameters.mRefinementLevel( Ii )( Ia ) );

                            aMaxRefinementPerPattern( Ik ) = std::max( aMaxRefinementPerPattern( Ik ), mParameters.mRefinementLevel( Ii )( Ia ) );
                        }
                    }
                }
            }
        }



        //--------------------------------------------------------------------------------------------------------------

        void Refinement_Mini_Performer::perform_refinement_old(
                std::shared_ptr<hmr::HMR>          aHMR,
                Cell< std::shared_ptr<Performer> > aPerformers,
                bool                               aSimultaneous)
        {
            moris::sint tMaxNumRefinements = get_max_refinement_level( aPerformers );

            moris::Matrix< DDSMat > tMeshIndices;
            moris::uint             tNumMeshes = 0;

            get_all_refinement_mesh_indices(
                    aPerformers,
                    tMeshIndices,
                    tNumMeshes );

            // sint tRefinementNumber = 0;
            sint tNumPerformers = aPerformers.size();

            for (uint Ii = 0; Ii < tNumMeshes; Ii++)
            {
                sint tMeshIndex = tMeshIndices( Ii );

                for (sint Ij = 0; Ij < tMaxNumRefinements; Ij++)
                {
                    // Create mesh //FIXME
                    std::shared_ptr<hmr::Mesh> tMesh = aHMR->create_mesh( tMeshIndex );

                    uint tLagrangeMeshPattern = tMesh->get_lagrange_mesh_pattern();

                    for (sint Ik = 0; Ik < tNumPerformers; Ik++)
                    {
                        // Queue refinement
                        queue_single_refinement(aHMR, tMesh, aPerformers( Ik ), Ij, tMeshIndex);
                    }

                    //refine
                    // Perform refinement and update index
                    if (true)
                    {
                        aHMR->perform_refinement( tLagrangeMeshPattern );
                        aHMR->update_refinement_pattern( tLagrangeMeshPattern );
                    }

                    //if (tPerformRefinement)
                    //{
                    //    aHMR->perform_refinement_based_on_working_pattern( 0, false );
                    //}
                }
            }

            for (uint Ii = 0; Ii < tNumMeshes; Ii++)
            {
                sint tMeshIndex = tMeshIndices( Ii );

                while ( true )
                {
                    // Create mesh //FIXME
                    std::shared_ptr<hmr::Mesh> tMesh = aHMR->create_mesh( tMeshIndex );

                    uint tLagrangeMeshPattern = tMesh->get_lagrange_mesh_pattern();

                    uint tNumQueuedElements = 0;

                    for (sint Ik = 0; Ik < tNumPerformers; Ik++)
                    {
                        // Queue refinement
                        tNumQueuedElements += queue_low_level_elements_for_refinement(aHMR, tMesh, aPerformers( Ik ), tMeshIndex);
                    }

                    if( tNumQueuedElements == 0 )
                    {
                        break;
                    }

                    //refine
                    // Perform refinement and update index
                    if (true)
                    {
                        aHMR->perform_refinement( tLagrangeMeshPattern );
                        aHMR->update_refinement_pattern( tLagrangeMeshPattern );
                    }

                    //if (tPerformRefinement)
                    //{
                    //    aHMR->perform_refinement_based_on_working_pattern( 0, false );
                    //}
                    break;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Refinement_Mini_Performer::queue_single_refinement(std::shared_ptr<hmr::HMR>  aHMR,
                std::shared_ptr<hmr::Mesh> aMesh,
                std::shared_ptr<Performer> aPerformer,
                sint                       aRefinementNumber,
                sint                       aMeshIndex )
        {
            //uint tLagrangeMeshPattern = aMesh->get_lagrange_mesh_pattern();

            // Loop over fields
            for (uint Ik = 0; Ik < aPerformer->get_num_refinement_fields(); Ik++)
            {
                const moris::Matrix< DDSMat > & tNumRefinements      = aPerformer->get_num_refinements( Ik );
                const moris::Matrix< DDSMat > & tLagrangeMeshIndices = aPerformer->get_refinement_mesh_indices( Ik );

                // loop over tLagrangeMeshIndices // if aMeshIndex put in queue
                for (uint Ii = 0; Ii < tLagrangeMeshIndices.numel(); Ii++)
                {
                    if( tLagrangeMeshIndices( Ii ) == aMeshIndex )
                    {
                        if( tNumRefinements( Ii ) > aRefinementNumber )
                        {
                            // Loop over nodes and get field values
                            Matrix<DDRMat> tFieldValues(aMesh->get_num_nodes(), 1);
                            for (uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++)
                            {
                                tFieldValues(tNodeIndex) = aPerformer->get_field_value(Ik, tNodeIndex, aMesh->get_node_coordinate(tNodeIndex));
                            }

                            // Put elements on queue and set flag for refinement
                            aHMR->based_on_field_put_elements_on_queue(tFieldValues, 0, aPerformer->get_refinement_function_index(Ik, aRefinementNumber));
                        }
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Refinement_Mini_Performer::queue_low_level_elements_for_refinement(
                std::shared_ptr<hmr::HMR>  aHMR,
                std::shared_ptr<hmr::Mesh> aMesh,
                std::shared_ptr<Performer> aPerformer,
                sint                       aMeshIndex )
        {
            uint tNumElements = 0;
            // Loop over fields
            for (uint Ik = 0; Ik < aPerformer->get_num_refinement_fields(); Ik++)
            {
                const moris::Matrix< DDSMat > & tLagrangeMeshIndices = aPerformer->get_refinement_mesh_indices( Ik );

                // loop over tLagrangeMeshIndices // if aMeshIndex put in queue
                for (uint Ii = 0; Ii < tLagrangeMeshIndices.numel(); Ii++)
                {
                    if( tLagrangeMeshIndices( Ii ) == aMeshIndex )
                    {
                        // Loop over nodes and get field values
                        Matrix<DDRMat> tFieldValues(aMesh->get_num_nodes(), 1);
                        for (uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++)
                        {
                            tFieldValues(tNodeIndex) = aPerformer->get_field_value(Ik, tNodeIndex, aMesh->get_node_coordinate(tNodeIndex));
                        }

                        // Put elements on queue and set flag for refinement //FIXME this is untested for a refinement function,
                        tNumElements += aHMR->based_on_field_put_low_level_elements_on_queue(tFieldValues, 0, aPerformer->get_refinement_function_index(Ik, 0));
                    }
                }
            }

            return tNumElements;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::sint Refinement_Mini_Performer::get_max_refinement_level( const Cell< std::shared_ptr<Performer> > & aPerformers )
        {
            sint tMaxNumRefinements = 0;

            for (uint Ia = 0; Ia < aPerformers.size(); Ia++)
            {
                uint tNumFields = aPerformers( Ia )->get_num_refinement_fields();

                // Loop over fields
                for (uint Ik = 0; Ik < tNumFields; Ik++)
                {
                    const moris::Matrix< DDSMat > & tNumRefinements = aPerformers( Ia )->get_num_refinements( Ik );

                    for (uint Ii = 0; Ii < tNumRefinements.numel(); Ii++)
                    {
                        tMaxNumRefinements = std::max( tMaxNumRefinements, tNumRefinements( Ii ));
                    }
                }
            }

            return tMaxNumRefinements;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Refinement_Mini_Performer::get_all_refinement_mesh_indices(
                const Cell< std::shared_ptr<Performer> > & aPerformers,
                moris::Matrix< DDSMat >                  & aAllPatternMap,
                moris::uint                              & aNumPattern )
        {
            moris::Matrix< DDSMat > tCombinedPattern( hmr::gNumberOfPatterns, 1 , -1 );

            aNumPattern = 0;

            for (uint Ia = 0; Ia < aPerformers.size(); Ia++ )
            {
                uint tNumFields = aPerformers( Ia )->get_num_refinement_fields();

                // Loop over fields
                for (uint Ik = 0; Ik < tNumFields; Ik++)
                {
                    const moris::Matrix< DDSMat > & tRefinementMeshIndex = aPerformers( Ia )->get_refinement_mesh_indices( Ik );

                    for (uint Ii = 0; Ii < tRefinementMeshIndex.numel(); Ii++)
                    {
                        if( tCombinedPattern( tRefinementMeshIndex( Ii ) ) == -1)
                        {
                            tCombinedPattern( tRefinementMeshIndex( Ii ) ) = 1;

                            aNumPattern++;
                        }
                    }
                }
            }

            aAllPatternMap.set_size( aNumPattern, 1 , -1 );

            uint tCounter = 0;

            for (uint Ia = 0; Ia < tCombinedPattern.numel(); Ia++ )
            {
                if( tCombinedPattern( Ia ) == 1 )
                {
                    aAllPatternMap( tCounter++ ) = Ia;
                }
            }
        }
    }
}

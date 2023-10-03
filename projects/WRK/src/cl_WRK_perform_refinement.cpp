/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_perform_refinement.cpp
 *
 */

#include "cl_WRK_Performer.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_WRK_perform_refinement.hpp"
#include "HMR_Globals.hpp"
#include "cl_HMR_Background_Element_Base.hpp"
#include <unordered_map>

// Logging package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "cl_MTK_Field.hpp"

#include "cl_Param_List.hpp"

namespace moris
{
    namespace wrk
    {

        Refinement_Mini_Performer::Refinement_Mini_Performer(
                ParameterList&                aParameterlist,
                std::shared_ptr< Library_IO > aLibrary )
                : mLibrary( aLibrary )
        {
            // set field names
            moris::Cell< std::string > tFieldNames;
            string_to_cell(
                    aParameterlist.get< std::string >( "field_names" ),
                    mParameters.mFieldNames );

            if ( !aParameterlist.get< std::string >( "levels_of_refinement" ).empty() )
            {
                // set refinement level
                string_to_cell_mat(
                        aParameterlist.get< std::string >( "levels_of_refinement" ),
                        mParameters.mRefinementLevel );
            }

            // set refinementpattern
            string_to_cell_mat(
                    aParameterlist.get< std::string >( "refinement_pattern" ),
                    mParameters.mRefinementPattern );

            if ( !aParameterlist.get< std::string >( "remeshing_copy_old_pattern_to_pattern" ).empty() )
            {
                // set refinement pattern
                string_to_cell_mat(
                        aParameterlist.get< std::string >( "remeshing_copy_old_pattern_to_pattern" ),
                        mParameters.mRefinemenCopytPatternToPattern );
            }

            string_to_cell(
                    aParameterlist.get< std::string >( "refinement_function_name" ),
                    mParameters.mRefinementFunctionName );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Refinement_Mini_Performer::perform_refinement(
                Cell< std::shared_ptr< mtk::Field > >& aFields,
                std::shared_ptr< hmr::HMR >            aHMR )
        {
            Tracer tTracer( "WRK", "Refinement Mini Performer", "Perform refinement" );
            // create field name to index map
            moris::map< std::string, moris_index > tFieldNameToIndexMap;

            for ( uint Ik = 0; Ik < aFields.size(); Ik++ )
            {
                tFieldNameToIndexMap[ aFields( Ik )->get_label() ] = Ik;
            }

            Cell< moris_index >                                      tPattern;
            moris::Cell< moris::Cell< std::string > >                tFieldNames;
            moris::Cell< moris::Cell< uint > >                       tRefinements;
            moris::Cell< sint >                                      tMaxRefinementPerLevel;
            moris::Cell< moris::Cell< hmr::Refinement_Function_2 > > tRefinementFunctions;

            this->prepare_input_for_refinement(
                    tPattern,
                    tFieldNames,
                    tRefinements,
                    tMaxRefinementPerLevel,
                    tRefinementFunctions );

            if ( ( mParameters.mRefinementFunctionName.size() > 0 ) && mLibrary != nullptr )
            {
                mParameters.mRefinementFunction = mLibrary->load_function< moris::hmr::Refinement_Function >( mParameters.mRefinementFunctionName( 0 ) );
            }

            for ( uint Ik = 0; Ik < tPattern.size(); Ik++ )
            {
                // get pattern
                uint tActivationPattern = tPattern( Ik );

                for ( uint Ii = 0; Ii < tMaxRefinementPerLevel.size(); Ii++ )
                {
                    // Mesh has changed after first refinement and therefore has to be rebuilt
                    // This will only work with analytic fields. When remeshing the L2 id done in the remesh mini performer
                    if ( Ii > 0 )
                    {
                        // Interpolation_Mesh_HMR * tInterpolationMesh = aHMR->create_interpolation_mesh( tOrder, tPattern, tPattern );
                        ////   mtk::Mesh_Pair

                        for ( uint Ia = 0; Ia < tFieldNames( Ik ).size(); Ia++ )
                        {
                            //  mtk::Field * tField = aFields( aFieldNameToIndexMap.find( tFieldNames( Ia ) ) );
                            //  tField->
                        }
                    }

                    for ( uint Ia = 0; Ia < tFieldNames( Ik ).size(); Ia++ )
                    {
                        std::shared_ptr< mtk::Field > tField = aFields( tFieldNameToIndexMap.find( tFieldNames( Ik )( Ia ) ) );

                        uint tOrder = tField->get_lagrange_order();

                        const Matrix< DDRMat >& tFieldValues = tField->get_values();

                        // Put elements on queue and set flag for refinement
                        aHMR->based_on_field_put_elements_on_queue(
                                tFieldValues,
                                tActivationPattern,
                                tOrder,
                                mParameters.mRefinementFunction );
                    }
                }

                aHMR->perform_refinement( tActivationPattern );
                aHMR->update_refinement_pattern( tActivationPattern );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Refinement_Mini_Performer::perform_refinement_2(
                Cell< std::shared_ptr< mtk::Field > >& aFields,
                std::shared_ptr< hmr::HMR >            aHMR )
        {
            Tracer tTracer( "WRK", "Refinement Mini Performer", "Perform advanced refinement" );
            // create field name to index map
            moris::map< std::string, moris_index > tFieldNameToIndexMap;

            for ( uint Ik = 0; Ik < aFields.size(); Ik++ )
            {
                tFieldNameToIndexMap[ aFields( Ik )->get_label() ] = Ik;
            }

            Cell< moris_index >                                      tPattern;
            moris::Cell< moris::Cell< std::string > >                tFieldNames;
            moris::Cell< moris::Cell< uint > >                       tRefinements;
            moris::Cell< sint >                                      tMaxRefinementPerLevel;
            moris::Cell< moris::Cell< hmr::Refinement_Function_2 > > tRefinementFunctions;

            this->prepare_input_for_refinement(
                    tPattern,
                    tFieldNames,
                    tRefinements,
                    tMaxRefinementPerLevel,
                    tRefinementFunctions );

            // get pattern to save
            uint tNumPatternToSave = mParameters.mRefinemenCopytPatternToPattern.size();

            moris::map< sint, sint > tReferencePattern;
            for ( uint Ik = 0; Ik < tNumPatternToSave; Ik++ )
            {
                tReferencePattern[ mParameters.mRefinemenCopytPatternToPattern( Ik )( 0 ) ] =
                        mParameters.mRefinemenCopytPatternToPattern( Ik )( 1 );
            }

            MORIS_ERROR( mParameters.mRefinementFunctionName.size() > 0, "Refinement function names not set" );
            MORIS_ERROR( mLibrary != nullptr, "mLibrary not set" );

            uint tWorkingPattern = aHMR->get_parameters()->get_working_pattern();    // mParameters->get_working_pattern();

            // create a map with ids
            std::unordered_map< moris_index, luint > tMap;

            moris::Cell< hmr::Background_Element_Base* > tBGElements;
            aHMR->get_database()->get_background_mesh()->collect_all_elements( tBGElements );

            for ( uint Ib = 0; Ib < tBGElements.size(); Ib++ )
            {
                auto tHMRId_2 = tBGElements( Ib )->get_hmr_id();

                tMap[ tHMRId_2 ] = Ib;
            }

            for ( uint Ik = 0; Ik < tPattern.size(); Ik++ )
            {

                aHMR->get_database()->get_background_mesh()->reset_pattern( tWorkingPattern );

                // get pattern
                uint tActivationPattern = tPattern( Ik );

                sint tThisRefernecePattern = tReferencePattern.find( tActivationPattern );

                for ( uint Ia = 0; Ia < tFieldNames( Ik ).size(); Ia++ )
                {
                    std::string                   tFieldName = tFieldNames( Ik )( Ia );
                    std::shared_ptr< mtk::Field > tField     = aFields( tFieldNameToIndexMap.find( tFieldName ) );

                    // get interpolation mesh from mesh pair
                    moris::mtk::Mesh* tSourceMesh = tField->get_mesh_pair().get_interpolation_mesh();

                    MORIS_ERROR( tSourceMesh != nullptr, "Source mesh of field %s is nullptr", tFieldName.c_str() );

                    uint tNumElements = tSourceMesh->get_num_elems();

                    for ( uint Ii = 0; Ii < tNumElements; Ii++ )
                    {
                        mtk::Cell& tCell = tSourceMesh->get_mtk_cell( Ii );

                        uint tMaxLevel = MORIS_UINT_MAX;

                        enum hmr::ElementalRefienmentIndicator tRefIndicator =
                                tRefinementFunctions( Ik )( Ia )( &tCell, tField, tActivationPattern, tMaxLevel );

                        hmr::Background_Element_Base* tBGElementOld =
                                reinterpret_cast< hmr::Element* >( tCell.get_base_cell() )->get_background_element();

                        luint tHMRId = tBGElementOld->get_hmr_id();

                        auto tIter = tMap.find( tHMRId );

                        moris_index tIndex = tIter->second;

                        hmr::Background_Element_Base* tBGElementNew = tBGElements( tIndex );

                        while ( !tBGElementNew->is_active( tThisRefernecePattern ) )
                        {
                            // jump to parent
                            tBGElementNew = tBGElementNew->get_parent();
                        }

                        // in the case that the BG level is equal the max level, hold this element
                        if ( tRefIndicator == hmr::ElementalRefienmentIndicator::REFINE && tBGElementNew->get_level() >= tMaxLevel )
                        {
                            tRefIndicator = hmr::ElementalRefienmentIndicator::HOLD;
                        }

                        switch ( tRefIndicator )
                        {
                            case hmr::ElementalRefienmentIndicator::REFINE:
                            {
                                tBGElementNew->set_refined_flag( tWorkingPattern );
                                break;
                            }
                            case hmr::ElementalRefienmentIndicator::HOLD:
                            {
                                if ( tBGElementNew->get_level() > 0 )
                                {
                                    tBGElementNew->get_parent()    //
                                            ->set_refined_flag( tWorkingPattern );
                                }

                                break;
                            }
                            case hmr::ElementalRefienmentIndicator::COARSEN:
                            {
                                if ( tBGElementNew->get_level() > 1 )
                                {
                                    tBGElementNew->get_parent()    //
                                            ->get_parent()         //
                                            ->set_refined_flag( tWorkingPattern );
                                }

                                break;
                            }
                            case hmr::ElementalRefienmentIndicator::DROP:
                            {
                                // do nothing
                                break;
                            }
                            default:
                            {
                                MORIS_ERROR( false,
                                        "Refinement_Mini_Performer::perform_refinement_2 " );
                            }
                        }
                    }
                }

                aHMR->perform_refinement_based_on_working_pattern( tActivationPattern );
                aHMR->update_refinement_pattern( tActivationPattern );

                aHMR->get_database()->get_background_mesh()->update_database();
                aHMR->get_database()->update_bspline_meshes();
                aHMR->get_database()->update_lagrange_meshes();
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Refinement_Mini_Performer::perform_refinement_based_on_working_pattern(
                Cell< std::shared_ptr< mtk::Field > >& aFields,
                std::shared_ptr< hmr::HMR >            aHMR )
        {
            // create field name to index map
            moris::map< std::string, moris_index > tFieldNameToIndexMap;

            for ( uint Ik = 0; Ik < aFields.size(); Ik++ )
            {
                tFieldNameToIndexMap[ aFields( Ik )->get_label() ] = Ik;
            }

            if ( ( mParameters.mRefinementFunctionName.size() > 0 ) && mLibrary != nullptr )
            {
                mParameters.mRefinementFunction = mLibrary->load_function< moris::hmr::Refinement_Function >( mParameters.mRefinementFunctionName( 0 ) );
            }

            MORIS_ERROR( mParameters.mRefinementFunction != nullptr, "Refinement function not set!" );

            for ( uint Ik = 0; Ik < mParameters.mRefinementPattern( 0 ).numel(); Ik++ )
            {
                // get pattern
                uint tActivationPattern = mParameters.mRefinementPattern( 0 )( Ik );

                for ( uint Ia = 0; Ia < mParameters.mFieldNames.size(); Ia++ )
                {
                    std::shared_ptr< mtk::Field > tField = aFields( tFieldNameToIndexMap.find( mParameters.mFieldNames( Ia ) ) );

                    uint tOrder = tField->get_lagrange_order();

                    const Matrix< DDRMat >& tFieldValues = tField->get_values();

                    // Put elements on queue and set flag for refinement
                    aHMR->based_on_field_flag_elements_on_working_pattern(
                            tFieldValues,
                            4,    // tActivationPattern
                            tOrder,
                            mParameters.mRefinementFunction );
                }

                aHMR->perform_refinement_based_on_working_pattern( tActivationPattern );
                aHMR->update_refinement_pattern( tActivationPattern );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        uint
        Refinement_Mini_Performer::perform_refinement_low_level_elements(
                Cell< std::shared_ptr< mtk::Field > >& aFields,
                std::shared_ptr< hmr::HMR >            aHMR )
        {
            Tracer tTracer( "WRK", "Refinement Mini Performer", "Perform refinement of low level elements" );
            uint   tNumQueuedElements = 0;

            // create field name to index map
            moris::map< std::string, moris_index > tFieldNameToIndexMap;

            for ( uint Ik = 0; Ik < aFields.size(); Ik++ )
            {
                tFieldNameToIndexMap[ aFields( Ik )->get_label() ] = Ik;
            }

            Cell< moris_index >                                      tPattern;
            moris::Cell< moris::Cell< std::string > >                tFieldNames;
            moris::Cell< moris::Cell< uint > >                       tRefinements;
            moris::Cell< sint >                                      tMaxRefinementPerLevel;
            moris::Cell< moris::Cell< hmr::Refinement_Function_2 > > tRefinementFunctions;

            this->prepare_input_for_refinement(
                    tPattern,
                    tFieldNames,
                    tRefinements,
                    tMaxRefinementPerLevel,
                    tRefinementFunctions );

            for ( uint Ik = 0; Ik < tPattern.size(); Ik++ )
            {
                // get pattern
                uint tActivationPattern = tPattern( Ik );

                for ( uint Ia = 0; Ia < tFieldNames( Ik ).size(); Ia++ )
                {
                    std::shared_ptr< mtk::Field > tField = aFields( tFieldNameToIndexMap.find( tFieldNames( Ik )( Ia ) ) );

                    uint tOrder = tField->get_lagrange_order();

                    const Matrix< DDRMat >& tFieldValues = tField->get_values();

                    // Put elements on queue and set flag for refinement
                    tNumQueuedElements += aHMR->based_on_field_put_low_level_elements_on_queue(
                            tFieldValues,
                            tActivationPattern,
                            tOrder,
                            mParameters.mRefinementFunction );
                }

                aHMR->perform_refinement( tActivationPattern );
                aHMR->update_refinement_pattern( tActivationPattern );
            }

            return tNumQueuedElements;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Refinement_Mini_Performer::prepare_input_for_refinement(
                Cell< moris_index >&                                      aPatternForRefinement,
                moris::Cell< moris::Cell< std::string > >&                aFieldsForRefinement,
                moris::Cell< moris::Cell< uint > >&                       aRefinements,
                moris::Cell< sint >&                                      aMaxRefinementPerPattern,
                moris::Cell< moris::Cell< hmr::Refinement_Function_2 > >& aRefinementFunctions )
        {

            // produce unique list of pattern which will be refined
            for ( uint Ik = 0; Ik < mParameters.mRefinementPattern.size(); Ik++ )
            {
                for ( uint Ii = 0; Ii < mParameters.mRefinementPattern( Ik ).numel(); Ii++ )
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
            aFieldsForRefinement.resize( tNumberOfRefinementPattern );
            aRefinements.resize( tNumberOfRefinementPattern );
            aRefinementFunctions.resize( tNumberOfRefinementPattern );
            aMaxRefinementPerPattern.resize( tNumberOfRefinementPattern, 0 );

            // create list with field pointers and refinements per pattern
            for ( uint Ik = 0; Ik < tNumberOfRefinementPattern; Ik++ )
            {
                moris_index tPattern = aPatternForRefinement( Ik );

                // loop over all fields and corresponding patterns. Find the pattern which corresponds to tPattern and put it in list.
                // This is kind of a brute force algorithm. however there will be only a few fields
                for ( uint Ii = 0; Ii < mParameters.mRefinementPattern.size(); Ii++ )
                {
                    for ( uint Ia = 0; Ia < mParameters.mRefinementPattern( Ii ).numel(); Ia++ )
                    {
                        if ( tPattern == mParameters.mRefinementPattern( Ii )( Ia ) )
                        {
                            aFieldsForRefinement( Ik ).push_back( mParameters.mFieldNames( Ii ) );
                            if ( mParameters.mRefinementLevel.size() > 0 )
                            {
                                aRefinements( Ik ).push_back( mParameters.mRefinementLevel( Ii )( Ia ) );
                                aMaxRefinementPerPattern( Ik ) = std::max( aMaxRefinementPerPattern( Ik ), mParameters.mRefinementLevel( Ii )( Ia ) );
                            }
                            if ( mParameters.mRefinementFunctionName.size() > 0 )
                            {
                                aRefinementFunctions( Ik ).push_back(
                                        mLibrary->load_function< moris::hmr::Refinement_Function_2 >( mParameters.mRefinementFunctionName( Ii ) ) );
                            }
                        }
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Refinement_Mini_Performer::perform_refinement_old(
                std::shared_ptr< hmr::HMR >          aHMR,
                Cell< std::shared_ptr< Performer > > aPerformers,
                bool                                 aSimultaneous )
        {
            Tracer tTracer( "WRK", "Refinement Mini Performer", "Perform GEN refinement" );

            // get the maximum number of refinements that need to be performed and therefore refinement iterations
            moris::sint tMaxNumRefinements = this->get_max_refinement_level( aPerformers );

            // TODO: get number of what kind of meshes? (Lagrange meshes?, B-spline meshes?, output meshes?, other?)
            moris::Matrix< DDSMat > tMeshIndices;
            uint                    tNumMeshes = 0;
            this->get_all_refinement_mesh_indices(
                    aPerformers,
                    tMeshIndices,
                    tNumMeshes );

            // get number of performers that request refinement
            moris_index tNumPerformers = aPerformers.size();

            // get access to the database and background mesh
            std::shared_ptr< hmr::Database > tDataBase       = aHMR->get_database();
            hmr::Background_Mesh_Base*       tBackgroundMesh = tDataBase->get_background_mesh();

            // TODO: loop over what kind of meshes?
            for ( uint iMesh = 0; iMesh < tNumMeshes; iMesh++ )
            {
                // TODO: what kind of mesh index?
                moris_index tMeshIndex = tMeshIndices( iMesh );

                // iteratively refine meshes
                for ( moris_index iRefinement = 0; iRefinement < tMaxNumRefinements; iRefinement++ )
                {
                    // Create mesh //FIXME
                    std::shared_ptr< hmr::Mesh > tMesh = aHMR->create_mesh( tMeshIndex );

                    uint tLagrangeMeshPattern = tMesh->get_lagrange_mesh_pattern();

                    for ( moris_index iPerformer = 0; iPerformer < tNumPerformers; iPerformer++ )
                    {
                        // Queue refinement
                        this->queue_single_refinement(
                                aHMR,
                                tMesh,
                                aPerformers( iPerformer ),
                                iRefinement,
                                tMeshIndex );
                    }

                    //  perform refinement on the geometric meshes (replacing elements with their children, where necessary)
                    aHMR->perform_refinement( tLagrangeMeshPattern );
                    aHMR->update_refinement_pattern( tLagrangeMeshPattern );

                    // TODO: ?
                    tBackgroundMesh->update_database();
                    tDataBase->update_bspline_meshes();
                    tDataBase->update_lagrange_meshes();

                }    // end for: each refinement step

                // TODO: ?
                tBackgroundMesh->clear_refinement_queue();

            }    // end for: each mesh // TODO: what kind of meshes?

            // refinement loop to ensure that all intersected elements are on same refinement level
            for ( uint iMesh = 0; iMesh < tNumMeshes; iMesh++ )
            {
                moris_index tMeshIndex = tMeshIndices( iMesh );

                bool tRefinedAllLowLevelElements = false;

                // FIXME: ...
                const uint tMaxLowLevelRefinementSteps = 20;

                // set to true for using low level element refinement
                for ( uint iLowLevelRefinement = 0; iLowLevelRefinement < tMaxLowLevelRefinementSteps; ++iLowLevelRefinement )
                {
                    // Create mesh //FIXME
                    std::shared_ptr< hmr::Mesh > tMesh = aHMR->create_mesh( tMeshIndex );

                    uint tLagrangeMeshPattern = tMesh->get_lagrange_mesh_pattern();

                    uint tNumQueuedElements = 0;

                    for ( moris_index iPerformer = 0; iPerformer < tNumPerformers; iPerformer++ )
                    {
                        // Queue refinement
                        tNumQueuedElements += this->queue_low_level_elements_for_refinement(
                                aHMR,
                                tMesh,
                                aPerformers( iPerformer ),
                                tMeshIndex );
                    }

                    // compute number of queued elements across all processors
                    tNumQueuedElements = sum_all( tNumQueuedElements );

                    MORIS_LOG_INFO( "Found %d elements for low level refinement.", tNumQueuedElements );

                    if ( tNumQueuedElements == 0 )
                    {
                        tRefinedAllLowLevelElements = true;
                        break;
                    }

                    // Perform refinement and update indices of elements
                    aHMR->perform_refinement( tLagrangeMeshPattern );
                    aHMR->update_refinement_pattern( tLagrangeMeshPattern );

                    // FIXME should be removed such that loop is continued until all elements are refined

                }    // end for: each lower level refinement

                // check that all low level elements were refined
                MORIS_ERROR(
                        tRefinedAllLowLevelElements,
                        "Refinement_Mini_Performer::perform_refinement_old() - "
                        "could not refine all low level elements in %d steps.",
                        tMaxLowLevelRefinementSteps );

            }    // end for: each mesh // TODO: what kind of meshes?

            // TODO: ?
            tDataBase->update_bspline_meshes();
            tDataBase->update_lagrange_meshes();

        }    // end function:

        //--------------------------------------------------------------------------------------------------------------

        void
        Refinement_Mini_Performer::queue_single_refinement(
                std::shared_ptr< hmr::HMR >  aHMR,
                std::shared_ptr< hmr::Mesh > aMesh,
                std::shared_ptr< Performer > aPerformer,
                sint                         aRefinementNumber,
                sint                         aMeshIndex )
        {
            // uint tLagrangeMeshPattern = aMesh->get_lagrange_mesh_pattern();

            // Loop over fields
            for ( uint iField = 0; iField < aPerformer->get_num_refinement_fields(); iField++ )
            {
                const moris::Matrix< DDSMat >& tNumRefinements      = aPerformer->get_num_refinements( iField );
                const moris::Matrix< DDSMat >& tLagrangeMeshIndices = aPerformer->get_refinement_mesh_indices( iField );

                // loop over tLagrangeMeshIndices // if aMeshIndex put in queue
                for ( uint iLagMesh = 0; iLagMesh < tLagrangeMeshIndices.numel(); iLagMesh++ )
                {
                    if ( tLagrangeMeshIndices( iLagMesh ) == aMeshIndex )
                    {
                        if ( tNumRefinements( iLagMesh ) > aRefinementNumber )
                        {
                            // Loop over nodes and get field values
                            Matrix< DDRMat > tFieldValues( aMesh->get_num_nodes(), 1 );
                            for ( uint tNodeIndex = 0; tNodeIndex < aMesh->get_num_nodes(); tNodeIndex++ )
                            {
                                Matrix< DDRMat > tNodeCoordinates = aMesh->get_node_coordinate( tNodeIndex );
                                tFieldValues( tNodeIndex )        = aPerformer->get_field_value( iField, tNodeIndex, tNodeCoordinates );
                            }

                            // Put elements on queue and set flag for refinement
                            sint tRefinementFunctionIndex = aPerformer->get_refinement_function_index( iField, aRefinementNumber );
                            aHMR->based_on_field_put_elements_on_queue( tFieldValues, aMeshIndex, tRefinementFunctionIndex );
                        }
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        uint
        Refinement_Mini_Performer::queue_low_level_elements_for_refinement(
                std::shared_ptr< hmr::HMR >  aHMR,
                std::shared_ptr< hmr::Mesh > aMesh,
                std::shared_ptr< Performer > aPerformer,
                sint                         aMeshIndex )
        {
            uint tNumElements = 0;

            // Loop over fields
            for ( uint iField = 0; iField < aPerformer->get_num_refinement_fields(); iField++ )
            {
                const moris::Matrix< DDSMat >& tLagrangeMeshIndices = aPerformer->get_refinement_mesh_indices( iField );

                // loop over tLagrangeMesh indices // if aMeshIndex put in queue
                for ( uint iLagMesh = 0; iLagMesh < tLagrangeMeshIndices.numel(); iLagMesh++ )
                {
                    if ( tLagrangeMeshIndices( iLagMesh ) == aMeshIndex )
                    {
                        // get the number of nodes on the mesh
                        uint tNumNodes = aMesh->get_num_nodes();

                        // Loop over nodes and get field values
                        Matrix< DDRMat > tFieldValues( tNumNodes, 1 );

                        for ( uint tNodeIndex = 0; tNodeIndex < tNumNodes; tNodeIndex++ )
                        {
                            Matrix< DDRMat > tNodeCoordinates = aMesh->get_node_coordinate( tNodeIndex );
                            tFieldValues( tNodeIndex ) =
                                    aPerformer->get_field_value(
                                            iField,
                                            tNodeIndex,
                                            tNodeCoordinates );
                        }

                        // Put elements on queue and set flag for refinement // FIXME: this is untested for a refinement function,
                        sint tRefineFuncIndex = aPerformer->get_refinement_function_index( iField, 0 );
                        tNumElements += aHMR->based_on_field_put_low_level_elements_on_queue(
                                tFieldValues,
                                aMeshIndex,
                                tRefineFuncIndex );
                    }
                }
            }

            return tNumElements;
        }

        //--------------------------------------------------------------------------------------------------------------

        moris::sint
        Refinement_Mini_Performer::get_max_refinement_level( const Cell< std::shared_ptr< Performer > >& aPerformers )
        {
            sint tMaxNumRefinements = 0;

            for ( uint iPerformer = 0; iPerformer < aPerformers.size(); iPerformer++ )
            {
                uint tNumFields = aPerformers( iPerformer )->get_num_refinement_fields();

                // Loop over fields
                for ( uint iField = 0; iField < tNumFields; iField++ )
                {
                    const moris::Matrix< DDSMat >& tNumRefinements = aPerformers( iPerformer )->get_num_refinements( iField );

                    for ( uint iRefinement = 0; iRefinement < tNumRefinements.numel(); iRefinement++ )
                    {
                        tMaxNumRefinements = std::max( tMaxNumRefinements, tNumRefinements( iRefinement ) );
                    }
                }
            }

            return tMaxNumRefinements;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Refinement_Mini_Performer::get_all_refinement_mesh_indices(
                const Cell< std::shared_ptr< Performer > >& aPerformers,
                moris::Matrix< DDSMat >&                    aAllPatternMap,
                moris::uint&                                aNumPattern )
        {
            moris::Matrix< DDSMat > tCombinedPattern( hmr::gNumberOfPatterns, 1, -1 );

            aNumPattern = 0;

            for ( uint iPerformer = 0; iPerformer < aPerformers.size(); iPerformer++ )
            {
                uint tNumFields = aPerformers( iPerformer )->get_num_refinement_fields();

                // Loop over fields
                for ( uint iRefineField = 0; iRefineField < tNumFields; iRefineField++ )
                {
                    const moris::Matrix< DDSMat >& tRefinementMeshIndices = aPerformers( iPerformer )->get_refinement_mesh_indices( iRefineField );

                    for ( uint iMesh = 0; iMesh < tRefinementMeshIndices.numel(); iMesh++ )
                    {
                        sint tRefinementMeshIndex = tRefinementMeshIndices( iMesh );

                        if ( tCombinedPattern( tRefinementMeshIndex ) == -1 )
                        {
                            tCombinedPattern( tRefinementMeshIndex ) = 1;

                            aNumPattern++;
                        }
                    }
                }
            }

            aAllPatternMap.set_size( aNumPattern, 1, -1 );

            uint tCounter = 0;

            for ( uint iPattern = 0; iPattern < tCombinedPattern.numel(); iPattern++ )
            {
                if ( tCombinedPattern( iPattern ) == 1 )
                {
                    aAllPatternMap( tCounter++ ) = iPattern;
                }
            }

        }    // end function: Refinement_Mini_Performer::get_all_refinement_mesh_indices()

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace wrk
}    // namespace moris

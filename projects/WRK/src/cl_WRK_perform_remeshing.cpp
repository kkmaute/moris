
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_File.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "HMR_Globals.hpp"

#include "cl_MTK_Field.hpp"
#include "cl_MTK_Field_Discrete.hpp"
#include "cl_MTK_Mapper.hpp"

#include "cl_WRK_perform_remeshing.hpp"
#include "cl_WRK_perform_refinement.hpp"
#include "cl_Param_List.hpp"

#include "fn_PRM_MORIS_GENERAL_Parameters.hpp"

#include <memory>


namespace moris
{
    namespace wrk
    {

        Remeshing_Mini_Performer::Remeshing_Mini_Performer(
                ParameterList                 & aParameterlist,
                std::shared_ptr< Library_IO >   aLibrary )
        : mLibrary( aLibrary )
        {
            mParameters.mMode = aParameterlist.get< std::string >( "mode" );

            moris::map< std::string, moris::uint > tModeMap;
            tModeMap["ab_initio"] = 0;
            tModeMap["former"] = 1;

            MORIS_ERROR( tModeMap.key_exists( mParameters.mMode ),
                    "Remeshing_Mini_Performer::Remeshing_Mini_Performer(), Mode name does not exist" );

            mParameters.mModeIndex = tModeMap.find( mParameters.mMode );

            if( mParameters.mModeIndex == 0 )
            {
                // Mode 0
                moris::Cell< std::string > tFieldNames;
                string_to_cell(
                        aParameterlist.get< std::string >( "remeshing_field_names" ),
                        tFieldNames );

                // set refinement level
                Cell< Matrix< DDSMat > > tRefinementLevel;
                string_to_cell_mat(
                        aParameterlist.get< std::string >( "remeshing_levels_of_refinement" ),
                        tRefinementLevel );

                // set refinementpattern
                Cell< Matrix< DDSMat > >  tRefinementPattern;
                string_to_cell_mat(
                        aParameterlist.get< std::string >( "remeshing_refinement_pattern" ),
                        tRefinementPattern );

                mParameters.mRefinementsFieldNames_0 = tFieldNames;
                mParameters.mRefinementsMode_0       = tRefinementLevel;
                mParameters.mRefinementPatternMode_0 = tRefinementPattern;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Remeshing_Mini_Performer::perform_remeshing(
                moris::Cell< std::shared_ptr< mtk::Field > >          aSourceFields,
                moris::Cell< std::shared_ptr< hmr::HMR > >          & aHMRPerformers,
                moris::Cell< std::shared_ptr< mtk::Mesh_Manager > > & aMTKPerformer,
                moris::Cell< std::shared_ptr< mtk::Field > >        & aNewFields )
        {
            // find first discrete field. // assumption: all fields are based on the same interpolation mesh
            // If no discrete field is found use first with index 0 TODO
            sint tFirstDiscreteFieldIndex = 0;
            for(uint Ik = 0; Ik < aSourceFields.size(); Ik++ )
            {
                if( aSourceFields( Ik )->get_field_is_discrete() )
                {
                    tFirstDiscreteFieldIndex = Ik;
                    break;
                }
            }

            // get mesh pair from discrete field. We assume all fields are based on the same mesh
            mtk::Mesh_Pair tMeshPairIn = aSourceFields( tFirstDiscreteFieldIndex )->get_mesh_pair();

            // get interpolation mesh from mesh pair
            moris::mtk::Mesh * tSourceMesh = tMeshPairIn.get_interpolation_mesh();

            // check if interpolation mesh is an HMR mesh
            MORIS_ERROR( tSourceMesh->get_mesh_type() == MeshType::HMR,
                    "Mapper::map_input_field_to_output_field() Source mesh is not and HMR mesh" );

            // get Lagrange and Discretization order
            uint tSourceLagrangeOrder     = aSourceFields( tFirstDiscreteFieldIndex )->get_lagrange_order();
            uint tDiscretizationOrder     = aSourceFields( tFirstDiscreteFieldIndex )->get_discretization_order();
            uint tDiscretizationMeshIndex = aSourceFields( tFirstDiscreteFieldIndex )->get_discretization_mesh_index();

            // extract pattern from mesh on which this field id based on
            // both Lagrange and discretization order if they are not the same
            Matrix< DDLUMat > tElementCounterPerLevelAndPattern;
            moris::Cell< Matrix< DDLUMat > > tElementPerPattern;

            hmr::File tFile;
            tFile.save_refinement_pattern(
                    tSourceMesh->get_HMR_lagrange_mesh(),
                    tDiscretizationMeshIndex,
                    tElementCounterPerLevelAndPattern,
                    tElementPerPattern );

            // copy HMR parameters from old HMR performer
            hmr::Parameters * tParameters = aHMRPerformers( 0 )->get_parameters();

            // unset parameter owning flag. Input HMR does not own this parameteropject anymore
            aHMRPerformers( 0 )->get_database()->unset_parameter_owning_flag();

            //FIXME adjust parameters here. eg no initial refinement

            // build new HMR performer with copied parameters
            std::shared_ptr< hmr::HMR > tHMRPerformerNew = std::make_shared< hmr::HMR >( tParameters );

            // refine pattern 5 and 6 with given pattern
            tHMRPerformerNew->get_database()->load_refinement_pattern(
                    tElementCounterPerLevelAndPattern,
                    tElementPerPattern);

            uint tNumPattern = tElementPerPattern.size();

            hmr::Interpolation_Mesh_HMR * tOldInterpolationMesh = nullptr;

            if( tNumPattern == 1 )
            {
                // create mesh based on pattern 5 and 6
                tOldInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
                        tHMRPerformerNew->get_database(),
                        tSourceLagrangeOrder,
                        5,
                        tDiscretizationOrder,
                        5); // order, Lagrange pattern, bspline pattern
            }
            else
            {
                // create mesh based on pattern 5 and 6
                tOldInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
                        tHMRPerformerNew->get_database(),
                        tSourceLagrangeOrder,
                        5,
                        tDiscretizationOrder,
                        6); // order, Lagrange pattern, bspline pattern
            }

            // Create  mesh pair
            mtk::Mesh_Pair tMeshPairOld(tOldInterpolationMesh, nullptr, true);

            // create list of fields
            Cell< std::shared_ptr< mtk::Field > > tOldFields;

            // transfer nodal values from old-HMR-Fields to new-HMR-fields.
            // these fields are based on the same lagrange mesh
            this->map_fields(
                    aSourceFields,
                    tOldFields,
                    tMeshPairOld,
                    0,
                    false ); // discretization mesh index is zero because mesh has only on discretization ( tOldInterpolationMesh )

            // use new-HMR-fields to perform refinement
            this->perform_refinement( tHMRPerformerNew, tOldFields );

            // overwrite old HMR performer with new HMR performer
            aHMRPerformers( 0 ) = tHMRPerformerNew;

            // create MTK performer - will be used for HMR mesh
            aMTKPerformer( 0 ) =std::make_shared< mtk::Mesh_Manager >();

            // Set performer to HMR
            aHMRPerformers( 0 )->set_performer( aMTKPerformer( 0 ) );

            // Set HMR performer to MTK performer
            aMTKPerformer( 0 )->set_performer( aHMRPerformers( 0 ) );

            //------------------------------------------------------------------------------

            //FIXME //FIXME
            uint tSourceBSplinePattern = tSourceMesh->get_HMR_lagrange_mesh() ->get_bspline_pattern( tDiscretizationMeshIndex );
            //uint tPattern =0;

            // create mesh for this pattern
            hmr::Interpolation_Mesh_HMR * tInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
                    tHMRPerformerNew->get_database(),
                    tSourceLagrangeOrder,
                    tSourceBSplinePattern,
                    tDiscretizationOrder,
                    tSourceBSplinePattern);

            mtk::Mesh_Pair tMeshPair(tInterpolationMesh, nullptr, true);

            this->map_fields(
                    tOldFields,
                    aNewFields,
                    tMeshPair,
                    0,
                    true ); //FIXME =  discretization meshindex

            //------------------------------------------------------------------------------

            // HMR finalize
            aHMRPerformers( 0 )->perform();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Remeshing_Mini_Performer::perform_refinement(
                std::shared_ptr< hmr::HMR >           aHMRPerformer,
                Cell< std::shared_ptr< mtk::Field > > aSourceFields )
        {
            // switch based on mode index
            switch ( mParameters.mModeIndex )
            {
                case 0 :
                {
                    this->perform_refinement_mode_0( aHMRPerformer, aSourceFields );

                    break;
                }
                case 1 :
                {
                    this->perform_refinement_mode_1( aHMRPerformer, aSourceFields );

                    break;
                }
                default:
                    MORIS_ERROR( false, "Remeshing_Mini_Performer::perform_refinement() - Refinement mode does not exist" );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Remeshing_Mini_Performer::perform_refinement_mode_0(
                std::shared_ptr< hmr::HMR >             aHMRPerformer,
                Cell< std::shared_ptr< mtk::Field > > & aSourceFields )
        {
            // perform initial uniform refinement
            aHMRPerformer->perform_initial_refinement();

            sint tFirstDiscreteFieldIndex = 0;
            for(uint Ik = 0; Ik < aSourceFields.size(); Ik++ )
            {
                if( aSourceFields( Ik )->get_field_is_discrete() )
                {
                    tFirstDiscreteFieldIndex = Ik;
                    break;
                }
            }

            // get input field order
            uint tLagrangeOrder = aSourceFields( tFirstDiscreteFieldIndex )->get_lagrange_order();
            uint tDiscretizationOrder = aSourceFields( tFirstDiscreteFieldIndex )->get_discretization_order();

            std::shared_ptr< hmr::Database > tHMRDatabase = aHMRPerformer->get_database();

            Cell< moris_index >                       tRefinementPattern;
            moris::Cell< moris::Cell< std::string > > tFieldNames;
            moris::Cell< moris::Cell< uint > >        tRefinements;
            moris::Cell< sint >                       tMaxRefinementPerLevel;

            this->prepare_input_for_refinement(
                    tRefinementPattern,
                    tFieldNames,
                    tRefinements,
                    tMaxRefinementPerLevel );

            uint tNumPattern = tRefinementPattern.size();

            // loop over pattern
            for( uint Ik = 0; Ik < tNumPattern; Ik++ )
            {
                // get pattern
                uint tPattern = tRefinementPattern( Ik );

                for( sint Ii = 0; Ii < tMaxRefinementPerLevel( Ik ); Ii++ )
                {
                    // create mesh for this pattern
                    hmr::Interpolation_Mesh_HMR * tInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
                            tHMRDatabase,
                            tLagrangeOrder,
                            tPattern,
                            tDiscretizationOrder,
                            tPattern);

                    mtk::Mesh_Pair tMeshPair(tInterpolationMesh, nullptr, true);

                    Cell< std::shared_ptr< mtk::Field > > aTargetFields;

                    this->map_fields(
                            aSourceFields,
                            aTargetFields,
                            tMeshPair,
                            tPattern,
                            true ); //FIXME tPattern = DiscretizationMeshiondex

                    // create refinement parameter list
                    moris::ParameterList tRefinementParameterlist;
                    this->create_refinement_input_list( tRefinementParameterlist, tPattern );

                    // create refinement mini performer and perform refinement
                    wrk::Refinement_Mini_Performer tRefinementMiniPerformer( tRefinementParameterlist );
                    tRefinementMiniPerformer.perform_refinement( aTargetFields, aHMRPerformer );

                    tHMRDatabase->get_background_mesh()->update_database();
                    tHMRDatabase->update_bspline_meshes();
                    tHMRDatabase->update_lagrange_meshes();
                }

                {
                    uint tCounter = 0;
                    while( true )
                    {
                        // create mesh for this pattern
                        hmr::Interpolation_Mesh_HMR * tInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
                                tHMRDatabase,
                                tLagrangeOrder,
                                tPattern,
                                tDiscretizationOrder,
                                tPattern);

                        mtk::Mesh_Pair tMeshPair(tInterpolationMesh, nullptr, true);

                        Cell< std::shared_ptr< mtk::Field > > aTargetFields;

                        this->map_fields(
                                aSourceFields,
                                aTargetFields,
                                tMeshPair,
                                tPattern,
                                true ); //FIXME tPattern = DiscretizationMeshiondex

                        // create refinement parameter list
                        moris::ParameterList tRefinementParameterlist;
                        this->create_refinement_input_list( tRefinementParameterlist, tPattern );

                        // create refinement mini performer and perform refinement
                        wrk::Refinement_Mini_Performer tRefinementMiniPerformer( tRefinementParameterlist );
                        uint tRefinedElements = tRefinementMiniPerformer.perform_refinement_low_level_elements( aTargetFields, aHMRPerformer );

                        uint tSumRefEle =sum_all( tRefinedElements );

                        if( tSumRefEle == 0)
                        {
                            break;
                        }
                        tHMRDatabase->get_background_mesh()->update_database();
                        tHMRDatabase->update_bspline_meshes();
                        tHMRDatabase->update_lagrange_meshes();

                        tCounter++;
                    }
                }
            }
        }


        //--------------------------------------------------------------------------------------------------------------
        void Remeshing_Mini_Performer::perform_refinement_mode_1(
                std::shared_ptr< hmr::HMR >           aHMRPerformer,
                Cell< std::shared_ptr< mtk::Field > > aSourceFields )
        {

        }
        //--------------------------------------------------------------------------------------------------------------

        void Remeshing_Mini_Performer::map_fields(
                Cell< std::shared_ptr< mtk::Field > > & aSourceFields,
                Cell< std::shared_ptr< mtk::Field > > & aTargetFields,
                mtk::Mesh_Pair                        & aMeshPair,
                uint                                    aDiscretizationMeshIndex,
                bool                                    aMapFields)
        {
            uint tNumFields = aSourceFields.size();

            aTargetFields.resize( tNumFields, nullptr );

            for( uint If = 0; If< tNumFields; If++ )
            {
                if( aSourceFields( If )->get_field_is_discrete() )
                {
                    // create field object for this mesh
                    aTargetFields( If )= std::make_shared< mtk::Field_Discrete >( aMeshPair, aDiscretizationMeshIndex );
                    aTargetFields( If )->set_label( aSourceFields( If )->get_label() );

                    if( aMapFields )
                    {
                        // create mapper and map input field to new field
                        mtk::Mapper tMapper;
                        tMapper.map_input_field_to_output_field( aSourceFields( If ).get(), aTargetFields( If ).get() );
                    }
                    else
                    {
                        aTargetFields( If )->unlock_field();
                        //aTargetFields( If )->set_nodal_values( aSourceFields( If )->get_nodal_values() );
                        aTargetFields( If )->set_coefficients( aSourceFields( If )->get_coefficients() );
                        aTargetFields( If )->compute_nodal_values();
                    }
                }
                else
                {
                    aTargetFields( If ) = aSourceFields( If );
                    aTargetFields( If )->unlock_field();
                    aTargetFields( If )->set_mesh_pair( aMeshPair );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Remeshing_Mini_Performer::prepare_input_for_refinement(
                Cell< moris_index >                       & aPatternForRefinement,
                moris::Cell< moris::Cell< std::string > > & aFieldsForRefinement,
                moris::Cell< moris::Cell< uint > >        & aRefinements,
                moris::Cell< sint >                       & aMaxRefinementPerPattern )
        {
            // produce unique list of pattern which will be refined
            for( uint Ik = 0; Ik< mParameters.mRefinementPatternMode_0.size(); Ik++ )
            {
                for( uint Ii = 0; Ii< mParameters.mRefinementPatternMode_0( Ik ).numel(); Ii++ )
                {
                    aPatternForRefinement.push_back( mParameters.mRefinementPatternMode_0( Ik )( Ii ) );
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
            aMaxRefinementPerPattern.resize( tNumberOfRefinementPattern, MORIS_SINT_MAX );

            // create list with field pointers and refinements per pattern
            for( uint Ik = 0; Ik< tNumberOfRefinementPattern; Ik++ )
            {
                moris_index tPattern = aPatternForRefinement( Ik );

                // loop over all fields and corresponding patterns. Find the pattern which corresponds to tPattern and put it in list.
                // This is kind of a brute force algorithm. however there will be only a few fields
                for( uint Ii = 0; Ii< mParameters.mRefinementPatternMode_0.size(); Ii++ )
                {
                    for( uint Ia = 0; Ia< mParameters.mRefinementPatternMode_0( Ii ).numel(); Ia++ )
                    {
                        if( tPattern == mParameters.mRefinementPatternMode_0( Ii )( Ia ) )
                        {
                            aFieldsForRefinement( Ik ).push_back( mParameters.mRefinementsFieldNames_0( Ii ) );

                            // aRefinements are not use tight now but implemented for future use
                            aRefinements( Ik )        .push_back( mParameters.mRefinementsMode_0( Ii )( Ia ) );

                            sint tRefPatt = mParameters.mRefinementsMode_0( Ii )( Ia );

                            MORIS_ERROR( tRefPatt == aMaxRefinementPerPattern( Ik ) || aMaxRefinementPerPattern( Ik ) == MORIS_SINT_MAX,
                                    "prepare_input_for_refinement(), This implementation is limited to one refinement level per pattern."
                                    "It can be extended if needed." );

                            aMaxRefinementPerPattern( Ik ) = tRefPatt;
                        }
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Remeshing_Mini_Performer::create_refinement_input_list(
                moris::ParameterList & aRefinementParameterlist,
                uint                   aPattern )
        {
            std::string tFieldNames = "";
            std::string tPattern    = "";
            std::string tRefinement = "";

            for ( uint Ik = 0; Ik < mParameters.mRefinementsFieldNames_0.size(); Ik++ )
            {
                tFieldNames = tFieldNames + mParameters.mRefinementsFieldNames_0( Ik ) + ",";
                tPattern    = tPattern + "1;";
                tRefinement = tRefinement + ios::stringify( aPattern ) + ";";
            }

            tFieldNames.pop_back();
            tPattern.pop_back();
            tRefinement.pop_back();

            prm::create_refinement_parameterlist( aRefinementParameterlist );
            aRefinementParameterlist.set( "field_names" , tFieldNames );
            aRefinementParameterlist.set( "levels_of_refinement" , tPattern );
            aRefinementParameterlist.set( "refinement_pattern" , tRefinement );
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}

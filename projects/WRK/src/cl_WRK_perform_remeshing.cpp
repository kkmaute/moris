
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
                moris::Cell< std::shared_ptr< mtk::Field > > const  & aSourceFields,
                moris::Cell< std::shared_ptr< hmr::HMR > >          & aHMRPerformers,
                moris::Cell< std::shared_ptr< mtk::Mesh_Manager > > & aMTKPerformer,
                moris::Cell< std::shared_ptr< mtk::Field > >        & aNewFields )
        {
            mtk::Mesh_Pair tMeshPairIn = aSourceFields( 1 )->get_mesh_pair();

            moris::mtk::Mesh * tSourceMesh = tMeshPairIn.get_interpolation_mesh();

            MORIS_ERROR( tSourceMesh->get_mesh_type() == MeshType::HMR,
                    "Mapper::map_input_field_to_output_field() Source mesh is not and HMR mesh" );

            std::shared_ptr< hmr::Database > tHMRDatabase = tSourceMesh->get_HMR_database();

            // get Lagrange and Discretization orer
            uint tSourceLagrangeOrder = aSourceFields( 1 )->get_lagrange_order();
            uint tDiscretizationOrder = aSourceFields( 1 )->get_discretization_order();

            // extract pattern from mesh on which this field id based on
            // both Lagrange and discretization order if they are not the same
            Matrix< DDLUMat > aElementCounterPerLevelAndPattern;
            moris::Cell< Matrix< DDLUMat > > aElementPerPattern;

            hmr::File tFile;
            tFile.save_refinement_pattern(
                    tSourceMesh->get_HMR_lagrange_mesh(),
                    aSourceFields( 1 )->get_discretization_mesh_index(),
                    aElementCounterPerLevelAndPattern,
                    aElementPerPattern );

            // copy HMR parameters from old HMR performer
            hmr::Parameters * tParameters = aHMRPerformers( 0 )->get_parameters();

            aHMRPerformers( 0 )->get_database()->unset_parameter_owning_flag();

            //FIXME adjust parameters here. eg no initial refinement

            // build new HMR performer with copied parameters
            std::shared_ptr< hmr::HMR > tHMRPerformerNew =
                    std::make_shared< hmr::HMR >( tParameters );

            // FIXME put this line into HMR / all of this in function map discrete fields.
            std::shared_ptr< hmr::Database > tHMRDatabaseNew = tHMRPerformerNew->get_database();

            // refine pattern 5 and 6 with given pattern
            tHMRDatabaseNew->load_refinement_pattern(
                    aElementCounterPerLevelAndPattern,
                    aElementPerPattern);

            // create mesh based on pattern 5 and 6
            hmr::Interpolation_Mesh_HMR * tOldInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
                    tHMRDatabaseNew,
                    tSourceLagrangeOrder,
                    5,
                    tDiscretizationOrder,
                    5); // order, Lagrange pattern, bspline pattern

            hmr::Integration_Mesh_HMR* tOldIntegrationMesh = new hmr::Integration_Mesh_HMR(
                    tSourceLagrangeOrder,
                    5,
                    tOldInterpolationMesh);

            // Create  mesh pair
            mtk::Mesh_Pair tMeshPairOld(tOldInterpolationMesh, tOldIntegrationMesh, true);

            // build field with mesh
            std::shared_ptr< mtk::Field > tFieldOld = std::make_shared< mtk::Field_Discrete >( tMeshPairOld );

            // copy values from input mesh to New/Old mesh ( this mesh is build based on the new HMR performer )
            tFieldOld->unlock_field();
            tFieldOld->set_nodal_values( aSourceFields( 1 )->get_nodal_values() );
            tFieldOld->set_label( aSourceFields( 1 )->get_label() );

            Cell< std::shared_ptr< mtk::Field > > tOldFields( 2 );
            tOldFields( 0 ) = aSourceFields( 0 );
            tOldFields( 1 ) = tFieldOld;


            this->perform_refinement( tHMRPerformerNew, tOldFields );

            aHMRPerformers( 0 ) = tHMRPerformerNew;

            // create MTK performer - will be used for HMR mesh
            aMTKPerformer( 0 ) =std::make_shared< mtk::Mesh_Manager >();

            // Set performer to HMR
            aHMRPerformers( 0 )->set_performer( aMTKPerformer( 0 ) );

            // Set HMR performer to MTK performer
            aMTKPerformer( 0 )->set_performer( aHMRPerformers( 0 ) );

            mtk::Field * tFieldOnPattern = nullptr;

            this->map_field(
                    aHMRPerformers( 0 ),
                    tFieldOld.get(),
                    tFieldOnPattern );

            // HMR finalize
            aHMRPerformers( 0 )->perform();

            // build field with mesh
            std::shared_ptr< mtk::Field > tFieldNew =
                    std::make_shared< mtk::Field_Discrete >( aMTKPerformer( 0 )->get_mesh_pair(0), 0 );
            tFieldNew->set_label( "Level_Set_Field" );

            tFieldNew->unlock_field();
            tFieldNew->set_coefficients( tFieldOnPattern->get_coefficients());
            //tFieldNew->unlock_field();
            //tFieldNew->set_nodal_values( tFieldOnPattern->get_nodal_values());


            //tFieldNew->save_field_to_exodus( "FieldNew.exo");

            aNewFields.resize( 1, tFieldNew );
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
                std::shared_ptr< hmr::HMR >           aHMRPerformer,
                Cell< std::shared_ptr< mtk::Field > > aSourceFields )
        {
            // perform initial uniform refinement
            aHMRPerformer->perform_initial_refinement();

            // get input field order
            uint tLagrangeOrder = aSourceFields( 1 )->get_lagrange_order();
            uint tDiscretizationOrder = aSourceFields( 1 )->get_discretization_order();

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

                    hmr::Integration_Mesh_HMR* tIntegrationMesh = new hmr::Integration_Mesh_HMR(
                            tLagrangeOrder,
                            tPattern,
                            tInterpolationMesh);

                    mtk::Mesh_Pair tMeshPair(tInterpolationMesh, tIntegrationMesh, true);

                    uint tNumFields = aSourceFields.size();

                    Cell< std::shared_ptr< mtk::Field > > tFields( tNumFields, nullptr );

                    for( uint If = 0; If< aSourceFields.size(); If++ )
                    {
                        if( aSourceFields( If )->get_field_is_discrete() )
                        {
                            // create field object for this mesh
                            tFields( If )= std::make_shared< mtk::Field_Discrete >( tMeshPair, tPattern );
                            tFields( If )->set_label( aSourceFields( If )->get_label() );

                            // create mapper and map input field to new field
                            mtk::Mapper tMapper;
                            tMapper.map_input_field_to_output_field( aSourceFields( If ).get(), tFields( If ).get() );
                        }
                        else
                        {
                            tFields( If ) = aSourceFields( If );
                            tFields( If )->unlock_field();
                            tFields( If )->set_mesh_pair( tMeshPair );
                        }
                    }

                    // create refinement parameter list
                    moris::ParameterList tRefinementParameterlist;
                    this->create_refinement_input_list( tRefinementParameterlist, tPattern );

                    // create refinement mini performer and perform refinement
                    wrk::Refinement_Mini_Performer tRefinementMiniPerformer( tRefinementParameterlist );
                    tRefinementMiniPerformer.perform_refinement( tFields, aHMRPerformer );

                    tHMRDatabase->get_background_mesh()->update_database();
                    tHMRDatabase->update_bspline_meshes();
                    tHMRDatabase->update_lagrange_meshes();
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
        void Remeshing_Mini_Performer::map_field(
                std::shared_ptr< hmr::HMR >     aHMRPerformer,
                mtk::Field                  *   aSourceField,
                mtk::Field                  * & aTargetField )
        {
            std::shared_ptr< hmr::Database > tHMRDatabase = aHMRPerformer->get_database();

            uint tLagrangeOrder = 1;                  //TODO
            uint tDiscretizationOrder = 1;            //TODO
            uint tPattern =0;

            // create mesh for this pattern
            hmr::Interpolation_Mesh_HMR * tInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
                    tHMRDatabase,
                    tLagrangeOrder,
                    tPattern,
                    tDiscretizationOrder,
                    tPattern);

            hmr::Integration_Mesh_HMR* tIntegrationMesh = new hmr::Integration_Mesh_HMR(
                    tLagrangeOrder,
                    tPattern,
                    tInterpolationMesh);

            mtk::Mesh_Pair tMeshPair(tInterpolationMesh, tIntegrationMesh, true);

            // create field object for this mesh
            mtk::Field_Discrete *  tFieldOnPattern = new mtk::Field_Discrete( tMeshPair, tPattern );
            tFieldOnPattern->set_label( "Level_Set_Field" );

            mtk::Mapper tMapper;
            // create mapper and map input field to new field
            tMapper.map_input_field_to_output_field( aSourceField, tFieldOnPattern );

            //tFieldOnPattern->save_field_to_exodus( "Field1111.exo");

            aTargetField = tFieldOnPattern;
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


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
                Matrix< DDUMat > tRefinementLevel;
                string_to_mat(
                        aParameterlist.get< std::string >( "remeshing_levels_of_refinement" ),
                        tRefinementLevel );

                Matrix< DDUMat > tRefinementPattern;
                string_to_mat(
                        aParameterlist.get< std::string >( "remeshing_refinement_pattern" ),
                        tRefinementPattern );

                mParameters.mRefinementsMode_0 = tRefinementLevel;
                mParameters.mRefinementPatternMode_0 = tRefinementPattern;
            }



        }

        //--------------------------------------------------------------------------------------------------------------

        void Remeshing_Mini_Performer::perform_remeshing(
                mtk::Field                                 * aSourceField,
                moris::Cell< std::shared_ptr< hmr::HMR > > & aHMRPerformers )
        {
            mtk::Mesh_Pair * tMeshPairIn = aSourceField->get_mesh_pair();

            moris::mtk::Mesh * tSourceMesh = tMeshPairIn->mInterpolationMesh;

            MORIS_ERROR( tSourceMesh->get_mesh_type() == MeshType::HMR,
                    "Mapper::map_input_field_to_output_field() Source mesh is not and HMR mesh" );

            std::shared_ptr< hmr::Database > tHMRDatabase = tSourceMesh->get_HMR_database();

            // get Lagrange and Discretization orer
            uint tSourceLagrangeOrder = aSourceField->get_lagrange_order();
            uint tDiscretizationOrder = aSourceField->get_discretization_order();

            // extract pattern from mesh on which this field id based on
            // both Lagrange and discretization order if they are not the same
            Matrix< DDLUMat > aElementCounterPerLevelAndPattern;
            moris::Cell< Matrix< DDLUMat > > aElementPerPattern;

            hmr::File tFile;
            tFile.save_refinement_pattern(
                    tSourceMesh->get_HMR_lagrange_mesh(),
                    aSourceField->get_discretization_mesh_index(),
                    aElementCounterPerLevelAndPattern,
                    aElementPerPattern );

            // copy HMR parameters from old HMR performer
            hmr::Parameters * tParameters = aHMRPerformers( 0 )->get_parameters();

            //FIXME adjust parameters here. eg no initial refinement

            // build new HMR performer with copied parameters
            std::shared_ptr< hmr::HMR > tHMRPerformerNew =
                    std::make_shared< hmr::HMR >( tParameters );

            // FIXME put this line into HMR
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

            // Create  mesh pair
            mtk::Mesh_Pair tMeshPairOld;
            tMeshPairOld.mInterpolationMesh = tOldInterpolationMesh;
            tMeshPairOld.mIsOwned   = true;

            // build field with mesh
            mtk::Field tFieldOld( &tMeshPairOld );

            // copy values from input mesh to New/Old mesh ( this mesh is build based on the new HMR performer )
            tFieldOld.unlock_field();
            tFieldOld.set_nodal_values( aSourceField->get_nodal_values() );

            this->perform_refinement( tHMRPerformerNew, &tFieldOld );

            tHMRPerformerNew->finalize();

            aHMRPerformers( 0 ) = tHMRPerformerNew;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Remeshing_Mini_Performer::perform_refinement(
                std::shared_ptr< hmr::HMR >   aHMRPerformer,
                mtk::Field                  * aSourceField )
        {
            // switch based on mode index
            switch ( mParameters.mModeIndex )
            {
                case 0 :
                {
                    this->perform_refinement_mode_0( aHMRPerformer, aSourceField );

                    break;
                }
                case 1 :
                {
                    this->perform_refinement_mode_1( aHMRPerformer, aSourceField );

                    break;
                }
                default:
                    MORIS_ERROR( false, "Remeshing_Mini_Performer::perform_refinement() - Refinement mode does not exist" );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Remeshing_Mini_Performer::perform_refinement_mode_0(
                std::shared_ptr< hmr::HMR >   aHMRPerformer,
                mtk::Field                  * aSourceField )
        {
            // perform initial uniform refinement
            aHMRPerformer->perform_initial_refinement();

            // get input field order
            uint tLagrangeOrder = aSourceField->get_lagrange_order();
            uint tDiscretizationOrder = aSourceField->get_discretization_order();

            std::shared_ptr< hmr::Database > tHMRDatabase = aHMRPerformer->get_database();

            // get refinement level and pattern
            Matrix< DDUMat > tRefinementLevel = mParameters.mRefinementsMode_0;
            Matrix< DDUMat > tRefinementPattern = mParameters.mRefinementPatternMode_0;

            uint tNumPattern = tRefinementPattern.numel();

            // loop over pattern
            for( uint Ik = 0; Ik < tNumPattern; Ik++ )
            {
                // get pattern
                uint tPattern = tRefinementPattern( Ik );

                for( uint Ii = 0; Ii < tRefinementLevel( Ik ); Ii++ )
                {

                    // create mesh for this pattern
                    hmr::Interpolation_Mesh_HMR * tInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
                            tHMRDatabase,
                            tLagrangeOrder,
                            tPattern,
                            tDiscretizationOrder,
                            tPattern);

                    mtk::Mesh_Pair tMeshPair;
                    tMeshPair.mInterpolationMesh = tInterpolationMesh;
                    tMeshPair.mIsOwned   = true;

                    // create field object for this mesh
                    mtk::Field_Discrete tFieldOnPattern( &tMeshPair, tPattern );
                    tFieldOnPattern.set_label( "Field_for_refinement" );

                    // create mapper and map input field to new field
                    mtk::Mapper tMapper;
                    tMapper.map_input_field_to_output_field( aSourceField, &tFieldOnPattern );

                    // create refinement parameter list
                    moris::ParameterList tRefinementParameterlist;
                    prm::create_refinement_parameterlist( tRefinementParameterlist );
                    tRefinementParameterlist.set( "field_names" , "Field_for_refinement" );
                    tRefinementParameterlist.set( "levels_of_refinement" , ios::stringify( 1 ) );
                    tRefinementParameterlist.set( "refinement_pattern" , ios::stringify( tPattern ) );

                    // put field pointer in cell
                    Cell< mtk::Field* > tFields( 1, &tFieldOnPattern );

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
                std::shared_ptr< hmr::HMR >   aHMRPerformer,
                mtk::Field                  * aSourceField )
        {

        }

    }
}

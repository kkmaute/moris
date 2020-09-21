
#include "cl_VIS_Output_Manager.hpp"

#include "cl_VIS_Factory.hpp"

#include "cl_MDL_Model.hpp"
#include "cl_MSI_Equation_Set.hpp"
#include "cl_MSI_Equation_Model.hpp"

#include "cl_FEM_Set.hpp"

#include "fn_Parsing_Tools.hpp"

// Logging package
#include "cl_Logger.hpp"

extern moris::Comm_Manager gMorisComm;
extern moris::Logger       gLogger;

namespace moris
{
    namespace vis
    {

        //-----------------------------------------------------------------------------------------------------------

        void Output_Manager::set_outputs(
                const uint                              aOutputIndex,
                const enum VIS_Mesh_Type                aMeshType,
                const std::string                     & aMeshPath,
                const std::string                     & aMeshName,
                const std::string                     & aTempPath,
                const std::string                     & aTempName,
                const moris::Cell< std::string >      & aBlockNames,
                const moris::Cell< std::string >      & aFieldNames,
                const moris::Cell< enum Field_Type >  & aFieldType,
                const moris::Cell< std::string >      & aQINames,
                const uint                              aSaveFrequency,
                const real                              aTimeOffset)
        {
            // create output data object
            vis::Output_Data tOutputData;

            // fill output data object
            tOutputData.mMeshIndex  = aOutputIndex;
            tOutputData.mMeshType   = aMeshType;
            tOutputData.mMeshName   = aMeshName;
            tOutputData.mMeshPath   = aMeshPath;
            tOutputData.mTempName   = aTempName;
            tOutputData.mTempPath   = aTempPath;
            tOutputData.mSetNames   = aBlockNames;
            tOutputData.mFieldNames = aFieldNames;
            tOutputData.mFieldType  = aFieldType;
            tOutputData.mQINames    = aQINames;

            tOutputData.mSaveFrequency = aSaveFrequency;
            tOutputData.mTimeOffset    = aTimeOffset;

            // resize list of output data objects
            uint tSize = mOutputData.size();
            uint OutputDataSize = std::max( tSize, aOutputIndex + 1 );

            mOutputData.resize( OutputDataSize );

            // assign output data object to list
            mOutputData( aOutputIndex ) = tOutputData;

            // resize mesh list
            mVisMesh.resize( mOutputData.size(), nullptr );

            mVisMeshCreatedAndOpen.resize( mOutputData.size(), false );
        }

        //-----------------------------------------------------------------------------------------------------------

        void Output_Manager::set_outputs( moris::ParameterList aParamterelist )
        {
            // create output data object
            vis::Output_Data tOutputData;

            // fill output data object
            tOutputData.mMeshIndex  = aParamterelist.get< moris::sint >( "Output_Index" );
            tOutputData.mMeshType   = static_cast< moris::vis::VIS_Mesh_Type >( aParamterelist.get< moris::uint >( "Mesh_Type" ) );

            tOutputData.mOutputPath = std::get< 0 >( aParamterelist.get< std::pair< std::string, std::string > >( "File_Name" ) );
            tOutputData.mMeshName   = std::get< 1 >( aParamterelist.get< std::pair< std::string, std::string > >( "File_Name" ) );

            // note: file path for temp file currently ignored
            tOutputData.mTempPath   = std::get< 0 >( aParamterelist.get< std::pair< std::string, std::string > >( "Temp_Name" ) );
            tOutputData.mTempName   = std::get< 1 >( aParamterelist.get< std::pair< std::string, std::string > >( "Temp_Name" ) );

            tOutputData.mSaveFrequency  = aParamterelist.get< moris::sint >( "Save_Frequency" );
            tOutputData.mTimeOffset     = aParamterelist.get< moris::real >( "Time_Offset" );

            moris::Cell< std::string > tSetNames;
            string_to_cell( aParamterelist.get< std::string >( "Set_Names" ), tSetNames );
            tOutputData.mSetNames   = tSetNames;

            moris::Cell< std::string > tFieldNames;
            string_to_cell(
                    aParamterelist.get< std::string >( "Field_Names" ),
                    tFieldNames );
            tOutputData.mFieldNames = tFieldNames;

            moris::Cell< enum vis::Field_Type > tFieldTypes;
            moris::map< std::string, enum vis::Field_Type > tFieldTypeMap = get_vis_field_type_map();
            string_to_cell(
                    aParamterelist.get< std::string >( "Field_Type" ) ,
                    tFieldTypes,
                    tFieldTypeMap );
            tOutputData.mFieldType  = tFieldTypes;

            // check that length of Field_Names and Field_Type are consistent
            MORIS_ERROR( tFieldNames.size() == tFieldTypes.size(),"Output_Manager::set_outputs - Number of Field Names and Field Types differ.");

            moris::Cell< std::string > tQINames;
            string_to_cell( aParamterelist.get< std::string >( "IQI_Names"), tQINames );
            tOutputData.mQINames = tQINames;

            // check that length of Field_Names and Field_Type are consistent
            MORIS_ERROR( tFieldNames.size() == tQINames.size(),"Output_Manager::set_outputs - Number of Field Names and QI Names differ.");

            // resize list of output data objects
            sint tSize = mOutputData.size();
            sint OutputDataSize = std::max( tSize, tOutputData.mMeshIndex + 1 );

            mOutputData.resize( OutputDataSize );

            // assign output data object to list
            mOutputData( tOutputData.mMeshIndex ) = tOutputData;

            // resize mesh list
            mVisMesh.resize( mOutputData.size(), nullptr );

            mVisMeshCreatedAndOpen.resize( mOutputData.size(), false );
        }

        //-----------------------------------------------------------------------------------------------------------

        void Output_Manager::create_visualization_mesh(
                const uint          aVisMeshIndex,
                mtk::Mesh_Manager * aMesh,
                const uint          aMeshPairIndex)
        {
            MORIS_ERROR( mOutputData( aVisMeshIndex ).mMeshIndex == ( sint )aVisMeshIndex,
                    "create_visualization_meshes(), Visualization mesh not set" );

            mMTKMesh = aMesh;
            mMTKMeshPairIndex = aMeshPairIndex;

            // create vis factory
            vis::VIS_Factory tVisFactory( aMesh, aMeshPairIndex );

            // create vis mesh
            mVisMesh( aVisMeshIndex ) = tVisFactory.create_visualization_mesh( mOutputData( aVisMeshIndex ) );

            // resize list of writers to list of outputs. memory allocation stays intact
            mWriter.resize( mOutputData.size(), nullptr );

            // create writer for this mesh
            mWriter( aVisMeshIndex ) = new moris::mtk::Writer_Exodus( mVisMesh( aVisMeshIndex ) );

        }

        //-----------------------------------------------------------------------------------------------------------

        void Output_Manager::set_visualization_sets(
                const uint                             aVisMeshIndex,
                std::shared_ptr< MSI::Equation_Model > aEquationModel )
        {
            // get number of requested sets
            uint tNumRequestedSets = mOutputData( aVisMeshIndex ).mSetNames.size();

            // get mtk set index to fem set index map
            //map< moris_index, moris_index > & tMeshSetToFemSetMap
            //map< std::pair< moris_index, bool >, moris_index > & tMeshSetToFemSetMap
            map< std::tuple< moris_index, bool, bool >, moris_index > & tMeshSetToFemSetMap =
                    aEquationModel->get_mesh_set_to_fem_set_index_map();

            // get equation sets
            moris::Cell< MSI::Equation_Set * > tEquationSets
            = aEquationModel->get_equation_sets();

            // get integration mesh
            mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;

            mMTKMesh->get_mesh_pair( mMTKMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

            // loop over equation sets.
            for( uint Ii = 0; Ii < tNumRequestedSets; Ii++ )
            {
                // get mtk set index
                moris_index tSetIndex = tIntegrationMesh->get_set_index_by_name( mOutputData( aVisMeshIndex ).mSetNames( Ii ) );

                if ( tMeshSetToFemSetMap.key_exists( std::make_tuple( tSetIndex, false, false ) ) )
                {
                    // find set index for this block index
                    moris_index tEquationSetIndex = tMeshSetToFemSetMap.find( std::make_tuple( tSetIndex, false, false ) );

                    // set vis set to fem set. +1 because 0 is reserved for fem
                    tEquationSets(
                            tEquationSetIndex )->set_visualization_set(
                                    aVisMeshIndex + 1,
                                    mVisMesh( aVisMeshIndex )->get_set_by_index( Ii ),
                                    mOnlyPrimary );
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void Output_Manager::write_mesh( const uint aVisMeshIndex )
        {
            // specify file path
            std::string tMeshFilePath = mOutputData( aVisMeshIndex ).mMeshPath;

            // get file name
            std::string tMeshFileName = mOutputData( aVisMeshIndex ).mMeshName;

            // specify file path for temporary file
            std::string tMeshTempPath = mOutputData( aVisMeshIndex ).mTempPath;

            // get file name of temporary file
            std::string tMeshTempName = mOutputData( aVisMeshIndex ).mTempName;

            // augment file name if time offset > 0
            if ( mOutputData( aVisMeshIndex ).mTimeOffset > 0)
            {
                // get optimization iteration
                uint tOptIter = gLogger.get_opt_iteration();

                // set name
                std::string tOptIterStrg = std::to_string(tOptIter);
                tMeshFileName += ".e-s." + std::string(4-tOptIterStrg.length(),'0') + tOptIterStrg;

                // determine time shift
                mTimeShift = tOptIter * mOutputData( aVisMeshIndex ).mTimeOffset;
            }

            std::string tMassage = "Writing " + tMeshFileName + " to " + tMeshFilePath +".";

            MORIS_LOG( tMassage.c_str() );

            // write mesh to file
            mWriter( aVisMeshIndex )->write_mesh( tMeshFilePath, tMeshFileName, tMeshTempPath, tMeshTempName );

            // add nodal elemental and global fields to mesh
            this->add_nodal_fields( aVisMeshIndex );
            this->add_elemetal_fields( aVisMeshIndex );
            this->add_global_fields( aVisMeshIndex );

            // write standard outputs like IDs and Indices to file
            this->write_mesh_indices( aVisMeshIndex );

            // reset field write counter
            mOutputData( aVisMeshIndex ).mFieldWriteCounter=0;
        }

        //-----------------------------------------------------------------------------------------------------------

        void Output_Manager::add_nodal_fields( const uint aVisMeshIndex )
        {
            // set list of nodal field names to input + 2
            moris::Cell<std::string> tNodalFieldNames( 2 + mOutputData( aVisMeshIndex ).mFieldNames.size() );

            // set standard field names
            tNodalFieldNames( 0 ) = "Mesh_Id";
            tNodalFieldNames( 1 ) = "Mesh_Index";

            uint tCounter = 2;

            // loop over field names and check if fields are nodal fields
            for( uint Ik = 0; Ik < mOutputData( aVisMeshIndex ).mFieldNames.size(); Ik++ )
            {
                if( mOutputData( aVisMeshIndex ).mFieldType( Ik ) == Field_Type::NODAL )
                {
                    tNodalFieldNames( tCounter++ ) = mOutputData( aVisMeshIndex ).mFieldNames( Ik );
                }
            }

            tNodalFieldNames.resize( tCounter );

            // pass nodal field names to writer
            mWriter( aVisMeshIndex )->set_nodal_fields( tNodalFieldNames );
        }

        //-----------------------------------------------------------------------------------------------------------

        void Output_Manager::add_elemetal_fields( const uint aVisMeshIndex )
        {
            // allocate cell for storing elemental field names; 3 default fields are added
            moris::Cell<std::string> tElementalFieldNames( 3 + mOutputData( aVisMeshIndex ).mFieldNames.size() );

            tElementalFieldNames( 0 ) = "Mesh_Id";
            tElementalFieldNames( 1 ) = "Mesh_Index";
            tElementalFieldNames( 2 ) = "Proc_Index";

            // set field counter to 3 to account for default fields
            uint tCounter = 3;

            // loop over field names and check if fields are elemental fields
            for( uint Ik = 0; Ik < mOutputData( aVisMeshIndex ).mFieldNames.size(); Ik++ )
            {
                if( mOutputData( aVisMeshIndex ).mFieldType( Ik ) == Field_Type::ELEMENTAL )
                {
                    tElementalFieldNames( tCounter++ ) = mOutputData( aVisMeshIndex ).mFieldNames( Ik );
                }
            }

            // trim cell of element field names
            tElementalFieldNames.resize( tCounter );

            // write field names to file
            mWriter( aVisMeshIndex )->set_elemental_fields( tElementalFieldNames );
        }

        //-----------------------------------------------------------------------------------------------------------

        void Output_Manager::add_global_fields( const uint aVisMeshIndex )
        {
            moris::Cell<std::string> tGlobalFieldNames( mOutputData( aVisMeshIndex ).mFieldNames.size() );

            uint tCounter = 0;

            // loop over field names and check if fields are global fields
            for( uint Ik = 0; Ik < mOutputData( aVisMeshIndex ).mFieldNames.size(); Ik++ )
            {
                if( mOutputData( aVisMeshIndex ).mFieldType( Ik ) == Field_Type::GLOBAL )
                {
                    tGlobalFieldNames( tCounter++ ) = mOutputData( aVisMeshIndex ).mFieldNames( Ik );
                }
            }

            tGlobalFieldNames.resize( tCounter );

            mWriter( aVisMeshIndex )->set_global_variables( tGlobalFieldNames );
        }

        //-----------------------------------------------------------------------------------------------------------

        void Output_Manager::write_mesh_indices( const uint aVisMeshIndex )     //FIXME
        {
            // get mesh set indices
            uint tRequestedSets = mOutputData( aVisMeshIndex ).mSetNames.size();

            // loop over mesh sets
            for( uint Ii = 0; Ii < tRequestedSets; Ii++ )
            {
                // get vis set by index
                moris::mtk::Set * tSet = mVisMesh( aVisMeshIndex )->get_set_by_index( Ii );

                bool tOnlyPrimaryCells = true ;

                switch( mOutputData( aVisMeshIndex ).mMeshType )
                {
                    case  vis::VIS_Mesh_Type::STANDARD:
                        tOnlyPrimaryCells = true ;
                        break;

                    case  vis::VIS_Mesh_Type::OVERLAPPING_INTERFACE:
                        tOnlyPrimaryCells = false;
                        break;

                    case  vis::VIS_Mesh_Type::FULL_DISCONTINOUS:
                        MORIS_ERROR( false, "create_visualization_mesh() - Mesh type FULL_DISCONTINOUS not implemented yet. " );
                        break;

                    default:
                        MORIS_ERROR( false, "create_visualization_mesh() - Mesh type not specified. " );
                        break;
                }

                // get number of cells on set
                uint tNumCells = tSet->get_num_cells_on_set( tOnlyPrimaryCells );

                // check whether number of cells > 0; otherwise skip remainder
                if (tNumCells == 0)
                {
                    continue;
                }

                // get cell indices on set
                moris::Matrix< DDSMat > tCellIndex = tSet->get_cell_inds_on_block( tOnlyPrimaryCells );

                // find the maximal index for resizing purposes
                sint tMaxIndex = tCellIndex.max();

                // create cell assembly map( index to position )
                Matrix< DDSMat > tCellAsseblyMap( tMaxIndex + 1, 1, -1 );

                // loop over cells and put them in map
                for( uint Ik = 0; Ik < tNumCells; Ik++ )
                {
                    tCellAsseblyMap( tCellIndex( Ik ) ) = Ik;
                }

                // create cell of index, id field, and proc indices
                moris::Cell< Matrix< DDRMat > > tIdIndex( 3 );
                tIdIndex( 0 ).set_size( tCellIndex.numel(), 1 );
                tIdIndex( 1 ).set_size( tCellIndex.numel(), 1 );
                tIdIndex( 2 ).set_size( tCellIndex.numel(), 1 );

                // get clusters from vis set
                moris::Cell< mtk::Cluster const* > tMeshClusterList = tSet->get_clusters_on_set();

                // loop over clusters and get ids and indices
                for( uint Ik = 0; Ik < tMeshClusterList.size(); Ik++ )
                {
                    // get primary cells
                    const moris::Cell<moris::mtk::Cell const *> & tPrimaryCells =
                            tMeshClusterList( Ik )->get_primary_cells_in_cluster();

                    // loop over primary cells
                    for( uint Ia = 0; Ia < tPrimaryCells.size(); Ia++ )
                    {
                        // get index of vis cell
                        moris_index tIndex = tPrimaryCells( Ia )->get_index();

                        moris_id tMeshId =
                                reinterpret_cast< const vis::Cell_Visualization* >( tPrimaryCells( Ia ) )->get_mesh_cell_id();

                        moris_index tMeshIndex =
                                reinterpret_cast< const vis::Cell_Visualization* >( tPrimaryCells( Ia ) )->get_mesh_cell_index();

                        tIdIndex( 0 )( tCellAsseblyMap( tIndex ) ) = tMeshId;
                        tIdIndex( 1 )( tCellAsseblyMap( tIndex ) ) = tMeshIndex;
                        tIdIndex( 2 )( tCellAsseblyMap( tIndex ) ) = par_rank();
                    }

                    const moris::Cell<moris::mtk::Cell const *> & tVoidCells = tMeshClusterList( Ik )->get_void_cells_in_cluster();

                    for( uint Ia = 0; Ia < tVoidCells.size(); Ia++ )
                    {
                        moris_index tIndex = tVoidCells( Ia )->get_index();

                        moris_id tMeshId =
                                reinterpret_cast< const vis::Cell_Visualization* >( tVoidCells( Ia ) )->get_mesh_cell_id();

                        moris_index tMeshIndex =
                                reinterpret_cast< const vis::Cell_Visualization* >( tVoidCells( Ia ) )->get_mesh_cell_index();

                        tIdIndex( 0 )( tCellAsseblyMap( tIndex ) ) = tMeshId;
                        tIdIndex( 1 )( tCellAsseblyMap( tIndex ) ) = tMeshIndex;
                        tIdIndex( 2 )( tCellAsseblyMap( tIndex ) ) = par_rank();
                    }
                }

                mWriter( aVisMeshIndex )->write_elemental_field( tSet->get_set_name() , "Mesh_Id",    tIdIndex( 0 ) );
                mWriter( aVisMeshIndex )->write_elemental_field( tSet->get_set_name() , "Mesh_Index", tIdIndex( 1 ) );
                mWriter( aVisMeshIndex )->write_elemental_field( tSet->get_set_name() , "Proc_Index", tIdIndex( 2 ) );
            }
        }

        //-----------------------------------------------------------------------------------------------------------
        void Output_Manager::write_field(
                const uint                             aVisMeshIndex,
                const real                             aTime,
                std::shared_ptr< MSI::Equation_Model > aEquationModel )
        {
            // increment field write counter
            mOutputData( aVisMeshIndex ).mFieldWriteCounter++;

            // write time to file
            mWriter( aVisMeshIndex )->set_time( aTime + mTimeShift );

            // get mesh set to fem set index map
            map< std::tuple< moris_index, bool, bool >, moris_index > & tMeshSetToFemSetMap =
                    aEquationModel->get_mesh_set_to_fem_set_index_map( );

            // get equation sets
            moris::Cell< MSI::Equation_Set * > tEquationSets = aEquationModel->get_equation_sets();

            // get integration mesh
            mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;

            mMTKMesh->get_mesh_pair( mMTKMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

            // number of fields in vis mesh
            uint tNumFields = mOutputData( aVisMeshIndex ).mFieldNames.size();

            // allocate vectors for names and values of global variables
            moris::Cell<std::string> tGlobalVariableNames(tNumFields);
            moris::Matrix<DDRMat>    tGlobalVarialbeValues(tNumFields,1,MORIS_REAL_MAX);
            uint tGlobalVariablesCounter = 0;

            // loop over all fields of this output object
            for( uint Ik = 0; Ik < tNumFields; Ik++ )
            {
                // get field name
                std::string tFieldName = mOutputData( aVisMeshIndex ).mFieldNames( Ik );

                // nodal and global field values
                Matrix< DDRMat > tNodalValues( mVisMesh( aVisMeshIndex )->get_num_nodes(), 1, std::numeric_limits<real>::quiet_NaN() );
                moris::real      tGlobalValue = 0.0;

                // loop over all blocks on this output object
                for( uint Ii = 0; Ii < mOutputData( aVisMeshIndex ).mSetNames.size(); Ii++ )
                {
                    moris_index tSetIndex = tIntegrationMesh->get_set_index_by_name( mOutputData( aVisMeshIndex ).mSetNames( Ii ) );

                    if ( tMeshSetToFemSetMap.key_exists( std::make_tuple( tSetIndex, false, false ) ) )
                    {
                        // find set index for this block index
                        moris_index tEquationSetIndex = tMeshSetToFemSetMap.find( std::make_tuple( tSetIndex, false, false ) );

                        // elemental field values
                        Matrix< DDRMat > tElementValues;

                        tEquationSets( tEquationSetIndex )->compute_quantity_of_interest(
                                aVisMeshIndex + 1,
                                &tElementValues,
                                &tNodalValues,
                                &tGlobalValue,
                                mOutputData( aVisMeshIndex ).mQINames( Ik ),
                                mOutputData( aVisMeshIndex ).mFieldType( Ik ) );

                        if( mOutputData( aVisMeshIndex ).mFieldType( Ik ) == Field_Type::ELEMENTAL )
                        {
                            mWriter( aVisMeshIndex )->write_elemental_field(
                                    mOutputData( aVisMeshIndex ).mSetNames( Ii ),
                                    tFieldName,
                                    tElementValues );
                        }
                    }
                }

                real tTotalGlobalValue;
                switch ( mOutputData( aVisMeshIndex ).mFieldType( Ik ) )
                {
                    case Field_Type::NODAL:
                        mWriter( aVisMeshIndex )->write_nodal_field( tFieldName, tNodalValues );
                        break;
                    case Field_Type::GLOBAL:
                        tTotalGlobalValue = sum_all(tGlobalValue);
                        tGlobalVariableNames( tGlobalVariablesCounter )  = tFieldName;
                        tGlobalVarialbeValues(tGlobalVariablesCounter,0) = tTotalGlobalValue;
                        tGlobalVariablesCounter++;
                        MORIS_LOG_INFO ("Global Variable: %s = %e",tFieldName.c_str(),tTotalGlobalValue);
                        break;
                    case Field_Type::ELEMENTAL:
                        // do nothing here - case is handled above
                        break;
                    default:
                        MORIS_ERROR(false,"undefined FieldType option\n");
                }

            }

            // write global variables
            if ( tGlobalVariablesCounter > 0 )
            {
                // trim list of global variable names as correct size needed by writing routine
                tGlobalVariableNames.resize(tGlobalVariablesCounter);

                // write global variables
                mWriter( aVisMeshIndex )->write_global_variables( tGlobalVariableNames, tGlobalVarialbeValues );
            }

            // check if a copy of the current mesh file should be created
            sint tFieldWriteCounter = mOutputData( aVisMeshIndex ).mFieldWriteCounter;
            sint tSaveFrequency     = mOutputData( aVisMeshIndex ).mSaveFrequency;

            if ( std::remainder( tFieldWriteCounter,tSaveFrequency) == 0 )
            {
                mWriter( aVisMeshIndex )->save_mesh( );
            }
        }

        //-----------------------------------------------------------------------------------------------------------
    }
}

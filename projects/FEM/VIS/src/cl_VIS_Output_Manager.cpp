
#include "cl_VIS_Output_Manager.hpp"

#include "cl_VIS_Factory.hpp"

#include "cl_MDL_Model.hpp"
#include "cl_MSI_Equation_Set.hpp"
#include "cl_MSI_Equation_Model.hpp"

#include "cl_FEM_Set.hpp"

#include "fn_Parsing_Tools.hpp"

// Logging package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

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

            // read and check mesh set names
            moris::Cell< std::string > tSetNames;
            string_to_cell(
                    aParamterelist.get< std::string >( "Set_Names" ),
                    tSetNames );

            MORIS_ERROR( tSetNames.size() > 0, "At least one mesh set name needs to be provided for Vis mesh\n");
            for (auto tName:tSetNames)
            {
                MORIS_ERROR( tName.length() > 0, "Empty strings for set names in Vis mesh are not allowed\n");
            }

            tOutputData.mSetNames   = tSetNames;

            // read and check field names
            moris::Cell< std::string > tFieldNames;
            string_to_cell(
                    aParamterelist.get< std::string >( "Field_Names" ),
                    tFieldNames );

            MORIS_ERROR( tFieldNames.size() > 0, "At least one field name needs to be provided for Vis mesh\n");
            for (auto tName:tFieldNames)
            {
                MORIS_ERROR( tName.length() > 0, "Empty strings for field names in Vis mesh are not allowed\n");
            }

            tOutputData.mFieldNames = tFieldNames;

            // read and check field types
            moris::Cell< enum vis::Field_Type > tFieldTypes;
            moris::map< std::string, enum vis::Field_Type > tFieldTypeMap = get_vis_field_type_map();
            string_to_cell(
                    aParamterelist.get< std::string >( "Field_Type" ) ,
                    tFieldTypes,
                    tFieldTypeMap );

            MORIS_ERROR( tFieldTypes.size() > 0, "At least one field type needs to be provided for Vis mesh\n");
            for (auto tName:tFieldNames)
            {
                MORIS_ERROR( tName.length() > 0, "Empty strings for field types in Vis mesh are not allowed\n");
            }

            tOutputData.mFieldType  = tFieldTypes;

            // check that length of Field_Names and Field_Type are consistent
            MORIS_ERROR( tFieldNames.size() == tFieldTypes.size(),"Output_Manager::set_outputs - Number of Field Names and Field Types differ.");

            // read and check IQI names
            moris::Cell< std::string > tQINames;
            string_to_cell(
                    aParamterelist.get< std::string >( "IQI_Names"),
                    tQINames );

            MORIS_ERROR( tQINames.size() > 0, "At least one IQI name needs to be provided for Vis mesh\n");
            for (auto tName:tQINames)
            {
                MORIS_ERROR( tName.length() > 0, "Empty strings for IQI name in Vis mesh are not allowed\n");
            }

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

        void Output_Manager::setup_vis_mesh_for_output(
                const uint                               aVisMeshIndex,
                std::shared_ptr< mtk::Mesh_Manager >     aMesh,
                const uint                               aMeshPairIndex,
                std::shared_ptr< MSI::Equation_Model >   aEquationModel)
        {
            if( mVisMeshCreatedAndOpen( aVisMeshIndex ) == false )
            {
                Tracer tTracer( "Output_Manager", "VisMesh", "CreateAndWriteVisMesh" );

                this->create_visualization_mesh( aVisMeshIndex, aMesh, aMeshPairIndex );

                this->set_visualization_sets( aVisMeshIndex, aEquationModel );

                this->write_mesh( aVisMeshIndex );

                mVisMeshCreatedAndOpen( aVisMeshIndex ) = true;
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void Output_Manager::create_visualization_mesh(
                const uint                           aVisMeshIndex,
                std::shared_ptr< mtk::Mesh_Manager > aMesh,
                const uint                           aMeshPairIndex)
        {
            Tracer tTracer( "Output_Manager", "VisMesh", "CreateVisMesh" );

            MORIS_ERROR( mOutputData( aVisMeshIndex ).mMeshIndex == ( sint )aVisMeshIndex,
                    "create_visualization_meshes(), Visualization mesh not set" );

            mMTKMesh          = aMesh;
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
            Tracer tTracer( "Output_Manager", "VisMesh", "CreateVisSets" );

            // get number of requested sets
            uint tNumRequestedSets = mOutputData( aVisMeshIndex ).mSetNames.size();

            // get mtk set index to fem set index map
            map< std::tuple< moris_index, bool, bool >, moris_index > & tMeshSetToFemSetMap =
                    aEquationModel->get_mesh_set_to_fem_set_index_map();

            // get equation sets
            moris::Cell< MSI::Equation_Set * > tEquationSets = aEquationModel->get_equation_sets();

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
            Tracer tTracer( "Output_Manager", "VisMesh", "WriteVisMesh" );

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

            // Log mesh writing message
            std::string tMessage = "Writing " + tMeshFileName + " to " + tMeshFilePath +".";
            MORIS_LOG( tMessage.c_str() );

            // write mesh to file
            mWriter( aVisMeshIndex )->write_mesh(
                    tMeshFilePath,
                    tMeshFileName,
                    tMeshTempPath,
                    tMeshTempName);

            // add nodal elemental and global fields to mesh
            this->add_nodal_fields( aVisMeshIndex );
            this->add_elemental_fields( aVisMeshIndex );
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
            tNodalFieldNames( 0 ) = "Node_Id";
            tNodalFieldNames( 1 ) = "Node_Index";

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

        void Output_Manager::add_elemental_fields( const uint aVisMeshIndex )
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

            // getting the vertices on the vis mesh
            moris::Cell< moris::mtk::Vertex const * > tVertices = mVisMesh( aVisMeshIndex )->get_all_vertices();

            // get number of verts
            uint tNumVerts = tVertices.size();

            // create cell of index, id field, and proc indices
            moris::Cell< Matrix< DDRMat > > tVertIdIndex( 2 );
            tVertIdIndex( 0 ).set_size( tNumVerts, 1 );
            tVertIdIndex( 1 ).set_size( tNumVerts, 1 );

            for ( uint iVert = 0; iVert < tNumVerts; iVert++ )
            {
                // get the vis index
                moris_index tIndex = tVertices( iVert )->get_index();

                // get the moris ID
                moris_id tVertId =
                        reinterpret_cast< const vis::Vertex_Visualization* >( tVertices( iVert ) )->get_integration_id();

                // get the moris index
                moris_index tVertIndex =
                        reinterpret_cast< const vis::Vertex_Visualization* >( tVertices( iVert ) )->get_integration_index();

                // assign moris ID and index
                tVertIdIndex( 0 )( tIndex ) = tVertId;
                tVertIdIndex( 1 )( tIndex ) = tVertIndex;
            }

            // write the moris vertex ids and indices
            mWriter( aVisMeshIndex )->write_nodal_field( "Node_Id",    tVertIdIndex( 0 ) );
            mWriter( aVisMeshIndex )->write_nodal_field( "Node_Index", tVertIdIndex( 1 ) );
        }

        //-----------------------------------------------------------------------------------------------------------

        void Output_Manager::write_field(
                const uint                             aVisMeshIndex,
                const real                             aTime,
                std::shared_ptr< MSI::Equation_Model > aEquationModel )
        {
            Tracer tTracer( "Output_Manager", "VisMesh", "WriteFields" );

            // number of fields in vis mesh
            uint tNumFields = mOutputData( aVisMeshIndex ).mFieldNames.size();

            // number of set names
            uint tNumSetNames = mOutputData( aVisMeshIndex ).mSetNames.size();

            // loop over field names and check if fields are global fields
            moris::Cell< std::string > tGlobalIQINames;
            moris::Cell< std::string > tGlobalFieldNames;
            moris::Cell< std::string > tNodalIQINames;
            moris::Cell< std::string > tNodalFieldNames;
            moris::Cell< std::string > tElementalIQINames;
            moris::Cell< std::string > tElementalFieldNames;

            for( uint Ik = 0; Ik < tNumFields; Ik++ )
            {
                // get the field type
                Field_Type tFieldType = mOutputData( aVisMeshIndex ).mFieldType( Ik );

                // get the iQI name
                std::string tIQIName = mOutputData( aVisMeshIndex ).mQINames( Ik );

                // get the field name
                std::string tFieldName = mOutputData( aVisMeshIndex ).mFieldNames( Ik );

                // switch on field type
                switch( tFieldType )
                {
                    case Field_Type::GLOBAL :
                        tGlobalIQINames.push_back( tIQIName );
                        tGlobalFieldNames.push_back( tFieldName );
                        break;
                    case Field_Type::NODAL :
                        tNodalIQINames.push_back( tIQIName );
                        tNodalFieldNames.push_back( tFieldName );
                        break;
                    case Field_Type::ELEMENTAL :
                        tElementalIQINames.push_back( tIQIName );
                        tElementalFieldNames.push_back( tFieldName );
                        break;
                    default:
                        MORIS_ERROR( false, "Unknown field type.");
                }
            }
            uint tNumGlobalIQIs    = tGlobalIQINames.size();
            uint tNumNodalIQIs     = tNodalIQINames.size();
            uint tNumElementalIQIs = tElementalIQINames.size();

            // increment field write counter
            mOutputData( aVisMeshIndex ).mFieldWriteCounter++;

            // write time to file
            mWriter( aVisMeshIndex )->set_time( aTime + mTimeShift );

            // get mesh set to fem set index map
            map< std::tuple< moris_index, bool, bool >, moris_index > & tMeshSetToFemSetMap =
                    aEquationModel->get_mesh_set_to_fem_set_index_map( );

            // get equation sets
            moris::Cell< MSI::Equation_Set * > & tEquationSets =
                    aEquationModel->get_equation_sets();

            // get integration mesh
            mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;

            mMTKMesh->get_mesh_pair( mMTKMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

            // init nodal field values
            Matrix< DDRMat > tNodalValues(
                    mVisMesh( aVisMeshIndex )->get_num_nodes(),
                    tNumNodalIQIs,
                    std::numeric_limits<real>::quiet_NaN() );

            // init global field values
            Matrix< DDRMat > tGlobalValues( 1, tNumGlobalIQIs, 0.0 );

            // loop over all blocks on this output object
            for( uint Ii = 0; Ii < tNumSetNames; Ii++ )
            {
                // get mesh set index from name
                moris_index tSetIndex = tIntegrationMesh->get_set_index_by_name(
                        mOutputData( aVisMeshIndex ).mSetNames( Ii ) );

                if ( tMeshSetToFemSetMap.key_exists( std::make_tuple( tSetIndex, false, false ) ) )
                {
                    // find set index for this block index
                    moris_index tEquationSetIndex =
                            tMeshSetToFemSetMap.find( std::make_tuple( tSetIndex, false, false ) );

                    // global values
                    tEquationSets( tEquationSetIndex )->compute_quantity_of_interest_global(
                            aVisMeshIndex + 1,
                            &tGlobalValues,
                            tGlobalIQINames );

                    // nodal values
                    tEquationSets( tEquationSetIndex )->compute_quantity_of_interest_nodal(
                            aVisMeshIndex + 1,
                            &tNodalValues,
                            tNodalIQINames );

                    // elemental field values
                    Matrix< DDRMat > tElementValues;

                    // elemental values
                    tEquationSets( tEquationSetIndex )->compute_quantity_of_interest_elemental(
                            aVisMeshIndex + 1,
                            &tElementValues,
                            tElementalIQINames );

                    for( uint iElemField = 0; iElemField < tNumElementalIQIs; iElemField++ )
                    {
                        // get the elemental field name
                        std::string tFieldName = tElementalFieldNames( iElemField );

                        // get elemental field values
                        Matrix< DDRMat > tFieldValues = tElementValues.get_column( iElemField );

                        // write elemental field
                        mWriter( aVisMeshIndex )->write_elemental_field(
                                mOutputData( aVisMeshIndex ).mSetNames( Ii ),
                                tFieldName,
                                tFieldValues );
                    }
                }
            }

            for( uint iNodalField = 0; iNodalField < tNumNodalIQIs; iNodalField++ )
            {
                // get the elemental field name
                std::string tFieldName = tNodalFieldNames( iNodalField );

                // get elemental field values
                Matrix< DDRMat > tFieldValues = tNodalValues.get_column( iNodalField );

                // write nodal field
                mWriter( aVisMeshIndex )->write_nodal_field( tFieldName, tFieldValues );
            }

            Matrix< DDRMat > tGlobalVariableValues( tNumGlobalIQIs, 1, MORIS_REAL_MAX );
            for( uint iGlobalField = 0; iGlobalField < tNumGlobalIQIs; iGlobalField++ )
            {
                real tTotalGlobalValue;

                // get the global field name
                std::string tFieldName = tGlobalFieldNames( iGlobalField );

                // get global field values
                tTotalGlobalValue = sum_all( tGlobalValues( iGlobalField ) );

                tGlobalVariableValues( iGlobalField ) = tTotalGlobalValue;
                
                MORIS_LOG_SPEC(tFieldName, tGlobalVariableValues( iGlobalField ));
                // MORIS_LOG_INFO ("Global Variable: %s = %e",tFieldName.c_str(), tGlobalVariableValues( iGlobalField ) );
            }

            // write global variables
            if( tNumGlobalIQIs > 0 )
            {
                mWriter( aVisMeshIndex )->write_global_variables( tGlobalFieldNames, tGlobalVariableValues );
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

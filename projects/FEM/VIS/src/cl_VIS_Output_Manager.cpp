/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_VIS_Output_Manager.cpp
 *
 */

#include "cl_VIS_Output_Manager.hpp"

#include <set>
#include "cl_MTK_Mesh_Manager.hpp"
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

        void
        Output_Manager::set_outputs(
                const uint                       aOutputIndex,
                const enum VIS_Mesh_Type         aMeshType,
                const std::string&               aMeshPath,
                const std::string&               aMeshName,
                const std::string&               aTempPath,
                const std::string&               aTempName,
                const Vector< std::string >&     aBlockNames,
                const Vector< std::string >&     aFieldNames,
                const Vector< enum Field_Type >& aFieldType,
                const Vector< std::string >&     aQINames,
                const uint                       aSaveFrequency,
                const real                       aTimeOffset )
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
            uint tSize          = mOutputData.size();
            uint OutputDataSize = std::max( tSize, aOutputIndex + 1 );

            mOutputData.resize( OutputDataSize );

            // assign output data object to list
            mOutputData( aOutputIndex ) = tOutputData;

            // resize mesh list
            mVisMesh.resize( mOutputData.size(), nullptr );

            mVisMeshCreatedAndOpen.resize( mOutputData.size(), false );
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Output_Manager::set_outputs( moris::ParameterList aParameterlist )
        {
            // create output data object
            vis::Output_Data tOutputData;

            // fill output data object
            tOutputData.mMeshIndex = aParameterlist.get< moris::sint >( "Output_Index" );
            tOutputData.mMeshType  = static_cast< moris::vis::VIS_Mesh_Type >( aParameterlist.get< moris::uint >( "Mesh_Type" ) );

            tOutputData.mOutputPath = std::get< 0 >( aParameterlist.get< std::pair< std::string, std::string > >( "File_Name" ) );
            tOutputData.mMeshName   = std::get< 1 >( aParameterlist.get< std::pair< std::string, std::string > >( "File_Name" ) );

            // note: file path for temp file currently ignored
            tOutputData.mTempPath = std::get< 0 >( aParameterlist.get< std::pair< std::string, std::string > >( "Temp_Name" ) );
            tOutputData.mTempName = std::get< 1 >( aParameterlist.get< std::pair< std::string, std::string > >( "Temp_Name" ) );

            tOutputData.mSaveFrequency = aParameterlist.get< moris::sint >( "Save_Frequency" );
            tOutputData.mTimeOffset    = aParameterlist.get< moris::real >( "Time_Offset" );

            // read and check mesh set names
            Vector< std::string > tSetNames;
            string_to_cell(
                    aParameterlist.get< std::string >( "Set_Names" ),
                    tSetNames );

            MORIS_ERROR( tSetNames.size() > 0,
                    "Output_Manager::set_outputs() - At least one mesh set name needs to be provided for Vis mesh." );

            for ( auto tName : tSetNames )
            {
                MORIS_ERROR( tName.length() > 0,
                        "Output_Manager::set_outputs() - Empty strings for set names in Vis mesh are not allowed." );
            }

            tOutputData.mSetNames = tSetNames;

            // read and check field names
            Vector< std::string > tFieldNames;
            string_to_cell(
                    aParameterlist.get< std::string >( "Field_Names" ),
                    tFieldNames );

            MORIS_ERROR( tFieldNames.size() > 0,
                    "Output_Manager::set_outputs() - At least one field name needs to be provided for Vis mesh." );

            for ( auto tName : tFieldNames )
            {
                MORIS_ERROR( tName.length() > 0,
                        "Output_Manager::set_outputs() - Empty strings for field names in Vis mesh are not allowed." );
            }

            tOutputData.mFieldNames = tFieldNames;

            // read and check field types
            Vector< enum vis::Field_Type >                  tFieldTypes;
            moris::map< std::string, enum vis::Field_Type > tFieldTypeMap = get_vis_field_type_map();
            string_to_cell(
                    aParameterlist.get< std::string >( "Field_Type" ),
                    tFieldTypes,
                    tFieldTypeMap );

            MORIS_ERROR( tFieldTypes.size() > 0,
                    "Output_Manager::set_outputs() - At least one field type needs to be provided for Vis mesh." );

            for ( auto tName : tFieldNames )
            {
                MORIS_ERROR( tName.length() > 0,
                        "Output_Manager::set_outputs() - Empty strings for field types in Vis mesh are not allowed." );
            }

            tOutputData.mFieldType = tFieldTypes;

            // check that length of Field_Names and Field_Type are consistent
            MORIS_ERROR( tFieldNames.size() == tFieldTypes.size(),
                    "Output_Manager::set_outputs() - Number of Field Names and Field Types differ." );

            // read and check IQI names
            Vector< std::string > tQINames;
            string_to_cell(
                    aParameterlist.get< std::string >( "IQI_Names" ),
                    tQINames );

            MORIS_ERROR( tQINames.size() > 0,
                    "Output_Manager::set_outputs() - At least one IQI name needs to be provided for Vis mesh\n" );

            for ( auto tName : tQINames )
            {
                MORIS_ERROR( tName.length() > 0,
                        "Output_Manager::set_outputs() - Empty strings for IQI name in Vis mesh are not allowed\n" );
            }

            tOutputData.mQINames = tQINames;

            // check that length of Field_Names and Field_Type are consistent
            MORIS_ERROR( tFieldNames.size() == tQINames.size(),
                    "Output_Manager::set_outputs() - Number of Field Names and QI Names differ." );

            // resize list of output data objects
            sint tSize          = mOutputData.size();
            sint OutputDataSize = std::max( tSize, tOutputData.mMeshIndex + 1 );

            mOutputData.resize( OutputDataSize );

            // assign output data object to list
            mOutputData( tOutputData.mMeshIndex ) = tOutputData;

            // resize mesh list
            mVisMesh.resize( mOutputData.size(), nullptr );

            mVisMeshCreatedAndOpen.resize( mOutputData.size(), false );
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Output_Manager::setup_vis_mesh_for_output(
                const uint                             aVisMeshIndex,
                std::shared_ptr< mtk::Mesh_Manager >   aMesh,
                const uint                             aMeshPairIndex,
                std::shared_ptr< MSI::Equation_Model > aEquationModel )
        {
            if ( mVisMeshCreatedAndOpen( aVisMeshIndex ) == false )
            {
                Tracer tTracer( "VIS", "Output Manager", "Setup VIS mesh" );

                this->create_visualization_mesh( aVisMeshIndex, aMesh, aMeshPairIndex );

                this->set_visualization_sets( aVisMeshIndex, aEquationModel );

                this->write_mesh( aVisMeshIndex );

                mVisMeshCreatedAndOpen( aVisMeshIndex ) = true;
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Output_Manager::create_visualization_mesh(
                const uint                           aVisMeshIndex,
                std::shared_ptr< mtk::Mesh_Manager > aMesh,
                const uint                           aMeshPairIndex )
        {
            Tracer tTracer( "VIS", "Output Manager", "Create VIS mesh" );

            MORIS_ERROR( mOutputData( aVisMeshIndex ).mMeshIndex == (sint)aVisMeshIndex,
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

        void
        Output_Manager::set_visualization_sets(
                const uint                             aVisMeshIndex,
                std::shared_ptr< MSI::Equation_Model > aEquationModel )
        {
            Tracer tTracer( "VIS", "Output Manager", "Create Sets" );

            // get number of requested sets
            uint tNumRequestedSets = mOutputData( aVisMeshIndex ).mSetNames.size();

            // get the names of the requested sets
            Vector< std::string > const & tMeshSetNames = mOutputData( aVisMeshIndex ).mSetNames;

            // get mtk set index to fem set index map
            map< std::tuple< moris_index, bool, bool >, moris_index >& tMeshSetToFemSetMap =
                    aEquationModel->get_mesh_set_to_fem_set_index_map();

            // get equation sets
            Vector< MSI::Equation_Set* > tEquationSets = aEquationModel->get_equation_sets();

            // get the interpolation and integration meshes
            mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
            mMTKMesh->get_mesh_pair( mMTKMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

            // loop over requested equation sets
            for ( uint iSet = 0; iSet < tNumRequestedSets; iSet++ )
            {
                // get the current set's name
                std::string tSetName = tMeshSetNames( iSet );

                // get FEM/MTK set index
                moris_index tFemSetIndex = tIntegrationMesh->get_set_index_by_name( tSetName );

                if ( tMeshSetToFemSetMap.key_exists( std::make_tuple( tFemSetIndex, false, false ) ) )
                {
                    // find set index for this block index
                    moris_index tEquationSetIndex = tMeshSetToFemSetMap.find( std::make_tuple( tFemSetIndex, false, false ) );

                    // get VIS set index
                    mtk::Set* tVisSet = mVisMesh( aVisMeshIndex )->get_set_by_name( tSetName );

                    // set vis set to fem set. +1 because 0 is reserved for fem
                    tEquationSets( tEquationSetIndex )->set_visualization_set( aVisMeshIndex, tVisSet, mOnlyPrimary );
                }

                // warn the user that the current set name (probably specified in the VIS parameter list) will be ignored as it does not exist
                else
                {
                    MORIS_LOG_WARNING(
                            "Ignoring set with name '%s' has been requested for output on VIS-mesh #%i "
                            "as it does not exist in the list of FEM-sets",
                            tMeshSetNames( iSet ).c_str(),
                            aVisMeshIndex );
                }
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Output_Manager::write_mesh( const uint aVisMeshIndex )
        {
            Tracer tTracer( "VIS", "Output Manager", "Write mesh to file" );

            // specify file path
            std::string tMeshFilePath = mOutputData( aVisMeshIndex ).mMeshPath;

            // get file name
            std::string tMeshFileName = mOutputData( aVisMeshIndex ).mMeshName;

            // specify file path for temporary file
            std::string tMeshTempPath = mOutputData( aVisMeshIndex ).mTempPath;

            // get file name of temporary file
            std::string tMeshTempName = mOutputData( aVisMeshIndex ).mTempName;

            // augment file name if time offset > 0
            if ( mOutputData( aVisMeshIndex ).mTimeOffset > 0 )
            {
                // get optimization iteration
                uint tOptIter = gLogger.get_opt_iteration();

                // set name
                std::string tOptIterStr = std::to_string( tOptIter );
                tMeshFileName += ".e-s." + std::string( 4 - tOptIterStr.length(), '0' ) + tOptIterStr;

                // determine time shift
                mTimeShift = tOptIter * mOutputData( aVisMeshIndex ).mTimeOffset;
            }

            // Log mesh writing message
            MORIS_LOG( "Writing %s to %s.", tMeshFileName.c_str(), tMeshFilePath.c_str() );

            // writes the exodus mesh to file (only the mesh as a geometric entity without any output data)
            mWriter( aVisMeshIndex )->write_mesh( tMeshFilePath, tMeshFileName, tMeshTempPath, tMeshTempName );

            // add nodal elemental and global fields to mesh (only tell the mesh they are there, it doesn't populate them yet)
            this->add_nodal_fields( aVisMeshIndex );
            this->add_elemental_fields( aVisMeshIndex );
            this->add_faceted_fields( aVisMeshIndex );
            this->add_global_fields( aVisMeshIndex );

            // write standard outputs like IDs and Indices to file
            this->write_mesh_indices( aVisMeshIndex );

            // reset field write counter
            mOutputData( aVisMeshIndex ).mFieldWriteCounter = 0;
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Output_Manager::add_nodal_fields( const uint aVisMeshIndex )
        {
            // get the field names
            Vector< std::string > const & tFieldNames = mOutputData( aVisMeshIndex ).mFieldNames;

            // set list of nodal field names to input + 2
            Vector< std::string > tNodalFieldNames( 2 + tFieldNames.size() );

            // set standard field names
            tNodalFieldNames( 0 ) = "Node_Id";
            tNodalFieldNames( 1 ) = "Node_Index";

            uint tCounter = 2;

            // loop over field names and check if fields are nodal fields
            for ( uint iField = 0; iField < tFieldNames.size(); iField++ )
            {
                if ( mOutputData( aVisMeshIndex ).mFieldType( iField ) == Field_Type::NODAL )
                {
                    tNodalFieldNames( tCounter++ ) = tFieldNames( iField );
                }
            }

            tNodalFieldNames.resize( tCounter );

            // pass nodal field names to writer
            mWriter( aVisMeshIndex )->set_nodal_fields( tNodalFieldNames );
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Output_Manager::add_elemental_fields( const uint aVisMeshIndex )
        {
            // get the field names
            const Vector< std::string > tFieldNames = mOutputData( aVisMeshIndex ).mFieldNames;
            uint                        tNumFields  = tFieldNames.size();

            // allocate cell for storing elemental field names; 3 default fields are added
            Vector< std::string > tElementalFieldNames( 3 + tNumFields );

            tElementalFieldNames( 0 ) = "Mesh_Id";
            tElementalFieldNames( 1 ) = "Mesh_Index";
            tElementalFieldNames( 2 ) = "Proc_Index";

            // set field counter to 3 to account for default fields
            uint tCounter = 3;

            // loop over field names and check if fields are elemental fields
            for ( uint iField = 0; iField < tNumFields; iField++ )
            {
                vis::Field_Type tOutputFieldType = mOutputData( aVisMeshIndex ).mFieldType( iField );

                if ( tOutputFieldType == Field_Type::ELEMENTAL_INT || tOutputFieldType == Field_Type::ELEMENTAL_AVG )
                {
                    tElementalFieldNames( tCounter++ ) = tFieldNames( iField );
                }
            }

            // trim cell of element field names
            tElementalFieldNames.resize( tCounter );

            // write field names to file
            mWriter( aVisMeshIndex )->set_elemental_fields( tElementalFieldNames );
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Output_Manager::add_faceted_fields( const uint aVisMeshIndex )
        {
            // allocate cell for storing elemental field names; 3 default fields are added
            Vector< std::string > tFacetedFieldNames( mOutputData( aVisMeshIndex ).mFieldNames.size() );

            // initialize counter
            uint tCounter = 0;

            // loop over field names and check if fields are elemental fields
            for ( uint iField = 0; iField < mOutputData( aVisMeshIndex ).mFieldNames.size(); iField++ )
            {
                vis::Field_Type tOutputFieldType = mOutputData( aVisMeshIndex ).mFieldType( iField );

                if ( tOutputFieldType == Field_Type::FACETED_INT || tOutputFieldType == Field_Type::FACETED_AVG )
                {
                    tFacetedFieldNames( tCounter++ ) = mOutputData( aVisMeshIndex ).mFieldNames( iField );
                }
            }

            // trim cell of element field names
            tFacetedFieldNames.resize( tCounter );

            // write field names to file
            mWriter( aVisMeshIndex )->set_side_set_fields( tFacetedFieldNames );
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Output_Manager::add_global_fields( const uint aVisMeshIndex )
        {
            Vector< std::string > tGlobalFieldNames( mOutputData( aVisMeshIndex ).mFieldNames.size() );

            uint tCounter = 0;

            // loop over field names and check if fields are global fields
            for ( uint iField = 0; iField < mOutputData( aVisMeshIndex ).mFieldNames.size(); iField++ )
            {
                if ( mOutputData( aVisMeshIndex ).mFieldType( iField ) == Field_Type::GLOBAL )
                {
                    tGlobalFieldNames( tCounter++ ) = mOutputData( aVisMeshIndex ).mFieldNames( iField );
                }
            }

            tGlobalFieldNames.resize( tCounter );

            mWriter( aVisMeshIndex )->set_global_variables( tGlobalFieldNames );
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Output_Manager::write_mesh_indices( const uint aVisMeshIndex )
        {
            // ---------------------------------
            // write nodal mesh info

            // getting the vertices on the vis mesh
            Vector< moris::mtk::Vertex const * > tVertices = mVisMesh( aVisMeshIndex )->get_all_vertices();

            // get number of verts
            uint tNumVerts = tVertices.size();

            // create cell of index, id field, and proc indices
            Vector< Matrix< DDRMat > > tVertIdIndex( 2 );
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
            mWriter( aVisMeshIndex )->write_nodal_field( "Node_Id", tVertIdIndex( 0 ) );
            mWriter( aVisMeshIndex )->write_nodal_field( "Node_Index", tVertIdIndex( 1 ) );

            // ---------------------------------
            // write elemental mesh info

            // get mesh set indices
            uint tRequestedSets = mOutputData( aVisMeshIndex ).mSetNames.size();

            // loop over mesh sets
            for ( uint iSet = 0; iSet < tRequestedSets; iSet++ )
            {
                // get vis set by index
                moris::mtk::Set* tSet = mVisMesh( aVisMeshIndex )->get_set_by_index( iSet );

                bool tOnlyPrimaryCells = true;

                switch ( mOutputData( aVisMeshIndex ).mMeshType )
                {
                    case vis::VIS_Mesh_Type::STANDARD:
                        tOnlyPrimaryCells = true;
                        break;

                    case vis::VIS_Mesh_Type::STANDARD_WITH_OVERLAP:
                        tOnlyPrimaryCells = false;
                        break;

                    case vis::VIS_Mesh_Type::FULL_DISCONTINUOUS:
                        tOnlyPrimaryCells = true;
                        break;

                    case vis::VIS_Mesh_Type::FULL_DISCONTINUOUS_WITH_OVERLAP:
                        tOnlyPrimaryCells = false;
                        break;

                    default:
                        MORIS_ERROR( false, "Output_Manager::create_visualization_mesh() - Mesh type not specified or unknown." );
                        break;
                }

                // don't write indices on facets
                if ( tSet->get_set_type() != mtk::SetType::BULK )
                {
                    continue;
                }

                // get number of cells on set
                uint tNumCells = tSet->get_num_cells_on_set( tOnlyPrimaryCells );

                // check whether number of cells > 0; otherwise skip remainder
                if ( tNumCells == 0 )
                {
                    continue;
                }

                // get cell indices on set
                moris::Matrix< DDSMat > tCellIndex = tSet->get_cell_inds_on_block( tOnlyPrimaryCells );

                // find the maximal index for resizing purposes
                sint tMaxIndex = tCellIndex.max();

                // create cell assembly map( index to position )
                Matrix< DDSMat > tCellAssemblyMap( tMaxIndex + 1, 1, -1 );

                // loop over cells and put them in map
                for ( uint iCell = 0; iCell < tNumCells; iCell++ )
                {
                    tCellAssemblyMap( tCellIndex( iCell ) ) = iCell;
                }

                // create cell of index, id field, and proc indices
                Vector< Matrix< DDRMat > > tIdIndex( 3 );
                tIdIndex( 0 ).set_size( tCellIndex.numel(), 1 );
                tIdIndex( 1 ).set_size( tCellIndex.numel(), 1 );
                tIdIndex( 2 ).set_size( tCellIndex.numel(), 1 );

                // get clusters from vis set
                Vector< mtk::Cluster const * > tMeshClusterList = tSet->get_clusters_on_set();

                // loop over clusters and get ids and indices
                for ( uint iCluster = 0; iCluster < tMeshClusterList.size(); iCluster++ )
                {
                    // get primary cells
                    const Vector< moris::mtk::Cell const * >& tPrimaryCells =
                            tMeshClusterList( iCluster )->get_primary_cells_in_cluster();

                    // loop over primary cells
                    for ( uint iPrimaryCell = 0; iPrimaryCell < tPrimaryCells.size(); iPrimaryCell++ )
                    {
                        // get index of vis cell
                        moris_index tIndex = tPrimaryCells( iPrimaryCell )->get_index();

                        moris_id tMeshId =
                                reinterpret_cast< const vis::Cell_Visualization* >( tPrimaryCells( iPrimaryCell ) )->get_mesh_cell_id();

                        moris_index tMeshIndex =
                                reinterpret_cast< const vis::Cell_Visualization* >( tPrimaryCells( iPrimaryCell ) )->get_mesh_cell_index();

                        tIdIndex( 0 )( tCellAssemblyMap( tIndex ) ) = tMeshId;
                        tIdIndex( 1 )( tCellAssemblyMap( tIndex ) ) = tMeshIndex;
                        tIdIndex( 2 )( tCellAssemblyMap( tIndex ) ) = par_rank();
                    }

                    const Vector< moris::mtk::Cell const * >& tVoidCells = tMeshClusterList( iCluster )->get_void_cells_in_cluster();

                    for ( uint iVoidCell = 0; iVoidCell < tVoidCells.size(); iVoidCell++ )
                    {
                        moris_index tIndex = tVoidCells( iVoidCell )->get_index();

                        moris_id tMeshId =
                                reinterpret_cast< const vis::Cell_Visualization* >( tVoidCells( iVoidCell ) )->get_mesh_cell_id();

                        moris_index tMeshIndex =
                                reinterpret_cast< const vis::Cell_Visualization* >( tVoidCells( iVoidCell ) )->get_mesh_cell_index();

                        tIdIndex( 0 )( tCellAssemblyMap( tIndex ) ) = tMeshId;
                        tIdIndex( 1 )( tCellAssemblyMap( tIndex ) ) = tMeshIndex;
                        tIdIndex( 2 )( tCellAssemblyMap( tIndex ) ) = par_rank();
                    }

                }    // end for: each cluster on set

                const std::string tSetName = tSet->get_set_name();
                mWriter( aVisMeshIndex )->write_elemental_field( tSetName, "Mesh_Id", tIdIndex( 0 ) );
                mWriter( aVisMeshIndex )->write_elemental_field( tSetName, "Mesh_Index", tIdIndex( 1 ) );
                mWriter( aVisMeshIndex )->write_elemental_field( tSetName, "Proc_Index", tIdIndex( 2 ) );

            }    // end for: each set
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Output_Manager::write_field(
                const uint                             aVisMeshIndex,
                const real                             aTime,
                std::shared_ptr< MSI::Equation_Model > aEquationModel )
        {
            // log/time this operation
            Tracer tTracer( "VIS", "Output Manager", "Write Fields" );

            // number of set names
            uint tNumSetNames = mOutputData( aVisMeshIndex ).mSetNames.size();

            // initialize lists of IQIs and their output field names
            Vector< Vector< std::string > > tIQINames;
            Vector< Vector< std::string > > tFieldNames;
            Vector< uint >                  tNumIQIsForFieldType;
            this->get_IQI_and_field_names(
                    aVisMeshIndex,
                    tIQINames,
                    tFieldNames,
                    tNumIQIsForFieldType );

            // increment field write counter (signal that another output has been performed)
            mOutputData( aVisMeshIndex ).mFieldWriteCounter++;

            // write time to file
            mWriter( aVisMeshIndex )->set_time( aTime + mTimeShift );

            // get mesh set to fem set index map
            map< std::tuple< moris_index, bool, bool >, moris_index >& tMeshSetToFemSetMap =
                    aEquationModel->get_mesh_set_to_fem_set_index_map();

            // get equation sets
            Vector< MSI::Equation_Set* >& tEquationSets = aEquationModel->get_equation_sets();

            // access the integration mesh for output
            mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
            mMTKMesh->get_mesh_pair( mMTKMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

            // get number of global and nodal IQIs
            uint tNumGlobalIQIs = tNumIQIsForFieldType( (uint)Field_Type::GLOBAL );
            uint tNumNodalIQIs  = tNumIQIsForFieldType( (uint)Field_Type::NODAL );

            // initialize list of nodal field values to be filled below
            Matrix< DDRMat > tNodalValues(
                    mVisMesh( aVisMeshIndex )->get_num_nodes(),
                    tNumNodalIQIs,
                    std::numeric_limits< real >::quiet_NaN() );

            // initialize list of global field values to be filled below
            Matrix< DDRMat > tGlobalValues( 1, tNumGlobalIQIs, 0.0 );

            // loop over all sets on this VIS mesh
            for ( uint iSet = 0; iSet < tNumSetNames; iSet++ )
            {
                // get the current set's name
                std::string tSetName = mOutputData( aVisMeshIndex ).mSetNames( iSet );

                // get mesh set index from name
                moris_index tFemSetIndex = tIntegrationMesh->get_set_index_by_name( tSetName );

                // only output sets that actually exist in equation model
                if ( !tMeshSetToFemSetMap.key_exists( std::make_tuple( tFemSetIndex, false, false ) ) )
                {
                    continue;
                }

                // find set index for this block index
                moris_index tEquationSetIndex =
                        tMeshSetToFemSetMap.find( std::make_tuple( tFemSetIndex, false, false ) );

                // get access to the FEM set
                MSI::Equation_Set* tFemSet = tEquationSets( tEquationSetIndex );

                // skip the set if it is empty on the current proc
                if ( tFemSet->get_num_equation_objects() == 0 )
                {
                    continue;
                }

                // compute the fields for the current set
                this->compute_fields_for_set(
                        aVisMeshIndex,
                        tFemSet,
                        tIQINames,
                        tFieldNames,
                        &tGlobalValues,
                        &tNodalValues );

            }    // end for: each set requested for output mesh

            for ( uint iNodalField = 0; iNodalField < tNumNodalIQIs; iNodalField++ )
            {
                // get the elemental field name
                std::string tFieldName = tFieldNames( (uint)Field_Type::NODAL )( iNodalField );

                // get elemental field values
                Matrix< DDRMat > tFieldValues = tNodalValues.get_column( iNodalField );

                // write nodal field
                mWriter( aVisMeshIndex )->write_nodal_field( tFieldName, tFieldValues );
            }

            Matrix< DDRMat > tGlobalVariableValues( tNumGlobalIQIs, 1, MORIS_REAL_MAX );
            for ( uint iGlobalField = 0; iGlobalField < tNumGlobalIQIs; iGlobalField++ )
            {
                // get the global field name
                std::string tFieldName = tFieldNames( (uint)Field_Type::GLOBAL )( iGlobalField );

                // get global field values
                real tTotalGlobalValue = sum_all( tGlobalValues( iGlobalField ) );

                // store global value
                tGlobalVariableValues( iGlobalField ) = tTotalGlobalValue;

                // write global values to console
                MORIS_LOG_SPEC( tFieldName, tGlobalVariableValues( iGlobalField ) );
            }

            // write global variables to exodus
            if ( tNumGlobalIQIs > 0 )
            {
                mWriter( aVisMeshIndex )->write_global_variables( tFieldNames( (uint)Field_Type::GLOBAL ), tGlobalVariableValues );
            }

            // check if a copy of the current mesh file should be created
            sint tFieldWriteCounter = mOutputData( aVisMeshIndex ).mFieldWriteCounter;
            sint tSaveFrequency     = mOutputData( aVisMeshIndex ).mSaveFrequency;

            if ( std::remainder( tFieldWriteCounter, tSaveFrequency ) == 0 )
            {
                mWriter( aVisMeshIndex )->save_mesh();
            }
        }

        //-----------------------------------------------------------------------------------------------------------

        void
        Output_Manager::get_IQI_and_field_names(
                const uint                       aVisMeshIndex,
                Vector< Vector< std::string > >& aIQINames,
                Vector< Vector< std::string > >& aFieldNames,
                Vector< uint >&                  aNumIQIsForFieldType )
        {
            // number of fields in vis mesh
            uint tNumFields = mOutputData( aVisMeshIndex ).mFieldNames.size();

            /* initialize lists of IQIs and their output field names
             * input: Field_Type converted to integer
             * output: list of IQIs and their field names for this field type */
            aIQINames.resize( (uint)Field_Type::END_ENUM );
            aFieldNames.resize( (uint)Field_Type::END_ENUM );
            aIQINames.reserve( 2 * tNumFields );
            aFieldNames.reserve( 2 * tNumFields );

            // loop over all output fields and sort them into lists by type (i.e. global, nodal, or elemental)
            for ( uint iField = 0; iField < tNumFields; iField++ )
            {
                // get the type of the current field
                uint tFieldType = (uint)mOutputData( aVisMeshIndex ).mFieldType( iField );

                // get the IQI name
                std::string tIQIName = mOutputData( aVisMeshIndex ).mQINames( iField );

                // get the field name
                std::string tFieldName = mOutputData( aVisMeshIndex ).mFieldNames( iField );

                // sort the IQIs and their corresponding field name into the lists
                aIQINames( tFieldType ).push_back( tIQIName );
                aFieldNames( tFieldType ).push_back( tFieldName );
            }

            // figure out how many fields there are of each type
            aNumIQIsForFieldType.resize( (uint)Field_Type::END_ENUM );

            for ( uint iFieldType = 0; iFieldType < (uint)Field_Type::END_ENUM; iFieldType++ )
            {
                aNumIQIsForFieldType( iFieldType ) = aIQINames( iFieldType ).size();
            }

        }    // end function: VIS::Output_Manager::get_IQI_and_field_names()

        //-----------------------------------------------------------------------------------------------------------

        void
        Output_Manager::compute_fields_for_set(
                const uint                              aVisMeshIndex,
                MSI::Equation_Set*                      aFemSet,
                Vector< Vector< std::string > > const & aIQINames,
                Vector< Vector< std::string > > const & aFieldNames,
                Matrix< DDRMat >*                       aGlobalFieldValues,
                Matrix< DDRMat >*                       aNodalFieldValues )
        {
            // compute fields for each type
            for ( uint iFieldType = 0; iFieldType < (uint)Field_Type::END_ENUM; iFieldType++ )
            {
                // get the requested IQIs for the current type
                Vector< std::string > const & tIQINamesForType   = aIQINames( iFieldType );
                Vector< std::string > const & tFieldNamesForType = aFieldNames( iFieldType );

                // get the number of fields for the current type
                uint tNumFieldsForType = tIQINamesForType.size();

                // skip this field type if there aren't any IQIs
                if ( tNumFieldsForType == 0 )
                {
                    continue;
                }

                // get the field type as an enum (for better readability)
                Field_Type tFieldType = (Field_Type)iFieldType;

                // perform check that we only output using supported Field types
                std::set< Field_Type > tUnsupportedFieldTypes = { Field_Type::UNDEFINED, Field_Type::NODAL_IP, Field_Type::END_ENUM };
                MORIS_ERROR( tUnsupportedFieldTypes.find( tFieldType ) == tUnsupportedFieldTypes.end(),
                        "VIS::Output_Manager::compute_fields_for_set() - Trying to output to unsupported field type." );

                // all supported elemental field types are some sort of elemental values
                std::set< Field_Type > tElementalFieldTypes = {
                    Field_Type::ELEMENTAL_AVG,
                    Field_Type::ELEMENTAL_INT,
                    Field_Type::FACETED_AVG,
                    Field_Type::FACETED_INT
                };

                // compute global fields
                if ( tFieldType == Field_Type::GLOBAL )
                {
                    aFemSet->compute_quantity_of_interest_global(
                            aVisMeshIndex,
                            aGlobalFieldValues,
                            tIQINamesForType );
                }

                // compute nodal fields, if there are any
                else if ( tFieldType == Field_Type::NODAL )
                {
                    aFemSet->compute_quantity_of_interest_nodal(
                            aVisMeshIndex,
                            aNodalFieldValues,
                            tIQINamesForType );
                }

                // compute elemental fields and output them immediately to avoid unnecessarily storing data
                else if ( tElementalFieldTypes.find( tFieldType ) != tElementalFieldTypes.end() )
                {
                    this->compute_and_write_elemental_fields_on_set(
                            aVisMeshIndex,
                            aFemSet,
                            tFieldType,
                            tIQINamesForType,
                            tFieldNamesForType );
                }

            }    // end for: each field type

        }    // end function: VIS::Output_Manager::compute_fields_for_set()

        //-----------------------------------------------------------------------------------------------------------

        void
        Output_Manager::compute_and_write_elemental_fields_on_set(
                const uint                    aVisMeshIndex,
                MSI::Equation_Set*            aFemSet,
                const Field_Type              aFieldType,
                Vector< std::string > const & aIQINamesForType,
                Vector< std::string > const & aFieldNamesForType )
        {
            // get the number of IQIs for the type
            uint tNumIQIsForType = aIQINamesForType.size();
            MORIS_ASSERT(
                    tNumIQIsForType > 0,
                    "VIS::Output_Manager::compute_and_write_elemental_fields_on_set() - "
                    "No IQIs to compute. This function should only be called if there are IQIs to be computed." );

            // sanity check the  input
#ifdef MORIS_HAVE_DEBUG
            std::set< Field_Type > tSupportedFieldTypes = {
                Field_Type::ELEMENTAL_AVG,
                Field_Type::ELEMENTAL_INT,
                Field_Type::FACETED_AVG,
                Field_Type::FACETED_INT
            };
            MORIS_ASSERT(
                    tSupportedFieldTypes.find( aFieldType ) != tSupportedFieldTypes.end(),
                    "VIS::Output_Manager::compute_and_write_elemental_fields_on_set() - "
                    "Trying to compute non-elemental field type which is not supported in this function." );
#endif

            // decode the field type into FACETED vs. ELEMENTAL and AVG vs. INT
            bool tIsFacetedFieldType  = ( aFieldType == Field_Type::FACETED_AVG || aFieldType == Field_Type::FACETED_INT );
            bool tIsAveragedFieldType = ( aFieldType == Field_Type::FACETED_AVG || aFieldType == Field_Type::ELEMENTAL_AVG );

            // get the current set's name
            const std::string tSetName = aFemSet->get_set_name();

            // get the set's type
            const fem::Element_Type tSetType = aFemSet->get_element_type();

            // check if set is a faceted set
            bool tIsFacetedSet = ( tSetType == fem::Element_Type::SIDESET )
                              || ( tSetType == fem::Element_Type::DOUBLE_SIDESET )
                              || ( tSetType == fem::Element_Type::NONCONFORMAL_SIDESET );

            // skip mis-matched sets, i.e. faceted fields are only outputted on side sets and elemental fields only on block sets
            if ( tIsFacetedFieldType != tIsFacetedSet )
            {
                return;
            }

            // initialize array of elemental values
            Matrix< DDRMat > tElementValues;

            // compute elemental values
            aFemSet->compute_quantity_of_interest_elemental(
                    aVisMeshIndex,
                    &tElementValues,
                    aIQINamesForType,
                    tIsAveragedFieldType );

            // write each elemental field to output immediately
            for ( uint iElemField = 0; iElemField < tNumIQIsForType; iElemField++ )
            {
                // get the elemental field name
                std::string tFieldName = aFieldNamesForType( iElemField );

                // get elemental field values
                Matrix< DDRMat > tFieldValues = tElementValues.get_column( iElemField );

                // write elemental field (write as facet or as elemental values to exodus depending on the set type)
                if ( tIsFacetedFieldType )
                {
                    mWriter( aVisMeshIndex )->write_side_set_field( tSetName, tFieldName, tFieldValues );
                }
                else
                {
                    mWriter( aVisMeshIndex )->write_elemental_field( tSetName, tFieldName, tFieldValues );
                }
            }

        }    // end function: VIS::Output_Manager::compute_and_write_elemental_fields_on_set()

        //-----------------------------------------------------------------------------------------------------------

    }    // namespace vis
}    // namespace moris

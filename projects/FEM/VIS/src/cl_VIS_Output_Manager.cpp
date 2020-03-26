
#include "cl_VIS_Output_Manager.hpp"

#include "cl_VIS_Factory.hpp"

#include "cl_MDL_Model.hpp"
#include "cl_MSI_Equation_Set.hpp"
#include "cl_MSI_Equation_Model.hpp"

#include "cl_FEM_Set.hpp"

#include "fn_Parsing_Tools.hpp"

extern moris::Comm_Manager gMorisComm;

namespace moris
{
    namespace vis
    {

//-----------------------------------------------------------------------------------------------------------

    void Output_Manager::set_outputs( const uint                              aOutputIndex,
                                      const enum VIS_Mesh_Type                aMeshType,
                                      const std::string                     & aMeshPath,
                                      const std::string                     & aMeshName,
                                      const moris::Cell< std::string >      & aBlockNames,
                                      const moris::Cell< std::string >      & aFieldNames,
                                      const moris::Cell< enum Field_Type >  & aFieldType,
                                      const moris::Cell< enum Output_Type > & aEnum )
    {
        // create output data object
        vis::Output_Data tOutputData;

        // fill output data object
        tOutputData.mMeshIndex  = aOutputIndex;
        tOutputData.mMeshType   = aMeshType;
        tOutputData.mMeshName   = aMeshName;
        tOutputData.mMeshPath   = aMeshPath;
        tOutputData.mSetNames   = aBlockNames;
        tOutputData.mFieldNames = aFieldNames;
        tOutputData.mFieldType  = aFieldType;
        tOutputData.mOutputType = aEnum;

        // resize list of output data objects
        uint tSize = mOutputData.size();
        uint OutputDataSize = std::max( tSize, aOutputIndex + 1 );

        mOutputData.resize( OutputDataSize );

        // assign output data object to list
        mOutputData( aOutputIndex ) = tOutputData;

        // resize mesh list
        mVisMesh.resize( mOutputData.size(), nullptr );
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

        moris::Cell< std::string > tSetNames;
        string_to_cell( aParamterelist.get< std::string >( "Set_Names" ), tSetNames );
        tOutputData.mSetNames   = tSetNames;

        moris::Cell< std::string > tFieldNames;
        string_to_cell( aParamterelist.get< std::string >( "Field_Names" ), tFieldNames );
        tOutputData.mFieldNames = tFieldNames;

        moris::Cell< enum vis::Field_Type > tFieldTypes;
        moris::map< std::string, enum vis::Field_Type > tFieldTypeMap = get_vis_field_type_map();
        string_to_cell( aParamterelist.get< std::string >( "Field_Type" ) ,
                        tFieldTypes,
                        tFieldTypeMap );
        tOutputData.mFieldType  = tFieldTypes;

        moris::Cell< enum vis::Output_Type > tOutputTypes;
        moris::map< std::string, enum vis::Output_Type > tOutputTypeMap = get_vis_output_type_map();
        string_to_cell( aParamterelist.get< std::string >( "Output_Type" ) ,
                        tOutputTypes,
                        tOutputTypeMap );
        tOutputData.mOutputType = tOutputTypes;

        // resize list of output data objects
        sint tSize = mOutputData.size();
        sint OutputDataSize = std::max( tSize, tOutputData.mMeshIndex + 1 );

        mOutputData.resize( OutputDataSize );

        // assign output data object to list
        mOutputData( tOutputData.mMeshIndex ) = tOutputData;

        // resize mesh list
        mVisMesh.resize( mOutputData.size(), nullptr );
    }

//-----------------------------------------------------------------------------------------------------------

    void Output_Manager::create_visualization_mesh( const uint                aVisMeshIndex,
                                                          mtk::Mesh_Manager * aMesh,
                                                    const uint                aMeshPairIndex)
    {
        MORIS_ERROR( mOutputData( aVisMeshIndex ).mMeshIndex == ( sint )aVisMeshIndex, "create_visualization_meshes(), Visualization mesh not set" );

        mMTKMesh = aMesh;
        mMTKMeshPairIndex = aMeshPairIndex;

        // create vis factory
        vis::VIS_Factory tVisFactory( aMesh, aMeshPairIndex );

        // create vis mesh
        mVisMesh( aVisMeshIndex ) = tVisFactory.create_visualization_mesh( mOutputData( aVisMeshIndex ) );

        // resize list of writers to list of outputs. memory allocation stays intact
        mWriter.resize( mOutputData.size(), nullptr );

        // create writer for this mesh
        mWriter( aVisMeshIndex ) = new Writer_Exodus( mVisMesh( aVisMeshIndex ) );
    }

//-----------------------------------------------------------------------------------------------------------

    void Output_Manager::set_visualization_sets( const uint                                   aVisMeshIndex,
                                                       std::shared_ptr< MSI::Equation_Model > aEquationModel )
    {
        // get number of requested sets
        uint tNumRequestedSets = mOutputData( aVisMeshIndex ).mSetNames.size();

        // get mtk set index to fem set index map
        map< moris_index, moris_index > & tMeshSetToFemSetMap
        = aEquationModel->get_mesh_set_to_fem_set_index_map();

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

            if ( tMeshSetToFemSetMap.key_exists( tSetIndex ) )
            {
                // find set index for this block index
                moris_index tEquationSetIndex = tMeshSetToFemSetMap.find( tSetIndex );

                // set vis set to fem set. +1 because 0 is reserved for fem
                tEquationSets( tEquationSetIndex )->set_visualization_set( aVisMeshIndex + 1,
                                                                           mVisMesh( aVisMeshIndex )->get_set_by_index( Ii ),
                                                                           mOnlyPrimary );
            }
        }
    }

//-----------------------------------------------------------------------------------------------------------

    void Output_Manager::write_mesh( const uint aVisMeshIndex,
                                     const real tTime )
    {
        // specify file path
        std::string tMeshFilePath = mOutputData( aVisMeshIndex ).mMeshPath;

        // get file name
        std::string tMeshFileName = mOutputData( aVisMeshIndex ).mMeshName;

        std::string tMassage = "Writing " + tMeshFileName + " to " + tMeshFilePath +". \n";

        MORIS_LOG( tMassage.c_str() );

        // write mesh to file
        mWriter( aVisMeshIndex )->write_mesh( tMeshFilePath, tMeshFileName );

        // add nodal elemental and global fields to mesh
        this->add_nodal_fields( aVisMeshIndex );
        this->add_elemetal_fields( aVisMeshIndex );
        this->add_global_fields( aVisMeshIndex );

        // write time to file
        mWriter( aVisMeshIndex )->set_time( tTime );

        // write standard outputs like IDs and Indices to file
        //this->write_mesh_indices( aVisMeshIndex );
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
        moris::Cell<std::string> tElementalFieldNames( 2 + mOutputData( aVisMeshIndex ).mFieldNames.size() );

        tElementalFieldNames( 0 ) = "Mesh_Id";
        tElementalFieldNames( 1 ) = "Mesh_Index";

        uint tCounter = 2;

        // loop over field names and check if fields are elemental fields
        for( uint Ik = 0; Ik < mOutputData( aVisMeshIndex ).mFieldNames.size(); Ik++ )
        {
            if( mOutputData( aVisMeshIndex ).mFieldType( Ik ) == Field_Type::ELEMENTAL )
            {
                tElementalFieldNames( tCounter++ ) = mOutputData( aVisMeshIndex ).mFieldNames( Ik );
            }
        }

        tElementalFieldNames.resize( tCounter );

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
                case ( vis::VIS_Mesh_Type::STANDARD ):
                    tOnlyPrimaryCells = true ;
                    break;

                case ( vis::VIS_Mesh_Type::OVERLAPPING_INTERFACE ):
                    tOnlyPrimaryCells = false;
                    break;

                case ( vis::VIS_Mesh_Type::FULL_DISCONTINOUS ):
                     MORIS_ERROR( false, "create_visualization_mesh() - Mesh type FULL_DISCONTINOUS not implemented yet. " );
                     break;

                default:
                    MORIS_ERROR( false, "create_visualization_mesh() - Mesh type not specified. " );
                    break;
            }

            // get number of cells on set
            uint tNumCells = tSet->get_num_cells_on_set( tOnlyPrimaryCells );

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

            // create cell of inde and id field
            moris::Cell< Matrix< DDRMat > > tIdIndex( 2 );
            tIdIndex( 0 ).set_size( tCellIndex.numel(), 1 );
            tIdIndex( 1 ).set_size( tCellIndex.numel(), 1 );

            // get clusters from vis set
            moris::Cell< mtk::Cluster const* > tMeshClusterList = tSet->get_clusters_on_set();

            // loop over clusters and get ids and indices
            for( uint Ik = 0; Ik < tMeshClusterList.size(); Ik++ )
            {
                // get primary cells
                const moris::Cell<moris::mtk::Cell const *> & tPrimaryCells = tMeshClusterList( Ik )->get_primary_cells_in_cluster();

                // loop over primary cells
                for( uint Ia = 0; Ia < tPrimaryCells.size(); Ia++ )
                {
                    // get index of vis cell
                    moris_index tIndex = tPrimaryCells( Ia )->get_index();

                    moris_id tMeshId =  reinterpret_cast< const vis::Cell_Visualization* >( tPrimaryCells( Ia ) )->get_mesh_cell_id();

                    moris_index tMeshIndex = reinterpret_cast< const vis::Cell_Visualization* >( tPrimaryCells( Ia ) )->get_mesh_cell_index();

                    tIdIndex( 0 )( tCellAsseblyMap( tIndex ) ) = tMeshId;
                    tIdIndex( 1 )( tCellAsseblyMap( tIndex ) ) = tMeshIndex;
                }

                const moris::Cell<moris::mtk::Cell const *> & tVoidCells = tMeshClusterList( Ik )->get_void_cells_in_cluster();

                for( uint Ia = 0; Ia < tVoidCells.size(); Ia++ )
                {
                    moris_index tIndex = tVoidCells( Ia )->get_index();

                    moris_id tMeshId = reinterpret_cast< const vis::Cell_Visualization* >( tVoidCells( Ia ) )->get_mesh_cell_id();

                    moris_index tMeshIndex = reinterpret_cast< const vis::Cell_Visualization* >( tVoidCells( Ia ) )->get_mesh_cell_index();

                    tIdIndex( 0 )( tCellAsseblyMap( tIndex ) ) = tMeshId;
                    tIdIndex( 1 )( tCellAsseblyMap( tIndex ) ) = tMeshIndex;
                }
            }

            mWriter( aVisMeshIndex )->write_elemental_field( tSet->get_set_name() , "Mesh_Id", tIdIndex( 0 ) );
            mWriter( aVisMeshIndex )->write_elemental_field( tSet->get_set_name() , "Mesh_Index", tIdIndex( 1 ) );
        }
    }

//-----------------------------------------------------------------------------------------------------------
    void Output_Manager::write_field( const uint                                   aVisMeshIndex,
                                            std::shared_ptr< MSI::Equation_Model > aEquationModel )
    {
        // get mesh set to fem set index map
        map< moris_index, moris_index > & tMeshSetToFemSetMap
        = aEquationModel->get_mesh_set_to_fem_set_index_map( );

        // get equation sets
        moris::Cell< MSI::Equation_Set * > tEquationSets = aEquationModel->get_equation_sets();

        // get integration mesh
        mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
        mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
        mMTKMesh->get_mesh_pair( mMTKMeshPairIndex, tInterpolationMesh, tIntegrationMesh );


        // loop over all fields of this output object
        for( uint Ik = 0; Ik < mOutputData( aVisMeshIndex ).mFieldNames.size(); Ik++ )
        {
            // get field name
            std::string tFieldName = mOutputData( aVisMeshIndex ).mFieldNames( Ik );

            // nodal and global field vaues
            Matrix< DDRMat > tNodalValues( mVisMesh( aVisMeshIndex )->get_num_nodes(), 1, 0.0 );   //FIXME replace default with NaN
            moris::real      tGlobalValue;

            // loop over all blocks on this output object
            for( uint Ii = 0; Ii < mOutputData( aVisMeshIndex ).mSetNames.size(); Ii++ )
            {
                moris_index tSetIndex = tIntegrationMesh->get_set_index_by_name( mOutputData( aVisMeshIndex ).mSetNames( Ii ) );

                if ( tMeshSetToFemSetMap.key_exists( tSetIndex ) )
                {
                    // find set index for this block index
                    moris_index tEquationSetIndex = tMeshSetToFemSetMap.find( tSetIndex );

                    // elemental field values
                    Matrix< DDRMat > tElementValues;

                    tEquationSets( tEquationSetIndex )->compute_quantity_of_interest( aVisMeshIndex + 1,
                                                                                      &tElementValues,
                                                                                      &tNodalValues,
                                                                                      &tGlobalValue,
                                                                                      mOutputData( aVisMeshIndex ).mOutputType( Ik ),
                                                                                      mOutputData( aVisMeshIndex ).mFieldType( Ik ) );

                    if( mOutputData( aVisMeshIndex ).mFieldType( Ik ) == Field_Type::ELEMENTAL )
                    {
                        mWriter( aVisMeshIndex )->write_elemental_field( mOutputData( aVisMeshIndex ).mSetNames( Ii ), tFieldName, tElementValues );
                    }
                }
            }

            if( mOutputData( aVisMeshIndex ).mFieldType( Ik ) == Field_Type::NODAL )
            {
                mWriter( aVisMeshIndex )->write_nodal_field( tFieldName, tNodalValues );
            }
            else if ( mOutputData( aVisMeshIndex ).mFieldType( Ik ) == Field_Type::GLOBAL )
            {
                mWriter( aVisMeshIndex )->write_global_variable( tFieldName, tGlobalValue );
            }
        }
    }

//-----------------------------------------------------------------------------------------------------------

    }
}

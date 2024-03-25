/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * main.cpp
 *
 */

#include "cl_Stopwatch.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "moris_typedefs.hpp" // COR/src

#include "cl_Vector.hpp"

#include "banner.hpp" // COR/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_save_matrix_to_binary_file.hpp"

#include "cl_MTK_Mesh_Factory.hpp"

#include "SDF_Tools.hpp"

#include "cl_SDF_State.hpp"
#include "cl_SDF_Arguments.hpp"
#include "cl_SDF_Parameters.hpp"
#include "cl_SDF_Mesh.hpp"
#include "cl_SDF_Object.hpp"
#include "cl_SDF_Core.hpp"
#include "cl_SDF_STK.hpp"
#include "cl_SDF_Field.hpp"
#include "cl_Logger.hpp" // MRS/IOS/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

using namespace moris;
using namespace sdf;

//------------------------------------------------------------------------------

//void
//save_exodus_mesh( const std::string & aMeshPath,
//------------------------------------------------------------------------------

void
unite_raycast(
        const uint                      & aNumberOfObjects,
        Vector< Matrix< DDRMat > > & aFieldValues,
        Vector< std::string >      & aFieldLabels )
{
    // get number of nodes
    uint tNumberOfNodes = aFieldValues( 0 ).length();

    Matrix< DDRMat > tEmpty;
    aFieldValues.push_back( tEmpty );
    Matrix< DDRMat > & tUnion = aFieldValues( aFieldLabels.size() );
    aFieldLabels.push_back( "Union_SDF" );

    tUnion.set_size( tNumberOfNodes, 1, 1.0 );

    // loop over all nodes
    for( uint k=0; k<tNumberOfNodes; ++k )
    {
        for( uint i=0; i<aNumberOfObjects; ++i )
        {
            // test sign
            if(  aFieldValues( i )( k ) < 0 )
            {
                tUnion( k ) = -1.0;
                break;
            }
        }
    }
}

void
unite_sdfs(
        const uint                      & aNumberOfObjects,
        Vector< Matrix< DDRMat > > & aFieldValues,
        Vector< std::string >      & aFieldLabels )
{
    // get number of nodes
    uint tNumberOfNodes = aFieldValues( 0 ).length();

    // get number of fields
    Matrix< DDRMat > tEmpty;

    aFieldValues.push_back( tEmpty );
    Matrix< DDRMat > & tUnion = aFieldValues( aFieldLabels.size() );
    aFieldLabels.push_back( "Union_SDF" );

    aFieldValues.push_back( tEmpty );
    Matrix< DDRMat > & tHasUnion = aFieldValues( aFieldLabels.size() );
    aFieldLabels.push_back( "Union_has_sdf" );

    tUnion.resize( tNumberOfNodes, 1 );
    tHasUnion.set_size( tNumberOfNodes, 1, 0.0 );

    real tEpsilon = 1e-9;

    // loop over all nodes
    for( uint k=0; k<tNumberOfNodes; ++k )
    {
        real tValue = 1e9;

        real tSign = 1.0;

        // loop over all objects
        for( uint i=0; i<aNumberOfObjects; ++i )
        {
            // test if SDF exists
            if ( aFieldValues( i+aNumberOfObjects )( k ) > 0 )
            {
                tValue = std::min( tValue, std::abs( aFieldValues( i )( k ) ) );
                tHasUnion( k ) = 1.0;
            }

            // test sign
            if( aFieldValues( i )( k ) < 0 && tSign > 0 )
            {
                tSign = -1.0;
            }
        }
        if( tHasUnion( k ) > 0.0 )
        {
            tUnion( k ) = tSign * tValue;
        }
        else
        {
            tUnion( k ) = tSign * tEpsilon;
        }
    }

    // broadcast min
    Matrix< DDRMat > tValues;
    comm_gather_and_broadcast( tUnion.min(), tValues );
    real tMinVal = tValues.min();

    // broadcasst max
    comm_gather_and_broadcast( tUnion.max(), tValues );
    real tMaxVal = tValues.max();

    // write dummy value into fields without sdf
    for( uint k=0; k<tNumberOfNodes; ++k )
    {
        if( ! tHasUnion( k ) )
        {
            if( tUnion( k ) < 0.0 )
            {
                tUnion( k ) = tMinVal;
            }
            else
            {
                tUnion( k ) = tMaxVal;
            }
        }
    }
}

void
perform_calculation(
        const Arguments   & aArguments,
        const bool aCalculateSDF )
{

    // step 1: load parameters
    Parameter_List           tGlobalParameters;
    Vector< Parameter_List > tObjectParameters;

    // load XML file
    load_sdf_parameter_list_from_xml(
            aArguments.get_parameter_path(),
            tGlobalParameters,
            tObjectParameters );

    // get verbose flag
    bool tVerbose = true;

    tic tTimer1;

    // step 2: create mesh objects
    mtk::Mesh * tMtkMesh = mtk::create_interpolation_mesh(
            mtk::MeshType::STK,
            aArguments.get_input_mesh_path(),
            nullptr,
            false );

    if( tVerbose )
    {
        real tElapsedTime = tTimer1.toc<moris::chronos::milliseconds>().wall;

            // print output
            std::fprintf( stdout,"%s Time for reading exodus file using MTK: %5.3f seconds.\n\n",
                    proc_string().c_str(),
                    ( double ) tElapsedTime / 1000 );
    }
    tic tTimer2;

    sdf::Mesh tMesh( tMtkMesh, tVerbose );
    if( tVerbose )
    {
        real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

        // print output
        std::fprintf( stdout,"%s Time for creating SDF Mesh: %5.3f seconds.\n\n",
                proc_string().c_str(),
                ( double ) tElapsedTime / 1000 );
    }

    // step 3: create output data

    uint tNumberOfObjects = tObjectParameters.size();

    // get number of nodes on mesh
    uint tNumberOfNodes    = tMesh.get_num_nodes();
    //uint tNumberOfElements = tMesh.get_num_elems();

    // create empty matrix
    Matrix< DDRMat > tEmpty( tNumberOfNodes, 1 );

    // create array for output data
    uint tNumberOfFields = 0;

    if( ! aCalculateSDF )
    {
        tNumberOfFields = tNumberOfObjects;
    }
    else
    {
        tNumberOfFields = 2*tNumberOfObjects;
    }

    Vector< Matrix< DDRMat > > tFieldValues( tNumberOfFields, tEmpty );
    Vector< std::string > tFieldLabels( tNumberOfFields, "" );

    // step 4: perform raycast

    // loop over all objects
    for( uint k=0; k<tNumberOfObjects; ++k )
    {
        // create SDF object
        sdf::Object tObject( tObjectParameters( k ).get< std::string >("stl_file") );

        // create core and set verbosity flag
        sdf::Core tCore( tMesh, tObject, tVerbose );

        // set parameters of core
        tCore.set_candidate_search_depth( tObjectParameters( k ).get< sint >("candidate_search_depth") );
        tCore.set_candidate_search_epsilon( tObjectParameters( k ).get< real >("candidate_search_epsilon") );

        // set label of output field
        tFieldLabels( k ) =  tObjectParameters( k ).get< std::string >("label");

        if( ! aCalculateSDF )
        {
            // perform raycast
            tCore.raycast_mesh();

            // create output matrix
            Matrix< DDRMat > tNodeIsInside( tNumberOfNodes, 1 );

            // loop over all nodes and write -1 if node is inside
            for( uint i=0; i<tNumberOfNodes; ++i )
            {
                if( tMesh.get_vertex( i )->is_inside() )
                {
                    tFieldValues( k )( i ) = -1.0;
                }
                else
                {
                    tFieldValues( k )( i ) = 1.0;
                }
            }
        }
        else
        {
            // perform raycast and sdf
            tCore.calculate_raycast_and_sdf( tFieldValues( k ) );

            // create flag if SDF exists
            uint j = tNumberOfObjects + k;
            for( uint i=0; i<tNumberOfNodes; ++i )
            {
                if( tMesh.get_vertex( i )->has_sdf() )
                {
                    tFieldValues( j )( i ) = 1.0;
                }
                else
                {
                    tFieldValues( j )( i ) = 0.0;
                }
            }

            // set label of sdf field
            tFieldLabels( j ) =  tFieldLabels( k ) + "_has_sdf";
        }
        // get path for output data and make it parallel
        std::string tBinaryPath = tObjectParameters( k ).get< std::string >( "output_values" );

        if( tBinaryPath.size() > 0 )
        {
            // save matrix to binary file
            save_matrix_to_binary_file( tFieldValues( k ), parallelize_path ( tBinaryPath ) );
        }

        // get HDF path
        std::string tHdf5Path = tObjectParameters( k ).get< std::string >( "output_hdf5" );

        if( tHdf5Path.size() > 0 )
        {
            // create field object
            Field tField( tFieldLabels( k ), tMesh, tFieldValues( k ) );

            // save field
            tField.save_field_to_hdf5( tHdf5Path );
        }
    }

    if( tNumberOfObjects > 1 )
    {
        if( aCalculateSDF )
        {
            unite_sdfs( tNumberOfObjects, tFieldValues, tFieldLabels );
        }
        else
        {
            unite_raycast( tNumberOfObjects, tFieldValues, tFieldLabels );
        }
    }

    // create output mesh
    std::string tOutputPath = aArguments.get_output_mesh_path();

    if( tOutputPath.size() > 0 )
    {
        sdf::STK tSTK( tMesh );
        tSTK.create_mesh_data( tFieldValues, tFieldLabels, aArguments.get_timestep() );

        tSTK.save_to_file( tOutputPath );
    }

    // tidy up memory
    delete tMtkMesh;
}

//------------------------------------------------------------------------------
int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );

//------------------------------------------------------------------------------

    // create arguments object
    Arguments tArguments( argc, argv );

//------------------------------------------------------------------------------

    switch ( tArguments.get_state() )
    {
        case( State::PRINT_USAGE ) :
        {
            // print system usage
            tArguments.print_usage();
            break;
        }
        case( State::PRINT_VERSION ) :
        {
            // print welcome banner and system information
            moris::print_banner( argc, argv );
            break;
        }
        case( State::PRINT_HELP ) :
        {
            // print help line and exit
            tArguments.print_help();
            break;
        }
        case( State::CALCULATE_RAYCAST ) :
        {
            // perform raycast
            perform_calculation( tArguments, false );
            break;
        }
        case( State::CALCULATE_RAYCAST_AND_SDF ) :
        {
            // perform raycast and calculate sdf
            perform_calculation( tArguments, true );
            break;
        }
        default:
        {
            // print system usage
            tArguments.print_usage();
            break;
        }
    }

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}


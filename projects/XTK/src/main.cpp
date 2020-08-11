/*
 * main.cpp
 *
 *  Created on: Mar 1, 2019
 *      Author: doble
 */

// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "HDF5_Tools.hpp"

// from linalg
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_print.hpp"

// XTK
#include "cl_XTK_Paramfile.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Output_Options.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "typedefs.hpp"
#include "cl_Logger.hpp" // MRS/IOS/src
#include "fn_compute_xtk_model_volumes.hpp"
#include "fn_compute_interface_surface_area.hpp"

#include "cl_GEN_Plane.hpp"
#include "cl_GEN_Sphere.hpp"
#include "cl_GEN_Superellipsoid.hpp"

// select namespaces
using namespace moris;
using namespace xtk;

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

moris_index
get_index_in_cell(
        Cell<std::string> & aLabels,
        std::string         aStr)
{
    auto  tIt = std::find(aLabels.begin(), aLabels.end(), aStr);
    MORIS_ERROR(tIt != aLabels.end(),"Radius not found in labels for sphere, please use r as the label");
    auto tPos = std::distance(aLabels.begin(), tIt);
    return tPos;
}

moris::ge::Geometry*
geometry_parse_factory(XTK_Problem_Params & aXTKProblemParams)
{
    moris::ge::Geometry* tGeometry = nullptr;
    if (aXTKProblemParams.mGeometryName == "sphere")
    {
        MORIS_ERROR( aXTKProblemParams.mRealGeomParams.size() == 4,"For a parsed sphere geometry there needs to be 4 parameters, r, xc, yc, zc");
        MORIS_ERROR( aXTKProblemParams.mRealGeomLabels.size() == 4,"For a parsed sphere geometry there needs to be 4 labels, r, xc, yc, zc");


        // find the radius in real parameters
        std::string tStr = "r";
        moris_index tPos = get_index_in_cell(aXTKProblemParams.mRealGeomLabels,tStr);
        moris::real tR =aXTKProblemParams.mRealGeomParams(tPos);

        // find the xc in real parameters
        tStr = "xc";
        tPos = get_index_in_cell(aXTKProblemParams.mRealGeomLabels,tStr);
        moris::real tXc =aXTKProblemParams.mRealGeomParams(tPos);

        // find the yc in real parameters
        tStr = "yc";
        tPos = get_index_in_cell(aXTKProblemParams.mRealGeomLabels,tStr);
        moris::real tYc = aXTKProblemParams.mRealGeomParams(tPos);

        // find the zc in real parameters
        tStr = "zc";
        tPos = get_index_in_cell(aXTKProblemParams.mRealGeomLabels,tStr);
        moris::real tZc =aXTKProblemParams.mRealGeomParams(tPos);

        tGeometry = new moris::ge::Sphere(tR,tXc,tYc,tZc);

    }
    else if (aXTKProblemParams.mGeometryName == "plane")
    {
        moris::Matrix<moris::DDRMat> tCenters(3,1);
        moris::Matrix<moris::DDRMat> tNormals(3,1);

        std::string tStr = "xc";
        moris_index tPos = get_index_in_cell(aXTKProblemParams.mRealGeomLabels,tStr);
        tCenters(0) =aXTKProblemParams.mRealGeomParams(tPos);

        // find the yc in real parameters
        tStr = "yc";
        tPos = get_index_in_cell(aXTKProblemParams.mRealGeomLabels,tStr);
        tCenters(1)= aXTKProblemParams.mRealGeomParams(tPos);

        // find the zc in real parameters
        tStr = "zc";
        tPos = get_index_in_cell(aXTKProblemParams.mRealGeomLabels,tStr);
        tCenters(2) =aXTKProblemParams.mRealGeomParams(tPos);

        // find the x normal
        tStr = "xnorm";
        tPos = get_index_in_cell(aXTKProblemParams.mRealGeomLabels,tStr);
        tNormals(0) = aXTKProblemParams.mRealGeomParams(tPos);

        tStr = "ynorm";
        tPos = get_index_in_cell(aXTKProblemParams.mRealGeomLabels,tStr);
        tNormals(1) = aXTKProblemParams.mRealGeomParams(tPos);

        tStr = "znorm";
        tPos = get_index_in_cell(aXTKProblemParams.mRealGeomLabels,tStr);
        tNormals(2) = aXTKProblemParams.mRealGeomParams(tPos);

        tGeometry = new moris::ge::Plane(tCenters(0), tCenters(1), tCenters(2), tNormals(0), tNormals(1), tNormals(2));
    }

    return tGeometry;

}

void run_xtk_problem(XTK_Problem_Params & aXTKProblemParams)
{
    // times to keep track of
    real tFullTime       = 0.0;
    real tMeshLoadTime   = 0.0;
    real tGeometryTime   = 0.0;
    real tDecompTime     = 0.0;
    real tSensTime       = 0.0;
    real tUnzipTime      = 0.0;
    real tGhostTime      = 0.0;
    real tEnrichmentTime = 0.0;
    real tOutputTime     = 0.0;
    real tWriteObjTime   = 0.0;
    real tWriteData      = 0.0;

    // print some stuff
    if(par_rank() == 0)
    {
        std::cout<<"-------------------------------------------------------------"<<std::endl;
        std::cout<<"| XTK PROBLEM INFO                                           |"<<std::endl;
        std::cout<<"-------------------------------------------------------------"<<std::endl;
        std::cout<<"Input Mesh File:  "<<aXTKProblemParams.mInputMeshFile<<std::endl;
        aXTKProblemParams.print_geom();
        aXTKProblemParams.print_operations();

        if(aXTKProblemParams.mExport)
        {
            std::cout<<"Output Mesh File: "<<aXTKProblemParams.mOutputMeshFile<< std::endl;
        }
        if(aXTKProblemParams.mOutputData)
        {
            std::cout<<"Output Data File: "<<aXTKProblemParams.mDataFile<< std::endl;
        }
        if(aXTKProblemParams.mWriteobj)
        {
            std::cout<<"Obj Surf File:    "<<aXTKProblemParams.mobjOutputFile<< std::endl;
        }

    }


    std::clock_t  tFullTimer = std::clock();


    // initialize mesh
    std::clock_t tOpTimer = std::clock();
    moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( aXTKProblemParams.mMeshType, aXTKProblemParams.mInputMeshFile );
    tMeshLoadTime = (std::clock() - tOpTimer)/(CLOCKS_PER_SEC/1000);

    // setup the geometry
    //TODO: support multiple geometries
    tOpTimer = std::clock();
    Cell<std::shared_ptr<moris::ge::Geometry>> tGeometry;
    tGeometry.resize(1);
    tGeometry(0) = std::shared_ptr<moris::ge::Geometry>(geometry_parse_factory(aXTKProblemParams));
    tGeometryTime = (std::clock() - tOpTimer)/(CLOCKS_PER_SEC/1000);

    // setup the geometry engine
    //TODO: support multiple geometries, and different phase tables
    //          Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    moris::ge::Phase_Table tPhaseTable (1,  "exp_base_2");
    moris::ge::Geometry_Engine tGeometryEngine(tGeometry, tPhaseTable, tMeshData);

    // setup the XTK model
    Model tXTKModel(3,tMeshData,&tGeometryEngine);
    tXTKModel.mVerbose  =  false;

    // decompose the mesh
    tOpTimer = std::clock();
    tXTKModel.decompose(aXTKProblemParams.mSubdivisionMethods);
    tDecompTime = (std::clock() - tOpTimer)/(CLOCKS_PER_SEC/1000);

    // compute sensitivity
    if(aXTKProblemParams.mComputeSens)
    {
        tOpTimer = std::clock();
        tXTKModel.communicate_interface_nodes();
        tSensTime = (std::clock() - tOpTimer)/(CLOCKS_PER_SEC/1000);
    }


    // compute sensitivity
    if(aXTKProblemParams.mUnzip)
    {
        tOpTimer = std::clock();
        tXTKModel.unzip_interface();
        tUnzipTime = (std::clock() - tOpTimer)/(CLOCKS_PER_SEC/1000);
    }

    // perform ghost stabilization
    if(aXTKProblemParams.mGhost)
    {
        tOpTimer = std::clock();
        tXTKModel.construct_face_oriented_ghost_penalization_cells();
        tGhostTime = (std::clock() - tOpTimer)/(CLOCKS_PER_SEC/1000);
    }

    // perform ghost stabilization
    if(aXTKProblemParams.mEnrich)
    {
        tOpTimer = std::clock();
        tXTKModel.perform_basis_enrichment(EntityRank::NODE);
        tEnrichmentTime = (std::clock() - tOpTimer)/(CLOCKS_PER_SEC/1000);
    }

    // export
    // TODO: support output option parsing
    if(aXTKProblemParams.mExport)
    {
        Output_Options tXTKOutputOptions;
        tXTKOutputOptions.mPackageDxDpSparsely = false;

        tOpTimer = std::clock();
        moris::mtk::Mesh* tOutputMTK = tXTKModel.get_output_mesh(tXTKOutputOptions);
        tOutputMTK->create_output_mesh(aXTKProblemParams.mOutputMeshFile);
        tOutputTime = (std::clock() - tOpTimer)/(CLOCKS_PER_SEC/1000);

        delete tOutputMTK;
    }

    if(aXTKProblemParams.mWriteobj)
    {
        tOpTimer = std::clock();
        moris::Cell<std::string> tBoundingSets = {"surface_1","surface_2","surface_3","surface_4","surface_5","surface_6"};


        tXTKModel.extract_surface_mesh_to_obj(aXTKProblemParams.mobjOutputFile,
                aXTKProblemParams.mPhaseForobj,
                tBoundingSets);
        tWriteObjTime = (std::clock() - tOpTimer)/(CLOCKS_PER_SEC/1000);
    }


    // stop full clock timer
    tFullTime = (std::clock() - tFullTimer)/(CLOCKS_PER_SEC/1000);

    if(aXTKProblemParams.mOutputData)
    {
        tOpTimer = std::clock();
        //TODO: abstract to support more phases
        moris::Matrix<moris::DDRMat> tNodeCoords = tXTKModel.get_background_mesh().get_all_node_coordinates_loc_inds();
        moris::real tParentPhase0Vol = compute_non_intersected_parent_element_volume_by_phase(0,tNodeCoords,tXTKModel);
        moris::real tParentPhase1Vol = compute_non_intersected_parent_element_volume_by_phase(1,tNodeCoords,tXTKModel);
        moris::real tChildPhase0Vol  = compute_child_element_volume_by_phase(0,tNodeCoords,tXTKModel);
        moris::real tChildPhase1Vol  = compute_child_element_volume_by_phase(1,tNodeCoords,tXTKModel);
        moris::real tMySurfaceArea   = compute_interface_surface_area(tNodeCoords,tXTKModel);

        moris::real tMyPhase0Volume = tParentPhase0Vol + tChildPhase0Vol;
        moris::real tMyPhase1Volume = tParentPhase1Vol + tChildPhase1Vol;

        // sum up volumes and surface areas globally
        // Collect all volumes
        moris::real tGlbVolume0 = sum_all(tMyPhase0Volume);

        // Collect all volumes
        moris::real tGlbVolume1 = sum_all(tMyPhase1Volume);

        moris::real tTotalVolume = tGlbVolume1 + tGlbVolume0;

        // Collect all surface areas
        moris::real tGlbSurf = sum_all(tMySurfaceArea);

        hid_t tHDFId = create_hdf5_file(aXTKProblemParams.mDataFile.c_str());
        herr_t tErr;
        save_scalar_to_hdf5_file(tHDFId,"full_time",tFullTime,tErr);
        save_scalar_to_hdf5_file(tHDFId,"mesh_in_time",tMeshLoadTime,tErr);
        save_scalar_to_hdf5_file(tHDFId,"geeom_time",tGeometryTime,tErr);
        save_scalar_to_hdf5_file(tHDFId,"decomp_time",tDecompTime,tErr);
        save_scalar_to_hdf5_file(tHDFId,"sens_time",tSensTime,tErr);
        save_scalar_to_hdf5_file(tHDFId,"unzip_time",tUnzipTime,tErr);
        save_scalar_to_hdf5_file(tHDFId,"ghost_time",tGhostTime,tErr);
        save_scalar_to_hdf5_file(tHDFId,"enrich_time",tEnrichmentTime,tErr);
        save_scalar_to_hdf5_file(tHDFId,"mesh_out_time",tOutputTime,tErr);
        save_scalar_to_hdf5_file(tHDFId,"surf_area",tGlbSurf,tErr);
        save_scalar_to_hdf5_file(tHDFId,"volume",tTotalVolume,tErr);
        save_scalar_to_hdf5_file(tHDFId,"par_size",par_size(),tErr);

        tWriteData = (std::clock() - tOpTimer)/(CLOCKS_PER_SEC/1000);
    }

    // print timing
    if(par_rank() == 0)
    {
        std::cout<<"-------------------------------------------------------------"<<std::endl;
        std::cout<<"| XTK TIMING INFO (note: Not all times displayed below)     |"<<std::endl;
        std::cout<<"-------------------------------------------------------------"<<std::endl;
        std::cout<<"Full XTK problem:    "<<tFullTime<<" ms."<<std::endl;
        std::cout<<"Mesh Loading:        "<<tMeshLoadTime<<" ms."<<std::endl;
        std::cout<<"Geometry Setup:      "<<tGeometryTime<<" ms."<<std::endl;
        std::cout<<"Decomposition:       "<<tDecompTime<<" ms."<<std::endl;
        std::cout<<"Sensitivity Compute: "<<tSensTime<<" ms."<<std::endl;
        std::cout<<"Interface Unzip:     "<<tUnzipTime<<" ms."<<std::endl;
        std::cout<<"Ghost Penalty:       "<<tGhostTime<<" ms."<<std::endl;
        std::cout<<"Basis Enrichment:    "<<tEnrichmentTime<<" ms."<<std::endl;
        std::cout<<"Mesh Output Time:    "<<tOutputTime<<" ms."<<std::endl;
        std::cout<<"Write Obj Time:      "<<tWriteObjTime<<" ms."<<std::endl;
        std::cout<<"Data Output Time:    "<<tWriteData<<" ms."<<std::endl;
    }
    delete tMeshData;
}

int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );

    if(par_rank() == 0)
    {
        std::cout<<"XML file = "<<argv[1]<<std::endl;
    }

    // initialize the xtk parameters
    Paramfile tParams(argv[1]);

    // retrieve the XTK problems
    Cell<XTK_Problem_Params> & tXTKProblems = tParams.get_xtk_problem_params();

    // run the problem
    run_xtk_problem(tXTKProblems(0));


    gMorisComm.finalize();

    return 0;

}

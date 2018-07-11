/*#include <catch.hpp>

#include "fn_hash_value.hpp" // ALG/src
#include "typedefs.hpp" // COR/src
#include "cl_Mat.hpp" // LNA/src
#include "GeUtilities.hpp" // GEN/src
#include "cl_GeTriangle.hpp" // GEN/src
#include "cl_GeSDFGen.hpp" // GEN/src

TEST_CASE("ge::SDFGen","[geomeng],[SDFGen]")
{
    moris::uint tModelDim = 3;    // ModelDimension
    moris::uint tPolynomial = 1;  // Polynomial degree of the basis functions

    moris::Mat<moris::uint> tNumElements ={{20,20,20}};

    // the path must be set correctly
    std::string tMeshFilePath = "./test/src/geomeng/objfiles/designbox.obj";

    // origin of hierarchical mesh
    moris::Mat<moris::real> tOrigin =  {{-4.88,3.37,-1.594}};

    // dimensions of hierarchical mesh
    moris::Mat<moris::real> tDimensions = {{6.36,2.76,3.588}};

    moris::Hierarchical_Mesh_Main tHMR(tModelDim,tPolynomial,tNumElements,tDimensions);

    tHMR.set_point_of_origin(tOrigin);

    // create SDF object
    ge::SDFGen tSDFGen(tHMR, tMeshFilePath);

    // some optional settings for the SDF Generator

    // candidate search depth (default is 0)
    //tSDFGen.set_candidate_search_depth(2);

    // diagonal epsilon (default is 0.01)
    //tSDFGen.set_buffer_epsilon(0.01);

    // activate heavytest (default is false)
    //tSDFGen.set_heavytest_switch(true);

    // set heavytest error tolerance (default is zero)
    //tSDFGen.set_heavytest_allowed_error(0.00001);

    // activate sweeping (default is false)
    //tSDFGen.set_sweeping_switch(true);

    // create mesh
    std::cout << "Creating Reference Mesh ...";
    tHMR.create_mesh();
    std::cout << " done." << std::endl;


    // calculate signed distance field
    std::cout << "Starting SDF ...";
    tSDFGen.calculate_sdf();
    std::cout << " done." << std::endl;


    // get SDF
    moris::Mat<moris::real> tSDF = tSDFGen.get_sdf();

    // get active nodes as Moris<uint>
    moris::Mat<moris::uint> tActiveNodes = tSDFGen.get_active_nodes();

    // get active nodes as BoostBitset (if needed)
    //moris::BoostBitset  tActiveFlags = tSDFGen.get_active_flags();

#ifdef MORIS_USE_ARMA
    REQUIRE(hash_value(tSDF)         ==  15334645698042653713UL);
    REQUIRE(hash_value(tActiveNodes) ==  6189272998541787136UL);
#endif
}*/

//
// Created by messe on 1/25/18.
//

/*
#include <catch.hpp>

#include "cl_GeTriangle.hpp" // GEN/src


// solutions generated with Kurt's MATLAB code
TEST_CASE(
        "ge::Triangle",
        "[geomeng],[Triangle]")
{
    moris::Mat<moris::real> tNodeCoords = {{1.050229216800883, 0.649117827347835, 1.088464523227452},
                                         {1.417028287272334,   1.192051358390624, -0.179967737969674},
                                         {1.334891429874816, 0.076610713185436, 0.309230566913958}};

    // create triangle
    ge::Triangle tTriangle(0,1,2);

    // assign coordinates
    tTriangle.set_node_coords(tNodeCoords);

    moris::Mat<moris::real> tPoint= {{0.885670154103523},
                                   {0.730289641563979},
                                   {1.121888622244782}};

    SECTION("ge::Triangle : center"){
        // the center of the triangle
        moris::Mat<moris::real> tCenter = {{0.929270522458724},
                                         {0.809703969231095},
                                         {0.573577569991403}};

        moris::Mat<moris::real> tError = tCenter - tTriangle.get_center();
        REQUIRE(ge::norm(tError) < MORIS_GE_EPSILON);
    }

    SECTION("ge::Triangle : normal"){
        // the normal of the triangle
        moris::Mat<moris::real> tNormal = {{-0.912893263456338},
                                         {-0.235837187416504},
                                         { 0.333176695714916}};
        moris::Mat<moris::real> tError = tNormal - tTriangle.get_normal();
        REQUIRE(ge::norm(tError) < MORIS_GE_EPSILON);
    }

    SECTION("ge::Triangle : hesse"){

        // the hesse distance
        moris::real tHesse = -0.848180427118634;
        moris::real tError = tHesse - tTriangle.get_hesse();
        REQUIRE(ge::abs(tError) < MORIS_GE_EPSILON);
    }

    SECTION("ge::Triangle : bounding box"){

        // minimum coordinate
        moris::Mat<moris::real> tMinCoordinate = {{ 0.649117827347835},
                                                {-0.179967737969674},
                                                { 0.076610713185436}};

        // maximum coordinate
        moris::Mat<moris::real> tMaxCoordinate = {{1.088464523227452},
                                                {1.417028287272334},
                                                {1.334891429874816}};

        moris::Mat<moris::real> tError(3,1);
        for (moris::uint i=0; i<3; ++i){
            tError(i) = tMinCoordinate(i) - tTriangle.get_min_coord(i);
        }
        REQUIRE(ge::norm(tError) < MORIS_GE_EPSILON);

        for (moris::uint i=0; i<3; ++i){
            tError(i) = tMaxCoordinate(i) - tTriangle.get_max_coord(i);
        }
        REQUIRE(ge::norm(tError) < MORIS_GE_EPSILON);
    }

    SECTION("ge::Triangle : check edges"){
        // projection in y-z plane
        REQUIRE(tTriangle.check_edge(0,0,tPoint) == true);
        REQUIRE(tTriangle.check_edge(1,0,tPoint) == false);
        REQUIRE(tTriangle.check_edge(2,0,tPoint) == true);

        // projection in x-z plane
        REQUIRE(tTriangle.check_edge(0,1,tPoint) == true);
        REQUIRE(tTriangle.check_edge(1,1,tPoint) == true);
        REQUIRE(tTriangle.check_edge(2,1,tPoint) == false);

        // projection in x-y plane
        REQUIRE(tTriangle.check_edge(0,2,tPoint) == true);
        REQUIRE(tTriangle.check_edge(1,2,tPoint) == true);
        REQUIRE(tTriangle.check_edge(2,2,tPoint) == true);

    }

    SECTION("ge::Triangle : intersection tests"){
        moris::Mat<moris::real> tDirection = {   {0.677108707530041},
                                               {0.677108707530041},
                                               {-0.506439712788708}};


        moris::Mat<moris::real> tIntersection = {{1.0582225162727209846334585509245},
                                               {0.90284200373317699340953594768469},
                                               {0.99282903732709029970246501366325}};

        moris::Mat<moris::real> tError = tIntersection - tTriangle.intersect_with_line(tPoint, tDirection);
        REQUIRE(ge::norm(tError) < MORIS_GE_EPSILON);

        moris::real tIntersectionYZ = 1.1499023579142249836729065160288;

        REQUIRE(ge::abs(tIntersectionYZ-tTriangle.intersect_with_coordinate_axis(tPoint, 0))
                < MORIS_GE_EPSILON);

        moris::real tIntersectionXZ = 1.7530961017725150876322698604616;
        REQUIRE(ge::abs(tIntersectionXZ-tTriangle.intersect_with_coordinate_axis(tPoint, 1))
                < MORIS_GE_EPSILON);

        moris::real tIntersectionXY = 0.39790101461988117245214654257444;
        REQUIRE(ge::abs(tIntersectionXY-tTriangle.intersect_with_coordinate_axis(tPoint, 2))
                < MORIS_GE_EPSILON);
    }

    SECTION("ge::Triangle : area") {
        moris::real tArea = 0.974220833569323;

        REQUIRE(ge::abs(tArea-tTriangle.get_area()) < 10*MORIS_GE_EPSILON);
    }

    SECTION("ge::Triangle : Distance Tests (stage 0)") {

        moris::Mat<moris::real> tProjection = {{-0.488600443461426},
                                             { 0.109257414478981},
                                             { 0.241215798847012}};

        moris::Mat<moris::real> tBarycentric = {{ 0.691336629346961},
                                              {-0.099792207295336},
                                              { 0.408455577948375}};

        moris::Mat<moris::real> tError =  tProjection-tTriangle.project_point_to_local_cartesian(tPoint);
        REQUIRE(ge::norm(tError) < MORIS_GE_EPSILON);


        tError = tBarycentric-tTriangle.get_barycentric_from_local_cartesian(tProjection);
        REQUIRE(ge::norm(tError) < MORIS_GE_EPSILON);

        moris::Mat<moris::real> tDistances = {{0.954058443607797},
                                            {0.262060528044059},
                                            {0.641160884691928}};

        for(moris::uint k=0; k<3; ++k){
            REQUIRE(ge::abs(tDistances(k)-tTriangle.distance_point_to_edge_in_local_cartesian(tProjection, k))
              < MORIS_GE_EPSILON);
        }
    }

    SECTION("ge::Triangle : Distance Tests (stage 1)") {
        moris::Mat<moris::real> tPoints = {{ 0.059234, 0.975223, 0.052644, 0.705246, 0.670342,
                                           0.455441, 1.663529, 0.338271, 1.382485, 1.403031,
                                           1.646998, 1.779594, 0.512393, 0.691528, 0.678431,
                                           1.164209, 0.823881, 1.839212, 0.957501, 1.783941,
                                           1.024442, 0.208111, 1.288813, 1.623248,-0.015943},
                                         {-0.113611, 0.088797, 0.975141, 1.195421, 1.682441,
                                           1.012004, 0.694506, 0.804913, 1.378310, 0.413346,
                                           0.208015, 0.130009, 1.783196, 1.433054, 0.773284,
                                           1.513817, 0.104125, 1.447990, 0.393489, 0.810006,
                                           1.232155, 1.255812, 0.866007, 0.752600, 0.015000},
                                         { 0.130866,-0.064282,-0.295901, 1.515404, 1.378770,
                                           0.737314, 1.395056,-0.325510,-0.079550, 0.211413,
                                           1.526104, 0.966354,-0.035633, 0.241918,-0.407221,
                                           0.381194, 0.951702, 1.501818, 1.012948, 0.834930,
                                           0.943918, 0.794787, -0.352999, 0.022700, -0.259247}};

        moris::Mat<moris::real> tDistances = {{0.966945530788054}, {0.323467402524833}, {0.735932134607854},
                                            {0.440445580657558}, {0.464652493703906}, {0.439399071526912},
                                            {0.745953762707671}, {0.604685436852981}, {0.765441776142965},
                                            {0.459681594500460}, {0.998665297040420}, {0.802059606964542},
                                            {0.617044825654042}, {0.209197658398687}, {0.551743281547435},
                                            {0.466112521766804}, {0.462261824170279}, {0.807042092863447},
                                            {0.304215233076093}, {0.723787800878556}, {0.063122262878454},
                                            {0.638982157601383}, {0.717369656022604}, {0.803599693257272},
                                            {1.101429307019322}};

        moris::Mat<moris::real> tPointK(3,1);
        moris::Mat<moris::real> tError(25, 1);
        for(moris::uint k=0; k<25; ++k){
            for(moris::uint i=0; i<3; ++i){
                tPointK(i) = tPoints(i,k);
            }
            tError(k) = tDistances(k) - tTriangle.get_distance_to_point(tPointK);
        }

        REQUIRE(ge::norm(tError) < 10*MORIS_GE_EPSILON);

    }

}

*/
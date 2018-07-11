/*
 * SDF_Core.cpp
 *
 *  Created on: Feb 28, 2018
 *      Author: messe
 */

// sadly, direct access not work with eigen: http://eigen.tuxfamily.org/bz/show_bug.cgi?id=329

#include "cl_Ge_SDF_Core.hpp"

// -----------------------------------------------------------------------------

void
ge::SDF_Core::calculate_raycast(
        const ge::SDF_Mesh_Data & aMeshData,
        ge::SDF_Data            & aData)
{
    moris::tic tTimer;

    aData.init_data_fields (aMeshData.get_number_of_nodes(),
                            aMeshData.get_number_of_elements());

    // perform voxelizing algorithm
    // z-direction:
    // intialization of results and correction list
    voxelize(aMeshData, aData, 2);
    if(aData.mUnsureNewNodesCount > 0)
    {
        // y-direction
        // has only to be done for those points,
        // that were not determinable in the first step
        voxelize(aMeshData, aData, 1);
        if(aData.mUnsureNewNodesCount > 0)
        {
            // finally, do it for x-direction
            // x-direction
            voxelize(aMeshData, aData, 0);
        }
    }
    // write error message if not all nodes were found
    if( aData.mUnsureNewNodesCount > 0 )
    {
        MORIS_LOG_ERROR << "Raycast could not detect all nodes\n";
    }

    // needed for elements on surface and in volume
    calculate_candidate_points_and_buffer_diagonal(aMeshData, aData);

    // stop the timer
    moris::real tElapsedTime
            = tTimer.toc<moris::chronos::milliseconds>().wall;

    // print elapsed time
    if(moris::par_size() == 1)
    {
        std::fprintf(stdout, "Time for ray cast              : %5.3f [sec]\n",
                     tElapsedTime/1000);
    }
    else
    {
        std::fprintf(stdout, "Proc % i - Time for ray cast              : %5.3f [sec]\n",
                     (int) moris::par_rank(), tElapsedTime/1000);
    }

}
// -----------------------------------------------------------------------------

void
ge::SDF_Core::calculate_raycast_and_sdf(
    const ge::SDF_Mesh_Data & aMeshData,
    ge::SDF_Data            & aData)
{
    calculate_raycast(aMeshData, aData);

    // start the tomer
    moris::tic tTimer;

    calculate_udf(aMeshData, aData);
    sign_udf(aMeshData, aData);

    // stop the timer
    moris::real tElapsedTime
            = tTimer.toc<moris::chronos::milliseconds>().wall;

    // print elapsed time
    if(moris::par_size() == 1)
    {
        std::fprintf(stdout, "Time for signed distance field : %5.3f [sec]\n",
                     tElapsedTime/1000);
    }
    else
    {
        std::fprintf(stdout, "Proc % i - Time for signed distance field : %5.3f [sec]\n",
                     (int) moris::par_rank(), tElapsedTime/1000);
    }

}

// -----------------------------------------------------------------------------

void
ge::SDF_Core::save_to_vtk(
        const ge::SDF_Mesh_Data & aMeshData,
        const ge::SDF_Data      & aData,
        const std::string       & aFilePath)
{

    // modify filename
    std::string tFilePath;
    if ( moris::par_size() > 1 )
    {
        auto tFileExt = aFilePath.substr(aFilePath.find_last_of("."),aFilePath.length());
        auto tBasePath = aFilePath.substr(0,aFilePath.find_last_of("."));
        tFilePath = tBasePath + "_" +  std::to_string(moris::par_rank()) + tFileExt;

    }
    else
    {
        tFilePath = aFilePath;
    }

    // open the file
    std::ofstream tFile(tFilePath, std::ios::binary);

    float tFValue = 0;
    int   tIValue = 0;
    float tFChar = 0;
    int   tIChar = 0;

    // write header
    tFile << "# vtk DataFile Version 3.0" << std::endl;
    tFile << "GO BUFFS!" << std::endl;
    tFile << "BINARY" << std::endl;
    // write node data
    tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
    moris::uint tNumberOfNodes = aMeshData.get_number_of_nodes();
    tFile << "POINTS " << tNumberOfNodes << " float"  << std::endl;
    for ( moris::uint k = 0; k < tNumberOfNodes; ++k )
    {
        moris::Mat< moris::real > tNodeCoords = aMeshData.get_node_coords(k);
        for ( moris::uint i = 0; i < 3; ++i )
        {
            tFValue = (float)  tNodeCoords(i);
            tFChar = ge::swap_byte_endian(tFValue);
            tFile.write( (char*) &tFChar, sizeof(float));
        }
    }
    tFile << std::endl;
    auto tNumberOfElements = aMeshData.get_number_of_elements();

    //moris::uint tN = aMeshData.mElementTopology.n_rows();
    moris::uint tN = (aMeshData.get_nodes_of_element(0)).length();

    moris::uint tCellType = 0;
    // map element topology
    moris::Mat<moris::uint> tIndex( tN, 1 );
    if( tN == 8 )
    {
        tIndex(0) = 0;
        tIndex(1) = 1;
        tIndex(2) = 3;
        tIndex(3) = 2;
        tIndex(4) = 4;
        tIndex(5) = 5;
        tIndex(6) = 7;
        tIndex(7) = 6;
        tCellType = 12;
    } else if( tN == 27 )
    {
        tIndex( 0) =  0;
        tIndex( 1) =  2;
        tIndex( 2) =  8;
        tIndex( 3) =  6;
        tIndex( 4) = 18;
        tIndex( 5) = 20;
        tIndex( 6) = 26;
        tIndex( 7) = 24;
        tIndex( 8) =  1;
        tIndex( 9) =  5;
        tIndex(10) =  7;
        tIndex(11) =  3;
        tIndex(12) = 19;
        tIndex(13) = 23;
        tIndex(14) = 25;
        tIndex(15) = 21;
        tIndex(16) =  9;
        tIndex(17) = 11;
        tIndex(18) = 17;
        tIndex(19) = 15;
        tIndex(20) = 12;
        tIndex(21) = 14;
        tIndex(22) = 10;
        tIndex(23) = 16;
        tIndex(24) =  4;
        tIndex(25) = 22;
        tIndex(26) = 13;
        tCellType = 29;
    }
    else
    {
        throw std::invalid_argument( "Unsupported Polynomial");
    }

    // write element topology
    tFile << "CELLS " << tNumberOfElements << " "
          << (tN+1)*tNumberOfElements  << std::endl;
    for (moris::uint k = 0; k < tNumberOfElements; ++k)
    {

        tIValue = (int)  tN;
        tIChar = ge::swap_byte_endian(tIValue);
        tFile.write( (char*) &tIChar, sizeof(int));

        moris::Mat< moris::uint > tNodesOfElement = aMeshData.get_nodes_of_element(k);
        MORIS_ASSERT(tNodesOfElement.length() == tN, "For VTK output, all elements must be of the same type.");

        for (moris::uint i=0; i< tN; ++i)
        {
            tIValue = (int) tNodesOfElement(tIndex(i));
            tIChar = ge::swap_byte_endian(tIValue);
            tFile.write((char *) &tIChar, sizeof(int));
        }

    }
    tFile << std::endl;
    tFile << "CELL_TYPES " << tNumberOfElements << std::endl;
    tIValue = (int) tCellType;
    tIChar = ge::swap_byte_endian(tIValue);
    for (moris::uint k = 0; k < tNumberOfElements; ++k)
    {
        tFile.write( (char*) &tIChar, sizeof(int));
    }
    // write element data
    tFile << "CELL_DATA " << tNumberOfElements << std::endl;
    tFile << "SCALARS ElementId int" << std::endl;
    tFile << "LOOKUP_TABLE default" << std::endl;

    for (moris::uint k = 0; k < tNumberOfElements; ++k)
    {
        tIValue = (int) aMeshData.get_global_element_id(k);
        tIChar = ge::swap_byte_endian(tIValue);
        tFile.write( (char*) &tIChar, sizeof(float));
    }

    moris::Mat< moris::uint > tElementFlags(tNumberOfElements, 1, 0 );
    for (moris::uint k = 0; k < aData.mLocalElementsInVolume.length(); ++k)
    {
        tElementFlags( aData.mLocalElementsInVolume ( k ) ) = 1;
    }

    tFile << "SCALARS ElementsInVolume int" << std::endl;
    tFile << "LOOKUP_TABLE default" << std::endl;
    for (moris::uint k = 0; k < tNumberOfElements; ++k)
    {
        tIValue = (int) tElementFlags( k );
        tIChar = ge::swap_byte_endian(tIValue);
        tFile.write( (char*) &tIChar, sizeof(float));
    }

    tElementFlags.fill(0);

    for (moris::uint k = 0; k < aData.mLocalElementsAtSurface.length(); ++k)
    {
        tElementFlags( aData.mLocalElementsAtSurface ( k ) ) = 1;
    }

    tFile << "SCALARS ElementsAtSurface int" << std::endl;
    tFile << "LOOKUP_TABLE default" << std::endl;
    for (moris::uint k = 0; k < tNumberOfElements; ++k)
    {
        tIValue = (int) tElementFlags( k );
        tIChar = ge::swap_byte_endian(tIValue);
        tFile.write( (char*) &tIChar, sizeof(float));
    }
    // write node data
    tFile << "POINT_DATA " << tNumberOfNodes << std::endl;
    tFile << "SCALARS NodeID int" << std::endl;
    tFile << "LOOKUP_TABLE default" << std::endl;
    for (moris::uint k = 0; k < tNumberOfNodes; ++k)
    {
        tIValue = (int) aMeshData.get_global_node_id(k);
        tIChar = ge::swap_byte_endian(tIValue);
        tFile.write( (char*) &tIChar, sizeof(float));
    }

    tFile << "SCALARS LocalNodeInd int" << std::endl;
    tFile << "LOOKUP_TABLE default" << std::endl;
    for (moris::uint k = 0; k < tNumberOfNodes; ++k)
    {
        tIValue = (int) k;
        tIChar = ge::swap_byte_endian(tIValue);
        tFile.write( (char*) &tIChar, sizeof(float));
    }

    tFile << "SCALARS NodeIsInside int" << std::endl;
    tFile << "LOOKUP_TABLE default" << std::endl;
    for (moris::uint k = 0; k < tNumberOfNodes; ++k)
    {
        if(aData.mLocalNodeInsideFlags.test(k))
        {
            tIValue = (int) 1;
        }
        else
        {
            tIValue = (int) 0;
        }

        tIChar = ge::swap_byte_endian(tIValue);
        tFile.write( (char*) &tIChar, sizeof(float));
    }

    tFile << "SCALARS NodeIsCandidate int" << std::endl;
    tFile << "LOOKUP_TABLE default" << std::endl;
    for (moris::uint k = 0; k < tNumberOfNodes; ++k)
    {
        if(aData.mLocalNodeCandidateFlags.test(k))
        {
            tIValue = (int) 1;
        }
        else
        {
            tIValue = (int) 0;
        }

        tIChar = ge::swap_byte_endian(tIValue);
        tFile.write( (char*) &tIChar, sizeof(float));
    }
    tFile << "SCALARS NodeHasSDF int" << std::endl;
    tFile << "LOOKUP_TABLE default" << std::endl;
    for (moris::uint k = 0; k < tNumberOfNodes; ++k)
    {
        if(aData.mLocalNodeHasSdfFlag.test(k))
        {
            tIValue = (int) 1;
        }
        else
        {
            tIValue = (int) 0;
        }

        tIChar = ge::swap_byte_endian(tIValue);
        tFile.write( (char*) &tIChar, sizeof(float));
    }

    // write node data
    tFile << "SCALARS SDF float" << std::endl;
    tFile << "LOOKUP_TABLE default" << std::endl;
    for (moris::uint k = 0; k < tNumberOfNodes; ++k)
    {
        tFValue = (float) aData.mLocalSDF(k);
        tFChar = ge::swap_byte_endian(tFValue);
        tFile.write( (char*) &tFChar, sizeof(float));
    }

    // close the output file
    tFile.close();
}

// =============================================================================
//  Private
// =============================================================================

void
ge::SDF_Core::voxelize(
        const ge::SDF_Mesh_Data& aMeshData,
        ge::SDF_Data           & aData,
        const moris::uint        aAxis)
{
    // node coordinate
    moris::Mat<moris::real> tNodeCoords(3,1);

    // reset unsure nodes counter
    aData.mUnsureNewNodesCount = 0;

    // loop over all nodes
    for ( moris::uint k=0; k<aData.mUnsureNodes.length(); ++k )
    {
        // copy node data
        tNodeCoords = aMeshData.get_node_coords( aData.mUnsureNodes( k ));

        // preselect triangles for intersection test
        if(aAxis == 0)
            preselect_triangles_x(aData, tNodeCoords);
        else if (aAxis == 1)
            preselect_triangles_y(aData, tNodeCoords);
        else
            preselect_triangles_z(aData, tNodeCoords);

        // from the candidate triangles, perform intersection
        if(aData.mCandidateTriangles.length() > 0)
        {
            intersect_triangles( aData, aAxis, tNodeCoords );

            // intersect ray with triangles and check if node is inside
            if(aData.mIntersectedTriangles.length() > 0)
            {
                intersect_ray_with_triangles( aData, aAxis, tNodeCoords );
                check_if_node_is_inside(aMeshData,
                                        aData,
                                        aAxis,
                                        aData.mUnsureNodes( k ),
                                        tNodeCoords );
            }

        }

    }

    // copy back unsure nodes and adapt vector length
    for (moris::uint k=0; k<aData.mUnsureNewNodesCount; ++k)
    {
        aData.mUnsureNodes(k) = aData.mUnsureNodesNew(k);
    }
    aData.mUnsureNodes.resize(aData.mUnsureNewNodesCount, 1);

}

// -----------------------------------------------------------------------------

/* @TODO check if we really need the if statements */
void
ge::SDF_Core::preselect_triangles_x(
        ge::SDF_Data                  & aData,
        const moris::Mat<moris::real> & aNodeCoords )
{
    // x: k = x, j = z, i = y
#ifdef MORIS_USE_ARMA

    // check bounding box in J-direction
    aData.mCandJ = arma::find(
            (aNodeCoords(2) - aData.mTriangleMinCoordsZ) %
            (aData.mTriangleMaxCoordsZ - aNodeCoords(2)) > -MORIS_GE_EPSILON );

    // check bounding box in I-direction
    aData.mCandI = arma::find(
            (aNodeCoords(1) - aData.mTriangleMinCoordsY.elem(aData.mCandJ)) %
            (aData.mTriangleMaxCoordsY.elem(aData.mCandJ) - aNodeCoords(1)) > -MORIS_GE_EPSILON );

    // help vector to be written in aData.mCandidateTriangles.data()
    aData.mCandK = aData.mCandJ.elem(aData.mCandI);
    // resize data object
    aData.mCandidateTriangles.resize(aData.mCandK.n_elem, 1);

    // link to current object
    arma::Mat<moris::uint> &tCand = aData.mCandidateTriangles.data();

    // write data
    tCand = arma::conv_to<arma::Mat<moris::uint> >::from(aData.mCandK);

#else
    // loop over all triangles in J-Direction
    moris::uint tCountJ = 0;
    for (moris::uint k = 0; k<aData.mNumberOfTriangles; ++k)
    {
        // check bounding box in J-direction
        if ( (aNodeCoords(2) - aData.mTriangleMinCoordsZ(k)) *
                             (aData.mTriangleMaxCoordsZ(k) - aNodeCoords(2)) > -MORIS_GE_EPSILON )
        {
            // remember this triangle
            aData.mCandJ(tCountJ) = k;

            // increment counter
            ++tCountJ;
        }
    }

    // counter for triangles
    moris::uint tCount = 0;

    // reset candidate size
    aData.mCandidateTriangles.resize(aData.mNumberOfTriangles, 1);

    // loop over remaining triangles in I-direction
    for (moris::uint k = 0; k<tCountJ; ++k)
    {
        // check bounding box in I-direction
        if((aNodeCoords(1) - aData.mTriangleMinCoordsY(aData.mCandJ(k)))*
        (aData.mTriangleMaxCoordsY(aData.mCandJ(k)) - aNodeCoords(1)) > -MORIS_GE_EPSILON )
        {
            aData.mCandidateTriangles(tCount) = aData.mCandJ(k);
            ++tCount;
        }
    }

    aData.mCandidateTriangles.resize(tCount, 1);
#endif
}

// -----------------------------------------------------------------------------
void
ge::SDF_Core::preselect_triangles_y(
        ge::SDF_Data                  & aData,
        const moris::Mat<moris::real> & aNodeCoords )
{
    // y: k = y, j = x, i = z
#ifdef MORIS_USE_ARMA
    // check bounding box in J-direction
    aData.mCandJ = arma::find(
            (aNodeCoords(0) - aData.mTriangleMinCoordsX) %
            (aData.mTriangleMaxCoordsX - aNodeCoords(0)) > -MORIS_GE_EPSILON );

    // check bounding box in I-direction
    aData.mCandI = arma::find(
            (aNodeCoords(2) - aData.mTriangleMinCoordsZ.elem(aData.mCandJ)) %
            (aData.mTriangleMaxCoordsZ.elem(aData.mCandJ) - aNodeCoords(2)) > -MORIS_GE_EPSILON );

    // help vector to be written in aData.mCandidateTriangles.data()
    aData.mCandK = aData.mCandJ.elem(aData.mCandI);

    // resize data object
    aData.mCandidateTriangles.resize(aData.mCandK.n_elem, 1);

    // link to current object
    arma::Mat<moris::uint> &tCand = aData.mCandidateTriangles.data();

    // write data
    tCand = arma::conv_to<arma::Mat<moris::uint> >::from(aData.mCandK);

#else
    // loop over all triangles in J-Direction
    moris::uint tCountJ = 0;
    for (moris::uint k = 0; k<aData.mNumberOfTriangles; ++k)
    {
        // check bounding box in J-direction
        if ((aNodeCoords(0) - aData.mTriangleMinCoordsX(k)) *
                             (aData.mTriangleMaxCoordsX(k) - aNodeCoords(0)) > -MORIS_GE_EPSILON )
        {
            // remember this triangle
            aData.mCandJ(tCountJ) = k;

            // increment counter
            ++tCountJ;
        }
    }

    // counter for triangles
    moris::uint tCount = 0;

    // reset candidate size
    aData.mCandidateTriangles.resize(aData.mNumberOfTriangles, 1);

    // loop over remaining triangles in I-direction
    for (moris::uint k = 0; k<tCountJ; ++k)
    {
        // check bounding box in I-direction
        if((aNodeCoords(2) - aData.mTriangleMinCoordsZ(aData.mCandJ(k)))*
        (aData.mTriangleMaxCoordsZ(aData.mCandJ(k)) - aNodeCoords(2)) > -MORIS_GE_EPSILON )
        {
            aData.mCandidateTriangles(tCount) = aData.mCandJ(k);
            ++tCount;
        }
    }

    aData.mCandidateTriangles.resize(tCount, 1);
#endif
}

// -----------------------------------------------------------------------------

void
ge::SDF_Core::preselect_triangles_z(
        ge::SDF_Data                  & aData,
        const moris::Mat<moris::real> & aNodeCoords)
{
    // z: k = z, j = y, i = x
#ifdef MORIS_USE_ARMA

    //moris::bool_t tNothingFound = true;

    // check bounding box in J-direction
    aData.mCandJ = arma::find(
            (aNodeCoords(1) - aData.mTriangleMinCoordsY) %
            (aData.mTriangleMaxCoordsY - aNodeCoords(1)) > -MORIS_GE_EPSILON );
    // check bounding box in I-direction
    aData.mCandI = arma::find(
            (aNodeCoords(0) - aData.mTriangleMinCoordsX.elem(aData.mCandJ)) %
            (aData.mTriangleMaxCoordsX.elem(aData.mCandJ) - aNodeCoords(0)) > -MORIS_GE_EPSILON );

    // help vector to be written in aData.mCandidateTriangles.data()
    aData.mCandK = aData.mCandJ.elem(aData.mCandI);

    // resize data object
    aData.mCandidateTriangles.resize(aData.mCandK.n_elem, 1);

    // link to current object
    arma::Mat<moris::uint> &tCand = aData.mCandidateTriangles.data();

    // write data
    tCand = arma::conv_to<arma::Mat<moris::uint> >::from(aData.mCandK);
#else
    // loop over all triangles in J-Direction
    moris::uint tCountJ = 0;
    for (moris::uint k = 0; k<aData.mNumberOfTriangles; ++k)
    {
        // check bounding box in J-direction
        if( (aNodeCoords(1) - aData.mTriangleMinCoordsY(k)) *
                             (aData.mTriangleMaxCoordsY(k) - aNodeCoords(1)) > -MORIS_GE_EPSILON )
        {
            // remember this triangle
            aData.mCandJ(tCountJ) = k;

            // increment counter
            ++tCountJ;
        }
    }

    // counter for triangles
    moris::uint tCount = 0;

    // reset candidate size
    aData.mCandidateTriangles.resize(aData.mNumberOfTriangles, 1);

    // loop over remaining triangles in I-direction
    for (moris::uint k = 0; k<tCountJ; ++k)
    {
        // check bounding box in I-direction
        if((aNodeCoords(0) - aData.mTriangleMinCoordsX(aData.mCandJ(k)))*
        (aData.mTriangleMaxCoordsX(aData.mCandJ(k)) - aNodeCoords(0)) > -MORIS_GE_EPSILON )
        {
            aData.mCandidateTriangles(tCount) = aData.mCandJ(k);
            ++tCount;
        }
    }

    aData.mCandidateTriangles.resize(tCount, 1);
#endif
}

// -----------------------------------------------------------------------------

void
ge::SDF_Core::intersect_triangles(
        ge::SDF_Data                  & aData,
        const moris::uint               aAxis,
        const moris::Mat<moris::real> & aNodeCoords)
{
    // counter for intersected triangles
    moris::uint tCount = aData.mCandidateTriangles.length();
    aData.mIntersectedTriangles.resize(tCount, 1);
    moris::uint tThisTriangle = 0;

    // reset counter
    tCount = 0;

    // loop over all candidate triangles, and check intersection criterion
    for (moris:: uint k=0; k < aData.mCandidateTriangles.length(); ++k)
    {
        tThisTriangle = aData.mCandidateTriangles(k);

        if ( aData.mTriangles(tThisTriangle).check_edge( 0, aAxis, aNodeCoords ))
        {
            if ( aData.mTriangles(tThisTriangle).check_edge( 1, aAxis, aNodeCoords ))
            {
                if ( aData.mTriangles(tThisTriangle).check_edge( 2, aAxis, aNodeCoords ))
                {
                    aData.mIntersectedTriangles(tCount) = tThisTriangle;
                    ++tCount;
                }
            }
        }
    }

    // adapt size of list
    aData.mIntersectedTriangles.resize(tCount, 1);
}

// -----------------------------------------------------------------------------

void
ge::SDF_Core::intersect_ray_with_triangles(
        ge::SDF_Data                  & aData,
        const moris::uint               aAxis,
        const moris::Mat<moris::real> & aNodeCoords)
{
    // resize length of vector
    aData.mCoordsK.resize(aData.mIntersectedTriangles.length(), 1);

    // counter for coordinates
    //moris::uint tCount = 0;

    // temporary variable for K-Coordinate
    moris::real tKCoordinate;

    moris::uint tCount = 0;
    //  Find the k-coordinate of the locations where the ray crosses each triangle:
    for (moris::uint k = 0; k < aData.mIntersectedTriangles.length(); ++k)
    {
        moris::uint tThisTriangle = aData.mIntersectedTriangles(k);

        bool tError;

        // calculate intersection coordinate
        aData.mTriangles(tThisTriangle).intersect_with_coordinate_axis(
                aNodeCoords,
                aAxis,
                tKCoordinate,
                tError );

        // error meant we would have divided by zero. This triangle is ignored
        // otherwose, the value is written into the result vector
        if ( ! tError )
        {
            aData.mCoordsK( tCount++ ) = ge::round( tKCoordinate / MORIS_GE_EPSILON )* MORIS_GE_EPSILON;
        }
    }

    // resize coordinate array
    aData.mCoordsK.resize( tCount, 1);

    if ( tCount > 0)
    {
        // find unique entries
        aData.mCoordsK = ge::unique(aData.mCoordsK);
    }

}
// -----------------------------------------------------------------------------

void
ge::SDF_Core::check_if_node_is_inside(
        const ge::SDF_Mesh_Data         & aMeshData,
        ge::SDF_Data                    & aData,
        const moris::uint                 aAxis,
        const moris::uint                 aLocalNodeInd,
        const moris::Mat<moris::real>   & aNodeCoords)
{

    moris::uint tNumCoordsK = aData.mCoordsK.length();

    moris::bool_t tNodeIsInside = false;

    // only even number of intersections is considered
    if(tNumCoordsK % 2 == 0)
    {
        // loop over k-Coordinates.
        for (moris::uint k=0; k< tNumCoordsK / 2 ; ++k)
        {
            tNodeIsInside = (aNodeCoords( aAxis ) > aData.mCoordsK( 2 * k )) &&
                            (aNodeCoords( aAxis ) < aData.mCoordsK( 2 * k + 1 ));

            // break the loop if inside
            if ( tNodeIsInside )
            {
                break;
            }
        }

        // set the inside flag of this node to the corresponding value
        if( tNodeIsInside )
        {
            aData.mLocalNodeInsideFlags.set(aLocalNodeInd);
        }
        else
        {
            aData.mLocalNodeInsideFlags.reset(aLocalNodeInd);
        }
    }
    else
    {
        // add to unsure node list
        aData.mUnsureNodesNew(aData.mUnsureNewNodesCount) = aLocalNodeInd;

        // increment unsure nodes counter
        ++aData.mUnsureNewNodesCount;
    }
}

// -----------------------------------------------------------------------------

void
ge::SDF_Core::calculate_candidate_points_and_buffer_diagonal(
        const SDF_Mesh_Data & aMeshData,
        SDF_Data            & aData)
{
    // counter for elements near surface
    moris::uint tSurfaceElementCounter = 0;

    // counter for elements in volume
    moris::uint tVolumeElementCounter = 0;

    moris::real tBufferDiagonal = 0;

    // search all elements for sign change
    for ( moris::uint e=0; e < aMeshData.get_number_of_elements(); ++e )
    {
        // get the nodes of element e
        moris::Mat< moris::uint > tNodesOfThisElement = aMeshData.get_nodes_of_element(e);

        // we check if there is a sign change within the element (e).
        // If so, we flag all nodes of this element as candidate flags.
        moris::bool_t tFirstSign = aData.mLocalNodeInsideFlags.test(tNodesOfThisElement( 0 ));
        moris::bool_t tSameSign = true;

        for ( moris::uint k=1; k<tNodesOfThisElement.length(); ++k )
        {
            // check if sign is the same
            tSameSign = aData.mLocalNodeInsideFlags.test( tNodesOfThisElement( k ) ) ==  tFirstSign;

            // if so, break the loop
            if( !tSameSign )
            {
                break;
            }
        }

        // if there is a sign change
        if( !tSameSign )
        {
            // flag this element as surface element
            aData.mLocalElementsAtSurface(tSurfaceElementCounter) = e;

            // increment counter
            ++tSurfaceElementCounter;

            // also flag all nodes of this element as candidates
            for (moris::uint k = 0; k < tNodesOfThisElement.length(); ++k)
            {
                aData.mLocalNodeCandidateFlags.set(tNodesOfThisElement(k));
            }

            // update buffer diagonal
            tBufferDiagonal = ge::max(tBufferDiagonal, get_diagonal_length_of_element(aMeshData, tNodesOfThisElement));
        }
        else if (tFirstSign)
        {
            // this element is in the volume
            aData.mLocalElementsInVolume(tVolumeElementCounter) = e;

            // increment counter
            ++tVolumeElementCounter;
        }
    }

    // adapt array sizes
    aData.mLocalElementsAtSurface.resize(tSurfaceElementCounter, 1);
    aData.mLocalElementsInVolume.resize(tVolumeElementCounter, 1);

    // addidional search depth
    if( aData.mSettings.mCandidateSearchDepth > 0 )
    {
        for (moris::uint d = 1; d <= aData.mSettings.mCandidateSearchDepth; ++d)
        {
            // copy flags
            aData.mLocalNodeCandidateFlagsOld = aData.mLocalNodeCandidateFlags;

            // loop over all elements
            for (moris::uint e = 0; e < aMeshData.get_number_of_elements(); ++e)
            {
                // get the nodes of element e
                moris::Mat< moris::uint > tNodesOfThisElement = aMeshData.get_nodes_of_element(e);

                // switch telling if one of the nodes of this element is flagged
                moris::bool_t tCandidateSwitch = false;

                // loop over all nodes of this element
                for ( moris::uint k=0; k<tNodesOfThisElement.length(); ++k)
                {
                    // check if any node of this cell is flagged
                    tCandidateSwitch = aData.mLocalNodeCandidateFlagsOld.test(tNodesOfThisElement( k ));

                    // if so, exit the loop
                    if ( tCandidateSwitch )
                    {
                        break;
                    }
                }

                if (tCandidateSwitch)
                {
                    // flag all nodes of this cell
                    for ( moris::uint k=0; k<tNodesOfThisElement.length(); ++k)
                    {
                        aData.mLocalNodeCandidateFlags.set(tNodesOfThisElement( k ));
                    }

                    // get local node numbers
                    moris::Mat< moris::uint > tNodesOfThisElement = aMeshData.get_nodes_of_element(e);

                    // update buffer diagonal
                    tBufferDiagonal = ge::max(tBufferDiagonal,
                                              get_diagonal_length_of_element(aMeshData, tNodesOfThisElement));
                }

            }
        }
    }

    // save buffer diagonal to data
    aData.mBufferDiagonal = tBufferDiagonal*( 1 + aData.mSettings.mBufferDiagonalEpsilon )
                                           *( 1 + aData.mSettings.mCandidateSearchDepth );
}
// -----------------------------------------------------------------------------

moris::real
ge::SDF_Core::get_diagonal_length_of_element (
        const ge::SDF_Mesh_Data       & aMeshData,
        const moris::Mat<moris::uint> & aNodesOfElement) const
{
    // reset min and max coordinate
    moris::Mat<moris::real> tMinCoordinate(3,1, MORIS_GE_HUGE);
    moris::Mat<moris::real> tMaxCoordinate(3,1,-MORIS_GE_HUGE);

    // loop over all nodes
    for(moris::uint k=0; k<aNodesOfElement.length(); ++k)
    {

        // get coordinates of this node
        moris::Mat< moris::real > tNodeCoords = aMeshData.get_node_coords( aNodesOfElement( k ) );

        // now update min and max value
        for(moris::uint i=0; i<3; ++i)
        {
            tMinCoordinate( i ) = ge::min(tMinCoordinate( i ), tNodeCoords( i ) );
            tMaxCoordinate( i ) = ge::max(tMaxCoordinate( i ), tNodeCoords( i ) );
        }
    }

    // calculate the diagonal
    return ge::point_to_point_distance( tMinCoordinate, tMaxCoordinate );
}

// -----------------------------------------------------------------------------

void
ge::SDF_Core::calculate_udf(
        const ge::SDF_Mesh_Data & aMeshData,
                   ge::SDF_Data & aData)
{

    // reset SDF value
    aData.mLocalSDF.fill(MORIS_GE_HUGE);

    // convert candidate bitset to matrix with indices (faster)
    aData.mLocalCandidateNodes = ge::BitsetToMat(aData.mLocalNodeCandidateFlags);

    // list of nodes to calculate
    moris::Mat< moris::uint > tNodesWithinTriangle;

    // coordinate of node
    moris::Mat< moris::real > tCoordsOfThisNode(3,1);
    // loop over all triangles
    for ( moris::uint t=0; t<aData.mNumberOfTriangles; ++t )
    {
        // ask triangle for local nodes
        tNodesWithinTriangle = get_nodes_within_triangle_bounding_box(aMeshData, aData, t);

        for ( moris::uint k=0; k < tNodesWithinTriangle.length(); ++k )
        {
            // get node number
            moris::uint tThisNode = tNodesWithinTriangle(k);

            // get node coordinates
            tCoordsOfThisNode = aMeshData.get_node_coords(tThisNode);

            // calculate triangle distance
            moris::real tDistance = aData.mTriangles( t ).get_distance_to_point( tCoordsOfThisNode );

            // set active not flag in order to remember that we have calculated this node (debugging only)
            aData.mLocalNodeHasSdfFlag.set( tThisNode );

            // check if distance is smaller then the last calculated
            if (tDistance < aData.mLocalSDF( tThisNode ))
            {
                // write data to array
                aData.mLocalSDF(tThisNode) = tDistance;
            }

        }
    }
}

// -----------------------------------------------------------------------------

moris::Mat<moris::uint>
ge::SDF_Core::get_nodes_within_triangle_bounding_box(
        const ge::SDF_Mesh_Data & aMeshData,
        ge::SDF_Data            & aData,
        const moris::uint         aTriangle)
{
    // reserve memory for max and min coordinate within triangle
    moris::Mat< moris::real > tMinPoint( 3, 1 );
    moris::Mat< moris::real > tMaxPoint( 3, 1 );

    // get min coordinate from mesh
    moris::Mat < moris::real > tMinMeshCoord = aMeshData.get_min_coord();

    // get max coordinate from mesh
    moris::Mat < moris::real > tMaxMeshCoord = aMeshData.get_max_coord();

    // calculate minimum node coordinates
    tMinPoint( 0 ) = ge::max(aData.mTriangleMinCoordsX(aTriangle)-aData.mBufferDiagonal,
                             tMinMeshCoord( 0 ));
    tMinPoint( 1 ) = ge::max(aData.mTriangleMinCoordsY(aTriangle)-aData.mBufferDiagonal,
                             tMinMeshCoord( 1 ));
    tMinPoint( 2 ) = ge::max(aData.mTriangleMinCoordsZ(aTriangle)-aData.mBufferDiagonal,
                             tMinMeshCoord( 2 ));

    // calculate maximum node coordinates
    tMaxPoint( 0 ) = ge::min(aData.mTriangleMaxCoordsX(aTriangle)+aData.mBufferDiagonal,
                             tMaxMeshCoord( 0 ));
    tMaxPoint( 1 ) = ge::min(aData.mTriangleMaxCoordsY(aTriangle)+aData.mBufferDiagonal,
                             tMaxMeshCoord( 1 ));
    tMaxPoint( 2 ) = ge::min(aData.mTriangleMaxCoordsZ(aTriangle)+aData.mBufferDiagonal,
                             tMaxMeshCoord( 2 ));
    // reserve memory for relevant nodes
    moris::Mat< moris::uint > aNodes( aData.mLocalCandidateNodes.length(), 1 );
    // counter for aNodes
    moris::uint tNodeCount = 0;
    
    // local indes for this node
    moris::uint tThisNode;

    // coordiantes of current node
    moris::Mat< moris::real > tCoordsOfThisNode(3,1);

    // switch to tell of node is within triangle vicinity
    moris::bool_t tNodeIsWithinTriangle;

    // loop over all  candidate nodes
    for (moris::uint k=0; k < aData.mLocalCandidateNodes.length(); ++k)
    {
        // copy node ID
        tThisNode = aData.mLocalCandidateNodes(k);
        // check if node is in triangle
        tNodeIsWithinTriangle = true;

        // get coordinate of this node
        tCoordsOfThisNode = aMeshData.get_node_coords(tThisNode);

        // loop over all dimensions
        for( moris::uint i=0; i<3; ++i )
        {
            tNodeIsWithinTriangle = tNodeIsWithinTriangle &&
                                    (tCoordsOfThisNode( i ) >= tMinPoint( i )) &&
                                    (tCoordsOfThisNode( i ) <= tMaxPoint( i ));
            // exit loop if false
            if( ! tNodeIsWithinTriangle )
            {
                break;
            }
        }
        // check of node is inside of the triangle
        if ( tNodeIsWithinTriangle )
        {
            // append triangle to list if node is inside
            aNodes(tNodeCount) = tThisNode;

            // increment counter
            ++tNodeCount;
        }

    }

    // correct size of candidate list
    aNodes.resize(tNodeCount, 1);

    // return list
    return aNodes;
}

// -----------------------------------------------------------------------------

void
ge::SDF_Core::sign_udf(
        const SDF_Mesh_Data & aMeshData,
              SDF_Data      & aData)
{
    // calculate minumum and maximum value of sdf
    moris::real tMinSDF = 0;
    moris::real tMaxSDF = 0;
    for (moris::uint k=0; k<aMeshData.get_number_of_nodes(); ++k)
    {
        // check if value was calculated
        if (aData.mLocalNodeHasSdfFlag.test(k))
        {
            if (aData.mLocalNodeInsideFlags.test(k))
            {
                // remember minimum value
                tMinSDF = ge::max(tMinSDF, aData.mLocalSDF(k));
            } else
            {
                // remember maximum value
                tMaxSDF = ge::max(tMaxSDF, aData.mLocalSDF(k));
            }
        }
    }

    // set pseudo value for tMaxSDF if domain on proc is empty
    if (tMaxSDF == 0)
    {
        tMaxSDF = 0.001;
    }
    // flip sign of minimal value
    tMinSDF *= -1;

    // add extra value to min and max value
    tMaxSDF += aData.mBufferDiagonal;
    tMinSDF -= aData.mBufferDiagonal;

    // loop over all nodes
    for (moris::uint k=0; k<aMeshData.get_number_of_nodes(); ++k)
    {
        // check if value was calculated
        if (aData.mLocalNodeHasSdfFlag.test(k))
        {
            // check if node is inside
            if (aData.mLocalNodeInsideFlags.test(k))
            {
                // invert sign of value
                aData.mLocalSDF(k) *= -1;
            }
        }
        else
        {
            if ( aData.mLocalNodeInsideFlags.test(k) )
            {
                // node is inside, set pseudo value
                aData.mLocalSDF(k) = tMinSDF;
            }
            else
            {
                // node is outside, set pseudo value
                aData.mLocalSDF(k) = tMaxSDF;
            }
        }

    }
}

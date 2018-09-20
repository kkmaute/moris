#include "cl_Ge_SDF_Triangle_File.hpp"

// =============================================================================
//  Public
// =============================================================================

ge::SDF_Triangle_File::SDF_Triangle_File(const std::string& aFilePath)
{

    // check the file extension
    auto tFileExt = aFilePath.substr(aFilePath.find_last_of(".")+1,aFilePath.length());

    if (tFileExt == "obj")
    {
        this->load_from_obj_file(aFilePath);
    }
    else if(tFileExt == "msh"){
        this->load_from_msh_file(aFilePath);
    }
    else
    {
        throw std::invalid_argument( "SDF_Triangle_File needs an obj or msh file for creation." );
    }
}

// -----------------------------------------------------------------------------

void ge::SDF_Triangle_File::save_to_file(const std::string& aFilePath)
{

    // check the file extension
    auto tFileExt = aFilePath.substr(aFilePath.find_last_of(".")+1,
                                     aFilePath.length());
    if (tFileExt == "msh")
    {
        save_to_msh(aFilePath);
    }
    else if(tFileExt == "vtk")
    {
        save_to_vtk(aFilePath);
    }
    else if(tFileExt == "obj")
    {
        save_to_obj(aFilePath);
    }
    else{
        throw std::invalid_argument( "SDF_Triangle_File can only save msh and vtk files." );
    }

}

// =============================================================================
//  Private
// =============================================================================

void ge::SDF_Triangle_File::load_ascii_to_buffer(const std::string& aFilePath,
                                            moris::Cell<std::string>& aBuffer)
{
    // try to open ascii file
    std::ifstream tAsciiFile(aFilePath);
    std::string tLine;

    // load file into buffer, otherwise throw error
    if(tAsciiFile)
    {

        // count number of lines
        moris::uint tBufferLength = 0;
        while (!tAsciiFile.eof())
        {
            std::getline(tAsciiFile, tLine);
            ++tBufferLength;
        }
        tAsciiFile.close();
        tAsciiFile.open(aFilePath);

        // load file into buffer
        aBuffer.reserve(tBufferLength);
        while (!tAsciiFile.eof())
        {
            std::getline(tAsciiFile, tLine);
            aBuffer.push_back(tLine);
        }
        tAsciiFile.close();
    }
    else
    {
        std::cerr << "Something went wrong while trying to load from " <<
                  aFilePath << "." << std::endl;
    }
}

// -----------------------------------------------------------------------------

void ge::SDF_Triangle_File::load_from_obj_file(const std::string& aFilePath)
{
    // temporary variable containing the current line in the ascii file
    std::string tLine;

    // temporary buffer which contains the lines
    moris::Cell<std::string> tBuffer;

    // load asci file
    load_ascii_to_buffer(aFilePath, tBuffer);

    // count buffer length
    moris::uint tBufferLength = tBuffer.size();
    // count number of nodes and triangles
    mNumberOfNodes = 0;
    mNumberOfTriangles = 0;

    // loop over all lines
    for (uint k=0; k<tBufferLength; ++k)
    {
        if( tBuffer( k ).substr( 0, 2 ) == "v ")
        {
            ++mNumberOfNodes;
        }
        else if ( tBuffer( k ).substr( 0, 2 ) == "f ")
        {
            ++mNumberOfTriangles;
        }
    }
    // set size of node array
    mNodeCoords.set_size(3, mNumberOfNodes);
    // set size of triangle array
    mTriangles.set_size(3, mNumberOfTriangles);

    // vertex counter
    uint v=0;

    // triangle counter
    uint t=0;

    // x, y and z-coordinates for node
    float tX = 0;
    float tY = 0;
    float tZ = 0;

    // temporary one-based Ids for triangle nodes 1, 2 and 3
    uint tNode0 = 0;
    uint tNode1 = 0;
    uint tNode2 = 0;

    // loop over all lines
    for (uint k=0; k<tBufferLength; ++k)
    {
        if( tBuffer( k ).substr( 0, 2 ) == "v ")
        {
            // read vertex (node)

            // read ascii data into coordinates
            std::sscanf(tBuffer(k).substr(2,tBuffer(k).length()).c_str(), "%f %f %f", &tX, &tY, &tZ);

            // mSettings.mMeshHighPass:
            // ignore very small values, which are usually caused by machine precision
            if(ge::abs(tX) > mSettings.mMeshHighPass)
                mNodeCoords(0,v) = ( moris::real ) tX;
            else
                mNodeCoords(0,v) = 0.0;

            if(ge::abs(tY) > mSettings.mMeshHighPass )
                mNodeCoords(1,v) = ( moris::real ) tY;
            else
                mNodeCoords(1,v) = 0.0;

            if(ge::abs(tZ) > mSettings.mMeshHighPass )
                mNodeCoords(2,v) = ( moris::real ) tZ;
            else
                mNodeCoords(2,v) = 0.0;

            ++v;
        }
        else if ( tBuffer( k ).substr( 0, 2 ) == "f ")
        {
            // read face (triangle)

            // read triangle topology
            std::sscanf(tBuffer(k).substr(2,tBuffer(k).length()).c_str(),
                    "%u %u %u", &tNode0, &tNode1, &tNode2);

            // save data into array
            mTriangles(0,t) = (moris::uint) tNode0-1;
            mTriangles(1,t) = (moris::uint) tNode1-1;
            mTriangles(2,t) = (moris::uint) tNode2-1;

            ++t;
        }
    }

    // tidy up
    tBuffer.clear();
}

// -----------------------------------------------------------------------------

void ge::SDF_Triangle_File::load_from_msh_file(const std::string& aFilePath){

    // load ascii file into buffer
    moris::Cell<std::string> tBuffer;
    load_ascii_to_buffer(aFilePath, tBuffer);
    moris::uint tBufferSize = tBuffer.size();

    // count nodes
    moris::bool_t run = true;
    moris::uint tFlag = 0;
    while(run)
    {
        run = tBuffer(tFlag) != "$Nodes";
        run = run && tFlag < tBufferSize;
        tFlag++;
    }
    mNumberOfNodes = (moris::uint) std::stoul(tBuffer(tFlag));
    tFlag++;
    moris::uint tNodeFlag = tFlag;

    // count Elements
    run = true;
    while(run)
    {
        run = tBuffer(tFlag) != "$Elements";
        run = run && tFlag < tBufferSize;
        tFlag++;
    }

    moris::uint tNumberOfElements = (moris::uint) std::stoul(tBuffer(tFlag));
    tFlag++;
    mNumberOfTriangles = 0;
    moris::uint tPivot = 0;
    moris::uint tOldPivot = 0;
    moris::uint tNumberOfQuadrangles = 0;

    // loop over all elements to determine data type
    for (moris::uint k = tFlag; k < tFlag+tNumberOfElements; ++k)
    {
        tPivot = tBuffer(k).find(" ");
        for (moris::uint i = 0; i < 1; ++i)
        {
            tOldPivot = tPivot+1;
            tPivot = tBuffer(k).find(" ", tOldPivot);
        }

        // count triangles and quadrangles
        if(tBuffer(k).substr(tOldPivot,tPivot-tOldPivot) == "2")
        {
            mNumberOfTriangles++;
        }
        else if(tBuffer(k).substr(tOldPivot,tPivot-tOldPivot) == "3")
        {
            tNumberOfQuadrangles++;
        }
    }
    std::size_t tElementFlag = tFlag;

    if(mSettings.mQuadSplitMode == 1 || mSettings.mQuadSplitMode == 2)
    {
        mNumberOfTriangles += 2*tNumberOfQuadrangles;
        mNodeCoords.set_size(3, mNumberOfNodes);
    }
    else if(mSettings.mQuadSplitMode == 4)
    {
        mNumberOfTriangles += 4 * tNumberOfQuadrangles;
        mNodeCoords.set_size(3, mNumberOfNodes + tNumberOfQuadrangles);
    } else {
        throw std::invalid_argument( "SDF_Triangle_File->Settings mQuadSplitMode must be either 1, 2 or 4" );
    }

    // read nodes
    moris::uint tNode = 0;
    tFlag = tNodeFlag;
    for (moris::uint k = tFlag; k < tFlag+mNumberOfNodes; ++k)
    {
        tPivot = tBuffer(k).find(" ");
        tNode = std::stoul(tBuffer(k).substr(0,tPivot))-1;
        for (moris::uint i = 0; i < 3; ++i)
        {
            tOldPivot = tPivot+1;
            tPivot = tBuffer(k).find(" ", tOldPivot);
            mNodeCoords(i,tNode) = std::stod(tBuffer(k).substr(tOldPivot,tPivot-tOldPivot));
            // highpass
            if(std::abs(mNodeCoords(i,tNode)) < mSettings.mMeshHighPass)
                mNodeCoords(i,tNode) = 0;
        }
    }

    // set size of Triangle array
    mTriangles.set_size(3, mNumberOfTriangles);

    // read triangles and quadrangles
    moris::uint tTriangle = 0;
    moris::uint tParam = 0;
    moris::Matrix< moris::DDUMat > tQuad(4,1);
    moris::Matrix< moris::DDRMat > tCenter(3,1);

    // loop over all lines containing element topology
    tFlag = tElementFlag;
    for (moris::uint k = tFlag; k < tFlag+tNumberOfElements; ++k)
    {
        tPivot = tBuffer(k).find(" ");
        for (moris::uint i = 0; i < 1; ++i){
            tOldPivot = tPivot+1;
            tPivot = tBuffer(k).find(" ", tOldPivot);
        }
        if(tBuffer(k).substr(tOldPivot,tPivot-tOldPivot) == "2")
        {
            // read triangle topology

            tOldPivot = tPivot+1;
            tPivot = tBuffer(k).find(" ", tOldPivot);
            tParam = std::stoi(tBuffer(k).substr(tOldPivot,tPivot-tOldPivot));
            for (moris::uint i = 0; i < tParam; ++i)
            {
                tOldPivot = tPivot+1;
                tPivot = tBuffer(k).find(" ", tOldPivot);
            }
            for (moris::uint i = 0; i < 3; ++i)
            {
                tOldPivot = tPivot+1;
                tPivot = tBuffer(k).find(" ", tOldPivot);
                tNode = std::stoi(tBuffer(k).substr(tOldPivot,tPivot-tOldPivot))-1;
                mTriangles(i,tTriangle) = tNode;
            }
            tTriangle++;
        }
        else if(tBuffer(k).substr(tOldPivot,tPivot-tOldPivot) == "3")
        {
            // read quadrangle topology and split it into triangles
            tOldPivot = tPivot+1;
            tPivot = tBuffer(k).find(" ", tOldPivot);
            tParam = std::stoi(tBuffer(k).substr(tOldPivot,tPivot-tOldPivot));

            for (moris::uint i = 0; i < tParam; ++i)
            {
                tOldPivot = tPivot+1;
                tPivot = tBuffer(k).find(" ", tOldPivot);
            }
            for (moris::uint i = 0; i < 4; ++i)
            {
                tOldPivot = tPivot+1;
                tPivot = tBuffer(k).find(" ", tOldPivot);
                tNode = std::stoi(tBuffer(k).substr(tOldPivot,tPivot-tOldPivot))-1;
                tQuad(i) = tNode;
            }

            if(mSettings.mQuadSplitMode == 1)
            {
                // split quad into two triangles
                mTriangles(0, tTriangle) = tQuad(0);
                mTriangles(1, tTriangle) = tQuad(1);
                mTriangles(2, tTriangle) = tQuad(2);
                tTriangle++;

                mTriangles(0, tTriangle) = tQuad(1);
                mTriangles(1, tTriangle) = tQuad(2);
                mTriangles(2, tTriangle) = tQuad(4);
                tTriangle++;
            }
            else if(mSettings.mQuadSplitMode == 2)
            {
                // split quad into two triangles
                mTriangles(0, tTriangle) = tQuad(0);
                mTriangles(1, tTriangle) = tQuad(1);
                mTriangles(2, tTriangle) = tQuad(3);
                tTriangle++;

                mTriangles(0, tTriangle) = tQuad(1);
                mTriangles(1, tTriangle) = tQuad(2);
                mTriangles(2, tTriangle) = tQuad(3);
                tTriangle++;
            }
            else if(mSettings.mQuadSplitMode == 4)
            {
                // calculate center node
                tCenter.fill(0);
                for (moris::uint i = 0; i < 4; ++i)
                {
                    for (moris::uint j = 0; j < 3; ++j)
                        tCenter(j) += mNodeCoords(j, tQuad(i));
                }
                for (moris::uint i = 0; i < 4; ++i)
                    tCenter(i) *= 0.25;

                for (moris::uint j = 0; j < 3; ++j)
                    if (std::abs(tCenter(j)) < mSettings.mMeshHighPass)
                        tCenter(j) = 0;

                // add center node to nodelist
                for (moris::uint j = 0; j < 3; ++j)
                    mNodeCoords(j, mNumberOfNodes) = tCenter(j);

                // split quadrangle into four triangles

                mTriangles(0, tTriangle) = tQuad(0);
                mTriangles(1, tTriangle) = tQuad(1);
                mTriangles(2, tTriangle) = mNumberOfNodes;
                tTriangle++;

                mTriangles(0, tTriangle) = tQuad(1);
                mTriangles(1, tTriangle) = tQuad(2);
                mTriangles(2, tTriangle) = mNumberOfNodes;
                tTriangle++;

                mTriangles(0, tTriangle) = tQuad(2);
                mTriangles(1, tTriangle) = tQuad(3);
                mTriangles(2, tTriangle) = mNumberOfNodes;
                tTriangle++;

                mTriangles(0, tTriangle) = tQuad(3);
                mTriangles(1, tTriangle) = tQuad(0);
                mTriangles(2, tTriangle) = mNumberOfNodes;
                tTriangle++;

                // increment node counter
                mNumberOfNodes++;
            }
            else
            {
                throw std::invalid_argument(
                        "SDF_Triangle_File->mSettings->mQuadSplitMode must be either 1, 2, or 4");
            }
        }
    }
}

// -----------------------------------------------------------------------------

void ge::SDF_Triangle_File::save_to_msh(const std::string& aFilePath)
{
    std::ofstream tFile(aFilePath);
    tFile << "$MeshFormat" << std::endl;
    tFile << "2.2 0 8" << std::endl;
    tFile << "$EndMeshFormat" << std::endl;
    tFile << "$Nodes" << std::endl;
    tFile << mNumberOfNodes << std::endl;
    for (moris::uint k = 0; k < mNumberOfNodes; ++k)
    {
        tFile << k+1 << " " << mNodeCoords(0,k)
              << " " << mNodeCoords(1,k)
              << " " << mNodeCoords(2,k) << std::endl;
    }
    tFile << "$EndNodes" << std::endl;
    tFile << "$Elements" << std::endl;
    tFile << mNumberOfTriangles << std::endl;
    for (moris::uint k = 0; k < mNumberOfTriangles; ++k)
    {
        tFile << k+1 << " 2 2 0 1 " << mTriangles(0,k)+1
              << " "  << mTriangles(1,k)+1
              << " "  << mTriangles(2,k)+1 << std::endl;
    }
    tFile << "$EndElements" << std::endl;
    tFile.close();
}

// -----------------------------------------------------------------------------

void ge::SDF_Triangle_File::save_to_vtk(const std::string& aFilePath)
{
    std::ofstream tFile(aFilePath);

    tFile << "# vtk DataFile Version 3.0" << std::endl;
    tFile << "GO BUFFS!" << std::endl;

    if(mSettings.mSaveVTKasASCII){
        tFile << "ASCII" << std::endl;
        tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
        tFile << "POINTS " << mNumberOfNodes << " float" << std::endl;
        for (moris::uint k = 0; k < mNumberOfNodes; ++k)
        {
            tFile << mNodeCoords(0,k) << " "
                  << mNodeCoords(1,k) << " "
                  << mNodeCoords(2,k) << std::endl;
        }

        tFile << "CELLS " << mNumberOfTriangles << " " << 4*mNumberOfTriangles
              << std::endl;
        for (moris::uint k = 0; k < mNumberOfTriangles; ++k)
        {
            tFile << "3 " << mTriangles(0,k) << " "  << mTriangles(1,k) << " "
                  << mTriangles(2,k) << std::endl;
        }
        tFile << "CELL_TYPES " << mNumberOfTriangles << std::endl;
        for (moris::uint k = 0; k < mNumberOfTriangles; ++k)
        {
            tFile << "5" << std::endl;
        }
    }
    else
    {
        float tFValue = 0;
        int tIValue = 0;
        float tFChar = 0;
        moris::uint tIChar = 0;

        tFile << "BINARY" << std::endl;
        tFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
        tFile << "POINTS " << mNumberOfNodes << " float" << std::endl;
        for (moris::uint k = 0; k < mNumberOfNodes; ++k)
        {
            for (moris::uint i = 0; i < 3; ++i)
            {
                tFValue = (float) mNodeCoords(i, k);
                tFChar = ge::swap_byte_endian(tFValue);
                tFile.write((char *) &tFChar, sizeof(float));
            }
        }
        tFile << std::endl;
        tFile << "CELLS " << mNumberOfTriangles << " "
              << 4 * mNumberOfTriangles << std::endl;
        for (moris::uint k = 0; k < mNumberOfTriangles; ++k)
        {
            tIValue = (moris::uint) 3;
            tIChar = ge::swap_byte_endian(tIValue);
            tFile.write((char *) &tIChar, sizeof(moris::uint));
            for (moris::uint i = 0; i < 3; ++i) {
                tIChar = ge::swap_byte_endian(mTriangles(i, k));
                tFile.write((char *) &tIChar, sizeof(moris::uint));
            }
        }
        tFile << std::endl;
        tFile << "CELL_TYPES " << mNumberOfTriangles << std::endl;
        tIValue = (moris::uint) 5;
        tIChar = ge::swap_byte_endian(tIValue);
        for (moris::uint k = 0; k < mNumberOfTriangles; ++k)
        {
            tFile.write((char *) &tIChar, sizeof(moris::uint));
        }
    }
    tFile.close();
}

// -----------------------------------------------------------------------------

void ge::SDF_Triangle_File::save_to_obj(const std::string& aFilePath)
{
    std::ofstream tFile(aFilePath);
    tFile << "# " << mNumberOfNodes << std::endl;
    tFile << "# " << mNumberOfTriangles << std::endl;
    for (moris::uint k = 0; k < mNumberOfNodes; ++k)
    {
        tFile << "v " << mNodeCoords(0,k) << " "
              << mNodeCoords(1,k) << " "
              << mNodeCoords(2,k) << std::endl;
    }
    for (moris::uint k = 0; k < mNumberOfTriangles; ++k)
    {
        tFile << "f " << mTriangles(0,k)+1 << " " << mTriangles(1,k)+1 << " "
              << mTriangles(2,k)+1 << std::endl;
    }
    tFile.close();
}

// -----------------------------------------------------------------------------

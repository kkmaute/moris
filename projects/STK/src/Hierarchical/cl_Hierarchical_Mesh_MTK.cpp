/*
 * cl_Hierarchical_Mesh_MTK.cpp
 *
 *  Created on: Jan 4, 2018
 *      Author: gleim
 */

#include "cl_Hierarchical_Mesh_MTK.hpp" // STK/src/Hierarchical
using namespace moris;

void
Hierarchical_Mesh_MTK::create_MTK_file(
        std::string const & OutputFileName,
        uint            & aModelDim,
        Mat<uint>       & aNumberOfElementsPerDirection,
        Mat<uint>       & aFetopo,
        Mat<real>       & aFeCoord,
        Mat<uint>       & aBasisProcOwner,
        Mat<uint>       & aElementListOnProc,
        Mat<uint>       & aLagrangeToBSplineMap,
        Mat<uint>       & aLagrangeListOnProc,
        Cell<Mat<uint>> & aSideSet)
{
    // Create MORIS mesh using MTK database
    MtkSetsInfo aMeshSets;

    // Save SideSet
    MtkSideSetsInfo tSideSetStruc;
    tSideSetStruc.ElemIdsAndSideOrds = &aSideSet;
    Cell< std::string > tSetsName( aSideSet.size() );

    //Create names for the sets
    for ( uint i = 0; i < aSideSet.size(); i++ )
    {
        tSetsName(i) = "tSideSet_" + std::to_string( i );
    }
    tSideSetStruc.SSetNames   = tSetsName;
    if( aSideSet.size() > 0 )
    {
        aMeshSets.SideSetsInfo   = &tSideSetStruc;
    }

    //Create additional element information and save it in fields for output
    Mat<real> tElementLevel( aElementListOnProc.length(), 1, 0 );
    Mat<real> tElementId( aElementListOnProc.length(), 1, 0 );
    for ( uint i = 0; i < aElementListOnProc.length(); i++ )
    {
        tElementLevel( i ) = mBaseElement.give_element_level( aElementListOnProc( i ), aModelDim, aNumberOfElementsPerDirection );
        tElementId( i ) = aElementListOnProc( i );
    }

    //Create additional nodal information and save it in fields for output
    Mat<real> tNodalID( aLagrangeListOnProc.length(), 1, 0 );
    for ( uint i = 0; i < aLagrangeListOnProc.length(); i++ )
    {
        tNodalID( i ) = (real)aLagrangeListOnProc( i );
    }

    //Save field data in cells for MTK
    Cell < Mat< real > > tFieldData      = { tElementLevel,
                                             tElementId,
                                             tNodalID};

    Cell < std::string > tFieldName      = { "Element level",
                                             "Element ID ",
                                             "Lagrange ID (Internal)"};

    Cell < enum EntityRank > tFieldRanks = { EntityRank::ELEMENT,
                                             EntityRank::ELEMENT,
                                             EntityRank::NODE};
    MtkFieldsInfo aFieldsInfo;
    aFieldsInfo.FieldsData = &tFieldData;
    aFieldsInfo.FieldsName = tFieldName;
    aFieldsInfo.FieldsRank = tFieldRanks;

    //Save mesh data for MTK
    MtkMeshData aMeshData;
    aMeshData.SpatialDim = &aModelDim;
    aMeshData.ElemConn   = &aFetopo;
    aMeshData.NodeCoords = &aFeCoord;

    if ( par_size() > 1 )
    {
        aMeshData.EntProcOwner  = &aBasisProcOwner;
    }
    aMeshData.LocaltoGlobalElemMap = &aElementListOnProc;
    aMeshData.LocaltoGlobalNodeMap = &aLagrangeToBSplineMap;
    aMeshData.SetsInfo             = &aMeshSets;
    aMeshData.FieldsInfo           = &aFieldsInfo;

    // Create MORIS mesh using MTK database
    mesh tMesh( MeshType::MTK, aMeshData );
    std::string tOutputFileName = OutputFileName;
    tMesh.create_output_mesh( tOutputFileName );
}

//-------------------------------------------------------------------------------

void
Hierarchical_Mesh_MTK::create_MTK_file(
        std::string const &     OutputFileName,
        uint &                  aModelDim,
        Mat<uint> &             aNumberOfElementsPerDirection,
        Mat<uint> &             aFetopo,
        Mat<real> &             aFeCoord,
        Mat<uint> &             aBasisProcOwner,
        Mat<uint> &             aElementListOnProc,
        Mat<uint> &             aLagrangeToBSplineMap,
        Mat<uint> &             aLagrangeListOnProc,
        Cell<Mat<uint>> &       aSideSet,
        Cell<Mat<real>> &       aFieldData,
        Cell<std::string> &     aFieldName,
        Cell<enum EntityRank> & aFieldRank)
{
    // Create MORIS mesh using MTK database
    ///////////////////////////////////////
    MtkSetsInfo aMeshSets;

    // Save SideSet
    MtkSideSetsInfo tSideSetStruc;
    tSideSetStruc.ElemIdsAndSideOrds = &aSideSet;
    Cell< std::string > tSetsName( aSideSet.size() );
    //Create names for the sets
    for ( uint i = 0; i < aSideSet.size(); i++ )
    {
        tSetsName( i ) = "tSideSet_" + std::to_string(i);
    }
    tSideSetStruc.SSetNames   = tSetsName;
    if( aSideSet.size() > 0 )
    {
        aMeshSets.SideSetsInfo   = &tSideSetStruc;
    }
    //Create additional element information and save it in fields for output
    Mat<real> tElementLevel( aElementListOnProc.length(), 1, 0 );
    Mat<real> tElementId( aElementListOnProc.length(), 1, 0 );
    for ( uint i = 0; i < aElementListOnProc.length(); i++ )
    {
        tElementLevel( i ) = mBaseElement.give_element_level( aElementListOnProc( i ), aModelDim, aNumberOfElementsPerDirection );
        tElementId( i ) = aElementListOnProc( i );
    }
    //Create additional nodal information and save it in fields for output
    Mat<real> tNodalID( aLagrangeListOnProc.length(), 1, 0 );
    for ( uint i = 0; i < aLagrangeListOnProc.length(); i++ )
    {
        tNodalID( i ) = (real)aLagrangeListOnProc( i );
    }

    //Save field data in cells for MTK
    Cell < Mat< real > > tFieldData      = { tElementLevel,
                                             tElementId,
                                             tNodalID};
    Cell < std::string > tFieldName      = { "Element level",
                                             "Element ID ",
                                             "Lagrange ID (Internal)"};

    Cell < enum EntityRank > tFieldRanks = { EntityRank::ELEMENT,
                                             EntityRank::ELEMENT,
                                             EntityRank::NODE};

    //Append data from input and these field from this function
    aFieldData.append( tFieldData );
    aFieldName.append( tFieldName );
    aFieldRank.append( tFieldRanks );
    MtkFieldsInfo aFieldsInfo;
    aFieldsInfo.FieldsData = &aFieldData;
    aFieldsInfo.FieldsName = aFieldName;
    aFieldsInfo.FieldsRank = aFieldRank;
    //Create mesh data for MTK
    MtkMeshData aMeshData;
    aMeshData.SpatialDim = &aModelDim;
    aMeshData.ElemConn   = &aFetopo;
    aMeshData.NodeCoords = &aFeCoord;
    if (par_size() > 1 )
    {
        aMeshData.EntProcOwner  = &aBasisProcOwner;
    }
    aMeshData.LocaltoGlobalElemMap = &aElementListOnProc;
    aMeshData.LocaltoGlobalNodeMap = &aLagrangeToBSplineMap;

    aMeshData.SetsInfo   = &aMeshSets;
    aMeshData.FieldsInfo  = &aFieldsInfo;
    // Create MORIS mesh using MTK database
    mesh tMesh( MeshType::MTK, aMeshData );
    std::string tOutputFileName = OutputFileName;
    tMesh.create_output_mesh( tOutputFileName );
}

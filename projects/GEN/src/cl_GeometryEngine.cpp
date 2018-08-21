
#include "cl_GeometryEngine.hpp"
moris::GeometryEngine::GeometryEngine()
{

}
// ----------------------------------------------------------------------------

moris::Cell<moris::GeometryEngine::Object>
moris::GeometryEngine::is_intersected(moris::Mat<moris::real> const   & aNodeCoords,
                                      moris::Mat<moris::uint> const   & aNodetoEntityConn,
                                      moris::uint                       aCheckType)
{
    // Call the intersection_check function, (virtual function in the parent class)
    moris::Cell<moris::GeometryEngine::Object>tGeometryObjects =
            intersection_check(aNodeCoords, aNodetoEntityConn,aCheckType);
    return tGeometryObjects;
}
// ----------------------------------------------------------------------------
// Geometry Object Functions

// ----------------------------------------------------------------------------

moris::GeometryEngine::Object::Object()
{
    mParentEntityIndex = UINT_MAX;
}

// ----------------------------------------------------------------------------

moris::GeometryEngine::Object::Object(moris::uint     aParentEntityIndex)
{
    mParentEntityIndex = aParentEntityIndex;
}

// ----------------------------------------------------------------------------

moris::GeometryEngine::Object::~Object()
{

}

// ----------------------------------------------------------------------------

void
moris::GeometryEngine::Object::set_parent_entity_index(moris::uint aEntityIndex)
{
    mParentEntityIndex = aEntityIndex;
}

// ----------------------------------------------------------------------------

moris::uint
moris::GeometryEngine::Object::get_parent_entity_index()
{
    return mParentEntityIndex;
}

// ----------------------------------------------------------------------------

void
moris::GeometryEngine::Object::set_interface_location(moris::Mat<moris::real>   aInterfaceCoords,
                                                      moris::Mat<moris::uint>   aRelevantNodes)
{
//    moris::uint tNumAdd  = aInterfaceCoords.n_rows();
//    moris::uint tNumInit = mInterfaceCoords.n_rows();
//
//    mInterfaceCoords.resize(tNumInit + tNumAdd,3);
//    mRelevantNodes.resize(tNumInit + tNumAdd,{{0}});
//
//    moris::uint j = 0;
//    for(moris::uint i = tNumInit; i<tNumAdd; i++)
//    {
//
//        mRelevantNodes(i)       = aRelevantNodes.row(j);
//        mInterfaceCoords.row(i) = aInterfaceCoords.row(j);
//        j++;
//    }
}

// ----------------------------------------------------------------------------

void
moris::GeometryEngine::Object::set_interface_lcl_coord(moris::real aLclCoord)
{
    mInterfaceLclCoords.set_size(1,1);
    mInterfaceLclCoords(0,0) = aLclCoord;
}

void
moris::GeometryEngine::Object::set_interface_glb_coord(moris::Mat<moris::real>  &aGlbCoord)
{
    mInterfaceGlbCoords = aGlbCoord;
}
void
moris::GeometryEngine::Object::set_dx_dp(moris::real &adxdp)
{
    mdxdp = adxdp;
}

moris::Mat< moris::real >
moris::GeometryEngine::Object::get_interface_lcl_coord()
{
    return mInterfaceLclCoords;
}

#include "cl_GeometryEngine_Levelset.hpp"

moris::LevelsetEngine::LevelsetEngine(
        moris::tools::Function* aGeometryFunction )
{
    mLevelsetFunction = aGeometryFunction;
}

moris::LevelsetEngine::~LevelsetEngine()
{
    delete mLevelsetFunction;
}

// -------------------------------------------------------------------------

moris::Cell<moris::GeometryEngine::Object>
moris::LevelsetEngine::intersection_check(moris::Mat<moris::real>   const & aNodeCoords,
                                          moris::Mat<moris::uint>   const & aNodetoEntityConn,
                                          moris::uint                       aCheckType)
{

    //Get information for loops
    moris::uint tNumEntities = aNodetoEntityConn.n_rows(); // Number of elements in the entire mesh

    //Initialize
    moris::uint tIntersectedCount = 0;    // Intersected element counter
    moris::Cell<moris::GeometryEngine::Object> tGeoObjects(tNumEntities);

    // Compute level set values for entire mesh (populates tMeshNodeVars)
    moris::Mat< moris::real > tNodalValues = compute_levelset_value( aNodeCoords);

    //Loop over elements and determine if the element has an intersection
    for ( moris::uint i = 0; i < tNumEntities; i++ )
    {

        //Populate the intersection flag of this element with a bool
        moris::Cell<moris::Mat<moris::real>> tIntersectionInfo = compute_intersection_info(aNodetoEntityConn.row(i), tNodalValues, aNodeCoords,aCheckType);

        //moris::cout << tElementNodeVars;
        if (( tIntersectionInfo(0)(0,0) == true ))
        {
            tGeoObjects( tIntersectedCount ).set_parent_entity_index( i );
            if(aCheckType == 1)
            {
//                tGeoObjects(tIntersectedCount).set_interface_lcl_coord(tIntersectionInfo(0,1));
                tGeoObjects(tIntersectedCount).set_interface_lcl_coord(tIntersectionInfo(0)(0,1));
            }

            tIntersectedCount++;
            continue;
        }
    }

    // resize
    tGeoObjects.resize(tIntersectedCount,moris::GeometryEngine::Object());
    return tGeoObjects;

}

// -------------------------------------------------------------------------

moris::Mat<moris::real>
moris::LevelsetEngine::compute_levelset_value( moris::Mat< moris::real > const & aNodeCoords )
{

    //Initialize
    moris::uint tNumNodes = aNodeCoords.n_rows();
    moris::Mat<moris::real>  tNodeCoords( 1, 3, 0 ); // Node coordinates
    moris::Mat<moris::real>  tNodeVars(tNumNodes,2,0);
    for ( moris::uint n = 0; n < tNumNodes; n++ )
    {
        if ( tNodeVars( n, 1 ) == 1 )
        {
            // The levelset value at this node has already been computed
            continue;

        }
        else
        {
            tNodeCoords.row(0) = aNodeCoords.row(n);
            // Compute node level set value at coordinate
            tNodeVars( n, 0 ) = mLevelsetFunction->evaluate_function(tNodeCoords );
            tNodeVars( n, 1 ) = 1;
            continue;
        }
    }

    return tNodeVars.col(0);
}

// -------------------------------------------------------------------------

moris::Cell<moris::Mat<moris::real>>
moris::LevelsetEngine::compute_intersection_info( moris::Mat<moris::uint>  const & aEntityNodeInds,
                                                  moris::Mat<moris::real>  const & aNodeVars,
                                                  moris::Mat<moris::real>  const & aNodeCoords,
                                                  moris::uint                      aCheckType)
{
    moris::Cell<moris::Mat<moris::real>> tIntersectionInfo(3);

    //Initialize
    moris::uint tNodeInd = 0;
    moris::real tMax     = 0;
    moris::real tMin     = 0;

    moris::uint tNumNodes = aEntityNodeInds.n_cols();
    moris::Mat<moris::real> tEntityNodeVars(tNumNodes, 1);
    moris::Mat<moris::real> tEntityCoords(tNumNodes,3);

    //Loop through nodes and get levelset values from precomputed list
    for ( moris::uint n = 0; n < tNumNodes; n++ )
    {   //Get node id n
        tNodeInd = aEntityNodeInds( 0, n );

        //Grab levelset value from nodal levelset value list (-1 to adjust for indexing)
        tEntityNodeVars( n,0 ) = aNodeVars( tNodeInd , 0 );
    }

    //Compute the max and minimum levelset value for the element
    tMax = tEntityNodeVars.max();
    tMin = tEntityNodeVars.min();

//    If there is a sign change in element node variables return true, else return false
    moris::Mat<moris::real> tIntersection(1,2,0); // Initialize as false

    MORIS_ASSERT( tMax >= tMin, "Maximum levelset value is less than minimum levelset value check above calc" );
    MORIS_ASSERT(tMax!=0 || tMin!=0,"Zero levelset value at node!");

    if ( ( tMax > 0 ) && ( tMin < 0 ))
    {
        MORIS_ASSERT(( tMax > 0 ) && ( tMin < 0 ),"Error should not enter here unless a sign change occurs make sure above if statement didnt change");
        tIntersection(0,0) = 1;

        if(aCheckType == 1)
        {
            moris::Mat<moris::real> zero = {{0}};
            tIntersection(0,1) = moris::Interpolation::linear_interpolation_value(tEntityNodeVars,zero)(0,0);

            for(moris::uint i = 0; i<tNumNodes; i++)
            {
                tEntityCoords.row(i) = aNodeCoords.row(aEntityNodeInds(0,i));
            }

              // Global coordinate
              tIntersectionInfo(1) = moris::Interpolation::linear_interpolation_location(tEntityCoords,moris::Mat<moris::real>(1,1,tIntersection(0,1)));

              // Sensitivity
              tIntersectionInfo(2) = mLevelsetFunction->evaluate_dxdp(tIntersectionInfo(1));
        }
    }

    tIntersectionInfo(0) = tIntersection;

    return tIntersectionInfo;

}
moris::Mat<moris::real>
moris::LevelsetEngine::locate_interface(moris::Mat<moris::uint>  const & aElementNodeIds,
                                        moris::Mat<moris::real>  const & aNodeVars)
{
    moris::Mat< moris::real > dummy;
    return dummy;
}
// -------------------------------------------------------------------------

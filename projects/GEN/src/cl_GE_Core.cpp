/*
 * cl_GE_Core.cpp
 *
 *  Created on: Feb 22, 2019
 *      Author: sonne
 */

#include "cl_GE_Core.hpp"

using namespace moris;
using namespace ge;

//------------------------------------------------------------------------------

void
GE_Core::set_geometry( std::shared_ptr< Geometry > & aGeomPointer,
                       real                          aThreshold )
{
    mListOfGeoms.push_back( aGeomPointer );
    mThresholds.push_back(aThreshold);
}

//------------------------------------------------------------------------------

std::shared_ptr< Geometry >
GE_Core::get_geometry_pointer( uint aWhichGeometry )
{
    return mListOfGeoms( aWhichGeometry );
}

//------------------------------------------------------------------------------

moris::Cell< uint >
GE_Core::flag_element_list_for_refinement( moris::Cell< mtk::Cell* > & aElementList ,
                                  moris::Cell< real >       & aConstant, // constants relative to LS function
                                  uint                      aWhichGeometry )
{
    uint mNumberOfElements = aElementList.size(); //number of elements
    moris::Cell< uint > tRefineFlags(mNumberOfElements); //flag list to be returned

    for(uint i=0; i<mNumberOfElements; i++ )
    {
        moris::uint j = aElementList(i)->get_number_of_vertices(); // number of nodes in the element

        Matrix< DDRMat > mFlag(j,1);
        for( uint k=0; k<j; k++ )
        {
            mFlag(k,0) = mListOfGeoms(aWhichGeometry)->get_field_val_at_coordinate( aElementList(i)->get_vertex_pointers()(k)->get_coords() , aConstant );
        }

        if ( mFlag.max() > mThresholds(aWhichGeometry) && mFlag.min() < mThresholds(aWhichGeometry) )
        {
            tRefineFlags(i) = 1;
        }
        else
        {
            tRefineFlags(i) = 0;
        }
    }
    return tRefineFlags;
}

//------------------------------------------------------------------------------

moris::Cell< uint >
GE_Core::check_for_intersection( moris::Cell< real > &   aConstant,
                        mtk::Mesh*          &   aMeshPointer,
                        uint                    aWhichGeometry)
{
    uint mNumElements = aMeshPointer->get_num_elems();
    moris::Cell< mtk::Cell* > aEleList( mNumElements );

    for(uint i=0; i<mNumElements; i++)
    {
        // loop through all elements and build the moris::Cell< mtk::Cell* >
        aEleList(i) = & aMeshPointer->get_mtk_cell(i);
    }
    moris::Cell< uint > tRefFlags(mNumElements);

    tRefFlags = this->flag_element_list_for_refinement( aEleList,
                                                        aConstant,
                                                        aWhichGeometry);
    return tRefFlags;
}

//------------------------------------------------------------------------------

Matrix< DDRMat >
GE_Core::get_edge_normal_for_straight_edge_quad4(    uint const & aElemGlobInd,
                                            uint const & aEdgeSideOrd,
                                            mtk::Mesh* & aMeshPointer   )
{
    Matrix< IdMat > tEdgesOnElem = aMeshPointer->get_edges_connected_to_element_glob_ids( aElemGlobInd );
    Matrix< IdMat > tNodesOnElem = aMeshPointer->get_nodes_connected_to_element_glob_ids( aElemGlobInd );
    MORIS_ASSERT( tEdgesOnElem.numel() == 4, "get_edge_normal_for_straight_edge_quad4() is only valid for a 2D Quad4 element");

    Matrix< DDRMat > tNodeCoord0( 1, 2 );
    Matrix< DDRMat > tNodeCoord1( 1, 2 );
    Matrix< DDRMat > tVec( 1, 2 );
    Matrix< DDRMat > tNormal( 1, 2 );

    switch ( aEdgeSideOrd )
    {
    case( 0 )   :
    {
        tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(0) - 1 );
        tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(1) - 1 );
        tVec = tNodeCoord0 - tNodeCoord1;
        tVec(0,0) = tVec(0,0)/norm(tVec);   tVec(0,1) = tVec(0,1)/norm(tVec);
        tNormal(0,0) = tVec(0,1);   tNormal(0,1) = tVec(0,0);
        break;
    }
    case( 1 )   :
    {
        tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(1) - 1 );
        tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(2) - 1 );
        tVec = tNodeCoord1 - tNodeCoord0;
        tVec(0,0) = tVec(0,0)/norm(tVec);   tVec(0,1) = tVec(0,1)/norm(tVec);
        tNormal(0,0) = tVec(0,1);   tNormal(0,1) = tVec(0,0);
        break;
    }
    case( 2 )   :
    {
        tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(2) - 1 );
        tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(3) - 1 );
        tVec = tNodeCoord0 - tNodeCoord1;
        tVec(0,0) = tVec(0,0)/norm(tVec);   tVec(0,1) = tVec(0,1)/norm(tVec);
        tNormal(0,0) = tVec(0,1);   tNormal(0,1) = tVec(0,0);
        break;
    }
    case( 3 )   :
    {
        tNodeCoord0 = aMeshPointer->get_node_coordinate( tNodesOnElem(3) - 1 );
        tNodeCoord1 = aMeshPointer->get_node_coordinate( tNodesOnElem(0) - 1 );
        tVec = tNodeCoord1 - tNodeCoord0;
        tVec(0,0) = tVec(0,0)/norm(tVec);   tVec(0,1) = tVec(0,1)/norm(tVec);
        tNormal(0,0) = tVec(0,1);   tNormal(0,1) = tVec(0,0);
        break;
    }
    default     :
    {
        MORIS_ERROR( false, "unknown side ordinal rank");
        break;
    }
    }

    return tNormal;
}

//------------------------------------------------------------------------------

void
GE_Core::find_cells_within_levelset(
              Cell< mtk::Cell * >      & aCells,
              Cell< mtk::Cell * >      & aCandidates,
              const  Matrix< DDRMat >  & aVertexValues,
              const  uint                aUpperBound )
{


    // make sure that the field is a scalar field
    MORIS_ASSERT( aVertexValues.n_cols() == 1,
            "find_cells_within_levelset() can only be performed on scalar fields" );

    // make sure that node values are calculated
    //MORIS_ASSERT( tVertexValues.length() == aScalarField->get_num_nodes(),
    //        "number of field values does not match number of vertices on block" );

    // initialize output cell
    aCells.resize( aCandidates.size(), nullptr );

    // initialize counter
    uint tCount = 0;

    // loop over all candidates
    for( mtk::Cell * tCell : aCandidates )
    {
        // get cell of vertex pointers
        Cell< mtk::Vertex * > tVertices = tCell->get_vertex_pointers();

        // get number of vertices on this element
        uint tNumberOfVertices = tVertices.size();

        // assign matrix with vertex values
        Matrix< DDRMat > tCellValues( tNumberOfVertices, 1 );

        // loop over all vertices and extract scalar field
        for( uint k=0; k<tNumberOfVertices; ++k )
        {
            // copy value from field into element local matrix
            tCellValues( k ) = aVertexValues( tVertices( k )->get_index() );
        }

        // test if cell is inside
        if(  tCellValues.max() <= aUpperBound )
        {
            // copy pointer to output
            aCells( tCount++ ) = tCell;
        }
    }

    // shrink output to fit
    aCells.resize( tCount );
}

//------------------------------------------------------------------------------

void
GE_Core::find_cells_intersected_by_levelset(
        Cell< mtk::Cell * >         & aCells,
        Cell< mtk::Cell * >         & aCandidates,
        const  Matrix< DDRMat >     & aVertexValues,
        const  real                   aLowerBound,
        const  real                   aUpperBound )
{
    // make sure that input makes sense
    MORIS_ASSERT( aLowerBound <= aUpperBound,
            "find_cells_intersected_by_levelset() : aLowerBound bound must be less or equal aUpperBound" );

    // make sure that the field is a scalar field
    MORIS_ASSERT( aVertexValues.n_cols() == 1,
            "find_cells_within_levelset() can only be performed on scalar fields" );

    // initialize output cell
    aCells.resize( aCandidates.size(), nullptr );

    // initialize counter
    uint tCount = 0;

    // loop over all candidates
    for( mtk::Cell * tCell : aCandidates )
    {
        // get cell of vertex pointers
        Cell< mtk::Vertex * > tVertices = tCell->get_vertex_pointers();

        // get number of vertices on this element
        uint tNumberOfVertices = tVertices.size();

        // assign matrix with vertex values
        Matrix< DDRMat > tCellValues( tNumberOfVertices, 1 );

        // loop over all vertices and extract scalar field
        for( uint k=0; k<tNumberOfVertices; ++k )
        {
            // copy value from field into element local matrix
            tCellValues( k ) = aVertexValues( tVertices( k )->get_index() );
        }

        // test if cell is inside
        if ( tCellValues.min() <= aUpperBound && tCellValues.max() >= aLowerBound )
        {
            // copy pointer to output
            aCells( tCount++ ) = tCell;
        }
    }

    // shrink output to fit
    aCells.resize( tCount );
}

//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

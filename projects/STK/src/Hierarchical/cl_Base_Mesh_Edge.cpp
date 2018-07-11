/*
 * cl_Base_Mesh_Edge.cpp
 *
 *  Created on: Feb 23, 2018
 *      Author: gleim
 */

#include "cl_Base_Mesh_Edge.hpp"
using namespace moris;

uint
Base_Mesh_Edge::give_number_of_edges_x(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    MORIS_ASSERT( aModelDim > 1, "Edge function is only available for 2D and 3D");
    // output variable
    uint tEdgeNumber=0;
    if( aModelDim == 2 )
    {
        // Computes the number of edges from level zero until level "aLevel"
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            tEdgeNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) )
                                 * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) + 1 );
        }
    }
    else if( aModelDim == 3 )
    {
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            tEdgeNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) )
                                 * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) + 1 )
                                 * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) + 1 );
        }
    }
    return tEdgeNumber;
}

//--------------------------------------------------------------------------------

uint
Base_Mesh_Edge::give_number_of_edges_y(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    MORIS_ASSERT( aModelDim > 1, "Edge function is only available for 2D and 3D");
    // output variable
    uint tEdgeNumber=0;
    if( aModelDim == 2 )
    {
        // Computes the number of edges from level zero until level "aLevel"
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            tEdgeNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) + 1 )
                                 * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) );
        }
    }
    else if( aModelDim == 3 )
    {
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            tEdgeNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) + 1 )
                                 * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) )
                                 * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) + 1 );
        }
    }
    return tEdgeNumber;
}

//--------------------------------------------------------------------------------

uint
Base_Mesh_Edge::give_number_of_edges_z(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    uint tEdgeNumber=0; // output variable
    MORIS_ASSERT( aModelDim == 3, " Edges in z-direction are only in 3D available");
    for( uint i = 0; i < aLevel + 1; i++ )
    {
        // Computes the number of edges from level zero until level "aLevel"
        tEdgeNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) + 1 )
                             * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) + 1 )
                             * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) );
    }
    return tEdgeNumber;
}

//--------------------------------------------------------------------------------

uint
Base_Mesh_Edge::give_number_of_edges(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    MORIS_ASSERT( aModelDim > 1, "Edge function is only available for 2D and 3D");
    // output variable
    uint tEdgeNumber=0;
    if( aModelDim == 2 )
    {
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            // Computes the number of edges from level zero until level "aLevel"
            tEdgeNumber +=  ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) )
                                  * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) + 1 )
                                  + ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) + 1 )
                                  * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) );
        }
    }
    else if( aModelDim == 3 )
    {
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            tEdgeNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) + 1 )
                                 * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) )
                                 * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) + 1 )
                                 + ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) + 1 )
                                 * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) + 1 )
                                 * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) )
                                 + ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) )
                                 * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) + 1 )
                                 * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) + 1 );
        }
    }
    return tEdgeNumber;
}

//--------------------------------------------------------------------------------

uint
Base_Mesh_Edge::give_edge_x_of_position(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aIJKPosition) const
{
    MORIS_ASSERT( aModelDim > 1, "Edge function is only available for 2D and 3D");
    uint tEdgeId = 0;
    // Last level is needed to compute the number of edges until the last level
    uint tLevel = aLevel - 1;
    // Position of the edge in y-direction
    uint tY_Position = 0;
    if( aModelDim == 2 )
    {
        tY_Position = aIJKPosition( 1 ) * aNumberOfElementsPerDirection( 0 );
        if( aLevel > 0 )
        {
            //Give number of edges in x-direction for the last level to have the initial number at position 0,0
            tEdgeId = give_number_of_edges_x( tLevel, aModelDim, aNumberOfElementsPerDirection );
            // Position of the edge in y-direction for a specific level
            tY_Position = aIJKPosition( 1 ) * aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel );
        }
        //Edge with the calculated y-position
        tEdgeId = tEdgeId + aIJKPosition( 0 ) + tY_Position;
    }
    else if( aModelDim == 3 )
    {
        tY_Position = aIJKPosition( 1 ) * aNumberOfElementsPerDirection( 0 )
                            + aIJKPosition( 2 ) * aNumberOfElementsPerDirection( 0 ) * ( 1 + aNumberOfElementsPerDirection( 1 ) );
        if( aLevel > 0 )
        {
            //Give number of edges in x-direction for the last level to have the initial number at position 0,0,0
            tEdgeId = give_number_of_edges_x(tLevel,aModelDim,aNumberOfElementsPerDirection);
            // Position of the edge in y-direction for a specific level
            tY_Position = aIJKPosition( 1 ) * aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel )
                                + aIJKPosition( 2 ) * aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel )
                                * ( 1 + aNumberOfElementsPerDirection( 1 ) * pow( 2, aLevel ) );
        }
        //Edge with the calculated y-position
        tEdgeId = tEdgeId + aIJKPosition( 0 ) + tY_Position;
    }
    return tEdgeId;
}

//--------------------------------------------------------------------------------

uint
Base_Mesh_Edge::give_edge_y_of_position(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aIJKPosition) const
{
    MORIS_ASSERT( aModelDim > 1, "Edge function is only available for 2D and 3D");
    uint tEdgeId = 0;
    // Last level is needed to compute the number of edges until the last level
    uint tLevel = aLevel - 1;
    // Position of the edge in y-direction
    uint tY_Position = 0;
    if( aModelDim == 2 )
    {
        tY_Position = aIJKPosition( 1 ) * ( aNumberOfElementsPerDirection( 0 ) + 1 );
        if( aLevel > 0 )
        {
            //Give number of edges in x-direction for the last level to have the initial number at position 0,0
            tEdgeId = give_number_of_edges_y(tLevel,aModelDim,aNumberOfElementsPerDirection);
            // Position of the edge in y-direction for a specific level
            tY_Position = aIJKPosition( 1 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) * pow(2,aLevel) );
        }
        //Edge with the calculated y-position
        tEdgeId = tEdgeId + aIJKPosition( 0 ) + tY_Position;
    }
    else if(aModelDim == 3)
    {
        tY_Position = aIJKPosition( 1 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) )
                            + aIJKPosition( 2 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) ) * aNumberOfElementsPerDirection( 1 );
        if( aLevel>0 )
        {
            //Give number of edges in x-direction for the last level to have the initial number at position 0,0,0
            tEdgeId = give_number_of_edges_y(tLevel,aModelDim,aNumberOfElementsPerDirection);
            // Position of the edge in y-direction for a specific level
            tY_Position = aIJKPosition( 1 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) * pow(2,aLevel) )
                                + aIJKPosition( 2 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) * pow(2,aLevel) )
                                * ( aNumberOfElementsPerDirection( 1 ) * pow(2,aLevel) );
        }
        //Edge with the calculated y-position
        tEdgeId = tEdgeId + aIJKPosition( 0 ) + tY_Position;
    }
    return tEdgeId;
}

//--------------------------------------------------------------------------------

uint
Base_Mesh_Edge::give_edge_z_of_position(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aIJKPosition) const
{
    MORIS_ASSERT( aModelDim == 3, " Edges in z-direction are only in 3D available" );
    uint tEdgeId = 0;
    // Last level is needed to compute the number of edges until the last level
    uint tLevel = aLevel - 1;
    // Position of the edge in y-direction
    uint tY_Position = aIJKPosition( 1 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) )
                             + aIJKPosition( 2 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) ) * ( 1 + aNumberOfElementsPerDirection( 1 ) );
    if( aLevel > 0 )
    {
        //Give number of edges in x-direction for the last level to have the initial number at position 0,0,0
        tEdgeId = give_number_of_edges_z(tLevel,aModelDim,aNumberOfElementsPerDirection);
        // Position of the edge in y-direction for a specific level
        tY_Position = aIJKPosition( 1 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) * pow(2,aLevel) )
                            + aIJKPosition( 2 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) * pow(2,aLevel) )
                            * ( 1 + aNumberOfElementsPerDirection( 1 ) * pow(2,aLevel) );
    }
    //Edge with the calculated y-position
    tEdgeId = tEdgeId + aIJKPosition( 0 ) + tY_Position;
    return tEdgeId;
}

//--------------------------------------------------------------------------------

uint
Base_Mesh_Edge::give_edge_level(
        uint const & aEdgeId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    MORIS_ASSERT( aModelDim > 1, "Edge function is only available for 2D and 3D");
    // output variable
    uint tEdgeLevel = 1;
    // temporary variable for the while loop
    uint tLevel = 1;
    //Compute the relation of the different levels by the power of the level
    uint tPowLevel = pow(2,tEdgeLevel-1);
    if( aModelDim == 2 )
    {
        // tComputes the number of edges for level 0
        uint tNumberOfEdges = (tPowLevel * aNumberOfElementsPerDirection( 0 ) ) * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) + 1 )
                                    + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 ) * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) );
        while( tLevel > 0 )
        {
            // Checks if it fits in the level
            tLevel = floor( (real)( aEdgeId + 1 ) / tNumberOfEdges );
            if( (real)( aEdgeId + 1 ) / tNumberOfEdges == 1 || tLevel < 1 )
            {
                tLevel = 0;
            }
            else
            {
                tEdgeLevel++;
                tPowLevel = pow(2,tEdgeLevel-1);
                // Increase the number of edges of the next level
                tNumberOfEdges += ( tPowLevel * aNumberOfElementsPerDirection( 0 ) )
                                       * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) + 1 )
                                       + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 )
                                       * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) );
            }
        }
        tEdgeLevel--;
    }
    else if( aModelDim == 3 )
    {
        // temporary variable for the number of edges
        uint tNumber_of_edges = ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 )
                                      * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) ) * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) + 1 )
                                      + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 ) * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) + 1 )
                                      * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) ) + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) )
                                      * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) + 1 ) * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) + 1 );
        while( tLevel > 0 )
        {
            tLevel = floor( (real)( aEdgeId + 1 ) / tNumber_of_edges );
            if( (real)( aEdgeId + 1 ) / tNumber_of_edges == 1 || tLevel < 1 )
            {
                tLevel = 0;
            }
            else
            {
                tEdgeLevel++;
                tPowLevel = pow(2,tEdgeLevel-1);
                tNumber_of_edges += ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 )  * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) )
                                          * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) + 1 )  + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 )
                                          * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) + 1 ) * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) )
                                          + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) ) * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) + 1 )
                                          * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) + 1);
            }
        }
        tEdgeLevel--;
    }
    return tEdgeLevel;
}

//--------------------------------------------------------------------------------

Mat<uint>
Base_Mesh_Edge::give_element_edges(
        uint const & aElementId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    MORIS_ASSERT( aModelDim > 1, "Edge function is only available for 2D and 3D");
    // Compute position of element
    Mat<uint> tEdgePosition = mBaseElement.give_position_of_element( aElementId, aModelDim, aNumberOfElementsPerDirection );
    // Compute the level of the element
    uint tElementLevel = mBaseElement.give_element_level( aElementId, aModelDim, aNumberOfElementsPerDirection );
    Mat<uint> tElementEdges;
    Mat<uint> aIJKPosition( aModelDim , 1 , 0 );
    if( aModelDim == 2)
    {
        tElementEdges.set_size(4,1,UINT_MAX);
        aIJKPosition( 0 ) = tEdgePosition( 0 );
        aIJKPosition( 1 ) = tEdgePosition( 1 );
        //Compute the first edge number in x direction with the position of the element id
        uint tEdgeXNumber1 = give_edge_x_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 );
        aIJKPosition( 1 ) = tEdgePosition( 1 ) + 1;
        //Compute the second edge number in x direction with the position of the element id
        uint tEdgeXNumber2 = give_edge_x_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 );
        aIJKPosition( 1 ) = tEdgePosition( 1 );
        //Compute the first edge number in y direction with the position of the element id
        uint tEdgeYNumber1 = give_edge_y_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 ) + 1;
        aIJKPosition( 1 ) = tEdgePosition( 1 );
        //Compute the second edge number in y direction with the position of the element id
        uint tEdgeYNumber2 = give_edge_y_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        uint tEdgeXNumber = give_number_of_edges_x( tElementLevel, aModelDim, aNumberOfElementsPerDirection );
        uint tEdgeYNumber_old = 0;
        if( tElementLevel > 0)
        {
            // Give number of edges from last level
            tEdgeYNumber_old = give_number_of_edges_y( tElementLevel - 1, aModelDim, aNumberOfElementsPerDirection );
        }
        tElementEdges( 0 ) = tEdgeXNumber1 + tEdgeYNumber_old;
        tElementEdges( 1 ) = tEdgeXNumber2 + tEdgeYNumber_old;
        tElementEdges( 2 ) = tEdgeYNumber1 + tEdgeXNumber;
        tElementEdges( 3 ) = tEdgeYNumber2 + tEdgeXNumber;
    }
    else if( aModelDim == 3)
    {
        tElementEdges.set_size(12,1,UINT_MAX);
        aIJKPosition( 0 ) = tEdgePosition( 0 );
        aIJKPosition( 1 ) = tEdgePosition( 1 );
        aIJKPosition( 2 ) = tEdgePosition( 2 );
        //Compute the first edge number in x direction with the position of the element id
        uint tEdgeXNumber1 = give_edge_x_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 );
        aIJKPosition( 1 ) = tEdgePosition( 1 ) + 1;
        aIJKPosition( 2 ) = tEdgePosition( 2 );
        //Compute the second edge number in x direction with the position of the element id
        uint tEdgeXNumber2 = give_edge_x_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 );
        aIJKPosition( 1 ) = tEdgePosition( 1 );
        aIJKPosition( 2 ) = tEdgePosition( 2 ) + 1;
        //Compute the third edge number in x direction with the position of the element id
        uint tEdgeXNumber3 = give_edge_x_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 );
        aIJKPosition( 1 ) = tEdgePosition( 1 ) + 1;
        aIJKPosition( 2 ) = tEdgePosition( 2 ) + 1;
        //Compute the fourth edge number in x direction with the position of the element id
        uint tEdgeXNumber4 = give_edge_x_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 );
        aIJKPosition( 1 ) = tEdgePosition( 1 );
        aIJKPosition( 2 ) = tEdgePosition( 2 );
        //Compute the first edge number in y direction with the position of the element id
        uint tEdgeYNumber1 = give_edge_y_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 ) + 1;
        aIJKPosition( 1 ) = tEdgePosition( 1 );
        aIJKPosition( 2 ) = tEdgePosition( 2 );
        //Compute the second edge number in y direction with the position of the element id
        uint tEdgeYNumber2 = give_edge_y_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 );
        aIJKPosition( 1 ) = tEdgePosition( 1 );
        aIJKPosition( 2 ) = tEdgePosition( 2 ) + 1;
        //Compute the third edge number in y direction with the position of the element id
        uint tEdgeYNumber3 = give_edge_y_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 ) + 1;
        aIJKPosition( 1 ) = tEdgePosition( 1 );
        aIJKPosition( 2 ) = tEdgePosition( 2 ) + 1;
        //Compute the fourth edge number in y direction with the position of the element id
        uint tEdgeYNumber4 = give_edge_y_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 );
        aIJKPosition( 1 ) = tEdgePosition( 1 );
        aIJKPosition( 2 ) = tEdgePosition( 2 );
        //Compute the first edge number in z direction with the position of the element id
        uint tEdgeZNumber1 = give_edge_z_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 ) + 1;
        aIJKPosition( 1 ) = tEdgePosition( 1 );
        aIJKPosition( 2 ) = tEdgePosition( 2 );
        //Compute the second edge number in z direction with the position of the element id
        uint tEdgeZNumber2 = give_edge_z_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 );
        aIJKPosition( 1 ) = tEdgePosition( 1 ) + 1;
        aIJKPosition( 2 ) = tEdgePosition( 2 );
        //Compute the third edge number in z direction with the position of the element id
        uint tEdgeZNumber3 = give_edge_z_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tEdgePosition( 0 ) + 1;
        aIJKPosition( 1 ) = tEdgePosition( 1 ) + 1;
        aIJKPosition( 2 ) = tEdgePosition( 2 );
        //Compute the fourth edge number in z direction with the position of the element id
        uint tEdgeZNumber4 = give_edge_z_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        uint tEdgeXNumber = give_number_of_edges_x( tElementLevel, aModelDim, aNumberOfElementsPerDirection );
        uint tEdgeYNumber = give_number_of_edges_y( tElementLevel, aModelDim, aNumberOfElementsPerDirection );
        uint tEdgeYNumber_old = 0;
        uint tEdgeZNumber_old = 0;
        if( tElementLevel > 0 )
        {
            // Give number of edges from last level
            tEdgeYNumber_old = give_number_of_edges_y( tElementLevel - 1, aModelDim, aNumberOfElementsPerDirection );
            tEdgeZNumber_old = give_number_of_edges_z( tElementLevel - 1, aModelDim, aNumberOfElementsPerDirection );
        }
        tElementEdges( 0 ) = tEdgeXNumber1 + tEdgeYNumber_old + tEdgeZNumber_old;
        tElementEdges( 1 ) = tEdgeXNumber2 + tEdgeYNumber_old + tEdgeZNumber_old;
        tElementEdges( 2 ) = tEdgeXNumber3 + tEdgeYNumber_old + tEdgeZNumber_old;
        tElementEdges( 3 ) = tEdgeXNumber4 + tEdgeYNumber_old + tEdgeZNumber_old;
        tElementEdges( 4 ) = tEdgeYNumber1 + tEdgeXNumber + tEdgeZNumber_old;
        tElementEdges( 5 ) = tEdgeYNumber2 + tEdgeXNumber + tEdgeZNumber_old;
        tElementEdges( 6 ) = tEdgeYNumber3 + tEdgeXNumber + tEdgeZNumber_old;
        tElementEdges( 7 ) = tEdgeYNumber4 + tEdgeXNumber + tEdgeZNumber_old;
        tElementEdges( 8 ) = tEdgeZNumber1 + tEdgeXNumber + tEdgeYNumber;
        tElementEdges( 9 ) = tEdgeZNumber2 + tEdgeXNumber + tEdgeYNumber;
        tElementEdges( 10 ) = tEdgeZNumber3 + tEdgeXNumber + tEdgeYNumber;
        tElementEdges( 11 ) = tEdgeZNumber4 + tEdgeXNumber + tEdgeYNumber;
    }
    return tElementEdges;
}

//--------------------------------------------------------------------------------

Mat<uint>
Base_Mesh_Edge::give_elements_of_edge(
        uint const & aEdgeId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    MORIS_ASSERT( aModelDim > 1, "Edge function is only available for 2D and 3D");
    // Compute the position of the edge
    Mat<uint> tEdgePosition = give_edge_position( aEdgeId, aModelDim, aNumberOfElementsPerDirection );
    // Determine the level of the edge
    uint tEdgeLevel = give_edge_level( aEdgeId, aModelDim, aNumberOfElementsPerDirection);
    Mat<uint> tElements;
    if( aModelDim == 2 )
    {
        tElements.set_size(2,1);
        Mat<uint> tElementPositon(2,1);
        if( tEdgePosition( 0 ) < UINT_MAX && tEdgePosition( 1 ) < UINT_MAX )
        {
            tElementPositon( 0 ) = tEdgePosition( 0 );
            tElementPositon( 1 ) = tEdgePosition( 1 )-1;
            // Compute, based on the positon of the edge the positon of the element
            tElements( 0 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tEdgePosition( 0 );
            tElementPositon( 1 ) = tEdgePosition( 1 );
            tElements( 1 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
        }
        else if( tEdgePosition( 2 ) < UINT_MAX && tEdgePosition( 3 ) < UINT_MAX )
        {
            tElementPositon( 0 ) = tEdgePosition( 2 )-1;
            tElementPositon( 1 ) = tEdgePosition( 3 );
            tElements( 0 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tEdgePosition( 2 );
            tElementPositon( 1 ) = tEdgePosition( 3 );
            tElements( 1 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
        }
    }
    else if( aModelDim == 3 )
    {
        tElements.set_size(4,1);
        Mat<uint> tElementPositon(3,1);
        if( tEdgePosition( 0 ) < UINT_MAX && tEdgePosition( 1 ) < UINT_MAX && tEdgePosition( 2 ) < UINT_MAX )
        {
            tElementPositon( 0 ) = tEdgePosition( 0 );
            tElementPositon( 1 ) = tEdgePosition( 1 )-1;
            tElementPositon( 2 ) = tEdgePosition( 2 )-1;
            tElements( 0 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tEdgePosition( 0 );
            tElementPositon( 1 ) = tEdgePosition( 1 );
            tElementPositon( 2 ) = tEdgePosition( 2 )-1;
            tElements( 1 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tEdgePosition( 0 );
            tElementPositon( 1 ) = tEdgePosition( 1 )-1;
            tElementPositon( 2 ) = tEdgePosition( 2 );
            tElements( 2 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tEdgePosition( 0 );
            tElementPositon( 1 ) = tEdgePosition( 1 );
            tElementPositon( 2 ) = tEdgePosition( 2 );
            tElements( 3 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
        }
        else if( tEdgePosition( 3 ) < UINT_MAX && tEdgePosition( 4 ) < UINT_MAX && tEdgePosition( 5 ) < UINT_MAX )
        {
            tElementPositon( 0 ) = tEdgePosition( 3 )-1;
            tElementPositon( 1 ) = tEdgePosition( 4 );
            tElementPositon( 2 ) = tEdgePosition( 5 )-1;
            tElements( 0 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tEdgePosition( 3 );
            tElementPositon( 1 ) = tEdgePosition( 4 );
            tElementPositon( 2 ) = tEdgePosition( 5 )-1;
            tElements( 1 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tEdgePosition( 3 )-1;
            tElementPositon( 1 ) = tEdgePosition( 4 );
            tElementPositon( 2 ) = tEdgePosition( 5 );
            tElements( 2 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tEdgePosition( 3 );
            tElementPositon( 1 ) = tEdgePosition( 4 );
            tElementPositon( 2 ) = tEdgePosition( 5 );
            tElements( 3 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
        }
        else if( tEdgePosition( 6 ) < UINT_MAX && tEdgePosition( 7 ) < UINT_MAX && tEdgePosition( 8 ) < UINT_MAX )
        {
            tElementPositon( 0 ) = tEdgePosition( 6 )-1;
            tElementPositon( 1 ) = tEdgePosition( 7 )-1;
            tElementPositon( 2 ) = tEdgePosition( 8 );
            tElements( 0 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tEdgePosition( 6 );
            tElementPositon( 1 ) = tEdgePosition( 7 )-1;
            tElementPositon( 2 ) = tEdgePosition( 8 );
            tElements( 1 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tEdgePosition( 6 )-1;
            tElementPositon( 1 ) = tEdgePosition( 7 );
            tElementPositon( 2 ) = tEdgePosition( 8 );
            tElements( 2 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tEdgePosition( 6 );
            tElementPositon( 1 ) = tEdgePosition( 7 );
            tElementPositon( 2 ) = tEdgePosition( 8 );
            tElements( 3 ) = mBaseElement.give_element_of_position( tEdgeLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
        }
    }
    return tElements;
}

//--------------------------------------------------------------------------------

Mat<uint>
Base_Mesh_Edge::give_edge_position(
        uint const & aEdgeId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    MORIS_ASSERT( aModelDim > 1, "Edge function is only available for 2D and 3D");
    // Compute the level of the edge
    uint tEdgeLevel = give_edge_level( aEdgeId, aModelDim, aNumberOfElementsPerDirection );
    //Compute the relation of the different levels by the power of the level
    uint tPowLevel = pow( 2, tEdgeLevel );
    Mat<uint> tEdgePosition;
    if( aModelDim == 2 )
    {
        tEdgePosition.set_size( 4, 1, UINT_MAX );
        uint tEdgeNumber_old = 0;
        if( tEdgeLevel > 0 )
        {
            // Number of edges until the last smaller level
            tEdgeNumber_old = give_number_of_edges( tEdgeLevel - 1, aModelDim, aNumberOfElementsPerDirection );
        }
        // Number of edges until the last smaller level plus number of edges in x-direction of the current level
        uint tEdgeNumbery = tEdgeNumber_old + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) )
                                  * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) + 1 );
        if( aEdgeId >= tEdgeNumber_old && aEdgeId < tEdgeNumbery )
        {
            // Calculate the y-position
            tEdgePosition( 1 ) = ceil( (real)( aEdgeId + 1 - tEdgeNumber_old )
                    / (real)( aNumberOfElementsPerDirection( 0 ) * tPowLevel ) );
            // Calculate the x-direction with the calculated y-position
            tEdgePosition( 0 ) = aEdgeId+1 - tEdgeNumber_old - ( tEdgePosition( 1 ) - 1 )
                                       * ( aNumberOfElementsPerDirection( 0 ) * tPowLevel );
            // -1 to get an indexed basis position, starting with zero
            tEdgePosition( 1 )--;
            tEdgePosition( 0 )--;
        }
        if( aEdgeId >= tEdgeNumbery)
        {
            tEdgePosition( 3 ) = ceil( (real)( aEdgeId + 1 - tEdgeNumbery )
                    / (real)( aNumberOfElementsPerDirection( 0 ) * tPowLevel+1 ) );
            tEdgePosition( 2 ) = aEdgeId+1 - tEdgeNumbery - ( tEdgePosition( 3 ) - 1 )
                                       * ( aNumberOfElementsPerDirection( 0 ) * tPowLevel + 1 );
            // -1 to get an indexed basis position, starting with zero
            tEdgePosition( 3 )--;
            tEdgePosition( 2 )--;
        }
    }
    else if( aModelDim == 3 )
    {
        tEdgePosition.set_size( 9 , 1 , UINT_MAX );
        uint tEdgeNumber_old = 0;
        if( tEdgeLevel > 0 )
        {
            // Number of edges until the last smaller level
            tEdgeNumber_old = give_number_of_edges(tEdgeLevel-1,aModelDim,aNumberOfElementsPerDirection);
        }
        // Number of edges until the last smaller level plus number of edges in x-direction of the current level
        uint tEdgeNumbery = tEdgeNumber_old + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) )
                                  * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) + 1 ) * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) + 1);
        // Number of edges until the last smaller level plus number of edges in x- and y-direction of the current level
        uint tEdgeNumberz = tEdgeNumber_old + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) )
                                  * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) + 1 ) * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) + 1 )
                                  + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 ) * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) )
                                  * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) + 1 );
        if( aEdgeId >= tEdgeNumber_old && aEdgeId < tEdgeNumbery )
        {
            // Calculate the z-position
            tEdgePosition( 2 ) = ceil( (real)( aEdgeId + 1 - tEdgeNumber_old )
                    / (real)( ( aNumberOfElementsPerDirection( 0 ) * tPowLevel ) * ( 1 + aNumberOfElementsPerDirection( 1 ) * tPowLevel ) ) );
            // Calculate the y-position
            tEdgePosition( 1 ) = ceil( (real)( aEdgeId + 1 - tEdgeNumber_old - ( tEdgePosition( 2 ) - 1 ) * ( (aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                    * (real)( 1 + aNumberOfElementsPerDirection( 1 ) * tPowLevel ) ) ) / (real)( aNumberOfElementsPerDirection( 0 ) * tPowLevel ) );
            // Calculate the x-direction with the calculated y-position
            tEdgePosition( 0 ) = aEdgeId + 1 - tEdgeNumber_old - ( tEdgePosition( 1 ) - 1 ) * ( aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                                       - ( tEdgePosition( 2 ) - 1 ) * ( ( aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                                               * ( 1 + aNumberOfElementsPerDirection( 1 ) * tPowLevel ) );
            // -1 to get an indexed basis position, starting with zero
            tEdgePosition( 2 )--;
            tEdgePosition( 1 )--;
            tEdgePosition( 0 )--;
        }
        else if( aEdgeId >= tEdgeNumbery && aEdgeId < tEdgeNumberz )
        {
            // Calculate the z-position
            tEdgePosition( 5 ) = ceil( (real)( aEdgeId + 1 - tEdgeNumbery ) / (real)( ( 1 + aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                    * ( aNumberOfElementsPerDirection( 1 ) * tPowLevel ) ) );
            // Calculate the y-position
            tEdgePosition( 4 ) = ceil( (real)( aEdgeId + 1 - tEdgeNumbery - ( tEdgePosition( 5 ) - 1 )
                    * ( ( 1 + aNumberOfElementsPerDirection( 0 ) * tPowLevel ) * ( aNumberOfElementsPerDirection( 1 ) * tPowLevel ) ) )
                    / (real)( 1 + aNumberOfElementsPerDirection( 0 ) * tPowLevel ) );
            // Calculate the x-direction with the calculated y-position
            tEdgePosition( 3 ) = aEdgeId + 1 - tEdgeNumbery - ( tEdgePosition( 4 ) - 1 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                                       - ( tEdgePosition( 5 ) - 1 ) * ( ( 1 + aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                                               * ( aNumberOfElementsPerDirection( 1 ) * tPowLevel ) );
            // -1 to get an indexed basis position, starting with zero
            tEdgePosition( 5 )--;
            tEdgePosition( 4 )--;
            tEdgePosition( 3 )--;
        }
        else if( aEdgeId >= tEdgeNumberz )
        {
            // Calculate the z-position
            tEdgePosition( 8 ) = ceil( (real)( aEdgeId + 1 - tEdgeNumberz ) / (real)( ( 1 + aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                    * ( 1 + aNumberOfElementsPerDirection( 1 ) * tPowLevel ) ));
            // Calculate the y-position
            tEdgePosition( 7 ) = ceil( (real)( aEdgeId + 1 - tEdgeNumberz - ( tEdgePosition( 8 ) - 1 ) * ( ( 1 + aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                    * ( 1 + aNumberOfElementsPerDirection( 1 ) * tPowLevel ) ) )
                    / (real)( 1 + aNumberOfElementsPerDirection( 0 ) * tPowLevel ) );
            // Calculate the x-direction with the calculated y-position
            tEdgePosition( 6 ) = aEdgeId+1 - tEdgeNumberz - ( tEdgePosition( 7 ) - 1 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                                       - ( tEdgePosition( 8 ) - 1 ) * ( ( 1+ aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                                               * ( 1 + aNumberOfElementsPerDirection( 1 ) * tPowLevel ) );
            // -1 to get an indexed basis position, starting with zero
            tEdgePosition( 8 )--;
            tEdgePosition( 7 )--;
            tEdgePosition( 6 )--;
        }
    }
    return tEdgePosition;
}

//--------------------------------------------------------------------------------

uint
Base_Mesh_Edge::give_edge_owner(
        uint const & aEdgeId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aProcRange,
        Mat<uint> const & aProcNeighbors) const
{
    MORIS_ASSERT( aModelDim > 1, "Edge function is only available for 2D and 3D");
    uint tProcRank = par_rank();
    // Compute the position of the edge
    Mat<uint> tEdgePosition = give_edge_position( aEdgeId, aModelDim, aNumberOfElementsPerDirection );
    // Compute the level of the edge
    uint tEdgeLevel = give_edge_level( aEdgeId, aModelDim, aNumberOfElementsPerDirection );
    uint tProcOwner = UINT_MAX;
    if( par_size() > 1 )
    {
        // Compute the range on the respective level
        Mat<uint> tProcRange( aProcRange.length() , 1 , 0 );
        tProcRange( 0 ) = aProcRange( 0 ) * pow( 2, tEdgeLevel );
        tProcRange( 1 ) = aProcRange( 1 ) * pow( 2, tEdgeLevel );
        tProcRange( 2 ) = aProcRange( 2 ) * pow( 2, tEdgeLevel );
        tProcRange( 3 ) = aProcRange( 3 ) * pow( 2, tEdgeLevel );
        if( aModelDim == 2 )
        {
            // Check if it is an edge in x-direction
            if( tEdgePosition( 0 ) < UINT_MAX && tEdgePosition( 1 ) < UINT_MAX )
            {
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 )
                        && tEdgePosition( 1 ) > tProcRange( 2 ) && tEdgePosition( 1 ) <= tProcRange( 3 ) ) // Inner edges
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 )
                        && tEdgePosition( 1 ) == tProcRange( 2 ) && aProcNeighbors( 0 ) == UINT_MAX ) // Bottom proc edge
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 )
                        && tEdgePosition( 1 ) == tProcRange( 2 ) &&  aProcNeighbors( 0 ) < UINT_MAX ) // Bottom proc edge
                    tProcOwner = aProcNeighbors( 0 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && aProcNeighbors( 1 ) == UINT_MAX ) // Bottom right proc edge
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) > tProcRange( 2 )
                        && tEdgePosition( 1 ) <= tProcRange( 3 ) && aProcNeighbors( 1 ) < UINT_MAX ) // Right outer edge
                    tProcOwner = aProcNeighbors( 1 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && aProcNeighbors( 15 ) < UINT_MAX ) // Bottom right outer edge
                    tProcOwner = aProcNeighbors( 15 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && aProcNeighbors( 15 ) == UINT_MAX && aProcNeighbors( 1 ) < UINT_MAX ) // Bottom right outer edge
                    tProcOwner = aProcNeighbors( 1 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && aProcNeighbors( 15 ) == UINT_MAX && aProcNeighbors( 1 ) == UINT_MAX ) // Bottom right outer edge
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 )
                        && tEdgePosition( 1 ) == tProcRange( 3 )+1 &&  aProcNeighbors( 2 ) < UINT_MAX ) // Top outer edges
                    tProcOwner = aProcNeighbors( 2 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )+1
                        &&  aProcNeighbors( 17 ) < UINT_MAX ) // Top right outer edge
                    tProcOwner = aProcNeighbors( 17 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 )
                        && tEdgePosition( 1 ) == tProcRange( 2 )-1 &&  aProcNeighbors( 0 ) < UINT_MAX ) // Top outer edges
                    tProcOwner = aProcNeighbors( 0 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )-1
                        &&  aProcNeighbors( 15 ) < UINT_MAX ) // Top right outer edge
                    tProcOwner = aProcNeighbors( 15 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) > tProcRange( 2 )
                        && tEdgePosition( 1 ) <= tProcRange( 3 ) && aProcNeighbors( 3 ) < UINT_MAX ) // Left outer edges
                    tProcOwner = aProcNeighbors( 3 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 )
                        && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer edge
                    tProcOwner = aProcNeighbors( 8 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 )
                        && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 3 ) < UINT_MAX ) // Bottom left outer edge
                    tProcOwner = aProcNeighbors( 3 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 3 )+1
                        &&  aProcNeighbors( 16 ) < UINT_MAX ) // Top left outer edge
                    tProcOwner = aProcNeighbors( 16 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 )-1
                        &&  aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer edge
                    tProcOwner = aProcNeighbors( 8 );
            }
            else if( tEdgePosition( 2 ) < UINT_MAX && tEdgePosition( 3 ) < UINT_MAX )
            {
                // Check if it is an edge in y-direction
                if( tEdgePosition( 2 ) > tProcRange( 0 ) && tEdgePosition( 2 ) <= tProcRange( 1 )
                        && tEdgePosition( 3 ) >= tProcRange( 2 ) && tEdgePosition( 3 ) < tProcRange( 3 ) ) // Inner edges
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 2 ) == tProcRange( 0 ) && tEdgePosition( 3 ) >= tProcRange( 2 )
                        && tEdgePosition( 3 ) < tProcRange( 3 ) && aProcNeighbors( 3 ) == UINT_MAX ) // Left proc edge
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 2 ) == tProcRange( 0 ) && tEdgePosition( 3 ) >= tProcRange( 2 )
                        && tEdgePosition( 3 ) < tProcRange( 3 ) && aProcNeighbors( 3 ) < UINT_MAX )  // Left proc edge
                    tProcOwner = aProcNeighbors( 3 );
                else if( tEdgePosition( 2 ) == tProcRange( 1 )+1 && tEdgePosition( 3 ) >= tProcRange( 2 )
                        && tEdgePosition( 3 ) < tProcRange( 3 ) && aProcNeighbors( 1 ) < UINT_MAX )  // Right outer edges
                    tProcOwner = aProcNeighbors( 1 );
                else if( tEdgePosition( 2 ) == tProcRange( 1 )+1 && tEdgePosition( 3 ) >= tProcRange( 2 )
                        && tEdgePosition( 3 ) < tProcRange( 3 ) && aProcNeighbors( 1 ) == UINT_MAX ) // Right outer edges
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 2 ) > tProcRange( 0 ) && tEdgePosition( 2 ) <= tProcRange( 1 )
                        && tEdgePosition( 3 ) == tProcRange( 3 ) && aProcNeighbors( 2 ) < UINT_MAX ) // Top outer edges
                    tProcOwner = aProcNeighbors( 2 );
                else if( tEdgePosition( 2 ) > tProcRange( 0 ) && tEdgePosition( 2 ) <= tProcRange( 1 )
                        && tEdgePosition( 3 ) == tProcRange( 3 ) && aProcNeighbors( 2 ) == UINT_MAX ) // Top outer edges
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 2 ) == tProcRange( 1 )+1 && tEdgePosition( 3 ) == tProcRange( 3 )
                        && aProcNeighbors( 17 ) < UINT_MAX ) // Top right outer edge
                    tProcOwner = aProcNeighbors( 17 );
                else if( tEdgePosition( 2 ) == tProcRange( 0 ) && tEdgePosition( 3 ) == tProcRange( 3 )
                        && aProcNeighbors( 16 ) < UINT_MAX ) // Top outer left edge
                    tProcOwner = aProcNeighbors( 16 );
                else if( tEdgePosition( 2 ) == tProcRange( 0 ) && tEdgePosition( 3 ) == tProcRange( 3 )
                        && aProcNeighbors( 16 ) == UINT_MAX && aProcNeighbors( 2 ) < UINT_MAX ) // Top outer left edge
                    tProcOwner = aProcNeighbors( 2 );
                else if( tEdgePosition( 2 ) > tProcRange( 0 ) && tEdgePosition( 2 ) <= tProcRange( 1 )
                        && tEdgePosition( 3 ) == tProcRange( 2 )-1 && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom outer edges
                    tProcOwner = aProcNeighbors( 0 );
                else if( tEdgePosition( 2 ) == tProcRange( 0 ) && tEdgePosition( 3 ) == tProcRange( 2 )-1
                        && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom outer edges
                    tProcOwner = aProcNeighbors( 8 );
                else if( tEdgePosition( 2 ) == tProcRange( 0 ) && tEdgePosition( 3 ) == tProcRange( 2 )-1
                        && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom outer edges
                    tProcOwner = aProcNeighbors( 0 );
                else if( tEdgePosition( 2 ) == tProcRange( 1 )+1 && tEdgePosition( 3 ) == tProcRange( 2 )-1
                        && aProcNeighbors( 15 ) < UINT_MAX ) // Bottom  right outer edge
                    tProcOwner = aProcNeighbors( 15 );
                else if( tEdgePosition( 2 ) == tProcRange( 0 )-1 && tEdgePosition( 3 ) >= tProcRange( 2 )
                        && tEdgePosition( 3 ) < tProcRange( 3 ) && aProcNeighbors( 3 ) < UINT_MAX )  // Left outer edges
                    tProcOwner = aProcNeighbors( 3 );
                else if( tEdgePosition( 2 ) == tProcRange( 0 )-1 && tEdgePosition( 3 ) == tProcRange( 3 )
                        && aProcNeighbors( 16 ) < UINT_MAX )  // Top left outer edge
                    tProcOwner = aProcNeighbors( 16 );
                else if( tEdgePosition( 2 ) == tProcRange( 0 )-1 && tEdgePosition( 3 ) == tProcRange( 2 )-1
                        && aProcNeighbors( 8 ) < UINT_MAX )  // Bottom left outer edge
                    tProcOwner = aProcNeighbors( 8 );
            }
        }
        else if( aModelDim == 3)
        {
            tProcRange( 4 ) = aProcRange( 4 ) * pow( 2, tEdgeLevel );
            tProcRange( 5 ) = aProcRange( 5 ) * pow( 2, tEdgeLevel );
            // Check if it is an edge in x-direction
            if( tEdgePosition( 0 ) < UINT_MAX && tEdgePosition( 1 ) < UINT_MAX  && tEdgePosition( 2 ) < UINT_MAX )
            {
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) > tProcRange( 2 )
                        && tEdgePosition( 1 ) <= tProcRange( 3 ) && tEdgePosition( 2 ) > tProcRange( 4 ) && tEdgePosition( 2 ) <= tProcRange( 5 ) ) // Inner edges
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && tEdgePosition( 2 ) > tProcRange( 4 ) && tEdgePosition( 2 ) <= tProcRange( 5 ) && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom edges
                    tProcOwner = aProcNeighbors( 0 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && tEdgePosition( 2 ) > tProcRange( 4 ) && tEdgePosition( 2 ) <= tProcRange( 5 ) && aProcNeighbors( 0 ) == UINT_MAX ) // Bottom edges
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) > tProcRange( 2 )
                        && tEdgePosition( 1 ) <= tProcRange( 3 ) && tEdgePosition( 2 ) == tProcRange( 4 ) && aProcNeighbors( 4 ) < UINT_MAX ) // Back edges
                    tProcOwner = aProcNeighbors( 4 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) > tProcRange( 2 )
                        && tEdgePosition( 1 ) <= tProcRange( 3 ) && tEdgePosition( 2 ) == tProcRange( 4 ) && aProcNeighbors( 4 ) == UINT_MAX ) // Back edges
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && tEdgePosition( 2 ) == tProcRange( 4 ) &&   aProcNeighbors( 9 ) < UINT_MAX ) // Back bottom edges
                    tProcOwner = aProcNeighbors( 9 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 )  && tEdgePosition( 1 ) == tProcRange( 2 )
                        && tEdgePosition( 2 ) == tProcRange( 4 ) &&   aProcNeighbors( 9 ) == UINT_MAX && aProcNeighbors( 0 ) < UINT_MAX )  // Back bottom edges
                    tProcOwner = aProcNeighbors( 0 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 )  && tEdgePosition( 1 ) == tProcRange( 2 )
                        && tEdgePosition( 2 ) == tProcRange( 4 ) &&   aProcNeighbors( 9 ) == UINT_MAX && aProcNeighbors( 4 ) < UINT_MAX )  // Back bottom edges
                    tProcOwner = aProcNeighbors( 4 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 )  && tEdgePosition( 1 ) == tProcRange( 2 )
                        && tEdgePosition( 2 ) == tProcRange( 4 ) && aProcNeighbors( 9 ) == UINT_MAX && aProcNeighbors( 4 ) == UINT_MAX
                        && aProcNeighbors( 0 ) == UINT_MAX )  // Back bottom edges
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) > tProcRange( 2 ) && tEdgePosition( 1 ) <= tProcRange( 3 )
                        && tEdgePosition( 2 ) > tProcRange( 4 ) && tEdgePosition( 2 ) <= tProcRange( 5 ) && aProcNeighbors( 1 ) < UINT_MAX ) // Right outer edges
                    tProcOwner = aProcNeighbors( 1 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) > tProcRange( 4 )
                        && tEdgePosition( 2 ) <= tProcRange( 5 ) && aProcNeighbors( 15 ) < UINT_MAX ) // Right bottom outer edges
                    tProcOwner = aProcNeighbors( 15 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) > tProcRange( 4 )
                        && tEdgePosition( 2 ) <= tProcRange( 5 ) && aProcNeighbors( 15 ) == UINT_MAX && aProcNeighbors( 1 ) < UINT_MAX ) // Right bottom outer edges
                    tProcOwner = aProcNeighbors( 1 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) > tProcRange( 2 ) && tEdgePosition( 1 ) <= tProcRange( 3 )
                        && tEdgePosition( 2 ) > tProcRange( 4 ) && tEdgePosition( 2 ) <= tProcRange( 5 ) && aProcNeighbors( 3 ) < UINT_MAX ) // Left outer edges
                    tProcOwner = aProcNeighbors( 3 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) > tProcRange( 4 )
                        && tEdgePosition( 2 ) <= tProcRange( 5 ) && aProcNeighbors( 8 ) < UINT_MAX ) // Left bottom outer edges
                    tProcOwner = aProcNeighbors( 8 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) > tProcRange( 4 )
                        && tEdgePosition( 2 ) <= tProcRange( 5 ) && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 3 ) < UINT_MAX ) // Left bottom outer edges
                    tProcOwner = aProcNeighbors( 3 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) > tProcRange( 2 ) && tEdgePosition( 1 ) <= tProcRange( 3 )
                        && tEdgePosition( 2 ) == tProcRange( 4 ) && aProcNeighbors( 11 ) < UINT_MAX ) // Back Right outer edges
                    tProcOwner = aProcNeighbors( 11 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) > tProcRange( 2 ) && tEdgePosition( 1 ) <= tProcRange( 3 )
                        && tEdgePosition( 2 ) == tProcRange( 4 ) && aProcNeighbors( 11 ) == UINT_MAX && aProcNeighbors( 1 ) < UINT_MAX ) // Back Right outer edges
                    tProcOwner = aProcNeighbors( 1 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) > tProcRange( 2 ) && tEdgePosition( 1 ) <= tProcRange( 3 )
                        && tEdgePosition( 2 ) == tProcRange( 4 ) && aProcNeighbors( 7 ) < UINT_MAX ) // Back left outer edges
                    tProcOwner = aProcNeighbors( 7 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) > tProcRange( 2 ) && tEdgePosition( 1 ) <= tProcRange( 3 )
                        && tEdgePosition( 2 ) == tProcRange( 4 ) && aProcNeighbors( 7 ) == UINT_MAX && aProcNeighbors( 3 ) < UINT_MAX ) // Back left outer edges
                    tProcOwner = aProcNeighbors( 3 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 10 ) < UINT_MAX ) // Back bottom right outer edges
                    tProcOwner = aProcNeighbors( 10 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 10 ) == UINT_MAX && aProcNeighbors( 15 ) < UINT_MAX ) // Back bottom right outer edges
                    tProcOwner = aProcNeighbors( 15 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 10 ) == UINT_MAX && aProcNeighbors( 15 ) == UINT_MAX && aProcNeighbors( 11 ) < UINT_MAX ) // Back bottom right outer edges
                    tProcOwner = aProcNeighbors( 11 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 10 ) == UINT_MAX && aProcNeighbors( 15 ) == UINT_MAX && aProcNeighbors( 11 ) == UINT_MAX
                        && aProcNeighbors( 1 ) < UINT_MAX ) // Back bottom right outer edges
                    tProcOwner = aProcNeighbors( 1 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 6 ) < UINT_MAX ) // Back bottom left outer edges
                    tProcOwner = aProcNeighbors( 6 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 8 ) < UINT_MAX ) // Back bottom left outer edges
                    tProcOwner = aProcNeighbors( 8 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 7 ) < UINT_MAX ) // Back bottom left outer edges
                    tProcOwner = aProcNeighbors( 7 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 7 ) == UINT_MAX
                        && aProcNeighbors( 3 ) < UINT_MAX ) // Back bottom left outer edges
                    tProcOwner = aProcNeighbors( 3 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )+1
                        && tEdgePosition( 2 ) > tProcRange( 4 ) && tEdgePosition( 2 ) <= tProcRange( 5 )  && aProcNeighbors( 2 ) < UINT_MAX ) // Top outer edges
                    tProcOwner = aProcNeighbors( 2 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )+1
                        && tEdgePosition( 2 ) == tProcRange( 4 ) && aProcNeighbors( 13 ) < UINT_MAX ) // Top back outer edges
                    tProcOwner = aProcNeighbors( 13 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )+1
                        && tEdgePosition( 2 ) == tProcRange( 4 ) && aProcNeighbors( 13 ) == UINT_MAX && aProcNeighbors( 2 ) < UINT_MAX ) // Top back outer edges
                    tProcOwner = aProcNeighbors( 2 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )-1
                        && tEdgePosition( 2 ) > tProcRange( 4 ) && tEdgePosition( 2 ) <= tProcRange( 5 ) && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom outer edges
                    tProcOwner = aProcNeighbors( 0 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )-1
                        && tEdgePosition( 2 ) == tProcRange( 4 ) && aProcNeighbors( 9 ) < UINT_MAX ) // Bottom back outer edges
                    tProcOwner = aProcNeighbors( 9 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )-1
                        && tEdgePosition( 2 ) == tProcRange( 4 ) && aProcNeighbors( 9 ) == UINT_MAX && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom back outer edges
                    tProcOwner = aProcNeighbors( 0 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) > tProcRange( 2 )
                        && tEdgePosition( 1 ) <= tProcRange( 3 ) && tEdgePosition( 2 ) == tProcRange( 5 )+1 && aProcNeighbors( 5 ) < UINT_MAX ) // Front outer edges
                    tProcOwner = aProcNeighbors( 5 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && tEdgePosition( 2 ) == tProcRange( 5 )+1 && aProcNeighbors( 19 ) < UINT_MAX ) // Front outer edges
                    tProcOwner = aProcNeighbors( 19 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && tEdgePosition( 2 ) == tProcRange( 5 )+1 && aProcNeighbors( 19 ) == UINT_MAX && aProcNeighbors( 5 ) < UINT_MAX ) // Front outer edges
                    tProcOwner = aProcNeighbors( 5 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) > tProcRange( 2 )
                        && tEdgePosition( 1 ) <= tProcRange( 3 ) && tEdgePosition( 2 ) == tProcRange( 4 )-1 && aProcNeighbors( 4 ) < UINT_MAX ) // Back outer edges
                    tProcOwner = aProcNeighbors( 4 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && tEdgePosition( 2 ) == tProcRange( 4 )-1 && aProcNeighbors( 9 ) < UINT_MAX ) // Back outer edges
                    tProcOwner = aProcNeighbors( 9 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && tEdgePosition( 2 ) == tProcRange( 4 )-1 && aProcNeighbors( 9 ) == UINT_MAX && aProcNeighbors( 4 ) < UINT_MAX ) // Back outer edges
                    tProcOwner = aProcNeighbors( 4 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )+1
                        && tEdgePosition( 2 ) == tProcRange( 5 )+1 && aProcNeighbors( 24 ) < UINT_MAX ) // Front top outer edges
                    tProcOwner = aProcNeighbors( 24 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )+1
                        && tEdgePosition( 2 ) == tProcRange( 4 )-1 && aProcNeighbors( 13 ) < UINT_MAX ) // Back top outer edges
                    tProcOwner = aProcNeighbors( 13 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )-1
                        && tEdgePosition( 2 ) == tProcRange( 5 )+1 && aProcNeighbors( 19 ) < UINT_MAX ) // Front bottom outer edges
                    tProcOwner = aProcNeighbors( 19 );
                else if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) < tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )-1
                        && tEdgePosition( 2 ) == tProcRange( 4 )-1 && aProcNeighbors( 9 ) < UINT_MAX ) // Back bottom outer edges
                    tProcOwner = aProcNeighbors( 9 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) > tProcRange( 2 ) && tEdgePosition( 1 ) <= tProcRange( 3 )
                        && tEdgePosition( 2 ) == tProcRange( 5 )+1 && aProcNeighbors( 22 ) < UINT_MAX ) // Front right outer edges
                    tProcOwner = aProcNeighbors( 22 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 20 ) < UINT_MAX ) // Front right bottom outer edge
                    tProcOwner = aProcNeighbors( 20 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 20 ) == UINT_MAX && aProcNeighbors( 22 ) < UINT_MAX ) // Front right bottom outer edge
                    tProcOwner = aProcNeighbors( 22 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) > tProcRange( 2 ) && tEdgePosition( 1 ) <= tProcRange( 3 )
                        && tEdgePosition( 2 ) == tProcRange( 5 )+1 && aProcNeighbors( 21 ) < UINT_MAX ) // Front left outer edges
                    tProcOwner = aProcNeighbors( 21 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 18 ) < UINT_MAX ) // Front left bottom outer edge
                    tProcOwner = aProcNeighbors( 18 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 18 ) == UINT_MAX && aProcNeighbors( 21 ) < UINT_MAX ) // Front left bottom outer edge
                    tProcOwner = aProcNeighbors( 21 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) > tProcRange( 2 ) && tEdgePosition( 1 ) <= tProcRange( 3 )
                        && tEdgePosition( 2 ) == tProcRange( 4 )-1 && aProcNeighbors( 11 ) < UINT_MAX ) // Back right outer edges
                    tProcOwner = aProcNeighbors( 11 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 10 ) < UINT_MAX ) // Back right bottom outer edge
                    tProcOwner = aProcNeighbors( 10 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 10 ) == UINT_MAX && aProcNeighbors( 11 ) < UINT_MAX ) // Back right bottom outer edge
                    tProcOwner = aProcNeighbors( 11 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) > tProcRange( 2 ) && tEdgePosition( 1 ) <= tProcRange( 3 )
                        && tEdgePosition( 2 ) == tProcRange( 4 )-1 && aProcNeighbors( 7 ) < UINT_MAX ) // Back left outer edges
                    tProcOwner = aProcNeighbors( 7 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) < UINT_MAX ) // Back left bottom outer edge
                    tProcOwner = aProcNeighbors( 6 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 ) && tEdgePosition( 2 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 7 ) < UINT_MAX ) // Back left bottom outer edge
                    tProcOwner = aProcNeighbors( 7 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )+1 && tEdgePosition( 2 ) > tProcRange( 4 )
                        && tEdgePosition( 2 ) <= tProcRange( 5 ) && aProcNeighbors( 17 ) < UINT_MAX ) // Top right outer edges
                    tProcOwner = aProcNeighbors( 17 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )+1 && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 14 ) < UINT_MAX ) // Top right back outer edge
                    tProcOwner = aProcNeighbors( 14 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )+1 && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 14 ) == UINT_MAX && aProcNeighbors( 17 ) < UINT_MAX ) // Top right back outer edge
                    tProcOwner = aProcNeighbors( 17 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )+1 && tEdgePosition( 2 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 14 ) < UINT_MAX ) // Back top right  outer edge
                    tProcOwner = aProcNeighbors( 14 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )+1 && tEdgePosition( 2 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 25 ) < UINT_MAX ) // Back top right outer edge
                    tProcOwner = aProcNeighbors( 25 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 3 )+1 && tEdgePosition( 2 ) > tProcRange( 4 )
                        && tEdgePosition( 2 ) <= tProcRange( 5 ) && aProcNeighbors( 16 ) < UINT_MAX ) // Top left outer edges
                    tProcOwner = aProcNeighbors( 16 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 3 )+1 && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 12 ) < UINT_MAX ) // Top left back outer edge
                    tProcOwner = aProcNeighbors( 12 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 3 )+1 && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 12 ) == UINT_MAX && aProcNeighbors( 16 ) < UINT_MAX ) // Top left back outer edge
                    tProcOwner = aProcNeighbors( 16 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 3 )+1 && tEdgePosition( 2 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 12 ) < UINT_MAX ) // Back top left outer edge
                    tProcOwner = aProcNeighbors( 12 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 3 )+1 && tEdgePosition( 2 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 23 ) < UINT_MAX ) // Back top left outer edge
                    tProcOwner = aProcNeighbors( 23 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )-1 && tEdgePosition( 2 ) > tProcRange( 4 )
                        && tEdgePosition( 2 ) <= tProcRange( 5 ) && aProcNeighbors( 15 ) < UINT_MAX ) // Bottom right outer edges
                    tProcOwner = aProcNeighbors( 15 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )-1 && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 10 ) < UINT_MAX ) // Bottom right back outer edge
                    tProcOwner = aProcNeighbors( 10 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )-1 && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 10 ) == UINT_MAX && aProcNeighbors( 15 ) < UINT_MAX ) // Bottom right back outer edge
                    tProcOwner = aProcNeighbors( 15 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )-1 && tEdgePosition( 2 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 10 ) < UINT_MAX ) // Back bottom right  outer edge
                    tProcOwner = aProcNeighbors( 10 );
                else if( tEdgePosition( 0 ) == tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )-1 && tEdgePosition( 2 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 20 ) < UINT_MAX ) // Back bottom right outer edge
                    tProcOwner = aProcNeighbors( 20 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 )-1 && tEdgePosition( 2 ) > tProcRange( 4 )
                        && tEdgePosition( 2 ) <= tProcRange( 5 ) && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer edges
                    tProcOwner = aProcNeighbors( 8 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 )-1 && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 6 ) < UINT_MAX ) // Bottom left back outer edge
                    tProcOwner = aProcNeighbors( 6 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 )-1 && tEdgePosition( 2 ) == tProcRange( 4 )
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left back outer edge
                    tProcOwner = aProcNeighbors( 8 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 )-1 && tEdgePosition( 2 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) < UINT_MAX ) // Back bottom left outer edge
                    tProcOwner = aProcNeighbors( 6 );
                else if( tEdgePosition( 0 ) == tProcRange( 0 )-1 && tEdgePosition( 1 ) == tProcRange( 2 )-1 && tEdgePosition( 2 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 18 ) < UINT_MAX ) // Front bottom left outer edge
                    tProcOwner = aProcNeighbors( 18 );
            }
            else if( tEdgePosition( 3 ) < UINT_MAX && tEdgePosition( 4 ) < UINT_MAX && tEdgePosition( 5 ) < UINT_MAX  )
            {
                // Check if it is an edge in y-direction
                if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) >= tProcRange( 2 )
                        && tEdgePosition( 4 ) < tProcRange( 3 ) && tEdgePosition( 5 ) > tProcRange( 4 ) && tEdgePosition( 5 ) <= tProcRange( 5 ) ) // Inner edges
                    tProcOwner =tProcRank;
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) > tProcRange( 4 ) && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 3 ) == UINT_MAX ) // Left edges
                    tProcOwner =tProcRank;
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) > tProcRange( 4 ) && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 3 ) < UINT_MAX ) // Left edges
                    tProcOwner =aProcNeighbors( 3 );
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) >= tProcRange( 2 )
                        && tEdgePosition( 4 ) < tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )  && aProcNeighbors( 4 ) == UINT_MAX ) // Back edges
                    tProcOwner =tProcRank;
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) >= tProcRange( 2 )
                        && tEdgePosition( 4 ) < tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 ) && aProcNeighbors( 4 ) < UINT_MAX ) // Back  edges
                    tProcOwner =aProcNeighbors( 4 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 ) && aProcNeighbors( 7 ) < UINT_MAX ) // Left back edges
                    tProcOwner = aProcNeighbors( 7 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 ) && aProcNeighbors( 7 ) == UINT_MAX && aProcNeighbors( 4 ) < UINT_MAX ) // Left back edges
                    tProcOwner = aProcNeighbors( 4 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 ) && aProcNeighbors( 7 ) == UINT_MAX && aProcNeighbors( 3 ) < UINT_MAX ) // Left back edges
                    tProcOwner = aProcNeighbors( 3 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 ) && aProcNeighbors( 7 ) == UINT_MAX && aProcNeighbors( 4 ) == UINT_MAX && aProcNeighbors( 3 ) == UINT_MAX) // Left back edges
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) > tProcRange( 4 ) && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 1 ) < UINT_MAX ) // Right outer edges
                    tProcOwner =aProcNeighbors( 1 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 ) && aProcNeighbors( 11 ) < UINT_MAX ) // Right outer edges
                    tProcOwner =aProcNeighbors( 11 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 ) && aProcNeighbors( 11 ) == UINT_MAX && aProcNeighbors( 1 ) < UINT_MAX ) // Right outer edges
                    tProcOwner =aProcNeighbors( 1 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) > tProcRange( 4 ) && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 3 ) < UINT_MAX ) // Left outer edges
                    tProcOwner =aProcNeighbors( 3 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 ) && aProcNeighbors( 7 ) < UINT_MAX ) // Left outer edges
                    tProcOwner =aProcNeighbors( 7 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 ) && aProcNeighbors( 7 ) == UINT_MAX && aProcNeighbors( 3 ) < UINT_MAX ) // Left outer edges
                    tProcOwner =aProcNeighbors( 3 );
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) == tProcRange( 3 )
                        && tEdgePosition( 5 ) > tProcRange( 4 ) && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 2 ) < UINT_MAX ) // Top outer edges
                    tProcOwner =aProcNeighbors( 2 );
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) == tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 ) && aProcNeighbors( 13 ) < UINT_MAX ) // Top back outer edges
                    tProcOwner =aProcNeighbors( 13 );
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) == tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 ) && aProcNeighbors( 13 ) == UINT_MAX && aProcNeighbors( 2 ) < UINT_MAX ) // Top back outer edges
                    tProcOwner =aProcNeighbors( 2 );
                else if( tEdgePosition( 3 ) >= tProcRange( 0 )-1 && tEdgePosition( 3 ) <= tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 3 )
                        && tEdgePosition( 5 ) > tProcRange( 4 ) && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 16 ) < UINT_MAX ) // Top left outer edges
                    tProcOwner =aProcNeighbors( 16 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) > tProcRange( 4 )
                        && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 16 ) == UINT_MAX && aProcNeighbors( 2 ) < UINT_MAX ) // Top left outer edges
                    tProcOwner =aProcNeighbors( 2 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 12 ) < UINT_MAX ) // Top left outer edges
                    tProcOwner =aProcNeighbors( 12 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 12 ) == UINT_MAX && aProcNeighbors( 13 ) < UINT_MAX ) // Top left outer edges
                    tProcOwner =aProcNeighbors( 13 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 12 ) == UINT_MAX && aProcNeighbors( 13 ) == UINT_MAX && aProcNeighbors( 16 ) < UINT_MAX ) // Top left outer edges
                    tProcOwner =aProcNeighbors( 16 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 12 ) == UINT_MAX && aProcNeighbors( 13 ) == UINT_MAX && aProcNeighbors( 16 ) == UINT_MAX && aProcNeighbors( 2 ) < UINT_MAX ) // Top left outer edges
                    tProcOwner =aProcNeighbors( 2 );
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1
                        && tEdgePosition( 5 ) > tProcRange( 4 ) && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom outer edges
                    tProcOwner =aProcNeighbors( 0 );
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1
                        && tEdgePosition( 5 ) == tProcRange( 4 ) && aProcNeighbors( 9 ) < UINT_MAX ) // Bottom back outer edges
                    tProcOwner =aProcNeighbors( 9 );
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1
                        && tEdgePosition( 5 ) == tProcRange( 4 ) && aProcNeighbors( 9 ) == UINT_MAX && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom back outer edges
                    tProcOwner =aProcNeighbors( 0 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) > tProcRange( 4 )
                        && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer edges
                    tProcOwner =aProcNeighbors( 8 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) > tProcRange( 4 )
                        && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom left outer edges
                    tProcOwner =aProcNeighbors( 0 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 6 ) < UINT_MAX ) // Bottom left outer edges
                    tProcOwner =aProcNeighbors( 6 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 9 ) < UINT_MAX ) // Bottom left outer edges
                    tProcOwner =aProcNeighbors( 9 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer edges
                    tProcOwner =aProcNeighbors( 8 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 9 ) == UINT_MAX && aProcNeighbors( 8 ) == UINT_MAX
                        && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom left outer edges
                    tProcOwner =aProcNeighbors( 0 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) > tProcRange( 4 )
                        && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer edges
                    tProcOwner =aProcNeighbors( 8 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 6 ) < UINT_MAX ) // Bottom left outer edges
                    tProcOwner =aProcNeighbors( 6 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer edges
                    tProcOwner =aProcNeighbors( 8 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) > tProcRange( 4 )
                        && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 2 ) < UINT_MAX && aProcNeighbors( 17 ) < UINT_MAX ) // Top right outer edges
                    tProcOwner =aProcNeighbors( 17 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 2 ) < UINT_MAX && aProcNeighbors( 14 ) < UINT_MAX ) // Top right back outer edges
                    tProcOwner =aProcNeighbors( 14 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 2 ) < UINT_MAX && aProcNeighbors( 14 ) == UINT_MAX && aProcNeighbors( 17 ) < UINT_MAX ) // Top right back outer edges
                    tProcOwner =aProcNeighbors( 17 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) > tProcRange( 4 )
                        && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 2 ) < UINT_MAX && aProcNeighbors( 16 ) < UINT_MAX ) // Top left outer edges
                    tProcOwner =aProcNeighbors( 16 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 2 ) < UINT_MAX && aProcNeighbors( 12 ) < UINT_MAX ) // Top left back outer edges
                    tProcOwner =aProcNeighbors( 12 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 2 ) < UINT_MAX && aProcNeighbors( 12 ) == UINT_MAX && aProcNeighbors( 16 ) < UINT_MAX ) // Top left back outer edges
                    tProcOwner =aProcNeighbors( 16 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) > tProcRange( 4 )
                        && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 15 ) < UINT_MAX ) // Bottom right outer edges
                    tProcOwner =aProcNeighbors( 15 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 10 ) < UINT_MAX ) // Bottom right back outer edges
                    tProcOwner =aProcNeighbors( 10 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 10 ) == UINT_MAX && aProcNeighbors( 15 ) < UINT_MAX ) // Bottom right back outer edges
                    tProcOwner =aProcNeighbors( 15 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) > tProcRange( 4 )
                        && tEdgePosition( 5 ) <= tProcRange( 5 ) && aProcNeighbors( 2 ) < UINT_MAX && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer edges
                    tProcOwner =aProcNeighbors( 8 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 2 ) < UINT_MAX && aProcNeighbors( 6 ) < UINT_MAX ) // Bottom left back outer edges
                    tProcOwner =aProcNeighbors( 6 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )
                        && aProcNeighbors( 2 ) < UINT_MAX && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left back outer edges
                    tProcOwner =aProcNeighbors( 8 );
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) >= tProcRange( 2 )
                        && tEdgePosition( 4 ) < tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 5 )+1 && aProcNeighbors( 5 ) < UINT_MAX ) // Front edges
                    tProcOwner =aProcNeighbors( 5 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 5 )+1 && aProcNeighbors( 21 ) < UINT_MAX ) // Front left edges
                    tProcOwner =aProcNeighbors( 21 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 5 )+1 && aProcNeighbors( 21 ) == UINT_MAX && aProcNeighbors( 5 ) < UINT_MAX ) // Front left edges
                    tProcOwner =aProcNeighbors( 5 );
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) == tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 5 )+1 && aProcNeighbors( 24 ) < UINT_MAX ) // Front top edges
                    tProcOwner =aProcNeighbors( 24 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 16 ) < UINT_MAX ) // Front top left edges
                    tProcOwner =aProcNeighbors( 16 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 16 ) == UINT_MAX && aProcNeighbors( 23 ) < UINT_MAX ) // Front top left edges
                    tProcOwner =aProcNeighbors( 23 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 16 ) == UINT_MAX && aProcNeighbors( 23 ) == UINT_MAX && aProcNeighbors( 24 ) < UINT_MAX ) // Front top left edges
                    tProcOwner =aProcNeighbors( 24 );
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1
                        && tEdgePosition( 5 ) == tProcRange( 5 )+1 && aProcNeighbors( 19 ) < UINT_MAX ) // Front bottom edges
                    tProcOwner =aProcNeighbors( 19 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 18 ) < UINT_MAX ) // Front bottom left edges
                    tProcOwner =aProcNeighbors( 18 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 18 ) == UINT_MAX && aProcNeighbors( 19 ) < UINT_MAX ) // Front bottom left edges
                    tProcOwner =aProcNeighbors( 19 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 5 )+1 && aProcNeighbors( 21 ) < UINT_MAX ) // Front left outer edge
                    tProcOwner =aProcNeighbors( 21 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 5 )+1 && aProcNeighbors( 22 ) < UINT_MAX ) // Front right outer edge
                    tProcOwner =aProcNeighbors( 22 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 23 ) < UINT_MAX ) // Front top left outer edge
                    tProcOwner =aProcNeighbors( 23 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 25 ) < UINT_MAX ) // Front top right outer edge
                    tProcOwner =aProcNeighbors( 25 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 18 ) < UINT_MAX ) // Front bottom left outer edge
                    tProcOwner =aProcNeighbors( 18 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 5 )+1
                        && aProcNeighbors( 20 ) < UINT_MAX ) // Front bottom right outer edge
                    tProcOwner =aProcNeighbors( 20 );
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) >= tProcRange( 2 )
                        && tEdgePosition( 4 ) < tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )-1 && aProcNeighbors( 4 ) < UINT_MAX ) // Back edges
                    tProcOwner =aProcNeighbors( 4 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 )-1 && aProcNeighbors( 7 ) < UINT_MAX ) // Back left edges
                    tProcOwner =aProcNeighbors( 7 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 )-1 && aProcNeighbors( 7 ) == UINT_MAX && aProcNeighbors( 4 ) < UINT_MAX ) // Back left edges
                    tProcOwner =aProcNeighbors( 4 );
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) == tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 )-1 && aProcNeighbors( 13 ) < UINT_MAX ) // Back top edges
                    tProcOwner =aProcNeighbors( 13 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 12 ) < UINT_MAX ) // Back top left edge
                    tProcOwner =aProcNeighbors( 12 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 12 ) == UINT_MAX && aProcNeighbors( 13 ) < UINT_MAX ) // Back top left edge
                    tProcOwner =aProcNeighbors( 13 );
                else if( tEdgePosition( 3 ) > tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1
                        && tEdgePosition( 5 ) == tProcRange( 4 )-1 && aProcNeighbors( 9 ) < UINT_MAX ) // Back bottom edges
                    tProcOwner =aProcNeighbors( 9 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) < UINT_MAX ) // Back bottom left edges
                    tProcOwner =aProcNeighbors( 6 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 9 ) < UINT_MAX ) // Back bottom left edges
                    tProcOwner =aProcNeighbors( 9 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 )-1 && aProcNeighbors( 7 ) < UINT_MAX ) // Back left outer edge
                    tProcOwner =aProcNeighbors( 7 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) < tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 )-1 && aProcNeighbors( 11 ) < UINT_MAX ) // Back right outer edge
                    tProcOwner =aProcNeighbors( 11 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 12 ) < UINT_MAX ) // Back top left outer edge
                    tProcOwner =aProcNeighbors( 12 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) == tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 14 ) < UINT_MAX ) // Back top right outer edge
                    tProcOwner =aProcNeighbors( 14 );
                else if( tEdgePosition( 3 ) == tProcRange( 0 )-1 && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) < UINT_MAX ) // Back bottom left outer edge
                    tProcOwner =aProcNeighbors( 6 );
                else if( tEdgePosition( 3 ) == tProcRange( 1 )+1 && tEdgePosition( 4 ) == tProcRange( 2 )-1 && tEdgePosition( 5 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 10 ) < UINT_MAX ) // Back bottom right outer edge
                    tProcOwner =aProcNeighbors( 10 );
            }
            else if( tEdgePosition( 6 ) < UINT_MAX && tEdgePosition( 7 ) < UINT_MAX && tEdgePosition( 8 ) < UINT_MAX  )
            {
                // Check if it is an edge in z-direction
                if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) > tProcRange( 2 )
                        && tEdgePosition( 7 ) <= tProcRange( 3 ) && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) < tProcRange( 5 ) ) // Inner edges
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 2 )
                        && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 0 ) == UINT_MAX ) // Bottom edges
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 2 )
                        && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom edges
                    tProcOwner = aProcNeighbors( 0 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) > tProcRange( 2 ) && tEdgePosition( 7 ) <= tProcRange( 3 )
                        && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 3 ) < UINT_MAX ) // Left edges
                    tProcOwner = aProcNeighbors( 3 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) > tProcRange( 2 ) && tEdgePosition( 7 ) <= tProcRange( 3 )
                        && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 3 ) == UINT_MAX ) // Left edges
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left edges
                    tProcOwner = aProcNeighbors( 8 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom left edges
                    tProcOwner = aProcNeighbors( 0 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 3 ) < UINT_MAX ) // Bottom left edges
                    tProcOwner = aProcNeighbors( 3 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 0 ) == UINT_MAX && aProcNeighbors( 3 ) == UINT_MAX ) // Bottom left edges
                    tProcOwner = tProcRank;
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 2 )-1
                        && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom outer edges
                    tProcOwner = aProcNeighbors( 0 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 )-1 && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer edges
                    tProcOwner = aProcNeighbors( 8 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 )-1 && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom left outer edges
                    tProcOwner = aProcNeighbors( 0 );
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 3 )+1
                        && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 2 ) < UINT_MAX ) // Top outer edges
                    tProcOwner = aProcNeighbors( 2 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 3 )+1 && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 16 ) < UINT_MAX ) // Top left outer edges
                    tProcOwner = aProcNeighbors( 16 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 3 )+1 && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 16 ) == UINT_MAX && aProcNeighbors( 2 ) < UINT_MAX ) // Top left outer edges
                    tProcOwner = aProcNeighbors( 2 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) > tProcRange( 2 ) && tEdgePosition( 7 ) <= tProcRange( 3 )
                        && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 3 ) < UINT_MAX ) // Left outer edges
                    tProcOwner = aProcNeighbors( 3 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 8 ) < UINT_MAX ) // Left outer edges
                    tProcOwner = aProcNeighbors( 8 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 3 ) < UINT_MAX ) // Left outer edges
                    tProcOwner = aProcNeighbors( 3 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) > tProcRange( 2 ) && tEdgePosition( 7 ) <= tProcRange( 3 )
                        && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 1 ) < UINT_MAX ) // Right outer edges
                    tProcOwner = aProcNeighbors( 1 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 15 ) < UINT_MAX ) // Right outer edges
                    tProcOwner = aProcNeighbors( 15 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 15 ) == UINT_MAX && aProcNeighbors( 1 ) < UINT_MAX ) // Right outer edges
                    tProcOwner = aProcNeighbors( 1 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) == tProcRange( 2 )-1 && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 8 ) < UINT_MAX ) // Left bottom outer edges
                    tProcOwner = aProcNeighbors( 8 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) == tProcRange( 2 )-1 && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 15 ) < UINT_MAX ) // Right bottom outer edges
                    tProcOwner = aProcNeighbors( 15 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) == tProcRange( 3 )+1 && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 16 ) < UINT_MAX ) // Left top outer edges
                    tProcOwner = aProcNeighbors( 16 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) == tProcRange( 3 )+1 && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) < tProcRange( 5 ) && aProcNeighbors( 17 ) < UINT_MAX ) // Right top outer edges
                    tProcOwner = aProcNeighbors( 17 );
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) > tProcRange( 2 )
                        && tEdgePosition( 7 ) <= tProcRange( 3 ) && tEdgePosition( 8 ) == tProcRange( 5 ) && aProcNeighbors( 5 ) < UINT_MAX ) // Front edges
                    tProcOwner = aProcNeighbors( 5 );
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 2 )
                        && tEdgePosition( 8 ) == tProcRange( 5 ) && aProcNeighbors( 19 ) < UINT_MAX ) // Front bottom edges
                    tProcOwner = aProcNeighbors( 19 );
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 2 )
                        && tEdgePosition( 8 ) == tProcRange( 5 ) && aProcNeighbors( 19 ) == UINT_MAX && aProcNeighbors( 5 ) < UINT_MAX ) // Front bottom edges
                    tProcOwner = aProcNeighbors( 5 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) > tProcRange( 2 ) && tEdgePosition( 7 ) <= tProcRange( 3 )
                        && tEdgePosition( 8 ) == tProcRange( 5 ) && aProcNeighbors( 21 ) < UINT_MAX ) // Front left edges
                    tProcOwner = aProcNeighbors( 21 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) > tProcRange( 2 ) && tEdgePosition( 7 ) <= tProcRange( 3 )
                        && tEdgePosition( 8 ) == tProcRange( 5 ) && aProcNeighbors( 21 ) == UINT_MAX && aProcNeighbors( 5 ) < UINT_MAX ) // Front left edges
                    tProcOwner = aProcNeighbors( 5 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 18 ) < UINT_MAX ) // Front bottom left edges
                    tProcOwner = aProcNeighbors( 18 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 18 ) == UINT_MAX && aProcNeighbors( 19 ) < UINT_MAX ) // Front bottom left edges
                    tProcOwner = aProcNeighbors( 19 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 18 ) == UINT_MAX && aProcNeighbors( 21 ) < UINT_MAX ) // Front bottom left edges
                    tProcOwner = aProcNeighbors( 21 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 18 ) == UINT_MAX && aProcNeighbors( 19 ) == UINT_MAX && aProcNeighbors( 21 ) == UINT_MAX ) // Front bottom left edges
                    tProcOwner = aProcNeighbors( 5 );
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 3 )+1
                        && tEdgePosition( 8 ) == tProcRange( 5 ) && aProcNeighbors( 24 ) < UINT_MAX ) // Front top edges
                    tProcOwner = aProcNeighbors( 24 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 3 )+1 && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 23 ) < UINT_MAX ) // Front top left edge
                    tProcOwner = aProcNeighbors( 23 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 3 )+1 && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 23 ) == UINT_MAX && aProcNeighbors( 24 ) < UINT_MAX ) // Front top left edge
                    tProcOwner = aProcNeighbors( 24 );
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 2 )-1
                        && tEdgePosition( 8 ) == tProcRange( 5 ) && aProcNeighbors( 19 ) < UINT_MAX ) // Front bottom edges
                    tProcOwner = aProcNeighbors( 19 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 )-1 && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 18 ) < UINT_MAX ) // Front bottom left edge
                    tProcOwner = aProcNeighbors( 18 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 )-1 && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 18 ) == UINT_MAX && aProcNeighbors( 19 ) < UINT_MAX ) // Front bottom left edge
                    tProcOwner = aProcNeighbors( 19 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) == tProcRange( 2 )-1 && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 18 ) < UINT_MAX ) // Front bottom left outer edge
                    tProcOwner = aProcNeighbors( 18 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) == tProcRange( 2 )-1 && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 20 ) < UINT_MAX ) // Front bottom right outer edge
                    tProcOwner = aProcNeighbors( 20 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) == tProcRange( 3 )+1 && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 23 ) < UINT_MAX ) // Front top left outer edge
                    tProcOwner = aProcNeighbors( 23 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) == tProcRange( 3 )+1 && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 25 ) < UINT_MAX ) // Front top right outer edge
                    tProcOwner = aProcNeighbors( 25 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) > tProcRange( 2 ) && tEdgePosition( 7 ) <= tProcRange( 3 )
                        && tEdgePosition( 8 ) == tProcRange( 5 ) && aProcNeighbors( 22 ) < UINT_MAX ) // Front right outer edges
                    tProcOwner = aProcNeighbors( 22 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 20 ) < UINT_MAX ) // Front right outer edges
                    tProcOwner = aProcNeighbors( 20 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 20 ) == UINT_MAX && aProcNeighbors( 22 ) < UINT_MAX ) // Front right outer edges
                    tProcOwner = aProcNeighbors( 22 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) > tProcRange( 2 ) && tEdgePosition( 7 ) <= tProcRange( 3 )
                        && tEdgePosition( 8 ) == tProcRange( 5 ) && aProcNeighbors( 21 ) < UINT_MAX ) // Front left outer edges
                    tProcOwner = aProcNeighbors( 21 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 18 ) < UINT_MAX ) // Front left outer edges
                    tProcOwner = aProcNeighbors( 18 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 5 )
                        && aProcNeighbors( 18 ) == UINT_MAX && aProcNeighbors( 21 ) < UINT_MAX ) // Front left outer edges
                    tProcOwner = aProcNeighbors( 21 );
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) > tProcRange( 2 )
                        && tEdgePosition( 7 ) <= tProcRange( 3 ) && tEdgePosition( 8 ) == tProcRange( 4 )-1 && aProcNeighbors( 4 ) < UINT_MAX ) // Back edges
                    tProcOwner = aProcNeighbors( 4 );
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 2 )
                        && tEdgePosition( 8 ) == tProcRange( 4 )-1 && aProcNeighbors( 9 ) < UINT_MAX ) // Back bottom edges
                    tProcOwner = aProcNeighbors( 9 );
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 2 )
                        && tEdgePosition( 8 ) == tProcRange( 4 )-1 && aProcNeighbors( 9 ) == UINT_MAX && aProcNeighbors( 4 ) < UINT_MAX ) // Back bottom edges
                    tProcOwner = aProcNeighbors( 4 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) > tProcRange( 2 ) && tEdgePosition( 7 ) <= tProcRange( 3 )
                        && tEdgePosition( 8 ) == tProcRange( 4 )-1 && aProcNeighbors( 7 ) < UINT_MAX ) // Back left edges
                    tProcOwner = aProcNeighbors( 7 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) > tProcRange( 2 ) && tEdgePosition( 7 ) <= tProcRange( 3 )
                        && tEdgePosition( 8 ) == tProcRange( 4 )-1 && aProcNeighbors( 7 ) == UINT_MAX && aProcNeighbors( 4 ) < UINT_MAX ) // Back left edges
                    tProcOwner = aProcNeighbors( 4 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) < UINT_MAX ) // Back bottom left edges
                    tProcOwner = aProcNeighbors( 6 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 9 ) < UINT_MAX ) // Back bottom left edges
                    tProcOwner = aProcNeighbors( 9 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 7 ) < UINT_MAX ) // Back bottom left edges
                    tProcOwner = aProcNeighbors( 7 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 9 ) == UINT_MAX && aProcNeighbors( 7 ) == UINT_MAX ) // Back bottom left edges
                    tProcOwner = aProcNeighbors( 4 );
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 3 )+1
                        && tEdgePosition( 8 ) == tProcRange( 4 )-1 && aProcNeighbors( 13 ) < UINT_MAX ) // Back top edges
                    tProcOwner = aProcNeighbors( 13 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 3 )+1 && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 12 ) < UINT_MAX ) // Back top left edge
                    tProcOwner = aProcNeighbors( 12 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 3 )+1 && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 12 ) == UINT_MAX && aProcNeighbors( 13 ) < UINT_MAX ) // Back top left edge
                    tProcOwner = aProcNeighbors( 13 );
                else if( tEdgePosition( 6 ) > tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 2 )-1
                        && tEdgePosition( 8 ) == tProcRange( 4 )-1 && aProcNeighbors( 9 ) < UINT_MAX ) // Back bottom edges
                    tProcOwner = aProcNeighbors( 9 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 )-1 && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) < UINT_MAX ) // Back bottom left edge
                    tProcOwner = aProcNeighbors( 6 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 7 ) == tProcRange( 2 )-1 && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 9 ) < UINT_MAX ) // Back bottom left edge
                    tProcOwner = aProcNeighbors( 9 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) == tProcRange( 2 )-1 && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) < UINT_MAX ) // Back bottom left outer edge
                    tProcOwner = aProcNeighbors( 6 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) == tProcRange( 2 )-1 && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 10 ) < UINT_MAX ) // Back bottom right outer edge
                    tProcOwner = aProcNeighbors( 10 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) == tProcRange( 3 )+1 && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 12 ) < UINT_MAX ) // Back top left outer edge
                    tProcOwner = aProcNeighbors( 12 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) == tProcRange( 3 )+1 && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 14 ) < UINT_MAX ) // Back top right outer edge
                    tProcOwner = aProcNeighbors( 14 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) > tProcRange( 2 ) && tEdgePosition( 7 ) <= tProcRange( 3 )
                        && tEdgePosition( 8 ) == tProcRange( 4 )-1 && aProcNeighbors( 11 ) < UINT_MAX ) // Back right outer edges
                    tProcOwner = aProcNeighbors( 11 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 10 ) < UINT_MAX ) // Back right outer edges
                    tProcOwner = aProcNeighbors( 10 );
                else if( tEdgePosition( 6 ) == tProcRange( 1 )+1 && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 20 ) == UINT_MAX && aProcNeighbors( 11 ) < UINT_MAX ) // Back right outer edges
                    tProcOwner = aProcNeighbors( 11 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) > tProcRange( 2 ) && tEdgePosition( 7 ) <= tProcRange( 3 )
                        && tEdgePosition( 8 ) == tProcRange( 4 )-1 && aProcNeighbors( 7 ) < UINT_MAX ) // Back left outer edges
                    tProcOwner = aProcNeighbors( 7 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) < UINT_MAX ) // Back left outer edges
                    tProcOwner = aProcNeighbors( 6 );
                else if( tEdgePosition( 6 ) == tProcRange( 0 )-1 && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) == tProcRange( 4 )-1
                        && aProcNeighbors( 6 ) == UINT_MAX && aProcNeighbors( 7 ) < UINT_MAX ) // Back left outer edges
                    tProcOwner = aProcNeighbors( 7 );
            }
        }
    }
    return tProcOwner;
}

//--------------------------------------------------------------------------------

Mat<uint>
Base_Mesh_Edge::give_edge_share(
        uint const & aEdgeId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aProcRange,
        Mat<uint> const & aProcNeighbors) const
{
    MORIS_ASSERT( aModelDim > 1, "Edge function is only available for 2D and 3D");
    uint tProcRank = par_rank();
    // Compute positon of the edge
    Mat<uint> tEdgePosition = give_edge_position( aEdgeId, aModelDim, aNumberOfElementsPerDirection );
    // Compute the level of the edge
    uint tEdgeLevel = give_edge_level( aEdgeId, aModelDim, aNumberOfElementsPerDirection );
    // Arbitrary number for possible sharing elements ( It will find with thie algorithm muliple times the same neighbors and unique will kick them out)
    Mat<uint> tProcShare( 40 , 1 , UINT_MAX );
    uint tVar = 0;
    if( par_size() > 1 )
    {
        // Compute the range on the respective level
        Mat<uint> tProcRange( aProcRange.length() , 1 , 0 );
        tProcRange( 0 ) = aProcRange( 0 ) * pow( 2, tEdgeLevel );
        tProcRange( 1 ) = aProcRange( 1 ) * pow( 2, tEdgeLevel );
        tProcRange( 2 ) = aProcRange( 2 ) * pow( 2, tEdgeLevel );
        tProcRange( 3 ) = aProcRange( 3 ) * pow( 2, tEdgeLevel );
        if( aModelDim == 2 )
        {
            // Check if it is an edge in x-direction
            if( tEdgePosition( 0 ) < UINT_MAX && tEdgePosition( 1 ) < UINT_MAX )
            {
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) <= tProcRange( 1 )
                        && tEdgePosition( 1 ) >= tProcRange( 2 ) && tEdgePosition( 1 ) <= tProcRange( 3 ) )
                {
                    tProcShare(tVar) =tProcRank;
                    tVar++;
                }
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) <= tProcRange( 1 )
                        && tEdgePosition( 1 ) == tProcRange( 2 ))
                {
                    tProcShare(tVar) =aProcNeighbors( 0 );
                    tVar++;
                }
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) <= tProcRange( 1 )
                        && tEdgePosition( 1 ) == tProcRange( 3 ))
                {
                    tProcShare(tVar) =aProcNeighbors( 2 );
                    tVar++;
                }
            }
            else if( tEdgePosition( 2 ) < UINT_MAX && tEdgePosition( 3 ) < UINT_MAX )
            {
                // Check if it is an edge in y-direction
                if( tEdgePosition( 2 ) >= tProcRange( 0 ) && tEdgePosition( 2 ) <= tProcRange( 1 )
                        && tEdgePosition( 3 ) >= tProcRange( 2 ) && tEdgePosition( 3 ) <= tProcRange( 3 ) )
                {
                    tProcShare(tVar) =tProcRank;
                    tVar++;
                }
                if( tEdgePosition( 2 ) == tProcRange( 0 )  && tEdgePosition( 3 ) >= tProcRange( 2 )
                        && tEdgePosition( 3 ) <= tProcRange( 3 ) )
                {
                    tProcShare(tVar) =aProcNeighbors( 3 );
                    tVar++;
                }
                if( tEdgePosition( 2 ) == tProcRange( 1 )  && tEdgePosition( 3 ) >= tProcRange( 2 )
                        && tEdgePosition( 3 ) <= tProcRange( 3 ) )
                {
                    tProcShare(tVar) =aProcNeighbors( 1 );
                    tVar++;
                }
            }
        }
        else if( aModelDim == 3)
        {
            // Compute the range on the respective level
            tProcRange( 4 ) = aProcRange( 4 ) * pow( 2, tEdgeLevel );
            tProcRange( 5 ) = aProcRange( 5 ) * pow( 2, tEdgeLevel );
            // Check if it is an edge in x-direction
            if( tEdgePosition( 0 ) < UINT_MAX && tEdgePosition( 1 ) < UINT_MAX  && tEdgePosition( 2 ) < UINT_MAX )
            {
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) <= tProcRange( 1 ) && tEdgePosition( 1 ) >= tProcRange( 2 )
                        && tEdgePosition( 1 ) <= tProcRange( 3 ) && tEdgePosition( 2 ) >= tProcRange( 4 ) && tEdgePosition( 2 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = tProcRank;
                    tVar++;
                }
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) <= tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && tEdgePosition( 2 ) >= tProcRange( 4 ) && tEdgePosition( 2 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 0 );
                    tVar++;
                }
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) <= tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )
                        && tEdgePosition( 2 ) >= tProcRange( 4 ) && tEdgePosition( 2 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 2 );
                    tVar++;
                }
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) <= tProcRange( 1 ) && tEdgePosition( 1 ) >= tProcRange( 2 )
                        && tEdgePosition( 1 ) <= tProcRange( 3 ) && tEdgePosition( 2 ) == tProcRange( 4 )  )
                {
                    tProcShare(tVar) = aProcNeighbors( 4 );
                    tVar++;
                }
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) <= tProcRange( 1 ) && tEdgePosition( 1 ) >= tProcRange( 2 )
                        && tEdgePosition( 1 ) <= tProcRange( 3 ) && tEdgePosition( 2 ) == tProcRange( 5 )  )
                {
                    tProcShare(tVar) = aProcNeighbors( 5 );
                    tVar++;
                }
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) <= tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && tEdgePosition( 2 ) == tProcRange( 4 )  )
                {
                    tProcShare(tVar) = aProcNeighbors( 9 );
                    tVar++;
                }
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) <= tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )
                        && tEdgePosition( 2 ) == tProcRange( 4 )  )
                {
                    tProcShare(tVar) = aProcNeighbors( 2 );
                    tVar++;
                }
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) <= tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 2 )
                        && tEdgePosition( 2 ) == tProcRange( 5 )  )
                {
                    tProcShare(tVar) = aProcNeighbors( 19 );
                    tVar++;
                }
                if( tEdgePosition( 0 ) >= tProcRange( 0 ) && tEdgePosition( 0 ) <= tProcRange( 1 ) && tEdgePosition( 1 ) == tProcRange( 3 )
                        && tEdgePosition( 2 ) == tProcRange( 5 )  )
                {
                    tProcShare(tVar) = aProcNeighbors( 24 );
                    tVar++;
                }
            }
            else if( tEdgePosition( 3 ) < UINT_MAX && tEdgePosition( 4 ) < UINT_MAX && tEdgePosition( 5 ) < UINT_MAX  )
            {
                // Check if it is an edge in y-direction
                if( tEdgePosition( 3 ) >= tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) >= tProcRange( 2 )
                        && tEdgePosition( 4 ) <= tProcRange( 3 ) && tEdgePosition( 5 ) >= tProcRange( 4 ) && tEdgePosition( 5 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = tProcRank;
                    tVar++;
                }
                if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) <= tProcRange( 3 )
                        && tEdgePosition( 5 ) >= tProcRange( 4 ) && tEdgePosition( 5 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 3 );
                    tVar++;
                }
                if( tEdgePosition( 3 ) == tProcRange( 1 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) <= tProcRange( 3 )
                        && tEdgePosition( 5 ) >= tProcRange( 4 ) && tEdgePosition( 5 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 1 );
                    tVar++;
                }
                if( tEdgePosition( 3 ) >= tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) >= tProcRange( 2 )
                        && tEdgePosition( 4 ) <= tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 4 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 4 );
                    tVar++;
                }
                if( tEdgePosition( 3 ) >= tProcRange( 0 ) && tEdgePosition( 3 ) <= tProcRange( 1 ) && tEdgePosition( 4 ) >= tProcRange( 2 )
                        && tEdgePosition( 4 ) <= tProcRange( 3 ) && tEdgePosition( 5 ) == tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 5 );
                    tVar++;
                }
                if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) <= tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 7 );
                    tVar++;
                }
                if( tEdgePosition( 3 ) == tProcRange( 1 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) <= tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 4 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 11 );
                    tVar++;
                }
                if( tEdgePosition( 3 ) == tProcRange( 0 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) <= tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 21 );
                    tVar++;
                }
                if( tEdgePosition( 3 ) == tProcRange( 1 ) && tEdgePosition( 4 ) >= tProcRange( 2 ) && tEdgePosition( 4 ) <= tProcRange( 3 )
                        && tEdgePosition( 5 ) == tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 22 );
                    tVar++;
                }
            }
            else if( tEdgePosition( 6 ) < UINT_MAX && tEdgePosition( 7 ) < UINT_MAX && tEdgePosition( 8 ) < UINT_MAX  )
            {
                // Check if it is an edge in z-direction
                if( tEdgePosition( 6 ) >= tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) >= tProcRange( 2 )
                        && tEdgePosition( 7 ) <= tProcRange( 3 ) && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = tProcRank;
                    tVar++;
                }
                if( tEdgePosition( 6 ) == tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) >= tProcRange( 2 )
                        && tEdgePosition( 7 ) <= tProcRange( 3 ) && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 3 );
                    tVar++;
                }
                if( tEdgePosition( 6 ) == tProcRange( 1 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) >= tProcRange( 2 )
                        && tEdgePosition( 7 ) <= tProcRange( 3 ) && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 1 );
                    tVar++;
                }
                if( tEdgePosition( 6 ) >= tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 2 )
                        && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 0 );
                    tVar++;
                }
                if( tEdgePosition( 6 ) >= tProcRange( 0 ) && tEdgePosition( 6 ) <= tProcRange( 1 ) && tEdgePosition( 7 ) == tProcRange( 3 )
                        && tEdgePosition( 8 ) >= tProcRange( 4 ) && tEdgePosition( 8 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 2 );
                    tVar++;
                }
                if( tEdgePosition( 6 ) == tProcRange( 0 )  && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 8 );
                    tVar++;
                }
                if( tEdgePosition( 6 ) == tProcRange( 1 )  && tEdgePosition( 7 ) == tProcRange( 2 ) && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 15 );
                    tVar++;
                }
                if( tEdgePosition( 6 ) == tProcRange( 0 )  && tEdgePosition( 7 ) == tProcRange( 3 ) && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 16 );
                    tVar++;
                }
                if( tEdgePosition( 6 ) == tProcRange( 1 )  && tEdgePosition( 7 ) == tProcRange( 3 ) && tEdgePosition( 8 ) >= tProcRange( 4 )
                        && tEdgePosition( 8 ) <= tProcRange( 5 ) )
                {
                    tProcShare(tVar) = aProcNeighbors( 17 );
                    tVar++;
                }
            }
        }
    }
    tProcShare.resize( tVar, 1 );
    tProcShare = unique( tProcShare );
    if( isempty( tProcShare ) == 0 && tProcShare( tProcShare.length() - 1 ) == UINT_MAX )
    {
        tProcShare.resize( tProcShare.length() - 1, 1 );
    }
    return tProcShare;
}

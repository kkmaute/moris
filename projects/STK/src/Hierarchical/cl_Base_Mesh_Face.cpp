/*
 * cl_Base_Mesh_Face.cpp
 *
 *  Created on: Feb 27, 2018
 *      Author: gleim
 */

#include "cl_Base_Mesh_Face.hpp"
using namespace moris;

/* @TODO simplify formula */
uint
Base_Mesh_Face::give_number_of_faces_x(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // output variable
    uint tFaceNumber=0;
    if( aModelDim == 2 )
    {
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            // Count faces on each level until level "aLevel"
            tFaceNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) + 1 )
                         * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) );
        }
    }
    else if( aModelDim == 3 )
    {
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            tFaceNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) + 1 )
                         * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) )
                         * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) );
        }
    }
    return tFaceNumber;
}

//--------------------------------------------------------------------------------

/* @TODO simplify formula */
uint
Base_Mesh_Face::give_number_of_faces_y(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // output variable
    uint tFaceNumber=0;
    if( aModelDim == 2 )
    {
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            // Count faces on each level until level "aLevel"
            tFaceNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) )
                         * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) + 1 );
        }
    }
    else if( aModelDim == 3 )
    {
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            tFaceNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) )
                         * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) + 1 )
                         * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) );
        }
    }
    return tFaceNumber;
}

//--------------------------------------------------------------------------------
/* @TODO simplify formula */
uint
Base_Mesh_Face::give_number_of_faces_z(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // output variable
    uint tFaceNumber=0;
    MORIS_ASSERT( aModelDim == 3, " Faces in z-direction are only in 3D available");
    for( uint i = 0; i < aLevel + 1; i++ )
    {
        // Count faces on each level until level "aLevel"
        tFaceNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) )
                     * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) )
                     * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) + 1 );
    }
    return tFaceNumber;
}

//--------------------------------------------------------------------------------

/* @TODO simplify formula */
uint
Base_Mesh_Face::give_number_of_faces(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // output variable
    uint tFaceNumber=0;
    if(aModelDim == 2)
    {
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            // Count all faces on each level until level "aLevel"
            tFaceNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) + 1 )
                         * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) )
                         + ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) )
                         * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) + 1 );
        }
    }
    else if( aModelDim == 3 )
    {
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            tFaceNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) + 1 )
                         * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) ) * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) )
                         + ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) ) * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) + 1 )
                         * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) ) + ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) )
                         * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) ) * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) + 1 );
        }
    }
    return tFaceNumber;
}

//--------------------------------------------------------------------------------

/* @TODO simplify formula */
uint
Base_Mesh_Face::give_face_level(
        uint const & aFaceId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // output variable
    uint Face_level=1;
    // temporary variable for the while loop
    uint tLevel=1;
    //Compute the relation of the different levels by the power of the level
    uint tPowLevel = 1; //pow( 2, Face_level - 1 );
    if( aModelDim == 2 )
    {
        // temporary variable for the number of faces
        uint tNumber_of_faces = ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 )
                              * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) )
                              + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) )
                              * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) + 1 );
        while( tLevel > 0 )
        {
            tLevel = floor( (real)( aFaceId + 1 ) / tNumber_of_faces );
            if( (real)( aFaceId + 1 ) / tNumber_of_faces == 1 || tLevel < 1 )
            {
                tLevel = 0;
            }
            else
            {
                Face_level++;
                tPowLevel = pow( 2, Face_level - 1 );
                tNumber_of_faces += ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 )
                                  * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) )
                                  + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) )
                                  * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) + 1 );
            }
        }
        Face_level--;
    }
    else if( aModelDim == 3 )
    {
        // temporary variable for the number of faces
        uint tNumber_of_faces = ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 )
                              * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) ) * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) )
                              + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) ) * ( tPowLevel *aNumberOfElementsPerDirection( 1 ) + 1 )
                              * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) ) + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) )
                              * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) ) * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) + 1 );
        while( tLevel > 0 )
        {
            tLevel = floor( (real)( aFaceId + 1 ) / tNumber_of_faces );
            if( (real)( aFaceId + 1 ) / tNumber_of_faces == 1 || tLevel < 1 )
            {
                tLevel = 0;
            }
            else
            {
                Face_level++;
                tPowLevel = pow( 2, Face_level - 1 );
                tNumber_of_faces += ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 ) * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) )
                                  * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) ) + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) )
                                  * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) + 1 ) * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) )
                                  + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) ) * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) )
                                  * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) + 1 );
            }
        }
        Face_level--;
    }
    return Face_level;
}

//--------------------------------------------------------------------------------

/* @TODO simplify formula */

uint
Base_Mesh_Face::give_face_x_of_position(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aIJKPosition) const
{
    uint tFaceId = 0;
    //Compute the level minus 1 to compute the number of faces until the last level to determine the face on the current level
    uint tLevel = aLevel - 1;
    // Position of the face in y-direction
    uint tY_Position = 0;
    if( aModelDim == 2 )
    {
        tY_Position = aIJKPosition( 1 ) * ( aNumberOfElementsPerDirection( 0 ) + 1 );
        if( aLevel>0 )
        {
            //Give number of faces in x-direction for the last level to have the initial number at position 0,0
            tFaceId = give_number_of_faces_x(tLevel,aModelDim,aNumberOfElementsPerDirection);
            // Position of the face in y-direction for a specific level
            tY_Position = aIJKPosition( 1 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) * pow(2,aLevel) );
        }
        //Face with the calculated y-position
        tFaceId = tFaceId + aIJKPosition( 0 ) + tY_Position;
    }
    else if( aModelDim == 3 )
    {
        tY_Position = aIJKPosition( 1 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) )
                    + aIJKPosition( 2 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) ) * aNumberOfElementsPerDirection( 1 );
        if( aLevel>0 )
        {
            //Give number of faces in x-direction for the last level to have the initial number at position 0,0,0
            tFaceId = give_number_of_faces_x(tLevel,aModelDim,aNumberOfElementsPerDirection);
            // Position of the face in y-direction for a specific level
            tY_Position = aIJKPosition( 1 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) * pow(2,aLevel) )
                        + aIJKPosition( 2 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) * pow(2,aLevel) )
                        * ( aNumberOfElementsPerDirection( 1 ) * pow(2,aLevel) );
        }
        tFaceId = tFaceId + aIJKPosition( 0 ) + tY_Position; //face with the calculated y-position
    }
    return tFaceId;
}

//--------------------------------------------------------------------------------

/* @TODO simplify formula */

uint
Base_Mesh_Face::give_face_y_of_position(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aIJKPosition) const
{
    uint tFaceId = 0;
    //Compute the level minus 1 to compute the number of faces until the last level to determine the face on the current level
    uint tLevel = aLevel - 1;
    // Position of the face in y-direction
    uint tY_Position = 0;
    if( aModelDim == 2 )
    {
        tY_Position = aIJKPosition( 1 ) * aNumberOfElementsPerDirection( 0 );
        if( aLevel>0 )
        {
            //Give number of faces in x-direction for the last level to have the initial number at position 0,0
            tFaceId = give_number_of_faces_y(tLevel,aModelDim,aNumberOfElementsPerDirection);
            // Position of the face in y-direction for a specific level
            tY_Position = aIJKPosition( 1 ) * ( aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel ) );
        }
        //face with the calculated y-position
        tFaceId = tFaceId + aIJKPosition( 0 ) + tY_Position;
    }
    else if( aModelDim == 3 )
    {
        tY_Position = aIJKPosition( 1 ) * aNumberOfElementsPerDirection( 0 ) + aIJKPosition( 2 )
                    * aNumberOfElementsPerDirection( 0 ) * ( 1 + aNumberOfElementsPerDirection( 1 ) );
        if( aLevel > 0 )
        {
            //Give number of faces in x-direction for the last level to have the initial number at position 0,0,0
            tFaceId = give_number_of_faces_y(tLevel,aModelDim,aNumberOfElementsPerDirection);
            // Position of the face in y-direction for a specific level
            tY_Position = aIJKPosition( 1 ) * ( aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel ) )
                        + aIJKPosition( 2 ) * ( aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel ) )
                        * ( 1 + aNumberOfElementsPerDirection( 1 ) * pow( 2, aLevel ) );
        }
        //face with the calculated y-position
        tFaceId = tFaceId + aIJKPosition( 0 ) + tY_Position;
    }
    return tFaceId;
}

//--------------------------------------------------------------------------------

/* @TODO simplify formula */

uint
Base_Mesh_Face::give_face_z_of_position(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aIJKPosition) const
{
    MORIS_ASSERT( aModelDim == 3, " Edges in z-direction are only in 3D available" );
    uint tFaceId = 0;
    //Compute the level minus 1 to compute the number of faces until the last level to determine the face on the current level
    uint tLevel = aLevel - 1;
    // Position of the face in y-direction
    uint tY_Position = 0;
    tY_Position = aIJKPosition( 1 ) * ( aNumberOfElementsPerDirection( 0 ) )
            + aIJKPosition( 2 ) * aNumberOfElementsPerDirection( 0 ) * aNumberOfElementsPerDirection( 1 );
    if( aLevel>0 )
    {
        //Give number of face in x-direction for the last level to have the initial number at position 0,0,0
        tFaceId = give_number_of_faces_z(tLevel,aModelDim,aNumberOfElementsPerDirection);
        // Position of the face in y-direction for a specific level
        tY_Position = aIJKPosition( 1 ) * ( aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel ) )
                    + aIJKPosition( 2 ) * ( aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel ) )
                    * ( aNumberOfElementsPerDirection( 1 ) * pow( 2, aLevel ) );
    }
    //Face with the calculated y-position
    tFaceId = tFaceId + aIJKPosition( 0 ) + tY_Position;
    return tFaceId;
}

//--------------------------------------------------------------------------------

/* @TODO create graphic for face numbering */
Mat<uint>
Base_Mesh_Face::give_element_faces(
        uint const & aElementId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // Compute the Element position
    Mat<uint> tElementPosition = mBaseElement.give_position_of_element( aElementId, aModelDim, aNumberOfElementsPerDirection );
    // Determine the level of the element
    uint tElementLevel = mBaseElement.give_element_level( aElementId, aModelDim, aNumberOfElementsPerDirection );
    Mat<uint> tElementFaces;
    // Position for the face
    Mat<uint> aIJKPosition( aModelDim , 1 , 0 );
    if( aModelDim == 2 )
    {
        tElementFaces.set_size( 4, 1, UINT_MAX );
        aIJKPosition( 0 ) = tElementPosition( 0 );
        aIJKPosition( 1 ) = tElementPosition( 1 );
        uint tFaceXNumber1 = give_face_x_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tElementPosition( 0 ) + 1;
        aIJKPosition( 1 ) = tElementPosition( 1 );
        uint tFaceXNumber2 = give_face_x_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tElementPosition( 0 );
        aIJKPosition( 1 ) = tElementPosition( 1 );
        uint tFaceYNumber1 = give_face_y_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tElementPosition( 0 );
        aIJKPosition( 1 ) = tElementPosition( 1 ) + 1;
        uint tFaceYNumber2 = give_face_y_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        uint tFaceXNumber = give_number_of_faces_x( tElementLevel, aModelDim, aNumberOfElementsPerDirection );
        uint tFaceYNumber_old = 0;
        if( tElementLevel > 0 )
        {
            // Give number of faces from last level
            tFaceYNumber_old = give_number_of_faces_y( tElementLevel - 1, aModelDim, aNumberOfElementsPerDirection );
        }
        tElementFaces( 0 ) = tFaceXNumber1 + tFaceYNumber_old;
        tElementFaces( 1 ) = tFaceXNumber2 + tFaceYNumber_old;
        tElementFaces( 2 ) = tFaceYNumber1 + tFaceXNumber;
        tElementFaces( 3 ) = tFaceYNumber2 + tFaceXNumber;
    }
    else if( aModelDim == 3 )
    {
        tElementFaces.set_size(6,1,UINT_MAX);
        aIJKPosition( 0 ) = tElementPosition( 0 );
        aIJKPosition( 1 ) = tElementPosition( 1 );
        aIJKPosition( 2 ) = tElementPosition( 2 );
        uint tFaceXNumber1 = give_face_x_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tElementPosition( 0 ) + 1;
        aIJKPosition( 1 ) = tElementPosition( 1 );
        aIJKPosition( 2 ) = tElementPosition( 2 );
        uint tFaceXNumber2 = give_face_x_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tElementPosition( 0 );
        aIJKPosition( 1 ) = tElementPosition( 1 );
        aIJKPosition( 2 ) = tElementPosition( 2 );
        uint tFaceYNumber1 = give_face_y_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tElementPosition( 0 );
        aIJKPosition( 1 ) = tElementPosition( 1 ) + 1;
        aIJKPosition( 2 ) = tElementPosition( 2 );
        uint tFaceYNumber2 = give_face_y_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tElementPosition( 0 );
        aIJKPosition( 1 ) = tElementPosition( 1 );
        aIJKPosition( 2 ) = tElementPosition( 2 );
        uint tFaceZNumber1 = give_face_z_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        aIJKPosition( 0 ) = tElementPosition( 0 );
        aIJKPosition( 1 ) = tElementPosition( 1 );
        aIJKPosition( 2 ) = tElementPosition( 2 ) + 1;
        uint tFaceZNumber2 = give_face_z_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
        uint tFaceXNumber = give_number_of_faces_x( tElementLevel, aModelDim, aNumberOfElementsPerDirection );
        uint tFaceYNumber = give_number_of_faces_y( tElementLevel, aModelDim, aNumberOfElementsPerDirection );
        uint tFaceYNumber_old = 0;
        uint tFaceZNumber_old = 0;
        if( tElementLevel > 0 )
        {
            // Give number of faces from last level
            tFaceYNumber_old = give_number_of_faces_y( tElementLevel - 1, aModelDim, aNumberOfElementsPerDirection );
            tFaceZNumber_old = give_number_of_faces_z( tElementLevel - 1, aModelDim, aNumberOfElementsPerDirection );
        }
        tElementFaces( 0 ) = tFaceXNumber1 + tFaceYNumber_old + tFaceZNumber_old;
        tElementFaces( 1 ) = tFaceXNumber2 + tFaceYNumber_old + tFaceZNumber_old;
        tElementFaces( 2 ) = tFaceYNumber1 + tFaceXNumber + tFaceZNumber_old;
        tElementFaces( 3 ) = tFaceYNumber2 + tFaceXNumber + tFaceZNumber_old;
        tElementFaces( 4 ) = tFaceZNumber1 + tFaceXNumber + tFaceYNumber;
        tElementFaces( 5 ) = tFaceZNumber2 + tFaceXNumber + tFaceYNumber;
    }
    return tElementFaces;
}

//--------------------------------------------------------------------------------

Mat<uint>
Base_Mesh_Face::give_elements_of_face(
        uint const & aFaceId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // Compute the position of the edge
    Mat<uint> tFacePosition = give_face_position( aFaceId, aModelDim, aNumberOfElementsPerDirection );
    // Determine the level of the edge
    uint tFaceLevel = give_face_level( aFaceId, aModelDim, aNumberOfElementsPerDirection );
    Mat<uint> tElements( 2, 1 , 0 );
    if( aModelDim == 2 )
    {
        Mat<uint> tElementPositon( 2, 1, 0 );
        if( tFacePosition( 0 ) < UINT_MAX && tFacePosition( 1 ) < UINT_MAX )
        {
            tElementPositon( 0 ) = tFacePosition( 0 ) - 1;
            tElementPositon( 1 ) = tFacePosition( 1 );
            // Compute, based on the positon of the Face the positon of the element
            tElements( 0 ) = mBaseElement.give_element_of_position( tFaceLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tFacePosition( 0 );
            tElementPositon( 1 ) = tFacePosition( 1 );
            tElements( 1 ) = mBaseElement.give_element_of_position( tFaceLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
        }
        else if( tFacePosition( 2 ) < UINT_MAX && tFacePosition( 3 ) < UINT_MAX )
        {
            tElementPositon( 0 ) = tFacePosition( 2 );
            tElementPositon( 1 ) = tFacePosition( 3 ) - 1;
            tElements( 0 ) = mBaseElement.give_element_of_position( tFaceLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tFacePosition( 2 );
            tElementPositon( 1 ) = tFacePosition( 3 );
            tElements( 1 ) = mBaseElement.give_element_of_position( tFaceLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
        }
    }
    else if( aModelDim == 3 )
    {
        Mat<uint> tElementPositon(3,1);
        if( tFacePosition( 0 ) < UINT_MAX && tFacePosition( 1 ) < UINT_MAX && tFacePosition( 2 ) < UINT_MAX )
        {
            tElementPositon( 0 ) = tFacePosition( 0 ) - 1;
            tElementPositon( 1 ) = tFacePosition( 1 );
            tElementPositon( 2 ) = tFacePosition( 2 );
            tElements( 0 ) = mBaseElement.give_element_of_position( tFaceLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tFacePosition( 0 );
            tElementPositon( 1 ) = tFacePosition( 1 );
            tElementPositon( 2 ) = tFacePosition( 2 );
            tElements( 1 ) = mBaseElement.give_element_of_position( tFaceLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );

        }
        else if( tFacePosition( 3 ) < UINT_MAX && tFacePosition( 4 ) < UINT_MAX && tFacePosition( 5 ) < UINT_MAX )
        {
            tElementPositon( 0 ) = tFacePosition( 3 );
            tElementPositon( 1 ) = tFacePosition( 4 ) - 1;
            tElementPositon( 2 ) = tFacePosition( 5 );
            tElements( 0 ) = mBaseElement.give_element_of_position( tFaceLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tFacePosition( 3 );
            tElementPositon( 1 ) = tFacePosition( 4 );
            tElementPositon( 2 ) = tFacePosition( 5 );
            tElements( 1 ) = mBaseElement.give_element_of_position( tFaceLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
        }
        else if( tFacePosition( 6 ) < UINT_MAX && tFacePosition( 7 ) < UINT_MAX && tFacePosition( 8 ) < UINT_MAX )
        {
            tElementPositon( 0 ) = tFacePosition( 6 );
            tElementPositon( 1 ) = tFacePosition( 7 );
            tElementPositon( 2 ) = tFacePosition( 8 ) - 1;
            tElements( 0 ) = mBaseElement.give_element_of_position( tFaceLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
            tElementPositon( 0 ) = tFacePosition( 6 );
            tElementPositon( 1 ) = tFacePosition( 7 );
            tElementPositon( 2 ) = tFacePosition( 8 );
            tElements( 1 ) = mBaseElement.give_element_of_position( tFaceLevel, aModelDim, aNumberOfElementsPerDirection, tElementPositon );
        }
    }
    return tElements;
}

//--------------------------------------------------------------------------------

/* @TODO simplify formula */
Mat<uint>
Base_Mesh_Face::give_face_position(
        uint const & aFaceId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    uint tFaceLevel = give_face_level( aFaceId, aModelDim, aNumberOfElementsPerDirection );
    Mat<uint> tFacePosition;
    //Compute the relation of the different levels by the power of the level
    uint tPowLevel = pow( 2, tFaceLevel );
    if( aModelDim == 2 )
    {
        tFacePosition.set_size( 4, 1, UINT_MAX );
        uint tFaceNumber_old = 0;
        if( tFaceLevel > 0 )
        {
            // Number of faces until the last smaller level
            tFaceNumber_old = give_number_of_faces( tFaceLevel - 1, aModelDim, aNumberOfElementsPerDirection );
        }
        // Number of faces until the last smaller level plus number of faces in x-direction of the current level
        uint tFaceNumbery = tFaceNumber_old + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 )
                          * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) );

        if( aFaceId >= tFaceNumber_old && aFaceId < tFaceNumbery )
        {
            // Calculate the y-position
            tFacePosition( 1 ) = ceil( (real)( aFaceId + 1 - tFaceNumber_old ) / (real)( aNumberOfElementsPerDirection( 0 ) * tPowLevel + 1 ) );
            // Calculate the x-direction with the calculated y-position
            tFacePosition( 0 ) = aFaceId + 1 - tFaceNumber_old - ( tFacePosition( 1 ) - 1 )
                               * ( aNumberOfElementsPerDirection( 0 ) * tPowLevel + 1 );
            //  - 1 to get an indexed basis position, starting with zero
            tFacePosition( 1 )--;
            tFacePosition( 0 )--;
        }
        if( aFaceId >= tFaceNumbery)
        {
            tFacePosition( 3 ) = ceil( (real)( aFaceId + 1 - tFaceNumbery ) / (real)( aNumberOfElementsPerDirection( 0 ) * tPowLevel ) );
            tFacePosition( 2 ) = aFaceId + 1 - tFaceNumbery - ( tFacePosition( 3 ) - 1 )
                               * ( aNumberOfElementsPerDirection( 0 ) * tPowLevel );
            //  - 1 to get an indexed basis position, starting with zero
            tFacePosition( 3 )--;
            tFacePosition( 2 )--;
        }
    }
    else if( aModelDim == 3 )
    {
        tFacePosition.set_size( 9, 1, UINT_MAX );
        uint tFaceNumber_old = 0;
        if( tFaceLevel > 0 )
        {
            // Number of faces until the last smaller level
            tFaceNumber_old = give_number_of_faces( tFaceLevel - 1, aModelDim, aNumberOfElementsPerDirection );
        }
        // Number of faces until the last smaller level plus number of faces in x-direction of the current level
        uint tFaceNumbery = tFaceNumber_old + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 )
                          * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) ) * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) );
        // Number of faces until the last smaller level plus number of faces in x- and y-direction of the current level
        uint tFaceNumberz = tFaceNumber_old + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) + 1 )
                          * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) ) * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) )
                          + ( tPowLevel * aNumberOfElementsPerDirection( 0 ) ) * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) + 1 )
                          * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) );

        if( aFaceId >= tFaceNumber_old && aFaceId < tFaceNumbery )
        {
            // Calculate the z-position
            tFacePosition( 2 ) = ceil( (real)( aFaceId + 1 - tFaceNumber_old ) / (real)( ( aNumberOfElementsPerDirection( 0 ) * tPowLevel + 1 )
                               * ( aNumberOfElementsPerDirection( 1 ) * tPowLevel ) ) );
            // Calculate the y-position
            tFacePosition( 1 ) = ceil( (real)( aFaceId + 1 - tFaceNumber_old - ( tFacePosition( 2 ) - 1 ) * ( ( 1 + aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                               * ( aNumberOfElementsPerDirection( 1 ) * tPowLevel ) ) )
                               / (real)( 1 + aNumberOfElementsPerDirection( 0 ) * tPowLevel ) );
            // Calculate the x-direction with the calculated y-position
            tFacePosition( 0 ) = aFaceId + 1 - tFaceNumber_old - ( tFacePosition( 1 ) - 1 ) * ( 1 + aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                               - ( tFacePosition( 2 ) - 1 ) * ( ( 1 + aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                               * ( aNumberOfElementsPerDirection( 1 ) * tPowLevel ) );
            //  - 1 to get an indexed basis position, starting with zero
            tFacePosition( 2 )--;
            tFacePosition( 1 )--;
            tFacePosition( 0 )--;
        }
        else if( aFaceId >= tFaceNumbery && aFaceId < tFaceNumberz )
        {
            // Calculate the z-position
            tFacePosition( 5 ) = ceil( (real)( aFaceId + 1 - tFaceNumbery ) / (real)( ( aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                               * ( 1 + aNumberOfElementsPerDirection( 1 ) * tPowLevel ) ) );
            // Calculate the y-position
            tFacePosition( 4 ) = ceil( (real)( aFaceId + 1 - tFaceNumbery - ( tFacePosition( 5 ) - 1 ) * ( ( aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                               * ( 1 + aNumberOfElementsPerDirection( 1 ) * tPowLevel ) ) )  / (real)( aNumberOfElementsPerDirection( 0 ) * tPowLevel ) );
            // Calculate the x-direction with the calculated y-position
            tFacePosition( 3 ) = aFaceId + 1 - tFaceNumbery - ( tFacePosition( 4 ) - 1 ) * ( aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                               - ( tFacePosition( 5 ) - 1 ) * ( ( aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                               * ( 1 + aNumberOfElementsPerDirection( 1 ) * tPowLevel ) );
            //  - 1 to get an indexed basis position, starting with zero
            tFacePosition( 5 )--;
            tFacePosition( 4 )--;
            tFacePosition( 3 )--;
        }
        else if( aFaceId >= tFaceNumberz )
        {
            // Calculate the z-position
            tFacePosition( 8 ) = ceil( (real)( aFaceId + 1 - tFaceNumberz ) / (real)( ( aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                               * ( aNumberOfElementsPerDirection( 1 ) * tPowLevel ) ) );
            // Calculate the y-position
            tFacePosition( 7 ) = ceil( (real)( aFaceId + 1 - tFaceNumberz - ( tFacePosition( 8 ) - 1 ) * ( ( aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                               * ( aNumberOfElementsPerDirection( 1 ) * tPowLevel ) ) ) / (real)( aNumberOfElementsPerDirection( 0 ) * tPowLevel ) );
            // Calculate the x-direction with the calculated y-position
            tFacePosition( 6 ) = aFaceId + 1 - tFaceNumberz - ( tFacePosition( 7 ) - 1 ) * ( aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                               - ( tFacePosition( 8 ) - 1 ) * ( ( aNumberOfElementsPerDirection( 0 ) * tPowLevel )
                               * ( aNumberOfElementsPerDirection( 1 ) * tPowLevel ) );
            //  - 1 to get an indexed basis position, starting with zero
            tFacePosition( 8 )--;
            tFacePosition( 7 )--;
            tFacePosition( 6 )--;
        }
    }
    return tFacePosition;
}

//--------------------------------------------------------------------------------

/* @TODO reorder if statements */
uint
Base_Mesh_Face::give_face_owner(uint const & aFaceId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aProcRange,
        Mat<uint> const & aProcNeighbors) const
{
    //Compute the position of the face
    Mat<uint> tFacePosition = give_face_position( aFaceId, aModelDim, aNumberOfElementsPerDirection );
    // Determine the level of the face
    uint tFaceLevel = give_face_level( aFaceId, aModelDim, aNumberOfElementsPerDirection );
    uint tProcOwner = UINT_MAX;
    uint tProcRank = par_rank();
    if( par_size() > 1 )
    {
        Mat<uint> tProcRange( aProcRange.length() ,1 ,0 );
        tProcRange( 0 ) = aProcRange( 0 ) * pow( 2, tFaceLevel );
        tProcRange( 1 ) = aProcRange( 1 ) * pow( 2, tFaceLevel );
        tProcRange( 2 ) = aProcRange( 2 ) * pow( 2, tFaceLevel );
        tProcRange( 3 ) = aProcRange( 3 ) * pow( 2, tFaceLevel );
        if( aModelDim == 2 )
        {
            // Check if it is an edge in x-direction
            if( tFacePosition( 0 ) < UINT_MAX && tFacePosition( 1 ) < UINT_MAX )
            {
                if( tFacePosition( 0 ) > tProcRange( 0 ) && tFacePosition( 0 ) <= tProcRange( 1 )
                        && tFacePosition( 1 ) >= tProcRange( 2 ) && tFacePosition( 1 ) < tProcRange( 3 ) ) // Inner faces
                    tProcOwner = tProcRank;
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) < tProcRange( 3 ) && aProcNeighbors( 3 ) == UINT_MAX  ) // Left proc faces
                    tProcOwner =tProcRank;
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) < tProcRange( 3 ) && aProcNeighbors( 3 ) < UINT_MAX  ) // Left proc faces
                    tProcOwner =aProcNeighbors( 3 );
                else if( tFacePosition( 0 ) == tProcRange( 1 ) + 1 && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) < tProcRange( 3 ) && aProcNeighbors( 1 ) < UINT_MAX  ) // Right outer faces
                    tProcOwner =aProcNeighbors( 1 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) - 1 && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) < tProcRange( 3 ) && aProcNeighbors( 3 ) < UINT_MAX  ) // Left outer faces
                    tProcOwner =aProcNeighbors( 3 );
                else if( tFacePosition( 0 ) > tProcRange( 0 ) && tFacePosition( 0 ) <= tProcRange( 1 )
                        && tFacePosition( 1 ) == tProcRange( 3 ) && aProcNeighbors( 2 ) < UINT_MAX ) // Top outer faces
                    tProcOwner = aProcNeighbors( 2 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) == tProcRange( 3 )
                        && aProcNeighbors( 16 ) < UINT_MAX ) // Top left outer face
                    tProcOwner = aProcNeighbors( 16 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) == tProcRange( 3 )
                        && aProcNeighbors( 16 ) == UINT_MAX && aProcNeighbors( 2 ) < UINT_MAX ) // Top left outer face
                    tProcOwner = aProcNeighbors( 2 );
                else if( tFacePosition( 0 ) == tProcRange( 1 ) + 1 && tFacePosition( 1 ) == tProcRange( 3 )
                        && aProcNeighbors( 17 ) < UINT_MAX ) // Top right outer face
                    tProcOwner = aProcNeighbors( 17 );
                else if( tFacePosition( 0 ) > tProcRange( 0 ) && tFacePosition( 0 ) <= tProcRange( 1 )
                        && tFacePosition( 1 ) == tProcRange( 2 ) - 1 && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom outer faces
                    tProcOwner = aProcNeighbors( 0 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) == tProcRange( 2 ) - 1
                        && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer face
                    tProcOwner = aProcNeighbors( 8 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) == tProcRange( 2 ) - 1
                        && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom left outer face
                    tProcOwner = aProcNeighbors( 0 );
                else if( tFacePosition( 0 ) == tProcRange( 1 ) + 1 && tFacePosition( 1 ) == tProcRange( 2 ) - 1
                        && aProcNeighbors( 15 ) < UINT_MAX ) // Bottom right outer face
                    tProcOwner = aProcNeighbors( 15 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) - 1 && tFacePosition( 1 ) == tProcRange( 3 )
                        && aProcNeighbors( 16 ) < UINT_MAX  ) // Top left outer face
                    tProcOwner =aProcNeighbors( 16 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) - 1 && tFacePosition( 1 ) == tProcRange( 2 ) - 1
                        && aProcNeighbors( 8 ) < UINT_MAX  ) // Bottom left outer face
                    tProcOwner =aProcNeighbors( 8 );
            }
            else if( tFacePosition( 2 ) < UINT_MAX && tFacePosition( 3 ) < UINT_MAX )
            {
                // Check if it is an edge in y-direction
                if( tFacePosition( 2 ) >= tProcRange( 0 ) && tFacePosition( 2 ) < tProcRange( 1 )
                        && tFacePosition( 3 ) > tProcRange( 2 ) && tFacePosition( 3 ) <= tProcRange( 3 ) ) // Inner faces
                    tProcOwner =tProcRank;
                else if( tFacePosition( 2 ) >= tProcRange( 0 ) && tFacePosition( 2 ) < tProcRange( 1 )
                        && tFacePosition( 3 ) == tProcRange( 2 ) && aProcNeighbors( 0 ) == UINT_MAX ) // Bottom proc faces
                    tProcOwner =tProcRank;
                else if( tFacePosition( 2 ) >= tProcRange( 0 ) && tFacePosition( 2 ) < tProcRange( 1 )
                        && tFacePosition( 3 ) == tProcRange( 2 ) && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom proc faces
                    tProcOwner =aProcNeighbors( 0 );
                else if( tFacePosition( 2 ) == tProcRange( 1 ) && tFacePosition( 3 ) > tProcRange( 2 )
                        && tFacePosition( 3 ) <= tProcRange( 3 ) && aProcNeighbors( 1 ) < UINT_MAX ) // Right outer faces
                    tProcOwner =aProcNeighbors( 1 );
                else if( tFacePosition( 2 ) == tProcRange( 0 ) - 1 && tFacePosition( 3 ) > tProcRange( 2 )
                        && tFacePosition( 3 ) <= tProcRange( 3 ) && aProcNeighbors( 3 ) < UINT_MAX ) // Left outer faces
                    tProcOwner =aProcNeighbors( 3 );
                else if( tFacePosition( 2 ) == tProcRange( 0 ) - 1 && tFacePosition( 3 ) == tProcRange( 2 )
                        && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer face
                    tProcOwner =aProcNeighbors( 8 );
                else if( tFacePosition( 2 ) == tProcRange( 0 ) - 1 && tFacePosition( 3 ) == tProcRange( 2 )
                        && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 3 ) < UINT_MAX ) // Bottom left outer face
                    tProcOwner =aProcNeighbors( 3 );
                else if( tFacePosition( 2 ) == tProcRange( 1 ) && tFacePosition( 3 ) == tProcRange( 2 )
                        && aProcNeighbors( 15 ) < UINT_MAX ) // Bottom right outer face
                    tProcOwner =aProcNeighbors( 15 );
                else if( tFacePosition( 2 ) == tProcRange( 1 ) && tFacePosition( 3 ) == tProcRange( 2 )
                        && aProcNeighbors( 15 ) == UINT_MAX && aProcNeighbors( 1 ) < UINT_MAX ) // Bottom right outer face
                    tProcOwner =aProcNeighbors( 1 );
                else if( tFacePosition( 2 ) >= tProcRange( 0 ) && tFacePosition( 2 ) < tProcRange( 1 )
                        && tFacePosition( 3 ) == tProcRange( 3 ) + 1 && aProcNeighbors( 2 ) < UINT_MAX ) // Top outer faces
                    tProcOwner =aProcNeighbors( 2 );
                else if( tFacePosition( 2 ) == tProcRange( 1 ) && tFacePosition( 3 ) == tProcRange( 3 ) + 1
                        && aProcNeighbors( 17 ) < UINT_MAX ) // Top right outer face
                    tProcOwner =aProcNeighbors( 17 );
                else if( tFacePosition( 2 ) >= tProcRange( 0 ) && tFacePosition( 2 ) < tProcRange( 1 )
                        && tFacePosition( 3 ) == tProcRange( 2 ) - 1 && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom outer faces
                    tProcOwner =aProcNeighbors( 0 );
                else if( tFacePosition( 2 ) == tProcRange( 1 ) && tFacePosition( 3 ) == tProcRange( 2 ) - 1
                        && aProcNeighbors( 15 ) < UINT_MAX ) // Bottom right outer face
                    tProcOwner =aProcNeighbors( 15 );
                else if( tFacePosition( 2 ) == tProcRange( 0 ) - 1 && tFacePosition( 3 ) == tProcRange( 3 ) + 1
                        && aProcNeighbors( 16 ) < UINT_MAX ) // Top left outer face
                    tProcOwner =aProcNeighbors( 16 );
                else if( tFacePosition( 2 ) == tProcRange( 0 ) - 1 && tFacePosition( 3 ) == tProcRange( 2 ) - 1
                        && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer face
                    tProcOwner =aProcNeighbors( 8 );
            }
        }
        else if( aModelDim == 3)
        {
            tProcRange( 4 ) = aProcRange( 4 ) * pow( 2, tFaceLevel );
            tProcRange( 5 ) = aProcRange( 5 ) * pow( 2, tFaceLevel );
            // Check if it is an edge in x-direction
            if( tFacePosition( 0 ) < UINT_MAX && tFacePosition( 1 ) < UINT_MAX  && tFacePosition( 2 ) < UINT_MAX )
            {
                if( tFacePosition( 0 ) > tProcRange( 0 ) && tFacePosition( 0 ) <= tProcRange( 1 )
                        && tFacePosition( 1 ) >= tProcRange( 2 ) && tFacePosition( 1 ) < tProcRange( 3 )
                        && tFacePosition( 2 ) >= tProcRange( 4 ) && tFacePosition( 2 ) < tProcRange( 5 ) ) // Inner faces
                    tProcOwner = tProcRank;
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) < tProcRange( 3 ) && tFacePosition( 2 ) >= tProcRange( 4 )
                        && tFacePosition( 2 ) < tProcRange( 5 ) && aProcNeighbors( 3 ) < UINT_MAX ) // left faces
                    tProcOwner = aProcNeighbors( 3 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) < tProcRange( 3 ) && tFacePosition( 2 ) >= tProcRange( 4 )
                        && tFacePosition( 2 ) < tProcRange( 5 ) && aProcNeighbors( 3 ) == UINT_MAX )// left faces
                    tProcOwner = tProcRank;
                else if( tFacePosition( 0 ) == tProcRange( 1 ) + 1 && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) < tProcRange( 3 ) && tFacePosition( 2 ) >= tProcRange( 4 )
                        && tFacePosition( 2 ) < tProcRange( 5 ) && aProcNeighbors( 1 ) < UINT_MAX ) // right outer faces
                    tProcOwner = aProcNeighbors( 1 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) - 1 && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) < tProcRange( 3 ) && tFacePosition( 2 ) >= tProcRange( 4 )
                        && tFacePosition( 2 ) < tProcRange( 5 ) && aProcNeighbors( 3 ) < UINT_MAX ) // left outer faces
                    tProcOwner = aProcNeighbors( 3 );
                else if( tFacePosition( 0 ) > tProcRange( 0 ) && tFacePosition( 0 ) <= tProcRange( 1 )
                        && tFacePosition( 1 ) == tProcRange( 3 ) && tFacePosition( 2 ) >= tProcRange( 4 )
                        && tFacePosition( 2 ) < tProcRange( 5 ) && aProcNeighbors( 2 ) < UINT_MAX ) // Top outer faces
                    tProcOwner = aProcNeighbors( 2 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) == tProcRange( 3 )
                        && tFacePosition( 2 ) >= tProcRange( 4 ) && tFacePosition( 2 ) < tProcRange( 5 )
                        && aProcNeighbors( 16 ) < UINT_MAX ) // Top left outer face
                    tProcOwner = aProcNeighbors( 16 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) == tProcRange( 3 )
                        && tFacePosition( 2 ) >= tProcRange( 4 ) && tFacePosition( 2 ) < tProcRange( 5 )
                        && aProcNeighbors( 16 ) == UINT_MAX && aProcNeighbors( 2 ) < UINT_MAX  ) // Top left outer face
                    tProcOwner = aProcNeighbors( 2 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) - 1 && tFacePosition( 1 ) == tProcRange( 3 )
                        && tFacePosition( 2 ) >= tProcRange( 4 ) && tFacePosition( 2 ) < tProcRange( 5 )
                        && aProcNeighbors( 16 ) < UINT_MAX ) // Top left outer face
                    tProcOwner = aProcNeighbors( 16 );
                else if( tFacePosition( 0 ) > tProcRange( 0 ) && tFacePosition( 0 ) <= tProcRange( 1 )
                        && tFacePosition( 1 ) == tProcRange( 2 ) - 1 && tFacePosition( 2 ) >= tProcRange( 4 )
                        && tFacePosition( 2 ) < tProcRange( 5 ) && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom outer faces
                    tProcOwner = aProcNeighbors( 0 );
                else if( tFacePosition( 0 ) >= tProcRange( 0 ) - 1 && tFacePosition( 0 ) <= tProcRange( 0 )
                        && tFacePosition( 1 ) == tProcRange( 2 ) - 1 && tFacePosition( 2 ) >= tProcRange( 4 )
                        && tFacePosition( 2 ) < tProcRange( 5 ) && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer face
                    tProcOwner = aProcNeighbors( 8 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 2 ) >= tProcRange( 4 ) && tFacePosition( 2 ) < tProcRange( 5 )
                        && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 0 ) < UINT_MAX  ) // Bottom left outer face
                    tProcOwner = aProcNeighbors( 0 );
                else if( tFacePosition( 0 ) == tProcRange( 1 ) + 1 && tFacePosition( 1 ) == tProcRange( 3 )
                        && tFacePosition( 2 ) >= tProcRange( 4 ) && tFacePosition( 2 ) < tProcRange( 5 )
                        && aProcNeighbors( 17 ) < UINT_MAX ) // Top right outer face
                    tProcOwner = aProcNeighbors( 17 );
                else if( tFacePosition( 0 ) > tProcRange( 0 ) && tFacePosition( 0 ) <= tProcRange( 1 )
                        && tFacePosition( 1 ) >= tProcRange( 2 ) && tFacePosition( 1 ) < tProcRange( 3 )
                        && tFacePosition( 2 ) && tFacePosition( 2 ) == tProcRange( 5 ) && aProcNeighbors( 5 ) < UINT_MAX ) // Front outer faces
                    tProcOwner = aProcNeighbors( 5 );
                else if( tFacePosition( 0 ) >= tProcRange( 0 ) - 1 && tFacePosition( 0 ) <= tProcRange( 0 )
                        && tFacePosition( 1 ) >= tProcRange( 2 ) && tFacePosition( 1 ) < tProcRange( 3 )
                        && tFacePosition( 2 ) == tProcRange( 5 ) && aProcNeighbors( 21 ) < UINT_MAX ) // Front left outer face
                    tProcOwner = aProcNeighbors( 21 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) < tProcRange( 3 ) && tFacePosition( 2 ) == tProcRange( 5 )
                        && aProcNeighbors( 21 ) == UINT_MAX && aProcNeighbors( 5 ) < UINT_MAX  ) // Front left outer face
                    tProcOwner = aProcNeighbors( 5 );
                else if( tFacePosition( 0 ) == tProcRange( 1 ) + 1 && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) < tProcRange( 3 ) && tFacePosition( 2 ) == tProcRange( 5 )
                        && aProcNeighbors( 22 ) < UINT_MAX ) // Front right outer face
                    tProcOwner = aProcNeighbors( 22 );
                else if( tFacePosition( 0 ) > tProcRange( 0 ) && tFacePosition( 0 ) <= tProcRange( 1 )
                        && tFacePosition( 1 ) == tProcRange( 3 ) && tFacePosition( 2 ) && tFacePosition( 2 ) == tProcRange( 5 )
                        && aProcNeighbors( 24 ) < UINT_MAX ) // Top front outer faces
                    tProcOwner = aProcNeighbors( 24 );
                else if( tFacePosition( 0 ) >= tProcRange( 0 ) - 1 && tFacePosition( 0 ) <= tProcRange( 0 )
                        && tFacePosition( 1 ) == tProcRange( 3 ) && tFacePosition( 2 ) == tProcRange( 5 )
                        && aProcNeighbors( 23 ) < UINT_MAX ) // Top front left outer face
                    tProcOwner = aProcNeighbors( 23 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) == tProcRange( 3 )
                        && tFacePosition( 2 ) == tProcRange( 5 ) && aProcNeighbors( 23 ) == UINT_MAX
                        && aProcNeighbors( 24 ) < UINT_MAX  ) // Top front left outer face
                    tProcOwner = aProcNeighbors( 24 );
                else if( tFacePosition( 0 ) == tProcRange( 1 ) + 1 && tFacePosition( 1 ) == tProcRange( 3 )
                        && tFacePosition( 2 ) == tProcRange( 5 ) && aProcNeighbors( 25 ) < UINT_MAX ) // Top front right outer face
                    tProcOwner = aProcNeighbors( 25 );
                else if( tFacePosition( 0 ) > tProcRange( 0 ) && tFacePosition( 0 ) <= tProcRange( 1 )
                        && tFacePosition( 1 ) >= tProcRange( 2 ) && tFacePosition( 1 ) < tProcRange( 3 )
                        && tFacePosition( 2 ) && tFacePosition( 2 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 4 ) < UINT_MAX ) // Back outer faces
                    tProcOwner = aProcNeighbors( 4 );
                else if( tFacePosition( 0 ) >= tProcRange( 0 ) - 1 && tFacePosition( 0 ) <= tProcRange( 0 )
                        && tFacePosition( 1 ) >= tProcRange( 2 ) && tFacePosition( 1 ) < tProcRange( 3 )
                        && tFacePosition( 2 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 7 ) < UINT_MAX ) // Back left outer face
                    tProcOwner = aProcNeighbors( 7 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) < tProcRange( 3 ) && tFacePosition( 2 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 7 ) == UINT_MAX && aProcNeighbors( 4 ) < UINT_MAX  ) // Back left outer face
                    tProcOwner = aProcNeighbors( 4 );
                else if( tFacePosition( 0 ) == tProcRange( 1 ) + 1 && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) < tProcRange( 3 ) && tFacePosition( 2 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 11 ) < UINT_MAX ) // Back right outer face
                    tProcOwner = aProcNeighbors( 11 );
                else if( tFacePosition( 0 ) > tProcRange( 0 ) && tFacePosition( 0 ) <= tProcRange( 1 )
                        && tFacePosition( 1 ) == tProcRange( 3 ) && tFacePosition( 2 )
                        && tFacePosition( 2 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 13 ) < UINT_MAX ) // Top back outer faces
                    tProcOwner = aProcNeighbors( 13 );
                else if(  tFacePosition( 0 ) >= tProcRange( 0 ) - 1 && tFacePosition( 0 ) <= tProcRange( 0 )
                        && tFacePosition( 1 ) == tProcRange( 3 ) && tFacePosition( 2 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 12 ) < UINT_MAX ) // Top back left outer face
                    tProcOwner = aProcNeighbors( 12 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) == tProcRange( 3 )
                        && tFacePosition( 2 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 12 ) == UINT_MAX
                        && aProcNeighbors( 13 ) < UINT_MAX  ) // Top back left outer face
                    tProcOwner = aProcNeighbors( 13 );
                else if( tFacePosition( 0 ) == tProcRange( 1 ) + 1 && tFacePosition( 1 ) == tProcRange( 3 )
                        && tFacePosition( 2 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 14 ) < UINT_MAX ) // Top back right outer face
                    tProcOwner = aProcNeighbors( 14 );
                else if( tFacePosition( 0 ) > tProcRange( 0 ) && tFacePosition( 0 ) <= tProcRange( 1 )
                        && tFacePosition( 1 ) == tProcRange( 2 ) - 1 && tFacePosition( 2 ) && tFacePosition( 2 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 9 ) < UINT_MAX ) // Bottom back outer faces
                    tProcOwner = aProcNeighbors( 9 );
                else if( tFacePosition( 0 ) >= tProcRange( 0 ) - 1 && tFacePosition( 0 ) <= tProcRange( 0 )
                        && tFacePosition( 1 ) == tProcRange( 2 ) - 1 && tFacePosition( 2 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 6 ) < UINT_MAX ) // Bottom back left outer face
                    tProcOwner = aProcNeighbors( 6 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 2 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 6 ) == UINT_MAX
                        && aProcNeighbors( 9 ) < UINT_MAX  ) // Bottom back left outer face
                    tProcOwner = aProcNeighbors( 9 );
                else if( tFacePosition( 0 ) == tProcRange( 1 ) + 1 && tFacePosition( 1 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 2 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 10 ) < UINT_MAX ) // Bottom back right outer face
                    tProcOwner = aProcNeighbors( 10 );
                else if( tFacePosition( 0 ) == tProcRange( 1 ) + 1 && tFacePosition( 1 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 2 ) >= tProcRange( 4 ) && tFacePosition( 2 ) < tProcRange( 5 )
                        && aProcNeighbors( 15 ) < UINT_MAX ) // Top right outer face
                    tProcOwner = aProcNeighbors( 15 );
                else if( tFacePosition( 0 ) > tProcRange( 0 ) && tFacePosition( 0 ) <= tProcRange( 1 )
                        && tFacePosition( 1 ) == tProcRange( 2 ) - 1 && tFacePosition( 2 ) && tFacePosition( 2 ) == tProcRange( 5 )
                        && aProcNeighbors( 19 ) < UINT_MAX ) // Bottom front outer faces
                    tProcOwner = aProcNeighbors( 19 );
                else if( tFacePosition( 0 ) >= tProcRange( 0 ) - 1 && tFacePosition( 0 ) <= tProcRange( 0 )
                        && tFacePosition( 1 ) == tProcRange( 2 ) - 1 && tFacePosition( 2 ) == tProcRange( 5 )
                        && aProcNeighbors( 18 ) < UINT_MAX ) // Bottom front left outer face
                    tProcOwner = aProcNeighbors( 18 );
                else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 2 ) == tProcRange( 5 ) && aProcNeighbors( 18 ) == UINT_MAX && aProcNeighbors( 19 ) < UINT_MAX  ) // Bottom front left outer face
                    tProcOwner = aProcNeighbors( 19 );
                else if( tFacePosition( 0 ) == tProcRange( 1 ) + 1 && tFacePosition( 1 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 2 ) == tProcRange( 5 ) && aProcNeighbors( 20 ) < UINT_MAX ) // Bottom front right outer face
                    tProcOwner = aProcNeighbors( 20 );
            }
            else if( tFacePosition( 3 ) < UINT_MAX && tFacePosition( 4 ) < UINT_MAX && tFacePosition( 5 ) < UINT_MAX  )// Check if it is an edge in y-direction
            {
                if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) > tProcRange( 2 ) && tFacePosition( 4 ) <= tProcRange( 3 )
                        && tFacePosition( 5 ) >= tProcRange( 4 ) && tFacePosition( 5 ) < tProcRange( 5 ) ) // Inner faces
                    tProcOwner = tProcRank;
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) == tProcRange( 2 ) && tFacePosition( 5 ) >= tProcRange( 4 )
                        && tFacePosition( 5 ) < tProcRange( 5 ) && aProcNeighbors( 0 ) == UINT_MAX ) // Bottom faces
                    tProcOwner = tProcRank;
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) == tProcRange( 2 ) && tFacePosition( 5 ) >= tProcRange( 4 )
                        && tFacePosition( 5 ) < tProcRange( 5 ) && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom faces
                    tProcOwner = aProcNeighbors( 0 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) > tProcRange( 2 )
                        && tFacePosition( 4 ) <= tProcRange( 3 ) && tFacePosition( 5 ) >= tProcRange( 4 )
                        && tFacePosition( 5 ) < tProcRange( 5 ) && aProcNeighbors( 1 ) < UINT_MAX ) // Right outer faces
                    tProcOwner = aProcNeighbors( 1 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) == tProcRange( 2 )
                        && tFacePosition( 5 ) >= tProcRange( 4 ) && tFacePosition( 5 ) < tProcRange( 5 )
                        && aProcNeighbors( 15 ) < UINT_MAX ) // Bottom right outer faces
                    tProcOwner = aProcNeighbors( 15 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) == tProcRange( 2 )
                        && tFacePosition( 5 ) >= tProcRange( 4 ) && tFacePosition( 5 ) < tProcRange( 5 )
                        && aProcNeighbors( 15 ) == UINT_MAX && aProcNeighbors( 1 ) < UINT_MAX ) // Bottom right outer faces
                    tProcOwner = aProcNeighbors( 1 );
                else if( tFacePosition( 3 ) == tProcRange( 0 ) - 1 && tFacePosition( 4 ) > tProcRange( 2 )
                        && tFacePosition( 4 ) <= tProcRange( 3 ) && tFacePosition( 5 ) >= tProcRange( 4 )
                        && tFacePosition( 5 ) < tProcRange( 5 ) && aProcNeighbors( 3 ) < UINT_MAX ) // Left outer faces
                    tProcOwner = aProcNeighbors( 3 );
                else if( tFacePosition( 3 ) == tProcRange( 0 ) - 1 && tFacePosition( 4 ) == tProcRange( 2 )
                        && tFacePosition( 5 ) >= tProcRange( 4 ) && tFacePosition( 5 ) < tProcRange( 5 )
                        && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer faces
                    tProcOwner = aProcNeighbors( 8 );
                else if( tFacePosition( 3 ) == tProcRange( 0 ) - 1 && tFacePosition( 4 ) == tProcRange( 2 )
                        && tFacePosition( 5 ) >= tProcRange( 4 ) && tFacePosition( 5 ) < tProcRange( 5 )
                        && aProcNeighbors( 8 ) == UINT_MAX && aProcNeighbors( 3 ) < UINT_MAX ) // Bottom left outer faces
                    tProcOwner = aProcNeighbors( 3 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) > tProcRange( 2 ) && tFacePosition( 4 ) <= tProcRange( 3 )
                        && tFacePosition( 5 ) == tProcRange( 5 ) && aProcNeighbors( 5 ) < UINT_MAX ) // Front outer faces
                    tProcOwner = aProcNeighbors( 5 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) == tProcRange( 2 ) && tFacePosition( 5 ) == tProcRange( 5 )
                        && aProcNeighbors( 19 ) < UINT_MAX ) // Front bottom outer faces
                    tProcOwner = aProcNeighbors( 19 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) == tProcRange( 2 ) && tFacePosition( 5 ) == tProcRange( 5 )
                        && aProcNeighbors( 19 ) == UINT_MAX && aProcNeighbors( 5 ) < UINT_MAX ) // Front bottom outer faces
                    tProcOwner = aProcNeighbors( 5 );
                else if( tFacePosition( 3 ) == tProcRange( 0 ) - 1 && tFacePosition( 4 ) > tProcRange( 2 )
                        && tFacePosition( 4 ) <= tProcRange( 3 ) && tFacePosition( 5 ) == tProcRange( 5 )
                        && aProcNeighbors( 21 ) < UINT_MAX ) // Front outer faces
                    tProcOwner = aProcNeighbors( 21 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) - 1 && tFacePosition( 4 ) == tProcRange( 2 )
                        && tFacePosition( 5 ) == tProcRange( 5 ) && aProcNeighbors( 18 ) < UINT_MAX ) // Front left bottom outer faces
                    tProcOwner = aProcNeighbors( 18 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) - 1 && tFacePosition( 4 ) == tProcRange( 2 )
                        && tFacePosition( 5 ) == tProcRange( 5 ) && aProcNeighbors( 18 ) == UINT_MAX
                        && aProcNeighbors( 21 ) < UINT_MAX ) // Front left bottom outer faces
                    tProcOwner = aProcNeighbors( 21 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) > tProcRange( 2 )
                        && tFacePosition( 4 ) <= tProcRange( 3 ) && tFacePosition( 5 ) == tProcRange( 5 )
                        && aProcNeighbors( 22 ) < UINT_MAX ) // Front outer faces
                    tProcOwner = aProcNeighbors( 22 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) == tProcRange( 2 )
                        && tFacePosition( 5 ) == tProcRange( 5 ) && aProcNeighbors( 20 ) < UINT_MAX ) // Front bottom outer faces
                    tProcOwner = aProcNeighbors( 20 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) == tProcRange( 2 )
                        && tFacePosition( 5 ) == tProcRange( 5 ) && aProcNeighbors( 20 ) == UINT_MAX
                        && aProcNeighbors( 22 ) < UINT_MAX ) // Front bottom outer faces
                    tProcOwner = aProcNeighbors( 22 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) == tProcRange( 3 ) + 1 && tFacePosition( 5 ) >= tProcRange( 4 )
                        && tFacePosition( 5 ) < tProcRange( 5 ) && aProcNeighbors( 2 ) < UINT_MAX ) // Top outer faces
                    tProcOwner = aProcNeighbors( 2 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) == tProcRange( 3 ) + 1
                        && tFacePosition( 5 ) >= tProcRange( 4 ) && tFacePosition( 5 ) < tProcRange( 5 )
                        && aProcNeighbors( 17 ) < UINT_MAX ) // Top right outer face
                    tProcOwner = aProcNeighbors( 17 );
                else if( tFacePosition( 3 ) == tProcRange( 0 ) - 1 && tFacePosition( 4 ) == tProcRange( 3 ) + 1
                        && tFacePosition( 5 ) >= tProcRange( 4 ) && tFacePosition( 5 ) < tProcRange( 5 )
                        && aProcNeighbors( 16 ) < UINT_MAX ) // Top left outer face
                    tProcOwner = aProcNeighbors( 16 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) == tProcRange( 3 ) + 1 && tFacePosition( 5 ) == tProcRange( 5 )
                        && aProcNeighbors( 24 ) < UINT_MAX ) // Top front outer faces
                    tProcOwner = aProcNeighbors( 24 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) == tProcRange( 3 ) + 1
                        && tFacePosition( 5 ) == tProcRange( 5 ) && aProcNeighbors( 25 ) < UINT_MAX ) // Top front right outer face
                    tProcOwner = aProcNeighbors( 25 );
                else if( tFacePosition( 3 ) == tProcRange( 0 ) - 1 && tFacePosition( 4 ) == tProcRange( 3 ) + 1
                        && tFacePosition( 5 ) == tProcRange( 5 ) && aProcNeighbors( 23 ) < UINT_MAX ) // Top front left outer face
                    tProcOwner = aProcNeighbors( 23 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) > tProcRange( 2 ) && tFacePosition( 4 ) <= tProcRange( 3 )
                        && tFacePosition( 5 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 4 ) < UINT_MAX ) // Back outer faces
                    tProcOwner = aProcNeighbors( 4 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) == tProcRange( 2 ) && tFacePosition( 5 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 9 ) < UINT_MAX ) // Back bottom outer faces
                    tProcOwner = aProcNeighbors( 9 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) == tProcRange( 2 ) && tFacePosition( 5 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 9 ) == UINT_MAX && aProcNeighbors( 4 ) < UINT_MAX ) // Back bottom outer faces
                    tProcOwner = aProcNeighbors( 4 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) > tProcRange( 2 )
                        && tFacePosition( 4 ) <= tProcRange( 3 ) && tFacePosition( 5 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 11 ) < UINT_MAX ) // Back outer faces
                    tProcOwner = aProcNeighbors( 11 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) == tProcRange( 2 )
                        && tFacePosition( 5 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 10 ) < UINT_MAX ) // Back bottom outer faces
                    tProcOwner = aProcNeighbors( 10 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) == tProcRange( 2 )
                        && tFacePosition( 5 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 10 ) == UINT_MAX
                        && aProcNeighbors( 11 ) < UINT_MAX ) // Back bottom outer faces
                    tProcOwner = aProcNeighbors( 11 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 5 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 10 ) < UINT_MAX ) // Back bottom right outer faces
                    tProcOwner = aProcNeighbors( 10 );
                else if( tFacePosition( 3 ) == tProcRange( 0 ) - 1 && tFacePosition( 4 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 5 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 6 ) < UINT_MAX ) // Back bottom left outer faces
                    tProcOwner = aProcNeighbors( 6 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) == tProcRange( 2 ) - 1 && tFacePosition( 5 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 9 ) < UINT_MAX ) // Back bottom outer faces
                    tProcOwner = aProcNeighbors( 9 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) == tProcRange( 3 ) + 1 && tFacePosition( 5 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 13 ) < UINT_MAX ) // Top back outer faces
                    tProcOwner = aProcNeighbors( 13 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) == tProcRange( 3 ) + 1
                        && tFacePosition( 5 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 14 ) < UINT_MAX ) // Top back right outer face
                    tProcOwner = aProcNeighbors( 14 );
                else if( tFacePosition( 3 ) == tProcRange( 0 ) - 1 && tFacePosition( 4 ) == tProcRange( 3 ) + 1
                        && tFacePosition( 5 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 12 ) < UINT_MAX ) // Top back left outer face
                    tProcOwner = aProcNeighbors( 12 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) == tProcRange( 2 ) - 1 && tFacePosition( 5 ) >= tProcRange( 4 )
                        && tFacePosition( 5 ) < tProcRange( 5 ) && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom outer faces
                    tProcOwner = aProcNeighbors( 0 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 5 ) >= tProcRange( 4 ) && tFacePosition( 5 ) < tProcRange( 5 )
                        && aProcNeighbors( 15 ) < UINT_MAX ) // Bottom right outer face
                    tProcOwner = aProcNeighbors( 15 );
                else if( tFacePosition( 3 ) == tProcRange( 0 ) - 1 && tFacePosition( 4 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 5 ) >= tProcRange( 4 ) && tFacePosition( 5 ) < tProcRange( 5 )
                        && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer face
                    tProcOwner = aProcNeighbors( 8 );
                else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) == tProcRange( 2 ) - 1 && tFacePosition( 5 ) == tProcRange( 5 )
                        && aProcNeighbors( 19 ) < UINT_MAX ) // Bottom front outer faces
                    tProcOwner = aProcNeighbors( 19 );
                else if( tFacePosition( 3 ) == tProcRange( 1 ) && tFacePosition( 4 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 5 ) == tProcRange( 5 ) && aProcNeighbors( 20 ) < UINT_MAX ) // Bottom front right outer face
                    tProcOwner = aProcNeighbors( 20 );
                else if( tFacePosition( 3 ) == tProcRange( 0 ) - 1 && tFacePosition( 4 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 5 ) == tProcRange( 5 ) && aProcNeighbors( 18 ) < UINT_MAX ) // Bottom front left outer face
                    tProcOwner = aProcNeighbors( 18 );
                else if( tFacePosition( 3 ) == tProcRange( 0 ) - 1 && tFacePosition( 4 ) > tProcRange( 2 )
                        && tFacePosition( 4 ) <= tProcRange( 3 ) && tFacePosition( 5 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 7 ) < UINT_MAX ) // Back left  outer faces
                    tProcOwner = aProcNeighbors( 7 );
                else if( tFacePosition( 3 ) == tProcRange( 0 ) - 1 && tFacePosition( 4 ) == tProcRange( 2 )
                        && tFacePosition( 5 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 6 ) < UINT_MAX ) // Back left  outer faces
                    tProcOwner = aProcNeighbors( 6 );
                else if( tFacePosition( 3 ) == tProcRange( 0 ) - 1 && tFacePosition( 4 ) == tProcRange( 2 )
                        && tFacePosition( 5 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 6 ) == UINT_MAX
                        && aProcNeighbors( 7 ) < UINT_MAX ) // Back left  outer faces
                    tProcOwner = aProcNeighbors( 7 );
            }
            else if( tFacePosition( 6 ) < UINT_MAX && tFacePosition( 7 ) < UINT_MAX && tFacePosition( 8 ) < UINT_MAX  )// Check if it is an edge in z-direction
            {
                if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) >= tProcRange( 2 ) && tFacePosition( 7 ) < tProcRange( 3 )
                        && tFacePosition( 8 ) > tProcRange( 4 ) && tFacePosition( 8 ) <= tProcRange( 5 ) ) // Inner faces
                    tProcOwner = tProcRank;
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) >= tProcRange( 2 ) && tFacePosition( 7 ) < tProcRange( 3 )
                        && tFacePosition( 8 ) == tProcRange( 4 ) && aProcNeighbors( 4 ) == UINT_MAX  ) // Back faces
                    tProcOwner = tProcRank;
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) >= tProcRange( 2 ) && tFacePosition( 7 ) < tProcRange( 3 )
                        && tFacePosition( 8 ) == tProcRange( 4 ) && aProcNeighbors( 4 ) < UINT_MAX  ) // Back faces
                    tProcOwner = aProcNeighbors( 4 );
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) >= tProcRange( 2 ) && tFacePosition( 7 ) < tProcRange( 3 )
                        && tFacePosition( 8 ) > tProcRange( 4 ) && tFacePosition( 8 ) <= tProcRange( 5 ) ) // Inner faces
                    tProcOwner = tProcRank;
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) >= tProcRange( 2 )
                        && tFacePosition( 7 ) < tProcRange( 3 ) && tFacePosition( 8 ) > tProcRange( 4 )
                        && tFacePosition( 8 ) <= tProcRange( 5 ) && aProcNeighbors( 1 ) < UINT_MAX ) // Right outer faces
                    tProcOwner = aProcNeighbors( 1 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) >= tProcRange( 2 )
                        && tFacePosition( 7 ) < tProcRange( 3 ) && tFacePosition( 8 ) == tProcRange( 4 )
                        && aProcNeighbors( 11 ) < UINT_MAX ) // Right bottom outer faces
                    tProcOwner = aProcNeighbors( 11 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) >= tProcRange( 2 )
                        && tFacePosition( 7 ) < tProcRange( 3 ) && tFacePosition( 8 ) == tProcRange( 4 )
                        && aProcNeighbors( 11 ) == UINT_MAX && aProcNeighbors( 1 ) < UINT_MAX ) // Right bottom outer faces
                    tProcOwner = aProcNeighbors( 1 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) >= tProcRange( 2 )
                        && tFacePosition( 7 ) < tProcRange( 3 ) && tFacePosition( 8 ) > tProcRange( 4 )
                        && tFacePosition( 8 ) <= tProcRange( 5 ) && aProcNeighbors( 3 ) < UINT_MAX ) // Left outer faces
                    tProcOwner = aProcNeighbors( 3 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) >= tProcRange( 2 )
                        && tFacePosition( 7 ) < tProcRange( 3 ) && tFacePosition( 8 ) == tProcRange( 4 )
                        && aProcNeighbors( 7 ) < UINT_MAX ) // Left bottom outer faces
                    tProcOwner = aProcNeighbors( 7 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) >= tProcRange( 2 )
                        && tFacePosition( 7 ) < tProcRange( 3 ) && tFacePosition( 8 ) == tProcRange( 4 )
                        && aProcNeighbors( 7 ) == UINT_MAX && aProcNeighbors( 3 ) < UINT_MAX ) // Left bottom outer faces
                    tProcOwner = aProcNeighbors( 3 );
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) == tProcRange( 3 ) && tFacePosition( 8 ) > tProcRange( 4 )
                        && tFacePosition( 8 ) <= tProcRange( 5 ) && aProcNeighbors( 2 ) < UINT_MAX ) // Top outer faces
                    tProcOwner = aProcNeighbors( 2 );
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) == tProcRange( 3 ) && tFacePosition( 8 ) == tProcRange( 4 )
                        && aProcNeighbors( 13 ) < UINT_MAX ) // Top bottom outer faces
                    tProcOwner = aProcNeighbors( 13 );
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) == tProcRange( 3 ) && tFacePosition( 8 ) == tProcRange( 4 )
                        && aProcNeighbors( 13 ) == UINT_MAX && aProcNeighbors( 2 ) < UINT_MAX ) // Top bottom outer faces
                    tProcOwner = aProcNeighbors( 2 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) == tProcRange( 3 )
                        && tFacePosition( 8 ) > tProcRange( 4 ) && tFacePosition( 8 ) <= tProcRange( 5 )
                        && aProcNeighbors( 17 ) < UINT_MAX ) // Top right outer faces
                    tProcOwner = aProcNeighbors( 17 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) == tProcRange( 3 )
                        && tFacePosition( 8 ) > tProcRange( 4 ) && tFacePosition( 8 ) <= tProcRange( 5 )
                        && aProcNeighbors( 16 ) < UINT_MAX ) // Top left outer faces
                    tProcOwner = aProcNeighbors( 16 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) == tProcRange( 3 )
                        && tFacePosition( 8 ) == tProcRange( 4 )  && aProcNeighbors( 12 ) < UINT_MAX ) // Top left back outer faces
                    tProcOwner = aProcNeighbors( 12 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) == tProcRange( 3 )
                        && tFacePosition( 8 ) == tProcRange( 4 )  && aProcNeighbors( 12 ) == UINT_MAX
                        && aProcNeighbors( 16 ) < UINT_MAX ) // Top left back outer faces
                    tProcOwner = aProcNeighbors( 16 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) == tProcRange( 3 )
                        && tFacePosition( 8 ) == tProcRange( 4 )  && aProcNeighbors( 14 ) < UINT_MAX ) // Top right bottom outer faces
                    tProcOwner = aProcNeighbors( 14 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) == tProcRange( 3 )
                        && tFacePosition( 8 ) == tProcRange( 4 )  && aProcNeighbors( 14 ) == UINT_MAX
                        && aProcNeighbors( 17 ) < UINT_MAX ) // Top right bottom outer faces
                    tProcOwner = aProcNeighbors( 17 );
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) >= tProcRange( 2 ) && tFacePosition( 7 ) < tProcRange( 3 )
                        && tFacePosition( 8 ) == tProcRange( 5 ) + 1 && aProcNeighbors( 5 ) < UINT_MAX) // Front outer faces
                    tProcOwner = aProcNeighbors( 5 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) >= tProcRange( 2 )
                        && tFacePosition( 7 ) < tProcRange( 3 ) && tFacePosition( 8 ) == tProcRange( 5 ) + 1
                        && aProcNeighbors( 22 ) < UINT_MAX) // Front right outer faces
                    tProcOwner = aProcNeighbors( 22 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) >= tProcRange( 2 )
                        && tFacePosition( 7 ) < tProcRange( 3 ) && tFacePosition( 8 ) == tProcRange( 5 ) + 1
                        && aProcNeighbors( 21 ) < UINT_MAX) // Front left outer faces
                    tProcOwner = aProcNeighbors( 21 );
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) == tProcRange( 3 ) && tFacePosition( 8 ) == tProcRange( 5 ) + 1
                        && aProcNeighbors( 24 ) < UINT_MAX) // Front top outer faces
                    tProcOwner = aProcNeighbors( 24 );
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) == tProcRange( 2 ) - 1 && tFacePosition( 8 ) == tProcRange( 5 ) + 1
                        && aProcNeighbors( 19 ) < UINT_MAX) // Front bottom outer faces
                    tProcOwner = aProcNeighbors( 19 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) == tProcRange( 3 )
                        && tFacePosition( 8 ) == tProcRange( 5 ) + 1 && aProcNeighbors( 25 ) < UINT_MAX) // Front top right outer face
                    tProcOwner = aProcNeighbors( 25 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) == tProcRange( 3 )
                        && tFacePosition( 8 ) == tProcRange( 5 ) + 1 && aProcNeighbors( 23 ) < UINT_MAX) // Front top left outer face
                    tProcOwner = aProcNeighbors( 23 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 8 ) == tProcRange( 5 ) + 1 && aProcNeighbors( 20 ) < UINT_MAX) // Front bottom right outer face
                    tProcOwner = aProcNeighbors( 20 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 8 ) == tProcRange( 5 ) + 1 && aProcNeighbors( 18 ) < UINT_MAX) // Front bottom left outer face
                    tProcOwner = aProcNeighbors( 18 );
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) >= tProcRange( 2 ) && tFacePosition( 7 ) < tProcRange( 3 )
                        && tFacePosition( 8 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 4 ) < UINT_MAX) // Back outer faces
                    tProcOwner = aProcNeighbors( 4 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) >= tProcRange( 2 )
                        && tFacePosition( 7 ) < tProcRange( 3 ) && tFacePosition( 8 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 11 ) < UINT_MAX) // Back Bottom right outer faces
                    tProcOwner = aProcNeighbors( 11 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) >= tProcRange( 2 )
                        && tFacePosition( 7 ) < tProcRange( 3 ) && tFacePosition( 8 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 7 ) < UINT_MAX) // Back Bottom left outer faces
                    tProcOwner = aProcNeighbors( 7 );
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) == tProcRange( 3 ) && tFacePosition( 8 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 13 ) < UINT_MAX) // Back Bottom top outer faces
                    tProcOwner = aProcNeighbors( 13 );
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) == tProcRange( 2 ) - 1 && tFacePosition( 8 ) == tProcRange( 4 ) - 1
                        && aProcNeighbors( 9 ) < UINT_MAX) // Back Bottom bottom outer faces
                    tProcOwner = aProcNeighbors( 9 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) == tProcRange( 3 )
                        && tFacePosition( 8 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 14 ) < UINT_MAX) // Back top right outer faces
                    tProcOwner = aProcNeighbors( 14 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 8 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 10 ) < UINT_MAX) // Back bottom right outer faces
                    tProcOwner = aProcNeighbors( 10 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) == tProcRange( 3 )
                        && tFacePosition( 8 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 12 ) < UINT_MAX) // Back top left outer faces
                    tProcOwner = aProcNeighbors( 12 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 8 ) == tProcRange( 4 ) - 1 && aProcNeighbors( 6 ) < UINT_MAX) // Back bottom left outer faces
                    tProcOwner = aProcNeighbors( 6 );
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) == tProcRange( 2 ) - 1 && tFacePosition( 8 ) > tProcRange( 4 )
                        && tFacePosition( 8 ) <= tProcRange( 5 ) && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom outer faces
                    tProcOwner = aProcNeighbors( 0 );
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) == tProcRange( 2 ) - 1 && tFacePosition( 8 ) == tProcRange( 4 )
                        && aProcNeighbors( 9 ) < UINT_MAX ) // Bottom bottom outer faces
                    tProcOwner = aProcNeighbors( 9 );
                else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) == tProcRange( 2 ) - 1 && tFacePosition( 8 ) == tProcRange( 4 )
                        && aProcNeighbors( 9 ) == UINT_MAX && aProcNeighbors( 0 ) < UINT_MAX ) // Bottom bottom outer faces
                    tProcOwner = aProcNeighbors( 0 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 8 ) > tProcRange( 4 ) && tFacePosition( 8 ) <= tProcRange( 5 )
                        && aProcNeighbors( 15 ) < UINT_MAX ) // Bottom right outer faces
                    tProcOwner = aProcNeighbors( 15 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 8 ) == tProcRange( 4 )  && aProcNeighbors( 10 ) < UINT_MAX ) // Bottom right bottom outer faces
                    tProcOwner = aProcNeighbors( 10 );
                else if( tFacePosition( 6 ) == tProcRange( 1 ) && tFacePosition( 7 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 8 ) == tProcRange( 4 )  && aProcNeighbors( 10 ) == UINT_MAX
                        && aProcNeighbors( 15 ) < UINT_MAX ) // Bottom right bottom outer faces
                    tProcOwner = aProcNeighbors( 15 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 8 ) > tProcRange( 4 ) && tFacePosition( 8 ) <= tProcRange( 5 )
                        && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left outer faces
                    tProcOwner = aProcNeighbors( 8 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 8 ) == tProcRange( 4 )  && aProcNeighbors( 6 ) < UINT_MAX ) // Bottom left bottom outer faces
                    tProcOwner = aProcNeighbors( 6 );
                else if( tFacePosition( 6 ) == tProcRange( 0 ) - 1 && tFacePosition( 7 ) == tProcRange( 2 ) - 1
                        && tFacePosition( 8 ) == tProcRange( 4 )  && aProcNeighbors( 6 ) == UINT_MAX
                        && aProcNeighbors( 8 ) < UINT_MAX ) // Bottom left bottom outer faces
                    tProcOwner = aProcNeighbors( 8 );
            }
        }
    }
    return tProcOwner;
}

//--------------------------------------------------------------------------------

Mat<uint>
Base_Mesh_Face::give_face_share(uint const & aFaceId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aProcRange,
        Mat<uint> const & aProcNeighbors) const
{
    //Compute the positon of the face
    Mat<uint> tFacePosition = give_face_position( aFaceId, aModelDim, aNumberOfElementsPerDirection );
    // Determine the level of the face
    uint tFaceLevel = give_face_level( aFaceId,aModelDim, aNumberOfElementsPerDirection );
    Mat<uint> tProcShare( 40, 1, UINT_MAX );
    uint tVar = 0;
    uint tProcRank = par_rank();

    if( par_size() > 1 )
    {
        Mat<uint> tProcRange( aProcRange.length() ,1 ,0 );
        tProcRange( 0 ) = aProcRange( 0 ) * pow( 2, tFaceLevel );
        tProcRange( 1 ) = aProcRange( 1 ) * pow( 2, tFaceLevel );
        tProcRange( 2 ) = aProcRange( 2 ) * pow( 2, tFaceLevel );
        tProcRange( 3 ) = aProcRange( 3 ) * pow( 2, tFaceLevel );

        if( aModelDim == 2 )
        {
            // Check if it is an edge in x-direction
            if( tFacePosition( 0 ) < UINT_MAX && tFacePosition( 1 ) < UINT_MAX )
            {
                if( tFacePosition( 0 ) >= tProcRange( 0 ) && tFacePosition( 0 ) <= tProcRange( 1 )
                        && tFacePosition( 1 ) >= tProcRange( 2 ) && tFacePosition( 1 ) <= tProcRange( 3 ) )
                {
                    tProcShare(tVar) =tProcRank;
                    tVar++;
                }
                if( tFacePosition( 0 ) == tProcRange( 0 )  && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) <= tProcRange( 3 ) )
                {
                    tProcShare(tVar) =aProcNeighbors( 3 );
                    tVar++;
                }
                if( tFacePosition( 0 ) == tProcRange( 1 )  && tFacePosition( 1 ) >= tProcRange( 2 )
                        && tFacePosition( 1 ) <= tProcRange( 3 ) )
                {
                    tProcShare(tVar) =aProcNeighbors( 1 );
                    tVar++;
                }
            }
            else if( tFacePosition( 2 ) < UINT_MAX && tFacePosition( 3 ) < UINT_MAX )
            {
                // Check if it is an edge in y-direction
                if( tFacePosition( 2 ) >= tProcRange( 0 ) && tFacePosition( 2 ) <= tProcRange( 1 )
                        && tFacePosition( 3 ) >= tProcRange( 2 ) && tFacePosition( 3 ) <= tProcRange( 3 ) )
                {
                    tProcShare(tVar) =tProcRank;
                    tVar++;
                }
                if( tFacePosition( 2 ) >= tProcRange( 0 ) && tFacePosition( 2 ) <= tProcRange( 1 )
                        && tFacePosition( 3 ) == tProcRange( 2 ) )
                {
                    tProcShare(tVar) =aProcNeighbors( 0 );
                    tVar++;
                }
                if( tFacePosition( 2 ) >= tProcRange( 0 ) && tFacePosition( 2 ) <= tProcRange( 1 )
                        && tFacePosition( 3 ) == tProcRange( 3 ) )
                {
                    tProcShare(tVar) =aProcNeighbors( 2 );
                    tVar++;
                }
            }
        }
        else if( aModelDim == 3)
        {
            tProcRange( 4 ) = aProcRange( 4 ) * pow( 2, tFaceLevel );
            tProcRange( 5 ) = aProcRange( 5 ) * pow( 2, tFaceLevel );
            // Check if it is an edge in x-direction
            if( tFacePosition( 0 ) < UINT_MAX && tFacePosition( 1 ) < UINT_MAX  && tFacePosition( 2 ) < UINT_MAX )
            {
                if( tFacePosition( 0 ) > tProcRange( 0 ) && tFacePosition( 0 ) < tProcRange( 1 )
                        && tFacePosition( 1 ) >= tProcRange( 2 ) && tFacePosition( 1 ) < tProcRange( 3 )
                        && tFacePosition( 2 ) >= tProcRange( 4 ) && tFacePosition( 2 ) < tProcRange( 5 ) )
                {
                    tProcShare(tVar) = tProcRank;
                    tVar++;
                }
                else
                {
                    if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) >= tProcRange( 2 )
                            && tFacePosition( 1 ) < tProcRange( 3 ) && tFacePosition( 2 ) >= tProcRange( 4 )
                            && tFacePosition( 2 ) < tProcRange( 5 ) && aProcNeighbors( 3 ) == UINT_MAX )
                    {
                        tProcShare(tVar) = tProcRank;
                        tVar++;
                    }
                    else if( tFacePosition( 0 ) == tProcRange( 0 ) && tFacePosition( 1 ) >= tProcRange( 2 )
                            && tFacePosition( 1 ) < tProcRange( 3 ) && tFacePosition( 2 ) >= tProcRange( 4 )
                            && tFacePosition( 2 ) < tProcRange( 5 ) && aProcNeighbors( 3 ) < UINT_MAX )
                    {
                        tProcShare(tVar) = tProcRank;
                        tVar++;
                        tProcShare(tVar) = aProcNeighbors( 3 );
                        tVar++;
                    }
                    else if( tFacePosition( 0 ) == tProcRange( 1 ) && tFacePosition( 1 ) >= tProcRange( 2 )
                            && tFacePosition( 1 ) < tProcRange( 3 ) && tFacePosition( 2 ) >= tProcRange( 4 )
                            && tFacePosition( 2 ) < tProcRange( 5 ) && aProcNeighbors( 1 ) == UINT_MAX )
                    {
                        tProcShare(tVar) = tProcRank;
                        tVar++;
                    }
                    else if( tFacePosition( 0 ) == tProcRange( 1 ) && tFacePosition( 1 ) >= tProcRange( 2 )
                            && tFacePosition( 1 ) < tProcRange( 3 ) && tFacePosition( 2 ) >= tProcRange( 4 )
                            && tFacePosition( 2 ) < tProcRange( 5 ) && aProcNeighbors( 1 ) < UINT_MAX )
                    {
                        tProcShare(tVar) = tProcRank;
                        tVar++;
                        tProcShare(tVar) = aProcNeighbors( 1 );
                        tVar++;
                    }
                }
            }
            else if( tFacePosition( 3 ) < UINT_MAX && tFacePosition( 4 ) < UINT_MAX && tFacePosition( 5 ) < UINT_MAX  )
            {
                // Check if it is an edge in y-direction
                if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                        && tFacePosition( 4 ) > tProcRange( 2 ) && tFacePosition( 4 ) < tProcRange( 3 )
                        && tFacePosition( 5 ) >= tProcRange( 4 ) && tFacePosition( 5 ) < tProcRange( 5 ) )
                {
                    tProcShare(tVar) = tProcRank;
                    tVar++;
                }
                else
                {
                    if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                            && tFacePosition( 4 ) == tProcRange( 2 ) && tFacePosition( 5 ) >= tProcRange( 4 )
                            && tFacePosition( 5 ) < tProcRange( 5 ) && aProcNeighbors( 0 ) == UINT_MAX )
                    {
                        tProcShare(tVar) = tProcRank;
                        tVar++;
                    }
                    else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                            && tFacePosition( 4 ) == tProcRange( 2 ) && tFacePosition( 5 ) >= tProcRange( 4 )
                            && tFacePosition( 5 ) < tProcRange( 5 ) && aProcNeighbors( 0 ) < UINT_MAX )
                    {
                        tProcShare(tVar) = tProcRank;
                        tVar++;
                        tProcShare(tVar) = aProcNeighbors( 0 );
                        tVar++;
                    }
                    else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                            && tFacePosition( 4 ) == tProcRange( 3 ) && tFacePosition( 5 ) >= tProcRange( 4 )
                            && tFacePosition( 5 ) < tProcRange( 5 ) && aProcNeighbors( 2 ) == UINT_MAX )
                    {
                        tProcShare(tVar) = tProcRank;
                        tVar++;
                    }
                    else if( tFacePosition( 3 ) >= tProcRange( 0 ) && tFacePosition( 3 ) < tProcRange( 1 )
                            && tFacePosition( 4 ) == tProcRange( 3 ) && tFacePosition( 5 ) >= tProcRange( 4 )
                            && tFacePosition( 5 ) < tProcRange( 5 ) && aProcNeighbors( 2 ) < UINT_MAX )
                    {
                        tProcShare(tVar) = tProcRank;
                        tVar++;
                        tProcShare(tVar) = aProcNeighbors( 2 );
                        tVar++;
                    }
                }
            }
            else if( tFacePosition( 6 ) < UINT_MAX && tFacePosition( 7 ) < UINT_MAX && tFacePosition( 8 ) < UINT_MAX  )
            {
                // Check if it is an edge in z-direction
                if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                        && tFacePosition( 7 ) >= tProcRange( 2 ) && tFacePosition( 7 ) < tProcRange( 3 )
                        && tFacePosition( 8 ) > tProcRange( 4 ) && tFacePosition( 8 ) < tProcRange( 5 ) )
                {
                    tProcShare(tVar) = tProcRank;
                    tVar++;
                }
                else
                {
                    if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                            && tFacePosition( 7 ) >= tProcRange( 2 ) && tFacePosition( 7 ) < tProcRange( 3 )
                            && tFacePosition( 8 ) == tProcRange( 4 ) && aProcNeighbors( 4 ) == UINT_MAX )
                    {
                        tProcShare(tVar) = tProcRank;
                        tVar++;
                    }
                    else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                            && tFacePosition( 7 ) >= tProcRange( 2 ) && tFacePosition( 7 ) < tProcRange( 3 )
                            && tFacePosition( 8 ) == tProcRange( 4 ) && aProcNeighbors( 4 ) < UINT_MAX )
                    {
                        tProcShare(tVar) = tProcRank;
                        tVar++;
                        tProcShare(tVar) = aProcNeighbors( 4 );
                        tVar++;
                    }
                    else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                            && tFacePosition( 7 ) >= tProcRange( 2 ) && tFacePosition( 7 ) < tProcRange( 3 )
                            && tFacePosition( 8 ) == tProcRange( 5 ) && aProcNeighbors( 5 ) == UINT_MAX )
                    {
                        tProcShare(tVar) = tProcRank;
                        tVar++;
                    }
                    else if( tFacePosition( 6 ) >= tProcRange( 0 ) && tFacePosition( 6 ) < tProcRange( 1 )
                            && tFacePosition( 7 ) >= tProcRange( 2 ) && tFacePosition( 7 ) < tProcRange( 3 )
                            && tFacePosition( 8 ) == tProcRange( 5 ) && aProcNeighbors( 5 ) < UINT_MAX )
                    {
                        tProcShare(tVar) = tProcRank;
                        tVar++;
                        tProcShare(tVar) = aProcNeighbors( 5 );
                        tVar++;
                    }
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


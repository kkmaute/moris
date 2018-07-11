/*
 * cl_Lagrange_Basis.hpp
 *
 *  Created on: Feb 21, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_LAGRANGE_BASIS_HPP_
#define SRC_MESH_CL_LAGRANGE_BASIS_HPP_

#include "algorithms.hpp"
#include "cl_Base_Mesh_Element.hpp"
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_BoostBitset.hpp" // CON/src

namespace moris
{

    class Lagrange_Basis
    {
    protected:

    public:
        //Create Object of BaseElement
        Base_Mesh_Element mBaseElement;
        /**
         * Lagrange_Basis constructor
         */
        Lagrange_Basis()
    {
    }

        /**
         * Lagrange_Basis destructor.
         */
        ~Lagrange_Basis() = default;

        /**
         * Provides the number of basis functions within all levels from level 0 until aLevel
         *
         * @param[in] aLevel              Level of the basis functions.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] basis_number        Number of basis functions within all Levels until aLevel.
         *
         */
        /* static uint
        give_number_of_basis(
                uint const & aLevel,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection); */
        inline static uint
        give_number_of_basis(
                uint const & aLevel,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection)
        {
            // Variable to calculate the number of basis functions from level 0 until level "aLevel"
            uint tBasisNumber = 0;
            if( aModelDim == 1 )
            {
                for( uint i = 0; i < aLevel + 1; i++ )
                {
                    tBasisNumber += pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 ;
                }
            }
            else if( aModelDim == 2 )
            {
                for( uint i = 0; i < aLevel + 1; i++ )
                {
                    tBasisNumber +=  ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                                                                                                                    * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1 );
                }
            }
            else if( aModelDim == 3 )
            {
                for( uint i = 0; i < aLevel + 1 ; i++ )
                {
                    tBasisNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                                                                                                                    * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1 )
                                                                                                                    * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) * aPolynomial + 1 );
                }
            }
            return tBasisNumber;
        }
        /**
         * Provides the neighbours of an basis
         *
         * @param[in] aBasisId              Basis function Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial         Polynomial degree
         * @param[in] aBuffer             Provides the number of layers around the basis
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Basis neighbour        Basis Id plus neighbor Ids (9 for aModelDim = 2, 27 for aModelDim = 3).\n
         *                                2D:...----------------------------------...3D:....-----plain 1----------------------.....-----plain 2-------------------------\n.....-----plain 3-------------------------\n
         *                                ......|Element(6) Element(7) Element(8)|..........|Element(6) Element(7) Element(8)|.....|Element(15) Element(16) Element(17)|\n.....|Element(24) Element(25) Element(26)|\n
         *                                ......|Element(3) Element(4) Element(5)|..........|Element(3) Element(4) Element(5)|.....|Element(12) Element(13) Element(14)|\n.....|Element(21) Element(22) Element(23)|\n
         *                                ......|Element(0) Element(1) Element(2)|..........|Element(0) Element(1) Element(2)|.....|Element(9) Element(10) Element(11) |\n.....|Element(18) Element(19) Element(20)|\n
         *                                ......|--------------------------------|..........|--------------------------------|.....|-----------------------------------|  .....|-----------------------------------|\n
         *                                2D: Element(4) = aElement, 3D: Element(13) = aElement;
         *
         */
        Mat<uint>
        give_neighbor_of_basis(uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                uint const & aBuffer,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the level of the basis function
         *
         * @param[in] aBasisId            Basis function Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] level        Level of the element from the hierarchical mesh.
         *
         */
        /* uint
        give_basis_level(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection) const; */
        inline uint
        give_basis_level(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection) const
        {
            // output variable
            uint tBasisLevel = 1;
            // temporary variable for the while loop
            uint tLevel = 1;
            //Compute the relation of the different levels by the power of the level
            uint tPowLevel = pow( 2, tBasisLevel - 1 );
            if( aModelDim == 1 )
            {
                // temporary variable for the number of basis per level
                uint tNumberOfBasis = tPowLevel * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1;
                while(tLevel>0)
                {
                    // Check if it fits in the level
                    tLevel = floor( (real)( aBasisId + 1 ) / tNumberOfBasis );
                    if( (real)( aBasisId + 1 ) / tNumberOfBasis == 1 || tLevel < 1 )
                    {
                        tLevel = 0;
                    }
                    else
                    {
                        tBasisLevel++;
                        // Count the number of basis functions on the current level
                        tPowLevel = pow( 2, tBasisLevel - 1 );
                        tNumberOfBasis += tPowLevel * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1;
                    }
                }
                tBasisLevel--;
            }
            else if( aModelDim == 2 )
            {
                // temporary variable for the number of basis per level
                uint tNumberOfBasis = ( tPowLevel * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                                                                                                    * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1 );
                while( tLevel > 0 )
                {
                    // Check if it fits in the level
                    tLevel = floor( (real)( aBasisId + 1 ) / tNumberOfBasis );
                    if( (real)( aBasisId + 1 ) / tNumberOfBasis == 1 || tLevel < 1 )
                    {
                        tLevel = 0;
                    }
                    else
                    {
                        tBasisLevel++;
                        // Count the number of basis functions on the current level
                        tPowLevel = pow( 2, tBasisLevel - 1 );
                        tNumberOfBasis += ( tPowLevel * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                                                                                                        * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1 );
                    }
                }
                tBasisLevel--;
            }
            else if( aModelDim == 3 )
            {
                // temporary variable for the number of elements per level
                uint tNumberOfBasis = ( tPowLevel * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                                                                                                    * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1 )
                                                                                                    * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) * aPolynomial + 1 );
                while( tLevel > 0 )
                {
                    // Check if it fits in the level
                    tLevel = floor( (real)( aBasisId + 1 ) / tNumberOfBasis );
                    if( (real)( aBasisId + 1 ) / tNumberOfBasis == 1 || tLevel < 1 )
                    {
                        tLevel = 0;
                    }
                    else
                    {
                        tBasisLevel++;
                        // Count the number of basis functions on the current level
                        tPowLevel = pow( 2, tBasisLevel - 1 );
                        tNumberOfBasis += ( tPowLevel * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                                                                                                        * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1 )
                                                                                                        * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) * aPolynomial + 1 );
                    }
                }
                tBasisLevel--;
            }
            return tBasisLevel;
        }

        /**
         * Provides the position of a basis function from a tensorial grid
         *
         * @param[in] aBasisId            Basis function Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] Basis_position        Position I,J,K.
         *
         */
        /* Mat<uint>
        give_position_of_basis(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection) const; */
        inline Mat<uint>
        give_position_of_basis(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection) const
        {
            // Basis positions in the direction x,y,z
            Mat<uint> tBasisPosition( aModelDim , 1 , 0 );
            // Compute the level of the basis function
            uint tBasisLevel=give_basis_level( aBasisId, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
            // temporary variable for the number of elements per level
            uint tNumberOfBasis = 0;
            if( aModelDim == 2 )
            {
                if( tBasisLevel != 0 )
                {
                    //Count basis functions until tBasisLevel-1
                    tNumberOfBasis = give_number_of_basis(tBasisLevel-1, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
                }

                tBasisPosition( 1 ) = ceil( (real)(aBasisId + 1 - tNumberOfBasis )
                        / (real)( ( pow( 2, tBasisLevel ) * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 ) ) );
                // Calculate the x-direction with the calculated y-position
                tBasisPosition( 0 ) = aBasisId +1 - tNumberOfBasis
                        + ( pow( 2, tBasisLevel ) * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                        - tBasisPosition( 1 ) * ( pow( 2, tBasisLevel ) * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 );
                // -1 to get an indexed basis position, starting with zero
                tBasisPosition( 1 ) -= 1;
                tBasisPosition( 0 ) -= 1;
            }
            else if( aModelDim == 3 )
            {
                if( tBasisLevel !=0 )
                {
                    //Count basis functions until tBasisLevel-1
                    tNumberOfBasis = give_number_of_basis( tBasisLevel-1, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
                }
                //Temporary variables with precalculation for simpler coding style
                uint tPow2 = pow( 2, tBasisLevel );

                uint tConst0 = tPow2 * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1;
                uint tConst1 = tPow2 * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1;
                // Calculate the position in z-direction
                tBasisPosition( 2 ) = ceil( (real)( aBasisId + 1 - tNumberOfBasis ) / ( (real)( tConst0 * tConst1 ) ) );
                // Calculate the y-direction with the calculated z-position
                tBasisPosition( 1 ) = ceil( (real)( aBasisId + 1 - tNumberOfBasis - ( tBasisPosition( 2 ) - 1 ) * tConst0 * tConst1 ) / (real)tConst0 );
                // Calculate the x-direction with the calculated y and z-position
                tBasisPosition( 0 ) = aBasisId + 1 - tNumberOfBasis + tConst0 * ( 1 - tBasisPosition( 1 ) + tConst1 * ( 1 - tBasisPosition( 2 ) ) );
                // -1 to get an indexed basis position, starting with zero
                tBasisPosition( 2 ) -= 1;
                tBasisPosition( 1 ) -= 1;
                tBasisPosition( 0 ) -= 1;
            }
            return tBasisPosition;
        }

        /**
         * Provides the basis function Id with the position i,j,k
         *
         * @param[in] aLevel             Level of the basis function Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] tBasis       Basis function Id.
         *
         */
        /* uint
        give_basis_of_position(
                uint const & aLevel,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aIJKPosition) const; */
        inline uint
        give_basis_of_position(
                uint const & aLevel,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aIJKPosition) const
        {
            uint tBasisId = 0;
            // Position of the Element in y-direction
            uint tY_Position = 0;
            if( aModelDim == 1 )
            {
                if( aLevel > 0 )
                {
                    //Give number of basis functions for the last level to have the initial number at position 0,0
                    tBasisId = give_number_of_basis( aLevel-1, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
                }
                //Basis function with the number of basis of lower level
                tBasisId += aIJKPosition( 0 );
            }
            else if( aModelDim == 2 )
            {
                tY_Position = aIJKPosition( 1 )*( aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 );
                if( aLevel>0 )
                {
                    //Give number of basis functions for the last level to have the initial number at position 0,0
                    tBasisId = give_number_of_basis( aLevel-1, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
                    // Position of the Element in y-direction for a specific level
                    tY_Position = aIJKPosition( 1 ) * ( 1 + aPolynomial * aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel ) );
                }
                //Basis function with the calculated y-position
                tBasisId += aIJKPosition( 0 ) + tY_Position;
            }
            else if( aModelDim == 3 )
            {
                tY_Position = aIJKPosition( 1 ) * ( aPolynomial * aNumberOfElementsPerDirection( 0 ) + 1 )
                                                                                                   + aIJKPosition( 2 ) * ( aPolynomial * aNumberOfElementsPerDirection( 0 ) + 1 )
                                                                                                   * ( aPolynomial * aNumberOfElementsPerDirection( 1 ) + 1 );
                if( aLevel > 0 )
                {
                    //Give number of basis functions for the last level to have the initial number at position 0,0,0
                    tBasisId = give_number_of_basis( aLevel-1, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
                    // Position of the Element in y-direction for a specific level
                    tY_Position = aIJKPosition( 1 ) * ( aPolynomial * aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel ) +1  )
                                                                                                                    + aIJKPosition( 2 ) * ( aPolynomial * aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel ) + 1 )
                                                                                                                    * ( aPolynomial * aNumberOfElementsPerDirection( 1 ) * pow( 2, aLevel ) + 1 );
                }
                tBasisId += aIJKPosition( 0 ) + tY_Position; //Basis function with the calculated y-position
            }
            return tBasisId;
        }

        /**
         * Provides the basis function Id of the parent (Works only for a linear polynomial degree) (Returns UINT_MAX if there is no parent)
         *
         * @param[in] aBasisId            Basis function number.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] tBasis            Basis function Id of the parent.
         *
         */
        uint
        give_basis_of_parent(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the elements, which have support with the basis function from a tensorial grid
         *
         * @param[in] aBasisId            Basis function Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] tElements        Elements with support of basis function Id.
         *
         */
        Mat<uint>
        give_element_of_basis(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the coordinates for a specific basis function Id
         *
         * @param[in] aBasisId                Basis function Id
         * @param[in] aModelDim                  Number of dimensions.
         * @param[in] aPolynomial           Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection          NumElements=[Number of elements in x,y,z direction].
         * @param[in] aModelDimensions           Dimensions of the whole domain
         * @param[in] aModelDimensions_Offset    Offset of the problem, which is embedded in the outer layer of elements (it is needed to determine the point where the first basis functions sits, Point of origin (default = 0,0,0))
         *
         * @param[out] tCoordinates       Coordinates in x,y,z direction.
         *
         */
        Mat<real>
        give_coordinate_from_basis(uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<real> const & aModelDimensions,
                Mat<real> const & aModelDimensions_Offset) const;

        /**
         * Defines owner for all basis function Ids from the vector "aNodalLocaltoGlobal"
         *
         * @param[in] aModelDim                  Number of dimensions.
         * @param[in] aPolynomial           Polynomial degree of the basis functions.
         * @param[in] aPolynomial           Max Polynomial degree of Design and Lagrange basis (when in doubt, use same as aPolynomial)
         * @param[in] aNodalLocaltoGlobal   A list of node Ids, which are sitting on a proc (Needed for MTK/STK)
         * @param[in] aProcNeighbours       A list of proc neighbors
         * @param[in] aDecomp               Vector with a element range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
         * @param[in] aNumberBasis          Number of basis functions for the current level
         * @param[in] aBasisActive          A bitset with active basis functions
         *
         * @param[out] tNodeProcs           A vector with ownership of the node Ids
         *
         */
        Mat<uint>
        give_basis_proc_owner(
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aNodalLocaltoGlobal,
                Mat<uint> const & aProcNeighbour,
                Mat<uint> const & aDecomp,
                uint const & aNumberBasis,
                BoostBitset const & aBasisActive) const;

        /**
         * Defines owner of a basis function Id
         *
         * @param[in] aModelDim                  Basis function Id
         * @param[in] aModelDim                  Number of dimensions.
         * @param[in] aPolynomial           Polynomial degree of the basis functions.
         * @param[in] aProcNeighbours       A list of proc neighbors
         * @param[in] aDecomp               Vector with a element range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
         * @param[in] aNumberBasis          Number of basis functions for the current level
         *
         * @param[out] tProcOwner           Owner of the basis function id
         *
         */
        uint
        give_basis_proc_owner(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aProcNeighbour,
                Mat<uint> const & aDecomp) const;

        /**
         * Defines owner of a basis function Id
         *
         * @param[in] aBasisId                  Basis function Id
         * @param[in] aModelDim                  Number of dimensions.
         * @param[in] aPolynomial           Polynomial degree of the basis functions.
         * @param[in] aProcNeighbours       A list of proc neighbors
         * @param[in] aDecomp               Vector with a element range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
         * @param[in] aNumberBasis          Number of basis functions for the current level
         *
         * @param[out] tProcShare           Who shares this basis function
         *
         */
        Mat<uint>
        give_basis_share(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aProcNeighbour,
                Mat<uint> const & aDecomp) const;

        /**
         * Provides for a list of lagrange basis function Id the parent basis ids (if possible)
         *
         * @param[in] aBasisList         List of basis functions
         * @param[in] aModelDim                  Number of dimensions.
         * @param[in] aPolynomial           Polynomial degree of the basis functions.
         *
         * @param[out] tBasisList        List of updated basis functions
         *
         */
        Mat<uint>
        give_parents_of_basis(
                Mat<uint> const & aBasisList,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection) const;
    };
}

#endif /* SRC_MESH_CL_LAGRANGE_BASIS_HPP_ */

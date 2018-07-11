/*
 * cl_HIERARCHICAL_MESH_Main.cpp
 *
 *  Created on: Dec 7, 2017
 *      Author: gleim
 */

#include "cl_Hierarchical_Mesh_Main.hpp" // STK/src/Hierarchical
#include "cl_Hierarchical_Mesh_Solver_Input.hpp"

using namespace moris;

void
/* @TODO assert that number of elements are even when using MPI */
/* @TODO put initial element initialization into a private function */
/* @TODO assert that mMeshData.DeleteDomain is not demanded for 1D */
Hierarchical_Mesh_Main::create_mesh()
{
    moris::tic create_hmr;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // General Mesh Settings ( Max Polynomial, Number of Elements per Direction on whole mesh
    //                         including invisible elements )
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // FIXME : this should be moved to check_settings
    if ( mMeshData.ModelDim == 3 )
    {
        MORIS_ASSERT( mMeshData.NumberOfElementsPerDirection( 2 ) != 0,
                "Something is wrong in the settings, Dim cannot be 3 if number of elements in z-direction is zero" );
    }
    if ( mMeshData.ModelDim == 2 )
    {
        MORIS_ASSERT( mMeshData.NumberOfElementsPerDirection( 2 ) == 0 ,
                "Something is wrong in the settings, Dim cannot be smaller then 3 if number of elements in z-direction is non-zero" );
    }

    //Take Max polynomial of topology and design to create the coarsest mesh
    mMeshData.MaxPolynomial = ( mMeshData.Polynomial < mMeshData.PolynomialDesign ) ? mMeshData.PolynomialDesign : mMeshData.Polynomial;

    // Number of elements in x-direction (B-splines need a layer of +2*MaxPolynomial) (The number of elements should always be even for MPI)
    uint nex = mMeshData.NumberOfElementsPerDirection( 0 ) + 2 * mMeshData.MaxPolynomial;
    uint ney = 0; // Number of elements in y-direction (B-splines need a layer of +2*MaxPolynomial) (The number of elements should always be even for MPI)
    if ( mMeshData.ModelDim > 1 )
    {
        // Number of elements in y-direction (B-splines need a layer of +2*MaxPolynomial) (The number of elements should always be even for MPI)
        ney = mMeshData.NumberOfElementsPerDirection( 1 ) + 2 * mMeshData.MaxPolynomial;
    }
    uint nez = 0;
    if ( mMeshData.ModelDim == 3 )
    {
        // Number of elements in z-direction (In 2D, nez is zero!) (The number of elements should always be even for MPI)
        nez = mMeshData.NumberOfElementsPerDirection( 2 ) + 2 * mMeshData.MaxPolynomial;
    }
    // Select the number of elements in each direction, x,y,z
    mMeshData.NumberOfElementsPerDirection( 0 ) = nex;
    mMeshData.NumberOfElementsPerDirection( 1 ) = ney;
    mMeshData.NumberOfElementsPerDirection( 2 ) = nez;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Make sure that max possible element ID does not exceed UINT_MAX
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // FIXME clean up finding maximum refinement level routine.
    uint tMaxRefinementLevel = ( mSettings.MaxDesignLevelOfRefinement < mSettings.MaxLevelSetLevelOfRefinement )
                                 ? mSettings.MaxLevelSetLevelOfRefinement : mSettings.MaxDesignLevelOfRefinement;
    uint tMaxRefinementObject = mElementData.MaxLevelOfObjectRefinement.max();
    tMaxRefinementLevel = ( tMaxRefinementLevel < tMaxRefinementObject ) ? tMaxRefinementObject : tMaxRefinementLevel;
    uint tNumberElementsTestOld = 0;
    for ( uint i = 0; i < tMaxRefinementLevel + 1; i++ )
    {
        uint tNumberElementsTest = mBaseElement .give_number_of_elements(
                i,
                mMeshData.ModelDim,
                mMeshData.NumberOfElementsPerDirection );

        MORIS_ERROR( tNumberElementsTest < UINT_MAX && tNumberElementsTestOld < tNumberElementsTest, "The mesh exceeds more then MAX_UINT Elements" );
        tNumberElementsTestOld = tNumberElementsTest;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Calculate dimensions of whole mesh
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //Variables, that get calculated automatically from the variables above
    mMeshData.Dimensions_Orig = mMeshData.Dimensions;
    // Length in x-direction
    real lengthorig    = mMeshData.Dimensions( 0 );
    // Length in y-direction
    real heightorig    = mMeshData.Dimensions( 1 );
    // Thickness of the rectangular
    real thicknessorig = mMeshData.Dimensions( 2 );

    // Total length in x-direction with layer for the b-splines
    real length = lengthorig + lengthorig / ( nex - 2 * mMeshData.MaxPolynomial )
            * 2 * mMeshData.MaxPolynomial;

    real height = heightorig;
    if ( ney > 0 )
    {
        // Total length in y-direction with layer for the b-splines
        height = heightorig + heightorig / ( ney - 2 * mMeshData.MaxPolynomial )
                * 2 * mMeshData.MaxPolynomial;
    }

    real thickness = thicknessorig;
    if ( nez > 0 )
    {
        // Total length in z-direction with layer for the b-splines
        thickness = thicknessorig + thicknessorig / ( nez - 2 * mMeshData.MaxPolynomial )
                * 2 * mMeshData.MaxPolynomial;
    }

    //Save data in struc
    mMeshData.Dimensions( 0 ) = length;
    mMeshData.Dimensions( 1 ) = height;
    mMeshData.Dimensions( 2 ) = thickness;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Lengths of element on Level 0
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    mMeshData.ElementLength( 0 ) = length / (real) nex;
    if ( ney > 0 )
    {
        mMeshData.ElementLength( 1 ) = height / (real) ney;
    }
    if ( nez > 0 )
    {
        mMeshData.ElementLength( 2 ) = thickness / (real) nez;
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Calculate Offset of point zero
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //Offset for the coordinates, such that the active elements starts with the point of origin
    //Calculates an offset, such that the first element starts with the x-coordinate x=PointOfOrigin( 0 );
    real offset_x = -lengthorig / ( nex - 2 * mMeshData.MaxPolynomial )
            * mMeshData.MaxPolynomial + mSettings.PointOfOrigin( 0 );
    real offset_y = 0.0;
    if ( ney > 0 )
    {
        //Calculates an offset, such that the first element starts with the y-coordinate y=PointOfOrigin( 1 );
        offset_y = -heightorig / ( ney - 2 * mMeshData.MaxPolynomial )
                * mMeshData.MaxPolynomial + mSettings.PointOfOrigin( 1 );
    }
    real offset_z = 0.0;
    if ( nez > 0 )
    {
        //Calculates an offset, such that the first element starts with the z-coordinate z=PointOfOrigin( 2 );
        offset_z = -thicknessorig / ( nez - 2 * mMeshData.MaxPolynomial )
                * mMeshData.MaxPolynomial + mSettings.PointOfOrigin( 2 );
    }
    //Save data in struc
    mMeshData.DimensionsOffset.set_size(3,1);
    mMeshData.DimensionsOffset( 0 ) = offset_x;
    mMeshData.DimensionsOffset( 1 ) = offset_y;
    mMeshData.DimensionsOffset( 2 ) = offset_z;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Calculate Maximum possible ID of elements and initialize bitsets
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //Computes the number of elements
    mMeshData.NumberElements = mBaseElement.give_number_of_elements(
            mMeshData.Level,
            mMeshData.ModelDim,
            mMeshData.NumberOfElementsPerDirection );

    mElementData.ElementActive.resize( mMeshData.NumberElements );
    mElementData.ElementRefined.resize( mMeshData.NumberElements );
    // Some of the elements will be more then once in the list and will be kicked out by unique

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Calculate number of nodes ( or basis ) in each direction on first level
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // calculate maximum node ID on level 0
    uint tMaxBasisID =  give_number_of_basis( 0 )  - 1;

    mMeshData.NumberOfBasisPerDirection = give_position_of_basis( tMaxBasisID  );
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Create initial initialization pattern and remove elements affected by DeleteDomain
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if ( par_size() == 1 )
    {
        //Strategy in serial (for timing purposes), is to activate all elements and deactive the outer layer.
        //In MPI, only elements on each proc will be activated and then a broadcast to all processors
        //Activate the inner elements of the domain
        //
        // note that elements on the edge as well as initially deactivated domains are always marked as refined
        uint tElement = 0;
        //Temporary variable for loop
        uint tVar = 0;
        Mat<uint> tPosition( mMeshData.ModelDim, 1, 0 );
        mProcData.ElementListOnProc.set_size( mMeshData.NumberElements, 1 );
        mProcData.ElementRefinedListOnProc.set_size( 2* mMeshData.NumberElements, 1 );
        if ( mMeshData.ModelDim == 1 )
        {
            // Active all elements
            for ( uint k = 0; k < mMeshData.NumberElements; k++ )
            {
                mElementData.ElementActive.set( k );
                mProcData.ElementListOnProc( k ) = k;
            }
            //Deactivate the outer layer, which is needed for the b-splines
            for ( uint i = 0; i < mMeshData.MaxPolynomial; i++ )
            {
                //Position of the element in a tensorial grid
                tPosition( 0 ) = i;
                tElement = mBaseElement.give_element_of_position(
                        mMeshData.Level,
                        mMeshData.ModelDim,
                        mMeshData.NumberOfElementsPerDirection,
                        tPosition );

                mElementData.ElementActive.reset( tElement );
                mProcData.ElementListOnProc( tElement ) = UINT_MAX;
                mElementData.ElementRefined.set( tElement );
                mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                tVar++;
            }
            for ( uint i = nex - mMeshData.MaxPolynomial; i < nex; i++ )
            {
                // Position of the element in a tensorial grid
                tPosition( 0 ) = i;

                tElement = mBaseElement.give_element_of_position(
                        mMeshData.Level,
                        mMeshData.ModelDim,
                        mMeshData.NumberOfElementsPerDirection,
                        tPosition );

                mElementData.ElementActive.reset( tElement );
                mProcData.ElementListOnProc( tElement ) = UINT_MAX;
                mElementData.ElementRefined.set( tElement );
                mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                tVar++;
            }
        }
        if ( mMeshData.ModelDim == 2 )
        {
            // Active all elements
            for ( uint k = 0; k < mMeshData.NumberElements; k++ )
            {
                mElementData.ElementActive.set( k );
                mProcData.ElementListOnProc( k ) = k;
            }
            //Deactivate the outer layer, which is needed for the B-Splines
            for ( uint j = 0; j < ney; j++ )
            {
                for ( uint i = 0; i < mMeshData.MaxPolynomial; i++ )
                {
                    // Position of the element in a tensor grid
                    tPosition( 0 ) = i;
                    tPosition( 1 ) = j;

                    tElement = mBaseElement.give_element_of_position(
                            mMeshData.Level,
                            mMeshData.ModelDim,
                            mMeshData.NumberOfElementsPerDirection,
                            tPosition );

                    mElementData.ElementActive.reset( tElement );
                    mProcData.ElementListOnProc( tElement ) = UINT_MAX;
                    mElementData.ElementRefined.set( tElement );
                    mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                    tVar++;
                }
                for ( uint i = nex - mMeshData.MaxPolynomial; i < nex; i++ )
                {
                    // Position of the element in a tensor grid
                    tPosition( 0 ) = i;
                    tPosition( 1 ) = j;
                    tElement = mBaseElement.give_element_of_position(
                            mMeshData.Level,
                            mMeshData.ModelDim,
                            mMeshData.NumberOfElementsPerDirection,
                            tPosition );

                    mElementData.ElementActive.reset( tElement );
                    mProcData.ElementListOnProc( tElement ) = UINT_MAX;
                    mElementData.ElementRefined.set( tElement );
                    mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                    tVar++;
                }
            }
            for ( uint i = 0; i < nex; i++ )
            {
                for ( uint j = 0; j < mMeshData.MaxPolynomial; j++ )
                {
                    // Position of the element in a tensor grid
                    tPosition( 0 ) = i;
                    tPosition( 1 ) = j;

                    tElement = mBaseElement.give_element_of_position(
                            mMeshData.Level,
                            mMeshData.ModelDim,
                            mMeshData.NumberOfElementsPerDirection,
                            tPosition );

                    (mElementData.ElementActive).reset(tElement);
                    mProcData.ElementListOnProc(tElement) = UINT_MAX;
                    mElementData.ElementRefined.set( tElement );
                    mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                    tVar++;
                }
                for ( uint j = ney-mMeshData.MaxPolynomial; j<ney; j++)
                {
                    // Position of the element in a tensor grid
                    tPosition( 0 ) = i;
                    tPosition( 1 ) = j;

                    tElement = mBaseElement.give_element_of_position(
                            mMeshData.Level,
                            mMeshData.ModelDim,
                            mMeshData.NumberOfElementsPerDirection,
                            tPosition );

                    mElementData.ElementActive.reset( tElement );
                    mProcData.ElementListOnProc( tElement ) = UINT_MAX;
                    mElementData.ElementRefined.set( tElement );
                    mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                    tVar++;
                }
            }
            // If an initial domain needs to be deactivated
            Mat<real> DeleteDomain = mMeshData.DeleteDomain;
            Mat<real> tMatDelDomain0 = DeleteDomain.row( 0 );
            real tDelDomain0 = sum( tMatDelDomain0 );
            Mat<real> tMatDelDomain1 = DeleteDomain.row( 1 );
            real tDelDomain1 = sum( tMatDelDomain1 );
            if (  tDelDomain0 != 0.0 ||  tDelDomain1 != 0.0 )
            {
                Mat<uint> tPossibleDeleteElements = unique( ( mProcData.ElementListOnProc ) );
                tPossibleDeleteElements.resize( tPossibleDeleteElements.length() - 1, 1 );
                Mat<uint> tBasis;
                Mat<real> tElementMiddleCoord;
                for ( uint k = 0; k < tPossibleDeleteElements.length(); k++ )
                {
                    tBasis = mHMRElement.give_basis_of_element(
                            tPossibleDeleteElements( k ),
                            mMeshData.ModelDim,
                            mMeshData.Polynomial,
                            mMeshData.NumberOfElementsPerDirection );

                    tElementMiddleCoord = mHMRElement.give_middlecoordinate_from_element(
                            tPossibleDeleteElements( k),
                            mMeshData.ModelDim,
                            mMeshData.NumberOfElementsPerDirection,
                            mMeshData.Dimensions,
                            mMeshData.DimensionsOffset );

                    if ( tElementMiddleCoord( 0 ) > DeleteDomain( 0, 0 ) &&  tElementMiddleCoord( 0 )  < DeleteDomain( 1, 0 ) &&
                            tElementMiddleCoord( 1 ) > DeleteDomain( 0, 1 ) &&  tElementMiddleCoord( 1 )  < DeleteDomain( 1, 1 ) )
                    {
                        mElementData.ElementActive.reset( tPossibleDeleteElements( k ) );
                        mProcData.ElementListOnProc( tPossibleDeleteElements( k ) ) = UINT_MAX;
                        mElementData.ElementRefined.set( tElement );
                        mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                        tVar++;
                    }
                }
            }
        }
        else if ( mMeshData.ModelDim == 3 )
        {
            // Active all elements
            for ( uint k = 0; k < mMeshData.NumberElements; k++ )
            {
                mElementData.ElementActive.set( k );
                mProcData.ElementListOnProc( k ) = k;
            }
            //Deactivate the outer layer, which is needed for the B-Splines
            for ( uint k = 0; k < nez; k++ )
            {
                for ( uint j = 0; j < ney; j++ )
                {
                    for ( uint i = 0; i < mMeshData.MaxPolynomial; i++ )
                    {
                        // Position of the element in a tensorial grid
                        tPosition( 0 ) = i;
                        tPosition( 1 ) = j;
                        tPosition( 2 ) = k;

                        tElement = mBaseElement.give_element_of_position(
                                mMeshData.Level,
                                mMeshData.ModelDim,
                                mMeshData.NumberOfElementsPerDirection,
                                tPosition );

                        mElementData.ElementActive.reset( tElement );
                        mProcData.ElementListOnProc( tElement ) = UINT_MAX;
                        mElementData.ElementRefined.set( tElement );
                        mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                        tVar++;
                    }
                    for ( uint i = 1; i <= mMeshData.MaxPolynomial; i++ )
                    {
                        // Position of the element in a tensorial grid
                        tPosition( 0 ) = nex - i;
                        tPosition( 1 ) = j;
                        tPosition( 2 ) = k;

                        tElement = mBaseElement.give_element_of_position(
                                mMeshData.Level,
                                mMeshData.ModelDim,
                                mMeshData.NumberOfElementsPerDirection,
                                tPosition );

                        mElementData.ElementActive.reset( tElement );
                        mProcData.ElementListOnProc( tElement ) = UINT_MAX;
                        mElementData.ElementRefined.set( tElement );
                        mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                        tVar++;
                    }
                }
            }
            for ( uint i = 0; i < nex; i++ )
            {
                for ( uint j = 0; j < ney; j++ )
                {
                    for ( uint k = 0; k < mMeshData.MaxPolynomial; k++ )
                    {
                        // Position of the element in a tensorial grid
                        tPosition( 0 ) = i;
                        tPosition( 1 ) = j;
                        tPosition( 2 ) = k;

                        tElement = mBaseElement.give_element_of_position(
                                mMeshData.Level,
                                mMeshData.ModelDim,
                                mMeshData.NumberOfElementsPerDirection,
                                tPosition );

                        mElementData.ElementActive.reset( tElement );
                        mProcData.ElementListOnProc( tElement ) = UINT_MAX;
                        mElementData.ElementRefined.set( tElement );
                        mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                        tVar++;
                    }
                    for ( uint k = 1; k <= mMeshData.MaxPolynomial; k++ )
                    {
                        // Position of the element in a tensorial grid
                        tPosition( 0 ) = i;
                        tPosition( 1 ) = j;
                        tPosition( 2 ) = nez - k;

                        tElement = mBaseElement.give_element_of_position(
                                mMeshData.Level,
                                mMeshData.ModelDim,
                                mMeshData.NumberOfElementsPerDirection,
                                tPosition );

                        mElementData.ElementActive.reset( tElement );
                        mProcData.ElementListOnProc( tElement ) = UINT_MAX;
                        mElementData.ElementRefined.set( tElement );
                        mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                        tVar++;
                    }
                }
            }
            for ( uint i = 0; i < nex; i++ )
            {
                for ( uint k = 0; k < nez; k++ )
                {
                    for ( uint j = 0; j < mMeshData.MaxPolynomial; j++ )
                    {
                        // Position of the element in a tensorial grid
                        tPosition( 0 ) = i;
                        tPosition( 1 ) = j;
                        tPosition( 2 ) = k;

                        tElement = mBaseElement.give_element_of_position(
                                mMeshData.Level,
                                mMeshData.ModelDim,
                                mMeshData.NumberOfElementsPerDirection,
                                tPosition );

                        mElementData.ElementActive.reset( tElement );
                        mProcData.ElementListOnProc( tElement ) = UINT_MAX;
                        mElementData.ElementRefined.set( tElement );
                        mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                        tVar++;
                    }
                    for ( uint j = 1; j <= mMeshData.MaxPolynomial; j++ )
                    {
                        // Position of the element in a tensorial grid
                        tPosition( 0 ) = i;
                        tPosition( 1 ) = ney - j;
                        tPosition( 2 ) = k;

                        tElement = mBaseElement.give_element_of_position(
                                mMeshData.Level,
                                mMeshData.ModelDim,
                                mMeshData.NumberOfElementsPerDirection,
                                tPosition );

                        mElementData.ElementActive.reset( tElement );
                        mProcData.ElementListOnProc( tElement ) = UINT_MAX;
                        mElementData.ElementRefined.set( tElement );
                        mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                        tVar++;
                    }
                }
            }
            // If an initial domain needs to be deactivated
            Mat<real> DeleteDomain = mMeshData.DeleteDomain;
            Mat<real> tMatDelDomain0 = DeleteDomain.row( 0 );
            real tDelDomain0 = sum( tMatDelDomain0 );
            Mat<real> tMatDelDomain1 = DeleteDomain.row( 1 );
            real tDelDomain1 = sum( tMatDelDomain1 );
            if ( tDelDomain0 != 0.0 || tDelDomain1 != 0.0 )
            {
                Mat<uint> tPossibleDeleteElements = unique( mProcData.ElementListOnProc );
                tPossibleDeleteElements.resize( tPossibleDeleteElements.length() - 1, 1 );
                Mat<uint> tBasis;
                Mat<real> tElementMiddleCoord;
                for ( uint k = 0; k < tPossibleDeleteElements.length(); k++ )
                {
                    //Compute basis functions, which have support in this element
                    tBasis = mHMRElement.give_basis_of_element( tPossibleDeleteElements( k ), mMeshData.ModelDim, mMeshData.Polynomial, mMeshData.NumberOfElementsPerDirection );
                    //Compute a coordinate in the middle of an element to know if this in the domain, which we want to deactivate
                    tElementMiddleCoord = mHMRElement.give_middlecoordinate_from_element( tPossibleDeleteElements( k ), mMeshData.ModelDim, mMeshData.NumberOfElementsPerDirection, mMeshData.Dimensions, mMeshData.DimensionsOffset );
                    if ( tElementMiddleCoord( 0 ) > DeleteDomain( 0, 0 ) &&  tElementMiddleCoord( 0 )  < DeleteDomain( 1, 0 ) &&
                            tElementMiddleCoord( 1 ) > DeleteDomain( 0, 1 ) &&  tElementMiddleCoord( 1 )  < DeleteDomain( 1, 1 ) &&
                            tElementMiddleCoord( 2 ) > DeleteDomain( 0, 2 ) &&  tElementMiddleCoord( 2 )  < DeleteDomain( 1, 2 ) )
                    {
                        mElementData.ElementActive.reset( tPossibleDeleteElements( k ) );
                        mProcData.ElementListOnProc( tPossibleDeleteElements( k ) ) = UINT_MAX;
                        mElementData.ElementRefined.set( tElement );
                        mProcData.ElementRefinedListOnProc( tVar ) = tElement;
                        tVar++;
                    }
                }
            }
        }
        mProcData.ElementListOnProc = unique( mProcData.ElementListOnProc );
        mProcData.ElementListOnProc.resize( mProcData.ElementListOnProc.length() - 1, 1 );
        mProcData.ElementRefinedListOnProc.resize( tVar, 1 );
        mProcData.ElementRefinedListOnProc = unique( mProcData.ElementRefinedListOnProc );
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Mesh Decomposition for parallel processing
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //If MPI, then a decomposition needs to be done (Decomposition is made by MPI itself with MPI_Coord)
    this->mesh_decomposition();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Parent-Child relation T-Matrices
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //Compute all necessary T-matrices with their relationship between an element and their children
    mElementData.TMatrixParentChildRelationFEM =
            mTMatrix.give_Tmatrices_of_childs(
                    mMeshData.Polynomial,
                    mMeshData.ModelDim );

    mElementData.TMatrixParentChildRelationDesign = mTMatrix.give_Tmatrices_of_childs(
            mMeshData.PolynomialDesign,
            mMeshData.ModelDim );

    /*mElementData.TMatrixParentChildRelationDesignProjection =
            mTMatrix.give_Tmatrices_of_childs(
                    mMeshData.PolynomialDesign,
                    mMeshData.ModelDim ); */

    mElementData.TProjectionMatrixBasis = give_projection_matrix_basis();
    mElementData.TProjectionMatrixDesignBasis = give_projection_matrix_design_basis();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Init Sidesets
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Set the side set and broadcast information for MPI
    this->set_sideset( mMeshData.SideSet );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Check for imposed symmetry conditions
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    this->check_and_init_symmetry();

    real tElapsed_create_hmr = create_hmr.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest )
    {
        std::fprintf(stdout,"Time for creating the mesh : %f [sec]\n",tElapsed_create_hmr);
    }
}

//----------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::mesh_decomposition()
{
    //Temporary variable for loop
    uint tVar = 0;
    // Calculate the decomposition for each proc with 6 numbers ( X_Start, X_End, Y_Start, Y_End, Z_Start, Z_End)
    mProcData.DecompElement.set_size( 6, 1, 0 );
    // Calculate the decomposition for each proc (with an additional layer of 1 element in each direction) with 6 numbers ( X_Start - 1, X_End + 1, Y_Start - 1, Y_End + 1, Z_Start - 1, Z_End + 1)
    mProcData.DecompElementAura.set_size( 6, 1, 0 );

    mProcData.DecompLagrangeNode.set_size( 6, 1 );

    if ( par_size() == 1 )
    {
        // Initial active elements need to be saved
        mProcData.ElementListOnProcInit = mProcData.ElementListOnProc;
        // Active elements with an auro need to be saved (Aura is in serial not there)
        mProcData.ElementListOnProcWithAura = mProcData.ElementListOnProc;
        //Decomposition in x-direction
        mProcData.DecompElement( 0 ) = mMeshData.MaxPolynomial;
        mProcData.DecompElement( 1 ) = mMeshData.NumberOfElementsPerDirection( 0 ) - ( mMeshData.MaxPolynomial + 1 ) ;
        if ( mMeshData.ModelDim > 1 )
        {
            //Decomposition in y-direction
            mProcData.DecompElement( 2 ) = mMeshData.MaxPolynomial;
            mProcData.DecompElement( 3 ) = mMeshData.NumberOfElementsPerDirection( 1 ) - ( mMeshData.MaxPolynomial + 1 );
        }
        if ( mMeshData.ModelDim == 3 )
        {
            //Decomposition in z-direction
            mProcData.DecompElement( 4 ) = mMeshData.MaxPolynomial;
            mProcData.DecompElement( 5 ) = mMeshData.NumberOfElementsPerDirection( 2 ) - ( mMeshData.MaxPolynomial + 1 );
        }
        /* @TODO check if this is correct */
        mProcData.DecompLagrangeNode(0) = mMeshData.MaxPolynomial;
        mProcData.DecompLagrangeNode(1) = mMeshData.MaxPolynomial
                                                  *(1+mMeshData.NumberOfElementsPerDirection(0));
        mProcData.DecompLagrangeNode(2) = mMeshData.MaxPolynomial;
        mProcData.DecompLagrangeNode(3) = mMeshData.MaxPolynomial
                                                  *(1+mMeshData.NumberOfElementsPerDirection(1));
        mProcData.DecompLagrangeNode(4) = mMeshData.MaxPolynomial;
        mProcData.DecompLagrangeNode(5) = mMeshData.MaxPolynomial
                                                  *(1+mMeshData.NumberOfElementsPerDirection(2));
    }
    else if ( par_size() > 1 )
    {
        // Number of processors in x, y and z-direction
        int tProcDimensions[3];
        tProcDimensions[0] = 0;
        tProcDimensions[1] = 0;
        tProcDimensions[2] = 0;
        // Creates a grid of processors
        MPI_Dims_create( par_size(), mMeshData.ModelDim, tProcDimensions );
        Mat<uint> tDimsCartProc( 3, 1, 0 );
        tDimsCartProc( 0, 0 ) = tProcDimensions[0];
        tDimsCartProc( 1, 0 ) = tProcDimensions[1];
        // Dummy is needed in 2D for the if clauses in proc_coordinates_neighbors
        tDimsCartProc( 2, 0 ) = 1;
        if ( mMeshData.ModelDim == 3 )
        {
            tDimsCartProc( 2, 0 ) = tProcDimensions[2];
        }
        // Calculate the coordinates of each proc and the Neighbors of each proc
        this->proc_coordinates_neighbors( tDimsCartProc );
        if ( ( mMeshData.NumberOfElementsPerDirection( 0 ) - 2 * mMeshData.MaxPolynomial ) % tDimsCartProc( 0, 0 ) != 0 )
        {
            MORIS_LOG_ERROR << "Number of Elements in x-direction must be a multiple of " << tDimsCartProc( 0, 0 )
                                                                                                    << ". Number of Elements in x-direciton = " << ( mMeshData.NumberOfElementsPerDirection( 0 ) - 2 * mMeshData.MaxPolynomial );
            exit( 1 );
        }
        if ( ( mMeshData.NumberOfElementsPerDirection( 1 ) - 2 * mMeshData.MaxPolynomial ) % tDimsCartProc( 1, 0 ) != 0 )
        {
            MORIS_LOG_ERROR << "Number of Elements in y-direction must be a multiple of " << tDimsCartProc( 1, 0 )
                                                                                                    << ". Number of Elements in y-direciton = " << ( mMeshData.NumberOfElementsPerDirection( 1 ) - 2 * mMeshData.MaxPolynomial );
            exit( 1 );
        }
        if ( mMeshData.ModelDim == 3 && ( mMeshData.NumberOfElementsPerDirection( 2 ) - 2 * mMeshData.MaxPolynomial ) % tDimsCartProc( 2, 0 ) != 0 )
        {
            MORIS_LOG_ERROR << "Number of Elements in z-direction must be a multiple of " << tDimsCartProc( 2, 0 )
                                                                                                    << ". Number of Elements in z-direciton = " << ( mMeshData.NumberOfElementsPerDirection( 2 ) - 2 * mMeshData.MaxPolynomial );
            exit( 1 );
        }
        //Compute decomposition in x-direction
        mProcData.DecompElement( 0 ) = mMeshData.MaxPolynomial + ( mMeshData.NumberOfElementsPerDirection( 0 ) - 2 * mMeshData.MaxPolynomial )
                                                                                                / tDimsCartProc( 0, 0 ) * mProcData.ProcCoord( 0 );
        mProcData.DecompElement( 1 ) = mMeshData.MaxPolynomial + ( mMeshData.NumberOfElementsPerDirection( 0 ) - 2 * mMeshData.MaxPolynomial )
                                                                                                        / tDimsCartProc( 0, 0 ) * ( mProcData.ProcCoord( 0 ) + 1 ) - 1;
        if ( mMeshData.ModelDim > 1 )
        {
            //Compute decomposition in y-direction
            mProcData.DecompElement( 2 ) = mMeshData.MaxPolynomial + ( mMeshData.NumberOfElementsPerDirection( 1 ) - 2 * mMeshData.MaxPolynomial )
                                                                                                    / tDimsCartProc( 1, 0 ) * mProcData.ProcCoord( 1 );
            mProcData.DecompElement( 3 ) = mMeshData.MaxPolynomial +  ( mMeshData.NumberOfElementsPerDirection( 1 ) - 2 * mMeshData.MaxPolynomial )
                                                                                                    / tDimsCartProc( 1, 0 ) * ( mProcData.ProcCoord( 1 ) + 1 ) - 1;
        }
        if ( mMeshData.ModelDim == 3 )
        {
            //Compute decomposition in z-direction
            mProcData.DecompElement( 4 ) = mMeshData.MaxPolynomial + ( mMeshData.NumberOfElementsPerDirection( 2 ) - 2 * mMeshData.MaxPolynomial )
                                                                                                    / tDimsCartProc( 2, 0 ) * mProcData.ProcCoord( 2 );
            mProcData.DecompElement( 5 ) =  mMeshData.MaxPolynomial + ( mMeshData.NumberOfElementsPerDirection( 2 ) - 2 * mMeshData.MaxPolynomial )
                                                                                                    / tDimsCartProc( 2, 0 ) * ( mProcData.ProcCoord( 2 ) + 1 ) - 1;
        }
        //Determine decomposition with aura elements
        mProcData.DecompElementAura( 0 ) = mProcData.DecompElement( 0 );
        mProcData.DecompElementAura( 1 ) = mProcData.DecompElement( 1 );
        mProcData.DecompElementAura( 2 ) = mProcData.DecompElement( 2 );
        mProcData.DecompElementAura( 3 ) = mProcData.DecompElement( 3 );
        mProcData.DecompElementAura( 4 ) = mProcData.DecompElement( 4 );
        mProcData.DecompElementAura( 5 ) = mProcData.DecompElement( 5 );

        if ( mProcData.ProcNeighbor( 1 ) != UINT_MAX )
        {
            mProcData.DecompElementAura( 1 ) += 1;
        }
        if ( mProcData.ProcNeighbor( 3 ) != UINT_MAX )
        {
            mProcData.DecompElementAura( 0 ) -=1;
        }
        if ( mProcData.ProcNeighbor( 0 ) != UINT_MAX )
        {
            mProcData.DecompElementAura( 2 ) -=1;
        }
        if ( mProcData.ProcNeighbor( 2 ) != UINT_MAX )
        {
            mProcData.DecompElementAura( 3 ) +=1;
        }
        if ( mMeshData.ModelDim == 3)
        {
            if ( mProcData.ProcNeighbor( 4 ) != UINT_MAX )
            {
                mProcData.DecompElementAura( 4 ) -=1;
            }
            if ( mProcData.ProcNeighbor( 5 ) != UINT_MAX )
            {
                mProcData.DecompElementAura( 5 ) +=1;
            }
        }

        //Create element list with aura on processor domain
        tVar = 0;
        if ( mMeshData.ModelDim == 2 )
        {
            mProcData.ElementListOnProcWithAura.set_size( ( mProcData.DecompElementAura( 1 ) - mProcData.DecompElementAura( 0 ) + 1 )
                    * ( mProcData.DecompElementAura( 3 ) - mProcData.DecompElementAura( 2 ) + 1 ), 1 );
            for ( uint j = mProcData.DecompElementAura( 2 ); j <= mProcData.DecompElementAura( 3 ); j++ )
            {
                for ( uint i = mProcData.DecompElementAura( 0 ); i <= mProcData.DecompElementAura( 1 ); i++ )
                {
                    if ( mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) / 2
                            + mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) * ( i - 1 ) > mMeshData.DeleteDomain( 0, 0 )
                            && mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) / 2
                            + mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) * ( i - 1 ) < mMeshData.DeleteDomain( 1, 0 )
                            && mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) / 2
                            + mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) * ( j - 1 ) > mMeshData.DeleteDomain( 0, 1 )
                            && mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) / 2
                            + mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) * ( j - 1 ) < mMeshData.DeleteDomain( 1, 1 ) )
                    {
                    }
                    else
                    {
                        mProcData.ElementListOnProcWithAura(tVar) = i + j * mMeshData.NumberOfElementsPerDirection( 0 );
                        tVar++;
                    }
                }
            }
        }

        else if ( mMeshData.ModelDim == 3 )
        {
            mProcData.ElementListOnProcWithAura.set_size( ( mProcData.DecompElementAura( 1 ) - mProcData.DecompElementAura( 0 ) + 1 )
                    * ( mProcData.DecompElementAura( 3 ) - mProcData.DecompElementAura( 2 ) + 1 )
                    * ( mProcData.DecompElementAura( 5 ) - mProcData.DecompElementAura( 4 ) + 1 ), 1 );
            for ( uint k = mProcData.DecompElementAura( 4 ); k <= mProcData.DecompElementAura( 5 ); k++ )
            {
                for ( uint j = mProcData.DecompElementAura( 2 ); j <= mProcData.DecompElementAura( 3 ); j++ )
                {
                    for ( uint i = mProcData.DecompElementAura( 0 ); i <= mProcData.DecompElementAura( 1 ); i++ )
                    {
                        if ( mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) / 2
                                + mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) * ( i - 1 ) > mMeshData.DeleteDomain( 0, 0 )
                                && mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) / 2
                                + mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) * ( i - 1 ) < mMeshData.DeleteDomain( 1, 0 )
                                && mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) / 2
                                + mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) * ( j - 1 ) > mMeshData.DeleteDomain( 0, 1 )
                                && mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) / 2
                                + mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) * ( j - 1 ) < mMeshData.DeleteDomain( 1, 1 )
                                && mMeshData.Dimensions( 2 ) / mMeshData.NumberOfElementsPerDirection( 2 ) / 2
                                + mMeshData.Dimensions( 2 ) / mMeshData.NumberOfElementsPerDirection( 2 ) * ( k - 1 ) > mMeshData.DeleteDomain( 0, 1 )
                                && mMeshData.Dimensions( 2 ) / mMeshData.NumberOfElementsPerDirection( 2 ) / 2
                                + mMeshData.Dimensions( 2 ) / mMeshData.NumberOfElementsPerDirection( 2 ) * ( k - 1 ) < mMeshData.DeleteDomain( 1, 2 ) )
                        {
                        }
                        else
                        {
                            mProcData.ElementListOnProcWithAura( tVar ) = i + j * mMeshData.NumberOfElementsPerDirection( 0 )
                                                                                                            + k * mMeshData.NumberOfElementsPerDirection( 0 ) * mMeshData.NumberOfElementsPerDirection( 1 );
                            tVar++;
                        }
                    }
                }
            }
        }

        mProcData.ElementListOnProcWithAura.resize( tVar, 1 );
        tVar = 0;
        //Create element list and activate bitset on processor domain
        if ( mMeshData.ModelDim == 2 )
        {
            mProcData.ElementListOnProc.set_size( ( mProcData.DecompElement( 1 ) - mProcData.DecompElement( 0 ) + 1 )
                    * ( mProcData.DecompElement( 3 ) - mProcData.DecompElement( 2 ) + 1 ), 1 );
            mProcData.ElementListOnProcInit.set_size( ( mProcData.DecompElement( 1 ) - mProcData.DecompElement( 0 ) + 1 )
                    * ( mProcData.DecompElement( 3 ) - mProcData.DecompElement( 2 ) + 1 ), 1 );
            for ( uint j = mProcData.DecompElement( 2 ); j <= mProcData.DecompElement( 3 ); j++ )
            {
                for ( uint i = mProcData.DecompElement( 0 ); i <= mProcData.DecompElement( 1 ); i++ )
                {
                    if ( mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) / 2
                            + mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) * ( i - 1 ) > mMeshData.DeleteDomain( 0, 0 )
                            && mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) / 2
                            + mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) * ( i - 1 ) < mMeshData.DeleteDomain( 1, 0 )
                            && mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) / 2
                            + mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) * ( j - 1 ) > mMeshData.DeleteDomain( 0, 1 )
                            && mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) / 2
                            + mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) * ( j - 1 ) < mMeshData.DeleteDomain( 1, 1 ) )
                    {
                    }
                    else
                    {
                        mProcData.ElementListOnProc( tVar ) = i + j * mMeshData.NumberOfElementsPerDirection( 0 );
                        //Activate the inner elements
                        mElementData.ElementActive.set( i + j * mMeshData.NumberOfElementsPerDirection( 0 ) );
                        tVar++;
                    }
                }
            }
        }
        else if ( mMeshData.ModelDim == 3 )
        {
            mProcData.ElementListOnProc.set_size( ( mProcData.DecompElement( 1 ) - mProcData.DecompElement( 0 ) + 1 )
                    * ( mProcData.DecompElement( 3 ) - mProcData.DecompElement( 2 ) + 1 )
                    * ( mProcData.DecompElement( 5 ) - mProcData.DecompElement( 4 ) + 1 ), 1 );
            mProcData.ElementListOnProcInit.set_size( ( mProcData.DecompElement( 1 ) - mProcData.DecompElement( 0 ) + 1 )
                    * ( mProcData.DecompElement( 3 ) - mProcData.DecompElement( 2 ) + 1 )
                    * ( mProcData.DecompElement( 5 ) - mProcData.DecompElement( 4 ) + 1 ), 1 );
            for ( uint k = mProcData.DecompElement( 4 ); k <= mProcData.DecompElement( 5 ); k++ )
            {
                for ( uint j = mProcData.DecompElement( 2 ); j <= mProcData.DecompElement( 3 ); j++ )
                {
                    for ( uint i = mProcData.DecompElement( 0 ); i <= mProcData.DecompElement( 1 ); i++ )
                    {
                        if ( mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) / 2
                                + mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) * ( i - 1 ) > mMeshData.DeleteDomain( 0, 0 )
                                && mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) / 2
                                + mMeshData.Dimensions( 0 ) / mMeshData.NumberOfElementsPerDirection( 0 ) * ( i - 1 ) < mMeshData.DeleteDomain( 1, 0 )
                                && mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) / 2
                                + mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) * ( j - 1 ) > mMeshData.DeleteDomain( 0, 1 )
                                && mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) / 2
                                + mMeshData.Dimensions( 1 ) / mMeshData.NumberOfElementsPerDirection( 1 ) * ( j - 1 ) < mMeshData.DeleteDomain( 1, 1 )
                                && mMeshData.Dimensions( 2 ) / mMeshData.NumberOfElementsPerDirection( 2 ) / 2
                                + mMeshData.Dimensions( 2 ) / mMeshData.NumberOfElementsPerDirection( 2 ) * ( k - 1 ) > mMeshData.DeleteDomain( 0, 1 )
                                && mMeshData.Dimensions( 2 ) / mMeshData.NumberOfElementsPerDirection( 2 ) / 2
                                + mMeshData.Dimensions( 2 ) / mMeshData.NumberOfElementsPerDirection( 2 ) * ( k -  1) < mMeshData.DeleteDomain( 1, 2 ) )
                        {
                        }
                        else
                        {
                            mProcData.ElementListOnProc( tVar ) = i + j * mMeshData.NumberOfElementsPerDirection( 0 )
                                                                                                            + k *mMeshData.NumberOfElementsPerDirection( 0 ) * mMeshData.NumberOfElementsPerDirection( 1 );
                            //Activate the inner elements
                            mElementData.ElementActive.set( i + j * mMeshData.NumberOfElementsPerDirection( 0 )
                                    + k * mMeshData.NumberOfElementsPerDirection( 0 ) * mMeshData.NumberOfElementsPerDirection( 1 ) );
                            tVar++;
                        }
                    }
                }
            }
        }
        mProcData.ElementListOnProc.resize( tVar, 1 );
        //Broadcast active elements, that each processor is aware of all active elements
        mMPI.broadcast_bitset_logical_or(mElementData.ElementActive);
        // Needed to subtract aura elements
        BoostBitset tDummyBitset = mElementData.ElementActive;
        for ( uint i = 0; i < mProcData.ElementListOnProc.length(); i++ )
        {
            // Deactivate elements from own proc to find aura elements
            tDummyBitset.reset( mProcData.ElementListOnProc( i ) );
        }
        //Create a list of aura elements
        mProcData.ElementListOnProcAura.set_size( mElementData.ElementActive.count(), 1, 0 );
        uint tVar = 0;
        for ( uint i = 0; i < mProcData.ElementListOnProcWithAura.length(); i++ )
        {
            if ( tDummyBitset.test( mProcData.ElementListOnProcWithAura( i ) ) == 1 )
            {
                // Save all elements, which are in the Aura
                mProcData.ElementListOnProcAura( tVar ) = mProcData.ElementListOnProcWithAura( i );
                tVar++;
            }
        }

        mProcData.ElementListOnProcAura.resize( tVar, 1 );
        // Calculate the Decomposition of the b-spline basis functions
        //Compute the basis functions of the first element on this processor
        Mat<uint> tBasis_of_element = mHMRElement.give_basis_of_element(
                mProcData.ElementListOnProc( 0 ),
                mMeshData.ModelDim,
                mMeshData.Polynomial,  // this is OK
                mMeshData.NumberOfElementsPerDirection );

        //Compute the position of the first basis function on this processor
        //
        Mat<uint> tBasis_position = mBasis.give_position_of_basis(
                tBasis_of_element( 0 ),
                mMeshData.ModelDim,
                mMeshData.Polynomial, // FIXME: should this be MaxPolynomial?
                mMeshData.NumberOfElementsPerDirection );
        mProcData.DecompBSpline.set_size( 6, 1, 0 );
        mProcData.DecompEdgeFace.set_size( 6, 1, 0 );
        //First processors at the left boundary (x-direction) takes the position directly, all other processors need to add the polynomial degree
        mProcData.DecompEdgeFace( 0 ) = tBasis_position( 0 );
        if( mProcData.ProcCoord( 0 ) == 0 )
        {
            mProcData.DecompBSpline( 0 ) = tBasis_position( 0 );
        }
        else
        {
            mProcData.DecompBSpline( 0 ) = tBasis_position( 0 ) + mMeshData.MaxPolynomial;
        }
        //First processors at the bottom boundary (y-direction) takes the position directly, all other processors need to add the polynomial degree
        mProcData.DecompEdgeFace( 2 ) = tBasis_position( 1 );
        if( mProcData.ProcCoord( 1 ) == 0 )
        {
            mProcData.DecompBSpline( 2 ) = tBasis_position( 1 );
        }
        else
        {
            mProcData.DecompBSpline( 2 ) = tBasis_position( 1 ) + mMeshData.MaxPolynomial;
        }
        if ( mMeshData.ModelDim == 3 )
        {
            //First processors at the back boundary (z-direction) takes the position directly, all other processors need to add the polynomial degree
            mProcData.DecompEdgeFace( 4 ) = tBasis_position( 2 );
            if( mProcData.ProcCoord( 2 ) == 0 )
            {
                mProcData.DecompBSpline( 4 ) = tBasis_position( 2 );
            }
            else
            {
                mProcData.DecompBSpline( 4 ) = tBasis_position( 2 ) + mMeshData.MaxPolynomial;
            }
        }

        //Compute the basis functions of the last element on this processor
        tBasis_of_element = mHMRElement.give_basis_of_element( mProcData.ElementListOnProc( mProcData.ElementListOnProc.length() - 1 ),
                mMeshData.ModelDim, mMeshData.Polynomial, mMeshData.NumberOfElementsPerDirection );
        //Compute the position of the last basis function on this processor
        tBasis_position = mBasis.give_position_of_basis( tBasis_of_element( tBasis_of_element.length() - 1 ), mMeshData.ModelDim, mMeshData.Polynomial, mMeshData.NumberOfElementsPerDirection );
        //First processors at the left boundary (x-direction) takes the position directly, all other processors need to add the polynomial degree
        (mProcData.DecompEdgeFace)( 1 ) = tBasis_position( 0 );
        if( mProcData.ProcCoord( 0 ) == 0 )
        {
            (mProcData.DecompBSpline)( 1 ) = tBasis_position( 0 );
        }
        else
        {
            (mProcData.DecompBSpline)( 1 ) = tBasis_position( 0 ) + mMeshData.MaxPolynomial;
        }
        //First processors at the bottom boundary (y-direction) takes the position directly, all other processors need to add the polynomial degree
        (mProcData.DecompEdgeFace)( 3 ) = tBasis_position( 1 );
        if( mProcData.ProcCoord( 1 ) == 0 )
        {
            (mProcData.DecompBSpline)( 3 ) = tBasis_position( 1 );
        }
        else
        {
            (mProcData.DecompBSpline)( 3 ) = tBasis_position( 1 ) + mMeshData.MaxPolynomial;
        }
        if ( mMeshData.ModelDim == 3 )
        {
            //First processors at the back boundary (z-direction) takes the position directly, all other processors need to add the polynomial degree
            (mProcData.DecompEdgeFace)( 5 ) = tBasis_position( 2 );
            if( mProcData.ProcCoord( 2 ) == 0 )
            {
                (mProcData.DecompBSpline)( 5 ) = tBasis_position( 2 );
            }
            else
            {
                (mProcData.DecompBSpline)( 5 ) = tBasis_position( 2 ) + mMeshData.MaxPolynomial;
            }
        }
        // Calculate the Decomposition of the lagrange basis functions

        //Compute the basis functions of the last element on this processor
        tBasis_of_element = mLagrangeElement.give_basis_of_element(
                mProcData.ElementListOnProc( mProcData.ElementListOnProc.length() - 1 ),
                mMeshData.ModelDim,
                mMeshData.Polynomial,
                mMeshData.NumberOfElementsPerDirection );

        //Compute the position of the last basis function on this processor
        tBasis_position = mLagrangeBasis.give_position_of_basis(
                tBasis_of_element( tBasis_of_element.length() - 1 ),
                mMeshData.ModelDim,
                mMeshData.Polynomial,
                mMeshData.NumberOfElementsPerDirection );

        if( mProcData.ProcCoord( 0 ) == 0 )
        {
            (mProcData.DecompLagrangeNode)( 1 ) = tBasis_position( 0 );
        }
        else
        {
            (mProcData.DecompLagrangeNode)( 1 ) = tBasis_position( 0 ) + 1;
        }
        if( mProcData.ProcCoord( 1 ) == 0 )
        {
            (mProcData.DecompLagrangeNode)( 3 ) = tBasis_position( 1 );
        }
        else
        {
            (mProcData.DecompLagrangeNode)( 3 ) = tBasis_position( 1 ) + 1;
        }
        if ( mMeshData.ModelDim == 3 )
        {
            if( mProcData.ProcCoord( 2 ) == 0 )
            {
                (mProcData.DecompLagrangeNode)( 5 ) = tBasis_position( 2 );
            }
            else
            {
                (mProcData.DecompLagrangeNode)( 5 ) = tBasis_position( 2 ) + 1;
            }
        }
    }
    mProcData.ElementListOnProcInit = mProcData.ElementListOnProc;
}

//--------------------------------------------------------------------------------
// FIXME remove hard coded MPI
void
Hierarchical_Mesh_Main::proc_coordinates_neighbors(
        Mat<uint> & aDimsCartProc)
{
    // Coordinates of the procs
    (mProcData.ProcCoord).set_size( 3, 1, 0 );
    // New communicator for cart and coordinates
    MPI_Comm new_comm;
    int dims[3];
    // Number of processors in x-direction
    dims[0] = aDimsCartProc( 0 );
    // Number of processors in y-direction
    dims[1] = aDimsCartProc( 1 );
    // Number of processors in z-direction
    dims[2] = aDimsCartProc( 2 );

    int periods[3];
    // No periodic boundary conditions are needed
    periods[0] = 0;
    periods[1] = 0;
    periods[2] = 0;
    // Ranks of the processors stay unaltered in the new communicator "new_comm"
    int reorder = 1;
    int my_rank;
    int coords[3];
    //Compute a cart of processors to define their position
    MPI_Cart_create( MPI_COMM_WORLD, mMeshData.ModelDim, dims, periods, reorder, &new_comm );
    //Define ranks
    MPI_Comm_rank( new_comm, &my_rank );
    //Define coordinates for the current processor
    MPI_Cart_coords( new_comm, my_rank, 3, coords );

    (mProcData.ProcCoord)( 0 ) = (uint)coords[0];
    (mProcData.ProcCoord)( 1 ) = (uint)coords[1];
    if ( mMeshData.ModelDim == 2 )
    {
        coords[2] = 0;
    }
    (mProcData.ProcCoord)( 2 ) = (uint)coords[2];

    // Calculate the Neighbors of "my_rank"
    //3D Representation
    // Layer on the back side of the rank element
    //|-----------------|
    //| N23 | N24 | N25 |
    //|-----------------|
    //| N21 | N5  | N22 |
    //|-----------------|
    //| N18 | N19 | N20 |
    //|-----------------|
    // Layer with rank
    //|------------------|
    //| N16 |  N2  | N17 |
    //|------------------|
    //| N3  | rank | N1  |
    //|------------------|
    //| N8  | N0   | N15 |
    //|------------------|
    // Layer in front of the rank element
    //|------------------|
    //| N12 | N13  | N14 |
    //|------------------|
    //| N7  | N4  | N11  |
    //|------------------|
    //| N6  | N9  | N10  |
    //|------------------|

    //in 2D
    // Layer with rank
    //|------------------|
    //| N16 |  N2  | N17 |
    //|------------------|
    //| N3  | rank | N1  |
    //|------------------|
    //| N8  | N0   | N15 |

    int rank;
    int error = 0;
    int coordsnew[3];
    // Neighbors of proc "my_rank"
    mProcData.ProcNeighbor.set_size(26,1,UINT_MAX);
    //------------------------------------------------------ Neighbor 0 --------------------
    // Neighbor at the bottom
    coordsnew[0] = coords[0];
    coordsnew[1] = coords[1] - 1;
    coordsnew[2] = coords[2];
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 0 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 1 --------------------
    // Neighbor at the right side
    coordsnew[0] = coords[0] + 1;
    coordsnew[1] = coords[1];
    coordsnew[2] = coords[2];
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 1 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 2 --------------------
    // Neighbor at the top side
    coordsnew[0] = coords[0];
    coordsnew[1] = coords[1] + 1;
    coordsnew[2] = coords[2];
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 2 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 3 --------------------
    // Neighbor at the left side
    coordsnew[0] = coords[0] - 1;
    coordsnew[1] = coords[1];
    coordsnew[2] = coords[2];
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 3 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 4 --------------------
    // Neighbor at the front side
    coordsnew[0] = coords[0];
    coordsnew[1] = coords[1];
    coordsnew[2] = coords[2] - 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 4 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 5 --------------------
    // Neighbor at the back side
    coordsnew[0] = coords[0];
    coordsnew[1] = coords[1];
    coordsnew[2] = coords[2] + 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 5 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 6 --------------------
    // Neighbor for the bottom left node
    coordsnew[0] = coords[0] - 1;
    coordsnew[1] = coords[1] - 1;
    coordsnew[2] = coords[2] - 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 6 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 7 --------------------
    // Neighbor at the left edge (front)
    coordsnew[0] = coords[0] - 1;
    coordsnew[1] = coords[1];
    coordsnew[2] = coords[2] - 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 7 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 8 --------------------
    // Neighbor at the left bottom edge
    coordsnew[0] = coords[0] - 1;
    coordsnew[1] = coords[1] - 1;
    coordsnew[2] = coords[2];
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 8 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 9 --------------------
    // Neighbor at the front bottom edge
    coordsnew[0] = coords[0];
    coordsnew[1] = coords[1] - 1;
    coordsnew[2] = coords[2] - 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor(9) = (uint)rank;
    //------------------------------------------------------ Neighbor 10 --------------------
    // Neighbor for the front bottom right node
    coordsnew[0] = coords[0] + 1;
    coordsnew[1] = coords[1] - 1;
    coordsnew[2] = coords[2] - 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 10 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 11 --------------------
    // Neighbor for the front top right edge
    coordsnew[0] = coords[0] + 1;
    coordsnew[1] = coords[1];
    coordsnew[2] = coords[2] - 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 11 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 12 --------------------
    // Neighbor for the front top left node
    coordsnew[0] = coords[0] - 1;
    coordsnew[1] = coords[1] + 1;
    coordsnew[2] = coords[2] - 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 12 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 13 --------------------
    // Neighbor for the front top  edge
    coordsnew[0] = coords[0];
    coordsnew[1] = coords[1] + 1;
    coordsnew[2] = coords[2] - 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 13 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 14 --------------------
    // Neighbor for the front top right node
    coordsnew[0] = coords[0] + 1;
    coordsnew[1] = coords[1] + 1;
    coordsnew[2] = coords[2] - 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 14 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 15 --------------------
    // Neighbor bottom right node
    coordsnew[0] = coords[0] + 1;
    coordsnew[1] = coords[1] - 1;
    coordsnew[2] = coords[2];
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 15 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 16 --------------------
    // Neighbor top left node
    coordsnew[0] = coords[0] - 1;
    coordsnew[1] = coords[1] + 1;
    coordsnew[2] = coords[2];
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 16 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 17 --------------------
    // Neighbor top right node
    coordsnew[0] = coords[0] + 1;
    coordsnew[1] = coords[1] + 1;
    coordsnew[2] = coords[2];
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 17 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 18 --------------------
    // Neighbor in the back side bottom left node
    coordsnew[0] = coords[0] - 1;
    coordsnew[1] = coords[1] - 1;
    coordsnew[2] = coords[2] + 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 18 ) = (uint)rank;

    //------------------------------------------------------ Neighbor 19 --------------------
    // Neighbor in the back side bottom edge
    coordsnew[0] = coords[0];
    coordsnew[1] = coords[1] - 1;
    coordsnew[2] = coords[2] + 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 19 ) = (uint)rank;

    //------------------------------------------------------ Neighbor 20 --------------------
    // Neighbor in the back side bottom right node
    coordsnew[0] = coords[0] + 1;
    coordsnew[1] = coords[1] - 1;
    coordsnew[2] = coords[2] + 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 20 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 21 --------------------
    // Neighbor in the back side left edge
    coordsnew[0] = coords[0] - 1;
    coordsnew[1] = coords[1];
    coordsnew[2] = coords[2] + 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if (error == 0 )
        mProcData.ProcNeighbor( 21 ) = (uint)rank;

    //------------------------------------------------------ Neighbor 22 --------------------
    // Neighbor in the back side right edge
    coordsnew[0] = coords[0] + 1;
    coordsnew[1] = coords[1];
    coordsnew[2] = coords[2] + 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 22 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 23 --------------------
    // Neighbor in the back side top left node
    coordsnew[0] = coords[0] - 1;
    coordsnew[1] = coords[1] + 1;
    coordsnew[2] = coords[2] + 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 23 ) = (uint)rank;
    //------------------------------------------------------ Neighbor 24 --------------------
    coordsnew[0] = coords[0];
    coordsnew[1] = coords[1] + 1;
    coordsnew[2] = coords[2] + 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 24 ) = (uint)rank;

    //------------------------------------------------------ Neighbor 25 --------------------
    coordsnew[0] = coords[0] + 1;
    coordsnew[1] = coords[1] + 1;
    coordsnew[2] = coords[2] + 1;
    error = 1;
    if ( coordsnew[0] >= 0 && coordsnew[1] >= 0 && coordsnew[2] >= 0
            && coordsnew[0] < dims[0] && coordsnew[1] < dims[1] && coordsnew[2] < dims[2] )
        error = MPI_Cart_rank( new_comm, coordsnew, &rank );
    if ( error == 0 )
        mProcData.ProcNeighbor( 25 ) = (uint)rank;

}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::L2_projection(
        const map<uint, real> & aNodalField,
        const real              aRhsScaling,
        const real              aRhsOffset)
{
    // tUsingSolvers = 0 for old solvers, tUsingSolvers = 1 for new solvers
    uint tUsingSolvers = 0;

    switch(tUsingSolvers)
    {
    case (0):
    {
        uint tLevel;
        Mat<uint> tBasis;
        Mat<real> tOldSol;
        uint tBasisInElement = pow( mMeshData.PolynomialDesign + 1, mMeshData.ModelDim );
        Mat<real> tNewSol( tBasisInElement, 1 );
        Mat<real> ElementRHS;
        Mat<real> ElementRHSDummy;
        Mat<real> ElementMass;
        Mat<uint> tBasisOfElement;
        uint tVar = 0;
        Mat<uint> tListOfChildren;
        Cell<Mat<real>> tTMatrixOfChildren;
        Mat<uint> tIdField;
        Mat<real> tTMatrix;

        // Initialize global right hand side vector
        Mat<real> GLBRHS( mBasisData.DesignBSplineActive.count(), 1, 0.0);

        // Initialize vector-based storage of sparse mass matrix
        moris::Mat< moris::uint > RowInd( (mFieldData.SparseSize), 1, 0 );
        moris::Mat< moris::uint > ColInd( (mFieldData.SparseSize), 1, 0 );
        moris::Mat< moris::real > Values( (mFieldData.SparseSize), 1, 0 );

        // Start timer
        moris::tic tBuildSystemTiming;

        // Build linear system
        for ( uint i = 0; i < (mElementData.ElementListActiveDesign).length(); i++)
        {
            //Check if I need to be coarsened or refined
            mTMatrix.give_Tmatrix_and_IdField_Design_L2Projection_Coarsening(
                    mElementData.ElementListActiveDesign( i ),
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.LevelLastStep,
                    mMeshData.LevelDesignLastStep,
                    mMeshData.NumberOfElementsPerDirection,
                    mElementData.ElementActiveDesignLastStep,
                    tTMatrixOfChildren,
                    tListOfChildren,
                    mElementData.TMatrixParentChildRelationDesign );

            //If i have children, from where i get the solution, i have to use the coarsening operator
            if ( tListOfChildren.length() >= pow( 2, mMeshData.ModelDim ) )
            {
                ElementRHS.set_size( tBasisInElement, 1, 0 );
                for ( uint j = 0; j < tListOfChildren.length(); j++ )
                {
                    tBasisOfElement = mHMRElement.give_basis_of_element( tListOfChildren( j ), mMeshData.ModelDim, mMeshData.PolynomialDesign, mMeshData.NumberOfElementsPerDirection );
                    tTMatrix = tTMatrixOfChildren(j);
                    tOldSol.set_size(tBasisOfElement.length(),1,0);
                    for ( uint k = 0; k < tBasisOfElement.length(); k++)
                    {
                        tOldSol(k) = aRhsScaling*aNodalField.find(tBasisOfElement(k))+aRhsOffset;
                    }
                    tNewSol         = tOldSol;
                    tLevel          = mBaseElement.give_element_level(tListOfChildren(j),mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
                    ElementMass     = MassMatrix_for_L2_projection(mMeshData.ModelDim,mMeshData.PolynomialDesign);
                    ElementRHSDummy = 1/pow(pow(2,tLevel),mMeshData.ModelDim) * tTMatrix * ElementMass * tNewSol;
                    ElementRHS      = ElementRHS + ElementRHSDummy;
                }
            }
            else
            {
                //If I get the solution from the same element or a parent element, I have to use the refinement operator with the TMatrix
                if ( mSettings.TruncatedBsplines == false)
                {
                    mTMatrix.give_Tmatrix_and_IdField(
                            mElementData.ElementListActiveDesign(i),
                            mMeshData.ModelDim,
                            mMeshData.PolynomialDesign,
                            mMeshData.NumberOfElementsPerDirection,
                            mBasisData.DesignBSplineActiveLastStep,
                            tTMatrix,
                            tIdField,
                            mElementData.TMatrixParentChildRelationDesign);
                }
                else
                {
                    mTMatrix.give_Truncated_Tmatrix_and_IdField(
                            mElementData.ElementListActiveDesign(i),
                            mMeshData.ModelDim,
                            mMeshData.PolynomialDesign,
                            mMeshData.NumberOfElementsPerDirection,
                            mBasisData.DesignBSplineActiveLastStep,
                            tTMatrix,
                            tIdField,
                            mElementData.TMatrixParentChildRelationDesign);
                }
                tBasis = mHMRElement.give_basis_of_element(
                        mElementData.ElementListActiveDesign(i),
                        mMeshData.ModelDim,
                        mMeshData.PolynomialDesign,
                        mMeshData.NumberOfElementsPerDirection);

                real tWeight;
                for ( uint j = 0; j < tTMatrix.n_cols(); j++)
                {
                    tNewSol(j) = 0.0;
                    tWeight = 0.0;
                    for ( uint k = 0; k < tTMatrix.n_rows(); k++)
                    {
                        tWeight    += tTMatrix(k,j);
                        tNewSol(j) += tTMatrix(k,j) * (aRhsScaling*aNodalField.find(tIdField(k))+aRhsOffset);
                    }
                    if ( mSettings.PerformNormalization )
                    {
                        tNewSol(j) /= tWeight;
                    }
                }
                tLevel = mBaseElement.give_element_level(
                        mElementData.ElementListActiveDesign(i),
                        mMeshData.ModelDim,
                        mMeshData.NumberOfElementsPerDirection);

                ElementRHS = RHS_for_L2_projection(
                        mMeshData.ModelDim,
                        mMeshData.PolynomialDesign,
                        tNewSol);

                ElementRHS = 1/pow(pow(2,tLevel),mMeshData.ModelDim) * ElementRHS;
            }

            //Create T-Matrix from current element
            if ( mSettings.TruncatedBsplines == false)
            {
                mTMatrix.give_Tmatrix_and_IdField(
                        mElementData.ElementListActiveDesign(i),
                        mMeshData.ModelDim,
                        mMeshData.PolynomialDesign,
                        mMeshData.NumberOfElementsPerDirection,
                        mBasisData.DesignBSplineActive,
                        tTMatrix,
                        tIdField,
                        mElementData.TMatrixParentChildRelationDesign);
            }
            else
            {
                mTMatrix.give_Truncated_Tmatrix_and_IdField(
                        mElementData.ElementListActiveDesign(i),
                        mMeshData.ModelDim,
                        mMeshData.PolynomialDesign,
                        mMeshData.NumberOfElementsPerDirection,
                        mBasisData.DesignBSplineActive,
                        tTMatrix,
                        tIdField,
                        mElementData.TMatrixParentChildRelationDesign);
            }

            //Use T-Matrix to transform to current mesh and use ID-Field for assembling
            ElementRHS = tTMatrix * ElementRHS;
            ElementMass = MassMatrix_for_L2_projection(
                    mMeshData.ModelDim,
                    mMeshData.PolynomialDesign);

            tLevel = mBaseElement.give_element_level(
                    mElementData.ElementListActiveDesign(i),
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);

            ElementMass = 1/pow(pow(2,tLevel),mMeshData.ModelDim) * tTMatrix * ElementMass * trans(tTMatrix);
            for ( uint j = 0; j< tIdField.length(); j++)
            {
                GLBRHS(mBasisData.DesignBSplineListMap.find(tIdField(j))) += ElementRHS(j);
                for ( uint k = 0; k < tIdField.length(); k++)
                {
                    RowInd(tVar) = mBasisData.DesignBSplineListMap.find(tIdField(j));
                    ColInd(tVar) = mBasisData.DesignBSplineListMap.find(tIdField(k));
                    Values(tVar) = ElementMass(j,k);
                    tVar++;
                }
            }
        }

        // Trim vectors storing mass matrix
        RowInd.resize(tVar,1);
        ColInd.resize(tVar,1);
        Values.resize(tVar,1);

        // Build sparse matrix
        Sp_Mat<real> GLBMass(RowInd,ColInd,Values,GLBRHS.length(),GLBRHS.length()); // Generate global matrix from vectors

        // Write time for building linear system
        real tElapsedTimeBuildSystem = tBuildSystemTiming.toc<moris::chronos::seconds>().wall;
        std::fprintf(stdout,"Time for building projection system : %f [sec]\n",tElapsedTimeBuildSystem);

        std::cout << "Size of matrix to solve " << GLBRHS.length() << std::endl;

        // Start timer for solving linear system
        moris::tic tSolveSystemTiming;

        // call solver and write values into BSplineCoeff
        mFieldData.BSplineCoeff = solve(GLBMass, GLBRHS , "bicgstab");

        //mFieldData.BSplineCoeff.print("mFieldData.BSplineCoeff");

        // Write time for solving linear system
        real tElapsedTimeSolveSystem = tSolveSystemTiming.toc<moris::chronos::seconds>().wall;
        std::fprintf(stdout,"Time for solving projection system : %f [sec]\n",tElapsedTimeSolveSystem);
        break;
    }
    case(1):
    {
        // create solver input object
        Hierarchical_Mesh_Solver_Input*  tSolverInput;
        tSolverInput = new Hierarchical_Mesh_Solver_Input(  aRhsScaling,
                                                            aRhsOffset,
                                                            this );
        // Start timer for solving linear system
        moris::tic tSolveSystemTiming;

        // create solver factory
        Solver_Factory  tSolFactory;

        // create solver object
        std::shared_ptr<Linear_Solver> tLin = tSolFactory.create_solver( tSolverInput );

        //tLin->set_param("max_its")   = 200;
        //tLin->set_param("solver_type")   = AZ_cg;

        // call solve
        tLin->solve_linear_system();

        // Write time for solving linear system
        real tElapsedTimeSolveSystem = tSolveSystemTiming.toc<moris::chronos::seconds>().wall;
        std::fprintf(stdout,"Time for building and solving projection system : %f [sec]\n",tElapsedTimeSolveSystem);

        // Get solution vector
        mFieldData.BSplineCoeff.set_size ( mBasisData.DesignBSplineActive.count(), 1 );
        tLin->get_solution(mFieldData.BSplineCoeff);

        delete tSolverInput;

        break;
    }
    default:
        MORIS_ASSERT( false, "No L2_projection solver specified" );
        break;
    }
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::RHS_for_L2_projection(
        const uint      & aDim,
        const uint      & aPolynomial,
        const Mat<real> & aNodalField,
        Mat< real >     & aElementRHS)
{
    aElementRHS.set_size( pow( aPolynomial + 1, aDim), 1 );

    if ( aDim == 2)
    {
        if ( aPolynomial == 1)
        {
            aElementRHS( 0 ) = 1.0/36.0*(4.0*aNodalField( 0 ) + 2.0*aNodalField( 1 ) + 2.0*aNodalField( 2 ) + 1.0*aNodalField( 3 ));
            aElementRHS( 1 ) = 1.0/36.0*(2.0*aNodalField( 0 ) + 4.0*aNodalField( 1 ) + 1.0*aNodalField( 2 ) + 2.0*aNodalField( 3 ));
            aElementRHS( 2 ) = 1.0/36.0*(2.0*aNodalField( 0 ) + 1.0*aNodalField( 1 ) + 4.0*aNodalField( 2 ) + 2.0*aNodalField( 3 ));
            aElementRHS( 3 ) = 1.0/36.0*(1.0*aNodalField( 0 ) + 2.0*aNodalField( 1 ) + 2.0*aNodalField( 2 ) + 4.0*aNodalField( 3 ));
        }
        else if ( aPolynomial == 2)
        {
            aElementRHS( 0 ) = 0.002500000000000*aNodalField( 0 ) + 0.005416666666667*aNodalField( 1 ) + 0.000416666666667*aNodalField( 2 ) + 0.005416666666667*aNodalField( 3 ) + 0.011736111111111*aNodalField( 4 ) + 0.000902777777778*aNodalField( 5 ) + 0.000416666666667*aNodalField( 6 ) + 0.000902777777778*aNodalField( 7 ) + 0.000069444444444*aNodalField( 8 );
            aElementRHS( 1 ) = 0.005416666666667*aNodalField( 0 ) + 0.022500000000000*aNodalField( 1 ) + 0.005416666666667*aNodalField( 2 ) + 0.011736111111111*aNodalField( 3 ) + 0.048750000000000*aNodalField( 4 ) + 0.011736111111111*aNodalField( 5 ) + 0.000902777777778*aNodalField( 6 ) + 0.003750000000000*aNodalField( 7 ) + 0.000902777777778*aNodalField( 8 );
            aElementRHS( 2 ) = 0.000416666666667*aNodalField( 0 ) + 0.005416666666667*aNodalField( 1 ) + 0.002500000000000*aNodalField( 2 ) + 0.000902777777778*aNodalField( 3 ) + 0.011736111111111*aNodalField( 4 ) + 0.005416666666667*aNodalField( 5 ) + 0.000069444444444*aNodalField( 6 ) + 0.000902777777778*aNodalField( 7 ) + 0.000416666666667*aNodalField( 8 );
            aElementRHS( 3 ) = 0.005416666666667*aNodalField( 0 ) + 0.011736111111111*aNodalField( 1 ) + 0.000902777777778*aNodalField( 2 ) + 0.022500000000000*aNodalField( 3 ) + 0.048750000000000*aNodalField( 4 ) + 0.003750000000000*aNodalField( 5 ) + 0.005416666666667*aNodalField( 6 ) + 0.011736111111111*aNodalField( 7 ) + 0.000902777777778*aNodalField( 8 );
            aElementRHS( 4 ) = 0.011736111111111*aNodalField( 0 ) + 0.048750000000000*aNodalField( 1 ) + 0.011736111111111*aNodalField( 2 ) + 0.048750000000000*aNodalField( 3 ) + 0.202500000000000*aNodalField( 4 ) + 0.048750000000000*aNodalField( 5 ) + 0.011736111111111*aNodalField( 6 ) + 0.048750000000000*aNodalField( 7 ) + 0.011736111111111*aNodalField( 8 );
            aElementRHS( 5 ) = 0.000902777777778*aNodalField( 0 ) + 0.011736111111111*aNodalField( 1 ) + 0.005416666666667*aNodalField( 2 ) + 0.003750000000000*aNodalField( 3 ) + 0.048750000000000*aNodalField( 4 ) + 0.022500000000000*aNodalField( 5 ) + 0.000902777777778*aNodalField( 6 ) + 0.011736111111111*aNodalField( 7 ) + 0.005416666666667*aNodalField( 8 );
            aElementRHS( 6 ) = 0.000416666666667*aNodalField( 0 ) + 0.000902777777778*aNodalField( 1 ) + 0.000069444444444*aNodalField( 2 ) + 0.005416666666667*aNodalField( 3 ) + 0.011736111111111*aNodalField( 4 ) + 0.000902777777778*aNodalField( 5 ) + 0.002500000000000*aNodalField( 6 ) + 0.005416666666667*aNodalField( 7 ) + 0.000416666666667*aNodalField( 8 );
            aElementRHS( 7 ) = 0.000902777777778*aNodalField( 0 ) + 0.003750000000000*aNodalField( 1 ) + 0.000902777777778*aNodalField( 2 ) + 0.011736111111111*aNodalField( 3 ) + 0.048750000000000*aNodalField( 4 ) + 0.011736111111111*aNodalField( 5 ) + 0.005416666666667*aNodalField( 6 ) + 0.022500000000000*aNodalField( 7 ) + 0.005416666666667*aNodalField( 8 );
            aElementRHS( 8 ) = 0.000069444444444*aNodalField( 0 ) + 0.000902777777778*aNodalField( 1 ) + 0.000416666666667*aNodalField( 2 ) + 0.000902777777778*aNodalField( 3 ) + 0.011736111111111*aNodalField( 4 ) + 0.005416666666667*aNodalField( 5 ) + 0.000416666666667*aNodalField( 6 ) + 0.005416666666667*aNodalField( 7 ) + 0.002500000000000*aNodalField( 8 );
        }
    }
    else if ( aDim == 3)
    {
        if ( aPolynomial == 1)
        {
            aElementRHS( 0 ) = 1.0/216.0*(8.0*aNodalField( 0 ) + 4.0*aNodalField( 1 ) + 4.0*aNodalField( 2 ) + 2.0*aNodalField( 3 ) + 4.0*aNodalField( 4 ) + 2.0*aNodalField( 5 ) + 2.0*aNodalField( 6 ) + 1.0*aNodalField( 7 ));
            aElementRHS( 1 ) = 1.0/216.0*(4.0*aNodalField( 0 ) + 8.0*aNodalField( 1 ) + 2.0*aNodalField( 2 ) + 4.0*aNodalField( 3 ) + 2.0*aNodalField( 4 ) + 4.0*aNodalField( 5 ) + 1.0*aNodalField( 6 ) + 2.0*aNodalField( 7 ));
            aElementRHS( 2 ) = 1.0/216.0*(4.0*aNodalField( 0 ) + 2.0*aNodalField( 1 ) + 8.0*aNodalField( 2 ) + 4.0*aNodalField( 3 ) + 2.0*aNodalField( 4 ) + 1.0*aNodalField( 5 ) + 4.0*aNodalField( 6 ) + 2.0*aNodalField( 7 ));
            aElementRHS( 3 ) = 1.0/216.0*(2.0*aNodalField( 0 ) + 4.0*aNodalField( 1 ) + 4.0*aNodalField( 2 ) + 8.0*aNodalField( 3 ) + 1.0*aNodalField( 4 ) + 2.0*aNodalField( 5 ) + 2.0*aNodalField( 6 ) + 4.0*aNodalField( 7 ));
            aElementRHS( 4 ) = 1.0/216.0*(4.0*aNodalField( 0 ) + 2.0*aNodalField( 1 ) + 2.0*aNodalField( 2 ) + 1.0*aNodalField( 3 ) + 8.0*aNodalField( 4 ) + 4.0*aNodalField( 5 ) + 4.0*aNodalField( 6 ) + 2.0*aNodalField( 7 ));
            aElementRHS( 5 ) = 1.0/216.0*(2.0*aNodalField( 0 ) + 4.0*aNodalField( 1 ) + 1.0*aNodalField( 2 ) + 2.0*aNodalField( 3 ) + 4.0*aNodalField( 4 ) + 8.0*aNodalField( 5 ) + 2.0*aNodalField( 6 ) + 4.0*aNodalField( 7 ));
            aElementRHS( 6 ) = 1.0/216.0*(2.0*aNodalField( 0 ) + 1.0*aNodalField( 1 ) + 4.0*aNodalField( 2 ) + 2.0*aNodalField( 3 ) + 4.0*aNodalField( 4 ) + 2.0*aNodalField( 5 ) + 8.0*aNodalField( 6 ) + 4.0*aNodalField( 7 ));
            aElementRHS( 7 ) = 1.0/216.0*(1.0*aNodalField( 0 ) + 2.0*aNodalField( 1 ) + 2.0*aNodalField( 2 ) + 4.0*aNodalField( 3 ) + 2.0*aNodalField( 4 ) + 4.0*aNodalField( 5 ) + 4.0*aNodalField( 6 ) + 8.0*aNodalField( 7 ));
        }
        else if ( aPolynomial == 2)
        {
            aElementRHS(  0 ) = 1.250000e-04*aNodalField( 0 ) + 2.708333e-04*aNodalField( 1 ) + 2.083333e-05*aNodalField( 2 ) + 2.708333e-04*aNodalField( 3 ) + 5.868056e-04*aNodalField( 4 ) + 4.513889e-05*aNodalField( 5 ) + 2.083333e-05*aNodalField( 6 ) + 4.513889e-05*aNodalField( 7 ) + 3.472222e-06*aNodalField( 8 ) + 2.708333e-04*aNodalField(9) + 5.868056e-04*aNodalField( 10 ) + 4.513889e-05*aNodalField( 11 ) + 5.868056e-04*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 9.780093e-05*aNodalField( 14 ) + 4.513889e-05*aNodalField( 15 ) + 9.780093e-05*aNodalField( 16 ) + 7.523148e-06*aNodalField( 17 ) + 2.083333e-05*aNodalField( 18 ) + 4.513889e-05*aNodalField( 19 ) + 3.472222e-06*aNodalField( 20 ) + 4.513889e-05*aNodalField( 21 ) + 9.780093e-05*aNodalField( 22 ) + 7.523148e-06*aNodalField( 23 ) + 3.472222e-06*aNodalField( 24 ) + 7.523148e-06*aNodalField( 25 ) + 5.787037e-07*aNodalField( 26 );
            aElementRHS(  1 ) = 2.708333e-04*aNodalField( 0 ) + 1.125000e-03*aNodalField( 1 ) + 2.708333e-04*aNodalField( 2 ) + 5.868056e-04*aNodalField( 3 ) + 2.437500e-03*aNodalField( 4 ) + 5.868056e-04*aNodalField( 5 ) + 4.513889e-05*aNodalField( 6 ) + 1.875000e-04*aNodalField( 7 ) + 4.513889e-05*aNodalField( 8 ) + 5.868056e-04*aNodalField(9) + 2.437500e-03*aNodalField( 10 ) + 5.868056e-04*aNodalField( 11 ) + 1.271412e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 1.271412e-03*aNodalField( 14 ) + 9.780093e-05*aNodalField( 15 ) + 4.062500e-04*aNodalField( 16 ) + 9.780093e-05*aNodalField( 17 ) + 4.513889e-05*aNodalField( 18 ) + 1.875000e-04*aNodalField( 19 ) + 4.513889e-05*aNodalField( 20 ) + 9.780093e-05*aNodalField( 21 ) + 4.062500e-04*aNodalField( 22 ) + 9.780093e-05*aNodalField( 23 ) + 7.523148e-06*aNodalField( 24 ) + 3.125000e-05*aNodalField( 25 ) + 7.523148e-06*aNodalField( 26 );
            aElementRHS(  2 ) = 2.083333e-05*aNodalField( 0 ) + 2.708333e-04*aNodalField( 1 ) + 1.250000e-04*aNodalField( 2 ) + 4.513889e-05*aNodalField( 3 ) + 5.868056e-04*aNodalField( 4 ) + 2.708333e-04*aNodalField( 5 ) + 3.472222e-06*aNodalField( 6 ) + 4.513889e-05*aNodalField( 7 ) + 2.083333e-05*aNodalField( 8 ) + 4.513889e-05*aNodalField(9) + 5.868056e-04*aNodalField( 10 ) + 2.708333e-04*aNodalField( 11 ) + 9.780093e-05*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 5.868056e-04*aNodalField( 14 ) + 7.523148e-06*aNodalField( 15 ) + 9.780093e-05*aNodalField( 16 ) + 4.513889e-05*aNodalField( 17 ) + 3.472222e-06*aNodalField( 18 ) + 4.513889e-05*aNodalField( 19 ) + 2.083333e-05*aNodalField( 20 ) + 7.523148e-06*aNodalField( 21 ) + 9.780093e-05*aNodalField( 22 ) + 4.513889e-05*aNodalField( 23 ) + 5.787037e-07*aNodalField( 24 ) + 7.523148e-06*aNodalField( 25 ) + 3.472222e-06*aNodalField( 26 );
            aElementRHS(  3 ) = 2.708333e-04*aNodalField( 0 ) + 5.868056e-04*aNodalField( 1 ) + 4.513889e-05*aNodalField( 2 ) + 1.125000e-03*aNodalField( 3 ) + 2.437500e-03*aNodalField( 4 ) + 1.875000e-04*aNodalField( 5 ) + 2.708333e-04*aNodalField( 6 ) + 5.868056e-04*aNodalField( 7 ) + 4.513889e-05*aNodalField( 8 ) + 5.868056e-04*aNodalField(9) + 1.271412e-03*aNodalField( 10 ) + 9.780093e-05*aNodalField( 11 ) + 2.437500e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 4.062500e-04*aNodalField( 14 ) + 5.868056e-04*aNodalField( 15 ) + 1.271412e-03*aNodalField( 16 ) + 9.780093e-05*aNodalField( 17 ) + 4.513889e-05*aNodalField( 18 ) + 9.780093e-05*aNodalField( 19 ) + 7.523148e-06*aNodalField( 20 ) + 1.875000e-04*aNodalField( 21 ) + 4.062500e-04*aNodalField( 22 ) + 3.125000e-05*aNodalField( 23 ) + 4.513889e-05*aNodalField( 24 ) + 9.780093e-05*aNodalField( 25 ) + 7.523148e-06*aNodalField( 26 );
            aElementRHS(  4 ) = 5.868056e-04*aNodalField( 0 ) + 2.437500e-03*aNodalField( 1 ) + 5.868056e-04*aNodalField( 2 ) + 2.437500e-03*aNodalField( 3 ) + 1.012500e-02*aNodalField( 4 ) + 2.437500e-03*aNodalField( 5 ) + 5.868056e-04*aNodalField( 6 ) + 2.437500e-03*aNodalField( 7 ) + 5.868056e-04*aNodalField( 8 ) + 1.271412e-03*aNodalField(9) + 5.281250e-03*aNodalField( 10 ) + 1.271412e-03*aNodalField( 11 ) + 5.281250e-03*aNodalField( 12 ) + 2.193750e-02*aNodalField( 13 ) + 5.281250e-03*aNodalField( 14 ) + 1.271412e-03*aNodalField( 15 ) + 5.281250e-03*aNodalField( 16 ) + 1.271412e-03*aNodalField( 17 ) + 9.780093e-05*aNodalField( 18 ) + 4.062500e-04*aNodalField( 19 ) + 9.780093e-05*aNodalField( 20 ) + 4.062500e-04*aNodalField( 21 ) + 1.687500e-03*aNodalField( 22 ) + 4.062500e-04*aNodalField( 23 ) + 9.780093e-05*aNodalField( 24 ) + 4.062500e-04*aNodalField( 25 ) + 9.780093e-05*aNodalField( 26 );
            aElementRHS(  5 ) = 4.513889e-05*aNodalField( 0 ) + 5.868056e-04*aNodalField( 1 ) + 2.708333e-04*aNodalField( 2 ) + 1.875000e-04*aNodalField( 3 ) + 2.437500e-03*aNodalField( 4 ) + 1.125000e-03*aNodalField( 5 ) + 4.513889e-05*aNodalField( 6 ) + 5.868056e-04*aNodalField( 7 ) + 2.708333e-04*aNodalField( 8 ) + 9.780093e-05*aNodalField(9) + 1.271412e-03*aNodalField( 10 ) + 5.868056e-04*aNodalField( 11 ) + 4.062500e-04*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 2.437500e-03*aNodalField( 14 ) + 9.780093e-05*aNodalField( 15 ) + 1.271412e-03*aNodalField( 16 ) + 5.868056e-04*aNodalField( 17 ) + 7.523148e-06*aNodalField( 18 ) + 9.780093e-05*aNodalField( 19 ) + 4.513889e-05*aNodalField( 20 ) + 3.125000e-05*aNodalField( 21 ) + 4.062500e-04*aNodalField( 22 ) + 1.875000e-04*aNodalField( 23 ) + 7.523148e-06*aNodalField( 24 ) + 9.780093e-05*aNodalField( 25 ) + 4.513889e-05*aNodalField( 26 );
            aElementRHS(  6 ) = 2.083333e-05*aNodalField( 0 ) + 4.513889e-05*aNodalField( 1 ) + 3.472222e-06*aNodalField( 2 ) + 2.708333e-04*aNodalField( 3 ) + 5.868056e-04*aNodalField( 4 ) + 4.513889e-05*aNodalField( 5 ) + 1.250000e-04*aNodalField( 6 ) + 2.708333e-04*aNodalField( 7 ) + 2.083333e-05*aNodalField( 8 ) + 4.513889e-05*aNodalField(9) + 9.780093e-05*aNodalField( 10 ) + 7.523148e-06*aNodalField( 11 ) + 5.868056e-04*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 9.780093e-05*aNodalField( 14 ) + 2.708333e-04*aNodalField( 15 ) + 5.868056e-04*aNodalField( 16 ) + 4.513889e-05*aNodalField( 17 ) + 3.472222e-06*aNodalField( 18 ) + 7.523148e-06*aNodalField( 19 ) + 5.787037e-07*aNodalField( 20 ) + 4.513889e-05*aNodalField( 21 ) + 9.780093e-05*aNodalField( 22 ) + 7.523148e-06*aNodalField( 23 ) + 2.083333e-05*aNodalField( 24 ) + 4.513889e-05*aNodalField( 25 ) + 3.472222e-06*aNodalField( 26 );
            aElementRHS(  7 ) = 4.513889e-05*aNodalField( 0 ) + 1.875000e-04*aNodalField( 1 ) + 4.513889e-05*aNodalField( 2 ) + 5.868056e-04*aNodalField( 3 ) + 2.437500e-03*aNodalField( 4 ) + 5.868056e-04*aNodalField( 5 ) + 2.708333e-04*aNodalField( 6 ) + 1.125000e-03*aNodalField( 7 ) + 2.708333e-04*aNodalField( 8 ) + 9.780093e-05*aNodalField(9) + 4.062500e-04*aNodalField( 10 ) + 9.780093e-05*aNodalField( 11 ) + 1.271412e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 1.271412e-03*aNodalField( 14 ) + 5.868056e-04*aNodalField( 15 ) + 2.437500e-03*aNodalField( 16 ) + 5.868056e-04*aNodalField( 17 ) + 7.523148e-06*aNodalField( 18 ) + 3.125000e-05*aNodalField( 19 ) + 7.523148e-06*aNodalField( 20 ) + 9.780093e-05*aNodalField( 21 ) + 4.062500e-04*aNodalField( 22 ) + 9.780093e-05*aNodalField( 23 ) + 4.513889e-05*aNodalField( 24 ) + 1.875000e-04*aNodalField( 25 ) + 4.513889e-05*aNodalField( 26 );
            aElementRHS(  8 ) = 3.472222e-06*aNodalField( 0 ) + 4.513889e-05*aNodalField( 1 ) + 2.083333e-05*aNodalField( 2 ) + 4.513889e-05*aNodalField( 3 ) + 5.868056e-04*aNodalField( 4 ) + 2.708333e-04*aNodalField( 5 ) + 2.083333e-05*aNodalField( 6 ) + 2.708333e-04*aNodalField( 7 ) + 1.250000e-04*aNodalField( 8 ) + 7.523148e-06*aNodalField(9) + 9.780093e-05*aNodalField( 10 ) + 4.513889e-05*aNodalField( 11 ) + 9.780093e-05*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 5.868056e-04*aNodalField( 14 ) + 4.513889e-05*aNodalField( 15 ) + 5.868056e-04*aNodalField( 16 ) + 2.708333e-04*aNodalField( 17 ) + 5.787037e-07*aNodalField( 18 ) + 7.523148e-06*aNodalField( 19 ) + 3.472222e-06*aNodalField( 20 ) + 7.523148e-06*aNodalField( 21 ) + 9.780093e-05*aNodalField( 22 ) + 4.513889e-05*aNodalField( 23 ) + 3.472222e-06*aNodalField( 24 ) + 4.513889e-05*aNodalField( 25 ) + 2.083333e-05*aNodalField( 26 );
            aElementRHS(  9 ) = 2.708333e-04*aNodalField( 0 ) + 5.868056e-04*aNodalField( 1 ) + 4.513889e-05*aNodalField( 2 ) + 5.868056e-04*aNodalField( 3 ) + 1.271412e-03*aNodalField( 4 ) + 9.780093e-05*aNodalField( 5 ) + 4.513889e-05*aNodalField( 6 ) + 9.780093e-05*aNodalField( 7 ) + 7.523148e-06*aNodalField( 8 ) + 1.125000e-03*aNodalField(9) + 2.437500e-03*aNodalField( 10 ) + 1.875000e-04*aNodalField( 11 ) + 2.437500e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 4.062500e-04*aNodalField( 14 ) + 1.875000e-04*aNodalField( 15 ) + 4.062500e-04*aNodalField( 16 ) + 3.125000e-05*aNodalField( 17 ) + 2.708333e-04*aNodalField( 18 ) + 5.868056e-04*aNodalField( 19 ) + 4.513889e-05*aNodalField( 20 ) + 5.868056e-04*aNodalField( 21 ) + 1.271412e-03*aNodalField( 22 ) + 9.780093e-05*aNodalField( 23 ) + 4.513889e-05*aNodalField( 24 ) + 9.780093e-05*aNodalField( 25 ) + 7.523148e-06*aNodalField( 26 );
            aElementRHS( 10 ) = 5.868056e-04*aNodalField( 0 ) + 2.437500e-03*aNodalField( 1 ) + 5.868056e-04*aNodalField( 2 ) + 1.271412e-03*aNodalField( 3 ) + 5.281250e-03*aNodalField( 4 ) + 1.271412e-03*aNodalField( 5 ) + 9.780093e-05*aNodalField( 6 ) + 4.062500e-04*aNodalField( 7 ) + 9.780093e-05*aNodalField( 8 ) + 2.437500e-03*aNodalField(9) + 1.012500e-02*aNodalField( 10 ) + 2.437500e-03*aNodalField( 11 ) + 5.281250e-03*aNodalField( 12 ) + 2.193750e-02*aNodalField( 13 ) + 5.281250e-03*aNodalField( 14 ) + 4.062500e-04*aNodalField( 15 ) + 1.687500e-03*aNodalField( 16 ) + 4.062500e-04*aNodalField( 17 ) + 5.868056e-04*aNodalField( 18 ) + 2.437500e-03*aNodalField( 19 ) + 5.868056e-04*aNodalField( 20 ) + 1.271412e-03*aNodalField( 21 ) + 5.281250e-03*aNodalField( 22 ) + 1.271412e-03*aNodalField( 23 ) + 9.780093e-05*aNodalField( 24 ) + 4.062500e-04*aNodalField( 25 ) + 9.780093e-05*aNodalField( 26 );
            aElementRHS( 11 ) = 4.513889e-05*aNodalField( 0 ) + 5.868056e-04*aNodalField( 1 ) + 2.708333e-04*aNodalField( 2 ) + 9.780093e-05*aNodalField( 3 ) + 1.271412e-03*aNodalField( 4 ) + 5.868056e-04*aNodalField( 5 ) + 7.523148e-06*aNodalField( 6 ) + 9.780093e-05*aNodalField( 7 ) + 4.513889e-05*aNodalField( 8 ) + 1.875000e-04*aNodalField(9) + 2.437500e-03*aNodalField( 10 ) + 1.125000e-03*aNodalField( 11 ) + 4.062500e-04*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 2.437500e-03*aNodalField( 14 ) + 3.125000e-05*aNodalField( 15 ) + 4.062500e-04*aNodalField( 16 ) + 1.875000e-04*aNodalField( 17 ) + 4.513889e-05*aNodalField( 18 ) + 5.868056e-04*aNodalField( 19 ) + 2.708333e-04*aNodalField( 20 ) + 9.780093e-05*aNodalField( 21 ) + 1.271412e-03*aNodalField( 22 ) + 5.868056e-04*aNodalField( 23 ) + 7.523148e-06*aNodalField( 24 ) + 9.780093e-05*aNodalField( 25 ) + 4.513889e-05*aNodalField( 26 );
            aElementRHS( 12 ) = 5.868056e-04*aNodalField( 0 ) + 1.271412e-03*aNodalField( 1 ) + 9.780093e-05*aNodalField( 2 ) + 2.437500e-03*aNodalField( 3 ) + 5.281250e-03*aNodalField( 4 ) + 4.062500e-04*aNodalField( 5 ) + 5.868056e-04*aNodalField( 6 ) + 1.271412e-03*aNodalField( 7 ) + 9.780093e-05*aNodalField( 8 ) + 2.437500e-03*aNodalField(9) + 5.281250e-03*aNodalField( 10 ) + 4.062500e-04*aNodalField( 11 ) + 1.012500e-02*aNodalField( 12 ) + 2.193750e-02*aNodalField( 13 ) + 1.687500e-03*aNodalField( 14 ) + 2.437500e-03*aNodalField( 15 ) + 5.281250e-03*aNodalField( 16 ) + 4.062500e-04*aNodalField( 17 ) + 5.868056e-04*aNodalField( 18 ) + 1.271412e-03*aNodalField( 19 ) + 9.780093e-05*aNodalField( 20 ) + 2.437500e-03*aNodalField( 21 ) + 5.281250e-03*aNodalField( 22 ) + 4.062500e-04*aNodalField( 23 ) + 5.868056e-04*aNodalField( 24 ) + 1.271412e-03*aNodalField( 25 ) + 9.780093e-05*aNodalField( 26 );
            aElementRHS( 13 ) = 1.271412e-03*aNodalField( 0 ) + 5.281250e-03*aNodalField( 1 ) + 1.271412e-03*aNodalField( 2 ) + 5.281250e-03*aNodalField( 3 ) + 2.193750e-02*aNodalField( 4 ) + 5.281250e-03*aNodalField( 5 ) + 1.271412e-03*aNodalField( 6 ) + 5.281250e-03*aNodalField( 7 ) + 1.271412e-03*aNodalField( 8 ) + 5.281250e-03*aNodalField(9) + 2.193750e-02*aNodalField( 10 ) + 5.281250e-03*aNodalField( 11 ) + 2.193750e-02*aNodalField( 12 ) + 9.112500e-02*aNodalField( 13 ) + 2.193750e-02*aNodalField( 14 ) + 5.281250e-03*aNodalField( 15 ) + 2.193750e-02*aNodalField( 16 ) + 5.281250e-03*aNodalField( 17 ) + 1.271412e-03*aNodalField( 18 ) + 5.281250e-03*aNodalField( 19 ) + 1.271412e-03*aNodalField( 20 ) + 5.281250e-03*aNodalField( 21 ) + 2.193750e-02*aNodalField( 22 ) + 5.281250e-03*aNodalField( 23 ) + 1.271412e-03*aNodalField( 24 ) + 5.281250e-03*aNodalField( 25 ) + 1.271412e-03*aNodalField( 26 );
            aElementRHS( 14 ) = 9.780093e-05*aNodalField( 0 ) + 1.271412e-03*aNodalField( 1 ) + 5.868056e-04*aNodalField( 2 ) + 4.062500e-04*aNodalField( 3 ) + 5.281250e-03*aNodalField( 4 ) + 2.437500e-03*aNodalField( 5 ) + 9.780093e-05*aNodalField( 6 ) + 1.271412e-03*aNodalField( 7 ) + 5.868056e-04*aNodalField( 8 ) + 4.062500e-04*aNodalField(9) + 5.281250e-03*aNodalField( 10 ) + 2.437500e-03*aNodalField( 11 ) + 1.687500e-03*aNodalField( 12 ) + 2.193750e-02*aNodalField( 13 ) + 1.012500e-02*aNodalField( 14 ) + 4.062500e-04*aNodalField( 15 ) + 5.281250e-03*aNodalField( 16 ) + 2.437500e-03*aNodalField( 17 ) + 9.780093e-05*aNodalField( 18 ) + 1.271412e-03*aNodalField( 19 ) + 5.868056e-04*aNodalField( 20 ) + 4.062500e-04*aNodalField( 21 ) + 5.281250e-03*aNodalField( 22 ) + 2.437500e-03*aNodalField( 23 ) + 9.780093e-05*aNodalField( 24 ) + 1.271412e-03*aNodalField( 25 ) + 5.868056e-04*aNodalField( 26 );
            aElementRHS( 15 ) = 4.513889e-05*aNodalField( 0 ) + 9.780093e-05*aNodalField( 1 ) + 7.523148e-06*aNodalField( 2 ) + 5.868056e-04*aNodalField( 3 ) + 1.271412e-03*aNodalField( 4 ) + 9.780093e-05*aNodalField( 5 ) + 2.708333e-04*aNodalField( 6 ) + 5.868056e-04*aNodalField( 7 ) + 4.513889e-05*aNodalField( 8 ) + 1.875000e-04*aNodalField(9) + 4.062500e-04*aNodalField( 10 ) + 3.125000e-05*aNodalField( 11 ) + 2.437500e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 4.062500e-04*aNodalField( 14 ) + 1.125000e-03*aNodalField( 15 ) + 2.437500e-03*aNodalField( 16 ) + 1.875000e-04*aNodalField( 17 ) + 4.513889e-05*aNodalField( 18 ) + 9.780093e-05*aNodalField( 19 ) + 7.523148e-06*aNodalField( 20 ) + 5.868056e-04*aNodalField( 21 ) + 1.271412e-03*aNodalField( 22 ) + 9.780093e-05*aNodalField( 23 ) + 2.708333e-04*aNodalField( 24 ) + 5.868056e-04*aNodalField( 25 ) + 4.513889e-05*aNodalField( 26 );
            aElementRHS( 16 ) = 9.780093e-05*aNodalField( 0 ) + 4.062500e-04*aNodalField( 1 ) + 9.780093e-05*aNodalField( 2 ) + 1.271412e-03*aNodalField( 3 ) + 5.281250e-03*aNodalField( 4 ) + 1.271412e-03*aNodalField( 5 ) + 5.868056e-04*aNodalField( 6 ) + 2.437500e-03*aNodalField( 7 ) + 5.868056e-04*aNodalField( 8 ) + 4.062500e-04*aNodalField(9) + 1.687500e-03*aNodalField( 10 ) + 4.062500e-04*aNodalField( 11 ) + 5.281250e-03*aNodalField( 12 ) + 2.193750e-02*aNodalField( 13 ) + 5.281250e-03*aNodalField( 14 ) + 2.437500e-03*aNodalField( 15 ) + 1.012500e-02*aNodalField( 16 ) + 2.437500e-03*aNodalField( 17 ) + 9.780093e-05*aNodalField( 18 ) + 4.062500e-04*aNodalField( 19 ) + 9.780093e-05*aNodalField( 20 ) + 1.271412e-03*aNodalField( 21 ) + 5.281250e-03*aNodalField( 22 ) + 1.271412e-03*aNodalField( 23 ) + 5.868056e-04*aNodalField( 24 ) + 2.437500e-03*aNodalField( 25 ) + 5.868056e-04*aNodalField( 26 );
            aElementRHS( 17 ) = 7.523148e-06*aNodalField( 0 ) + 9.780093e-05*aNodalField( 1 ) + 4.513889e-05*aNodalField( 2 ) + 9.780093e-05*aNodalField( 3 ) + 1.271412e-03*aNodalField( 4 ) + 5.868056e-04*aNodalField( 5 ) + 4.513889e-05*aNodalField( 6 ) + 5.868056e-04*aNodalField( 7 ) + 2.708333e-04*aNodalField( 8 ) + 3.125000e-05*aNodalField(9) + 4.062500e-04*aNodalField( 10 ) + 1.875000e-04*aNodalField( 11 ) + 4.062500e-04*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 2.437500e-03*aNodalField( 14 ) + 1.875000e-04*aNodalField( 15 ) + 2.437500e-03*aNodalField( 16 ) + 1.125000e-03*aNodalField( 17 ) + 7.523148e-06*aNodalField( 18 ) + 9.780093e-05*aNodalField( 19 ) + 4.513889e-05*aNodalField( 20 ) + 9.780093e-05*aNodalField( 21 ) + 1.271412e-03*aNodalField( 22 ) + 5.868056e-04*aNodalField( 23 ) + 4.513889e-05*aNodalField( 24 ) + 5.868056e-04*aNodalField( 25 ) + 2.708333e-04*aNodalField( 26 );
            aElementRHS( 18 ) = 2.083333e-05*aNodalField( 0 ) + 4.513889e-05*aNodalField( 1 ) + 3.472222e-06*aNodalField( 2 ) + 4.513889e-05*aNodalField( 3 ) + 9.780093e-05*aNodalField( 4 ) + 7.523148e-06*aNodalField( 5 ) + 3.472222e-06*aNodalField( 6 ) + 7.523148e-06*aNodalField( 7 ) + 5.787037e-07*aNodalField( 8 ) + 2.708333e-04*aNodalField(9) + 5.868056e-04*aNodalField( 10 ) + 4.513889e-05*aNodalField( 11 ) + 5.868056e-04*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 9.780093e-05*aNodalField( 14 ) + 4.513889e-05*aNodalField( 15 ) + 9.780093e-05*aNodalField( 16 ) + 7.523148e-06*aNodalField( 17 ) + 1.250000e-04*aNodalField( 18 ) + 2.708333e-04*aNodalField( 19 ) + 2.083333e-05*aNodalField( 20 ) + 2.708333e-04*aNodalField( 21 ) + 5.868056e-04*aNodalField( 22 ) + 4.513889e-05*aNodalField( 23 ) + 2.083333e-05*aNodalField( 24 ) + 4.513889e-05*aNodalField( 25 ) + 3.472222e-06*aNodalField( 26 );
            aElementRHS( 19 ) = 4.513889e-05*aNodalField( 0 ) + 1.875000e-04*aNodalField( 1 ) + 4.513889e-05*aNodalField( 2 ) + 9.780093e-05*aNodalField( 3 ) + 4.062500e-04*aNodalField( 4 ) + 9.780093e-05*aNodalField( 5 ) + 7.523148e-06*aNodalField( 6 ) + 3.125000e-05*aNodalField( 7 ) + 7.523148e-06*aNodalField( 8 ) + 5.868056e-04*aNodalField(9) + 2.437500e-03*aNodalField( 10 ) + 5.868056e-04*aNodalField( 11 ) + 1.271412e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 1.271412e-03*aNodalField( 14 ) + 9.780093e-05*aNodalField( 15 ) + 4.062500e-04*aNodalField( 16 ) + 9.780093e-05*aNodalField( 17 ) + 2.708333e-04*aNodalField( 18 ) + 1.125000e-03*aNodalField( 19 ) + 2.708333e-04*aNodalField( 20 ) + 5.868056e-04*aNodalField( 21 ) + 2.437500e-03*aNodalField( 22 ) + 5.868056e-04*aNodalField( 23 ) + 4.513889e-05*aNodalField( 24 ) + 1.875000e-04*aNodalField( 25 ) + 4.513889e-05*aNodalField( 26 );
            aElementRHS( 20 ) = 3.472222e-06*aNodalField( 0 ) + 4.513889e-05*aNodalField( 1 ) + 2.083333e-05*aNodalField( 2 ) + 7.523148e-06*aNodalField( 3 ) + 9.780093e-05*aNodalField( 4 ) + 4.513889e-05*aNodalField( 5 ) + 5.787037e-07*aNodalField( 6 ) + 7.523148e-06*aNodalField( 7 ) + 3.472222e-06*aNodalField( 8 ) + 4.513889e-05*aNodalField(9) + 5.868056e-04*aNodalField( 10 ) + 2.708333e-04*aNodalField( 11 ) + 9.780093e-05*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 5.868056e-04*aNodalField( 14 ) + 7.523148e-06*aNodalField( 15 ) + 9.780093e-05*aNodalField( 16 ) + 4.513889e-05*aNodalField( 17 ) + 2.083333e-05*aNodalField( 18 ) + 2.708333e-04*aNodalField( 19 ) + 1.250000e-04*aNodalField( 20 ) + 4.513889e-05*aNodalField( 21 ) + 5.868056e-04*aNodalField( 22 ) + 2.708333e-04*aNodalField( 23 ) + 3.472222e-06*aNodalField( 24 ) + 4.513889e-05*aNodalField( 25 ) + 2.083333e-05*aNodalField( 26 );
            aElementRHS( 21 ) = 4.513889e-05*aNodalField( 0 ) + 9.780093e-05*aNodalField( 1 ) + 7.523148e-06*aNodalField( 2 ) + 1.875000e-04*aNodalField( 3 ) + 4.062500e-04*aNodalField( 4 ) + 3.125000e-05*aNodalField( 5 ) + 4.513889e-05*aNodalField( 6 ) + 9.780093e-05*aNodalField( 7 ) + 7.523148e-06*aNodalField( 8 ) + 5.868056e-04*aNodalField(9) + 1.271412e-03*aNodalField( 10 ) + 9.780093e-05*aNodalField( 11 ) + 2.437500e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 4.062500e-04*aNodalField( 14 ) + 5.868056e-04*aNodalField( 15 ) + 1.271412e-03*aNodalField( 16 ) + 9.780093e-05*aNodalField( 17 ) + 2.708333e-04*aNodalField( 18 ) + 5.868056e-04*aNodalField( 19 ) + 4.513889e-05*aNodalField( 20 ) + 1.125000e-03*aNodalField( 21 ) + 2.437500e-03*aNodalField( 22 ) + 1.875000e-04*aNodalField( 23 ) + 2.708333e-04*aNodalField( 24 ) + 5.868056e-04*aNodalField( 25 ) + 4.513889e-05*aNodalField( 26 );
            aElementRHS( 22 ) = 9.780093e-05*aNodalField( 0 ) + 4.062500e-04*aNodalField( 1 ) + 9.780093e-05*aNodalField( 2 ) + 4.062500e-04*aNodalField( 3 ) + 1.687500e-03*aNodalField( 4 ) + 4.062500e-04*aNodalField( 5 ) + 9.780093e-05*aNodalField( 6 ) + 4.062500e-04*aNodalField( 7 ) + 9.780093e-05*aNodalField( 8 ) + 1.271412e-03*aNodalField(9) + 5.281250e-03*aNodalField( 10 ) + 1.271412e-03*aNodalField( 11 ) + 5.281250e-03*aNodalField( 12 ) + 2.193750e-02*aNodalField( 13 ) + 5.281250e-03*aNodalField( 14 ) + 1.271412e-03*aNodalField( 15 ) + 5.281250e-03*aNodalField( 16 ) + 1.271412e-03*aNodalField( 17 ) + 5.868056e-04*aNodalField( 18 ) + 2.437500e-03*aNodalField( 19 ) + 5.868056e-04*aNodalField( 20 ) + 2.437500e-03*aNodalField( 21 ) + 1.012500e-02*aNodalField( 22 ) + 2.437500e-03*aNodalField( 23 ) + 5.868056e-04*aNodalField( 24 ) + 2.437500e-03*aNodalField( 25 ) + 5.868056e-04*aNodalField( 26 );
            aElementRHS( 23 ) = 7.523148e-06*aNodalField( 0 ) + 9.780093e-05*aNodalField( 1 ) + 4.513889e-05*aNodalField( 2 ) + 3.125000e-05*aNodalField( 3 ) + 4.062500e-04*aNodalField( 4 ) + 1.875000e-04*aNodalField( 5 ) + 7.523148e-06*aNodalField( 6 ) + 9.780093e-05*aNodalField( 7 ) + 4.513889e-05*aNodalField( 8 ) + 9.780093e-05*aNodalField(9) + 1.271412e-03*aNodalField( 10 ) + 5.868056e-04*aNodalField( 11 ) + 4.062500e-04*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 2.437500e-03*aNodalField( 14 ) + 9.780093e-05*aNodalField( 15 ) + 1.271412e-03*aNodalField( 16 ) + 5.868056e-04*aNodalField( 17 ) + 4.513889e-05*aNodalField( 18 ) + 5.868056e-04*aNodalField( 19 ) + 2.708333e-04*aNodalField( 20 ) + 1.875000e-04*aNodalField( 21 ) + 2.437500e-03*aNodalField( 22 ) + 1.125000e-03*aNodalField( 23 ) + 4.513889e-05*aNodalField( 24 ) + 5.868056e-04*aNodalField( 25 ) + 2.708333e-04*aNodalField( 26 );
            aElementRHS( 24 ) = 3.472222e-06*aNodalField( 0 ) + 7.523148e-06*aNodalField( 1 ) + 5.787037e-07*aNodalField( 2 ) + 4.513889e-05*aNodalField( 3 ) + 9.780093e-05*aNodalField( 4 ) + 7.523148e-06*aNodalField( 5 ) + 2.083333e-05*aNodalField( 6 ) + 4.513889e-05*aNodalField( 7 ) + 3.472222e-06*aNodalField( 8 ) + 4.513889e-05*aNodalField(9) + 9.780093e-05*aNodalField( 10 ) + 7.523148e-06*aNodalField( 11 ) + 5.868056e-04*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 9.780093e-05*aNodalField( 14 ) + 2.708333e-04*aNodalField( 15 ) + 5.868056e-04*aNodalField( 16 ) + 4.513889e-05*aNodalField( 17 ) + 2.083333e-05*aNodalField( 18 ) + 4.513889e-05*aNodalField( 19 ) + 3.472222e-06*aNodalField( 20 ) + 2.708333e-04*aNodalField( 21 ) + 5.868056e-04*aNodalField( 22 ) + 4.513889e-05*aNodalField( 23 ) + 1.250000e-04*aNodalField( 24 ) + 2.708333e-04*aNodalField( 25 ) + 2.083333e-05*aNodalField( 26 );
            aElementRHS( 25 ) = 7.523148e-06*aNodalField( 0 ) + 3.125000e-05*aNodalField( 1 ) + 7.523148e-06*aNodalField( 2 ) + 9.780093e-05*aNodalField( 3 ) + 4.062500e-04*aNodalField( 4 ) + 9.780093e-05*aNodalField( 5 ) + 4.513889e-05*aNodalField( 6 ) + 1.875000e-04*aNodalField( 7 ) + 4.513889e-05*aNodalField( 8 ) + 9.780093e-05*aNodalField(9) + 4.062500e-04*aNodalField( 10 ) + 9.780093e-05*aNodalField( 11 ) + 1.271412e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 1.271412e-03*aNodalField( 14 ) + 5.868056e-04*aNodalField( 15 ) + 2.437500e-03*aNodalField( 16 ) + 5.868056e-04*aNodalField( 17 ) + 4.513889e-05*aNodalField( 18 ) + 1.875000e-04*aNodalField( 19 ) + 4.513889e-05*aNodalField( 20 ) + 5.868056e-04*aNodalField( 21 ) + 2.437500e-03*aNodalField( 22 ) + 5.868056e-04*aNodalField( 23 ) + 2.708333e-04*aNodalField( 24 ) + 1.125000e-03*aNodalField( 25 ) + 2.708333e-04*aNodalField( 26 );
            aElementRHS( 26 ) = 5.787037e-07*aNodalField( 0 ) + 7.523148e-06*aNodalField( 1 ) + 3.472222e-06*aNodalField( 2 ) + 7.523148e-06*aNodalField( 3 ) + 9.780093e-05*aNodalField( 4 ) + 4.513889e-05*aNodalField( 5 ) + 3.472222e-06*aNodalField( 6 ) + 4.513889e-05*aNodalField( 7 ) + 2.083333e-05*aNodalField( 8 ) + 7.523148e-06*aNodalField(9) + 9.780093e-05*aNodalField( 10 ) + 4.513889e-05*aNodalField( 11 ) + 9.780093e-05*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 5.868056e-04*aNodalField( 14 ) + 4.513889e-05*aNodalField( 15 ) + 5.868056e-04*aNodalField( 16 ) + 2.708333e-04*aNodalField( 17 ) + 3.472222e-06*aNodalField( 18 ) + 4.513889e-05*aNodalField( 19 ) + 2.083333e-05*aNodalField( 20 ) + 4.513889e-05*aNodalField( 21 ) + 5.868056e-04*aNodalField( 22 ) + 2.708333e-04*aNodalField( 23 ) + 2.083333e-05*aNodalField( 24 ) + 2.708333e-04*aNodalField( 25 ) + 1.250000e-04*aNodalField( 26 );

        }
    }
}

//-------------------------------------------------------------------------------------------------------------------------
Mat<real>
Hierarchical_Mesh_Main::RHS_for_L2_projection(
        const uint      & aDim,
        const uint      & aPolynomial,
        const Mat<real> & aNodalField)
{
    Mat<real> Element_RHS(pow(aPolynomial + 1,aDim),1);
    if ( aDim == 2)
    {
        if ( aPolynomial == 1)
        {
            Element_RHS( 0 ) = 1.0/36.0*(4.0*aNodalField( 0 ) + 2.0*aNodalField( 1 ) + 2.0*aNodalField( 2 ) + 1.0*aNodalField( 3 ));
            Element_RHS( 1 ) = 1.0/36.0*(2.0*aNodalField( 0 ) + 4.0*aNodalField( 1 ) + 1.0*aNodalField( 2 ) + 2.0*aNodalField( 3 ));
            Element_RHS( 2 ) = 1.0/36.0*(2.0*aNodalField( 0 ) + 1.0*aNodalField( 1 ) + 4.0*aNodalField( 2 ) + 2.0*aNodalField( 3 ));
            Element_RHS( 3 ) = 1.0/36.0*(1.0*aNodalField( 0 ) + 2.0*aNodalField( 1 ) + 2.0*aNodalField( 2 ) + 4.0*aNodalField( 3 ));
        }
        else if ( aPolynomial == 2)
        {
            Element_RHS( 0 ) = 0.002500000000000*aNodalField( 0 ) + 0.005416666666667*aNodalField( 1 ) + 0.000416666666667*aNodalField( 2 ) + 0.005416666666667*aNodalField( 3 ) + 0.011736111111111*aNodalField( 4 ) + 0.000902777777778*aNodalField( 5 ) + 0.000416666666667*aNodalField( 6 ) + 0.000902777777778*aNodalField( 7 ) + 0.000069444444444*aNodalField( 8 );
            Element_RHS( 1 ) = 0.005416666666667*aNodalField( 0 ) + 0.022500000000000*aNodalField( 1 ) + 0.005416666666667*aNodalField( 2 ) + 0.011736111111111*aNodalField( 3 ) + 0.048750000000000*aNodalField( 4 ) + 0.011736111111111*aNodalField( 5 ) + 0.000902777777778*aNodalField( 6 ) + 0.003750000000000*aNodalField( 7 ) + 0.000902777777778*aNodalField( 8 );
            Element_RHS( 2 ) = 0.000416666666667*aNodalField( 0 ) + 0.005416666666667*aNodalField( 1 ) + 0.002500000000000*aNodalField( 2 ) + 0.000902777777778*aNodalField( 3 ) + 0.011736111111111*aNodalField( 4 ) + 0.005416666666667*aNodalField( 5 ) + 0.000069444444444*aNodalField( 6 ) + 0.000902777777778*aNodalField( 7 ) + 0.000416666666667*aNodalField( 8 );
            Element_RHS( 3 ) = 0.005416666666667*aNodalField( 0 ) + 0.011736111111111*aNodalField( 1 ) + 0.000902777777778*aNodalField( 2 ) + 0.022500000000000*aNodalField( 3 ) + 0.048750000000000*aNodalField( 4 ) + 0.003750000000000*aNodalField( 5 ) + 0.005416666666667*aNodalField( 6 ) + 0.011736111111111*aNodalField( 7 ) + 0.000902777777778*aNodalField( 8 );
            Element_RHS( 4 ) = 0.011736111111111*aNodalField( 0 ) + 0.048750000000000*aNodalField( 1 ) + 0.011736111111111*aNodalField( 2 ) + 0.048750000000000*aNodalField( 3 ) + 0.202500000000000*aNodalField( 4 ) + 0.048750000000000*aNodalField( 5 ) + 0.011736111111111*aNodalField( 6 ) + 0.048750000000000*aNodalField( 7 ) + 0.011736111111111*aNodalField( 8 );
            Element_RHS( 5 ) = 0.000902777777778*aNodalField( 0 ) + 0.011736111111111*aNodalField( 1 ) + 0.005416666666667*aNodalField( 2 ) + 0.003750000000000*aNodalField( 3 ) + 0.048750000000000*aNodalField( 4 ) + 0.022500000000000*aNodalField( 5 ) + 0.000902777777778*aNodalField( 6 ) + 0.011736111111111*aNodalField( 7 ) + 0.005416666666667*aNodalField( 8 );
            Element_RHS( 6 ) = 0.000416666666667*aNodalField( 0 ) + 0.000902777777778*aNodalField( 1 ) + 0.000069444444444*aNodalField( 2 ) + 0.005416666666667*aNodalField( 3 ) + 0.011736111111111*aNodalField( 4 ) + 0.000902777777778*aNodalField( 5 ) + 0.002500000000000*aNodalField( 6 ) + 0.005416666666667*aNodalField( 7 ) + 0.000416666666667*aNodalField( 8 );
            Element_RHS( 7 ) = 0.000902777777778*aNodalField( 0 ) + 0.003750000000000*aNodalField( 1 ) + 0.000902777777778*aNodalField( 2 ) + 0.011736111111111*aNodalField( 3 ) + 0.048750000000000*aNodalField( 4 ) + 0.011736111111111*aNodalField( 5 ) + 0.005416666666667*aNodalField( 6 ) + 0.022500000000000*aNodalField( 7 ) + 0.005416666666667*aNodalField( 8 );
            Element_RHS( 8 ) = 0.000069444444444*aNodalField( 0 ) + 0.000902777777778*aNodalField( 1 ) + 0.000416666666667*aNodalField( 2 ) + 0.000902777777778*aNodalField( 3 ) + 0.011736111111111*aNodalField( 4 ) + 0.005416666666667*aNodalField( 5 ) + 0.000416666666667*aNodalField( 6 ) + 0.005416666666667*aNodalField( 7 ) + 0.002500000000000*aNodalField( 8 );
        }
    }
    else if ( aDim == 3)
    {
        if ( aPolynomial == 1)
        {
            Element_RHS( 0 ) = 1.0/216.0*(8.0*aNodalField( 0 ) + 4.0*aNodalField( 1 ) + 4.0*aNodalField( 2 ) + 2.0*aNodalField( 3 ) + 4.0*aNodalField( 4 ) + 2.0*aNodalField( 5 ) + 2.0*aNodalField( 6 ) + 1.0*aNodalField( 7 ));
            Element_RHS( 1 ) = 1.0/216.0*(4.0*aNodalField( 0 ) + 8.0*aNodalField( 1 ) + 2.0*aNodalField( 2 ) + 4.0*aNodalField( 3 ) + 2.0*aNodalField( 4 ) + 4.0*aNodalField( 5 ) + 1.0*aNodalField( 6 ) + 2.0*aNodalField( 7 ));
            Element_RHS( 2 ) = 1.0/216.0*(4.0*aNodalField( 0 ) + 2.0*aNodalField( 1 ) + 8.0*aNodalField( 2 ) + 4.0*aNodalField( 3 ) + 2.0*aNodalField( 4 ) + 1.0*aNodalField( 5 ) + 4.0*aNodalField( 6 ) + 2.0*aNodalField( 7 ));
            Element_RHS( 3 ) = 1.0/216.0*(2.0*aNodalField( 0 ) + 4.0*aNodalField( 1 ) + 4.0*aNodalField( 2 ) + 8.0*aNodalField( 3 ) + 1.0*aNodalField( 4 ) + 2.0*aNodalField( 5 ) + 2.0*aNodalField( 6 ) + 4.0*aNodalField( 7 ));
            Element_RHS( 4 ) = 1.0/216.0*(4.0*aNodalField( 0 ) + 2.0*aNodalField( 1 ) + 2.0*aNodalField( 2 ) + 1.0*aNodalField( 3 ) + 8.0*aNodalField( 4 ) + 4.0*aNodalField( 5 ) + 4.0*aNodalField( 6 ) + 2.0*aNodalField( 7 ));
            Element_RHS( 5 ) = 1.0/216.0*(2.0*aNodalField( 0 ) + 4.0*aNodalField( 1 ) + 1.0*aNodalField( 2 ) + 2.0*aNodalField( 3 ) + 4.0*aNodalField( 4 ) + 8.0*aNodalField( 5 ) + 2.0*aNodalField( 6 ) + 4.0*aNodalField( 7 ));
            Element_RHS( 6 ) = 1.0/216.0*(2.0*aNodalField( 0 ) + 1.0*aNodalField( 1 ) + 4.0*aNodalField( 2 ) + 2.0*aNodalField( 3 ) + 4.0*aNodalField( 4 ) + 2.0*aNodalField( 5 ) + 8.0*aNodalField( 6 ) + 4.0*aNodalField( 7 ));
            Element_RHS( 7 ) = 1.0/216.0*(1.0*aNodalField( 0 ) + 2.0*aNodalField( 1 ) + 2.0*aNodalField( 2 ) + 4.0*aNodalField( 3 ) + 2.0*aNodalField( 4 ) + 4.0*aNodalField( 5 ) + 4.0*aNodalField( 6 ) + 8.0*aNodalField( 7 ));
        }
        else if ( aPolynomial == 2)
        {
            Element_RHS( 0 ) = 1.250000e-04*aNodalField( 0 ) + 2.708333e-04*aNodalField( 1 ) + 2.083333e-05*aNodalField( 2 ) + 2.708333e-04*aNodalField( 3 ) + 5.868056e-04*aNodalField( 4 ) + 4.513889e-05*aNodalField( 5 ) + 2.083333e-05*aNodalField( 6 ) + 4.513889e-05*aNodalField( 7 ) + 3.472222e-06*aNodalField( 8 ) + 2.708333e-04*aNodalField(9) + 5.868056e-04*aNodalField( 10 ) + 4.513889e-05*aNodalField( 11 ) + 5.868056e-04*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 9.780093e-05*aNodalField( 14 ) + 4.513889e-05*aNodalField( 15 ) + 9.780093e-05*aNodalField( 16 ) + 7.523148e-06*aNodalField( 17 ) + 2.083333e-05*aNodalField( 18 ) + 4.513889e-05*aNodalField( 19 ) + 3.472222e-06*aNodalField( 20 ) + 4.513889e-05*aNodalField( 21 ) + 9.780093e-05*aNodalField( 22 ) + 7.523148e-06*aNodalField( 23 ) + 3.472222e-06*aNodalField( 24 ) + 7.523148e-06*aNodalField( 25 ) + 5.787037e-07*aNodalField( 26 );
            Element_RHS( 1 ) = 2.708333e-04*aNodalField( 0 ) + 1.125000e-03*aNodalField( 1 ) + 2.708333e-04*aNodalField( 2 ) + 5.868056e-04*aNodalField( 3 ) + 2.437500e-03*aNodalField( 4 ) + 5.868056e-04*aNodalField( 5 ) + 4.513889e-05*aNodalField( 6 ) + 1.875000e-04*aNodalField( 7 ) + 4.513889e-05*aNodalField( 8 ) + 5.868056e-04*aNodalField(9) + 2.437500e-03*aNodalField( 10 ) + 5.868056e-04*aNodalField( 11 ) + 1.271412e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 1.271412e-03*aNodalField( 14 ) + 9.780093e-05*aNodalField( 15 ) + 4.062500e-04*aNodalField( 16 ) + 9.780093e-05*aNodalField( 17 ) + 4.513889e-05*aNodalField( 18 ) + 1.875000e-04*aNodalField( 19 ) + 4.513889e-05*aNodalField( 20 ) + 9.780093e-05*aNodalField( 21 ) + 4.062500e-04*aNodalField( 22 ) + 9.780093e-05*aNodalField( 23 ) + 7.523148e-06*aNodalField( 24 ) + 3.125000e-05*aNodalField( 25 ) + 7.523148e-06*aNodalField( 26 );
            Element_RHS( 2 ) = 2.083333e-05*aNodalField( 0 ) + 2.708333e-04*aNodalField( 1 ) + 1.250000e-04*aNodalField( 2 ) + 4.513889e-05*aNodalField( 3 ) + 5.868056e-04*aNodalField( 4 ) + 2.708333e-04*aNodalField( 5 ) + 3.472222e-06*aNodalField( 6 ) + 4.513889e-05*aNodalField( 7 ) + 2.083333e-05*aNodalField( 8 ) + 4.513889e-05*aNodalField(9) + 5.868056e-04*aNodalField( 10 ) + 2.708333e-04*aNodalField( 11 ) + 9.780093e-05*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 5.868056e-04*aNodalField( 14 ) + 7.523148e-06*aNodalField( 15 ) + 9.780093e-05*aNodalField( 16 ) + 4.513889e-05*aNodalField( 17 ) + 3.472222e-06*aNodalField( 18 ) + 4.513889e-05*aNodalField( 19 ) + 2.083333e-05*aNodalField( 20 ) + 7.523148e-06*aNodalField( 21 ) + 9.780093e-05*aNodalField( 22 ) + 4.513889e-05*aNodalField( 23 ) + 5.787037e-07*aNodalField( 24 ) + 7.523148e-06*aNodalField( 25 ) + 3.472222e-06*aNodalField( 26 );
            Element_RHS( 3 ) = 2.708333e-04*aNodalField( 0 ) + 5.868056e-04*aNodalField( 1 ) + 4.513889e-05*aNodalField( 2 ) + 1.125000e-03*aNodalField( 3 ) + 2.437500e-03*aNodalField( 4 ) + 1.875000e-04*aNodalField( 5 ) + 2.708333e-04*aNodalField( 6 ) + 5.868056e-04*aNodalField( 7 ) + 4.513889e-05*aNodalField( 8 ) + 5.868056e-04*aNodalField(9) + 1.271412e-03*aNodalField( 10 ) + 9.780093e-05*aNodalField( 11 ) + 2.437500e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 4.062500e-04*aNodalField( 14 ) + 5.868056e-04*aNodalField( 15 ) + 1.271412e-03*aNodalField( 16 ) + 9.780093e-05*aNodalField( 17 ) + 4.513889e-05*aNodalField( 18 ) + 9.780093e-05*aNodalField( 19 ) + 7.523148e-06*aNodalField( 20 ) + 1.875000e-04*aNodalField( 21 ) + 4.062500e-04*aNodalField( 22 ) + 3.125000e-05*aNodalField( 23 ) + 4.513889e-05*aNodalField( 24 ) + 9.780093e-05*aNodalField( 25 ) + 7.523148e-06*aNodalField( 26 );
            Element_RHS( 4 ) = 5.868056e-04*aNodalField( 0 ) + 2.437500e-03*aNodalField( 1 ) + 5.868056e-04*aNodalField( 2 ) + 2.437500e-03*aNodalField( 3 ) + 1.012500e-02*aNodalField( 4 ) + 2.437500e-03*aNodalField( 5 ) + 5.868056e-04*aNodalField( 6 ) + 2.437500e-03*aNodalField( 7 ) + 5.868056e-04*aNodalField( 8 ) + 1.271412e-03*aNodalField(9) + 5.281250e-03*aNodalField( 10 ) + 1.271412e-03*aNodalField( 11 ) + 5.281250e-03*aNodalField( 12 ) + 2.193750e-02*aNodalField( 13 ) + 5.281250e-03*aNodalField( 14 ) + 1.271412e-03*aNodalField( 15 ) + 5.281250e-03*aNodalField( 16 ) + 1.271412e-03*aNodalField( 17 ) + 9.780093e-05*aNodalField( 18 ) + 4.062500e-04*aNodalField( 19 ) + 9.780093e-05*aNodalField( 20 ) + 4.062500e-04*aNodalField( 21 ) + 1.687500e-03*aNodalField( 22 ) + 4.062500e-04*aNodalField( 23 ) + 9.780093e-05*aNodalField( 24 ) + 4.062500e-04*aNodalField( 25 ) + 9.780093e-05*aNodalField( 26 );
            Element_RHS( 5 ) = 4.513889e-05*aNodalField( 0 ) + 5.868056e-04*aNodalField( 1 ) + 2.708333e-04*aNodalField( 2 ) + 1.875000e-04*aNodalField( 3 ) + 2.437500e-03*aNodalField( 4 ) + 1.125000e-03*aNodalField( 5 ) + 4.513889e-05*aNodalField( 6 ) + 5.868056e-04*aNodalField( 7 ) + 2.708333e-04*aNodalField( 8 ) + 9.780093e-05*aNodalField(9) + 1.271412e-03*aNodalField( 10 ) + 5.868056e-04*aNodalField( 11 ) + 4.062500e-04*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 2.437500e-03*aNodalField( 14 ) + 9.780093e-05*aNodalField( 15 ) + 1.271412e-03*aNodalField( 16 ) + 5.868056e-04*aNodalField( 17 ) + 7.523148e-06*aNodalField( 18 ) + 9.780093e-05*aNodalField( 19 ) + 4.513889e-05*aNodalField( 20 ) + 3.125000e-05*aNodalField( 21 ) + 4.062500e-04*aNodalField( 22 ) + 1.875000e-04*aNodalField( 23 ) + 7.523148e-06*aNodalField( 24 ) + 9.780093e-05*aNodalField( 25 ) + 4.513889e-05*aNodalField( 26 );
            Element_RHS( 6 ) = 2.083333e-05*aNodalField( 0 ) + 4.513889e-05*aNodalField( 1 ) + 3.472222e-06*aNodalField( 2 ) + 2.708333e-04*aNodalField( 3 ) + 5.868056e-04*aNodalField( 4 ) + 4.513889e-05*aNodalField( 5 ) + 1.250000e-04*aNodalField( 6 ) + 2.708333e-04*aNodalField( 7 ) + 2.083333e-05*aNodalField( 8 ) + 4.513889e-05*aNodalField(9) + 9.780093e-05*aNodalField( 10 ) + 7.523148e-06*aNodalField( 11 ) + 5.868056e-04*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 9.780093e-05*aNodalField( 14 ) + 2.708333e-04*aNodalField( 15 ) + 5.868056e-04*aNodalField( 16 ) + 4.513889e-05*aNodalField( 17 ) + 3.472222e-06*aNodalField( 18 ) + 7.523148e-06*aNodalField( 19 ) + 5.787037e-07*aNodalField( 20 ) + 4.513889e-05*aNodalField( 21 ) + 9.780093e-05*aNodalField( 22 ) + 7.523148e-06*aNodalField( 23 ) + 2.083333e-05*aNodalField( 24 ) + 4.513889e-05*aNodalField( 25 ) + 3.472222e-06*aNodalField( 26 );
            Element_RHS( 7 ) = 4.513889e-05*aNodalField( 0 ) + 1.875000e-04*aNodalField( 1 ) + 4.513889e-05*aNodalField( 2 ) + 5.868056e-04*aNodalField( 3 ) + 2.437500e-03*aNodalField( 4 ) + 5.868056e-04*aNodalField( 5 ) + 2.708333e-04*aNodalField( 6 ) + 1.125000e-03*aNodalField( 7 ) + 2.708333e-04*aNodalField( 8 ) + 9.780093e-05*aNodalField(9) + 4.062500e-04*aNodalField( 10 ) + 9.780093e-05*aNodalField( 11 ) + 1.271412e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 1.271412e-03*aNodalField( 14 ) + 5.868056e-04*aNodalField( 15 ) + 2.437500e-03*aNodalField( 16 ) + 5.868056e-04*aNodalField( 17 ) + 7.523148e-06*aNodalField( 18 ) + 3.125000e-05*aNodalField( 19 ) + 7.523148e-06*aNodalField( 20 ) + 9.780093e-05*aNodalField( 21 ) + 4.062500e-04*aNodalField( 22 ) + 9.780093e-05*aNodalField( 23 ) + 4.513889e-05*aNodalField( 24 ) + 1.875000e-04*aNodalField( 25 ) + 4.513889e-05*aNodalField( 26 );
            Element_RHS( 8 ) = 3.472222e-06*aNodalField( 0 ) + 4.513889e-05*aNodalField( 1 ) + 2.083333e-05*aNodalField( 2 ) + 4.513889e-05*aNodalField( 3 ) + 5.868056e-04*aNodalField( 4 ) + 2.708333e-04*aNodalField( 5 ) + 2.083333e-05*aNodalField( 6 ) + 2.708333e-04*aNodalField( 7 ) + 1.250000e-04*aNodalField( 8 ) + 7.523148e-06*aNodalField(9) + 9.780093e-05*aNodalField( 10 ) + 4.513889e-05*aNodalField( 11 ) + 9.780093e-05*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 5.868056e-04*aNodalField( 14 ) + 4.513889e-05*aNodalField( 15 ) + 5.868056e-04*aNodalField( 16 ) + 2.708333e-04*aNodalField( 17 ) + 5.787037e-07*aNodalField( 18 ) + 7.523148e-06*aNodalField( 19 ) + 3.472222e-06*aNodalField( 20 ) + 7.523148e-06*aNodalField( 21 ) + 9.780093e-05*aNodalField( 22 ) + 4.513889e-05*aNodalField( 23 ) + 3.472222e-06*aNodalField( 24 ) + 4.513889e-05*aNodalField( 25 ) + 2.083333e-05*aNodalField( 26 );
            Element_RHS(9) = 2.708333e-04*aNodalField( 0 ) + 5.868056e-04*aNodalField( 1 ) + 4.513889e-05*aNodalField( 2 ) + 5.868056e-04*aNodalField( 3 ) + 1.271412e-03*aNodalField( 4 ) + 9.780093e-05*aNodalField( 5 ) + 4.513889e-05*aNodalField( 6 ) + 9.780093e-05*aNodalField( 7 ) + 7.523148e-06*aNodalField( 8 ) + 1.125000e-03*aNodalField(9) + 2.437500e-03*aNodalField( 10 ) + 1.875000e-04*aNodalField( 11 ) + 2.437500e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 4.062500e-04*aNodalField( 14 ) + 1.875000e-04*aNodalField( 15 ) + 4.062500e-04*aNodalField( 16 ) + 3.125000e-05*aNodalField( 17 ) + 2.708333e-04*aNodalField( 18 ) + 5.868056e-04*aNodalField( 19 ) + 4.513889e-05*aNodalField( 20 ) + 5.868056e-04*aNodalField( 21 ) + 1.271412e-03*aNodalField( 22 ) + 9.780093e-05*aNodalField( 23 ) + 4.513889e-05*aNodalField( 24 ) + 9.780093e-05*aNodalField( 25 ) + 7.523148e-06*aNodalField( 26 );
            Element_RHS( 10 ) = 5.868056e-04*aNodalField( 0 ) + 2.437500e-03*aNodalField( 1 ) + 5.868056e-04*aNodalField( 2 ) + 1.271412e-03*aNodalField( 3 ) + 5.281250e-03*aNodalField( 4 ) + 1.271412e-03*aNodalField( 5 ) + 9.780093e-05*aNodalField( 6 ) + 4.062500e-04*aNodalField( 7 ) + 9.780093e-05*aNodalField( 8 ) + 2.437500e-03*aNodalField(9) + 1.012500e-02*aNodalField( 10 ) + 2.437500e-03*aNodalField( 11 ) + 5.281250e-03*aNodalField( 12 ) + 2.193750e-02*aNodalField( 13 ) + 5.281250e-03*aNodalField( 14 ) + 4.062500e-04*aNodalField( 15 ) + 1.687500e-03*aNodalField( 16 ) + 4.062500e-04*aNodalField( 17 ) + 5.868056e-04*aNodalField( 18 ) + 2.437500e-03*aNodalField( 19 ) + 5.868056e-04*aNodalField( 20 ) + 1.271412e-03*aNodalField( 21 ) + 5.281250e-03*aNodalField( 22 ) + 1.271412e-03*aNodalField( 23 ) + 9.780093e-05*aNodalField( 24 ) + 4.062500e-04*aNodalField( 25 ) + 9.780093e-05*aNodalField( 26 );
            Element_RHS( 11 ) = 4.513889e-05*aNodalField( 0 ) + 5.868056e-04*aNodalField( 1 ) + 2.708333e-04*aNodalField( 2 ) + 9.780093e-05*aNodalField( 3 ) + 1.271412e-03*aNodalField( 4 ) + 5.868056e-04*aNodalField( 5 ) + 7.523148e-06*aNodalField( 6 ) + 9.780093e-05*aNodalField( 7 ) + 4.513889e-05*aNodalField( 8 ) + 1.875000e-04*aNodalField(9) + 2.437500e-03*aNodalField( 10 ) + 1.125000e-03*aNodalField( 11 ) + 4.062500e-04*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 2.437500e-03*aNodalField( 14 ) + 3.125000e-05*aNodalField( 15 ) + 4.062500e-04*aNodalField( 16 ) + 1.875000e-04*aNodalField( 17 ) + 4.513889e-05*aNodalField( 18 ) + 5.868056e-04*aNodalField( 19 ) + 2.708333e-04*aNodalField( 20 ) + 9.780093e-05*aNodalField( 21 ) + 1.271412e-03*aNodalField( 22 ) + 5.868056e-04*aNodalField( 23 ) + 7.523148e-06*aNodalField( 24 ) + 9.780093e-05*aNodalField( 25 ) + 4.513889e-05*aNodalField( 26 );
            Element_RHS( 12 ) = 5.868056e-04*aNodalField( 0 ) + 1.271412e-03*aNodalField( 1 ) + 9.780093e-05*aNodalField( 2 ) + 2.437500e-03*aNodalField( 3 ) + 5.281250e-03*aNodalField( 4 ) + 4.062500e-04*aNodalField( 5 ) + 5.868056e-04*aNodalField( 6 ) + 1.271412e-03*aNodalField( 7 ) + 9.780093e-05*aNodalField( 8 ) + 2.437500e-03*aNodalField(9) + 5.281250e-03*aNodalField( 10 ) + 4.062500e-04*aNodalField( 11 ) + 1.012500e-02*aNodalField( 12 ) + 2.193750e-02*aNodalField( 13 ) + 1.687500e-03*aNodalField( 14 ) + 2.437500e-03*aNodalField( 15 ) + 5.281250e-03*aNodalField( 16 ) + 4.062500e-04*aNodalField( 17 ) + 5.868056e-04*aNodalField( 18 ) + 1.271412e-03*aNodalField( 19 ) + 9.780093e-05*aNodalField( 20 ) + 2.437500e-03*aNodalField( 21 ) + 5.281250e-03*aNodalField( 22 ) + 4.062500e-04*aNodalField( 23 ) + 5.868056e-04*aNodalField( 24 ) + 1.271412e-03*aNodalField( 25 ) + 9.780093e-05*aNodalField( 26 );
            Element_RHS( 13 ) = 1.271412e-03*aNodalField( 0 ) + 5.281250e-03*aNodalField( 1 ) + 1.271412e-03*aNodalField( 2 ) + 5.281250e-03*aNodalField( 3 ) + 2.193750e-02*aNodalField( 4 ) + 5.281250e-03*aNodalField( 5 ) + 1.271412e-03*aNodalField( 6 ) + 5.281250e-03*aNodalField( 7 ) + 1.271412e-03*aNodalField( 8 ) + 5.281250e-03*aNodalField(9) + 2.193750e-02*aNodalField( 10 ) + 5.281250e-03*aNodalField( 11 ) + 2.193750e-02*aNodalField( 12 ) + 9.112500e-02*aNodalField( 13 ) + 2.193750e-02*aNodalField( 14 ) + 5.281250e-03*aNodalField( 15 ) + 2.193750e-02*aNodalField( 16 ) + 5.281250e-03*aNodalField( 17 ) + 1.271412e-03*aNodalField( 18 ) + 5.281250e-03*aNodalField( 19 ) + 1.271412e-03*aNodalField( 20 ) + 5.281250e-03*aNodalField( 21 ) + 2.193750e-02*aNodalField( 22 ) + 5.281250e-03*aNodalField( 23 ) + 1.271412e-03*aNodalField( 24 ) + 5.281250e-03*aNodalField( 25 ) + 1.271412e-03*aNodalField( 26 );
            Element_RHS( 14 ) = 9.780093e-05*aNodalField( 0 ) + 1.271412e-03*aNodalField( 1 ) + 5.868056e-04*aNodalField( 2 ) + 4.062500e-04*aNodalField( 3 ) + 5.281250e-03*aNodalField( 4 ) + 2.437500e-03*aNodalField( 5 ) + 9.780093e-05*aNodalField( 6 ) + 1.271412e-03*aNodalField( 7 ) + 5.868056e-04*aNodalField( 8 ) + 4.062500e-04*aNodalField(9) + 5.281250e-03*aNodalField( 10 ) + 2.437500e-03*aNodalField( 11 ) + 1.687500e-03*aNodalField( 12 ) + 2.193750e-02*aNodalField( 13 ) + 1.012500e-02*aNodalField( 14 ) + 4.062500e-04*aNodalField( 15 ) + 5.281250e-03*aNodalField( 16 ) + 2.437500e-03*aNodalField( 17 ) + 9.780093e-05*aNodalField( 18 ) + 1.271412e-03*aNodalField( 19 ) + 5.868056e-04*aNodalField( 20 ) + 4.062500e-04*aNodalField( 21 ) + 5.281250e-03*aNodalField( 22 ) + 2.437500e-03*aNodalField( 23 ) + 9.780093e-05*aNodalField( 24 ) + 1.271412e-03*aNodalField( 25 ) + 5.868056e-04*aNodalField( 26 );
            Element_RHS( 15 ) = 4.513889e-05*aNodalField( 0 ) + 9.780093e-05*aNodalField( 1 ) + 7.523148e-06*aNodalField( 2 ) + 5.868056e-04*aNodalField( 3 ) + 1.271412e-03*aNodalField( 4 ) + 9.780093e-05*aNodalField( 5 ) + 2.708333e-04*aNodalField( 6 ) + 5.868056e-04*aNodalField( 7 ) + 4.513889e-05*aNodalField( 8 ) + 1.875000e-04*aNodalField(9) + 4.062500e-04*aNodalField( 10 ) + 3.125000e-05*aNodalField( 11 ) + 2.437500e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 4.062500e-04*aNodalField( 14 ) + 1.125000e-03*aNodalField( 15 ) + 2.437500e-03*aNodalField( 16 ) + 1.875000e-04*aNodalField( 17 ) + 4.513889e-05*aNodalField( 18 ) + 9.780093e-05*aNodalField( 19 ) + 7.523148e-06*aNodalField( 20 ) + 5.868056e-04*aNodalField( 21 ) + 1.271412e-03*aNodalField( 22 ) + 9.780093e-05*aNodalField( 23 ) + 2.708333e-04*aNodalField( 24 ) + 5.868056e-04*aNodalField( 25 ) + 4.513889e-05*aNodalField( 26 );
            Element_RHS( 16 ) = 9.780093e-05*aNodalField( 0 ) + 4.062500e-04*aNodalField( 1 ) + 9.780093e-05*aNodalField( 2 ) + 1.271412e-03*aNodalField( 3 ) + 5.281250e-03*aNodalField( 4 ) + 1.271412e-03*aNodalField( 5 ) + 5.868056e-04*aNodalField( 6 ) + 2.437500e-03*aNodalField( 7 ) + 5.868056e-04*aNodalField( 8 ) + 4.062500e-04*aNodalField(9) + 1.687500e-03*aNodalField( 10 ) + 4.062500e-04*aNodalField( 11 ) + 5.281250e-03*aNodalField( 12 ) + 2.193750e-02*aNodalField( 13 ) + 5.281250e-03*aNodalField( 14 ) + 2.437500e-03*aNodalField( 15 ) + 1.012500e-02*aNodalField( 16 ) + 2.437500e-03*aNodalField( 17 ) + 9.780093e-05*aNodalField( 18 ) + 4.062500e-04*aNodalField( 19 ) + 9.780093e-05*aNodalField( 20 ) + 1.271412e-03*aNodalField( 21 ) + 5.281250e-03*aNodalField( 22 ) + 1.271412e-03*aNodalField( 23 ) + 5.868056e-04*aNodalField( 24 ) + 2.437500e-03*aNodalField( 25 ) + 5.868056e-04*aNodalField( 26 );
            Element_RHS( 17 ) = 7.523148e-06*aNodalField( 0 ) + 9.780093e-05*aNodalField( 1 ) + 4.513889e-05*aNodalField( 2 ) + 9.780093e-05*aNodalField( 3 ) + 1.271412e-03*aNodalField( 4 ) + 5.868056e-04*aNodalField( 5 ) + 4.513889e-05*aNodalField( 6 ) + 5.868056e-04*aNodalField( 7 ) + 2.708333e-04*aNodalField( 8 ) + 3.125000e-05*aNodalField(9) + 4.062500e-04*aNodalField( 10 ) + 1.875000e-04*aNodalField( 11 ) + 4.062500e-04*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 2.437500e-03*aNodalField( 14 ) + 1.875000e-04*aNodalField( 15 ) + 2.437500e-03*aNodalField( 16 ) + 1.125000e-03*aNodalField( 17 ) + 7.523148e-06*aNodalField( 18 ) + 9.780093e-05*aNodalField( 19 ) + 4.513889e-05*aNodalField( 20 ) + 9.780093e-05*aNodalField( 21 ) + 1.271412e-03*aNodalField( 22 ) + 5.868056e-04*aNodalField( 23 ) + 4.513889e-05*aNodalField( 24 ) + 5.868056e-04*aNodalField( 25 ) + 2.708333e-04*aNodalField( 26 );
            Element_RHS( 18 ) = 2.083333e-05*aNodalField( 0 ) + 4.513889e-05*aNodalField( 1 ) + 3.472222e-06*aNodalField( 2 ) + 4.513889e-05*aNodalField( 3 ) + 9.780093e-05*aNodalField( 4 ) + 7.523148e-06*aNodalField( 5 ) + 3.472222e-06*aNodalField( 6 ) + 7.523148e-06*aNodalField( 7 ) + 5.787037e-07*aNodalField( 8 ) + 2.708333e-04*aNodalField(9) + 5.868056e-04*aNodalField( 10 ) + 4.513889e-05*aNodalField( 11 ) + 5.868056e-04*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 9.780093e-05*aNodalField( 14 ) + 4.513889e-05*aNodalField( 15 ) + 9.780093e-05*aNodalField( 16 ) + 7.523148e-06*aNodalField( 17 ) + 1.250000e-04*aNodalField( 18 ) + 2.708333e-04*aNodalField( 19 ) + 2.083333e-05*aNodalField( 20 ) + 2.708333e-04*aNodalField( 21 ) + 5.868056e-04*aNodalField( 22 ) + 4.513889e-05*aNodalField( 23 ) + 2.083333e-05*aNodalField( 24 ) + 4.513889e-05*aNodalField( 25 ) + 3.472222e-06*aNodalField( 26 );
            Element_RHS( 19 ) = 4.513889e-05*aNodalField( 0 ) + 1.875000e-04*aNodalField( 1 ) + 4.513889e-05*aNodalField( 2 ) + 9.780093e-05*aNodalField( 3 ) + 4.062500e-04*aNodalField( 4 ) + 9.780093e-05*aNodalField( 5 ) + 7.523148e-06*aNodalField( 6 ) + 3.125000e-05*aNodalField( 7 ) + 7.523148e-06*aNodalField( 8 ) + 5.868056e-04*aNodalField(9) + 2.437500e-03*aNodalField( 10 ) + 5.868056e-04*aNodalField( 11 ) + 1.271412e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 1.271412e-03*aNodalField( 14 ) + 9.780093e-05*aNodalField( 15 ) + 4.062500e-04*aNodalField( 16 ) + 9.780093e-05*aNodalField( 17 ) + 2.708333e-04*aNodalField( 18 ) + 1.125000e-03*aNodalField( 19 ) + 2.708333e-04*aNodalField( 20 ) + 5.868056e-04*aNodalField( 21 ) + 2.437500e-03*aNodalField( 22 ) + 5.868056e-04*aNodalField( 23 ) + 4.513889e-05*aNodalField( 24 ) + 1.875000e-04*aNodalField( 25 ) + 4.513889e-05*aNodalField( 26 );
            Element_RHS( 20 ) = 3.472222e-06*aNodalField( 0 ) + 4.513889e-05*aNodalField( 1 ) + 2.083333e-05*aNodalField( 2 ) + 7.523148e-06*aNodalField( 3 ) + 9.780093e-05*aNodalField( 4 ) + 4.513889e-05*aNodalField( 5 ) + 5.787037e-07*aNodalField( 6 ) + 7.523148e-06*aNodalField( 7 ) + 3.472222e-06*aNodalField( 8 ) + 4.513889e-05*aNodalField(9) + 5.868056e-04*aNodalField( 10 ) + 2.708333e-04*aNodalField( 11 ) + 9.780093e-05*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 5.868056e-04*aNodalField( 14 ) + 7.523148e-06*aNodalField( 15 ) + 9.780093e-05*aNodalField( 16 ) + 4.513889e-05*aNodalField( 17 ) + 2.083333e-05*aNodalField( 18 ) + 2.708333e-04*aNodalField( 19 ) + 1.250000e-04*aNodalField( 20 ) + 4.513889e-05*aNodalField( 21 ) + 5.868056e-04*aNodalField( 22 ) + 2.708333e-04*aNodalField( 23 ) + 3.472222e-06*aNodalField( 24 ) + 4.513889e-05*aNodalField( 25 ) + 2.083333e-05*aNodalField( 26 );
            Element_RHS( 21 ) = 4.513889e-05*aNodalField( 0 ) + 9.780093e-05*aNodalField( 1 ) + 7.523148e-06*aNodalField( 2 ) + 1.875000e-04*aNodalField( 3 ) + 4.062500e-04*aNodalField( 4 ) + 3.125000e-05*aNodalField( 5 ) + 4.513889e-05*aNodalField( 6 ) + 9.780093e-05*aNodalField( 7 ) + 7.523148e-06*aNodalField( 8 ) + 5.868056e-04*aNodalField(9) + 1.271412e-03*aNodalField( 10 ) + 9.780093e-05*aNodalField( 11 ) + 2.437500e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 4.062500e-04*aNodalField( 14 ) + 5.868056e-04*aNodalField( 15 ) + 1.271412e-03*aNodalField( 16 ) + 9.780093e-05*aNodalField( 17 ) + 2.708333e-04*aNodalField( 18 ) + 5.868056e-04*aNodalField( 19 ) + 4.513889e-05*aNodalField( 20 ) + 1.125000e-03*aNodalField( 21 ) + 2.437500e-03*aNodalField( 22 ) + 1.875000e-04*aNodalField( 23 ) + 2.708333e-04*aNodalField( 24 ) + 5.868056e-04*aNodalField( 25 ) + 4.513889e-05*aNodalField( 26 );
            Element_RHS( 22 ) = 9.780093e-05*aNodalField( 0 ) + 4.062500e-04*aNodalField( 1 ) + 9.780093e-05*aNodalField( 2 ) + 4.062500e-04*aNodalField( 3 ) + 1.687500e-03*aNodalField( 4 ) + 4.062500e-04*aNodalField( 5 ) + 9.780093e-05*aNodalField( 6 ) + 4.062500e-04*aNodalField( 7 ) + 9.780093e-05*aNodalField( 8 ) + 1.271412e-03*aNodalField(9) + 5.281250e-03*aNodalField( 10 ) + 1.271412e-03*aNodalField( 11 ) + 5.281250e-03*aNodalField( 12 ) + 2.193750e-02*aNodalField( 13 ) + 5.281250e-03*aNodalField( 14 ) + 1.271412e-03*aNodalField( 15 ) + 5.281250e-03*aNodalField( 16 ) + 1.271412e-03*aNodalField( 17 ) + 5.868056e-04*aNodalField( 18 ) + 2.437500e-03*aNodalField( 19 ) + 5.868056e-04*aNodalField( 20 ) + 2.437500e-03*aNodalField( 21 ) + 1.012500e-02*aNodalField( 22 ) + 2.437500e-03*aNodalField( 23 ) + 5.868056e-04*aNodalField( 24 ) + 2.437500e-03*aNodalField( 25 ) + 5.868056e-04*aNodalField( 26 );
            Element_RHS( 23 ) = 7.523148e-06*aNodalField( 0 ) + 9.780093e-05*aNodalField( 1 ) + 4.513889e-05*aNodalField( 2 ) + 3.125000e-05*aNodalField( 3 ) + 4.062500e-04*aNodalField( 4 ) + 1.875000e-04*aNodalField( 5 ) + 7.523148e-06*aNodalField( 6 ) + 9.780093e-05*aNodalField( 7 ) + 4.513889e-05*aNodalField( 8 ) + 9.780093e-05*aNodalField(9) + 1.271412e-03*aNodalField( 10 ) + 5.868056e-04*aNodalField( 11 ) + 4.062500e-04*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 2.437500e-03*aNodalField( 14 ) + 9.780093e-05*aNodalField( 15 ) + 1.271412e-03*aNodalField( 16 ) + 5.868056e-04*aNodalField( 17 ) + 4.513889e-05*aNodalField( 18 ) + 5.868056e-04*aNodalField( 19 ) + 2.708333e-04*aNodalField( 20 ) + 1.875000e-04*aNodalField( 21 ) + 2.437500e-03*aNodalField( 22 ) + 1.125000e-03*aNodalField( 23 ) + 4.513889e-05*aNodalField( 24 ) + 5.868056e-04*aNodalField( 25 ) + 2.708333e-04*aNodalField( 26 );
            Element_RHS( 24 ) = 3.472222e-06*aNodalField( 0 ) + 7.523148e-06*aNodalField( 1 ) + 5.787037e-07*aNodalField( 2 ) + 4.513889e-05*aNodalField( 3 ) + 9.780093e-05*aNodalField( 4 ) + 7.523148e-06*aNodalField( 5 ) + 2.083333e-05*aNodalField( 6 ) + 4.513889e-05*aNodalField( 7 ) + 3.472222e-06*aNodalField( 8 ) + 4.513889e-05*aNodalField(9) + 9.780093e-05*aNodalField( 10 ) + 7.523148e-06*aNodalField( 11 ) + 5.868056e-04*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 9.780093e-05*aNodalField( 14 ) + 2.708333e-04*aNodalField( 15 ) + 5.868056e-04*aNodalField( 16 ) + 4.513889e-05*aNodalField( 17 ) + 2.083333e-05*aNodalField( 18 ) + 4.513889e-05*aNodalField( 19 ) + 3.472222e-06*aNodalField( 20 ) + 2.708333e-04*aNodalField( 21 ) + 5.868056e-04*aNodalField( 22 ) + 4.513889e-05*aNodalField( 23 ) + 1.250000e-04*aNodalField( 24 ) + 2.708333e-04*aNodalField( 25 ) + 2.083333e-05*aNodalField( 26 );
            Element_RHS( 25 ) = 7.523148e-06*aNodalField( 0 ) + 3.125000e-05*aNodalField( 1 ) + 7.523148e-06*aNodalField( 2 ) + 9.780093e-05*aNodalField( 3 ) + 4.062500e-04*aNodalField( 4 ) + 9.780093e-05*aNodalField( 5 ) + 4.513889e-05*aNodalField( 6 ) + 1.875000e-04*aNodalField( 7 ) + 4.513889e-05*aNodalField( 8 ) + 9.780093e-05*aNodalField(9) + 4.062500e-04*aNodalField( 10 ) + 9.780093e-05*aNodalField( 11 ) + 1.271412e-03*aNodalField( 12 ) + 5.281250e-03*aNodalField( 13 ) + 1.271412e-03*aNodalField( 14 ) + 5.868056e-04*aNodalField( 15 ) + 2.437500e-03*aNodalField( 16 ) + 5.868056e-04*aNodalField( 17 ) + 4.513889e-05*aNodalField( 18 ) + 1.875000e-04*aNodalField( 19 ) + 4.513889e-05*aNodalField( 20 ) + 5.868056e-04*aNodalField( 21 ) + 2.437500e-03*aNodalField( 22 ) + 5.868056e-04*aNodalField( 23 ) + 2.708333e-04*aNodalField( 24 ) + 1.125000e-03*aNodalField( 25 ) + 2.708333e-04*aNodalField( 26 );
            Element_RHS( 26 ) = 5.787037e-07*aNodalField( 0 ) + 7.523148e-06*aNodalField( 1 ) + 3.472222e-06*aNodalField( 2 ) + 7.523148e-06*aNodalField( 3 ) + 9.780093e-05*aNodalField( 4 ) + 4.513889e-05*aNodalField( 5 ) + 3.472222e-06*aNodalField( 6 ) + 4.513889e-05*aNodalField( 7 ) + 2.083333e-05*aNodalField( 8 ) + 7.523148e-06*aNodalField(9) + 9.780093e-05*aNodalField( 10 ) + 4.513889e-05*aNodalField( 11 ) + 9.780093e-05*aNodalField( 12 ) + 1.271412e-03*aNodalField( 13 ) + 5.868056e-04*aNodalField( 14 ) + 4.513889e-05*aNodalField( 15 ) + 5.868056e-04*aNodalField( 16 ) + 2.708333e-04*aNodalField( 17 ) + 3.472222e-06*aNodalField( 18 ) + 4.513889e-05*aNodalField( 19 ) + 2.083333e-05*aNodalField( 20 ) + 4.513889e-05*aNodalField( 21 ) + 5.868056e-04*aNodalField( 22 ) + 2.708333e-04*aNodalField( 23 ) + 2.083333e-05*aNodalField( 24 ) + 2.708333e-04*aNodalField( 25 ) + 1.250000e-04*aNodalField( 26 );

        }
    }
    return Element_RHS;
}

void
Hierarchical_Mesh_Main::MassMatrix_for_L2_projection(
        const uint & aDim,
        const uint & aPolynomial,
        Mat<real> &  aElementMass
)
{
    // set size of output matrix
    aElementMass.set_size(
            pow(aPolynomial + 1,aDim),
            pow(aPolynomial + 1,aDim));

    // FIXME tidy this mess up, add comments and write matrices in row major order!
    if ( aDim == 2)
    {
        if ( aPolynomial == 1)
        {
            aElementMass(0,0) = 1.0/9.0; aElementMass(0,1) = 1.0/18.0; aElementMass(0,2) = 1.0/18.0; aElementMass(0,3) = 1.0/36.0;
            aElementMass(1,0) = 1.0/18.0; aElementMass(1,1) = 1.0/9.0; aElementMass(1,2) = 1.0/36.0; aElementMass(1,3) = 1.0/18.0;
            aElementMass(2,0) = 1.0/18.0; aElementMass(2,1) = 1.0/36.0; aElementMass(2,2) = 1.0/9.0; aElementMass(2,3) = 1.0/18.0;
            aElementMass(3,0) = 1.0/36.0; aElementMass(3,1) = 1.0/18.0; aElementMass(3,2) = 1.0/18.0; aElementMass(3,3) = 1.0/9.0;
        }
        else if ( aPolynomial == 2)
        {
            aElementMass(0,0) = 0.002500000000000;  aElementMass(0,1) = 0.005416666666667;  aElementMass(0,2) = 0.000416666666667;  aElementMass(0,3) = 0.005416666666667;  aElementMass(0,4) = 0.011736111111111;  aElementMass(0,5) = 0.000902777777778;  aElementMass(0,6) = 0.000416666666667;  aElementMass(0,7) = 0.000902777777778;  aElementMass(0,8) = 0.000069444444444;
            aElementMass(1,0) = 0.005416666666667;  aElementMass(1,1) = 0.022500000000000;  aElementMass(1,2) = 0.005416666666667;  aElementMass(1,3) = 0.011736111111111;  aElementMass(1,4) = 0.048750000000000;  aElementMass(1,5) = 0.011736111111111;  aElementMass(1,6) = 0.000902777777778;  aElementMass(1,7) = 0.003750000000000;  aElementMass(1,8) = 0.000902777777778;
            aElementMass(2,0) = 0.000416666666667;  aElementMass(2,1) = 0.005416666666667;  aElementMass(2,2) = 0.002500000000000;  aElementMass(2,3) = 0.000902777777778;  aElementMass(2,4) = 0.011736111111111;  aElementMass(2,5) = 0.005416666666667;  aElementMass(2,6) = 0.000069444444444;  aElementMass(2,7) = 0.000902777777778;  aElementMass(2,8) = 0.000416666666667;
            aElementMass(3,0) = 0.005416666666667;  aElementMass(3,1) = 0.011736111111111;  aElementMass(3,2) = 0.000902777777778;  aElementMass(3,3) = 0.022500000000000;  aElementMass(3,4) = 0.048750000000000;  aElementMass(3,5) = 0.003750000000000;  aElementMass(3,6) = 0.005416666666667;  aElementMass(3,7) = 0.011736111111111;  aElementMass(3,8) = 0.000902777777778;
            aElementMass(4,0) = 0.011736111111111;  aElementMass(4,1) = 0.048750000000000;  aElementMass(4,2) = 0.011736111111111;  aElementMass(4,3) = 0.048750000000000;  aElementMass(4,4) = 0.202500000000000;  aElementMass(4,5) = 0.048750000000000;  aElementMass(4,6) = 0.011736111111111;  aElementMass(4,7) = 0.048750000000000;  aElementMass(4,8) = 0.011736111111111;
            aElementMass(5,0) = 0.000902777777778;  aElementMass(5,1) = 0.011736111111111;  aElementMass(5,2) = 0.005416666666667;  aElementMass(5,3) = 0.003750000000000;  aElementMass(5,4) = 0.048750000000000;  aElementMass(5,5) = 0.022500000000000;  aElementMass(5,6) = 0.000902777777778;  aElementMass(5,7) = 0.011736111111111;  aElementMass(5,8) = 0.005416666666667;
            aElementMass(6,0) = 0.000416666666667;  aElementMass(6,1) = 0.000902777777778;  aElementMass(6,2) = 0.000069444444444;  aElementMass(6,3) = 0.005416666666667;  aElementMass(6,4) = 0.011736111111111;  aElementMass(6,5) = 0.000902777777778;  aElementMass(6,6) = 0.002500000000000;  aElementMass(6,7) = 0.005416666666667;  aElementMass(6,8) = 0.000416666666667;
            aElementMass(7,0) = 0.000902777777778;  aElementMass(7,1) = 0.003750000000000;  aElementMass(7,2) = 0.000902777777778;  aElementMass(7,3) = 0.011736111111111;  aElementMass(7,4) = 0.048750000000000;  aElementMass(7,5) = 0.011736111111111;  aElementMass(7,6) = 0.005416666666667;  aElementMass(7,7) = 0.022500000000000;  aElementMass(7,8) = 0.005416666666667;
            aElementMass(8,0) = 0.000069444444444;  aElementMass(8,1) = 0.000902777777778;  aElementMass(8,2) = 0.000416666666667;  aElementMass(8,3) = 0.000902777777778;  aElementMass(8,4) = 0.011736111111111;  aElementMass(8,5) = 0.005416666666667;  aElementMass(8,6) = 0.000416666666667;  aElementMass(8,7) = 0.005416666666667;  aElementMass(8,8) = 0.002500000000000;
        }
    }
    else if ( aDim == 3)
    {
        if ( aPolynomial == 1)
        {
            aElementMass(0,0) = 1.0/216.0*8.0; aElementMass(0,1) = 1.0/216.0*4.0; aElementMass(0,2) = 1.0/216.0*4.0; aElementMass(0,3) = 1.0/216.0*2.0; aElementMass(0,4) = 1.0/216.0*4.0; aElementMass(0,5) = 1.0/216.0*2.0; aElementMass(0,6) = 1.0/216.0*2.0; aElementMass(0,7) = 1.0/216.0*1.0;
            aElementMass(1,0) = 1.0/216.0*4.0; aElementMass(1,1) = 1.0/216.0*8.0; aElementMass(1,2) = 1.0/216.0*2.0; aElementMass(1,3) = 1.0/216.0*4.0; aElementMass(1,4) = 1.0/216.0*2.0; aElementMass(1,5) = 1.0/216.0*4.0; aElementMass(1,6) = 1.0/216.0*1.0; aElementMass(1,7) = 1.0/216.0*2.0;
            aElementMass(2,0) = 1.0/216.0*4.0; aElementMass(2,1) = 1.0/216.0*2.0; aElementMass(2,2) = 1.0/216.0*8.0; aElementMass(2,3) = 1.0/216.0*4.0; aElementMass(2,4) = 1.0/216.0*2.0; aElementMass(2,5) = 1.0/216.0*1.0; aElementMass(2,6) = 1.0/216.0*4.0; aElementMass(2,7) = 1.0/216.0*2.0;
            aElementMass(3,0) = 1.0/216.0*2.0; aElementMass(3,1) = 1.0/216.0*4.0; aElementMass(3,2) = 1.0/216.0*4.0; aElementMass(3,3) = 1.0/216.0*8.0; aElementMass(3,4) = 1.0/216.0*1.0; aElementMass(3,5) = 1.0/216.0*2.0; aElementMass(3,6) = 1.0/216.0*2.0; aElementMass(3,7) = 1.0/216.0*4.0;
            aElementMass(4,0) = 1.0/216.0*4.0; aElementMass(4,1) = 1.0/216.0*2.0; aElementMass(4,2) = 1.0/216.0*2.0; aElementMass(4,3) = 1.0/216.0*1.0; aElementMass(4,4) = 1.0/216.0*8.0; aElementMass(4,5) = 1.0/216.0*4.0; aElementMass(4,6) = 1.0/216.0*4.0; aElementMass(4,7) = 1.0/216.0*2.0;
            aElementMass(5,0) = 1.0/216.0*2.0; aElementMass(5,1) = 1.0/216.0*4.0; aElementMass(5,2) = 1.0/216.0*1.0; aElementMass(5,3) = 1.0/216.0*2.0; aElementMass(5,4) = 1.0/216.0*4.0; aElementMass(5,5) = 1.0/216.0*8.0; aElementMass(5,6) = 1.0/216.0*2.0; aElementMass(5,7) = 1.0/216.0*4.0;
            aElementMass(6,0) = 1.0/216.0*2.0; aElementMass(6,1) = 1.0/216.0*1.0; aElementMass(6,2) = 1.0/216.0*4.0; aElementMass(6,3) = 1.0/216.0*2.0; aElementMass(6,4) = 1.0/216.0*4.0; aElementMass(6,5) = 1.0/216.0*2.0; aElementMass(6,6) = 1.0/216.0*8.0; aElementMass(6,7) = 1.0/216.0*4.0;
            aElementMass(7,0) = 1.0/216.0*1.0; aElementMass(7,1) = 1.0/216.0*2.0; aElementMass(7,2) = 1.0/216.0*2.0; aElementMass(7,3) = 1.0/216.0*4.0; aElementMass(7,4) = 1.0/216.0*2.0; aElementMass(7,5) = 1.0/216.0*4.0; aElementMass(7,6) = 1.0/216.0*4.0; aElementMass(7,7) = 1.0/216.0*8.0;
        }
        else if ( aPolynomial == 2)
        {
            aElementMass(0,0) = 1.250000e-04; aElementMass(0,1) = 2.708333e-04; aElementMass(0,2) = 2.083333e-05; aElementMass(0,3) = 2.708333e-04; aElementMass(0,4) = 5.868056e-04; aElementMass(0,5) = 4.513889e-05; aElementMass(0,6) = 2.083333e-05; aElementMass(0,7) = 4.513889e-05; aElementMass(0,8) = 3.472222e-06; aElementMass(0,9) = 2.708333e-04; aElementMass(0,10) = 5.868056e-04; aElementMass(0,11) = 4.513889e-05; aElementMass(0,12) = 5.868056e-04; aElementMass(0,13) = 1.271412e-03; aElementMass(0,14) = 9.780093e-05; aElementMass(0,15) = 4.513889e-05; aElementMass(0,16) = 9.780093e-05; aElementMass(0,17) = 7.523148e-06; aElementMass(0,18) = 2.083333e-05; aElementMass(0,19) = 4.513889e-05; aElementMass(0,20) = 3.472222e-06; aElementMass(0,21) = 4.513889e-05; aElementMass(0,22) = 9.780093e-05; aElementMass(0,23) = 7.523148e-06; aElementMass(0,24) = 3.472222e-06; aElementMass(0,25) = 7.523148e-06; aElementMass(0,26) = 5.787037e-07;
            aElementMass(1,0) = 2.708333e-04; aElementMass(1,1) = 1.125000e-03; aElementMass(1,2) = 2.708333e-04; aElementMass(1,3) = 5.868056e-04; aElementMass(1,4) = 2.437500e-03; aElementMass(1,5) = 5.868056e-04; aElementMass(1,6) = 4.513889e-05; aElementMass(1,7) = 1.875000e-04; aElementMass(1,8) = 4.513889e-05; aElementMass(1,9) = 5.868056e-04; aElementMass(1,10) = 2.437500e-03; aElementMass(1,11) = 5.868056e-04; aElementMass(1,12) = 1.271412e-03; aElementMass(1,13) = 5.281250e-03; aElementMass(1,14) = 1.271412e-03; aElementMass(1,15) = 9.780093e-05; aElementMass(1,16) = 4.062500e-04; aElementMass(1,17) = 9.780093e-05; aElementMass(1,18) = 4.513889e-05; aElementMass(1,19) = 1.875000e-04; aElementMass(1,20) = 4.513889e-05; aElementMass(1,21) = 9.780093e-05; aElementMass(1,22) = 4.062500e-04; aElementMass(1,23) = 9.780093e-05; aElementMass(1,24) = 7.523148e-06; aElementMass(1,25) = 3.125000e-05; aElementMass(1,26) = 7.523148e-06;
            aElementMass(2,0) = 2.083333e-05; aElementMass(2,1) = 2.708333e-04; aElementMass(2,2) = 1.250000e-04; aElementMass(2,3) = 4.513889e-05; aElementMass(2,4) = 5.868056e-04; aElementMass(2,5) = 2.708333e-04; aElementMass(2,6) = 3.472222e-06; aElementMass(2,7) = 4.513889e-05; aElementMass(2,8) = 2.083333e-05; aElementMass(2,9) = 4.513889e-05; aElementMass(2,10) = 5.868056e-04; aElementMass(2,11) = 2.708333e-04; aElementMass(2,12) = 9.780093e-05; aElementMass(2,13) = 1.271412e-03; aElementMass(2,14) = 5.868056e-04; aElementMass(2,15) = 7.523148e-06; aElementMass(2,16) = 9.780093e-05; aElementMass(2,17) = 4.513889e-05; aElementMass(2,18) = 3.472222e-06; aElementMass(2,19) = 4.513889e-05; aElementMass(2,20) = 2.083333e-05; aElementMass(2,21) = 7.523148e-06; aElementMass(2,22) = 9.780093e-05; aElementMass(2,23) = 4.513889e-05; aElementMass(2,24) = 5.787037e-07; aElementMass(2,25) = 7.523148e-06; aElementMass(2,26) = 3.472222e-06;
            aElementMass(3,0) = 2.708333e-04; aElementMass(3,1) = 5.868056e-04; aElementMass(3,2) = 4.513889e-05; aElementMass(3,3) = 1.125000e-03; aElementMass(3,4) = 2.437500e-03; aElementMass(3,5) = 1.875000e-04; aElementMass(3,6) = 2.708333e-04; aElementMass(3,7) = 5.868056e-04; aElementMass(3,8) = 4.513889e-05; aElementMass(3,9) = 5.868056e-04; aElementMass(3,10) = 1.271412e-03; aElementMass(3,11) = 9.780093e-05; aElementMass(3,12) = 2.437500e-03; aElementMass(3,13) = 5.281250e-03; aElementMass(3,14) = 4.062500e-04; aElementMass(3,15) = 5.868056e-04; aElementMass(3,16) = 1.271412e-03; aElementMass(3,17) = 9.780093e-05; aElementMass(3,18) = 4.513889e-05; aElementMass(3,19) = 9.780093e-05; aElementMass(3,20) = 7.523148e-06; aElementMass(3,21) = 1.875000e-04; aElementMass(3,22) = 4.062500e-04; aElementMass(3,23) = 3.125000e-05; aElementMass(3,24) = 4.513889e-05; aElementMass(3,25) = 9.780093e-05; aElementMass(3,26) = 7.523148e-06;
            aElementMass(4,0) = 5.868056e-04; aElementMass(4,1) = 2.437500e-03; aElementMass(4,2) = 5.868056e-04; aElementMass(4,3) = 2.437500e-03; aElementMass(4,4) = 1.012500e-02; aElementMass(4,5) = 2.437500e-03; aElementMass(4,6) = 5.868056e-04; aElementMass(4,7) = 2.437500e-03; aElementMass(4,8) = 5.868056e-04; aElementMass(4,9) = 1.271412e-03; aElementMass(4,10) = 5.281250e-03; aElementMass(4,11) = 1.271412e-03; aElementMass(4,12) = 5.281250e-03; aElementMass(4,13) = 2.193750e-02; aElementMass(4,14) = 5.281250e-03; aElementMass(4,15) = 1.271412e-03; aElementMass(4,16) = 5.281250e-03; aElementMass(4,17) = 1.271412e-03; aElementMass(4,18) = 9.780093e-05; aElementMass(4,19) = 4.062500e-04; aElementMass(4,20) = 9.780093e-05; aElementMass(4,21) = 4.062500e-04; aElementMass(4,22) = 1.687500e-03; aElementMass(4,23) = 4.062500e-04; aElementMass(4,24) = 9.780093e-05; aElementMass(4,25) = 4.062500e-04; aElementMass(4,26) = 9.780093e-05;
            aElementMass(5,0) = 4.513889e-05; aElementMass(5,1) = 5.868056e-04; aElementMass(5,2) = 2.708333e-04; aElementMass(5,3) = 1.875000e-04; aElementMass(5,4) = 2.437500e-03; aElementMass(5,5) = 1.125000e-03; aElementMass(5,6) = 4.513889e-05; aElementMass(5,7) = 5.868056e-04; aElementMass(5,8) = 2.708333e-04; aElementMass(5,9) = 9.780093e-05; aElementMass(5,10) = 1.271412e-03; aElementMass(5,11) = 5.868056e-04; aElementMass(5,12) = 4.062500e-04; aElementMass(5,13) = 5.281250e-03; aElementMass(5,14) = 2.437500e-03; aElementMass(5,15) = 9.780093e-05; aElementMass(5,16) = 1.271412e-03; aElementMass(5,17) = 5.868056e-04; aElementMass(5,18) = 7.523148e-06; aElementMass(5,19) = 9.780093e-05; aElementMass(5,20) = 4.513889e-05; aElementMass(5,21) = 3.125000e-05; aElementMass(5,22) = 4.062500e-04; aElementMass(5,23) = 1.875000e-04; aElementMass(5,24) = 7.523148e-06; aElementMass(5,25) = 9.780093e-05; aElementMass(5,26) = 4.513889e-05;
            aElementMass(6,0) = 2.083333e-05; aElementMass(6,1) = 4.513889e-05; aElementMass(6,2) = 3.472222e-06; aElementMass(6,3) = 2.708333e-04; aElementMass(6,4) = 5.868056e-04; aElementMass(6,5) = 4.513889e-05; aElementMass(6,6) = 1.250000e-04; aElementMass(6,7) = 2.708333e-04; aElementMass(6,8) = 2.083333e-05; aElementMass(6,9) = 4.513889e-05; aElementMass(6,10) = 9.780093e-05; aElementMass(6,11) = 7.523148e-06; aElementMass(6,12) = 5.868056e-04; aElementMass(6,13) = 1.271412e-03; aElementMass(6,14) = 9.780093e-05; aElementMass(6,15) = 2.708333e-04; aElementMass(6,16) = 5.868056e-04; aElementMass(6,17) = 4.513889e-05; aElementMass(6,18) = 3.472222e-06; aElementMass(6,19) = 7.523148e-06; aElementMass(6,20) = 5.787037e-07; aElementMass(6,21) = 4.513889e-05; aElementMass(6,22) = 9.780093e-05; aElementMass(6,23) = 7.523148e-06; aElementMass(6,24) = 2.083333e-05; aElementMass(6,25) = 4.513889e-05; aElementMass(6,26) = 3.472222e-06;
            aElementMass(7,0) = 4.513889e-05; aElementMass(7,1) = 1.875000e-04; aElementMass(7,2) = 4.513889e-05; aElementMass(7,3) = 5.868056e-04; aElementMass(7,4) = 2.437500e-03; aElementMass(7,5) = 5.868056e-04; aElementMass(7,6) = 2.708333e-04; aElementMass(7,7) = 1.125000e-03; aElementMass(7,8) = 2.708333e-04; aElementMass(7,9) = 9.780093e-05; aElementMass(7,10) = 4.062500e-04; aElementMass(7,11) = 9.780093e-05; aElementMass(7,12) = 1.271412e-03; aElementMass(7,13) = 5.281250e-03; aElementMass(7,14) = 1.271412e-03; aElementMass(7,15) = 5.868056e-04; aElementMass(7,16) = 2.437500e-03; aElementMass(7,17) = 5.868056e-04; aElementMass(7,18) = 7.523148e-06; aElementMass(7,19) = 3.125000e-05; aElementMass(7,20) = 7.523148e-06; aElementMass(7,21) = 9.780093e-05; aElementMass(7,22) = 4.062500e-04; aElementMass(7,23) = 9.780093e-05; aElementMass(7,24) = 4.513889e-05; aElementMass(7,25) = 1.875000e-04; aElementMass(7,26) = 4.513889e-05;
            aElementMass(8,0) = 3.472222e-06; aElementMass(8,1) = 4.513889e-05; aElementMass(8,2) = 2.083333e-05; aElementMass(8,3) = 4.513889e-05; aElementMass(8,4) = 5.868056e-04; aElementMass(8,5) = 2.708333e-04; aElementMass(8,6) = 2.083333e-05; aElementMass(8,7) = 2.708333e-04; aElementMass(8,8) = 1.250000e-04; aElementMass(8,9) = 7.523148e-06; aElementMass(8,10) = 9.780093e-05; aElementMass(8,11) = 4.513889e-05; aElementMass(8,12) = 9.780093e-05; aElementMass(8,13) = 1.271412e-03; aElementMass(8,14) = 5.868056e-04; aElementMass(8,15) = 4.513889e-05; aElementMass(8,16) = 5.868056e-04; aElementMass(8,17) = 2.708333e-04; aElementMass(8,18) = 5.787037e-07; aElementMass(8,19) = 7.523148e-06; aElementMass(8,20) = 3.472222e-06; aElementMass(8,21) = 7.523148e-06; aElementMass(8,22) = 9.780093e-05; aElementMass(8,23) = 4.513889e-05; aElementMass(8,24) = 3.472222e-06; aElementMass(8,25) = 4.513889e-05; aElementMass(8,26) = 2.083333e-05;
            aElementMass(9,0) = 2.708333e-04; aElementMass(9,1) = 5.868056e-04; aElementMass(9,2) = 4.513889e-05; aElementMass(9,3) = 5.868056e-04; aElementMass(9,4) = 1.271412e-03; aElementMass(9,5) = 9.780093e-05; aElementMass(9,6) = 4.513889e-05; aElementMass(9,7) = 9.780093e-05; aElementMass(9,8) = 7.523148e-06; aElementMass(9,9) = 1.125000e-03; aElementMass(9,10) = 2.437500e-03; aElementMass(9,11) = 1.875000e-04; aElementMass(9,12) = 2.437500e-03; aElementMass(9,13) = 5.281250e-03; aElementMass(9,14) = 4.062500e-04; aElementMass(9,15) = 1.875000e-04; aElementMass(9,16) = 4.062500e-04; aElementMass(9,17) = 3.125000e-05; aElementMass(9,18) = 2.708333e-04; aElementMass(9,19) = 5.868056e-04; aElementMass(9,20) = 4.513889e-05; aElementMass(9,21) = 5.868056e-04; aElementMass(9,22) = 1.271412e-03; aElementMass(9,23) = 9.780093e-05; aElementMass(9,24) = 4.513889e-05; aElementMass(9,25) = 9.780093e-05; aElementMass(9,26) = 7.523148e-06;
            aElementMass(10,0) = 5.868056e-04; aElementMass(10,1) = 2.437500e-03; aElementMass(10,2) = 5.868056e-04; aElementMass(10,3) = 1.271412e-03; aElementMass(10,4) = 5.281250e-03; aElementMass(10,5) = 1.271412e-03; aElementMass(10,6) = 9.780093e-05; aElementMass(10,7) = 4.062500e-04; aElementMass(10,8) = 9.780093e-05; aElementMass(10,9) = 2.437500e-03; aElementMass(10,10) = 1.012500e-02; aElementMass(10,11) = 2.437500e-03; aElementMass(10,12) = 5.281250e-03; aElementMass(10,13) = 2.193750e-02; aElementMass(10,14) = 5.281250e-03; aElementMass(10,15) = 4.062500e-04; aElementMass(10,16) = 1.687500e-03; aElementMass(10,17) = 4.062500e-04; aElementMass(10,18) = 5.868056e-04; aElementMass(10,19) = 2.437500e-03; aElementMass(10,20) = 5.868056e-04; aElementMass(10,21) = 1.271412e-03; aElementMass(10,22) = 5.281250e-03; aElementMass(10,23) = 1.271412e-03; aElementMass(10,24) = 9.780093e-05; aElementMass(10,25) = 4.062500e-04; aElementMass(10,26) = 9.780093e-05;
            aElementMass(11,0) = 4.513889e-05; aElementMass(11,1) = 5.868056e-04; aElementMass(11,2) = 2.708333e-04; aElementMass(11,3) = 9.780093e-05; aElementMass(11,4) = 1.271412e-03; aElementMass(11,5) = 5.868056e-04; aElementMass(11,6) = 7.523148e-06; aElementMass(11,7) = 9.780093e-05; aElementMass(11,8) = 4.513889e-05; aElementMass(11,9) = 1.875000e-04; aElementMass(11,10) = 2.437500e-03; aElementMass(11,11) = 1.125000e-03; aElementMass(11,12) = 4.062500e-04; aElementMass(11,13) = 5.281250e-03; aElementMass(11,14) = 2.437500e-03; aElementMass(11,15) = 3.125000e-05; aElementMass(11,16) = 4.062500e-04; aElementMass(11,17) = 1.875000e-04; aElementMass(11,18) = 4.513889e-05; aElementMass(11,19) = 5.868056e-04; aElementMass(11,20) = 2.708333e-04; aElementMass(11,21) = 9.780093e-05; aElementMass(11,22) = 1.271412e-03; aElementMass(11,23) = 5.868056e-04; aElementMass(11,24) = 7.523148e-06; aElementMass(11,25) = 9.780093e-05; aElementMass(11,26) = 4.513889e-05;
            aElementMass(12,0) = 5.868056e-04; aElementMass(12,1) = 1.271412e-03; aElementMass(12,2) = 9.780093e-05; aElementMass(12,3) = 2.437500e-03; aElementMass(12,4) = 5.281250e-03; aElementMass(12,5) = 4.062500e-04; aElementMass(12,6) = 5.868056e-04; aElementMass(12,7) = 1.271412e-03; aElementMass(12,8) = 9.780093e-05; aElementMass(12,9) = 2.437500e-03; aElementMass(12,10) = 5.281250e-03; aElementMass(12,11) = 4.062500e-04; aElementMass(12,12) = 1.012500e-02; aElementMass(12,13) = 2.193750e-02; aElementMass(12,14) = 1.687500e-03; aElementMass(12,15) = 2.437500e-03; aElementMass(12,16) = 5.281250e-03; aElementMass(12,17) = 4.062500e-04; aElementMass(12,18) = 5.868056e-04; aElementMass(12,19) = 1.271412e-03; aElementMass(12,20) = 9.780093e-05; aElementMass(12,21) = 2.437500e-03; aElementMass(12,22) = 5.281250e-03; aElementMass(12,23) = 4.062500e-04; aElementMass(12,24) = 5.868056e-04; aElementMass(12,25) = 1.271412e-03; aElementMass(12,26) = 9.780093e-05;
            aElementMass(13,0) = 1.271412e-03; aElementMass(13,1) = 5.281250e-03; aElementMass(13,2) = 1.271412e-03; aElementMass(13,3) = 5.281250e-03; aElementMass(13,4) = 2.193750e-02; aElementMass(13,5) = 5.281250e-03; aElementMass(13,6) = 1.271412e-03; aElementMass(13,7) = 5.281250e-03; aElementMass(13,8) = 1.271412e-03; aElementMass(13,9) = 5.281250e-03; aElementMass(13,10) = 2.193750e-02; aElementMass(13,11) = 5.281250e-03; aElementMass(13,12) = 2.193750e-02; aElementMass(13,13) = 9.112500e-02; aElementMass(13,14) = 2.193750e-02; aElementMass(13,15) = 5.281250e-03; aElementMass(13,16) = 2.193750e-02; aElementMass(13,17) = 5.281250e-03; aElementMass(13,18) = 1.271412e-03; aElementMass(13,19) = 5.281250e-03; aElementMass(13,20) = 1.271412e-03; aElementMass(13,21) = 5.281250e-03; aElementMass(13,22) = 2.193750e-02; aElementMass(13,23) = 5.281250e-03; aElementMass(13,24) = 1.271412e-03; aElementMass(13,25) = 5.281250e-03; aElementMass(13,26) = 1.271412e-03;
            aElementMass(14,0) = 9.780093e-05; aElementMass(14,1) = 1.271412e-03; aElementMass(14,2) = 5.868056e-04; aElementMass(14,3) = 4.062500e-04; aElementMass(14,4) = 5.281250e-03; aElementMass(14,5) = 2.437500e-03; aElementMass(14,6) = 9.780093e-05; aElementMass(14,7) = 1.271412e-03; aElementMass(14,8) = 5.868056e-04; aElementMass(14,9) = 4.062500e-04; aElementMass(14,10) = 5.281250e-03; aElementMass(14,11) = 2.437500e-03; aElementMass(14,12) = 1.687500e-03; aElementMass(14,13) = 2.193750e-02; aElementMass(14,14) = 1.012500e-02; aElementMass(14,15) = 4.062500e-04; aElementMass(14,16) = 5.281250e-03; aElementMass(14,17) = 2.437500e-03; aElementMass(14,18) = 9.780093e-05; aElementMass(14,19) = 1.271412e-03; aElementMass(14,20) = 5.868056e-04; aElementMass(14,21) = 4.062500e-04; aElementMass(14,22) = 5.281250e-03; aElementMass(14,23) = 2.437500e-03; aElementMass(14,24) = 9.780093e-05; aElementMass(14,25) = 1.271412e-03; aElementMass(14,26) = 5.868056e-04;
            aElementMass(15,0) = 4.513889e-05; aElementMass(15,1) = 9.780093e-05; aElementMass(15,2) = 7.523148e-06; aElementMass(15,3) = 5.868056e-04; aElementMass(15,4) = 1.271412e-03; aElementMass(15,5) = 9.780093e-05; aElementMass(15,6) = 2.708333e-04; aElementMass(15,7) = 5.868056e-04; aElementMass(15,8) = 4.513889e-05; aElementMass(15,9) = 1.875000e-04; aElementMass(15,10) = 4.062500e-04; aElementMass(15,11) = 3.125000e-05; aElementMass(15,12) = 2.437500e-03; aElementMass(15,13) = 5.281250e-03; aElementMass(15,14) = 4.062500e-04; aElementMass(15,15) = 1.125000e-03; aElementMass(15,16) = 2.437500e-03; aElementMass(15,17) = 1.875000e-04; aElementMass(15,18) = 4.513889e-05; aElementMass(15,19) = 9.780093e-05; aElementMass(15,20) = 7.523148e-06; aElementMass(15,21) = 5.868056e-04; aElementMass(15,22) = 1.271412e-03; aElementMass(15,23) = 9.780093e-05; aElementMass(15,24) = 2.708333e-04; aElementMass(15,25) = 5.868056e-04; aElementMass(15,26) = 4.513889e-05;
            aElementMass(16,0) = 9.780093e-05; aElementMass(16,1) = 4.062500e-04; aElementMass(16,2) = 9.780093e-05; aElementMass(16,3) = 1.271412e-03; aElementMass(16,4) = 5.281250e-03; aElementMass(16,5) = 1.271412e-03; aElementMass(16,6) = 5.868056e-04; aElementMass(16,7) = 2.437500e-03; aElementMass(16,8) = 5.868056e-04; aElementMass(16,9) = 4.062500e-04; aElementMass(16,10) = 1.687500e-03; aElementMass(16,11) = 4.062500e-04; aElementMass(16,12) = 5.281250e-03; aElementMass(16,13) = 2.193750e-02; aElementMass(16,14) = 5.281250e-03; aElementMass(16,15) = 2.437500e-03; aElementMass(16,16) = 1.012500e-02; aElementMass(16,17) = 2.437500e-03; aElementMass(16,18) = 9.780093e-05; aElementMass(16,19) = 4.062500e-04; aElementMass(16,20) = 9.780093e-05; aElementMass(16,21) = 1.271412e-03; aElementMass(16,22) = 5.281250e-03; aElementMass(16,23) = 1.271412e-03; aElementMass(16,24) = 5.868056e-04; aElementMass(16,25) = 2.437500e-03; aElementMass(16,26) = 5.868056e-04;
            aElementMass(17,0) = 7.523148e-06; aElementMass(17,1) = 9.780093e-05; aElementMass(17,2) = 4.513889e-05; aElementMass(17,3) = 9.780093e-05; aElementMass(17,4) = 1.271412e-03; aElementMass(17,5) = 5.868056e-04; aElementMass(17,6) = 4.513889e-05; aElementMass(17,7) = 5.868056e-04; aElementMass(17,8) = 2.708333e-04; aElementMass(17,9) = 3.125000e-05; aElementMass(17,10) = 4.062500e-04; aElementMass(17,11) = 1.875000e-04; aElementMass(17,12) = 4.062500e-04; aElementMass(17,13) = 5.281250e-03; aElementMass(17,14) = 2.437500e-03; aElementMass(17,15) = 1.875000e-04; aElementMass(17,16) = 2.437500e-03; aElementMass(17,17) = 1.125000e-03; aElementMass(17,18) = 7.523148e-06; aElementMass(17,19) = 9.780093e-05; aElementMass(17,20) = 4.513889e-05; aElementMass(17,21) = 9.780093e-05; aElementMass(17,22) = 1.271412e-03; aElementMass(17,23) = 5.868056e-04; aElementMass(17,24) = 4.513889e-05; aElementMass(17,25) = 5.868056e-04; aElementMass(17,26) = 2.708333e-04;
            aElementMass(18,0) = 2.083333e-05; aElementMass(18,1) = 4.513889e-05; aElementMass(18,2) = 3.472222e-06; aElementMass(18,3) = 4.513889e-05; aElementMass(18,4) = 9.780093e-05; aElementMass(18,5) = 7.523148e-06; aElementMass(18,6) = 3.472222e-06; aElementMass(18,7) = 7.523148e-06; aElementMass(18,8) = 5.787037e-07; aElementMass(18,9) = 2.708333e-04; aElementMass(18,10) = 5.868056e-04; aElementMass(18,11) = 4.513889e-05; aElementMass(18,12) = 5.868056e-04; aElementMass(18,13) = 1.271412e-03; aElementMass(18,14) = 9.780093e-05; aElementMass(18,15) = 4.513889e-05; aElementMass(18,16) = 9.780093e-05; aElementMass(18,17) = 7.523148e-06; aElementMass(18,18) = 1.250000e-04; aElementMass(18,19) = 2.708333e-04; aElementMass(18,20) = 2.083333e-05; aElementMass(18,21) = 2.708333e-04; aElementMass(18,22) = 5.868056e-04; aElementMass(18,23) = 4.513889e-05; aElementMass(18,24) = 2.083333e-05; aElementMass(18,25) = 4.513889e-05; aElementMass(18,26) = 3.472222e-06;
            aElementMass(19,0) = 4.513889e-05; aElementMass(19,1) = 1.875000e-04; aElementMass(19,2) = 4.513889e-05; aElementMass(19,3) = 9.780093e-05; aElementMass(19,4) = 4.062500e-04; aElementMass(19,5) = 9.780093e-05; aElementMass(19,6) = 7.523148e-06; aElementMass(19,7) = 3.125000e-05; aElementMass(19,8) = 7.523148e-06; aElementMass(19,9) = 5.868056e-04; aElementMass(19,10) = 2.437500e-03; aElementMass(19,11) = 5.868056e-04; aElementMass(19,12) = 1.271412e-03; aElementMass(19,13) = 5.281250e-03; aElementMass(19,14) = 1.271412e-03; aElementMass(19,15) = 9.780093e-05; aElementMass(19,16) = 4.062500e-04; aElementMass(19,17) = 9.780093e-05; aElementMass(19,18) = 2.708333e-04; aElementMass(19,19) = 1.125000e-03; aElementMass(19,20) = 2.708333e-04; aElementMass(19,21) = 5.868056e-04; aElementMass(19,22) = 2.437500e-03; aElementMass(19,23) = 5.868056e-04; aElementMass(19,24) = 4.513889e-05; aElementMass(19,25) = 1.875000e-04; aElementMass(19,26) = 4.513889e-05;
            aElementMass(20,0) = 3.472222e-06; aElementMass(20,1) = 4.513889e-05; aElementMass(20,2) = 2.083333e-05; aElementMass(20,3) = 7.523148e-06; aElementMass(20,4) = 9.780093e-05; aElementMass(20,5) = 4.513889e-05; aElementMass(20,6) = 5.787037e-07; aElementMass(20,7) = 7.523148e-06; aElementMass(20,8) = 3.472222e-06; aElementMass(20,9) = 4.513889e-05; aElementMass(20,10) = 5.868056e-04; aElementMass(20,11) = 2.708333e-04; aElementMass(20,12) = 9.780093e-05; aElementMass(20,13) = 1.271412e-03; aElementMass(20,14) = 5.868056e-04; aElementMass(20,15) = 7.523148e-06; aElementMass(20,16) = 9.780093e-05; aElementMass(20,17) = 4.513889e-05; aElementMass(20,18) = 2.083333e-05; aElementMass(20,19) = 2.708333e-04; aElementMass(20,20) = 1.250000e-04; aElementMass(20,21) = 4.513889e-05; aElementMass(20,22) = 5.868056e-04; aElementMass(20,23) = 2.708333e-04; aElementMass(20,24) = 3.472222e-06; aElementMass(20,25) = 4.513889e-05; aElementMass(20,26) = 2.083333e-05;
            aElementMass(21,0) = 4.513889e-05; aElementMass(21,1) = 9.780093e-05; aElementMass(21,2) = 7.523148e-06; aElementMass(21,3) = 1.875000e-04; aElementMass(21,4) = 4.062500e-04; aElementMass(21,5) = 3.125000e-05; aElementMass(21,6) = 4.513889e-05; aElementMass(21,7) = 9.780093e-05; aElementMass(21,8) = 7.523148e-06; aElementMass(21,9) = 5.868056e-04; aElementMass(21,10) = 1.271412e-03; aElementMass(21,11) = 9.780093e-05; aElementMass(21,12) = 2.437500e-03; aElementMass(21,13) = 5.281250e-03; aElementMass(21,14) = 4.062500e-04; aElementMass(21,15) = 5.868056e-04; aElementMass(21,16) = 1.271412e-03; aElementMass(21,17) = 9.780093e-05; aElementMass(21,18) = 2.708333e-04; aElementMass(21,19) = 5.868056e-04; aElementMass(21,20) = 4.513889e-05; aElementMass(21,21) = 1.125000e-03; aElementMass(21,22) = 2.437500e-03; aElementMass(21,23) = 1.875000e-04; aElementMass(21,24) = 2.708333e-04; aElementMass(21,25) = 5.868056e-04; aElementMass(21,26) = 4.513889e-05;
            aElementMass(22,0) = 9.780093e-05; aElementMass(22,1) = 4.062500e-04; aElementMass(22,2) = 9.780093e-05; aElementMass(22,3) = 4.062500e-04; aElementMass(22,4) = 1.687500e-03; aElementMass(22,5) = 4.062500e-04; aElementMass(22,6) = 9.780093e-05; aElementMass(22,7) = 4.062500e-04; aElementMass(22,8) = 9.780093e-05; aElementMass(22,9) = 1.271412e-03; aElementMass(22,10) = 5.281250e-03; aElementMass(22,11) = 1.271412e-03; aElementMass(22,12) = 5.281250e-03; aElementMass(22,13) = 2.193750e-02; aElementMass(22,14) = 5.281250e-03; aElementMass(22,15) = 1.271412e-03; aElementMass(22,16) = 5.281250e-03; aElementMass(22,17) = 1.271412e-03; aElementMass(22,18) = 5.868056e-04; aElementMass(22,19) = 2.437500e-03; aElementMass(22,20) = 5.868056e-04; aElementMass(22,21) = 2.437500e-03; aElementMass(22,22) = 1.012500e-02; aElementMass(22,23) = 2.437500e-03; aElementMass(22,24) = 5.868056e-04; aElementMass(22,25) = 2.437500e-03; aElementMass(22,26) = 5.868056e-04;
            aElementMass(23,0) = 7.523148e-06; aElementMass(23,1) = 9.780093e-05; aElementMass(23,2) = 4.513889e-05; aElementMass(23,3) = 3.125000e-05; aElementMass(23,4) = 4.062500e-04; aElementMass(23,5) = 1.875000e-04; aElementMass(23,6) = 7.523148e-06; aElementMass(23,7) = 9.780093e-05; aElementMass(23,8) = 4.513889e-05; aElementMass(23,9) = 9.780093e-05; aElementMass(23,10) = 1.271412e-03; aElementMass(23,11) = 5.868056e-04; aElementMass(23,12) = 4.062500e-04; aElementMass(23,13) = 5.281250e-03; aElementMass(23,14) = 2.437500e-03; aElementMass(23,15) = 9.780093e-05; aElementMass(23,16) = 1.271412e-03; aElementMass(23,17) = 5.868056e-04; aElementMass(23,18) = 4.513889e-05; aElementMass(23,19) = 5.868056e-04; aElementMass(23,20) = 2.708333e-04; aElementMass(23,21) = 1.875000e-04; aElementMass(23,22) = 2.437500e-03; aElementMass(23,23) = 1.125000e-03; aElementMass(23,24) = 4.513889e-05; aElementMass(23,25) = 5.868056e-04; aElementMass(23,26) = 2.708333e-04;
            aElementMass(24,0) = 3.472222e-06; aElementMass(24,1) = 7.523148e-06; aElementMass(24,2) = 5.787037e-07; aElementMass(24,3) = 4.513889e-05; aElementMass(24,4) = 9.780093e-05; aElementMass(24,5) = 7.523148e-06; aElementMass(24,6) = 2.083333e-05; aElementMass(24,7) = 4.513889e-05; aElementMass(24,8) = 3.472222e-06; aElementMass(24,9) = 4.513889e-05; aElementMass(24,10) = 9.780093e-05; aElementMass(24,11) = 7.523148e-06; aElementMass(24,12) = 5.868056e-04; aElementMass(24,13) = 1.271412e-03; aElementMass(24,14) = 9.780093e-05; aElementMass(24,15) = 2.708333e-04; aElementMass(24,16) = 5.868056e-04; aElementMass(24,17) = 4.513889e-05; aElementMass(24,18) = 2.083333e-05; aElementMass(24,19) = 4.513889e-05; aElementMass(24,20) = 3.472222e-06; aElementMass(24,21) = 2.708333e-04; aElementMass(24,22) = 5.868056e-04; aElementMass(24,23) = 4.513889e-05; aElementMass(24,24) = 1.250000e-04; aElementMass(24,25) = 2.708333e-04; aElementMass(24,26) = 2.083333e-05;
            aElementMass(25,0) = 7.523148e-06; aElementMass(25,1) = 3.125000e-05; aElementMass(25,2) = 7.523148e-06; aElementMass(25,3) = 9.780093e-05; aElementMass(25,4) = 4.062500e-04; aElementMass(25,5) = 9.780093e-05; aElementMass(25,6) = 4.513889e-05; aElementMass(25,7) = 1.875000e-04; aElementMass(25,8) = 4.513889e-05; aElementMass(25,9) = 9.780093e-05; aElementMass(25,10) = 4.062500e-04; aElementMass(25,11) = 9.780093e-05; aElementMass(25,12) = 1.271412e-03; aElementMass(25,13) = 5.281250e-03; aElementMass(25,14) = 1.271412e-03; aElementMass(25,15) = 5.868056e-04; aElementMass(25,16) = 2.437500e-03; aElementMass(25,17) = 5.868056e-04; aElementMass(25,18) = 4.513889e-05; aElementMass(25,19) = 1.875000e-04; aElementMass(25,20) = 4.513889e-05; aElementMass(25,21) = 5.868056e-04; aElementMass(25,22) = 2.437500e-03; aElementMass(25,23) = 5.868056e-04; aElementMass(25,24) = 2.708333e-04; aElementMass(25,25) = 1.125000e-03; aElementMass(25,26) = 2.708333e-04;
            aElementMass(26,0) = 5.787037e-07; aElementMass(26,1) = 7.523148e-06; aElementMass(26,2) = 3.472222e-06; aElementMass(26,3) = 7.523148e-06; aElementMass(26,4) = 9.780093e-05; aElementMass(26,5) = 4.513889e-05; aElementMass(26,6) = 3.472222e-06; aElementMass(26,7) = 4.513889e-05; aElementMass(26,8) = 2.083333e-05; aElementMass(26,9) = 7.523148e-06; aElementMass(26,10) = 9.780093e-05; aElementMass(26,11) = 4.513889e-05; aElementMass(26,12) = 9.780093e-05; aElementMass(26,13) = 1.271412e-03; aElementMass(26,14) = 5.868056e-04; aElementMass(26,15) = 4.513889e-05; aElementMass(26,16) = 5.868056e-04; aElementMass(26,17) = 2.708333e-04; aElementMass(26,18) = 3.472222e-06; aElementMass(26,19) = 4.513889e-05; aElementMass(26,20) = 2.083333e-05; aElementMass(26,21) = 4.513889e-05; aElementMass(26,22) = 5.868056e-04; aElementMass(26,23) = 2.708333e-04; aElementMass(26,24) = 2.083333e-05; aElementMass(26,25) = 2.708333e-04; aElementMass(26,26) = 1.250000e-04;
        }
    }
}
//--------------------------------------------------------------------------------

Mat<real>
Hierarchical_Mesh_Main::MassMatrix_for_L2_projection(
        uint & aDim,
        uint & aPolynomial)
{
    Mat<real> Element_Mass(pow(aPolynomial + 1,aDim),pow(aPolynomial + 1,aDim));
    if ( aDim == 2)
    {
        if ( aPolynomial == 1)
        {
            Element_Mass(0,0) = 1.0/9.0; Element_Mass(0,1) = 1.0/18.0; Element_Mass(0,2) = 1.0/18.0; Element_Mass(0,3) = 1.0/36.0;
            Element_Mass(1,0) = 1.0/18.0; Element_Mass(1,1) = 1.0/9.0; Element_Mass(1,2) = 1.0/36.0; Element_Mass(1,3) = 1.0/18.0;
            Element_Mass(2,0) = 1.0/18.0; Element_Mass(2,1) = 1.0/36.0; Element_Mass(2,2) = 1.0/9.0; Element_Mass(2,3) = 1.0/18.0;
            Element_Mass(3,0) = 1.0/36.0; Element_Mass(3,1) = 1.0/18.0; Element_Mass(3,2) = 1.0/18.0; Element_Mass(3,3) = 1.0/9.0;
        }
        else if ( aPolynomial == 2)
        {
            Element_Mass(0,0) = 0.002500000000000;  Element_Mass(0,1) = 0.005416666666667;  Element_Mass(0,2) = 0.000416666666667;  Element_Mass(0,3) = 0.005416666666667;  Element_Mass(0,4) = 0.011736111111111;  Element_Mass(0,5) = 0.000902777777778;  Element_Mass(0,6) = 0.000416666666667;  Element_Mass(0,7) = 0.000902777777778;  Element_Mass(0,8) = 0.000069444444444;
            Element_Mass(1,0) = 0.005416666666667;  Element_Mass(1,1) = 0.022500000000000;  Element_Mass(1,2) = 0.005416666666667;  Element_Mass(1,3) = 0.011736111111111;  Element_Mass(1,4) = 0.048750000000000;  Element_Mass(1,5) = 0.011736111111111;  Element_Mass(1,6) = 0.000902777777778;  Element_Mass(1,7) = 0.003750000000000;  Element_Mass(1,8) = 0.000902777777778;
            Element_Mass(2,0) = 0.000416666666667;  Element_Mass(2,1) = 0.005416666666667;  Element_Mass(2,2) = 0.002500000000000;  Element_Mass(2,3) = 0.000902777777778;  Element_Mass(2,4) = 0.011736111111111;  Element_Mass(2,5) = 0.005416666666667;  Element_Mass(2,6) = 0.000069444444444;  Element_Mass(2,7) = 0.000902777777778;  Element_Mass(2,8) = 0.000416666666667;
            Element_Mass(3,0) = 0.005416666666667;  Element_Mass(3,1) = 0.011736111111111;  Element_Mass(3,2) = 0.000902777777778;  Element_Mass(3,3) = 0.022500000000000;  Element_Mass(3,4) = 0.048750000000000;  Element_Mass(3,5) = 0.003750000000000;  Element_Mass(3,6) = 0.005416666666667;  Element_Mass(3,7) = 0.011736111111111;  Element_Mass(3,8) = 0.000902777777778;
            Element_Mass(4,0) = 0.011736111111111;  Element_Mass(4,1) = 0.048750000000000;  Element_Mass(4,2) = 0.011736111111111;  Element_Mass(4,3) = 0.048750000000000;  Element_Mass(4,4) = 0.202500000000000;  Element_Mass(4,5) = 0.048750000000000;  Element_Mass(4,6) = 0.011736111111111;  Element_Mass(4,7) = 0.048750000000000;  Element_Mass(4,8) = 0.011736111111111;
            Element_Mass(5,0) = 0.000902777777778;  Element_Mass(5,1) = 0.011736111111111;  Element_Mass(5,2) = 0.005416666666667;  Element_Mass(5,3) = 0.003750000000000;  Element_Mass(5,4) = 0.048750000000000;  Element_Mass(5,5) = 0.022500000000000;  Element_Mass(5,6) = 0.000902777777778;  Element_Mass(5,7) = 0.011736111111111;  Element_Mass(5,8) = 0.005416666666667;
            Element_Mass(6,0) = 0.000416666666667;  Element_Mass(6,1) = 0.000902777777778;  Element_Mass(6,2) = 0.000069444444444;  Element_Mass(6,3) = 0.005416666666667;  Element_Mass(6,4) = 0.011736111111111;  Element_Mass(6,5) = 0.000902777777778;  Element_Mass(6,6) = 0.002500000000000;  Element_Mass(6,7) = 0.005416666666667;  Element_Mass(6,8) = 0.000416666666667;
            Element_Mass(7,0) = 0.000902777777778;  Element_Mass(7,1) = 0.003750000000000;  Element_Mass(7,2) = 0.000902777777778;  Element_Mass(7,3) = 0.011736111111111;  Element_Mass(7,4) = 0.048750000000000;  Element_Mass(7,5) = 0.011736111111111;  Element_Mass(7,6) = 0.005416666666667;  Element_Mass(7,7) = 0.022500000000000;  Element_Mass(7,8) = 0.005416666666667;
            Element_Mass(8,0) = 0.000069444444444;  Element_Mass(8,1) = 0.000902777777778;  Element_Mass(8,2) = 0.000416666666667;  Element_Mass(8,3) = 0.000902777777778;  Element_Mass(8,4) = 0.011736111111111;  Element_Mass(8,5) = 0.005416666666667;  Element_Mass(8,6) = 0.000416666666667;  Element_Mass(8,7) = 0.005416666666667;  Element_Mass(8,8) = 0.002500000000000;
        }
    }
    else if ( aDim == 3)
    {
        if ( aPolynomial == 1)
        {
            Element_Mass(0,0) = 1.0/216.0*8.0; Element_Mass(0,1) = 1.0/216.0*4.0; Element_Mass(0,2) = 1.0/216.0*4.0; Element_Mass(0,3) = 1.0/216.0*2.0; Element_Mass(0,4) = 1.0/216.0*4.0; Element_Mass(0,5) = 1.0/216.0*2.0; Element_Mass(0,6) = 1.0/216.0*2.0; Element_Mass(0,7) = 1.0/216.0*1.0;
            Element_Mass(1,0) = 1.0/216.0*4.0; Element_Mass(1,1) = 1.0/216.0*8.0; Element_Mass(1,2) = 1.0/216.0*2.0; Element_Mass(1,3) = 1.0/216.0*4.0; Element_Mass(1,4) = 1.0/216.0*2.0; Element_Mass(1,5) = 1.0/216.0*4.0; Element_Mass(1,6) = 1.0/216.0*1.0; Element_Mass(1,7) = 1.0/216.0*2.0;
            Element_Mass(2,0) = 1.0/216.0*4.0; Element_Mass(2,1) = 1.0/216.0*2.0; Element_Mass(2,2) = 1.0/216.0*8.0; Element_Mass(2,3) = 1.0/216.0*4.0; Element_Mass(2,4) = 1.0/216.0*2.0; Element_Mass(2,5) = 1.0/216.0*1.0; Element_Mass(2,6) = 1.0/216.0*4.0; Element_Mass(2,7) = 1.0/216.0*2.0;
            Element_Mass(3,0) = 1.0/216.0*2.0; Element_Mass(3,1) = 1.0/216.0*4.0; Element_Mass(3,2) = 1.0/216.0*4.0; Element_Mass(3,3) = 1.0/216.0*8.0; Element_Mass(3,4) = 1.0/216.0*1.0; Element_Mass(3,5) = 1.0/216.0*2.0; Element_Mass(3,6) = 1.0/216.0*2.0; Element_Mass(3,7) = 1.0/216.0*4.0;
            Element_Mass(4,0) = 1.0/216.0*4.0; Element_Mass(4,1) = 1.0/216.0*2.0; Element_Mass(4,2) = 1.0/216.0*2.0; Element_Mass(4,3) = 1.0/216.0*1.0; Element_Mass(4,4) = 1.0/216.0*8.0; Element_Mass(4,5) = 1.0/216.0*4.0; Element_Mass(4,6) = 1.0/216.0*4.0; Element_Mass(4,7) = 1.0/216.0*2.0;
            Element_Mass(5,0) = 1.0/216.0*2.0; Element_Mass(5,1) = 1.0/216.0*4.0; Element_Mass(5,2) = 1.0/216.0*1.0; Element_Mass(5,3) = 1.0/216.0*2.0; Element_Mass(5,4) = 1.0/216.0*4.0; Element_Mass(5,5) = 1.0/216.0*8.0; Element_Mass(5,6) = 1.0/216.0*2.0; Element_Mass(5,7) = 1.0/216.0*4.0;
            Element_Mass(6,0) = 1.0/216.0*2.0; Element_Mass(6,1) = 1.0/216.0*1.0; Element_Mass(6,2) = 1.0/216.0*4.0; Element_Mass(6,3) = 1.0/216.0*2.0; Element_Mass(6,4) = 1.0/216.0*4.0; Element_Mass(6,5) = 1.0/216.0*2.0; Element_Mass(6,6) = 1.0/216.0*8.0; Element_Mass(6,7) = 1.0/216.0*4.0;
            Element_Mass(7,0) = 1.0/216.0*1.0; Element_Mass(7,1) = 1.0/216.0*2.0; Element_Mass(7,2) = 1.0/216.0*2.0; Element_Mass(7,3) = 1.0/216.0*4.0; Element_Mass(7,4) = 1.0/216.0*2.0; Element_Mass(7,5) = 1.0/216.0*4.0; Element_Mass(7,6) = 1.0/216.0*4.0; Element_Mass(7,7) = 1.0/216.0*8.0;
        }
        else if ( aPolynomial == 2)
        {
            Element_Mass(0,0) = 1.250000e-04; Element_Mass(0,1) = 2.708333e-04; Element_Mass(0,2) = 2.083333e-05; Element_Mass(0,3) = 2.708333e-04; Element_Mass(0,4) = 5.868056e-04; Element_Mass(0,5) = 4.513889e-05; Element_Mass(0,6) = 2.083333e-05; Element_Mass(0,7) = 4.513889e-05; Element_Mass(0,8) = 3.472222e-06; Element_Mass(0,9) = 2.708333e-04; Element_Mass(0,10) = 5.868056e-04; Element_Mass(0,11) = 4.513889e-05; Element_Mass(0,12) = 5.868056e-04; Element_Mass(0,13) = 1.271412e-03; Element_Mass(0,14) = 9.780093e-05; Element_Mass(0,15) = 4.513889e-05; Element_Mass(0,16) = 9.780093e-05; Element_Mass(0,17) = 7.523148e-06; Element_Mass(0,18) = 2.083333e-05; Element_Mass(0,19) = 4.513889e-05; Element_Mass(0,20) = 3.472222e-06; Element_Mass(0,21) = 4.513889e-05; Element_Mass(0,22) = 9.780093e-05; Element_Mass(0,23) = 7.523148e-06; Element_Mass(0,24) = 3.472222e-06; Element_Mass(0,25) = 7.523148e-06; Element_Mass(0,26) = 5.787037e-07;
            Element_Mass(1,0) = 2.708333e-04; Element_Mass(1,1) = 1.125000e-03; Element_Mass(1,2) = 2.708333e-04; Element_Mass(1,3) = 5.868056e-04; Element_Mass(1,4) = 2.437500e-03; Element_Mass(1,5) = 5.868056e-04; Element_Mass(1,6) = 4.513889e-05; Element_Mass(1,7) = 1.875000e-04; Element_Mass(1,8) = 4.513889e-05; Element_Mass(1,9) = 5.868056e-04; Element_Mass(1,10) = 2.437500e-03; Element_Mass(1,11) = 5.868056e-04; Element_Mass(1,12) = 1.271412e-03; Element_Mass(1,13) = 5.281250e-03; Element_Mass(1,14) = 1.271412e-03; Element_Mass(1,15) = 9.780093e-05; Element_Mass(1,16) = 4.062500e-04; Element_Mass(1,17) = 9.780093e-05; Element_Mass(1,18) = 4.513889e-05; Element_Mass(1,19) = 1.875000e-04; Element_Mass(1,20) = 4.513889e-05; Element_Mass(1,21) = 9.780093e-05; Element_Mass(1,22) = 4.062500e-04; Element_Mass(1,23) = 9.780093e-05; Element_Mass(1,24) = 7.523148e-06; Element_Mass(1,25) = 3.125000e-05; Element_Mass(1,26) = 7.523148e-06;
            Element_Mass(2,0) = 2.083333e-05; Element_Mass(2,1) = 2.708333e-04; Element_Mass(2,2) = 1.250000e-04; Element_Mass(2,3) = 4.513889e-05; Element_Mass(2,4) = 5.868056e-04; Element_Mass(2,5) = 2.708333e-04; Element_Mass(2,6) = 3.472222e-06; Element_Mass(2,7) = 4.513889e-05; Element_Mass(2,8) = 2.083333e-05; Element_Mass(2,9) = 4.513889e-05; Element_Mass(2,10) = 5.868056e-04; Element_Mass(2,11) = 2.708333e-04; Element_Mass(2,12) = 9.780093e-05; Element_Mass(2,13) = 1.271412e-03; Element_Mass(2,14) = 5.868056e-04; Element_Mass(2,15) = 7.523148e-06; Element_Mass(2,16) = 9.780093e-05; Element_Mass(2,17) = 4.513889e-05; Element_Mass(2,18) = 3.472222e-06; Element_Mass(2,19) = 4.513889e-05; Element_Mass(2,20) = 2.083333e-05; Element_Mass(2,21) = 7.523148e-06; Element_Mass(2,22) = 9.780093e-05; Element_Mass(2,23) = 4.513889e-05; Element_Mass(2,24) = 5.787037e-07; Element_Mass(2,25) = 7.523148e-06; Element_Mass(2,26) = 3.472222e-06;
            Element_Mass(3,0) = 2.708333e-04; Element_Mass(3,1) = 5.868056e-04; Element_Mass(3,2) = 4.513889e-05; Element_Mass(3,3) = 1.125000e-03; Element_Mass(3,4) = 2.437500e-03; Element_Mass(3,5) = 1.875000e-04; Element_Mass(3,6) = 2.708333e-04; Element_Mass(3,7) = 5.868056e-04; Element_Mass(3,8) = 4.513889e-05; Element_Mass(3,9) = 5.868056e-04; Element_Mass(3,10) = 1.271412e-03; Element_Mass(3,11) = 9.780093e-05; Element_Mass(3,12) = 2.437500e-03; Element_Mass(3,13) = 5.281250e-03; Element_Mass(3,14) = 4.062500e-04; Element_Mass(3,15) = 5.868056e-04; Element_Mass(3,16) = 1.271412e-03; Element_Mass(3,17) = 9.780093e-05; Element_Mass(3,18) = 4.513889e-05; Element_Mass(3,19) = 9.780093e-05; Element_Mass(3,20) = 7.523148e-06; Element_Mass(3,21) = 1.875000e-04; Element_Mass(3,22) = 4.062500e-04; Element_Mass(3,23) = 3.125000e-05; Element_Mass(3,24) = 4.513889e-05; Element_Mass(3,25) = 9.780093e-05; Element_Mass(3,26) = 7.523148e-06;
            Element_Mass(4,0) = 5.868056e-04; Element_Mass(4,1) = 2.437500e-03; Element_Mass(4,2) = 5.868056e-04; Element_Mass(4,3) = 2.437500e-03; Element_Mass(4,4) = 1.012500e-02; Element_Mass(4,5) = 2.437500e-03; Element_Mass(4,6) = 5.868056e-04; Element_Mass(4,7) = 2.437500e-03; Element_Mass(4,8) = 5.868056e-04; Element_Mass(4,9) = 1.271412e-03; Element_Mass(4,10) = 5.281250e-03; Element_Mass(4,11) = 1.271412e-03; Element_Mass(4,12) = 5.281250e-03; Element_Mass(4,13) = 2.193750e-02; Element_Mass(4,14) = 5.281250e-03; Element_Mass(4,15) = 1.271412e-03; Element_Mass(4,16) = 5.281250e-03; Element_Mass(4,17) = 1.271412e-03; Element_Mass(4,18) = 9.780093e-05; Element_Mass(4,19) = 4.062500e-04; Element_Mass(4,20) = 9.780093e-05; Element_Mass(4,21) = 4.062500e-04; Element_Mass(4,22) = 1.687500e-03; Element_Mass(4,23) = 4.062500e-04; Element_Mass(4,24) = 9.780093e-05; Element_Mass(4,25) = 4.062500e-04; Element_Mass(4,26) = 9.780093e-05;
            Element_Mass(5,0) = 4.513889e-05; Element_Mass(5,1) = 5.868056e-04; Element_Mass(5,2) = 2.708333e-04; Element_Mass(5,3) = 1.875000e-04; Element_Mass(5,4) = 2.437500e-03; Element_Mass(5,5) = 1.125000e-03; Element_Mass(5,6) = 4.513889e-05; Element_Mass(5,7) = 5.868056e-04; Element_Mass(5,8) = 2.708333e-04; Element_Mass(5,9) = 9.780093e-05; Element_Mass(5,10) = 1.271412e-03; Element_Mass(5,11) = 5.868056e-04; Element_Mass(5,12) = 4.062500e-04; Element_Mass(5,13) = 5.281250e-03; Element_Mass(5,14) = 2.437500e-03; Element_Mass(5,15) = 9.780093e-05; Element_Mass(5,16) = 1.271412e-03; Element_Mass(5,17) = 5.868056e-04; Element_Mass(5,18) = 7.523148e-06; Element_Mass(5,19) = 9.780093e-05; Element_Mass(5,20) = 4.513889e-05; Element_Mass(5,21) = 3.125000e-05; Element_Mass(5,22) = 4.062500e-04; Element_Mass(5,23) = 1.875000e-04; Element_Mass(5,24) = 7.523148e-06; Element_Mass(5,25) = 9.780093e-05; Element_Mass(5,26) = 4.513889e-05;
            Element_Mass(6,0) = 2.083333e-05; Element_Mass(6,1) = 4.513889e-05; Element_Mass(6,2) = 3.472222e-06; Element_Mass(6,3) = 2.708333e-04; Element_Mass(6,4) = 5.868056e-04; Element_Mass(6,5) = 4.513889e-05; Element_Mass(6,6) = 1.250000e-04; Element_Mass(6,7) = 2.708333e-04; Element_Mass(6,8) = 2.083333e-05; Element_Mass(6,9) = 4.513889e-05; Element_Mass(6,10) = 9.780093e-05; Element_Mass(6,11) = 7.523148e-06; Element_Mass(6,12) = 5.868056e-04; Element_Mass(6,13) = 1.271412e-03; Element_Mass(6,14) = 9.780093e-05; Element_Mass(6,15) = 2.708333e-04; Element_Mass(6,16) = 5.868056e-04; Element_Mass(6,17) = 4.513889e-05; Element_Mass(6,18) = 3.472222e-06; Element_Mass(6,19) = 7.523148e-06; Element_Mass(6,20) = 5.787037e-07; Element_Mass(6,21) = 4.513889e-05; Element_Mass(6,22) = 9.780093e-05; Element_Mass(6,23) = 7.523148e-06; Element_Mass(6,24) = 2.083333e-05; Element_Mass(6,25) = 4.513889e-05; Element_Mass(6,26) = 3.472222e-06;
            Element_Mass(7,0) = 4.513889e-05; Element_Mass(7,1) = 1.875000e-04; Element_Mass(7,2) = 4.513889e-05; Element_Mass(7,3) = 5.868056e-04; Element_Mass(7,4) = 2.437500e-03; Element_Mass(7,5) = 5.868056e-04; Element_Mass(7,6) = 2.708333e-04; Element_Mass(7,7) = 1.125000e-03; Element_Mass(7,8) = 2.708333e-04; Element_Mass(7,9) = 9.780093e-05; Element_Mass(7,10) = 4.062500e-04; Element_Mass(7,11) = 9.780093e-05; Element_Mass(7,12) = 1.271412e-03; Element_Mass(7,13) = 5.281250e-03; Element_Mass(7,14) = 1.271412e-03; Element_Mass(7,15) = 5.868056e-04; Element_Mass(7,16) = 2.437500e-03; Element_Mass(7,17) = 5.868056e-04; Element_Mass(7,18) = 7.523148e-06; Element_Mass(7,19) = 3.125000e-05; Element_Mass(7,20) = 7.523148e-06; Element_Mass(7,21) = 9.780093e-05; Element_Mass(7,22) = 4.062500e-04; Element_Mass(7,23) = 9.780093e-05; Element_Mass(7,24) = 4.513889e-05; Element_Mass(7,25) = 1.875000e-04; Element_Mass(7,26) = 4.513889e-05;
            Element_Mass(8,0) = 3.472222e-06; Element_Mass(8,1) = 4.513889e-05; Element_Mass(8,2) = 2.083333e-05; Element_Mass(8,3) = 4.513889e-05; Element_Mass(8,4) = 5.868056e-04; Element_Mass(8,5) = 2.708333e-04; Element_Mass(8,6) = 2.083333e-05; Element_Mass(8,7) = 2.708333e-04; Element_Mass(8,8) = 1.250000e-04; Element_Mass(8,9) = 7.523148e-06; Element_Mass(8,10) = 9.780093e-05; Element_Mass(8,11) = 4.513889e-05; Element_Mass(8,12) = 9.780093e-05; Element_Mass(8,13) = 1.271412e-03; Element_Mass(8,14) = 5.868056e-04; Element_Mass(8,15) = 4.513889e-05; Element_Mass(8,16) = 5.868056e-04; Element_Mass(8,17) = 2.708333e-04; Element_Mass(8,18) = 5.787037e-07; Element_Mass(8,19) = 7.523148e-06; Element_Mass(8,20) = 3.472222e-06; Element_Mass(8,21) = 7.523148e-06; Element_Mass(8,22) = 9.780093e-05; Element_Mass(8,23) = 4.513889e-05; Element_Mass(8,24) = 3.472222e-06; Element_Mass(8,25) = 4.513889e-05; Element_Mass(8,26) = 2.083333e-05;
            Element_Mass(9,0) = 2.708333e-04; Element_Mass(9,1) = 5.868056e-04; Element_Mass(9,2) = 4.513889e-05; Element_Mass(9,3) = 5.868056e-04; Element_Mass(9,4) = 1.271412e-03; Element_Mass(9,5) = 9.780093e-05; Element_Mass(9,6) = 4.513889e-05; Element_Mass(9,7) = 9.780093e-05; Element_Mass(9,8) = 7.523148e-06; Element_Mass(9,9) = 1.125000e-03; Element_Mass(9,10) = 2.437500e-03; Element_Mass(9,11) = 1.875000e-04; Element_Mass(9,12) = 2.437500e-03; Element_Mass(9,13) = 5.281250e-03; Element_Mass(9,14) = 4.062500e-04; Element_Mass(9,15) = 1.875000e-04; Element_Mass(9,16) = 4.062500e-04; Element_Mass(9,17) = 3.125000e-05; Element_Mass(9,18) = 2.708333e-04; Element_Mass(9,19) = 5.868056e-04; Element_Mass(9,20) = 4.513889e-05; Element_Mass(9,21) = 5.868056e-04; Element_Mass(9,22) = 1.271412e-03; Element_Mass(9,23) = 9.780093e-05; Element_Mass(9,24) = 4.513889e-05; Element_Mass(9,25) = 9.780093e-05; Element_Mass(9,26) = 7.523148e-06;
            Element_Mass(10,0) = 5.868056e-04; Element_Mass(10,1) = 2.437500e-03; Element_Mass(10,2) = 5.868056e-04; Element_Mass(10,3) = 1.271412e-03; Element_Mass(10,4) = 5.281250e-03; Element_Mass(10,5) = 1.271412e-03; Element_Mass(10,6) = 9.780093e-05; Element_Mass(10,7) = 4.062500e-04; Element_Mass(10,8) = 9.780093e-05; Element_Mass(10,9) = 2.437500e-03; Element_Mass(10,10) = 1.012500e-02; Element_Mass(10,11) = 2.437500e-03; Element_Mass(10,12) = 5.281250e-03; Element_Mass(10,13) = 2.193750e-02; Element_Mass(10,14) = 5.281250e-03; Element_Mass(10,15) = 4.062500e-04; Element_Mass(10,16) = 1.687500e-03; Element_Mass(10,17) = 4.062500e-04; Element_Mass(10,18) = 5.868056e-04; Element_Mass(10,19) = 2.437500e-03; Element_Mass(10,20) = 5.868056e-04; Element_Mass(10,21) = 1.271412e-03; Element_Mass(10,22) = 5.281250e-03; Element_Mass(10,23) = 1.271412e-03; Element_Mass(10,24) = 9.780093e-05; Element_Mass(10,25) = 4.062500e-04; Element_Mass(10,26) = 9.780093e-05;
            Element_Mass(11,0) = 4.513889e-05; Element_Mass(11,1) = 5.868056e-04; Element_Mass(11,2) = 2.708333e-04; Element_Mass(11,3) = 9.780093e-05; Element_Mass(11,4) = 1.271412e-03; Element_Mass(11,5) = 5.868056e-04; Element_Mass(11,6) = 7.523148e-06; Element_Mass(11,7) = 9.780093e-05; Element_Mass(11,8) = 4.513889e-05; Element_Mass(11,9) = 1.875000e-04; Element_Mass(11,10) = 2.437500e-03; Element_Mass(11,11) = 1.125000e-03; Element_Mass(11,12) = 4.062500e-04; Element_Mass(11,13) = 5.281250e-03; Element_Mass(11,14) = 2.437500e-03; Element_Mass(11,15) = 3.125000e-05; Element_Mass(11,16) = 4.062500e-04; Element_Mass(11,17) = 1.875000e-04; Element_Mass(11,18) = 4.513889e-05; Element_Mass(11,19) = 5.868056e-04; Element_Mass(11,20) = 2.708333e-04; Element_Mass(11,21) = 9.780093e-05; Element_Mass(11,22) = 1.271412e-03; Element_Mass(11,23) = 5.868056e-04; Element_Mass(11,24) = 7.523148e-06; Element_Mass(11,25) = 9.780093e-05; Element_Mass(11,26) = 4.513889e-05;
            Element_Mass(12,0) = 5.868056e-04; Element_Mass(12,1) = 1.271412e-03; Element_Mass(12,2) = 9.780093e-05; Element_Mass(12,3) = 2.437500e-03; Element_Mass(12,4) = 5.281250e-03; Element_Mass(12,5) = 4.062500e-04; Element_Mass(12,6) = 5.868056e-04; Element_Mass(12,7) = 1.271412e-03; Element_Mass(12,8) = 9.780093e-05; Element_Mass(12,9) = 2.437500e-03; Element_Mass(12,10) = 5.281250e-03; Element_Mass(12,11) = 4.062500e-04; Element_Mass(12,12) = 1.012500e-02; Element_Mass(12,13) = 2.193750e-02; Element_Mass(12,14) = 1.687500e-03; Element_Mass(12,15) = 2.437500e-03; Element_Mass(12,16) = 5.281250e-03; Element_Mass(12,17) = 4.062500e-04; Element_Mass(12,18) = 5.868056e-04; Element_Mass(12,19) = 1.271412e-03; Element_Mass(12,20) = 9.780093e-05; Element_Mass(12,21) = 2.437500e-03; Element_Mass(12,22) = 5.281250e-03; Element_Mass(12,23) = 4.062500e-04; Element_Mass(12,24) = 5.868056e-04; Element_Mass(12,25) = 1.271412e-03; Element_Mass(12,26) = 9.780093e-05;
            Element_Mass(13,0) = 1.271412e-03; Element_Mass(13,1) = 5.281250e-03; Element_Mass(13,2) = 1.271412e-03; Element_Mass(13,3) = 5.281250e-03; Element_Mass(13,4) = 2.193750e-02; Element_Mass(13,5) = 5.281250e-03; Element_Mass(13,6) = 1.271412e-03; Element_Mass(13,7) = 5.281250e-03; Element_Mass(13,8) = 1.271412e-03; Element_Mass(13,9) = 5.281250e-03; Element_Mass(13,10) = 2.193750e-02; Element_Mass(13,11) = 5.281250e-03; Element_Mass(13,12) = 2.193750e-02; Element_Mass(13,13) = 9.112500e-02; Element_Mass(13,14) = 2.193750e-02; Element_Mass(13,15) = 5.281250e-03; Element_Mass(13,16) = 2.193750e-02; Element_Mass(13,17) = 5.281250e-03; Element_Mass(13,18) = 1.271412e-03; Element_Mass(13,19) = 5.281250e-03; Element_Mass(13,20) = 1.271412e-03; Element_Mass(13,21) = 5.281250e-03; Element_Mass(13,22) = 2.193750e-02; Element_Mass(13,23) = 5.281250e-03; Element_Mass(13,24) = 1.271412e-03; Element_Mass(13,25) = 5.281250e-03; Element_Mass(13,26) = 1.271412e-03;
            Element_Mass(14,0) = 9.780093e-05; Element_Mass(14,1) = 1.271412e-03; Element_Mass(14,2) = 5.868056e-04; Element_Mass(14,3) = 4.062500e-04; Element_Mass(14,4) = 5.281250e-03; Element_Mass(14,5) = 2.437500e-03; Element_Mass(14,6) = 9.780093e-05; Element_Mass(14,7) = 1.271412e-03; Element_Mass(14,8) = 5.868056e-04; Element_Mass(14,9) = 4.062500e-04; Element_Mass(14,10) = 5.281250e-03; Element_Mass(14,11) = 2.437500e-03; Element_Mass(14,12) = 1.687500e-03; Element_Mass(14,13) = 2.193750e-02; Element_Mass(14,14) = 1.012500e-02; Element_Mass(14,15) = 4.062500e-04; Element_Mass(14,16) = 5.281250e-03; Element_Mass(14,17) = 2.437500e-03; Element_Mass(14,18) = 9.780093e-05; Element_Mass(14,19) = 1.271412e-03; Element_Mass(14,20) = 5.868056e-04; Element_Mass(14,21) = 4.062500e-04; Element_Mass(14,22) = 5.281250e-03; Element_Mass(14,23) = 2.437500e-03; Element_Mass(14,24) = 9.780093e-05; Element_Mass(14,25) = 1.271412e-03; Element_Mass(14,26) = 5.868056e-04;
            Element_Mass(15,0) = 4.513889e-05; Element_Mass(15,1) = 9.780093e-05; Element_Mass(15,2) = 7.523148e-06; Element_Mass(15,3) = 5.868056e-04; Element_Mass(15,4) = 1.271412e-03; Element_Mass(15,5) = 9.780093e-05; Element_Mass(15,6) = 2.708333e-04; Element_Mass(15,7) = 5.868056e-04; Element_Mass(15,8) = 4.513889e-05; Element_Mass(15,9) = 1.875000e-04; Element_Mass(15,10) = 4.062500e-04; Element_Mass(15,11) = 3.125000e-05; Element_Mass(15,12) = 2.437500e-03; Element_Mass(15,13) = 5.281250e-03; Element_Mass(15,14) = 4.062500e-04; Element_Mass(15,15) = 1.125000e-03; Element_Mass(15,16) = 2.437500e-03; Element_Mass(15,17) = 1.875000e-04; Element_Mass(15,18) = 4.513889e-05; Element_Mass(15,19) = 9.780093e-05; Element_Mass(15,20) = 7.523148e-06; Element_Mass(15,21) = 5.868056e-04; Element_Mass(15,22) = 1.271412e-03; Element_Mass(15,23) = 9.780093e-05; Element_Mass(15,24) = 2.708333e-04; Element_Mass(15,25) = 5.868056e-04; Element_Mass(15,26) = 4.513889e-05;
            Element_Mass(16,0) = 9.780093e-05; Element_Mass(16,1) = 4.062500e-04; Element_Mass(16,2) = 9.780093e-05; Element_Mass(16,3) = 1.271412e-03; Element_Mass(16,4) = 5.281250e-03; Element_Mass(16,5) = 1.271412e-03; Element_Mass(16,6) = 5.868056e-04; Element_Mass(16,7) = 2.437500e-03; Element_Mass(16,8) = 5.868056e-04; Element_Mass(16,9) = 4.062500e-04; Element_Mass(16,10) = 1.687500e-03; Element_Mass(16,11) = 4.062500e-04; Element_Mass(16,12) = 5.281250e-03; Element_Mass(16,13) = 2.193750e-02; Element_Mass(16,14) = 5.281250e-03; Element_Mass(16,15) = 2.437500e-03; Element_Mass(16,16) = 1.012500e-02; Element_Mass(16,17) = 2.437500e-03; Element_Mass(16,18) = 9.780093e-05; Element_Mass(16,19) = 4.062500e-04; Element_Mass(16,20) = 9.780093e-05; Element_Mass(16,21) = 1.271412e-03; Element_Mass(16,22) = 5.281250e-03; Element_Mass(16,23) = 1.271412e-03; Element_Mass(16,24) = 5.868056e-04; Element_Mass(16,25) = 2.437500e-03; Element_Mass(16,26) = 5.868056e-04;
            Element_Mass(17,0) = 7.523148e-06; Element_Mass(17,1) = 9.780093e-05; Element_Mass(17,2) = 4.513889e-05; Element_Mass(17,3) = 9.780093e-05; Element_Mass(17,4) = 1.271412e-03; Element_Mass(17,5) = 5.868056e-04; Element_Mass(17,6) = 4.513889e-05; Element_Mass(17,7) = 5.868056e-04; Element_Mass(17,8) = 2.708333e-04; Element_Mass(17,9) = 3.125000e-05; Element_Mass(17,10) = 4.062500e-04; Element_Mass(17,11) = 1.875000e-04; Element_Mass(17,12) = 4.062500e-04; Element_Mass(17,13) = 5.281250e-03; Element_Mass(17,14) = 2.437500e-03; Element_Mass(17,15) = 1.875000e-04; Element_Mass(17,16) = 2.437500e-03; Element_Mass(17,17) = 1.125000e-03; Element_Mass(17,18) = 7.523148e-06; Element_Mass(17,19) = 9.780093e-05; Element_Mass(17,20) = 4.513889e-05; Element_Mass(17,21) = 9.780093e-05; Element_Mass(17,22) = 1.271412e-03; Element_Mass(17,23) = 5.868056e-04; Element_Mass(17,24) = 4.513889e-05; Element_Mass(17,25) = 5.868056e-04; Element_Mass(17,26) = 2.708333e-04;
            Element_Mass(18,0) = 2.083333e-05; Element_Mass(18,1) = 4.513889e-05; Element_Mass(18,2) = 3.472222e-06; Element_Mass(18,3) = 4.513889e-05; Element_Mass(18,4) = 9.780093e-05; Element_Mass(18,5) = 7.523148e-06; Element_Mass(18,6) = 3.472222e-06; Element_Mass(18,7) = 7.523148e-06; Element_Mass(18,8) = 5.787037e-07; Element_Mass(18,9) = 2.708333e-04; Element_Mass(18,10) = 5.868056e-04; Element_Mass(18,11) = 4.513889e-05; Element_Mass(18,12) = 5.868056e-04; Element_Mass(18,13) = 1.271412e-03; Element_Mass(18,14) = 9.780093e-05; Element_Mass(18,15) = 4.513889e-05; Element_Mass(18,16) = 9.780093e-05; Element_Mass(18,17) = 7.523148e-06; Element_Mass(18,18) = 1.250000e-04; Element_Mass(18,19) = 2.708333e-04; Element_Mass(18,20) = 2.083333e-05; Element_Mass(18,21) = 2.708333e-04; Element_Mass(18,22) = 5.868056e-04; Element_Mass(18,23) = 4.513889e-05; Element_Mass(18,24) = 2.083333e-05; Element_Mass(18,25) = 4.513889e-05; Element_Mass(18,26) = 3.472222e-06;
            Element_Mass(19,0) = 4.513889e-05; Element_Mass(19,1) = 1.875000e-04; Element_Mass(19,2) = 4.513889e-05; Element_Mass(19,3) = 9.780093e-05; Element_Mass(19,4) = 4.062500e-04; Element_Mass(19,5) = 9.780093e-05; Element_Mass(19,6) = 7.523148e-06; Element_Mass(19,7) = 3.125000e-05; Element_Mass(19,8) = 7.523148e-06; Element_Mass(19,9) = 5.868056e-04; Element_Mass(19,10) = 2.437500e-03; Element_Mass(19,11) = 5.868056e-04; Element_Mass(19,12) = 1.271412e-03; Element_Mass(19,13) = 5.281250e-03; Element_Mass(19,14) = 1.271412e-03; Element_Mass(19,15) = 9.780093e-05; Element_Mass(19,16) = 4.062500e-04; Element_Mass(19,17) = 9.780093e-05; Element_Mass(19,18) = 2.708333e-04; Element_Mass(19,19) = 1.125000e-03; Element_Mass(19,20) = 2.708333e-04; Element_Mass(19,21) = 5.868056e-04; Element_Mass(19,22) = 2.437500e-03; Element_Mass(19,23) = 5.868056e-04; Element_Mass(19,24) = 4.513889e-05; Element_Mass(19,25) = 1.875000e-04; Element_Mass(19,26) = 4.513889e-05;
            Element_Mass(20,0) = 3.472222e-06; Element_Mass(20,1) = 4.513889e-05; Element_Mass(20,2) = 2.083333e-05; Element_Mass(20,3) = 7.523148e-06; Element_Mass(20,4) = 9.780093e-05; Element_Mass(20,5) = 4.513889e-05; Element_Mass(20,6) = 5.787037e-07; Element_Mass(20,7) = 7.523148e-06; Element_Mass(20,8) = 3.472222e-06; Element_Mass(20,9) = 4.513889e-05; Element_Mass(20,10) = 5.868056e-04; Element_Mass(20,11) = 2.708333e-04; Element_Mass(20,12) = 9.780093e-05; Element_Mass(20,13) = 1.271412e-03; Element_Mass(20,14) = 5.868056e-04; Element_Mass(20,15) = 7.523148e-06; Element_Mass(20,16) = 9.780093e-05; Element_Mass(20,17) = 4.513889e-05; Element_Mass(20,18) = 2.083333e-05; Element_Mass(20,19) = 2.708333e-04; Element_Mass(20,20) = 1.250000e-04; Element_Mass(20,21) = 4.513889e-05; Element_Mass(20,22) = 5.868056e-04; Element_Mass(20,23) = 2.708333e-04; Element_Mass(20,24) = 3.472222e-06; Element_Mass(20,25) = 4.513889e-05; Element_Mass(20,26) = 2.083333e-05;
            Element_Mass(21,0) = 4.513889e-05; Element_Mass(21,1) = 9.780093e-05; Element_Mass(21,2) = 7.523148e-06; Element_Mass(21,3) = 1.875000e-04; Element_Mass(21,4) = 4.062500e-04; Element_Mass(21,5) = 3.125000e-05; Element_Mass(21,6) = 4.513889e-05; Element_Mass(21,7) = 9.780093e-05; Element_Mass(21,8) = 7.523148e-06; Element_Mass(21,9) = 5.868056e-04; Element_Mass(21,10) = 1.271412e-03; Element_Mass(21,11) = 9.780093e-05; Element_Mass(21,12) = 2.437500e-03; Element_Mass(21,13) = 5.281250e-03; Element_Mass(21,14) = 4.062500e-04; Element_Mass(21,15) = 5.868056e-04; Element_Mass(21,16) = 1.271412e-03; Element_Mass(21,17) = 9.780093e-05; Element_Mass(21,18) = 2.708333e-04; Element_Mass(21,19) = 5.868056e-04; Element_Mass(21,20) = 4.513889e-05; Element_Mass(21,21) = 1.125000e-03; Element_Mass(21,22) = 2.437500e-03; Element_Mass(21,23) = 1.875000e-04; Element_Mass(21,24) = 2.708333e-04; Element_Mass(21,25) = 5.868056e-04; Element_Mass(21,26) = 4.513889e-05;
            Element_Mass(22,0) = 9.780093e-05; Element_Mass(22,1) = 4.062500e-04; Element_Mass(22,2) = 9.780093e-05; Element_Mass(22,3) = 4.062500e-04; Element_Mass(22,4) = 1.687500e-03; Element_Mass(22,5) = 4.062500e-04; Element_Mass(22,6) = 9.780093e-05; Element_Mass(22,7) = 4.062500e-04; Element_Mass(22,8) = 9.780093e-05; Element_Mass(22,9) = 1.271412e-03; Element_Mass(22,10) = 5.281250e-03; Element_Mass(22,11) = 1.271412e-03; Element_Mass(22,12) = 5.281250e-03; Element_Mass(22,13) = 2.193750e-02; Element_Mass(22,14) = 5.281250e-03; Element_Mass(22,15) = 1.271412e-03; Element_Mass(22,16) = 5.281250e-03; Element_Mass(22,17) = 1.271412e-03; Element_Mass(22,18) = 5.868056e-04; Element_Mass(22,19) = 2.437500e-03; Element_Mass(22,20) = 5.868056e-04; Element_Mass(22,21) = 2.437500e-03; Element_Mass(22,22) = 1.012500e-02; Element_Mass(22,23) = 2.437500e-03; Element_Mass(22,24) = 5.868056e-04; Element_Mass(22,25) = 2.437500e-03; Element_Mass(22,26) = 5.868056e-04;
            Element_Mass(23,0) = 7.523148e-06; Element_Mass(23,1) = 9.780093e-05; Element_Mass(23,2) = 4.513889e-05; Element_Mass(23,3) = 3.125000e-05; Element_Mass(23,4) = 4.062500e-04; Element_Mass(23,5) = 1.875000e-04; Element_Mass(23,6) = 7.523148e-06; Element_Mass(23,7) = 9.780093e-05; Element_Mass(23,8) = 4.513889e-05; Element_Mass(23,9) = 9.780093e-05; Element_Mass(23,10) = 1.271412e-03; Element_Mass(23,11) = 5.868056e-04; Element_Mass(23,12) = 4.062500e-04; Element_Mass(23,13) = 5.281250e-03; Element_Mass(23,14) = 2.437500e-03; Element_Mass(23,15) = 9.780093e-05; Element_Mass(23,16) = 1.271412e-03; Element_Mass(23,17) = 5.868056e-04; Element_Mass(23,18) = 4.513889e-05; Element_Mass(23,19) = 5.868056e-04; Element_Mass(23,20) = 2.708333e-04; Element_Mass(23,21) = 1.875000e-04; Element_Mass(23,22) = 2.437500e-03; Element_Mass(23,23) = 1.125000e-03; Element_Mass(23,24) = 4.513889e-05; Element_Mass(23,25) = 5.868056e-04; Element_Mass(23,26) = 2.708333e-04;
            Element_Mass(24,0) = 3.472222e-06; Element_Mass(24,1) = 7.523148e-06; Element_Mass(24,2) = 5.787037e-07; Element_Mass(24,3) = 4.513889e-05; Element_Mass(24,4) = 9.780093e-05; Element_Mass(24,5) = 7.523148e-06; Element_Mass(24,6) = 2.083333e-05; Element_Mass(24,7) = 4.513889e-05; Element_Mass(24,8) = 3.472222e-06; Element_Mass(24,9) = 4.513889e-05; Element_Mass(24,10) = 9.780093e-05; Element_Mass(24,11) = 7.523148e-06; Element_Mass(24,12) = 5.868056e-04; Element_Mass(24,13) = 1.271412e-03; Element_Mass(24,14) = 9.780093e-05; Element_Mass(24,15) = 2.708333e-04; Element_Mass(24,16) = 5.868056e-04; Element_Mass(24,17) = 4.513889e-05; Element_Mass(24,18) = 2.083333e-05; Element_Mass(24,19) = 4.513889e-05; Element_Mass(24,20) = 3.472222e-06; Element_Mass(24,21) = 2.708333e-04; Element_Mass(24,22) = 5.868056e-04; Element_Mass(24,23) = 4.513889e-05; Element_Mass(24,24) = 1.250000e-04; Element_Mass(24,25) = 2.708333e-04; Element_Mass(24,26) = 2.083333e-05;
            Element_Mass(25,0) = 7.523148e-06; Element_Mass(25,1) = 3.125000e-05; Element_Mass(25,2) = 7.523148e-06; Element_Mass(25,3) = 9.780093e-05; Element_Mass(25,4) = 4.062500e-04; Element_Mass(25,5) = 9.780093e-05; Element_Mass(25,6) = 4.513889e-05; Element_Mass(25,7) = 1.875000e-04; Element_Mass(25,8) = 4.513889e-05; Element_Mass(25,9) = 9.780093e-05; Element_Mass(25,10) = 4.062500e-04; Element_Mass(25,11) = 9.780093e-05; Element_Mass(25,12) = 1.271412e-03; Element_Mass(25,13) = 5.281250e-03; Element_Mass(25,14) = 1.271412e-03; Element_Mass(25,15) = 5.868056e-04; Element_Mass(25,16) = 2.437500e-03; Element_Mass(25,17) = 5.868056e-04; Element_Mass(25,18) = 4.513889e-05; Element_Mass(25,19) = 1.875000e-04; Element_Mass(25,20) = 4.513889e-05; Element_Mass(25,21) = 5.868056e-04; Element_Mass(25,22) = 2.437500e-03; Element_Mass(25,23) = 5.868056e-04; Element_Mass(25,24) = 2.708333e-04; Element_Mass(25,25) = 1.125000e-03; Element_Mass(25,26) = 2.708333e-04;
            Element_Mass(26,0) = 5.787037e-07; Element_Mass(26,1) = 7.523148e-06; Element_Mass(26,2) = 3.472222e-06; Element_Mass(26,3) = 7.523148e-06; Element_Mass(26,4) = 9.780093e-05; Element_Mass(26,5) = 4.513889e-05; Element_Mass(26,6) = 3.472222e-06; Element_Mass(26,7) = 4.513889e-05; Element_Mass(26,8) = 2.083333e-05; Element_Mass(26,9) = 7.523148e-06; Element_Mass(26,10) = 9.780093e-05; Element_Mass(26,11) = 4.513889e-05; Element_Mass(26,12) = 9.780093e-05; Element_Mass(26,13) = 1.271412e-03; Element_Mass(26,14) = 5.868056e-04; Element_Mass(26,15) = 4.513889e-05; Element_Mass(26,16) = 5.868056e-04; Element_Mass(26,17) = 2.708333e-04; Element_Mass(26,18) = 3.472222e-06; Element_Mass(26,19) = 4.513889e-05; Element_Mass(26,20) = 2.083333e-05; Element_Mass(26,21) = 4.513889e-05; Element_Mass(26,22) = 5.868056e-04; Element_Mass(26,23) = 2.708333e-04; Element_Mass(26,24) = 2.083333e-05; Element_Mass(26,25) = 2.708333e-04; Element_Mass(26,26) = 1.250000e-04;
        }
    }
    return Element_Mass;
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::create_mesh_data_new()
{
    moris::tic  create_mesh_data;

    //Compute number of basis of an element
    uint tNumBasisOfElement = pow( mMeshData.Polynomial + 1, mMeshData.ModelDim );

    // global element IDs for elements on this proc
    mProcData.ElementListOnProc = mRefinement.give_active_elements(mMeshData.ModelDim,
                                                                             mMeshData.Polynomial,
                                                                             mMeshData.NumberOfElementsPerDirection,
                                                                             mMeshData.Level,
                                                                             mProcData.ElementListOnProcInit,
                                                                             mElementData.ElementActive);

    // number of elements on this proc
    uint tSizeActiveElements = mProcData.ElementListOnProc.length();

    // this field will contain the element topology with respect to lagrange nodes in global non-consecutive IDs
    mElementData.FeTopo.set_size( tSizeActiveElements, tNumBasisOfElement );

    // will contain a unique list of Lagrange IDs
    // Local: consecutive index of lagrange element
    // Global: Non-Consecutive ID of lagrange basis
    mProcData.LagrangeListOnProc.set_size( tSizeActiveElements * tNumBasisOfElement, 1, UINT_MAX);

    // Change ordering to the classical order of FE connectivity (For ParaView)
    Mat<uint> tOrder = give_vector_for_reorder();
    real tWeight;
    uint tMaxIdField = 0;
    uint tMaxIdFieldDesign = 0;
    mFieldData.SparseSize = 0; // Define the Sparse size for L2Projection
    Mat<uint> tIdField;
    Mat<uint> tIdFieldDesign;
    Mat<real> tTMatrix;
    Mat<real> tTMatrixDesign;
    Mat<uint> tLagrangeBasisOfElement;
    // If Design field is not activated, create a dummy with zeros
    if ( mBasisData.DesignBSplineActive.count() == 0)
    {
        mBasisData.DesignBSplineActive = mBasisData.BSplineActive;
    }
    // If Design basis are on a lower level active, resize it to the FEM level to have the same size of the bitset,
    // because the function for the T-Matrix wants to check everything
    if ( mBasisData.DesignBSplineActive.size() < mBasisData.BSplineActive.size())
    {
        mBasisData.DesignBSplineActive.resize(mBasisData.BSplineActive.size());
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // in this block, we loop over all active elements and
    //    * create the element topology with respect to (non-consecutive) Lagrange basis
    //    * create a unique list of all active Lagrange basis
    //    * calculate the T-Matrices and ID fields and save them column-wise in temporary arrays for later processing
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // temporary variables for Loop
    Cell< Mat< real > > tTmatrixNodalBased( tSizeActiveElements * tNumBasisOfElement );
    Cell< Mat< uint > > tIdFieldNodalBased( tSizeActiveElements * tNumBasisOfElement );
    Cell< Mat< real > > tTmatrixNodalBasedDesign( tSizeActiveElements * tNumBasisOfElement );
    Cell< Mat< uint > > tIdFieldNodalBasedDesign( tSizeActiveElements * tNumBasisOfElement );

    // Temporary variable for loop, counting global Lagrange basis
    uint tVar = 0;

    // counter for B-Splines
    uint tBsplineCount = 0;

    // loop over all elements on proc
    for ( uint i = 0; i<tSizeActiveElements; i++ )
    {
        //Compute lagrange basis functions of that element
        tLagrangeBasisOfElement =  give_active_lagrange_basis_of_element( mProcData.ElementListOnProc( i ) );

        //Compute the T-matrix for FEM
        this->give_Tmatrix_and_IdField( mProcData.ElementListOnProc( i ) );

        //Map the Bspline basis to lagrange nodes (projection matrix was calculated during create_mesh)
        mElementData.TMatrix = mElementData.TMatrix * mElementData.TProjectionMatrixBasis;

        //Compute the T-matrix for Design
        this->give_Tmatrix_and_IdField_DesignVariables_new( mProcData.ElementListOnProc( i ) );

        //Map the Bspline basis to lagrange nodes
        mElementData.TMatrixDesign = mElementData.TMatrixDesign * mElementData.TProjectionMatrixDesignBasis;

        //Compute the entries for the sparsity matrix
        mFieldData.SparseSize += pow( mElementData.IdFieldDesign.length(), 2 );

        // Save biggest ID field for static field in MTK
        if (  mElementData.IdField.length() > tMaxIdField)
            tMaxIdField =  mElementData.IdField.length();

        // Save biggest ID field of design variables for static field in MTK
        if ( mElementData.IdFieldDesign.length() > tMaxIdFieldDesign)
            tMaxIdFieldDesign = mElementData.IdFieldDesign.length();

        // Loop over lagrange basis functions of this element
        for ( uint j = 0; j < tNumBasisOfElement; j++ )
        {
            // save global Lagrange ID into topology, respecting the ParaView node order
            mElementData.FeTopo( i, tOrder( j ) ) = tLagrangeBasisOfElement( j );

            // remember this Lagrange basis (field is made unique later)
            mProcData.LagrangeListOnProc( tVar ) = tLagrangeBasisOfElement( j );

            // ID field for this basis containing the B-Spline Basis number
            tIdFieldNodalBased( tVar ) = mElementData.IdField;

            //Count the number of total entries in the cell
            tBsplineCount += mElementData.IdField.length();

            // save T-Matrix column
            tTmatrixNodalBased( tVar ) = mElementData.TMatrix.col( j );

            // ID field for this basis containing the B-Spline Basis number (design field)
            tIdFieldNodalBasedDesign( tVar ) = mElementData.IdFieldDesign;

            // save T-matrix column (design field)
            tTmatrixNodalBasedDesign( tVar ) = mElementData.TMatrixDesign.col( j );
            tVar++;
        }
    }

    // Find unique nodes
    Mat<uint> tUnique_list = find_unique( mProcData.LagrangeListOnProc );

    // Make Lagrange basis list unique
    mProcData.LagrangeListOnProc = unique( mProcData.LagrangeListOnProc );

    // This data field is a legacy from the old create_mesh_data function. If lagrange basis are used,
    // it is identical to LagrangeListOnProc
    //mBasisData.NodalLocaltoGlobal = mProcData.LagrangeListOnProc;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // in this block, we loop over all B-Splines that we create a unique list of B-Spline basis
    // by copying all entries from the Cell tIdFieldNodalBased into mBasisData.BsplineBasisList,
    // and making it unique
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Temporary variable to compute the length of an elemental Id field
    uint tCellIdFieldLength = 0;

    // Temporary variable, which is needed to copy data from rows
    uint tCopyLower = 0;
    uint tCopyUpper = 0;

    // reserve memory for B-Spline List
    mBasisData.BsplineBasisList.set_size( tBsplineCount, 1, UINT_MAX);

    // loop over all entries in the cell tIdFieldNodalBased, and copy data
    for ( uint i = 0; i < tIdFieldNodalBased.size(); i++ )
    {
        //Compute the upper point for copying data
        tCellIdFieldLength = tIdFieldNodalBased( i ).length();
        tCopyUpper = tCopyLower + tCellIdFieldLength;
        mBasisData.BsplineBasisList.rows( tCopyLower, tCopyUpper - 1 )
            = tIdFieldNodalBased( i ).rows( 0, tCellIdFieldLength - 1 );

        //Lower point is now the old upper point
        tCopyLower = tCopyUpper;
    }

    // make list unique
    mBasisData.BsplineBasisList = unique( mBasisData.BsplineBasisList );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // in this block, we save the node coordinates for the active Lagrange basis in consecutive order
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    //List of coordinates
    mElementData.FeCoord.set_size( mProcData.LagrangeListOnProc.length(), mMeshData.ModelDim );

    // temporary vector to store coordinate information
    Mat<real> tCoordinates;

    // loop over all active Lagrange basis
    for ( uint i = 0; i < mProcData.LagrangeListOnProc.length(); i++ )
    {
        // Calculate coordinate
        tCoordinates =  give_coordinate_from_lagrange_basis( mProcData.LagrangeListOnProc( i ) );

        // save coordinate into FeCoord
        mElementData.FeCoord(i,0) = tCoordinates( 0 );
        mElementData.FeCoord(i,1) = tCoordinates( 1 );
        if ( mMeshData.ModelDim == 3)
            mElementData.FeCoord(i,2) = tCoordinates( 2 );
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // in this block, we create the T-Matrices
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Tmatrix entries are saved in a Nodal field + Number of supported basis functions
    mBasisData.TMatrixNodalField.set_size(mProcData.LagrangeListOnProc.length(),tMaxIdField + 1,0);

    // ID Field entries are saved in a Nodal field + Number of supported basis functions
    mBasisData.IdFieldNodalField.set_size(mProcData.LagrangeListOnProc.length(),tMaxIdField + 1,0);

    // Tmatrix entries are saved in a Nodal field + Number of supported basis functions
    mBasisData.TMatrixDesignNodalField.set_size(mProcData.LagrangeListOnProc.length(),tMaxIdFieldDesign + 1,0);

    // ID Field entries are saved in a Nodal field + Number of supported basis functions
    mBasisData.IdFieldDesignNodalField.set_size(mProcData.LagrangeListOnProc.length(),tMaxIdFieldDesign + 1,0);
    //mBasisData.IdFieldDesignNodalList.set_size(mProcData.LagrangeListOnProc.length()*tMaxIdFieldDesign,1,0);

    // B-Spline counter for IdFieldDesignNodalList
    uint tVarDesignNodalList = 0;

    // we will resize TMatrixNodalField and IdFieldNodalField later, so safe the max number of B-Spline basis functions
    tMaxIdField = 0;

    // same for TMatrixDesignNodalField,IdFieldDesignNodalField
    tMaxIdFieldDesign = 0;

    // loop over all active Lagrange basis
    for ( uint i = 0; i < mProcData.LagrangeListOnProc.length(); i++)
    {
        // we start with the B-Spline basis counter=1 (0 contains the number of unique B-Spline bases
        tVar = 1;

        // get the ID of the T-Matrix for each node:
        // using tUnique_list, we pick the correct entry form Lagrange Basis i in consecutive ordering
        // we copy the B-Spline indices into the ID field,
        // and the T-Matrix row into the corresponding row of the Lagrange basis
        for ( uint j = 0; j< (tIdFieldNodalBased(tUnique_list(i))).length(); j++)
        {
            if ( (tTmatrixNodalBased(tUnique_list(i)))(j) > 0.0)
            {
                    // copy entries from T-Matrix into corresponding row of Lagrange Basis
                    mBasisData.TMatrixNodalField(i,tVar) = (tTmatrixNodalBased(tUnique_list(i)))(j);
                    mBasisData.IdFieldNodalField(i,tVar) = (tIdFieldNodalBased(tUnique_list(i)))(j);

                    // increment B-Spline counter
                    tVar++;
            }
        }

        // update maximum number of B-Splines
        if ( tVar > tMaxIdField)
            tMaxIdField = tVar;

        // safe number if B-Splines in first row of T-Matrix and ID field
        mBasisData.TMatrixNodalField(i,0) = tVar - 1;
        mBasisData.IdFieldNodalField(i,0) = tVar - 1;

        // now do the same for the B-Splines in the design field, start counting with 1
        tVar = 1;

        // weight for T-Matrix averaging. Tobi says this weighting functionality is deprecated
        // (mSettings.PerformNormalization is set to false by default)
        tWeight = 0.0;

        // using tUnique_list, we pick the correct entry form Lagrange Basis i in consecutive ordering
        // we copy the B-Spline indices into the design ID field,
        for ( uint j = 0; j< (tIdFieldNodalBasedDesign(tUnique_list(i))).length(); j++)
        {
            if ( (tTmatrixNodalBasedDesign(tUnique_list(i)))(j) > 0.0)
            {

                // add value to the weight for the T-Matrix
                tWeight += (tTmatrixNodalBasedDesign(tUnique_list(i)))(j);

                // copy entries from T-Matrix into corresponding row of Lagrange Basis
                mBasisData.TMatrixDesignNodalField(i,tVar) = (tTmatrixNodalBasedDesign(tUnique_list(i)))(j);
                mBasisData.IdFieldDesignNodalField(i,tVar) = (tIdFieldNodalBasedDesign(tUnique_list(i)))(j);

                // increment counter for T-Matrix and ID field
                tVar++;

                //Save all Ids from the design variables in a list (for outputting purposes)
                //mBasisData.IdFieldDesignNodalList(tVarDesignNodalList) = (tIdFieldNodalBasedDesign(tUnique_list(i)))(j);

                // increment counter for IdFieldDesignNodalList
                tVarDesignNodalList++;
            }
        }

        // perform the normalization if PerformNormalization is set (off by default)
        if ( mSettings.PerformNormalization )
        {
            mBasisData.TMatrixDesignNodalField.row(i) /= tWeight;
        }

        // increment max number if B-Spline basis
        if ( tVar > tMaxIdFieldDesign)
                tMaxIdFieldDesign = tVar;

        // save number of B-Spline basis for Lagrange Basis (i) in first column
        mBasisData.TMatrixDesignNodalField(i,0) = tVar - 1;
        mBasisData.IdFieldDesignNodalField(i,0) = tVar - 1;
    } // end of loop over all Lagrange basis
    // adapt size of matrices. Remember that we started counting with tVar=1,
    // since first column contains number of B-Splines
    mBasisData.TMatrixNodalField.resize(mProcData.LagrangeListOnProc.length(),tMaxIdField);
    mBasisData.IdFieldNodalField.resize(mProcData.LagrangeListOnProc.length(),tMaxIdField);
    mBasisData.TMatrixDesignNodalField.resize(mProcData.LagrangeListOnProc.length(),tMaxIdFieldDesign);
    mBasisData.IdFieldDesignNodalField.resize(mProcData.LagrangeListOnProc.length(),tMaxIdFieldDesign);

    // adapt B-Spline list for design field, and make array unique
    //mBasisData.IdFieldDesignNodalList.resize(tVarDesignNodalList,1);
    //mBasisData.IdFieldDesignNodalList = unique(mBasisData.IdFieldDesignNodalList);
    // create maps for node numbering
    if (mBasisData.LagrangeToBSplineMap.length() != mProcData.LagrangeListOnProc.length())
    {
        create_maps();
    }

    //Use a filter if filter radius is larger then zero
    if ( mSettings.FilterRadius > 0.0 )
    {
        moris::tic updateidandtmatrixdesign;
        mLagrangeFilter.Update_IDandTMatrix_design(
                mMeshData.ModelDim,
                mMeshData.PolynomialDesign,
                mMeshData.NumberOfElementsPerDirection,
                mMeshData.Level,
                mSettings.FilterRadius,
                mElementData.ElementActiveDesign,
                mBasisData.DesignBSplineActive,
                mBasisData.DesignLagrangeBasisActive,
                mBasisData.NumberOfActiveDesignLagrangeBasis,
                //mBasisData.DesignBSplineListMap,
                mMeshData.Dimensions,
                mMeshData.DimensionsOffset,
                mMeshData.Dimensions_Orig,
                mBasisData.IdFieldDesignNodalField,
                mBasisData.TMatrixDesignNodalField,
                mSettings.FilterLevelWeight,
                mSettings.PerformNormalization,
                mSettings.PointOfOrigin,
                mProcData.LagrangeListOnProc
        );

        real tElapsedupdateidandtmatrixdesign = updateidandtmatrixdesign.toc<moris::chronos::seconds>().wall;
        std::fprintf(stdout,"Time for using the filter: %f [sec]\n",tElapsedupdateidandtmatrixdesign);
    }

    //Determine owner if MPI
    if ( par_size() > 1 )
    {
        mBasisData.LagrangeNodeProcOwner.set_size( mProcData.LagrangeListOnProc.length(), 1, UINT_MAX );
        for ( uint i = 0; i < mProcData.LagrangeListOnProc.length(); i++ )
        {
            mBasisData.LagrangeNodeProcOwner( i ) = give_lagrange_basis_owner( mProcData.LagrangeListOnProc( i ) );
        }
    }

    real tElapsedcreate_mesh_data = create_mesh_data.toc<moris::chronos::seconds>().wall;
    if ( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for creating mesh data : %f [sec]\n",tElapsedcreate_mesh_data);
}

//--------------------------------------------------------------------------------

/* @TODO consider renaming this function to give_basis_vicinity_for_bounding_box */
Mat<uint>
Hierarchical_Mesh_Main::give_basis_vicinity_for_triangle(
        Mat<real> const & aCoordinate0,
        Mat<real> const & aCoordinate1)
{
    uint tElement0 = mHMRElement.give_element_for_coordinate_on_level_zero(mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection,mMeshData.Dimensions,mMeshData.DimensionsOffset,mSettings.PointOfOrigin,aCoordinate0);
    uint tElement1 = mHMRElement.give_element_for_coordinate_on_level_zero(mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection,mMeshData.Dimensions,mMeshData.DimensionsOffset,mSettings.PointOfOrigin,aCoordinate1);
    Mat<uint> tElementPos0 = mBaseElement.give_position_of_element(tElement0,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
    Mat<uint> tElementPos1 = mBaseElement.give_position_of_element(tElement1,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
    uint tVar = 0;
    Mat<uint> tBasisOfElement;
    Mat<uint> tListOfBasis;
    Mat<uint> tListOfPossibleElements;
    Mat<uint> tIJKPosition(mMeshData.ModelDim,1,0);
    Mat<uint> tActiveElements;
    if ( mMeshData.ModelDim == 2)
    {
        tListOfPossibleElements.set_size((tElementPos1( 0 ) + 1-tElementPos0( 0 ))*(tElementPos1( 1 ) + 1-tElementPos0( 1 )),1,UINT_MAX);
        for ( uint i = tElementPos0( 0 ); i <= tElementPos1( 0 ); i++)
        {
            for ( uint j = tElementPos0( 1 ); j <= tElementPos1( 1 ); j++)
            {
                tIJKPosition( 0 ) = i; tIJKPosition( 1 ) = j;
                tListOfPossibleElements(tVar) = mBaseElement.give_element_of_position(0,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection,tIJKPosition); // Check for elements on level zero
                tVar++;
            }
        }
    }
    else if ( mMeshData.ModelDim == 3)
    {
        tListOfPossibleElements.set_size((tElementPos1( 0 ) + 1-tElementPos0( 0 ))*(tElementPos1( 1 ) + 1-tElementPos0( 1 ))*(tElementPos1( 2 ) + 1-tElementPos0( 2 )),1,UINT_MAX);
        for ( uint i = tElementPos0( 0 ); i <= tElementPos1( 0 ); i++)
        {
            for ( uint j = tElementPos0( 1 ); j <= tElementPos1( 1 ); j++)
            {
                for ( uint k = tElementPos0( 2 ); k <= tElementPos1( 2 ); k++)
                {
                    tIJKPosition( 0 ) = i; tIJKPosition( 1 ) = j; tIJKPosition( 2 ) = k;
                    tListOfPossibleElements(tVar) = mBaseElement.give_element_of_position(0,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection,tIJKPosition);  // Check for elements on level zero
                    tVar++;
                }
            }
        }
    }
    tListOfPossibleElements = unique(tListOfPossibleElements);

    tActiveElements.set_size(tListOfPossibleElements.length()*pow(pow(2,mMeshData.ModelDim),(mMeshData.Level + 1)),1,0);
    uint tVarc = tListOfPossibleElements.length();
    tListOfPossibleElements.resize(tListOfPossibleElements.length()*pow(pow(2,mMeshData.ModelDim),(mMeshData.Level + 1)),1);
    Mat<uint> tChildren;
    tVar = 0;
    uint tVarb = 0;
    while( tListOfPossibleElements(tVarb) > 0)
    {
        if ( mElementData.ElementActive.test(tListOfPossibleElements(tVarb)) == 1 )
        {
            tActiveElements(tVar) = tListOfPossibleElements(tVarb);
            tVar++;
        }
        else
        {
            tChildren = mBaseElement.give_children_of_element(tListOfPossibleElements(tVarb),mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
            for ( uint i = 0; i < tChildren.length(); i++)
            {
                tListOfPossibleElements(tVarc) = tChildren(i);
                tVarc++;
            }
        }
        tVarb++;
    }
    tActiveElements.resize(tVar,1);

    tListOfBasis.set_size(tActiveElements.length()*pow(mMeshData.Polynomial + 1,mMeshData.ModelDim),1,UINT_MAX);
    for( uint i = 0; i < tActiveElements.length(); i++)
    {
        tBasisOfElement = mHMRElement.give_basis_of_element(tActiveElements(i),mMeshData.ModelDim,mMeshData.Polynomial,mMeshData.NumberOfElementsPerDirection);
        tListOfBasis.rows(i*tBasisOfElement.length(),(i + 1)*(tBasisOfElement.length()) - 1) = tBasisOfElement.rows(0,tBasisOfElement.length() - 1);
    }
    tListOfBasis = unique(tListOfBasis);
    return tListOfBasis;
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::floodfill_for_STL(
        Mat<real> & aNodalField)
{
    Mat<uint> tElementListOnProc = mRefinement.give_active_elements(mMeshData.ModelDim,mMeshData.Polynomial,mMeshData.NumberOfElementsPerDirection,mMeshData.Level,mProcData.ElementListOnProcInit,mElementData.ElementActive);
    uint tVarb = 0;
    uint tVarc = 1;
    real tVal = 1.0;
    //    uint tSwitch = 0;
    Mat<uint> tElementNeighbor;
    Mat<uint> tBasisOfElement;
    Mat<uint> tPossibleElements;
    BoostBitset tActiveElements = mElementData.ElementActive;
    //    BoostBitset tActiveBasis = mBasisData.BSplineActive;
    Mat<uint> tCheckNeighbors;
    real tSum;
    Mat<real> tElementValues(pow(mMeshData.Polynomial + 1,mMeshData.ModelDim),1,0);
    if ( mMeshData.ModelDim == 2)
    {
        Mat<uint> tCheckdirections2D = {{0, 3},{0, 1},{3, 2},{1, 2}};
        tCheckNeighbors = tCheckdirections2D;
    }
    else if ( mMeshData.ModelDim == 3)
    {
        Mat<uint> tCheckdirections3D = {{0, 3, 4},{0, 1, 4},{3, 2, 4},{1, 2, 4}, {0, 3, 5},{0, 1, 5},{3, 2, 5},{1, 2, 5}};
        tCheckNeighbors = tCheckdirections3D;
    }
    while( tActiveElements.count() > 0)
    {
        tPossibleElements.set_size(tActiveElements.count() + 1,1,0);
        tPossibleElements( 0 ) = tActiveElements.find_first(); // Start with the first active element;
        tActiveElements.reset(tPossibleElements( 0 ));
        tVarb = 0;
        tVarc = 1;
        while( tPossibleElements(tVarb) > 0)
        {
            //            aNodalField.print("aNodalField");
            tBasisOfElement = mHMRElement.give_basis_of_element(tPossibleElements(tVarb),mMeshData.ModelDim,mMeshData.Polynomial,mMeshData.NumberOfElementsPerDirection);
            tElementNeighbor = mHMRElement.give_active_face_neighbor_of_element(tPossibleElements(tVarb),mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection,mMeshData.Level,mElementData.ElementActive);
            //            tElementNeighbor.print("tElementNeighbor");
            //            tCheckNeighbors.print("tCheckNeighbors");
            for ( uint i = 0; i < tBasisOfElement.length(); i++)
            {
                tElementValues(i) = aNodalField(tBasisOfElement(i)) ;
            }
            tSum = sum(tElementValues);
            if ( tSum == 0.0 )
            {
                tVal = 1;
            }
            else if ( tElementValues.max() > 0.0 )
            {
                tVal = 1;
            }
            else if ( tElementValues.min() < 0.0 )
            {
                tVal =  - 1;
            }

            for ( uint i = 0; i < tBasisOfElement.length(); i++)
            {
                //                tActiveBasis.reset(tBasisOfElement(i));
                if ( aNodalField(tBasisOfElement(i)) == 0.0 )
                {
                    aNodalField(tBasisOfElement(i)) = tVal;
                    for ( uint j = 0; j < tElementNeighbor.n_rows(); j++)
                    {
                        //                        std::cout << " tElementNeighbor(j,1) " << tElementNeighbor(j,1) << " tCheckNeighbors(i,0) " << tCheckNeighbors(i,0) << " tElementNeighbor(j,0) " << tElementNeighbor(j,0) << " tActiveElements.test(tElementNeighbor(j,0))  " << tActiveElements.test(tElementNeighbor(j,0))  << std::endl;
                        if ( tElementNeighbor(j,1) == tCheckNeighbors(i,0) && tActiveElements.test(tElementNeighbor(j,0)) == 1)
                        {
                            tPossibleElements(tVarc) = tElementNeighbor(j,0);
                            tVarc++;
                            tActiveElements.reset(tElementNeighbor(j,0));
                        }
                        //                        std::cout << " tElementNeighbor(j,1) " << tElementNeighbor(j,1) << " tCheckNeighbors(i,1) " << tCheckNeighbors(i,1) << " tElementNeighbor(j,0) " << tElementNeighbor(j,0) << " tActiveElements.test(tElementNeighbor(j,0))  " << tActiveElements.test(tElementNeighbor(j,0))  << std::endl;
                        if ( tElementNeighbor(j,1) == tCheckNeighbors(i,1) && tActiveElements.test(tElementNeighbor(j,0)) == 1)
                        {
                            tPossibleElements(tVarc) = tElementNeighbor(j,0);
                            tVarc++;
                            tActiveElements.reset(tElementNeighbor(j,0));
                        }
                        if ( mMeshData.ModelDim == 3 && tElementNeighbor(j,1) == tCheckNeighbors(i,2) && tActiveElements.test(tElementNeighbor(j,0)) == 1)
                        {
                            tPossibleElements(tVarc) = tElementNeighbor(j,0);
                            tVarc++;
                            tActiveElements.reset(tElementNeighbor(j,0));
                        }
                    }
                }
            }
            tVarb++;
        }
    }
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::copy_parser_information(
        Hierarchical_Mesh_Main& aHierarchicalMeshIn,
        Hierarchical_Mesh_Main& aHierarchicalMeshOut)
{

    // copy settings struct
    aHierarchicalMeshOut.mSettings = aHierarchicalMeshIn.mSettings;

    /* aHierarchicalMeshOut.mSettings.Refinement                   = aHierarchicalMeshIn.mSettings.Refinement;
    aHierarchicalMeshOut.mSettings.BufferElements               = aHierarchicalMeshIn.mSettings.BufferElements;
    aHierarchicalMeshOut.mSettings.FeatureResolution            = aHierarchicalMeshIn.mSettings.FeatureResolution;
    aHierarchicalMeshOut.mSettings.MaxDesignLevelOfRefinement   = aHierarchicalMeshIn.mSettings.MaxDesignLevelOfRefinement;
    aHierarchicalMeshOut.mSettings.MaxLevelSetLevelOfRefinement = aHierarchicalMeshIn.mSettings.MaxLevelSetLevelOfRefinement;
    aHierarchicalMeshOut.mSettings.FilterRadius                 = aHierarchicalMeshIn.mSettings.FilterRadius;
    aHierarchicalMeshOut.mSettings.ElementLowerBound            = aHierarchicalMeshIn.mSettings.ElementLowerBound;
    aHierarchicalMeshOut.mSettings.ElementIncrementBound        = aHierarchicalMeshIn.mSettings.ElementIncrementBound;
    aHierarchicalMeshOut.mSettings.ElementUpperBound            = aHierarchicalMeshIn.mSettings.ElementUpperBound;
    aHierarchicalMeshOut.mSettings.FeatureLowerBound            = aHierarchicalMeshIn.mSettings.FeatureLowerBound;
    aHierarchicalMeshOut.mSettings.MinimumRefinement            = aHierarchicalMeshIn.mSettings.MinimumRefinement;
    aHierarchicalMeshOut.mSettings.PerformNormalization         = aHierarchicalMeshIn.mSettings.PerformNormalization;
    aHierarchicalMeshOut.mSettings.FilterLevelWeight            = aHierarchicalMeshIn.mSettings.FilterLevelWeight;
    aHierarchicalMeshOut.mSettings.AdaptiveBufferLayer          = aHierarchicalMeshIn.mSettings.AdaptiveBufferLayer;
    aHierarchicalMeshOut.mSettings.Staircasebuffer              = aHierarchicalMeshIn.mSettings.Staircasebuffer;
    aHierarchicalMeshOut.mSettings.TruncatedBsplines            = aHierarchicalMeshIn.mSettings.TruncatedBsplines;
    aHierarchicalMeshOut.mSettings.LvlSetMethod                 = aHierarchicalMeshIn.mSettings.LvlSetMethod;
    aHierarchicalMeshOut.mSettings.ProjectionFlag               = aHierarchicalMeshIn.mSettings.ProjectionFlag;
    aHierarchicalMeshOut.mSettings.SimpExp                      = aHierarchicalMeshIn.mSettings.SimpExp;
    aHierarchicalMeshOut.mSettings.ElementEdgeWidth             = aHierarchicalMeshIn.mSettings.ElementEdgeWidth;
    aHierarchicalMeshOut.mSettings.LSthresh                     = aHierarchicalMeshIn.mSettings.LSthresh;
    aHierarchicalMeshOut.mSettings.LSscale                      = aHierarchicalMeshIn.mSettings.LSscale;
    aHierarchicalMeshOut.mSettings.DensityShift                 = aHierarchicalMeshIn.mSettings.DensityShift;
    aHierarchicalMeshOut.mSettings.Prjbeta                      = aHierarchicalMeshIn.mSettings.Prjbeta;
    aHierarchicalMeshOut.mSettings.SDFField                               = aHierarchicalMeshIn.mSettings.SDFField;
    aHierarchicalMeshOut.mSettings.SDFCandidateSearchDepth                = aHierarchicalMeshIn.mSettings.SDFCandidateSearchDepth;
    aHierarchicalMeshOut.mSettings.SDFObjectFiles                         = aHierarchicalMeshIn.mSettings.SDFObjectFiles;
    aHierarchicalMeshOut.mSettings.MaxSurfaceRefinement         = aHierarchicalMeshIn.mSettings.MaxSurfaceRefinement;
    aHierarchicalMeshOut.mSettings.MaxVolumeRefinement          = aHierarchicalMeshIn.mSettings.MaxVolumeRefinement;
    aHierarchicalMeshOut.mSettings.DesignVariablesOnInitialMesh = aHierarchicalMeshIn.mSettings.DesignVariablesOnInitialMesh;
    aHierarchicalMeshOut.mSettings.RefineBaseOnExoFile                    = aHierarchicalMeshIn.mSettings.RefineBaseOnExoFile;
    aHierarchicalMeshOut.mSettings.PointOfOrigin                  = aHierarchicalMeshIn.mSettings.PointOfOrigin;
    aHierarchicalMeshOut.mSettings.InputOutputFilesForMorisAreBinary              = aHierarchicalMeshIn.mSettings.InputOutputFilesForMorisAreBinary;
    aHierarchicalMeshOut.mSettings.SDFAlpha                               = aHierarchicalMeshIn.mSettings.SDFAlpha;
    aHierarchicalMeshOut.mSettings.SDFBeta                                = aHierarchicalMeshIn.mSettings.SDFBeta;
    aHierarchicalMeshOut.mSettings.UseSymmetry                            = aHierarchicalMeshIn.mSettings.UseSymmetry;
    aHierarchicalMeshOut.mSettings.SymmetryPlane                          = aHierarchicalMeshIn.mSettings.SymmetryPlane; */

    // resize input vectors according to given number of object files
    aHierarchicalMeshOut.mSettings.SDFAlpha.resize(aHierarchicalMeshOut.mSettings.SDFObjectFiles.size(), 1);
    aHierarchicalMeshOut.mSettings.SDFBeta.resize(aHierarchicalMeshOut.mSettings.SDFObjectFiles.size(), 1);

    //Set data in private struct
    //Receive number of elements per direction from dummy HMR
    Mat<uint> tNumberOfElementsPerDirection = aHierarchicalMeshIn.give_set_number_elements_per_direction();
    //set number of elements per direction for HMR
    aHierarchicalMeshOut.set_number_elements_per_direction( tNumberOfElementsPerDirection );

    //Receive domain dimensions from dummy HMR
    Mat<real> tDomainDimensions = aHierarchicalMeshIn.give_domain_dimensions();
    //set domain dimensions for HMR
    aHierarchicalMeshOut.set_domain_dimensions( tDomainDimensions );

    //Receive side set from dummy HMR
    Mat<uint> tSideSet = aHierarchicalMeshIn.give_side_set();
    //Set side set for HMR
    aHierarchicalMeshOut.set_side_set( tSideSet );

    //Receive domain which needs to be deactivated from the beginning from dummy HMR
    Mat<real> tDeactiveDomain = aHierarchicalMeshIn.give_deactive_domain();
    //Set domain which needs to be deactivated from the beginning for HMR
    aHierarchicalMeshOut.set_deactive_domain( tDeactiveDomain );

    //Activate timings for several functions
    aHierarchicalMeshOut.set_timing_tests();

    // print parameters to screen
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.Refinement                   = %d\n",aHierarchicalMeshOut.mSettings.Refinement);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.BufferElements               = %d\n",aHierarchicalMeshOut.mSettings.BufferElements);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.FeatureResolution            = %d\n",aHierarchicalMeshOut.mSettings.FeatureResolution);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.MaxDesignLevelOfRefinement   = %d\n",aHierarchicalMeshOut.mSettings.MaxDesignLevelOfRefinement);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.MaxLevelSetLevelOfRefinement = %d\n",aHierarchicalMeshOut.mSettings.MaxLevelSetLevelOfRefinement);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.MaxSurfaceRefinement         = %d\n",aHierarchicalMeshOut.mSettings.MaxSurfaceRefinement);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.MaxVolumeRefinement          = %d\n",aHierarchicalMeshOut.mSettings.MaxVolumeRefinement);

    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.FilterRadius                 = %f\n",aHierarchicalMeshOut.mSettings.FilterRadius);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.DensityLowerBound            = %f\n",aHierarchicalMeshOut.mSettings.ElementLowerBound);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.DensityIncrementBound        = %f\n",aHierarchicalMeshOut.mSettings.ElementIncrementBound);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.DensityUpperBound            = %f\n",aHierarchicalMeshOut.mSettings.ElementUpperBound);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.FeatureLowerBound            = %f\n",aHierarchicalMeshOut.mSettings.FeatureLowerBound);

    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.MinimumRefinement            = %d\n",aHierarchicalMeshOut.mSettings.MinimumRefinement);

    std::fprintf(stdout,"    aTensorMeshInputData.nex                                    = %d\n",aHierarchicalMeshOut.mMeshData.NumberOfElementsPerDirection( 0 ) );
    std::fprintf(stdout,"    aTensorMeshInputData.ney                                    = %d\n",aHierarchicalMeshOut.mMeshData.NumberOfElementsPerDirection( 1 ) );
    std::fprintf(stdout,"    aTensorMeshInputData.nez                                    = %d\n",aHierarchicalMeshOut.mMeshData.NumberOfElementsPerDirection( 2 ) );

    std::fprintf(stdout,"    aTensorMeshInputData.lengthorig                             = %f\n",aHierarchicalMeshOut.mMeshData.Dimensions( 0 ) );
    std::fprintf(stdout,"    aTensorMeshInputData.heightorig                             = %f\n",aHierarchicalMeshOut.mMeshData.Dimensions( 1 ) );
    std::fprintf(stdout,"    aTensorMeshInputData.thicknessorig                          = %f\n",aHierarchicalMeshOut.mMeshData.Dimensions( 2 ) );

    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.Elementwidth                 = %f\n",aHierarchicalMeshOut.mSettings.ElementEdgeWidth);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.InitialElementwidth          = %f\n",aHierarchicalMeshOut.mSettings.InitialElementEdgeWidth);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.sLSthresh                    = %f\n",aHierarchicalMeshOut.mSettings.LSthresh);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.sLSscale                     = %f\n",aHierarchicalMeshOut.mSettings.LSscale);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.DensityShift                 = %f\n",aHierarchicalMeshOut.mSettings.DensityShift);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.SimpExponent                 = %f\n",aHierarchicalMeshOut.mSettings.SimpExp);
    std::fprintf(stdout,"    aHierarchicalMesh.mElementList.Prjbeta                      = %f\n",aHierarchicalMeshOut.mSettings.Prjbeta);

    std::fprintf(stdout,"    aTensorMeshInputData.tVec                                 = ");
    for (uint i=0;i<5;++i) fprintf(stdout,"%d,",aHierarchicalMeshOut.mMeshData.SideSet( i ) );
    fprintf(stdout,"%d\n",aHierarchicalMeshOut.mMeshData.SideSet( 5 ) );

    std::fprintf(stdout,"    aTensorMeshInputData.remPt1                               = ");
    for (uint i=0;i<2;++i) fprintf(stdout,"%3.1f,",aHierarchicalMeshOut.mMeshData.DeleteDomain( 0, i ) );
    fprintf(stdout,"%3.1f\n",aHierarchicalMeshOut.mMeshData.DeleteDomain( 0, 2 ) );

    std::fprintf(stdout,"    aTensorMeshInputData.remPt2                               = ");
    for (uint i=0;i<2;++i) fprintf(stdout,"%3.1f,",aHierarchicalMeshOut.mMeshData.DeleteDomain( 1, i ) );
    fprintf(stdout,"%3.1f\n",aHierarchicalMeshOut.mMeshData.DeleteDomain( 1, 2 ) );
    if ( aHierarchicalMeshOut.mSettings.SDFAlpha.length() > 0 )
    {
        std::fprintf(stdout,"    aHierarchicalMesh.mSettings.SDFAlpha                      = ");
        fprintf(stdout,"%3.1f", aHierarchicalMeshOut.mSettings.SDFAlpha(0) );
        for( uint k=1; k<aHierarchicalMeshOut.mSettings.SDFAlpha.length(); ++k )
        {
            fprintf(stdout,",%3.1f", aHierarchicalMeshOut.mSettings.SDFAlpha(k) );
        }
        fprintf(stdout,"\n");
    }

    if ( aHierarchicalMeshOut.mSettings.SDFBeta.length() > 0 )
    {
        std::fprintf(stdout,"    aHierarchicalMesh.mSettings.SDFBeta                       = ");
        fprintf(stdout,"%3.1f", aHierarchicalMeshOut.mSettings.SDFBeta(0) );
        for( uint k=1; k<aHierarchicalMeshOut.mSettings.SDFBeta.length(); ++k )
        {
            fprintf(stdout,",%3.1f", aHierarchicalMeshOut.mSettings.SDFBeta(k) );
        }
        fprintf(stdout,"\n");
    }

    std::fprintf(stdout,"    aTensorMeshInputData.PointOfOrigin                        = ");
    for (uint i=0;i<2;++i) fprintf(stdout,"%3.1f,",aHierarchicalMeshOut.mSettings.PointOfOrigin( i ) );
    fprintf(stdout,"%3.1f\n",aHierarchicalMeshOut.mSettings.PointOfOrigin( 2 ) );
    if ( aHierarchicalMeshOut.mSettings.ProjectionFlag )
    {
        std::fprintf(stdout,"    Projection turned on\n");
    }
    else
    {
        std::fprintf(stdout,"    Projection turned off\n");
    }

    if (  aHierarchicalMeshOut.mSettings.PerformNormalization )
    {
        std::fprintf(stdout,"    Normalization turned on\n");
    }
    else
    {
        std::fprintf(stdout,"    Normalization turned off\n");
    }

    if (  aHierarchicalMeshOut.mSettings.FilterLevelWeight )
    {
        std::fprintf(stdout,"    Level weight for filter turned on\n");
    }
    else
    {
        std::fprintf(stdout,"    Level weight for filter turned off\n");
    }

    if (  aHierarchicalMeshOut.mSettings.AdaptiveBufferLayer )
    {
        std::fprintf(stdout,"    Adaptive buffer layer turned on\n");
    }
    else
    {
        std::fprintf(stdout,"    Adaptive buffer layer turned off\n");
    }

    if (  aHierarchicalMeshOut.mSettings.Staircasebuffer )
    {
        std::fprintf(stdout,"    Staircase buffer layer turned on\n");
    }
    else
    {
        std::fprintf(stdout,"    Staircase buffer layer turned off\n");
    }
    if (  aHierarchicalMeshOut.mSettings.TruncatedBsplines )
    {
        std::fprintf(stdout,"    Truncated B-splines turned on\n");
    }
    else
    {
        std::fprintf(stdout,"    Truncated B-splines turned off\n");
    }
    if (  aHierarchicalMeshOut.mSettings.LvlSetMethod )
    {
        std::fprintf(stdout,"    Criteria for level set field is turned on\n");
    }
    else
    {
        std::fprintf(stdout,"    Criteria for level set field is turned off\n");
    }

    if (  aHierarchicalMeshOut.mSettings.SDFField )
    {
        std::fprintf(stdout,"    SDF is turned on\n");
    }
    else
    {
        std::fprintf(stdout,"    SDF is turned off\n");
    }

    if (  aHierarchicalMeshOut.mSettings.DesignVariablesOnInitialMesh )
    {
        std::fprintf(stdout,"    Design variables active on initial mesh is turned on\n");
    }
    else
    {
        std::fprintf(stdout,"    Design variables active on initial mesh is turned off\n");
    }

    if (  aHierarchicalMeshOut.mSettings.RefineBaseOnExoFile )
    {
        std::fprintf(stdout,"    Read data from exodus file is turned on\n");
    }
    else
    {
        std::fprintf(stdout,"    Read data from exodus is turned off\n");
    }

    if (  aHierarchicalMeshOut.mSettings.InputOutputFilesForMorisAreAscii  )
    {
        std::fprintf(stdout,"    MORIS files are written as ASCII.\n");
    }
    else
    {
        std::fprintf(stdout,"    MORIS files are written as binary.\n");
    }

    if (  aHierarchicalMeshOut.mSettings.OutputFilesForFemdocAreBinary  )
    {
        std::fprintf(stdout,"    FEMDOC files are written as binary.\n");
    }
    else
    {
        std::fprintf(stdout,"    FEMDOC files are written as ASCII.\n");
    }
}

//------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::give_Tmatrix_and_IdField_Reorder(
        uint const & aElementId)
{
    if( mSettings.TruncatedBsplines == false)
    {
        mTMatrix.give_Tmatrix_and_IdField_Reorder(
                aElementId,mMeshData.ModelDim,mMeshData.Polynomial,
                mMeshData.NumberOfElementsPerDirection,mBasisData.BSplineActive,
                mElementData.TMatrix,mElementData.IdField,
                mElementData.TMatrixParentChildRelationFEM);
    }
    else
    {
        mTMatrix.give_Truncated_Tmatrix_and_IdField_Reorder(
                aElementId,mMeshData.ModelDim,mMeshData.Polynomial,
                mMeshData.NumberOfElementsPerDirection,mBasisData.BSplineActive,
                mElementData.TMatrix,mElementData.IdField,
                mElementData.TMatrixParentChildRelationFEM);
    }
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::give_Tmatrix_and_IdField(
        uint const & aElementId)
{
    if( mSettings.TruncatedBsplines == false)
    {
        mTMatrix.give_Tmatrix_and_IdField(
                aElementId,mMeshData.ModelDim,mMeshData.Polynomial,
                mMeshData.NumberOfElementsPerDirection,mBasisData.BSplineActive,
                mElementData.TMatrix,mElementData.IdField,mElementData.TMatrixParentChildRelationFEM);
    }
    else
    {
        mTMatrix.give_Truncated_Tmatrix_and_IdField(
                aElementId,mMeshData.ModelDim,mMeshData.Polynomial,
                mMeshData.NumberOfElementsPerDirection,mBasisData.BSplineActive,
                mElementData.TMatrix,mElementData.IdField,mElementData.TMatrixParentChildRelationFEM);
    }
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::give_Tmatrix_and_IdField_DesignVariables(
        uint const & aElementId)
{
    if( mSettings.TruncatedBsplines == false)
    {
        mTMatrix.give_Tmatrix_and_IdField_DesignToFEMProjection(
                aElementId,
                mMeshData.ModelDim,
                mMeshData.PolynomialDesign,
                mMeshData.NumberOfElementsPerDirection,
                mBasisData.DesignBSplineActive,
                mElementData.TMatrixDesign,
                mElementData.IdFieldDesign,
                mElementData.TMatrixParentChildRelationDesign);
    }
    else
    {
        mTMatrix.give_Truncated_Tmatrix_and_IdField_DesignToFEMProjection(
                aElementId,
                mMeshData.ModelDim,
                mMeshData.PolynomialDesign,
                mMeshData.NumberOfElementsPerDirection,
                mBasisData.DesignBSplineActive,
                mElementData.TMatrixDesign,
                mElementData.IdFieldDesign,
                mElementData.TMatrixParentChildRelationDesign);
    }
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::give_Tmatrix_and_IdField_DesignVariables_new(
        uint const & aElementId)
{
    if( mSettings.TruncatedBsplines == false)
    {
        mTMatrix.give_Tmatrix_and_IdField(
                aElementId,
                mMeshData.ModelDim,
                mMeshData.Polynomial,
                mMeshData.NumberOfElementsPerDirection,
                mBasisData.DesignBSplineActive,
                mElementData.TMatrixDesign,
                mElementData.IdFieldDesign,
                mElementData.TMatrixParentChildRelationFEM);
    }
    else
    {
        mTMatrix.give_Truncated_Tmatrix_and_IdField(
                aElementId,
                mMeshData.ModelDim,
                mMeshData.Polynomial,
                mMeshData.NumberOfElementsPerDirection,
                mBasisData.DesignBSplineActive,
                mElementData.TMatrixDesign,
                mElementData.IdFieldDesign,
                mElementData.TMatrixParentChildRelationFEM);
    }
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::give_Tmatrix_and_IdField_DesignVariables_LastStep(
        uint const & aElementId)
{
    if( mSettings.TruncatedBsplines == false)
    {
        mTMatrix.give_Tmatrix_and_IdField_DesignToFEMProjection(
                aElementId,
                mMeshData.ModelDim,
                mMeshData.PolynomialDesign,
                mMeshData.NumberOfElementsPerDirection,
                mBasisData.DesignBSplineActiveLastStep,
                mElementData.TMatrixDesign,
                mElementData.IdFieldDesign,
                mElementData.TMatrixParentChildRelationDesign);
    }
    else
    {
        mTMatrix.give_Truncated_Tmatrix_and_IdField_DesignToFEMProjection(
                aElementId,
                mMeshData.ModelDim,
                mMeshData.PolynomialDesign,
                mMeshData.NumberOfElementsPerDirection,
                mBasisData.DesignBSplineActiveLastStep,
                mElementData.TMatrixDesign,
                mElementData.IdFieldDesign,
                mElementData.TMatrixParentChildRelationDesign);
    }
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::give_Tmatrix_and_IdField_Specific_NaturalCoord(
        const uint      & aElementId,
        const Mat<real> & aNaturalCoordinate)
{
    if( mSettings.TruncatedBsplines == false)
    {
        mTMatrix.give_Tmatrix_and_IdField_SpecificPoint(
                aElementId,
                mMeshData.ModelDim,
                mMeshData.Polynomial,
                mMeshData.NumberOfElementsPerDirection,
                mBasisData.BSplineActive,
                mElementData.TMatrix,
                mElementData.IdField,
                aNaturalCoordinate,
                mElementData.TMatrixParentChildRelationFEM);
    }
    else
    {
        mTMatrix.give_Truncated_Tmatrix_and_IdField_SpecificPoint(
                aElementId,
                mMeshData.ModelDim,
                mMeshData.Polynomial,
                mMeshData.NumberOfElementsPerDirection,
                mBasisData.BSplineActive,
                mElementData.TMatrix,
                mElementData.IdField,
                aNaturalCoordinate,
                mElementData.TMatrixParentChildRelationFEM);
    }
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::give_Tmatrix_and_IdField_Specific_PhysicalCoord(
        Mat<real> const & aPhysicalCoordinate)
{
    uint tElementId = this->give_active_element_for_coordinate(aPhysicalCoordinate); // Determine from the physical coordinate the element ID
    uint tElementLevel = this->give_element_level(tElementId); // Determine the level of the element
    Mat<uint> tBasisOfElement = mHMRElement.give_basis_of_element(tElementId,mMeshData.ModelDim,1,mMeshData.NumberOfElementsPerDirection); // Determine the basis functions of the element in a linear approach
    Mat<real> tElementCoord1 = this->give_coordinate_from_basis(tBasisOfElement(0)); // Take the first basis function and calculate the coordinates (Corner coordinate of the element bottom left)
    Mat<real> tNaturalCoordinate(1,mMeshData.ModelDim,0);
    for(uint i = 0; i < mMeshData.ModelDim; i++)
    {
        // Determine the natural coordinate within the element (0 <= xi_i <= 1)
        tNaturalCoordinate(i) = (aPhysicalCoordinate(i)-tElementCoord1(i))/(mMeshData.ElementLength(i)*pow(0.5,tElementLevel));
    }
    if( mSettings.TruncatedBsplines == false)
    {
        mTMatrix.Hierarchical_Mesh_TMatrix::give_Tmatrix_and_IdField_SpecificPoint(
                tElementId,
                mMeshData.ModelDim,
                mMeshData.Polynomial,
                mMeshData.NumberOfElementsPerDirection,
                mBasisData.BSplineActive,
                mElementData.TMatrix,
                mElementData.IdField,
                tNaturalCoordinate,
                mElementData.TMatrixParentChildRelationFEM);
    }
    else
    {
        mTMatrix.Hierarchical_Mesh_TMatrix::give_Truncated_Tmatrix_and_IdField_SpecificPoint(
                tElementId,
                mMeshData.ModelDim,
                mMeshData.Polynomial,
                mMeshData.NumberOfElementsPerDirection,
                mBasisData.BSplineActive,
                mElementData.TMatrix,
                mElementData.IdField,
                tNaturalCoordinate,
                mElementData.TMatrixParentChildRelationFEM);
    }
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Main::give_IdField(
        uint const & aElementId)
{
    Mat<uint> tIdField;

    mTMatrix.give_IdField(
            aElementId,
            mMeshData.ModelDim,
            mMeshData.Polynomial,
            mMeshData.NumberOfElementsPerDirection,
            mBasisData.BSplineActive,
            tIdField);

    return tIdField;
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Main::give_Id_DesignField(
        uint const & aElementId)
{
    Mat<uint> tIdField;

    mTMatrix.give_IdField(
            aElementId,
            mMeshData.ModelDim,
            mMeshData.Polynomial,
            mMeshData.NumberOfElementsPerDirection,
            mBasisData.DesignBSplineActive,
            tIdField);

    return tIdField;
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::hierarchical_element_refinement()
{
    MORIS_ASSERT( mElementData.DeactivateElement.n_cols() <= 1, " mElementData.DeactivateElement must be a row vector");

    // Initialize level
    uint tLevel = 0;

    // the timer for the refinement procedure
    moris::tic  elementrefinement;
    if ( mElementData.DeactivateElement.length() > 0 )
    {

        if ( mSettings.UseSymmetry )
        {
            make_deactivate_list_symmetric();
        }

        // if LevelLastStep > 0, check Element active/passive bitsets, and resize if necessary
        uint tNumberElements;
        uint tMaxElement = mElementData.DeactivateElement.max();

        tLevel = give_element_level(tMaxElement);

        if( mMeshData.LevelLastStep > 0 )
        {
            // calculate the maximum possible list of elements
            uint tNumberElementsA = mBaseElement.give_number_of_elements(mMeshData.LevelLastStep,
                                                                   mMeshData.ModelDim,
                                                                   mMeshData.NumberOfElementsPerDirection);

            uint tNumberElementsB =  give_number_of_lagrange_basis(tLevel);

            // get maximum
            tNumberElements = (tNumberElementsA > tNumberElementsB) ? (tNumberElementsA) : (tNumberElementsB);
        }
        else
        {
            tNumberElements = give_number_of_lagrange_basis(tLevel);
        }

        // check size of active bitset
        if (mElementData.ElementActive.size() < tNumberElements)
        {
            mElementData.ElementActive.resize(tNumberElements);
        }
        // check size of passive bitset
        if (mElementData.ElementRefined.size() < tNumberElements)
        {
            mElementData.ElementRefined.resize(tNumberElements);
        }

        // Deactivate elements to broadcast this message
        for( uint i=0; i< mElementData.DeactivateElement.length(); i++ )
        {
            mElementData.ElementActive.reset( mElementData.DeactivateElement( i ) ); // Deactivate elements

        }

        mMPI.broadcast_bitset_logical_and( mElementData.ElementActive );

        mRefinement.hierarchical_element_refinement(mMeshData.ModelDim,
                                                    mMeshData.Polynomial,
                                                    mMeshData.NumberOfElementsPerDirection,
                                                    mElementData.DeactivateElement,
                                                    mElementData.ElementActive,
                                                    mElementData.ElementRefined);

        tLevel = mBaseElement.give_element_level(tMaxElement,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection)+1; // Give current level
    }

    mMPI.gather_value_and_bcast_max(tLevel);
    if ((mMeshData.Level)< tLevel)
    {
        mMeshData.Level = tLevel;
        mMeshData.NumberElements = mBaseElement.give_number_of_elements(
                tLevel,
                mMeshData.ModelDim,
                mMeshData.NumberOfElementsPerDirection);

        mBasisData.NumberBasis = mBasis.give_number_of_basis(
                tLevel,
                mMeshData.ModelDim,
                mMeshData.Polynomial,
                mMeshData.NumberOfElementsPerDirection);

        // Increase size of active and passive bitset
        if( mElementData.ElementActive.size() < mMeshData.NumberElements )
        {
            mElementData.ElementActive.resize(mMeshData.NumberElements);
            mElementData.ElementRefined.resize(mMeshData.NumberElements);
        }
        // Increase size of active and passive bitset
        if (mBasisData.BSplineActive.size() < mBasisData.NumberBasis )
        {
            mBasisData.BSplineActive.resize( mBasisData.NumberBasis );
            mBasisData.BSplineRefined.resize( mBasisData.NumberBasis );
        }
        this->update_passive_outer_layer();
    }

    // this line is correct!
    mMPI.broadcast_bitset_logical_or( mElementData.ElementActive );
    mMPI.broadcast_bitset_logical_or( mElementData.ElementRefined );

    // An element can not be both refined and active.
    // if the refine flag is set, the element is deactivated
    for( moris::uint k=0; k<mElementData.ElementRefined.size(); ++k )
    {
        if( mElementData.ElementRefined.test(k) )
        {
            mElementData.ElementActive.reset(k);
        }
    }

    real tElapsedelementrefinement = elementrefinement.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for de/activating elements : %f [sec]\n",tElapsedelementrefinement);

    mElementData.DeactivateElement.set_size( 0, 0 );

    // update internal element information
    mProcData.ElementListOnProc = give_element_on_proc();
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::activate_basisfunction()
{
    moris::tic  activatebasis;

    mProcData.ElementListOnProc = mRefinement.give_active_elements(
            mMeshData.ModelDim,
            mMeshData.Polynomial,
            mMeshData.NumberOfElementsPerDirection,
            mMeshData.Level,
            mProcData.ElementListOnProcInit,
            mElementData.ElementActive);

    mRefinement.activate_basisfunction(
            mMeshData.ModelDim,
            mMeshData.Polynomial,
            mMeshData.NumberOfElementsPerDirection,
            mBasisData.BSplineActive,
            mElementData.ElementActive,
            mProcData.ElementListOnProc,
            mBasisData.BSplineActiveList);

    mBasisData.NumberBasis = mBasis.give_number_of_basis(
            mMeshData.Level,
            mMeshData.ModelDim,
            mMeshData.Polynomial,
            mMeshData.NumberOfElementsPerDirection);

    if (mBasisData.BSplineActive.size() < mBasisData.NumberBasis )
    {
        mBasisData.BSplineActive.resize( mBasisData.NumberBasis );
        mBasisData.BSplineRefined.resize( mBasisData.NumberBasis );
    }

    mMPI.broadcast_bitset_logical_or(mBasisData.BSplineActive);

    real tElapsedactivatebasis = activatebasis.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for activating basis : %f [sec]\n",tElapsedactivatebasis);
}

//--------------------------------------------------------------------------------

// FIXME: experimental function, turned out to be buggy.
void
Hierarchical_Mesh_Main::activate_basisfunction_new()
{
    MORIS_ASSERT( false, "Please do not use activate_basisfunction_new. This funciton is not roubust yet.");
    moris::tic  activatebasis;
    mProcData.ElementListOnProc = give_active_elements();
    mProcData.ElementRefinedListOnProc = give_passive_elements();

    mRefinement.activate_basisfunction_new(
            mMeshData.ModelDim,
            mMeshData.Polynomial,
            mMeshData.NumberOfElementsPerDirection,
            mMeshData.Level,
            mBasisData.BSplineActive,
            mBasisData.BSplineRefined,
            mElementData.ElementActive,
            mElementData.ElementRefined,
            mProcData.ElementListOnProc,
            mProcData.ElementRefinedListOnProc,
            mBasisData.BSplineActiveList);

    mBasisData.NumberBasis = mBasis.give_number_of_basis(
            mMeshData.Level,
            mMeshData.ModelDim,
            mMeshData.Polynomial,
            mMeshData.NumberOfElementsPerDirection);

    mMPI.broadcast_bitset_logical_or(mBasisData.BSplineActive);
    real tElapsedactivatebasis = activatebasis.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for activating basis : %f [sec]\n",tElapsedactivatebasis);
}

//--------------------------------------------------------------------------------
void
Hierarchical_Mesh_Main::activate_design_lagrange_basis()
{
    // create timer object
    moris::tic tTimer;

    uint tNumberOfBasis = give_number_of_lagrange_basis( mMeshData.Level );
    // set size ofLagrange basis
    mBasisData.DesignLagrangeBasisActive.resize(tNumberOfBasis);

    // set all values in basis to false
    mBasisData.DesignLagrangeBasisActive.reset();

    // loop over all elements on proc
    for ( uint e=0; e< mProcData.ElementListOnProc.length(); ++e )
    {
        uint tThisElement = mProcData.ElementListOnProc ( e ) ;

        // get all basis of this element
        Mat<uint> tLagrangeBasisOfElement = give_lagrange_basis_of_element( tThisElement );

        // loop over all basis of this element
        for ( uint b=0; b< tLagrangeBasisOfElement.length(); ++b )
        {
            uint tThisBasis = tLagrangeBasisOfElement( b );

            uint tLevelOfThisBasis = give_basis_level ( tThisBasis );

            // get position of this node
            Mat<uint> tIJK = give_position_of_basis( tThisBasis );

            bool tParentExists = true;

            // look up
            while ( tParentExists && tLevelOfThisBasis > 0 )
            {
                for (uint i=0; i< mMeshData.ModelDim; ++i )
                {
                    // check if position is divisable by two
                    tParentExists = tParentExists && ( tIJK( i ) % 2 == 0 ) ;

                    // break the loop if we have no parent
                    if (! tParentExists)
                        break;
                }

                // use basis from other level if basis exists
                if ( tParentExists )
                {
                    // take position from basis on upper level
                    for (uint i=0; i< mMeshData.ModelDim; ++i)
                    {
                        tIJK( i ) /= 2;
                    }

                    // decrement level
                    --tLevelOfThisBasis;
                }
            }

            // get new basis number
            tThisBasis = give_basis_of_position( tLevelOfThisBasis, tIJK );
            // activate this basis ( including hanging nodes )
            mBasisData.DesignLagrangeBasisActive.set( tThisBasis );
        }
    }

    // boradcast logical or
    mMPI.broadcast_bitset_logical_or(mBasisData.DesignLagrangeBasisActive);

    mBasisData.NumberOfActiveDesignLagrangeBasis = mBasisData.DesignLagrangeBasisActive.count();

    // stop timer
    real tElapsedTime = tTimer.toc<moris::chronos::seconds>().wall;
    // print timing information
    if( mSettings.TimingTest == true )
            std::fprintf(stdout,"Time for activating Lagrange basis : %f [sec]\n",tElapsedTime);
}
//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::activate_basisdesignfunction()
{
    moris::tic  activatedesignbasis;
    mElementData.ElementListActiveDesign = mRefinement.give_active_elements(
            mMeshData.ModelDim,
            mMeshData.Polynomial,
            mMeshData.NumberOfElementsPerDirection,
            mMeshData.Level,
            mProcData.ElementListOnProcInit,
            mElementData.ElementActive);

    mElementData.ElementActiveDesign = mElementData.ElementActive;

    mRefinement.activate_basisfunction(
            mMeshData.ModelDim,
            mMeshData.PolynomialDesign,
            mMeshData.NumberOfElementsPerDirection,
            mBasisData.DesignBSplineActive,
            mElementData.ElementActive,
            mElementData.ElementListActiveDesign,
            mBasisData.DesignBSplineActiveList);

    mMPI.broadcast_bitset_logical_or(mBasisData.DesignBSplineActive);
    real tElapsedactivatedesignbasis = activatedesignbasis.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for activating design basis : %f [sec]\n",tElapsedactivatedesignbasis);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::refinement_criteria_nodal_field(
        map<uint, real> & aNodalField,
        const uint & aMaxRefinementLevel,
        const real & aLowerBound,
        const real & aUpperBound)

{
    std::printf("Nodal field refinement: %i elements initially to be refined\n",(uint)mElementData.DeactivateElement.length());

    //Check for refinement if max level of refinement is higher then 0
    if( aMaxRefinementLevel > 0 )
    {
        Mat<uint> tDeactivateElement;
        mRefinement.refinement_criteria_nodal_field(
                mMeshData.ModelDim,
                mMeshData.PolynomialDesign,
                mMeshData.NumberOfElementsPerDirection,
                mProcData.ElementListOnProcLastStepDesign,
                tDeactivateElement,
                mElementData.ElementActive,
                mBasisData.DesignBSplineActiveLastStep,
                aNodalField,
                mSettings.TruncatedBsplines,
                aLowerBound,
                aUpperBound,
                aMaxRefinementLevel,
                mSettings.BufferElements,
                mSettings.AdaptiveBufferLayer,
                mSettings.Staircasebuffer);
        
        std::printf("Nodal field refinement: %i elements flagged for refinement\n",(uint)tDeactivateElement.length());

        if( isempty(mElementData.DeactivateElement)) // If nothing is in deactivate elements, it is equal, otherwise resize the list and add the new entries
        {
            mElementData.DeactivateElement = tDeactivateElement;
        }
        else if ( tDeactivateElement.length() > 0 )
        {
            uint tDeactivateElementlength = mElementData.DeactivateElement.length();
            mElementData.DeactivateElement.resize(tDeactivateElementlength+tDeactivateElement.length(),1);
            mElementData.DeactivateElement.rows(tDeactivateElementlength,tDeactivateElementlength+tDeactivateElement.length()-1) = tDeactivateElement.rows(0,tDeactivateElement.length()-1);
            mElementData.DeactivateElement = unique(mElementData.DeactivateElement);
        }
        else
        {
            std::printf("Nodal field refinement: no additional elements found for refinement\n");
        }
    }
    else
    {
        std::printf("Nodal field refinement: no elements flagged for refinement as aMaxRefinementLevel = 0\n");
    }

    std::printf("Nodal field refinement: %i elements to be refined\n",(uint)mElementData.DeactivateElement.length());
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::refinement_criteria_element_field(
        map<uint, real> & aElementField,
        map<uint, real> & aNodalField,
        const uint & aMaxRefinementLevel,
        const real & aLowerBound,
        const real & aUpperBound)
{
    std::printf("Elemental field refinement: Max ref level %u, lower bound %f, upper bound %f\n",  aMaxRefinementLevel, aLowerBound, aUpperBound);
    std::printf("Elemental field refinement: %i elements initially to be refined\n",(uint)mElementData.DeactivateElement.length());

    //Check for refinement if max level of refinement is higher then 0
    if( aMaxRefinementLevel > 0 )
    {
        Mat<uint> tDeactivateElement;
        mRefinement.refinement_criteria_element_field(
                mMeshData.ModelDim,
                mMeshData.PolynomialDesign,
                mMeshData.NumberOfElementsPerDirection,
                mProcData.ElementListOnProcLastStepDesign,
                tDeactivateElement,
                aElementField,
                mElementData.ElementActive,
                mBasisData.DesignBSplineActiveLastStep,
                aNodalField,
                mSettings.TruncatedBsplines,
                aLowerBound,
                mSettings.ElementIncrementBound,
                aUpperBound,
                aMaxRefinementLevel,
                mSettings.BufferElements,
                mSettings.AdaptiveBufferLayer,
                mSettings.Staircasebuffer);

        std::printf("Elemental field refinement: %i elements flagged for refinement\n",(uint)tDeactivateElement.length());

        //Check the size of the element active bitset. IF this is changed, change the passive to the same size
        if( mElementData.ElementActive.size() > mElementData.ElementRefined.size() )
        {
            mElementData.ElementRefined.resize( mElementData.ElementActive.size() );
        }
        
        // If nothing is in deactivate elements, it is equal, otherwise resize the list and add the new entries
        if( isempty(mElementData.DeactivateElement)) 
        {
            mElementData.DeactivateElement = tDeactivateElement;
        }
        else if ( tDeactivateElement.length() > 0 )
        {
            uint tDeactivateElementlength = mElementData.DeactivateElement.length();

            mElementData.DeactivateElement.resize(tDeactivateElementlength+tDeactivateElement.length(),1);
            mElementData.DeactivateElement.rows(tDeactivateElementlength,tDeactivateElementlength+tDeactivateElement.length()-1) = tDeactivateElement.rows(0,tDeactivateElement.length()-1);
            mElementData.DeactivateElement = unique(mElementData.DeactivateElement);
        }
        else
        {            
            std::printf("Elemental field refinement: no additional elements found for refinement\n");
        }
    }
    else
    {
        std::printf("Elemental field refinement: no elements flagged for refinement as aMaxRefinementLevel = 0\n");
    }

    std::printf("Elemental field refinement: %i elements to be refined\n",(uint)mElementData.DeactivateElement.length());
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::find_elements_of_object(
        uint & aMaxLevelOfRefinement,
        const Mat<uint> & aCandidateElementsToDeactivate)
{
    std::printf("SDF refinement refinement: %i elements initially to be refined\n",(uint)mElementData.DeactivateElement.length());

    //Check for refinement if max level of refinement is higher then 0
    if( aMaxLevelOfRefinement > 0 )
    {
        //In the first step, the elements of the last step does not exist, so the initial active elements are the elements from the last step
        if( mElementData.ElementActiveLastStep.count() == 0 )
        {
            mElementData.ElementActiveLastStep = mElementData.ElementActive;
        }
        Mat<uint> tElementsToDeactivate;

        if( aCandidateElementsToDeactivate.length() > 0)
        {
            mRefinement.find_elements_of_object(
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    aMaxLevelOfRefinement,
                    mMeshData.NumberOfElementsPerDirection,
                    mElementData.ElementActiveLastStep,
                    mSettings.BufferElements,
                    mMeshData.Level,
                    mSettings.AdaptiveBufferLayer,
                    mSettings.Staircasebuffer,
                    aCandidateElementsToDeactivate,
                    tElementsToDeactivate);

            std::printf("SDF refinement: %i elements flagged for refinement\n",(uint)tElementsToDeactivate.length());

            // If nothing is in deactivate elements, it is equal, otherwise resize the list and add the new entries
            if (isempty(mElementData.DeactivateElement))
            {
                mElementData.DeactivateElement = tElementsToDeactivate;
            } else
            {
                uint tDeactivateElementlength = mElementData.DeactivateElement.length();

                mElementData.DeactivateElement.resize(tDeactivateElementlength + tElementsToDeactivate.length(),
                                                             1);
                mElementData.DeactivateElement.rows(tDeactivateElementlength,
                                                           tDeactivateElementlength + tElementsToDeactivate.length() -
                                                           1) = tElementsToDeactivate.rows(0,
                                                                                           tElementsToDeactivate.length() -
                                                                                           1);
                mElementData.DeactivateElement = unique(mElementData.DeactivateElement);
            }
        }
    }
    else
    {
        std::printf("SDF refinement refinement: no elements flagged for refinement as MaxDesignLevelOfRefinement = 0\n");
    }

    std::printf("SDF refinement refinement: %i elements to be refined\n",(uint)mElementData.DeactivateElement.length());
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::find_elements_in_stencil()
{
    /*if( mSettings.Refinement )
    {
        Mat<uint> tElementListOnProc = mRefinement.give_active_elements(mMeshData.ModelDim,mMeshData.Polynomial,mMeshData.NumberOfElementsPerDirection,mMeshData.Level,mProcData.ElementListOnProcInit,mElementData.ElementActive);
        mElementData.DeactivateElement = mRefinement.find_elements_in_stencil(
                mMeshData.ModelDim,
                mMeshData.PolynomialDesign,
                mMeshData.NumberOfElementsPerDirection,
                mElementData.ElementActive,
                tElementListOnProc,
                mFieldData.ElementDensity,
                mSettings.FeatureLowerBound,
                mSettings.FeatureResolution);
    }*/

    // call stencil only if feature size resolution is set
    if ( mSettings.FeatureResolution > 0 )
    {
        // perform search for elements in stencil
        Mat<uint> tElementsInStencil = mRefinement.find_elements_in_stencil(
                        mMeshData.ModelDim,
                        mMeshData.PolynomialDesign,
                        mMeshData.NumberOfElementsPerDirection,
                        mElementData.ElementActiveDesignLastStep,
                        mProcData.ElementListOnProcLastStepDesign,
                        mFieldData.ElementDensity,
                        mSettings.FeatureLowerBound,
                        mSettings.FeatureResolution,
                        mSettings.BufferElements,
                        mSettings.AdaptiveBufferLayer,
                        mSettings.Staircasebuffer );

        // count number of elements
        uint tNumberOfElementsInStencil = tElementsInStencil.length();

        // if there are some, add them to refinement list
        if( tNumberOfElementsInStencil > 0 )
        {

            std::printf("Found %u elements in stencil and add them to refinement list\n", tNumberOfElementsInStencil );

            uint tNumberOfElements = mElementData.DeactivateElement.length();

            if ( tNumberOfElements == 0 )
            {
                mElementData.DeactivateElement = tElementsInStencil ;
            }
            else
            {
                mElementData.DeactivateElement.resize( tNumberOfElements + tNumberOfElementsInStencil, 1);
                mElementData.DeactivateElement.rows(
                        tNumberOfElements,
                        tNumberOfElements + tNumberOfElementsInStencil - 1) = tElementsInStencil.rows(0, tNumberOfElementsInStencil - 1);
            }
        }
        else
        {
            std::printf("No elements found in stencil. \n");
        }
    }
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::update_basis_values_on_design_elements(
        map<uint, real> & aNodalField)
{
    moris::tic updatevalues;

    mRefinement.update_basis_values_on_design_elements(
            mMeshData.ModelDim,
            mMeshData.PolynomialDesign,
            mMeshData.NumberOfElementsPerDirection,
            mProcData.ElementListOnProcLastStepDesign,
            mSettings.TruncatedBsplines,
            mBasisData.DesignBSplineActiveLastStep,
            aNodalField,
            mElementData.TMatrixParentChildRelationDesign);

    real tElapsedupdatevalues = updatevalues.toc<moris::chronos::seconds>().wall;

    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for updating values on design elements: %f [sec]\n",tElapsedupdatevalues);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::set_nodeset(
        Mat<uint> & aNodeSet)
{
    mMeshData.NodeSet = mSets.set_nodeset(
            mMeshData.ModelDim,
            mMeshData.Polynomial,
            mMeshData.NumberOfElementsPerDirection,
            aNodeSet,mProcData.ElementListOnProc);

    Mat<uint> SetCol = (mMeshData.NodeSet ).col(2);
    mMeshData.NumSets(0) = SetCol.max();
    mMPI.gather_value_and_bcast_max(mMeshData.NumSets(0));
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::set_sideset(
        Mat<uint> & aSideSet)
{
    mMeshData.NumSets(1) = aSideSet.max();
    mMeshData.SideSet = mSets.set_sideset(
            mMeshData.ModelDim,
            mMeshData.Polynomial,
            mMeshData.NumberOfElementsPerDirection,
            mProcData.DecompElement,
            aSideSet);
    mMPI.gather_value_and_bcast_max(mMeshData.NumSets(1));
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::create_MTK_file(
        std::string const & aOutputFileName)
{
    mMeshData.SideSetMTK = mSets.update_sideset(mMeshData.ModelDim,
                                                       mMeshData.NumberOfElementsPerDirection,
                                                       mMeshData.Level,
                                                       mMeshData.NumSets,
                                                       mMeshData.SideSet,
                                                       mElementData.ElementActive);
    mMTK.create_MTK_file(aOutputFileName,
                         mMeshData.ModelDim,
                         mMeshData.NumberOfElementsPerDirection,
                         mElementData.FeTopoMTK,
                         mElementData.FeCoord,
                         mBasisData.LagrangeNodeProcOwner,
                         mProcData.ElementListOnProc,
                         //mBasisData.NodalLocaltoGlobal,
                         mBasisData.LagrangeToBSplineMap,
                         mProcData.LagrangeListOnProc,
                         mMeshData.SideSetMTK);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::create_MTK_file(
        std::string const & aOutputFileName,
        Cell<Mat<real>> & aFieldData,
        Cell<std::string> & aFieldName,
        Cell<enum EntityRank> & aFieldRank)
{

    mMeshData.SideSetMTK = mSets.update_sideset(mMeshData.ModelDim,
                                                       mMeshData.NumberOfElementsPerDirection,
                                                       mMeshData.Level,
                                                       mMeshData.NumSets,
                                                       mMeshData.SideSet,
                                                       mElementData.ElementActive);
    mMTK.create_MTK_file(aOutputFileName,
                         mMeshData.ModelDim,
                         mMeshData.NumberOfElementsPerDirection,
                         mElementData.FeTopoMTK,
                         mElementData.FeCoord,
                         mBasisData.LagrangeNodeProcOwner,
                         mProcData.ElementListOnProc,
                         //mBasisData.NodalLocaltoGlobal,
                         mBasisData.LagrangeToBSplineMap,
                         mProcData.LagrangeListOnProc,
                         mMeshData.SideSetMTK,
                         aFieldData,
                         aFieldName,
                         aFieldRank);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::create_MTK_file_last_step(
        std::string const & aOutputFileName)
{
    moris::tic  createoldmesh;

    //Copy data to a dummy and copy data from last steps
    mElementData.ElementActiveDummy    = mElementData.ElementActive;
    mElementData.ElementActive         = mElementData.ElementActiveLastStep;
    mBasisData.BSplineActiveDummy      = mBasisData.BSplineActive;
    mBasisData.BSplineActive           = mBasisData.BSplineActiveLastStep;
    mElementData.ElementActiveDesign   = mElementData.ElementActiveDesignLastStep;
    mBasisData.DesignBSplineActive       = mBasisData.DesignBSplineActiveLastStep;

    // Save the filter radius, it's not needed to create the filter for the solution of the last step
    real tFilterRadiusDummy                 = mSettings.FilterRadius;
    mSettings.FilterRadius = 0.0;

    // Save the old level to restore it after the creation of the mesh
    uint tLevelDummy = mMeshData.Level;
    mMeshData.Level = mMeshData.LevelLastStep;

    //Create mesh data with data from last step
    this->create_mesh_data_new();

    //If we read from an exo file, ADV and PDV have no sense and they are declared with zeros
    /*if( mSettings.RefineBaseOnExoFile == true )
    {
        mFieldData.NodalADV.set_size( mFieldData.NodalLevelSet.length(), 1, 0 );
        mFieldData.NodalPDV.set_size( mFieldData.NodalLevelSet.length(), 1, 0 );
    }*/

    //Create output files
    Mat<real> NodLvlSetField_output((mElementData.FeCoord).n_rows(),1,0);
    Mat<real> NodADVField_output((mElementData.FeCoord).n_rows(),1,0);
    Mat<real> NodPDVField_output((mElementData.FeCoord).n_rows(),1,0);

    for(uint i  = 0; i < mBasisData.IdFieldDesignNodalField.n_rows(); i++)
    {
        for(uint j = 0; j < mBasisData.IdFieldDesignNodalField(i,0); j++)
        {
            NodLvlSetField_output(i) += mFieldData.NodalLevelSet.find(mBasisData.IdFieldDesignNodalField(i,j+1)) * mBasisData.TMatrixDesignNodalField(i,j+1);
            NodADVField_output(i)    += mFieldData.NodalADV.find(mBasisData.IdFieldDesignNodalField(i,j+1))    * mBasisData.TMatrixDesignNodalField(i,j+1);
            NodPDVField_output(i)    += mFieldData.NodalPDV.find(mBasisData.IdFieldDesignNodalField(i,j+1))    * mBasisData.TMatrixDesignNodalField(i,j+1);
        }
    }

    // create maps needed for reorder and MTK output
    if (mBasisData.LagrangeToBSplineMap.length() != NodLvlSetField_output.length())
    {
        create_maps();
    }

    Mat<real> ElemSetField_output((mProcData.ElementListOnProc).length(),1,UINT_MAX);
    Mat<real> ElemDenistyField_output((mProcData.ElementListOnProc).length(),1,UINT_MAX);

    for(uint i=0; i<(mProcData.ElementListOnProc).length(); i++)
    {
        ElemSetField_output(i) = mBaseElement.give_element_level(mProcData.ElementListOnProc(i),mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
        ElemDenistyField_output(i) = mFieldData.ElementDensity.find(mProcData.ElementListOnProc(i));
    }

    //Create Cells with field data for MTK
    Cell < Mat< real > > aFieldData      = { NodLvlSetField_output, ElemSetField_output, NodADVField_output, NodPDVField_output, ElemDenistyField_output};
    Cell < std::string > aFieldName      = { "NodalLevelSet", "ElementLevel", " Nodal_ADV ", " Nodal_PDV ", "ElementDensity"};
    Cell < enum EntityRank > aFieldRanks = { EntityRank::NODE, EntityRank::ELEMENT,  EntityRank::NODE,  EntityRank::NODE, EntityRank::ELEMENT};

    // FIXME make this output print first density
    // add SDF density to output, if SDF was activated
    /* if ( mSettings.SDFField )
    {
        // reserve memory four output data
        Mat<real>  tSDFDensity_output ( mProcData.ElementListOnProc.length(), 1);

        // loop over all elements and read data from map
        for( uint k=0; k<mProcData.ElementListOnProc.length(); ++k )
        {
            tSDFDensity_output( k ) = mFieldData.ElementVolumeDensitySDF.find( mProcData.ElementListOnProc( k ) );
        }

        // add vector to field data list
        aFieldData.push_back( tSDFDensity_output );

        // add name
        aFieldName.push_back( "SDF_Volume_Density" );

        // add rank
        aFieldRanks.push_back( EntityRank::ELEMENT );
    } */

    //Create the output with MTK
    this->create_MTK_file(aOutputFileName,aFieldData,aFieldName,aFieldRanks);

    //Copy data back from dummy
    mSettings.FilterRadius  = tFilterRadiusDummy; // Save back the solution of the filter radius from the dummy variable
    mElementData.ElementActive = mElementData.ElementActiveDummy;
    mBasisData.BSplineActive   = mBasisData.BSplineActiveDummy;
    mMeshData.Level         = tLevelDummy;

    real tElapsedcreateoldmesh = createoldmesh.toc<moris::chronos::seconds>().wall;

    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for creating file for last step : %f [sec]\n",tElapsedcreateoldmesh);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::create_MTK_file_with_adv(
        std::string const & aOutputFileName)
{
    moris::tic  createnewmesh;

    // rearaange array
    Cell < Mat< real > > aFieldData      = { mFieldData.LagrangeADVField };
    Cell < std::string > aFieldName      = { "NodalADV"};
    Cell < enum EntityRank > aFieldRanks = { EntityRank::NODE};
    this->create_MTK_file(aOutputFileName,aFieldData,aFieldName,aFieldRanks);
    real tElapsedcreatenewmesh = createnewmesh.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for creating file for new step : %f [sec]\n",tElapsedcreatenewmesh);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::create_MTK_file_with_adv_and_obj(
        std::string const & aOutputFileName )
{
    moris::tic  createnewmesh;
    if (mBasisData.LagrangeToBSplineMap.length() != mProcData.LagrangeListOnProc.length())
    {
        create_maps();
    }

    Cell < Mat< real > > aFieldData ;
    Cell < std::string > aFieldName ;
    Cell < enum EntityRank > aFieldRanks ;

    // can only write Nodal Abstract Design Variable field, if L2 projection was performed
    if ( mSettings.ProjectionFlag )
    {
        aFieldData.push_back ( mFieldData.LagrangeADVField );
        aFieldName.push_back ( "NodalADV" );
        aFieldRanks.push_back( EntityRank::NODE );
    }

    //Create sdf field data from computed sdf field
    if( mFieldData.ObjectSDF.size() > 0 )
    {
        uint tNumberOfFields;

        tNumberOfFields = 2*mFieldData.ObjectSDF.size();

        Cell<Mat<real> > aSDFFieldData( mFieldData.ObjectSDF.size() );
        Cell<std::string> aSDFFieldName( tNumberOfFields );
        Cell<enum EntityRank> aSDFFieldRanks( tNumberOfFields );
        Cell< Mat<real> > tObjectSDF;
        //Create SDF field, based on data from SDF module
        //Mat<real> NodalSDFField_output;

        for (uint j = 0; j < mFieldData.ObjectSDF.size(); j++)
        {
            aSDFFieldRanks(j) = EntityRank::NODE;
            aSDFFieldName(j)  = "NP_PROP_" + std::to_string(j + 1);
            tObjectSDF.push_back( mFieldData.ObjectSDF( j ) ) ;

        }
        aFieldData.append(tObjectSDF);
        uint tOffset = mFieldData.ObjectSDF.size();
        Mat<real> tSDFFlags;
        for( uint j = 0 ; j < mFieldData.ObjectSDF.size(); j++ )
        {
            tSDFFlags.set_size( mFieldData.ObjectSDF( j ).length(), 1, 0 );
            for ( uint i = 0; i < mFieldData.ObjectSDF( j ).length(); i++ )
            {
                if( mFieldData.ObjectFlagSDF( j ).test( i ) )
                {
                    tSDFFlags( i ) = 1;
                }
            }
            aSDFFieldData ( j )           =  tSDFFlags ;
            aSDFFieldRanks( j + tOffset)  = EntityRank::NODE;
            aSDFFieldName ( j + tOffset ) = "HAS_NP_PROP_" + std::to_string( j + 1 );
        }

        aFieldData.append(aSDFFieldData);
        aFieldName.append(aSDFFieldName);
        aFieldRanks.append(aSDFFieldRanks);

    }
    //Save data in MTK file
    this->create_MTK_file( aOutputFileName, aFieldData, aFieldName, aFieldRanks );
    real tElapsedcreatenewmesh = createnewmesh.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for creating file for new step : %f [sec]\n",tElapsedcreatenewmesh);
}
//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_element_dof_connectivity_ascii()
{
    moris::tic  save_dofconnectivity;
    mOutput.output_element_dof_connectivity_ascii(
            mMeshData.ModelDim,
            mMeshData.Polynomial,
            mMeshData.NumberOfElementsPerDirection,
            mElementData.FeTopoMTK,
            mProcData.ElementListOnProc,
            mBasisData.BSplineActive,
            mSettings.Refinement,
            mSettings.TruncatedBsplines,
            mElementData.TMatrixParentChildRelationFEM,
            mBasisData.BsplineBasisList,
            mBasisData.BSplineMap );

    real tElapsedsave_dofconnectivity = save_dofconnectivity.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving elemental dof connectivity : %f [sec]\n",tElapsedsave_dofconnectivity);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_element_dof_connectivity_binary()
{
    moris::tic  save_dofconnectivity;

    mOutput.output_element_dof_connectivity_binary(
            mMeshData.ModelDim,
            mMeshData.Polynomial,
            mMeshData.NumberOfElementsPerDirection,
            mElementData.FeTopoMTK,
            mProcData.ElementListOnProc,
            mBasisData.BSplineActive,
            mSettings.Refinement,
            mSettings.TruncatedBsplines,
            mElementData.TMatrixParentChildRelationFEM,
            mBasisData.BsplineBasisList,
            mBasisData.BSplineMap );

    real tElapsedsave_dofconnectivity = save_dofconnectivity.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving elemental dof connectivity : %f [sec]\n",tElapsedsave_dofconnectivity);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_node_dof_connectivity_ascii()
{
    moris::tic save_dofconnectivity_new;

    mOutput.output_node_dof_connectivity_ascii(
            mBasisData.LagrangeToBSplineMap,
            mBasisData.IdFieldNodalField,
            mBasisData.TMatrixNodalField,
            mBasisData.BSplineMap,
            mSettings.Refinement);

    real tElapsedsave_dofconnectivity_new = save_dofconnectivity_new.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving nodal dof connectivity : %f [sec]\n",tElapsedsave_dofconnectivity_new);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_node_dof_connectivity_binary()
{
    moris::tic save_dofconnectivity_new;

    mOutput.output_node_dof_connectivity_binary(
            mBasisData.LagrangeToBSplineMap,
            mBasisData.IdFieldNodalField,
            mBasisData.TMatrixNodalField,
            mBasisData.BSplineMap,
            mSettings.Refinement);

    real tElapsedsave_dofconnectivity_new = save_dofconnectivity_new.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving nodal dof connectivity : %f [sec]\n",tElapsedsave_dofconnectivity_new);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_design_variables_for_moris_ascii()
{
    moris::tic  save_designvariables;

    mOutput.output_design_variables_for_moris_ascii(
            mBasisData.LagrangeToBSplineMap,
            mBasisData.IdFieldDesignNodalField,
            mBasisData.TMatrixDesignNodalField,
            mBasisData.DesignBSplineListMap,
            mSettings.Refinement,
            mSettings.UseSymmetry
    );

    real tElapsedsave_designvariables= save_designvariables.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving nodal design variables for MORIS: %f [sec]\n",tElapsedsave_designvariables);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_design_variables_for_femdoc_ascii()
{
    moris::tic  save_designvariables;

    mOutput.output_design_variables_for_femdoc_ascii(
            mBasisData.LagrangeToBSplineMap,
            mBasisData.IdFieldDesignNodalField,
            mBasisData.TMatrixDesignNodalField,
            mBasisData.DesignBSplineListMap,
            mSettings.Refinement,
            mSettings.UseSymmetry
    );

    real tElapsedsave_designvariables= save_designvariables.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving nodal design variables for FEMDOC: %f [sec]\n",tElapsedsave_designvariables);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_design_variables_for_moris_binary()
{
    moris::tic  save_designvariables;

    mOutput.output_design_variables_for_moris_binary(mBasisData.LagrangeToBSplineMap,
                                   mBasisData.IdFieldDesignNodalField,
                                   mBasisData.TMatrixDesignNodalField,
                                   mBasisData.DesignBSplineListMap,
                                   mSettings.Refinement,
                                   mSettings.UseSymmetry );

    real tElapsedsave_designvariables= save_designvariables.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving nodal design variables for MORIS: %f [sec]\n",tElapsedsave_designvariables);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_design_variables_for_femdoc_binary()
{
    moris::tic  save_designvariables;

    mOutput.output_design_variables_for_femdoc_binary(
            mBasisData.LagrangeToBSplineMap,
            mBasisData.IdFieldDesignNodalField,
            mBasisData.TMatrixDesignNodalField,
            mBasisData.DesignBSplineListMap,
            mSettings.Refinement,
            mSettings.UseSymmetry
    );

    real tElapsedsave_designvariables= save_designvariables.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving nodal design variables for FEMDOC: %f [sec]\n",tElapsedsave_designvariables);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_active_design_elements()
{
    moris::tic  save_active_design_elements;
    mOutput.output_active_design_elements(mElementData.ElementListActiveDesign,
                                          mSettings.Refinement);
    real tElapsedsave_active_design_elements = save_active_design_elements.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving active design elements: %f [sec]\n",tElapsedsave_active_design_elements);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_active_fem_elements()
{
    moris::tic  save_active_fem_elements;
    mOutput.save_active_fem_elements(mSettings.Refinement,mProcData.ElementListOnProc);
    real tElapsedsave_active_fem_elements = save_active_fem_elements.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving active FEM elements: %f [sec]\n",tElapsedsave_active_fem_elements);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_active_design_basis()
{
    moris::tic  save_active_design_basis;
    mOutput.save_active_design_basis(mSettings.Refinement,mBasisData.DesignBSplineActiveList);
    real tElapsedsave_active_design_basis = save_active_design_basis.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving active design basis : %f [sec]\n",tElapsedsave_active_design_basis);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_active_basis()
{
    moris::tic  save_active_basis;
    mOutput.save_active_basis(mSettings.Refinement,mBasisData.BSplineActiveList);
    real tElapsedsave_active_basis = save_active_basis.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving active basis : %f [sec]\n",tElapsedsave_active_basis);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_coordinate_list()
{
    moris::tic  save_coordinate_list;
    mOutput.save_coordinate_list(mSettings.Refinement,mBasisData.NodalLocaltoGlobal);
    real tElapsedsave_coordinate_list = save_coordinate_list.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving the coord list : %f [sec]\n",tElapsedsave_coordinate_list);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::save_all_data_files_for_continuation()
{
    moris::tic  savedata;

    this->save_coordinate_list();
    if ( mSettings.OutputFilesForFemdocAreBinary )
    {
        this->save_element_dof_connectivity_binary();
        this->save_node_dof_connectivity_binary();
        this->save_design_variables_for_femdoc_binary();
    }
    else
    {
        this->save_element_dof_connectivity_ascii();
        this->save_node_dof_connectivity_ascii();
        this->save_design_variables_for_femdoc_ascii();
    }
    
    if ( mSettings.InputOutputFilesForMorisAreAscii )
    {
        this->save_design_variables_for_moris_ascii();
    }
    else
    {
        this->save_design_variables_for_moris_binary();
    }

    this->save_active_design_elements();
    this->save_active_fem_elements();
    this->save_active_design_basis();
    this->save_active_basis();

    real tElapsedcsavedata = savedata.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for saving data in files : %f [sec]\n",tElapsedcsavedata);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::read_active_fem_elements()
{
    moris::tic  read_active_fem_elements;

    mInput.read_active_fem_elements(
            mProcData.ElementListOnProcLastStepFEM,
            mMeshData.LevelLastStep,
            mMeshData.ModelDim,
            mMeshData.NumberOfElementsPerDirection,
            mElementData.ElementActiveLastStep);

    real tElapsedread_active_fem_elements = read_active_fem_elements.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for reading active FEM elements : %f [sec]\n",tElapsedread_active_fem_elements);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::read_active_design_elements()
{
    moris::tic  read_active_design_elements;

    mInput.read_active_design_elements(
            mMeshData.ModelDim,
            mMeshData.NumberOfElementsPerDirection,
            mMeshData.LevelDesignLastStep,
            mProcData.ElementListOnProcLastStepDesign,
            mElementData.ElementActiveDesignLastStep);

    real tElapsedread_active_design_elements = read_active_design_elements.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for reading active design elements : %f [sec]\n",tElapsedread_active_design_elements);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::read_active_design_basis()
{
    moris::tic  read_active_design_basis;

    mInput.read_active_design_basis(
            mMeshData.ModelDim,
            mMeshData.PolynomialDesign,
            mMeshData.NumberOfElementsPerDirection,
            mBasisData.DesignBSplineActiveListLastStep,
            mBasisData.DesignBSplineActiveLastStep);

    real tElapsedread_active_design_basis = read_active_design_basis.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for reading active design basis : %f [sec]\n",tElapsedread_active_design_basis);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::read_active_basis()
{
    moris::tic  read_active_basis;

    mInput.read_active_basis(
            mMeshData.ModelDim,
            mMeshData.PolynomialDesign,
            mMeshData.NumberOfElementsPerDirection,
            mBasisData.BSplineActiveListLastStep,
            mBasisData.BSplineActiveLastStep);

    real tElapsedread_active_basis = read_active_basis.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for reading active basis : %f [sec]\n",tElapsedread_active_basis);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::read_design_variables_ascii()
{
    moris::tic  read_design_variables_ascii;

    mInput.read_design_variables_ascii(
            mBasisData.IdFieldDesignNodalField,
            mBasisData.TMatrixDesignNodalField);

    real tElapsedread_design_variables_ascii = read_design_variables_ascii.toc<moris::chronos::seconds>().wall;

    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for reading design variables : %f [sec]\n",tElapsedread_design_variables_ascii);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::read_coordinate_list()
{
    moris::tic  read_coordinate_list;
    mInput.read_coordinate_list(mBasisData.CoordinateList);
    real tElapsedread_coordinate_list = read_coordinate_list.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for reading coordinate list : %f [sec]\n",tElapsedread_coordinate_list);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::read_absdesvariables_file()
{
    moris::tic  read_absdesvariables_file;
    mInput.read_absdesvariables_file(
            mMeshData.ModelDim,
            mMeshData.PolynomialDesign,
            mMeshData.NumberOfElementsPerDirection,
            mMeshData.NumberOfBasisPerDirection,
            mBasisData.DesignBSplineActiveListLastStep,
            mFieldData.NodalADV,
            mFieldData.NodalPDV,
            mFieldData.NodalLevelSet,
            mFieldData.NodalSDFs,
            mFieldData.NodalFieldExists,
            mBasisData.IdFieldDesignNodalField,
            mBasisData.TMatrixDesignNodalField,
            mSettings.SimpExp,
            mSettings.ElementEdgeWidth,
            mSettings.InitialElementEdgeWidth,
            mSettings.LSscale,
            mSettings.LSthresh,
            mSettings.DensityShift,
            mSettings.Prjbeta,
            mSettings.SDFAlpha,
            mSettings.SDFBeta,
            ( ! mSettings.InputOutputFilesForMorisAreAscii ),
            mSettings.UseSymmetry,
            mSettings.SymmetryPlane,
            mMeshData.BasisSymmetryIndex );

    real tElapsedread_absdesvariables_file = read_absdesvariables_file.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for reading AbsDesVariables : %f [sec]\n",tElapsedread_absdesvariables_file);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::read_all_data_files_for_continuation()
{
    moris::tic read_data;

    this->read_active_fem_elements();
    this->read_active_design_elements();
    this->read_active_design_basis();
    this->read_active_basis();

    // read SDF field files
    if ( mSettings.SDFField )
    {
        this->read_sdf_for_active_design_basis();
    }

    //Read absdesVariables only, if no exofile is read
    if( mSettings.RefineBaseOnExoFile == false )
    {
        this->read_absdesvariables_file();
    }

    this->read_coordinate_list();

    real tElapsed_read_data = read_data.toc<moris::chronos::seconds>().wall;

    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for reading last step : %f [sec]\n",tElapsed_read_data);
}

//--------------------------------------------------------------------------------
void
Hierarchical_Mesh_Main::read_sdf_for_active_design_basis()
{
    // get number of SDF fields
    uint tNumberOfSDFs = mSettings.SDFObjectFiles.size();

    // clear nodal field map
    mFieldData.NodalSDFs.clear();

    // loop over all basis
    for ( uint k=0; k<tNumberOfSDFs; ++k )
    {
        // create temporary map object
        map< uint, real > tMap;

        // create name for input file
        std::string tFilePath = "SDF_" +  std::to_string(k) +".dat";

        // read field from input file
        mInput.read_field_from_file(
                tFilePath,
                mBasisData.BSplineActiveLastStep,
                mSettings.UseSymmetry,
                mBasisData.DesignBSplineActiveListLastStep,
                mMeshData.ModelDim,
                mMeshData.PolynomialDesign,
                mSettings.SymmetryPlane,
                mMeshData.BasisSymmetryIndex,
                mMeshData.NumberOfElementsPerDirection,
                mMeshData.NumberOfBasisPerDirection,
                tMap );

        // append map to field
        mFieldData.NodalSDFs.push_back( std::move( tMap ) );
    }
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::broadcast_bitset_logical_and(
        BoostBitset & aBitset)
{
    moris::tic sendrecvdeactivateionbitset;
    mMPI.broadcast_bitset_logical_and(aBitset);
    real tElapsedsendrecvdeactivateionbitset = sendrecvdeactivateionbitset.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for MPI and deactivation of the bitsets : %f [sec]\n",tElapsedsendrecvdeactivateionbitset);
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::broadcast_bitset_logical_or(BoostBitset & aBitset)
{
    moris::tic sendrecvdeactivateionbitset;
    mMPI.broadcast_bitset_logical_or(aBitset);
    real tElapsedsendrecvdeactivateionbitset = sendrecvdeactivateionbitset.toc<moris::chronos::seconds>().wall;
    if( mSettings.TimingTest == true )
        std::fprintf(stdout,"Time for MPI and activation of the bitsets : %f [sec]\n",tElapsedsendrecvdeactivateionbitset);
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Main::get_entities_owned_and_shared_by_current_proc(
        enum EntityRank   aEntityRank) const
{
    switch ( aEntityRank )
    {
        case( EntityRank::ELEMENT ):
                {
            return mProcData.ElementListOnProc;
                }
        case( EntityRank::NODE ):
                {
            //Number of nodes per Element
            uint tNumberBasisPerElement = pow( mMeshData.Polynomial + 1, mMeshData.ModelDim );
            //Create a list to get all active nodes
            Mat<uint> tBasis( mProcData.ElementListOnProc.length() * tNumberBasisPerElement, 1, 0 );
            Mat<uint> tBasisOfElement;
            //Set all active nodes of the loop in bitset and then count active ones
            for ( uint i = 0; i < mProcData.ElementListOnProc.length(); i++ )
            {
                tBasisOfElement = give_active_lagrange_basis_of_element( mProcData.ElementListOnProc( i ) );
                tBasis.rows(i*tNumberBasisPerElement, ( i + 1 ) * tNumberBasisPerElement - 1 ) = tBasisOfElement.rows( 0, tNumberBasisPerElement - 1 );
            }
            //Only unique is needed
            tBasis = unique( tBasis );
            return tBasis;
                }
        case( EntityRank::EDGE ):
                {
            MORIS_ASSERT( mMeshData.ModelDim >= 2 && mMeshData.ModelDim <=3, " Edges are only in 2D and 3D available");
            //Determine the number of edges per element
            uint tNumberElementEdges = 0;
            if ( mMeshData.ModelDim == 2 )
            {
                tNumberElementEdges = 4;
            }
            else if ( mMeshData.ModelDim == 3 )
            {
                tNumberElementEdges = 12;
            }
            Mat<uint> tElementEdges;
            //Collect all possible edges
            Mat<uint> tPossibleEdges( mProcData.ElementListOnProc.length() * tNumberElementEdges, 1 );
            for ( uint i = 0; i < mProcData.ElementListOnProc.length(); i++ )
            {
                tElementEdges = give_element_edges( mProcData.ElementListOnProc( i ) );
                tPossibleEdges.rows( i * tElementEdges.length(), ( i + 1 ) * tElementEdges.length() - 1 )
                                        = tElementEdges.rows( 0, tElementEdges.length() - 1 );
            }
            //Only unique is needed
            tPossibleEdges = unique( tPossibleEdges );
            return tPossibleEdges;
                }
        case( EntityRank::FACE ):
                {
            MORIS_ASSERT( mMeshData.ModelDim >= 2 && mMeshData.ModelDim <=3, " Faces are only in 2D and 3D available");
            Mat<uint> tElementFaces;
            //Collect all possible faces
            Mat<uint> tPossibleFaces( mProcData.ElementListOnProc.length() * 2 * mMeshData.ModelDim, 1 );
            for ( uint i = 0; i < mProcData.ElementListOnProc.length(); i++ )
            {
                tElementFaces = mBaseFace.give_element_faces( mProcData.ElementListOnProc(i),
                        mMeshData.ModelDim, mMeshData.NumberOfElementsPerDirection );
                tPossibleFaces.rows( i * tElementFaces.length(), ( i + 1 ) * tElementFaces.length() - 1 )
                                        = tElementFaces.rows( 0, tElementFaces.length() - 1 );
            }
            //Only unique is needed
            tPossibleFaces = unique( tPossibleFaces );
            return tPossibleFaces;
                }
        default:
        {
            MORIS_LOG_ERROR << "Specified mesh type not supported by MORIS";
            return Mat< uint >( 0, 0 );
        }
    }
}

//--------------------------------------------------------------------------------
Mat<uint>
Hierarchical_Mesh_Main::get_bspline_entities_owned_and_shared_by_current_proc(
        enum EntityRank   aEntityRank) const
{

    // refer to correct element list depending on mesh mode
    /* const Mat< uint > *tElementsOnProc = mBasisData.MeshImplementationUsesOutputMesh ?
            ( & mElementData.ElementListActiveDesign ) :
            ( & mProcData.ElementListOnProcLastStepDesign ); */

    switch ( aEntityRank )
    {
        case( EntityRank::ELEMENT ):
        {
            return mProcData.ElementListOnProcLastStepDesign;
        }
        case( EntityRank::NODE ):
        {
            // get level of max element
            uint tLevel = mBaseElement.give_element_level(
                    mProcData.ElementListOnProcLastStepDesign.max(),
                    mMeshData.ModelDim,
                    mMeshData.NumberOfElementsPerDirection);

            // get max number of nodes
            uint tNumberOfBasis = mBasis.give_number_of_basis(
                    tLevel,
                    mMeshData.ModelDim,
                    mMeshData.PolynomialDesign,
                    mMeshData.NumberOfElementsPerDirection);

            // initialize bitset
            BoostBitset tBasisFlags( tNumberOfBasis );

            // loop over all elements
            for ( uint k=0; k< mProcData.ElementListOnProcLastStepDesign.length(); ++k)
            {
                // get basis of element
                Mat< uint > tBasisOfElement = give_basis_of_element(
                        mProcData.ElementListOnProcLastStepDesign( k ) );

                // flag basis
                for ( uint i=0; i < tBasisOfElement.length(); ++i )
                {
                    tBasisFlags.set( tBasisOfElement ( i ) );
                }
            }

            // create outout node list
            Mat< uint > aBasisList( tBasisFlags.count(), 1);

            // reset counter
            uint tCount = 0;

            // loop over all basis
            for (uint j=0; j<tNumberOfBasis; ++j)
            {
               if ( tBasisFlags.test( j ) )
               {
                   // add basis to output
                   aBasisList( tCount ) = j;

                   // increment counter
                   ++tCount;
               }
            }

            return aBasisList;
        }
        default:
        {
            MORIS_LOG_ERROR << "Specified mesh type not supported by MORIS";
            return Mat< uint >( 0, 0 );
        }
    }
}

//--------------------------------------------------------------------------------

uint
Hierarchical_Mesh_Main::get_num_entities_owned_and_shared_by_current_proc(
        enum EntityRank   aEntityRank) const
{
    switch ( aEntityRank )
    {
        case( EntityRank::ELEMENT ):
                {
            return get_entities_owned_and_shared_by_current_proc( EntityRank::ELEMENT ).length();
                }
        case( EntityRank::NODE ):
                {
            return get_entities_owned_and_shared_by_current_proc( EntityRank::NODE ).length();
                }
        case( EntityRank::EDGE ):
                {
            return get_entities_owned_and_shared_by_current_proc( EntityRank::EDGE ).length();
                }
        case( EntityRank::FACE ):
                {
            return get_entities_owned_and_shared_by_current_proc( EntityRank::FACE ).length();
                }
        default:
        {
            MORIS_LOG_ERROR << "Specified mesh type not supported by MORIS";
            return 0;
        }
    }
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Main::get_entities_universal(
        enum EntityRank   aEntityRank) const
{
    switch ( aEntityRank )
    {
        case( EntityRank::ELEMENT ):
                {
            return this -> give_active_elements_with_aura();
                }
        case( EntityRank::NODE ):
                {
            //Determine the active elements with the aura
            Mat<uint> ElementListOnProcWithAura = this -> give_active_elements_with_aura();
            //Number of nodes per Element
            uint tNumberBasisPerElement = pow( mMeshData.Polynomial + 1, mMeshData.ModelDim );
            //Create a list to get all active nodes
            Mat<uint> tActiveBasisWithAura( ElementListOnProcWithAura.length() * tNumberBasisPerElement, 1, 0 );
            Mat<uint> tBasisOfElement;
            //Set all active nodes of the loop in bitset and then count active ones
            for ( uint i = 0; i < ElementListOnProcWithAura.length(); i++ )
            {
                tBasisOfElement = give_active_lagrange_basis_of_element( ElementListOnProcWithAura( i ) );
                tActiveBasisWithAura.rows(i*tNumberBasisPerElement, ( i + 1 ) * tNumberBasisPerElement - 1 ) = tBasisOfElement.rows( 0, tNumberBasisPerElement - 1 );
            }
            tActiveBasisWithAura = unique( tActiveBasisWithAura );
            return tActiveBasisWithAura;
                }
        case( EntityRank::EDGE ):
                {
            MORIS_ASSERT( mMeshData.ModelDim >= 2 && mMeshData.ModelDim <=3, " Edges are only in 2D and 3D available");
            //Determine first the number of edges per element
            uint tNumberElementEdges = 0;
            if( mMeshData.ModelDim == 2 )
            {
                tNumberElementEdges = 4;
            }
            else if( mMeshData.ModelDim == 3 )
            {
                tNumberElementEdges = 12;
            }
            Mat<uint> tElementListWithAura = this -> give_active_elements_with_aura();
            Mat<uint> tElementEdges;
            //Collect all possible edges
            Mat<uint> tPossibleEdges( tElementListWithAura.length() * tNumberElementEdges, 1 );
            for ( uint i = 0; i < tElementListWithAura.length(); i++ )
            {
                tElementEdges = give_element_edges( tElementListWithAura( i ) );
                tPossibleEdges.rows( i * tElementEdges.length(), ( i + 1  ) * tElementEdges.length() - 1 )
                                        = tElementEdges.rows( 0, tElementEdges.length() - 1 );
            }
            //Only unique edges are needed
            tPossibleEdges = unique( tPossibleEdges );
            return tPossibleEdges;
                }
        case( EntityRank::FACE ):
                {
            MORIS_ASSERT( mMeshData.ModelDim >= 2 && mMeshData.ModelDim <=3, " Faces are only in 2D and 3D available");
            //Determine active elements iwth aura
            Mat<uint> tElementListWithAura = this -> give_active_elements_with_aura();
            Mat<uint> tElementFaces;
            //Collect all possible faces
            Mat<uint> tPossibleFaces( tElementListWithAura.length() * 2 * mMeshData.ModelDim, 1 );
            for ( uint i = 0; i < tElementListWithAura.length(); i++ )
            {
                tElementFaces = mBaseFace.give_element_faces( tElementListWithAura(i),
                        mMeshData.ModelDim, mMeshData.NumberOfElementsPerDirection );
                tPossibleFaces.rows( i * tElementFaces.length(), ( i + 1 ) * tElementFaces.length() - 1 )
                                        = tElementFaces.rows( 0, tElementFaces.length() - 1 );
            }
            //Only unique faces are needed
            tPossibleFaces = unique(tPossibleFaces);
            return tPossibleFaces;
                }
        default:
        {
            MORIS_LOG_ERROR << "Specified mesh type not supported by MORIS";
            return Mat< uint >( 0, 0 );
        }
    }
}

//--------------------------------------------------------------------------------

uint
Hierarchical_Mesh_Main::get_num_entities_universal(
        enum EntityRank   aEntityRank) const
{
    switch ( aEntityRank )
    {
        case( EntityRank::ELEMENT ):
                {
            return get_entities_universal( EntityRank::ELEMENT ).length();
                }
        case( EntityRank::NODE ):
                {
            return get_entities_universal( EntityRank::NODE ).length();
                }
        case( EntityRank::EDGE ):
                {
            return get_entities_universal( EntityRank::EDGE ).length();
                }
        case( EntityRank::FACE ):
                {
            return get_entities_universal( EntityRank::FACE ).length();
                }
        default:
        {
            MORIS_LOG_ERROR << "Specified mesh type not supported by MORIS";
            return 0;
        }
    }
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Main::get_entities_owned_current_proc(
        enum EntityRank   aEntityRank) const
{
    switch ( aEntityRank )
    {
        case( EntityRank::ELEMENT ):
                {
            return mProcData.ElementListOnProc;
                }
        case( EntityRank::NODE ):
        {
            //Get locally owned and globally shared and determine then the own basis functions
            Mat<uint> tPossibleBasis = get_entities_owned_and_shared_by_current_proc(EntityRank::NODE);
            Mat<uint> tLocalOwnedBasis( tPossibleBasis.length(), 1 );
            uint tVar = 0;
            for ( uint i = 0; i < tPossibleBasis.length(); i++ )
            {
                if ( give_lagrange_basis_owner( tPossibleBasis( i ) ) == par_rank() )
                {
                    tLocalOwnedBasis( tVar ) = tPossibleBasis( i );
                    tVar++;
                }
            }
            tLocalOwnedBasis.resize( tVar, 1 );
            return tLocalOwnedBasis;
        }
        case( EntityRank::EDGE ):
        {
            //Get locally owned and globally shared and determine then the own edges
            Mat<uint> tPossibleEdges = get_entities_owned_and_shared_by_current_proc( EntityRank::EDGE );
            Mat<uint> tLocalOwnedEdges( tPossibleEdges.length(), 1 );
            uint tVar = 0;
            for ( uint i = 0; i < tPossibleEdges.length(); i++ )
            {
                if( give_edge_owner( tPossibleEdges( i ) ) == par_rank() )
                {
                    tLocalOwnedEdges( tVar ) = tPossibleEdges( i );
                    tVar++;
                }
            }
            tLocalOwnedEdges.resize( tVar, 1 );
            return tLocalOwnedEdges;
        }
        case( EntityRank::FACE ):
        {
            //Get locally owned and globally shared and determine then the own faces
            Mat<uint> tPossibleFaces = get_entities_owned_and_shared_by_current_proc( EntityRank::FACE );
            Mat<uint> tLocalOwnedFaces( tPossibleFaces.length(), 1 );
            uint tVar = 0;
            for ( uint i = 0; i < tPossibleFaces.length(); i++ )
            {
                if( give_face_owner( tPossibleFaces( i ) ) == par_rank() )
                {
                    tLocalOwnedFaces( tVar ) = tPossibleFaces( i );
                    tVar++;
                }
            }
            tLocalOwnedFaces.resize( tVar, 1 );
            return tLocalOwnedFaces;
        }
        default:
        {
            MORIS_LOG_ERROR << "Specified mesh type not supported by MORIS";
            return Mat< uint >( 0, 0 );
        }
    }
}

//--------------------------------------------------------------------------------

uint
Hierarchical_Mesh_Main::get_num_entities_owned_current_proc(
        enum EntityRank   aEntityRank) const
{
    switch ( aEntityRank )
    {
        case( EntityRank::ELEMENT ):
                {
            return get_entities_owned_current_proc( EntityRank::ELEMENT ).length();
                }
        case( EntityRank::NODE ):
                {
            return get_entities_owned_current_proc( EntityRank::NODE ).length();
                }
        case( EntityRank::EDGE ):
                {
            return get_entities_owned_current_proc( EntityRank::EDGE ).length();
                }
        case(EntityRank::FACE):
                {
            return get_entities_owned_current_proc( EntityRank::FACE ).length();
                }
        default:
        {
            MORIS_LOG_ERROR << "Specified mesh type not supported by MORIS";
            return 0;
        }
    }
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Main::get_entities_in_aura( enum EntityRank   aEntityRank) const
{
    switch (aEntityRank)
    {
        case(EntityRank::ELEMENT):
                                                                                             {
            return mProcData.ElementListOnProcAura;
                                                                                             }
        case(EntityRank::NODE):
                                                                                             {
            Mat<uint> tNodes_universal = get_entities_universal(EntityRank::NODE);
            BoostBitset tAuraBasis(tNodes_universal.max()+1);
            for(uint i = 0; i < tNodes_universal.length(); i++)
            {
                tAuraBasis.set(tNodes_universal(i));
            }
            Mat<uint> tNodes_owned_and_shared = get_entities_owned_and_shared_by_current_proc(EntityRank::NODE);
            for(uint i = 0; i < tNodes_owned_and_shared.length(); i++)
            {
                tAuraBasis.reset(tNodes_owned_and_shared(i));
            }
            uint tVar = 0;
            Mat<uint> tNodesInAura(tAuraBasis.count(),1);
            for(uint i = 0; i < tNodes_universal.length(); i++)
            {
                if( tAuraBasis.test(tNodes_universal(i)) == 1)
                {
                    tNodesInAura(tVar) = tNodes_universal(i);
                    tVar++;
                }
            }
            return tNodesInAura;
                                                                                             }
        case(EntityRank::EDGE):
                                                                                          {
            Mat<uint> tEdges_universal = get_entities_universal(EntityRank::EDGE);
            BoostBitset tAuraEdges(tEdges_universal.max()+1);
            for(uint i = 0; i < tEdges_universal.length(); i++)
            {
                tAuraEdges.set(tEdges_universal(i));
            }
            Mat<uint> tEdges_owned_and_shared = get_entities_owned_and_shared_by_current_proc(EntityRank::EDGE);
            for(uint i = 0; i < tEdges_owned_and_shared.length(); i++)
            {
                tAuraEdges.reset(tEdges_owned_and_shared(i));
            }
            uint tVar = 0;
            Mat<uint> tEdgesInAura(tAuraEdges.count(),1);
            for(uint i = 0; i < tEdges_universal.length(); i++)
            {
                if( tAuraEdges.test(tEdges_universal(i)) == 1)
                {
                    tEdgesInAura(tVar) = tEdges_universal(i);
                    tVar++;
                }
            }
            return tEdgesInAura;
                                                                                          }
        case(EntityRank::FACE):
                                                                       {
            Mat<uint> tFaces_universal = get_entities_universal(EntityRank::FACE);
            BoostBitset tAuraFaces(tFaces_universal.max()+1);
            for(uint i = 0; i < tFaces_universal.length(); i++)
            {
                tAuraFaces.set(tFaces_universal(i));
            }
            Mat<uint> tFaces_owned_and_shared = get_entities_owned_and_shared_by_current_proc(EntityRank::FACE);
            for(uint i = 0; i < tFaces_owned_and_shared.length(); i++)
            {
                tAuraFaces.reset(tFaces_owned_and_shared(i));
            }
            uint tVar = 0;
            Mat<uint> tFacesInAura(tAuraFaces.count(),1);
            for(uint i = 0; i < tFaces_universal.length(); i++)
            {
                if( tAuraFaces.test(tFaces_universal(i)) == 1)
                {
                    tFacesInAura(tVar) = tFaces_universal(i);
                    tVar++;
                }
            }
            return tFacesInAura;
                                                                       }
        default:
        {
            MORIS_LOG_ERROR << "Specified mesh type not supported by MORIS";
            return Mat< uint >( 0, 0 );
        }
    }
}

//--------------------------------------------------------------------------------

uint
Hierarchical_Mesh_Main::get_num_entities_in_aura( enum EntityRank   aEntityRank) const
{
    switch (aEntityRank)
    {
        case(EntityRank::ELEMENT):
                                                                                               {
            return get_entities_in_aura(EntityRank::ELEMENT).length();
                                                                                               }
        case(EntityRank::NODE):
                                                                                               {
            return get_entities_in_aura(EntityRank::NODE).length();
                                                                                               }
        case(EntityRank::EDGE):
                                                                                            {
            return get_entities_in_aura(EntityRank::EDGE).length();
                                                                                            }
        case(EntityRank::FACE):
                                                                         {
            return get_entities_in_aura(EntityRank::FACE).length();
                                                                         }
        default:
        {
            MORIS_LOG_ERROR << "Specified mesh type not supported by MORIS";
            return 0;
        }
    }
}

//--------------------------------------------------------------------------------

uint
Hierarchical_Mesh_Main::parallel_owner_rank_by_entity_id(
        uint aEntityID,
        enum EntityRank aEntityRank) const
{
    switch (aEntityRank)
    {
        case(EntityRank::NODE):
                                                                 {
            return give_basis_owner( aEntityID );
                                                                 }
        case(EntityRank::EDGE):
                                                                 {
            return give_edge_owner( aEntityID );
                                                                 }
        case(EntityRank::FACE):
                                                                 {
            return give_face_owner( aEntityID );
                                                                 }
        default:
        {
            MORIS_LOG_ERROR << "Specified mesh type not supported by MORIS";
            return 0;
        }
    }
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Main::get_procs_sharing_entity_by_id(
        uint aEntityID,
        enum EntityRank aEntityRank) const
{
    switch (aEntityRank)
    {
        case(EntityRank::NODE):
                                                                 {
            return give_basis_share( aEntityID );
                                                                 }
        case(EntityRank::EDGE):
                                                                 {
            return give_edge_share( aEntityID );
                                                                 }
        case(EntityRank::FACE):
                                                                 {
            return give_face_share( aEntityID );
                                                                 }
        default:
        {
            MORIS_LOG_ERROR << "Specified mesh type not supported by MORIS";
            return Mat< uint >( 0, 0 );
        }
    }
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Main::give_active_elements_of_edge(uint const & aEdgeId) const
{
    //Determine elements of edge
    Mat<uint> tElementsOfEdge =  mBaseEdge.give_elements_of_edge(aEdgeId,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
    tElementsOfEdge.print("tElementsOfEdge");
    //Temporary variable for loop
    uint tVar = 0;
    //Active elements of edge
    Mat<uint> tActiveElementsOfEdge( tElementsOfEdge.length(), 1, 0 );
    //Which element is it of the list
    Mat<uint> tWhichElementIsActive( tElementsOfEdge.length(), 1, 0 );
    //List of elements, which will finally determined
    Mat<uint> tActiveElementsOfEdgeFinal;
    //Check which elements are active
    for ( uint i = 0; i < tElementsOfEdge.length(); i++ )
    {
        if( mElementData.ElementActive.test( tElementsOfEdge( i ) ) == 1 )
        {
            tActiveElementsOfEdge( tVar ) = tElementsOfEdge( i );
            tWhichElementIsActive( tVar ) = i;
            tVar++;
        }
    }
    //If not all elements are active, then find the active neighbor elements
    if( tVar < tElementsOfEdge.length() )
    {
        //Check from the first active element all active neighbors
        Mat<uint> tAllActiveElementNeighbors;
        if( mMeshData.ModelDim == 2 )
        {
            //Use face neighbors to find elements, which are connected to an edge
            tAllActiveElementNeighbors = this->give_active_face_neighbor_of_element( tActiveElementsOfEdge( 0 ) );
        }
        else if( mMeshData.ModelDim == 3 )
        {
            //Use all neighbors (also diagonal neighbors), to find elements, which are connected to an ede
            //3D uses here internally a different numbering for the ordinals
            tAllActiveElementNeighbors = this->give_active_neighbor_of_element( tActiveElementsOfEdge( 0 ) );
        }
        tAllActiveElementNeighbors.print("tAllActiveElementNeighbors");
        // Compute the position of the edge
        Mat<uint> tEdgePosition = mBaseEdge.give_edge_position( aEdgeId, mMeshData.ModelDim, mMeshData.NumberOfElementsPerDirection );
        //Extract the right active element neighbors with the position of the edge
        tActiveElementsOfEdgeFinal.set_size( tAllActiveElementNeighbors.n_rows(), 1, 0 );
        tActiveElementsOfEdgeFinal( 0 ) = tActiveElementsOfEdge( 0 );
        tVar = 1;

        if( mMeshData.ModelDim == 2 )
        {
            //If edge is in x-direction
            if( tEdgePosition( 0 ) < UINT_MAX && tEdgePosition( 1 ) < UINT_MAX )
            {
                if( tWhichElementIsActive( 0 ) == 0 )
                {
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 2 )
                        {
                            tActiveElementsOfEdgeFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
                else if( tWhichElementIsActive( 0 ) == 1 )
                {
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 0 )
                        {
                            tActiveElementsOfEdgeFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
            }
            else if( tEdgePosition( 2 ) < UINT_MAX && tEdgePosition( 3 ) < UINT_MAX )
            {
                //If edge is in y-direction
                if( tWhichElementIsActive( 0 ) == 0 )
                {
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 1 )
                        {
                            tActiveElementsOfEdgeFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
                else if( tWhichElementIsActive( 0 ) == 1 )
                {
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 3 )
                        {
                            tActiveElementsOfEdgeFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
            }
        }
        else if( mMeshData.ModelDim == 3 )
        {
            MORIS_ASSERT( mMeshData.ModelDim < 2, " Edge to Element connectivity is right now not available in 3D " );
            //If edge is in x-direction
            if( tEdgePosition( 0 ) < UINT_MAX && tEdgePosition( 1 ) < UINT_MAX && tEdgePosition( 2 ) < UINT_MAX )
            {

            }
            else if( tEdgePosition( 3 ) < UINT_MAX && tEdgePosition( 4 ) < UINT_MAX && tEdgePosition( 5 ) < UINT_MAX )
            {
                //If edge is in y-direction

            }
            else if( tEdgePosition( 6 ) < UINT_MAX && tEdgePosition( 7 ) < UINT_MAX && tEdgePosition( 8 ) < UINT_MAX )
            {

            }
        }
    }

    tActiveElementsOfEdgeFinal.resize( tVar, 1 );
    return tActiveElementsOfEdgeFinal;
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Main::give_active_elements_of_face(uint const & aFaceId) const
{
    //Determine elements of Face
    Mat<uint> tElementsOfFace =  mBaseFace.give_elements_of_face(aFaceId,mMeshData.ModelDim,mMeshData.NumberOfElementsPerDirection);
    //Temporary variable for loop
    uint tVar = 0;
    //Active elements of Face
    Mat<uint> tActiveElementsOfFace( tElementsOfFace.length(), 1, 0 );
    //Which element is it of the list
    Mat<uint> tWhichElementIsActive( tElementsOfFace.length(), 1, 0 );
    //List of elements, which will finally determined
    Mat<uint> tActiveElementsOfFaceFinal;
    //Check which elements are active
    for ( uint i = 0; i < tElementsOfFace.length(); i++ )
    {
        if( mElementData.ElementActive.test( tElementsOfFace( i ) ) == 1 )
        {
            tActiveElementsOfFace( tVar ) = tElementsOfFace( i );
            tWhichElementIsActive( tVar ) = i;
            tVar++;
        }
    }
    //If not all elements are active, then find the active neighbor elements
    if( tVar < tElementsOfFace.length() )
    {
        //Check from the first active element all active neighbors
        Mat<uint> tAllActiveElementNeighbors = this->give_active_face_neighbor_of_element( tActiveElementsOfFace( 0 ) );
        // Compute the position of the Face
        Mat<uint> tFacePosition = mBaseFace.give_face_position( aFaceId, mMeshData.ModelDim, mMeshData.NumberOfElementsPerDirection );
        //Extract the right active element neighbors with the position of the Face
        tActiveElementsOfFaceFinal.set_size( tAllActiveElementNeighbors.n_rows(), 1, 0 );
        tActiveElementsOfFaceFinal( 0 ) = tActiveElementsOfFace( 0 );
        tVar = 1;
        //Check, which active neighbor elements are connected to the face
        if( mMeshData.ModelDim == 2 )
        {
            //Check if face in x or y-direction
            if( tFacePosition( 0 ) < UINT_MAX && tFacePosition( 1 ) < UINT_MAX )
            {
                //If first element is active, then take neighbors from right
                if( tWhichElementIsActive( 0 ) == 0 )
                {
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 1 )
                        {
                            tActiveElementsOfFaceFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
                else if( tWhichElementIsActive( 0 ) == 1 )
                {
                    //If second element is active, then take neighbors from left
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 3 )
                        {
                            tActiveElementsOfFaceFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
            }
            else if( tFacePosition( 2 ) < UINT_MAX && tFacePosition( 3 ) < UINT_MAX )
            {
                //If first element is active, then take neighbors from top
                if( tWhichElementIsActive( 0 ) == 0 )
                {
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 2 )
                        {
                            tActiveElementsOfFaceFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
                else if( tWhichElementIsActive( 0 ) == 1 )
                {
                    //If first element is active, then take neighbors from bottom
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 0 )
                        {
                            tActiveElementsOfFaceFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
            }
        }
        else if( mMeshData.ModelDim == 3 )
        {
            if( tFacePosition( 0 ) < UINT_MAX && tFacePosition( 1 ) < UINT_MAX && tFacePosition( 2 ) < UINT_MAX )
            {
                //If first element is active, then take neighbors from right
                if( tWhichElementIsActive( 0 ) == 0 )
                {
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 1 )
                        {
                            tActiveElementsOfFaceFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
                else if( tWhichElementIsActive( 0 ) == 1 )
                {
                    //If first element is active, then take neighbors from left
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 3 )
                        {
                            tActiveElementsOfFaceFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
            }
            else if( tFacePosition( 3 ) < UINT_MAX && tFacePosition( 4 ) < UINT_MAX && tFacePosition( 5 ) < UINT_MAX )
            {
                //If first element is active, then take neighbors from top
                if( tWhichElementIsActive( 0 ) == 0 )
                {
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 2 )
                        {
                            tActiveElementsOfFaceFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
                else if( tWhichElementIsActive( 0 ) == 1 )
                {
                    //If first element is active, then take neighbors from bottom
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 0 )
                        {
                            tActiveElementsOfFaceFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
            }
            else if( tFacePosition( 6 ) < UINT_MAX && tFacePosition( 7 ) < UINT_MAX && tFacePosition( 8 ) < UINT_MAX )
            {
                //If first element is active, then take neighbors from front
                if( tWhichElementIsActive( 0 ) == 0 )
                {
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 5 )
                        {
                            tActiveElementsOfFaceFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
                else if( tWhichElementIsActive( 0 ) == 1 )
                {
                    //If first element is active, then take neighbors from back
                    for ( uint i = 0; i < tAllActiveElementNeighbors.n_rows(); i++)
                    {
                        if( tAllActiveElementNeighbors( i, 1 ) == 4 )
                        {
                            tActiveElementsOfFaceFinal( tVar ) =  tAllActiveElementNeighbors( i, 0 );
                            tVar++;
                        }
                    }
                }
            }
        }
    }
    tActiveElementsOfFaceFinal.resize( tVar, 1 );
    return tActiveElementsOfFaceFinal;
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Main::give_parents_of_lagrange_basis(Mat<uint> const & aBasisList ) const
{
    //Level of basis
    uint tLevelOfBasis = 0;
    //"Parent" of basis
    uint tPossibleParentOfBasis = UINT_MAX;
    //Copy first the possible basis and then check if one of them has a "parent"
    Mat<uint> tParentOfBasis = aBasisList;
    //Determine the lowest possible basis id (Check for "parents" )
    Mat<uint> tLowestBasis;
    //Loop over all possible basis
    for ( uint i = 0; i < tParentOfBasis.length(); i++ )
    {
        tLevelOfBasis = give_lagrange_basis_level( tParentOfBasis( i ) );
        tPossibleParentOfBasis = tParentOfBasis( i );
        while( tLevelOfBasis >  0 && tPossibleParentOfBasis < UINT_MAX )
        {
            //Parent of Basis is there
            tParentOfBasis( i ) = tPossibleParentOfBasis;
            //Check for next parent
            tPossibleParentOfBasis = give_lagrange_basis_of_parent( tPossibleParentOfBasis );
        }
    }
    return tParentOfBasis;
}

//--------------------------------------------------------------------------------

/* Mat<uint>
Hierarchical_Mesh_Main::give_local_parents_of_lagrange_basis(Mat<uint> const & aBasisList ) const
{
    //Level of basis
    uint tLevelOfBasis = 0;
    //"Parent" of basis
    uint tPossibleParentOfBasis = UINT_MAX;
    //Copy first the possible basis and then check if one of them has a "parent"
    Mat<uint> tParentOfBasis = aBasisList;
    //Determine the lowest possible basis id (Check for "parents" )
    Mat<uint> tLowestBasis;
    //Loop over all possible basis
    for ( uint i = 0; i < tParentOfBasis.length(); i++ )
    {
        tLevelOfBasis = give_local_lagrange_basis_level( tParentOfBasis( i ) );
        tPossibleParentOfBasis = tParentOfBasis( i );
        while( tLevelOfBasis >  0 && tPossibleParentOfBasis < UINT_MAX )
        {
            //Parent of Basis is there
            tParentOfBasis( i ) = tPossibleParentOfBasis;
            //Check for next parent
            tPossibleParentOfBasis = give_local_lagrange_basis_of_parent( tPossibleParentOfBasis );
        }
    }
    return tParentOfBasis;
} */

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Main::give_active_lagrange_basis_of_element(uint const & aElementId) const
{
    //Determine the possible lagrange basis functions
    Mat<uint> tPossibleBasis = give_lagrange_basis_of_element( aElementId );
    //Determine parent basis functions, if possible.
    return give_parents_of_lagrange_basis( tPossibleBasis );
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Main::give_edge_active_lagrange_nodes(uint const & aEdgeId) const
{
    //Determine the possible lagrange basis functions
    Mat<uint> tPossibleBasis = give_edge_lagrange_nodes( aEdgeId );
    //Determine parent basis functions, if possible.
    return give_parents_of_lagrange_basis( tPossibleBasis );
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Main::give_face_active_lagrange_basis(uint const & aFaceId) const
{
    //Determine the possible lagrange basis functions
    Mat<uint> tPossibleBasis = give_face_lagrange_basis( aFaceId );
    //Determine parent basis functions, if possible.
    return give_parents_of_lagrange_basis( tPossibleBasis );
}

//--------------------------------------------------------------------------------

/* Mat<uint>
Hierarchical_Mesh_Main::give_local_active_lagrange_basis_of_element(uint const & aElementId) const
{
    //Determine the possible lagrange basis functions
    Mat<uint> tPossibleBasis = give_local_lagrange_basis_of_element( aElementId );
    //Determine parent basis functions, if possible.
    return give_local_parents_of_lagrange_basis( tPossibleBasis );
}*/

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::read_data_from_exo_file(
        std::string & aGFile)
{
    //Load file and read data from file
    mesh tExoMesh( MeshType::MTK, aGFile );
    Mat<real> tNodalLvlsetField = tExoMesh.get_field_values( EntityRank::NODE, "nodlevelset" );
    Mat<real> tDensity = tExoMesh.get_intersected_data_field_set(EntityRank::ELEMENT, "matpropstrucdensity", "block_2");
    //Save nodal data in struct at the right position (active nodes)
    Mat<real> tNodalLvlsetFieldOld( mBasisData.CoordinateList.max() + 1, 1, 0);
    for ( uint i = 0; i < mBasisData.CoordinateList.length(); i++ )
    {
        tNodalLvlsetFieldOld( mBasisData.CoordinateList( i ) ) = tNodalLvlsetField( i );
    }
    //Map solution to Bspines (Least square fit)
    //Matrix of transpose(T) * T
    Mat<real> tMatrixM;
    //Right hand side for the least square fit
    Mat<real> tRightHandSide;
    //Solution of the least suare fit
    Mat<real> tBasisSolution;
    //Lagrange basis of an element
    Mat<uint> tLagrangeBasisOfElement;
    //Solution vector of elemental lagrange nodes
    Mat<real> tLagrangeSolution;
    //Compute size of vector
    //mFieldData.NodalLevelSet.set_size( mBasisData.DesignBSplineActiveListLastStep.max() + 1, 1, 0);
    for ( uint i = 0; i < mProcData.ElementListOnProcLastStepDesign.length(); i++ )
    {
        //Compute Lagrange basis of the element
        tLagrangeBasisOfElement = this->give_active_lagrange_basis_of_element( mProcData.ElementListOnProcLastStepDesign( i ) );
        //Get the elemental lagrange solution
        tLagrangeSolution.set_size( tLagrangeBasisOfElement.length(), 1, 0);
        for ( uint j = 0; j < tLagrangeBasisOfElement.length(); j++ )
        {
            tLagrangeSolution( j ) = tNodalLvlsetFieldOld( tLagrangeBasisOfElement( j ) );
        }
        //Compute the Tmatrix and Id field
         this->give_Tmatrix_and_IdField_DesignVariables_LastStep( mProcData.ElementListOnProcLastStepDesign( i ) );
         //Create the Mass matrix with trans(T) * T
        tMatrixM = trans( mElementData.TMatrixDesign ) * mElementData.TMatrixDesign;
        //Create right hand side with trans(T) * solutionvector
        tRightHandSide =  trans( mElementData.TMatrixDesign ) * tLagrangeSolution;
        //Solve the Bspline coefficients
        tBasisSolution = solve( tMatrixM, tRightHandSide, "superlu");
        //Put the solutions in the bspline solution vector
        for ( uint j = 0; j < mElementData.IdFieldDesign.length(); j++ )
        {
            mFieldData.NodalLevelSet[  mElementData.IdFieldDesign( j ) ] = tBasisSolution( j );
        }
    }

    //Save element data in struct at the right position (active elements)
    //mFieldData.ElementDensity.set_size( mProcData.ElementListOnProcLastStepFEM.max() + 1, 1, 0);
    for ( uint i = 0; i < mProcData.ElementListOnProcLastStepFEM.length(); i++ )
    {
        mFieldData.ElementDensity[ mProcData.ElementListOnProcLastStepFEM( i ) ] = tDensity( i );
    }
}

//--------------------------------------------------------------------------------

// rearranges the Lagrange DOFs in the order MTK needs.
Mat<uint>
Hierarchical_Mesh_Main::give_vector_for_reorder() const
{
    Mat<uint> tOrder( pow( mMeshData.Polynomial + 1, mMeshData.ModelDim ), 1, 0);
    if ( mMeshData.ModelDim == 2)
    {
        if( mMeshData.Polynomial == 1 )
        {
            tOrder( 0 ) = 0;
            tOrder( 1 ) = 1;
            tOrder( 2 ) = 3;
            tOrder( 3 ) = 2;
        }
        else if( mMeshData.Polynomial == 2 )
        {
            tOrder( 0 ) = 0;
            tOrder( 1 ) = 4;
            tOrder( 2 ) = 1;
            tOrder( 3 ) = 7;
            tOrder( 4 ) = 8;
            tOrder( 5 ) = 5;
            tOrder( 6 ) = 3;
            tOrder( 7 ) = 6;
            tOrder( 8 ) = 2;
        }
    }
    else if ( mMeshData.ModelDim == 3)
    {
        if( mMeshData.Polynomial == 1 )
        {
            tOrder( 0 ) = 0;
            tOrder( 1 ) = 1;
            tOrder( 2 ) = 3;
            tOrder( 3 ) = 2;
            tOrder( 4 ) = 4;
            tOrder( 5 ) = 5;
            tOrder( 6 ) = 7;
            tOrder( 7 ) = 6;
        }
        else if( mMeshData.Polynomial == 2 )
        {
            tOrder( 0 ) = 1;
            tOrder( 1 ) = 13;
            tOrder( 2 ) = 5;
            tOrder( 3 ) = 9;
            tOrder( 4 ) = 24;
            tOrder( 5 ) = 17;
            tOrder( 6 ) = 2;
            tOrder( 7 ) = 14;
            tOrder( 8 ) = 6;
            tOrder( 9 ) = 8;
            tOrder( 10 ) = 25;
            tOrder( 11 ) = 16;
            tOrder( 12 ) = 21;
            tOrder( 13 ) = 20;
            tOrder( 14 ) = 22;
            tOrder( 15 ) = 10;
            tOrder( 16 ) = 26;
            tOrder( 17 ) = 18;
            tOrder( 18 ) = 0;
            tOrder( 19 ) = 12;
            tOrder( 20 ) = 4;
            tOrder( 21 ) = 11;
            tOrder( 22 ) = 23;
            tOrder( 23 ) = 19;
            tOrder( 24 ) = 3;
            tOrder( 25 ) = 15;
            tOrder( 26 ) = 7;
        }
    }
    return tOrder;
}

//--------------------------------------------------------------------------------

uint
Hierarchical_Mesh_Main::convert_lagrange_id_to_bspline_id(
        const uint aMaxLevel,
        const uint aGlobalLagrangeID) const
{
    // get the i-j-k position from the global element ID
    Mat<uint> tIJK = give_position_of_lagrange_basis(aGlobalLagrangeID);

    // calculate level of basis
    uint tLevel = give_lagrange_basis_level(aGlobalLagrangeID);

    bool_t tBSplineIsActive = false;
    uint aBSplineID = UINT_MAX;

    while ( ! tBSplineIsActive && tLevel <= aMaxLevel )
    {
        // get the id of the current B-Spline Basis
        aBSplineID = give_basis_of_position( tLevel, tIJK );

        // test if B-Spline is active
        if ( aBSplineID < mBasisData.BSplineActive.size() )
        {
            tBSplineIsActive = mBasisData.BSplineActive.test(aBSplineID);
        }
        else
        {
            tBSplineIsActive = false;
        }
        if ( ! tBSplineIsActive )
        {
            // increase level
            ++tLevel;

            // multiply everything for next level
            for(uint k=0; k<mMeshData.ModelDim; ++k)
            {
                tIJK(k) *= 2;
            }
        }
    }

    // set BSpline ID to UINT_MAX if node is hanging
    if ( tLevel > aMaxLevel)
    {
        aBSplineID = UINT_MAX;
    }

    return aBSplineID;
}

//--------------------------------------------------------------------------------

uint
Hierarchical_Mesh_Main::convert_lagrange_id_to_design_bspline_id(
        const uint aMaxLevel,
        const uint aGlobalLagrangeID) const

{
    // get the i-j-k position from the global element ID
    Mat<uint> tIJK = give_position_of_lagrange_basis(aGlobalLagrangeID);

    // calculate level of basis
    uint tLevel = give_lagrange_basis_level(aGlobalLagrangeID);

    bool_t tBSplineIsActive = false;
    uint aBSplineID = UINT_MAX;

    while ( ! tBSplineIsActive && tLevel <= aMaxLevel )
    {
        // get the id of the current B-Spline Basis
        aBSplineID = give_basis_of_position( tLevel, tIJK );

        // test if B-Spline is active
        tBSplineIsActive = mBasisData.BSplineActive.test(aBSplineID);

        if ( ! tBSplineIsActive )
        {
            // increase level
            ++tLevel;

            // multiply everything for next level
            for(uint k=0; k<mMeshData.ModelDim; ++k)
            {
                tIJK(k) *= 2;
            }
        }
    }

    // in contrast to convert_lagrange_id_to_bspline_id,
    // the basis is not overwritten if it is a hanging node

    return aBSplineID;
}

//--------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::create_maps()
{

    // std::fprintf(stdout, " Proc %u : creating maps ...\n", (uint) par_rank()) ;

    moris::tic tTimer;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Step 0: Perform consistency of B-Spline and Element activation
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // check for debugging purpose
    check_active_elements();
    check_active_basis();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Step 1: Create consecutive Lagrange node numbering
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    uint tNumberOfBasis =  mBasis.give_number_of_basis(
            mMeshData.Level + mMeshData.MaxPolynomial-1,
            mMeshData.ModelDim,
            mMeshData.MaxPolynomial,
            mMeshData.NumberOfElementsPerDirection);

    // create Lagrange bitset
    mBasisData.LagrangeActive.resize( tNumberOfBasis );

    // reset Lagrange bitset
    mBasisData.LagrangeActive.reset();

    // set Lagrange bitset
    for(uint i=0; i < mProcData.LagrangeListOnProc.length(); ++i)
    {
        mBasisData.LagrangeActive.set( mProcData.LagrangeListOnProc( i ) );
    }

    // broadcast Lagrange Bitset
    broadcast_bitset_logical_and( mBasisData.LagrangeActive );
    // create map for consecutive numbering

    // counter for basis
    uint tCount = 0;
    map< uint, uint > tLagrangeIdToGlobalInd;
    for(uint i=0; i < mBasisData.LagrangeActive.size(); ++i)
    {
        if ( mBasisData.LagrangeActive.test( i ) )
        {
            tLagrangeIdToGlobalInd[ i ] = ++tCount;
        }
    }

    // assign size for field NodalLocaltoGlobal
    // FIXME: make this one work for parallel
    mBasisData.NodalLocaltoGlobal.set_size( mProcData.LagrangeListOnProc.length(), 1);

    for(uint i=0; i < mProcData.LagrangeListOnProc.length(); i++)
    {
        mBasisData.NodalLocaltoGlobal(i) = i+1;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Step 2: Create Map for Node-Dof connectivity
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // create basis map for B-Splines
    mBasisData.BSplineMap.clear();
    tCount = 0;

    for(uint i=0; i < mBasisData.BSplineActive.size(); ++i)
    {
        if ( mBasisData.BSplineActive.test( i ) )
        {

            mBasisData.BSplineMap[ i ] = ++tCount;
        }
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Step 3: Create Map for Design variables
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // FIXME: make sure that this also works for parallel
    Mat<uint> tDesignBSplineList(mBasisData.NodalLocaltoGlobal.length()*mBasisData.IdFieldDesignNodalField.n_cols(),1,UINT_MAX);

    tCount = 0;
    for(uint i = 0; i<mBasisData.NodalLocaltoGlobal.length(); ++i)
    {
        for(uint j = 0; j< mBasisData.IdFieldDesignNodalField(i,0); j++)
        {
            tDesignBSplineList( tCount ) = mBasisData.IdFieldDesignNodalField(i,j+1);
            ++tCount;
        }
    }

    tDesignBSplineList = unique(tDesignBSplineList);
    tDesignBSplineList.resize(tDesignBSplineList.length()-1,1);

    // mBasisData.DesignBSplineListMap.set_size(tDesignBSplineList.max()+1,1,0);
    mBasisData.DesignBSplineListMap.clear();

    if ( mSettings.UseSymmetry )
    {
        tCount = 0;

        // loop over all non-mirror basis
        for(uint i = 0; i<tDesignBSplineList.length(); ++i)
        {
            if ( basis_is_symmetry_master (  tDesignBSplineList(i) ) )
            {
                mBasisData.DesignBSplineListMap[ tDesignBSplineList(i) ] = tCount++;
            }
        }

        // reset master counter
        mMeshData.NumberOfSymmetryMasterBasis = 0;
        // loop over all mirrored basis
        for(uint i = 0; i<tDesignBSplineList.length(); ++i)
        {
            if ( ! basis_is_symmetry_master (  tDesignBSplineList(i) ) )
            {
                uint tMaster = give_symmetry_master_of_basis( tDesignBSplineList(i) );
                mBasisData.DesignBSplineListMap[ tDesignBSplineList(i) ]
                    = mBasisData.DesignBSplineListMap.find( tMaster );
            }
            else
            {
                // increment master counter
                ++mMeshData.NumberOfSymmetryMasterBasis;
            }
        }
    }
    else
    {
        for(uint i = 0; i<tDesignBSplineList.length(); ++i)
        {
            mBasisData.DesignBSplineListMap[ tDesignBSplineList(i) ] = i;
        }

        // all basis are master if there is no symmerty
        mMeshData.NumberOfSymmetryMasterBasis = tDesignBSplineList.length();
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Step 4a: Create temporary Lagrange list map
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // FIXME is max level of elements > max level of basis?
    // get maximum level of mesh
    //uint tMaxLevel = give_basis_level( mBasisData.LagrangeActive.size()-1 );
    uint tMaxLevel = mMeshData.Level + mMeshData.MaxPolynomial-1;
    // assign memory for map
    uint tNumberOfNodes = mBasisData.LagrangeActive.count();

    Mat<uint> tLagrangeToBSplineMap( tNumberOfNodes, 1, UINT_MAX );

    // loop over all basis and create map
    tCount = 0;
    for (uint k=0; k< mBasisData.LagrangeActive.size(); ++k)
    {
        if ( mBasisData.LagrangeActive.test( k ) )
        {
            tLagrangeToBSplineMap( tCount++ ) = convert_lagrange_id_to_bspline_id( tMaxLevel, k );
        }
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Step 4b: Create final Lagrange list map
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // count not-hanging nodes
    tCount = 0;
    for( uint i = 0; i<tNumberOfNodes; ++i )
    {
        if (tLagrangeToBSplineMap(i) < UINT_MAX)
        {
            ++tCount;
        }
    }

    // save number of non hanging nodes
    mBasisData.NumberOfNonHangingNodes = tCount;

    // create map for Lagrange basis
    mBasisData.LagrangeToBSplineMap.set_size( tNumberOfNodes, 1, UINT_MAX );

    for(uint i = 0; i< tNumberOfNodes; ++i)
    {
        // check if node is not hanging
        if ( tLagrangeToBSplineMap(i) < UINT_MAX )
        {
            mBasisData.LagrangeToBSplineMap(i) = mBasisData.BSplineMap.find(tLagrangeToBSplineMap(i));
            //std::cout << i+1 << " " << mBasisData.LagrangeToBSplineMap(i) << std::endl;
        }
        else
        {
            // node is hanging
            ++tCount;
            mBasisData.LagrangeToBSplineMap(i) = tCount;
            //std::cout << i+1 << " " << mBasisData.LagrangeToBSplineMap(i) << std::endl;
        }
    }

    // save number of hanging nodes
    mBasisData.NumberOfHangingNodes = tNumberOfNodes - mBasisData.NumberOfNonHangingNodes;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Step 5: Create element topology
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Do not overwrite Fetopo. Create new field
    mElementData.FeTopoMTK.set_size(
            mElementData.FeTopo.n_rows(), mElementData.FeTopo.n_cols());

    for(uint j = 0; j < mElementData.FeTopo.n_cols(); j++)
    {
        for(uint i = 0; i < mElementData.FeTopo.n_rows(); i++)
        {
            mElementData.FeTopoMTK(i,j) =  mBasisData.LagrangeToBSplineMap(
                    tLagrangeIdToGlobalInd.find(mElementData.FeTopo(i,j))-1);
        }
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Step 6: Assert that B-Spline to Lagrange Map and FE-Topology are consistent
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // make sure that element connectivity is correct
    check_dof_consistency();

    // create a bitset for the map
    moris::BoostBitset tMapFlags( mBasisData.LagrangeToBSplineMap.length() );

    // std::cout << "max basis " <<  mBasisData.LagrangeToBSplineMap.max() << std::endl;
    // std::cout << "length    " <<  mBasisData.LagrangeToBSplineMap.length() << std::endl;

    // flag bitset
    for( uint k = 0; k < mBasisData.LagrangeToBSplineMap.length(); ++k )
    {
        tMapFlags.set( mBasisData.LagrangeToBSplineMap( k ) - 1 );
    }
    // create a bitset for the topology
    moris::BoostBitset tTopoFlags( mBasisData.LagrangeToBSplineMap.length() );

    // flag bitset
    for( uint j = 0; j < mElementData.FeTopoMTK.n_cols(); ++j )
    {
        for( uint i = 0; i < mElementData.FeTopoMTK.n_rows(); ++i )
        {
            tTopoFlags.set( mElementData.FeTopoMTK( i, j ) - 1 );
        }
    }

    //MORIS_ASSERT( mBasisData.LagrangeToBSplineMap.length() == tMapFlags.count(), "B-Spline to Lagrange Map not consistent" );
    //MORIS_ASSERT( mBasisData.LagrangeToBSplineMap.length() == tTopoFlags.count(), "FE-Topology not consistent" );
    // FIXME: do this error message the proper way
    if ( mBasisData.LagrangeToBSplineMap.length() != tMapFlags.count()  )
    {
        std::cout <<  "Error: B-Spline to Lagrange Map not consistent" << std::endl;
        exit(1);
    }

    // FIXME: do this error message the proper way
    if ( mBasisData.LagrangeToBSplineMap.length() != tTopoFlags.count() )
    {
        std::cout <<  "FE-Topology not consistent" << std::endl;
        exit(1);
    }
    real tElapsedTime = tTimer.toc<moris::chronos::seconds>().wall;
    // print timing information
    if( mSettings.TimingTest == true )
        std::fprintf( stdout,"Time for creating Maps: %f [sec]\n", tElapsedTime );

}

// -------------------------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::check_active_basis()
{
    // loop over all possible basis
    for ( uint k = 0; k < mBasisData.BSplineActive.size(); ++k )
    {
        // check if basis if active
        if ( mBasisData.BSplineActive.test( k ) )
        {
            // get level of basis
            uint tLevel = give_basis_level( k ) ;

            // get position of this basis
            Mat < uint > tPosition = give_position_of_basis( k );

            while ( tLevel > 0 )
            {

                // check if there is a basis above
                bool tUpperBasisExist = true;
                for ( uint i=0; i<mMeshData.ModelDim; ++i )
                {
                    tUpperBasisExist = tUpperBasisExist && ( tPosition(i) % 2 == 0 ) ;
                }

                // go one basis up if there exists one
                if ( tUpperBasisExist )
                {
                    for ( uint i=0; i<mMeshData.ModelDim; ++i )
                    {
                        tPosition(i) /= 2;
                    }

                    // go one level up
                    --tLevel;

                    // get the number of this basis
                    uint tBasis = give_basis_of_position( tLevel, tPosition );

                    /*if( mBasisData.BSplineActive.test( tBasis ) )
                                std::cout << "bad basis " << tBasis << std::endl;*/

                    // check if this basis is active
                    //MORIS_ASSERT( ! mBasisData.BSplineActive.test( tBasis ),
                    //              " Error: Duplicate basis activated " );

                    // FIXME: do this error message the proper way
                    if ( mBasisData.BSplineActive.test( tBasis ) )
                    {
                        std::cout <<  "Error: Duplicate basis activated" << std::endl;
                        exit(1);
                    }

                }
                else
                {
                    // continue with next basis k in loop
                    break;
                }
            }
        }
    }
}

// -------------------------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::check_active_elements()
{
    for (uint e=0; e < mProcData.ElementListOnProc.length(); ++e)
    {
        // get ID of element

        uint tElementID = mProcData.ElementListOnProc( e );

        // get level of element
        uint tLevel = give_element_level( tElementID );

        while ( tLevel > 0 )
        {
            // get parent of element
            tElementID = give_parent_of_element( tElementID );

            // MORIS_ASSERT( ! mElementData.ElementActive.test( tElementID ),
            //                                 " Error: Element parent activated " );

            // FIXME: do this error message the proper way

            if ( mElementData.ElementActive.test( tElementID ) )
            {
                std::cout <<  "Error: Element parent activated: "<< e << " " << mProcData.ElementListOnProc( e ) << std::endl;
                exit(1);
            }
            // this one does not work yet
            //MORIS_ASSERT( mElementData.ElementRefined.test( tElementID ),
            //                                 " Error: Element parent not refined " );

            // calculate level of element for next loop
            tLevel = give_element_level( tElementID );
        }
    }
}

// -------------------------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::calculate_ADV_on_new_field()
{

    // assign memory
    mFieldData.LagrangeADVField.set_size( mElementData.FeCoord.n_rows(), 1, 0 );

    // check if project ADV or first SDF
    if ( ! mSettings.SDFField  || mSettings.Refinement == 1 )
    {
        // call L2 Projection for NodalADV defined on last bspline basis mesh
        L2_projection( mFieldData.NodalADV );
    }
    else
    {
        // compute scaling and offset for adv field needed for proper reconstruction of level set field in femdoc
        real tScaling = -1.0/(mSettings.LSscale*mSettings.InitialElementEdgeWidth);
        real tOffset  = mSettings.LSthresh;

        // call L2 Projection SDF-0 defined on last bspline basis mesh
        L2_projection( mFieldData.NodalSDFs(0), tScaling, tOffset );
    }

    // create Lagrange field for MTK output
    mFieldData.LagrangeADVField.set_size( mElementData.FeCoord.n_rows(), 1, 0.0);

    for(uint i=0; i<mElementData.FeCoord.n_rows(); i++)
    {
        for(uint j = 0; j < mBasisData.IdFieldDesignNodalField(i,0); j++)
        {
            mFieldData.LagrangeADVField(i) += mBasisData.TMatrixDesignNodalField(i,j+1)*
                (mFieldData.BSplineCoeff)(mBasisData.DesignBSplineListMap.find(mBasisData.IdFieldDesignNodalField(i,j+1)));
        }
    }

    /* FILE* fid=std::fopen("AbsDesVariables_new.dat","wb");
    double * buffer = (double*) alloca(mFieldData.BSplineCoeff.length()*sizeof(double));
    for ( uint i = 0; i < mFieldData.BSplineCoeff.length(); i++)
        buffer[i] = mFieldData.BSplineCoeff(i);
    std::fwrite(buffer,sizeof(double),mFieldData.BSplineCoeff.length(),fid);
    std::fclose(fid); */

    if ( mSettings.UseSymmetry )
    {
        real tError = 0;

        // make sure that B-Spline Coeffs are close to zero for all slave entries
        for ( uint k=mMeshData.NumberOfSymmetryMasterBasis; k<mFieldData.BSplineCoeff.length(); ++k )
        {
            tError += pow( mFieldData.BSplineCoeff( k ), 2 );
        }
        tError = pow( tError, 0.5 );
        MORIS_ASSERT( tError < 1e-6, "something went wrong with symmetry and L2 projection");

        // resize output vector
        mFieldData.BSplineCoeff.resize( mMeshData.NumberOfSymmetryMasterBasis, 1 );
    }
    save_vector_to_binary_file( "AbsDesVariables_new.dat", mFieldData.BSplineCoeff );
}

// -------------------------------------------------------------------------------------------------
void
Hierarchical_Mesh_Main::create_nodal_fields_from_sdf()
{
    // get basis list for element topology
    Mat< uint > tBasisList = get_bspline_entities_owned_and_shared_by_current_proc( EntityRank::NODE );

    // get number of SDF fields
    uint tNumberOfSDFs = mFieldData.ObjectSDF.size();

    // clear cell
    mFieldData.NodalSDFs.clear();

    // loop over all SDFs
    for ( uint k=0; k<tNumberOfSDFs; ++k )
    {
        // crete temporary map
        map< uint, real > tMap;

        // create map for SDF k
        for ( uint i=0; i < mFieldData.ObjectSDF( k ).length(); ++i )
        {
            tMap[ tBasisList( i ) ] = mFieldData.ObjectSDF( k )( i ) ;
        }

        // add map to cell
        mFieldData.NodalSDFs.push_back( std::move( tMap ) );
    }
}

// -------------------------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::create_field_maps_from_sdf(
        const Mat< real > & aSignOfSDF,
        const Mat< uint > & aElementsAtSurface,
        const Mat< uint > & aElementsInVolume )
{
    // get basis list for element topology
    Mat< uint > tBasisList = get_bspline_entities_owned_and_shared_by_current_proc( EntityRank::NODE );

    // clear sdf map
    mFieldData.NodalSignOfSDF.clear();

    // transform moris::Mat to moris::map
    for ( uint i=0; i < aSignOfSDF.length(); ++i )
    {
        mFieldData.NodalSignOfSDF[ tBasisList( i ) ] = aSignOfSDF( i );
    }

    // initialize surface density map
    mFieldData.ElementSurfaceDensitySDF.clear();

    // loop over all active elements
    for ( uint e=0; e<mProcData.ElementListOnProcLastStepDesign.length(); ++e )
    {
        mFieldData.ElementSurfaceDensitySDF[ mProcData.ElementListOnProcLastStepDesign( e ) ] = 0.0;
    }

    // write volume density=1 into elements from volume
    for ( uint k=0; k< aElementsAtSurface.length(); ++k )
    {
        mFieldData.ElementSurfaceDensitySDF[ aElementsAtSurface( k ) ] = 1.0;
    }

    // initialize volume density map
    mFieldData.ElementVolumeDensitySDF.clear();

    // loop over all active elements
    for ( uint e=0; e<mProcData.ElementListOnProcLastStepDesign.length(); ++e )
    {
        mFieldData.ElementVolumeDensitySDF[ mProcData.ElementListOnProcLastStepDesign( e ) ] = 0.0;
    }

    // write volume density=1 into elements from volume
    for ( uint k=0; k< aElementsInVolume.length(); ++k )
    {
        mFieldData.ElementVolumeDensitySDF[ aElementsInVolume( k ) ] = 1.0;
    }
}

// -------------------------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::init_sdf_bitset()
{
    // get number of SDF objects
    uint tNumberOfObjects = mSettings.SDFObjectFiles.size();

    // set size of bitset
    mFieldData.RefineAgainstObjectSDF.resize( tNumberOfObjects );

    // loop over all objects and turn them on
    for ( uint k=0; k<tNumberOfObjects; ++k )
    {
        mFieldData.RefineAgainstObjectSDF.set( k );
    }

    // deactivate first object if refinement switch is set
    if ( mSettings.Refinement )
    {
        mFieldData.RefineAgainstObjectSDF.reset( 0 );
    }
}

// -------------------------------------------------------------------------------------------------
void
Hierarchical_Mesh_Main::check_and_init_symmetry()
{
    // check if symmetry switch is set
    if ( mSettings.UseSymmetry )
    {

        MORIS_ASSERT( mMeshData.MaxPolynomial == 1, "Symmetry functionality not tested for polynomials > 1" );

        MORIS_ASSERT(  mMeshData.NumberOfElementsPerDirection( mSettings.SymmetryPlane ) % 2 == 0,
                "number of elements in symmetry direction must be even when symmetry feature is used" );

        // cout symmetry axis
        std::fprintf( stdout, "Select Symmetry Plane %u\n", mSettings.SymmetryPlane );
        MORIS_ASSERT( mSettings.SymmetryPlane <= mMeshData.ModelDim,
                      "Selected symmetry plane must be within model dimension");

        // calculate index for element symmetry
        mMeshData.ElementSymmetryIndex
        =  ceil ( (real)  mMeshData.NumberOfElementsPerDirection( mSettings.SymmetryPlane ) / 2 );
        std::fprintf( stdout, "Element Symmetry index  %u\n", mMeshData.ElementSymmetryIndex  );

        // calculate index for basis symmetry
        mMeshData.BasisSymmetryIndex
            =  ceil ( (real)  mMeshData.NumberOfBasisPerDirection( mSettings.SymmetryPlane ) / 2 );

        // calculate offset for symmetry
        //mMeshData.SymmetryOffset = (uint) ( (real) mMeshData.NumberOfBasisPerDirection( mSettings.SymmetryPlane ) ) % 2;

        std::fprintf( stdout, "Basis Symmetry index  %u\n", mMeshData.BasisSymmetryIndex  );

    }
}

// -------------------------------------------------------------------------------------------------

uint
Hierarchical_Mesh_Main::give_symmetry_master_of_basis( const uint& aBasisID ) const
{
    // check if symmetry switch is set
    if ( mSettings.UseSymmetry )
    {

     return mBasis.give_symmetry_master_of_basis(
                aBasisID,
                mMeshData.ModelDim,
                mMeshData.Polynomial,  // FIXME: check if need MaxPolynomial here
                mSettings.SymmetryPlane,
                mMeshData.BasisSymmetryIndex,
                mMeshData.NumberOfElementsPerDirection,
                mMeshData.NumberOfBasisPerDirection );

    }
    else
    {
        return aBasisID;
    }
}
// -------------------------------------------------------------------------------------------------

bool
Hierarchical_Mesh_Main::basis_is_symmetry_master(  const uint& aBasisID ) const
{
    return mBasis.basis_is_symmetry_master(
            aBasisID,
            mMeshData.ModelDim,
            mMeshData.Polynomial,
            mSettings.SymmetryPlane,
            mMeshData.BasisSymmetryIndex,
            mMeshData.NumberOfElementsPerDirection );
}

// -------------------------------------------------------------------------------------------------
uint
Hierarchical_Mesh_Main::give_symmetry_master_of_element( const uint& aElementID ) const
{
    // check if symmetry switch is set
    if ( mSettings.UseSymmetry )
    {
      // calculate position of basis
      Mat< uint > tIJK = give_position_of_element( aElementID );

      // calculate level of basis
      uint tLevel = give_element_level( aElementID );

      // calculate scaling factor for level
      uint tScale = pow( 2, tLevel );

      if ( tIJK( mSettings.SymmetryPlane ) >= tScale*mMeshData.ElementSymmetryIndex )
      {
          // reflect position at defined plane
          tIJK( mSettings.SymmetryPlane ) =  tScale*
                  mMeshData.NumberOfElementsPerDirection( mSettings.SymmetryPlane ) - 1
                   - tIJK( mSettings.SymmetryPlane );

          // return reflected basis
          return give_element_of_position( tLevel, tIJK );
      }
      else
      {
          // ID of this element
          return aElementID;
      }
    }
    else
    {
        return aElementID;
    }
}

// -------------------------------------------------------------------------------------------------
uint
Hierarchical_Mesh_Main::give_symmetry_slave_of_element( const uint& aElementID ) const
{
    // check if symmetry switch is set
    if ( mSettings.UseSymmetry )
    {
      // calculate position of basis
      Mat< uint > tIJK = give_position_of_element( aElementID );

      // calculate level of basis
      uint tLevel = give_element_level( aElementID );

      // calculate scaling factor for level
      uint tScale = pow( 2, tLevel );

      if ( tIJK( mSettings.SymmetryPlane ) < tScale*mMeshData.ElementSymmetryIndex )
      {
          // reflect position at defined plane
          tIJK( mSettings.SymmetryPlane ) =  tScale*
                  mMeshData.NumberOfElementsPerDirection( mSettings.SymmetryPlane ) - 1
                   - tIJK( mSettings.SymmetryPlane );

          // return reflected basis
          return give_element_of_position( tLevel, tIJK );
      }
      else
      {
          // ID of this element
          return aElementID;
      }
    }
    else
    {
        return aElementID;
    }
}

// -------------------------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::make_deactivate_list_symmetric()
{
    // number of candidates
    uint tNumberOfCandidates = mElementData.DeactivateElement.length();

    std::fprintf( stdout, "Elements to check for symmetry condition %u.\n", tNumberOfCandidates ) ;
    if ( tNumberOfCandidates > 0)
    {
        // get max number of elements
        uint tMaxNumberOfElements = mElementData.ElementActive.size();

        // create candidate bitset from deactivate list
        BoostBitset tElementList( tMaxNumberOfElements );

        // loop over all candidates
        for ( uint k=0; k< tNumberOfCandidates; ++k )
        {
            tElementList.set( mElementData.DeactivateElement( k ) );
        }

        // new list of candidates
        Mat< uint > tCandidates( 2*tNumberOfCandidates, 1 );

        // counter for candidates
        uint tCount = 0;

        // loop over all candidates
        for ( uint k=0; k< tNumberOfCandidates; ++k )
        {
            // FIXME: here, we determine master and slave elements of each entry
            //        we are doing a bit too much work here, and this can be done in a more
            //        elegany way
            // assume that this is a slave
            uint tThisElement = mElementData.DeactivateElement( k );

            // get master of element
            uint tMaster = give_symmetry_master_of_element( tThisElement );
            if ( tThisElement == tMaster )
            {
                // if this is a master, get slave
                uint tSlave = give_symmetry_slave_of_element( tThisElement );

                // now check if slave is on list
                if ( ! tElementList.test( tSlave ) )
                {
                    // add master to list
                    tCandidates( tCount++ ) = tSlave;
                }
            }
            else
            {
                // this is a slave

                // now check if master is on list
                if ( ! tElementList.test( tMaster ) )
                {
                    // add master to list
                    tCandidates( tCount++ ) = tMaster;
                }

            }
        }

        if ( tCount > 0)
        {
            // resize list
            tCandidates.resize( tCount, 1);

            // make result unique
            tCandidates = unique( tCandidates );

            tCount = tCandidates.length();

            std::fprintf( stdout, "Elements refined due to symmetry condition %u.\n", tCount ) ;

            // increase deactivate list
            mElementData.DeactivateElement.resize( tNumberOfCandidates + tCount, 1);

            // append candidates to list
            for( uint k=0; k<tCount; ++k)
            {
                mElementData.DeactivateElement( tNumberOfCandidates + k ) = tCandidates( k );
            }
        }
        else
        {
            std::fprintf( stdout, "No element needs to be refined due to symmetry enforcement.\n" );
        }
    }

}

// -------------------------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::check_dof_consistency()
{

    bool verbose = false;

    std::fprintf( stdout, "Checking DOF consistency of %u elements ... ",  ( uint ) mProcData.ElementListOnProc.length() );

    // temporary map
    map< uint, uint > tMap;

    for( uint k=0; k<mBasisData.LagrangeToBSplineMap.length(); ++k )
    {
        tMap[ mBasisData.LagrangeToBSplineMap( k ) ] = k;
    }

    // loop over all elements
    for ( uint e=0; e<mProcData.ElementListOnProc.length(); ++e )
    {

        // Compute number of basis of an element
        uint tNumBasisOfElement = pow( mMeshData.Polynomial + 1, mMeshData.ModelDim );

        // List of Lagrange basis
        Mat< uint > tLagrangeBasis( tNumBasisOfElement, 1);

        // get ID of element
        uint tThisElement = mProcData.ElementListOnProc( e );

        Mat< real > tTMatrix;
        Mat< uint > tMorisIdField;

        // uncomment the next line if you want to print information for a specific element
        //verbose = tThisElement == 438548;

        // calculate T-Matrix for this element
        if( mSettings.TruncatedBsplines )
        {
            mTMatrix.give_Truncated_Tmatrix_and_IdField_Reorder(
                    tThisElement,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection,
                    mBasisData.BSplineActive,
                    tTMatrix,
                    tMorisIdField,
                    mElementData.TMatrixParentChildRelationFEM);
        }
        else
        {
            mTMatrix.give_Tmatrix_and_IdField_Reorder(
                    tThisElement,
                    mMeshData.ModelDim,
                    mMeshData.Polynomial,
                    mMeshData.NumberOfElementsPerDirection,
                    mBasisData.BSplineActive,
                    tTMatrix,
                    tMorisIdField,
                    mElementData.TMatrixParentChildRelationFEM );
        }
        if ( verbose )
        {
            std::cout << "Element ID " << tThisElement << std::endl;
        }

        Mat< uint > tFemdocIdField( tMorisIdField.length(), 1);
        for( uint k=0; k<tMorisIdField.length(); ++k )
        {
            tFemdocIdField( k ) = mBasisData.BSplineMap.find( tMorisIdField( k ) );
        }

        if ( verbose )
        {
            std::cout << "B-Spline List (MORIS - FEMDOC) " << std::endl;
            for( uint k=0; k<tMorisIdField.length(); ++k )
            {
                std::cout << k << " : " << tMorisIdField( k ) << " " <<  tFemdocIdField( k ) <<   std::endl;
            }
            std::cout << std::endl;
            std::cout << "Lagrange List (MORIS - FEMDOC) " << std::endl;
            for( uint k=0; k<tNumBasisOfElement; ++k )
            {
                std::cout << k << " : " << mElementData.FeTopo( e, k ) << " " << mElementData.FeTopoMTK( e, k ) <<   std::endl;
            }
            std::cout << std::endl;
            std::cout << "Weights: " << std::endl;
        }

        // get Lagrange basis of this element
        for( uint k=0; k<tNumBasisOfElement; ++k )
        {

            // index of Lagrange basis in Memory
            uint tFemdocLagrangeIndex = tMap.find( mElementData.FeTopoMTK( e, k ) ) ;
            real tSumOfWeights = 0;
            if ( verbose )
            {
                std::cout << std::endl;
                std::cout << " Node " << mElementData.FeTopoMTK( e, k ) << std::endl;
            }

            for( uint i=0; i< mBasisData.IdFieldNodalField(tFemdocLagrangeIndex , 0 ); ++i )
            {
                // get ID of B-Spline for Femdoc
                uint tFemdocBSplineID = mBasisData.BSplineMap.find( mBasisData.IdFieldNodalField( tFemdocLagrangeIndex, i+1 ) );

                // find index of B-Spline in element
                uint tIndex = UINT_MAX;
                for( uint j=0; j<tFemdocIdField.length(); ++j )
                {
                    if (tFemdocBSplineID == tFemdocIdField( j ) )
                    {
                        tIndex = j;
                        break;
                    }
                }

                // error if id is not found
                if ( tIndex == UINT_MAX )
                {
                    std::cout << "Error: Inconsistent DOF-Connectivity detected in element " << tThisElement << std::endl;
                    std::cout << "       Could not find DOF number " << k << " : " << tFemdocBSplineID << std::endl;

                    tFemdocIdField.print("ID_Field");
                    tTMatrix.print("T_Matrix");

                    exit(1);

                }

                // check that T-Matrix is consistent
                if ( std::abs(  tTMatrix( tIndex, k ) - mBasisData.TMatrixNodalField( tFemdocLagrangeIndex, i+1 ) ) > 1e-6 && mSettings.TruncatedBsplines )
                {
                    std::cout << "Error: Inconsistent DOF-Connectivity weight detected in element " << tThisElement << " at Lagrange Node " << mElementData.FeTopoMTK( e, k ) << std::endl;

                    tFemdocIdField.print("ID Field");
                    tTMatrix.print("T-Matrix");

                    exit(1);
                }
                tSumOfWeights += tTMatrix( tIndex, k );
                if ( verbose )
                {
                    std::cout <<  tFemdocBSplineID << " " << i  << " " << tIndex << " " <<  tTMatrix( tIndex, k  ) << " " << mBasisData.TMatrixNodalField( tFemdocLagrangeIndex, i+1 ) << " " << std::endl;
                }

            }

            // test for unity
            if ( mSettings.Refinement )
            {
                if ( std::abs( tSumOfWeights - 1) > 1E-6 )
                {
                    std::cout << "Error: partition of unity not fulfilled at Lagrange Node  "<< mElementData.FeTopoMTK( e, k )  << std::endl;
                    exit(1);
                }
            }
        }

    }

    std::fprintf( stdout, " passed.\n");
}

// -------------------------------------------------------------------------------------------------

void
Hierarchical_Mesh_Main::activate_symmetry_design_bitset()
{

    // get number of SDFs
    uint tNumberOfSDFs = mFieldData.ObjectSDF.size();

    // test if we are symmetric

    if ( mSettings.UseSymmetry && tNumberOfSDFs > 0 )
    {
        // get number of basis
        /* uint tNumberOfActiveBasis = mFieldData.ObjectSDF( 0 ).length();

        MORIS_ASSERT( tNumberOfActiveBasis == mBasisData.DesignBSplineActive.count(),
                                     "something went wrong while making SDF symmetric: number of active basis does not match." );
        */
        // initialize counters
        uint tCount = 0;

        // get number of
        // loop over all basis
        uint tNumberOfBasis =  mBasisData.DesignBSplineActive.size();

        // reset symmetric field
        mBasisData.DesignBSplineActiveSymmetric.clear();

        // assign memory for symmetric field
        mBasisData.DesignBSplineActiveSymmetric.resize( tNumberOfBasis );

        for( uint k=0; k<tNumberOfBasis; ++k )
        {
            // test if this basis is active
            if (  mBasisData.DesignBSplineActive.test( k ) )
            {
                // test if basis is symmetric
                if (  basis_is_symmetry_master( k ) )
                {

                    // flag this basis
                    mBasisData.DesignBSplineActiveSymmetric.set( k );

                    // increment counter for symmetric basus
                    ++tCount;
                }
            }

        }

        MORIS_ASSERT( mMeshData.NumberOfSymmetryMasterBasis == tCount,
                             "something went wrong while making SDF symmetric" );

    }

}

// -------------------------------------------------------------------------------------------------

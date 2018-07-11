/*
 * cl_Hierarchical_Mesh_TMatrix.cpp
 *
 *  Created on: Dec 7, 2017
 *      Author: gleim
 */

#include "cl_Hierarchical_Mesh_TMatrix.hpp" // STK/src/Hierarchical
using namespace moris;

void
Hierarchical_Mesh_TMatrix::give_Tmatrix_and_IdField_Reorder(
        uint const & aElementId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        BoostBitset const & aBasisActive,
        Mat<real> & aTMatrix,
        Mat<uint> & aIdField,
        Cell<Mat<real>> const & aTMatrixParentChild)
{
    //Number of basis function within an element
    uint tNumberBasisInElement = pow( aPolynomial + 1, aModelDim );
    // Change ordering to the classical order of FE connectivity (for FemDoc and Paraview)
    Mat<uint> tOrder = give_vector_for_reorder( aModelDim, aPolynomial );
    // Level of the element
    uint tLevel = mBaseElement.give_element_level( aElementId, aModelDim, aNumberOfElementsPerDirection ) + 1;
    //Tmatrix of the current element
    Mat<real> tT = eye( tNumberBasisInElement, tNumberBasisInElement );
    //Tmatrix of the parent element
    aTMatrix.set_size( tLevel * tNumberBasisInElement, tNumberBasisInElement, 0 );
    //IdField of the parent Element
    aIdField.set_size( tLevel * tNumberBasisInElement, 1, 0 );
    //Basis of the curent element
    Mat<uint> tBasis = mHMRElement.give_basis_of_element( aElementId, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    //Temporary variable for a loop
    uint tVar = 0;
    // Help matrix to reorder everything in a classical sense (counter clockwise numbering)
    Mat<real> tMatrixDummy;
    for ( uint i = 0; i < tNumberBasisInElement; i++ )
    {
        if ( aBasisActive.test( tBasis( tOrder( i ) ) ) == 1 )
        {
            //Add the row of the Tmatrix of the active basis functions of the curent element
            tMatrixDummy = tT.row( tOrder( i ) );
            for ( uint j = 0; j < tNumberBasisInElement; j++ )
            {
                aTMatrix( tVar, j ) =  tMatrixDummy( tOrder( j ) );
            }
            //Add the row of the IdField of the active basis functions of the current element
            aIdField( tVar ) = tBasis( tOrder( i ) );
            tVar++;
        }
    }
    //Asking for T-matrix on level l-1 (one coarser level)
    tLevel--;
    uint tChildId = aElementId;
    Mat<uint> tParentChildRelation;
    Mat<uint> tChildren( pow( 2, aModelDim ), 1 );
    while(tLevel>0){
        // Who is the parent of the child and which children is it (First entry is parent ID and second entry ist the child number)
        tParentChildRelation = mBaseElement.give_parent_child_realation_of_element( tChildId, aModelDim, aNumberOfElementsPerDirection );
        //Basis of the parent element
        tBasis = mHMRElement.give_basis_of_element( tParentChildRelation(0), aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        //Create T-matrix for the parent element
        tT = aTMatrixParentChild( tParentChildRelation( 1 ) ) * tT;
        for ( uint i = 0; i < tNumberBasisInElement; i++ )
        {
            if ( aBasisActive.test( tBasis( tOrder( i ) ) ) == 1 )
            {
                //Add the row of the Tmatrix of the active basis functions of the curent element
                tMatrixDummy = tT.row( tOrder( i ) );
                for ( uint j = 0; j < tNumberBasisInElement; j++ )
                {
                    aTMatrix( tVar, j ) =  tMatrixDummy( tOrder( j ) );
                }
                //Add the row of the IdField of the active basis functions of the parents
                aIdField( tVar ) = tBasis( tOrder( i ) );
                tVar++;
            }
        }
        tLevel--;
        tChildId = tParentChildRelation( 0 );
    }
    //Resize the Tmatrix and the IdField
    aTMatrix.resize( tVar, tNumberBasisInElement );
    aIdField.resize( tVar, 1 );
}

//-------------------------------------------------------------------------------

void
Hierarchical_Mesh_TMatrix::give_Truncated_Tmatrix_and_IdField_Reorder(
        uint const & aElementId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        BoostBitset const & aBasisActive,
        Mat<real> & aTMatrix,
        Mat<uint> & aIdField,
        Cell<Mat<real>> const & aTMatrixParentChild)
{
    //Number of basis function within an element
    uint tNumberBasisInElement = pow( aPolynomial + 1, aModelDim );
    // Change ordering to the classical order of FE connectivity (for FemDoc and Paraview)
    Mat<uint> tOrder = give_vector_for_reorder( aModelDim, aPolynomial );
    // Level of the element
    uint tLevel = mBaseElement.give_element_level(aElementId,aModelDim,aNumberOfElementsPerDirection)+1;
    //Tmatrix of the curent element
    Mat<real> tT = eye( tNumberBasisInElement, tNumberBasisInElement );
    //Tmatrix of the curent element
    Mat<real> tTeye = eye( tNumberBasisInElement, tNumberBasisInElement );
    //Tmatrix of the parent element
    aTMatrix.set_size( tLevel * tNumberBasisInElement, tNumberBasisInElement, 0 );
    //IdField of the parent Element
    aIdField.set_size( tLevel * tNumberBasisInElement,1 ,0 );
    //Basis of the curent element
    Mat<uint> tBasis = mHMRElement.give_basis_of_element( aElementId, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    //Temporary variable for a loop
    uint tVar = 0;
    uint tVarb = 0;
    uint tVarc = 0;
    uint tVarcell = 0;
    //Stores the tMatrix of a cell of the current child in a matrix
    Mat<real> tTmatrix_child;
    // Stores a active basis functions for each level with the eye matrix
    Cell<Mat<real>> tMatrixActiveBasisEye( tLevel );
    //  Stores a T-matrices for each level
    Cell<Mat<real>> tTMatrixOnLevel( tLevel );
    Mat<real> tMatrixActiveBasis( tLevel * tNumberBasisInElement, tNumberBasisInElement, 0 );
    Mat<real> tTMatrixOfLevel( tLevel * tNumberBasisInElement, tNumberBasisInElement, 0 );
    Mat<real> tTMatrix( tLevel * tNumberBasisInElement, tNumberBasisInElement, 0 );
    Mat<real> tTMatrixDummy( tLevel * tNumberBasisInElement, tNumberBasisInElement, 0 );
    // Help matrix to calculate the truncated T-Matrix
    Mat<real> tTruncatedTMatrixOfLevel;
    // Help matrix to calculate the truncated T-Matrix
    Mat<real> tMatrixActiveBasisDummy;
    // Help matrix to reorder everything in a classical sense (counter clockwise numbering)
    Mat<real> tHelpMata;
    Mat<real> tHelpMatb;
    Mat<real> tHelpMatc;
    for ( uint i = 0; i < tNumberBasisInElement; i++ )
    {
        if ( aBasisActive.test( tBasis( tOrder( i ) ) ) == 1 )
        {
            //Add the row of the Tmatrix of the active basis functions of the curent element
            tHelpMata = tT.row( tOrder( i ) );
            for ( uint j = 0; j < tNumberBasisInElement; j++ )
            {
                aTMatrix( tVar, j ) =  tHelpMata( tOrder( j ) );
                //Add the row of the Tmatrix of the active basis functions of the curent element
                tMatrixActiveBasis( tVar, j ) = tHelpMata( tOrder( j ) );
                //Add the row of the Tmatrix of the active basis functions of the curent element
                tTMatrixDummy( tVar, j ) = tHelpMata( tOrder( j ) );
            }
            //Add the row of the IdField of the active basis functions of the curent element
            aIdField( tVar ) = tBasis( tOrder( i ) );
            tVar++;
        }
    }
    // Save curent size of aIdField
    tVarb = tVar;
    // Save curent size of aTMatrix
    tVarc = tVar;
    tTMatrixDummy.resize( tVar, tNumberBasisInElement );
    tMatrixActiveBasis.resize( tVar, tNumberBasisInElement );
    // Save active basis functions with a one for their position
    tMatrixActiveBasisEye(tVarcell) = tMatrixActiveBasis;
    // Save the TMatrix of the current level
    tTMatrixOnLevel(tVarcell) = tTMatrixDummy;
    tVarcell++;
    //Asking for T-matrix on level l-1 (one coarser level)
    tLevel--;
    uint tChildId = aElementId;
    Mat<uint> tParentChildRelation;
    Mat<uint> tChildren( pow( 2, aModelDim ), 1 );
    while(tLevel>0){
        // Who is the parent of the child and which children is it (First entry is parent ID and second entry ist the child number)
        tParentChildRelation = mBaseElement.give_parent_child_realation_of_element( tChildId, aModelDim, aNumberOfElementsPerDirection );
        //Basis of the parent element
        tBasis = mHMRElement.give_basis_of_element( tParentChildRelation( 0 ), aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        // Which parent child relation is it
        tTmatrix_child = aTMatrixParentChild( tParentChildRelation( 1 ) );
        //Create T-matrix for the parent element
        tT=tTmatrix_child*tT;
        tVar = 0;
        tTMatrix.set_size( tNumberBasisInElement, tNumberBasisInElement, 0 );
        tTMatrixOfLevel.set_size( tNumberBasisInElement, tNumberBasisInElement, 0 );
        tMatrixActiveBasis.set_size( tNumberBasisInElement, tNumberBasisInElement, 0 );
        for ( uint i = 0; i < tNumberBasisInElement; i++ )
        {
            if ( aBasisActive.test( tBasis( tOrder( i ) ) ) == 1 )
            {
                //Add the row of the Tmatrix of the active basis functions of the parents
                tHelpMata = tT.row( tOrder( i ) );
                //Add the row of the Tmatrix of the active basis functions of the parents
                tHelpMatb = tTmatrix_child.row( tOrder( i ) );
                // Add active basis functions with a one for their position
                tHelpMatc = tTeye.row( tOrder( i ) );
                for ( uint j = 0; j < tNumberBasisInElement; j++ )
                {
                    //Add the row of the Tmatrix of the active basis functions of the parents
                    tTMatrix( tVar, j ) = tHelpMata( tOrder( j ) );
                    //Add the row of the Tmatrix of the active basis functions of the parents
                    tTMatrixOfLevel( tVar, j ) = tHelpMatb( tOrder( j ) );
                    // Add active basis functions with a one for their position
                    tMatrixActiveBasis( tVar, j ) = tHelpMatc( tOrder( j ) );
                }
                //Add the row of the IdField of the active basis functions of the parents
                aIdField( tVarb ) = tBasis( tOrder( i ) );
                tVar++;
                tVarb++;
            }
        }
        tLevel--;
        tChildId = tParentChildRelation( 0 );
        tTMatrix.resize( tVar, tNumberBasisInElement );
        tMatrixActiveBasis.resize( tVar, tNumberBasisInElement );
        tTMatrixOfLevel.resize( tVar, tNumberBasisInElement );
        tMatrixActiveBasisDummy = tMatrixActiveBasisEye( tVarcell - 1 );
        tTMatrixDummy = tTMatrixOnLevel( tVarcell - 1 );
        // Take the T-Matrix of the current level and subrtract parts from the children level, which have an influence
        tTruncatedTMatrixOfLevel = tTMatrix - ( ( tTMatrixOfLevel * trans( tMatrixActiveBasisDummy ) ) * tTMatrixDummy );
        tTMatrixOnLevel( tVarcell ) = tTMatrix;
        tMatrixActiveBasisEye( tVarcell ) = tMatrixActiveBasis;
        tVarcell++;
        for ( uint i = 0; i < tVarcell - 2; i++ )
        {
            tMatrixActiveBasis = tMatrixActiveBasisEye( i );
            tTMatrix = tTMatrixOnLevel( i );
            // Take the T-Matrix of the current level and subtract parts from all finer levels, which have an influence
            tTruncatedTMatrixOfLevel = tTruncatedTMatrixOfLevel - (tTruncatedTMatrixOfLevel * trans( tMatrixActiveBasis ) ) * tMatrixActiveBasis;
        }
        if ( tVar > 0 )
            aTMatrix.rows( tVarc, tVarb - 1 ) = tTruncatedTMatrixOfLevel.rows( 0, tVar - 1 );
        tVarc = tVarb;
    }
    //Resize the Tmatrix and the IdField
    aTMatrix.resize( tVarb, tNumberBasisInElement);
    aIdField.resize( tVarb, 1);
}

//-------------------------------------------------------------------------------

void
Hierarchical_Mesh_TMatrix::give_IdField(
        uint const        & aElementId,
        uint const        & aModelDim,
        uint const        & aPolynomial,
        Mat<uint> const   & aNumberOfElementsPerDirection,
        BoostBitset const & aBasisActive,
        Mat<uint>         & aIdField)
{
    //Number of basis function within an element
    uint tNumberBasisInElement = pow( aPolynomial + 1, aModelDim );
    // Level of the element
    uint tLevel = mBaseElement.give_element_level( aElementId, aModelDim, aNumberOfElementsPerDirection ) + 1;
    //IdField of the parent Element
    aIdField.set_size( tLevel * tNumberBasisInElement, 1, 0 );
    //Basis of the current element
    Mat<uint> tBasis = mHMRElement.give_basis_of_element( aElementId, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    //Temporary variable for a loop
    uint tVar = 0;
    for ( uint i = 0; i < tNumberBasisInElement; i++ )
    {
        if ( aBasisActive.test( tBasis( i ) ) == 1 )
        {
            //Add the row of the IdField of the active basis functions of the curent element
            aIdField(tVar) = tBasis( i );
            tVar++;
        }
    }
    //Asking for T-matrix on level l-1 (one coarser level)
    tLevel--;
    uint tChildId=aElementId;
    uint tParentId;
    while( tLevel > 0 )
    {
        // Who is the parent of the child and which children is it
        tParentId = mBaseElement.give_parent_of_element( tChildId, aModelDim, aNumberOfElementsPerDirection );
        //Basis of the parent element
        tBasis = mHMRElement.give_basis_of_element( tParentId, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        for ( uint i = 0; i < tNumberBasisInElement; i++ )
        {
            if ( aBasisActive.test( tBasis( i ) ) == 1 )
            {
                //Add the row of the IdField of the active basis functions of the parents
                aIdField( tVar ) = tBasis( i );
                tVar++;
            }
        }
        tLevel--;
        tChildId=tParentId;
    }
    //Resize the IdField
    aIdField.resize( tVar, 1 );
}

//-------------------------------------------------------------------------------

void
Hierarchical_Mesh_TMatrix::give_Tmatrix_and_IdField(
        uint const & aElementId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        BoostBitset const & aBasisActive,
        Mat<real> & aTMatrix,
        Mat<uint> & aIdField,
        Cell<Mat<real>> const & aTMatrixParentChild)
{
    //Number of basis function within an element
    uint tNumberBasisInElement = pow( aPolynomial + 1, aModelDim );
    // Level of the element
    uint tLevel = mBaseElement.give_element_level( aElementId, aModelDim, aNumberOfElementsPerDirection ) + 1;
    //Tmatrix of the current element
    Mat<real> tTdummy = eye( tNumberBasisInElement, tNumberBasisInElement );
    //Tmatrix of the parent element
    aTMatrix.set_size( tLevel * tNumberBasisInElement, tNumberBasisInElement, 0 );
    //IdField of the parent Element
    aIdField.set_size( tLevel * tNumberBasisInElement, 1, 0 );
    //Basis of the current element
    Mat<uint> tBasis = mHMRElement.give_basis_of_element( aElementId, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    uint tVar = 0; //Temporary variable for a loop

    for ( uint i = 0; i < tNumberBasisInElement; i++ )
    {
        if ( aBasisActive.test( tBasis( i ) ) == 1 )
        {
            //Add the row of the Tmatrix of the active basis functions of the current element
            aTMatrix.row( tVar ) = tTdummy.row( i );
            //Add the row of the IdField of the active basis functions of the current element
            aIdField( tVar ) = tBasis( i );
            tVar++;
        }
    }
    //Asking for T-matrix on level l-1 (one coarser level)
    tLevel--;
    uint tChildId=aElementId;
    Mat<uint> tParentChildRelation;
    Mat<uint> tChildren( pow( 2, aModelDim ), 1 );
    while( tLevel > 0 )
    {
        // Who is the parent of the child and which children is it (First entry is parent ID and second entry ist the child number)
        tParentChildRelation = mBaseElement.give_parent_child_realation_of_element( tChildId, aModelDim, aNumberOfElementsPerDirection );
        //Basis of the parent element
        tBasis = mHMRElement.give_basis_of_element( tParentChildRelation( 0 ), aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        //Create T-matrix for the parent element
        tTdummy = aTMatrixParentChild( tParentChildRelation( 1 ) ) * tTdummy;
        for ( uint i = 0; i < tNumberBasisInElement; i++ )
        {
            if ( aBasisActive.test( tBasis( i ) ) == 1 )
            {
                //Add the row of the Tmatrix of the active basis functions of the curent element
                aTMatrix.row( tVar ) = tTdummy.row( i );
                //Add the row of the IdField of the active basis functions of the parents
                aIdField( tVar ) = tBasis( i );
                tVar++;
            }
        }
        tLevel--;
        tChildId = tParentChildRelation( 0 );
    }
    //Resize the Tmatrix and the IdField
    aTMatrix.resize( tVar, tNumberBasisInElement );
    aIdField.resize( tVar, 1 );
}

//-------------------------------------------------------------------------------

void
Hierarchical_Mesh_TMatrix::give_Truncated_Tmatrix_and_IdField(
        uint const & aElementId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        BoostBitset const & aBasisActive,
        Mat<real> & aTMatrix,
        Mat<uint> & aIdField,
        Cell<Mat<real>> const & aTMatrixParentChild)
{
    //Number of basis function within an element
    uint tNumberBasisInElement = pow( aPolynomial + 1, aModelDim );
    // Level of the element
    uint tLevel = mBaseElement.give_element_level( aElementId, aModelDim, aNumberOfElementsPerDirection ) + 1;
    //Tmatrix of the curent element
    Mat<real> tTdummy = eye( tNumberBasisInElement, tNumberBasisInElement );
    //Tmatrix for reduction of dependent parts from upper levels
    Mat<real> tTeye = eye( tNumberBasisInElement, tNumberBasisInElement );
    //Tmatrix of the parent element
    aTMatrix.set_size( tLevel * tNumberBasisInElement, tNumberBasisInElement, 0 );
    //IdField of the parent Element
    aIdField.set_size( tLevel * tNumberBasisInElement, 1, 0 );
    //Basis of the curent element
    Mat<uint> tBasis = mHMRElement.give_basis_of_element( aElementId, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    //Temporary variable for a loop
    uint tVar = 0;
    uint tVarb = 0;
    uint tVarc = 0;
    uint tVarcell = 0;
    //Stores the tMatrix of a cell of the current child in a matrix
    Mat<real> tTmatrix_child;
    // Stores a active basis functions for each level with the eye matrix
    Cell<Mat<real>> tMatrixActiveBasisEye( tLevel );
    //  Stores a T-matrices for each level
    Cell<Mat<real>> tTMatrixOnLevel( tLevel );
    Mat<real> tMatrixActiveBasis( tLevel * tNumberBasisInElement, tNumberBasisInElement, 0 );
    Mat<real> tTMatrixOfLevel( tLevel * tNumberBasisInElement, tNumberBasisInElement, 0 );
    Mat<real> tTMatrix( tLevel * tNumberBasisInElement, tNumberBasisInElement, 0 );
    Mat<real> tTMatrixDummy( tLevel * tNumberBasisInElement, tNumberBasisInElement, 0 );
    // Help matrix to calculate the truncated T-Matrix
    Mat<real> tTruncatedTMatrixOfLevel;
    Mat<real> tMatrixActiveBasisDummy;
    for (uint i = 0; i<tNumberBasisInElement; i++)
    {
        if ( aBasisActive.test( tBasis( i ) ) == 1 )
        {
            //Add the row of the Tmatrix of the active basis functions of the current element
            aTMatrix.row( tVar ) = tTdummy.row( i );
            //Add the row of the Tmatrix of the active basis functions of the current element
            tMatrixActiveBasis.row( tVar ) = tTdummy.row( i );
            //Add the row of the Tmatrix of the active basis functions of the current element
            tTMatrixDummy.row( tVar ) = tTdummy.row( i );
            //Add the row of the IdField of the active basis functions of the current element
            aIdField( tVar ) = tBasis( i );
            tVar++;
        }
    }
    // Save curent size of aIdField
    tVarb = tVar;
    tVarc = tVar;
    tTMatrixDummy.resize( tVar, tNumberBasisInElement );
    tMatrixActiveBasis.resize( tVar, tNumberBasisInElement );
    // Save active basis functions with a one for their position
    tMatrixActiveBasisEye( tVarcell ) = tMatrixActiveBasis;
    // Save the TMatrix of the current level
    tTMatrixOnLevel( tVarcell ) = tTMatrixDummy;
    tVarcell++;
    tLevel--; //Asking for T-matrix on level l-1 (one coarser level)
    uint tChildId = aElementId;
    Mat<uint> tParentChildRelation;
    Mat<uint> tChildren( pow( 2, aModelDim ), 1 );
    while( tLevel > 0 )
    {
        // Who is the parent of the child and which children is it (First entry is parent ID and second entry ist the child number)
        tParentChildRelation = mBaseElement.give_parent_child_realation_of_element( tChildId, aModelDim, aNumberOfElementsPerDirection );
        //Basis of the parent element
        tBasis = mHMRElement.give_basis_of_element( tParentChildRelation( 0 ), aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        tTmatrix_child = aTMatrixParentChild( tParentChildRelation( 1 ) );
        //Create T-matrix for the parent element
        tTdummy = tTmatrix_child * tTdummy;
        tVar = 0;
        tTMatrix.set_size( tNumberBasisInElement, tNumberBasisInElement, 0 );
        tTMatrixOfLevel.set_size( tNumberBasisInElement, tNumberBasisInElement, 0 );
        tMatrixActiveBasis.set_size( tNumberBasisInElement, tNumberBasisInElement, 0 );
        for ( uint i = 0; i < tNumberBasisInElement; i++ )
        {
            if ( aBasisActive.test( tBasis( i ) ) == 1 )
            {
                //Add the row of the Tmatrix of the active basis functions of the parents
                tTMatrix.row( tVar ) = tTdummy.row( i );
                //Add the row of the IdField of the active basis functions of the parents
                aIdField( tVarb ) = tBasis( i );
                //Add the row of the Tmatrix of the active basis functions of the parents
                tTMatrixOfLevel.row( tVar ) = tTmatrix_child.row( i );
                // Add active basis functions with a one for their position
                tMatrixActiveBasis.row( tVar ) = tTeye.row( i );
                tVar++;
                tVarb++;
            }
        }
        tLevel--;
        tChildId = tParentChildRelation( 0 );
        tTMatrix.resize( tVar, tNumberBasisInElement );
        tMatrixActiveBasis.resize( tVar, tNumberBasisInElement );
        tTMatrixOfLevel.resize( tVar, tNumberBasisInElement );
        tMatrixActiveBasisDummy = tMatrixActiveBasisEye( tVarcell - 1 );
        tTMatrixDummy = tTMatrixOnLevel( tVarcell - 1 );
        // Take the T-Matrix of the current level and subrtract parts from the children level, which have an influence
        tTruncatedTMatrixOfLevel = tTMatrix - ( ( tTMatrixOfLevel * trans( tMatrixActiveBasisDummy ) ) * tTMatrixDummy );
        tTMatrixOnLevel( tVarcell ) = tTMatrix;
        tMatrixActiveBasisEye( tVarcell ) = tMatrixActiveBasis;
        tVarcell++;
        for ( uint i = 0; i < tVarcell - 2; i++ )
        {
            tMatrixActiveBasis = tMatrixActiveBasisEye( i );
            // Take the T-Matrix of the current level and subtract parts from all finer levels, which have an influence
            tTruncatedTMatrixOfLevel = tTruncatedTMatrixOfLevel - ( tTruncatedTMatrixOfLevel * trans( tMatrixActiveBasis ) ) * tMatrixActiveBasis;
        }
        if ( tVar > 0 )
            aTMatrix.rows( tVarc, tVarb - 1 ) = tTruncatedTMatrixOfLevel.rows( 0, tVar - 1 );
        tVarc = tVarb;
    }
    //Resize the Tmatrix and the IdField
    aTMatrix.resize( tVarb, tNumberBasisInElement );
    aIdField.resize( tVarb, 1 );
}

//-------------------------------------------------------------------------------

void
Hierarchical_Mesh_TMatrix::give_Tmatrix_and_IdField_DesignToFEMProjection(
        uint const & aElementId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        BoostBitset const & aBasisActive,
        Mat<real> & aTMatrix,
        Mat<uint> & aIdField,
        Cell<Mat<real>> const & aTMatrixParentChild)
{
    this->give_Tmatrix_and_IdField( aElementId, aModelDim, aPolynomial, aNumberOfElementsPerDirection, aBasisActive, aTMatrix, aIdField, aTMatrixParentChild );
    //Projection to the FEM mesh
    Mat<real> tT_Project =  give_projection_matrix( aModelDim, aPolynomial );
    aTMatrix = aTMatrix * tT_Project;
}

//-------------------------------------------------------------------------------

void
Hierarchical_Mesh_TMatrix::give_Truncated_Tmatrix_and_IdField_DesignToFEMProjection(
        uint const & aElementId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        BoostBitset const & aBasisActive,
        Mat<real> & aTMatrix,
        Mat<uint> & aIdField,
        Cell<Mat<real>> const & aTMatrixParentChild)
{
    this->give_Truncated_Tmatrix_and_IdField( aElementId, aModelDim, aPolynomial, aNumberOfElementsPerDirection, aBasisActive, aTMatrix, aIdField, aTMatrixParentChild );
    //Projection to the FEM mesh
    Mat<real> tT_Project =  give_projection_matrix( aModelDim, aPolynomial );
    // Project one field to another one ( Design to FEM )
    aTMatrix = aTMatrix * tT_Project;
}

//-------------------------------------------------------------------------------

void
Hierarchical_Mesh_TMatrix::give_Tmatrix_and_IdField_SpecificPoint(
        uint const & aElementId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        BoostBitset const & aBasisActive,
        Mat<real> & aTMatrix,
        Mat<uint> & aIdField,
        Mat<real> const & aNaturalCoordinate,
        Cell<Mat<real>> const & aTMatrixParentChild)
{
    this->give_Tmatrix_and_IdField( aElementId, aModelDim, aPolynomial, aNumberOfElementsPerDirection, aBasisActive, aTMatrix, aIdField, aTMatrixParentChild );
    //Calculate the Shape functions for the specific point
    Mat<real> tT_Project;
    //Create temporary TMatrix for calculation of a list of natural coordinates
    Mat<real> tTMatrixTemp( aTMatrix.n_rows(), aNaturalCoordinate.n_rows() );
    for ( uint i = 0; i < aNaturalCoordinate.n_rows(); i++ )
    {
        if ( aModelDim == 1 )
        {
            tT_Project = Bspline::build_spline_uniform_1d( aNaturalCoordinate( i, 0 ), aPolynomial );
        }
        else
        {
            tT_Project = Bspline::build_spline_uniform_nd( aNaturalCoordinate.row( i ), aPolynomial, aModelDim );
        }
        tTMatrixTemp.col( i ) = aTMatrix * tT_Project;
    }
    // Save new matrix to argument
    aTMatrix = tTMatrixTemp;
}

//-------------------------------------------------------------------------------

void
Hierarchical_Mesh_TMatrix::give_Truncated_Tmatrix_and_IdField_SpecificPoint(
        uint const & aElementNumber,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        BoostBitset const & aBasisActive,
        Mat<real> & aTMatrix,
        Mat<uint> & aIdField,
        Mat<real> const & aNaturalCoordinate,
        Cell<Mat<real>> const & aTMatrixParentChild)
{
    this->give_Truncated_Tmatrix_and_IdField( aElementNumber, aModelDim, aPolynomial, aNumberOfElementsPerDirection, aBasisActive, aTMatrix, aIdField, aTMatrixParentChild );
    //Calculate the Shape functions for the specific point
    Mat<real> tT_Project;
    //Create temporary TMatrix for calculation of a list of natural coordinates
    Mat<real> tTMatrixTemp( aTMatrix.n_rows(), aNaturalCoordinate.n_rows() );
    for ( uint i = 0; i < aNaturalCoordinate.n_rows(); i++ )
    {
        if ( aModelDim == 1 )
        {
            tT_Project = Bspline::build_spline_uniform_1d( aNaturalCoordinate( i, 0 ), aPolynomial );
        }
        else
        {
            tT_Project = Bspline::build_spline_uniform_nd( aNaturalCoordinate.row( i ), aPolynomial, aModelDim );
        }
        tTMatrixTemp.col( i ) = aTMatrix * tT_Project;
    }
    // Save new matrix to argument
    aTMatrix = tTMatrixTemp;
}

//-------------------------------------------------------------------------------

void
Hierarchical_Mesh_TMatrix::give_Tmatrix_and_IdField_Design_L2Projection_Coarsening(
        uint const & aElementId,
        uint const & aModelDim,
        uint const & aPolynomial,
        uint const & aLevelFem,
        uint const & aLevelDesign,
        Mat<uint> const & aNumberOfElementsPerDirection,
        BoostBitset const & aElementActiveLastStep,
        Cell<Mat<real>> & aTMatrixOfChildren,
        Mat<uint> & aListOfChildren,
        Cell<Mat<real>> const & aTMatrixParentChild)
{
    //Number of basis function within an element
    uint tNumberBasisInElement = pow( aPolynomial + 1, aModelDim );
    //Number of possible elements
    uint tPowPowElements = pow(pow(2,aModelDim),(aLevelFem)+1);
    //Tmatrix of the current element
    Mat<real> tT = eye( tNumberBasisInElement, tNumberBasisInElement );
    Cell<Mat<real>> tTforchildren( tPowPowElements );
    Mat<uint> tPossibleChildrenList( tPowPowElements, 1, 0 );
    Mat<uint> tActiveChildrenList( tPowPowElements, 1, 0 );
    Mat<uint> tListActiveChildren( tPowPowElements, 1, 0 );
    //Temporary variable for a loop
    uint tVara = 0;
    uint tVarb = 1;
    uint tVarc = 1;
    uint tVard = 0;
    Mat<uint> tChildren( pow( 2, aModelDim ), 1 );
    tPossibleChildrenList( 0 ) = aElementId;
    tTforchildren( 0 ) = tT;
    uint tLevel = mBaseElement.give_element_level( aElementId, aModelDim, aNumberOfElementsPerDirection );
    while( tPossibleChildrenList( tVara ) > 0 && tLevel < aLevelDesign )
    {
        // Give the children of the parent element
        tChildren = mBaseElement.give_children_of_element( tPossibleChildrenList( tVara ), aModelDim, aNumberOfElementsPerDirection );
        for ( uint i = 0; i < tChildren.length(); i++ )
        {
            tLevel = mBaseElement.give_element_level( tChildren( i ), aModelDim, aNumberOfElementsPerDirection );
            Mat<real> tTmatrix_child = aTMatrixParentChild( i );
            tT = tTforchildren( tVara );
            //Create T-matrix for the parent element
            tT = tTmatrix_child * tT;
            tTforchildren( tVarc ) = tT;
            if ( aElementActiveLastStep.test( tChildren( i ) ) == 1 )
            {
                tActiveChildrenList( tVard ) = tVarc;
                tListActiveChildren( tVard ) = tChildren( i );
                tVard++;
            }
            else if ( tLevel < aLevelDesign )
            {
                tPossibleChildrenList( tVarb ) =  tChildren( i );
                tVarb++;
            }
            tVarc++;
        }
        tVara++;
    }
    tActiveChildrenList.resize( tVard, 1 );
    tListActiveChildren.resize( tVard, 1 );
    Cell<Mat<real>> tTMatrixActiveChildren( tVard );
    for ( uint i = 0; i < tVard; i++ )
    {
        tTMatrixActiveChildren( i ) = tTforchildren( tActiveChildrenList( i ) );
    }
    aListOfChildren = tListActiveChildren;
    aTMatrixOfChildren = tTMatrixActiveChildren;
}

//-------------------------------------------------------------------------------

Mat<real>
Hierarchical_Mesh_TMatrix::give_projection_matrix(
        uint const & aModelDim,
        uint const & aPolynomial)
{
    // Projection matrix
    Mat<real> tT_Project;
    if ( aModelDim == 2 )
    {
        if ( aPolynomial == 1 )
        {
            Mat<real> tTtemp = {{1, 0, 0, 0},{0, 1, 0 ,0},{0, 0, 1, 0},{0, 0, 0, 1}};
            tT_Project = tTtemp;
        }
        else if ( aPolynomial == 2 )
        {
            Mat<real> tTtemp = {
                    {0.2500,         0,         0,         0},
                    {0.2500,    0.2500,         0,         0},
                    {     0,    0.2500,         0,         0},
                    {0.2500,         0,    0.2500,         0},
                    {0.2500,    0.2500,    0.2500,    0.2500},
                    {     0,    0.2500,         0,    0.2500},
                    {     0,         0,    0.2500,         0},
                    {     0,         0,    0.2500,    0.2500},
                    {     0,         0,         0,    0.2500}};
            tT_Project = tTtemp;
        }
        else if ( aPolynomial == 3 )
        {
            Mat<real> tTtemp = {
                    {0.0278,         0,         0,         0},{0.1111,    0.0278,         0,         0},{0.0278,    0.1111,         0,         0},{     0,    0.0278,         0,         0},{0.1111,         0,    0.0278,         0},
                    {0.4444,    0.1111,    0.1111,    0.0278},{0.1111,    0.4444,    0.0278,    0.1111},{     0,    0.1111,         0,    0.0278},{0.0278,         0,    0.1111,         0},{0.1111,    0.0278,    0.4444,    0.1111},
                    {0.0278,    0.1111,    0.1111,    0.4444},{     0,    0.0278,         0,    0.1111},{     0,         0,    0.0278,         0},{     0,         0,    0.1111,    0.0278},{     0,         0,    0.0278,    0.1111},{     0,         0,         0,    0.0278}};
            tT_Project = tTtemp;
        }
        else
        {
            MORIS_LOG_ERROR << " Polynomial degree is not implemented";
        }
    }
    else if ( aModelDim == 3 )
    {
        if ( aPolynomial == 1 )
        {
            Mat<real> tTtemp = {
                    { 1,     0,     0,     0,     0,     0,     0,     0},
                    {0,     1,     0,     0,     0,     0,     0,     0},
                    {0,     0,     1,     0,     0,     0,     0,     0},
                    {0,     0,     0,     1,     0,     0,     0,     0},
                    {0,     0,     0,     0,     1,     0,     0,     0},
                    {0,     0,     0,     0,     0,     1,     0,     0},
                    {0,     0,     0,     0,     0,     0,     1,     0},
                    {0,     0,     0,     0,     0,     0,     0,     1}};
            tT_Project = tTtemp;
        }
        else if ( aPolynomial == 2 )
        {
            Mat<real> tTtemp = {
                    {0.1250,         0,         0,         0,         0,         0,         0,         0},{0.1250,    0.1250,         0,         0,         0,         0,         0,         0},{     0,    0.1250,         0,         0,         0,         0,         0,         0},
                    {0.1250,         0,    0.1250,         0,         0,         0,         0,         0},{0.1250,    0.1250,    0.1250,    0.1250,         0,         0,         0,         0},{     0,    0.1250,         0,    0.1250,         0,         0,         0,         0},
                    {     0,         0,    0.1250,         0,         0,         0,         0,         0},{     0,         0,    0.1250,    0.1250,         0,         0,         0,         0},{     0,         0,         0,    0.1250,         0,         0,         0,         0},
                    {0.1250,         0,         0,         0,    0.1250,         0,         0,         0},{0.1250,    0.1250,         0,         0,    0.1250,    0.1250,         0,         0},{     0,    0.1250,         0,         0,         0,    0.1250,         0,         0},
                    {0.1250,         0,    0.1250,         0,    0.1250,         0,    0.1250,         0},{0.1250,    0.1250,    0.1250,    0.1250,    0.1250,    0.1250,    0.1250,    0.1250},{     0,    0.1250,         0,    0.1250,         0,    0.1250,         0,    0.1250},{     0,         0,    0.1250,         0,         0,         0,    0.1250,         0},
                    {     0,         0,    0.1250,    0.1250,         0,         0,    0.1250,    0.1250},{     0,         0,         0,    0.1250,         0,         0,         0,    0.1250},{     0,         0,         0,         0,    0.1250,         0,         0,         0},{     0,         0,         0,         0,    0.1250,    0.1250,         0,         0},
                    {     0,         0,         0,         0,         0,    0.1250,         0,         0},{     0,         0,         0,         0,    0.1250,         0,    0.1250,         0},{     0,         0,         0,         0,    0.1250,    0.1250,    0.1250,    0.1250},{     0,         0,         0,         0,         0,    0.1250,         0,    0.1250},
                    {     0,         0,         0,         0,         0,         0,    0.1250,         0},{     0,         0,         0,         0,         0,         0,    0.1250,    0.1250},{     0,         0,         0,         0,         0,         0,         0,    0.1250}};
            tT_Project = tTtemp;
        }
        else if ( aPolynomial == 3 )
        {
            Mat<real> tTtemp = {
                    {0.004629629629630,                   0,                   0,                   0,                   0,                   0,                   0,                   0 },{0.018518518518519,   0.004629629629630,                   0,                   0,                   0,                   0,                   0,                   0},{0.004629629629630,   0.018518518518519,                   0,                   0,                   0,                   0,                   0,                   0},
                    {                0,   0.004629629629630,                   0,                   0,                   0,                   0,                   0,                   0},{0.018518518518519,                   0,   0.004629629629630,                   0,                   0,                   0,                   0,                   0},{0.074074074074074,   0.018518518518519,   0.018518518518519,   0.004629629629630,                   0,                   0,                   0,                   0},{0.018518518518519,   0.074074074074074,   0.004629629629630,   0.018518518518519,                   0,                   0,                   0,                   0},
                    {                0,   0.018518518518519,                   0,   0.004629629629630,                   0,                   0,                   0,                   0},{0.004629629629630,                   0,   0.018518518518519,                   0,                   0,                   0,                   0,                   0},{0.018518518518519,   0.004629629629630,   0.074074074074074,   0.018518518518519,                   0,                   0,                   0,                   0},{0.004629629629630,   0.018518518518519,   0.018518518518519,   0.074074074074074,                   0,                   0,                   0,                   0},
                    {                0,   0.004629629629630,                   0,   0.018518518518519,                   0,                   0,                   0,                   0},{                0,                   0,   0.004629629629630,                   0,                   0,                   0,                   0,                   0},{                0,                   0,   0.018518518518519,   0.004629629629630,                   0,                   0,                   0,                   0},{                0,                   0,   0.004629629629630,   0.018518518518519,                   0,                   0,                   0,                   0},
                    {                0,                   0,                   0,   0.004629629629630,                   0,                   0,                   0,                   0},{0.018518518518519,                   0,                   0,                   0,   0.004629629629630,                   0,                   0,                   0},{0.074074074074074,   0.018518518518519,                   0,                   0,   0.018518518518519,   0.004629629629630,                   0,                   0},{0.018518518518519,   0.074074074074074,                   0,                   0,   0.004629629629630,   0.018518518518519,                   0,                   0},
                    {                0,   0.018518518518519,                   0,                   0,                   0,   0.004629629629630,                   0,                   0},{0.074074074074074,                   0,   0.018518518518519,                   0,   0.018518518518519,                   0,   0.004629629629630,                   0},{0.296296296296296,   0.074074074074074,   0.074074074074074,   0.018518518518519,   0.074074074074074,   0.018518518518519,   0.018518518518519,   0.004629629629630},{0.074074074074074,   0.296296296296296,   0.018518518518519,   0.074074074074074,   0.018518518518519,   0.074074074074074,   0.004629629629630,   0.018518518518519},
                    {                0,   0.074074074074074,                   0,   0.018518518518519,                   0,   0.018518518518519,                   0,   0.004629629629630},{0.018518518518519,                   0,   0.074074074074074,                   0,   0.004629629629630,                   0,   0.018518518518519,                   0},{0.074074074074074,   0.018518518518519,   0.296296296296296,   0.074074074074074,   0.018518518518519,   0.004629629629630,   0.074074074074074,   0.018518518518519},
                    {0.018518518518519,   0.074074074074074,   0.074074074074074,   0.296296296296296,   0.004629629629630,   0.018518518518519,   0.018518518518519,   0.074074074074074},{                0,   0.018518518518519,                   0,   0.074074074074074,                   0,   0.004629629629630,                   0,   0.018518518518519},{                0,                   0,   0.018518518518519,                   0,                   0,                   0,   0.004629629629630,                   0},
                    {                0,                   0,   0.074074074074074,   0.018518518518519,                   0,                   0,   0.018518518518519,   0.004629629629630},{                0,                   0,   0.018518518518519,   0.074074074074074,                   0,                   0,   0.004629629629630,   0.018518518518519},{                0,                   0,                   0,   0.018518518518519,                   0,                   0,                   0,   0.004629629629630},{0.004629629629630,                   0,                   0,                   0,   0.018518518518519,                   0,                   0,                   0},
                    {0.018518518518519,   0.004629629629630,                   0,                   0,   0.074074074074074,   0.018518518518519,                   0,                   0},{0.004629629629630,   0.018518518518519,                   0,                   0,   0.018518518518519,   0.074074074074074,                   0,                   0},{                0,   0.004629629629630,                   0,                   0,                   0,   0.018518518518519,                   0,                   0},{0.018518518518519,                   0,   0.004629629629630,                   0,   0.074074074074074,                   0,   0.018518518518519,                   0},
                    {0.074074074074074,   0.018518518518519,   0.018518518518519,   0.004629629629630,   0.296296296296296,   0.074074074074074,   0.074074074074074,   0.018518518518519},{0.018518518518519,   0.074074074074074,   0.004629629629630,   0.018518518518519,   0.074074074074074,   0.296296296296296,   0.018518518518519,   0.074074074074074},{                0,   0.018518518518519,                   0,   0.004629629629630,                   0,   0.074074074074074,                   0,   0.018518518518519},{0.004629629629630,                   0,   0.018518518518519,                   0,   0.018518518518519,                   0,   0.074074074074074,                   0},
                    {0.018518518518519,   0.004629629629630,   0.074074074074074,   0.018518518518519,   0.074074074074074,   0.018518518518519,   0.296296296296296,   0.074074074074074},{0.004629629629630,   0.018518518518519,   0.018518518518519,   0.074074074074074,   0.018518518518519,   0.074074074074074,   0.074074074074074,   0.296296296296296},{                0,   0.004629629629630,                   0,   0.018518518518519,                   0,   0.018518518518519,                   0,   0.074074074074074},{                0,                   0,   0.004629629629630,                   0,                   0,                   0,   0.018518518518519,                   0},
                    {                0,                   0,   0.018518518518519,   0.004629629629630,                   0,                   0,   0.074074074074074,   0.018518518518519},{                0,                   0,   0.004629629629630,   0.018518518518519,                   0,                   0,   0.018518518518519,   0.074074074074074},{                0,                   0,                   0,   0.004629629629630,                   0,                   0,                   0,   0.018518518518519},{                0,                   0,                   0,                   0,   0.004629629629630,                   0,                   0,                   0},
                    {                0,                   0,                   0,                   0,   0.018518518518519,   0.004629629629630,                   0,                   0},{                0,                   0,                   0,                   0,   0.004629629629630,   0.018518518518519,                   0,                   0},{                0,                   0,                   0,                   0,                   0,   0.004629629629630,                   0,                   0},{                0,                   0,                   0,                   0,   0.018518518518519,                   0,   0.004629629629630,                   0},
                    {                0,                   0,                   0,                   0,   0.074074074074074,   0.018518518518519,   0.018518518518519,   0.004629629629630},{                0,                   0,                   0,                   0,   0.018518518518519,   0.074074074074074,   0.004629629629630,   0.018518518518519},{                0,                   0,                   0,                   0,                   0,   0.018518518518519,                   0,   0.004629629629630},{                0,                   0,                   0,                   0,   0.004629629629630,                   0,   0.018518518518519,                   0},
                    {                0,                   0,                   0,                   0,   0.018518518518519,   0.004629629629630,   0.074074074074074,   0.018518518518519},{                0,                   0,                   0,                   0,   0.004629629629630,   0.018518518518519,   0.018518518518519,   0.074074074074074},{                0,                   0,                   0,                   0,                   0,   0.004629629629630,                   0,   0.018518518518519},{                0,                   0,                   0,                   0,                   0,                   0,   0.004629629629630,                   0},
                    {                0,                   0,                   0,                   0,                   0,                   0,   0.018518518518519,   0.004629629629630},{                0,                   0,                   0,                   0,                   0,                   0,   0.004629629629630,   0.018518518518519},{                0,                   0,                   0,                   0,                   0,                   0,                   0,   0.004629629629630}};
            tT_Project = tTtemp;
        }
        else
        {
            MORIS_LOG_ERROR << " Polynomial degree is not implemented";
        }
    }
    return tT_Project;
}

//-------------------------------------------------------------------------------

Mat<real>
Hierarchical_Mesh_TMatrix::give_projection_matrix_new(
        uint const & aModelDim,
        uint const & aPolynomialLagrange,
        uint const & aPolynomial) const
{
    //Compute number of basis of an element
    uint tNumBasisOfElement = pow( aPolynomial + 1, aModelDim );
    //Compute number of design basis of an element
    uint tNumLagrangeBasisOfElement = pow( aPolynomialLagrange + 1, aModelDim );
    //Create a Matrix to project the solution from a Bspline basis to a lagrange
    Mat<real> tProjectTMatrix( tNumBasisOfElement, tNumLagrangeBasisOfElement );
    //Compute for the lagrange nodes the respective Bspline functions
    real tNumberLagrangeBasis1D  = aPolynomialLagrange + 1.0;
    if ( aModelDim == 2 )
    {
        Mat<real> tBspline;
        //Natural coordinate
        Mat<real> tNaturalCoordinate( aModelDim, 1, 0);
        //Temporary variable for loop
        uint tVar = 0;
        for ( real j = 0; j < tNumberLagrangeBasis1D; j++ )
        {
            for ( real i = 0; i < tNumberLagrangeBasis1D; i++ )
            {
                //Compute natural coordinates
                tNaturalCoordinate(0) = i / aPolynomialLagrange;
                tNaturalCoordinate(1) = j / aPolynomialLagrange;
                //Compute Bspline
                tBspline = Bspline::build_spline_uniform_nd( tNaturalCoordinate, aPolynomial, aModelDim );
                tProjectTMatrix.col( tVar ) = tBspline.col( 0 );
                tVar++;
            }
        }
    }
    else if ( aModelDim == 3 )
    {
        Mat<real> tBspline;
        //Natural coordinate
        Mat<real> tNaturalCoordinate( aModelDim, 1, 0);
        //Temporary variable for loop
        uint tVar = 0;
        for ( real k = 0; k < tNumberLagrangeBasis1D; k++ )
        {
            for ( real j = 0; j < tNumberLagrangeBasis1D; j++ )
            {
                for ( real i = 0; i < tNumberLagrangeBasis1D; i++ )
                {
                    //Compute natural coordinates
                    tNaturalCoordinate(0) = i / aPolynomialLagrange;
                    tNaturalCoordinate(1) = j / aPolynomialLagrange;
                    tNaturalCoordinate(2) = k / aPolynomialLagrange;
                    //Compute Bspline
                    tBspline = Bspline::build_spline_uniform_nd( tNaturalCoordinate, aPolynomial, aModelDim );
                    tProjectTMatrix.col( tVar ) = tBspline.col( 0 );
                    tVar++;
                }
            }
        }
    }
    return tProjectTMatrix;
}

//-------------------------------------------------------------------------------

Cell<Mat<real>>
Hierarchical_Mesh_TMatrix::give_Tmatrices_of_childs(
        uint const & aPolynomial,
        uint const & aModelDim)
{
    // Temporary variables for a loop
    uint tVar1=-1;
    uint tVar2=0;
    // Cell to store all the T-matrices of the childs
    Cell<Mat<real>> tT_child( pow( 2, aModelDim ) );
    // Left T-matrix of the parent
    Mat<real> tChild_left( aPolynomial + 1, aPolynomial + 1, 0 );
    // Left T-matrix of the parent
    Mat<real> tChild_right( aPolynomial + 1, aPolynomial + 1, 0 );
    //Number of basis function within an element
    uint tNumberBasisInElement = pow( aPolynomial + 1, aModelDim );
    if ( aPolynomial == 1 )
    {
        tChild_left( 0, 0 )  = 1;   tChild_left( 0, 1 )  = 0.5;
        tChild_left( 1, 0 )  = 0;   tChild_left( 1, 1 )  = 0.5;
        tChild_right( 0, 0 ) = 0.5; tChild_right( 0, 1 ) = 0;
        tChild_right( 1, 0 ) = 0.5; tChild_right( 1, 1 ) = 1;
    }
    else if ( aPolynomial==2 )
    {
        tChild_left( 0, 0 )  = 0.75; tChild_left( 0, 1 )  = 0.25; tChild_left( 0, 2 )  = 0;
        tChild_left( 1, 0 )  = 0.25; tChild_left( 1, 1 )  = 0.75; tChild_left( 1, 2 )  = 0.75;
        tChild_left( 2, 0 )  = 0;    tChild_left( 2, 1 )  = 0;    tChild_left( 2, 2 )  = 0.25;
        tChild_right( 0, 0 ) = 0.25; tChild_right( 0, 1 ) = 0;    tChild_right( 0, 2 ) = 0;
        tChild_right( 1, 0 ) = 0.75; tChild_right( 1, 1 ) = 0.75; tChild_right( 1, 2 ) = 0.25;
        tChild_right( 2, 0 ) = 0;    tChild_right( 2, 1 ) = 0.25; tChild_right( 2, 2 ) = 0.75;
    }
    else if ( aPolynomial==3 )
    {
        tChild_left( 0, 0 )  = 0.5;   tChild_left( 0, 1 )  = 0.125; tChild_left( 0, 2 )  = 0;     tChild_left( 0, 3 )  = 0;
        tChild_left( 1, 0 )  = 0.5;   tChild_left( 1, 1 )  = 0.75;  tChild_left( 1, 2 )  = 0.5;   tChild_left( 1, 3 )  = 0.125;
        tChild_left( 2, 0 )  = 0;     tChild_left( 2, 1 )  = 0.125; tChild_left( 2, 2 )  = 0.5;   tChild_left( 2, 3 )  = 0.75;
        tChild_left( 3, 0 )  = 0;     tChild_left( 3, 1 )  = 0;     tChild_left( 3, 2 )  = 0;     tChild_left( 3, 3 )  = 0.125;
        tChild_right( 0, 0 ) = 0.125; tChild_right( 0, 1 ) = 0;     tChild_right( 0, 2 ) = 0;     tChild_right( 0, 3 ) = 0;
        tChild_right( 1, 0 ) = 0.75;  tChild_right( 1, 1 ) = 0.5;   tChild_right( 1, 2 ) = 0.125; tChild_right( 1, 3 ) = 0;
        tChild_right( 2, 0 ) = 0.125; tChild_right( 2, 1 ) = 0.5;   tChild_right( 2, 2 ) = 0.75;  tChild_right( 2, 3 ) = 0.5;
        tChild_right( 3, 0 ) = 0;     tChild_right( 3, 1 ) = 0;     tChild_right( 3, 2 ) = 0.125; tChild_right( 3, 3 ) = 0.5;
    }
    else
    {
        MORIS_LOG_ERROR << "T-matrix for this polynomial degree is not available";
    }
    // Sace the left and right child in the cell tT_child
    if ( aModelDim == 1 )
    {
        tT_child( 0 ) = tChild_left;
        tT_child( 1 ) = tChild_right;
    }
    // T-matrix of the four childs are calculated automatically for each polynomial degree
    if ( aModelDim == 2 )
    {
        if ( aPolynomial == 1 )
        {
            Mat<real> tChild1 ={
                    {1.0000, 0.5000, 0.5000, 0.2500},
                    {0,      0.5000, 0,      0.2500},
                    {0,      0,      0.5000, 0.2500},
                    {0,      0,      0,      0.2500} };
            Mat<real> tChild2 ={
                    {0.5000, 0,      0.2500, 0     },
                    {0.5000, 1.0000, 0.2500, 0.5000},
                    {0,      0,      0.2500, 0     },
                    {0,      0,      0.2500, 0.5000} };
            Mat<real> tChild3 ={
                    {0.5000, 0.2500, 0,      0     },
                    {0,      0.2500, 0,      0     },
                    {0.5000, 0.2500, 1.0000, 0.5000},
                    {0,      0.2500, 0,      0.5000} };
            Mat<real> tChild4 ={
                    {0.2500, 0,      0,      0     },
                    {0.2500, 0.5000, 0,      0     },
                    {0.2500, 0,      0.5000, 0     },
                    {0.2500, 0.5000, 0.5000, 1.0000} };
            tT_child( 0 ) = tChild1;
            tT_child( 1 ) = tChild2;
            tT_child( 2 ) = tChild3;
            tT_child( 3 ) = tChild4;
        }
        else if ( aPolynomial == 2 )
        {
            Mat<real> tChild1 ={
                    {0.562500,0.187500,0.000000,0.187500,0.062500,0.000000,0.000000,0.000000,0.000000},
                    {0.187500,0.562500,0.562500,0.062500,0.187500,0.187500,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.187500,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000},
                    {0.187500,0.062500,0.000000,0.562500,0.187500,0.000000,0.562500,0.187500,0.000000},
                    {0.062500,0.187500,0.187500,0.187500,0.562500,0.562500,0.187500,0.562500,0.562500},
                    {0.000000,0.000000,0.062500,0.000000,0.000000,0.187500,0.000000,0.000000,0.187500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.187500,0.062500,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.062500,0.187500,0.187500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.062500}};
            Mat<real> tChild2 ={
                    {0.187500,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.562500,0.562500,0.187500,0.187500,0.187500,0.062500,0.000000,0.000000,0.000000},
                    {0.000000,0.187500,0.562500,0.000000,0.062500,0.187500,0.000000,0.000000,0.000000},
                    {0.062500,0.000000,0.000000,0.187500,0.000000,0.000000,0.187500,0.000000,0.000000},
                    {0.187500,0.187500,0.062500,0.562500,0.562500,0.187500,0.562500,0.562500,0.187500},
                    {0.000000,0.062500,0.187500,0.000000,0.187500,0.562500,0.000000,0.187500,0.562500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.187500,0.187500,0.062500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.062500,0.187500}};
            Mat<real> tChild3 ={
                    {0.187500,0.062500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.062500,0.187500,0.187500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.562500,0.187500,0.000000,0.562500,0.187500,0.000000,0.187500,0.062500,0.000000},
                    {0.187500,0.562500,0.562500,0.187500,0.562500,0.562500,0.062500,0.187500,0.187500},
                    {0.000000,0.000000,0.187500,0.000000,0.000000,0.187500,0.000000,0.000000,0.062500},
                    {0.000000,0.000000,0.000000,0.187500,0.062500,0.000000,0.562500,0.187500,0.000000},
                    {0.000000,0.000000,0.000000,0.062500,0.187500,0.187500,0.187500,0.562500,0.562500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.187500}};
            Mat<real> tChild4 ={
                    {0.062500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.187500,0.187500,0.062500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.062500,0.187500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.187500,0.000000,0.000000,0.187500,0.000000,0.000000,0.062500,0.000000,0.000000},
                    {0.562500,0.562500,0.187500,0.562500,0.562500,0.187500,0.187500,0.187500,0.062500},
                    {0.000000,0.187500,0.562500,0.000000,0.187500,0.562500,0.000000,0.062500,0.187500},
                    {0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.187500,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.187500,0.187500,0.062500,0.562500,0.562500,0.187500},
                    {0.000000,0.000000,0.000000,0.000000,0.062500,0.187500,0.000000,0.187500,0.562500}};
            tT_child( 0 ) = tChild1;
            tT_child( 1 ) = tChild2;
            tT_child( 2 ) = tChild3;
            tT_child( 3 ) = tChild4;
        }
        else if ( aPolynomial == 3 )
        {
            Mat<real> tChild1 = {
                    {0.250000,0.062500,0.000000,0.000000,0.062500,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.250000,0.375000,0.250000,0.062500,0.062500,0.093750,0.062500,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.062500,0.250000,0.375000,0.000000,0.015625,0.062500,0.093750,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.250000,0.062500,0.000000,0.000000,0.375000,0.093750,0.000000,0.000000,0.250000,0.062500,0.000000,0.000000,0.062500,0.015625,0.000000,0.000000},
                    {0.250000,0.375000,0.250000,0.062500,0.375000,0.562500,0.375000,0.093750,0.250000,0.375000,0.250000,0.062500,0.062500,0.093750,0.062500,0.015625},
                    {0.000000,0.062500,0.250000,0.375000,0.000000,0.093750,0.375000,0.562500,0.000000,0.062500,0.250000,0.375000,0.000000,0.015625,0.062500,0.093750},
                    {0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.093750,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.015625},
                    {0.000000,0.000000,0.000000,0.000000,0.062500,0.015625,0.000000,0.000000,0.250000,0.062500,0.000000,0.000000,0.375000,0.093750,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.062500,0.093750,0.062500,0.015625,0.250000,0.375000,0.250000,0.062500,0.375000,0.562500,0.375000,0.093750},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.062500,0.093750,0.000000,0.062500,0.250000,0.375000,0.000000,0.093750,0.375000,0.562500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.093750},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.062500,0.015625,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.062500,0.093750,0.062500,0.015625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.062500,0.093750},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625}};
            Mat<real> tChild2 = {
                    {0.062500,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.375000,0.250000,0.062500,0.000000,0.093750,0.062500,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.062500,0.250000,0.375000,0.250000,0.015625,0.062500,0.093750,0.062500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.062500,0.250000,0.000000,0.000000,0.015625,0.062500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.062500,0.000000,0.000000,0.000000,0.093750,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000},
                    {0.375000,0.250000,0.062500,0.000000,0.562500,0.375000,0.093750,0.000000,0.375000,0.250000,0.062500,0.000000,0.093750,0.062500,0.015625,0.000000},
                    {0.062500,0.250000,0.375000,0.250000,0.093750,0.375000,0.562500,0.375000,0.062500,0.250000,0.375000,0.250000,0.015625,0.062500,0.093750,0.062500},
                    {0.000000,0.000000,0.062500,0.250000,0.000000,0.000000,0.093750,0.375000,0.000000,0.000000,0.062500,0.250000,0.000000,0.000000,0.015625,0.062500},
                    {0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.093750,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.093750,0.062500,0.015625,0.000000,0.375000,0.250000,0.062500,0.000000,0.562500,0.375000,0.093750,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.015625,0.062500,0.093750,0.062500,0.062500,0.250000,0.375000,0.250000,0.093750,0.375000,0.562500,0.375000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.062500,0.000000,0.000000,0.062500,0.250000,0.000000,0.000000,0.093750,0.375000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.093750,0.062500,0.015625,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.062500,0.093750,0.062500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.062500}};
            Mat<real> tChild3 ={
                    {0.062500,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.062500,0.093750,0.062500,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.015625,0.062500,0.093750,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.375000,0.093750,0.000000,0.000000,0.250000,0.062500,0.000000,0.000000,0.062500,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.375000,0.562500,0.375000,0.093750,0.250000,0.375000,0.250000,0.062500,0.062500,0.093750,0.062500,0.015625,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.093750,0.375000,0.562500,0.000000,0.062500,0.250000,0.375000,0.000000,0.015625,0.062500,0.093750,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.093750,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000},
                    {0.062500,0.015625,0.000000,0.000000,0.250000,0.062500,0.000000,0.000000,0.375000,0.093750,0.000000,0.000000,0.250000,0.062500,0.000000,0.000000},
                    {0.062500,0.093750,0.062500,0.015625,0.250000,0.375000,0.250000,0.062500,0.375000,0.562500,0.375000,0.093750,0.250000,0.375000,0.250000,0.062500},
                    {0.000000,0.015625,0.062500,0.093750,0.000000,0.062500,0.250000,0.375000,0.000000,0.093750,0.375000,0.562500,0.000000,0.062500,0.250000,0.375000},
                    {0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.093750,0.000000,0.000000,0.000000,0.062500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.062500,0.015625,0.000000,0.000000,0.250000,0.062500,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.062500,0.093750,0.062500,0.015625,0.250000,0.375000,0.250000,0.062500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.062500,0.093750,0.000000,0.062500,0.250000,0.375000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.062500}};
            Mat<real> tChild4 = {
                    {0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.093750,0.062500,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.0000000},
                    {0.015625,0.062500,0.093750,0.062500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.015625,0.062500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.093750,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.562500,0.375000,0.093750,0.000000,0.375000,0.250000,0.062500,0.000000,0.093750,0.062500,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.093750,0.375000,0.562500,0.375000,0.062500,0.250000,0.375000,0.250000,0.015625,0.062500,0.093750,0.062500,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.093750,0.375000,0.000000,0.000000,0.062500,0.250000,0.000000,0.000000,0.015625,0.062500,0.000000,0.000000,0.000000,0.000000},
                    {0.015625,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.093750,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000},
                    {0.093750,0.062500,0.015625,0.000000,0.375000,0.250000,0.062500,0.000000,0.562500,0.375000,0.093750,0.000000,0.375000,0.250000,0.062500,0.000000},
                    {0.015625,0.062500,0.093750,0.062500,0.062500,0.250000,0.375000,0.250000,0.093750,0.375000,0.562500,0.375000,0.062500,0.250000,0.375000,0.250000},
                    {0.000000,0.000000,0.015625,0.062500,0.000000,0.000000,0.062500,0.250000,0.000000,0.000000,0.093750,0.375000,0.000000,0.000000,0.062500,0.250000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.093750,0.062500,0.015625,0.000000,0.375000,0.250000,0.062500,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.062500,0.093750,0.062500,0.062500,0.250000,0.375000,0.250000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.062500,0.000000,0.000000,0.062500,0.250000}};
            tT_child( 0 ) = tChild1;
            tT_child( 1 ) = tChild2;
            tT_child( 2 ) = tChild3;
            tT_child( 3 ) = tChild4;
        }
        else if ( aPolynomial >= 4 )
        {
            // Initialize the Childs
            Mat<real> tChild1( tNumberBasisInElement, tNumberBasisInElement, 0 );
            Mat<real> tChild2( tNumberBasisInElement, tNumberBasisInElement, 0 );
            Mat<real> tChild3( tNumberBasisInElement, tNumberBasisInElement, 0 );
            Mat<real> tChild4( tNumberBasisInElement, tNumberBasisInElement, 0 );

            for ( uint k = 0; k < aPolynomial + 1; k++ )
            {
                for ( uint l = 0; l < aPolynomial + 1; l++ )
                {
                    tVar1++;
                    tVar2=0;
                    for ( uint i = 0; i < aPolynomial + 1; i++ )
                    {
                        for ( uint j = 0; j < aPolynomial + 1; j++ )
                        {
                            // Calculate the Childs in 2D with the help of the Childs of 1D
                            tChild1( tVar1, tVar2 ) =  tChild_left( k, i ) *  tChild_left( l, j );
                            tChild2( tVar1, tVar2 ) =  tChild_left( k, i ) * tChild_right( l, j );
                            tChild3( tVar1, tVar2 ) = tChild_right( k, i ) *  tChild_left( l, j );
                            tChild4( tVar1, tVar2 ) = tChild_right( k, i ) * tChild_right( l, j );
                            tVar2++;
                        }
                    }
                }
            }
            //Save the four Childs in the Cell
            tT_child( 0 ) = tChild1;
            tT_child( 1 ) = tChild2;
            tT_child( 2 ) = tChild3;
            tT_child( 3 ) = tChild4;
        }
    }
    // T-matrix of the eight childs are calculated automatically for each polynomial degree
    if ( aModelDim == 3 )
    {
        if ( aPolynomial == 1 )
        {
            Mat<real> tChild1 = {
                    {1.000000,0.500000,0.500000,0.250000,0.500000,0.250000,0.250000,0.125000},
                    {0.000000,0.500000,0.000000,0.250000,0.000000,0.250000,0.000000,0.125000},
                    {0.000000,0.000000,0.500000,0.250000,0.000000,0.000000,0.250000,0.125000},
                    {0.000000,0.000000,0.000000,0.250000,0.000000,0.000000,0.000000,0.125000},
                    {0.000000,0.000000,0.000000,0.000000,0.500000,0.250000,0.250000,0.125000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.250000,0.000000,0.125000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.250000,0.125000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.125000}};
            Mat<real> tChild2 =  {
                    {0.500000,0.000000,0.250000,0.000000,0.250000,0.000000,0.125000,0.000000},
                    {0.500000,1.000000,0.250000,0.500000,0.250000,0.500000,0.125000,0.250000},
                    {0.000000,0.000000,0.250000,0.000000,0.000000,0.000000,0.125000,0.000000},
                    {0.000000,0.000000,0.250000,0.500000,0.000000,0.000000,0.125000,0.250000},
                    {0.000000,0.000000,0.000000,0.000000,0.250000,0.000000,0.125000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.250000,0.500000,0.125000,0.250000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.125000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.125000,0.250000}};
            Mat<real> tChild3 = {
                    {0.500000,0.250000,0.000000,0.000000,0.250000,0.125000,0.000000,0.000000},
                    {0.000000,0.250000,0.000000,0.000000,0.000000,0.125000,0.000000,0.000000},
                    {0.500000,0.250000,1.000000,0.500000,0.250000,0.125000,0.500000,0.250000},
                    {0.000000,0.250000,0.000000,0.500000,0.000000,0.125000,0.000000,0.250000},
                    {0.000000,0.000000,0.000000,0.000000,0.250000,0.125000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.125000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.250000,0.125000,0.500000,0.250000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.125000,0.000000,0.250000}};
            Mat<real> tChild4 = {
                    {0.250000,0.000000,0.000000,0.000000,0.125000,0.000000,0.000000,0.000000},
                    {0.250000,0.500000,0.000000,0.000000,0.125000,0.250000,0.000000,0.000000},
                    {0.250000,0.000000,0.500000,0.000000,0.125000,0.000000,0.250000,0.000000},
                    {0.250000,0.500000,0.500000,1.000000,0.125000,0.250000,0.250000,0.500000},
                    {0.000000,0.000000,0.000000,0.000000,0.125000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.125000,0.250000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.125000,0.000000,0.250000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.125000,0.250000,0.250000,0.500000}};
            Mat<real> tChild5 = {
                    {0.500000,0.250000,0.250000,0.125000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.250000,0.000000,0.125000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.250000,0.125000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.125000,0.000000,0.000000,0.000000,0.000000},
                    {0.500000,0.250000,0.250000,0.125000,1.000000,0.500000,0.500000,0.250000},
                    {0.000000,0.250000,0.000000,0.125000,0.000000,0.500000,0.000000,0.250000},
                    {0.000000,0.000000,0.250000,0.125000,0.000000,0.000000,0.500000,0.250000},
                    {0.000000,0.000000,0.000000,0.125000,0.000000,0.000000,0.000000,0.250000,}};
            Mat<real> tChild6 = {
                    {0.250000,0.000000,0.125000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.250000,0.500000,0.125000,0.250000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.125000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.125000,0.250000,0.000000,0.000000,0.000000,0.000000},
                    {0.250000,0.000000,0.125000,0.000000,0.500000,0.000000,0.250000,0.000000},
                    {0.250000,0.500000,0.125000,0.250000,0.500000,1.000000,0.250000,0.500000},
                    {0.000000,0.000000,0.125000,0.000000,0.000000,0.000000,0.250000,0.000000},
                    {0.000000,0.000000,0.125000,0.250000,0.000000,0.000000,0.250000,0.500000}};
            Mat<real> tChild7 = {
                    {0.250000,0.125000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.125000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.250000,0.125000,0.500000,0.250000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.125000,0.000000,0.250000,0.000000,0.000000,0.000000,0.000000},
                    {0.250000,0.125000,0.000000,0.000000,0.500000,0.250000,0.000000,0.000000},
                    {0.000000,0.125000,0.000000,0.000000,0.000000,0.250000,0.000000,0.000000},
                    {0.250000,0.125000,0.500000,0.250000,0.500000,0.250000,1.000000,0.500000},
                    {0.000000,0.125000,0.000000,0.250000,0.000000,0.250000,0.000000,0.500000,}};
            Mat<real> tChild8 = {
                    {0.125000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.125000,0.250000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.125000,0.000000,0.250000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.125000,0.250000,0.250000,0.500000,0.000000,0.000000,0.000000,0.000000},
                    {0.125000,0.000000,0.000000,0.000000,0.250000,0.000000,0.000000,0.000000},
                    {0.125000,0.250000,0.000000,0.000000,0.250000,0.500000,0.000000,0.000000},
                    {0.125000,0.000000,0.250000,0.000000,0.250000,0.000000,0.500000,0.000000},
                    {0.125000,0.250000,0.250000,0.500000,0.250000,0.500000,0.500000,1.000000,}};
            tT_child( 0 ) = tChild1;
            tT_child( 1 ) = tChild2;
            tT_child( 2 ) = tChild3;
            tT_child( 3 ) = tChild4;
            tT_child( 4 ) = tChild5;
            tT_child( 5 ) = tChild6;
            tT_child( 6 ) = tChild7;
            tT_child( 7 ) = tChild8;
        }
        else if ( aPolynomial == 2 )
        {
            Mat<real> tChild1 = {
                    {0.421875,0.140625,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.421875,0.421875,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.046875,0.000000,0.421875,0.140625,0.000000,0.421875,0.140625,0.000000,0.046875,0.015625,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.140625,0.140625,0.140625,0.421875,0.421875,0.140625,0.421875,0.421875,0.015625,0.046875,0.046875,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.046875,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.421875,0.140625,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.421875,0.140625,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.140625,0.140625,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.140625,0.421875,0.421875,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.140625,0.421875,0.421875,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000},
                    {0.046875,0.015625,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.421875,0.140625,0.000000,0.421875,0.140625,0.000000,0.140625,0.046875,0.000000,0.421875,0.140625,0.000000,0.421875,0.140625,0.000000},
                    {0.015625,0.046875,0.046875,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.140625,0.421875,0.421875,0.140625,0.421875,0.421875,0.046875,0.140625,0.140625,0.140625,0.421875,0.421875,0.140625,0.421875,0.421875},
                    {0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625}};
            Mat<real> tChild2 = {
                    {0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.421875,0.421875,0.140625,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.140625,0.421875,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.140625,0.046875,0.421875,0.421875,0.140625,0.421875,0.421875,0.140625,0.046875,0.046875,0.015625,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.046875,0.140625,0.000000,0.140625,0.421875,0.000000,0.140625,0.421875,0.000000,0.015625,0.046875,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.140625,0.046875,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.421875,0.421875,0.140625,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.421875,0.421875,0.140625,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000},
                    {0.000000,0.046875,0.140625,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.140625,0.421875,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.140625,0.421875,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000},
                    {0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000},
                    {0.046875,0.046875,0.015625,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.421875,0.421875,0.140625,0.421875,0.421875,0.140625,0.140625,0.140625,0.046875,0.421875,0.421875,0.140625,0.421875,0.421875,0.140625},
                    {0.000000,0.015625,0.046875,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.140625,0.421875,0.000000,0.140625,0.421875,0.000000,0.046875,0.140625,0.000000,0.140625,0.421875,0.000000,0.140625,0.421875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875}};
            Mat<real> tChild3 = {
                    {0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.421875,0.140625,0.000000,0.421875,0.140625,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.421875,0.421875,0.140625,0.421875,0.421875,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.421875,0.140625,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.140625,0.421875,0.421875,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.046875,0.015625,0.000000,0.421875,0.140625,0.000000,0.421875,0.140625,0.000000,0.140625,0.046875,0.000000,0.421875,0.140625,0.000000,0.421875,0.140625,0.000000,0.140625,0.046875,0.000000},
                    {0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.015625,0.046875,0.046875,0.140625,0.421875,0.421875,0.140625,0.421875,0.421875,0.046875,0.140625,0.140625,0.140625,0.421875,0.421875,0.140625,0.421875,0.421875,0.046875,0.140625,0.140625},
                    {0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875},
                    {0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.421875,0.140625,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.421875,0.140625,0.000000},
                    {0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.140625,0.421875,0.421875,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.140625,0.421875,0.421875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.046875,0.015625,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.015625,0.046875,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.140625,0.046875,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.046875,0.140625,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875}};
            Mat<real> tChild4 = {
                    {0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.421875,0.421875,0.140625,0.421875,0.421875,0.140625,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.140625,0.421875,0.000000,0.140625,0.421875,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.421875,0.421875,0.140625,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.140625,0.421875,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000},
                    {0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.046875,0.046875,0.015625,0.421875,0.421875,0.140625,0.421875,0.421875,0.140625,0.140625,0.140625,0.046875,0.421875,0.421875,0.140625,0.421875,0.421875,0.140625,0.140625,0.140625,0.046875},
                    {0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.015625,0.046875,0.000000,0.140625,0.421875,0.000000,0.140625,0.421875,0.000000,0.046875,0.140625,0.000000,0.140625,0.421875,0.000000,0.140625,0.421875,0.000000,0.046875,0.140625},
                    {0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.421875,0.421875,0.140625,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.421875,0.421875,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.140625,0.421875,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.140625,0.421875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.046875,0.046875,0.015625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.015625,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.140625,0.140625,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.046875,0.140625}};
            Mat<real> tChild5 = {
                    {0.140625,0.046875,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.140625,0.140625,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.015625,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.015625,0.046875,0.046875,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.421875,0.140625,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.421875,0.140625,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.421875,0.421875,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.140625,0.421875,0.421875,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000},
                    {0.140625,0.046875,0.000000,0.421875,0.140625,0.000000,0.421875,0.140625,0.000000,0.140625,0.046875,0.000000,0.421875,0.140625,0.000000,0.421875,0.140625,0.000000,0.046875,0.015625,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000},
                    {0.046875,0.140625,0.140625,0.140625,0.421875,0.421875,0.140625,0.421875,0.421875,0.046875,0.140625,0.140625,0.140625,0.421875,0.421875,0.140625,0.421875,0.421875,0.015625,0.046875,0.046875,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625},
                    {0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.421875,0.140625,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.140625,0.421875,0.421875,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.421875,0.140625,0.000000,0.421875,0.140625,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.140625,0.421875,0.421875,0.140625,0.421875,0.421875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875}};
            Mat<real> tChild6 = {
                    {0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.140625,0.046875,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.046875,0.140625,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.046875,0.015625,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.015625,0.046875,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.421875,0.421875,0.140625,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.421875,0.421875,0.140625,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000},
                    {0.000000,0.140625,0.421875,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.140625,0.421875,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000},
                    {0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000},
                    {0.140625,0.140625,0.046875,0.421875,0.421875,0.140625,0.421875,0.421875,0.140625,0.140625,0.140625,0.046875,0.421875,0.421875,0.140625,0.421875,0.421875,0.140625,0.046875,0.046875,0.015625,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875},
                    {0.000000,0.046875,0.140625,0.000000,0.140625,0.421875,0.000000,0.140625,0.421875,0.000000,0.046875,0.140625,0.000000,0.140625,0.421875,0.000000,0.140625,0.421875,0.000000,0.015625,0.046875,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.421875,0.421875,0.140625,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.140625,0.421875,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.421875,0.421875,0.140625,0.421875,0.421875,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.140625,0.421875,0.000000,0.140625,0.421875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625}};
            Mat<real> tChild7 = {
                    {0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.421875,0.140625,0.000000,0.421875,0.140625,0.000000,0.140625,0.046875,0.000000,0.421875,0.140625,0.000000,0.421875,0.140625,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.046875,0.015625,0.000000},
                    {0.140625,0.421875,0.421875,0.140625,0.421875,0.421875,0.046875,0.140625,0.140625,0.140625,0.421875,0.421875,0.140625,0.421875,0.421875,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.015625,0.046875,0.046875},
                    {0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625},
                    {0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.421875,0.140625,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.421875,0.140625,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.140625,0.046875,0.000000},
                    {0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.140625,0.421875,0.421875,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.140625,0.421875,0.421875,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.046875,0.140625,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.140625,0.046875,0.000000,0.046875,0.015625,0.000000,0.421875,0.140625,0.000000,0.421875,0.140625,0.000000,0.140625,0.046875,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.015625,0.046875,0.046875,0.140625,0.421875,0.421875,0.140625,0.421875,0.421875,0.046875,0.140625,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.421875,0.140625,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.046875,0.046875,0.140625,0.140625,0.000000,0.000000,0.000000,0.046875,0.140625,0.140625,0.140625,0.421875,0.421875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625}};
            Mat<real> tChild8 = {
                    {0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000},
                    {0.421875,0.421875,0.140625,0.421875,0.421875,0.140625,0.140625,0.140625,0.046875,0.421875,0.421875,0.140625,0.421875,0.421875,0.140625,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.046875,0.046875,0.015625},
                    {0.000000,0.140625,0.421875,0.000000,0.140625,0.421875,0.000000,0.046875,0.140625,0.000000,0.140625,0.421875,0.000000,0.140625,0.421875,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.015625,0.046875},
                    {0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.421875,0.421875,0.140625,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.421875,0.421875,0.140625,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.140625,0.140625,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.140625,0.421875,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.140625,0.421875,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.046875,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.140625,0.000000,0.000000,0.140625,0.000000,0.000000,0.046875,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.140625,0.140625,0.046875,0.046875,0.046875,0.015625,0.421875,0.421875,0.140625,0.421875,0.421875,0.140625,0.140625,0.140625,0.046875},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.046875,0.140625,0.000000,0.015625,0.046875,0.000000,0.140625,0.421875,0.000000,0.140625,0.421875,0.000000,0.046875,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.140625,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.046875,0.015625,0.140625,0.140625,0.046875,0.000000,0.000000,0.000000,0.140625,0.140625,0.046875,0.421875,0.421875,0.140625},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.046875,0.000000,0.046875,0.140625,0.000000,0.000000,0.000000,0.000000,0.046875,0.140625,0.000000,0.140625,0.421875}};
            tT_child( 0 ) = tChild1;
            tT_child( 1 ) = tChild2;
            tT_child( 2 ) = tChild3;
            tT_child( 3 ) = tChild4;
            tT_child( 4 ) = tChild5;
            tT_child( 5 ) = tChild6;
            tT_child( 6 ) = tChild7;
            tT_child( 7 ) = tChild8;
        }
        else if ( aPolynomial >= 3 )
        {
            // Initialize the Childs
            Mat<real> tChild1( tNumberBasisInElement, tNumberBasisInElement, 0 );
            Mat<real> tChild2( tNumberBasisInElement, tNumberBasisInElement, 0 );
            Mat<real> tChild3( tNumberBasisInElement, tNumberBasisInElement, 0 );
            Mat<real> tChild4( tNumberBasisInElement, tNumberBasisInElement, 0 );
            Mat<real> tChild5( tNumberBasisInElement, tNumberBasisInElement, 0 );
            Mat<real> tChild6( tNumberBasisInElement, tNumberBasisInElement, 0 );
            Mat<real> tChild7( tNumberBasisInElement, tNumberBasisInElement, 0 );
            Mat<real> tChild8( tNumberBasisInElement, tNumberBasisInElement, 0 );

            for ( uint k = 0; k < aPolynomial + 1; k++ )
            {
                for ( uint l = 0; l < aPolynomial + 1; l++ )
                {
                    for ( uint m = 0; m < aPolynomial + 1; m++ )
                    {
                        tVar1++;
                        tVar2=0;
                        for (uint i = 0; i < aPolynomial + 1; i++ )
                        {
                            for (uint j = 0; j < aPolynomial + 1; j++ )
                            {
                                for (uint n = 0; n < aPolynomial + 1; n++ )
                                {
                                    // Calculate the Childs in 3D with the help of the Childs of 1D
                                    tChild1( tVar1, tVar2 ) =  tChild_left( m, n ) *  tChild_left( l, j ) *  tChild_left( k, i );
                                    tChild2( tVar1, tVar2 ) = tChild_right( m, n ) *  tChild_left( l, j ) *  tChild_left( k, i );
                                    tChild3( tVar1, tVar2 ) =  tChild_left( m, n ) * tChild_right( l, j ) *  tChild_left( k, i );
                                    tChild4( tVar1, tVar2 ) = tChild_right( m, n ) * tChild_right( l, j ) *  tChild_left( k, i );
                                    tChild5( tVar1, tVar2 ) =  tChild_left( m, n ) *  tChild_left( l, j ) * tChild_right( k, i );
                                    tChild6( tVar1, tVar2 ) = tChild_right( m, n ) *  tChild_left( l, j ) * tChild_right( k, i );
                                    tChild7( tVar1, tVar2 ) =  tChild_left( m, n ) * tChild_right( l, j ) * tChild_right( k, i );
                                    tChild8( tVar1, tVar2 ) = tChild_right( m, n ) * tChild_right( l, j ) * tChild_right( k, i );
                                    tVar2++;
                                }
                            }
                        }
                    }
                }
            }
            //Save the four Childs in the Cell
            tT_child( 0 ) = tChild1;
            tT_child( 1 ) = tChild2;
            tT_child( 2 ) = tChild3;
            tT_child( 3 ) = tChild4;
            tT_child( 4 ) = tChild5;
            tT_child( 5 ) = tChild6;
            tT_child( 6 ) = tChild7;
            tT_child( 7 ) = tChild8;
        }
    }
    return tT_child;
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_TMatrix::give_vector_for_reorder(
        uint const & aModelDim,
        uint const & aPolynomial) const
{
    Mat<uint> tOrder( pow( aPolynomial + 1, aModelDim ), 1, 0);
    if ( aModelDim == 2)
    {
        if( aPolynomial == 1 )
        {
            tOrder( 0 ) = 0;
            tOrder( 1 ) = 1;
            tOrder( 2 ) = 3;
            tOrder( 3 ) = 2;
        }
        else if( aPolynomial == 2 )
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
    else if ( aModelDim == 3)
    {
        if( aPolynomial == 1 )
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
        else if( aPolynomial == 2 )
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

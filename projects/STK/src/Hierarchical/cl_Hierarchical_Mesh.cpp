/*
 * cl_Hierarchical_Mesh.cpp
 *
 *  Created on: May 10, 2017
 *      Author: gleim
 */
#include "cl_Hierarchical_Mesh.hpp" // STK/src/Hierarchical
using namespace moris;
Mat<uint>
moris::Hierarchical_Mesh::give_vector_entries(Mat<uint> & aMat, Mat<uint> & bMat)
{
    Mat<uint> cMat(bMat.length(),1,0);
    for(uint i = 0; i< bMat.length(); i++)
    {
        cMat(i) = aMat(bMat(i));
    }
    return cMat;
}

Mat<real>
moris::Hierarchical_Mesh::give_vector_entries(Mat<real> & aMat, Mat<uint> & bMat)
{
    Mat<real> cMat(bMat.length(),1,0);
    for(uint i = 0; i< bMat.length(); i++)
    {
        cMat(i) = aMat(bMat(i));
    }
    return cMat;
}

void moris::Hierarchical_Mesh::give_Tmatrix_and_id_field(uint & aElementNumber,
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    Mat<uint> tOrder(pow((aElementList.Polynomial)+1,(aElementList.Dim)),1); // Change ordering to the classical order of FE connectivity (for FemDoc and Paraview)
    if( aElementList.Dim == 2)
    {
        tOrder(0) = 0;  tOrder(1) = 1;  tOrder(2) = 3; tOrder(3) = 2;
    }
    else if( aElementList.Dim == 3)
    {
        tOrder(0) = 0;  tOrder(1) = 1;  tOrder(2) = 3; tOrder(3) = 2; tOrder(4) = 4; tOrder(5) = 5; tOrder(6) = 7; tOrder(7) = 6;
    }
    uint tLevel = give_element_level(aElementNumber,aElementList.Dim,aElementList.NumElements)+1;// Level of the element
    Cell<Mat<real>> ttT_child = give_Tmatrices_of_childs(aElementList.Polynomial,aElementList.Dim); // Give the T-matrices of the childs of the element "aElementNumber"
    Mat<real> tT = eye( pow(aElementList.Polynomial+1,aElementList.Dim), pow(aElementList.Polynomial+1,aElementList.Dim) ); //Tmatrix of the curent element
    Mat<real> tTMatrix(tLevel*pow(aElementList.Polynomial+1,aElementList.Dim),pow(aElementList.Polynomial+1,aElementList.Dim),0); //Tmatrix of the parent element
    Mat<uint> tIdField(tLevel*pow(aElementList.Polynomial+1,aElementList.Dim),1,0); //IdField of the parent Element
    Mat<uint> tBasis = give_basis_of_element(aElementNumber,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements); //Basis of the curent element
    uint tVar = 0; //Temporary variable for a loop
    Mat<real> tMatrixDummy;

    for(uint i = 0; i<pow(aElementList.Polynomial+1,aElementList.Dim); i++)
    {
        if( (aBasisList.BasisActive).test(tBasis(tOrder(i))) == 1 )
        {
            tMatrixDummy = tT.row(tOrder(i)); //Add the row of the Tmatrix of the active basis functions of the curent element
            for(uint j = 0; j<pow(aElementList.Polynomial+1,aElementList.Dim); j++)
                tTMatrix(tVar,j) =  tMatrixDummy(tOrder(j));
            tIdField(tVar) = tBasis(tOrder(i)); //Add the row of the IdField of the active basis functions of the curent element
            tVar++;
        }
    }
    tLevel = tLevel - 1;     //    //Asking for T-matrix on level l-1 (one coarser level)
    uint tChild=aElementNumber;
    uint tParent;
    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    while(tLevel>0){
        tParent = give_parent_of_element(tChild,aElementList.Dim,aElementList.NumElements); // Who is the parent of the child
        tChildren = give_children_of_element(tParent,aElementList.Dim,aElementList.NumElements); // Give the children of the parent element
        tBasis = give_basis_of_element(tParent,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements); //Basis of the parent element
        for(uint i = 0; i<pow(aElementList.Polynomial+1,aElementList.Dim); i++)
        {
            if(tChildren(i)==tChild) //Give me the T-matrix of the child for the parent
            {
                Mat<real> tTmatrix_child = ttT_child(i);
                tT=tTmatrix_child*tT; //Create T-matrix for the parent element
                break;
            }
        }
        for(uint i = 0; i<pow(aElementList.Polynomial+1,aElementList.Dim); i++)
        {
            if( (aBasisList.BasisActive).test(tBasis(tOrder(i))) == 1 )
            {
                tMatrixDummy = tT.row(tOrder(i)); //Add the row of the Tmatrix of the active basis functions of the curent element
                for(uint j = 0; j<pow(aElementList.Polynomial+1,aElementList.Dim); j++)
                    tTMatrix(tVar,j) =  tMatrixDummy(tOrder(j));
                tIdField(tVar) = tBasis(tOrder(i)); //Add the row of the IdField of the active basis functions of the parents
                tVar++;
            }
        }
        tLevel = tLevel - 1;
        tChild=tParent;
    }
    //Resize the Tmatrix and the IdField
    tTMatrix.resize(tVar,pow(aElementList.Polynomial+1,aElementList.Dim));
    tIdField.resize(tVar,1);
    // Save data in a struc
    aElementList.TMatrix = tTMatrix;
    aElementList.IdField = tIdField;
}

void moris::Hierarchical_Mesh::give_Truncated_Tmatrix_and_id_field(uint & aElementNumber,
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    Mat<uint> tOrder(pow((aElementList.Polynomial)+1,(aElementList.Dim)),1); // Change ordering to the classical order of FE connectivity (for FemDoc and Paraview)
    if( aElementList.Dim == 2)
    {
        tOrder(0) = 0;  tOrder(1) = 1;  tOrder(2) = 3; tOrder(3) = 2;
    }
    else if( aElementList.Dim == 3)
    {
        tOrder(0) = 0;  tOrder(1) = 1;  tOrder(2) = 3; tOrder(3) = 2; tOrder(4) = 4; tOrder(5) = 5; tOrder(6) = 7; tOrder(7) = 6;
    }
    uint tLevel = give_element_level(aElementNumber,aElementList.Dim,aElementList.NumElements)+1;// Level of the element
    Cell<Mat<real>> ttT_child = give_Tmatrices_of_childs(aElementList.Polynomial,aElementList.Dim); // Give the T-matrices of the childs of the element "aElementNumber"
    Mat<real> tT = eye( pow(aElementList.Polynomial+1,aElementList.Dim), pow(aElementList.Polynomial+1,aElementList.Dim) ); //Tmatrix of the curent element
    Mat<real> tTeye = eye( pow(aElementList.Polynomial+1,aElementList.Dim), pow(aElementList.Polynomial+1,aElementList.Dim) ); //Tmatrix of the curent element
    Mat<real> tTMatrix(tLevel*pow(aElementList.Polynomial+1,aElementList.Dim),pow(aElementList.Polynomial+1,aElementList.Dim),0); //Tmatrix of the parent element
    Mat<real> aTMatrix(tLevel*pow(aElementList.Polynomial+1,aElementList.Dim),pow(aElementList.Polynomial+1,aElementList.Dim),0); //Tmatrix of the parent element
    Mat<uint> tIdField(tLevel*pow(aElementList.Polynomial+1,aElementList.Dim),1,0); //IdField of the parent Element
    Mat<uint> tBasis = give_basis_of_element(aElementNumber,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements); //Basis of the curent element
    uint tVar = 0; uint tVarb = 0; uint tVarc = 0; uint tVarcell = 0; //Temporary variable for a loop
    Mat<real> tMatrixDummy; Mat<real> tTmatrix_child_dummy; Mat<real> tTeyeDummy;
    Cell<Mat<real>> tMatrixActiveBasisEye(tLevel); // Stores a active basis functions for each level with the eye matrix
    Cell<Mat<real>> tTMatrixOnLevel(tLevel);//  Stores a T-matrices for each level
    Mat<real> tMatrixActiveBasis(tLevel*pow(aElementList.Polynomial+1,aElementList.Dim),pow(aElementList.Polynomial+1,aElementList.Dim),0);
    Mat<real> tTMatrixOfLevel(tLevel*pow(aElementList.Polynomial+1,aElementList.Dim),pow(aElementList.Polynomial+1,aElementList.Dim),0);
    Mat<real> tTMatrixDummy(tLevel*pow(aElementList.Polynomial+1,aElementList.Dim),pow(aElementList.Polynomial+1,aElementList.Dim),0);
    Mat<real> tTruncatedTMatrixOfLevel;
    Mat<real> tMatrixActiveBasisDummy;
    Mat<real> tTmatrix_child;

    for(uint i = 0; i<pow(aElementList.Polynomial+1,aElementList.Dim); i++)
    {
        if( (aBasisList.BasisActive).test(tBasis(tOrder(i))) == 1 )
        {
            tMatrixDummy = tT.row(tOrder(i)); //Add the row of the Tmatrix of the active basis functions of the curent element
            for(uint j = 0; j<pow(aElementList.Polynomial+1,aElementList.Dim); j++)
            {
                aTMatrix(tVar,j) =  tMatrixDummy(tOrder(j));
                tTMatrix(tVar,j) =  tMatrixDummy(tOrder(j));
                tMatrixActiveBasis(tVar,j) = tMatrixDummy(tOrder(j)); //Add the row of the Tmatrix of the active basis functions of the curent element
                tTMatrixDummy(tVar,j) = tMatrixDummy(tOrder(j)); //Add the row of the Tmatrix of the active basis functions of the curent element
            }
            tIdField(tVar) = tBasis(tOrder(i)); //Add the row of the IdField of the active basis functions of the curent element
            tVar++;
        }
    }

    tVarb = tVar; // Save curent size of aIdField
    tVarc = tVar; // Save curent size of aTMatrix
    tTMatrixDummy.resize(tVar,pow(aElementList.Polynomial+1,aElementList.Dim));
    tMatrixActiveBasis.resize(tVar,pow(aElementList.Polynomial+1,aElementList.Dim));
    tMatrixActiveBasisEye(tVarcell) = tMatrixActiveBasis; // Save active basis functions with a one for their position
    tTMatrixOnLevel(tVarcell) = tTMatrixDummy; // Save the TMatrix of the current level
    tVarcell++;

    //    aTMatrix.print("aTMatrix");
    //    tTMatrix.print("tTMatrix");
    //    tMatrixActiveBasis.print("tMatrixActiveBasis");
    //    tTMatrixDummy.print("tTMatrixDummy");
    tLevel = tLevel - 1;     //    //Asking for T-matrix on level l-1 (one coarser level)
    uint tChild=aElementNumber;
    uint tParent;
    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    while(tLevel>0){
        tParent = give_parent_of_element(tChild,aElementList.Dim,aElementList.NumElements); // Who is the parent of the child
        tChildren = give_children_of_element(tParent,aElementList.Dim,aElementList.NumElements); // Give the children of the parent element
        tBasis = give_basis_of_element(tParent,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements); //Basis of the parent element
        for(uint i = 0; i<pow(aElementList.Polynomial+1,aElementList.Dim); i++)
        {
            if(tChildren(i)==tChild) //Give me the T-matrix of the child for the parent
            {
                tTmatrix_child = ttT_child(i);
                //                tTmatrix_child.print("tTmatrix_child");
                tT=tTmatrix_child*tT; //Create T-matrix for the parent element
                //                tT.print("tT");
                break;
            }
        }
        tTMatrix.set_size(pow(aElementList.Polynomial+1,aElementList.Dim),pow(aElementList.Polynomial+1,aElementList.Dim),0);
        tTMatrixOfLevel.set_size(pow(aElementList.Polynomial+1,aElementList.Dim),pow(aElementList.Polynomial+1,aElementList.Dim),0);
        tMatrixActiveBasis.set_size(pow(aElementList.Polynomial+1,aElementList.Dim),pow(aElementList.Polynomial+1,aElementList.Dim),0);
        tVar = 0;
        for(uint i = 0; i<pow(aElementList.Polynomial+1,aElementList.Dim); i++)
        {
            if( (aBasisList.BasisActive).test(tBasis(tOrder(i))) == 1 )
            {
                tMatrixDummy = tT.row(tOrder(i)); //Add the row of the Tmatrix of the active basis functions of the curent element
                tTmatrix_child_dummy = tTmatrix_child.row(tOrder(i)); //Add the row of the Tmatrix of the active basis functions of the parents
                tTeyeDummy = tTeye.row(tOrder(i));  // Add active basis functions with a one for their position
                for(uint j = 0; j<pow(aElementList.Polynomial+1,aElementList.Dim); j++)
                {
                    tTMatrix(tVar,j) =  tMatrixDummy(tOrder(j));
                    tMatrixActiveBasis(tVar,j) = tTeyeDummy(tOrder(j));
                    tTMatrixOfLevel(tVar,j) = tTmatrix_child_dummy(tOrder(j));

                }
                tIdField(tVarb) = tBasis(tOrder(i)); //Add the row of the IdField of the active basis functions of the parents
                tVar++;
                tVarb++;
            }
        }
        tLevel = tLevel - 1;
        tChild=tParent;

        tTMatrix.resize(tVar,pow(aElementList.Polynomial+1,aElementList.Dim));
        tMatrixActiveBasis.resize(tVar,pow(aElementList.Polynomial+1,aElementList.Dim));
        tTMatrixOfLevel.resize(tVar,pow(aElementList.Polynomial+1,aElementList.Dim));
        tMatrixActiveBasisDummy = tMatrixActiveBasisEye(tVarcell-1);
        tTMatrixDummy = tTMatrixOnLevel(tVarcell-1);

        //        tTMatrixOfLevel.print("tTMatrixOfLevel");
        //        tMatrixActiveBasisDummy.print("tMatrixActiveBasisDummy");
        //        tTMatrixDummy.print("tTMatrixDummy");
        tTruncatedTMatrixOfLevel = tTMatrix - (( tTMatrixOfLevel * trans(tMatrixActiveBasisDummy) ) * tTMatrixDummy); // Take the T-Matrix of the current level and subrtract parts from the children level, which have an influence
        //        tTruncatedTMatrixOfLevel = tTMatrix - (( tTMatrixOfLevel * trans(tMatrixActiveBasisDummy) ) * tMatrixActiveBasisDummy); // Take the T-Matrix of the current level and subrtract parts from the children level, which have an influence

        //        tTruncatedTMatrixOfLevel.print("tTruncatedTMatrixOfLevel");

        tTMatrixOnLevel(tVarcell) = tTruncatedTMatrixOfLevel;
        tMatrixActiveBasisEye(tVarcell) = tMatrixActiveBasis;
        tVarcell++;
        for(uint i = 0; i < tVarcell-2; i++)
        {
            tMatrixActiveBasis = tMatrixActiveBasisEye(i);
            tTMatrix = tTMatrixOnLevel(i);

            //            tTruncatedTMatrixOfLevel.print("tTruncatedTMatrixOfLevel");
            //            tMatrixActiveBasis.print("tMatrixActiveBasis");
            //            tTMatrix.print("tTMatrix");

            //            tTruncatedTMatrixOfLevel = tTruncatedTMatrixOfLevel - (tTruncatedTMatrixOfLevel * trans(tMatrixActiveBasis)) * tTMatrix; // Take the T-Matrix of the current level and subtract parts from all finer levels, which have an influence
            tTruncatedTMatrixOfLevel = tTruncatedTMatrixOfLevel - (tTruncatedTMatrixOfLevel * trans(tMatrixActiveBasis)) * tMatrixActiveBasis; // Take the T-Matrix of the current level and subtract parts from all finer levels, which have an influence

            //            tTruncatedTMatrixOfLevel.print("tTruncatedTMatrixOfLevel");
        }
        if( tVar > 0)
            aTMatrix.rows(tVarc,tVarb-1) = tTruncatedTMatrixOfLevel.rows(0,tVar-1);
        tVarc = tVarb;
    }
    //Resize the Tmatrix and the IdField
    aTMatrix.resize(tVarb,pow(aElementList.Polynomial+1,aElementList.Dim));
    tIdField.resize(tVarb,1);
    // Save data in a struc
    aElementList.TMatrix = aTMatrix;
    aElementList.IdField = tIdField;
}

void moris::Hierarchical_Mesh::give_Tmatrix_and_id_field_for_designvariables(uint & aElementNumber,
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tLevel = give_element_level(aElementNumber,aElementList.Dim,aElementList.NumElements)+1;// Level of the element
    Cell<Mat<real>> ttT_child = give_Tmatrices_of_childs(aElementList.PolynomialDesign,aElementList.Dim); // Give the T-matrices of the childs of the element "aElementNumber"
    Mat<real> tT = eye( pow(aElementList.PolynomialDesign+1,aElementList.Dim), pow(aElementList.PolynomialDesign+1,aElementList.Dim) ); //Tmatrix of the curent element
    Mat<real> tT_Project;//(pow(aElementList.PolynomialDesign+1,aElementList.Dim), pow(aElementList.PolynomialDesign+1,aElementList.Dim),0); //Projection to the FEM mesh
    Mat<real> tTMatrix(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0); //Tmatrix of the parent element
    Mat<uint> tIdField(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),1,0); //IdField of the parent Element
    Mat<uint> tBasis = give_basis_of_element(aElementNumber,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements); //Basis of the curent element
    uint tVar = 0; //Temporary variable for a loop
    tT_Project =  give_projection_matrix(aElementList);

    for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
    {
        if( (aBasisList.DesignBSplineActive).test(tBasis(i)) == 1 )
        {
            tTMatrix.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the curent element
            tIdField(tVar) = tBasis(i); //Add the row of the IdField of the active basis functions of the curent element
            tVar++;
        }
    }
    tLevel = tLevel - 1;
    //    //Asking for T-matrix on level l-1 (one coarser level)
    uint tChild=aElementNumber;
    uint tParent;
    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    while(tLevel>0){
        tParent = give_parent_of_element(tChild,aElementList.Dim,aElementList.NumElements); // Who is the parent of the child
        tChildren = give_children_of_element(tParent,aElementList.Dim,aElementList.NumElements); // Give the children of the parent element
        tBasis = give_basis_of_element(tParent,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements); //Basis of the parent element
        for(uint i = 0; i<tChildren.length(); i++)
        {
            if(tChildren(i)==tChild) //Give me the T-matrix of the child for the parent
            {
                Mat<real> tTmatrix_child = ttT_child(i);
                tT=tTmatrix_child*tT; //Create T-matrix for the parent element
                break;
            }
        }
        for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
        {
            if( (aBasisList.DesignBSplineActive).test(tBasis(i)) == 1 )
            {
                tTMatrix.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the parents
                tIdField(tVar) = tBasis(i); //Add the row of the IdField of the active basis functions of the parents
                tVar++;
            }
        }
        tLevel = tLevel - 1;
        tChild=tParent;
    }
    //Resize the Tmatrix and the IdField
    tTMatrix.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
    tIdField.resize(tVar,1);
    // Save data in a struc and multiply t-matrix with projection matrix
    aElementList.TMatrixDesign = tTMatrix * tT_Project;
    aElementList.IdFieldDesign = tIdField;
}

void moris::Hierarchical_Mesh::give_Truncated_Tmatrix_and_id_field_for_designvariables(uint & aElementNumber,
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tLevel = give_element_level(aElementNumber,aElementList.Dim,aElementList.NumElements)+1;// Level of the element
    Cell<Mat<real>> ttT_child = give_Tmatrices_of_childs(aElementList.PolynomialDesign,aElementList.Dim); // Give the T-matrices of the childs of the element "aElementNumber"
    Mat<real> tT = eye( pow(aElementList.PolynomialDesign+1,aElementList.Dim), pow(aElementList.PolynomialDesign+1,aElementList.Dim) ); //Tmatrix of the curent element
    Mat<real> tTeye = eye( pow(aElementList.PolynomialDesign+1,aElementList.Dim), pow(aElementList.PolynomialDesign+1,aElementList.Dim) ); //Tmatrix of the curent element
    Mat<real> tT_Project;//(pow(aElementList.PolynomialDesign+1,aElementList.Dim), pow(aElementList.PolynomialDesign+1,aElementList.Dim),0); //Projection to the FEM mesh
    Mat<real> aTMatrix(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0); //Tmatrix of the parent element
    Mat<uint> tIdField(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),1,0); //IdField of the parent Element
    Mat<uint> tBasis = give_basis_of_element(aElementNumber,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements); //Basis of the curent element
    uint tVar = 0; uint tVarb = 0; uint tVarc = 0; uint tVarcell = 0; //Temporary variable for a loop
    Mat<real> tTmatrix_child; //Stores the tMatrix of a cell of the current child in a matrix
    Cell<Mat<real>> tMatrixActiveBasisEye(tLevel); // Stores a active basis functions for each level with the eye matrix
    Cell<Mat<real>> tTMatrixOnLevel(tLevel);//  Stores a T-matrices for each level
    Mat<real> tMatrixActiveBasis(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
    Mat<real> tTMatrixOfLevel(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
    Mat<real> tTMatrix(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
    Mat<real> tTMatrixDummy(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
    Mat<real> tTruncatedTMatrixOfLevel;
    Mat<real> tMatrixActiveBasisDummy;
    tT_Project =  give_projection_matrix(aElementList);

    for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
    {
        if( (aBasisList.DesignBSplineActive).test(tBasis(i)) == 1 )
        {
            aTMatrix.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the curent element
            tMatrixActiveBasis.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the curent element
            tTMatrixDummy.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the curent element
            tIdField(tVar) = tBasis(i); //Add the row of the IdField of the active basis functions of the curent element
            tVar++;
        }
    }
    tVarb = tVar; // Save curent size of aIdField
    tVarc = tVar; // Save curent size of aTMatrix
    tTMatrixDummy.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
    tMatrixActiveBasis.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
    tMatrixActiveBasisEye(tVarcell) = tMatrixActiveBasis; // Save active basis functions with a one for their position
    tTMatrixOnLevel(tVarcell) = tTMatrixDummy; // Save the TMatrix of the current level
    tVarcell++;
    tLevel = tLevel - 1;
    //    //Asking for T-matrix on level l-1 (one coarser level)
    uint tChild=aElementNumber;
    uint tParent;
    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    while(tLevel>0){
        tParent = give_parent_of_element(tChild,aElementList.Dim,aElementList.NumElements); // Who is the parent of the child
        tChildren = give_children_of_element(tParent,aElementList.Dim,aElementList.NumElements); // Give the children of the parent element
        tBasis = give_basis_of_element(tParent,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements); //Basis of the parent element
        for(uint i = 0; i<tChildren.length(); i++)
        {
            if(tChildren(i)==tChild) //Give me the T-matrix of the child for the parent
            {
                tTmatrix_child = ttT_child(i);
                tT=tTmatrix_child*tT; //Create T-matrix for the parent element
                break;
            }
        }
        tVar = 0;
        tTMatrix.set_size(pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
        tTMatrixOfLevel.set_size(pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
        tMatrixActiveBasis.set_size(pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
        for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
        {
            if( (aBasisList.DesignBSplineActive).test(tBasis(i)) == 1 )
            {
                tTMatrix.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the parents
                tIdField(tVarb) = tBasis(i); //Add the row of the IdField of the active basis functions of the parents
                tTMatrixOfLevel.row(tVar) = tTmatrix_child.row(i); //Add the row of the Tmatrix of the active basis functions of the parents
                tMatrixActiveBasis.row(tVar) = tTeye.row(i);  // Add active basis functions with a one for their position
                tVar++;
                tVarb++;
            }
        }
        tLevel = tLevel - 1;
        tChild=tParent;

        tTMatrix.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
        tMatrixActiveBasis.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
        tTMatrixOfLevel.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
        tMatrixActiveBasisDummy = tMatrixActiveBasisEye(tVarcell-1);
        tTMatrixDummy = tTMatrixOnLevel(tVarcell-1);

        tTruncatedTMatrixOfLevel = tTMatrix - (( tTMatrixOfLevel * trans(tMatrixActiveBasisDummy) ) * tTMatrixDummy); // Take the T-Matrix of the current level and subrtract parts from the children level, which have an influence

        tTMatrixOnLevel(tVarcell) = tTruncatedTMatrixOfLevel;
        tMatrixActiveBasisEye(tVarcell) = tMatrixActiveBasis;
        tVarcell++;
        for(uint i = 0; i < tVarcell-2; i++)
        {
            tMatrixActiveBasis = tMatrixActiveBasisEye(i);
            tTMatrix = tTMatrixOnLevel(i);
            //            tTruncatedTMatrixOfLevel = tTruncatedTMatrixOfLevel - (tTruncatedTMatrixOfLevel * trans(tMatrixActiveBasis)) * tTMatrix; // Take the T-Matrix of the current level and subtract parts from all finer levels, which have an influence
            tTruncatedTMatrixOfLevel = tTruncatedTMatrixOfLevel - (tTruncatedTMatrixOfLevel * trans(tMatrixActiveBasis)) * tMatrixActiveBasis; // Take the T-Matrix of the current level and subtract parts from all finer levels, which have an influence

        }
        if( tVar > 0)
            aTMatrix.rows(tVarc,tVarb-1) = tTruncatedTMatrixOfLevel.rows(0,tVar-1);
        tVarc = tVarb;

    }
    //Resize the Tmatrix and the IdField
    aTMatrix.resize(tVarb,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
    tIdField.resize(tVarb,1);
    // Save data in a struc and multiply t-matrix with projection matrix
    aElementList.TMatrixDesign = aTMatrix * tT_Project;
    aElementList.IdFieldDesign = tIdField;
}

Mat<real>
moris::Hierarchical_Mesh::give_projection_matrix(
        ElementList_struc & aElementList)
{
    Mat<real> tT_Project; // Projection matrix
    if( aElementList.Dim == 2)
    {
        if( aElementList.PolynomialDesign == 1)
        {
            Mat<real> tTtemp = {{1, 0, 0, 0},{0, 1, 0 ,0},{0, 0, 1, 0},{0, 0, 0, 1}};
            tT_Project = tTtemp;
        }
        else if( aElementList.PolynomialDesign == 2)
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
        else if( aElementList.PolynomialDesign == 3)
        {
            Mat<real> tTtemp = {{    0.0278,         0,         0,         0},{0.1111,    0.0278,         0,         0},{0.0278,    0.1111,         0,         0},{     0,    0.0278,         0,         0},{0.1111,         0,    0.0278,         0},
                    {0.4444,    0.1111,    0.1111,    0.0278},{0.1111,    0.4444,    0.0278,    0.1111},{     0,    0.1111,         0,    0.0278},{0.0278,         0,    0.1111,         0},{0.1111,    0.0278,    0.4444,    0.1111},
                    {0.0278,    0.1111,    0.1111,    0.4444},{     0,    0.0278,         0,    0.1111},{     0,         0,    0.0278,         0},{     0,         0,    0.1111,    0.0278},{     0,         0,    0.0278,    0.1111},{     0,         0,         0,    0.0278}};
            tT_Project = tTtemp;
        }
        else
        {
            MORIS_LOG_ERROR << " Polynomial degree is not implemented";
        }
    }
    else if( aElementList.Dim == 3)
    {
        if( aElementList.PolynomialDesign == 1)
        {
            Mat<real> tTtemp = {{ 1,     0,     0,     0,     0,     0,     0,     0},
                    {0,     1,     0,     0,     0,     0,     0,     0},
                    {0,     0,     1,     0,     0,     0,     0,     0},
                    {0,     0,     0,     1,     0,     0,     0,     0},
                    {0,     0,     0,     0,     1,     0,     0,     0},
                    {0,     0,     0,     0,     0,     1,     0,     0},
                    {0,     0,     0,     0,     0,     0,     1,     0},
                    {0,     0,     0,     0,     0,     0,     0,     1}};
            tT_Project = tTtemp;
        }
        else if( aElementList.PolynomialDesign == 2)
        {
            Mat<real> tTtemp = {{    0.1250,         0,         0,         0,         0,         0,         0,         0},{0.1250,    0.1250,         0,         0,         0,         0,         0,         0},{     0,    0.1250,         0,         0,         0,         0,         0,         0},
                    {0.1250,         0,    0.1250,         0,         0,         0,         0,         0},{0.1250,    0.1250,    0.1250,    0.1250,         0,         0,         0,         0},{     0,    0.1250,         0,    0.1250,         0,         0,         0,         0},
                    {     0,         0,    0.1250,         0,         0,         0,         0,         0},{     0,         0,    0.1250,    0.1250,         0,         0,         0,         0},{     0,         0,         0,    0.1250,         0,         0,         0,         0},
                    {0.1250,         0,         0,         0,    0.1250,         0,         0,         0},{0.1250,    0.1250,         0,         0,    0.1250,    0.1250,         0,         0},{     0,    0.1250,         0,         0,         0,    0.1250,         0,         0},
                    {0.1250,         0,    0.1250,         0,    0.1250,         0,    0.1250,         0},{0.1250,    0.1250,    0.1250,    0.1250,    0.1250,    0.1250,    0.1250,    0.1250},{     0,    0.1250,         0,    0.1250,         0,    0.1250,         0,    0.1250},{     0,         0,    0.1250,         0,         0,         0,    0.1250,         0},
                    {     0,         0,    0.1250,    0.1250,         0,         0,    0.1250,    0.1250},{     0,         0,         0,    0.1250,         0,         0,         0,    0.1250},{     0,         0,         0,         0,    0.1250,         0,         0,         0},{     0,         0,         0,         0,    0.1250,    0.1250,         0,         0},
                    {     0,         0,         0,         0,         0,    0.1250,         0,         0},{     0,         0,         0,         0,    0.1250,         0,    0.1250,         0},{     0,         0,         0,         0,    0.1250,    0.1250,    0.1250,    0.1250},{     0,         0,         0,         0,         0,    0.1250,         0,    0.1250},
                    {     0,         0,         0,         0,         0,         0,    0.1250,         0},{     0,         0,         0,         0,         0,         0,    0.1250,    0.1250},{     0,         0,         0,         0,         0,         0,         0,    0.1250}};
            tT_Project = tTtemp;
        }
        else if( aElementList.PolynomialDesign == 3)
        {
            Mat<real> tTtemp = {{0.004629629629630,                   0,                   0,                   0,                   0,                   0,                   0,                   0 },{0.018518518518519,   0.004629629629630,                   0,                   0,                   0,                   0,                   0,                   0},{0.004629629629630,   0.018518518518519,                   0,                   0,                   0,                   0,                   0,                   0},
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

void moris::Hierarchical_Mesh::give_Tmatrix_and_id_field_for_designvariables_projection(uint & aElementNumber,
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tLevel = give_element_level(aElementNumber,aElementList.Dim,aElementList.NumElements)+1;// Level of the element
    Cell<Mat<real>> ttT_child = give_Tmatrices_of_childs(aElementList.PolynomialDesign,aElementList.Dim); // Give the T-matrices of the childs of the element "aElementNumber"
    Mat<real> tT = eye( pow(aElementList.PolynomialDesign+1,aElementList.Dim), pow(aElementList.PolynomialDesign+1,aElementList.Dim) ); //Tmatrix of the curent element
    Mat<real> tTMatrix(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0); //Tmatrix of the parent element
    Mat<uint> tIdField(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),1,0); //IdField of the parent Element
    Mat<uint> tBasis = give_basis_of_element(aElementNumber,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements); //Basis of the curent element
    uint tVar = 0; //Temporary variable for a loop

    for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
    {
        if( (aBasisList.DesignBSplineActiveLastStep).test(tBasis((i))) == 1 )
        {
            tTMatrix.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the curent element
            tIdField(tVar) = tBasis((i)); //Add the row of the IdField of the active basis functions of the curent element
            tVar++;
        }
    }
    tLevel = tLevel - 1; //Asking for T-matrix on level l-1 (one coarser level)
    uint tChild=aElementNumber;
    uint tParent;
    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    while(tLevel>0){
        tParent = give_parent_of_element(tChild,aElementList.Dim,aElementList.NumElements); // Who is the parent of the child
        tChildren = give_children_of_element(tParent,aElementList.Dim,aElementList.NumElements); // Give the children of the parent element
        tBasis = give_basis_of_element(tParent,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements); //Basis of the parent element
        for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
        {
            if(tChildren(i)==tChild) //Give me the T-matrix of the child for the parent
            {
                Mat<real> tTmatrix_child = ttT_child(i);
                tT=tTmatrix_child*tT; //Create T-matrix for the parent element
                break;
            }
        }
        for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
        {
            if( (aBasisList.DesignBSplineActiveLastStep).test(tBasis((i))) == 1 )
            {
                tTMatrix.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the parents
                tIdField(tVar) = tBasis((i)); //Add the row of the IdField of the active basis functions of the parents
                tVar++;
            }
        }
        tLevel = tLevel - 1;
        tChild=tParent;
    }
    //Resize the Tmatrix and the IdField
    tTMatrix.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
    tIdField.resize(tVar,1);
    // Save data in a struc
    aElementList.TMatrix = tTMatrix;
    aElementList.IdField = tIdField;
}

void moris::Hierarchical_Mesh::give_Truncated_Tmatrix_and_id_field_for_designvariables_projection(uint & aElementNumber,
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tLevel = give_element_level(aElementNumber,aElementList.Dim,aElementList.NumElements)+1;// Level of the element
    Cell<Mat<real>> ttT_child = give_Tmatrices_of_childs(aElementList.PolynomialDesign,aElementList.Dim); // Give the T-matrices of the childs of the element "aElementNumber"
    Mat<real> tT = eye( pow(aElementList.PolynomialDesign+1,aElementList.Dim), pow(aElementList.PolynomialDesign+1,aElementList.Dim) ); //Tmatrix of the curent element
    Mat<real> tTeye = eye( pow(aElementList.PolynomialDesign+1,aElementList.Dim), pow(aElementList.PolynomialDesign+1,aElementList.Dim) ); //Tmatrix of the curent element
    Mat<real> aTMatrix(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0); //Tmatrix of the parent element
    Mat<uint> tIdField(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),1,0); //IdField of the parent Element
    Mat<uint> tBasis = give_basis_of_element(aElementNumber,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements); //Basis of the curent element
    uint tVar = 0; uint tVarb = 0; uint tVarc = 0; uint tVarcell = 0; //Temporary variable for a loop
    Mat<real> tTmatrix_child; //Stores the tMatrix of a cell of the current child in a matrix
    Cell<Mat<real>> tMatrixActiveBasisEye(tLevel); // Stores a active basis functions for each level with the eye matrix
    Cell<Mat<real>> tTMatrixOnLevel(tLevel);//  Stores a T-matrices for each level
    Mat<real> tMatrixActiveBasis(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
    Mat<real> tTMatrixOfLevel(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
    Mat<real> tTMatrix(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
    Mat<real> tTMatrixDummy(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
    Mat<real> tTruncatedTMatrixOfLevel;
    Mat<real> tMatrixActiveBasisDummy;

    for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
    {
        if( (aBasisList.DesignBSplineActiveLastStep).test(tBasis((i))) == 1 )
        {
            aTMatrix.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the curent element
            tMatrixActiveBasis.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the curent element
            tTMatrixDummy.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the curent element
            tIdField(tVar) = tBasis((i)); //Add the row of the IdField of the active basis functions of the curent element
            tVar++;
        }
    }
    tVarb = tVar; // Save curent size of aIdField
    tVarc = tVar; // Save curent size of aTMatrix
    tTMatrixDummy.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
    tMatrixActiveBasis.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
    tMatrixActiveBasisEye(tVarcell) = tMatrixActiveBasis; // Save active basis functions with a one for their position
    tTMatrixOnLevel(tVarcell) = tTMatrixDummy; // Save the TMatrix of the current level
    tVarcell++;

    tLevel = tLevel - 1; //Asking for T-matrix on level l-1 (one coarser level)
    uint tChild=aElementNumber;
    uint tParent;
    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    while(tLevel>0){
        tParent = give_parent_of_element(tChild,aElementList.Dim,aElementList.NumElements); // Who is the parent of the child
        tChildren = give_children_of_element(tParent,aElementList.Dim,aElementList.NumElements); // Give the children of the parent element
        tBasis = give_basis_of_element(tParent,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements); //Basis of the parent element
        for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
        {
            if(tChildren(i)==tChild) //Give me the T-matrix of the child for the parent
            {
                tTmatrix_child = ttT_child(i);
                tT=tTmatrix_child*tT; //Create T-matrix for the parent element
                break;
            }
        }
        tVar = 0;
        tTMatrix.set_size(pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
        tTMatrixOfLevel.set_size(pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
        tMatrixActiveBasis.set_size(pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
        for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
        {
            if( (aBasisList.DesignBSplineActiveLastStep).test(tBasis((i))) == 1 )
            {
                tTMatrix.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the parents
                tIdField(tVarb) = tBasis(i); //Add the row of the IdField of the active basis functions of the parents
                tTMatrixOfLevel.row(tVar) = tTmatrix_child.row(i); //Add the row of the Tmatrix of the active basis functions of the parents
                tMatrixActiveBasis.row(tVar) = tTeye.row(i);  // Add active basis functions with a one for their position
                tVar++;
                tVarb++;
            }
        }
        tLevel = tLevel - 1;
        tChild=tParent;

        tTMatrix.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
        tMatrixActiveBasis.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
        tTMatrixOfLevel.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
        tMatrixActiveBasisDummy = tMatrixActiveBasisEye(tVarcell-1);
        tTMatrixDummy = tTMatrixOnLevel(tVarcell-1);

        tTruncatedTMatrixOfLevel = tTMatrix - (( tTMatrixOfLevel * trans(tMatrixActiveBasisDummy) ) * tTMatrixDummy); // Take the T-Matrix of the current level and subrtract parts from the children level, which have an influence

        tTMatrixOnLevel(tVarcell) = tTruncatedTMatrixOfLevel;
        tMatrixActiveBasisEye(tVarcell) = tMatrixActiveBasis;
        tVarcell++;
        for(uint i = 0; i < tVarcell-2; i++)
        {
            tMatrixActiveBasis = tMatrixActiveBasisEye(i);
            tTMatrix = tTMatrixOnLevel(i);
            //            tTruncatedTMatrixOfLevel = tTruncatedTMatrixOfLevel - (tTruncatedTMatrixOfLevel * trans(tMatrixActiveBasis)) * tTMatrix; // Take the T-Matrix of the current level and subtract parts from all finer levels, which have an influence
            tTruncatedTMatrixOfLevel = tTruncatedTMatrixOfLevel - (tTruncatedTMatrixOfLevel * trans(tMatrixActiveBasis)) * tMatrixActiveBasis; // Take the T-Matrix of the current level and subtract parts from all finer levels, which have an influence

        }
        if( tVar > 0)
            aTMatrix.rows(tVarc,tVarb-1) = tTruncatedTMatrixOfLevel.rows(0,tVar-1);
        tVarc = tVarb;
    }
    //Resize the Tmatrix and the IdField
    aTMatrix.resize(tVarb,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
    tIdField.resize(tVarb,1);
    // Save data in a struc
    aElementList.TMatrix = aTMatrix;
    aElementList.IdField = tIdField;
}
void moris::Hierarchical_Mesh::give_Tmatrix_and_id_field_for_designvariables_projection_coarsening(uint & aElementNumber,
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    Cell<Mat<real>> ttT_child = give_Tmatrices_of_childs(aElementList.PolynomialDesign,aElementList.Dim); // Give the T-matrices of the childs of the element "aElementNumber"
    Mat<real> tT = eye( pow(aElementList.PolynomialDesign+1,aElementList.Dim), pow(aElementList.PolynomialDesign+1,aElementList.Dim) ); //Tmatrix of the curent element
    Cell<Mat<real>> tTforchildren(pow(pow(2,aElementList.Dim),(aElementList.LevelLastStep)+1));
    Mat<uint> tPossibleChildrenList(pow(pow(2,aElementList.Dim),(aElementList.LevelLastStep)+1),1,0);
    Mat<uint> tActiveChildrenList(pow(pow(2,aElementList.Dim),(aElementList.LevelLastStep)+1),1,0);
    Mat<uint> tListActiveChildren(pow(pow(2,aElementList.Dim),(aElementList.LevelLastStep)+1),1,0);
    uint tVara = 0, tVarb = 1, tVarc = 1, tVard = 0;//Temporary variable for a loop
    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    tPossibleChildrenList(0) = aElementNumber;
    tTforchildren(0) = tT;
    uint tLevel = give_element_level(aElementNumber,aElementList.Dim,aElementList.NumElements);
    while( tPossibleChildrenList(tVara) > 0 && tLevel < (aElementList.LevelDesignLastStep))
    {
        tChildren = give_children_of_element(tPossibleChildrenList(tVara),aElementList.Dim,aElementList.NumElements); // Give the children of the parent element
        for(uint i = 0; i < tChildren.length(); i++)
        {
            tLevel = give_element_level(tChildren(i),aElementList.Dim,aElementList.NumElements);
            Mat<real> tTmatrix_child = ttT_child(i);
            tT = tTforchildren(tVara);
            tT=tTmatrix_child*tT; //Create T-matrix for the parent element
            tTforchildren(tVarc) = tT;
            if( (aElementList.ElementActiveLastStep).test(tChildren(i)) == 1 )
            {
                tActiveChildrenList(tVard) = tVarc;
                tListActiveChildren(tVard) = tChildren(i);
                tVard++;
            }
            else if( tLevel < (aElementList.LevelDesignLastStep) )
            {
                tPossibleChildrenList(tVarb) =  tChildren(i);
                tVarb++;
            }
            tVarc++;
        }
        tVara++;
    }
    tActiveChildrenList.resize(tVard,1);
    tListActiveChildren.resize(tVard,1);
    Cell<Mat<real>> tTMatrixActiveChildren(tVard);
    for(uint i = 0; i < tVard; i++)
    {
        tTMatrixActiveChildren(i) = tTforchildren(tActiveChildrenList(i));
    }
    (aElementList.ListOfChildren) = tListActiveChildren;
    (aElementList.tTMatrixOfChildren) = tTMatrixActiveChildren;
}

void moris::Hierarchical_Mesh::give_Truncated_Tmatrix_and_id_field_for_designvariables_projection_coarsening(uint & aElementNumber,
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    Cell<Mat<real>> ttT_child = give_Tmatrices_of_childs(aElementList.PolynomialDesign,aElementList.Dim); // Give the T-matrices of the childs of the element "aElementNumber"
    Mat<real> tT = eye( pow(aElementList.PolynomialDesign+1,aElementList.Dim), pow(aElementList.PolynomialDesign+1,aElementList.Dim) ); //Tmatrix of the curent element
    Cell<Mat<real>> tTforchildren(pow(pow(2,aElementList.Dim),(aElementList.LevelLastStep)+1));
    Mat<uint> tPossibleChildrenList(pow(pow(2,aElementList.Dim),(aElementList.LevelLastStep)+1),1,0);
    Mat<uint> tActiveChildrenList(pow(pow(2,aElementList.Dim),(aElementList.LevelLastStep)+1),1,0);
    Mat<uint> tListActiveChildren(pow(pow(2,aElementList.Dim),(aElementList.LevelLastStep)+1),1,0);
    uint tVara = 0, tVarb = 1, tVarc = 1, tVard = 0;//Temporary variable for a loop
    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    tPossibleChildrenList(0) = aElementNumber;
    tTforchildren(0) = tT;
    uint tLevel = give_element_level(aElementNumber,aElementList.Dim,aElementList.NumElements);
    while( tPossibleChildrenList(tVara) > 0 && tLevel < (aElementList.LevelDesignLastStep))
    {
        tChildren = give_children_of_element(tPossibleChildrenList(tVara),aElementList.Dim,aElementList.NumElements); // Give the children of the parent element
        for(uint i = 0; i < tChildren.length(); i++)
        {
            tLevel = give_element_level(tChildren(i),aElementList.Dim,aElementList.NumElements);
            Mat<real> tTmatrix_child = ttT_child(i);
            tT = tTforchildren(tVara);
            tT=tTmatrix_child*tT; //Create T-matrix for the parent element
            tTforchildren(tVarc) = tT;
            if( (aElementList.ElementActiveLastStep).test(tChildren(i)) == 1 )
            {
                tActiveChildrenList(tVard) = tVarc;
                tListActiveChildren(tVard) = tChildren(i);
                tVard++;
            }
            else if( tLevel < (aElementList.LevelDesignLastStep) )
            {
                tPossibleChildrenList(tVarb) =  tChildren(i);
                tVarb++;
            }
            tVarc++;
        }
        tVara++;
    }
    tActiveChildrenList.resize(tVard,1);
    tListActiveChildren.resize(tVard,1);
    Cell<Mat<real>> tTMatrixActiveChildren(tVard);
    for(uint i = 0; i < tVard; i++)
    {
        tTMatrixActiveChildren(i) = tTforchildren(tActiveChildrenList(i));
    }
    (aElementList.ListOfChildren) = tListActiveChildren;
    (aElementList.tTMatrixOfChildren) = tTMatrixActiveChildren;
}

void moris::Hierarchical_Mesh::give_Tmatrix_and_id_field_design(uint & aElementNumber,
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tLevel = give_element_level(aElementNumber,aElementList.Dim,aElementList.NumElements)+1;// Level of the element
    Cell<Mat<real>> ttT_child = give_Tmatrices_of_childs(aElementList.PolynomialDesign,aElementList.Dim); // Give the T-matrices of the childs of the element "aElementNumber"
    Mat<real> tT = eye( pow(aElementList.PolynomialDesign+1,aElementList.Dim), pow(aElementList.PolynomialDesign+1,aElementList.Dim) ); //Tmatrix of the curent element
    Mat<real> tTMatrix(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0); //Tmatrix of the parent element
    Mat<uint> tIdField(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),1,0); //IdField of the parent Element
    Mat<uint> tBasis = give_basis_of_element(aElementNumber,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements); //Basis of the curent element
    uint tVar = 0; //Temporary variable for a loop

    for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
    {
        if( (aBasisList.DesignBSplineActive).test(tBasis((i))) == 1 )
        {
            tTMatrix.row(tVar) = tT.row((i)); //Add the row of the Tmatrix of the active basis functions of the curent element
            tIdField(tVar) = tBasis((i)); //Add the row of the IdField of the active basis functions of the curent element
            tVar++;
        }
    }
    tLevel = tLevel - 1;
    //    //Asking for T-matrix on level l-1 (one coarser level)
    uint tChild=aElementNumber;
    uint tParent;
    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    while(tLevel>0){
        tParent = give_parent_of_element(tChild,aElementList.Dim,aElementList.NumElements); // Who is the parent of the child
        tChildren = give_children_of_element(tParent,aElementList.Dim,aElementList.NumElements); // Give the children of the parent element
        tBasis = give_basis_of_element(tParent,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements); //Basis of the parent element
        for(uint i = 0; i<tChildren.length(); i++)
        {
            if(tChildren(i)==tChild) //Give me the T-matrix of the child for the parent
            {
                Mat<real> tTmatrix_child = ttT_child(i);
                tT=tTmatrix_child*tT; //Create T-matrix for the parent element
                break;
            }
        }
        for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
        {
            if( (aBasisList.DesignBSplineActive).test(tBasis((i))) == 1 )
            {
                tTMatrix.row(tVar) = tT.row((i)); //Add the row of the Tmatrix of the active basis functions of the parents
                tIdField(tVar) = tBasis((i)); //Add the row of the IdField of the active basis functions of the parents
                tVar++;
            }
        }
        tLevel = tLevel - 1;
        tChild=tParent;
    }
    //Resize the Tmatrix and the IdField
    tTMatrix.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
    tIdField.resize(tVar,1);
    // Save data in a struc
    aElementList.TMatrix = tTMatrix;
    aElementList.IdField = tIdField;
}

void moris::Hierarchical_Mesh::give_Truncated_Tmatrix_and_id_field_design(uint & aElementNumber,
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tLevel = give_element_level(aElementNumber,aElementList.Dim,aElementList.NumElements)+1;// Level of the element
    Cell<Mat<real>> ttT_child = give_Tmatrices_of_childs(aElementList.PolynomialDesign,aElementList.Dim); // Give the T-matrices of the childs of the element "aElementNumber"
    Mat<real> tT = eye( pow(aElementList.PolynomialDesign+1,aElementList.Dim), pow(aElementList.PolynomialDesign+1,aElementList.Dim) ); //Tmatrix of the curent element
    Mat<real> tTeye = eye( pow(aElementList.PolynomialDesign+1,aElementList.Dim), pow(aElementList.PolynomialDesign+1,aElementList.Dim) ); //Tmatrix of the curent element
    Mat<real> aTMatrix(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0); //Tmatrix of the parent element
    Mat<uint> tIdField(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),1,0); //IdField of the parent Element
    Mat<uint> tBasis = give_basis_of_element(aElementNumber,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements); //Basis of the curent element
    uint tVar = 0; uint tVarb = 0; uint tVarc = 0; uint tVarcell = 0; //Temporary variable for a loop
    Mat<real> tTmatrix_child; //Stores the tMatrix of a cell of the current child in a matrix
    Cell<Mat<real>> tMatrixActiveBasisEye(tLevel); // Stores a active basis functions for each level with the eye matrix
    Cell<Mat<real>> tTMatrixOnLevel(tLevel);//  Stores a T-matrices for each level
    Mat<real> tMatrixActiveBasis(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
    Mat<real> tTMatrixOfLevel(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
    Mat<real> tTMatrix(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
    Mat<real> tTMatrixDummy(tLevel*pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
    Mat<real> tTruncatedTMatrixOfLevel;
    Mat<real> tMatrixActiveBasisDummy;

    for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
    {
        if( (aBasisList.DesignBSplineActive).test(tBasis((i))) == 1 )
        {
            aTMatrix.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the curent element
            tIdField(tVar) = tBasis(i); //Add the row of the IdField of the active basis functions of the curent element
            tMatrixActiveBasis.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the curent element
            tTMatrixDummy.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the curent element
            tVar++;
        }
    }
    tVarb = tVar; // Save curent size of aIdField
    tVarc = tVar; // Save curent size of aTMatrix
    tTMatrixDummy.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
    tMatrixActiveBasis.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
    tMatrixActiveBasisEye(tVarcell) = tMatrixActiveBasis; // Save active basis functions with a one for their position
    tTMatrixOnLevel(tVarcell) = tTMatrixDummy; // Save the TMatrix of the current level
    tVarcell++;

    tLevel = tLevel - 1;
    //    //Asking for T-matrix on level l-1 (one coarser level)
    uint tChild=aElementNumber;
    uint tParent;
    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    while(tLevel>0){
        tParent = give_parent_of_element(tChild,aElementList.Dim,aElementList.NumElements); // Who is the parent of the child
        tChildren = give_children_of_element(tParent,aElementList.Dim,aElementList.NumElements); // Give the children of the parent element
        tBasis = give_basis_of_element(tParent,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements); //Basis of the parent element
        for(uint i = 0; i<tChildren.length(); i++)
        {
            if(tChildren(i)==tChild) //Give me the T-matrix of the child for the parent
            {
                tTmatrix_child = ttT_child(i);
                tT=tTmatrix_child*tT; //Create T-matrix for the parent element
                break;
            }
        }
        tVar = 0;
        tTMatrix.set_size(pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
        tTMatrixOfLevel.set_size(pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
        tMatrixActiveBasis.set_size(pow(aElementList.PolynomialDesign+1,aElementList.Dim),pow(aElementList.PolynomialDesign+1,aElementList.Dim),0);
        for(uint i = 0; i<pow(aElementList.PolynomialDesign+1,aElementList.Dim); i++)
        {
            if( (aBasisList.DesignBSplineActive).test(tBasis((i))) == 1 )
            {
                tTMatrix.row(tVar) = tT.row(i); //Add the row of the Tmatrix of the active basis functions of the parents
                tIdField(tVarb) = tBasis(i); //Add the row of the IdField of the active basis functions of the parents
                tTMatrixOfLevel.row(tVar) = tTmatrix_child.row(i); //Add the row of the Tmatrix of the active basis functions of the parents
                tMatrixActiveBasis.row(tVar) = tTeye.row(i);  // Add active basis functions with a one for their position
                tVar++;
                tVarb++;
            }
        }
        tLevel = tLevel - 1;
        tChild=tParent;

        tTMatrix.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
        tMatrixActiveBasis.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
        tTMatrixOfLevel.resize(tVar,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
        tMatrixActiveBasisDummy = tMatrixActiveBasisEye(tVarcell-1);
        tTMatrixDummy = tTMatrixOnLevel(tVarcell-1);

        tTruncatedTMatrixOfLevel = tTMatrix - (( tTMatrixOfLevel * trans(tMatrixActiveBasisDummy) ) * tTMatrixDummy); // Take the T-Matrix of the current level and subrtract parts from the children level, which have an influence

        tTMatrixOnLevel(tVarcell) = tTruncatedTMatrixOfLevel;
        tMatrixActiveBasisEye(tVarcell) = tMatrixActiveBasis;
        tVarcell++;
        for(uint i = 0; i < tVarcell-2; i++)
        {
            tMatrixActiveBasis = tMatrixActiveBasisEye(i);
            tTMatrix = tTMatrixOnLevel(i);
            //            tTruncatedTMatrixOfLevel = tTruncatedTMatrixOfLevel - (tTruncatedTMatrixOfLevel * trans(tMatrixActiveBasis)) * tTMatrix; // Take the T-Matrix of the current level and subtract parts from all finer levels, which have an influence
            tTruncatedTMatrixOfLevel = tTruncatedTMatrixOfLevel - (tTruncatedTMatrixOfLevel * trans(tMatrixActiveBasis)) * tMatrixActiveBasis; // Take the T-Matrix of the current level and subtract parts from all finer levels, which have an influence

        }
        if( tVar > 0)
            aTMatrix.rows(tVarc,tVarb-1) = tTruncatedTMatrixOfLevel.rows(0,tVar-1);
        tVarc = tVarb;
    }
    //Resize the Tmatrix and the IdField
    aTMatrix.resize(tVarb,pow(aElementList.PolynomialDesign+1,aElementList.Dim));
    tIdField.resize(tVarb,1);
    // Save data in a struc
    aElementList.TMatrix = aTMatrix;
    aElementList.IdField = tIdField;
}

uint
moris::Hierarchical_Mesh::give_element_of_position(uint & aLevel,
        uint & aDim,
        Mat<uint> & aNumElements,
        Mat<uint> & aIJKPosition)
{
    uint tElement=0;
    uint tLevel = aLevel - 1;
    uint tY_Position; // Position of the Element in y-direction
    if( aDim == 2 )
    {
        tY_Position = aIJKPosition(1)*aNumElements(0);
        if( aLevel>0 )
        {
            tElement = give_number_of_elements(tLevel,aDim,aNumElements);//Give number of elements for the last level to have the initial number at position 0,0
            tY_Position = aIJKPosition(1)*aNumElements(0)*pow(2,aLevel);  // Position of the Element in y-direction for a specific level
        }
        tElement = tElement + aIJKPosition(0) + tY_Position;
    }
    else if(aDim == 3)
    {
        tY_Position = aIJKPosition(1)*aNumElements(0)+aIJKPosition(2)*aNumElements(0)*aNumElements(1);
        if( aLevel>0 )
        {
            tElement = give_number_of_elements(tLevel,aDim,aNumElements);//Give number of elements for the last level to have the initial number at position 0,0,0
            tY_Position = aIJKPosition(1)*aNumElements(0)*pow(2,aLevel)+aIJKPosition(2)*aNumElements(0)*pow(2,aLevel)*aNumElements(1)*pow(2,aLevel);  // Position of the Element in y-direction for a specific level
        }
        tElement = tElement + aIJKPosition(0) + tY_Position;
    }
    return tElement;
}

Mat<uint>
moris::Hierarchical_Mesh::give_position_of_element(uint & aElement,
        uint & aDim,
        Mat<uint> & aNumElements)
{
    Mat<uint> tElement_position(aDim,1); // Element positions in the direction x,y,z
    uint tElement_level=give_element_level(aElement,aDim,aNumElements);
    uint tNumber_of_elemente=0; // temporary variable for the number of elements per level
    if( aDim == 2 )
    {
        if(tElement_level !=0)
        {
            tElement_level -= 1; // Count the elements until the level tElement_level-1
            tNumber_of_elemente = give_number_of_elements(tElement_level,aDim,aNumElements);
            tElement_level += 1; // +1 to give the original level of elment aElement
        }
        tElement_position(1) = ceil((aElement+1-tNumber_of_elemente)/(pow(2,tElement_level)*aNumElements(0))); // Position is id based ( starting with 1,2,...)
        tElement_position(0) = aElement+1 - tNumber_of_elemente + pow(2,tElement_level)*aNumElements(0) - tElement_position(1)*pow(2,tElement_level)*aNumElements(0); // Position is id based ( starting with 1,2,...)
        tElement_position(1) -= 1; // Position is index based ( starting with 0,1,2,...)
        tElement_position(0) -= 1; // Position is index based ( starting with 0,1,2,...)

    }
    else if(aDim == 3)
    {
        if(tElement_level !=0)
        {
            tElement_level -= 1; // Count the elements until the level tElement_level-1
            tNumber_of_elemente = give_number_of_elements(tElement_level,aDim,aNumElements);
            tElement_level += 1; // +1 to give the original level of elment aElement
        }
        tElement_position(2) = ceil((aElement+1-tNumber_of_elemente)/(pow(2,tElement_level)*aNumElements(0)*pow(2,tElement_level)*aNumElements(1))); // Position is id based ( starting with 1,2,...)
        tElement_position(1) = ceil((aElement+1-tNumber_of_elemente-(tElement_position(2)-1)*pow(2,tElement_level)*aNumElements(0)*pow(2,tElement_level)*aNumElements(1))/(pow(2,tElement_level)*aNumElements(0))); // Position is id based ( starting with 1,2,...)
        tElement_position(0) = aElement+1 - tNumber_of_elemente + pow(2,tElement_level)*aNumElements(0) - tElement_position(1)*pow(2,tElement_level)*aNumElements(0) + pow(2,tElement_level)*aNumElements(0)*pow(2,tElement_level)*aNumElements(1) - tElement_position(2)*pow(2,tElement_level)*aNumElements(0)*pow(2,tElement_level)*aNumElements(1); // Position is id based ( starting with 1,2,...)
        tElement_position(2) -= 1; // Position is index based ( starting with 0,1,2,...)
        tElement_position(1) -= 1; // Position is index based ( starting with 0,1,2,...)
        tElement_position(0) -= 1; // Position is index based ( starting with 0,1,2,...)
    }
    return tElement_position;
}

uint
moris::Hierarchical_Mesh::give_element_level(uint & aElement,
        uint & aDim,
        Mat<uint> & aNumElements)
{
    uint tElement_level=1; // output variable
    uint tLevel=1; // temporary variable for the while loop
    if( aDim == 2 )
    {
        uint tNumber_of_elemente=pow(2,tElement_level-1)*aNumElements(0)*pow(2,tElement_level-1)*aNumElements(1); // temporary variable for the number of elements per level
        while(tLevel>0)
        {
            tLevel=floor((real)(aElement+1)/tNumber_of_elemente);
            if( (real)(aElement+1)/tNumber_of_elemente == 1 || tLevel<1)
            {
                tLevel = 0;
            }
            else
            {
                tElement_level += 1;
                tNumber_of_elemente += pow(2,tElement_level-1)*aNumElements(0)*pow(2,tElement_level-1)*aNumElements(1);
            }
        }
        tElement_level -= 1;
    }
    else if(aDim == 3)
    {
        uint tNumber_of_elemente=pow(2,tElement_level-1)*aNumElements(0)*pow(2,tElement_level-1)*aNumElements(1)*pow(2,tElement_level-1)*aNumElements(2); // temporary variable for the number of elements per level
        while(tLevel>0)
        {
            tLevel=floor((real)(aElement+1)/tNumber_of_elemente);
            if( (real)(aElement+1)/tNumber_of_elemente == 1 || tLevel<1)
            {
                tLevel = 0;
            }
            else
            {
                tElement_level += 1;
                tNumber_of_elemente += pow(2,tElement_level-1)*aNumElements(0)*pow(2,tElement_level-1)*aNumElements(1)*pow(2,tElement_level-1)*aNumElements(2);
            }
        }
        tElement_level -= 1;
    }
    return tElement_level;
}

uint
moris::Hierarchical_Mesh::give_number_of_elements(uint & aLevel,
        uint & aDim,
        Mat<uint> & aNumElements)
{
    uint tElement_number=0; // output variable
    if( aDim == 2 )
    {
        for(uint i = 0; i<aLevel+1; i++)
        {
            tElement_number += pow(2,i)*aNumElements(0)*pow(2,i)*aNumElements(1);
        }
    }
    else if(aDim == 3)
    {
        for(uint i = 0; i<aLevel+1; i++)
        {
            tElement_number += pow(2,i)*aNumElements(0)*pow(2,i)*aNumElements(1)*pow(2,i)*aNumElements(2);
        }
    }
    return tElement_number;
}

uint
moris::Hierarchical_Mesh::give_parent_of_element(uint & aElement,
        uint & aDim,
        Mat<uint> & aNumElements)
{
    Mat<uint> tElement_position=give_position_of_element(aElement,aDim,aNumElements);
    uint tElement_level=give_element_level(aElement,aDim,aNumElements);
    uint tParent = UINT_MAX;
    uint tParent_level = 0;
    if(tElement_level>0)
    {
        tParent_level = tElement_level - 1;
        Mat<uint> tParent_position(aDim,1); // Position of the Parent element
        if( aDim == 2 )
        {
            if(((tElement_position(0)+1)% 2) == 0 && ((tElement_position(1)+1) % 2) == 0)
            {
                tParent_position(0) = (tElement_position(0)+1)/2-1;
                tParent_position(1) = (tElement_position(1)+1)/2-1;
            }
            else if(((tElement_position(0)+1) % 2) == 0 && ((tElement_position(1)+1) % 2) != 0)
            {
                tParent_position(0) = (tElement_position(0)+1)/2-1;
                tParent_position(1) = (tElement_position(1)+1+1)/2-1;
            }
            else if(((tElement_position(0)+1) % 2) != 0 && ((tElement_position(1)+1) % 2) == 0)
            {
                tParent_position(0) = (tElement_position(0)+1+1)/2-1;
                tParent_position(1) = (tElement_position(1)+1)/2-1;
            }
            else if(((tElement_position(0)+1) % 2) != 0 && ((tElement_position(1)+1) % 2) != 0)
            {
                tParent_position(0) = (tElement_position(0)+1+1)/2-1;
                tParent_position(1) = (tElement_position(1)+1+1)/2-1;
            }
        }
        else if(aDim == 3)
        {
            if(((tElement_position(0)+1) % 2) == 0 && ((tElement_position(1)+1) % 2) == 0 && ((tElement_position(2)+1) % 2) == 0)
            {
                tParent_position(0) = (tElement_position(0)+1)/2-1;
                tParent_position(1) = (tElement_position(1)+1)/2-1;
                tParent_position(2) = (tElement_position(2)+1)/2-1;
            }
            else if(((tElement_position(0)+1) % 2) == 0 && ((tElement_position(1)+1) % 2) != 0 && ((tElement_position(2)+1) % 2) == 0)
            {
                tParent_position(0) = (tElement_position(0)+1)/2-1;
                tParent_position(1) = (tElement_position(1)+1+1)/2-1;
                tParent_position(2) = (tElement_position(2)+1)/2-1;
            }
            else if(((tElement_position(0)+1) % 2) != 0 && ((tElement_position(1)+1) % 2) == 0 && ((tElement_position(2)+1) % 2) == 0)
            {
                tParent_position(0) = (tElement_position(0)+1+1)/2-1;
                tParent_position(1) = (tElement_position(1)+1)/2-1;
                tParent_position(2) = (tElement_position(2)+1)/2-1;
            }
            else if(((tElement_position(0)+1) % 2) != 0 && ((tElement_position(1)+1) % 2) != 0 && ((tElement_position(2)+1) % 2) == 0)
            {
                tParent_position(0) = (tElement_position(0)+1+1)/2-1;
                tParent_position(1) = (tElement_position(1)+1+1)/2-1;
                tParent_position(2) = (tElement_position(2)+1)/2-1;
            }
            else if(((tElement_position(0)+1) % 2) == 0 && ((tElement_position(1)+1) % 2) == 0 && ((tElement_position(2)+1) % 2) != 0)
            {
                tParent_position(0) = (tElement_position(0)+1)/2-1;
                tParent_position(1) = (tElement_position(1)+1)/2-1;
                tParent_position(2) = (tElement_position(2)+1+1)/2-1;
            }
            else if(((tElement_position(0)+1) % 2) == 0 && ((tElement_position(1)+1) % 2) != 0 && ((tElement_position(2)+1) % 2) != 0)
            {
                tParent_position(0) = (tElement_position(0)+1)/2-1;
                tParent_position(1) = (tElement_position(1)+1+1)/2-1;
                tParent_position(2) = (tElement_position(2)+1+1)/2-1;
            }
            else if(((tElement_position(0)+1) % 2) != 0 && ((tElement_position(1)+1) % 2) == 0 && ((tElement_position(2)+1) % 2) != 0)
            {
                tParent_position(0) = (tElement_position(0)+1+1)/2-1;
                tParent_position(1) = (tElement_position(1)+1)/2-1;
                tParent_position(2) = (tElement_position(2)+1+1)/2-1;
            }
            else if(((tElement_position(0)+1) % 2) != 0 && ((tElement_position(1)+1) % 2) != 0 && ((tElement_position(2)+1) % 2) != 0)
            {
                tParent_position(0) = (tElement_position(0)+1+1)/2-1;
                tParent_position(1) = (tElement_position(1)+1+1)/2-1;
                tParent_position(2) = (tElement_position(2)+1+1)/2-1;
            }
        }
        tParent=give_element_of_position(tParent_level,aDim,aNumElements,tParent_position);
    }
    return tParent;
}

uint
moris::Hierarchical_Mesh::give_parent_of_level_x(uint & aElement,
        uint & aDim,
        Mat<uint> & aNumElements,
        uint & aWhichLevel)
{
    uint tParent = give_parent_of_element(aElement,aDim,aNumElements);
    uint tLevel;
    if( tParent < UINT_MAX)
    {
        tParent = aElement;
        tLevel = give_element_level(tParent,aDim,aNumElements);
        while( tLevel > aWhichLevel)
        {
            tParent = give_parent_of_element(tParent,aDim,aNumElements);
            tLevel = give_element_level(tParent,aDim,aNumElements);
        }
    }
    return tParent;
}

moris::Mat<uint>
moris::Hierarchical_Mesh::give_children_of_element(uint & aElement,
        uint & aDim,
        Mat<uint> & aNumElements)
{
    Mat<uint> tElement_position=give_position_of_element(aElement,aDim,aNumElements);
    uint tChildren_level=give_element_level(aElement,aDim,aNumElements)+1;
    Mat<uint> tChildren(pow(2,aDim),1); // Vector with all children of the element aElement
    Mat<uint> tChildren_position(aDim,1); // Position of the children with respect to the element aElement
    if(aDim == 2)
    {
        tChildren_position(0) = 2 * (tElement_position(0)+1) - 1 - 1;   tChildren_position(1) = 2 * (tElement_position(1)+1) - 1 - 1;
        (tChildren(0))=give_element_of_position(tChildren_level,aDim,aNumElements,tChildren_position);
        tChildren_position(0) = 2 * (tElement_position(0)+1) - 1;       tChildren_position(1) = 2 * (tElement_position(1)+1) - 1 - 1;
        (tChildren(1))=give_element_of_position(tChildren_level,aDim,aNumElements,tChildren_position);
        tChildren_position(0) = 2 * (tElement_position(0)+1) - 1 - 1;   tChildren_position(1) = 2 * (tElement_position(1)+1) - 1;
        (tChildren(2))=give_element_of_position(tChildren_level,aDim,aNumElements,tChildren_position);
        tChildren_position(0) = 2 * (tElement_position(0)+1) - 1;       tChildren_position(1) = 2 * (tElement_position(1)+1) - 1;
        (tChildren(3))=give_element_of_position(tChildren_level,aDim,aNumElements,tChildren_position);
    }
    else if(aDim == 3)
    {
        tChildren_position(0) = 2 * (tElement_position(0)+1) - 1 - 1;   tChildren_position(1) = 2 * (tElement_position(1)+1) - 1 - 1;   tChildren_position(2) = 2 * (tElement_position(2)+1) - 1 - 1;
        (tChildren(0))=give_element_of_position(tChildren_level,aDim,aNumElements,tChildren_position);
        tChildren_position(0) = 2 * (tElement_position(0)+1) - 1;       tChildren_position(1) = 2 * (tElement_position(1)+1) - 1 - 1;   tChildren_position(2) = 2 * (tElement_position(2)+1) - 1 - 1;
        (tChildren(1))=give_element_of_position(tChildren_level,aDim,aNumElements,tChildren_position);
        tChildren_position(0) = 2 * (tElement_position(0)+1) - 1 - 1;   tChildren_position(1) = 2 * (tElement_position(1)+1) - 1;       tChildren_position(2) = 2 * (tElement_position(2)+1) - 1 - 1;
        (tChildren(2))=give_element_of_position(tChildren_level,aDim,aNumElements,tChildren_position);
        tChildren_position(0) = 2 * (tElement_position(0)+1) - 1;       tChildren_position(1) = 2 * (tElement_position(1)+1) - 1;       tChildren_position(2) = 2 * (tElement_position(2)+1) - 1 - 1;
        (tChildren(3))=give_element_of_position(tChildren_level,aDim,aNumElements,tChildren_position);
        tChildren_position(0) = 2 * (tElement_position(0)+1) - 1 - 1;   tChildren_position(1) = 2 * (tElement_position(1)+1) - 1 - 1;   tChildren_position(2) = 2 * (tElement_position(2)+1) - 1;
        (tChildren(4))=give_element_of_position(tChildren_level,aDim,aNumElements,tChildren_position);
        tChildren_position(0) = 2 * (tElement_position(0)+1) - 1;       tChildren_position(1) = 2 * (tElement_position(1)+1) - 1 - 1;   tChildren_position(2) = 2 * (tElement_position(2)+1) - 1;
        (tChildren(5))=give_element_of_position(tChildren_level,aDim,aNumElements,tChildren_position);
        tChildren_position(0) = 2 * (tElement_position(0)+1) - 1 - 1;   tChildren_position(1) = 2 * (tElement_position(1)+1) - 1;       tChildren_position(2) = 2 * (tElement_position(2)+1) - 1;
        (tChildren(6))=give_element_of_position(tChildren_level,aDim,aNumElements,tChildren_position);
        tChildren_position(0) = 2 * (tElement_position(0)+1) - 1;       tChildren_position(1) = 2 * (tElement_position(1)+1) - 1;       tChildren_position(2) = 2 * (tElement_position(2)+1) - 1;
        (tChildren(7))=give_element_of_position(tChildren_level,aDim,aNumElements,tChildren_position);
    }
    return tChildren;
}
moris::Mat<uint>
moris::Hierarchical_Mesh::give_neighbour_of_element(uint & aElement,
        uint & aDim,
        uint & aBuffer,
        Mat<uint> & aNumElements)
{
    Mat<uint> tElement_position=give_position_of_element(aElement,aDim,aNumElements);
    uint tElement_level=give_element_level(aElement,aDim,aNumElements);
    Mat<uint> Element_neighbour(1,pow(1+2*aBuffer,aDim),UINT_MAX); // Element plus neighbours (9 for aDim = 2, 27 for aDim = 3)
    Mat<uint> tNeighbour_position(1,aDim); // Position of a neighbour element;
    if(aDim == 2)
    {
        if( aBuffer == 1)
        {
            tNeighbour_position(0) = tElement_position(0) - 1;   tNeighbour_position(1) = tElement_position(1) - 1;
            (Element_neighbour(0))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0);       tNeighbour_position(1) = tElement_position(1) - 1;
            (Element_neighbour(1))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) + 1;   tNeighbour_position(1) = tElement_position(1) - 1;
            (Element_neighbour(2))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) - 1;   tNeighbour_position(1) = tElement_position(1);
            (Element_neighbour(3))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0);       tNeighbour_position(1) = tElement_position(1);
            (Element_neighbour(4))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) + 1;   tNeighbour_position(1) = tElement_position(1);
            (Element_neighbour(5))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) - 1;   tNeighbour_position(1) = tElement_position(1) + 1;
            (Element_neighbour(6))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0);       tNeighbour_position(1) = tElement_position(1) + 1;
            (Element_neighbour(7))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) + 1;   tNeighbour_position(1) = tElement_position(1) + 1;
            (Element_neighbour(8))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
        }
        else if( aBuffer > 1) // If a higher buffer layer is needed, a loop generates the neighbours
        {
            uint tVar = 0;
            Mat<uint> aIJKPosition(aDim,1,0);
            for(int j = -((int)aBuffer); j <= ((int)aBuffer); j++ )
            {
                for(int i = -((int)aBuffer); i <= ((int)aBuffer); i++ )
                {
                    if ( tElement_position(0)+i >= 1*pow(2,tElement_level) && tElement_position(0)+i < aNumElements(0)*pow(2,tElement_level)-1*pow(2,tElement_level) && tElement_position(1)+j >= 1*pow(2,tElement_level) && tElement_position(1)+j < aNumElements(1)*pow(2,tElement_level)-1*pow(2,tElement_level) )
                    {
                        aIJKPosition(0) = tElement_position(0) + i;    aIJKPosition(1) = tElement_position(1) + j;
                        Element_neighbour(tVar)=give_element_of_position(tElement_level,aDim,aNumElements,aIJKPosition);
                    }
                    tVar++;
                }
            }
        }
    }
    else if(aDim == 3)
    {
        if( aBuffer == 1)
        {
            tNeighbour_position(0) = tElement_position(0) - 1;   tNeighbour_position(1) = tElement_position(1) - 1;   tNeighbour_position(2) = tElement_position(2) - 1;
            (Element_neighbour(0))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0);       tNeighbour_position(1) = tElement_position(1) - 1;   tNeighbour_position(2) = tElement_position(2) - 1;
            (Element_neighbour(1))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) + 1;   tNeighbour_position(1) = tElement_position(1) - 1;   tNeighbour_position(2) = tElement_position(2) - 1;
            (Element_neighbour(2))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) - 1;   tNeighbour_position(1) = tElement_position(1);       tNeighbour_position(2) = tElement_position(2) - 1;
            (Element_neighbour(3))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0);       tNeighbour_position(1) = tElement_position(1);       tNeighbour_position(2) = tElement_position(2) - 1;
            (Element_neighbour(4))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) + 1;   tNeighbour_position(1) = tElement_position(1);       tNeighbour_position(2) = tElement_position(2) - 1;
            (Element_neighbour(5))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) - 1;   tNeighbour_position(1) = tElement_position(1) + 1;   tNeighbour_position(2) = tElement_position(2) - 1;
            (Element_neighbour(6))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0);       tNeighbour_position(1) = tElement_position(1) + 1;   tNeighbour_position(2) = tElement_position(2) - 1;
            (Element_neighbour(7))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) + 1;   tNeighbour_position(1) = tElement_position(1) + 1;   tNeighbour_position(2) = tElement_position(2) - 1;
            (Element_neighbour(8))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);

            tNeighbour_position(0) = tElement_position(0) - 1;   tNeighbour_position(1) = tElement_position(1) - 1;   tNeighbour_position(2) = tElement_position(2);
            (Element_neighbour(9))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0);       tNeighbour_position(1) = tElement_position(1) - 1;   tNeighbour_position(2) = tElement_position(2);
            (Element_neighbour(10))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) + 1;   tNeighbour_position(1) = tElement_position(1) - 1;   tNeighbour_position(2) = tElement_position(2);
            (Element_neighbour(11))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) - 1;   tNeighbour_position(1) = tElement_position(1);       tNeighbour_position(2) = tElement_position(2);
            (Element_neighbour(12))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0);       tNeighbour_position(1) = tElement_position(1);       tNeighbour_position(2) = tElement_position(2);
            (Element_neighbour(13))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) + 1;   tNeighbour_position(1) = tElement_position(1);       tNeighbour_position(2) = tElement_position(2);
            (Element_neighbour(14))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) - 1;   tNeighbour_position(1) = tElement_position(1) + 1;   tNeighbour_position(2) = tElement_position(2);
            (Element_neighbour(15))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0);       tNeighbour_position(1) = tElement_position(1) + 1;   tNeighbour_position(2) = tElement_position(2);
            (Element_neighbour(16))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) + 1;   tNeighbour_position(1) = tElement_position(1) + 1;   tNeighbour_position(2) = tElement_position(2);
            (Element_neighbour(17))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);

            tNeighbour_position(0) = tElement_position(0) - 1;   tNeighbour_position(1) = tElement_position(1) - 1;   tNeighbour_position(2) = tElement_position(2) + 1;
            (Element_neighbour(18))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0);       tNeighbour_position(1) = tElement_position(1) - 1;   tNeighbour_position(2) = tElement_position(2) + 1;
            (Element_neighbour(19))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) + 1;   tNeighbour_position(1) = tElement_position(1) - 1;   tNeighbour_position(2) = tElement_position(2) + 1;
            (Element_neighbour(20))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) - 1;   tNeighbour_position(1) = tElement_position(1);       tNeighbour_position(2) = tElement_position(2) + 1;
            (Element_neighbour(21))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0);       tNeighbour_position(1) = tElement_position(1);       tNeighbour_position(2) = tElement_position(2) + 1;
            (Element_neighbour(22))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) + 1;   tNeighbour_position(1) = tElement_position(1);       tNeighbour_position(2) = tElement_position(2) + 1;
            (Element_neighbour(23))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) - 1;   tNeighbour_position(1) = tElement_position(1) + 1;   tNeighbour_position(2) = tElement_position(2) + 1;
            (Element_neighbour(24))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0);       tNeighbour_position(1) = tElement_position(1) + 1;   tNeighbour_position(2) = tElement_position(2) + 1;
            (Element_neighbour(25))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tElement_position(0) + 1;   tNeighbour_position(1) = tElement_position(1) + 1;   tNeighbour_position(2) = tElement_position(2) + 1;
            (Element_neighbour(26))=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
        }
        else if( aBuffer > 1) // If a higher buffer layer is needed, a loop generates the neighbours
        {
            uint tVar = 0;
            Mat<uint> aIJKPosition(aDim,1,0);
            for(int k = -((int)aBuffer); k <= ((int)aBuffer); k++ )
            {
                for(int j = -((int)aBuffer); j <= ((int)aBuffer); j++ )
                {
                    for(int i = -((int)aBuffer); i <= ((int)aBuffer); i++ )
                    {
                        if ( tElement_position(0)+i >= 1*pow(2,tElement_level) && tElement_position(0)+i < aNumElements(0)*pow(2,tElement_level)-1*pow(2,tElement_level) && tElement_position(1)+j >= 1*pow(2,tElement_level) && tElement_position(1)+j < aNumElements(1)*pow(2,tElement_level)-1*pow(2,tElement_level) && tElement_position(2)+k >= 1*pow(2,tElement_level) && tElement_position(2)+k < aNumElements(2)*pow(2,tElement_level)-1*pow(2,tElement_level) )
                        {
                            aIJKPosition(0) = tElement_position(0) + i; aIJKPosition(1) = tElement_position(1) + j; aIJKPosition(2) = tElement_position(2) + k;
                            Element_neighbour(tVar)=give_element_of_position(tElement_level,aDim,aNumElements,aIJKPosition);
                        }
                        tVar++;
                    }
                }
            }
        }
    }
    return Element_neighbour;
}

Mat<uint>
moris::Hierarchical_Mesh::give_active_neighbour_of_element(uint & aElement,
        ElementList_struc & aElementList)
{
    uint tBuffer = 1;
    Mat<uint> tNeighbour = give_neighbour_of_element(aElement,aElementList.Dim,tBuffer,aElementList.NumElements);
    Mat<uint> tPossibleNeighbours;
    Mat<uint> tActiveNeighbours(2*aElementList.Dim*pow(2,aElementList.Level),1,UINT_MAX);
    uint tVar = 0, tVara = 0, tVarb = 1;
    uint tParent = 0;
    uint tLevel = 0;
    Mat<uint> tChildren;
    Mat<uint> tNeighbours;
    Mat<uint> tChildNeighbours;
    uint tSwitch = 0;
    if( aElementList.Dim == 2)
    {
        Mat<uint> tNeighbours2D = {{1, 3, 5, 7 }}; // Neighbours in 2D (positions extracted from give_neighbour_of_element)
        tNeighbours = tNeighbours2D;
        Mat<uint> tChildNeighbours2D = {{2, 3}, {1, 3}, {0, 2}, {0, 1}}; // Possible childrens for each neighbour
        tChildNeighbours = tChildNeighbours2D;
    }
    else if( aElementList.Dim == 3)
    {
        Mat<uint> tNeighbours3D = {{4, 10, 12, 14, 16, 22}};  // Neighbours in 2D (positions extracted from give_neighbour_of_element)
        tNeighbours = tNeighbours3D;
        Mat<uint> tChildNeighbours3D = {{4, 5, 6, 7}, {2, 3, 6, 7}, {1, 3, 5, 7}, {0, 2, 4, 6}, {0, 1, 4, 5}, {0, 1, 2, 3}}; // Possible childrens for each neighbour
        tChildNeighbours = tChildNeighbours3D;
    }
    if( aElementList.Dim == 2)
    {
        for(uint i = 0; i < tNeighbours.length(); i++)
        {
            if( (aElementList.ElementActive).test(tNeighbour(tNeighbours(i))) == 1 )
            {
                tActiveNeighbours(tVar) = tNeighbour(tNeighbours(i));
                tVar++;
            }
            else
            {
                tPossibleNeighbours.set_size(2*aElementList.Dim*pow(2,aElementList.Level),1,0);
                tPossibleNeighbours(0) =  tNeighbour(tNeighbours(i));
                tVara = 0;
                tVarb = 1;
                tSwitch = 0;
                while( tPossibleNeighbours(tVara) > 0 )
                {
                    tParent = give_parent_of_element(tPossibleNeighbours(tVara),aElementList.Dim,aElementList.NumElements);
                    if( tParent < UINT_MAX && (aElementList.ElementActive).test(tParent) == 1 )
                    {
                        tActiveNeighbours(tVar) = tParent;
                        tVar++;
                        tSwitch = 1; // If parent is found, no need to check childrens
                    }
                    else if(  tParent < UINT_MAX )
                    {
                        tPossibleNeighbours(tVarb) =  tParent;
                        tVarb++;
                    }
                    tVara++;
                }
                tPossibleNeighbours.set_size(2*aElementList.Dim*pow(2,aElementList.Level),1,0);
                tPossibleNeighbours(0) =  tNeighbour(tNeighbours(i));
                tVara = 0;
                tVarb = 1;
                while( tPossibleNeighbours(tVara) > 0 && tSwitch == 0 )
                {
                    tChildren = give_children_of_element(tPossibleNeighbours(tVara),aElementList.Dim,aElementList.NumElements); // Give the children of the parent element
                    tLevel = give_element_level(tChildren(0),aElementList.Dim,aElementList.NumElements);
                    if( (aElementList.ElementActive).test(tChildren(tChildNeighbours(i,0))) == 1 )
                    {
                        tActiveNeighbours(tVar) = tChildren(tChildNeighbours(i,0));
                        tVar++;
                    }
                    else if( tLevel < aElementList.Level )
                    {
                        tPossibleNeighbours(tVarb) =  tChildren(tChildNeighbours(i,0));
                        tVarb++;
                    }
                    if( (aElementList.ElementActive).test(tChildren(tChildNeighbours(i,1))) == 1 )
                    {
                        tActiveNeighbours(tVar) = tChildren(tChildNeighbours(i,1));
                        tVar++;
                    }
                    else if( tLevel < aElementList.Level )
                    {
                        tPossibleNeighbours(tVarb) =  tChildren(tChildNeighbours(i,1));
                        tVarb++;
                    }
                    tVara++;
                }
            }
        }
    }
    tActiveNeighbours.resize(tVar,1);
    return tActiveNeighbours;
}

moris::Mat<uint>
moris::Hierarchical_Mesh::give_neighbour_of_basis(uint & aBasis,
        uint & aDim,
        uint & aPolynomial,
        uint & aBuffer,
        Mat<uint> & aNumElements)
{
    Mat<uint> tBasis_position=give_position_of_basis(aBasis,aDim,aPolynomial,aNumElements);
    uint tBasis_level=give_basis_level(aBasis,aDim,aPolynomial,aNumElements);
    Mat<uint> Basis_neighbour(1,pow(1+2*aBuffer,aDim),UINT_MAX); // Element plus neighbours (9 for aDim = 2, 27 for aDim = 3)
    Mat<uint> tNeighbour_position(1,aDim); // Position of a neighbour element;
    if(aDim == 2)
    {
        if( aBuffer == 1)
        {
            tNeighbour_position(0) = tBasis_position(0) - 1;   tNeighbour_position(1) = tBasis_position(1) - 1;
            Basis_neighbour(0) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0);       tNeighbour_position(1) = tBasis_position(1) - 1;
            Basis_neighbour(1) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) + 1;   tNeighbour_position(1) = tBasis_position(1) - 1;
            Basis_neighbour(2) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) - 1;   tNeighbour_position(1) = tBasis_position(1);
            Basis_neighbour(3) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0);       tNeighbour_position(1) = tBasis_position(1);
            Basis_neighbour(4) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) + 1;   tNeighbour_position(1) = tBasis_position(1);
            Basis_neighbour(5) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) - 1;   tNeighbour_position(1) = tBasis_position(1) + 1;
            Basis_neighbour(6) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0);       tNeighbour_position(1) = tBasis_position(1) + 1;
            Basis_neighbour(7) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) + 1;   tNeighbour_position(1) = tBasis_position(1) + 1;
            Basis_neighbour(8) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
        }
        else if( aBuffer > 1) // If a higher buffer layer is needed, a loop generates the neighbours
        {
            //            uint tVar = 0;
            //            Mat<uint> aIJKPosition(aDim,1,0);
            //            for(int j = -((int)aBuffer); j <= ((int)aBuffer); j++ )
            //            {
            //                for(int i = -((int)aBuffer); i <= ((int)aBuffer); i++ )
            //                {
            //                    if ( tBasis_position(0)+i >= 1*pow(2,tElement_level) && tBasis_position(0)+i < aNumElements(0)*pow(2,tElement_level)-1*pow(2,tElement_level) && tBasis_position(1)+j >= 1*pow(2,tElement_level) && tBasis_position(1)+j < aNumElements(1)*pow(2,tElement_level)-1*pow(2,tElement_level) )
            //                    {
            //                        tNeighbour_position(0) = tBasis_position(0) + i;    tNeighbour_position(1) = tBasis_position(1) + j;
            //                        Basis_neighbour(tVar)=give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
            //                    }
            //                    tVar++;
            //                }
            //            }
        }
    }
    else if(aDim == 3)
    {
        if( aBuffer == 1)
        {
            tNeighbour_position(0) = tBasis_position(0) - 1;   tNeighbour_position(1) = tBasis_position(1) - 1;   tNeighbour_position(2) = tBasis_position(2) - 1;
            Basis_neighbour(0) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0);       tNeighbour_position(1) = tBasis_position(1) - 1;   tNeighbour_position(2) = tBasis_position(2) - 1;
            Basis_neighbour(1) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) + 1;   tNeighbour_position(1) = tBasis_position(1) - 1;   tNeighbour_position(2) = tBasis_position(2) - 1;
            Basis_neighbour(2) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) - 1;   tNeighbour_position(1) = tBasis_position(1);       tNeighbour_position(2) = tBasis_position(2) - 1;
            Basis_neighbour(3) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0);       tNeighbour_position(1) = tBasis_position(1);       tNeighbour_position(2) = tBasis_position(2) - 1;
            Basis_neighbour(4) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) + 1;   tNeighbour_position(1) = tBasis_position(1);       tNeighbour_position(2) = tBasis_position(2) - 1;
            Basis_neighbour(5) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) - 1;   tNeighbour_position(1) = tBasis_position(1) + 1;   tNeighbour_position(2) = tBasis_position(2) - 1;
            Basis_neighbour(6) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0);       tNeighbour_position(1) = tBasis_position(1) + 1;   tNeighbour_position(2) = tBasis_position(2) - 1;
            Basis_neighbour(7) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) + 1;   tNeighbour_position(1) = tBasis_position(1) + 1;   tNeighbour_position(2) = tBasis_position(2) - 1;
            Basis_neighbour(8) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);

            tNeighbour_position(0) = tBasis_position(0) - 1;   tNeighbour_position(1) = tBasis_position(1) - 1;   tNeighbour_position(2) = tBasis_position(2);
            Basis_neighbour(9) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0);       tNeighbour_position(1) = tBasis_position(1) - 1;   tNeighbour_position(2) = tBasis_position(2);
            Basis_neighbour(10) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) + 1;   tNeighbour_position(1) = tBasis_position(1) - 1;   tNeighbour_position(2) = tBasis_position(2);
            Basis_neighbour(11) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) - 1;   tNeighbour_position(1) = tBasis_position(1);       tNeighbour_position(2) = tBasis_position(2);
            Basis_neighbour(12) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0);       tNeighbour_position(1) = tBasis_position(1);       tNeighbour_position(2) = tBasis_position(2);
            Basis_neighbour(13) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) + 1;   tNeighbour_position(1) = tBasis_position(1);       tNeighbour_position(2) = tBasis_position(2);
            Basis_neighbour(14) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) - 1;   tNeighbour_position(1) = tBasis_position(1) + 1;   tNeighbour_position(2) = tBasis_position(2);
            Basis_neighbour(15) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0);       tNeighbour_position(1) = tBasis_position(1) + 1;   tNeighbour_position(2) = tBasis_position(2);
            Basis_neighbour(16) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) + 1;   tNeighbour_position(1) = tBasis_position(1) + 1;   tNeighbour_position(2) = tBasis_position(2);
            Basis_neighbour(17) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);

            tNeighbour_position(0) = tBasis_position(0) - 1;   tNeighbour_position(1) = tBasis_position(1) - 1;   tNeighbour_position(2) = tBasis_position(2) + 1;
            Basis_neighbour(18) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0);       tNeighbour_position(1) = tBasis_position(1) - 1;   tNeighbour_position(2) = tBasis_position(2) + 1;
            Basis_neighbour(19) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) + 1;   tNeighbour_position(1) = tBasis_position(1) - 1;   tNeighbour_position(2) = tBasis_position(2) + 1;
            Basis_neighbour(20) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) - 1;   tNeighbour_position(1) = tBasis_position(1);       tNeighbour_position(2) = tBasis_position(2) + 1;
            Basis_neighbour(21) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0);       tNeighbour_position(1) = tBasis_position(1);       tNeighbour_position(2) = tBasis_position(2) + 1;
            Basis_neighbour(22) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) + 1;   tNeighbour_position(1) = tBasis_position(1);       tNeighbour_position(2) = tBasis_position(2) + 1;
            Basis_neighbour(23) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) - 1;   tNeighbour_position(1) = tBasis_position(1) + 1;   tNeighbour_position(2) = tBasis_position(2) + 1;
            Basis_neighbour(24) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0);       tNeighbour_position(1) = tBasis_position(1) + 1;   tNeighbour_position(2) = tBasis_position(2) + 1;
            Basis_neighbour(25) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
            tNeighbour_position(0) = tBasis_position(0) + 1;   tNeighbour_position(1) = tBasis_position(1) + 1;   tNeighbour_position(2) = tBasis_position(2) + 1;
            Basis_neighbour(26) = give_basis_of_position(tBasis_level,aDim,aPolynomial,aNumElements,tNeighbour_position);
        }
        else if( aBuffer > 1) // If a higher buffer layer is needed, a loop generates the neighbours
        {
            //            uint tVar = 0;
            //            Mat<uint> aIJKPosition(aDim,1,0);
            //            for(int k = -((int)aBuffer); k <= ((int)aBuffer); k++ )
            //            {
            //                for(int j = -((int)aBuffer); j <= ((int)aBuffer); j++ )
            //                {
            //                    for(int i = -((int)aBuffer); i <= ((int)aBuffer); i++ )
            //                    {
            //                        if ( tBasis_position(0)+i >= 1*pow(2,tElement_level) && tBasis_position(0)+i < aNumElements(0)*pow(2,tElement_level)-1*pow(2,tElement_level) && tBasis_position(1)+j >= 1*pow(2,tElement_level) && tBasis_position(1)+j < aNumElements(1)*pow(2,tElement_level)-1*pow(2,tElement_level) && tBasis_position(2)+k >= 1*pow(2,tElement_level) && tBasis_position(2)+k < aNumElements(2)*pow(2,tElement_level)-1*pow(2,tElement_level) )
            //                        {
            //                            aIJKPosition(0) = tBasis_position(0) + i; aIJKPosition(1) = tBasis_position(1) + j; aIJKPosition(2) = tBasis_position(2) + k;
            //                            Basis_neighbour(tVar)=give_element_of_position(tElement_level,aDim,aNumElements,aIJKPosition);
            //                        }
            //                        tVar++;
            //                    }
            //                }
            //            }
        }
    }
    return Basis_neighbour;
}

moris::Mat<uint>
moris::Hierarchical_Mesh::give_neighbour_of_element_for_floodfill(uint & aElement,
        uint & aDim,
        uint & aBuffer,
        Mat<uint> & aNumElements)
{
    Mat<uint> tNeighbourElements = give_neighbour_of_element(aElement,aDim,aBuffer,aNumElements);
    Mat<uint> tNeighbourElementsMod(2*aDim,1,UINT_MAX);
    if( aDim == 2)
    {
        tNeighbourElementsMod(0) = tNeighbourElements(1);
        tNeighbourElementsMod(1) = tNeighbourElements(3);
        tNeighbourElementsMod(2) = tNeighbourElements(5);
        tNeighbourElementsMod(3) = tNeighbourElements(7);
    }
    else if( aDim == 3)
    {
        tNeighbourElementsMod(0) = tNeighbourElements(10);
        tNeighbourElementsMod(1) = tNeighbourElements(12);
        tNeighbourElementsMod(2) = tNeighbourElements(14);
        tNeighbourElementsMod(3) = tNeighbourElements(16);
        tNeighbourElementsMod(4) = tNeighbourElements(4);
        tNeighbourElementsMod(5) = tNeighbourElements(22);
    }
    return tNeighbourElementsMod;
}

moris::Mat<uint>
moris::Hierarchical_Mesh::give_neighbour_stencil_of_element(uint & aFeatureResolution,
        uint & aElement,
        uint & aDim,
        uint & aPolynomial,
        Mat<uint> & aNumElements,
        BoostBitset & ElementActive)
{
    Mat<uint> tElement_position=give_position_of_element(aElement,aDim,aNumElements);
    uint tElement_level=give_element_level(aElement,aDim,aNumElements);
    Mat<uint> tNeighbour_position(aDim,1,UINT_MAX); // Position of a neighbour element;
    Mat<uint> Element_neighbour;
    uint tVar = 0; // Temporary variable for loop
    uint tElement;
    int i = 0, j = 0, k = 0;
    if(aDim == 2)
    {
        Element_neighbour.set_size(4,2*aFeatureResolution+3,0);
        tVar = 0;
        j = 0;
        for(i = -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++)
        {
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < aNumElements(0)*pow(2,tElement_level)-aPolynomial &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < aNumElements(1)*pow(2,tElement_level)-aPolynomial)
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(0,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(0,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        i = 0;
        for(j = -((int)aFeatureResolution+1); j<=((int)aFeatureResolution+1); j++)
        {
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < aNumElements(0)*pow(2,tElement_level)-aPolynomial &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < aNumElements(1)*pow(2,tElement_level)-aPolynomial)
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(1,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(1,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        for(j = -((int)aFeatureResolution+1); j<=((int)aFeatureResolution+1); j++)
        {
            i = j;
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < aNumElements(0)*pow(2,tElement_level)-aPolynomial &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < aNumElements(1)*pow(2,tElement_level)-aPolynomial)
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(2,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(2,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        for(j = -((int)aFeatureResolution+1); j<=((int)aFeatureResolution+1); j++)
        {
            i = -j;
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < aNumElements(0)*pow(2,tElement_level)-aPolynomial &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < aNumElements(1)*pow(2,tElement_level)-aPolynomial)
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(3,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(3,tVar) = tElement;
                }
            }
            tVar++;
        }
    }
    else if(aDim == 3)
    {
        Element_neighbour.set_size(12,2*aFeatureResolution+3,0);
        tVar = 0;
        k = 0;
        j = 0;
        for(i = -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++)
        {
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < aNumElements(0)*pow(2,tElement_level)-aPolynomial &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < aNumElements(1)*pow(2,tElement_level)-aPolynomial &&
                    tElement_position(2)+k >= aPolynomial*pow(2,tElement_level) && tElement_position(2)+k < aNumElements(2)*pow(2,tElement_level)-aPolynomial)
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;  tNeighbour_position(2) = tElement_position(2) + k;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(0,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(0,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        k = 0;
        i = 0;
        for(j = -((int)aFeatureResolution+1); j<=((int)aFeatureResolution+1); j++)
        {
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < (aNumElements(0)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < (aNumElements(1)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(2)+k >= aPolynomial*pow(2,tElement_level) && tElement_position(2)+k < (aNumElements(2)-aPolynomial)*pow(2,tElement_level))
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;  tNeighbour_position(2) = tElement_position(2) + k;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(1,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(1,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        j = 0;
        i = 0;
        for(k = -((int)aFeatureResolution+1); k<=((int)aFeatureResolution+1); k++)
        {
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < (aNumElements(0)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < (aNumElements(1)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(2)+k >= aPolynomial*pow(2,tElement_level) && tElement_position(2)+k < (aNumElements(2)-aPolynomial)*pow(2,tElement_level))
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;  tNeighbour_position(2) = tElement_position(2) + k;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(2,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(2,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        k = 0;
        for(i = -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++)
        {
            j = i;
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < (aNumElements(0)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < (aNumElements(1)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(2)+k >= aPolynomial*pow(2,tElement_level) && tElement_position(2)+k < (aNumElements(2)-aPolynomial)*pow(2,tElement_level))
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;  tNeighbour_position(2) = tElement_position(2) + k;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(3,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(3,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        k = 0;
        for(i = -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++)
        {
            j = -i;
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < (aNumElements(0)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < (aNumElements(1)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(2)+k >= aPolynomial*pow(2,tElement_level) && tElement_position(2)+k < (aNumElements(2)-aPolynomial)*pow(2,tElement_level))
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;  tNeighbour_position(2) = tElement_position(2) + k;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(4,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(4,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        i = 0;
        for(j = -((int)aFeatureResolution+1); j<=((int)aFeatureResolution+1); j++)
        {
            k = j;
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < (aNumElements(0)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < (aNumElements(1)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(2)+k >= aPolynomial*pow(2,tElement_level) && tElement_position(2)+k < (aNumElements(2)-aPolynomial)*pow(2,tElement_level))
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;  tNeighbour_position(2) = tElement_position(2) + k;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(5,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(5,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        i = 0;
        for(j = -((int)aFeatureResolution+1); j<=((int)aFeatureResolution+1); j++)
        {
            k = -j;
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < (aNumElements(0)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < (aNumElements(1)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(2)+k >= aPolynomial*pow(2,tElement_level) && tElement_position(2)+k < (aNumElements(2)-aPolynomial)*pow(2,tElement_level))
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;  tNeighbour_position(2) = tElement_position(2) + k;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(6,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(6,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        j = 0;
        for(i= -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++)
        {
            k = i;
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < (aNumElements(0)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < (aNumElements(1)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(2)+k >= aPolynomial*pow(2,tElement_level) && tElement_position(2)+k < (aNumElements(2)-aPolynomial)*pow(2,tElement_level))
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;  tNeighbour_position(2) = tElement_position(2) + k;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(7,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(7,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        j = 0;
        for(i= -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++)
        {
            k = -i;
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < (aNumElements(0)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < (aNumElements(1)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(2)+k >= aPolynomial*pow(2,tElement_level) && tElement_position(2)+k < (aNumElements(2)-aPolynomial)*pow(2,tElement_level))
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;  tNeighbour_position(2) = tElement_position(2) + k;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(8,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(8,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        for(i= -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++)
        {
            k = i;
            j = i;
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < (aNumElements(0)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < (aNumElements(1)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(2)+k >= aPolynomial*pow(2,tElement_level) && tElement_position(2)+k < (aNumElements(2)-aPolynomial)*pow(2,tElement_level))
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;  tNeighbour_position(2) = tElement_position(2) + k;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(9,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(9,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        for(i= -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++)
        {
            k = i;
            j = -i;
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < (aNumElements(0)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < (aNumElements(1)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(2)+k >= aPolynomial*pow(2,tElement_level) && tElement_position(2)+k < (aNumElements(2)-aPolynomial)*pow(2,tElement_level))
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;  tNeighbour_position(2) = tElement_position(2) + k;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(10,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(10,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        for(i= -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++)
        {
            k = -i;
            j = i;
            if( tElement_position(0)+i >= aPolynomial*pow(2,tElement_level) && tElement_position(0)+i < (aNumElements(0)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(1)+j >= aPolynomial*pow(2,tElement_level) && tElement_position(1)+j < (aNumElements(1)-aPolynomial)*pow(2,tElement_level) &&
                    tElement_position(2)+k >= aPolynomial*pow(2,tElement_level) && tElement_position(2)+k < (aNumElements(2)-aPolynomial)*pow(2,tElement_level))
            {
                tNeighbour_position(0) = tElement_position(0) + i;   tNeighbour_position(1) = tElement_position(1) + j;  tNeighbour_position(2) = tElement_position(2) + k;
                tElement = give_element_of_position(tElement_level,aDim,aNumElements,tNeighbour_position);
                Element_neighbour(11,tVar) = 0;
                if( ElementActive.test(tElement) == 1)
                {
                    Element_neighbour(11,tVar) = tElement;
                }
            }
            tVar++;
        }
    }
    return Element_neighbour;
}

uint
moris::Hierarchical_Mesh::give_basis_level(uint & aBasis,
        uint & aDim,
        uint & aPolynomial,
        Mat<uint> & aNumElements)
{
    uint Basis_level=1; // output variable
    uint tLevel=1; // temporary variable for the while loop
    if( aDim == 2 )
    {
        uint tNumber_of_basis=(pow(2,Basis_level-1)*aNumElements(0)+aPolynomial)*(pow(2,Basis_level-1)*aNumElements(1)+aPolynomial); // temporary variable for the number of elements per level
        while(tLevel>0)
        {
            tLevel=floor((real)(aBasis+1)/tNumber_of_basis);
            if( (real)(aBasis+1)/tNumber_of_basis == 1 || tLevel<1)
            {
                tLevel = 0;
            }
            else
            {
                Basis_level += 1;
                tNumber_of_basis += (pow(2,Basis_level-1)*aNumElements(0)+aPolynomial)*(pow(2,Basis_level-1)*aNumElements(1)+aPolynomial);
            }
        }
        Basis_level -= 1;
    }
    else if(aDim == 3)
    {
        uint tNumber_of_basis=(pow(2,Basis_level-1)*aNumElements(0)+aPolynomial)*(pow(2,Basis_level-1)*aNumElements(1)+aPolynomial)*(pow(2,Basis_level-1)*aNumElements(2)+aPolynomial); // temporary variable for the number of elements per level
        while(tLevel>0)
        {
            tLevel=floor((real)(aBasis+1)/tNumber_of_basis);
            if( (real)(aBasis+1)/tNumber_of_basis == 1 || tLevel<1)
            {
                tLevel = 0;
            }
            else
            {
                Basis_level += 1;
                tNumber_of_basis += (pow(2,Basis_level-1)*aNumElements(0)+aPolynomial)*(pow(2,Basis_level-1)*aNumElements(1)+aPolynomial)*(pow(2,Basis_level-1)*aNumElements(2)+aPolynomial);
            }
        }
        Basis_level -= 1;
    }
    return Basis_level;
}

uint
moris::Hierarchical_Mesh::give_number_of_basis(uint & aLevel,
        uint & aDim,
        uint & aPolynomial,
        Mat<uint> & aNumElements)
{
    uint Basis_number=0; // output variable
    if( aDim == 2 )
    {
        for(uint i = 0; i<aLevel+1; i++)
        {
            Basis_number += (pow(2,i)*aNumElements(0)+aPolynomial)*(pow(2,i)*aNumElements(1)+aPolynomial);
        }
    }
    else if(aDim == 3)
    {
        for(uint i = 0; i<aLevel+1; i++)
        {
            Basis_number += (pow(2,i)*aNumElements(0)+aPolynomial)*(pow(2,i)*aNumElements(1)+aPolynomial)*(pow(2,i)*aNumElements(2)+aPolynomial);
        }
    }
    return Basis_number;
}

moris::Mat<uint>
moris::Hierarchical_Mesh::give_position_of_basis(uint & aBasis,
        uint & aDim,
        uint & aPolynomial,
        Mat<uint> & aNumElements)
{
    Mat<uint> Basis_position(1,aDim); // Element positions in the direction x,y,z
    uint tBasis_level=give_basis_level(aBasis,aDim,aPolynomial,aNumElements);
    uint tNumber_of_basis=0; // temporary variable for the number of elements per level
    if( aDim == 2 )
    {
        if(tBasis_level !=0)
        {
            tBasis_level -= 1; //Count basis functions until tBasis_level-1
            tNumber_of_basis = give_number_of_basis(tBasis_level,aDim,aPolynomial,aNumElements);
            tBasis_level += 1; //+1 to give the original level of basis aBasis
        }
        Basis_position(1) = ceil((aBasis+1-tNumber_of_basis)/((pow(2,tBasis_level)*aNumElements(0)+aPolynomial)));
        Basis_position(0) = aBasis+1 - tNumber_of_basis + (pow(2,tBasis_level)*aNumElements(0)+aPolynomial) - Basis_position(1)*(pow(2,tBasis_level)*aNumElements(0)+aPolynomial); // Calculate the x-direction with the calculated y-position
        Basis_position(1) -= 1; // -1 to get an indexed basis position, starting with zero
        Basis_position(0) -= 1; // -1 to get an indexed basis position, starting with zero
    }
    else if(aDim == 3)
    {
        if(tBasis_level !=0)
        {
            tBasis_level -= 1; //Count basis functions until tBasis_level-1
            tNumber_of_basis = give_number_of_basis(tBasis_level,aDim,aPolynomial,aNumElements);
            tBasis_level += 1; //+1 to give the original level of basis aBasis
        }
        Basis_position(2) = ceil((aBasis+1-tNumber_of_basis)/((pow(2,tBasis_level)*aNumElements(0)+aPolynomial)*(pow(2,tBasis_level)*aNumElements(1)+aPolynomial))); // Calculate the position in z-direction
        Basis_position(1) = ceil((aBasis+1-tNumber_of_basis-(Basis_position(2)-1)*(pow(2,tBasis_level)*aNumElements(0)+aPolynomial)*(pow(2,tBasis_level)*aNumElements(1)+aPolynomial))/((pow(2,tBasis_level)*aNumElements(0)+aPolynomial))); // Calculate the y-direction with the calculated z-position
        Basis_position(0) = aBasis+1 - tNumber_of_basis + (pow(2,tBasis_level)*aNumElements(0)+aPolynomial) - Basis_position(1)*(pow(2,tBasis_level)*aNumElements(0)+aPolynomial) + (pow(2,tBasis_level)*aNumElements(0)+aPolynomial)*(pow(2,tBasis_level)*aNumElements(1)+aPolynomial) - Basis_position(2)*(pow(2,tBasis_level)*aNumElements(0)+aPolynomial)*(pow(2,tBasis_level)*aNumElements(1)+aPolynomial); // Calculate the x-direction with the calculated y and z-position
        Basis_position(2) -= 1; // -1 to get an indexed basis position, starting with zero
        Basis_position(1) -= 1; // -1 to get an indexed basis position, starting with zero
        Basis_position(0) -= 1; // -1 to get an indexed basis position, starting with zero
    }
    return Basis_position;
}

uint
moris::Hierarchical_Mesh::give_basis_of_position(uint & aLevel,
        uint & aDim,
        uint & aPolynomial,
        Mat<uint> & aNumElements,
        Mat<uint> & aIJKPosition)
{
    uint tBasis = 0;
    uint tLevel = aLevel - 1;
    uint tY_Position = 0; // Position of the Element in y-direction
    if( aDim == 2 )
    {
        tY_Position = aIJKPosition(1)*(aNumElements(0)+aPolynomial);
        if( aLevel>0 )
        {
            tBasis = give_number_of_basis(tLevel,aDim,aPolynomial,aNumElements); //Give number of basis functions for the last level to have the initial number at position 0,0
            tY_Position = aIJKPosition(1)*(aPolynomial+aNumElements(0)*pow(2,aLevel));  // Position of the Element in y-direction for a specific level
        }
        tBasis = tBasis + aIJKPosition(0) + tY_Position; //Basis function with the calculated y-position
    }
    else if(aDim == 3)
    {
        tY_Position = aIJKPosition(1)*(aPolynomial+aNumElements(0))+aIJKPosition(2)*(aPolynomial+aNumElements(0))*(aPolynomial+aNumElements(1));
        if( aLevel>0 )
        {
            tBasis = give_number_of_basis(tLevel,aDim,aPolynomial,aNumElements); //Give number of basis functions for the last level to have the initial number at position 0,0,0
            tY_Position = aIJKPosition(1)*(aPolynomial+aNumElements(0)*pow(2,aLevel))+aIJKPosition(2)*(aPolynomial+aNumElements(0)*pow(2,aLevel))*(aPolynomial+aNumElements(1)*pow(2,aLevel));  // Position of the Element in y-direction for a specific level
        }
        tBasis = tBasis + aIJKPosition(0) + tY_Position; //Basis function with the calculated y-position
    }
    return tBasis;
}

uint
moris::Hierarchical_Mesh::give_basis_of_parent(uint & aBasis,
        uint & aDim,
        uint & aPolynomial,
        Mat<uint> & aNumElements)
{
    Mat<uint> tIJKPosition = give_position_of_basis(aBasis,aDim,aPolynomial,aNumElements);
    uint tBasis = UINT_MAX;
    uint tLevel = give_basis_level(aBasis,aDim,aPolynomial,aNumElements);
    if( tLevel > 0)
    {
        tLevel -= 1; // Parent basis function is on a lower level!
        if( aDim == 2)
        {
            if( tIJKPosition(0) %2 == 0 && tIJKPosition(1) %2 == 0) // Check if a Position has a even number, otherwise there is no parent existent
            {
                tIJKPosition(0) /=2; tIJKPosition(1) /=2;
                tBasis = give_basis_of_position(tLevel,aDim,aPolynomial,aNumElements,tIJKPosition);
            }
        }
        else if( aDim == 3)
        {
            if( tIJKPosition(0) %2 == 0 && tIJKPosition(1) %2 == 0 && tIJKPosition(2) %2== 0) // Check if a Position has a even number, otherwise there is no parent existent
            {
                tIJKPosition(0) /=2; tIJKPosition(1) /=2; tIJKPosition(2) /=2; // Divide by two, to get the position on the coarser level
                tBasis = give_basis_of_position(tLevel,aDim,aPolynomial,aNumElements,tIJKPosition);
            }
        }
    }
    return tBasis;
}

moris::Mat<uint>
moris::Hierarchical_Mesh::give_basis_of_element(uint & aElement,
        uint & aDim,
        uint & aPolynomial,
        Mat<uint> & aNumElements)
{
    Mat<uint> tBasis(pow(aPolynomial+1,aDim),1); // Basis functions of an element
    uint tElement_level=give_element_level(aElement,aDim,aNumElements);
    Mat<uint> tElement_position=give_position_of_element(aElement,aDim,aNumElements);
    uint tNumber_of_basis = 0;
    uint tVar = 0; //Temporary variable for the for loop
    if(tElement_level !=0)
    {
        tElement_level -= 1; //Count basis functions until tBasis_level-1
        tNumber_of_basis = give_number_of_basis(tElement_level,aDim,aPolynomial,aNumElements);
        tElement_level += 1; //+1 to give the original level of basis aBasis
    }
    if( aDim == 2 )
    {
        if(aPolynomial == 1)
        {
            tBasis(0) = tElement_position(0) + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(1) = tElement_position(0) + 1 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(2) = tElement_position(0) + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(3) = tElement_position(0) + 1 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
        }
        else if(aPolynomial == 2)
        {
            tBasis(0) = tElement_position(0) +     tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(1) = tElement_position(0) + 1 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(2) = tElement_position(0) + 2 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(3) = tElement_position(0) +     (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(4) = tElement_position(0) + 1 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(5) = tElement_position(0) + 2 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(6) = tElement_position(0) +     (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(7) = tElement_position(0) + 1 + (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(8) = tElement_position(0) + 2 + (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
        }
        else if(aPolynomial == 3)
        {
            tBasis(0) = tElement_position(0) +     tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(1) = tElement_position(0) + 1 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(2) = tElement_position(0) + 2 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(3) = tElement_position(0) + 3 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(4) = tElement_position(0) +     (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(5) = tElement_position(0) + 1 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(6) = tElement_position(0) + 2 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(7) = tElement_position(0) + 3 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(8) = tElement_position(0) +     (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(9) =  tElement_position(0) + 1 + (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(10) = tElement_position(0) + 2 + (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(11) = tElement_position(0) + 3 + (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(12) = tElement_position(0) +     (tElement_position(1)+3)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(13) = tElement_position(0) + 1 + (tElement_position(1)+3)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(14) = tElement_position(0) + 2 + (tElement_position(1)+3)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
            tBasis(15) = tElement_position(0) + 3 + (tElement_position(1)+3)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
        }
        else // The basis functions will be calculated, if higher order polynomials are needed
        {
            for(uint i = 0; i<aPolynomial+1; i++)
            {
                for(uint j = 0; j<aPolynomial+1; j++)
                {
                    tBasis(tVar) = (tElement_position(0)+j) + (tElement_position(1)+i)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)+tNumber_of_basis;
                    tVar++;
                }
            }
        }
    }
    else if(aDim == 3)
    {
        if(aPolynomial == 1)
        {
            tBasis(0) = tElement_position(0) +     tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) +     tElement_position(2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(1) = tElement_position(0) + 1 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) +     tElement_position(2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(2) = tElement_position(0) +     (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + tElement_position(2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(3) = tElement_position(0) + 1 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + tElement_position(2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(4) = tElement_position(0) +     tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) +     (tElement_position(2)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(5) = tElement_position(0) + 1 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) +     (tElement_position(2)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(6) = tElement_position(0) +     (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(7) = tElement_position(0) + 1 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
        }
        else if(aPolynomial == 2)
        {
            tBasis(0) = tElement_position(0) +     tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) +     tElement_position(2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(1) = tElement_position(0) + 1 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) +     tElement_position(2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(2) = tElement_position(0) + 2 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) +     tElement_position(2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(3) = tElement_position(0) +     (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + tElement_position(2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(4) = tElement_position(0) + 1 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + tElement_position(2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(5) = tElement_position(0) + 2 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + tElement_position(2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(6) = tElement_position(0) +     (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + tElement_position(2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(7) = tElement_position(0) + 1 + (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + tElement_position(2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(8) = tElement_position(0) + 2 + (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + tElement_position(2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(9) =  tElement_position(0) +     tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) +     (tElement_position(2)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(10) = tElement_position(0) + 1 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) +     (tElement_position(2)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(11) = tElement_position(0) + 2 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) +     (tElement_position(2)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(12) = tElement_position(0) +     (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(13) = tElement_position(0) + 1 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(14) = tElement_position(0) + 2 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(15) = tElement_position(0) +     (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(16) = tElement_position(0) + 1 + (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(17) = tElement_position(0) + 2 + (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(18) = tElement_position(0) +     tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) +     (tElement_position(2)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(19) = tElement_position(0) + 1 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) +     (tElement_position(2)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(20) = tElement_position(0) + 2 + tElement_position(1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) +     (tElement_position(2)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(21) = tElement_position(0) +     (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(22) = tElement_position(0) + 1 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(23) = tElement_position(0) + 2 + (tElement_position(1)+1)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(24) = tElement_position(0) +     (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(25) = tElement_position(0) + 1 + (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
            tBasis(26) = tElement_position(0) + 2 + (tElement_position(1)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+2)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial) + tNumber_of_basis;
        }
        else // The basis functions will be calculated, if higher order polynomials are needed
        {
            for(uint k = 0; k<aPolynomial+1; k++)
            {
                for(uint i = 0; i<aPolynomial+1; i++)
                {
                    for(uint j = 0; j<aPolynomial+1; j++)
                    {
                        tBasis(tVar) = (tElement_position(0)+j) + (tElement_position(1)+i)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial) + (tElement_position(2)+k)*(pow(2,tElement_level)*aNumElements(0)+aPolynomial)*(pow(2,tElement_level)*aNumElements(1)+aPolynomial)+tNumber_of_basis;
                        tVar++;
                    }
                }
            }
        }
    }
    return tBasis;
}

moris::Mat<uint>
moris::Hierarchical_Mesh::give_element_of_basis(uint & aBasis,
        uint & aDim,
        uint & aPolynomial,
        Mat<uint> & aNumElements)
{
    uint tBasis_level=give_basis_level(aBasis,aDim,aPolynomial,aNumElements);
    Mat<uint> tBasis_position=give_position_of_basis(aBasis,aDim,aPolynomial,aNumElements);
    Mat<uint> tElements(pow((aPolynomial+1),aDim),1,UINT_MAX); // Elements, which have support with the basis function aBasis
    Mat<uint> tElement_position(aDim,1); // Position of the element with respect to the basis function
    if( aDim == 2 )
    {// Elements are hardcoded for aPolynomial = 1 & 2. aPolynomial > 2 is generated by a loop
        if(aPolynomial == 1)
        {
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);   tElement_position(1) = tBasis_position(1);
                (tElements(3))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1);
                (tElements(2))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if((int)(tBasis_position(0))>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);   tElement_position(1) = tBasis_position(1)-1;
                (tElements(1))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if((int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1)-1;
                (tElements(0))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
        }
        if(aPolynomial == 2)
        {
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);   tElement_position(1) = tBasis_position(1);
                (tElements(8))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1);
                (tElements(7))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-2)>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(0)-2) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-2;   tElement_position(1) = tBasis_position(1);
                (tElements(6))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);   tElement_position(1) = tBasis_position(1)-1;
                (tElements(5))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0) - 1;   tElement_position(1) = tBasis_position(1)-1;
                (tElements(4))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-2)>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(0)-2) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0) - 2;   tElement_position(1) = tBasis_position(1)-1;
                (tElements(3))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1)-2)>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1)-2) < (int)aNumElements(1)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);   tElement_position(1) = tBasis_position(1)-2;
                (tElements(2))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1)-2)>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1)-2) < (int)aNumElements(1)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1)-2;
                (tElements(1))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-2)>=0 && (int)(tBasis_position(1)-2)>=0 && (int)(tBasis_position(0)-2) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1)-2) < (int)aNumElements(1)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-2;   tElement_position(1) = tBasis_position(1)-2;
                (tElements(0))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
        }
        else if( aPolynomial > 2) // The elements will be calculated, if higher order polynomials are needed
        {
            uint tVar=0; // Temporary variable for the for loops
            for(int j=-aPolynomial; j<1; j++)
            {
                for(int i=-aPolynomial; i<1; i++)
                {
                    if((int)(tBasis_position(0) + i)>=0 && (int)(tBasis_position(1) + j)>=0 && (int)(tBasis_position(0) + i) < (int)aNumElements(0)*pow(2,tBasis_level) && (int)(tBasis_position(1) + j) < (int)aNumElements(1)*pow(2,tBasis_level) )
                    {
                        tElement_position(0) = tBasis_position(0) + i;   tElement_position(1) = tBasis_position(1) + j;
                        (tElements(tVar))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
                        tVar++;
                    }
                }
            }
        }
    }
    else if(aDim == 3)
    {// Elements are hardcoded for aPolynomial = 1 & 2. aPolynomial > 2 is generated by a loop
        if(aPolynomial == 1)
        {
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(2))>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);     tElement_position(1) = tBasis_position(1);     tElement_position(2) = tBasis_position(2);
                (tElements(7))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(2))>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1);     tElement_position(2) = tBasis_position(2);
                (tElements(6))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(2))>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);     tElement_position(1) = tBasis_position(1)-1;   tElement_position(2) = tBasis_position(2);
                (tElements(5))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(2))>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1)-1;   tElement_position(2) = tBasis_position(2);
                (tElements(4))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(2)-1)>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-1) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);     tElement_position(1) = tBasis_position(1);     tElement_position(2) = tBasis_position(2)-1;
                (tElements(3))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(2)-1)>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-1) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1);     tElement_position(2) = tBasis_position(2)-1;
                (tElements(2))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(2)-1)>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-1) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);     tElement_position(1) = tBasis_position(1)-1;   tElement_position(2) = tBasis_position(2)-1;
                (tElements(1))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(2)-1)>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-1) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1)-1;   tElement_position(2) = tBasis_position(2)-1;
                (tElements(0))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
        }
        if(aPolynomial == 2)
        {
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(2))>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);     tElement_position(1) = tBasis_position(1);     tElement_position(2) = tBasis_position(2);
                (tElements(26))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(2))>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1);     tElement_position(2) = tBasis_position(2);
                (tElements(25))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-2)>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(2))>=0 && (int)(tBasis_position(0)-2) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-2;   tElement_position(1) = tBasis_position(1);     tElement_position(2) = tBasis_position(2);
                (tElements(24))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(2))>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);     tElement_position(1) = tBasis_position(1)-1;   tElement_position(2) = tBasis_position(2);
                (tElements(23))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(2))>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1)-1;   tElement_position(2) = tBasis_position(2);
                (tElements(22))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-2)>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(2))>=0 && (int)(tBasis_position(0)-2) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-2;   tElement_position(1) = tBasis_position(1)-1;   tElement_position(2) = tBasis_position(2);
                (tElements(21))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1)-2)>=0 && (int)(tBasis_position(2))>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-2) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);     tElement_position(1) = tBasis_position(1)-2;   tElement_position(2) = tBasis_position(2);
                (tElements(20))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1)-2)>=0 && (int)(tBasis_position(2))>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-2) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1)-2;   tElement_position(2) = tBasis_position(2);
                (tElements(19))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-2)>=0 && (int)(tBasis_position(1)-2)>=0 && (int)(tBasis_position(2))>=0 && (int)(tBasis_position(0)-2) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-2) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-2;   tElement_position(1) = tBasis_position(1)-2;   tElement_position(2) = tBasis_position(2);
                (tElements(18))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(2)-1)>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-1) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);     tElement_position(1) = tBasis_position(1);     tElement_position(2) = tBasis_position(2) - 1;
                (tElements(17))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(2)-1)>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-1) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1);     tElement_position(2) = tBasis_position(2) - 1;
                (tElements(16))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-2)>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(2)-1)>=0 && (int)(tBasis_position(0)-2) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-1) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-2;   tElement_position(1) = tBasis_position(1);     tElement_position(2) = tBasis_position(2) - 1;
                (tElements(15))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(2)-1)>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-1) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);     tElement_position(1) = tBasis_position(1)-1;   tElement_position(2) = tBasis_position(2) - 1;
                (tElements(14))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(2)-1)>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-1) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1)-1;   tElement_position(2) = tBasis_position(2) - 1;
                (tElements(13))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-2)>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(2)-1)>=0 && (int)(tBasis_position(0)-2) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-1) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-2;   tElement_position(1) = tBasis_position(1)-1;   tElement_position(2) = tBasis_position(2) - 1;
                (tElements(12))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1)-2)>=0 && (int)(tBasis_position(2)-1)>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-2) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-1) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);     tElement_position(1) = tBasis_position(1)-2;   tElement_position(2) = tBasis_position(2) - 1;
                (tElements(11))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1)-2)>=0 && (int)(tBasis_position(2)-1)>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-2) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-1) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1)-2;   tElement_position(2) = tBasis_position(2) - 1;
                (tElements(10))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-2)>=0 && (int)(tBasis_position(1)-2)>=0 && (int)(tBasis_position(2)-1)>=0 && (int)(tBasis_position(0)-2) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-2) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-1) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-2;   tElement_position(1) = tBasis_position(1)-2;   tElement_position(2) = tBasis_position(2) - 1;
                (tElements(9))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(2)-2)>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-2) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);     tElement_position(1) = tBasis_position(1);     tElement_position(2) = tBasis_position(2) - 2;
                (tElements(8))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(2)-2)>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-2) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1);     tElement_position(2) = tBasis_position(2) - 2;
                (tElements(7))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-2)>=0 && (int)(tBasis_position(1))>=0 && (int)(tBasis_position(2)-2)>=0 && (int)(tBasis_position(0)-2) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-2) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-2;   tElement_position(1) = tBasis_position(1);     tElement_position(2) = tBasis_position(2) - 2;
                (tElements(6))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(2)-2)>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-2) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);     tElement_position(1) = tBasis_position(1)-1;   tElement_position(2) = tBasis_position(2) - 2;
                (tElements(5))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(2)-2)>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-2) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1)-1;   tElement_position(2) = tBasis_position(2) - 2;
                (tElements(4))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-2)>=0 && (int)(tBasis_position(1)-1)>=0 && (int)(tBasis_position(2)-2)>=0 && (int)(tBasis_position(0)-2) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-1) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-2) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-2;   tElement_position(1) = tBasis_position(1)-1;   tElement_position(2) = tBasis_position(2) - 2;
                (tElements(3))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0))>=0 && (int)(tBasis_position(1)-2)>=0 && (int)(tBasis_position(2)-2)>=0 && (int)(tBasis_position(0)) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-2) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-2) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0);     tElement_position(1) = tBasis_position(1)-2;   tElement_position(2) = tBasis_position(2) - 2;
                (tElements(2))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-1)>=0 && (int)(tBasis_position(1)-2)>=0 && (int)(tBasis_position(2)-2)>=0 && (int)(tBasis_position(0)-1) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-2) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-2) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-1;   tElement_position(1) = tBasis_position(1)-2;   tElement_position(2) = tBasis_position(2) - 2;
                (tElements(1))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
            if( (int)(tBasis_position(0)-2)>=0 && (int)(tBasis_position(1)-2)>=0 && (int)(tBasis_position(2)-2)>=0 && (int)(tBasis_position(0)-2) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1)-2) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2)-2) < (int)aNumElements(2)*pow(2,tBasis_level) )
            {
                tElement_position(0) = tBasis_position(0)-2;   tElement_position(1) = tBasis_position(1)-2;   tElement_position(2) = tBasis_position(2) - 2;
                (tElements(0))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
            }
        }
        else if( aPolynomial > 2) // The elements will be calculated, if higher order polynomials are needed
        {
            uint tVar=0; // Temporary variable for the for loops
            for(int k=-aPolynomial; k<1; k++)
            {
                for(int j=-aPolynomial; j<1; j++)
                {
                    for(int i=-aPolynomial; i<1; i++)
                    {
                        if((int)(tBasis_position(0) + i)>=0 && (int)(tBasis_position(1) + j)>=0 && (int)(tBasis_position(2) + k)>=0 && (int)(tBasis_position(0) + i) < (int)aNumElements(0)*pow(2,tBasis_level)  && (int)(tBasis_position(1) + j) < (int)aNumElements(1)*pow(2,tBasis_level)  && (int)(tBasis_position(2) + k) < (int)aNumElements(2)*pow(2,tBasis_level) )
                        {
                            tElement_position(0) = tBasis_position(0) + i;   tElement_position(1) = tBasis_position(1) + j;   tElement_position(2) = tBasis_position(2) + k;
                            (tElements(tVar))=give_element_of_position(tBasis_level,aDim,aNumElements,tElement_position);
                            tVar++;
                        }
                    }
                }
            }
        }
    }
    return tElements;
}

void moris::Hierarchical_Mesh::give_active_elements(
        ElementList_struc & aElementList)
{
    uint tMaxNumElemLevel = 0;
    for(uint level = 0; level<aElementList.Level+1; level ++)
    {
        tMaxNumElemLevel += pow(pow(2,aElementList.Dim),level);
    }
    Mat<uint> tPossibleElementListOnProc((aElementList.ElementListOnProcInit).length()*tMaxNumElemLevel,1,UINT_MAX);
    tMaxNumElemLevel -= pow(pow(2,aElementList.Dim),aElementList.Level);
    tPossibleElementListOnProc.rows(0,(aElementList.ElementListOnProcInit).length()-1) = (aElementList.ElementListOnProcInit).rows(0,(aElementList.ElementListOnProcInit).length()-1);
    uint tVar = (aElementList.ElementListOnProcInit).length();
    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    Mat<uint> tChildrenDummy(pow(2,aElementList.Dim),1);
    for(uint i = 0; i<(aElementList.ElementListOnProcInit).length()*tMaxNumElemLevel; i++)
    {
        tChildren = give_children_of_element(tPossibleElementListOnProc(i),aElementList.Dim,aElementList.NumElements);
        tPossibleElementListOnProc.rows(tVar+i*pow(2,aElementList.Dim),tVar+(i+1)*pow(2,aElementList.Dim)-1) = tChildren.rows(0,tChildren.length()-1);
    }
    tPossibleElementListOnProc = unique(tPossibleElementListOnProc);
    if( (aElementList.ElementActive).size() < (tPossibleElementListOnProc.max()+1)  )
        (aElementList.ElementActive).resize(tPossibleElementListOnProc.max()+1);
    (aElementList.ElementListOnProc).set_size(tPossibleElementListOnProc.length(),1);
    tVar = 0;
    for(uint i = 0; i<tPossibleElementListOnProc.length(); i++)
    {
        if( (aElementList.ElementActive).test(tPossibleElementListOnProc(i)) == 1)
        {
            (aElementList.ElementListOnProc)(tVar) = tPossibleElementListOnProc(i);
            tVar++;
        }
    }
    (aElementList.ElementListOnProc).resize(tVar,1);
    (aElementList.ElementListOnProc) = sort((aElementList.ElementListOnProc));
}

void moris::Hierarchical_Mesh::give_active_elements_with_aura(
        ElementList_struc & aElementList)
{
    uint tMaxNumElemLevel = 0;
    for(uint level = 0; level<aElementList.Level+1; level ++)
    {
        tMaxNumElemLevel += pow(pow(2,aElementList.Dim),level);
    }
    Mat<uint> tPossibleElementListOnProc((aElementList.ElementListOnProcAura).length()*tMaxNumElemLevel,1,UINT_MAX);
    tMaxNumElemLevel -= pow(pow(2,aElementList.Dim),aElementList.Level);
    tPossibleElementListOnProc.rows(0,(aElementList.ElementListOnProcAura).length()-1) = (aElementList.ElementListOnProcAura).rows(0,(aElementList.ElementListOnProcAura).length()-1);
    uint tVar = (aElementList.ElementListOnProcAura).length();
    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    for(uint i = 0; i<(aElementList.ElementListOnProcAura).length()*tMaxNumElemLevel; i++)
    {
        tChildren = give_children_of_element(tPossibleElementListOnProc(i),aElementList.Dim,aElementList.NumElements);
        tPossibleElementListOnProc.rows(tVar+i*pow(2,aElementList.Dim),tVar+(i+1)*pow(2,aElementList.Dim)-1) = tChildren.rows(0,tChildren.length()-1);
    }
    tPossibleElementListOnProc = unique(tPossibleElementListOnProc);
    if( (aElementList.ElementActive).size() < (tPossibleElementListOnProc.max()+1) )
        (aElementList.ElementActive).resize(tPossibleElementListOnProc.max()+1);
    (aElementList.ElementListOnProc).set_size(tPossibleElementListOnProc.length(),1);
    tVar = 0;
    for(uint i = 0; i<tPossibleElementListOnProc.length(); i++)
    {
        if( (aElementList.ElementActive).test(tPossibleElementListOnProc(i)) == 1)
        {
            (aElementList.ElementListOnProc)(tVar) = tPossibleElementListOnProc(i);
            tVar++;
        }
    }
    (aElementList.ElementListOnProc).resize(tVar,1);
    (aElementList.ElementListOnProc) = sort((aElementList.ElementListOnProc));
}

Mat<real>
moris::Hierarchical_Mesh::give_coordinate_from_basis(uint & aBasis,
        uint & aPolynomial,
        ElementList_struc & aElementList)
{
    Mat<real> tCoordinates(1,aElementList.Dim,0);
    uint tBasis_level = give_basis_level(aBasis,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
    // Get the Position of the basis function
    Mat<uint> tBasis_position = give_position_of_basis(aBasis,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
    //Calculate the Coordinates of the basis function with respect to the dimensions and the offset
    tCoordinates(0) = (tBasis_position(0)*(aElementList.Dimensions)(0))/((aElementList.NumElements)(0)*pow(2,tBasis_level)) + (aElementList.Dimensions_Offset)(0) - (aPolynomial-1)*(aElementList.Dimensions)(0)/((aElementList.NumElements)(0)*pow(2,tBasis_level)*2); //X-coordinates + offset (if needed);
    tCoordinates(1) = (tBasis_position(1)*(aElementList.Dimensions)(1))/((aElementList.NumElements)(1)*pow(2,tBasis_level)) + (aElementList.Dimensions_Offset)(1) - (aPolynomial-1)*(aElementList.Dimensions)(1)/((aElementList.NumElements)(1)*pow(2,tBasis_level)*2); //Y-Coordinates + offset (if needed);
    if (aElementList.Dim == 3)
        tCoordinates(2) = (tBasis_position(2)*(aElementList.Dimensions)(2))/((aElementList.NumElements)(2)*pow(2,tBasis_level)) + (aElementList.Dimensions_Offset)(2) - (aPolynomial-1)*(aElementList.Dimensions)(2)/((aElementList.NumElements)(2)*pow(2,tBasis_level)*2); //Z-Coordinates + offset (if needed)

    return tCoordinates;
}

Mat<real>
moris::Hierarchical_Mesh::give_middlecoordinate_from_element(
        uint & aElement,
        ElementList_struc & aElementList)
{
    Mat<real> tCoordinates(1,aElementList.Dim,0);
    uint tElementlevel = give_basis_level(aElement,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
    Mat<uint> tElement_position = give_position_of_element(aElement,aElementList.Dim,aElementList.NumElements);     // Get the Position of the element
    //Calculate the Coordinates in the middle  of the element with respect to the dimensions and the offset
    tCoordinates(0) = (tElement_position(0)*(aElementList.Dimensions)(0))/((aElementList.NumElements)(0)*pow(2,tElementlevel)) + (aElementList.Dimensions_Offset)(0) + (aElementList.Dimensions)(0)/((aElementList.NumElements)(0)*pow(2,tElementlevel)*2); //X-coordinates + offset (if needed);
    tCoordinates(1) = (tElement_position(1)*(aElementList.Dimensions)(1))/((aElementList.NumElements)(1)*pow(2,tElementlevel)) + (aElementList.Dimensions_Offset)(1) + (aElementList.Dimensions)(1)/((aElementList.NumElements)(1)*pow(2,tElementlevel)*2); //Y-Coordinates + offset (if needed);
    if (aElementList.Dim == 3)
        tCoordinates(2) = (tElement_position(2)*(aElementList.Dimensions)(2))/((aElementList.NumElements)(2)*pow(2,tElementlevel)) + (aElementList.Dimensions_Offset)(2) + (aElementList.Dimensions)(2)/((aElementList.NumElements)(2)*pow(2,tElementlevel)*2); //Z-Coordinates + offset (if needed)
    return tCoordinates;
}

void moris::Hierarchical_Mesh::give_deactivated_elements(
        ElementList_struc & aElementList)
{
    //    uint tProcSize = par_size();
    //        uint tProcRank = par_rank();
    uint tVar = 0; // Temporary variable for loop
    uint tLevel = 0;
    Mat<uint> tNeighbour(pow(3,aElementList.Dim),1); //Neighbours from an element
    Mat<uint> tBasis(pow(aElementList.Polynomial+1,aElementList.Dim),1); //Basis functions
    if ( isempty(aElementList.DeactivateElement) == 0)
    {
        uint tMaxElement = (aElementList.DeactivateElement).max();
        Mat<uint> tChildren = give_children_of_element(tMaxElement,aElementList.Dim,aElementList.NumElements);
        tMaxElement = tChildren.max();
        tLevel = give_element_level(tMaxElement,aElementList.Dim,aElementList.NumElements);
        uint tNumberElements = give_number_of_elements(tLevel,aElementList.Dim,aElementList.NumElements);
        if( (aElementList.ElementActive).size() < tNumberElements )
            (aElementList.ElementActive).resize(tNumberElements);

        if( aElementList.RefinementLevel == 1 && aElementList.MapDeactivateElementsToCoarsest == 1)
        {
            uint tParent;
            Mat<uint> tParentElements((aElementList.DeactivateElement).length(),1,UINT_MAX);
            for(uint i=0; i<(aElementList.DeactivateElement).length(); i++)
            {
                tParent = give_parent_of_level_x((aElementList.DeactivateElement)(i),aElementList.Dim,aElementList.NumElements,tLevel);
                if( tParent != UINT_MAX)
                {
                    tParentElements(i) = tParent;
                    (aElementList.ElementActive).reset(tParent);
                }
                else
                {
                    tParentElements(i) = (aElementList.DeactivateElement)(i);
                    (aElementList.ElementActive).reset((aElementList.DeactivateElement)(i));
                }
            }
            if( aElementList.BufferElements > 0)
            {
                Mat<uint> tNeighbour(pow(3,aElementList.Dim),1);
                (aElementList.DeactivateElement).set_size(tParentElements.length()*tNeighbour.length()*aElementList.BufferElements,1,UINT_MAX);
                (aElementList.DeactivateElement).rows(0,tParentElements.length()-1) = tParentElements.rows(0,tParentElements.length()-1);
                tVar = tParentElements.length();
                for(uint i = 0; i<aElementList.BufferElements; i++)
                {
                    for(uint j = 0; j<tParentElements.length(); j++)
                    {
                        tNeighbour = give_neighbour_of_element(tParentElements(j),aElementList.Dim,aElementList.BufferElements,aElementList.NumElements);
                        for(uint k = 0; k<tNeighbour.length(); k++)
                        {
                            if( (aElementList.ElementActive).test(tNeighbour(k)) == 1 )
                            {
                                (aElementList.DeactivateElement)(tVar) = tNeighbour(k);
                                (aElementList.ElementActive).reset(tNeighbour(k));
                                tVar++;
                            }
                        }
                    }
                }
                (aElementList.DeactivateElement).resize(tVar,1);
            }
            (aElementList.DeactivateElement)= unique(tParentElements);
            (aElementList.DeactivateElementInit) = (aElementList.DeactivateElement);
        }

        for(uint i=0; i<(aElementList.DeactivateElement).length(); i++)
        {
            (aElementList.ElementActive).reset((aElementList.DeactivateElement)(i)); // Deactivate elements
        }
        tVar = (aElementList.DeactivateElement).length();
        (aElementList.DeactivateElement).resize(tVar+(aElementList.ElementActive).count(),1);
        (aElementList.DeactivateElementRefLvl).resize(tVar+(aElementList.ElementActive).count(),1);
    }

    // Broadcast the deactivated elements to all procs for the first refinement level (It is needed to find dependencies at the border between procs)
    if( aElementList.RefinementLevel == 1)
        this->sendrecv_deactivation_and_bcast_element_data(aElementList);

    if ( isempty(aElementList.DeactivateElement) == 0)
    {
        this->give_active_elements(aElementList);
        //Check for linar dependencies
        if( aElementList.Dim == 2 && isempty(aElementList.ElementListOnProc) == 0 )
        {
            for(uint i=0; i<(aElementList.ElementListOnProc).length(); i++)
            {
                if( (aElementList.ElementActive).test((aElementList.ElementListOnProc)(i)) == 1 )
                {
                    uint tBuffer = 1; // Only one layer of neighbour elements are needed
                    tNeighbour = give_neighbour_of_element((aElementList.ElementListOnProc)(i),aElementList.Dim,tBuffer,aElementList.NumElements);
                    tBasis = give_basis_of_element((aElementList.ElementListOnProc)(i),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
                    if( (aElementList.ElementActive).test(tNeighbour(3)) == 0 && (aElementList.ElementActive).test(tNeighbour(5)) == 0 && (aElementList.ElementActive).test(tNeighbour(4)) == 1)
                    {
                        (aElementList.DeactivateElement)(tVar) = tNeighbour(4);
                        (aElementList.ElementActive).reset(tNeighbour(4));
                        (aElementList.DeactivateElementRefLvl)(tVar) = 0;
                        tVar++;
                    }
                    if( (aElementList.ElementActive).test(tNeighbour(1)) == 0 && (aElementList.ElementActive).test(tNeighbour(7)) == 0 && (aElementList.ElementActive).test(tNeighbour(4)) == 1)
                    {
                        (aElementList.DeactivateElement)(tVar) = tNeighbour(4);
                        (aElementList.ElementActive).reset(tNeighbour(4));
                        (aElementList.DeactivateElementRefLvl)(tVar) = 0;
                        tVar++;
                    }
                }
            }
        }
        else if( aElementList.Dim == 3 && isempty(aElementList.ElementListOnProc) == 0)
        {
            for(uint i=0; i<(aElementList.ElementListOnProc).length(); i++)
            {
                if( (aElementList.ElementActive).test((aElementList.ElementListOnProc)(i)) == 1 )
                {
                    uint tBuffer = 1; // Only one layer of neighbour elements are needed
                    tNeighbour = give_neighbour_of_element((aElementList.ElementListOnProc)(i),aElementList.Dim,tBuffer,aElementList.NumElements);
                    tBasis = give_basis_of_element((aElementList.ElementListOnProc)(i),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
                    if( (aElementList.ElementActive).test(tNeighbour(12)) == 0 && (aElementList.ElementActive).test(tNeighbour(14)) == 0 && (aElementList.ElementActive).test(tNeighbour(13)) == 1)
                    {
                        (aElementList.DeactivateElement)(tVar) = tNeighbour(13);
                        (aElementList.ElementActive).reset(tNeighbour(13));
                        (aElementList.DeactivateElementRefLvl)(tVar) = 0;
                        tVar++;
                    }
                    if( (aElementList.ElementActive).test(tNeighbour(10)) == 0 && (aElementList.ElementActive).test(tNeighbour(16)) == 0 && (aElementList.ElementActive).test(tNeighbour(13)) == 1)
                    {
                        (aElementList.DeactivateElement)(tVar) = tNeighbour(13);
                        (aElementList.ElementActive).reset(tNeighbour(13));
                        (aElementList.DeactivateElementRefLvl)(tVar) = 0;
                        tVar++;
                    }
                    if( (aElementList.ElementActive).test(tNeighbour(4)) == 0 && (aElementList.ElementActive).test(tNeighbour(22)) == 0 && (aElementList.ElementActive).test(tNeighbour(13)) == 1)
                    {
                        (aElementList.DeactivateElement)(tVar) = tNeighbour(13);
                        (aElementList.ElementActive).reset(tNeighbour(13));
                        (aElementList.DeactivateElementRefLvl)(tVar) = 0;
                        tVar++;
                    }
                }
            }
        }
    }
    if ( isempty(aElementList.DeactivateElement) == 0)
    {
        (aElementList.DeactivateElement).resize(tVar,1);
        (aElementList.DeactivateElementRefLvl).resize(tVar,1);
    }
}

void moris::Hierarchical_Mesh::hierarchical_element_refinement(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    (aElementList.DeactivateElement) = unique((aElementList.DeactivateElement));
    uint tMaxElement = (aElementList.DeactivateElement).max();
    Mat<uint> tChildren = give_children_of_element(tMaxElement,aElementList.Dim,aElementList.NumElements);
    tMaxElement = tChildren.max();
    uint tLevel = give_element_level(tMaxElement,aElementList.Dim,aElementList.NumElements);
    uint tNumberElements = give_number_of_elements(tLevel,aElementList.Dim,aElementList.NumElements);
    if( (aElementList.ElementActive).size() < tNumberElements )
        (aElementList.ElementActive).resize(tNumberElements);

    if ( aElementList.Dim == 2)
    {
        for(uint i=0; i<(aElementList.DeactivateElement).length(); i++)
        {
            tChildren=give_children_of_element((aElementList.DeactivateElement)(i),aElementList.Dim,aElementList.NumElements); //Children of the deactivated element
            (aElementList.ElementActive).set(tChildren(0));
            (aElementList.ElementActive).set(tChildren(1));
            (aElementList.ElementActive).set(tChildren(2));
            (aElementList.ElementActive).set(tChildren(3));
            (aElementList.ElementActive).reset((aElementList.DeactivateElement)(i));
        }
    }
    else if ( aElementList.Dim == 3)
    {
        for(uint i=0; i<(aElementList.DeactivateElement).length(); i++)
        {
            tChildren=give_children_of_element((aElementList.DeactivateElement)(i),aElementList.Dim,aElementList.NumElements); //Children of the deactivated element
            (aElementList.ElementActive).set(tChildren(0));
            (aElementList.ElementActive).set(tChildren(1));
            (aElementList.ElementActive).set(tChildren(2));
            (aElementList.ElementActive).set(tChildren(3));
            (aElementList.ElementActive).set(tChildren(4));
            (aElementList.ElementActive).set(tChildren(5));
            (aElementList.ElementActive).set(tChildren(6));
            (aElementList.ElementActive).set(tChildren(7));
            (aElementList.ElementActive).reset((aElementList.DeactivateElement)(i));
        }
    }
    tMaxElement = (aElementList.DeactivateElement).max();
    tLevel = give_element_level(tMaxElement,aElementList.Dim,aElementList.NumElements)+1; // Give current level
    if ((aElementList.Level)< tLevel)
    {
        (aElementList.Level) = tLevel;
        (aElementList.NumberElements) = give_number_of_elements(tLevel,aElementList.Dim,aElementList.NumElements);
        (aBasisList.NumberBasis) = give_number_of_basis(tLevel,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
    }
}

void moris::Hierarchical_Mesh::hierarchical_basisfunction_refinement(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    this->give_active_elements_with_aura(aElementList); // One layer of elements from the neighbor procs are taken into account
    Mat<uint> tBasis_of_element(pow(aElementList.Polynomial+1,aElementList.Dim),1);
    Mat<uint> tPossibleBasisFunction((aElementList.ElementListOnProc).length()*pow(aElementList.Polynomial+1,aElementList.Dim),1);
    uint tVarb = 0; // Temporary variable for loop
    for(uint i=0; i< (aElementList.ElementListOnProc).length(); i++)
    {
        tBasis_of_element=give_basis_of_element((aElementList.ElementListOnProc)(i),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
        tPossibleBasisFunction.rows(tVarb*pow(aElementList.Polynomial+1,aElementList.Dim),(tVarb+1)*pow(aElementList.Polynomial+1,aElementList.Dim)-1) = tBasis_of_element.rows(0,pow(aElementList.Polynomial+1,aElementList.Dim)-1);
        tVarb++;
    }
    tPossibleBasisFunction=unique(tPossibleBasisFunction);
    (aBasisList.BasisActiveList).set_size(tPossibleBasisFunction.length(),1,0);
    uint tVarActiveBasis = 0;
    uint tMaxBasis = tPossibleBasisFunction.max();
    uint tLevel = give_basis_level(tMaxBasis,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
    uint tBasisTotal = give_number_of_basis(tLevel,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
    (aBasisList.BasisActive).resize(tBasisTotal);
    (aBasisList.BasisActive).reset();
    uint tParent = 0;
    uint tBasis_level;
    uint tVar = 0; // Temporary variable for loop
    Mat<uint> tWhichElementsActive(pow(aElementList.Polynomial+1,aElementList.Dim),1);
    Mat<uint> tElement_of_basis(pow(aElementList.Polynomial+1,aElementList.Dim),1);
    Mat<uint> tBasisOfElement;
    Mat<uint> tListOfDeactivatedBasis(pow(aElementList.Polynomial+1,aElementList.Dim)*pow(aElementList.Polynomial+1,aElementList.Dim),1);
    Mat<uint> tCountsOfNumbers; // Provides the counts of the number of values
    uint tHelp; // Temporary variable for the loop
    uint tSwitch = 0;
    uint tFurtherLevel;
    uint tVarTemp = 0;
    uint tParentOfParent;
    uint tSwitchOfParent;

    for(uint i=0; i<tPossibleBasisFunction.length(); i++) // Update basis functions from the list of possible basis functions
    {
        tBasis_level = give_basis_level(tPossibleBasisFunction(i),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
        if( tBasis_level == 0 ) // Check the level of the element i
        {
            tElement_of_basis=give_element_of_basis(tPossibleBasisFunction(i),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
            for(uint j=0; j<tElement_of_basis.length(); j++) // Update basis functions on each level, except level 0
            {
                tWhichElementsActive(j) = (aElementList.ElementActive).test(tElement_of_basis(j));
            }
            if ( sum( tWhichElementsActive ) > 0 )
            {
                (aBasisList.BasisActive).set(tPossibleBasisFunction(i));
                (aBasisList.BasisActiveList)(tVarActiveBasis) = tPossibleBasisFunction(i);
                tVarActiveBasis++;
            }
        }
        else
        {
            break;
        }
    }

    for(uint i=0; i<tPossibleBasisFunction.length(); i++) // Update basis functions from the list of possible basis functions
    {
        tSwitch = 0;
        tFurtherLevel = 0;
        tBasis_level = give_basis_level(tPossibleBasisFunction(i),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
        if( tBasis_level > 0) // Check the level of the element i
        {
            tElement_of_basis=give_element_of_basis(tPossibleBasisFunction(i),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
            for(uint j=0; j<tElement_of_basis.length(); j++) // Update basis functions on each level, except level 0
            {
                tWhichElementsActive(j) = (aElementList.ElementActive).test(tElement_of_basis(j));
            }
            if ( sum( tWhichElementsActive ) > 0 )
            {
                tListOfDeactivatedBasis = tListOfDeactivatedBasis*0;
                while( tSwitch < 1)
                {
                    tSwitchOfParent = 0;
                    tVar = 0; // set to zero for the loop below
                    for(uint j=0; j<tElement_of_basis.length(); j++) // Update basis functions on each level, except level 0
                    {
                        tParent=give_parent_of_element(tElement_of_basis(j),aElementList.Dim,aElementList.NumElements); // Parent of element
                        for(uint k = 0; k<tFurtherLevel; k++)
                        {
                            tParent=give_parent_of_element(tParent,aElementList.Dim,aElementList.NumElements); // Parent of element
                        }
                        tParentOfParent =give_parent_of_element(tParent,aElementList.Dim,aElementList.NumElements); // Parent of parent element
                        if( tParent == UINT_MAX )
                        {
                            tSwitch = 1;
                            tVar = tVarTemp;
                            break;
                        }
                        tBasisOfElement = give_basis_of_element(tParent,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements); // Get Basis functions of parent element
                        for(uint k=0; k<pow(aElementList.Polynomial+1,aElementList.Dim); k++) // Update basis functions on each level, except level 0
                        {
                            if( (aBasisList.BasisActive).test(tBasisOfElement(k)) == 0 && (aElementList.ElementActive).test(tParent) == 0 )
                            {
                                tListOfDeactivatedBasis(tVar) = tBasisOfElement(k);
                                tVar++;
                            }
                            else
                            {
                                tSwitch = 1;
                            }
                        }
                        if( tParentOfParent != UINT_MAX && (aElementList.ElementActive).test(tParentOfParent) == 1)
                        {
                            tSwitchOfParent = 1;
                        }
                    }
                    tFurtherLevel++;
                    tVarTemp = tVar;
                    if( tSwitchOfParent == 1)
                    {
                        tSwitch = 0;
                    }
                }
                if( tVar > 0 && tSwitch == 1)
                {
                    Mat<uint> ReducedListOfDeactivatedBasis = tListOfDeactivatedBasis.rows(0,tVar-1);
                    if (unique(ReducedListOfDeactivatedBasis).length() == 1)
                    {
                        ReducedListOfDeactivatedBasis.resize(ReducedListOfDeactivatedBasis.length()+1,1);
                    }
                    tCountsOfNumbers = histc(ReducedListOfDeactivatedBasis,unique(ReducedListOfDeactivatedBasis)); // Provides the counts of the number of values
                    Mat<uint> tHelpMat = (tCountsOfNumbers == pow(aElementList.Polynomial+1,aElementList.Dim));
                    tHelp = sum(tHelpMat);
                    if ( tHelp > 0) // If all Elements with the support of the curent basis function  are active or passive and at least one element is active, then the basis function is active
                    {
                        (aBasisList.BasisActive).set(tPossibleBasisFunction(i));
                        (aBasisList.BasisActiveList)(tVarActiveBasis) = tPossibleBasisFunction(i);
                        tVarActiveBasis++;
                    }
                }
            }
        }
    }
    (aBasisList.BasisActiveList).resize(tVarActiveBasis,1);
}

void moris::Hierarchical_Mesh::hierarchical_basisfunction_designvariables_refinement(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    //    this->give_active_elements(aElementList);
    this->give_active_elements_with_aura(aElementList); // One layer of elements from the neighbor procs are taken into account
    Mat<uint> tBasis_of_element(pow(aElementList.PolynomialDesign+1,aElementList.Dim),1);
    Mat<uint> tPossibleBasisFunction((aElementList.ElementListOnProc).length()*pow(aElementList.PolynomialDesign+1,aElementList.Dim),1);
    uint tVarb = 0; // Temporary variable for loop

    for(uint i=0; i< (aElementList.ElementListOnProc).length(); i++)
    {
        tBasis_of_element=give_basis_of_element((aElementList.ElementListOnProc)(i),aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
        tPossibleBasisFunction.rows(tVarb*pow(aElementList.PolynomialDesign+1,aElementList.Dim),(tVarb+1)*pow(aElementList.PolynomialDesign+1,aElementList.Dim)-1) = tBasis_of_element.rows(0,pow(aElementList.PolynomialDesign+1,aElementList.Dim)-1);
        tVarb++;
    }
    (aElementList.ElementListActiveDesign) = (aElementList.ElementListOnProc);
    tPossibleBasisFunction=unique(tPossibleBasisFunction);
    (aBasisList.BasisActiveDesignList).set_size(tPossibleBasisFunction.length(),1,0);
    uint tMaxBasis = tPossibleBasisFunction.max();
    uint tLevel = give_basis_level(tMaxBasis,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
    uint tBasisTotal = give_number_of_basis(tLevel,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
    (aBasisList.DesignBSplineActive).resize(tBasisTotal);
    (aBasisList.DesignBSplineActive).reset();
    uint tParent = 0;
    uint tBasis_level;
    uint tVar = 0; // Temporary variable for loop
    Mat<uint> tWhichElementsActive(pow(aElementList.PolynomialDesign+1,aElementList.Dim),1);
    Mat<uint> tElement_of_basis(pow(aElementList.PolynomialDesign+1,aElementList.Dim),1);
    Mat<uint> tBasisOfElement;
    Mat<uint> tListOfDeactivatedBasis(pow(aElementList.PolynomialDesign+1,aElementList.Dim)*pow(aElementList.PolynomialDesign+1,aElementList.Dim),1);
    Mat<uint> tCountsOfNumbers; // Provides the counts of the number of values
    uint tHelp; // Temporary variable for the loop
    uint tSwitch = 0;
    uint tFurtherLevel;
    uint tVarTemp = 0;
    uint tParentOfParent;
    uint tSwitchOfParent;
    uint tVarActiveBasis = 0;
    for(uint i=0; i<tPossibleBasisFunction.length(); i++) // Update basis functions from the list of possible basis functions
    {
        //        if( tPossibleBasisFunction(i) == 105)
        //            std::cout << " tPossibleBasisFunction(i) " << tPossibleBasisFunction(i) << std::endl;
        //        if( tPossibleBasisFunction(i) == 106)
        //            std::cout << " tPossibleBasisFunction(i) " << tPossibleBasisFunction(i) << std::endl;
        //        if( tPossibleBasisFunction(i) == 138)
        //            std::cout << " tPossibleBasisFunction(i) " << tPossibleBasisFunction(i) << std::endl;
        //        if( tPossibleBasisFunction(i) == 139)
        //            std::cout << " tPossibleBasisFunction(i) " << tPossibleBasisFunction(i) << std::endl;

        tBasis_level = give_basis_level(tPossibleBasisFunction(i),aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
        if( tBasis_level == 0 ) // Check the level of the element i
        {
            tElement_of_basis=give_element_of_basis(tPossibleBasisFunction(i),aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
            for(uint j=0; j<tElement_of_basis.length(); j++) // Update basis functions on each level, except level 0
            {
                tWhichElementsActive(j) = (aElementList.ElementActive).test(tElement_of_basis(j));
            }
            if ( sum( tWhichElementsActive ) > 0 )
            {
                (aBasisList.DesignBSplineActive).set(tPossibleBasisFunction(i));
                (aBasisList.BasisActiveDesignList)(tVarActiveBasis) = tPossibleBasisFunction(i);
                tVarActiveBasis++;
            }
        }
        else
        {
            break;
        }
    }
    for(uint i=0; i<tPossibleBasisFunction.length(); i++) // Update basis functions from the list of possible basis functions
    {
        //        if( tPossibleBasisFunction(i) == 831)
        //            std::cout << " tPossibleBasisFunction(i) " << tPossibleBasisFunction(i) << std::endl;
        //        if( tPossibleBasisFunction(i) == 832)
        //            std::cout << " tPossibleBasisFunction(i) " << tPossibleBasisFunction(i) << std::endl;
        //        if( tPossibleBasisFunction(i) == 896)
        //            std::cout << " tPossibleBasisFunction(i) " << tPossibleBasisFunction(i) << std::endl;
        //        if( tPossibleBasisFunction(i) == 897)
        //            std::cout << " tPossibleBasisFunction(i) " << tPossibleBasisFunction(i) << std::endl;
        tSwitch = 0;
        tFurtherLevel = 0;
        tBasis_level = give_basis_level(tPossibleBasisFunction(i),aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
        if( tBasis_level > 0) // Check the level of the element i
        {
            tElement_of_basis=give_element_of_basis(tPossibleBasisFunction(i),aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
            for(uint j=0; j<tElement_of_basis.length(); j++) // Update basis functions on each level, except level 0
            {
                tWhichElementsActive(j) = (aElementList.ElementActive).test(tElement_of_basis(j));
            }
            if ( sum( tWhichElementsActive ) > 0 )
            {
                tListOfDeactivatedBasis = tListOfDeactivatedBasis*0;
                while( tSwitch < 1)
                {
                    tSwitchOfParent = 0;
                    tVar = 0; // set to zero for the loop below
                    for(uint j=0; j<tElement_of_basis.length(); j++) // Update basis functions on each level, except level 0
                    {
                        tParent=give_parent_of_element(tElement_of_basis(j),aElementList.Dim,aElementList.NumElements); // Parent of element
                        for(uint k = 0; k<tFurtherLevel; k++)
                        {
                            tParent=give_parent_of_element(tParent,aElementList.Dim,aElementList.NumElements); // Parent of element
                        }
                        tParentOfParent =give_parent_of_element(tParent,aElementList.Dim,aElementList.NumElements); // Parent of parent element
                        if( tParent == UINT_MAX )
                        {
                            tSwitch = 1;
                            tVar = tVarTemp;
                            break;
                        }
                        tBasisOfElement = give_basis_of_element(tParent,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements); // Get Basis functions of parent element
                        for(uint k=0; k<pow(aElementList.PolynomialDesign+1,aElementList.Dim); k++) // Update basis functions on each level, except level 0
                        {
                            if( (aBasisList.DesignBSplineActive).test(tBasisOfElement(k)) == 0 && (aElementList.ElementActive).test(tParent) == 0 )
                            {
                                tListOfDeactivatedBasis(tVar) = tBasisOfElement(k);
                                tVar++;
                            }
                            else
                            {
                                tSwitch = 1;
                            }
                        }

                        if( tParentOfParent != UINT_MAX && (aElementList.ElementActive).test(tParentOfParent) == 1)
                        {
                            tSwitchOfParent = 1;
                        }
                    }
                    tFurtherLevel++;
                    tVarTemp = tVar;
                    if( tSwitchOfParent == 1)
                    {
                        tSwitch = 0;
                    }
                }
                if( tVar > 0 && tSwitch == 1)
                {
                    Mat<uint> ReducedListOfDeactivatedBasis = tListOfDeactivatedBasis.rows(0,tVar-1);
                    if (unique(ReducedListOfDeactivatedBasis).length() == 1)
                    {
                        ReducedListOfDeactivatedBasis.resize(ReducedListOfDeactivatedBasis.length()+1,1);
                    }
                    tCountsOfNumbers = histc(ReducedListOfDeactivatedBasis,unique(ReducedListOfDeactivatedBasis)); // Provides the counts of the number of values
                    Mat<uint> tHelpMat = (tCountsOfNumbers == pow(aElementList.PolynomialDesign+1,aElementList.Dim));
                    tHelp = sum(tHelpMat);
                    if ( tHelp > 0) // If all Elements with the support of the curent basis function  are active or passive and at least one element is active, then the basis function is active
                    {
                        (aBasisList.DesignBSplineActive).set(tPossibleBasisFunction(i));
                        (aBasisList.BasisActiveDesignList)(tVarActiveBasis) = tPossibleBasisFunction(i);
                        tVarActiveBasis++;
                    }
                }
            }
        }
    }
    (aBasisList.BasisActiveDesignList).resize(tVarActiveBasis,1);
}

void moris::Hierarchical_Mesh::hierarhical_mesh_refinement(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    this->sendrecv_refinement_and_bcast_element_data(aElementList);

    if(aElementList.Refinement > 0)
    {
        if(aElementList.Refinement == 1)
        {
            this->Find_elements_in_stencil(aElementList,aBasisList);
        }
        this->give_deactivated_elements(aElementList); // Deactivate Elements and basis functions
        this->sendrecv_deactivation_and_bcast_element_data(aElementList);

        if ( isempty(aElementList.DeactivateElement) == 0)
            this->hierarchical_element_refinement(aElementList,aBasisList);
    }

    this->gather_value_and_bcast_max(aElementList.Level);
    (aElementList.NumberElements) = this->give_number_of_elements(aElementList.Level,aElementList.Dim,aElementList.NumElements);
    (aBasisList.NumberBasis) = this->give_number_of_basis(aElementList.Level,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
    this->sendrecv_refinement_and_bcast_element_data(aElementList);

    this->hierarchical_basisfunction_refinement(aElementList,aBasisList);
}

void moris::Hierarchical_Mesh::create_data_for_mesh(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    this->update_side_node_set(aElementList,aBasisList);
    this->create_mesh_data(aElementList,aBasisList);
    this->give_node_proc_owner(aElementList,aBasisList);
}

void moris::Hierarchical_Mesh::create_dofconnectivity_outputfile(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();

    // Write Dof connectivity for FEMDOC
    Mat<uint> tOrder(pow((aElementList.Polynomial)+1,(aElementList.Dim)),1); // Change ordering to the classical order of FE connectivity
    uint tVar = 0;
    Mat<real> tHelp;
    Mat<uint> tHelpList;
    std::ofstream DofNodeConnecitivity;
    if(  aElementList.Refinement > 0)
    {
        if( tProcSize == 1)
        {
            DofNodeConnecitivity.open ("ElementalDofNodeConnecitivity.data");
        }
        else
        {
            DofNodeConnecitivity.open ("ElementalDofNodeConnecitivity.data_" + std::to_string(tProcRank));
        }
    }
    else
    {
        if( tProcSize == 1)
        {
            DofNodeConnecitivity.open ("ElementalDofNodeConnecitivity.data.org");
        }
        else
        {
            DofNodeConnecitivity.open ("ElementalDofNodeConnecitivity.data.org_" + std::to_string(tProcRank));
        }
    }

    // Create a consecutive numbering
    (aElementList.ConsecutiveNumbering).set_size((aElementList.NodalLocaltoGlobal).max()+1,1,UINT_MAX);

    for(uint i=0; i < (aElementList.NodalLocaltoGlobal).length(); i++)
        (aElementList.ConsecutiveNumbering)((aElementList.NodalLocaltoGlobal)(i)) = i+1;

    for(uint i=0; i < (aElementList.NodalLocaltoGlobal).length(); i++)
        (aElementList.NodalLocaltoGlobal)(i) = i+1;

    for(uint i = 0; i < (aElementList.Fetopo).n_rows(); i++)
    {
        for(uint j = 0; j < (aElementList.Fetopo).n_cols(); j++)
        {
            (aElementList.Fetopo)(i,j) =  (aElementList.ConsecutiveNumbering)((aElementList.Fetopo)(i,j));
        }
    }
    // End of creating the consecutive numbering
    for(uint i=0; i<(aElementList.ElemLocaltoGlobal).length(); i++)
    {

        if( aElementList.TruncatedBsplines == false)
        {
            this->give_Tmatrix_and_id_field((aElementList.ElementListOnProc)(i),aElementList,aBasisList);
        }
        else
        {
            this->give_Truncated_Tmatrix_and_id_field((aElementList.ElementListOnProc)(i),aElementList,aBasisList);
        }

        Mat<real> tTMatrix = (aElementList.TMatrix);
        tHelpList.set_size(tTMatrix.n_rows(),1,0);
        tVar = 0;
        for(uint j=0; j<(aElementList.IdField).length(); j++)
        {
            tHelp = tTMatrix.row(j);
            if( sum( tHelp ) != 0.0 )
            {
                tHelpList(tVar) = j;
                tVar++;
            }
        }
        tHelpList.resize(tVar,1);
        DofNodeConnecitivity << (aElementList.ElemLocaltoGlobal)(i) << "\n";
        DofNodeConnecitivity <<  tHelpList.length() << " " << pow((aElementList.Polynomial+1),aElementList.Dim) << "\n";

        MORIS_ASSERT( (aElementList.IdField).length() >= pow((aElementList.Polynomial+1),aElementList.Dim), "AbsDesVariables.dat and ActiveDesignBasis.data have not the same number of entries ." );

        for(uint j=0; j<tHelpList.length(); j++)
        {
            for(uint k=0; k<pow((aElementList.Polynomial+1),aElementList.Dim); k++)
            {
                DofNodeConnecitivity << tTMatrix(tHelpList(j),k) << " ";
            }
            DofNodeConnecitivity << "\n";
        }
        for(uint j=0; j<tHelpList.length(); j++)
        {
            DofNodeConnecitivity <<  (aElementList.ConsecutiveNumbering)((aElementList.IdField)(tHelpList(j))) << "\n";
        }
    }
    DofNodeConnecitivity.close();
}

void moris::Hierarchical_Mesh::create_designvariables_outputfile(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    uint tVar = 0;
    Mat<uint> DesignBSplineList((aElementList.NodalLocaltoGlobal).length()*(aElementList.IdFieldFieldDesign).n_cols(),1,UINT_MAX);
    for(uint i = 0; i<(aElementList.NodalLocaltoGlobal).length(); i++) //(BasisList.BasisActive).size()
    {
        for(uint j = 0; j< (aElementList.IdFieldFieldDesign)(i,0); j++)
        {
            DesignBSplineList(tVar) = (aElementList.IdFieldFieldDesign)(i,j+1);
            tVar++;
        }
    }
    DesignBSplineList = unique(DesignBSplineList);
    DesignBSplineList.resize(DesignBSplineList.length()-1,1);

    (aBasisList.DesignBSplineListMap).set_size(DesignBSplineList.max()+1,1,0);
    for(uint i = 0; i<DesignBSplineList.length(); i++) //(BasisList.BasisActive).size()
    {
        (aBasisList.DesignBSplineListMap)(DesignBSplineList(i)) = i;
    }

    // Write Dof connectivity for FEMDOC
    std::ofstream DesignVariables;
    if(  aElementList.Refinement > 0)
    {
        if( tProcSize == 1)
        {
            DesignVariables.open ("DesignVariables.data");
        }
        else
        {
            DesignVariables.open ("DesignVariables.data_" + std::to_string(tProcRank));
        }
    }
    else
    {
        if( tProcSize == 1)
        {
            DesignVariables.open ("DesignVariables.data.org");
        }
        else
        {
            DesignVariables.open ("DesignVariables.data.org_" + std::to_string(tProcRank));
        }
    }
    for(uint i=0; i<(aElementList.NodalLocaltoGlobal).length(); i++)
    {
        DesignVariables << i+1 << " " << 2*(aElementList.TMatrixFieldDesign)(i,0) << "\n";
        for(uint j = 0; j< (aElementList.IdFieldFieldDesign)(i,0); j++)
        {
            DesignVariables << (aBasisList.DesignBSplineListMap)((aElementList.IdFieldFieldDesign)(i,j+1)) << " ";
        }
        for(uint j = 0; j< (aElementList.IdFieldFieldDesign)(i,0); j++)
        {
            DesignVariables << (aElementList.TMatrixFieldDesign)(i,j+1) << " ";
        }
        DesignVariables << "\n";
    }
    DesignVariables.close();
    // Write Dof connectivity for FEMDOC
    std::ofstream DesignVariablesMoris;
    if(  aElementList.Refinement > 0)
    {
        if( tProcSize == 1)
        {
            DesignVariablesMoris.open ("DesignVariables_for_moris_new.data");
        }
        else
        {
            DesignVariablesMoris.open ("DesignVariables_for_moris_new.data_" + std::to_string(tProcRank));
        }
    }
    else
    {
        if( tProcSize == 1)
        {
            DesignVariablesMoris.open ("DesignVariables_for_moris.data.org");
        }
        else
        {
            DesignVariablesMoris.open ("DesignVariables_for_moris.data.org_" + std::to_string(tProcRank));
        }
    }
    DesignVariablesMoris << (aElementList.IdFieldFieldDesign).n_rows() << "\n";
    DesignVariablesMoris << (aElementList.IdFieldFieldDesign).n_cols() << "\n";
    for(uint i=0; i<(aElementList.NodalLocaltoGlobal).length(); i++)
    {
        DesignVariablesMoris << (aElementList.NodalLocaltoGlobal)(i) << "\n";
        DesignVariablesMoris << 2*(aElementList.TMatrixFieldDesign)(i,0) << "\n";
        for(uint j = 0; j< (aElementList.IdFieldFieldDesign)(i,0); j++)
        {
            DesignVariablesMoris << (aBasisList.DesignBSplineListMap)((aElementList.IdFieldFieldDesign)(i,j+1)) << "\n";
        }
        for(uint j = 0; j< (aElementList.IdFieldFieldDesign)(i,0); j++)
        {
            DesignVariablesMoris << (aElementList.TMatrixFieldDesign)(i,j+1) << "\n";
        }
    }
    DesignVariablesMoris.close();
}

void moris::Hierarchical_Mesh::save_active_design_elements_in_file(
        ElementList_struc & aElementList)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    // Save active elements
    std::ofstream ActiveElements;
    if(  aElementList.Refinement > 0)
    {
        if( tProcSize == 1)
        {
            ActiveElements.open ("ActiveDesignElements_new.data");
        }
        else
        {
            ActiveElements.open ("ActiveDesignElements_new.data_" + std::to_string(tProcRank));
        }
    }
    else
    {
        if( tProcSize == 1)
        {
            ActiveElements.open ("ActiveDesignElements.data.org");
        }
        else
        {
            ActiveElements.open ("ActiveDesignElements_new.data.org_" + std::to_string(tProcRank));
        }
    }
    ActiveElements << (aElementList.ElementListActiveDesign).length()  << "\n";
    for(uint i=0; i<(aElementList.ElementListActiveDesign).length(); i++)
    {
        ActiveElements << (aElementList.ElementListActiveDesign)(i)  << "\n";
    }
}

void moris::Hierarchical_Mesh::save_active_fem_elements_in_file(
        ElementList_struc & aElementList)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    // Save active elements

    // Save active elements
    std::ofstream ActiveFEMElements;
    if(  aElementList.Refinement > 0)
    {
        if( tProcSize == 1)
        {
            ActiveFEMElements.open ("ActiveFEMElements_new.data");
        }
        else
        {
            ActiveFEMElements.open ("ActiveFEMElements_new.data_" + std::to_string(tProcRank));
        }
    }
    else
    {
        if( tProcSize == 1)
        {
            ActiveFEMElements.open ("ActiveFEMElements.data.org");
        }
        else
        {
            ActiveFEMElements.open ("ActiveFEMElements.data.org_" + std::to_string(tProcRank));
        }
    }
    ActiveFEMElements << (aElementList.ElementListOnProc).length()  << "\n";
    for(uint i=0; i<(aElementList.ElementListOnProc).length(); i++)
    {
        ActiveFEMElements << (aElementList.ElementListOnProc)(i)  << "\n";
    }
    ActiveFEMElements.close();
}

void moris::Hierarchical_Mesh::save_active_design_basis_in_file(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    // Save active design basis
    std::ofstream ActiveDesignBasis;
    if(  aElementList.Refinement > 0)
    {
        if( tProcSize == 1)
        {
            ActiveDesignBasis.open ("ActiveDesignBasis_new.data");
        }
        else
        {
            ActiveDesignBasis.open ("ActiveDesignBasis_new.data_" + std::to_string(tProcRank));
        }
    }
    else
    {
        if( tProcSize == 1)
        {
            ActiveDesignBasis.open ("ActiveDesignBasis.data.org");
        }
        else
        {
            ActiveDesignBasis.open ("ActiveDesignBasis.data.org_" + std::to_string(tProcRank));
        }
    }
    ActiveDesignBasis << (aBasisList.BasisActiveDesignList).length()  << "\n";
    for(uint i=0; i<(aBasisList.BasisActiveDesignList).length(); i++)
    {
        ActiveDesignBasis << (aBasisList.BasisActiveDesignList)(i)  << "\n";
    }
    ActiveDesignBasis.close();
}

void moris::Hierarchical_Mesh::save_active_basis_in_file(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    // Save active design basis
    std::ofstream ActiveBasis;
    if(  aElementList.Refinement > 0)
    {
        if( tProcSize == 1)
        {
            ActiveBasis.open ("ActiveBasis_new.data");
        }
        else
        {
            ActiveBasis.open ("ActiveBasis_new.data_" + std::to_string(tProcRank));
        }
    }
    else
    {
        if( tProcSize == 1)
        {
            ActiveBasis.open ("ActiveBasis.data.org");
        }
        else
        {
            ActiveBasis.open ("ActiveBasis.data.org_" + std::to_string(tProcRank));
        }
    }
    ActiveBasis << (aBasisList.BasisActiveList).length()  << "\n";
    for(uint i=0; i<(aBasisList.BasisActiveList).length(); i++)
    {
        ActiveBasis << (aBasisList.BasisActiveList)(i)  << "\n";
    }
    ActiveBasis.close();
}

void moris::Hierarchical_Mesh::save_coordlist_in_file(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    // Save active design basis
    std::ofstream CoordList;
    if(  aElementList.Refinement > 0)
    {
        if( tProcSize == 1)
        {
            CoordList.open ("CoordList_new.data");
        }
        else
        {
            CoordList.open ("CoordList_new.data_" + std::to_string(tProcRank));
        }
    }
    else
    {
        if( tProcSize == 1)
        {
            CoordList.open ("CoordList.data.org");
        }
        else
        {
            CoordList.open ("CoordList.data.org_" + std::to_string(tProcRank));
        }
    }
    CoordList << (aElementList.NodalLocaltoGlobalExist).length()  << "\n";
    for(uint i=0; i<(aElementList.NodalLocaltoGlobalExist).length(); i++)
    {
        CoordList << (aElementList.NodalLocaltoGlobalExist)(i)  << "\n";
    }
    CoordList.close();
}

void moris::Hierarchical_Mesh::save_basis_to_element_support(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    // Save active design basis
    std::ofstream BasisToElement;
    if(  aElementList.Refinement > 0)
    {
        if( tProcSize == 1)
        {
            BasisToElement.open ("BasisToElementSupport_new.data");
        }
        else
        {
            BasisToElement.open ("BasisToElementSupport_new.data_" + std::to_string(tProcRank));
        }
    }
    else
    {
        if( tProcSize == 1)
        {
            BasisToElement.open ("BasisToElementSupport.data.org");
        }
        else
        {
            BasisToElement.open ("BasisToElementSupport.data.org_" + std::to_string(tProcRank));
        }
    }
    Mat<uint> tBasisToElement((aElementList.IdFieldField).max()+1,1000,0);
    for(uint i = 0; i < (aElementList.ElemLocaltoGlobal).length(); i++)
    {
        if( aElementList.TruncatedBsplines == false)
        {
            this->give_Tmatrix_and_id_field((aElementList.ElemLocaltoGlobal)(i),aElementList,aBasisList);
        }
        else
        {
            this->give_Truncated_Tmatrix_and_id_field((aElementList.ElemLocaltoGlobal)(i),aElementList,aBasisList);
        }
        for(uint j = 0; j <  (aElementList.IdField).length(); j++)
        {
            tBasisToElement((aElementList.IdField)(j),0)++;
            tBasisToElement((aElementList.IdField)(j),tBasisToElement((aElementList.IdField)(j),0)) = (aElementList.ElemLocaltoGlobal)(i);
        }
    }
    Mat<uint> tCol = tBasisToElement.col(0);
    uint tMaxCol = tCol.max();
    tBasisToElement.resize((aElementList.IdFieldField).max()+1,tMaxCol+1);
    Mat<uint> tBasisToElementList(aElementList.TrueNodalLocaltoGlobalSize,tMaxCol+1); // Number of real nodes without dummy nodes (which are needed to create the mesh)
    for(uint i = 0; i < aElementList.TrueNodalLocaltoGlobalSize; i++)
    {
        tBasisToElementList.row(i) = tBasisToElement.row((aElementList.NodalLocaltoGlobal)(i));
    }
    BasisToElement << tBasisToElementList.n_rows()  << "\n";
    for(uint i = 0; i < tBasisToElementList.n_rows(); i++)
    {
        for(uint j = 0; j <= tBasisToElementList(i,0); j++)
        {
            BasisToElement << tBasisToElementList(i,j)  << " ";
        }
        BasisToElement  << "\n";
    }
    BasisToElement.close();
}

void moris::Hierarchical_Mesh::update_side_node_set(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    Cell<Mat<uint>> tSetsEntIds((aElementList.NumSets)(0)); // Initilize Cell for mesh output
    Cell<std::string > tSetsName((aElementList.NumSets)(0)); // Initilize Cell for mesh output
    uint tVarCell = 0;
    if( (aElementList.NumSets)(0) > 0)
    {
        if( (aBasisList.NodeSet).n_rows() > 0)
        {
            // Update NodeSet
            uint tHelp = (aBasisList.NodeSet)((aBasisList.NodeSet).n_rows()-1,0); // Which max level is in the NodeSet list
            uint Length = (aBasisList.NodeSet).n_rows(); // How many NodeSets exsists
            if ( tHelp < aElementList.Level )
            {
                (aBasisList.NodeSet).resize((aBasisList.NodeSet).n_rows()*(aElementList.Level+1),3);
                Mat<uint> tElement_of_basis(pow(aElementList.Polynomial+1,aElementList.Dim),1);
                Mat<uint> tChildren(pow(2,aElementList.Dim),1);
                Mat<uint> tBasis_of_children(pow(aElementList.Polynomial+1,aElementList.Dim),1);
                Mat<uint> tBasis_of_element(pow(aElementList.Polynomial+1,aElementList.Dim),1);
                Mat<uint> tHelpMat; // Temporary variable for loop
                uint WhichBasis; // Temporary variable for loop
                uint tVar = Length; // Temporary variable for loop
                for(uint i=0; i<Length*(aElementList.Level); i++)
                {
                    WhichBasis = (aBasisList.NodeSet)(i,1);
                    tElement_of_basis = give_element_of_basis(WhichBasis,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
                    tChildren = give_children_of_element(tElement_of_basis(0),aElementList.Dim,aElementList.NumElements);
                    tBasis_of_element= give_basis_of_element(tElement_of_basis(0),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
                    tHelpMat = (tBasis_of_element == WhichBasis);
                    tHelpMat = find(tHelpMat);
                    tHelp = tChildren(tHelpMat(0));
                    tBasis_of_children = give_basis_of_element(tHelp,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
                    tHelp = tBasis_of_children(tHelpMat(0));
                    (aBasisList.NodeSet)(tVar,0) = (aBasisList.NodeSet)(i,0)+1;
                    (aBasisList.NodeSet)(tVar,1) = tHelp;
                    (aBasisList.NodeSet)(tVar,2) = (aBasisList.NodeSet)(i,2);
                    tVar++;
                }
                (aBasisList.NodeSet).resize(tVar,3);
            }
        }
        Mat<uint> tMatNodeSet;
        for(uint i = 0; i<(aElementList.NumSets)(0); i++ )
        {
            uint tVar = 0;
            if( (aBasisList.NodeSet).n_rows() > 0 )
            {
                tMatNodeSet.set_size((aBasisList.NodeSet).n_rows(),1);
                for( uint j = 0; j < (aBasisList.NodeSet).n_rows(); j++)
                {
                    if( (aBasisList.BasisActive).test((aBasisList.NodeSet)(j,1)) == 1 && (aBasisList.NodeSet)(j,2) == i+1)
                    {
                        tMatNodeSet(tVar) = (aBasisList.NodeSet)(j,1)+1;
                        tVar++;
                    }
                }
            }
            if( tVar > 0)
            {
                tMatNodeSet.resize(tVar,1);
                tSetsEntIds(tVarCell) = tMatNodeSet;
                tSetsName(tVarCell) = "tNodeset_" + std::to_string(i+1);
                tVarCell++;
            }
            else
            {
                tSetsName(tVarCell) = "tNodeset_" + std::to_string(i+1);
                tVarCell++;
            }
        }
        (aElementList.NEntIds) = tSetsEntIds;
        (aElementList.NSetNames) = tSetsName;
    }
    //Update SideSet
    tSetsEntIds.resize((aElementList.NumSets)(1));
    tSetsName.resize((aElementList.NumSets)(1));
    tVarCell = 0;
    if( (aElementList.NumSets)(1) > 0)
    {
        if( (aElementList.SideSet).n_rows() > 0 )
        {
            uint tHelp = (aElementList.SideSet)((aElementList.SideSet).n_rows()-1,0); // Which max level is in the SideSet list
            uint Length = (aElementList.SideSet).n_rows(); // How many Side sets exsists
            uint tVar = Length; // Temporary variable for loop
            Mat<uint> tChildren(1,pow(2,aElementList.Dim),0);
            if ( tHelp < aElementList.Level )
            {
                if( aElementList.Dim == 2)
                {
                    tHelp = 0;
                    for(uint level = 0; level<= aElementList.Level; level++)
                    {
                        tHelp += pow(2,level); // 1 Sideset is always divided in 2 side sets in the children
                    }
                    (aElementList.SideSet).resize(Length*tHelp,6); // Level, Element, 4 Side sets = 6 cols
                    for(uint i=0; i<Length*(tHelp-pow(2,aElementList.Level)); i++)
                    {
                        if ( (aElementList.SideSet)(i,0) >= aElementList.Level)
                        {
                            break;
                        }
                        tChildren = give_children_of_element((aElementList.SideSet)(i,1),aElementList.Dim,aElementList.NumElements);
                        if( (aElementList.SideSet)(i,2+0) > 0 || (aElementList.SideSet)(i,2+3) > 0)
                        {
                            (aElementList.SideSet)(tVar,0) =  (aElementList.SideSet)(i,0) + 1; // Children level = Parent level + 1
                            (aElementList.SideSet)(tVar,1) = tChildren(0);
                            (aElementList.SideSet)(tVar,2+0) = (aElementList.SideSet)(i,2+0); // +2, because col(0) = Level, col(1) = Element, col(2:5) = SideSet(1:4)
                            (aElementList.SideSet)(tVar,2+1) = 0;
                            (aElementList.SideSet)(tVar,2+2) = 0;
                            (aElementList.SideSet)(tVar,2+3) = (aElementList.SideSet)(i,2+3);
                            tVar++;
                        }
                        if( (aElementList.SideSet)(i,2+0) > 0 || (aElementList.SideSet)(i,2+1) > 0)
                        {
                            (aElementList.SideSet)(tVar,0) =  (aElementList.SideSet)(i,0) + 1;
                            (aElementList.SideSet)(tVar,1) = tChildren(1);
                            (aElementList.SideSet)(tVar,2+0) = (aElementList.SideSet)(i,2+0);
                            (aElementList.SideSet)(tVar,2+1) = (aElementList.SideSet)(i,2+1);
                            (aElementList.SideSet)(tVar,2+2) = 0;
                            (aElementList.SideSet)(tVar,2+3) = 0;
                            tVar++;
                        }
                        if( (aElementList.SideSet)(i,2+2) > 0 || (aElementList.SideSet)(i,2+3) > 0)
                        {
                            (aElementList.SideSet)(tVar,0) =  (aElementList.SideSet)(i,0) + 1;
                            (aElementList.SideSet)(tVar,1) = tChildren(2);
                            (aElementList.SideSet)(tVar,2+0) = 0;
                            (aElementList.SideSet)(tVar,2+1) = 0;
                            (aElementList.SideSet)(tVar,2+2) = (aElementList.SideSet)(i,2+2);
                            (aElementList.SideSet)(tVar,2+3) = (aElementList.SideSet)(i,2+3);
                            tVar++;
                        }
                        if( (aElementList.SideSet)(i,2+1) > 0 || (aElementList.SideSet)(i,2+2) > 0)
                        {
                            (aElementList.SideSet)(tVar,0) =  (aElementList.SideSet)(i,0) + 1;
                            (aElementList.SideSet)(tVar,1) = tChildren(3);
                            (aElementList.SideSet)(tVar,2+0) = 0;
                            (aElementList.SideSet)(tVar,2+1) = (aElementList.SideSet)(i,2+1);
                            (aElementList.SideSet)(tVar,2+2) = (aElementList.SideSet)(i,2+2);
                            (aElementList.SideSet)(tVar,2+3) = 0;
                            tVar++;
                        }
                    }
                    (aElementList.SideSet).resize(tVar,6);
                }
                else if( aElementList.Dim == 3 )
                {
                    tHelp = 0;
                    for(uint level = 0; level<= aElementList.Level; level++)
                    {
                        tHelp += pow(4,level);  // 1 Sideset is always divided in 4 side sets in the children (1 Face of the cube)
                    }
                    (aElementList.SideSet).resize(Length*tHelp,8); // Level, Element, 6 Side sets = 8 cols
                    for(uint i=0; i<Length*(tHelp-pow(4,aElementList.Level)); i++)
                    {
                        if ( (aElementList.SideSet)(i,0) >= aElementList.Level)
                        {
                            break;
                        }
                        tChildren = give_children_of_element((aElementList.SideSet)(i,1),aElementList.Dim,aElementList.NumElements);
                        if( (aElementList.SideSet)(i,2+0) > 0 || (aElementList.SideSet)(i,2+3) > 0 || (aElementList.SideSet)(i,2+4) > 0)
                        {
                            (aElementList.SideSet)(tVar,0) = (aElementList.SideSet)(i,0) + 1; // Children level = Parent level + 1
                            (aElementList.SideSet)(tVar,1) = tChildren(0);
                            (aElementList.SideSet)(tVar,2+0) = (aElementList.SideSet)(i,2+0); // +2, because col(0) = Level, col(1) = Element, col(2:7) = SideSet(1:6)
                            (aElementList.SideSet)(tVar,2+1) = 0;
                            (aElementList.SideSet)(tVar,2+2) = 0;
                            (aElementList.SideSet)(tVar,2+3) = (aElementList.SideSet)(i,2+3);
                            (aElementList.SideSet)(tVar,2+4) = (aElementList.SideSet)(i,2+4);
                            (aElementList.SideSet)(tVar,2+5) = 0;
                            tVar++;
                        }
                        if( (aElementList.SideSet)(i,2+0) > 0 || (aElementList.SideSet)(i,2+1) > 0 || (aElementList.SideSet)(i,2+4) > 0)
                        {
                            (aElementList.SideSet)(tVar,0) = (aElementList.SideSet)(i,0) + 1;
                            (aElementList.SideSet)(tVar,1) = tChildren(1);
                            (aElementList.SideSet)(tVar,2+0) = (aElementList.SideSet)(i,2+0);
                            (aElementList.SideSet)(tVar,2+1) = (aElementList.SideSet)(i,2+1);
                            (aElementList.SideSet)(tVar,2+2) = 0;
                            (aElementList.SideSet)(tVar,2+3) = 0;
                            (aElementList.SideSet)(tVar,2+4) = (aElementList.SideSet)(i,2+4);
                            (aElementList.SideSet)(tVar,2+5) = 0;
                            tVar++;
                        }
                        if( (aElementList.SideSet)(i,2+2) > 0 || (aElementList.SideSet)(i,2+3) > 0 || (aElementList.SideSet)(i,2+4) > 0)
                        {
                            (aElementList.SideSet)(tVar,0) = (aElementList.SideSet)(i,0) + 1;
                            (aElementList.SideSet)(tVar,1) = tChildren(2);
                            (aElementList.SideSet)(tVar,2+0) = 0;
                            (aElementList.SideSet)(tVar,2+1) = 0;
                            (aElementList.SideSet)(tVar,2+2) = (aElementList.SideSet)(i,2+2);
                            (aElementList.SideSet)(tVar,2+3) = (aElementList.SideSet)(i,2+3);
                            (aElementList.SideSet)(tVar,2+4) = (aElementList.SideSet)(i,2+4);
                            (aElementList.SideSet)(tVar,2+5) = 0;
                            tVar++;
                        }
                        if( (aElementList.SideSet)(i,2+1) > 0 || (aElementList.SideSet)(i,2+2) > 0 || (aElementList.SideSet)(i,2+4) > 0)
                        {
                            (aElementList.SideSet)(tVar,0) = (aElementList.SideSet)(i,0) + 1;
                            (aElementList.SideSet)(tVar,1) = tChildren(3);
                            (aElementList.SideSet)(tVar,2+0) = 0;
                            (aElementList.SideSet)(tVar,2+1) = (aElementList.SideSet)(i,2+1);
                            (aElementList.SideSet)(tVar,2+2) = (aElementList.SideSet)(i,2+2);
                            (aElementList.SideSet)(tVar,2+3) = 0;
                            (aElementList.SideSet)(tVar,2+4) = (aElementList.SideSet)(i,2+4);
                            (aElementList.SideSet)(tVar,2+5) = 0;
                            tVar++;
                        }
                        if( (aElementList.SideSet)(i,2+0) > 0 || (aElementList.SideSet)(i,2+3) > 0 || (aElementList.SideSet)(i,2+5) > 0)
                        {
                            (aElementList.SideSet)(tVar,0) = (aElementList.SideSet)(i,0) + 1;
                            (aElementList.SideSet)(tVar,1) = tChildren(4);
                            (aElementList.SideSet)(tVar,2+0) = (aElementList.SideSet)(i,2+0);
                            (aElementList.SideSet)(tVar,2+1) = 0;
                            (aElementList.SideSet)(tVar,2+2) = 0;
                            (aElementList.SideSet)(tVar,2+3) = (aElementList.SideSet)(i,2+3);
                            (aElementList.SideSet)(tVar,2+4) = 0;
                            (aElementList.SideSet)(tVar,2+5) = (aElementList.SideSet)(i,2+5);
                            tVar++;
                        }
                        if( (aElementList.SideSet)(i,2+0) > 0 || (aElementList.SideSet)(i,2+1) > 0 || (aElementList.SideSet)(i,2+5) > 0)
                        {
                            (aElementList.SideSet)(tVar,0) = (aElementList.SideSet)(i,0) + 1;
                            (aElementList.SideSet)(tVar,1) = tChildren(5);
                            (aElementList.SideSet)(tVar,2+0) = (aElementList.SideSet)(i,2+0);
                            (aElementList.SideSet)(tVar,2+1) = (aElementList.SideSet)(i,2+1);
                            (aElementList.SideSet)(tVar,2+2) = 0;
                            (aElementList.SideSet)(tVar,2+3) = 0;
                            (aElementList.SideSet)(tVar,2+4) = 0;
                            (aElementList.SideSet)(tVar,2+5) = (aElementList.SideSet)(i,2+5);
                            tVar++;
                        }

                        if( (aElementList.SideSet)(i,2+2) > 0 || (aElementList.SideSet)(i,2+3) > 0 || (aElementList.SideSet)(i,2+5) > 0)
                        {
                            (aElementList.SideSet)(tVar,0) = (aElementList.SideSet)(i,0) + 1;
                            (aElementList.SideSet)(tVar,1) = tChildren(6);
                            (aElementList.SideSet)(tVar,2+0) = 0;
                            (aElementList.SideSet)(tVar,2+1) = 0;
                            (aElementList.SideSet)(tVar,2+2) = (aElementList.SideSet)(i,2+2);
                            (aElementList.SideSet)(tVar,2+3) = (aElementList.SideSet)(i,2+3);
                            (aElementList.SideSet)(tVar,2+4) = 0;
                            (aElementList.SideSet)(tVar,2+5) = (aElementList.SideSet)(i,2+5);
                            tVar++;
                        }
                        if( (aElementList.SideSet)(i,2+1) > 0 || (aElementList.SideSet)(i,2+2) > 0 || (aElementList.SideSet)(i,2+5) > 0)
                        {
                            (aElementList.SideSet)(tVar,0) = (aElementList.SideSet)(i,0) + 1;
                            (aElementList.SideSet)(tVar,1) = tChildren(7);
                            (aElementList.SideSet)(tVar,2+0) = 0;
                            (aElementList.SideSet)(tVar,2+1) = (aElementList.SideSet)(i,2+1);
                            (aElementList.SideSet)(tVar,2+2) = (aElementList.SideSet)(i,2+2);
                            (aElementList.SideSet)(tVar,2+3) = 0;
                            (aElementList.SideSet)(tVar,2+4) = 0;
                            (aElementList.SideSet)(tVar,2+5) = (aElementList.SideSet)(i,2+5);
                            tVar++;
                        }
                    }
                    (aElementList.SideSet).resize(tVar,8);
                }
            }
        }
        Mat<uint> tMatSideSet;
        for(uint i = 0; i<(aElementList.NumSets)(1); i++ )
        {
            uint tVar = 0;
            if( (aElementList.SideSet).n_rows() > 0 )
            {
                tMatSideSet.set_size((aElementList.SideSet).n_rows(),2);
                for(uint j = 0; j < (aElementList.SideSet).n_rows(); j++)
                {
                    if( (aElementList.SideSet)(j,0) > aElementList.Level) // Check only side sets, which are in the range of the highest refinement level
                        break;
                    if( (aElementList.ElementActive).test((aElementList.SideSet)(j,1)) == 1 && (aElementList.SideSet)(j,2) == i+1)
                    {
                        tMatSideSet(tVar,0) = (aElementList.SideSet)(j,1);
                        tMatSideSet(tVar,1) = 0;
                        tVar++;
                    }
                    if( (aElementList.ElementActive).test((aElementList.SideSet)(j,1)) == 1 && (aElementList.SideSet)(j,3) == i+1)
                    {
                        tMatSideSet(tVar,0) = (aElementList.SideSet)(j,1);
                        tMatSideSet(tVar,1) = 1;
                        tVar++;
                    }
                    if( (aElementList.ElementActive).test((aElementList.SideSet)(j,1)) == 1 && (aElementList.SideSet)(j,4) == i+1)
                    {
                        tMatSideSet(tVar,0) = (aElementList.SideSet)(j,1);
                        tMatSideSet(tVar,1) = 2;
                        tVar++;
                    }
                    if( (aElementList.ElementActive).test((aElementList.SideSet)(j,1)) == 1 && (aElementList.SideSet)(j,5) == i+1)
                    {
                        tMatSideSet(tVar,0) = (aElementList.SideSet)(j,1);
                        tMatSideSet(tVar,1) = 3;
                        tVar++;
                    }
                    if( aElementList.Dim == 3)
                    {
                        if( (aElementList.ElementActive).test((aElementList.SideSet)(j,1)) == 1 && (aElementList.SideSet)(j,6) == i+1)
                        {
                            tMatSideSet(tVar,0) = (aElementList.SideSet)(j,1);
                            tMatSideSet(tVar,1) = 4;
                            tVar++;
                        }
                        if( (aElementList.ElementActive).test((aElementList.SideSet)(j,1)) == 1 && (aElementList.SideSet)(j,7) == i+1)
                        {
                            tMatSideSet(tVar,0) = (aElementList.SideSet)(j,1);
                            tMatSideSet(tVar,1) = 5;
                            tVar++;
                        }
                    }
                }
            }
            if( tVar > 0)
            {
                tMatSideSet.resize(tVar,2);
                tSetsEntIds(tVarCell) = tMatSideSet;
                tSetsName(tVarCell) = "tSideSet_" + std::to_string(i);
                tVarCell++;
            }
            else
            {
                tSetsName(tVarCell) = "tSideSet_" + std::to_string(i);
                tVarCell++;
            }
        }
        (aElementList.ElemIdsAndSideOrds) = tSetsEntIds;
        (aElementList.SSetNames) = tSetsName;
    }
}

void moris::Hierarchical_Mesh::give_node_proc_owner(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    uint tBasis;
    uint testval = 0;
    uint tBasis_level;
    Mat<uint> tNodeBoundary(2,1);
    Mat<uint> tBasisPositionDummy(1,aElementList.Dim,0);
    Mat<uint> tBasisPosition(1,3,0);
    if ( tProcSize > 1)
    {
        (aBasisList.NodeProcs).set_size((aElementList.NodalLocaltoGlobal).length(),1,tProcRank);
        if ( tProcRank > 0)
        {
            for(uint i = 0; i<(aElementList.NodalLocaltoGlobal).length(); i++)
            {
                if( ((aElementList.NodalLocaltoGlobal)(i) - 1) < aBasisList.NumberBasis )
                {
                    tBasis = (aElementList.NodalLocaltoGlobal)(i) - 1;
                    tBasisPositionDummy = give_position_of_basis(tBasis,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
                    if( aElementList.Dim == 2)
                    {
                        tBasisPosition(0) =tBasisPositionDummy(0) ; tBasisPosition(1) = tBasisPositionDummy(1) ;tBasisPosition(2) = 0;
                    }
                    else if( aElementList.Dim == 3)
                    {
                        tBasisPosition = tBasisPositionDummy;
                    }
                    tBasis_level = give_basis_level(tBasis,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);

                    if((aBasisList.BasisActive).test(tBasis) == 1 )
                    {
                        if ( tBasisPosition(0) == aBasisList.Decomp(0)*pow(2,tBasis_level) && tBasisPosition(1) == aBasisList.Decomp(2)*pow(2,tBasis_level) && tBasisPosition(2) == aBasisList.Decomp(4)*pow(2,tBasis_level) && (aElementList.ProcNeighbour)(6) < UINT_MAX && (aElementList.ProcNeighbour)(6) < tProcRank) // Bottom left node
                        {
                            (aBasisList.NodeProcs)(i) = (aElementList.ProcNeighbour)(6);
                            testval++;
                        }
                        else if ( tBasisPosition(0) == aBasisList.Decomp(0)*pow(2,tBasis_level) && tBasisPosition(2) == aBasisList.Decomp(4)*pow(2,tBasis_level) && (aElementList.ProcNeighbour)(7) < UINT_MAX && (aElementList.ProcNeighbour)(7) < tProcRank) // Bottom left node
                        {
                            (aBasisList.NodeProcs)(i) = (aElementList.ProcNeighbour)(7);
                            testval++;
                        }
                        else if ( tBasisPosition(0) == aBasisList.Decomp(0)*pow(2,tBasis_level) && tBasisPosition(1) == aBasisList.Decomp(2)*pow(2,tBasis_level) && (aElementList.ProcNeighbour)(8) < UINT_MAX && (aElementList.ProcNeighbour)(8) < tProcRank) // Bottom left node
                        {
                            (aBasisList.NodeProcs)(i) = (aElementList.ProcNeighbour)(8);
                            testval++;
                        }
                        else if ( tBasisPosition(1) == aBasisList.Decomp(2)*pow(2,tBasis_level) && tBasisPosition(2) == aBasisList.Decomp(4)*pow(2,tBasis_level) && (aElementList.ProcNeighbour)(9) < UINT_MAX && (aElementList.ProcNeighbour)(9) < tProcRank) // Bottom left node
                        {
                            (aBasisList.NodeProcs)(i) = (aElementList.ProcNeighbour)(9);
                            testval++;
                        }
                        else if ( tBasisPosition(1) == aBasisList.Decomp(2)*pow(2,tBasis_level) && (aElementList.ProcNeighbour)(0) < UINT_MAX && (aElementList.ProcNeighbour)(0) < tProcRank) // Bottom
                        {
                            (aBasisList.NodeProcs)(i) = (aElementList.ProcNeighbour)(0);
                            testval++;
                        }
                        else if ( tBasisPosition(0) == aBasisList.Decomp(1)*pow(2,tBasis_level) && (aElementList.ProcNeighbour)(1) < UINT_MAX && (aElementList.ProcNeighbour)(1) < tProcRank) // Right
                        {
                            (aBasisList.NodeProcs)(i) = (aElementList.ProcNeighbour)(1);
                            testval++;
                        }
                        else if ( tBasisPosition(1) == aBasisList.Decomp(3)*pow(2,tBasis_level) && (aElementList.ProcNeighbour)(2) < UINT_MAX && (aElementList.ProcNeighbour)(2) < tProcRank) // Top
                        {
                            (aBasisList.NodeProcs)(i) = (aElementList.ProcNeighbour)(2);
                            testval++;
                        }
                        else if ( tBasisPosition(0) == aBasisList.Decomp(0)*pow(2,tBasis_level) && (aElementList.ProcNeighbour)(3) < UINT_MAX && (aElementList.ProcNeighbour)(3) < tProcRank) // Left (funkt)
                        {
                            (aBasisList.NodeProcs)(i) = (aElementList.ProcNeighbour)(3);
                            testval++;
                        }
                        else if ( tBasisPosition(2) == aBasisList.Decomp(4)*pow(2,tBasis_level) && (aElementList.ProcNeighbour)(4) < UINT_MAX && (aElementList.ProcNeighbour)(4) < tProcRank) // Front
                        {
                            (aBasisList.NodeProcs)(i) = (aElementList.ProcNeighbour)(4);
                            testval++;
                        }
                        else if ( tBasisPosition(2) == aBasisList.Decomp(5)*pow(2,tBasis_level) && (aElementList.ProcNeighbour)(5) < UINT_MAX && (aElementList.ProcNeighbour)(5) < tProcRank) // Back
                        {
                            (aBasisList.NodeProcs)(i) = (aElementList.ProcNeighbour)(5);
                            testval++;
                        }
                    }
                }
            }
        }
    }
}
moris::Mat<uint>
moris::Hierarchical_Mesh::remove_elements(ElementList_struc & aElementList)
{
    uint tNELevel = give_element_level(aElementList.NumberElements,aElementList.Dim,aElementList.NumElements);
    uint tNRemoveElements  = 0;
    uint tVar = 0; // Temporary variable for a loop
    Mat<uint> Position(aElementList.Dim,1,0); // Vector for position of element
    for(uint level=0; level<tNELevel; level++)
    {
        if( aElementList.Dim == 2)
        {
            tNRemoveElements += (2*(aElementList.NumElements)(0)+2*(aElementList.NumElements)(1))*pow(4,level);
        }
        else if( aElementList.Dim == 3)
        {
            tNRemoveElements += ((aElementList.NumElements)(0)*(aElementList.NumElements)(1)*2 + (aElementList.NumElements)(0)*(aElementList.NumElements)(2)*2 + (aElementList.NumElements)(1)*(aElementList.NumElements)(2)*2)*pow(8,level);
        }
    }
    Mat<uint> tRemoveElements(tNRemoveElements,1,0);
    for(uint level=0; level<tNELevel; level++)
    {
        if( aElementList.Dim == 2)
        {
            for(uint i=0; i<aElementList.Polynomial*pow(2,level); i++)
            {
                for(uint j=0; j<(aElementList.NumElements)(1)*pow(2,level); j++)
                {
                    Position(0) = i; Position(1) = j;
                    tRemoveElements(tVar) = give_element_of_position(level,aElementList.Dim,aElementList.NumElements,Position);
                    tVar++;
                }
            }
            for(uint i=((aElementList.NumElements)(0)*pow(2,level)-(aElementList.Polynomial)*pow(2,level)); i<(aElementList.NumElements)(0)*pow(2,level); i++)
            {
                for(uint j=0; j<(aElementList.NumElements)(1)*pow(2,level); j++)
                {
                    Position(0) = i; Position(1) = j;
                    tRemoveElements(tVar) = give_element_of_position(level,aElementList.Dim,aElementList.NumElements,Position);
                    tVar++;
                }
            }
            for(uint j=0; j<(aElementList.Polynomial)*pow(2,level); j++)
            {
                for(uint i=0; i<(aElementList.NumElements)(0)*pow(2,level); i++)
                {
                    Position(0) = i; Position(1) = j;
                    tRemoveElements(tVar) = give_element_of_position(level,aElementList.Dim,aElementList.NumElements,Position);
                    tVar++;
                }
            }
            for(uint j=((aElementList.NumElements)(1)*pow(2,level)-(aElementList.Polynomial)*pow(2,level)); j<(aElementList.NumElements)(1)*pow(2,level); j++)
            {
                for(uint i=0; i<(aElementList.NumElements)(0)*pow(2,level); i++)
                {
                    Position(0) = i; Position(1) = j;
                    tRemoveElements(tVar) = give_element_of_position(level,aElementList.Dim,aElementList.NumElements,Position);
                    tVar++;
                }
            }
        }
        else if( aElementList.Dim == 3)
        {
            for(uint i=0; i<(aElementList.Polynomial)*pow(2,level); i++)
            {
                for(uint j=0; j<(aElementList.NumElements)(1)*pow(2,level); j++)
                {
                    for(uint k=0; k<(aElementList.NumElements)(2)*pow(2,level); k++)
                    {
                        Position(0) = i; Position(1) = j; Position(2) = k;
                        tRemoveElements(tVar) = give_element_of_position(level,aElementList.Dim,aElementList.NumElements,Position);
                        tVar++;
                    }
                }
            }
            for(uint i=((aElementList.NumElements)(0)*pow(2,level)-(aElementList.Polynomial)*pow(2,level)); i<(aElementList.NumElements)(0)*pow(2,level); i++)
            {
                for(uint j=0; j<(aElementList.NumElements)(1)*pow(2,level); j++)
                {
                    for(uint k=0; k<(aElementList.NumElements)(2)*pow(2,level); k++)
                    {
                        Position(0) = i; Position(1) = j; Position(2) = k;
                        tRemoveElements(tVar) = give_element_of_position(level,aElementList.Dim,aElementList.NumElements,Position);
                        tVar++;
                    }
                }
            }
            for(uint j=0; j<(aElementList.Polynomial)*pow(2,level); j++)
            {
                for(uint i=0; i<(aElementList.NumElements)(0)*pow(2,level); i++)
                {
                    for(uint k=0; k<(aElementList.NumElements)(2)*pow(2,level); k++)
                    {
                        Position(0) = i; Position(1) = j; Position(2) = k;
                        tRemoveElements(tVar) = give_element_of_position(level,aElementList.Dim,(aElementList.NumElements),Position);
                        tVar++;
                    }
                }
            }
            for(uint j=((aElementList.NumElements)(1)*pow(2,level)-(aElementList.Polynomial)*pow(2,level)); j<(aElementList.NumElements)(1)*pow(2,level); j++)
            {
                for(uint i=0; i<(aElementList.NumElements)(0)*pow(2,level); i++)
                {
                    for(uint k=0; k<(aElementList.NumElements)(2)*pow(2,level); k++)
                    {
                        Position(0) = i; Position(1) = j; Position(2) = k;
                        tRemoveElements(tVar) = give_element_of_position(level,aElementList.Dim,(aElementList.NumElements),Position);
                        tVar++;
                    }
                }
            }
            for(uint k=0; k<(aElementList.Polynomial)*pow(2,level); k++)
            {
                for(uint i=0; i<(aElementList.NumElements)(0)*pow(2,level); i++)
                {
                    for(uint j=0; j<(aElementList.NumElements)(1)*pow(2,level); j++)
                    {
                        Position(0) = i; Position(1) = j; Position(2) = k;
                        tRemoveElements(tVar) = give_element_of_position(level,aElementList.Dim,aElementList.NumElements,Position);
                        tVar++;
                    }
                }
            }
            for(uint k=((aElementList.NumElements)(2)*pow(2,level)-(aElementList.Polynomial)*pow(2,level)); k<(aElementList.NumElements)(2)*pow(2,level); k++)
            {
                for(uint i=0; i<(aElementList.NumElements)(0)*pow(2,level); i++)
                {
                    for(uint j=0; j<(aElementList.NumElements)(1)*pow(2,level); j++)
                    {
                        Position(0) = i; Position(1) = j; Position(2) = k;
                        tRemoveElements(tVar) = give_element_of_position(level,aElementList.Dim,aElementList.NumElements,Position);
                        tVar++;
                    }
                }
            }
        }
    }
    tRemoveElements = unique(tRemoveElements);
    return tRemoveElements;
}

Cell<Mat<real>>
moris::Hierarchical_Mesh::give_Tmatrices_of_childs(uint & aPolynomialDegree,
        uint & aDim)
{
    uint tVar1=-1,tVar2=0; // Temporary variables for a loop
    Cell<Mat<real>> tT_child(pow(2,aDim));// Cell to store all the T-matrices of the childs
    Mat<real> tChild_left(aPolynomialDegree+1,aPolynomialDegree+1,0);// Left T-matrix of the parent
    Mat<real> tChild_right(aPolynomialDegree+1,aPolynomialDegree+1,0);// Left T-matrix of the parent

    if(aPolynomialDegree==1)
    {
        tChild_left(0,0) = 1; tChild_left(0,1) = 0.5; tChild_left(1,0) = 0; tChild_left(1,1) = 0.5;
        tChild_right(0,0) = 0.5; tChild_right(0,1) = 0; tChild_right(1,0) = 0.5; tChild_right(1,1) = 1;
    }
    else if(aPolynomialDegree==2)
    {
        tChild_left(0,0) = 0.75; tChild_left(0,1) = 0.25; tChild_left(0,2) = 0;
        tChild_left(1,0) = 0.25; tChild_left(1,1) = 0.75; tChild_left(1,2) = 0.75;
        tChild_left(2,0) = 0;    tChild_left(2,1) = 0;   tChild_left(2,2) = 0.25;
        tChild_right(0,0) = 0.25; tChild_right(0,1) = 0;   tChild_right(0,2) = 0;
        tChild_right(1,0) = 0.27; tChild_right(1,1) = 0.75; tChild_right(1,2) = 0.25;
        tChild_right(2,0) = 0;    tChild_right(2,1) = 0.25; tChild_right(2,2) = 0.75;
    }
    else if(aPolynomialDegree==3)
    {
        tChild_left(0,0) = 0.5; tChild_left(0,1) = 0.125; tChild_left(0,2) = 0;   tChild_left(0,3) = 0;
        tChild_left(1,0) = 0.5; tChild_left(1,1) = 0.75; tChild_left(1,2) = 0.5; tChild_left(1,3) = 0.125;
        tChild_left(2,0) = 0;   tChild_left(2,1) = 0.125; tChild_left(2,2) = 0.5; tChild_left(2,3) = 0.75;
        tChild_left(3,0) = 0;   tChild_left(3,1) = 0;    tChild_left(3,2) = 0;   tChild_left(3,3) = 0.125;
        tChild_right(0,0) = 0.125; tChild_right(0,1) = 0;  tChild_right(0,2) = 0;     tChild_left(0,3) = 0;
        tChild_right(1,0) = 0.75;  tChild_right(1,1) = 0.5; tChild_right(1,2) = 0.125; tChild_right(1,3) = 0;
        tChild_right(2,0) = 0.125; tChild_right(2,1) = 0.5; tChild_right(2,2) = 0.75;  tChild_right(2,3) = 0.5;
        tChild_right(3,0) = 0;     tChild_right(3,1) = 0;  tChild_right(3,2) = 0.125; tChild_right(3,3) = 0.5;
    }
    else
    {
        MORIS_LOG_ERROR << "T-matrix for this polynomial degree is not available";
    }

    // Sace the left and right child in the cell tT_child
    if(aDim==1)
    {
        tT_child(0) = tChild_left;
        tT_child(1) = tChild_right;
    }

    // T-matrix of the four childs are calculated automatically for each polynomial degree
    if(aDim==2)
    {
        if(aPolynomialDegree==1)
        {
            Mat<real> tChild1 ={{     1.0000,    0.5000,    0.5000,   0.2500},
                    {0,    0.5000,         0,    0.2500},
                    {0,         0,    0.5000,    0.2500},
                    {0,         0,         0,    0.2500}};
            Mat<real> tChild2 ={{    0.5000,         0,    0.2500,         0},
                    {0.5000,    1.0000,    0.2500,    0.5000},
                    {    0 ,        0,    0.2500,         0},
                    {    0 ,        0,    0.2500,    0.5000}};
            Mat<real> tChild3 ={{    0.5000,    0.2500,         0,         0},
                    {0,    0.2500,         0,         0},
                    {0.5000,    0.2500,    1.0000,    0.5000},
                    {    0,    0.2500,         0,    0.5000}};
            Mat<real> tChild4 ={{     0.2500,         0 ,        0,         0},
                    {0.2500,    0.5000,         0,         0},
                    {0.2500,         0,    0.5000,         0},
                    {0.2500,    0.5000,    0.5000,    1.0000}};
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
        }
        else if(aPolynomialDegree==2)
        {
            Mat<real> tChild1 ={{0.562500,0.187500,0.000000,0.187500,0.062500,0.000000,0.000000,0.000000,0.000000},
                    {0.187500,0.562500,0.562500,0.062500,0.187500,0.187500,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.187500,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000},
                    {0.187500,0.062500,0.000000,0.562500,0.187500,0.000000,0.562500,0.187500,0.000000},
                    {0.062500,0.187500,0.187500,0.187500,0.562500,0.562500,0.187500,0.562500,0.562500},
                    {0.000000,0.000000,0.062500,0.000000,0.000000,0.187500,0.000000,0.000000,0.187500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.187500,0.062500,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.062500,0.187500,0.187500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.062500}};
            Mat<real> tChild2 ={{0.187500,0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.562500,0.562500,0.187500,0.187500,0.187500,0.062500,0.000000,0.000000,0.000000},
                    {0.000000,0.187500,0.562500,0.000000,0.062500,0.187500,0.000000,0.000000,0.000000},
                    {0.062500,0.000000,0.000000,0.187500,0.000000,0.000000,0.187500,0.000000,0.000000},
                    {0.187500,0.187500,0.062500,0.562500,0.562500,0.187500,0.562500,0.562500,0.187500},
                    {0.000000,0.062500,0.187500,0.000000,0.187500,0.562500,0.000000,0.187500,0.562500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.187500,0.187500,0.062500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.062500,0.187500}};
            Mat<real> tChild3 ={{0.187500,0.062500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.062500,0.187500,0.187500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.062500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.562500,0.187500,0.000000,0.562500,0.187500,0.000000,0.187500,0.062500,0.000000},
                    {0.187500,0.562500,0.562500,0.187500,0.562500,0.562500,0.062500,0.187500,0.187500},
                    {0.000000,0.000000,0.187500,0.000000,0.000000,0.187500,0.000000,0.000000,0.062500},
                    {0.000000,0.000000,0.000000,0.187500,0.062500,0.000000,0.562500,0.187500,0.000000},
                    {0.000000,0.000000,0.000000,0.062500,0.187500,0.187500,0.187500,0.562500,0.562500},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.187500}};
            Mat<real> tChild4 ={{0.062500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.187500,0.187500,0.062500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.062500,0.187500,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.187500,0.000000,0.000000,0.187500,0.000000,0.000000,0.062500,0.000000,0.000000},
                    {0.562500,0.562500,0.187500,0.562500,0.562500,0.187500,0.187500,0.187500,0.062500},
                    {0.000000,0.187500,0.562500,0.000000,0.187500,0.562500,0.000000,0.062500,0.187500},
                    {0.000000,0.000000,0.000000,0.062500,0.000000,0.000000,0.187500,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.187500,0.187500,0.062500,0.562500,0.562500,0.187500},
                    {0.000000,0.000000,0.000000,0.000000,0.062500,0.187500,0.000000,0.187500,0.562500}};
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
        }
        else if(aPolynomialDegree==3)
        {
            Mat<real> tChild1 = {{0.250000,0.062500,0.000000,0.000000,0.062500,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
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
            Mat<real> tChild2 = {{0.062500,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
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
            Mat<real> tChild3 ={{0.062500,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
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
            Mat<real> tChild4 = {{0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
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
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
        }
        else if(aPolynomialDegree>=4)
        {
            // Initialize the Childs
            Mat<real> tChild1(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild2(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild3(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild4(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);

            for(moris::uint k = 0; k<aPolynomialDegree+1; k++)
            {
                for(moris::uint l = 0; l<aPolynomialDegree+1; l++)
                {
                    tVar1++;
                    tVar2=0;
                    for(moris::uint i = 0; i<aPolynomialDegree+1; i++)
                    {
                        for(moris::uint j = 0; j<aPolynomialDegree+1; j++)
                        {
                            // Calculate the Childs in 2D with the help of the Childs of 1D
                            tChild1(tVar1,tVar2) = tChild_left(k,i)*tChild_left(l,j);
                            tChild2(tVar1,tVar2) = tChild_left(k,i)*tChild_right(l,j);
                            tChild3(tVar1,tVar2) = tChild_right(k,i)*tChild_left(l,j);
                            tChild4(tVar1,tVar2) = tChild_right(k,i)*tChild_right(l,j);
                            tVar2++;
                        }
                    }
                }
            }
            //Save the four Childs in the Cell
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
        }
    }
    // T-matrix of the eight childs are calculated automatically for each polynomial degree
    if(aDim==3)
    {
        if(aPolynomialDegree==1)
        {
            Mat<real> tChild1 = {{1.000000,0.500000,0.500000,0.250000,0.500000,0.250000,0.250000,0.125000},
                    {0.000000,0.500000,0.000000,0.250000,0.000000,0.250000,0.000000,0.125000},
                    {0.000000,0.000000,0.500000,0.250000,0.000000,0.000000,0.250000,0.125000},
                    {0.000000,0.000000,0.000000,0.250000,0.000000,0.000000,0.000000,0.125000},
                    {0.000000,0.000000,0.000000,0.000000,0.500000,0.250000,0.250000,0.125000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.250000,0.000000,0.125000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.250000,0.125000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.125000}};
            Mat<real> tChild2 =  {{0.500000,0.000000,0.250000,0.000000,0.250000,0.000000,0.125000,0.000000},
                    {0.500000,1.000000,0.250000,0.500000,0.250000,0.500000,0.125000,0.250000},
                    {0.000000,0.000000,0.250000,0.000000,0.000000,0.000000,0.125000,0.000000},
                    {0.000000,0.000000,0.250000,0.500000,0.000000,0.000000,0.125000,0.250000},
                    {0.000000,0.000000,0.000000,0.000000,0.250000,0.000000,0.125000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.250000,0.500000,0.125000,0.250000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.125000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.125000,0.250000}};
            Mat<real> tChild3 = {{0.500000,0.250000,0.000000,0.000000,0.250000,0.125000,0.000000,0.000000},
                    {0.000000,0.250000,0.000000,0.000000,0.000000,0.125000,0.000000,0.000000},
                    {0.500000,0.250000,1.000000,0.500000,0.250000,0.125000,0.500000,0.250000},
                    {0.000000,0.250000,0.000000,0.500000,0.000000,0.125000,0.000000,0.250000},
                    {0.000000,0.000000,0.000000,0.000000,0.250000,0.125000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.125000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.250000,0.125000,0.500000,0.250000},
                    {0.000000,0.000000,0.000000,0.000000,0.000000,0.125000,0.000000,0.250000}};
            Mat<real> tChild4 = {{0.250000,0.000000,0.000000,0.000000,0.125000,0.000000,0.000000,0.000000},
                    {0.250000,0.500000,0.000000,0.000000,0.125000,0.250000,0.000000,0.000000},
                    {0.250000,0.000000,0.500000,0.000000,0.125000,0.000000,0.250000,0.000000},
                    {0.250000,0.500000,0.500000,1.000000,0.125000,0.250000,0.250000,0.500000},
                    {0.000000,0.000000,0.000000,0.000000,0.125000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.125000,0.250000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.125000,0.000000,0.250000,0.000000},
                    {0.000000,0.000000,0.000000,0.000000,0.125000,0.250000,0.250000,0.500000}};
            Mat<real> tChild5 = {{0.500000,0.250000,0.250000,0.125000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.250000,0.000000,0.125000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.250000,0.125000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.000000,0.125000,0.000000,0.000000,0.000000,0.000000},
                    {0.500000,0.250000,0.250000,0.125000,1.000000,0.500000,0.500000,0.250000},
                    {0.000000,0.250000,0.000000,0.125000,0.000000,0.500000,0.000000,0.250000},
                    {0.000000,0.000000,0.250000,0.125000,0.000000,0.000000,0.500000,0.250000},
                    {0.000000,0.000000,0.000000,0.125000,0.000000,0.000000,0.000000,0.250000,}};
            Mat<real> tChild6 = {{0.250000,0.000000,0.125000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.250000,0.500000,0.125000,0.250000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.125000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.000000,0.125000,0.250000,0.000000,0.000000,0.000000,0.000000},
                    {0.250000,0.000000,0.125000,0.000000,0.500000,0.000000,0.250000,0.000000},
                    {0.250000,0.500000,0.125000,0.250000,0.500000,1.000000,0.250000,0.500000},
                    {0.000000,0.000000,0.125000,0.000000,0.000000,0.000000,0.250000,0.000000},
                    {0.000000,0.000000,0.125000,0.250000,0.000000,0.000000,0.250000,0.500000}};
            Mat<real> tChild7 = {{0.250000,0.125000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.125000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.250000,0.125000,0.500000,0.250000,0.000000,0.000000,0.000000,0.000000},
                    {0.000000,0.125000,0.000000,0.250000,0.000000,0.000000,0.000000,0.000000},
                    {0.250000,0.125000,0.000000,0.000000,0.500000,0.250000,0.000000,0.000000},
                    {0.000000,0.125000,0.000000,0.000000,0.000000,0.250000,0.000000,0.000000},
                    {0.250000,0.125000,0.500000,0.250000,0.500000,0.250000,1.000000,0.500000},
                    {0.000000,0.125000,0.000000,0.250000,0.000000,0.250000,0.000000,0.500000,}};
            Mat<real> tChild8 = {{0.125000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.125000,0.250000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.125000,0.000000,0.250000,0.000000,0.000000,0.000000,0.000000,0.000000},
                    {0.125000,0.250000,0.250000,0.500000,0.000000,0.000000,0.000000,0.000000},
                    {0.125000,0.000000,0.000000,0.000000,0.250000,0.000000,0.000000,0.000000},
                    {0.125000,0.250000,0.000000,0.000000,0.250000,0.500000,0.000000,0.000000},
                    {0.125000,0.000000,0.250000,0.000000,0.250000,0.000000,0.500000,0.000000},
                    {0.125000,0.250000,0.250000,0.500000,0.250000,0.500000,0.500000,1.000000,}};
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
            tT_child(4) = tChild5;
            tT_child(5) = tChild6;
            tT_child(6) = tChild7;
            tT_child(7) = tChild8;
        }
        else if(aPolynomialDegree==2)
        {
            Mat<real> tChild1 = {{0.421875,0.140625,0.000000,0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.140625,0.046875,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
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
            Mat<real> tChild2 = {{0.140625,0.000000,0.000000,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
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
            Mat<real> tChild3 = {{0.140625,0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
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
            Mat<real> tChild4 = {{0.046875,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
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
            Mat<real> tChild5 = {{0.140625,0.046875,0.000000,0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
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
            Mat<real> tChild6 = {{0.046875,0.000000,0.000000,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
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
            Mat<real> tChild7 = {{0.046875,0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
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
            Mat<real> tChild8 = {{0.015625,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000},
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
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
            tT_child(4) = tChild5;
            tT_child(5) = tChild6;
            tT_child(6) = tChild7;
            tT_child(7) = tChild8;
        }
        else if(aPolynomialDegree>=3)
        {
            // Initialize the Childs
            Mat<real> tChild1(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild2(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild3(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild4(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild5(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild6(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild7(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild8(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);

            for(uint k = 0; k<aPolynomialDegree+1; k++)
            {
                for(uint l = 0; l<aPolynomialDegree+1; l++)
                {
                    for(uint m = 0; m<aPolynomialDegree+1; m++)
                    {
                        tVar1++;
                        tVar2=0;
                        for(uint i = 0; i<aPolynomialDegree+1; i++)
                        {
                            for(uint j = 0; j<aPolynomialDegree+1; j++)
                            {
                                for(uint n = 0; n<aPolynomialDegree+1; n++)
                                {
                                    // Calculate the Childs in 3D with the help of the Childs of 1D
                                    tChild1(tVar1,tVar2) = tChild_left(m,n)*tChild_left(l,j)*tChild_left(k,i);
                                    tChild2(tVar1,tVar2) = tChild_right(m,n)*tChild_left(l,j)*tChild_left(k,i);
                                    tChild3(tVar1,tVar2) = tChild_left(m,n)*tChild_right(l,j)*tChild_left(k,i);
                                    tChild4(tVar1,tVar2) = tChild_right(m,n)*tChild_right(l,j)*tChild_left(k,i);
                                    tChild5(tVar1,tVar2) = tChild_left(m,n)*tChild_left(l,j)*tChild_right(k,i);
                                    tChild6(tVar1,tVar2) = tChild_right(m,n)*tChild_left(l,j)*tChild_right(k,i);
                                    tChild7(tVar1,tVar2) = tChild_left(m,n)*tChild_right(l,j)*tChild_right(k,i);
                                    tChild8(tVar1,tVar2) = tChild_right(m,n)*tChild_right(l,j)*tChild_right(k,i);
                                    tVar2++;
                                }
                            }
                        }
                    }
                }
            }
            //Save the four Childs in the Cell
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
            tT_child(4) = tChild5;
            tT_child(5) = tChild6;
            tT_child(6) = tChild7;
            tT_child(7) = tChild8;
        }
    }
    return tT_child;
}

Cell<Mat<real>>
moris::Hierarchical_Mesh::give_Tmatrices_of_childs_for_design(uint & aPolynomialDegree,
        uint & aDim)
{
    uint tVar1=-1,tVar2=0; // Temporary variables for a loop
    Cell<Mat<real>> tT_child(pow(2,aDim));// Cell to store all the T-matrices of the childs
    Mat<real> tChild_left(aPolynomialDegree+1,aPolynomialDegree+1,0);// Left T-matrix of the parent
    Mat<real> tChild_right(aPolynomialDegree+1,aPolynomialDegree+1,0);// Left T-matrix of the parent

    if(aPolynomialDegree==1)
    {
        tChild_left(0,0) = 1; tChild_left(0,1) = 0.5; tChild_left(1,0) = 0; tChild_left(1,1) = 0.5;
        tChild_right(0,0) = 0.5; tChild_right(0,1) = 0; tChild_right(1,0) = 0.5; tChild_right(1,1) = 1;
    }
    else if(aPolynomialDegree==2)
    {
        tChild_left(0,0) = 0.75; tChild_left(0,1) = 0.25; tChild_left(0,2) = 0;
        tChild_left(1,0) = 0.25; tChild_left(1,1) = 0.75; tChild_left(1,2) = 0.75;
        tChild_left(2,0) = 0;    tChild_left(2,1) = 0;   tChild_left(2,2) = 0.25;
        tChild_right(0,0) = 0.25; tChild_right(0,1) = 0;   tChild_right(0,2) = 0;
        tChild_right(1,0) = 0.27; tChild_right(1,1) = 0.75; tChild_right(1,2) = 0.25;
        tChild_right(2,0) = 0;    tChild_right(2,1) = 0.25; tChild_right(2,2) = 0.75;
    }
    else if(aPolynomialDegree==3)
    {
        tChild_left(0,0) = 0.5; tChild_left(0,1) = 0.125; tChild_left(0,2) = 0;   tChild_left(0,3) = 0;
        tChild_left(1,0) = 0.5; tChild_left(1,1) = 0.75; tChild_left(1,2) = 0.5; tChild_left(1,3) = 0.125;
        tChild_left(2,0) = 0;   tChild_left(2,1) = 0.125; tChild_left(2,2) = 0.5; tChild_left(2,3) = 0.75;
        tChild_left(3,0) = 0;   tChild_left(3,1) = 0;    tChild_left(3,2) = 0;   tChild_left(3,3) = 0.125;
        tChild_right(0,0) = 0.125; tChild_right(0,1) = 0;  tChild_right(0,2) = 0;     tChild_left(0,3) = 0;
        tChild_right(1,0) = 0.75;  tChild_right(1,1) = 0.5; tChild_right(1,2) = 0.125; tChild_right(1,3) = 0;
        tChild_right(2,0) = 0.125; tChild_right(2,1) = 0.5; tChild_right(2,2) = 0.75;  tChild_right(2,3) = 0.5;
        tChild_right(3,0) = 0;     tChild_right(3,1) = 0;  tChild_right(3,2) = 0.125; tChild_right(3,3) = 0.5;
    }
    else
    {
        MORIS_LOG_ERROR << "T-matrix for this polynomial degree is not available";
    }

    // Sace the left and right child in the cell tT_child
    if(aDim==1)
    {
        tT_child(0) = tChild_left;
        tT_child(1) = tChild_right;
    }

    // T-matrix of the four childs are calculated automatically for each polynomial degree
    if(aDim==2)
    {
        if(aPolynomialDegree==1)
        {
            Mat<real> tChild1 ={{    1.0000,    0.5000,    0.2500,     0.5000},
                    {0,    0.5000,    0.2500,         0},
                    {0,         0,    0.2500,         0},
                    {0,         0,    0.2500,    0.5000}};
            Mat<real> tChild2 ={{    0.5000,         0 ,        0,    0.2500},
                    {0.5000,    1.0000,    0.5000,    0.2500},
                    {0,         0,    0.5000,    0.2500},
                    {0,         0,         0,    0.2500}};
            Mat<real> tChild3 ={{    0.5000,    0.2500,         0,         0},
                    {0,    0.2500,         0,         0},
                    {0,    0.2500,    0.5000,         0},
                    {0.5000,    0.2500,    0.5000,    1.0000}};
            Mat<real> tChild4 ={{    0.2500,         0,         0,         0},
                    {0.2500,    0.5000,         0,         0},
                    {0.2500,    0.5000,    1.0000,    0.5000},
                    {0.2500,         0,         0,    0.5000}};
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
        }
        else if(aPolynomialDegree==2)
        {
            Mat<real> tChild1 ={{0.5625, 0.1875, 0, 0.1875, 0.0625, 0, 0, 0, 0},{0.1875, 0.5625, 0.5625, 0.0625, 0.1875, 0.1875, 0, 0, 0},{0, 0, 0.1875, 0, 0, 0.0625, 0, 0, 0},
                    {0.1875, 0.0625, 0, 0.5625, 0.1875, 0, 0.5625, 0.1875, 0},{0.0625, 0.1875, 0.1875, 0.1875, 0.5625, 0.5625, 0.1875, 0.5625, 0.5625},
                    {0, 0, 0.0625, 0, 0, 0.1875, 0, 0, 0.1875},{0, 0, 0, 0, 0, 0, 0.1875, 0.0625, 0},{0, 0, 0, 0, 0, 0, 0.0625, 0.1875, 0.1875},{0, 0, 0, 0, 0, 0, 0, 0, 0.0625}};
            Mat<real> tChild2 ={{0.1875, 0, 0, 0.0625, 0, 0, 0, 0, 0},{0.5625, 0.5625, 0.1875, 0.1875, 0.1875, 0.0625, 0, 0, 0},{0, 0.1875, 0.5625, 0, 0.0625, 0.1875, 0, 0, 0},
                    {0.0625, 0, 0, 0.1875, 0, 0, 0.1875, 0, 0},{0.1875, 0.1875, 0.0625, 0.5625, 0.5625, 0.1875, 0.5625, 0.5625, 0.1875},{0, 0.0625, 0.1875, 0, 0.1875, 0.5625, 0, 0.1875, 0.5625},
                    {0, 0, 0, 0, 0, 0, 0.0625, 0, 0},{0, 0, 0, 0, 0, 0, 0.1875, 0.1875, 0.0625},{0, 0, 0, 0, 0, 0, 0, 0.0625, 0.1875}};
            Mat<real> tChild3 ={{0.1875, 0.0625, 0, 0, 0, 0, 0, 0, 0},{0.0625, 0.1875, 0.1875, 0, 0, 0, 0, 0, 0},{0, 0, 0.0625, 0, 0, 0, 0, 0, 0},
                    {0.5625, 0.1875, 0, 0.5625, 0.1875, 0, 0.1875, 0.0625, 0},{0.1875, 0.5625, 0.5625, 0.1875, 0.5625, 0.5625, 0.0625, 0.1875, 0.1875},{0, 0, 0.1875, 0, 0, 0.1875, 0, 0, 0.0625},
                    {0, 0, 0, 0.1875, 0.0625, 0, 0.5625, 0.1875, 0},{0, 0, 0, 0.0625, 0.1875, 0.1875, 0.1875, 0.5625, 0.5625},{0, 0, 0, 0, 0, 0.0625, 0, 0, 0.1875}};
            Mat<real> tChild4 ={{0.0625, 0, 0, 0, 0, 0, 0, 0, 0},{0.1875, 0.1875, 0.0625, 0, 0, 0, 0, 0, 0},{0, 0.0625, 0.1875, 0, 0, 0, 0, 0, 0},
                    {0.1875, 0, 0, 0.1875, 0, 0, 0.0625, 0, 0},{0.5625, 0.5625, 0.1875, 0.5625, 0.5625, 0.1875, 0.1875, 0.1875, 0.0625},{0, 0.1875, 0.5625, 0, 0.1875, 0.5625, 0, 0.0625, 0.1875},
                    {0, 0, 0, 0.0625, 0, 0, 0.1875, 0, 0},{0, 0, 0, 0.1875, 0.1875, 0.0625, 0.5625, 0.5625, 0.1875},{0, 0, 0, 0, 0.0625, 0.1875, 0, 0.1875, 0.5625}};
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
        }
        else if(aPolynomialDegree==3)
        {
            Mat<real> tChild1 ={{0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625}};
            Mat<real> tChild2 ={{0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500}};
            Mat<real> tChild3 ={{0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500}};
            Mat<real> tChild4 ={{0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.093750, 0.015625, 0.000000, 0.093750, 0.562500, 0.093750, 0.000000, 0.015625, 0.093750, 0.015625},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.062500, 0.000000, 0.000000, 0.375000, 0.375000, 0.000000, 0.000000, 0.062500, 0.062500},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.062500, 0.375000, 0.062500, 0.000000, 0.062500, 0.375000, 0.062500},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.250000, 0.250000, 0.000000, 0.000000, 0.250000, 0.250000}};
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
        }
        else if(aPolynomialDegree>=4)
        {
            // Initialize the Childs
            Mat<real> tChild1(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild2(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild3(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild4(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);

            for(moris::uint k = 0; k<aPolynomialDegree+1; k++)
            {
                for(moris::uint l = 0; l<aPolynomialDegree+1; l++)
                {
                    tVar1++;
                    tVar2=0;
                    for(moris::uint i = 0; i<aPolynomialDegree+1; i++)
                    {
                        for(moris::uint j = 0; j<aPolynomialDegree+1; j++)
                        {
                            // Calculate the Childs in 2D with the help of the Childs of 1D
                            tChild1(tVar1,tVar2) = tChild_left(k,i)*tChild_left(l,j);
                            tChild2(tVar1,tVar2) = tChild_left(k,i)*tChild_right(l,j);
                            tChild3(tVar1,tVar2) = tChild_right(k,i)*tChild_left(l,j);
                            tChild4(tVar1,tVar2) = tChild_right(k,i)*tChild_right(l,j);
                            tVar2++;
                        }
                    }
                }
            }
            //Save the four Childs in the Cell
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
        }
    }
    // T-matrix of the eight childs are calculated automatically for each polynomial degree
    if(aDim==3)
    {
        if(aPolynomialDegree==1)
        {
            Mat<real> tChild1 = {{             1.0000,    0.5000,    0.5000,    0.2500,    0.5000,    0.2500,    0.2500,    0.1250},
                    {         0,    0.5000,         0,    0.2500,         0,    0.2500,         0,    0.1250},
                    {         0,         0,    0.5000,    0.2500,         0,         0,    0.2500,    0.1250},
                    {         0,         0,         0,    0.2500,         0,         0,         0,    0.1250},
                    {         0,         0,         0,         0,    0.5000,    0.2500,    0.2500,    0.1250},
                    {         0,         0,         0,         0,         0,    0.2500,         0,    0.1250},
                    {         0,         0,         0,         0,         0,         0,    0.2500,    0.1250},
                    {         0,         0,         0,         0,         0,         0,         0,    0.1250}};
            Mat<real> tChild2 =  {{        0.5000,         0,    0.2500,         0,    0.2500,         0,    0.1250,         0},
                    {    0.5000,    1.0000,    0.2500,    0.5000,    0.2500,    0.5000,    0.1250,    0.2500},
                    {         0,         0,    0.2500,         0,         0,         0,    0.1250,         0},
                    {         0,         0,    0.2500,    0.5000,         0,         0,    0.1250,    0.2500},
                    {         0,         0,         0,         0,    0.2500,         0,    0.1250,         0},
                    {         0,         0,         0,         0,    0.2500,    0.5000,    0.1250,    0.2500},
                    {         0,         0,         0,         0,         0,         0,    0.1250,         0},
                    {         0,         0,         0,         0,         0,         0,    0.1250,    0.2500}};
            Mat<real> tChild3 = {{       0.5000,    0.2500,         0,         0,    0.2500,    0.1250,         0,         0},
                    {         0,    0.2500,         0,         0,         0,    0.1250,         0,         0},
                    {    0.5000,    0.2500,    1.0000,    0.5000,    0.2500,    0.1250,    0.5000,    0.2500},
                    {         0,    0.2500,         0,    0.5000,         0,    0.1250,         0,    0.2500},
                    {         0,         0,         0,         0,    0.2500,    0.1250,         0,         0},
                    {         0,         0,         0,         0,         0,    0.1250,         0,         0},
                    {         0,         0,         0,         0,    0.2500,    0.1250,    0.5000,    0.2500},
                    {         0,         0,         0,         0,         0,    0.1250,         0,    0.2500}};
            Mat<real> tChild4 = {{    0.2500,         0,         0,         0,    0.1250,         0,         0,         0},
                    {    0.2500,    0.5000,         0,         0,    0.1250,    0.2500,         0,         0},
                    {    0.2500,         0,    0.5000,         0,    0.1250,         0,    0.2500,         0},
                    {    0.2500,    0.5000,    0.5000,    1.0000,    0.1250,    0.2500,    0.2500,    0.5000},
                    {         0,         0,         0,         0,    0.1250,         0,         0,         0},
                    {         0,         0,         0,         0,    0.1250,    0.2500,         0,         0},
                    {         0,         0,         0,         0,    0.1250,         0,    0.2500,         0},
                    {         0,         0,         0,         0,    0.1250,    0.2500,    0.2500,    0.5000}};
            Mat<real> tChild5 = {{    0.5000,    0.2500,    0.2500,    0.1250,         0,         0,         0,         0},
                    {         0,    0.2500,         0,    0.1250,         0,         0,         0,         0                  },
                    {         0,         0,    0.2500,    0.1250,         0,         0,         0,         0                  },
                    {         0,         0,         0,    0.1250,         0,         0,         0,         0                  },
                    {    0.5000,    0.2500,    0.2500,    0.1250,    1.0000,    0.5000,    0.5000,    0.2500                  },
                    {         0,    0.2500,         0,    0.1250,         0,    0.5000,         0,    0.2500                  },
                    {         0,         0,    0.2500,    0.1250,         0,         0,    0.5000,    0.2500                  },
                    {         0,         0,         0,    0.1250,         0,         0,         0,    0.2500}};
            Mat<real> tChild6 = {{    0.2500,         0,    0.1250,         0,         0,         0,         0,         0},
                    {    0.2500,    0.5000,    0.1250,    0.2500,         0,         0,         0,         0               },
                    {         0,         0,    0.1250,         0,         0,         0,         0,         0               },
                    {         0,         0,    0.1250,    0.2500,         0,         0,         0,         0               },
                    {    0.2500,         0,    0.1250,         0,    0.5000,         0,    0.2500,         0               },
                    {    0.2500,    0.5000,    0.1250,    0.2500,    0.5000,    1.0000,    0.2500,    0.5000               },
                    {         0,         0,    0.1250,         0,         0,         0,    0.2500,         0               },
                    {         0,         0,    0.1250,    0.2500,         0,         0,    0.2500,    0.5000}};
            Mat<real> tChild7 = {{    0.2500,    0.1250,         0,         0,         0,         0,         0,         0},
                    {         0,    0.1250,         0,         0,         0,         0,         0,         0                   },
                    {    0.2500,    0.1250,    0.5000,    0.2500,         0,         0,         0,         0                   },
                    {         0,    0.1250,         0,    0.2500,         0,         0,         0,         0                   },
                    {    0.2500,    0.1250,         0,         0,    0.5000,    0.2500,         0,         0                   },
                    {         0,    0.1250,         0,         0,         0,    0.2500,         0,         0                   },
                    {    0.2500,    0.1250,    0.5000,    0.2500,    0.5000,    0.2500,    1.0000,    0.5000                   },
                    {         0,    0.1250,         0,    0.2500,         0,    0.2500,         0,    0.5000}};
            Mat<real> tChild8 = {{     0.1250,         0,         0,         0,         0,         0,         0,         0},
                    {    0.1250,    0.2500,         0,         0,         0,         0,         0,         0                },
                    {    0.1250,         0,    0.2500,         0,         0,         0,         0,         0                },
                    {    0.1250,    0.2500,    0.2500,    0.5000,         0,         0,         0,         0                },
                    {    0.1250,         0,         0,         0,    0.2500,         0,         0,         0                },
                    {    0.1250,    0.2500,         0,         0,    0.2500,    0.5000,         0,         0                },
                    {    0.1250,         0,    0.2500,         0,    0.2500,         0,    0.5000,         0                },
                    {    0.1250,    0.2500,    0.2500,    0.5000,    0.2500,    0.5000,    0.5000,    1.0000}};
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
            tT_child(4) = tChild5;
            tT_child(5) = tChild6;
            tT_child(6) = tChild7;
            tT_child(7) = tChild8;
        }
        else if(aPolynomialDegree==2)
        {
            Mat<real> tChild1 = {{0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625}};
            Mat<real> tChild2 = {{0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875}};
            Mat<real> tChild3 = {{0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875}};
            Mat<real> tChild4 = {{0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625}};
            Mat<real> tChild5 = {{0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875}};
            Mat<real> tChild6 = {{0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625}};
            Mat<real> tChild7 = {{0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625}};
            Mat<real> tChild8 = {{0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000, 0.000000, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.046875, 0.015625, 0.000000, 0.000000, 0.000000, 0.000000, 0.421875, 0.140625, 0.000000, 0.140625, 0.046875},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.015625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.421875, 0.000000, 0.046875, 0.140625},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875, 0.000000},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.015625, 0.000000, 0.140625, 0.046875, 0.000000, 0.000000, 0.000000, 0.000000, 0.140625, 0.046875, 0.000000, 0.421875, 0.140625},
                    {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.015625, 0.046875, 0.000000, 0.046875, 0.140625, 0.000000, 0.000000, 0.000000, 0.000000, 0.046875, 0.140625, 0.000000, 0.140625, 0.421875}};
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
            tT_child(4) = tChild5;
            tT_child(5) = tChild6;
            tT_child(6) = tChild7;
            tT_child(7) = tChild8;
        }
        else if(aPolynomialDegree>=3)
        {
            // Initialize the Childs
            Mat<real> tChild1(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild2(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild3(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild4(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild5(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild6(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild7(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);
            Mat<real> tChild8(pow(aPolynomialDegree+1,aDim),pow(aPolynomialDegree+1,aDim),0);

            for(uint k = 0; k<aPolynomialDegree+1; k++)
            {
                for(uint l = 0; l<aPolynomialDegree+1; l++)
                {
                    for(uint m = 0; m<aPolynomialDegree+1; m++)
                    {
                        tVar1++;
                        tVar2=0;
                        for(uint i = 0; i<aPolynomialDegree+1; i++)
                        {
                            for(uint j = 0; j<aPolynomialDegree+1; j++)
                            {
                                for(uint n = 0; n<aPolynomialDegree+1; n++)
                                {
                                    // Calculate the Childs in 3D with the help of the Childs of 1D
                                    tChild1(tVar1,tVar2) = tChild_left(m,n)*tChild_left(l,j)*tChild_left(k,i);
                                    tChild2(tVar1,tVar2) = tChild_right(m,n)*tChild_left(l,j)*tChild_left(k,i);
                                    tChild3(tVar1,tVar2) = tChild_left(m,n)*tChild_right(l,j)*tChild_left(k,i);
                                    tChild4(tVar1,tVar2) = tChild_right(m,n)*tChild_right(l,j)*tChild_left(k,i);
                                    tChild5(tVar1,tVar2) = tChild_left(m,n)*tChild_left(l,j)*tChild_right(k,i);
                                    tChild6(tVar1,tVar2) = tChild_right(m,n)*tChild_left(l,j)*tChild_right(k,i);
                                    tChild7(tVar1,tVar2) = tChild_left(m,n)*tChild_right(l,j)*tChild_right(k,i);
                                    tChild8(tVar1,tVar2) = tChild_right(m,n)*tChild_right(l,j)*tChild_right(k,i);
                                    tVar2++;
                                }
                            }
                        }
                    }
                }
            }
            //Save the four Childs in the Cell
            tT_child(0) = tChild1;
            tT_child(1) = tChild2;
            tT_child(2) = tChild3;
            tT_child(3) = tChild4;
            tT_child(4) = tChild5;
            tT_child(5) = tChild6;
            tT_child(6) = tChild7;
            tT_child(7) = tChild8;
        }
    }
    return tT_child;
}
void moris::Hierarchical_Mesh::create_mesh_data(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tProcRank = par_rank();
    uint tVar = 0, tVarb = 0, tVarc = 0, tVard = 0;//Temporary variable for loop
    this->give_active_elements(aElementList);
    uint tSizeActiveElements = (aElementList.ElementListOnProc).length();
    (aElementList.Fetopo).set_size(tSizeActiveElements,pow((aElementList.Polynomial)+1,(aElementList.Dim)));
    (aElementList.BlockSetOutput).set_size(tSizeActiveElements,1);
    (aElementList.ElemLocaltoGlobal).resize(tSizeActiveElements,1);
    (aElementList.NodalLocaltoGlobal).set_size(tSizeActiveElements*pow((aElementList.Polynomial)+1,(aElementList.Dim)),1,UINT_MAX);
    (aElementList.NodalLocaltoGlobalExist).set_size(tSizeActiveElements*pow((aElementList.Polynomial)+1,(aElementList.Dim)),1,UINT_MAX);
    (aElementList.NodalLocaltoGlobalMap).set_size(tSizeActiveElements*pow((aElementList.Polynomial)+1,(aElementList.Dim)),2,UINT_MAX); // Internal map for the Dummy Points, which are created for Paraview
    Mat<uint> DummyNodalLocaltoGlobal(tSizeActiveElements*pow((aElementList.Polynomial)+1,(aElementList.Dim)),1,UINT_MAX);
    Mat<uint> tBasis_of_element(pow((aElementList.Polynomial)+1,(aElementList.Dim)),1);
    Mat<uint> tDesignBSpline_of_element(pow((aElementList.PolynomialDesign)+1,(aElementList.Dim)),1);
    Mat<uint> tOrder(pow((aElementList.Polynomial)+1,(aElementList.Dim)),1); // Change ordering to the classical order of FE connectivity
    Mat<real> tCoordinates(1,aElementList.Dim);
    Cell<Mat<real>> tTmatrixNodalBased(tSizeActiveElements*pow((aElementList.Polynomial)+1,(aElementList.Dim)));
    Cell<Mat<uint>> tIdFieldNodalBased(tSizeActiveElements*pow((aElementList.Polynomial)+1,(aElementList.Dim)));
    Cell<Mat<real>> tTmatrixNodalBasedDesign(tSizeActiveElements*pow((aElementList.PolynomialDesign)+1,(aElementList.Dim)));
    Cell<Mat<uint>> tIdFieldNodalBasedDesign(tSizeActiveElements*pow((aElementList.PolynomialDesign)+1,(aElementList.Dim)));
    uint tParent;
    uint tDummyNum = aBasisList.NumberBasis*(1+(uint)tProcRank);
    uint tPos;
    real tWeight;
    //    uint tElementLevel;
    uint tMaxIdField = 0;
    uint tMaxIdFieldDesign = 0;
    (aElementList.SparseSize) = 0;

    if( (aBasisList.DesignBSplineActive).size() < (aBasisList.BasisActive).size() )
        (aBasisList.DesignBSplineActive).resize((aBasisList.BasisActive).size());

    if( aElementList.Dim == 2)
    {
        tOrder(0) = 0;  tOrder(1) = 1;  tOrder(2) = 3; tOrder(3) = 2;
    }
    else if( aElementList.Dim == 3)
    {
        tOrder(0) = 0;  tOrder(1) = 1;  tOrder(2) = 3; tOrder(3) = 2; tOrder(4) = 4; tOrder(5) = 5; tOrder(6) = 7; tOrder(7) = 6;
    }
    for(uint i = 0; i<tSizeActiveElements; i++ )
    {
        //                std::cout << " (aElementList.ElementListOnProc)(i) " << (aElementList.ElementListOnProc)(i) << std::endl;
        //                if( i == 139)
        //                    std::cout << " i " << i << std::endl;
        tBasis_of_element = trans(give_basis_of_element((aElementList.ElementListOnProc)(i),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements));
        tDesignBSpline_of_element = trans(give_basis_of_element((aElementList.ElementListOnProc)(i),aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements));
        //        tElementLevel = give_element_level((aElementList.ElementListOnProc)(i),aElementList.Dim,aElementList.NumElements);
        (aElementList.BlockSetOutput)(tVar,0) = 0;//tElementLevel; // Blocks are seperated for each level

        //        if( (aElementList.ElementListOnProc)(i) == 4484 )
        //            std::cout << " (aElementList.ElementListOnProc)(i) " << (aElementList.ElementListOnProc)(i) << std::endl;

        if( aElementList.TruncatedBsplines == false)
        {
            this->give_Tmatrix_and_id_field((aElementList.ElementListOnProc)(i),aElementList,aBasisList);
            this->give_Tmatrix_and_id_field_for_designvariables((aElementList.ElementListOnProc)(i),aElementList,aBasisList);
        }
        else
        {
            this->give_Truncated_Tmatrix_and_id_field((aElementList.ElementListOnProc)(i),aElementList,aBasisList);
            this->give_Truncated_Tmatrix_and_id_field_for_designvariables((aElementList.ElementListOnProc)(i),aElementList,aBasisList);
            //           (aElementList.TMatrix).print("(aElementList.TMatrix)");
            //            (aElementList.IdField).print("(aElementList.IdField)");
            //            (aElementList.TMatrixDesign).print("(aElementList.TMatrixDesign)");
            //            (aElementList.IdFieldDesign).print("(aElementList.IdFieldDesign)");
        }
        (aElementList.SparseSize) += pow((aElementList.IdFieldDesign).length(),2);
        //        (aElementList.IdField).print("(aElementList.IdField) ");
        //        (aElementList.TMatrixDesign).print("(aElementList.TMatrixDesign)");
        //        }
        if( (aElementList.IdField).length() > tMaxIdField)
            tMaxIdField = (aElementList.IdField).length(); // Save biggest ID field for static field in MTK

        if( (aElementList.IdFieldDesign).length() > tMaxIdFieldDesign)
            tMaxIdFieldDesign = (aElementList.IdFieldDesign).length(); // Save biggest ID field of design variables for static field in MTK

        //                (aElementList.IdFieldDesign).print("(aElementList.IdFieldDesign) " + std::to_string((aElementList.ElementListOnProc)(i)));
        //        (aElementList.TMatrixDesign).print("(aElementList.TMatrixDesign) " + std::to_string((aElementList.ElementListOnProc)(i)));
        for(uint j = 0; j< pow((aElementList.Polynomial)+1,(aElementList.Dim)); j++)
        {
            if ( (aBasisList.BasisActive).test(tBasis_of_element(j)) == 1 )
            {
                (aElementList.Fetopo)(tVar,tOrder(j)) = tBasis_of_element(j);
                (aElementList.NodalLocaltoGlobal)(tVarb) = tBasis_of_element(j);
                tIdFieldNodalBased(tVarb) = (aElementList.IdField);
                tTmatrixNodalBased(tVarb) = (aElementList.TMatrix).col((j)); //(aElementList.TMatrix).col(tOrder(j));
                tIdFieldNodalBasedDesign(tVarb) = (aElementList.IdFieldDesign);
                tTmatrixNodalBasedDesign(tVarb) = (aElementList.TMatrixDesign).col((j));
                tVarb++;
            }
            else
            {
                tParent = tBasis_of_element(j);
                tCoordinates =  give_coordinate_from_basis(tParent,aElementList.Polynomial,aElementList);
                while( tParent < UINT_MAX )
                {
                    tParent = give_basis_of_parent(tParent,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
                    if( tParent < UINT_MAX && (aBasisList.BasisActive).test(tParent) == 1 )
                        break;
                }
                if ( tParent != UINT_MAX && (aBasisList.BasisActive).test(tParent) == 1 )
                {
                    (aElementList.Fetopo)(tVar,tOrder(j)) = tParent;
                    (aElementList.NodalLocaltoGlobal)(tVarb) = tParent;
                    tIdFieldNodalBased(tVarb) = (aElementList.IdField);
                    tTmatrixNodalBased(tVarb) =  (aElementList.TMatrix).col((j)); //(aElementList.TMatrix).col(tOrder(j));
                    tIdFieldNodalBasedDesign(tVarb) = (aElementList.IdFieldDesign);
                    tTmatrixNodalBasedDesign(tVarb) = (aElementList.TMatrixDesign).col((j));
                    tVarb++;
                }
                else
                {
                    (aElementList.Fetopo)(tVar,tOrder(j)) = tDummyNum;
                    (aElementList.NodalLocaltoGlobal)(tVarb) = tDummyNum;
                    (aElementList.NodalLocaltoGlobalMap)(tVard,0) = tDummyNum;
                    (aElementList.NodalLocaltoGlobalMap)(tVard,1) = tBasis_of_element(j);
                    tVard++;
                    tIdFieldNodalBased(tVarb) = (aElementList.IdField);
                    tTmatrixNodalBased(tVarb) =  (aElementList.TMatrix).col((j)); //(aElementList.TMatrix).col(tOrder(j));
                    tIdFieldNodalBasedDesign(tVarb) = (aElementList.IdFieldDesign);
                    tTmatrixNodalBasedDesign(tVarb) = (aElementList.TMatrixDesign).col((j));
                    DummyNodalLocaltoGlobal(tVarc) = tBasis_of_element(j);
                    tVarc++;
                    tDummyNum++;
                    tVarb++;
                }
            }
        }
        (aElementList.ElemLocaltoGlobal)(tVar) = (aElementList.ElementListOnProc)(i);
        tVar++;
    }
    DummyNodalLocaltoGlobal.resize(tVarc,1);
    (aElementList.NodalLocaltoGlobalMap).resize(tVard,2);
    Mat<uint> tUnique_list = find_unique((aElementList.NodalLocaltoGlobal));
    (aElementList.NodalLocaltoGlobal) = unique((aElementList.NodalLocaltoGlobal));
    aElementList.NodalLocaltoGlobalExist = aElementList.NodalLocaltoGlobal;
    //Calculate the coordiantes for the respective basis functions
    (aElementList.ControlPoints).set_size((aElementList.NodalLocaltoGlobal).length(),aElementList.Dim);
    (aElementList.TMatrixField).set_size((aElementList.NodalLocaltoGlobal).length(),tMaxIdField+1,0); // Tmatrix entries are saved in a Nodal field + Number of supported basis functions
    (aElementList.IdFieldField).set_size((aElementList.NodalLocaltoGlobal).length(),tMaxIdField+1,0); // ID Field entries are saved in a Nodal field + Number of supported basis functions
    (aElementList.TMatrixFieldDesign).set_size((aElementList.NodalLocaltoGlobal).length(),tMaxIdFieldDesign+1,0); // Tmatrix entries are saved in a Nodal field + Number of supported basis functions
    (aElementList.IdFieldFieldDesign).set_size((aElementList.NodalLocaltoGlobal).length(),tMaxIdFieldDesign+1,0); // ID Field entries are saved in a Nodal field + Number of supported basis functions
    tMaxIdField = 0;
    tMaxIdFieldDesign = 0;
    for(uint i = 0; i<(aElementList.NodalLocaltoGlobal).length(); i++)
    {
        tCoordinates =  give_coordinate_from_basis((aElementList.NodalLocaltoGlobal)(i),aElementList.Polynomial,aElementList);
        (aElementList.ControlPoints).row(i) = tCoordinates.row(0);
        //        (aElementList.NodalLocaltoGlobal)(i) += 1;
        tVar = 1;
        for(uint j = 0; j< (tIdFieldNodalBased(tUnique_list(i))).length(); j++)
        {
            //            std::cout << " (aElementList.NodalLocaltoGlobal)(i) " << (aElementList.NodalLocaltoGlobal)(i) << " tUnique_list(i) " << tUnique_list(i) << " (tTmatrixNodalBased(tUnique_list(i)))(j) " << (tTmatrixNodalBased(tUnique_list(i)))(j) << " (tIdFieldNodalBased(tUnique_list(i)))(j) " << (tIdFieldNodalBased(tUnique_list(i)))(j) << std::endl;
            //            ((tTmatrixNodalBased(tUnique_list(i)))).print("(tTmatrixNodalBased(tUnique_list(i)))");
            //            (tIdFieldNodalBased(tUnique_list(i))).print("(tIdFieldNodalBased(tUnique_list(i)))");
            //            if( i == 3070)
            //            {
            //                        std::cout << " (aElementList.NodalLocaltoGlobal)(i) " << (aElementList.NodalLocaltoGlobal)(i) << " tUnique_list(i) " << tUnique_list(i) << " (tTmatrixNodalBased(tUnique_list(i)))(j) " << (tTmatrixNodalBased(tUnique_list(i)))(j) << " (tIdFieldNodalBased(tUnique_list(i)))(j) " << (tIdFieldNodalBased(tUnique_list(i)))(j) << std::endl;
            //            ((tTmatrixNodalBased(tUnique_list(i)))).print("(tTmatrixNodalBased(tUnique_list(i)))");
            //            (tIdFieldNodalBased(tUnique_list(i))).print("(tIdFieldNodalBased(tUnique_list(i)))");
            //            }
            if( (tTmatrixNodalBased(tUnique_list(i)))(j) > 0.0)
            {

                (aElementList.TMatrixField)(i,tVar) = (tTmatrixNodalBased(tUnique_list(i)))(j);
                (aElementList.IdFieldField)(i,tVar) = (tIdFieldNodalBased(tUnique_list(i)))(j);
                tVar++;
            }
        }
        if( tVar > tMaxIdField)
            tMaxIdField = tVar;
        (aElementList.TMatrixField)(i,0) = tVar-1;
        (aElementList.IdFieldField)(i,0) = tVar-1;
        tVar = 1;
        //        (tIdFieldNodalBasedDesign(tUnique_list(i))).print("(tIdFieldNodalBasedDesign(tUnique_list(i))) " + std::to_string(i));
        //        (tTmatrixNodalBasedDesign(tUnique_list(i))).print("(tTmatrixNodalBasedDesign(tUnique_list(i))) " + std::to_string(i));
        tWeight = 0.0;

        for(uint j = 0; j< (tIdFieldNodalBasedDesign(tUnique_list(i))).length(); j++)
        {
            if( (tTmatrixNodalBasedDesign(tUnique_list(i)))(j) > 0.0)
            {
                (aElementList.TMatrixFieldDesign)(i,tVar) = (tTmatrixNodalBasedDesign(tUnique_list(i)))(j);
                tWeight += (tTmatrixNodalBasedDesign(tUnique_list(i)))(j);
                (aElementList.IdFieldFieldDesign)(i,tVar) = (tIdFieldNodalBasedDesign(tUnique_list(i)))(j);
                tVar++;
            }
        }
        if ( aElementList.PerformNormalization )
        {
            (aElementList.TMatrixFieldDesign).row(i) /= tWeight;
        }
        if( tVar > tMaxIdFieldDesign)
            tMaxIdFieldDesign = tVar;
        (aElementList.TMatrixFieldDesign)(i,0) = tVar-1;
        (aElementList.IdFieldFieldDesign)(i,0) = tVar-1;

    }
    (aElementList.TMatrixField).resize((aElementList.NodalLocaltoGlobal).length(),tMaxIdField);
    (aElementList.IdFieldField).resize((aElementList.NodalLocaltoGlobal).length(),tMaxIdField);
    (aElementList.TMatrixFieldDesign).resize((aElementList.NodalLocaltoGlobal).length(),tMaxIdFieldDesign);
    (aElementList.IdFieldFieldDesign).resize((aElementList.NodalLocaltoGlobal).length(),tMaxIdFieldDesign);

    if( aElementList.FilterRadius > 0.0)
    {
        moris::tic updateidandtmatrixdesign;
        this->Update_IDandTMatrix_design(aElementList,aBasisList);
        real tElapsedupdateidandtmatrixdesign = updateidandtmatrixdesign.toc<moris::chronos::seconds>().wall;
        std::fprintf(stdout,"Time for using the filter: %f [sec]\n",tElapsedupdateidandtmatrixdesign);
    }

    aElementList.TrueNodalLocaltoGlobalSize = (aElementList.NodalLocaltoGlobal).n_rows()-tVarc;
    for(uint i = 0; i<tVarc; i++)
    {
        tCoordinates =  give_coordinate_from_basis(DummyNodalLocaltoGlobal(i),aElementList.Polynomial,aElementList);
        tPos = (aElementList.NodalLocaltoGlobal).n_rows()-tVarc+i;
        (aElementList.NodalLocaltoGlobalExist)(tPos) = DummyNodalLocaltoGlobal(i);
        (aElementList.ControlPoints).row(tPos) = tCoordinates.row(0);
    }
    Cell< std::string > BSetNames((aElementList.Level)+1);
    for(uint i = 0; i<((aElementList.Level)+1); i++ )
    {
        BSetNames(i) =  "tBlockset_" + std::to_string(i);
    }
    (aElementList.BSetNames) = BSetNames;
}

void moris::Hierarchical_Mesh::deactivate_element_list(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    Mat<real> Params(aElementList.Dim,1,0);
    Mat<uint> tBasis_of_element(pow(aElementList.Polynomial+1,aElementList.Dim),1,0);
    Mat<real> tLsg(pow(aElementList.Polynomial+1,aElementList.Dim),1,0);
    (aElementList.DeactivateElement).set_size((aElementList.ElementActive).count(),1,0);
    uint tVar = 0;
    for(uint i = 0; i<(aElementList.ElementListOnProc).n_rows(); i++)
    {
        if ((aElementList.ElementActive).test((aElementList.ElementListOnProc)(i)) == 1)
        {
            tBasis_of_element = give_basis_of_element((aElementList.ElementListOnProc)(i),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
            for(uint j = 0; j<pow(aElementList.Polynomial+1,aElementList.Dim); j++)
            {
                Params =  give_coordinate_from_basis(tBasis_of_element(j),aElementList.Polynomial,aElementList);
                tLsg(j) = 0.0;//mySpiralFunc(Params);
            }
            if ( tLsg.max() > 0 && tLsg.min() < 0 )
            {
                (aElementList.DeactivateElement)(tVar) = (aElementList.ElementListOnProc)(i);
                tVar++;
            }
        }
    }
    (aElementList.DeactivateElement).resize(tVar,1);
    aElementList.DeactivateElementInit = aElementList.DeactivateElement;
    (aElementList.DeactivateElementRefLvl).set_size((aElementList.DeactivateElement).length(),1,aElementList.Refinement);
}

void moris::Hierarchical_Mesh::sendrecv_deactivation_and_bcast_element_data(
        ElementList_struc & aElementList)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    if(tProcSize >1)
    {
        MPI_Status status;
        boost::dynamic_bitset<> tBitset1((aElementList.ElementActive).size());
        std::string tBuffer;
        tBitset1 = &(aElementList.ElementActive);
        boost::to_string(tBitset1,tBuffer);
        if(tProcRank == 0)
        {
            char *buf = new char[tBuffer.length()];
            for(uint i = 1; i<tProcSize; i++ )
            {
                MPI_Recv(buf,tBuffer.length(),MPI_CHAR,i,i,MPI_COMM_WORLD, &status);
                std::string tBuffer2(buf,tBuffer.length());
                boost::dynamic_bitset<> tBitset2(tBuffer2);
                tBitset1 = (tBitset1 & tBitset2); // Compare data and save on tBitsetAlpha all bits, which have a one
            }
            delete [] buf;
        }
        if(tProcRank > 0)
        {
            MPI_Send(tBuffer.c_str(),tBuffer.length(),MPI_CHAR,0,tProcRank,MPI_COMM_WORLD);
        }

        boost::to_string(tBitset1,tBuffer);
        broadcast(tBuffer);
        boost::dynamic_bitset<> tBitset3(tBuffer);
        (aElementList.ElementActive) = tBitset3;
    }
}

void moris::Hierarchical_Mesh::sendrecv_refinement_and_bcast_element_data(
        ElementList_struc & aElementList)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    if(tProcSize >1)
    {
        MPI_Status status;
        boost::dynamic_bitset<> tBitset1;
        std::string tBuffer;
        tBitset1 = &(aElementList.ElementActive);
        boost::to_string(tBitset1,tBuffer);
        if(tProcRank == 0)
        {
            char *buf = new char[tBuffer.length()];
            for(uint i = 1; i<tProcSize; i++ )
            {
                MPI_Recv(buf,tBuffer.length(),MPI_CHAR,i,i,MPI_COMM_WORLD, &status);
                std::string tBuffer2(buf,tBuffer.length());
                boost::dynamic_bitset<> tBitset2(tBuffer2);
                tBitset1 = (tBitset1 | tBitset2); // Compare data and save on tBitsetAlpha all bits, which have a one
            }
            delete [] buf;
        }
        if(tProcRank > 0)
        {
            MPI_Send(tBuffer.c_str(),tBuffer.length(),MPI_CHAR,0,tProcRank,MPI_COMM_WORLD);
        }

        boost::to_string(tBitset1,tBuffer);
        broadcast(tBuffer);
        boost::dynamic_bitset<> tBitset3(tBuffer);
        (aElementList.ElementActive) = tBitset3;
    }
}

void moris::Hierarchical_Mesh::sendrecv_refinement_and_bcast_basis_data(
        BasisList_struc & aBasisList)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    if(tProcSize >1)
    {
        MPI_Status status;
        boost::dynamic_bitset<> tBitset1;
        std::string tBuffer;
        tBitset1 = &(aBasisList.BasisActive);
        boost::to_string(tBitset1,tBuffer);
        if(tProcRank == 0)
        {
            char *buf = new char[tBuffer.length()];
            for(uint i = 1; i<tProcSize; i++ )
            {
                MPI_Recv(buf,tBuffer.length(),MPI_CHAR,i,i,MPI_COMM_WORLD, &status);
                std::string tBuffer2(buf,tBuffer.length());
                boost::dynamic_bitset<> tBitset2(tBuffer2);
                tBitset1 = (tBitset1 | tBitset2); // Compare data and save on tBitsetAlpha all bits, which have a one
            }
            delete [] buf;
        }
        if(tProcRank > 0)
        {
            MPI_Send(tBuffer.c_str(),tBuffer.length(),MPI_CHAR,0,tProcRank,MPI_COMM_WORLD);
        }
        boost::to_string(tBitset1,tBuffer);
        broadcast(tBuffer);
        boost::dynamic_bitset<> tBitset3(tBuffer);
        (aBasisList.BasisActive) = tBitset3;
    }
}

void moris::Hierarchical_Mesh::sendrecv_deactivation_and_bcast_basis_data(
        BasisList_struc & aBasisList)
{
    uint tProcSize = par_size();
    uint tProcRank = par_rank();
    if(tProcSize >1)
    {
        MPI_Status status;
        boost::dynamic_bitset<> tBitset1;
        std::string tBuffer;
        tBitset1 = &(aBasisList.BasisActive);
        boost::to_string(tBitset1,tBuffer);
        if(tProcRank == 0)
        {
            char *buf = new char[tBuffer.length()];
            for(uint i = 1; i<tProcSize; i++ )
            {
                MPI_Recv(buf,tBuffer.length(),MPI_CHAR,i,i,MPI_COMM_WORLD, &status);
                std::string tBuffer2(buf,tBuffer.length());
                boost::dynamic_bitset<> tBitset2(tBuffer2);
                tBitset1 = (tBitset1 & tBitset2); // Compare data and save on tBitsetAlpha all bits, which have a one
            }
            delete [] buf;
        }
        if(tProcRank > 0)
        {
            MPI_Send(tBuffer.c_str(),tBuffer.length(),MPI_CHAR,0,tProcRank,MPI_COMM_WORLD);
        }
        boost::to_string(tBitset1,tBuffer);
        broadcast(tBuffer);
        boost::dynamic_bitset<> tBitset3(tBuffer);
        (aBasisList.BasisActive) = tBitset3;
    }
}
void moris::Hierarchical_Mesh::gather_refinementlevel_and_bcast_data(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    int tProcSize = (int)par_size();
    uint tProcRank = par_rank();
    if(tProcSize >1)
    {
        int level = (int)(aElementList.Level);
        int *recvbuf;
        recvbuf = (int *)malloc(tProcSize*1*sizeof(int));
        uint tMaxLevel = 0;
        int tMaxLevell = 0;
        if(tProcRank == 0)
        {
            MPI_Gather(&level,1,MPI_INT,recvbuf,1,MPI_INT,0,MPI_COMM_WORLD);
            tMaxLevell = *std::max_element(recvbuf, recvbuf+tProcSize);
            tMaxLevel = (uint)tMaxLevell;
        }
        if(tProcRank > 0)
        {
            MPI_Gather(&level,1,MPI_INT,recvbuf,1,MPI_INT,0,MPI_COMM_WORLD);
        }
        broadcast(tMaxLevel);
        (aElementList.Level) = tMaxLevel;
    }
    (aElementList.NumberElements) = give_number_of_elements(aElementList.Level,aElementList.Dim,aElementList.NumElements);
    (aBasisList.NumberBasis) = give_number_of_basis(aElementList.Level,aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);

    if( (aElementList.ElementActive).size() < (aElementList.NumberElements))
        (aElementList.ElementActive).resize((aElementList.NumberElements));

    if( (aBasisList.BasisActive).size() < (aBasisList.NumberBasis))
        (aBasisList.BasisActive).resize((aBasisList.NumberBasis));

    if( (aElementList.ElementActiveDesign).size() < (aElementList.NumberElements))
        (aElementList.ElementActiveDesign).resize((aElementList.NumberElements));

    if( (aBasisList.DesignBSplineActive).size() < (aBasisList.NumberBasis))
        (aBasisList.DesignBSplineActive).resize((aBasisList.NumberBasis));
}

void moris::Hierarchical_Mesh::gather_value_and_bcast_max(
        uint & aMessage)
{
    int tProcSize = (int)par_size();
    uint tProcRank = par_rank();
    if(tProcSize >1)
    {
        int aMessageInt = (int)aMessage;
        int *recvbuf;
        recvbuf = (int *)malloc(tProcSize*1*sizeof(int));
        int tMaxValuee = 0;
        if(tProcRank == 0)
        {
            MPI_Gather(&aMessageInt,1,MPI_INT,recvbuf,1,MPI_INT,0,MPI_COMM_WORLD);
            tMaxValuee = *std::max_element(recvbuf, recvbuf+tProcSize);
            aMessage = (uint)tMaxValuee;
        }
        if(tProcRank > 0)
        {
            MPI_Gather(&aMessageInt,1,MPI_INT,recvbuf,1,MPI_INT,0,MPI_COMM_WORLD);
        }
        broadcast(aMessage);
    }
}
//void moris::Hierarchical_Mesh::sendrecv_refinement_of_neighbours_and_bcast_element_data(
//        ElementList_struc & aElementList)
//{
//    mpi::communicator world;
//    uint tProcSize = par_size();
//    if(tProcSize >1)
//    {
//        MPI_Status status;
//        MPI_Datatype tBitset;
//        MPI_Type_contiguous(1,MPI_DOUBLE,&tBitset);
//        MPI_Type_commit(&tBitset);
//        moris::MPITools tWorld;
//        uint tSizeOfDatatype = ceil((aElementList.ElementActive).size()/64.0); // MPI_DOUBLE = 8 bytes = 64bits
//        std::bitset<HMR_BITSET_MAXSIZE> tBitsetAlpha;
//        std::bitset<HMR_BITSET_MAXSIZE> tBitsetBeta;
//        tBitsetAlpha = &(aElementList.ElementActive);
//        Mat<uint> tProcList = unique((aElementList.ProcNeighbour));
//        if( tProcList(tProcList.length()-1) == UINT_MAX)
//            tProcList.resize(tProcList.length()-1,1);
//        for(uint i = 0; i<tProcList.length(); i++)
//        {
//            MPI_Sendrecv(&tBitsetAlpha,tSizeOfDatatype,tBitset,tProcList(i),201,&tBitsetBeta,tSizeOfDatatype,tBitset,tProcList(i),201,world,&status);
//            tBitsetAlpha = (tBitsetAlpha | tBitsetBeta); // Compare data and save on tBitsetAlpha all bits, which have a one
//        }
//        (aElementList.ElementActive) = tBitsetAlpha;
//    }
//}

//void moris::Hierarchical_Mesh::sendrecv_deactivation_of_neighbours_and_bcast_element_data(
//        ElementList_struc & aElementList)
//{
//    mpi::communicator world;
//    uint tProcSize = par_size();
//    if(tProcSize >1)
//    {
//        MPI_Status status;
//        MPI_Datatype tBitset;
//        MPI_Type_contiguous(1,MPI_DOUBLE,&tBitset);
//        MPI_Type_commit(&tBitset);
//        moris::MPITools tWorld;
//        uint tSizeOfDatatype = ceil((aElementList.ElementActive).size()/64.0); // MPI_DOUBLE = 8 bytes = 64bits
//        std::bitset<HMR_BITSET_MAXSIZE> tBitsetAlpha;
//        std::bitset<HMR_BITSET_MAXSIZE> tBitsetBeta;
//        tBitsetAlpha = &(aElementList.ElementActive);
//        for(uint i = 0; i<(aElementList.ProcNeighbour).length(); i++)
//        {
//            if( (aElementList.ProcNeighbour)(i) != UINT_MAX )
//            {
//                MPI_Sendrecv(&tBitsetAlpha,tSizeOfDatatype,tBitset,(aElementList.ProcNeighbour)(i),201,&tBitsetBeta,tSizeOfDatatype,tBitset,(aElementList.ProcNeighbour)(i),201,world,&status);
//                tBitsetAlpha = (tBitsetAlpha & tBitsetBeta); // Compare data and save on tBitsetAlpha all bits, which have a one
//            }
//        }
//        (aElementList.ElementActive) = tBitsetAlpha;
//    }
//}

void moris::Hierarchical_Mesh::Initial_mesh(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    //Activate basis functions
    this->hierarchical_basisfunction_refinement(aElementList,aBasisList);

    //Update all side and node sets
    this->update_side_node_set(aElementList,aBasisList);

    //Create the data for MTK
    this->create_mesh_data(aElementList,aBasisList);

    this->give_node_proc_owner(aElementList,aBasisList);
}

void moris::Hierarchical_Mesh::Find_elements_in_stencil(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    this->give_active_elements(aElementList);
    Mat<uint> tElementStencil;
    Mat<uint> tElementStencilDensity;
    //    Mat<uint> tElementStencilDensityReal;
    Mat<uint> tElementNeighbours;
    Mat<uint> tRefineElement((aElementList.ElementListOnProc).length(),1,0);
    Mat<uint> tRefineElementLvL((aElementList.ElementListOnProc).length(),1,0);
    uint tVar = 0;
    uint tSum = 0;
    Mat<uint> tBasis;
    Mat<uint> tChildList;
    //    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    //        uint tVara;
    //        uint tMaxNumElemLevel;
    //        uint tParent;
    //        uint tSizeDensityField = (aElementList.Density).length();
    //        uint tNumberElements;
    //        uint tMaxChild;
    //        uint tLevel;

    for(uint e = 0; e < (aElementList.ElementListOnProc).length(); e++)
    {
        if( (aElementList.Density)((aElementList.ElementListOnProc)(e)) > aElementList.FeatureLowerBound )
        {
            tElementStencil = this->give_neighbour_stencil_of_element(aElementList.FeatureResolution,(aElementList.ElementListOnProc)(e),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements,aElementList.ElementActive);
            //tElementStencil.print("tElementStencil");
            //            std::cout << " e " << e << " (aElementList.ElementListOnProc)(e) " << (aElementList.ElementListOnProc)(e)  << "  (aElementList.Density)((aElementList.ElementListOnProc)(e)) " << (aElementList.Density)((aElementList.ElementListOnProc)(e)) <<std::endl;
            Mat<uint> Pos  = give_position_of_element((aElementList.ElementListOnProc)(e),aElementList.Dim,aElementList.NumElements);
            //Pos.print("Pos");
            tElementStencilDensity.set_size(tElementStencil.n_rows(),tElementStencil.n_cols(),0);
            //            tElementStencilDensityReal.set_size(tElementStencil.n_rows(),tElementStencil.n_cols(),0);
            for(uint i = 0; i < tElementStencil.n_rows(); i++)
            {
                for(uint j = 0; j < tElementStencil.n_cols(); j++)
                {
                    if( tElementStencil(i,j) > 0 && tElementStencil(i,j) < UINT_MAX )
                    {
                        if( (aElementList.Density)(tElementStencil(i,j)) < aElementList.FeatureLowerBound )
                        {
                            tElementStencilDensity(i,j) = 0;
                            //                            tElementStencilDensityReal(i,j) = (aElementList.Density)(tElementStencil(i,j)) ;
                        }
                        else
                        {
                            tElementStencilDensity(i,j) = 1;
                            //                            tElementStencilDensityReal(i,j) = (aElementList.Density)(tElementStencil(i,j)) ;
                        }
                    }
                }
            }
            tElementNeighbours.set_size(1,tElementStencil.n_cols());

            //if ( e ==  341)
            //std::cout << " e " << e << " (aElementList.ElementListOnProc)(e) " << (aElementList.ElementListOnProc)(e) << " tVar " << tVar << std::endl;
            //            tElementStencil.print("tElementStencil");
            //            tElementStencilDensity.print("tElementStencilDensity");
            //            tElementStencilDensityReal.print("tElementStencilDensityReal");
            Mat<uint> tCount(tElementStencil.n_rows(),1,0);
            for(uint i = 0; i < tElementStencil.n_rows(); i++)
            {
                tSum = 1; // Count number of Elements with a density higher then DensityLowerBound
                tElementNeighbours.row(0) = tElementStencil.row(i);

                for(uint j = aElementList.FeatureResolution+2; j < tElementStencil.n_cols(); j++)
                {
                    if( tElementStencilDensity(i,j) == 1 || tElementStencil(i,j) == UINT_MAX )
                    {
                        tSum++;
                    }
                    else
                    {
                        break;
                    }
                }
                for(uint j = aElementList.FeatureResolution; j > 0; j--)
                {
                    if( tElementStencilDensity(i,j) == 1 || tElementStencil(i,j) == UINT_MAX )
                    {
                        tSum++;
                    }
                    else
                    {
                        break;
                    }
                }
                tCount(i) = tSum;
            }
            if( tCount.min() < aElementList.FeatureResolution)
            {
                tRefineElement(tVar) = (aElementList.ElementListOnProc)(e);
                tRefineElementLvL(tVar) = (uint)( log2( aElementList.FeatureResolution - tCount.min()));
                //                                if( tRefineElementLvL(tVar) > 0)
                //                                {
                //                                    tMaxNumElemLevel = 0;
                //                                    for(uint level = 1; level<=tRefineElementLvL(tVar) ; level ++)
                //                                        tMaxNumElemLevel += pow(pow(2,aElementList.Dim),level);
                //                                    tChildList.set_size(tMaxNumElemLevel,1,UINT_MAX);
                //                                    tMaxNumElemLevel -= pow(pow(2,aElementList.Dim),tRefineElementLvL(tVar));
                //                                    tChildren = give_children_of_element((aElementList.ElementListOnProc)(e),aElementList.Dim,aElementList.NumElements);
                //                                    tChildList.rows(0,tChildren.length()-1) = tChildren.rows(0,tChildren.length()-1);
                //                                    tVara = tChildren.length();
                //                                    for(uint k = 0; k<tMaxNumElemLevel; k++)
                //                                    {
                //                                        tChildren = give_children_of_element(tChildList(k),aElementList.Dim,aElementList.NumElements);
                //                                        tChildList.rows(tVara+k*pow(2,aElementList.Dim),tVara+(k+1)*pow(2,aElementList.Dim)-1) = tChildren.rows(0,tChildren.length()-1);
                //                                    }
                //                                    tChildList = sort(tChildList);
                //                                    if( tSizeDensityField < tChildList.max())
                //                                    {
                //                                        tMaxChild = tChildList.max();
                //                                        tLevel =   give_element_level(tMaxChild,aElementList.Dim,aElementList.NumElements);
                //                                        tNumberElements = give_number_of_elements(tLevel,aElementList.Dim,aElementList.NumElements);
                //                                        (aElementList.Density).resize(tNumberElements,1);
                //                                        tSizeDensityField = tNumberElements;
                //                                    }
                //                                    for(uint k = 0; k < tChildList.length(); k++)
                //                                    {
                //                                        tParent = give_parent_of_element(tChildList(k),aElementList.Dim,aElementList.NumElements);
                //                                        (aElementList.Density)(tChildList(k)) = (aElementList.Density)(tParent);
                //                                    }
                //                                }
                tVar++;
            }
        }
    }
    tRefineElement.resize(tVar,1);
    tRefineElementLvL.resize(tVar,1);
    //Find parent elements on level 0
    (aElementList.DeactivateElement).set_size(tVar,1,0);
    (aElementList.DeactivateElement) = tRefineElement;
    (aElementList.DeactivateElementRefLvl).set_size(tVar,1,0);
    (aElementList.DeactivateElementRefLvl) = tRefineElementLvL;

    aElementList.DeactivateElementInit = aElementList.DeactivateElement;
}

void moris::Hierarchical_Mesh::find_elements_for_refinement(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tMaxElement = (aElementList.ElemLocaltoGlobalLastStepFEM).max();
    Mat<uint> tBasisOfElement = give_basis_of_element(tMaxElement,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
    uint tMaxBasis = tBasisOfElement.max();
    uint tVar = 0;
    //    real tWeight = 0.0;
    Mat<uint> tChildren;
    uint tLevel = give_element_level((aElementList.ElemLocaltoGlobalLastStepFEM)((aElementList.ElemLocaltoGlobalLastStepFEM).length()-1),aElementList.Dim,aElementList.NumElements);
    tLevel++; // Add one level for the density of the children
    uint tNumberElements = give_number_of_elements(tLevel,aElementList.Dim,aElementList.NumElements);
    uint tMaxRefinementElement;
    (aElementList.Density).resize(tNumberElements,1);
    (aElementList.ElementActive).reset();
    uint tSizeBitset = (aElementList.ElementActiveDummy).size();
    (aElementList.ElementActive).resize(tSizeBitset);
    aElementList.ElementActive = aElementList.ElementActiveDummy;

    if( (aElementList.ElementActive).size() < tNumberElements )
        (aElementList.ElementActive).resize(tNumberElements);

    BoostBitset tElementActiveDummy = aElementList.ElementActive;

    if( (aBasisList.NodalPDV).length() < tMaxBasis)
        (aBasisList.NodalPDV).resize(tMaxBasis,1);
    if( (aBasisList.NodalADVNewField).length() < tMaxBasis)
        (aBasisList.NodalADVNewField).resize(tMaxBasis,1);

    Mat<real> tKnotVector(2*(aElementList.PolynomialDesign+1),1,0);
    for(uint i = 0; i < tKnotVector.length(); i++)
        tKnotVector(i) = -(real)aElementList.PolynomialDesign + (real)i;
    Mat<real> tXi(aElementList.Dim,1,0.5);
    Mat<uint> tElement_flag(aElementList.Dim,1,0);
    Mat<real> Bspline = moris::Bspline::build_spline_uniform_nd(tXi,tKnotVector,aElementList.PolynomialDesign,tElement_flag,aElementList.Dim);
    Mat<real> tNodalLvLSet(pow(aElementList.PolynomialDesign+1,aElementList.Dim),1); // Level set field of an element
    Mat<uint> tRefineElement(2*(aElementList.ElemLocaltoGlobalLastStepFEM).length(),1,0); // Possible are number of elements on proc * number of levels of children times number of elements, which have support with. Two times, because the density or level set can define an element
    // Check for element refinement. Check elemental density threshold and, if enabled, intersected elements of a level set field
    for(uint i = 0; i<(aElementList.ElemLocaltoGlobalLastStepFEM).length(); i++)
    {
        if( aElementList.TruncatedBsplines == false)
        {
            this->give_Tmatrix_and_id_field_for_designvariables_projection((aElementList.ElemLocaltoGlobalLastStepFEM)(i),aElementList,aBasisList);
        }
        else
        {
            this->give_Truncated_Tmatrix_and_id_field_for_designvariables_projection((aElementList.ElemLocaltoGlobalLastStepFEM)(i),aElementList,aBasisList);
        }
        tBasisOfElement = give_basis_of_element((aElementList.ElemLocaltoGlobalLastStepFEM)(i),aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
        for(uint j = 0; j < (aElementList.TMatrix).n_cols(); j++)
        {
            (aBasisList.NodalPDV)(tBasisOfElement(j)) = 0.0;
            (aBasisList.NodalADVNewField)(tBasisOfElement(j)) = 0.0;
            tNodalLvLSet(j) = 0.0;
            for(uint k = 0; k < (aElementList.TMatrix).n_rows(); k++)
            {
                (aBasisList.NodalPDV)(tBasisOfElement(j)) += (aElementList.TMatrix)(k,j) * (aBasisList.NodalField)((aElementList.IdField)(k));
                (aBasisList.NodalADVNewField)(tBasisOfElement(j)) += (aElementList.TMatrix)(k,j) * (aBasisList.NodalADV)((aElementList.IdField)(k));
                //                std::cout << " (aBasisList.NodalADVNewField)(tBasisOfElement(j)) " << (aBasisList.NodalADVNewField)(tBasisOfElement(j)) << " (aBasisList.NodalADV)((aElementList.IdField)(k)) " << (aBasisList.NodalADV)((aElementList.IdField)(k)) << " (aBasisList.NodalField)((aElementList.IdField)(k)) " << (aBasisList.NodalField)((aElementList.IdField)(k)) << " (aBasisList.NodalPDV)(tBasisOfElement(j)) " << (aBasisList.NodalPDV)(tBasisOfElement(j)) <<  std::endl;
                tNodalLvLSet(j) += (aElementList.TMatrix)(k,j) * (aBasisList.NodalLvlSet)((aElementList.IdField)(k));
            }
            if(  (aBasisList.NodalPDV)(tBasisOfElement(j)) > 1)
                (aBasisList.NodalPDV)(tBasisOfElement(j)) = 1;
            else if(  (aBasisList.NodalPDV)(tBasisOfElement(j)) < 0)
                (aBasisList.NodalPDV)(tBasisOfElement(j)) = 0;
        }
        (aElementList.Density)((aElementList.ElemLocaltoGlobalLastStepFEM)(i)) = 0.0;
        for( uint k = 0; k < Bspline.length(); k++)
        {
            (aElementList.Density)((aElementList.ElemLocaltoGlobalLastStepFEM)(i)) += Bspline(k) * (aBasisList.NodalPDV)(tBasisOfElement(k));
        }
        tLevel = give_element_level((aElementList.ElemLocaltoGlobalLastStepFEM)(i),aElementList.Dim,aElementList.NumElements);
        //The element needs to be refined, if the density of an element is within a threshold
        if( (aElementList.Density)((aElementList.ElemLocaltoGlobalLastStepFEM)(i)) >= (aElementList.DensityLowerBound+aElementList.DensityIncrementBound*tLevel) &&
                (aElementList.Density)((aElementList.ElemLocaltoGlobalLastStepFEM)(i)) <= aElementList.DensityUpperBound    )
        {
            tChildren = give_children_of_element((aElementList.ElemLocaltoGlobalLastStepFEM)(i),aElementList.Dim,aElementList.NumElements);
            for(uint j = 0; j < tChildren.length(); j++)
                (aElementList.Density)(tChildren(j)) = (aElementList.Density)((aElementList.ElemLocaltoGlobalLastStepFEM)(i));

            tLevel = give_element_level((aElementList.ElemLocaltoGlobalLastStepFEM)(i),aElementList.Dim,aElementList.NumElements);
            if( tLevel < aElementList.MaxDesignLevelOfRefinement)
            {
                tRefineElement(tVar) = (aElementList.ElemLocaltoGlobalLastStepFEM)(i);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = give_children_of_element(tRefineElement(tVar),aElementList.Dim,aElementList.NumElements);
                for(uint l = 0; l < tChildren.length(); l++)
                    tElementActiveDummy.set(tChildren(l));
                tVar++;
            }
            else
            {
                tLevel = aElementList.MaxDesignLevelOfRefinement-1;
                tRefineElement(tVar) = give_parent_of_level_x((aElementList.ElemLocaltoGlobalLastStepFEM)(i),aElementList.Dim,aElementList.NumElements,tLevel);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = give_children_of_element(tRefineElement(tVar),aElementList.Dim,aElementList.NumElements);
                for(uint l = 0; l < tChildren.length(); l++)
                    tElementActiveDummy.set(tChildren(l));
                tVar++;
            }
        }
        // If this is an intersected element, this element needs to be refined
        if( aElementList.LvlSetMethod == true && tNodalLvLSet.min() < 0 && tNodalLvLSet.max() > 0)
        {
            tLevel = give_element_level((aElementList.ElemLocaltoGlobalLastStepFEM)(i),aElementList.Dim,aElementList.NumElements);
            if( tLevel < aElementList.MaxLevelSetLevelOfRefinement)
            {
                tRefineElement(tVar) = (aElementList.ElemLocaltoGlobalLastStepFEM)(i);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = give_children_of_element(tRefineElement(tVar),aElementList.Dim,aElementList.NumElements);
                for(uint l = 0; l < tChildren.length(); l++)
                    tElementActiveDummy.set(tChildren(l));
                tVar++;
            }
            else
            {
                tLevel = aElementList.MaxLevelSetLevelOfRefinement-1;
                tRefineElement(tVar) = give_parent_of_level_x((aElementList.ElemLocaltoGlobalLastStepFEM)(i),aElementList.Dim,aElementList.NumElements,tLevel);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = give_children_of_element(tRefineElement(tVar),aElementList.Dim,aElementList.NumElements);
                for(uint l = 0; l < tChildren.length(); l++)
                    tElementActiveDummy.set(tChildren(l));
                tVar++;
            }
        }
    }
    tRefineElement.resize(tVar,1);
    tRefineElement = unique(tRefineElement);

    // Parent elements need also be on the refinement list
    if( tVar > 0)
    {
        tVar = tRefineElement.length();
        uint tLengthRefinement = tRefineElement.length();
        tRefineElement.resize(tRefineElement.length()*(tLevel+1),1);
        if( tLevel > 0)
        {
            for(uint i = 0; i < tLengthRefinement; i++ )
            {
                tLevel = give_element_level(tRefineElement(i),aElementList.Dim,aElementList.NumElements);
                for(uint j = 0; j < tLevel; j++)
                {
                    tRefineElement(tVar) = give_parent_of_level_x(tRefineElement(i),aElementList.Dim,aElementList.NumElements,j);
                    tElementActiveDummy.reset(tRefineElement(tVar));
                    tChildren = give_children_of_element(tRefineElement(tVar),aElementList.Dim,aElementList.NumElements);
                    for(uint l = 0; l < tChildren.length(); l++)
                        tElementActiveDummy.set(tChildren(l));
                    tVar++;
                }
            }
        }
        tRefineElement.resize(tVar,1);
        tRefineElement = unique(tRefineElement);
    }

    if( tVar > 0 )
    {
        tMaxRefinementElement = tRefineElement.max();
        tLevel = give_element_level(tMaxRefinementElement,aElementList.Dim,aElementList.NumElements) + 1;
        tNumberElements = give_number_of_elements(tLevel,aElementList.Dim,aElementList.NumElements);
        if( (aElementList.ElementActive).size() < tNumberElements )
        {
            (aElementList.ElementActive).resize(tNumberElements);
            tElementActiveDummy.resize(tNumberElements);
        }
    }

    // Add a buffer layer of elements
    if( tVar > 0 && aElementList.BufferElements > 0)
    {
        Mat<uint> tChildren;
        Mat<uint> tNeighbour(pow(1+2*aElementList.BufferElements*pow(2,aElementList.Level),aElementList.Dim),1,0); // Highest possible number of neigbours
        tVar = tRefineElement.length(); // Starting point to save new elements, which need to be deactivated
        uint tLengthRefinement = tVar; // Temporary variable for length of loop
        uint tBuffer = 0;
        tRefineElement.resize(tRefineElement.length()*(1+tNeighbour.length()),1); // Resize list of elements to the maximum possible list
        for(uint j = 0; j < tLengthRefinement; j++)
        {
            tLevel = give_element_level(tRefineElement(j),aElementList.Dim,aElementList.NumElements);
            if( aElementList.AdaptiveBufferLayer == true)
            {
                if( aElementList.Staircasebuffer == true)
                {
                    tBuffer = aElementList.BufferElements*(1+tLevel);
                }
                else
                {
                    tBuffer = aElementList.BufferElements*pow(2,tLevel);
                }
            }
            else
            {
                tBuffer = aElementList.BufferElements;
            }
            tNeighbour = give_neighbour_of_element(tRefineElement(j),aElementList.Dim,tBuffer,aElementList.NumElements);
            for(uint k = 0; k<tNeighbour.length(); k++)
            {
                if( tNeighbour(k) < tElementActiveDummy.size() &&  tElementActiveDummy.test(tNeighbour(k)) == 1 )
                {
                    tRefineElement(tVar) = tNeighbour(k);
                    tElementActiveDummy.reset(tNeighbour(k));
                    tChildren = give_children_of_element(tNeighbour(k),aElementList.Dim,aElementList.NumElements);
                    for(uint l = 0; l < tChildren.length(); l++)
                        tElementActiveDummy.set(tChildren(l));
                    tVar++;
                }
            }
        }
        tRefineElement.resize(tVar,1);
        tRefineElement = unique(tRefineElement);
    }

    // Parent elements need also be on the refinement list
    if( tVar > 0)
    {
        tMaxRefinementElement = tRefineElement.max();
        tLevel = give_element_level(tMaxRefinementElement,aElementList.Dim,aElementList.NumElements);
        tVar = tRefineElement.length();
        uint tLengthRefinement = tRefineElement.length();
        tRefineElement.resize(tRefineElement.length()*(tLevel+1),1);
        if( tLevel > 0)
        {
            for(uint i = 0; i < tLengthRefinement; i++ )
            {
                tLevel = give_element_level(tRefineElement(i),aElementList.Dim,aElementList.NumElements);
                for(uint j = 0; j < tLevel; j++)
                {
                    tRefineElement(tVar) = give_parent_of_level_x(tRefineElement(i),aElementList.Dim,aElementList.NumElements,j);
                    tElementActiveDummy.reset(tRefineElement(tVar));
                    tChildren = give_children_of_element(tRefineElement(tVar),aElementList.Dim,aElementList.NumElements);
                    for(uint l = 0; l < tChildren.length(); l++)
                        tElementActiveDummy.set(tChildren(l));
                    tVar++;
                }
            }
        }
        tRefineElement.resize(tVar,1);
        tRefineElement = unique(tRefineElement);
    }

    (aElementList.DeactivateElement) = tRefineElement;
    aElementList.DeactivateElementInit = aElementList.DeactivateElement;
}

Cell<Cell<Mat<uint>>>
moris::Hierarchical_Mesh::Floodfill_for_elementfield(
        Mat<uint> & aElementField,
        Mat<uint> & aListOfElement,
        ElementList_struc & aElementList)
{
    Mat<uint> tElementOfObject;
    Mat<uint> tHelp;
    Mat<uint> tWhereInVec;
    Mat<uint> tElementNeighbour;
    BoostBitset tActiveElements(aListOfElement.max()+1);
    Cell<Cell<Mat<uint>>> tObjects(aElementField.max()+1);
    uint iObj = 0, jObj = 0, iCount = 0;
    uint tBuffer = 1; // One layer of neighbour elements will be checked.
    uint tNumberOfUsedCells = 0;
    for(uint icell = 0; icell<= aElementField.max(); icell++)
    {
        Cell<Mat<uint>> tObjectsForPhase(10); // Cell of objects, which includes all objects for a phase. Default (random) =  10
        tHelp = ( aElementField == icell); // Find elements, which have the same colour (number)
        tWhereInVec = find( tHelp ); // Where are the colours (numbers)
        tActiveElements.reset();
        for(uint i = 0; i < tWhereInVec.length(); i++)
            tActiveElements.set(aListOfElement(tWhereInVec(i))); // Activate all elements of the respective phase
        while( tActiveElements.count() > 0) // Search for elements until all elements are deactivated
        {
            iObj = 1;
            iCount = 0;
            tElementOfObject.set_size(tActiveElements.count(),1,UINT_MAX); // Number of possible elements
            tElementOfObject(iCount) = tActiveElements.find_first(); // Start with the first active element
            tActiveElements.reset(tElementOfObject(iCount));
            while( tActiveElements.count() > 0 && tElementOfObject(iCount) < UINT_MAX) // Check elements, until all elements are deactivated or it founds no new neighbours
            {
                tElementNeighbour = give_neighbour_of_element_for_floodfill(tElementOfObject(iCount),aElementList.Dim,tBuffer,aElementList.NumElements);
                for(uint i = 0; i < tElementNeighbour.length(); i++)
                {
                    if( tElementNeighbour(i) < tActiveElements.size() && tActiveElements.test(tElementNeighbour(i)) == 1 ) // Check neighbours if they have the same phase and not checked yet
                    {
                        tElementOfObject(iObj) = tElementNeighbour(i);
                        tActiveElements.reset(tElementNeighbour(i));
                        iObj++;
                    }
                }
                iCount++;
            }
            tElementOfObject.resize(iObj,1); // Size the list of elements of one object
            if( tObjectsForPhase.size() < jObj+1) // If more objects then default, then resize!
                tObjectsForPhase.resize(jObj+1);
            tObjectsForPhase(jObj) = tElementOfObject;
            jObj++;
        }
        tObjectsForPhase.resize(jObj);
        tObjects(icell) = tObjectsForPhase; // Save objects of a phase in a cell
        tNumberOfUsedCells = icell+1;
    }
    if( tNumberOfUsedCells < tObjects.size() )
        tObjects.resize(tNumberOfUsedCells);
    return tObjects;
}

Cell<Cell<Mat<uint>>>
moris::Hierarchical_Mesh::Floodfill_for_levelset(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    Mat<uint> tBasisOfObject;
    Mat<uint> tElementOfObject;
    Mat<uint> tHelp;
    Mat<uint> tWhereInVec;
    Mat<uint> tBasisNeighbour;
    Mat<uint> tElementOfBasis;
    BoostBitset tActiveBasis(aBasisList.NumberBasis+1);
    Cell<Cell<Mat<uint>>> tObjects((aElementList.ElementField).max()+1);
    uint iObj = 0, jObj = 0, iCount = 0;
    uint tVar = 0;
    uint tBuffer = 1; // One layer of neighbour elements will be checked.
    uint tNumberOfUsedCells = 0;
    Mat<uint> tCheckCrossBasis, tDiagBasis;
    if( aElementList.Dim == 2)
    {
        Mat<uint> tHelpa = {{1, 3, 5, 7}}; // Basis functions in the directions of the coordinate system
        tCheckCrossBasis = tHelpa;
        Mat<uint> tHelpb = {{0, 2, 6, 8}}; // Basis functions in the diagonal direction of the coordinate system
        tDiagBasis = tHelpb;
    }
    else if( aElementList.Dim == 3)
    {
        Mat<uint> tHelpa = {{4, 10, 12, 14, 16, 22}}; // Basis functions in the directions of the coordinate system
        tCheckCrossBasis = tHelpa;
        Mat<uint> tHelpb = {{0, 1, 2, 3, 5, 6, 7, 8, 9, 11, 15, 17, 18, 19, 20, 21, 23, 24, 25, 26}}; // Basis functions in the diagonal direction of the coordinate system
        tDiagBasis = tHelpb;
    }
    for(uint icell = 0; icell<= (aElementList.ElementField).max(); icell++)
    {
        Cell<Mat<uint>> tObjectsForPhase(10); // Cell of objects, which includes all objects for a phase. Default (random) =  10
        tHelp = ( (aBasisList.NodalField) == icell); // Find elements, which have the same colour (number)
        tWhereInVec = find( tHelp ); // Where are the colours (numbers)
        tActiveBasis.reset();
        for(uint i = 0; i < tWhereInVec.length(); i++)
            tActiveBasis.set(tWhereInVec(i)); // Activate all elements of the respective phase
        while( tActiveBasis.count() > 0) // Search for elements until all elements are deactivated
        {
            iObj = 1;
            iCount = 0;
            tBasisOfObject.set_size(tActiveBasis.count(),1,UINT_MAX); // Number of possible elements
            tBasisOfObject(iCount) = tActiveBasis.find_first(); // Start with the first active element
            tActiveBasis.reset(tBasisOfObject(iCount));
            while( tActiveBasis.count() > 0 && tBasisOfObject(iCount) < UINT_MAX) // Check elements, until all elements are deactivated or it founds no new neighbours
            {
                tElementOfBasis = give_element_of_basis(tBasisOfObject(iCount),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
                tBasisNeighbour = give_neighbour_of_basis(tBasisOfObject(iCount),aElementList.Dim,aElementList.Polynomial,tBuffer,aElementList.NumElements);
                for(uint i = 0; i < tCheckCrossBasis.length(); i++)
                {
                    if( tBasisNeighbour(tCheckCrossBasis(i)) < tActiveBasis.size() && tActiveBasis.test(tBasisNeighbour(tCheckCrossBasis(i))) == 1 ) // Check neighbours if they have the same phase and not checked yet
                    {
                        tBasisOfObject(iObj) = tBasisNeighbour(tCheckCrossBasis(i));
                        tActiveBasis.reset(tBasisNeighbour(tCheckCrossBasis(i)));
                        iObj++;
                    }
                }
                for(uint i = 0; i < tDiagBasis.length(); i++)
                {
                    if( tBasisNeighbour(tDiagBasis(i)) < tActiveBasis.size() && tActiveBasis.test(tBasisNeighbour(tDiagBasis(i))) == 1 && aElementList.ElementField(tElementOfBasis(i)) == icell ) // Check neighbours if they have the same phase and not checked yet and for the diagonal, check if the same colur is within the element
                    {
                        tBasisOfObject(iObj) = tBasisNeighbour(tDiagBasis(i));
                        tActiveBasis.reset(tBasisNeighbour(tDiagBasis(i)));
                        iObj++;
                    }
                }
                iCount++;
            }
            tBasisOfObject.resize(iObj,1); // Size the list of elements of one object
            if( tObjectsForPhase.size() < jObj+1) // If more objects then default, then resize!
                tObjectsForPhase.resize(jObj+1);

            tElementOfObject.set_size(tBasisOfObject.length()*pow(aElementList.Polynomial+1,aElementList.Dim),1,0); // Find active elements, which have support with the basis functions
            tVar = 0;
            for(uint i = 0; i < tBasisOfObject.length(); i++)
            {
                tElementOfBasis = give_element_of_basis(tBasisOfObject(i),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
                for(uint j = 0; j < tElementOfBasis.length(); j++)
                {
                    if( (aElementList.ElementActive).test(tElementOfBasis(j)) == 1 )
                    {
                        tElementOfObject(tVar) = tElementOfBasis(j);
                        tVar++;
                    }
                }
            }
            tElementOfObject.resize(tVar,1);
            tElementOfObject = unique( tElementOfObject );
            tObjectsForPhase(jObj) = tElementOfObject;
            jObj++;
        }
        tObjectsForPhase.resize(jObj);
        tObjects(icell) = tObjectsForPhase; // Save objects of a phase in a cell
        tNumberOfUsedCells = icell+1;
    }
    if( tNumberOfUsedCells < tObjects.size() )
        tObjects.resize(tNumberOfUsedCells);
    return tObjects;
}
Mat<real>
moris::Hierarchical_Mesh::RHS_for_L2_projection(
        uint & aDim,
        uint & aPolynomial,
        Mat<real> & aNodalField)
{
    Mat<real> Element_RHS(pow(aPolynomial+1,aDim),1);
    if( aDim == 2)
    {
        if( aPolynomial == 1)
        {
            Element_RHS(0) = 1.0/36.0*(4.0*aNodalField(0) + 2.0*aNodalField(1) + 2.0*aNodalField(2) + 1.0*aNodalField(3));
            Element_RHS(1) = 1.0/36.0*(2.0*aNodalField(0) + 4.0*aNodalField(1) + 1.0*aNodalField(2) + 2.0*aNodalField(3));
            Element_RHS(2) = 1.0/36.0*(2.0*aNodalField(0) + 1.0*aNodalField(1) + 4.0*aNodalField(2) + 2.0*aNodalField(3));
            Element_RHS(3) = 1.0/36.0*(1.0*aNodalField(0) + 2.0*aNodalField(1) + 2.0*aNodalField(2) + 4.0*aNodalField(3));
        }
        else if( aPolynomial == 2)
        {
            Element_RHS(0) = 0.002500000000000*aNodalField(0) + 0.005416666666667*aNodalField(1) + 0.000416666666667*aNodalField(2) + 0.005416666666667*aNodalField(3) + 0.011736111111111*aNodalField(4) + 0.000902777777778*aNodalField(5) + 0.000416666666667*aNodalField(6) + 0.000902777777778*aNodalField(7) + 0.000069444444444*aNodalField(8);
            Element_RHS(1) = 0.005416666666667*aNodalField(0) + 0.022500000000000*aNodalField(1) + 0.005416666666667*aNodalField(2) + 0.011736111111111*aNodalField(3) + 0.048750000000000*aNodalField(4) + 0.011736111111111*aNodalField(5) + 0.000902777777778*aNodalField(6) + 0.003750000000000*aNodalField(7) + 0.000902777777778*aNodalField(8);
            Element_RHS(2) = 0.000416666666667*aNodalField(0) + 0.005416666666667*aNodalField(1) + 0.002500000000000*aNodalField(2) + 0.000902777777778*aNodalField(3) + 0.011736111111111*aNodalField(4) + 0.005416666666667*aNodalField(5) + 0.000069444444444*aNodalField(6) + 0.000902777777778*aNodalField(7) + 0.000416666666667*aNodalField(8);
            Element_RHS(3) = 0.005416666666667*aNodalField(0) + 0.011736111111111*aNodalField(1) + 0.000902777777778*aNodalField(2) + 0.022500000000000*aNodalField(3) + 0.048750000000000*aNodalField(4) + 0.003750000000000*aNodalField(5) + 0.005416666666667*aNodalField(6) + 0.011736111111111*aNodalField(7) + 0.000902777777778*aNodalField(8);
            Element_RHS(4) = 0.011736111111111*aNodalField(0) + 0.048750000000000*aNodalField(1) + 0.011736111111111*aNodalField(2) + 0.048750000000000*aNodalField(3) + 0.202500000000000*aNodalField(4) + 0.048750000000000*aNodalField(5) + 0.011736111111111*aNodalField(6) + 0.048750000000000*aNodalField(7) + 0.011736111111111*aNodalField(8);
            Element_RHS(5) = 0.000902777777778*aNodalField(0) + 0.011736111111111*aNodalField(1) + 0.005416666666667*aNodalField(2) + 0.003750000000000*aNodalField(3) + 0.048750000000000*aNodalField(4) + 0.022500000000000*aNodalField(5) + 0.000902777777778*aNodalField(6) + 0.011736111111111*aNodalField(7) + 0.005416666666667*aNodalField(8);
            Element_RHS(6) = 0.000416666666667*aNodalField(0) + 0.000902777777778*aNodalField(1) + 0.000069444444444*aNodalField(2) + 0.005416666666667*aNodalField(3) + 0.011736111111111*aNodalField(4) + 0.000902777777778*aNodalField(5) + 0.002500000000000*aNodalField(6) + 0.005416666666667*aNodalField(7) + 0.000416666666667*aNodalField(8);
            Element_RHS(7) = 0.000902777777778*aNodalField(0) + 0.003750000000000*aNodalField(1) + 0.000902777777778*aNodalField(2) + 0.011736111111111*aNodalField(3) + 0.048750000000000*aNodalField(4) + 0.011736111111111*aNodalField(5) + 0.005416666666667*aNodalField(6) + 0.022500000000000*aNodalField(7) + 0.005416666666667*aNodalField(8);
            Element_RHS(8) = 0.000069444444444*aNodalField(0) + 0.000902777777778*aNodalField(1) + 0.000416666666667*aNodalField(2) + 0.000902777777778*aNodalField(3) + 0.011736111111111*aNodalField(4) + 0.005416666666667*aNodalField(5) + 0.000416666666667*aNodalField(6) + 0.005416666666667*aNodalField(7) + 0.002500000000000*aNodalField(8);
        }
    }
    else if( aDim == 3)
    {
        if( aPolynomial == 1)
        {
            Element_RHS(0) = 1.0/216.0*(8.0*aNodalField(0) + 4.0*aNodalField(1) + 4.0*aNodalField(2) + 2.0*aNodalField(3) + 4.0*aNodalField(4) + 2.0*aNodalField(5) + 2.0*aNodalField(6) + 1.0*aNodalField(7));
            Element_RHS(1) = 1.0/216.0*(4.0*aNodalField(0) + 8.0*aNodalField(1) + 2.0*aNodalField(2) + 4.0*aNodalField(3) + 2.0*aNodalField(4) + 4.0*aNodalField(5) + 1.0*aNodalField(6) + 2.0*aNodalField(7));
            Element_RHS(2) = 1.0/216.0*(4.0*aNodalField(0) + 2.0*aNodalField(1) + 8.0*aNodalField(2) + 4.0*aNodalField(3) + 2.0*aNodalField(4) + 1.0*aNodalField(5) + 4.0*aNodalField(6) + 2.0*aNodalField(7));
            Element_RHS(3) = 1.0/216.0*(2.0*aNodalField(0) + 4.0*aNodalField(1) + 4.0*aNodalField(2) + 8.0*aNodalField(3) + 1.0*aNodalField(4) + 2.0*aNodalField(5) + 2.0*aNodalField(6) + 4.0*aNodalField(7));
            Element_RHS(4) = 1.0/216.0*(4.0*aNodalField(0) + 2.0*aNodalField(1) + 2.0*aNodalField(2) + 1.0*aNodalField(3) + 8.0*aNodalField(4) + 4.0*aNodalField(5) + 4.0*aNodalField(6) + 2.0*aNodalField(7));
            Element_RHS(5) = 1.0/216.0*(2.0*aNodalField(0) + 4.0*aNodalField(1) + 1.0*aNodalField(2) + 2.0*aNodalField(3) + 4.0*aNodalField(4) + 8.0*aNodalField(5) + 2.0*aNodalField(6) + 4.0*aNodalField(7));
            Element_RHS(6) = 1.0/216.0*(2.0*aNodalField(0) + 1.0*aNodalField(1) + 4.0*aNodalField(2) + 2.0*aNodalField(3) + 4.0*aNodalField(4) + 2.0*aNodalField(5) + 8.0*aNodalField(6) + 4.0*aNodalField(7));
            Element_RHS(7) = 1.0/216.0*(1.0*aNodalField(0) + 2.0*aNodalField(1) + 2.0*aNodalField(2) + 4.0*aNodalField(3) + 2.0*aNodalField(4) + 4.0*aNodalField(5) + 4.0*aNodalField(6) + 8.0*aNodalField(7));
        }
        else if( aPolynomial == 2)
        {
            Element_RHS(0) = 1.250000e-04*aNodalField(0) + 2.708333e-04*aNodalField(1) + 2.083333e-05*aNodalField(2) + 2.708333e-04*aNodalField(3) + 5.868056e-04*aNodalField(4) + 4.513889e-05*aNodalField(5) + 2.083333e-05*aNodalField(6) + 4.513889e-05*aNodalField(7) + 3.472222e-06*aNodalField(8) + 2.708333e-04*aNodalField(9) + 5.868056e-04*aNodalField(10) + 4.513889e-05*aNodalField(11) + 5.868056e-04*aNodalField(12) + 1.271412e-03*aNodalField(13) + 9.780093e-05*aNodalField(14) + 4.513889e-05*aNodalField(15) + 9.780093e-05*aNodalField(16) + 7.523148e-06*aNodalField(17) + 2.083333e-05*aNodalField(18) + 4.513889e-05*aNodalField(19) + 3.472222e-06*aNodalField(20) + 4.513889e-05*aNodalField(21) + 9.780093e-05*aNodalField(22) + 7.523148e-06*aNodalField(23) + 3.472222e-06*aNodalField(24) + 7.523148e-06*aNodalField(25) + 5.787037e-07*aNodalField(26);
            Element_RHS(1) = 2.708333e-04*aNodalField(0) + 1.125000e-03*aNodalField(1) + 2.708333e-04*aNodalField(2) + 5.868056e-04*aNodalField(3) + 2.437500e-03*aNodalField(4) + 5.868056e-04*aNodalField(5) + 4.513889e-05*aNodalField(6) + 1.875000e-04*aNodalField(7) + 4.513889e-05*aNodalField(8) + 5.868056e-04*aNodalField(9) + 2.437500e-03*aNodalField(10) + 5.868056e-04*aNodalField(11) + 1.271412e-03*aNodalField(12) + 5.281250e-03*aNodalField(13) + 1.271412e-03*aNodalField(14) + 9.780093e-05*aNodalField(15) + 4.062500e-04*aNodalField(16) + 9.780093e-05*aNodalField(17) + 4.513889e-05*aNodalField(18) + 1.875000e-04*aNodalField(19) + 4.513889e-05*aNodalField(20) + 9.780093e-05*aNodalField(21) + 4.062500e-04*aNodalField(22) + 9.780093e-05*aNodalField(23) + 7.523148e-06*aNodalField(24) + 3.125000e-05*aNodalField(25) + 7.523148e-06*aNodalField(26);
            Element_RHS(2) = 2.083333e-05*aNodalField(0) + 2.708333e-04*aNodalField(1) + 1.250000e-04*aNodalField(2) + 4.513889e-05*aNodalField(3) + 5.868056e-04*aNodalField(4) + 2.708333e-04*aNodalField(5) + 3.472222e-06*aNodalField(6) + 4.513889e-05*aNodalField(7) + 2.083333e-05*aNodalField(8) + 4.513889e-05*aNodalField(9) + 5.868056e-04*aNodalField(10) + 2.708333e-04*aNodalField(11) + 9.780093e-05*aNodalField(12) + 1.271412e-03*aNodalField(13) + 5.868056e-04*aNodalField(14) + 7.523148e-06*aNodalField(15) + 9.780093e-05*aNodalField(16) + 4.513889e-05*aNodalField(17) + 3.472222e-06*aNodalField(18) + 4.513889e-05*aNodalField(19) + 2.083333e-05*aNodalField(20) + 7.523148e-06*aNodalField(21) + 9.780093e-05*aNodalField(22) + 4.513889e-05*aNodalField(23) + 5.787037e-07*aNodalField(24) + 7.523148e-06*aNodalField(25) + 3.472222e-06*aNodalField(26);
            Element_RHS(3) = 2.708333e-04*aNodalField(0) + 5.868056e-04*aNodalField(1) + 4.513889e-05*aNodalField(2) + 1.125000e-03*aNodalField(3) + 2.437500e-03*aNodalField(4) + 1.875000e-04*aNodalField(5) + 2.708333e-04*aNodalField(6) + 5.868056e-04*aNodalField(7) + 4.513889e-05*aNodalField(8) + 5.868056e-04*aNodalField(9) + 1.271412e-03*aNodalField(10) + 9.780093e-05*aNodalField(11) + 2.437500e-03*aNodalField(12) + 5.281250e-03*aNodalField(13) + 4.062500e-04*aNodalField(14) + 5.868056e-04*aNodalField(15) + 1.271412e-03*aNodalField(16) + 9.780093e-05*aNodalField(17) + 4.513889e-05*aNodalField(18) + 9.780093e-05*aNodalField(19) + 7.523148e-06*aNodalField(20) + 1.875000e-04*aNodalField(21) + 4.062500e-04*aNodalField(22) + 3.125000e-05*aNodalField(23) + 4.513889e-05*aNodalField(24) + 9.780093e-05*aNodalField(25) + 7.523148e-06*aNodalField(26);
            Element_RHS(4) = 5.868056e-04*aNodalField(0) + 2.437500e-03*aNodalField(1) + 5.868056e-04*aNodalField(2) + 2.437500e-03*aNodalField(3) + 1.012500e-02*aNodalField(4) + 2.437500e-03*aNodalField(5) + 5.868056e-04*aNodalField(6) + 2.437500e-03*aNodalField(7) + 5.868056e-04*aNodalField(8) + 1.271412e-03*aNodalField(9) + 5.281250e-03*aNodalField(10) + 1.271412e-03*aNodalField(11) + 5.281250e-03*aNodalField(12) + 2.193750e-02*aNodalField(13) + 5.281250e-03*aNodalField(14) + 1.271412e-03*aNodalField(15) + 5.281250e-03*aNodalField(16) + 1.271412e-03*aNodalField(17) + 9.780093e-05*aNodalField(18) + 4.062500e-04*aNodalField(19) + 9.780093e-05*aNodalField(20) + 4.062500e-04*aNodalField(21) + 1.687500e-03*aNodalField(22) + 4.062500e-04*aNodalField(23) + 9.780093e-05*aNodalField(24) + 4.062500e-04*aNodalField(25) + 9.780093e-05*aNodalField(26);
            Element_RHS(5) = 4.513889e-05*aNodalField(0) + 5.868056e-04*aNodalField(1) + 2.708333e-04*aNodalField(2) + 1.875000e-04*aNodalField(3) + 2.437500e-03*aNodalField(4) + 1.125000e-03*aNodalField(5) + 4.513889e-05*aNodalField(6) + 5.868056e-04*aNodalField(7) + 2.708333e-04*aNodalField(8) + 9.780093e-05*aNodalField(9) + 1.271412e-03*aNodalField(10) + 5.868056e-04*aNodalField(11) + 4.062500e-04*aNodalField(12) + 5.281250e-03*aNodalField(13) + 2.437500e-03*aNodalField(14) + 9.780093e-05*aNodalField(15) + 1.271412e-03*aNodalField(16) + 5.868056e-04*aNodalField(17) + 7.523148e-06*aNodalField(18) + 9.780093e-05*aNodalField(19) + 4.513889e-05*aNodalField(20) + 3.125000e-05*aNodalField(21) + 4.062500e-04*aNodalField(22) + 1.875000e-04*aNodalField(23) + 7.523148e-06*aNodalField(24) + 9.780093e-05*aNodalField(25) + 4.513889e-05*aNodalField(26);
            Element_RHS(6) = 2.083333e-05*aNodalField(0) + 4.513889e-05*aNodalField(1) + 3.472222e-06*aNodalField(2) + 2.708333e-04*aNodalField(3) + 5.868056e-04*aNodalField(4) + 4.513889e-05*aNodalField(5) + 1.250000e-04*aNodalField(6) + 2.708333e-04*aNodalField(7) + 2.083333e-05*aNodalField(8) + 4.513889e-05*aNodalField(9) + 9.780093e-05*aNodalField(10) + 7.523148e-06*aNodalField(11) + 5.868056e-04*aNodalField(12) + 1.271412e-03*aNodalField(13) + 9.780093e-05*aNodalField(14) + 2.708333e-04*aNodalField(15) + 5.868056e-04*aNodalField(16) + 4.513889e-05*aNodalField(17) + 3.472222e-06*aNodalField(18) + 7.523148e-06*aNodalField(19) + 5.787037e-07*aNodalField(20) + 4.513889e-05*aNodalField(21) + 9.780093e-05*aNodalField(22) + 7.523148e-06*aNodalField(23) + 2.083333e-05*aNodalField(24) + 4.513889e-05*aNodalField(25) + 3.472222e-06*aNodalField(26);
            Element_RHS(7) = 4.513889e-05*aNodalField(0) + 1.875000e-04*aNodalField(1) + 4.513889e-05*aNodalField(2) + 5.868056e-04*aNodalField(3) + 2.437500e-03*aNodalField(4) + 5.868056e-04*aNodalField(5) + 2.708333e-04*aNodalField(6) + 1.125000e-03*aNodalField(7) + 2.708333e-04*aNodalField(8) + 9.780093e-05*aNodalField(9) + 4.062500e-04*aNodalField(10) + 9.780093e-05*aNodalField(11) + 1.271412e-03*aNodalField(12) + 5.281250e-03*aNodalField(13) + 1.271412e-03*aNodalField(14) + 5.868056e-04*aNodalField(15) + 2.437500e-03*aNodalField(16) + 5.868056e-04*aNodalField(17) + 7.523148e-06*aNodalField(18) + 3.125000e-05*aNodalField(19) + 7.523148e-06*aNodalField(20) + 9.780093e-05*aNodalField(21) + 4.062500e-04*aNodalField(22) + 9.780093e-05*aNodalField(23) + 4.513889e-05*aNodalField(24) + 1.875000e-04*aNodalField(25) + 4.513889e-05*aNodalField(26);
            Element_RHS(8) = 3.472222e-06*aNodalField(0) + 4.513889e-05*aNodalField(1) + 2.083333e-05*aNodalField(2) + 4.513889e-05*aNodalField(3) + 5.868056e-04*aNodalField(4) + 2.708333e-04*aNodalField(5) + 2.083333e-05*aNodalField(6) + 2.708333e-04*aNodalField(7) + 1.250000e-04*aNodalField(8) + 7.523148e-06*aNodalField(9) + 9.780093e-05*aNodalField(10) + 4.513889e-05*aNodalField(11) + 9.780093e-05*aNodalField(12) + 1.271412e-03*aNodalField(13) + 5.868056e-04*aNodalField(14) + 4.513889e-05*aNodalField(15) + 5.868056e-04*aNodalField(16) + 2.708333e-04*aNodalField(17) + 5.787037e-07*aNodalField(18) + 7.523148e-06*aNodalField(19) + 3.472222e-06*aNodalField(20) + 7.523148e-06*aNodalField(21) + 9.780093e-05*aNodalField(22) + 4.513889e-05*aNodalField(23) + 3.472222e-06*aNodalField(24) + 4.513889e-05*aNodalField(25) + 2.083333e-05*aNodalField(26);
            Element_RHS(9) = 2.708333e-04*aNodalField(0) + 5.868056e-04*aNodalField(1) + 4.513889e-05*aNodalField(2) + 5.868056e-04*aNodalField(3) + 1.271412e-03*aNodalField(4) + 9.780093e-05*aNodalField(5) + 4.513889e-05*aNodalField(6) + 9.780093e-05*aNodalField(7) + 7.523148e-06*aNodalField(8) + 1.125000e-03*aNodalField(9) + 2.437500e-03*aNodalField(10) + 1.875000e-04*aNodalField(11) + 2.437500e-03*aNodalField(12) + 5.281250e-03*aNodalField(13) + 4.062500e-04*aNodalField(14) + 1.875000e-04*aNodalField(15) + 4.062500e-04*aNodalField(16) + 3.125000e-05*aNodalField(17) + 2.708333e-04*aNodalField(18) + 5.868056e-04*aNodalField(19) + 4.513889e-05*aNodalField(20) + 5.868056e-04*aNodalField(21) + 1.271412e-03*aNodalField(22) + 9.780093e-05*aNodalField(23) + 4.513889e-05*aNodalField(24) + 9.780093e-05*aNodalField(25) + 7.523148e-06*aNodalField(26);
            Element_RHS(10) = 5.868056e-04*aNodalField(0) + 2.437500e-03*aNodalField(1) + 5.868056e-04*aNodalField(2) + 1.271412e-03*aNodalField(3) + 5.281250e-03*aNodalField(4) + 1.271412e-03*aNodalField(5) + 9.780093e-05*aNodalField(6) + 4.062500e-04*aNodalField(7) + 9.780093e-05*aNodalField(8) + 2.437500e-03*aNodalField(9) + 1.012500e-02*aNodalField(10) + 2.437500e-03*aNodalField(11) + 5.281250e-03*aNodalField(12) + 2.193750e-02*aNodalField(13) + 5.281250e-03*aNodalField(14) + 4.062500e-04*aNodalField(15) + 1.687500e-03*aNodalField(16) + 4.062500e-04*aNodalField(17) + 5.868056e-04*aNodalField(18) + 2.437500e-03*aNodalField(19) + 5.868056e-04*aNodalField(20) + 1.271412e-03*aNodalField(21) + 5.281250e-03*aNodalField(22) + 1.271412e-03*aNodalField(23) + 9.780093e-05*aNodalField(24) + 4.062500e-04*aNodalField(25) + 9.780093e-05*aNodalField(26);
            Element_RHS(11) = 4.513889e-05*aNodalField(0) + 5.868056e-04*aNodalField(1) + 2.708333e-04*aNodalField(2) + 9.780093e-05*aNodalField(3) + 1.271412e-03*aNodalField(4) + 5.868056e-04*aNodalField(5) + 7.523148e-06*aNodalField(6) + 9.780093e-05*aNodalField(7) + 4.513889e-05*aNodalField(8) + 1.875000e-04*aNodalField(9) + 2.437500e-03*aNodalField(10) + 1.125000e-03*aNodalField(11) + 4.062500e-04*aNodalField(12) + 5.281250e-03*aNodalField(13) + 2.437500e-03*aNodalField(14) + 3.125000e-05*aNodalField(15) + 4.062500e-04*aNodalField(16) + 1.875000e-04*aNodalField(17) + 4.513889e-05*aNodalField(18) + 5.868056e-04*aNodalField(19) + 2.708333e-04*aNodalField(20) + 9.780093e-05*aNodalField(21) + 1.271412e-03*aNodalField(22) + 5.868056e-04*aNodalField(23) + 7.523148e-06*aNodalField(24) + 9.780093e-05*aNodalField(25) + 4.513889e-05*aNodalField(26);
            Element_RHS(12) = 5.868056e-04*aNodalField(0) + 1.271412e-03*aNodalField(1) + 9.780093e-05*aNodalField(2) + 2.437500e-03*aNodalField(3) + 5.281250e-03*aNodalField(4) + 4.062500e-04*aNodalField(5) + 5.868056e-04*aNodalField(6) + 1.271412e-03*aNodalField(7) + 9.780093e-05*aNodalField(8) + 2.437500e-03*aNodalField(9) + 5.281250e-03*aNodalField(10) + 4.062500e-04*aNodalField(11) + 1.012500e-02*aNodalField(12) + 2.193750e-02*aNodalField(13) + 1.687500e-03*aNodalField(14) + 2.437500e-03*aNodalField(15) + 5.281250e-03*aNodalField(16) + 4.062500e-04*aNodalField(17) + 5.868056e-04*aNodalField(18) + 1.271412e-03*aNodalField(19) + 9.780093e-05*aNodalField(20) + 2.437500e-03*aNodalField(21) + 5.281250e-03*aNodalField(22) + 4.062500e-04*aNodalField(23) + 5.868056e-04*aNodalField(24) + 1.271412e-03*aNodalField(25) + 9.780093e-05*aNodalField(26);
            Element_RHS(13) = 1.271412e-03*aNodalField(0) + 5.281250e-03*aNodalField(1) + 1.271412e-03*aNodalField(2) + 5.281250e-03*aNodalField(3) + 2.193750e-02*aNodalField(4) + 5.281250e-03*aNodalField(5) + 1.271412e-03*aNodalField(6) + 5.281250e-03*aNodalField(7) + 1.271412e-03*aNodalField(8) + 5.281250e-03*aNodalField(9) + 2.193750e-02*aNodalField(10) + 5.281250e-03*aNodalField(11) + 2.193750e-02*aNodalField(12) + 9.112500e-02*aNodalField(13) + 2.193750e-02*aNodalField(14) + 5.281250e-03*aNodalField(15) + 2.193750e-02*aNodalField(16) + 5.281250e-03*aNodalField(17) + 1.271412e-03*aNodalField(18) + 5.281250e-03*aNodalField(19) + 1.271412e-03*aNodalField(20) + 5.281250e-03*aNodalField(21) + 2.193750e-02*aNodalField(22) + 5.281250e-03*aNodalField(23) + 1.271412e-03*aNodalField(24) + 5.281250e-03*aNodalField(25) + 1.271412e-03*aNodalField(26);
            Element_RHS(14) = 9.780093e-05*aNodalField(0) + 1.271412e-03*aNodalField(1) + 5.868056e-04*aNodalField(2) + 4.062500e-04*aNodalField(3) + 5.281250e-03*aNodalField(4) + 2.437500e-03*aNodalField(5) + 9.780093e-05*aNodalField(6) + 1.271412e-03*aNodalField(7) + 5.868056e-04*aNodalField(8) + 4.062500e-04*aNodalField(9) + 5.281250e-03*aNodalField(10) + 2.437500e-03*aNodalField(11) + 1.687500e-03*aNodalField(12) + 2.193750e-02*aNodalField(13) + 1.012500e-02*aNodalField(14) + 4.062500e-04*aNodalField(15) + 5.281250e-03*aNodalField(16) + 2.437500e-03*aNodalField(17) + 9.780093e-05*aNodalField(18) + 1.271412e-03*aNodalField(19) + 5.868056e-04*aNodalField(20) + 4.062500e-04*aNodalField(21) + 5.281250e-03*aNodalField(22) + 2.437500e-03*aNodalField(23) + 9.780093e-05*aNodalField(24) + 1.271412e-03*aNodalField(25) + 5.868056e-04*aNodalField(26);
            Element_RHS(15) = 4.513889e-05*aNodalField(0) + 9.780093e-05*aNodalField(1) + 7.523148e-06*aNodalField(2) + 5.868056e-04*aNodalField(3) + 1.271412e-03*aNodalField(4) + 9.780093e-05*aNodalField(5) + 2.708333e-04*aNodalField(6) + 5.868056e-04*aNodalField(7) + 4.513889e-05*aNodalField(8) + 1.875000e-04*aNodalField(9) + 4.062500e-04*aNodalField(10) + 3.125000e-05*aNodalField(11) + 2.437500e-03*aNodalField(12) + 5.281250e-03*aNodalField(13) + 4.062500e-04*aNodalField(14) + 1.125000e-03*aNodalField(15) + 2.437500e-03*aNodalField(16) + 1.875000e-04*aNodalField(17) + 4.513889e-05*aNodalField(18) + 9.780093e-05*aNodalField(19) + 7.523148e-06*aNodalField(20) + 5.868056e-04*aNodalField(21) + 1.271412e-03*aNodalField(22) + 9.780093e-05*aNodalField(23) + 2.708333e-04*aNodalField(24) + 5.868056e-04*aNodalField(25) + 4.513889e-05*aNodalField(26);
            Element_RHS(16) = 9.780093e-05*aNodalField(0) + 4.062500e-04*aNodalField(1) + 9.780093e-05*aNodalField(2) + 1.271412e-03*aNodalField(3) + 5.281250e-03*aNodalField(4) + 1.271412e-03*aNodalField(5) + 5.868056e-04*aNodalField(6) + 2.437500e-03*aNodalField(7) + 5.868056e-04*aNodalField(8) + 4.062500e-04*aNodalField(9) + 1.687500e-03*aNodalField(10) + 4.062500e-04*aNodalField(11) + 5.281250e-03*aNodalField(12) + 2.193750e-02*aNodalField(13) + 5.281250e-03*aNodalField(14) + 2.437500e-03*aNodalField(15) + 1.012500e-02*aNodalField(16) + 2.437500e-03*aNodalField(17) + 9.780093e-05*aNodalField(18) + 4.062500e-04*aNodalField(19) + 9.780093e-05*aNodalField(20) + 1.271412e-03*aNodalField(21) + 5.281250e-03*aNodalField(22) + 1.271412e-03*aNodalField(23) + 5.868056e-04*aNodalField(24) + 2.437500e-03*aNodalField(25) + 5.868056e-04*aNodalField(26);
            Element_RHS(17) = 7.523148e-06*aNodalField(0) + 9.780093e-05*aNodalField(1) + 4.513889e-05*aNodalField(2) + 9.780093e-05*aNodalField(3) + 1.271412e-03*aNodalField(4) + 5.868056e-04*aNodalField(5) + 4.513889e-05*aNodalField(6) + 5.868056e-04*aNodalField(7) + 2.708333e-04*aNodalField(8) + 3.125000e-05*aNodalField(9) + 4.062500e-04*aNodalField(10) + 1.875000e-04*aNodalField(11) + 4.062500e-04*aNodalField(12) + 5.281250e-03*aNodalField(13) + 2.437500e-03*aNodalField(14) + 1.875000e-04*aNodalField(15) + 2.437500e-03*aNodalField(16) + 1.125000e-03*aNodalField(17) + 7.523148e-06*aNodalField(18) + 9.780093e-05*aNodalField(19) + 4.513889e-05*aNodalField(20) + 9.780093e-05*aNodalField(21) + 1.271412e-03*aNodalField(22) + 5.868056e-04*aNodalField(23) + 4.513889e-05*aNodalField(24) + 5.868056e-04*aNodalField(25) + 2.708333e-04*aNodalField(26);
            Element_RHS(18) = 2.083333e-05*aNodalField(0) + 4.513889e-05*aNodalField(1) + 3.472222e-06*aNodalField(2) + 4.513889e-05*aNodalField(3) + 9.780093e-05*aNodalField(4) + 7.523148e-06*aNodalField(5) + 3.472222e-06*aNodalField(6) + 7.523148e-06*aNodalField(7) + 5.787037e-07*aNodalField(8) + 2.708333e-04*aNodalField(9) + 5.868056e-04*aNodalField(10) + 4.513889e-05*aNodalField(11) + 5.868056e-04*aNodalField(12) + 1.271412e-03*aNodalField(13) + 9.780093e-05*aNodalField(14) + 4.513889e-05*aNodalField(15) + 9.780093e-05*aNodalField(16) + 7.523148e-06*aNodalField(17) + 1.250000e-04*aNodalField(18) + 2.708333e-04*aNodalField(19) + 2.083333e-05*aNodalField(20) + 2.708333e-04*aNodalField(21) + 5.868056e-04*aNodalField(22) + 4.513889e-05*aNodalField(23) + 2.083333e-05*aNodalField(24) + 4.513889e-05*aNodalField(25) + 3.472222e-06*aNodalField(26);
            Element_RHS(19) = 4.513889e-05*aNodalField(0) + 1.875000e-04*aNodalField(1) + 4.513889e-05*aNodalField(2) + 9.780093e-05*aNodalField(3) + 4.062500e-04*aNodalField(4) + 9.780093e-05*aNodalField(5) + 7.523148e-06*aNodalField(6) + 3.125000e-05*aNodalField(7) + 7.523148e-06*aNodalField(8) + 5.868056e-04*aNodalField(9) + 2.437500e-03*aNodalField(10) + 5.868056e-04*aNodalField(11) + 1.271412e-03*aNodalField(12) + 5.281250e-03*aNodalField(13) + 1.271412e-03*aNodalField(14) + 9.780093e-05*aNodalField(15) + 4.062500e-04*aNodalField(16) + 9.780093e-05*aNodalField(17) + 2.708333e-04*aNodalField(18) + 1.125000e-03*aNodalField(19) + 2.708333e-04*aNodalField(20) + 5.868056e-04*aNodalField(21) + 2.437500e-03*aNodalField(22) + 5.868056e-04*aNodalField(23) + 4.513889e-05*aNodalField(24) + 1.875000e-04*aNodalField(25) + 4.513889e-05*aNodalField(26);
            Element_RHS(20) = 3.472222e-06*aNodalField(0) + 4.513889e-05*aNodalField(1) + 2.083333e-05*aNodalField(2) + 7.523148e-06*aNodalField(3) + 9.780093e-05*aNodalField(4) + 4.513889e-05*aNodalField(5) + 5.787037e-07*aNodalField(6) + 7.523148e-06*aNodalField(7) + 3.472222e-06*aNodalField(8) + 4.513889e-05*aNodalField(9) + 5.868056e-04*aNodalField(10) + 2.708333e-04*aNodalField(11) + 9.780093e-05*aNodalField(12) + 1.271412e-03*aNodalField(13) + 5.868056e-04*aNodalField(14) + 7.523148e-06*aNodalField(15) + 9.780093e-05*aNodalField(16) + 4.513889e-05*aNodalField(17) + 2.083333e-05*aNodalField(18) + 2.708333e-04*aNodalField(19) + 1.250000e-04*aNodalField(20) + 4.513889e-05*aNodalField(21) + 5.868056e-04*aNodalField(22) + 2.708333e-04*aNodalField(23) + 3.472222e-06*aNodalField(24) + 4.513889e-05*aNodalField(25) + 2.083333e-05*aNodalField(26);
            Element_RHS(21) = 4.513889e-05*aNodalField(0) + 9.780093e-05*aNodalField(1) + 7.523148e-06*aNodalField(2) + 1.875000e-04*aNodalField(3) + 4.062500e-04*aNodalField(4) + 3.125000e-05*aNodalField(5) + 4.513889e-05*aNodalField(6) + 9.780093e-05*aNodalField(7) + 7.523148e-06*aNodalField(8) + 5.868056e-04*aNodalField(9) + 1.271412e-03*aNodalField(10) + 9.780093e-05*aNodalField(11) + 2.437500e-03*aNodalField(12) + 5.281250e-03*aNodalField(13) + 4.062500e-04*aNodalField(14) + 5.868056e-04*aNodalField(15) + 1.271412e-03*aNodalField(16) + 9.780093e-05*aNodalField(17) + 2.708333e-04*aNodalField(18) + 5.868056e-04*aNodalField(19) + 4.513889e-05*aNodalField(20) + 1.125000e-03*aNodalField(21) + 2.437500e-03*aNodalField(22) + 1.875000e-04*aNodalField(23) + 2.708333e-04*aNodalField(24) + 5.868056e-04*aNodalField(25) + 4.513889e-05*aNodalField(26);
            Element_RHS(22) = 9.780093e-05*aNodalField(0) + 4.062500e-04*aNodalField(1) + 9.780093e-05*aNodalField(2) + 4.062500e-04*aNodalField(3) + 1.687500e-03*aNodalField(4) + 4.062500e-04*aNodalField(5) + 9.780093e-05*aNodalField(6) + 4.062500e-04*aNodalField(7) + 9.780093e-05*aNodalField(8) + 1.271412e-03*aNodalField(9) + 5.281250e-03*aNodalField(10) + 1.271412e-03*aNodalField(11) + 5.281250e-03*aNodalField(12) + 2.193750e-02*aNodalField(13) + 5.281250e-03*aNodalField(14) + 1.271412e-03*aNodalField(15) + 5.281250e-03*aNodalField(16) + 1.271412e-03*aNodalField(17) + 5.868056e-04*aNodalField(18) + 2.437500e-03*aNodalField(19) + 5.868056e-04*aNodalField(20) + 2.437500e-03*aNodalField(21) + 1.012500e-02*aNodalField(22) + 2.437500e-03*aNodalField(23) + 5.868056e-04*aNodalField(24) + 2.437500e-03*aNodalField(25) + 5.868056e-04*aNodalField(26);
            Element_RHS(23) = 7.523148e-06*aNodalField(0) + 9.780093e-05*aNodalField(1) + 4.513889e-05*aNodalField(2) + 3.125000e-05*aNodalField(3) + 4.062500e-04*aNodalField(4) + 1.875000e-04*aNodalField(5) + 7.523148e-06*aNodalField(6) + 9.780093e-05*aNodalField(7) + 4.513889e-05*aNodalField(8) + 9.780093e-05*aNodalField(9) + 1.271412e-03*aNodalField(10) + 5.868056e-04*aNodalField(11) + 4.062500e-04*aNodalField(12) + 5.281250e-03*aNodalField(13) + 2.437500e-03*aNodalField(14) + 9.780093e-05*aNodalField(15) + 1.271412e-03*aNodalField(16) + 5.868056e-04*aNodalField(17) + 4.513889e-05*aNodalField(18) + 5.868056e-04*aNodalField(19) + 2.708333e-04*aNodalField(20) + 1.875000e-04*aNodalField(21) + 2.437500e-03*aNodalField(22) + 1.125000e-03*aNodalField(23) + 4.513889e-05*aNodalField(24) + 5.868056e-04*aNodalField(25) + 2.708333e-04*aNodalField(26);
            Element_RHS(24) = 3.472222e-06*aNodalField(0) + 7.523148e-06*aNodalField(1) + 5.787037e-07*aNodalField(2) + 4.513889e-05*aNodalField(3) + 9.780093e-05*aNodalField(4) + 7.523148e-06*aNodalField(5) + 2.083333e-05*aNodalField(6) + 4.513889e-05*aNodalField(7) + 3.472222e-06*aNodalField(8) + 4.513889e-05*aNodalField(9) + 9.780093e-05*aNodalField(10) + 7.523148e-06*aNodalField(11) + 5.868056e-04*aNodalField(12) + 1.271412e-03*aNodalField(13) + 9.780093e-05*aNodalField(14) + 2.708333e-04*aNodalField(15) + 5.868056e-04*aNodalField(16) + 4.513889e-05*aNodalField(17) + 2.083333e-05*aNodalField(18) + 4.513889e-05*aNodalField(19) + 3.472222e-06*aNodalField(20) + 2.708333e-04*aNodalField(21) + 5.868056e-04*aNodalField(22) + 4.513889e-05*aNodalField(23) + 1.250000e-04*aNodalField(24) + 2.708333e-04*aNodalField(25) + 2.083333e-05*aNodalField(26);
            Element_RHS(25) = 7.523148e-06*aNodalField(0) + 3.125000e-05*aNodalField(1) + 7.523148e-06*aNodalField(2) + 9.780093e-05*aNodalField(3) + 4.062500e-04*aNodalField(4) + 9.780093e-05*aNodalField(5) + 4.513889e-05*aNodalField(6) + 1.875000e-04*aNodalField(7) + 4.513889e-05*aNodalField(8) + 9.780093e-05*aNodalField(9) + 4.062500e-04*aNodalField(10) + 9.780093e-05*aNodalField(11) + 1.271412e-03*aNodalField(12) + 5.281250e-03*aNodalField(13) + 1.271412e-03*aNodalField(14) + 5.868056e-04*aNodalField(15) + 2.437500e-03*aNodalField(16) + 5.868056e-04*aNodalField(17) + 4.513889e-05*aNodalField(18) + 1.875000e-04*aNodalField(19) + 4.513889e-05*aNodalField(20) + 5.868056e-04*aNodalField(21) + 2.437500e-03*aNodalField(22) + 5.868056e-04*aNodalField(23) + 2.708333e-04*aNodalField(24) + 1.125000e-03*aNodalField(25) + 2.708333e-04*aNodalField(26);
            Element_RHS(26) = 5.787037e-07*aNodalField(0) + 7.523148e-06*aNodalField(1) + 3.472222e-06*aNodalField(2) + 7.523148e-06*aNodalField(3) + 9.780093e-05*aNodalField(4) + 4.513889e-05*aNodalField(5) + 3.472222e-06*aNodalField(6) + 4.513889e-05*aNodalField(7) + 2.083333e-05*aNodalField(8) + 7.523148e-06*aNodalField(9) + 9.780093e-05*aNodalField(10) + 4.513889e-05*aNodalField(11) + 9.780093e-05*aNodalField(12) + 1.271412e-03*aNodalField(13) + 5.868056e-04*aNodalField(14) + 4.513889e-05*aNodalField(15) + 5.868056e-04*aNodalField(16) + 2.708333e-04*aNodalField(17) + 3.472222e-06*aNodalField(18) + 4.513889e-05*aNodalField(19) + 2.083333e-05*aNodalField(20) + 4.513889e-05*aNodalField(21) + 5.868056e-04*aNodalField(22) + 2.708333e-04*aNodalField(23) + 2.083333e-05*aNodalField(24) + 2.708333e-04*aNodalField(25) + 1.250000e-04*aNodalField(26);

        }
    }
    return Element_RHS;
}
Mat<real>
moris::Hierarchical_Mesh::MassMatrix_for_L2_projection(
        uint & aDim,
        uint & aPolynomial)
{
    Mat<real> Element_Mass(pow(aPolynomial+1,aDim),pow(aPolynomial+1,aDim));
    if( aDim == 2)
    {
        if( aPolynomial == 1)
        {
            Element_Mass(0,0) = 1.0/9.0; Element_Mass(0,1) = 1.0/18.0; Element_Mass(0,2) = 1.0/18.0; Element_Mass(0,3) = 1.0/36.0;
            Element_Mass(1,0) = 1.0/18.0; Element_Mass(1,1) = 1.0/9.0; Element_Mass(1,2) = 1.0/36.0; Element_Mass(1,3) = 1.0/18.0;
            Element_Mass(2,0) = 1.0/18.0; Element_Mass(2,1) = 1.0/36.0; Element_Mass(2,2) = 1.0/9.0; Element_Mass(2,3) = 1.0/18.0;
            Element_Mass(3,0) = 1.0/36.0; Element_Mass(3,1) = 1.0/18.0; Element_Mass(3,2) = 1.0/18.0; Element_Mass(3,3) = 1.0/9.0;
        }
        else if( aPolynomial == 2)
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
    else if( aDim == 3)
    {
        if( aPolynomial == 1)
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
        else if( aPolynomial == 2)
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
void moris::Hierarchical_Mesh::Filter_for_smoothing(uint & aBasis,
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    Mat<uint> tBasisPosition = this->give_position_of_basis(aBasis,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
    Mat<real> tBasisCoordinate = this->give_coordinate_from_basis(aBasis,aElementList.PolynomialDesign,aElementList);
    uint taBasisLevel = this->give_basis_level(aBasis,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
    Mat<uint> tElementOfBasis = this->give_element_of_basis(aBasis,aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
    uint tWhichLevel = 0;
    uint tParent = UINT_MAX;
    uint tElementOfaBasis = UINT_MAX; // Save the Element in which the Basis function live.
    for(uint i = 0; i < tElementOfBasis.length(); i++)
    {
        if( (aElementList.ElementActiveDesign).test(tElementOfBasis(i)) == 1)
        {
            tParent = give_parent_of_level_x(tElementOfBasis(i),aElementList.Dim,aElementList.NumElements,tWhichLevel);
            tElementOfaBasis = tElementOfBasis(i);
            if( tParent == UINT_MAX)
                tParent = tElementOfBasis(i);
            break;
        }
    }
    uint tNumberOfElements  = this->give_number_of_elements(aElementList.Level,aElementList.Dim,aElementList.NumElements);
    BoostBitset Element_dummy(aElementList.NumberElements+1); // Dummy bitset to remember the tested elements in the filter radius
    Mat<uint> tPossibleElementList(tNumberOfElements+1,1,0);
    tPossibleElementList(0) = tParent; // Set first element for searching
    tPossibleElementList(1) = tElementOfaBasis; // Save also active element, which have support with the basis function "aBasis"
    Element_dummy.set(tParent);
    uint tl = 0;
    uint tVar = 2;
    Mat<uint> tElementPosition;
    uint tElementLevel, tPossibleElementListLevel;
    Mat<real> tElementMiddleCoordinate;
    real tDistance;
    Mat<uint> tChildren;
    // Search for elements and childrens, until no new one can found
    while( tPossibleElementList(tl) > 0)
    {
        //        std::cout << " tPossibleElementList(tl) " << tPossibleElementList(tl) << std::endl;
        uint tBuffer = 1; // Only one layer of neighbour elements are needed
        Mat<uint> tElement_neighbour = this->give_neighbour_of_element(tPossibleElementList(tl),aElementList.Dim,tBuffer,aElementList.NumElements);
        tPossibleElementListLevel = this->give_element_level(tPossibleElementList(tl),aElementList.Dim,aElementList.NumElements);
        for( uint i = 0; i < tElement_neighbour.length(); i++)
        {
            tElementLevel = this->give_element_level(tElement_neighbour(i),aElementList.Dim,aElementList.NumElements);
            if( tElement_neighbour(i) > 0 && tElement_neighbour(i) <= aElementList.NumberElements && Element_dummy.test(tElement_neighbour(i)) == 0 && tPossibleElementListLevel == tElementLevel)
            {
                //                std::cout << " tElement_neighbour(i) " << tElement_neighbour(i) << std::endl;
                //                tElementPosition = this->give_position_of_element(tElement_neighbour(i),aElementList.Dim,aElementList.NumElements);

                //                if( tElementPosition(0) >= aElementList.Polynomial*pow(2,tElementLevel) && tElementPosition(0) < (aElementList.NumElements(0)-aElementList.Polynomial)*pow(2,tElementLevel)
                //                        && tElementPosition(1) >= aElementList.Polynomial*pow(2,tElementLevel) && tElementPosition(1) < (aElementList.NumElements(1)-aElementList.Polynomial)*pow(2,tElementLevel) )
                //                {
                tElementMiddleCoordinate = give_middlecoordinate_from_element(tElement_neighbour(i),aElementList);
                //                tElementMiddleCoordinate.print("tElementMiddleCoordinate");
                //                tBasisCoordinate.print("tBasisCoordinate");
                //                tElement_neighbour.print("tElement_neighbour");
                tElementMiddleCoordinate = tBasisCoordinate - tElementMiddleCoordinate;
                tDistance = tElementMiddleCoordinate.norm();
                if( (aElementList.FilterRadius - tDistance) > 0.0 )
                {
                    tPossibleElementList(tVar) = tElement_neighbour(i);
                    Element_dummy.set(tElement_neighbour(i));
                    tVar++;
                }
                //                }
            }
        }
        if( tPossibleElementListLevel < aElementList.Level && (aElementList.ElementActiveDesign).test(tPossibleElementList(tl)) == 0)
        {
            tChildren = give_children_of_element(tPossibleElementList(tl),aElementList.Dim,aElementList.NumElements);
            for( uint i = 0; i < tChildren.length(); i++)
            {
                tElementMiddleCoordinate = give_middlecoordinate_from_element(tChildren(i),aElementList);
                tElementMiddleCoordinate = tBasisCoordinate - tElementMiddleCoordinate;
                tDistance = tElementMiddleCoordinate.norm();
                if( (aElementList.FilterRadius - tDistance) > 0.0 )
                {
                    tPossibleElementList(tVar) = tChildren(i);
                    Element_dummy.set(tChildren(i));
                    tVar++;
                }
            }
        }
        tl++;
    }
    tPossibleElementList.resize(tVar,1);
    tPossibleElementList = unique(tPossibleElementList);
    if( (aElementList.ElementActiveDesign).size() < (tPossibleElementList.max()+1)  )
        (aElementList.ElementActiveDesign).resize(tPossibleElementList.max()+1);
    Mat<uint> tPossibleBasis(tVar*pow(aElementList.PolynomialDesign+1,aElementList.Dim),1,0);
    Mat<uint> tBasis;
    uint tBasisLevel;
    Mat<real> tBasisCoordinates;
    tVar = 0;
    //Save all basis design functions of each element
    for(uint i = 0; i < tPossibleElementList.length(); i++)
    {
        if( (aElementList.ElementActiveDesign).test(tPossibleElementList(i)) == 1 )
        {
            tBasis = give_basis_of_element(tPossibleElementList(i),aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
            tPossibleBasis.rows(tVar*tBasis.length(),(tVar+1)*tBasis.length()-1) = tBasis.rows(0,tBasis.length()-1);
            tVar++;
        }
    }
    tPossibleBasis.resize(tVar*tBasis.length(),1);
    tPossibleBasis = unique(tPossibleBasis);
    if( (aBasisList.DesignBSplineActive).size() < (tPossibleBasis.max()+1)  )
        (aBasisList.DesignBSplineActive).resize(tPossibleBasis.max()+1);
    (aBasisList.BasisFilterList).set_size(tPossibleBasis.length(),3,0);
    tVar = 0;
    //Search for active basis design functions and check the radius
    for( uint i = 0; i < tPossibleBasis.length(); i++)
    {
        if( (aBasisList.DesignBSplineActive).test(tPossibleBasis(i)) == 1  )
        {
            //            std::cout << " tPossibleBasis(i) " << tPossibleBasis(i) << std::endl;
            tBasisLevel = give_basis_level(tPossibleBasis(i),aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
            tBasisCoordinates = give_coordinate_from_basis(tPossibleBasis(i),aElementList.PolynomialDesign,aElementList);
            //            tBasisCoordinates.print("tBasisCoordinates");
            tBasisCoordinates = tBasisCoordinate - tBasisCoordinates;
            tDistance = tBasisCoordinates.norm();
            if( (aElementList.FilterRadius - tDistance) > 0.0 )
            {
                (aBasisList.BasisFilterList)(tVar,0) = tPossibleBasis(i);
                (aBasisList.BasisFilterList)(tVar,1) = aElementList.FilterRadius - tDistance;
                (aBasisList.BasisFilterList)(tVar,2) = ((real)taBasisLevel+1.0)/((real)tBasisLevel+1.0);
                tVar++;
            }
        }
    }
    (aBasisList.BasisFilterList).resize(tVar,3);
}

void moris::Hierarchical_Mesh::Update_IDandTMatrix_design(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    Mat<uint> tOldIdList;
    Mat<real> tBasisFiler;
    Mat<real> tNewTMatrix;
    Mat<uint> tNewIdField;
    Cell<Mat<real>> tNewTMatrixCell((aElementList.IdFieldFieldDesign).n_rows());
    Cell<Mat<uint>> tNewIdFieldCell((aElementList.IdFieldFieldDesign).n_rows());
    real tWeight;
    Mat<uint> tFindUniqueList;
    Mat<real> tHelp;
    Mat<real> tHelpa, tHelpb;
    Mat<uint> tHelpc;
    Mat<uint> tUniqueMap;
    uint tCopy1 = 0, tCopy2 = 0; // Length for beginning and end is needed
    uint tMaxIDField = 0;
    for(uint i=0; i<(aElementList.NodalLocaltoGlobal).length(); i++)
    {
        //        if( i == 3070)
        //        std::cout << " (aElementList.NodalLocaltoGlobal) " << (aElementList.NodalLocaltoGlobal) << std::endl;
        tBasisFiler.set_size((aElementList.Level+1)*1000*aElementList.PolynomialDesign,3,0);
        tCopy1 = 0;
        tCopy2 = 0;
        tOldIdList.set_size((aElementList.IdFieldFieldDesign)(i,0),1,0);
        //        ((aElementList.IdFieldFieldDesign).row(i)).print(" ((aElementList.IdFieldFieldDesign).row(i))");
        for(uint j = 0; j< (aElementList.IdFieldFieldDesign)(i,0); j++)
        {
            tOldIdList(j) = (aElementList.IdFieldFieldDesign)(i,j+1);

            this->Filter_for_smoothing((aElementList.IdFieldFieldDesign)(i,j+1),aElementList,aBasisList);
            tWeight = 0.0;
            tHelpa = (aBasisList.BasisFilterList).col(1);
            if( aElementList.FilterLevelWeight == true)
            {
                tHelpb = (aBasisList.BasisFilterList).col(2);
            }
            else
            {
                tHelpb.set_size(tHelpa.length(),1,1);
            }
            tWeight = dot(tHelpa,tHelpb);
            (aBasisList.BasisFilterList).col(1) /= tWeight;
            (aBasisList.BasisFilterList).col(1) *= (aElementList.TMatrixFieldDesign)(i,j+1);
            tCopy2 += (aBasisList.BasisFilterList).n_rows();

            if( tCopy2 > tBasisFiler.n_rows() )
                tBasisFiler.resize((tCopy1+tCopy2),3);

            tBasisFiler.rows(tCopy1,tCopy2-1) = (aBasisList.BasisFilterList).rows(0,(aBasisList.BasisFilterList).n_rows()-1);
            tCopy1 += (aBasisList.BasisFilterList).n_rows();
        }
        tBasisFiler.resize(tCopy2,3);
        tHelpc.set_size(tBasisFiler.n_rows(),1);
        for(uint j = 0; j < tBasisFiler.n_rows(); j++)
            tHelpc(j) = tBasisFiler(j,0);
        tFindUniqueList = unique(tHelpc);
        tUniqueMap.set_size(tFindUniqueList.max()+1,1,0);
        for(uint j = 0; j < tFindUniqueList.length(); j++)
            tUniqueMap(tFindUniqueList(j)) = j+1;
        tNewTMatrix.set_size(1,tFindUniqueList.length()+1,0);
        tNewIdField.set_size(1,tFindUniqueList.length()+1,0);
        for(uint j = 0; j < tBasisFiler.n_rows(); j++)
        {
            tNewIdField(tUniqueMap(tBasisFiler(j,0))) = tBasisFiler(j,0);
            if( aElementList.FilterLevelWeight == true)
            {
                tNewTMatrix(tUniqueMap(tBasisFiler(j,0))) += tBasisFiler(j,1)*tBasisFiler(j,2);
            }
            else
            {
                tNewTMatrix(tUniqueMap(tBasisFiler(j,0))) += tBasisFiler(j,1);
            }
        }
        if( aElementList.PerformNormalization == true)
        {
            tWeight= sum(tNewTMatrix);
            tNewTMatrix = tNewTMatrix / tWeight;
        }
        tNewTMatrix(0) = tFindUniqueList.length();
        tNewIdField(0) = tFindUniqueList.length();
        if ( tMaxIDField < tNewIdField.length())
            tMaxIDField = tNewIdField.length();

        tNewTMatrixCell(i) = tNewTMatrix;
        tNewIdFieldCell(i) = tNewIdField;
    }
    (aElementList.IdFieldFieldDesign).set_size((aElementList.NodalLocaltoGlobal).length(),tMaxIDField,0);
    (aElementList.TMatrixFieldDesign).set_size((aElementList.NodalLocaltoGlobal).length(),tMaxIDField,0);
    for(uint i=0; i<(aElementList.NodalLocaltoGlobal).length(); i++)
    {
        tHelpc =  tNewIdFieldCell(i).row(0);
        tHelpc.resize(1,tMaxIDField);
        tHelpa =  tNewTMatrixCell(i).row(0);
        tHelpa.resize(1,tMaxIDField);
        (aElementList.IdFieldFieldDesign).row(i) = tHelpc.row(0);
        (aElementList.TMatrixFieldDesign).row(i) = tHelpa.row(0);
    }
}

void moris::Hierarchical_Mesh::L2_projection(
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    uint tLevel;
    Mat<uint> tBasis;
    Mat<real> tOldSol;
    Mat<real> tNewSol(pow(aElementList.PolynomialDesign+1,aElementList.Dim),1);
    Mat<real> tNewSolPDV(pow(aElementList.PolynomialDesign+1,aElementList.Dim),1);
    Mat<real> ElementRHS;
    Mat<real> ElementRHSDummy;
    Mat<real> ElementMass;
    Mat<real> tTmatrix;
    Mat<uint> tBasisOfElement;
    uint tVar = 0;

    if( (aBasisList.DesignBSplineActiveLastStep).size() < aBasisList.NumberBasis)
        (aBasisList.DesignBSplineActiveLastStep).resize(aBasisList.NumberBasis);

    Mat<real> GLBRHS((aBasisList.DesignBSplineListMap).max()+1, 1, 0);
    moris::Mat< moris::uint > RowInd( (aElementList.SparseSize), 1, 0 );
    moris::Mat< moris::uint > ColInd( (aElementList.SparseSize), 1, 0 );
    moris::Mat< moris::real > Values( (aElementList.SparseSize), 1, 0 );

    moris::tic tBuildSystemTiming;

    for(uint i = 0; i < (aElementList.ElementListActiveDesign).length(); i++)
    {
        //        if( (aElementList.ElementListActiveDesign)(i) == 165 || (aElementList.ElementListActiveDesign)(i) == 197 || (aElementList.ElementListActiveDesign)(i) == 198 || (aElementList.ElementListActiveDesign)(i) == 170 )
        //            std::cout << "(ElementList.ElementListActiveDesign)(i) " << (aElementList.ElementListActiveDesign)(i) << std::endl;

        //        std::cout << "(ElementList.ElementListActiveDesign)(i) " << (aElementList.ElementListActiveDesign)(i) << std::endl;
        //                                if( (aElementList.ElementListActiveDesign)(i) == 2436)
        //                                    std::cout << "(ElementList.ElemLocaltoGlobal)(i) " << (aElementList.ElementListActiveDesign)(i) << std::endl;
        //                                if( (aElementList.ElementListActiveDesign)(i) == 780)
        //                                    std::cout << "(ElementList.ElemLocaltoGlobal)(i) " << (aElementList.ElementListActiveDesign)(i) << std::endl;

        this->give_Tmatrix_and_id_field_for_designvariables_projection_coarsening((aElementList.ElementListActiveDesign)(i),aElementList,aBasisList);
        //                        (aElementList.IdField).print("(ElementList.IdField) old");
        //                        (aElementList.TMatrix).print("(ElementList.TMatrix) old");
        if( (aElementList.ListOfChildren).length() >= pow(2,aElementList.Dim) )
        {
            ElementRHS.set_size(pow(aElementList.PolynomialDesign+1,aElementList.Dim),1,0);
            for(uint j = 0; j< (aElementList.ListOfChildren).length(); j++)
            {
                tBasisOfElement = give_basis_of_element((aElementList.ListOfChildren)(j),aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
                //                (aElementList.IdField) = (aElementList.IdFieldOfChildren)(j);
                //                (aElementList.TMatrix) = (aElementList.TMatrixOfChildren)(j);
                tTmatrix = (aElementList.tTMatrixOfChildren)(j);
                //                                                (aElementList.IdField).print("(ElementList.IdField)");
                //                                                (aElementList.TMatrix).print(" (ElementList.TMatrix) ");
                tOldSol.set_size(tBasisOfElement.length(),1,0);
                for(uint k = 0; k < tBasisOfElement.length(); k++)
                {
                    //                  tOldSol(k) = (aBasisList.NodalPDV)(tBasisOfElement(k));
                    tOldSol(k) = (aBasisList.NodalADVNewField)(tBasisOfElement(k));
                }
                //                tOldSol.print("tOldSol");
                //                tNewSol = trans((ElementList.TMatrix))*tOldSol;
                tNewSol = tOldSol;

                //                                tNewSol.print("tNewSol");
                //                tTmatrix.print("tTmatrix");
                tLevel = give_element_level((aElementList.ListOfChildren)(j),aElementList.Dim,aElementList.NumElements);
                ElementMass = MassMatrix_for_L2_projection(aElementList.Dim,aElementList.PolynomialDesign);
                ElementRHSDummy = 1/pow(pow(2,tLevel),aElementList.Dim) * tTmatrix * ElementMass * tNewSol;
                ElementRHS = ElementRHS + ElementRHSDummy;
            }
        }
        else
        {
            if( aElementList.TruncatedBsplines == false)
            {
                this->give_Tmatrix_and_id_field_for_designvariables_projection((aElementList.ElementListActiveDesign)(i),aElementList,aBasisList);
            }
            else
            {
                this->give_Truncated_Tmatrix_and_id_field_for_designvariables_projection((aElementList.ElementListActiveDesign)(i),aElementList,aBasisList);
            }
            tBasis = give_basis_of_element((aElementList.ElementListActiveDesign)(i),aElementList.Dim,aElementList.PolynomialDesign,aElementList.NumElements);
            real tWeight;
            for(uint j = 0; j < (aElementList.TMatrix).n_cols(); j++)
            {
                tNewSol(j) = 0.0;
                tNewSolPDV(j) = 0.0;
                tWeight = 0.0;
                for(uint k = 0; k < (aElementList.TMatrix).n_rows(); k++)
                {
                    tWeight += (aElementList.TMatrix)(k,j);
                    tNewSolPDV(j) += (aElementList.TMatrix)(k,j) * (aBasisList.NodalPDV)((aElementList.IdField)(k));
                    tNewSol(j) += (aElementList.TMatrix)(k,j) * (aBasisList.NodalADVNewField)((aElementList.IdField)(k));
                }
                if ( aElementList.PerformNormalization )
                {
                    tNewSol(j) /= tWeight;
                }
            }
            //                        tNewSol = trans((ElementList.TMatrix))*tOldSol;
            //                                                            tNewSol.print("tNewSol");
            //                                                            tNewSolPDV.print("tNewSolPDV");
            tLevel = give_element_level((aElementList.ElementListActiveDesign)(i),aElementList.Dim,aElementList.NumElements);
            ElementRHS = RHS_for_L2_projection(aElementList.Dim,aElementList.PolynomialDesign,tNewSol);
            ElementRHS = 1/pow(pow(2,tLevel),aElementList.Dim)  * ElementRHS;
        }

        //                    ElementRHS.print("ElementRHS");
        if( aElementList.TruncatedBsplines == false)
        {
            this->give_Tmatrix_and_id_field_design((aElementList.ElementListActiveDesign)(i),aElementList,aBasisList);
        }
        else
        {
            this->give_Truncated_Tmatrix_and_id_field_design((aElementList.ElementListActiveDesign)(i),aElementList,aBasisList);
        }
        //        (aElementList.IdField).print("(ElementList.IdField) new");
        //        (aElementList.TMatrix).print("(ElementList.TMatrix) new");
        //        ElementRHS.print("ElementRHS");
        ElementRHS = (aElementList.TMatrix) * ElementRHS;

        //        ElementRHS.print("ElementRHS");

        ElementMass = MassMatrix_for_L2_projection(aElementList.Dim,aElementList.PolynomialDesign);

        tLevel = give_element_level((aElementList.ElementListActiveDesign)(i),aElementList.Dim,aElementList.NumElements);
        ElementMass = 1/pow(pow(2,tLevel),aElementList.Dim)  * (aElementList.TMatrix) * ElementMass * trans((aElementList.TMatrix));
        //        ElementMass.print("ElementMass");
        //                                            ElementRHS.print("ElementRHS");
        for(uint j = 0; j< (aElementList.IdField).length(); j++)
        {
            //                        std::cout << " (ElementList.IdField)(j) " << (aElementList.IdField)(j) << " (BasisList.DesignBSplineListMap)((ElementList.IdField)(j)) " << (aBasisList.DesignBSplineListMap)((aElementList.IdField)(j)) << std::endl;
            GLBRHS((aBasisList.DesignBSplineListMap)((aElementList.IdField)(j))) += ElementRHS(j);
            for(uint k = 0; k < (aElementList.IdField).length(); k++)
            {
                //                std::cout << " (ElementList.IdField)(j) " << (aElementList.IdField)(j) << " (BasisList.DesignBSplineListMap)((ElementList.IdField)(j)) " << (aBasisList.DesignBSplineListMap)((aElementList.IdField)(j)) << " (BasisList.DesignBSplineListMap)((ElementList.IdField)(k)) " << (aBasisList.DesignBSplineListMap)((aElementList.IdField)(k)) << std::endl;
                RowInd(tVar) = (aBasisList.DesignBSplineListMap)((aElementList.IdField)(j));
                ColInd(tVar) = (aBasisList.DesignBSplineListMap)((aElementList.IdField)(k));
                Values(tVar) = ElementMass(j,k);
                tVar++;
            }
        }
    }
    RowInd.resize(tVar,1);
    ColInd.resize(tVar,1);
    Values.resize(tVar,1);
    Sp_Mat<real> GLBMass(RowInd,ColInd,Values,GLBRHS.length(),GLBRHS.length()); // Generate global matrix from vectors

    real tElapsedTimeBuildSystem = tBuildSystemTiming.toc<moris::chronos::seconds>().wall;
    std::fprintf(stdout,"Time for building projection system : %f [sec]\n",tElapsedTimeBuildSystem);

    moris::tic tSolveSystemTiming;
    Mat<real> GLBBsplineCoeff = solve(GLBMass,GLBRHS , "superlu");

    real tElapsedTimeSolveSystem = tSolveSystemTiming.toc<moris::chronos::seconds>().wall;
    std::fprintf(stdout,"Time for solving projection system : %f [sec]\n",tElapsedTimeSolveSystem);

    (aElementList.Coeff) = GLBBsplineCoeff;

    FILE* fid=std::fopen("AbsDesVariables_new.dat","wb");
    double * buffer = (double*) alloca(GLBBsplineCoeff.length()*sizeof(double));
    for(uint i = 0; i < GLBBsplineCoeff.length(); i++)
        buffer[i] = GLBBsplineCoeff(i);
    std::fwrite(buffer,sizeof(double),GLBBsplineCoeff.length(),fid);
    std::fclose(fid);

}

void moris::Hierarchical_Mesh::initial_refinement(
        ElementList_struc & aElementList)
{
    uint tMaxNumElemLevel = 0;
    for(uint level = 0; level<aElementList.InitialRefinement; level ++)
    {
        tMaxNumElemLevel += pow(pow(2,aElementList.Dim),level);
    }
    Mat<uint> tDeactiveElements((aElementList.ElementListOnProcInit).length()*tMaxNumElemLevel,1,UINT_MAX);
    tMaxNumElemLevel -= pow(pow(2,aElementList.Dim),aElementList.InitialRefinement-1);
    tDeactiveElements.rows(0,(aElementList.ElementListOnProcInit).length()-1) = (aElementList.ElementListOnProcInit).rows(0,(aElementList.ElementListOnProcInit).length()-1);
    uint tVar = (aElementList.ElementListOnProcInit).length();
    Mat<uint> tChildren(pow(2,aElementList.Dim),1);
    Mat<uint> tChildrenDummy(pow(2,aElementList.Dim),1);
    for(uint i = 0; i<(aElementList.ElementListOnProcInit).length()*tMaxNumElemLevel; i++)
    {
        tChildren = give_children_of_element(tDeactiveElements(i),aElementList.Dim,aElementList.NumElements);
        tDeactiveElements.rows(tVar+i*pow(2,aElementList.Dim),tVar+(i+1)*pow(2,aElementList.Dim)-1) = tChildren.rows(0,tChildren.length()-1);
    }
    (aElementList.DeactivateElement) = unique(tDeactiveElements);
}

void
moris::Hierarchical_Mesh::floodfill_for_STL(
        Mat<real> aNodalField,
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    this->give_active_elements(aElementList); // One layer of elements from the neighbor procs are taken into account
    uint tVarb = 0;
    uint tVarc = 1;
    real tVal = 1.0;
    //    uint tSwitch = 0;
    Mat<uint> tElementNeighbor;
    Mat<uint> tBasisOfElement;
    Mat<uint> tPossibleElements;
    BoostBitset tActiveElements = aElementList.ElementActive;
    //    BoostBitset tActiveBasis = mBasisData.BasisActive;
    Mat<uint> tCheckNeighbors;
    real tSum;
    Mat<real> tElementValues(pow(aElementList.Polynomial+1,aElementList.Dim),1,0);
    if( aElementList.Dim == 2)
    {
        Mat<uint> tCheckdirections2D = {{0, 3},{0, 1},{3, 2},{1, 2}};
        tCheckNeighbors = tCheckdirections2D;
    }
    else if( aElementList.Dim == 3)
    {
        Mat<uint> tCheckdirections3D = {{0, 3, 4},{0, 1, 4},{3, 2, 4},{1, 2, 4}, {0, 3, 5},{0, 1, 5},{3, 2, 5},{1, 2, 5}};
        tCheckNeighbors = tCheckdirections3D;
    }
    while( tActiveElements.count() > 0)
    {
        tPossibleElements.set_size(tActiveElements.count()+1,1,0);
        tPossibleElements(0) = tActiveElements.find_first(); // Start with the first active element;
        tActiveElements.reset(tPossibleElements(0));
        tVarb = 0;
        tVarc = 1;
        while( tPossibleElements(tVarb) > 0)
        {
            //            aNodalField.print("aNodalField");
            tBasisOfElement = give_basis_of_element(tPossibleElements(tVarb),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
            tElementNeighbor = give_active_neighbor_of_element(tPossibleElements(tVarb),aElementList.Dim,aElementList.NumElements,aElementList.Level,aElementList.ElementActive);
            //            tElementNeighbor.print("tElementNeighbor");
            //            tCheckNeighbors.print("tCheckNeighbors");
            for(uint i = 0; i < tBasisOfElement.length(); i++)
            {
                tElementValues(i) = aNodalField(tBasisOfElement(i)) ;
            }
            tSum = sum(tElementValues);
            if( tSum == 0.0 )
            {
                tVal = 1;
            }
            else if( tElementValues.max() > 0.0 )
            {
                tVal = 1;
            }
            else if( tElementValues.min() < 0.0 )
            {
                tVal = -1;
            }

            for(uint i = 0; i < tBasisOfElement.length(); i++)
            {
                //                tActiveBasis.reset(tBasisOfElement(i));
                if( aNodalField(tBasisOfElement(i)) == 0.0 )
                {
                    aNodalField(tBasisOfElement(i)) = tVal;
                    for(uint j = 0; j < tElementNeighbor.n_rows(); j++)
                    {
                        //                        std::cout << " tElementNeighbor(j,1) " << tElementNeighbor(j,1) << " tCheckNeighbors(i,0) " << tCheckNeighbors(i,0) << " tElementNeighbor(j,0) " << tElementNeighbor(j,0) << " tActiveElements.test(tElementNeighbor(j,0))  " << tActiveElements.test(tElementNeighbor(j,0))  << std::endl;
                        if( tElementNeighbor(j,1) == tCheckNeighbors(i,0) && tActiveElements.test(tElementNeighbor(j,0)) == 1)
                        {
                            tPossibleElements(tVarc) = tElementNeighbor(j,0);
                            tVarc++;
                            tActiveElements.reset(tElementNeighbor(j,0));
                        }
                        //                        std::cout << " tElementNeighbor(j,1) " << tElementNeighbor(j,1) << " tCheckNeighbors(i,1) " << tCheckNeighbors(i,1) << " tElementNeighbor(j,0) " << tElementNeighbor(j,0) << " tActiveElements.test(tElementNeighbor(j,0))  " << tActiveElements.test(tElementNeighbor(j,0))  << std::endl;
                        if( tElementNeighbor(j,1) == tCheckNeighbors(i,1) && tActiveElements.test(tElementNeighbor(j,0)) == 1)
                        {
                            tPossibleElements(tVarc) = tElementNeighbor(j,0);
                            tVarc++;
                            tActiveElements.reset(tElementNeighbor(j,0));
                        }
                        if( aElementList.Dim == 3 && tElementNeighbor(j,1) == tCheckNeighbors(i,2) && tActiveElements.test(tElementNeighbor(j,0)) == 1)
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

Mat<uint>
moris::Hierarchical_Mesh::give_active_neighbor_of_element(
        uint & aElement,
        uint & aDim,
        Mat<uint> & aNumElements,
        uint & aLevel,
        BoostBitset & aElementActive)
{
    uint tBuffer = 1;
    Mat<uint> tNeighbor = give_neighbour_of_element(aElement,aDim,tBuffer,aNumElements);
    Mat<uint> tPossibleNeighbors;
    Mat<uint> tActiveNeighbors(2*aDim*pow(pow(2,aLevel),aDim),2,UINT_MAX);
    uint tVar = 0, tVara = 0, tVarb = 1;
    uint tParent = 0;
    uint tLevel = 0;
    Mat<uint> tChildren;
    Mat<uint> tNeighbors;
    Mat<uint> tNeighborsOrdinal;
    Mat<uint> tChildNeighbors;
    Mat<uint> tChildNeighborsOrdinal;
    uint tSwitch = 0;
    if( aDim == 2)
    {
        Mat<uint> tNeighbors2D = {{1, 3, 5, 7 }}; // Neighbors in 2D (positions extracted from give_neighbor_of_element)
        tNeighbors = tNeighbors2D;
        Mat<uint> tNeighborsOrdinal2D = {{0, 3, 1, 2 }}; // Neighbors in 2D (positions extracted from give_neighbor_of_element)
        tNeighborsOrdinal = tNeighborsOrdinal2D;
        Mat<uint> tChildNeighbors2D = {{2, 3}, {1, 3}, {0, 2}, {0, 1}}; // Possible childrens for each neighbor
        tChildNeighbors = tChildNeighbors2D;
        Mat<uint> tChildNeighborsOrdinal2D = {{0, 0}, {3, 3}, {1, 1}, {2, 2}}; // Possible childrens for each neighbor
        tChildNeighborsOrdinal = tChildNeighborsOrdinal2D;
    }
    else if( aDim == 3)
    {
        Mat<uint> tNeighbors3D = {{4, 10, 12, 14, 16, 22}};  // Neighbors in 2D (positions extracted from give_neighbor_of_element)
        tNeighbors = tNeighbors3D;
        Mat<uint> tNeighborsOrdinal3D = {{4, 0, 3, 1, 2, 5}};  // Neighbors in 2D (positions extracted from give_neighbor_of_element)
        tNeighborsOrdinal = tNeighborsOrdinal3D;
        Mat<uint> tChildNeighbors3D = {{4, 5, 6, 7}, {2, 3, 6, 7}, {1, 3, 5, 7}, {0, 2, 4, 6}, {0, 1, 4, 5}, {0, 1, 2, 3}}; // Possible childrens for each neighbor
        tChildNeighbors = tChildNeighbors3D;
        Mat<uint> tChildNeighborsOrdinal3D = {{4, 4, 4, 4}, {0, 0, 0, 0}, {3, 3, 3, 3}, {1, 1, 1, 1}, {2, 2, 2, 2}, {5, 5, 5, 5}}; // Possible childrens for each neighbor
        tChildNeighborsOrdinal = tChildNeighborsOrdinal3D;
    }
    for(uint i = 0; i < tNeighbors.length(); i++)
    {
        if( aElementActive.test(tNeighbor(tNeighbors(i))) == 1 )
        {
            tActiveNeighbors(tVar,0) = tNeighbor(tNeighbors(i));
            tActiveNeighbors(tVar,1) = tNeighborsOrdinal(i);
            tVar++;
        }
        else
        {
            tLevel = give_element_level(tNeighbor(tNeighbors(i)),aDim,aNumElements);
            if( tLevel+1 <= aLevel)
            {
                tPossibleNeighbors.set_size(2*aDim*pow(2,aLevel),1,0);
                tPossibleNeighbors(0) =  tNeighbor(tNeighbors(i));
                tVara = 0;
                tVarb = 1;
                tSwitch = 0;
                while( tPossibleNeighbors(tVara) > 0 )
                {
                    tParent = give_parent_of_element(tPossibleNeighbors(tVara),aDim,aNumElements);
                    if( tParent < UINT_MAX && aElementActive.test(tParent) == 1 )
                    {
                        tActiveNeighbors(tVar,0) = tParent;
                        tActiveNeighbors(tVar,1) = tNeighborsOrdinal(i);
                        tVar++;
                        tSwitch = 1; // If parent is found, no need to check childrens
                    }
                    else if(  tParent < UINT_MAX )
                    {
                        tPossibleNeighbors(tVarb) =  tParent;
                        tVarb++;
                    }
                    tVara++;
                }
                tPossibleNeighbors.set_size(2*aDim*pow(2,aLevel),1,0);
                tPossibleNeighbors(0) =  tNeighbor(tNeighbors(i));
                tVara = 0;
                tVarb = 1;
                while( tPossibleNeighbors(tVara) > 0 && tSwitch == 0 )
                {
                    tChildren = give_children_of_element(tPossibleNeighbors(tVara),aDim,aNumElements); // Give the children of the parent element
                    tLevel = give_element_level(tChildren(0),aDim,aNumElements);
                    if( aElementActive.test(tChildren(tChildNeighbors(i,0))) == 1 )
                    {
                        tActiveNeighbors(tVar,0) = tChildren(tChildNeighbors(i,0));
                        tActiveNeighbors(tVar,1) = tChildNeighborsOrdinal(i,0);
                        tVar++;
                    }
                    else if( tLevel < aLevel )
                    {
                        tPossibleNeighbors(tVarb) =  tChildren(tChildNeighbors(i,0));
                        tVarb++;
                    }
                    if( aElementActive.test(tChildren(tChildNeighbors(i,1))) == 1 )
                    {
                        tActiveNeighbors(tVar,0) = tChildren(tChildNeighbors(i,1));
                        tActiveNeighbors(tVar,1) = tChildNeighborsOrdinal(i,1);
                        tVar++;
                    }
                    else if( tLevel < aLevel )
                    {
                        tPossibleNeighbors(tVarb) =  tChildren(tChildNeighbors(i,1));
                        tVarb++;
                    }
                    if( aDim == 3)
                    {
                        if( aElementActive.test(tChildren(tChildNeighbors(i,2))) == 1 )
                        {
                            tActiveNeighbors(tVar,0) = tChildren(tChildNeighbors(i,2));
                            tActiveNeighbors(tVar,1) = tChildNeighborsOrdinal(i,2);
                            tVar++;
                        }
                        else if( tLevel < aLevel )
                        {
                            tPossibleNeighbors(tVarb) =  tChildren(tChildNeighbors(i,2));
                            tVarb++;
                        }
                        if( aElementActive.test(tChildren(tChildNeighbors(i,3))) == 1 )
                        {
                            tActiveNeighbors(tVar,0) = tChildren(tChildNeighbors(i,3));
                            tActiveNeighbors(tVar,1) = tChildNeighborsOrdinal(i,3);
                            tVar++;
                        }
                        else if( tLevel < aLevel )
                        {
                            tPossibleNeighbors(tVarb) =  tChildren(tChildNeighbors(i,3));
                            tVarb++;
                        }
                    }
                    tVara++;
                }
            }
        }
    }
    tActiveNeighbors.resize(tVar,2);
    return tActiveNeighbors;
}

void moris::Hierarchical_Mesh::find_elements_of_object(
        Mat<real> aNodalField,
        Mat<uint> & tRefineElement,
        uint & MaxLevelOfRefinement,
        ElementList_struc & aElementList,
        BasisList_struc & aBasisList)
{
    Mat<uint> tBasis;
    uint tLevel;
    uint tVar = 0;
    uint tMaxRefinementElement;

    Mat<real> tElementSolution(pow(aElementList.Polynomial+1,aElementList.Dim),1,0);
    Mat<uint> tChildren;
    tRefineElement.set_size((aElementList.ElementListOnProc).length(),1,0);
    uint tNumberElements = give_number_of_elements(MaxLevelOfRefinement,aElementList.Dim,aElementList.NumElements);
    if( aElementList.ElementActive.size() < tNumberElements)
        aElementList.ElementActive.resize(tNumberElements);

    BoostBitset tElementActiveDummy = aElementList.ElementActive;
    for(uint i = 0; i < (aElementList.ElementListOnProc).length(); i++)
    {
//                tBasis = give_basis_of_element((aElementList.ElementListOnProc)(i),aElementList.Dim,aElementList.Polynomial,aElementList.NumElements);
        this->give_Tmatrix_and_id_field((aElementList.ElementListOnProc)(i),aElementList,aBasisList);
        tBasis = aElementList.IdField;
        tElementSolution.set_size(tBasis.length(),1,0);
        for(uint j = 0; j < tBasis.length(); j++)
        {
            tElementSolution(j) = aNodalField(tBasis(j));
        }
        if( tElementSolution.max() > 0.0 && tElementSolution.min() < 0.0 )
        {
            tLevel = give_element_level((aElementList.ElementListOnProc)(i),aElementList.Dim,aElementList.NumElements);
            if( tLevel < MaxLevelOfRefinement)
            {
                tRefineElement(tVar) = (aElementList.ElementListOnProc)(i);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = give_children_of_element(tRefineElement(tVar),aElementList.Dim,aElementList.NumElements);
                for(uint l = 0; l < tChildren.length(); l++)
                    tElementActiveDummy.set(tChildren(l));
                tVar++;
            }
            else
            {
                tLevel = MaxLevelOfRefinement-1;
                tRefineElement(tVar) = give_parent_of_level_x((aElementList.ElementListOnProc)(i),aElementList.Dim,aElementList.NumElements,tLevel);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = give_children_of_element(tRefineElement(tVar),aElementList.Dim,aElementList.NumElements);
                for(uint l = 0; l < tChildren.length(); l++)
                    tElementActiveDummy.set(tChildren(l));
                tVar++;
            }
        }
    }
    tRefineElement.resize(tVar,1);

    // Parent elements need also be on the refinement list
    if( tVar > 0)
    {
        tVar = tRefineElement.length();
        uint tLengthRefinement = tRefineElement.length();
        tRefineElement.resize(tRefineElement.length()*(tLevel+1),1);
        if( tLevel > 0)
        {
            for(uint i = 0; i < tLengthRefinement; i++ )
            {
                tLevel = give_element_level(tRefineElement(i),aElementList.Dim,aElementList.NumElements);
                for(uint j = 0; j < tLevel; j++)
                {
                    tRefineElement(tVar) = give_parent_of_level_x(tRefineElement(i),aElementList.Dim,aElementList.NumElements,j);
                    tElementActiveDummy.reset(tRefineElement(tVar));
                    tChildren = give_children_of_element(tRefineElement(tVar),aElementList.Dim,aElementList.NumElements);
                    for(uint l = 0; l < tChildren.length(); l++)
                        tElementActiveDummy.set(tChildren(l));
                    tVar++;
                }
            }
        }
        tRefineElement.resize(tVar,1);
        tRefineElement = unique(tRefineElement);
    }

    if( tVar > 0 )
    {
        tMaxRefinementElement = tRefineElement.max();
        tLevel = give_element_level(tMaxRefinementElement,aElementList.Dim,aElementList.NumElements) + 1;
        tNumberElements = give_number_of_elements(tLevel,aElementList.Dim,aElementList.NumElements);
        if( (aElementList.ElementActive).size() < tNumberElements )
        {
            (aElementList.ElementActive).resize(tNumberElements);
            tElementActiveDummy.resize(tNumberElements);
        }
    }

    // Add a buffer layer of elements
    if( tVar > 0 && aElementList.BufferElements > 0)
    {
        Mat<uint> tChildren;
        Mat<uint> tNeighbour(pow(1+2*aElementList.BufferElements*pow(2,aElementList.Level),aElementList.Dim),1,0); // Highest possible number of neigbours
        tVar = tRefineElement.length(); // Starting point to save new elements, which need to be deactivated
        uint tLengthRefinement = tVar; // Temporary variable for length of loop
        uint tBuffer = 0;
        tRefineElement.resize(tRefineElement.length()*(1+tNeighbour.length()),1); // Resize list of elements to the maximum possible list
        for(uint j = 0; j < tLengthRefinement; j++)
        {
            tLevel = give_element_level(tRefineElement(j),aElementList.Dim,aElementList.NumElements);
            if( aElementList.AdaptiveBufferLayer == true)
            {
                if( aElementList.Staircasebuffer == true)
                {
                    tBuffer = aElementList.BufferElements*(1+tLevel);
                }
                else
                {
                    tBuffer = aElementList.BufferElements*pow(2,tLevel);
                }
            }
            else
            {
                tBuffer = aElementList.BufferElements;
            }
            tNeighbour = give_neighbour_of_element(tRefineElement(j),aElementList.Dim,tBuffer,aElementList.NumElements);
            for(uint k = 0; k<tNeighbour.length(); k++)
            {
                if( tNeighbour(k) < tElementActiveDummy.size() &&  tElementActiveDummy.test(tNeighbour(k)) == 1 )
                {
                    tRefineElement(tVar) = tNeighbour(k);
                    tElementActiveDummy.reset(tNeighbour(k));
                    tChildren = give_children_of_element(tNeighbour(k),aElementList.Dim,aElementList.NumElements);
                    for(uint l = 0; l < tChildren.length(); l++)
                        tElementActiveDummy.set(tChildren(l));
                    tVar++;
                }
            }
        }
        tRefineElement.resize(tVar,1);
        tRefineElement = unique(tRefineElement);
    }

    // Parent elements need also be on the refinement list
    if( tVar > 0)
    {
        tMaxRefinementElement = tRefineElement.max();
        tLevel = give_element_level(tMaxRefinementElement,aElementList.Dim,aElementList.NumElements);
        tVar = tRefineElement.length();
        uint tLengthRefinement = tRefineElement.length();
        tRefineElement.resize(tRefineElement.length()*(tLevel+1),1);
        if( tLevel > 0)
        {
            for(uint i = 0; i < tLengthRefinement; i++ )
            {
                tLevel = give_element_level(tRefineElement(i),aElementList.Dim,aElementList.NumElements);
                for(uint j = 0; j < tLevel; j++)
                {
                    tRefineElement(tVar) = give_parent_of_level_x(tRefineElement(i),aElementList.Dim,aElementList.NumElements,j);
                    tElementActiveDummy.reset(tRefineElement(tVar));
                    tChildren = give_children_of_element(tRefineElement(tVar),aElementList.Dim,aElementList.NumElements);
                    for(uint l = 0; l < tChildren.length(); l++)
                        tElementActiveDummy.set(tChildren(l));
                    tVar++;
                }
            }
        }
        tRefineElement.resize(tVar,1);
        tRefineElement = unique(tRefineElement);
    }
}


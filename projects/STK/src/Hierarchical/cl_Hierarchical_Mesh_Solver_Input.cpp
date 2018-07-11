/*
 * cl_Hierarchical_Mesh_Solver_Input.cpp
 *
 *  Created on: Apr 9, 2018
 *      Author: schmidt
 */
#include "cl_Hierarchical_Mesh_Solver_Input.hpp"
//#include "cl_Hierarchical_Mesh_Main.hpp"

using namespace moris;

Hierarchical_Mesh_Solver_Input::Hierarchical_Mesh_Solver_Input( moris::real aRhsScaling,
                                                                moris::real  aRhsOffset,
                                                                Hierarchical_Mesh_Main*   aHMR ) : mHMR( aHMR )
{
    mRhsScaling = aRhsScaling;
    mRhsOffset  = aRhsOffset;
}

Hierarchical_Mesh_Solver_Input::~Hierarchical_Mesh_Solver_Input()
{
}

uint Hierarchical_Mesh_Solver_Input::get_num_my_dofs()
{
    // FIXME : this only works in serial at the moment
    MORIS_ASSERT(  par_size() == 1 ,
            " get_num_my_dofs does not work in parallel at the moment.");

    return mHMR->mBasisData.DesignBSplineActiveList.length();
}

uint Hierarchical_Mesh_Solver_Input::get_num_element_dof()
{
    return 0;
}

//--------------------------------------------------------------------------------

Mat <int>
Hierarchical_Mesh_Solver_Input::get_my_local_global_map()
{

    // FIXME : this only works in serial at the moment
    MORIS_ASSERT(  par_size() == 1 ,
            " get_my_local_global_map does not work in parallel at the moment.");

    // create a mat of length of active bspline basis (= number of non-hanging nodes )
    Mat <int> aMat(mHMR->mBasisData.NumberOfNonHangingNodes, 1 );
    for (uint k=0; k<mHMR->mBasisData.NumberOfNonHangingNodes; ++k)
    {
        aMat(k) = k;
    }
    return aMat;
}

//--------------------------------------------------------------------------------

uint
Hierarchical_Mesh_Solver_Input::get_num_my_elements()
{
    return mHMR->mElementData.ElementListActiveDesign.length();
}

//--------------------------------------------------------------------------------

void Hierarchical_Mesh_Solver_Input::get_element_matrix(
        const uint  & aMyElementInd,
        Mat< real > & aElementMatrix)
{
    // DOFs of this element (not needed in this function)
    Mat<uint> tIdField;

    // T-Matrix of this element
    Mat<real> tTMatrix;

    // create T-Matrix from current element
    if ( mHMR->mSettings.TruncatedBsplines == false)
    {
        // get normal T-Matrix
        mHMR->mTMatrix.give_Tmatrix_and_IdField(
                mHMR->mElementData.ElementListActiveDesign( aMyElementInd ),
                mHMR->mMeshData.ModelDim,
                mHMR->mMeshData.PolynomialDesign,
                mHMR->mMeshData.NumberOfElementsPerDirection,
                mHMR->mBasisData.DesignBSplineActive,
                tTMatrix,
                tIdField,
                mHMR->mElementData.TMatrixParentChildRelationDesign
                );
    }
    else
    {
        // get truncated T-Matrix
        mHMR->mTMatrix.give_Truncated_Tmatrix_and_IdField(
                mHMR->mElementData.ElementListActiveDesign( aMyElementInd ),
                mHMR->mMeshData.ModelDim,
                mHMR->mMeshData.PolynomialDesign,
                mHMR->mMeshData.NumberOfElementsPerDirection,
                mHMR->mBasisData.DesignBSplineActive,
                tTMatrix,
                tIdField,
                mHMR->mElementData.TMatrixParentChildRelationDesign);
    }

    // get mass matrix for this element
    mHMR->MassMatrix_for_L2_projection(
            mHMR->mMeshData.ModelDim,
            mHMR->mMeshData.PolynomialDesign,
            aElementMatrix );

    // get level of this element
    uint tLevel = mHMR->mBaseElement.give_element_level(
            mHMR->mElementData.ElementListActiveDesign( aMyElementInd ),
            mHMR->mMeshData.ModelDim,
            mHMR->mMeshData.NumberOfElementsPerDirection);

    // calculate element matrix for output
    aElementMatrix= 1/pow(pow(2,tLevel),mHMR->mMeshData.ModelDim)
            * tTMatrix * aElementMatrix * trans(tTMatrix);

}

//--------------------------------------------------------------------------------

void Hierarchical_Mesh_Solver_Input::get_element_topology(
        const uint & aMyElementInd,
        Mat< int > & aElementTopology)
{
    // create temporary array
    Mat< uint > tIdField;

    // call Id Field from T-Matrix module
    mHMR->mTMatrix.give_IdField(
            mHMR->mElementData.ElementListActiveDesign(aMyElementInd),
            mHMR->mMeshData.ModelDim,
            mHMR->mMeshData.PolynomialDesign,
            mHMR->mMeshData.NumberOfElementsPerDirection,
            mHMR->mBasisData.DesignBSplineActive,
            tIdField);

    // assign size for output matrix
    aElementTopology.set_size(tIdField.length(), 1);

    // convert to consecutive basis numbering
    for (uint k=0; k< tIdField.length(); ++k)
    {
        aElementTopology(k) = mHMR->mBasisData.DesignBSplineListMap.find( tIdField(k) );
    }
}

//--------------------------------------------------------------------------------

Mat< uint > Hierarchical_Mesh_Solver_Input::get_constr_dof()
{
    Mat <uint> aMat;
    return aMat;
}

//--------------------------------------------------------------------------------

void Hierarchical_Mesh_Solver_Input::get_element_rhs(
         const uint            & aMyElementInd,
         Mat< real >           & aElementRHS )
 {
//     moris::real              aRhsScaling = 1.0;
//     moris::real              aRhsOffset  = 0.0;

    Mat<real> tTMatrix;
    uint      tLevel;
    Mat<uint> tIdField;

     // number of DOFs of this element
     uint tNumberOfBasisInElement = pow(
             mHMR->mMeshData.PolynomialDesign + 1,
             mHMR->mMeshData.ModelDim );

     // list containing IDs of children of this element
     Mat< uint > tListOfChildren;

     // cell containinc T-Matrices of each child of this element
     Cell< Mat< real > > tTMatrixOfChildren;

     // check if this element needs to be coarsened or refined
     mHMR->mTMatrix.give_Tmatrix_and_IdField_Design_L2Projection_Coarsening(
             mHMR->mElementData.ElementListActiveDesign( aMyElementInd ),
             mHMR->mMeshData.ModelDim,
             mHMR->mMeshData.Polynomial,
             mHMR->mMeshData.LevelLastStep,
             mHMR->mMeshData.LevelDesignLastStep,
             mHMR->mMeshData.NumberOfElementsPerDirection,
             mHMR->mElementData.ElementActiveDesignLastStep,
             tTMatrixOfChildren,
             tListOfChildren,
             mHMR->mElementData.TMatrixParentChildRelationDesign );

     // mass matrix of this Element
     Mat<real> tElementMass;

     // call mass matrix for this element
     mHMR->MassMatrix_for_L2_projection(
              mHMR->mMeshData.ModelDim,
              mHMR->mMeshData.PolynomialDesign,
              tElementMass );

     // if this element has children, the coarsening operator needs to be used
     if ( tListOfChildren.length() >= pow( 2,  mHMR->mMeshData.ModelDim ) )
     {
         // set memory for RHS
         aElementRHS.set_size( tNumberOfBasisInElement , 1, 0 );

         // loop over all children
         for ( uint j = 0; j < tListOfChildren.length(); ++j )
         {
             // get IDs of DOFs of this child
             Mat< uint > tBSplinesOfThisChild =
                            mHMR->mHMRElement.give_basis_of_element(
                            tListOfChildren( j ),
                            mHMR->mMeshData.ModelDim,
                            mHMR->mMeshData.PolynomialDesign,
                            mHMR->mMeshData.NumberOfElementsPerDirection );

             // matrix to contain B-Spline DOFs of current child
             Mat< real > tBSplineDOFsOfThisChild ( tBSplinesOfThisChild.length(), 1, 0 );

             // loop over all DOFs of this child
             for ( uint k = 0; k < tBSplinesOfThisChild.length(); ++k )
             {
                 // copy DOFs into temporary vector
                 tBSplineDOFsOfThisChild( k ) =
                               mRhsScaling * mHMR->mFieldData.NodalADV.find( tBSplinesOfThisChild( k ) + mRhsOffset );
             }

             // get level of this child
             uint tLevel =  mHMR->mBaseElement.give_element_level(
                     tListOfChildren( j ),
                     mHMR->mMeshData.ModelDim,
                     mHMR->mMeshData.NumberOfElementsPerDirection );

             // add values to RHS
             aElementRHS = aElementRHS + 1/pow(pow(2,tLevel),mHMR->mMeshData.ModelDim)
                     * tTMatrixOfChildren( j ) * tElementMass * tBSplineDOFsOfThisChild;
         }
     }
     else
     {  // if this element has no children
         // get T-Matrix of this element
         if ( mHMR->mSettings.TruncatedBsplines == false)
         {
             // normal T-Matrix
             mHMR->mTMatrix.give_Tmatrix_and_IdField(
                     mHMR->mElementData.ElementListActiveDesign( aMyElementInd ),
                     mHMR->mMeshData.ModelDim,
                     mHMR->mMeshData.PolynomialDesign,
                     mHMR->mMeshData.NumberOfElementsPerDirection,
                     mHMR->mBasisData.DesignBSplineActiveLastStep,
                     tTMatrix,
                     tIdField,
                     mHMR->mElementData.TMatrixParentChildRelationDesign);
         }
         else
         {
             // truncated T-Matrix
             mHMR->mTMatrix.give_Truncated_Tmatrix_and_IdField(
                     mHMR->mElementData.ElementListActiveDesign( aMyElementInd ),
                     mHMR->mMeshData.ModelDim,
                     mHMR->mMeshData.PolynomialDesign,
                     mHMR->mMeshData.NumberOfElementsPerDirection,
                     mHMR->mBasisData.DesignBSplineActiveLastStep,
                     tTMatrix,
                     tIdField,
                     mHMR->mElementData.TMatrixParentChildRelationDesign );
         }
         
         // matrix to contain B-Spline DOFs of current child
         Mat< real > tBSplineDOFsOfThisElement ( tNumberOfBasisInElement, 1, 0 );

         // loop over all columns of T-Matrix
         for ( uint j = 0; j < tTMatrix.n_cols(); j++)
         {
             // weighting factor for RHS
             real tWeight=0.0;

             // loop over all columns of this T-Matrix
             for ( uint k = 0; k < tTMatrix.n_rows(); k++)
             {
                 // add weight to factor
                 tWeight    += tTMatrix( k, j );

                 // copy B-Spline DOF
                 tBSplineDOFsOfThisElement( j ) += tTMatrix( k, j )
                         * (mRhsScaling * mHMR->mFieldData.NodalADV.find( tIdField( k ) ) + mRhsOffset);
             }

             // use weight to normalize B-Spline DOF
             if ( mHMR->mSettings.PerformNormalization )
             {
                 tBSplineDOFsOfThisElement(j) /= tWeight;
             }

         }

         tLevel = mHMR->mBaseElement.give_element_level(
                 mHMR->mElementData.ElementListActiveDesign( aMyElementInd ),
                 mHMR->mMeshData.ModelDim,
                 mHMR->mMeshData.NumberOfElementsPerDirection);

         // get untransformed part of RHS
         mHMR->RHS_for_L2_projection(
                 mHMR->mMeshData.ModelDim,
                 mHMR->mMeshData.PolynomialDesign,
                 tBSplineDOFsOfThisElement,
                 aElementRHS );

         aElementRHS = 1/pow(pow(2,tLevel),mHMR->mMeshData.ModelDim) * aElementRHS;

         //Create T-Matrix from current element
         if ( mHMR->mSettings.TruncatedBsplines == false)
         {
             mHMR->mTMatrix.give_Tmatrix_and_IdField(
                     mHMR->mElementData.ElementListActiveDesign( aMyElementInd ),
                     mHMR->mMeshData.ModelDim,
                     mHMR->mMeshData.PolynomialDesign,
                     mHMR->mMeshData.NumberOfElementsPerDirection,
                     mHMR->mBasisData.DesignBSplineActive,
                     tTMatrix,
                     tIdField,
                     mHMR->mElementData.TMatrixParentChildRelationDesign);
         }
         else
         {
              mHMR->mTMatrix.give_Truncated_Tmatrix_and_IdField(
                     mHMR->mElementData.ElementListActiveDesign( aMyElementInd ),
                     mHMR->mMeshData.ModelDim,
                     mHMR->mMeshData.PolynomialDesign,
                     mHMR->mMeshData.NumberOfElementsPerDirection,
                     mHMR->mBasisData.DesignBSplineActive,
                     tTMatrix,
                     tIdField,
                     mHMR->mElementData.TMatrixParentChildRelationDesign);
         }

         // multiply RHS with T-Matrix
         aElementRHS = tTMatrix * aElementRHS ;
     }
 }


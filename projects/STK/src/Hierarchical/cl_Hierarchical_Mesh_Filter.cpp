/*
 * cl_Hierarchical_Mesh_Filter.cpp
 *
 *  Created on: Jan 25, 2018
 *      Author: gleim
 */

#include "cl_Hierarchical_Mesh_Filter.hpp" // STK/src/Hierarchical
#include "fn_norm.hpp"
using namespace moris;

void Hierarchical_Mesh_Filter::Update_IDandTMatrix_design(
        uint const & aDim,
        uint const & aPolynomialDesign,
        Mat<uint> const & aNumElements,
        uint const & aLevel,
        real const & aFilterRadius,
        BoostBitset & aElementActiveDesign,
        BoostBitset & aDesignBSplineActive,
        Mat<real> const & aDimensions,
        Mat<real> const & aDimensions_Offset,
        Mat<real> const & aDimensionsOriginal,
        Mat<uint> & aIdFieldDesignNodalField,
        Mat<real> & aTMatrixDesignNodalField,
        bool const & aFilterLevelWeight,
        bool const & aPerformNormalization,
        Mat<real> const & aPointOfOrigin,
        Mat<uint> & IdFieldDesignNodalList)
{
    //    Mat<uint> tOldIdList;
    Mat<real> tBasisFiler;
    Mat<real> tNewTMatrix;
    Mat<uint> tNewIdField;
    Cell<Mat<real>> tNewTMatrixCell(aIdFieldDesignNodalField.n_rows());
    Cell<Mat<uint>> tNewIdFieldCell(aIdFieldDesignNodalField.n_rows());
    real tWeight;
    Mat<uint> tFindUniqueList;
    Mat<real> tHelp;
    Mat<real> tHelpa, tHelpb;
    Mat<uint> tHelpc;
    Mat<uint> tUniqueMap;
    uint tCopy1 = 0, tCopy2 = 0; // Length for beginning and end is needed
    uint tMaxIDField = 0;
    Mat<real> tBasisFilterList;
    // Serial uses a global numbering and saving for the BasisFilterCell (faster in serial), MPI uses a local numbering
    if( par_size() == 1)
    {
        Cell<Mat<real>> tBasisFilterCell(IdFieldDesignNodalList.max()+1);
        for(uint j = 0; j< IdFieldDesignNodalList.length(); j++)
        {
            // tBasisFilterList = Filter_for_smoothing(aIdFieldDesignNodalField(i,j+1),aDim,aPolynomialDesign,aNumElements,aLevel,aFilterRadius,aElementActiveDesign,aDesignBSplineActive,aDimensions,aDimensions_Offset); // Slower then the new filter
            tBasisFilterList = Filter_for_smoothing_new(IdFieldDesignNodalList(j),aDim,aPolynomialDesign,aNumElements,aLevel,aFilterRadius,aDesignBSplineActive,aDimensions,aDimensions_Offset,aDimensionsOriginal,aPointOfOrigin);
            tBasisFilterCell(IdFieldDesignNodalList(j)) = tBasisFilterList;
        }
        for(uint i=0; i < aIdFieldDesignNodalField.n_rows(); i++)
        {
            tCopy1 = 0;
            tCopy2 = 0;
            for(uint j = 0; j< aIdFieldDesignNodalField(i,0); j++)
            {
                tCopy2 += (tBasisFilterCell(aIdFieldDesignNodalField(i,j+1))).n_rows();
            }
            tBasisFiler.set_size(tCopy2,3);
            tCopy1 = 0;
            tCopy2 = 0;
            for(uint j = 0; j< aIdFieldDesignNodalField(i,0); j++)
            {
                tBasisFilterList = tBasisFilterCell(aIdFieldDesignNodalField(i,j+1));
                tWeight = 0.0;
                tHelpa = tBasisFilterList.col(1);
                if( aFilterLevelWeight == true)
                {
                    tHelpb = tBasisFilterList.col(2);
                }
                else
                {
                    tHelpb.set_size(tHelpa.length(),1,1);
                }
                tWeight = dot(tHelpa,tHelpb);
                tBasisFilterList.col(1) /= tWeight;
                tBasisFilterList.col(1) *= aTMatrixDesignNodalField(i,j+1);
                tCopy2 += tBasisFilterList.n_rows();
                tBasisFiler.rows(tCopy1,tCopy2-1) = tBasisFilterList.rows(0,tBasisFilterList.n_rows()-1);
                tCopy1 += tBasisFilterList.n_rows();
            }
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
                if( aFilterLevelWeight == true)
                {
                    tNewTMatrix(tUniqueMap(tBasisFiler(j,0))) += tBasisFiler(j,1)*tBasisFiler(j,2);
                }
                else
                {
                    tNewTMatrix(tUniqueMap(tBasisFiler(j,0))) += tBasisFiler(j,1);
                }
            }
            if( aPerformNormalization == true)
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
    }
    else
    {
        for(uint i=0; i < aIdFieldDesignNodalField.n_rows(); i++)
        {
//            std::cout << " i " << i << " aIdFieldDesignNodalField(i,0) " << aIdFieldDesignNodalField(i,0) << std::endl;
            Cell<Mat<real>> tBasisFilterCell(aIdFieldDesignNodalField(i,0));
            tCopy1 = 0;
            tCopy2 = 0;
            for(uint j = 0; j< aIdFieldDesignNodalField(i,0); j++)
            {
//                std::cout << " j " << j << " aIdFieldDesignNodalField(i,j+1) " << aIdFieldDesignNodalField(i,j+1) << std::endl;
                //tBasisFilterList = Filter_for_smoothing(aIdFieldDesignNodalField(i,j+1),aDim,aPolynomialDesign,aNumElements,aLevel,aFilterRadius,aElementActiveDesign,aDesignBSplineActive,aDimensions,aDimensions_Offset);// Slower then the new filter
                tBasisFilterList = Filter_for_smoothing_new(aIdFieldDesignNodalField(i,j+1),aDim,aPolynomialDesign,aNumElements,aLevel,aFilterRadius,aDesignBSplineActive,aDimensions,aDimensions_Offset,aDimensionsOriginal,aPointOfOrigin);
                tWeight = 0.0;
                tHelpa = tBasisFilterList.col(1);
                if( aFilterLevelWeight == true)
                {
                    tHelpb = tBasisFilterList.col(2);
                }
                else
                {
                    tHelpb.set_size(tHelpa.length(),1,1);
                }
                tWeight = dot(tHelpa,tHelpb);
                tBasisFilterList.col(1) /= tWeight;
                tBasisFilterList.col(1) *= aTMatrixDesignNodalField(i,j+1);
                tCopy2 += tBasisFilterList.n_rows();
                tBasisFilterCell(j) = tBasisFilterList;
            }
            tBasisFiler.set_size(tCopy2,3);
            tCopy1 = 0;
            tCopy2 = 0;
            for(uint j = 0; j< aIdFieldDesignNodalField(i,0); j++)
            {
                tCopy2 += (tBasisFilterCell(j)).n_rows();
                tBasisFiler.rows(tCopy1,tCopy2-1) = (tBasisFilterCell(j)).rows(0,(tBasisFilterCell(j)).n_rows()-1);
                tCopy1 += (tBasisFilterCell(j)).n_rows();
            }
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
                if( aFilterLevelWeight == true)
                {
                    tNewTMatrix(tUniqueMap(tBasisFiler(j,0))) += tBasisFiler(j,1)*tBasisFiler(j,2);
                }
                else
                {
                    tNewTMatrix(tUniqueMap(tBasisFiler(j,0))) += tBasisFiler(j,1);
                }
            }
            if( aPerformNormalization == true)
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
    }
    aIdFieldDesignNodalField.set_size(aIdFieldDesignNodalField.n_rows(),tMaxIDField,0);
    aTMatrixDesignNodalField.set_size(aIdFieldDesignNodalField.n_rows(),tMaxIDField,0);
    for(uint i=0; i < aIdFieldDesignNodalField.n_rows(); i++)
    {
        tHelpc =  tNewIdFieldCell(i).row(0);
        tHelpc.resize(1,tMaxIDField);
        tHelpa =  tNewTMatrixCell(i).row(0);
        tHelpa.resize(1,tMaxIDField);
        aIdFieldDesignNodalField.row(i) = tHelpc.row(0);
        aTMatrixDesignNodalField.row(i) = tHelpa.row(0);
    }
}

Mat<real>
Hierarchical_Mesh_Filter::Filter_for_smoothing(
        uint const & aBasis,
        uint const & aDim,
        uint const & aPolynomialDesign,
        Mat<uint> const & aNumElements,
        uint const & aLevel,
        real const & aFilterRadius,
        BoostBitset & aElementActiveDesign,
        BoostBitset & aDesignBSplineActive,
        Mat<real> const & aDimensions,
        Mat<real> const & aDimensions_Offset)
{
    Mat<uint> tBasisPosition = mBasis.give_position_of_basis(aBasis,aDim,aPolynomialDesign,aNumElements);
    Mat<real> tBasisCoordinate = mBasis.give_coordinate_from_basis(aBasis,aDim,aPolynomialDesign,aNumElements,aDimensions,aDimensions_Offset);
    uint taBasisLevel = mBasis.give_basis_level(aBasis,aDim,aPolynomialDesign,aNumElements);
    Mat<uint> tElementOfBasis = mBasis.give_element_of_basis(aBasis,aDim,aPolynomialDesign,aNumElements);
    uint tWhichLevel = 0;
    uint tParent = UINT_MAX;
    uint tElementOfaBasis = UINT_MAX; // Save the Element in which the Basis function live.
    for(uint i = 0; i < tElementOfBasis.length(); i++)
    {
        if( aElementActiveDesign.test(tElementOfBasis(i)) == 1)
        {
            tParent = mBaseElement.give_parent_of_level_x(tElementOfBasis(i),aDim,aNumElements,tWhichLevel);
            tElementOfaBasis = tElementOfBasis(i);
            if( tParent == UINT_MAX)
                tParent = tElementOfBasis(i);
            break;
        }
    }
    uint tNumberElements = mBaseElement.give_number_of_elements(aLevel,aDim,aNumElements);
    BoostBitset Element_dummy(tNumberElements+1); // Dummy bitset to remember the tested elements in the filter radius
    Mat<uint> tPossibleElementList(aElementActiveDesign.count()+1,1,0);
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
        Mat<uint> tElement_neighbour = mBaseElement.give_neighbor_of_element(tPossibleElementList(tl),aDim,tBuffer,aNumElements);
        tPossibleElementListLevel = mBaseElement.give_element_level(tPossibleElementList(tl),aDim,aNumElements);
        for( uint i = 0; i < tElement_neighbour.length(); i++)
        {
            tElementLevel = mBaseElement.give_element_level(tElement_neighbour(i),aDim,aNumElements);
            if( tElement_neighbour(i) > 0 && tElement_neighbour(i) < Element_dummy.size() && Element_dummy.test(tElement_neighbour(i)) == 0 && tPossibleElementListLevel == tElementLevel)
            {
                //                std::cout << " tElement_neighbour(i) " << tElement_neighbour(i) << std::endl;
                //                tElementPosition = this->give_position_of_element(tElement_neighbour(i),aDim,aNumElements);

                //                if( tElementPosition(0) >= aElementList.Polynomial*pow(2,tElementLevel) && tElementPosition(0) < (aNumElements(0)-aElementList.Polynomial)*pow(2,tElementLevel)
                //                        && tElementPosition(1) >= aElementList.Polynomial*pow(2,tElementLevel) && tElementPosition(1) < (aNumElements(1)-aElementList.Polynomial)*pow(2,tElementLevel) )
                //                {
                tElementMiddleCoordinate = mHMRElement.give_middlecoordinate_from_element(tElement_neighbour(i),aDim,aNumElements,aDimensions,aDimensions_Offset);
                //                tElementMiddleCoordinate.print("tElementMiddleCoordinate");
                //                tBasisCoordinate.print("tBasisCoordinate");
                //                tElement_neighbour.print("tElement_neighbour");
                tElementMiddleCoordinate = tBasisCoordinate - tElementMiddleCoordinate;
                tDistance = norm(tElementMiddleCoordinate);
                if( (aFilterRadius - tDistance) > 0.0 )
                {
                    tPossibleElementList(tVar) = tElement_neighbour(i);
                    Element_dummy.set(tElement_neighbour(i));
                    tVar++;
                }
                //                }
            }
        }
        if( tPossibleElementListLevel < aLevel && aElementActiveDesign.test(tPossibleElementList(tl)) == 0)
        {
            tChildren = mBaseElement.give_children_of_element(tPossibleElementList(tl),aDim,aNumElements);
            for( uint i = 0; i < tChildren.length(); i++)
            {
                tElementMiddleCoordinate = mHMRElement.give_middlecoordinate_from_element(tChildren(i),aDim,aNumElements,aDimensions,aDimensions_Offset);
                tElementMiddleCoordinate = tBasisCoordinate - tElementMiddleCoordinate;
                tDistance = norm(tElementMiddleCoordinate);
                if( (aFilterRadius - tDistance) > 0.0 )
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
    if( aElementActiveDesign.size() < (tPossibleElementList.max()+1)  )
        aElementActiveDesign.resize(tPossibleElementList.max()+1);
    Mat<uint> tPossibleBasis(tVar*pow(aPolynomialDesign+1,aDim),1,0);
    Mat<uint> tBasis;
    uint tBasisLevel;
    Mat<real> tBasisCoordinates;
    tVar = 0;
    //Save all basis design functions of each element
    for(uint i = 0; i < tPossibleElementList.length(); i++)
    {
        if( aElementActiveDesign.test(tPossibleElementList(i)) == 1 )
        {
            tBasis = mHMRElement.give_basis_of_element(tPossibleElementList(i),aDim,aPolynomialDesign,aNumElements);
            tPossibleBasis.rows(tVar*tBasis.length(),(tVar+1)*tBasis.length()-1) = tBasis.rows(0,tBasis.length()-1);
            tVar++;
        }
    }
    tPossibleBasis.resize(tVar*tBasis.length(),1);
    tPossibleBasis = unique(tPossibleBasis);
    if( aDesignBSplineActive.size() < (tPossibleBasis.max()+1)  )
        aDesignBSplineActive.resize(tPossibleBasis.max()+1);
    Mat<real> tBasisFilterList(tPossibleBasis.length(),3,0);
    tVar = 0;
    //Search for active basis design functions and check the radius
    for( uint i = 0; i < tPossibleBasis.length(); i++)
    {
        if( aDesignBSplineActive.test(tPossibleBasis(i)) == 1  )
        {
            //            std::cout << " tPossibleBasis(i) " << tPossibleBasis(i) << std::endl;
            tBasisLevel = mBasis.give_basis_level(tPossibleBasis(i),aDim,aPolynomialDesign,aNumElements);
            tBasisCoordinates = mBasis.give_coordinate_from_basis(tPossibleBasis(i),aDim,aPolynomialDesign,aNumElements,aDimensions,aDimensions_Offset);
            //            tBasisCoordinates.print("tBasisCoordinates");
            tBasisCoordinates = tBasisCoordinate - tBasisCoordinates;
            tDistance = norm(tBasisCoordinates);
            if( (aFilterRadius - tDistance) > 0.0 )
            {
                tBasisFilterList(tVar,0) = tPossibleBasis(i);
                tBasisFilterList(tVar,1) = aFilterRadius - tDistance;
                tBasisFilterList(tVar,2) = ((real)taBasisLevel+1.0)/((real)tBasisLevel+1.0);
                tVar++;
            }
        }
    }
    tBasisFilterList.resize(tVar,3);
    return tBasisFilterList;
}

Mat<real>
Hierarchical_Mesh_Filter::Filter_for_smoothing_new(
        uint const & aBasis,
        uint const & aDim,
        uint const & aPolynomialDesign,
        Mat<uint> const & aNumElements,
        uint const & aLevel,
        real const & aFilterRadius,
        BoostBitset & aDesignBSplineActive,
        Mat<real> const & aDimensions,
        Mat<real> const & aDimensions_Offset,
        Mat<real> const & aDimensionsOriginal,
        Mat<real> const & aPointOfOrigin)
{
    Mat<uint> tBasisPosition = mBasis.give_position_of_basis(aBasis,aDim,aPolynomialDesign,aNumElements);
    Mat<real> tBasisCoordinate = mBasis.give_coordinate_from_basis(aBasis,aDim,aPolynomialDesign,aNumElements,aDimensions,aDimensions_Offset);
    uint taBasisLevel = mBasis.give_basis_level(aBasis,aDim,aPolynomialDesign,aNumElements);
    Mat<real> tCoord0 = tBasisCoordinate;
    for(uint i = 0; i < tCoord0.length(); i++)
        tCoord0(i) -= aFilterRadius;

    if( tCoord0(0) < aPointOfOrigin(0) )
        tCoord0(0) = aPointOfOrigin(0)+0.000000000001;
    if( tCoord0(1) < aPointOfOrigin(1) )
        tCoord0(1) = aPointOfOrigin(1)+0.000000000001;
    if( aDim == 3 &&  tCoord0(2) < aPointOfOrigin(2) )
        tCoord0(2) = aPointOfOrigin(2)+0.000000000001;
    if( tCoord0(0) > (aPointOfOrigin(0)+aDimensionsOriginal(0)) )
        tCoord0(0) = aPointOfOrigin(0)+aDimensionsOriginal(0)-0.000000000001;
    if( tCoord0(1) > (aPointOfOrigin(1)+aDimensionsOriginal(1)) )
        tCoord0(1) = aPointOfOrigin(1)+aDimensionsOriginal(1)-0.000000000001;
    if( aDim == 3 &&  tCoord0(2) > (aPointOfOrigin(2)+aDimensionsOriginal(2)) )
        tCoord0(2) = aPointOfOrigin(2)+aDimensionsOriginal(2)-0.000000000001;

    Mat<real> tCoord1 = tBasisCoordinate;
    for(uint i = 0; i < tCoord0.length(); i++)
        tCoord1(i) += aFilterRadius;

    if( tCoord1(0) < aPointOfOrigin(0) )
        tCoord1(0) = aPointOfOrigin(0)+0.000000000001;
    if( tCoord1(1) < aPointOfOrigin(1) )
        tCoord1(1) = aPointOfOrigin(1)+0.000000000001;
    if( aDim == 3 &&  tCoord1(2) < aPointOfOrigin(2) )
        tCoord1(2) = aPointOfOrigin(2)+0.000000000001;
    if( tCoord1(0) > (aPointOfOrigin(0)+aDimensionsOriginal(0)) )
        tCoord1(0) = aPointOfOrigin(0)+aDimensionsOriginal(0)-0.000000000001;
    if( tCoord1(1) > (aPointOfOrigin(1)+aDimensionsOriginal(1)) )
        tCoord1(1) = aPointOfOrigin(1)+aDimensionsOriginal(1)-0.000000000001;
    if( aDim == 3 &&  tCoord1(2) > (aPointOfOrigin(2)+aDimensionsOriginal(2)) )
        tCoord1(2) = aPointOfOrigin(2)+aDimensionsOriginal(2)-0.000000000001;

    uint tElement0 = mHMRElement.give_element_for_coordinate_on_level_zero(aDim,aNumElements,aDimensions,aDimensions_Offset,aPointOfOrigin,tCoord0);
    uint tElement1 = mHMRElement.give_element_for_coordinate_on_level_zero(aDim,aNumElements,aDimensions,aDimensions_Offset,aPointOfOrigin,tCoord1);
    Mat<uint> tElementPos0 = mBaseElement.give_position_of_element(tElement0,aDim,aNumElements);
    Mat<uint> tElementPos1 = mBaseElement.give_position_of_element(tElement1,aDim,aNumElements);
    uint tVar = 0;
    Mat<uint> tBasisOfElement;
    Mat<uint> tListOfBasis;
    Mat<uint> tListOfPossibleElements;
    Mat<uint> tIJKPosition(aDim,1,0);
    if( aDim == 2)
    {
        tListOfPossibleElements.set_size((tElementPos1(0)+1-tElementPos0(0))*(tElementPos1(1)+1-tElementPos0(1)),1,UINT_MAX);
        for(uint i = tElementPos0(0); i <= tElementPos1(0); i++)
        {
            for(uint j = tElementPos0(1); j <= tElementPos1(1); j++)
            {
                tIJKPosition(0) = i; tIJKPosition(1) = j;
                tListOfPossibleElements(tVar) = mBaseElement.give_element_of_position(0,aDim,aNumElements,tIJKPosition);
                tVar++;
            }
        }
    }
    else if( aDim == 3)
    {
        tListOfPossibleElements.set_size((tElementPos1(0)+1-tElementPos0(0))*(tElementPos1(1)+1-tElementPos0(1))*(tElementPos1(2)+1-tElementPos0(2)),1,UINT_MAX);
        for(uint i = tElementPos0(0); i <= tElementPos1(0); i++)
        {
            for(uint j = tElementPos0(1); j <= tElementPos1(1); j++)
            {
                for(uint k = tElementPos0(2); k <= tElementPos1(2); k++)
                {
                    tIJKPosition(0) = i; tIJKPosition(1) = j; tIJKPosition(2) = k;
                    tListOfPossibleElements(tVar) = mBaseElement.give_element_of_position(0,aDim,aNumElements,tIJKPosition);
                    tVar++;
                }
            }
        }
    }
    tVar = tListOfPossibleElements.length();
    uint tVarc = tVar;
    uint tMaxNumElemLevel = 0;
    for(uint level = 0; level<aLevel+1; level ++)
    {
        tMaxNumElemLevel += pow(pow(2,aDim),level);
    }
    tListOfPossibleElements.resize(tListOfPossibleElements.length()*tMaxNumElemLevel,1);
    tMaxNumElemLevel -= pow(pow(2,aDim),aLevel);
    Mat<uint> tChildren;
    for(uint i = 0; i<tVarc*tMaxNumElemLevel; i++)
    {
        tChildren = mBaseElement.give_children_of_element(tListOfPossibleElements(i),aDim,aNumElements);
        tListOfPossibleElements.rows(tVar+i*pow(2,aDim),tVar+(i+1)*pow(2,aDim)-1) = tChildren.rows(0,tChildren.length()-1);
    }
    tListOfBasis.set_size(tListOfPossibleElements.length()*pow(aPolynomialDesign+1,aDim),1,UINT_MAX);
    for( uint i = 0; i < tListOfPossibleElements.length(); i++)
    {
        tBasisOfElement = mHMRElement.give_basis_of_element(tListOfPossibleElements(i),aDim,aPolynomialDesign,aNumElements);
        tListOfBasis.rows(i*tBasisOfElement.length(),(i+1)*(tBasisOfElement.length())-1) = tBasisOfElement.rows(0,tBasisOfElement.length()-1);
    }
    tListOfBasis = unique(tListOfBasis);
    Mat<real> tBasisFilterList(tListOfBasis.length(),3,0);
    tVar = 0;
    real tDistance;
    uint tBasisLevel;
    Mat<real> tBasisCoordinates;
    //Search for active basis design functions and check the radius
    MORIS_ASSERT( aDesignBSplineActive.size() > tListOfBasis.max(), "Wrong basis function is found, mayebe outside the domain ?");
    for( uint i = 0; i < tListOfBasis.length(); i++)
    {
        if( aDesignBSplineActive.test(tListOfBasis(i)) == 1  )
        {
            //            std::cout << " tPossibleBasis(i) " << tPossibleBasis(i) << std::endl;
            tBasisLevel = mBasis.give_basis_level(tListOfBasis(i),aDim,aPolynomialDesign,aNumElements);
            tBasisCoordinates = mBasis.give_coordinate_from_basis(tListOfBasis(i),aDim,aPolynomialDesign,aNumElements,aDimensions,aDimensions_Offset);
            //            tBasisCoordinates.print("tBasisCoordinates");
            tBasisCoordinates = tBasisCoordinate - tBasisCoordinates;
            tDistance = norm(tBasisCoordinates);
            if( (aFilterRadius - tDistance) > 0.0 )
            {
                tBasisFilterList(tVar,0) = tListOfBasis(i);
                tBasisFilterList(tVar,1) = aFilterRadius - tDistance;
                tBasisFilterList(tVar,2) = ((real)taBasisLevel+1.0)/((real)tBasisLevel+1.0);
                tVar++;
            }
        }
    }
    tBasisFilterList.resize(tVar,3);
    return tBasisFilterList;
}

/*
 * cl_Hierarchical_Mesh.hpp
 *
 *  Created on: May 10, 2017
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_HPP_

// MORIS library header files.
#include <fstream>
#include "algorithms.hpp"
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_Cell.hpp" // CON/src
#include "fn_unique.hpp" // LNA/src
#include "fn_find_unique.hpp" // LNA/src
#include "fn_find.hpp" // LNA/src
#include "fn_sum.hpp" // LNA/src
#include "fn_reshape.hpp" // LNA/src
#include "fn_histc.hpp" // LNA/src
#include "fn_isempty.hpp" // LNA/src
#include "fn_dot.hpp" // LNA/src
#include "cl_Bitset.hpp" // CON/src
#include "chronos.hpp"
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_Bspline.hpp" // MOD/src
#include "cl_BoostBitset.hpp" // CON/src
#include <boost/dynamic_bitset.hpp>

namespace moris
{

    // Argument list for for the elements of the B-Splines
    struct ElementList_struc
    {
        Mat< real > ControlPoints;
        Mat< uint > SideSet;
        Mat<uint> BlockSet;
        Mat<uint> BlockSetOutput;
        uint Dim = 2;
        uint Polynomial = 1;
        uint PolynomialDesign = 1;
        uint MaxPolynomial = 1;
        uint Level = 0;
        uint LevelLastStep = 0;
        uint LevelDesignLastStep = 0;
        uint NumberElements = 0;
        uint RefinementLevel=1;
        uint BufferElements = 1;
        uint FeatureResolution = 3;
        uint Refinement = 0;
        uint TrueNodalLocaltoGlobalSize;
        uint MaxDesignLevelOfRefinement = 2;
        uint MaxLevelSetLevelOfRefinement = 2;
        uint MaxSurfaceRefinement = 0;
        uint MaxVolumeRefinement = 0;
        uint InitialRefinement = 0;
        real FilterRadius = 0.4;
        real DensityLowerBound = 0.5;
        real DensityIncrementBound = 0.0;
        real DensityUpperBound = 1.0;
        real FeatureLowerBound = 0.5;
        uint DesignVariablesOn = 1;
        bool SDFField = false;
        bool PerformNormalization = false;
        bool LvlSetMethod = false;
        bool FilterLevelWeight = false;
        bool AdaptiveBufferLayer = false;
        bool Staircasebuffer = false;
        bool TruncatedBsplines = false;
        uint MapDeactivateElementsToCoarsest = 0;
        uint SparseSize = 0;
        Mat<real> DeleteDomain;
        Mat<uint> ElementListOnProc;
        Mat<uint> ElementListOnProcInit;
        Mat<uint> ElementListOnProcAura;
        Mat<uint> ElementListActiveDesign;
        Mat<real> Dimensions;
        Mat<real> Dimensions_Offset;
        Mat<uint> NumElements;
        BoostBitset ElementActive;
        BoostBitset ElementActiveDummy;
        BoostBitset ElementActiveDesign;
        BoostBitset ElementActiveLastStep;
        BoostBitset ElementActiveFEMLastStep;
        Mat<uint> DeactivateElement;
        Mat<uint> DeactivateElementRefLvl;
        Mat<uint> DeactivateElementInit;
        Mat<uint> DeactivateElementList;
        Mat<real> Density;
        Mat<real> ElementField;
        Mat<uint> Fetopo;
        Mat<uint> ElemLocaltoGlobal;
        Mat<uint> ElemLocaltoGlobalLastStep;
        Mat<uint> ElemLocaltoGlobalLastStepFEM;
        Mat<uint> NodalLocaltoGlobal;
        Mat<uint> NodalLocaltoGlobalMap;
        Mat<uint> NodalLocaltoGlobalExist;
        Mat<uint> Decomp;
        Mat<uint> DecompAura;
        Mat<uint> ProcNeighbour;
        Mat<uint> ProcNeighbourSendRecv;
        Mat<uint> MPICoord;
        Mat<uint> NumSets; //First row are number of node sets, second row for number of side sets
        Cell<std::string > BSetNames;
        Cell< Mat< uint > > ElemIdsAndSideOrds;
        Cell< std::string > SSetNames;
        Cell< Mat< uint > > NEntIds;
        Cell< std::string > NSetNames;
        Mat< real > TMatrix;
        Mat< uint > IdField;
        Mat<real> TMatrixField;
        Mat<uint> IdFieldField;
        Mat< real > TMatrixDesign;
        Mat< uint > IdFieldDesign;
        Mat<real> TMatrixFieldDesign;
        Mat<uint> IdFieldFieldDesign;
        Mat<uint> ListOfChildren;
        Cell<Mat<uint>> IdFieldOfChildren;
        Cell<Mat<real>> TMatrixOfChildren;
        Cell<Mat<real>> tTMatrixOfChildren;
        Mat<real> Coeff;
        Mat<uint> ConsecutiveNumbering;
        real Elementwidth = 0.033333;
        real sLSthresh = 0.5;
        real Prjbeta = 0.01;
        real DensityShift = 0.0;
        real SimpExponent = 3.0;
        real sPrjbeta = 0.001;
        uint MaxLevelOfDesignboxRefinement =2;
        uint MaxLevelOfCargoboxRefinement =2;
        uint MaxLevelOfBoltsRefinement =2;
        Mat<uint> MaxLevelOfObjectRefinement= {{2},{2},{2}};
    };

    // Argument list for the basis functions of the B-Splines
    struct BasisList_struc
    {
        BoostBitset BasisActive;
        BoostBitset DesignBSplineActive;
        BoostBitset DesignBSplineActiveLastStep;
        BoostBitset BasisActivetLastStep;
        BoostBitset BasisActiveDummy;
        Mat<uint> BasisActiveDesignList;
        Mat<uint> BasisActiveDesignListLastStep;
        Mat<uint> BasisActiveListLastStep;
        Mat<uint> BasisActiveList;
        Mat<uint> DesignBSplineListMap;
        Mat<uint> NodeSet;
        Mat<real> BasisFilterList;
        uint NumberBasis;
        Mat<real> NodalField;
        Mat<real> NodalLvlSet;
        Mat<real> NodalADV;
        Mat<real> NodalADVNewField;
        Mat<real> NodalPDV;
        Mat<real> NodalLvLSetField;
        Mat<real> NodalSDFField;
        Mat<uint> NodeProcs;
        Mat<uint> Decomp;
        Mat<real> PointOfOrigin = {{0.0},{0.0},{0.0}};
    };

    // Argument list for for the elements of the B-Splines
    struct TMatrix_and_IdField_struc
    {
        Mat< real > TMatrix;
        Mat< uint > IdField;
    };

    class Hierarchical_Mesh
    {
    protected:

    public:
        ElementList_struc mElementList; // Element function specific argument list
        BasisList_struc   mBasisList;   // Basis function specific argument list
        TMatrix_and_IdField_struc   mTMatrix_and_IdField;   // TMatrix and IdField specific argument list
        /**
         * Hierarchical_Mesh constructor
         */
        Hierarchical_Mesh()
        {
        }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Hierarchical_Mesh() = default;
        static Mat<uint>
        give_vector_entries(Mat<uint> & aMat, Mat<uint> & bMat);

        static Mat<real>
        give_vector_entries(Mat<real> & aMat, Mat<uint> & bMat);

        /**
         * Creates the T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,....
         *
         * @param[in] aElementNumber      Element number.
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] Tmatrix_and_Id     Data is stored in aElementList
         *
         */
        void give_Tmatrix_and_id_field(uint & aElementNumber,
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);
        void give_Truncated_Tmatrix_and_id_field(uint & aElementNumber,
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Creates the T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... for the design variables
         *
         * @param[in] aElementNumber      Element number.
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] Tmatrix_and_Id     Data is stored in aElementList
         *
         */
        void give_Tmatrix_and_id_field_for_designvariables(uint & aElementNumber,
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);
        void give_Truncated_Tmatrix_and_id_field_for_designvariables(uint & aElementNumber,
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Give the projection matrix for the design variables to the FEM mesh
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         *
         * @param[out] tProject     Projection matrix
         *
         */
        static Mat<real>
        give_projection_matrix(
                ElementList_struc & aElementList);

        /**
         * Creates the T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... for the L2-projection of the design variables from the old to the new mesh
         *
         * @param[in] aElementNumber      Element number.
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] Tmatrix_and_Id     Data is stored in aElementList
         *
         */
        void give_Tmatrix_and_id_field_for_designvariables_projection(uint & aElementNumber,
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);
        void give_Truncated_Tmatrix_and_id_field_for_designvariables_projection(uint & aElementNumber,
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Creates the T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... for the projection of the design variables from the old to the new mesh. It is needed for the projection from a fine to a coarse
         *
         * @param[in] aElementNumber      Element number.
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] Tmatrix_and_Id     Data is stored in aElementList
         *
         */
        void give_Tmatrix_and_id_field_for_designvariables_projection_coarsening(uint & aElementNumber,
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);
        void give_Truncated_Tmatrix_and_id_field_for_designvariables_projection_coarsening(uint & aElementNumber,
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Creates the T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... for the design variables without projection to fem mesh
         *
         * @param[in] aElementNumber      Element number.
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] Tmatrix_and_Id     Data is stored in aElementList
         *
         */
        void give_Tmatrix_and_id_field_design(uint & aElementNumber,
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);
        void give_Truncated_Tmatrix_and_id_field_design(uint & aElementNumber,
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);
        /**
         * Creates the T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,....
         *
         * @param[in] aElementNumber      Element number.
         * @param[in] aPolynomial         Polynomial degree of the b-spline.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Tmatrix            The T-matrix of the respective element.
         *
         */
        static Mat<real>
        give_Tmatrix(uint & aElementNumber,
                uint & aPolynomial,
                uint & aDim,
                Mat<uint> & aNumElements);

        /**
         * Provides the element number with the position i,j,k
         *
         * @param[in] aLevel              Level of the element.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aPolynomialDegree   Polynomial degree of the b-spline.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] Bspline            Element number.
         *
         */
        static uint
        give_element_of_position(uint & aLevel,
                uint & aDim,
                Mat<uint> & aNumElements,
                Mat<uint> & aIJKPosition);

        /**
         * Provides the position of an element from a tensorial grid
         *
         * @param[in] aElement            Element number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Element_position        Position I,J,K.
         *
         */
        static Mat<uint>
        give_position_of_element(uint & aElement,
                uint & aDim,
                Mat<uint> & aNumElements);

        /**
         * Provides the level of the element from the hierarchical mesh.
         *
         * @param[in] aElement            Element number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] level        Level of the element from the hierarchical mesh.
         *
         */
        static uint
        give_element_level(uint & aElement,
                uint & aDim,
                Mat<uint> & aNumElements);

        /**
         * Provides the parent of an element
         *
         * @param[in] aElement            Element number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Parent        Parent of an element.
         *
         */
        static uint
        give_parent_of_element(uint & aElement,
                uint & aDim,
                Mat<uint> & aNumElements);

        /**
         * Provides the parent of an element for a specific level
         *
         * @param[in] aElement            Element number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Parent        Parent of an element.
         *
         */
        static uint
        give_parent_of_level_x(uint & aElement,
                uint & aDim,
                Mat<uint> & aNumElements,
                uint & aWhichLevel);
        /**
         * Provides the number of elements within all levels from level 0 until aLevel
         *
         * @param[in] aLevel              Level of the element.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] element_number        Number of elements within all Levels until aLevel.
         *
         */
        static uint
        give_number_of_elements(uint & aLevel,
                uint & aDim,
                Mat<uint> & aNumElements);

        /**
         * Provides the children of an element
         *
         * @param[in] aElement            Element number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Children        Childrens of an element (4 children for Dim=2 and 8 children for Dim=3).
         *
         */
        static Mat<uint>
        give_children_of_element(uint & aElement,
                uint & aDim,
                Mat<uint> & aNumElements);

        /**
         * Provides the neighbours of an element
         *
         * @param[in] aElement            Element number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aBuffer             Provides the number of layers around the element
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Element_neighbour        Element plus neighbours (9 for aDim = 2, 27 for aDim = 3).\n
         *                                2D:...----------------------------------...3D:....-----plain 1----------------------.....-----plain 2-------------------------\n.....-----plain 3-------------------------\n
         *                                ......|Element(6) Element(7) Element(8)|..........|Element(6) Element(7) Element(8)|.....|Element(15) Element(16) Element(17)|\n.....|Element(24) Element(25) Element(26)|\n
         *                                ......|Element(3) Element(4) Element(5)|..........|Element(3) Element(4) Element(5)|.....|Element(12) Element(13) Element(14)|\n.....|Element(21) Element(22) Element(23)|\n
         *                                ......|Element(0) Element(1) Element(2)|..........|Element(0) Element(1) Element(2)|.....|Element(9) Element(10) Element(11) |\n.....|Element(18) Element(19) Element(20)|\n
         *                                ......|--------------------------------|..........|--------------------------------|.....|-----------------------------------|  .....|-----------------------------------|\n
         *                                2D: Element(4) = aElement, 3D: Element(13) = aElement;
         *
         */
        static Mat<uint>
        give_neighbour_of_element(uint & aElement,
                uint & aDim,
                uint & aBuffer,
                Mat<uint> & aNumElements);

        static Mat<uint>
        give_active_neighbour_of_element(uint & aElement,
                ElementList_struc & aElementList);

        /**
         * Provides the neighbours of an basis
         *
         * @param[in] aBasis              Basis number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aPolynomial         Polynomial degree
         * @param[in] aBuffer             Provides the number of layers around the element
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Basis neighbour        Basis plus neighbours (9 for aDim = 2, 27 for aDim = 3).\n
         *                                2D:...----------------------------------...3D:....-----plain 1----------------------.....-----plain 2-------------------------\n.....-----plain 3-------------------------\n
         *                                ......|Element(6) Element(7) Element(8)|..........|Element(6) Element(7) Element(8)|.....|Element(15) Element(16) Element(17)|\n.....|Element(24) Element(25) Element(26)|\n
         *                                ......|Element(3) Element(4) Element(5)|..........|Element(3) Element(4) Element(5)|.....|Element(12) Element(13) Element(14)|\n.....|Element(21) Element(22) Element(23)|\n
         *                                ......|Element(0) Element(1) Element(2)|..........|Element(0) Element(1) Element(2)|.....|Element(9) Element(10) Element(11) |\n.....|Element(18) Element(19) Element(20)|\n
         *                                ......|--------------------------------|..........|--------------------------------|.....|-----------------------------------|  .....|-----------------------------------|\n
         *                                2D: Element(4) = aElement, 3D: Element(13) = aElement;
         *
         */
        static Mat<uint>
        give_neighbour_of_basis(uint & aBasis,
                uint & aDim,
                uint & aPolynomial,
                uint & aBuffer,
                Mat<uint> & aNumElements);

        /**
         * Provides the neighbours of an element for a flood fill algorithm. This function uses the function "give_neighbour_of_element" and extracts only the nessecary elements
         *
         * @param[in] aElement            Element number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aBuffer             Provides the number of layers around the element
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Element_neighbour        Element plus neighbours (4 for aDim = 2, 6 for aDim = 3).
         *
         */
        static Mat<uint>
        give_neighbour_of_element_for_floodfill(uint & aElement,
                uint & aDim,
                uint & aBuffer,
                Mat<uint> & aNumElements);

        /**
         * Provides a stencil of neighbours of an element
         *
         * @param[in] aFeatureResolution    Size of a feature resolution
         * @param[in] aElement            Element number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         * * @param[in] ElementActive       Bitset with active elements
         *
         * @param[out] Element_stencil     Gives a matrix with neighbour elements. In each row is one specific direction ( 2D: 2 axes + 2 diagonals, 3D: 3 axes + 9 diagonals)
         *
         */
        static Mat<uint>
        give_neighbour_stencil_of_element(uint & aFeatureResolution,
                uint & aElement,
                uint & aDim,
                uint & aPolynomial,
                Mat<uint> & aNumElements,
                BoostBitset & ElementActive);

        /**
         * Provides the level of the basis function from the hierarchical mesh.
         *
         * @param[in] aBasis            Basis function number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] level        Level of the element from the hierarchical mesh.
         *
         */
        static uint
        give_basis_level(uint & aBasis,
                uint & aDim,
                uint & aPolynomial,
                Mat<uint> & aNumElements);

        /**
         * Provides the number of basis functions within all levels from level 0 until aLevel
         *
         * @param[in] aLevel              Level of the basis functions.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] basis_number        Number of basis functions within all Levels until aLevel.
         *
         */
        static uint
        give_number_of_basis(uint & aLevel,
                uint & aDim,
                uint & aPolynomial,
                Mat<uint> & aNumElements);

        /**
         * Provides the position of a basis function from a tensorial grid
         *
         * @param[in] aBasis            Basis function number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Basis_position        Position I,J,K.
         *
         */
        static Mat<uint>
        give_position_of_basis(uint & aElement,
                uint & aDim,
                uint & aPolynomial,
                Mat<uint> & aNumElements);

        /**
         * Provides the basis function number with the position i,j,k
         *
         * @param[in] aBasis            Basis function number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] basis_number        Number of basis function.
         *
         */
        static uint
        give_basis_of_position(uint & aLevel,
                uint & aDim,
                uint & aPolynomial,
                Mat<uint> & aNumElements,
                Mat<uint> & aIJKPosition);

        /**
         * Provides the basis function of the parent (Works only for a linear polynomial degree)
         *
         * @param[in] aBasis            Basis function number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] basis_number        Number of basis function of the parent.
         *
         */
        static uint
        give_basis_of_parent(uint & aBasis,
                uint & aDim,
                uint & aPolynomial,
                Mat<uint> & aNumElements);

        /**
         * Provides the basis function of the child
         *
         * @param[in] aBasis            Basis function number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] basis_number        Number of basis function of the parent.
         *
         */
        static uint
        give_basis_of_child(uint & aBasis,
                uint & aDim,
                uint & aPolynomial,
                Mat<uint> & aNumElements);

        /**
         * Provides the basis functions of an element from a tensorial grid
         *
         * @param[in] aLevel              Level of the basis functions.
         * @param[in] aElement            Element number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] basis        Basis functions of an element.
         *
         */
        static Mat<uint>
        give_basis_of_element(uint & aElement,
                uint & aDim,
                uint & aPolynomial,
                Mat<uint> & aNumElements);

        /**
         * Provides the elements, which have support with the basis function from a tensorial grid
         *
         * @param[in] aBasis            Basis number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Elements        Elements with support of basis function aBasis.
         *
         */
        static Mat<uint>
        give_element_of_basis(uint & aBasis,
                uint & aDim,
                uint & aPolynomial,
                Mat<uint> & aNumElements);

        /**
         * Provides a list of active elements
         *
         * @param[in] aElementList        Struc with all information about the elements.
         *
         * @param[out] aElementList       Updates the vector "ListActiveElements".
         *
         */
        void give_active_elements(
                ElementList_struc & aElementList);

        /**
         * Provides a list of active elements with an aura (1 element layer around the cpu)
         *
         * @param[in] aElementList        Struc with all information about the elements.
         *
         * @param[out] aElementList       Updates the vector "ListActiveElements".
         *
         */
        void give_active_elements_with_aura(
                ElementList_struc & aElementList);

        /**
         * Provides the coordinates for a specific basis function
         *
         * @param[in] aBasis        A specific basis function.
         * @param[in] aBasis        A specific basis function.
         *@param[in] aElementList        Struc with all information about the elements.
         *
         * @param[out] coordinates       Coordinates in x,y,z direction.
         *
         */
        static Mat<real>
        give_coordinate_from_basis(uint & aBasis,
                uint & aPolynomial,
                ElementList_struc & aElementList);
        /**
         * Provides the coordinate in the middle of an element
         *
         * @param[in] aElement        A specific element.
         *@param[in] aElementList        Struc with all information about the elements.
         *@param[in]
         * @param[out] coordinates       Coordinates in x,y,z direction.
         *
         */
        static Mat<real>
        give_middlecoordinate_from_element(
                uint & aElement,
                ElementList_struc & aElementList);
        /**
         * Provides a list of deactivated elements
         *
         * @param[in] aDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aElementActive      Sparsity matrix with list of active/deactive elements.
         * @param[in] aBasisActive        Sparsity matrix with list of active/deactive basis functions.
         * @param[in] aDeactivateElement        Vector of all elements, which need to be deactivated.
         *
         * @param[out] aDeactivateElement       Updated list of deactivated elements (linear dependencies, ...) .
         *
         */
        void give_deactivated_elements(
                ElementList_struc & aElementList);

        /**
         * Activates or deactivates elements, inheret the side sets, block sets and calculates the coordiantes of the children
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         * @param[in] aDeactivateElement        Vector of all elements, which need to be deactivated.
         *
         * @param[out] aElementList       Updates the elements: active/deactive, block sets, side sets, coordinates .
         *
         */
        void hierarchical_element_refinement(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Activates or deactivates basis functions, inheret the node sets
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] aBasisList       Updates the basis functions: active/deactive, node sets
         *
         */
        void hierarchical_basisfunction_refinement(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Activates or deactivates basis functions for the design variables, inheret the node sets
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] aBasisList       Updates the basis functions: active/deactive, node sets
         *
         */
        void hierarchical_basisfunction_designvariables_refinement(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);
        /**
         * Deactivate/Activate elements on finer leves and activate basis functions until a certain level of refinement
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] aBasisList       Updates the basis functions: active/deactive, node sets
         * @param[in] aElementList      Updates the elements: (see top: ElementList_struc)
         *
         */
        void hierarhical_mesh_refinement(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Deactivate/Activate elements on finer leves and activate design basis functions until a certain level of refinement
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] aBasisList       Updates the design basis functions: active/deactive, node sets
         * @param[in] aElementList      Updates the elements: (see top: ElementList_struc)
         *
         */
        void hierarhical_design_mesh_refinement(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Create all the data, which is necessarry for the mesh class
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] aBasisList       Updates the basis functions: active/deactive, node sets
         * @param[in] aElementList      Updates the elements: (see top: ElementList_struc)
         *
         */
        void create_data_for_mesh(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Create the dof connectivity list for femdoc
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] file               Creates a file with the dof-connectivity information
         *
         */
        void create_dofconnectivity_outputfile(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Create the design variables in a file for femdoc
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] file           Creates a file with nodal information (Which design variables interpolates in the fem nodes)
         *
         */
        void create_designvariables_outputfile(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Save the active design elements in a file
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] file           Creates a file and saves all active design elements
         *
         */
        void save_active_design_elements_in_file(
                ElementList_struc & aElementList);

        /**
         * Save the active fem elements in a file
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] file           Creates a file and saves all active fem elements
         *
         */
        void save_active_fem_elements_in_file(
                ElementList_struc & aElementList);

        /**
         * Save the active design basis functions
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] file           Creates a file and saves all active design basis functions
         *
         */
        void save_active_design_basis_in_file(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Save the active  basis functions
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] file           Creates a file and saves all active basis functions
         *
         */
        void save_active_basis_in_file(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        void save_coordlist_in_file(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);
        /**
         * Save a list with the information of the support of the basis functions within an element. The number of rows are the basis functions and the columns have the element ID's
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] file           Creates a file and saves all active design basis functions
         *
         */
        void save_basis_to_element_support(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Updates the SideSet, NodeSet for the active elements/ basis functions
         * (Update of block sets and coordinates are created in 'create_mesh_data')
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] aBasisList       Update of node sets
         * @param[out] aElementList       Update of side set
         *
         */
        void update_side_node_set(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Activates or deactivates basis functions, inheret the node sets
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] aBasisList       Updates the basis functions: active/deactive, node sets
         *
         */
        void give_node_proc_owner(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Provides a list of elements, which needs to be removed (Elements on the layer for the b-splines)
         *
         * @param[in] aDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aElementActive      Sparsity matrix with list of active/deactive elements.
         *
         * @param[out] RemoveElements       A list of elements, which needs to be removed (Elements on the layer for the b-splines).
         *
         */
        static moris::Mat<uint>
        remove_elements(ElementList_struc & aElementList);

        /**
         * Gives the T-matrices of all childs for 1D, 2D or 3D, for the polynomial degrees p = 1, 2, 3
         *
         * @param[in] aPolynomialDegree   Polynomial degree of the b-splines for each direction p_1 = p_2 = p_3.
         * @param[in] aDim                Number of dimensions.
         *
         * @param[out] T_child            Gives the Tmatrices of all childs in a Cell. The order is:\n
         *                                2D:...------------------... 3D:...-----plain 1------.....-----plain 2------\n
         *                                ......|Child 3..Child 4|..........|Child 3..Child 4|.....|Child 7..Child 8|\n
         *                                ......|................|..........|................|.....|................|\n
         *                                ......|Child 1..Child 2|..........|Child 1..Child 2|.....|Child 5..Child 6|\n
         *                                ......|----------------|..........|----------------|.....|----------------|
         *
         */
        static Cell<Mat<real>>
        give_Tmatrices_of_childs(uint & aPolynomialDegree,
                uint & aDim);
        static Cell<Mat<real>>
        give_Tmatrices_of_childs_for_design(uint & aPolynomialDegree,
                uint & aDim);

        /**
         * Creates the mesh data for MTK (Create element connectivity, coordinates and block sets)
         *
         * @param[in] aElementList      Struc with element information (see top: ElementList_struc)
         * @param[in] aBasisList         Struc with basis function information (see top: BasisList_struc)
         *
         * @param[out] aElementList     New list for the Element connectivity, which considers only active basis functions and hanging nodes
         *
         */
        void create_mesh_data(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Creates a vector of deactivated elements (based on a nodal level set field)
         *
         * @param[in] aElementList       ElementList.DeactivateElement
         * @param[in] aBasisList        Information about the basis functions
         *
         * @param[out] aElementList       ElementList.DeactivateElements
         *
         */
        void deactivate_element_list(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * All bitsets are send to proc 0, which compares data and save on tBitsetAlpha only those, where both vectors have at the same position a one. Afterwards this tBitsetAlpha is broadcasted to all procs
         *
         * @param[in] aElementList       Struc which contains the bitset "ElementActive"
         *
         * @param[out] aElementList      Struc which contains the bitset "ElementActive" (New data due to broadcast)
         *
         */
        void sendrecv_deactivation_and_bcast_element_data(
                ElementList_struc & aElementList);

        /**
         * All bitsets are send to proc 0, which compares data and save on tBitsetAlpha all bits with a one. Afterwards this tBitsetAlpha is broadcasted to all procs
         *
         * @param[in] aElementList       Struc which contains the bitset "ElementActive"
         *
         * @param[out] aElementList      Struc which contains the bitset "ElementActive" (New data due to broadcast)
         *
         */
        void sendrecv_refinement_and_bcast_element_data(
                ElementList_struc & aElementList);

        /**
         * Broadcast a Message to all Procs
         *
         * @param[in] aMessage      A bitset with the active elements
         *
         * @param[out] aMessage     The message of Proc 0
         *
         */
//        void broadcast(BoostBitset & aMessage);

        /**
         * All bitsets are send to proc 0, which compares data and save on tBitsetAlpha all bits with a one. Afterwards this tBitsetAlpha is broadcasted to all procs
         *
         * @param[in] aBasisList       Struc which contains the bitset "BasisActive"
         *
         * @param[out] aBasisList      Struc which contains the bitset "BasisActive" (New data due to broadcast)
         *
         */
        void sendrecv_refinement_and_bcast_basis_data(
                BasisList_struc & aBasisList);

        /**
         * All bitsets are send to proc 0, which compares data and save on tBitsetAlpha only those, where both vectors have at the same position a one. Afterwards this tBitsetAlpha is broadcasted to all procs
         *
         * @param[in] aElementList       Struc which contains the bitset "BasisActive"
         *
         * @param[out] aElementList      Struc which contains the bitset "BasisActive" (New data due to broadcast)
         *
         */
        void sendrecv_deactivation_and_bcast_basis_data(
                BasisList_struc & aBasisList);

        /**
         * Gather all max element levels and bcast the max level, that every processor is aware of the highes level. Furthermore, the max number of elements and basis functions is updated
         *
         * @param[in] aElementList       aElementList.Level
         *
         * @param[out] aElementList      aElementList.Level (updated)
         *
         */
        void gather_refinementlevel_and_bcast_data(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Gather a message to proc 0, take the max of the proc values and broadcast the max.
         *
         * @param[in] aElementList       aElementList.Level
         *
         * @param[out] aElementList      aElementList.Level (updated)
         *
         */
        void gather_value_and_bcast_max(
                uint & aMessage);

//        void sendrecv_refinement_of_neighbours_and_bcast_element_data(
//                ElementList_struc & aElementList);

//        void sendrecv_deactivation_of_neighbours_and_bcast_element_data(
//                ElementList_struc & aElementList);

        void sendrecvring_refinement_and_bcast_element_data(
                ElementList_struc & aElementList);

        void sendrecvring_deactivation_and_bcast_element_data(
                ElementList_struc & aElementList);

        // Create the initial mesh and output
        void Initial_mesh(        ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Find Elements in a stencil, which need to be refined
         *
         * @param[in] aElementList       ElementList.Stencilsize
         * @param[in] aBasisList       Data of BasisList
         *
         * @param[out] aElementList       ElementList.DeactivateElements
         *
         */
        void Find_elements_in_stencil(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /* Find Elements, which have support with a basis function above a threshold, which need to be refined. And the second criteria are intersected Elements ( Level set zero contour)
         *
         * @param[in] aElementList     Data of ElementList
         * @param[in] aBasisList       Data of BasisList
         *
         * @param[out] aElementList       ElementList.DeactivateElements
         *
         */
        void find_elements_for_refinement(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);
        /**
         * Creates a List of all objects in all phases for a density field
         *
         * @param[in] aElementList      A vector with a list of elements
         * @param[in] aElementField     A vector list with the element field information: colours (instead of colours, uint's are used from 0 to ... )
         * @param[in] aElementList       Struc which contains the information about the mesh
         *
         * @param[out] ObjectData      Obtains all the objects in all the phases
         *
         */
        static Cell<Cell<Mat<uint>>>
        Floodfill_for_elementfield(
                Mat<uint> & aElementField,
                Mat<uint> & aListOfElement,
                ElementList_struc & aElementList);

        /**
         * Creates a List of all objects in all phases for a level set field
         *
         * @param[in] aElementField      A cell of elements with respect to the phases
         * @param[in] aElementList       Struc which contains the information about the mesh
         * @param[in] aBasisList       Data of BasisList
         *
         * @param[out] ObjectData      Obtains all the objects in all the phases
         *
         */
        static Cell<Cell<Mat<uint>>>
        Floodfill_for_levelset(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Provides the classical rhight hand side, which is needed for a L2 projeciton
         *
         * @param[in] aDim               Number of Dimensions
         * @param[in] aPolynomial        Which polynomial degree
         *
         * @param[out] RHS              RHS vector for a "load" on a unit element
         *
         */
        static Mat<real>
        RHS_for_L2_projection(
                uint & aDim,
                uint & aPolynomial,
                Mat<real> & aNodalField);

        /**
         * Provides the classical mass matrix, which is needed for a L2 projeciton
         *
         * @param[in] aDim               Number of Dimensions
         * @param[in] aPolynomial        Which polynomial degree
         *
         * @param[out] Mass              Mass matrix for a unit element
         *
         */
        static Mat<real>
        MassMatrix_for_L2_projection(
                uint & aDim,
                uint & aPolynomial);

        /**
         * Provides a list of basis functions, which are in a specified circle
         *
         * @param[in] aBasis             Basis function, which is used to search for neighbours
         * @param[in] aElementList       Struc which contains the information about the mesh
         * @param[in] aBasisList       Struc which contains the information about the basis functions
         *
         * @param[out] Mat             Output is a matrix with three columns. First column are a list of elements, second and third column are the distance and the ratio level
         *
         */
        void Filter_for_smoothing(uint & aBasis,
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Provides a list of basis functions, which are in a specified circle
         *
         * @param[in] aElementList       Struc which contains the information about the mesh
         * @param[in] aBasisList       Struc which contains the information about the basis functions
         *
         * @param[out] ID and T-Matrix   Updates with the help of the function "Filter_for_smoothing" the ID's and T-Matrix
         *
         */
        void Update_IDandTMatrix_design(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Creates a L2 projection for the design variables from the old mesh to the refined mesh of the design variables
         *
         * @param[in] aBasisList       Data of BasisList
         * @param[in] aElementList       Struc which contains the information about the mesh
         *
         * @param[out] AbsDersVariables.dat    It creates an output file for femdoc
         *
         */
        void L2_projection(
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        /**
         * Creates a list of deactivated elements, for the initial refinement
         *
         * @param[in] aElementList       Struc which contains the information about the mesh
         *
         * @param[out] DeactivateElements   List of elements, which gets deactivated
         *
         */
        void initial_refinement(
                ElementList_struc & aElementList);

        /**
         * Provides the ownership of an Element
         *
         * @param[in] aElement            Element number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aProcRange        Vector with a range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
         *
         * @param[out] Rank        Proc rank
         *
         */
        static uint
        give_element_owner(uint & aElement,
                uint & aDim,
                Mat<uint> & aNumElements,
                Mat<uint> & aProcRange,
                Mat<uint> & aProcNeighbours);

        /**
         * Provides the a list, which share an element (owner plus aura procs)
         *
         * @param[in] aElement            Element number.
         * @param[in] aDim                Number of dimensions.
         * @param[in] aNumElements        NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aProcRange        Vector with a range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
         *
         * @param[out] Share        Vector with procs, which share the element
         *
         */
        static Mat<uint>
        give_element_share(uint & aElement,
                uint & aDim,
                Mat<uint> & aNumElements,
                Mat<uint> & aProcRange,
                Mat<uint> & aProcNeighbours);

        void
        floodfill_for_STL(
                Mat<real> aNodalField,
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);

        static Mat<uint>
        give_active_neighbor_of_element(
                uint & aElement,
                uint & aDim,
                Mat<uint> & aNumElements,
                uint & aLevel,
                BoostBitset & aElementActive);

        void find_elements_of_object(
                Mat<real> aNodalField,
                Mat<uint> & tRefineElement,
                uint & MaxLevelOfRefinement,
                ElementList_struc & aElementList,
                BasisList_struc & aBasisList);
    };

}
#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_HPP_ */
